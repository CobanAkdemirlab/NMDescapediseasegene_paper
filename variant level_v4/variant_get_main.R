#!/usr/bin/env Rscript
# =============================================================================
# variant_get_main.R
#
# Build the four variant sets used for the variant-level feature comparison:
#
#   INPUT   fs disease gene list   (hgnc_symbol)
#           snv disease gene list  (hgnc_symbol)
#
#   OUTPUT  fs  disease variants   ClinVar P/LP frameshift NMD-escape
#           snv disease variants   ClinVar P/LP stop-gain  NMD-escape
#           fs  control variants   gnomAD frameshift NMD-escape, P/LP removed
#           snv control variants   gnomAD stop-gain  NMD-escape, P/LP removed
#
# Controls are gnomAD variants on the SAME canonical transcripts as the
# disease variants, with anything P/LP in ClinVar removed. That makes the
# disease/control contrast a within-gene contrast, which is what the
# (1 | gene) random intercept in the downstream model is there to exploit.
#
# Merged from: get_fs_variant_new.R, get_snv_variant_new.R,
#              get_snv_control_variant_new.R, get_gnomAD_control.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(biomaRt)
  library(S4Vectors)
})

# ---------------------------------------------------------------------------
# 0. Configuration -- edit paths here, nothing below should need touching
# ---------------------------------------------------------------------------

CFG <- list(

  # --- input gene lists (one column: hgnc_symbol) --------------------------
  fs_gene_file  = "fs_can_AD_acat_FDR0.20_all.txt",
  snv_gene_file = "snv_can_ADrestricted_bh_FDR0.20_all.txt",

  # --- ClinVar aenmd result objects (GRanges .rds) -------------------------
  # fs objects are NOT pre-filtered for the canonical escape rule
  fs_plp_rds    = "fs.rds",
  fs_vus_rds    = "fs_vus.rds",
  fs_benign_rds = "fs_benign.rds",
  # snv objects ARE already escape-filtered ("_can_filtered")
  snv_plp_rds    = "snv_plp_ptc_nmdesc_can_filtered20260201.rds",
  snv_vus_rds    = "snv_vus_ptc_nmdesc_can_filtered20260201.rds",
  snv_benign_rds = "snv_benign_ptc_nmdesc_can_filtered20260201.rds",

  # --- gnomAD PTC table ----------------------------------------------------
  gnomad_file = "ptc_can_NMD_df.csv",

  # --- optional: length-matched control GENE lists -------------------------
  # Leave NULL to skip. Supplying these adds two extra output files of
  # gnomAD variants in control genes (the gene-level control design).
  fs_control_gene_file  = NULL,
  snv_control_gene_file = NULL,

  # --- output names read by variant_compare_main.R -------------------------
  # These four are written verbatim, so the downstream script finds them by the
  # names it reads. Every other set keeps the dated "<stem>_<tag>.csv" form.
  out_names = c(
    snv_disease_variants_clinvar_plp = "snv_variants20260201_plp_dbh_clinvar.csv",
    fs_disease_variants_clinvar_plp  = "fs_variants20260201_plp_acat_clinvar.csv",
    snv_control_variants_gnomad      = "gnomad_snv_filtered_acat_0831.csv",
    fs_control_variants_gnomad       = "gnomad_fs_filtered_bh_0831.csv"
  ),

  # --- run options ---------------------------------------------------------
  out_dir       = "out",
  tag           = format(Sys.Date(), "%Y%m%d"),
  write_vus     = TRUE,   # also emit VUS / benign ClinVar sets
  biomart_cache = "biomart_cache.rds",  # NULL to disable
  ensembl_host  = "https://www.ensembl.org"
)

# Optional repo path helper. Falls back to plain relative paths if absent,
# so the script still runs standalone.
.p <- c("gene level_v3/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../lib/paths.R", "../../../../lib/paths.R")
.p <- .p[file.exists(.p)]
if (length(.p)) source(.p[1]) else {
  data_file <- function(x, ...) x
  out_file  <- function(x, ...) file.path(CFG$out_dir, x)
}
rm(.p)
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Helpers
# ---------------------------------------------------------------------------

msg <- function(...) cat(sprintf(...), "\n", sep = "")

read_gene_list <- function(path) {
  if (is.null(path)) return(NULL)
  df <- readr::read_csv(data_file(path), show_col_types = FALSE)
  col <- if ("hgnc_symbol" %in% names(df)) "hgnc_symbol" else names(df)[1]
  unique(stats::na.omit(df[[col]]))
}

# One biomaRt connection, reused; results memoised to disk so reruns are fast
# and a flaky Ensembl endpoint doesn't cost you the whole run.
.mart <- NULL
get_mart <- function() {
  if (is.null(.mart)) {
    # "genes" is the current biomart name; the configured host is the fallback
    .mart <<- tryCatch(
      biomaRt::useEnsembl("genes", dataset = "hsapiens_gene_ensembl"),
      error = function(e)
        biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
                            host = CFG$ensembl_host))
  }
  .mart
}

.cache <- if (!is.null(CFG$biomart_cache) && file.exists(CFG$biomart_cache)) {
  readRDS(CFG$biomart_cache)
} else list()

cached_bm <- function(attributes, filters, values) {
  values <- sort(unique(as.character(values)))
  if (!length(values)) {
    return(stats::setNames(
      data.frame(matrix(character(0), ncol = length(attributes))), attributes))
  }
  k <- digest_key(attributes, filters, values)
  if (!is.null(.cache[[k]])) return(.cache[[k]])
  res <- biomaRt::getBM(attributes = attributes, filters = filters,
                        values = values, mart = get_mart())
  .cache[[k]] <<- res
  if (!is.null(CFG$biomart_cache)) saveRDS(.cache, CFG$biomart_cache)
  res
}

digest_key <- function(...) {
  paste(vapply(list(...), function(x) paste(x, collapse = ","), character(1)),
        collapse = "||")
}

# Canonical transcript per gene symbol, with UniProt accession attached.
canonical_tx <- function(genes, label = "") {
  bm <- cached_bm(
    attributes = c("hgnc_symbol", "ensembl_gene_id", "ensembl_transcript_id",
                   "transcript_is_canonical", "uniprotswissprot"),
    filters = "hgnc_symbol", values = genes)
  tx <- bm %>%
    dplyr::filter(!is.na(transcript_is_canonical), transcript_is_canonical == 1) %>%
    dplyr::group_by(ensembl_transcript_id) %>%
    # a transcript can return several UniProt rows; keep one, prefer non-blank
    # (nzchar(NA) is TRUE, so test for NA explicitly)
    dplyr::arrange(dplyr::desc(!is.na(uniprotswissprot) & nzchar(uniprotswissprot)),
                   .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(hgnc_symbol, ensembl_gene_id, ensembl_transcript_id,
                  uniprotswissprot)
  missing <- setdiff(genes, tx$hgnc_symbol)
  msg("[%s] %d genes in -> %d canonical transcripts%s", label, length(genes),
      nrow(tx),
      if (length(missing)) sprintf("; %d unmapped: %s", length(missing),
                                   paste(head(missing, 10), collapse = ", ")) else "")
  tx
}

# Flatten an aenmd GRanges into transcript / key, optionally applying the
# canonical NMD-escape rule (last or penultimate exon).
aenmd_to_df <- function(gr, apply_escape_rule = NA, label = "") {
  md  <- S4Vectors::mcols(gr)
  res <- md[["res_aenmd"]]
  has_flags <- all(c("is_last", "is_penultimate") %in% names(res))
  if (is.na(apply_escape_rule)) apply_escape_rule <- has_flags
  if (apply_escape_rule && !has_flags)
    stop(sprintf("[%s] escape rule requested but is_last/is_penultimate absent", label))

  df <- data.frame(
    transcript = as.character(res[["transcript"]]),
    key        = as.character(md[["key"]]),
    stringsAsFactors = FALSE
  )
  if (apply_escape_rule) {
    keep <- which(as.logical(res[["is_last"]]) | as.logical(res[["is_penultimate"]]))
    msg("[%s] escape rule: %d / %d variants kept", label, length(keep), nrow(df))
    df <- df[keep, , drop = FALSE]
  } else {
    msg("[%s] escape rule already applied upstream; %d variants", label, nrow(df))
  }
  df
}

# Split "1:026780752|C|A" into position / ref / alt without relying on greedy
# regex backtracking.
parse_key <- function(x) {
  p <- strsplit(as.character(x), "|", fixed = TRUE)
  n <- lengths(p)
  if (any(n != 3L))
    warning(sprintf("%d keys are not chr:pos|ref|alt", sum(n != 3L)))
  data.frame(
    locus = vapply(p, function(z) if (length(z) >= 1) z[1] else NA_character_, ""),
    ref   = vapply(p, function(z) if (length(z) >= 2) z[2] else NA_character_, ""),
    alt   = vapply(p, function(z) if (length(z) >= 3) z[3] else NA_character_, ""),
    stringsAsFactors = FALSE
  )
}

# Canonical form for cross-source matching: drop "chr", unpad the position,
# uppercase the alleles. ClinVar keys and gnomAD ids are both chr:pos|ref|alt
# but zero-padding is not guaranteed to agree.
norm_key <- function(x) {
  pr <- parse_key(x)
  loc <- sub("^chr", "", pr$locus, ignore.case = TRUE)
  chr <- sub(":.*$", "", loc)
  pos <- sub("^.*:", "", loc)
  pos <- sub("^0+(?=\\d)", "", pos, perl = TRUE)
  paste0(chr, ":", pos, "|", toupper(pr$ref), "|", toupper(pr$alt))
}

# True frameshift only: net length change not a multiple of 3.
drop_inframe <- function(df, key_col, label = "") {
  pr <- parse_key(df[[key_col]])
  d  <- abs(nchar(pr$ref) - nchar(pr$alt))
  keep <- d %% 3 != 0
  msg("[%s] in-frame removed: %d / %d variants kept", label, sum(keep), nrow(df))
  df[keep, , drop = FALSE]
}

# Attach hgnc_symbol + uniprotswissprot from a canonical_tx() table.
annotate_tx <- function(df, tx) {
  dplyr::left_join(df, tx, by = c("transcript" = "ensembl_transcript_id"))
}

write_out <- function(df, stem) {
  if (is.null(df)) {
    msg("  SKIPPED %s (input missing)", stem)
    return(invisible(NULL))
  }
  nm <- if (!is.null(CFG$out_names) && stem %in% names(CFG$out_names))
          CFG$out_names[[stem]] else sprintf("%s_%s.csv", stem, CFG$tag)
  path <- out_file(nm)
  utils::write.csv(df, path, row.names = FALSE)
  msg("  wrote %-46s %6d variants  %4d transcripts",
      basename(path), nrow(df), dplyr::n_distinct(df$transcript))
  invisible(path)
}

# ---------------------------------------------------------------------------
# 2. Gene lists and canonical transcripts
# ---------------------------------------------------------------------------

msg("\n=== 1. gene lists -> canonical transcripts ===")
fs_gene  <- read_gene_list(CFG$fs_gene_file)
snv_gene <- read_gene_list(CFG$snv_gene_file)

fs_tx  <- canonical_tx(fs_gene,  "fs")
snv_tx <- canonical_tx(snv_gene, "snv")

# ---------------------------------------------------------------------------
# 3. ClinVar disease variants
# ---------------------------------------------------------------------------

msg("\n=== 2. ClinVar disease variants ===")

clinvar_set <- function(rds, tx, label, escape_rule = NA, frameshift = FALSE) {
  path <- data_file(rds)
  if (!file.exists(path)) {
    msg("[%s] SKIPPED, file not found: %s", label, path)
    return(NULL)
  }
  df <- aenmd_to_df(readRDS(path), apply_escape_rule = escape_rule, label = label)
  df <- df[df$transcript %in% tx$ensembl_transcript_id, , drop = FALSE]
  msg("[%s] on canonical transcripts of gene list: %d variants", label, nrow(df))
  if (frameshift) df <- drop_inframe(df, "key", label)
  annotate_tx(df, tx) %>% dplyr::distinct()
}

fs_plp  <- clinvar_set(CFG$fs_plp_rds,  fs_tx,  "fs P/LP",  escape_rule = TRUE, frameshift = TRUE)
snv_plp <- clinvar_set(CFG$snv_plp_rds, snv_tx, "snv P/LP", escape_rule = NA)

if (CFG$write_vus) {
  fs_vus     <- clinvar_set(CFG$fs_vus_rds,     fs_tx,  "fs VUS",     TRUE, TRUE)
  fs_benign  <- clinvar_set(CFG$fs_benign_rds,  fs_tx,  "fs benign",  TRUE, TRUE)
  snv_vus    <- clinvar_set(CFG$snv_vus_rds,    snv_tx, "snv VUS",    NA)
  snv_benign <- clinvar_set(CFG$snv_benign_rds, snv_tx, "snv benign", NA)
}

# ---------------------------------------------------------------------------
# 4. gnomAD control variants (same transcripts, ClinVar P/LP removed)
# ---------------------------------------------------------------------------

msg("\n=== 3. gnomAD control variants ===")
gnomad <- readr::read_csv(data_file(CFG$gnomad_file), show_col_types = FALSE)
stopifnot(all(c("transcript", "id", "type") %in% names(gnomad)))
msg("gnomAD table: %d variants, %d transcripts",
    nrow(gnomad), dplyr::n_distinct(gnomad$transcript))

exclude_plp <- function(df, plp, label) {
  if (is.null(plp) || !nrow(plp)) {
    msg("[%s] no ClinVar P/LP set available -- nothing excluded", label)
    return(df)
  }
  hit <- norm_key(df$id) %in% norm_key(plp$key)
  msg("[%s] ClinVar P/LP overlap removed: %d of %d", label, sum(hit), nrow(df))
  if (sum(hit) == 0)
    warning(sprintf("[%s] zero P/LP overlap -- check that ClinVar keys and ",
                    label), "gnomAD ids use the same chr:pos|ref|alt convention")
  df[!hit, , drop = FALSE]
}

# stop-gain controls: type == "snv"
snv_ctrl <- gnomad %>%
  dplyr::filter(transcript %in% snv_tx$ensembl_transcript_id, type == "snv") %>%
  exclude_plp(snv_plp, "snv control") %>%
  annotate_tx(snv_tx)

# frameshift controls: everything that is not an snv, then in-frame removed
fs_ctrl <- gnomad %>%
  dplyr::filter(transcript %in% fs_tx$ensembl_transcript_id, type != "snv") %>%
  exclude_plp(fs_plp, "fs control") %>%
  drop_inframe("id", "fs control") %>%
  annotate_tx(fs_tx)

# ---------------------------------------------------------------------------
# 5. Optional: gnomAD variants in length-matched control GENES
# ---------------------------------------------------------------------------

control_gene_set <- function(gene_file, type_filter, frameshift, label) {
  genes <- read_gene_list(gene_file)
  if (is.null(genes)) return(NULL)
  tx <- canonical_tx(genes, label)
  out <- gnomad %>%
    dplyr::filter(transcript %in% tx$ensembl_transcript_id) %>%
    dplyr::filter(if (type_filter == "snv") type == "snv" else type != "snv")
  if (frameshift) out <- drop_inframe(out, "id", label)
  annotate_tx(out, tx)
}

if (!is.null(CFG$snv_control_gene_file) || !is.null(CFG$fs_control_gene_file)) {
  msg("\n=== 4. gnomAD variants in control genes ===")
  snv_ctrl_gene <- control_gene_set(CFG$snv_control_gene_file, "snv", FALSE, "snv control genes")
  fs_ctrl_gene  <- control_gene_set(CFG$fs_control_gene_file,  "fs",  TRUE,  "fs control genes")
}

# ---------------------------------------------------------------------------
# 6. Write outputs
# ---------------------------------------------------------------------------

msg("\n=== 5. writing outputs to %s ===", CFG$out_dir)
write_out(fs_plp,   "fs_disease_variants_clinvar_plp")
write_out(snv_plp,  "snv_disease_variants_clinvar_plp")
write_out(fs_ctrl,  "fs_control_variants_gnomad")
write_out(snv_ctrl, "snv_control_variants_gnomad")

if (CFG$write_vus) {
  for (nm in c("fs_vus", "fs_benign", "snv_vus", "snv_benign")) {
    obj <- get0(nm, ifnotfound = NULL)
    if (!is.null(obj)) write_out(obj, sprintf("%s_variants_clinvar", nm))
  }
}
if (exists("snv_ctrl_gene") && !is.null(snv_ctrl_gene))
  write_out(snv_ctrl_gene, "snv_control_gene_variants_gnomad")
if (exists("fs_ctrl_gene") && !is.null(fs_ctrl_gene))
  write_out(fs_ctrl_gene, "fs_control_gene_variants_gnomad")

# ---------------------------------------------------------------------------
# 7. Per-gene summary (disease vs control counts)
# ---------------------------------------------------------------------------

msg("\n=== 6. per-gene summary ===")

gene_summary <- function(disease, control, tx, dis_col, ctl_col) {
  d <- disease %>% dplyr::group_by(hgnc_symbol) %>%
    dplyr::summarise(!!dis_col := dplyr::n_distinct(key), .groups = "drop")
  c_ <- control %>% dplyr::group_by(hgnc_symbol) %>%
    dplyr::summarise(!!ctl_col := dplyr::n_distinct(id), .groups = "drop")
  tx %>%
    dplyr::select(hgnc_symbol) %>% dplyr::distinct() %>%
    dplyr::left_join(d, by = "hgnc_symbol") %>%
    dplyr::left_join(c_, by = "hgnc_symbol") %>%
    tidyr::replace_na(stats::setNames(list(0L, 0L), c(dis_col, ctl_col))) %>%
    dplyr::arrange(dplyr::desc(.data[[dis_col]]))
}

if (requireNamespace("tidyr", quietly = TRUE) &&
    !is.null(snv_plp) && !is.null(fs_plp)) {
  snv_summary <- gene_summary(snv_plp, snv_ctrl, snv_tx, "n_clinvar_snv", "n_gnomad_snv")
  fs_summary  <- gene_summary(fs_plp,  fs_ctrl,  fs_tx,  "n_clinvar_fs",  "n_gnomad_fs")

  utils::write.csv(snv_summary, out_file(sprintf("snv_gene_summary_%s.csv", CFG$tag)), row.names = FALSE)
  utils::write.csv(fs_summary,  out_file(sprintf("fs_gene_summary_%s.csv",  CFG$tag)),  row.names = FALSE)

  # Genes with no variants on one side contribute nothing to a within-gene
  # contrast -- worth knowing before fitting (1 | gene).
  msg("snv: %d genes with 0 disease variants, %d with 0 control variants",
      sum(snv_summary$n_clinvar_snv == 0), sum(snv_summary$n_gnomad_snv == 0))
  msg("fs:  %d genes with 0 disease variants, %d with 0 control variants",
      sum(fs_summary$n_clinvar_fs == 0), sum(fs_summary$n_gnomad_fs == 0))
  msg("snv: %d genes have both; fs: %d genes have both",
      sum(snv_summary$n_clinvar_snv > 0 & snv_summary$n_gnomad_snv > 0),
      sum(fs_summary$n_clinvar_fs  > 0 & fs_summary$n_gnomad_fs  > 0))
  print(utils::head(snv_summary))
  print(utils::head(fs_summary))
} else {
  msg("per-gene summary skipped (tidyr missing, or a ClinVar set is unavailable)")
}

msg("\ndone.")
