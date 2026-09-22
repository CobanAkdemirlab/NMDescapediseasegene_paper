#!/usr/bin/env Rscript
# =============================================================================
# NMDesc -- run_analysis.R
# -----------------------------------------------------------------------------
# One entry point for the seven-step pipeline:
#
#   1  ClinVar VCF   -> SNV NMDesc-escape disease genes      gene_get_main.R
#   2  FS disease genes (per-transcript ACAT/BH enrichment)  fs_transcript_level.R
#                                                            bind_result_dbh.R
#   3  SNV + FS control gene lists                           get_snv_control_gene.R
#                                                            get_fs_control_gene.R
#   4  SNV + FS disease variant lists (ClinVar P/LP)         variant_get_main.R
#   5  SNV + FS control variant lists (gnomAD)               variant_get_main.R
#   6  gene-level comparison (CDS-matched pairs)             gene_compare_main.R
#   7  variant-level comparison (GLM / GLMM / Bayesian)      variant_compare_main.R
#
# Usage
#   Rscript run_analysis.R                 # all seven steps
#   Rscript run_analysis.R 4 5 6 7         # a subset, in the order given
#   Rscript run_analysis.R --dry-run       # preflight only, run nothing
#   NMDESC_DATA=/path/to/ship Rscript run_analysis.R
#   NMDESC_OUT=/path/to/results Rscript run_analysis.R
#
# Design notes
#   * All steps run in the GLOBAL environment and in sequence, because the
#     stage scripts pass objects to each other through globals (step 3 reads
#     `mart`, `v_ch` and `omim_AD_symbols` that step 1 leaves behind) and call
#     source() themselves, which lands in globalenv.
#   * Several stage scripts read and write bare relative filenames. Each step
#     therefore runs inside a per-step scratch directory under .run/ that is
#     seeded with the bare-relative inputs it needs (resolved through
#     data_file()) and drained into out_dir() afterwards, so the next step
#     finds them. No stage script is edited to achieve this.
#   * A step whose preflight fails does NOT abort the run. If its declared
#     outputs already exist in the data roots it is marked `cached` and the
#     chain continues from those files; otherwise it is marked `blocked` with
#     the reason, and later steps get their chance.
# =============================================================================

# ---- 0. repository root -----------------------------------------------------
.find_repo <- function() {
  p <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NA)
  if (!is.na(p) && dir.exists(file.path(p, "gene level_v5"))) return(p)
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) {
    p <- dirname(normalizePath(f[1]))
    if (dir.exists(file.path(p, "gene level_v5"))) return(p)
  }
  for (cand in c(getwd(), dirname(getwd()))) {
    cand <- path.expand(cand)
    if (dir.exists(file.path(cand, "gene level_v5"))) return(cand)
  }
  stop("cannot locate the repository root (no 'gene level_v5' directory found)",
       call. = FALSE)
}
REPO <- normalizePath(if (exists("REPO")) path.expand(REPO) else .find_repo(),
                      mustWork = TRUE)

# paths.R supplies data_file(), data_root(), out_dir() and out_file().
source(file.path(REPO, "gene level_v5/lib/paths.R"))

# A top-level lib/ so the stage scripts' "../../lib/paths.R" and "lib/paths.R"
# probes resolve from the scratch directories.
if (!file.exists(file.path(REPO, "lib")))
  try(file.symlink(file.path(REPO, "gene level_v5/lib"), file.path(REPO, "lib")),
      silent = TRUE)

RUNDIR <- file.path(REPO, ".run")
dir.create(RUNDIR, showWarnings = FALSE, recursive = TRUE)
# Repo-shaped links so a script running from .run/stepN still finds its helpers
for (d in c("gene level_v5", "variant level_v5", "lib")) {
  src <- file.path(REPO, d); dst <- file.path(RUNDIR, d)
  if (dir.exists(src) && !file.exists(dst)) try(file.symlink(src, dst), silent = TRUE)
}

FIGDIR <- file.path(REPO, "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. step table ----------------------------------------------------------
# needs_pkgs  : hard requirements; a missing one blocks the step
# needs_data  : input basenames that must resolve through data_file()
# needs_objs  : globals an earlier step must have left in the session
# stage_files : filenames a script reads or writes outside the path helpers;
#               copied into the scratch directory before the step and drained
#               back after. Empty for every step, since all of them resolve
#               through data_file() / out_file(); the hook stays available for
#               a script that does not.
# produces    : outputs later steps read by name; if all are already on disk a
#               blocked step is reported as `cached` instead
STEPS <- list(
  list(id = 1, label = "SNV disease genes  (ClinVar VCF -> NMDesc enrichment)",
       scripts     = "gene_get_main.R",
       needs_pkgs  = c("aenmd", "aenmd.data.ensdb.v105", "GenomicRanges",
                       "GenomicFeatures", "VariantAnnotation", "AnnotationDbi",
                       "BSgenome.Hsapiens.UCSC.hg38", "txdbmaker", "biomaRt",
                       "enrichR", "dplyr", "readr", "data.table"),
       needs_data  = c("clinvar_20260201.vcf.gz", "variant_summary.txt",
                       "omim_AD_symbols.csv"),
       stage_files = character(0),
       # clinvar_20260201_nmd.rds and snv_plp_ptc20260201.rds are what let
       # steps 2 and 3 run outside this session.
       produces    = c("clinvar_20260201_nmd.rds", "v_ch20260201.csv",
                       "snv_plp_ptc20260201.rds",
                       "snv_can_ADrestricted_bh_FDR0.20_all.txt")),

  list(id = 2, label = "FS disease genes  (per-transcript ACAT / BH enrichment)",
       scripts     = c("fs_transcript_level.R", "bind_result_dbh.R"),
       needs_pkgs  = c("biomaRt", "dplyr", "stringr", "readr"),
       needs_data  = c("clinvar_20260201_nmd.rds", "v_ch20260201.csv",
                       "omim_AD_symbols.csv"),
       stage_files = character(0),
       produces    = c("fs.rds", "PTC_info20260201_region.csv",
                       "fs_can_AD_acat_FDR0.20_all.txt")),

  list(id = 3, label = "Control gene lists  (SNV + FS, AD-restricted)",
       scripts     = c("get_snv_control_gene.R", "get_fs_control_gene.R"),
       needs_pkgs  = c("biomaRt", "dplyr", "stringr", "readr"),
       # Both scripts resolve their own inputs, so no session objects are
       # required: the control pool and the exclusion set are read from the two
       # files step 1 writes.
       needs_data  = c("v_ch20260201.csv", "fs.rds", "omim_AD_symbols.csv",
                       "snv_plp_ptc20260201.rds",
                       "snv_plp_ptc_nmdesc_can20260201.rds"),
       stage_files = character(0),
       produces    = c("snv_control_genes_AD.csv", "fs_control_genes_AD.csv")),

  list(id = 4, label = "Disease variant lists  (ClinVar P/LP, SNV + FS)",
       scripts     = "variant_get_main.R",
       needs_pkgs  = c("biomaRt", "dplyr", "readr", "S4Vectors", "tidyr"),
       needs_data  = c("fs.rds", "snv_plp_ptc_nmdesc_can_filtered20260201.rds",
                       "ptc_can_NMD_df.csv",
                       "snv_can_ADrestricted_bh_FDR0.20_all.txt",
                       "fs_can_AD_acat_FDR0.20_all.txt"),
       stage_files = character(0),
       produces    = c("snv_variants20260201_plp_dbh_clinvar.csv",
                       "fs_variants20260201_plp_acat_clinvar.csv")),

  list(id = 5, label = "Control variant lists  (gnomAD, SNV + FS)",
       scripts     = "variant_get_main.R",
       needs_pkgs  = c("biomaRt", "dplyr", "readr", "S4Vectors", "tidyr"),
       needs_data  = c("ptc_can_NMD_df.csv", "snv_control_genes_AD.csv",
                       "fs_control_genes_AD.csv"),
       stage_files = character(0),
       # Step 5 is step 4's script with the two control-gene inputs switched on.
       presets     = quote({
         CFG$snv_control_gene_file <- "snv_control_genes_AD.csv"
         CFG$fs_control_gene_file  <- "fs_control_genes_AD.csv"
       }),
       produces    = c("gnomad_snv_filtered_acat_0831.csv",
                       "gnomad_fs_filtered_bh_0831.csv")),

  list(id = 6, label = "Gene-level comparison  (CDS-matched pairs)",
       scripts     = "gene_compare_main.R",
       needs_pkgs  = c("biomaRt", "dplyr", "tidyr", "readr", "ggplot2",
                       "patchwork", "readxl"),
       needs_data  = c("PTC_info20260201_region.csv", "omim_AD_symbols.csv",
                       "snv_can_ADrestricted_bh_FDR0.20_all.txt",
                       "fs_can_AD_acat_FDR0.20_all.txt",
                       "snv_control_genes_AD.csv", "fs_control_genes_AD.csv",
                       "human (1).txt", "gnomad.v2.1.1.lof_metrics.by_gene.txt",
                       "NIHMS1818854-supplement-2(A).csv",
                       "NIHMS1818854-supplement-2(B).csv",
                       "Copy of NIHMS1818854-supplement-2.xls"),
       stage_files = character(0),
       produces    = c("gene_all_0826_matched.csv", "results_matched_paired.csv",
                       "result_gene_level_matched.pdf")),

  list(id = 7, label = "Variant-level comparison  (GLM / GLMM / Bayesian)",
       scripts     = "variant_compare_main.R",
       needs_pkgs  = c("biomaRt", "dplyr", "tidyr", "readr", "ggplot2",
                       "patchwork", "lme4", "broom", "broom.mixed"),
       needs_data  = c("snv_variants20260201_plp_dbh_clinvar.csv",
                       "fs_variants20260201_plp_acat_clinvar.csv",
                       "gnomad_snv_filtered_acat_0831.csv",
                       "gnomad_fs_filtered_bh_0831.csv",
                       "human (1).txt"),
       stage_files = character(0),
       produces    = c("mixed_effect_fdr.csv", "unadjusted_flags_fdr.csv",
                       "result_variant_mixed_effect.pdf"))
)

# ---- 2. helpers -------------------------------------------------------------
rule <- function(ch = "=") cat(strrep(ch, 78), "\n", sep = "")

locate <- function(name) {
  hits <- list.files(REPO, pattern = paste0("^", basename(name), "$"),
                     recursive = TRUE, full.names = TRUE)
  hits <- hits[!grepl("/(backup|[.]run)/", hits)]
  if (length(hits)) hits[1] else NA_character_
}

have_pkg <- function(p) vapply(p, requireNamespace, logical(1), quietly = TRUE)

resolve <- function(names) vapply(names, function(n) {
  p <- tryCatch(data_file(n, must = FALSE), error = function(e) NA_character_)
  if (is.na(p)) NA_character_ else p
}, "")

preflight <- function(st) {
  pk <- if (length(st$needs_pkgs)) st$needs_pkgs[!have_pkg(st$needs_pkgs)] else character(0)
  dt <- if (length(st$needs_data)) names(which(is.na(resolve(st$needs_data)))) else character(0)
  ob <- if (length(st$needs_objs))
          st$needs_objs[!vapply(st$needs_objs, exists, logical(1), envir = globalenv())]
        else character(0)
  list(pkgs = pk, data = dt, objs = ob,
       ok = length(pk) == 0 && length(dt) == 0 && length(ob) == 0)
}

# scratch dir: repo-shaped, seeded with the step's bare-relative inputs
make_scratch <- function(st) {
  d <- file.path(RUNDIR, sprintf("step%d", st$id))
  unlink(d, recursive = TRUE); dir.create(d, recursive = TRUE)
  for (sub in c("gene level_v5", "variant level_v5", "lib")) {
    src <- file.path(REPO, sub)
    if (dir.exists(src)) try(file.symlink(src, file.path(d, sub)), silent = TRUE)
  }
  # the running script's own directory, file by file, so SCRIPT_DIR-relative
  # source() calls and sibling helpers resolve
  for (s in st$scripts) {
    p <- locate(s); if (is.na(p)) next
    for (f in list.files(dirname(p), full.names = TRUE)) {
      dst <- file.path(d, basename(f))
      if (!file.exists(dst)) try(file.symlink(normalizePath(f), dst), silent = TRUE)
    }
  }
  for (n in st$stage_files) {
    p <- tryCatch(data_file(n, must = FALSE), error = function(e) NA_character_)
    if (!is.na(p)) file.copy(p, file.path(d, basename(n)), overwrite = TRUE)
  }
  d
}

# move anything the step created in the scratch dir into out_dir(), which sits
# under a data root, so data_file() finds it for later steps
drain_scratch <- function(d, seeded) {
  fs <- list.files(d, full.names = TRUE, recursive = FALSE)
  fs <- fs[!dir.exists(fs) & !basename(fs) %in% basename(seeded)]
  fs <- fs[!vapply(fs, function(f) nzchar(Sys.readlink(f)), logical(1))]
  if (!length(fs)) return(character(0))
  tgt <- out_dir(); moved <- character(0)
  for (f in fs) if (file.copy(f, file.path(tgt, basename(f)), overwrite = TRUE))
    moved <- c(moved, basename(f))
  moved
}

STATUS <- list()
note <- function(st, state, detail = "", files = character(0))
  STATUS[[length(STATUS) + 1]] <<- list(id = st$id, label = st$label,
                                        state = state, detail = detail,
                                        n_files = length(files))

run_step <- function(st, dry = FALSE) {
  cat("\n"); rule()
  cat(sprintf("step %d  %s\n", st$id, st$label))
  rule()

  pf <- preflight(st)
  if (length(pf$pkgs)) cat("  missing packages :", paste(pf$pkgs, collapse = ", "), "\n")
  if (length(pf$data)) cat("  missing inputs   :", paste(pf$data, collapse = ", "), "\n")
  if (length(pf$objs)) cat("  missing objects  :", paste(pf$objs, collapse = ", "),
                           "\n                     (an earlier step in this session must define them)\n")
  if (pf$ok) cat("  preflight        : ok\n")

  # How much of this step's output is already on disk decides whether a blocked
  # step can be skipped over: all of it -> cached, some -> partial, none -> blocked.
  have <- if (length(st$produces)) resolve(st$produces) else character(0)
  n_have <- sum(!is.na(have))

  if (!pf$ok) {
    reason <- paste(c(pf$pkgs, pf$data, pf$objs), collapse = ", ")
    if (n_have == length(st$produces) && n_have > 0) {
      cat("  -> skipped: outputs already present, the chain continues from them\n")
      for (p in have) cat("       ", sub(path.expand("~"), "~", p), "\n")
      note(st, "cached", reason)
    } else if (n_have > 0) {
      cat(sprintf("  -> PARTIAL: %d of %d outputs on disk, the rest are missing\n",
                  n_have, length(st$produces)))
      for (i in seq_along(have))
        cat(sprintf("       %-6s %s\n", if (is.na(have[i])) "absent" else "ok",
                    if (is.na(have[i])) st$produces[i]
                    else sub(path.expand("~"), "~", have[i])))
      note(st, "partial", reason)
    } else {
      cat("  -> BLOCKED: preflight unmet and no cached outputs\n")
      note(st, "blocked", reason)
    }
    return(invisible())
  }
  if (dry) { cat("  -> dry run, not executed\n"); note(st, "dry-run"); return(invisible()) }

  scratch <- make_scratch(st)
  wd <- getwd(); res <- "ok"
  for (s in st$scripts) {
    p <- locate(s)
    if (is.na(p)) { res <- paste("script not found:", s); break }
    cat("  running", s, "\n")
    r <- tryCatch({
      setwd(scratch)
      if (!is.null(st$presets)) eval(st$presets, envir = globalenv())
      sys.source(normalizePath(p), envir = globalenv())
      "ok"
    }, error = function(e) paste("failed:", conditionMessage(e)))
    setwd(wd)
    if (r != "ok") { res <- paste0(s, ": ", r); break }
  }
  setwd(wd)
  moved <- drain_scratch(scratch, st$stage_files)
  if (length(moved))
    cat("  outputs moved to", sub(path.expand("~"), "~", out_dir()), ":",
        paste(moved, collapse = ", "), "\n")
  cat("  ->", res, "\n")
  note(st, if (res == "ok") "ok" else "failed", if (res == "ok") "" else res, moved)
}

# ---- 3. drive ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
dry  <- any(args %in% c("--dry-run", "-n"))
sel  <- suppressWarnings(as.integer(args[grepl("^[0-9]+$", args)]))
sel  <- sel[!is.na(sel)]
if (!length(sel)) sel <- 1:7

rule()
cat("NMDesc pipeline\n")
cat("  repo      :", sub(path.expand("~"), "~", REPO), "\n")
cat("  data roots:", paste(sub(path.expand("~"), "~", DATA_ROOTS[dir.exists(DATA_ROOTS)]),
                           collapse = "\n               "), "\n")
cat("  out_dir() :", sub(path.expand("~"), "~", out_dir()), "\n")
cat("  figures   :", sub(path.expand("~"), "~", FIGDIR), "\n")
cat("  steps     :", paste(sel, collapse = " "), if (dry) "(dry run)" else "", "\n")
rule()

for (i in sel) {
  st <- Filter(function(s) s$id == i, STEPS)
  if (!length(st)) { cat("\nno such step:", i, "\n"); next }
  run_step(st[[1]], dry = dry)
}

# ---- 4. summary -------------------------------------------------------------
cat("\n"); rule(); cat("summary\n"); rule()
for (s in STATUS)
  cat(sprintf("  step %d  %-8s  %s%s\n", s$id, s$state, s$label,
              if (nzchar(s$detail)) paste0("\n            ", s$detail) else ""))

cat("\n"); rule(); cat("generated files\n"); rule()
dirs <- unique(c(FIGDIR, tryCatch(out_dir(), error = function(e) NULL),
                 list.files(RUNDIR, pattern = "^step[0-9]+$", full.names = TRUE)))
dirs <- dirs[dir.exists(dirs)]
for (d in dirs) {
  fs <- list.files(d, pattern = "[.](pdf|png|csv|rds|txt)$", full.names = TRUE)
  fs <- fs[!vapply(fs, function(f) nzchar(Sys.readlink(f)), logical(1))]
  cat("\n ", sub(path.expand("~"), "~", d), "\n")
  if (!length(fs)) { cat("    no files\n"); next }
  fs <- fs[order(file.mtime(fs), decreasing = TRUE)]
  for (f in head(fs, 25)) cat(sprintf("   %9.0f KB  %s\n", file.size(f) / 1024, basename(f)))
  if (length(fs) > 25) cat(sprintf("   ... %d more\n", length(fs) - 25))
}
