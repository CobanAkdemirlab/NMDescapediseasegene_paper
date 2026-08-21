# =============================================================================
# gene_main.R
# -----------------------------------------------------------------------------
# Gene-level comparison of NMD-escape (NMDesc) disease genes vs. controls.
#
# Structure
#   3.   gene_level comparision
#   3.0  input data      – libraries, config, functions, load gene_all
#   3.1  add gene level annotation
#   3.2  see if cds length and NMDesc region length could influence the result
# =============================================================================


# #############################################################################
# 3. gene_level comparision
# #############################################################################


# =============================================================================
# 3.0 input data
# =============================================================================

# ---- Libraries -------------------------------------------------------------
suppressPackageStartupMessages({
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(scales)
  library(broom)
  library(data.table)   # fread
  library(biomaRt)      # getBM / useEnsembl
  library(STRINGdb)
  library(igraph)
  library(GenomicRanges)
  library(IRanges)
})

# ---- Resolve namespace conflicts -------------------------------------------
## biomaRt/AnnotationDbi, STRINGdb, GenomicRanges and IRanges all export their
## own select()/filter()/rename()/etc., which mask dplyr's after they load.
## Force dplyr's verbs to win for the rest of the script. (GenomicRanges/IRanges
## calls that we DO want are already fully qualified, e.g. GenomicRanges::reduce.)
select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
slice     <- dplyr::slice
count     <- dplyr::count
first     <- dplyr::first
setdiff   <- dplyr::setdiff
union     <- dplyr::union
intersect <- dplyr::intersect

# ---- Configuration ---------------------------------------------------------
## All external inputs live here. Point these at your local copies.
CONFIG <- list(
  
  # --- Inputs to build gene_all (section 3.0) ---
  ptc_info            = "PTC_info20260201_region.csv",
  omim_ad_symbols     = "omim_AD_symbols.csv",
  snv_gene_list       = "~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/lists/FDR0.05/snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201_NMDesc_wald_enriched_can.txt",
  fs_gene_list        = "fs_can_AD_FDR0.05_wald_gene.csv",
  snv_control_list    = "snv_control_genes_AD.csv",
  fs_control_list     = "fs_control_genes_AD.csv",
  
  # --- Feature-annotation inputs (section 3.1) ---
  path_touni          = "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/others/NIHMS1818854-supplement-2(A).csv",
  path_motif          = "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/others/NIHMS1818854-supplement-2(B).csv",
  path_lcs            = "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/others/Copy of NIHMS1818854-supplement-2.xls",
  ppi_file            = "~/Downloads/human (1) (2).txt",
  gtex_path           = "/Users/jxu14/Downloads/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct",
  lof_metrics_path    = "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/others/gnomad.v2.1.1.lof_metrics.by_gene.txt",
  
  # --- Shared plotting ---
  group_colors = c(
    "snv"         = "#1f77b4",
    "snv_control" = "#aec7e8",
    "fs"          = "#ff7f0e",
    "fs_control"  = "#ffbb78"
  ),
  comparisons = list(c("fs", "fs_control"), c("snv", "snv_control"))
)

# =============================================================================
# 3.0 (cont.) function definitions
#   Defined here so section 3.1 can call them. Each is self-contained and
#   writes its own CSV / figure outputs.
# =============================================================================

# fn: get_gene_info: fetch all required info for one gene list -----------
get_gene_info <- function(
    hgnc_symbols   = NULL,       # provide either hgnc_symbols OR ensembl_tx_ids
    ensembl_tx_ids = NULL,
    mart,
    PTC_combined   = NULL,       # required only for fs / fs_control groups
    group_type     = c("snv", "fs", "snv_control", "fs_control")
) {
  group_type       <- match.arg(group_type)
  filter_by_symbol <- !is.null(hgnc_symbols)
  filter_val       <- if (filter_by_symbol) unique(hgnc_symbols) else unique(ensembl_tx_ids)
  filter_name      <- if (filter_by_symbol) "hgnc_symbol" else "ensembl_transcript_id"
  
  # 1. CDS sequence (canonical only, non-empty)
  cds_df <- getBM(
    attributes = c("hgnc_symbol", "ensembl_transcript_id",
                   "transcript_is_canonical", "coding"),
    filters    = filter_name,
    values     = filter_val,
    mart       = mart
  ) %>%
    filter(transcript_is_canonical == 1, !is.na(coding), coding != "")
  
  # 2. NMD-escape region
  if (group_type %in% c("snv", "snv_control")) {
    
    # Exon-level CDS coordinates
    exon_info <- getBM(
      attributes = c("ensembl_transcript_id", "rank", "cds_start", "cds_end"),
      filters    = filter_name,
      values     = filter_val,
      mart       = mart
    )
    
    # Canonical transcript map
    tx_info <- getBM(
      attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_is_canonical"),
      filters    = filter_name,
      values     = filter_val,
      mart       = mart
    ) %>%
      filter(transcript_is_canonical == 1) %>%
      distinct(hgnc_symbol, ensembl_transcript_id, .keep_all = TRUE)
    
    exon_info <- exon_info %>%
      inner_join(tx_info, by = "ensembl_transcript_id") %>%
      filter(!is.na(rank), !is.na(cds_start), !is.na(cds_end))
    
    # NMD-escape region = last exon + 50 bp upstream of penultimate exon boundary
    nmdesc_info <- exon_info %>%
      group_by(ensembl_transcript_id) %>%
      summarise(
        nmdesc_end       = max(cds_end),
        last_exon_length = cds_end[which.max(rank)] - cds_start[which.max(rank)] + 1,
        nmdesc_start     = max(cds_end) - 50 - last_exon_length,
        .groups = "drop"
      )
    
    nmdesc_df <- cds_df %>%
      inner_join(nmdesc_info, by = "ensembl_transcript_id") %>%
      mutate(
        NMD_region_start = nmdesc_start,
        NMD_region_end   = nmdesc_end,
        nmdesc_cds       = mapply(substr, coding, nmdesc_start, nmdesc_end)
      ) %>%
      select(-any_of(c("last_exon_length", "nmdesc_start", "nmdesc_end")))
    
  } else {
    # fs / fs_control: use median PTC region across plus1 / plus2 alleles
    stopifnot(!is.null(PTC_combined),
              all(c("transcript", "median_can_region_start",
                    "median_can_region_end") %in% names(PTC_combined)))
    
    nmdesc_df <- cds_df %>%
      inner_join(PTC_combined, by = c("ensembl_transcript_id" = "transcript")) %>%
      mutate(
        NMD_region_start = median_can_region_start,
        NMD_region_end   = median_can_region_end,
        nmdesc_cds       = mapply(substr, coding,
                                  median_can_region_start, median_can_region_end)
      ) %>%
      select(-median_can_region_start, -median_can_region_end)
  }
  
  nmdesc_df %>%
    mutate(
      transcript_is_canonical = as.character(transcript_is_canonical),
      group                   = group_type
    )
}


# fn: build_gene_all -----------------------------------------------------
build_gene_all <- function(
    snv_gene,
    fs_gene,
    snv_control_gene_AD,   # vector of transcript IDs
    fs_control_gene_AD,    # vector of transcript IDs
    PTC_info,
    mart,
    omim_ad_symbols_path = "omim_AD_symbols.csv",
    output_csv           = "gene_all.csv"
) {
  # Pre-compute PTC median region (shared by fs and fs_control)
  PTC_combined <- PTC_info %>%
    group_by(transcript) %>%
    summarise(
      median_can_region_start = round(median(can_region_start, na.rm = TRUE)),
      median_can_region_end   = round(median(can_region_end,   na.rm = TRUE)),
      .groups = "drop"
    )
  
  # Fetch info for each group
  snv_df         <- get_gene_info(hgnc_symbols   = snv_gene,            mart = mart, group_type = "snv")
  fs_df          <- get_gene_info(hgnc_symbols   = fs_gene,             mart = mart, group_type = "fs",         PTC_combined = PTC_combined)
  snv_control_df <- get_gene_info(ensembl_tx_ids = snv_control_gene_AD, mart = mart, group_type = "snv_control")
  fs_control_df  <- get_gene_info(ensembl_tx_ids = fs_control_gene_AD,  mart = mart, group_type = "fs_control", PTC_combined = PTC_combined)
  
  # Combine
  gene_all <- bind_rows(snv_df, fs_df, snv_control_df, fs_control_df) %>%
    mutate(NMDesc_region_length = NMD_region_end - NMD_region_start + 1,
           row_id = row_number())
  
  # Append CDS length
  cds_len <- getBM(
    attributes = c("ensembl_transcript_id", "cds_length"),
    filters    = "ensembl_transcript_id",
    values     = unique(gene_all$ensembl_transcript_id),
    mart       = mart
  )
  gene_all <- gene_all %>% left_join(cds_len, by = "ensembl_transcript_id")
  
  # Append UniProt Swiss-Prot ID
  uniprot <- getBM(
    attributes = c("ensembl_transcript_id", "uniprotswissprot"),
    filters    = "ensembl_transcript_id",
    values     = unique(gene_all$ensembl_transcript_id),
    mart       = mart
  )
  gene_all <- gene_all %>%
    left_join(uniprot, by = "ensembl_transcript_id") %>%
    rename(uniprot = uniprotswissprot)
  
  # Keep AD genes only
  omim_AD_symbols <- read.csv(omim_ad_symbols_path, header = FALSE)$V1
  gene_all <- gene_all %>% filter(hgnc_symbol %in% omim_AD_symbols)
  
  write.csv(gene_all, output_csv, row.names = FALSE)
  message("Saved: ", output_csv)
  invisible(gene_all)
}


# fn: calculate_ppi_degree_centrality ------------------------------------
calculate_ppi_degree_centrality <- function(
    gene_input,
    output_csv      = "centrality_data.csv",
    output_fig      = "centrality_plot.png",
    string_version  = "11.5",
    species         = 9606,
    score_threshold = 400
) {
  # STRING setup
  string_db <- STRINGdb$new(
    version         = string_version,
    species         = species,
    score_threshold = score_threshold,
    input_directory = ""
  )
  
  # Map genes to STRING IDs. string_db$map() keeps the input columns
  # (including `group`) and adds STRING_id, so no second join is needed.
  mapped_genes <- string_db$map(gene_input, "hgnc_symbol", removeUnmappedRows = TRUE)
  
  # Interactions -> undirected network
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  ppi_network <- graph_from_data_frame(
    data.frame(from = interactions$from, to = interactions$to),
    directed = FALSE
  )
  ppi_network <- simplify(ppi_network)
  
  # Degree centrality
  degree_scores <- degree(ppi_network, normalized = TRUE)
  degree_df     <- data.frame(STRING_id = names(degree_scores), Degree = degree_scores)
  
  # Merge with group labels
  mapped_genes_degree <- mapped_genes %>%
    dplyr::select(STRING_id, hgnc_symbol, group) %>%
    dplyr::distinct()
  
  centrality_data <- left_join(degree_df, mapped_genes_degree, by = "STRING_id") %>%
    dplyr::filter(!is.na(group))
  
  # Plot
  p <- ggplot(centrality_data, aes(x = group, y = Degree, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    stat_compare_means(
      method      = "t.test",
      comparisons = list(c("fs", "fs_control"), c("snv", "snv_control")),
      label       = "p.signif"
    ) +
    scale_y_continuous(trans = "log10") +
    scale_fill_manual(values = CONFIG$group_colors) +
    labs(title = "Degree Centrality of Genes in STRING PPI Network",
         y = "Degree Centrality (log10)", x = "Gene Category") +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.title.y    = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )
  
  write.csv(centrality_data, output_csv, row.names = FALSE)
  message("CSV saved to: ", output_csv)
  ggsave(output_fig, plot = p, width = 7, height = 6, dpi = 300)
  message("Figure saved to: ", output_fig)
  
  invisible(list(plot = p, data = centrality_data))
}


# fn: plot_gc_content ----------------------------------------------------
plot_gc_content <- function(
    gene_all,
    comparisons = CONFIG$comparisons,
    output_csv  = "gc_content.csv",
    output_fig  = "gc_content.png"
) {
  get_gc_content <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(NA_real_)
    bases <- strsplit(toupper(sequence), "")[[1]]
    round(sum(bases %in% c("G", "C")) / length(bases) * 100, 2)
  }
  
  gene_all$gc_content        <- sapply(gene_all$coding,     get_gc_content)
  gene_all$nmdesc_gc_content <- sapply(gene_all$nmdesc_cds, get_gc_content)
  
  make_plot <- function(y_var, y_label, title) {
    ggplot(gene_all, aes(x = group, y = .data[[y_var]], fill = group)) +
      geom_boxplot(width = 0.7, outlier.shape = 16, outlier.size = 1.5) +
      stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.format") +
      scale_y_continuous(labels = function(x) paste0(x, "%")) +
      scale_fill_manual(values = CONFIG$group_colors) +
      labs(x = "", y = y_label, title = title) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  }
  
  p1 <- make_plot("gc_content",        "GC content (%)", "Overall Gene GC Content")
  p2 <- make_plot("nmdesc_gc_content", "GC content (%)", "GC Content in NMD-escape Regions")
  combined <- ggarrange(p1, p2, ncol = 2, nrow = 1)
  
  write.csv(
    gene_all[, c("hgnc_symbol", "ensembl_transcript_id", "group",
                 "gc_content", "nmdesc_gc_content")],
    output_csv, row.names = FALSE
  )
  message("CSV saved to: ", output_csv)
  ggsave(output_fig, plot = combined, width = 10, height = 6, dpi = 300)
  message("Figure saved to: ", output_fig)
  
  invisible(list(plot = combined, data = gene_all))
}


# fn: plot_repeat_content ------------------------------------------------
plot_repeat_content <- function(
    gene_all,
    comparisons = CONFIG$comparisons,
    output_csv  = "repeat_content.csv",
    output_fig  = "repeat_content.png"
) {
  detect_repeats <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(0)
    matches <- str_extract_all(sequence, "([ATGC]{1,6})\\1+")[[1]]
    if (length(matches) == 0) return(0)
    sum(nchar(matches))
  }
  repeat_fraction <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(NA_real_)
    detect_repeats(sequence) / nchar(sequence)
  }
  detect_homopolymer <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(0)
    matches <- str_extract_all(sequence, "([ATGC])\\1{3,}")[[1]]
    if (length(matches) == 0) return(0)
    sum(nchar(matches))
  }
  homopolymer_fraction <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(NA_real_)
    detect_homopolymer(sequence) / nchar(sequence)
  }
  
  gene_all$repeat_fraction             <- sapply(gene_all$coding,     repeat_fraction)
  gene_all$nmdesc_repeat_fraction      <- sapply(gene_all$nmdesc_cds, repeat_fraction)
  gene_all$homopolymer_fraction        <- sapply(gene_all$coding,     homopolymer_fraction)
  gene_all$nmdesc_homopolymer_fraction <- sapply(gene_all$nmdesc_cds, homopolymer_fraction)
  
  make_plot <- function(y_var, y_label, title) {
    ggplot(gene_all, aes(x = group, y = .data[[y_var]], fill = group)) +
      geom_boxplot(width = 0.7, outlier.shape = 16, outlier.size = 1.5) +
      stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.format") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      scale_fill_manual(values = CONFIG$group_colors) +
      labs(x = "", y = y_label, title = title) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  }
  
  p1 <- make_plot("repeat_fraction",             "Repeat fraction",      "Overall Repeat Content")
  p2 <- make_plot("nmdesc_repeat_fraction",      "Repeat fraction",      "NMD-escape Repeat Content")
  p3 <- make_plot("homopolymer_fraction",        "Homopolymer fraction", "Overall Homopolymer Content")
  p4 <- make_plot("nmdesc_homopolymer_fraction", "Homopolymer fraction", "NMD-escape Homopolymer Content")
  combined <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
  
  out_cols <- c("hgnc_symbol", "ensembl_transcript_id", "group",
                "repeat_fraction", "nmdesc_repeat_fraction",
                "homopolymer_fraction", "nmdesc_homopolymer_fraction")
  write.csv(gene_all[, out_cols], output_csv, row.names = FALSE)
  message("CSV saved to: ", output_csv)
  ggsave(output_fig, plot = combined, width = 12, height = 10, dpi = 300)
  message("Figure saved to: ", output_fig)
  
  invisible(list(plot = combined, data = gene_all))
}


# fn: annotate_motif_flags -----------------------------------------------
annotate_motif_flags <- function(
    gene_all,
    path_touni,
    path_motif,
    path_LCS,
    mart,
    output_motif_csv = "gene_motif_flags.csv",
    output_lcs_csv   = "gene_LCS_flags.csv"
) {
  get_max_from_cell <- function(x) {
    if (is.na(x) || x == "") return(NA_real_)
    nums <- str_extract_all(x, "\\d+")[[1]]
    if (length(nums) == 0) return(NA_real_)
    max(as.numeric(nums), na.rm = TRUE)
  }
  parse_max <- function(df) {
    df %>% mutate(across(-1, ~ vapply(as.character(.x), get_max_from_cell, numeric(1))))
  }
  add_uniprot <- function(df, touni) {
    df %>% mutate(uniprot = touni$`Uniprot ID`[match(Protein, touni$Protein)])
  }
  flag_cols <- function(merged_df, col_map) {
    for (flag_name in names(col_map)) {
      src_col <- col_map[[flag_name]]
      if (!src_col %in% colnames(merged_df)) {
        warning("Column '", src_col, "' not found. Available: ",
                paste(colnames(merged_df), collapse = ", "))
        next
      }
      merged_df[[flag_name]] <- merged_df[[src_col]] * 3 >= merged_df$NMD_region_start
    }
    merged_df
  }
  
  # Load external data
  touni     <- read_csv(path_touni, show_col_types = FALSE)
  motif_doc <- read_csv(path_motif, show_col_types = FALSE)
  LCS_doc   <- read_excel(path_LCS, sheet = "G")
  
  motif_max <- motif_doc %>% parse_max() %>% add_uniprot(touni)
  LCS_max   <- LCS_doc   %>% parse_max() %>% add_uniprot(touni)
  
  # Replace empty uniprot strings with NA
  gene_all <- gene_all %>% mutate(uniprot = na_if(trimws(uniprot), ""))
  
  # Re-fetch uniprot for genes where it is missing
  missing_genes <- gene_all %>% filter(is.na(uniprot)) %>% pull(hgnc_symbol) %>% unique()
  if (length(missing_genes) > 0) {
    message(length(missing_genes), " genes missing UniProt ID — re-fetching from BioMart")
    uni_map <- getBM(
      attributes = c("hgnc_symbol", "uniprotswissprot"),
      filters    = "hgnc_symbol",
      values     = missing_genes,
      mart       = mart
    ) %>%
      filter(uniprotswissprot != "") %>%
      distinct(hgnc_symbol, .keep_all = TRUE)
    
    gene_all <- gene_all %>%
      left_join(uni_map, by = "hgnc_symbol") %>%
      mutate(uniprot = coalesce(uniprot, uniprotswissprot)) %>%
      select(-uniprotswissprot)
  }
  message("Genes with UniProt after re-fetch: ",
          sum(!is.na(gene_all$uniprot)), " / ", nrow(gene_all))
  
  # Motif flags
  motif_merged <- merge(motif_max, gene_all, by = "uniprot") %>%
    flag_cols(list(
      gene_protein_flag = "Protein Features",
      gene_domains_flag = "Domains",
      gene_slim_flag    = "SLiMs",
      gene_morf_flag    = "MORFs",
      gene_ptm_flag     = "PTMs",
      gene_nls_flag     = "NLSs/NESs"
    ))
  
  # LCS flags
  lcs_merged <- merge(LCS_max, gene_all, by = "uniprot") %>%
    flag_cols(list(gene_LCS_flag = "LCSs"))
  
  # Merge flags back to gene_all
  motif_flags <- c("gene_protein_flag", "gene_domains_flag", "gene_slim_flag",
                   "gene_morf_flag", "gene_ptm_flag", "gene_nls_flag")
  for (flag in motif_flags) {
    gene_all[[flag]] <- motif_merged[[flag]][match(gene_all$uniprot, motif_merged$uniprot)]
  }
  gene_all$gene_LCS_flag <- lcs_merged$gene_LCS_flag[match(gene_all$uniprot, lcs_merged$uniprot)]
  
  write.csv(motif_merged, output_motif_csv, row.names = FALSE)
  write.csv(lcs_merged,   output_lcs_csv,   row.names = FALSE)
  message("Motif flags saved to: ", output_motif_csv)
  message("LCS flags saved to:   ", output_lcs_csv)
  
  invisible(list(gene_all = gene_all, motif_merged = motif_merged, lcs_merged = lcs_merged))
}


# fn: run_pfam_overlap_analysis ------------------------------------------
run_pfam_overlap_analysis <- function(gene_all, ensembl, output_prefix = "pfam_overlap") {
  
  # 1. PFAM domains
  pfam_domains <- getBM(
    attributes = c("ensembl_transcript_id", "pfam", "pfam_start", "pfam_end"),
    filters    = "ensembl_transcript_id",
    values     = unique(gene_all$ensembl_transcript_id),
    mart       = ensembl
  )
  
  # 2. Merge PFAM info, convert aa coordinates to bp
  gene_all_pfam <- gene_all %>%
    left_join(pfam_domains, by = "ensembl_transcript_id") %>%
    mutate(pfam_start_bp = pfam_start * 3 - 2,
           pfam_end_bp   = pfam_end * 3)
  
  # 3. Per-record overlap with the NMD-escape region
  gene_all_pfam <- gene_all_pfam %>%
    mutate(
      overlap_start      = pmax(NMD_region_start, pfam_start_bp),
      overlap_end        = pmin(NMD_region_end, pfam_end_bp),
      overlap_valid      = !is.na(overlap_start) & !is.na(overlap_end) & (overlap_start <= overlap_end),
      overlap_length_raw = ifelse(overlap_valid, overlap_end - overlap_start + 1, 0)
    )
  
  # 4. Reduce intervals per row_id to avoid double counting
  idx_list <- split(seq_len(nrow(gene_all_pfam)), gene_all_pfam$row_id)
  
  pfam_overlap_summary <- lapply(idx_list, function(idx) {
    df          <- gene_all_pfam[idx, , drop = FALSE]
    nm_len      <- df$NMDesc_region_length[1]
    this_row_id <- df$row_id[1]
    df_valid    <- df[df$overlap_valid %in% TRUE, , drop = FALSE]
    
    if (nrow(df_valid) == 0 || is.na(nm_len) || nm_len <= 0) {
      return(data.frame(
        row_id                = this_row_id,
        pfam_overlap_length   = 0,
        pfam_overlap_flag     = 0,
        pfam_overlap_fraction = 0,
        n_overlapping_pfam    = 0,
        stringsAsFactors      = FALSE
      ))
    }
    
    gr <- GenomicRanges::GRanges(
      seqnames = rep("x", nrow(df_valid)),
      ranges   = IRanges::IRanges(start = as.numeric(df_valid$overlap_start),
                                  end   = as.numeric(df_valid$overlap_end))
    )
    gr_red        <- GenomicRanges::reduce(gr)
    total_overlap <- sum(IRanges::width(gr_red))
    
    data.frame(
      row_id                = this_row_id,
      pfam_overlap_length   = total_overlap,
      pfam_overlap_flag     = as.integer(total_overlap > 20),
      pfam_overlap_fraction = total_overlap / nm_len,
      n_overlapping_pfam    = nrow(df_valid),
      stringsAsFactors      = FALSE
    )
  })
  pfam_overlap_summary <- dplyr::bind_rows(pfam_overlap_summary)
  
  # 5. Merge back
  gene_all <- gene_all %>% left_join(pfam_overlap_summary, by = "row_id")
  
  # 6. Sanity checks
  print(table(gene_all$pfam_overlap_flag, useNA = "ifany"))
  print(summary(gene_all$pfam_overlap_fraction))
  print(summary(gene_all$pfam_overlap_length))
  
  # 7. Group order / colors
  gene_all$group <- factor(gene_all$group,
                           levels = c("snv", "snv_control", "fs", "fs_control"))
  group_colors2 <- c("snv" = "#2ca02c", "snv_control" = "#98df8a",
                     "fs"  = "#1f77b4", "fs_control"  = "#aec7e8")
  comparisons <- list(c("snv", "snv_control"), c("fs", "fs_control"))
  
  # 8. Fraction
  p1 <- ggplot(gene_all, aes(x = group, y = pfam_overlap_fraction, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
    scale_fill_manual(values = group_colors2, guide = "none") +
    labs(title = "PFAM overlap fraction across variant groups",
         x = NULL, y = "PFAM overlap fraction") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold")) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.format")
  
  # 9. Length
  p2 <- ggplot(gene_all, aes(x = group, y = pfam_overlap_length, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
    scale_fill_manual(values = group_colors2, guide = "none") +
    labs(title = "PFAM overlap length across variant groups",
         x = NULL, y = "PFAM overlap length (bp)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold")) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.format")
  
  # 10. Flag proportion
  pfam_flag_summary <- gene_all %>%
    group_by(group) %>%
    summarise(n = n(),
              n_overlap    = sum(pfam_overlap_flag == 1, na.rm = TRUE),
              prop_overlap = n_overlap / n,
              .groups = "drop")
  
  p3 <- ggplot(pfam_flag_summary, aes(x = group, y = prop_overlap, fill = group)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = group_colors2, guide = "none") +
    labs(title = "Proportion of variants with PFAM overlap > 20 bp",
         x = NULL, y = "Proportion") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold")) +
    geom_text(aes(label = paste0(n_overlap, "/", n)), vjust = -0.3, size = 4)
  
  print(p1); print(p2); print(p3)
  
  # 11. Fisher tests
  for (grp_pair in list(c("snv", "snv_control"), c("fs", "fs_control"))) {
    cat("\n--- Fisher test:", grp_pair[1], "vs", grp_pair[2], "---\n")
    tab <- table(gene_all$group[gene_all$group %in% grp_pair],
                 gene_all$pfam_overlap_flag[gene_all$group %in% grp_pair])
    print(tab); print(fisher.test(tab))
  }
  
  # 12. Save
  write.csv(gene_all, paste0(output_prefix, "_gene_all.csv"), row.names = FALSE)
  ggsave(paste0(output_prefix, "_fraction_plot.png"), p1, width = 8, height = 6, dpi = 300)
  ggsave(paste0(output_prefix, "_length_plot.png"),   p2, width = 8, height = 6, dpi = 300)
  ggsave(paste0(output_prefix, "_flag_barplot.png"),  p3, width = 8, height = 6, dpi = 300)
  
  invisible(list(gene_all = gene_all,
                 pfam_overlap_summary = pfam_overlap_summary,
                 plots = list(p1 = p1, p2 = p2, p3 = p3)))
}


# fn: run_ppi_overlap_analysis -------------------------------------------
run_ppi_overlap_analysis <- function(gene_all, ppi_file_path, output_prefix = "ppi_overlap") {
  
  # 1. Read PPI network data
  human_ppi <- fread(ppi_file_path, data.table = FALSE)
  
  # 2. Helper: parse a "[a, b, c]" string into a numeric vector
  convert_to_c <- function(x) {
    x <- gsub("\\[|\\]", "", x)
    if (x == "") return(c())
    as.numeric(trimws(unlist(strsplit(x, ","))))
  }
  
  # 3. Flag genes whose interface residues fall in the NMD-escape region
  gene_all$matched_uniprot <- 0
  for (i in seq_len(nrow(gene_all))) {
    uid         <- gene_all$uniprot[i]
    aft_NMD_ind <- gene_all$NMD_region_start[i]
    if (is.na(uid) || uid == "") next
    if (is.na(aft_NMD_ind)) next
    
    re_1 <- unlist(lapply(human_ppi$interface_residues1[human_ppi$uniprot1 == uid], convert_to_c)) * 3 - 2
    re_2 <- unlist(lapply(human_ppi$interface_residues2[human_ppi$uniprot2 == uid], convert_to_c)) * 3 - 2
    if (length(re_1) == 0 && length(re_2) == 0) next
    
    if (any(re_1 >= aft_NMD_ind, na.rm = TRUE) || any(re_2 >= aft_NMD_ind, na.rm = TRUE)) {
      gene_all$matched_uniprot[i] <- 1
    }
  }
  
  # 4. Rename to ppi_overlap
  gene_all <- gene_all %>% rename(ppi_overlap = matched_uniprot)
  
  # 5. Summary
  ppi_summary <- gene_all %>%
    group_by(group) %>%
    summarise(total_genes   = n(),
              matched_genes = sum(ppi_overlap),
              percentage    = matched_genes / total_genes * 100,
              .groups = "drop")
  print(ppi_summary)
  
  # 6. Group order / colors
  gene_all$group <- factor(gene_all$group, levels = c("snv", "snv_control", "fs", "fs_control"))
  group_colors <- c("snv" = "#2ca02c", "snv_control" = "#98df8a",
                    "fs"  = "#1f77b4", "fs_control"  = "#aec7e8")
  ppi_summary$group <- factor(ppi_summary$group, levels = c("snv", "snv_control", "fs", "fs_control"))
  
  # 7. Bar plot
  p <- ggplot(ppi_summary, aes(x = group, y = percentage, fill = group)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = group_colors, guide = "none") +
    labs(title = "Percentage of genes with PPI network overlap", x = NULL, y = "Percentage (%)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold")) +
    geom_text(aes(label = paste0(matched_genes, "/", total_genes)), vjust = -0.3, size = 4)
  print(p)
  
  # 8. Fisher tests
  for (grp_pair in list(c("snv", "snv_control"), c("fs", "fs_control"))) {
    cat("\n--- Fisher test:", grp_pair[1], "vs", grp_pair[2], "---\n")
    tab <- table(gene_all$group[gene_all$group %in% grp_pair],
                 gene_all$ppi_overlap[gene_all$group %in% grp_pair])
    print(tab); print(fisher.test(tab))
  }
  
  # 9. Save
  write.csv(gene_all,    paste0(output_prefix, "_gene_all.csv"), row.names = FALSE)
  write.csv(ppi_summary, paste0(output_prefix, "_summary.csv"),  row.names = FALSE)
  ggsave(paste0(output_prefix, "_barplot.png"), p, width = 8, height = 6, dpi = 300)
  
  invisible(list(gene_all = gene_all, ppi_summary = ppi_summary, plot = p))
}


# fn: run_tau_analysis ---------------------------------------------------
run_tau_analysis <- function(gene_all, gtex_path, output_prefix = "tau") {
  
  # 1. GTEx median TPM
  gtex <- read_tsv(gtex_path, skip = 2)
  gtex$Name <- sub("\\..*", "", gtex$Name)
  
  gene_expr <- gtex %>%
    dplyr::select(-tidyselect::any_of("Name")) %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    dplyr::rename(Gene = Description)
  
  write.csv(gene_expr, paste0(output_prefix, "_gene_matrix.csv"), row.names = FALSE)
  
  # 2. Build expression matrix
  gene.df     <- as.data.frame(gene_expr)
  gene.matrix <- as.matrix(gene.df[, -1])
  rownames(gene.matrix) <- gene.df$Gene
  
  # 3. Tau (tissue specificity)
  tau <- function(x) {
    if (any(is.na(x))) stop("NA values found.")
    if (any(x < 0))    stop("Negative values found. Maybe data is log-transformed?")
    sum(1 - x / max(x)) / (length(x) - 1)
  }
  tau_scores <- apply(gene.matrix, 1, tau)
  
  # 4. Attach tau to gene_all
  tau_all <- gene_all %>%
    filter(hgnc_symbol %in% names(tau_scores)) %>%
    mutate(tau = tau_scores[hgnc_symbol])
  
  # 5. Group order / colors
  tau_all$group <- factor(tau_all$group, levels = c("fs", "fs_control", "snv", "snv_control"))
  group_colors <- c("fs" = "#1f77b4", "fs_control" = "#aec7e8",
                    "snv" = "#2ca02c", "snv_control" = "#98df8a")
  comparisons <- list(c("fs", "fs_control"), c("snv", "snv_control"))
  
  # 6. Violin plot
  p <- ggplot(tau_all, aes(x = group, y = tau, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") +
    scale_fill_manual(values = group_colors) +
    scale_x_discrete(labels = c("fs" = "fs", "fs_control" = "fs_Control",
                                "snv" = "Nonsense", "snv_control" = "Nonsense_Control")) +
    theme_minimal() +
    labs(title = "Tissue specificity (Tau) across gene categories",
         y = "Tissue Specificity (Tau)", x = "Gene group") +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          legend.position = "none",
          plot.title   = element_text(hjust = 0.5, face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA))
  print(p)
  
  # 7. Save
  write.csv(tau_all, paste0(output_prefix, "_all.csv"), row.names = FALSE)
  ggsave(paste0(output_prefix, "_violin_plot.png"), p, width = 8, height = 6, dpi = 300)
  
  invisible(list(tau_all = tau_all, tau_scores = tau_scores, plot = p))
}


# fn: plot_gene_level_features -------------------------------------------
plot_gene_level_features <- function(
    gene_all,
    lof_metrics_path,
    ensembl     = NULL,
    out_dir     = ".",
    prefix      = "gene_level",
    comparisons = CONFIG$comparisons
) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  group_colors <- c("fs" = "#ff7f0e", "fs_control" = "#ffbb78",
                    "snv" = "#2ca02c", "snv_control" = "#98df8a")
  
  # 1. Load and merge LOF metrics; bin pLI / LOEUF
  Lof_metrics <- read.delim(lof_metrics_path)
  pli_all <- merge(gene_all, Lof_metrics, by.x = "hgnc_symbol", by.y = "gene") %>%
    dplyr::select(hgnc_symbol, pLI, oe_lof_upper, group) %>%
    distinct() %>%
    mutate(
      pli_cat   = cut(pLI,          breaks = c(0, 0.35, 0.66, 1), labels = c("Low", "Medium", "High"), include.lowest = TRUE),
      loeuf_cat = cut(oe_lof_upper, breaks = c(0, 0.2, 0.6, 2),   labels = c("Low", "Medium", "High"), include.lowest = TRUE)
    )
  write_csv(pli_all, file.path(out_dir, paste0(prefix, "_pli_loeuf_category.csv")))
  
  base_box <- function(data, y, title, ylab, log_y = FALSE) {
    g <- ggplot(data, aes(x = group, y = .data[[y]], fill = group)) +
      geom_boxplot(width = 0.1, color = "black") +
      theme_minimal() +
      labs(title = title, y = ylab, x = "Group") +
      theme(axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
            axis.title.x = element_text(size = 12, face = "bold"),
            plot.title   = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none") +
      scale_fill_manual(values = group_colors) +
      stat_compare_means(comparisons = comparisons, method = "wilcox.test")
    if (log_y) g <- g + scale_y_log10()
    g
  }
  
  # 2-4. pLI / LOEUF / CDS length
  pli_plot   <- base_box(pli_all,  "pLI",          "pLI distribution by gene group",   "pLI")
  loeuf_plot <- base_box(pli_all,  "oe_lof_upper", "LOEUF distribution by gene group", "LOEUF")
  cds_plot   <- base_box(gene_all, "cds_length",   "CDS length distribution by gene group", "CDS length, bp", log_y = TRUE)
  
  ggsave(file.path(out_dir, paste0(prefix, "_pli_distribution.pdf")),        plot = pli_plot,   width = 8, height = 5)
  ggsave(file.path(out_dir, paste0(prefix, "_loeuf_distribution.pdf")),      plot = loeuf_plot, width = 8, height = 5)
  ggsave(file.path(out_dir, paste0(prefix, "_cds_length_distribution.pdf")), plot = cds_plot,   width = 8, height = 5)
  
  # 5. Exon count (optional, needs ensembl)
  exon_plot <- NULL; exon_counts <- NULL
  if (!is.null(ensembl)) {
    exons <- getBM(
      attributes = c("ensembl_transcript_id", "cds_start", "cds_end", "rank", "strand"),
      filters    = "ensembl_transcript_id",
      values     = unique(gene_all$ensembl_transcript_id),
      mart       = ensembl
    )
    exon_counts <- exons %>%
      inner_join(gene_all, by = "ensembl_transcript_id") %>%
      group_by(ensembl_transcript_id, hgnc_symbol, strand, group) %>%
      summarise(exon_num = max(rank, na.rm = TRUE), .groups = "drop") %>%
      distinct()
    write_csv(exon_counts, file.path(out_dir, paste0(prefix, "_exon_counts.csv")))
    
    exon_plot <- base_box(exon_counts, "exon_num", "Exon count by gene group", "Exon number", log_y = TRUE)
    ggsave(file.path(out_dir, paste0(prefix, "_exon_count_distribution.pdf")), plot = exon_plot, width = 8, height = 5)
  }
  
  list(pli_all = pli_all, exon_counts = exon_counts,
       pli_plot = pli_plot, loeuf_plot = loeuf_plot,
       cds_plot = cds_plot, exon_plot = exon_plot)
}


# =============================================================================
# 3.0 (cont.) load gene_all
#   Run build_gene_all() once to (re)generate gene_all.csv, then read it back.
#   If gene_all.csv already exists you can skip straight to read.csv().
# =============================================================================

# PTC_info <- read.csv(CONFIG$ptc_info)
# gene_all <- build_gene_all(
#   snv_gene            = read_csv(CONFIG$snv_gene_list, col_names = FALSE),
#   fs_gene             = read_csv(CONFIG$fs_gene_list),
#   snv_control_gene_AD = read_csv(CONFIG$snv_control_list),
#   fs_control_gene_AD  = read_csv(CONFIG$fs_control_list),
#   PTC_info            = PTC_info,
#   mart                = ensembl,
#   omim_ad_symbols_path = CONFIG$omim_ad_symbols,
#   output_csv          = "gene_all.csv"
# )

gene_all <- read.csv("gene_all.csv")


# =============================================================================
# 3.1 add gene level annotation
#   Each call writes its own CSV/figure; those CSVs are re-read in 3.2.
# =============================================================================

calculate_ppi_degree_centrality(
  gene_all,
  output_csv = "wald_ppi_degree_centrality_results.csv"
)

plot_gc_content(
  gene_all   = gene_all,
  output_csv = "gc_content.csv",
  output_fig = "gc_content.png"
)

plot_repeat_content(
  gene_all   = gene_all,
  output_csv = "repeat_content.csv",
  output_fig = "repeat_content.png"
)

annotate_motif_flags(
  gene_all         = gene_all,
  path_touni       = CONFIG$path_touni,
  path_motif       = CONFIG$path_motif,
  path_LCS         = CONFIG$path_lcs,
  mart             = ensembl,
  output_motif_csv = "gene_motif_flags.csv",
  output_lcs_csv   = "gene_LCS_flags.csv"
)

run_pfam_overlap_analysis(
  gene_all      = gene_all,
  ensembl       = ensembl,
  output_prefix = "pfam_overlap"
)

run_ppi_overlap_analysis(
  gene_all      = gene_all,
  ppi_file_path = CONFIG$ppi_file,
  output_prefix = "ppi_overlap"
)

run_tau_analysis(
  gene_all      = gene_all,
  gtex_path     = CONFIG$gtex_path,
  output_prefix = "tau"
)

plot_gene_level_features(
  gene_all         = gene_all,
  lof_metrics_path = CONFIG$lof_metrics_path,
  ensembl          = NULL,
  out_dir          = ".",
  prefix           = "gene_level"
)


# =============================================================================
# 3.2 see if cds length and NMDesc region length could influence the result
# =============================================================================

## 3.2.1 Read the per-feature outputs from 3.1
wald_ppi_degree_centrality <- read_csv("wald_ppi_degree_centrality_results.csv")
gc_content                 <- read_csv("gc_content.csv")
repeat_content             <- read_csv("repeat_content.csv")
gene_motif_flags           <- read_csv("gene_motif_flags.csv")
gene_LCS_flags             <- read_csv("gene_LCS_flags.csv")
pfam_overlap               <- read_csv("pfam_overlap_gene_all.csv")
ppi_overlap                <- read_csv("ppi_overlap_gene_all.csv")
tau_results                <- read_csv("tau_gene_matrix.csv")   # NOTE: expression matrix, not merged below
gene_level                 <- read_csv("gene_level_pli_loeuf_category.csv")

## 3.2.2 Merge everything into one data frame
gene_all_merged <- gene_all %>%
  left_join(wald_ppi_degree_centrality %>% select(hgnc_symbol, group, Degree),
            by = c("hgnc_symbol", "group")) %>%
  left_join(gc_content %>% select(ensembl_transcript_id, group, gc_content, nmdesc_gc_content),
            by = c("ensembl_transcript_id", "group")) %>%
  left_join(repeat_content %>% select(ensembl_transcript_id, group, repeat_fraction, nmdesc_repeat_fraction,
                                      homopolymer_fraction, nmdesc_homopolymer_fraction),
            by = c("ensembl_transcript_id", "group")) %>%
  left_join(gene_motif_flags %>% select(ensembl_transcript_id, group, gene_protein_flag, gene_domains_flag,
                                        gene_slim_flag, gene_morf_flag, gene_ptm_flag, gene_nls_flag),
            by = c("ensembl_transcript_id", "group")) %>%
  left_join(gene_LCS_flags %>% select(ensembl_transcript_id, group, gene_LCS_flag),
            by = c("ensembl_transcript_id", "group")) %>%
  left_join(pfam_overlap %>% select(ensembl_transcript_id, group, pfam_overlap_length,
                                    pfam_overlap_flag, pfam_overlap_fraction, n_overlapping_pfam),
            by = c("ensembl_transcript_id", "group")) %>%
  left_join(ppi_overlap %>% select(ensembl_transcript_id, group, ppi_overlap),
            by = c("ensembl_transcript_id", "group")) %>%
  left_join(gene_level %>% select(hgnc_symbol, group, pLI, oe_lof_upper, pli_cat, loeuf_cat),
            by = c("hgnc_symbol", "group"))

## 3.2.3 Build model_data and save
model_data <- gene_all_merged %>%
  mutate(is_nmdesc = if_else(group %in% c("fs_control", "snv_control"), 0L, 1L))

write.csv(model_data, "gene_model_data.csv", row.names = FALSE)
model_data <- read_csv("gene_model_data.csv")

## 3.2.4 Preliminary regression: PPI degree centrality +/- length covariates
model1_snv <- glm(is_nmdesc ~ Degree,
                  data = model_data[model_data$group %in% c("snv", "snv_control"), ], family = binomial)
model2_snv <- glm(is_nmdesc ~ cds_length + NMDesc_region_length + Degree,
                  data = model_data[model_data$group %in% c("snv", "snv_control"), ], family = binomial)
model1_fs  <- glm(is_nmdesc ~ Degree,
                  data = model_data[model_data$group %in% c("fs", "fs_control"), ], family = binomial)
model2_fs  <- glm(is_nmdesc ~ cds_length + NMDesc_region_length + Degree,
                  data = model_data[model_data$group %in% c("fs", "fs_control"), ], family = binomial)

tidy(model1_snv, exponentiate = TRUE, conf.int = TRUE) %>% arrange(p.value)
tidy(model2_snv, exponentiate = TRUE, conf.int = TRUE) %>% arrange(p.value)
tidy(model1_fs,  exponentiate = TRUE, conf.int = TRUE) %>% arrange(p.value)
tidy(model2_fs,  exponentiate = TRUE, conf.int = TRUE) %>% arrange(p.value)

## 3.2.5 Prepare data for the flag-level association scan
model_data <- model_data %>%
  mutate(
    is_disease = as.integer(group %in% c("snv", "fs")),
    gene_set   = case_when(
      group %in% c("snv", "snv_control") ~ "SNV",
      group %in% c("fs",  "fs_control")  ~ "FS"
    )
  )

flag_cols <- c("gene_protein_flag", "gene_domains_flag", "gene_slim_flag",
               "gene_morf_flag", "gene_ptm_flag", "gene_nls_flag",
               "gene_LCS_flag", "pfam_overlap_flag", "ppi_overlap")

# confounders <- "cds_length + NMDesc_region_length"
confounders <- "cds_length + NMDesc_region_length + gc_content"

## 3.2.6 Drop flags without enough variation in a given gene set
check_flags <- function(data, gene_set_label) {
  data %>%
    filter(gene_set == gene_set_label) %>%
    select(all_of(flag_cols)) %>%
    mutate(across(everything(), as.integer)) %>%
    summarise(across(everything(), ~ n_distinct(na.omit(.)))) %>%
    pivot_longer(everything(), names_to = "flag", values_to = "n_unique") %>%
    filter(n_unique < 2) %>%
    pull(flag)
}

skip_snv <- check_flags(model_data, "SNV")
skip_fs  <- check_flags(model_data, "FS")
message("Skipping in SNV: ", if (length(skip_snv)) paste(skip_snv, collapse = ", ") else "none")
message("Skipping in FS: ",  if (length(skip_fs))  paste(skip_fs,  collapse = ", ") else "none")

## 3.2.7 Odds-ratio helpers (unadjusted + adjusted logistic models per flag)
extract_or <- function(fit, coef_name) {
  ci <- tryCatch(confint(fit)[coef_name, ], error = function(e) c(NA, NA))
  data.frame(
    flag    = coef_name,
    OR      = exp(coef(fit)[coef_name]),
    OR_low  = exp(ci[1]),
    OR_high = exp(ci[2]),
    p_value = summary(fit)$coefficients[coef_name, "Pr(>|z|)"],
    row.names = NULL
  )
}

run_unadjusted <- function(data, gene_set_label, skip_flags = character(0)) {
  active_flags <- setdiff(flag_cols, skip_flags)
  lapply(active_flags, function(col) {
    df <- data %>%
      select(is_disease, all_of(col)) %>%
      filter(!is.na(.data[[col]])) %>%
      mutate(across(all_of(col), as.integer))
    if (nrow(df) == 0 || length(unique(df[[col]])) < 2 || length(unique(df$is_disease)) < 2) {
      message(sprintf("Skipping %s in %s: no variation", col, gene_set_label)); return(NULL)
    }
    fit <- tryCatch(glm(as.formula(paste("is_disease ~", col)), data = df, family = binomial),
                    error = function(e) { message(sprintf("Model failed for %s in %s: %s", col, gene_set_label, e$message)); NULL })
    if (is.null(fit)) return(NULL)
    extract_or(fit, col)
  }) %>%
    bind_rows() %>%
    mutate(gene_set = gene_set_label, model = "Unadjusted")
}

run_adjusted <- function(data, gene_set_label, skip_flags = character(0)) {
  active_flags <- setdiff(flag_cols, skip_flags)
  lapply(active_flags, function(col) {
    df <- data %>%
      select(is_disease, all_of(col), cds_length, NMDesc_region_length, gc_content) %>%
      filter(!is.na(.data[[col]]), !is.na(cds_length), !is.na(NMDesc_region_length), !is.na(gc_content)) %>%
      mutate(across(all_of(col), as.integer))
    if (nrow(df) == 0 || length(unique(df[[col]])) < 2 || length(unique(df$is_disease)) < 2) {
      message(sprintf("Skipping %s in %s: no variation", col, gene_set_label)); return(NULL)
    }
    fit <- tryCatch(glm(as.formula(paste("is_disease ~", col, "+", confounders)), data = df, family = binomial),
                    error = function(e) { message(sprintf("Model failed for %s in %s: %s", col, gene_set_label, e$message)); NULL })
    if (is.null(fit)) return(NULL)
    extract_or(fit, col)
  }) %>%
    bind_rows() %>%
    mutate(gene_set = gene_set_label, model = "Adjusted")
}

unadj_snv <- run_unadjusted(filter(model_data, gene_set == "SNV"), "SNV", skip_snv)
unadj_fs  <- run_unadjusted(filter(model_data, gene_set == "FS"),  "FS",  skip_fs)
adj_snv   <- run_adjusted(filter(model_data, gene_set == "SNV"),   "SNV", skip_snv)
adj_fs    <- run_adjusted(filter(model_data, gene_set == "FS"),    "FS",  skip_fs)

## 3.2.8 Combine, BH-correct within (gene_set x model)
combined <- bind_rows(unadj_snv, unadj_fs, adj_snv, adj_fs) %>%
  group_by(gene_set, model) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig   = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  ) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = c("Unadjusted", "Adjusted")))

## 3.2.9 Forest plot of feature -> disease associations
forest_plot <- ggplot(combined, aes(x = OR, y = reorder(flag, OR), color = sig, shape = gene_set)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.2,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  facet_wrap(~ model) +
  scale_color_manual(values = c("***" = "red", "**" = "orange", "*" = "gold", "ns" = "grey60")) +
  scale_shape_manual(values = c("SNV" = 16, "FS" = 17)) +
  scale_x_log10() +
  labs(
    title    = "Gene-level feature association with disease (model_data)",
    subtitle = "Adjusted for cds_length + NMDesc_region_length + GC content",
    x = "Odds Ratio (log scale)", y = NULL,
    color = "Adj. p", shape = "Gene set"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(forest_plot)
ggsave("gene_level_feature_forest.png", forest_plot, width = 10, height = 6, dpi = 300)


# =============================================================================
# 3.3 CDS length and NMDesc region length matched analysis
#   Study genes are matched 1:1 to controls on CDS length + NMDesc region
#   length (pairs supplied in matching_log). For each gene-level feature flag,
#   the study vs. matched-control difference is tested WITHIN pairs using
#   McNemar's test (exact binomial when there are few discordant pairs),
#   then summarised as a 4-panel bar figure + a summary table.
# =============================================================================

# --- 3.3.0 build matching_log from gene_all -------------------------------
## For each disease gene (snv / fs) find the single nearest control gene
## (snv_control / fs_control) on CDS length + NMDesc region length. Distance is
## Euclidean on z-scored variables (pooled study+control mean/sd), so both
## variables contribute on a common scale. Matching is 1:1 greedy WITHOUT
## replacement by default (each control used once) — appropriate for the paired
## McNemar design in 3.3.1. Set with_replacement = TRUE to allow reuse, or a
## finite caliper to forbid poor matches.
build_matching_log <- function(
    gene_all,
    output_csv = "matching_log.csv",
    strata = list(
      list(study = "snv", control = "snv_control", label = "STOPGAIN"),
      list(study = "fs",  control = "fs_control",  label = "FS")
    ),
    match_vars       = c("cds_length", "NMDesc_region_length"),
    with_replacement = FALSE,
    caliper          = Inf,     # max z-scored distance allowed (Inf = no limit)
    seed             = 1
) {
  set.seed(seed)
  
  # one row per gene per group, with the matching variables (complete cases)
  pool <- gene_all %>%
    select(hgnc_symbol, group, all_of(match_vars)) %>%
    filter(if_all(all_of(match_vars), ~ !is.na(.))) %>%
    group_by(hgnc_symbol, group) %>%
    summarise(across(all_of(match_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  all_pairs <- list()
  
  for (st in strata) {
    study_df <- pool %>% filter(group == st$study)   %>% distinct(hgnc_symbol, .keep_all = TRUE)
    ctrl_df  <- pool %>% filter(group == st$control) %>% distinct(hgnc_symbol, .keep_all = TRUE)
    
    if (nrow(study_df) == 0 || nrow(ctrl_df) == 0) {
      message(sprintf("Stratum %s: missing %s genes — skipped",
                      st$label, if (nrow(study_df) == 0) "study" else "control"))
      next
    }
    
    # z-score each variable on the pooled (study + control) distribution
    pooled <- bind_rows(study_df, ctrl_df)
    ctr <- vapply(match_vars, function(v) mean(pooled[[v]]), numeric(1))
    scl <- vapply(match_vars, function(v) {
      s <- sd(pooled[[v]]); if (is.na(s) || s == 0) 1 else s
    }, numeric(1))
    z <- function(df) sweep(sweep(as.matrix(df[, match_vars, drop = FALSE]),
                                  2, ctr, "-"), 2, scl, "/")
    study_z <- z(study_df)
    ctrl_z  <- z(ctrl_df)
    
    available <- rep(TRUE, nrow(ctrl_df))
    n_matched <- 0L
    
    for (i in seq_len(nrow(study_df))) {
      d <- sqrt(colSums((t(ctrl_z) - study_z[i, ])^2))   # dist to every control
      if (!with_replacement) d[!available] <- Inf
      j <- which.min(d)
      if (length(j) == 0 || !is.finite(d[j]) || d[j] > caliper) next  # no match
      if (!with_replacement) available[j] <- FALSE
      n_matched <- n_matched + 1L
      
      row <- data.frame(
        study_gene = study_df$hgnc_symbol[i],
        ctrl_gene  = ctrl_df$hgnc_symbol[j],
        group      = st$label,
        match_dist = d[j],
        stringsAsFactors = FALSE
      )
      for (v in match_vars) {
        row[[paste0("study_", v)]] <- study_df[[v]][i]
        row[[paste0("ctrl_",  v)]] <- ctrl_df[[v]][j]
      }
      all_pairs[[length(all_pairs) + 1]] <- row
    }
    message(sprintf("Stratum %s: matched %d / %d study genes",
                    st$label, n_matched, nrow(study_df)))
  }
  
  matching_log <- bind_rows(all_pairs)
  write_csv(matching_log, output_csv)
  message("Saved: ", output_csv, " (", nrow(matching_log), " matched pairs)")
  invisible(matching_log)
}

build_matching_log(gene_all = gene_all, output_csv = "matching_log.csv")


# --- 3.3.1 matched McNemar / exact-test analysis --------------------------
run_matched_length_analysis <- function(
    gene_all,
    matching_log_path,
    flag_cols = c("ppi_overlap", "pfam_overlap_flag", "gene_nls_flag",
                  "gene_ptm_flag", "gene_slim_flag", "gene_morf_flag",
                  "gene_protein_flag", "gene_domains_flag", "gene_LCS_flag"),
    output_prefix = "matched_length_enrichment",
    exact_cutoff  = 25
) {
  
  # ---- nested helpers ------------------------------------------------------
  coerce_binary <- function(x) {
    if (is.logical(x))   return(as.numeric(x))
    if (is.character(x)) return(case_when(
      x %in% c("TRUE", "True", "true", "1", "yes", "Yes") ~ 1,
      x %in% c("FALSE", "False", "false", "0", "no", "No") ~ 0,
      TRUE ~ suppressWarnings(as.numeric(x))
    ))
    suppressWarnings(as.numeric(x))
  }
  
  format_p <- function(p) {
    if (is.na(p))     return("p = NA")
    if (p < 0.001)    return("p < 0.001")
    paste0("p = ", sprintf("%.3f", p))
  }
  
  # Positive counts / percentages for a feature, per match type + its control
  summarise_matched_feature <- function(data, feature, match_type_now) {
    study_col <- paste0(feature, "_study")
    ctrl_col  <- paste0(feature, "_ctrl")
    d <- data %>%
      filter(match_type == match_type_now,
             !is.na(.data[[study_col]]), !is.na(.data[[ctrl_col]])) %>%
      mutate(study = ifelse(.data[[study_col]] >= 1, 1, 0),
             ctrl  = ifelse(.data[[ctrl_col]]  >= 1, 1, 0))
    ctrl_label <- ifelse(match_type_now == "Frameshift", "FS Control", "Nonsense Control")
    tibble(
      group_label = c(match_type_now, ctrl_label),
      positive_n  = c(sum(d$study == 1), sum(d$ctrl == 1)),
      total_n     = c(nrow(d), nrow(d)),
      percent     = c(mean(d$study == 1) * 100, mean(d$ctrl == 1) * 100)
    )
  }
  
  # Within-pair test: exact binomial on discordant pairs, else McNemar
  get_mcnemar_p <- function(data, feature, match_type_now) {
    study_col <- paste0(feature, "_study")
    ctrl_col  <- paste0(feature, "_ctrl")
    d <- data %>%
      filter(match_type == match_type_now,
             !is.na(.data[[study_col]]), !is.na(.data[[ctrl_col]])) %>%
      mutate(study = ifelse(.data[[study_col]] >= 1, 1, 0),
             ctrl  = ifelse(.data[[ctrl_col]]  >= 1, 1, 0))
    if (nrow(d) == 0) return(NA_real_)
    
    tab <- table(factor(d$study, levels = c(0, 1)),
                 factor(d$ctrl,  levels = c(0, 1)))
    b <- unname(tab["0", "1"])   # control-only positive
    c <- unname(tab["1", "0"])   # study-only positive
    discordant_n <- b + c
    if (discordant_n == 0) return(NA_real_)
    
    tryCatch(
      if (discordant_n < exact_cutoff) {
        binom.test(x = c, n = discordant_n, p = 0.5, alternative = "two.sided")$p.value
      } else {
        mcnemar.test(tab, correct = TRUE)$p.value
      },
      error = function(e) NA_real_
    )
  }
  
  # One paired bar panel for a single feature
  make_matched_plot <- function(data, feature, title, panel_letter) {
    summary_df <- bind_rows(
      summarise_matched_feature(data, feature, "Frameshift"),
      summarise_matched_feature(data, feature, "Nonsense")
    ) %>%
      mutate(group_label = factor(
        group_label,
        levels = c("Frameshift", "FS Control", "Nonsense", "Nonsense Control")
      ))
    
    p_fs  <- get_mcnemar_p(data, feature, "Frameshift")
    p_snv <- get_mcnemar_p(data, feature, "Nonsense")
    
    ymax     <- max(summary_df$percent, na.rm = TRUE)
    ylim_top <- max(100, ymax + 18)
    
    group_colors <- c(
      "Frameshift"       = "#D07A3A",
      "FS Control"       = "#E2C7B5",
      "Nonsense"         = "#8ABA67",
      "Nonsense Control" = "#C9D6BE"
    )
    
    ggplot(summary_df, aes(x = group_label, y = percent, fill = group_label)) +
      geom_col(width = 0.62, color = "grey40", alpha = 0.95) +
      geom_text(aes(label = sprintf("%.1f%%", percent)),
                vjust = -0.55, size = 3.6, fontface = "bold") +
      geom_label(aes(y = pmax(percent * 0.5, 6),
                     label = paste0("n=", positive_n, "/", total_n)),
                 size = 3.0, label.size = 0, fill = "white", alpha = 0.8) +
      annotate("segment", x = 1, xend = 2, y = ymax + 8,   yend = ymax + 8,   linewidth = 0.6) +
      annotate("segment", x = 1, xend = 1, y = ymax + 6.5, yend = ymax + 8,   linewidth = 0.6) +
      annotate("segment", x = 2, xend = 2, y = ymax + 6.5, yend = ymax + 8,   linewidth = 0.6) +
      annotate("text",    x = 1.5, y = ymax + 10.5, label = format_p(p_fs),  fontface = "bold", size = 3.4) +
      annotate("segment", x = 3, xend = 4, y = ymax + 8,   yend = ymax + 8,   linewidth = 0.6) +
      annotate("segment", x = 3, xend = 3, y = ymax + 6.5, yend = ymax + 8,   linewidth = 0.6) +
      annotate("segment", x = 4, xend = 4, y = ymax + 6.5, yend = ymax + 8,   linewidth = 0.6) +
      annotate("text",    x = 3.5, y = ymax + 10.5, label = format_p(p_snv), fontface = "bold", size = 3.4) +
      scale_fill_manual(values = group_colors) +
      scale_y_continuous(limits = c(0, ylim_top), expand = expansion(mult = c(0, 0.02))) +
      labs(title = title, x = "Matched Gene Category", y = "Percentage (%)") +
      theme_bw(base_size = 13) +
      theme(
        legend.position    = "none",
        plot.title         = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.title.x       = element_text(face = "bold"),
        axis.title.y       = element_text(face = "bold"),
        axis.text.x        = element_text(size = 11),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin        = margin(12, 12, 12, 12)
      ) +
      annotate("text", x = -Inf, y = Inf, label = panel_letter,
               hjust = -0.8, vjust = 1.5, size = 20 / .pt, fontface = "bold")
  }
  
  # ---- 1. read matching log, standardise match type ------------------------
  match_log <- read_csv(matching_log_path, show_col_types = FALSE)
  miss_cols <- setdiff(c("study_gene", "ctrl_gene", "group"), colnames(match_log))
  if (length(miss_cols) > 0)
    stop("matching_log is missing columns: ", paste(miss_cols, collapse = ", "))
  
  match_log <- match_log %>%
    mutate(
      group      = toupper(group),
      match_type = case_when(
        group == "FS"       ~ "Frameshift",
        group == "STOPGAIN" ~ "Nonsense",
        TRUE                ~ NA_character_
      )
    ) %>%
    filter(!is.na(match_type))
  
  # ---- 2. keep flags present in gene_all, coerce to 0/1 --------------------
  flag_cols <- intersect(flag_cols, colnames(gene_all))
  if (length(flag_cols) == 0)
    stop("None of the requested flag columns were found in gene_all")
  for (cc in flag_cols) gene_all[[cc]] <- coerce_binary(gene_all[[cc]])
  
  # ---- 3. collapse to one row per gene (max across rows for each flag) -----
  gene_flags <- gene_all %>%
    select(hgnc_symbol, all_of(flag_cols)) %>%
    group_by(hgnc_symbol) %>%
    summarise(across(all_of(flag_cols), ~ {
      x <- suppressWarnings(as.numeric(.x))
      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    }), .groups = "drop")
  
  # ---- 4. attach flags to the study and control gene of each pair ----------
  study_flags <- gene_flags %>%
    rename(study_gene = hgnc_symbol) %>%
    rename_with(~ paste0(.x, "_study"), -study_gene)
  ctrl_flags <- gene_flags %>%
    rename(ctrl_gene = hgnc_symbol) %>%
    rename_with(~ paste0(.x, "_ctrl"), -ctrl_gene)
  
  matched_df <- match_log %>%
    left_join(study_flags, by = "study_gene") %>%
    left_join(ctrl_flags,  by = "ctrl_gene")
  
  # ---- 5. features to plot (4 panels) --------------------------------------
  feature_map <- tibble(
    feature = c("ppi_overlap", "pfam_overlap_flag", "gene_ptm_flag", "gene_LCS_flag"),
    title   = c("PPI Residue Enrichment", "Pfam Domain Enrichment",
                "Post-Translational Modification Enrichment",
                "Low Complexity Sequence Enrichment"),
    panel   = c("A", "B", "C", "D")
  ) %>%
    filter(feature %in% flag_cols)
  if (nrow(feature_map) == 0)
    stop("No features in feature_map were found in gene_all")
  
  # ---- 6. build panels + combined figure -----------------------------------
  plot_list <- lapply(seq_len(nrow(feature_map)), function(i)
    make_matched_plot(matched_df, feature_map$feature[i],
                      feature_map$title[i], feature_map$panel[i]))
  names(plot_list) <- feature_map$panel
  
  combined_fig <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title    = "Matched gene-level feature enrichment",
      subtitle = "Study genes matched to controls on CDS length + NMDesc region length",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 20, hjust = 0),
        plot.subtitle = element_text(size = 12, hjust = 0)
      )
    )
  
  # ---- 7. summary table + within-pair p-values -----------------------------
  summary_table <- bind_rows(lapply(feature_map$feature, function(feat) {
    this_title <- feature_map$title[match(feat, feature_map$feature)]
    bind_rows(
      summarise_matched_feature(matched_df, feat, "Frameshift") %>%
        mutate(comparison = "FS vs matched control"),
      summarise_matched_feature(matched_df, feat, "Nonsense") %>%
        mutate(comparison = "STOPGAIN vs matched control")
    ) %>%
      mutate(feature = feat, title = this_title, .before = 1)
  }))
  
  p_table <- tibble(
    feature    = feature_map$feature,
    title      = feature_map$title,
    p_fs       = vapply(feature_map$feature,
                        function(f) get_mcnemar_p(matched_df, f, "Frameshift"), numeric(1)),
    p_stopgain = vapply(feature_map$feature,
                        function(f) get_mcnemar_p(matched_df, f, "Nonsense"), numeric(1))
  )
  summary_table <- summary_table %>% left_join(p_table, by = c("feature", "title"))
  
  # ---- 8. save outputs -----------------------------------------------------
  ggsave(paste0(output_prefix, "_4panels.pdf"), combined_fig, width = 12, height = 10, dpi = 300)
  ggsave(paste0(output_prefix, "_4panels.png"), combined_fig, width = 12, height = 10, dpi = 300)
  write_csv(summary_table, paste0(output_prefix, "_summary.csv"))
  message("Saved: ", output_prefix, "_4panels.pdf / .png / _summary.csv")
  
  invisible(list(matched_df = matched_df, summary_table = summary_table, figure = combined_fig))
}

matched_length_results <- run_matched_length_analysis(
  gene_all          = model_data,          # annotated table from 3.2, has the flags
  matching_log_path = "matching_log.csv",
  output_prefix     = "matched_length_enrichment"
)
