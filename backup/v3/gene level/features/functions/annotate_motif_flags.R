annotate_motif_flags <- function(
    gene_all,
    path_touni,
    path_motif,
    path_LCS,
    mart,
    output_motif_csv = "gene_motif_flags.csv",
    output_lcs_csv   = "gene_LCS_flags.csv"
) {
  library(stringr)
  library(dplyr)
  library(readxl)
  library(readr)
  
  # ── Helpers ───────────────────────────────────────────────────────────────
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
  
  # ── Load external data ────────────────────────────────────────────────────
  touni     <- read_csv(path_touni, show_col_types = FALSE)
  motif_doc <- read_csv(path_motif, show_col_types = FALSE)
  LCS_doc   <- read_excel(path_LCS, sheet = "G")
  
  motif_max <- motif_doc %>% parse_max() %>% add_uniprot(touni)
  LCS_max   <- LCS_doc   %>% parse_max() %>% add_uniprot(touni)
  
  # ── Fix 1: replace empty uniprot strings with NA in gene_all ─────────────
  gene_all <- gene_all %>%
    mutate(uniprot = na_if(trimws(uniprot), ""))
  
  # ── Fix 2: re-fetch uniprot for genes where it's missing ─────────────────
  missing_genes <- gene_all %>%
    filter(is.na(uniprot)) %>%
    pull(hgnc_symbol) %>%
    unique()
  
  if (length(missing_genes) > 0) {
    message(length(missing_genes), " genes missing UniProt ID — re-fetching from BioMart")
    uni_map <- getBM(
      attributes = c("hgnc_symbol", "uniprotswissprot"),
      filters    = "hgnc_symbol",
      values     = missing_genes,
      mart       = mart
    ) %>%
      filter(uniprotswissprot != "") %>%          # drop empty BioMart returns
      distinct(hgnc_symbol, .keep_all = TRUE)     # one UniProt per gene
    
    gene_all <- gene_all %>%
      left_join(uni_map, by = "hgnc_symbol") %>%
      mutate(uniprot = coalesce(uniprot, uniprotswissprot)) %>%
      select(-uniprotswissprot)
  }
  
  message("Genes with UniProt after re-fetch: ",
          sum(!is.na(gene_all$uniprot)), " / ", nrow(gene_all))
  
  # ── Motif flags ───────────────────────────────────────────────────────────
  motif_merged <- merge(motif_max, gene_all, by = "uniprot") %>%
    flag_cols(list(
      gene_protein_flag = "Protein Features",
      gene_domains_flag = "Domains",
      gene_slim_flag    = "SLiMs",
      gene_morf_flag    = "MORFs",
      gene_ptm_flag     = "PTMs",
      gene_nls_flag     = "NLSs/NESs"       # Fix 1: corrected column name
    ))
  
  # ── LCS flags ─────────────────────────────────────────────────────────────
  lcs_merged <- merge(LCS_max, gene_all, by = "uniprot") %>%
    flag_cols(list(gene_LCS_flag = "LCSs"))
  
  # ── Merge flags back to gene_all ──────────────────────────────────────────
  motif_flags <- c("gene_protein_flag", "gene_domains_flag", "gene_slim_flag",
                   "gene_morf_flag",    "gene_ptm_flag",     "gene_nls_flag")
  for (flag in motif_flags) {
    gene_all[[flag]] <- motif_merged[[flag]][match(gene_all$uniprot, motif_merged$uniprot)]
  }
  gene_all$gene_LCS_flag <- lcs_merged$gene_LCS_flag[match(gene_all$uniprot, lcs_merged$uniprot)]
  
  # ── Save ──────────────────────────────────────────────────────────────────
  write.csv(motif_merged, output_motif_csv, row.names = FALSE)
  write.csv(lcs_merged,   output_lcs_csv,   row.names = FALSE)
  message("Motif flags saved to: ", output_motif_csv)
  message("LCS flags saved to:   ", output_lcs_csv)
  
  invisible(list(gene_all = gene_all, motif_merged = motif_merged, lcs_merged = lcs_merged))
}
