library(stringr)
library(dplyr)
library(readxl)
library(readr)
library(biomaRt)

#' Add motif and LCS flags to a variants data frame
#'
#' Mirrors the original script exactly:
#'   1. Uses biomaRt to map ensembl_transcript_id -> uniprotswissprot
#'   2. Loads touni, motif, and LCS supplement files
#'   3. Reduces multi-value cells to their maximum
#'   4. Merges with variants and computes per-variant boolean flags
#'
#' @param variants_all2  Data frame of variants; must contain columns
#'                       `ensembl_transcript_id` and `cds_mutation_loc`.
#' @param mart           A biomaRt Mart object (from useMart() or useEnsembl()).
#' @param touni_path     Path to supplement A CSV  (protein -> UniProt mapping).
#' @param motif_path     Path to supplement B CSV  (motif feature data).
#' @param lcs_path       Path to supplement Excel  (LCS data).
#' @param lcs_sheet      Sheet name or index for LCS data (default: "G").
#'
#' @return `variants_all2` with seven additional logical flag columns:
#'   `variant_protein_flag`, `variant_domain_flag`, `variant_slim_flag`,
#'   `variant_morf_flag`, `variant_ptm_flag`, `variant_nls_flag`,
#'   `variant_LCS_flag`.
variants_motif <- function(variants_all2,
                           mart,
                           touni_path,
                           motif_path,
                           lcs_path,
                           lcs_sheet = "G") {
  
  # ── helper: largest integer in a cell string ───────────────────────────────
  get_max_from_cell <- function(x) {
    if (is.na(x) || x == "") return(NA_real_)
    nums <- str_extract_all(x, "\\d+")[[1]]
    if (length(nums) == 0) return(NA_real_)
    max(as.numeric(nums), na.rm = TRUE)
  }
  
  # ── step 1: map transcript IDs -> uniprotswissprot (mirrors original getBM) ─
  tx_ids <- unique(variants_all2$ensembl_transcript_id)
  
  bm_raw <- getBM(
    attributes = c("ensembl_transcript_id", "uniprotswissprot"),
    filters    = "ensembl_transcript_id",
    values     = tx_ids,
    mart       = mart
  )
  
  variants_all2$uniprotswissprot <- bm_raw$uniprotswissprot[
    match(variants_all2$ensembl_transcript_id, bm_raw$ensembl_transcript_id)
  ]
  
  # ── step 2: load supplement files (mirrors original read_csv / read_excel) ──
  touni     <- read_csv(touni_path,  show_col_types = FALSE)
  motif_doc <- read_csv(motif_path,  show_col_types = FALSE)
  LCS_doc   <- read_excel(lcs_path,  sheet = lcs_sheet)
  
  # ── step 3: reduce cells to max (mirrors original mutate/across block) ──────
  motif_max <- motif_doc %>%
    mutate(across(-1, ~ vapply(as.character(.x), get_max_from_cell, numeric(1))))
  
  LCS_max <- LCS_doc %>%
    mutate(across(-1, ~ vapply(as.character(.x), get_max_from_cell, numeric(1))))
  
  # ── step 4: attach uniprot IDs to motif/LCS tables ────────────────────────
  motif_max$uniprot <- touni$`Uniprot ID`[match(motif_max$Protein, touni$Protein)]
  LCS_max$uniprot   <- touni$`Uniprot ID`[match(LCS_max$Protein,   touni$Protein)]
  
  # ── step 5: add uniprot alias to variants (mirrors variants_all$uniprot = ...) 
  variants_all2$uniprot <- variants_all2$uniprotswissprot
  
  # ── step 6: merge and compute motif flags (mirrors original merge + flag lines)
  motif_max3 <- merge(motif_max, variants_all2, by = "uniprot")
  
  motif_max3$variant_protein_feature_flag <- motif_max3$`Protein Features` * 3 >= motif_max3$cds_mutation_loc
  motif_max3$variant_domains_flag         <- motif_max3$`Domains`           * 3 >= motif_max3$cds_mutation_loc
  motif_max3$variant_slim_flag            <- motif_max3$`SLiMs`             * 3 >= motif_max3$cds_mutation_loc
  motif_max3$variant_morf_flag            <- motif_max3$MORFs               * 3 >= motif_max3$cds_mutation_loc
  motif_max3$variant_ptm_flag             <- motif_max3$`PTMs`              * 3 >= motif_max3$cds_mutation_loc
  motif_max3$variant_nls_flag             <- motif_max3$`NLSs`              * 3 >= motif_max3$cds_mutation_loc
  
  # ── step 7: merge and compute LCS flag ────────────────────────────────────
  LCS_max3 <- merge(LCS_max, variants_all2, by = "uniprot")
  LCS_max3$variant_LCS_flag <- LCS_max3$`LCSs` * 3 >= LCS_max3$cds_mutation_loc
  
  # ── step 8: map flags back to variants_all2 (mirrors original match lines) ─
  variants_all2$variant_protein_flag <- motif_max3$variant_protein_feature_flag[match(variants_all2$uniprot, motif_max3$uniprot)]
  variants_all2$variant_domain_flag  <- motif_max3$variant_domains_flag[match(variants_all2$uniprot, motif_max3$uniprot)]
  variants_all2$variant_slim_flag    <- motif_max3$variant_slim_flag[match(variants_all2$uniprot, motif_max3$uniprot)]
  variants_all2$variant_morf_flag    <- motif_max3$variant_morf_flag[match(variants_all2$uniprot, motif_max3$uniprot)]
  variants_all2$variant_ptm_flag     <- motif_max3$variant_ptm_flag[match(variants_all2$uniprot, motif_max3$uniprot)]
  variants_all2$variant_nls_flag     <- motif_max3$variant_nls_flag[match(variants_all2$uniprot, motif_max3$uniprot)]
  variants_all2$variant_LCS_flag     <- LCS_max3$variant_LCS_flag[match(variants_all2$uniprot, LCS_max3$uniprot)]
  
  variants_all2
}