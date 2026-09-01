# ==============================================================================
# 4. VARIANT-LEVEL COMPARISON
# ------------------------------------------------------------------------------
#   4.0  Config, data loading, annotation pipeline, and shared helper functions
#   4.1  Unmatched analysis
#   4.2  Mixed-effect model
#   4.3  Hierarchical Bayesian model
# ==============================================================================


# ==============================================================================
# 4.0  CONFIG, DATA, AND SHARED HELPER FUNCTIONS
# ==============================================================================

library(tidyverse)
library(biomaRt)
library(broom)
library(stringr)
library(dplyr)
library(readxl)
library(readr)
library(ggplot2)
library(scales)
library(ggpubr)

# ------------------------------------------------------------------------------
# 4.0.1  Paths and constants
# ------------------------------------------------------------------------------

.p <- c("gene level_v3/lib/paths.R", "../gene level_v3/lib/paths.R",
        "../../gene level_v3/lib/paths.R")
.p <- .p[file.exists(.p)]
source(.p[1]); rm(.p)

CLINVAR_DIR <- data_root("clinvar")
OTHERS_DIR  <- file.path(CLINVAR_DIR, "others")
DATA_DIR    <- dirname(CLINVAR_DIR)
OUT_DIR     <- out_dir("plots/pfam_ppi_analysis")
# -----------------------------------------------------------------------------

FDR_METHOD <- "BH"  
FLAG_COLS_ALL <- c(
  "variant_ppi_overlap", "ptc_after_max_pfam_end", "ptc_before_max_pfam_end",
  "variant_protein_flag", "variant_domain_flag", "variant_slim_flag",
  "variant_morf_flag", "variant_ptm_flag", "variant_nls_flag", "variant_LCS_flag"
)

# One flag per feature: ptc_after and ptc_before are complementary sides of the
# same continuous variable, so only ptc_before enters the analysis.
FLAG_COLS <- setdiff(FLAG_COLS_ALL, "ptc_after_max_pfam_end")

CONT_COL         <- "dist_to_cds_end_log"
CONFOUNDERS      <- c("cds_length")
CONFOUNDERS_FULL <- c("cds_length", "NMDesc_region_length", "GC_Content")

# Flags for the gene-matched panel: one Pfam flag, one PPI flag, plus motifs
CONT_PREDICTORS <- c("dist_to_cds_end_log")

GENE_MATCHED_FLAGS <- c(
  "variant_ppi_overlap",       # PPI: the single overlap flag
  "ptc_before_max_pfam_end",   # Pfam: keep the "domain disrupted" direction
  "variant_protein_flag", "variant_domain_flag", "variant_slim_flag",
  "variant_morf_flag", "variant_ptm_flag", "variant_nls_flag", "variant_LCS_flag",
  "dist_to_cds_end_log"        
)

# Publication-ready display names
FLAG_LABELS <- c(
  variant_ppi_overlap     = "PPI interface overlap",
  ptc_before_max_pfam_end = "Pfam domain disruption",
  variant_protein_flag    = "Protein feature",
  variant_domain_flag     = "Structured domain",
  variant_slim_flag       = "Short linear motif (SLiM)",
  variant_morf_flag       = "Molecular recognition feature (MoRF)",
  variant_ptm_flag        = "Post-translational modification (PTM)",
  variant_nls_flag        = "Nuclear localization signal (NLS)",
  variant_LCS_flag        = "Low-complexity segment",
  dist_to_cds_end_log     = "Distance to CDS end (log)"
)

setwd(CLINVAR_DIR)

#4.0 load data and add annotations
#4.0.1 adds cds_mutation_loc, dist_to_cds_end, etc
source("new_create_fasta_functions.R")
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
snv_variants = read.csv('snv_variants20260201_plp_dbh_clinvar.csv')
snv_dis <- create_fasta(snv_variants, output_dir = "snv_disease_fasta_output")

snv_control_variants = read.csv('gnomad_snv_filtered_acat_0831.csv')
gnomad_snv_filtered <- snv_control_variants
gnomad_snv_filtered$key <- gnomad_snv_filtered$id
gnomad_snv_variants <- gnomad_snv_filtered[,  c("transcript", "key")]
gnomad_snv_dis <- create_fasta(gnomad_snv_variants, output_dir = "snv_control_fasta_output")

fs_variants = read.csv('fs_variants20260201_plp_acat_clinvar.csv')
fs_dis <- create_fasta(fs_variants, output_dir = "fs_disease_fasta_output")

fs_control_variants = read.csv('gnomad_fs_filtered_bh_0831.csv')
gnomad_fs_filtered  <- fs_control_variants
gnomad_fs_variants  <- gnomad_fs_filtered[,  c("transcript", "key")]
gnomad_fs_dis <- create_fasta(gnomad_fs_variants, output_dir = "fs_control_fasta_output")

variants_all1 <- bind_rows(
  fs_dis %>% mutate(group = "fs_disease"),
  snv_dis %>% mutate(group = "snv_disease"),
  gnomad_fs_dis %>% mutate(group = "fs_control"),
  gnomad_snv_dis %>% mutate(group = "snv_control")
)
variants_all1$ensembl_transcript_id = sapply(strsplit(variants_all1$Variant_Key, "\\|"), function(x) x[1])
#add NMDesc region, coding sequence info from gene_all
all_variants = bind_rows(fs_variants, snv_variants, gnomad_fs_variants, gnomad_snv_variants)
#get transcript id from all_variants
variants_all1$ensembl_transcript_id = all_variants$transcript[match(variants_all1$Variant_Key, all_variants$key)]

#4.0.2 merges variant-level info with gene-level info
#run combine_gene.R get snv_nmdesc_df
snv_all <- snv_variants %>%
  left_join(
    snv_dis %>% #snv_dis is the output of create_fasta
      dplyr::select(Variant_Key, ptc_pos, cds_mutation_loc),
    by = c("key" = "Variant_Key")
  ) %>%
  left_join(
    snv_nmdesc_df %>% #snv_nmdesc_df contains gene level info
      dplyr::select(
        ensembl_transcript_id,
        can_region_start,
        can_region_end,
        snv_nmdesc_cds,
        coding
      ),
    by = c("transcript" = "ensembl_transcript_id")
  ) %>%
  #can_region_start and can_region_end as NMD_region_start and NMD_region_end
  rename(
    NMD_region_start = can_region_start,
    NMD_region_end = can_region_end,
    cds_ptc_loc = ptc_pos,
    nmdesc_cds = snv_nmdesc_cds,
  ) %>%
  dplyr::select(-any_of(c("genomic_loc", "PTC_pos"))) # remove the original genomic_pos column since we have renamed it to PTC_genomic_pos

gnomad_snv_all <- gnomad_snv_filtered %>%
  left_join(
    gnomad_snv_dis %>%
      dplyr::select(Variant_Key, ptc_pos, cds_mutation_loc),
    by = c("key" = "Variant_Key")
  ) %>%
  left_join(
    snv_nmdesc_df %>%
      dplyr::select(
        ensembl_transcript_id,
        can_region_start,
        can_region_end,
        snv_nmdesc_cds,
        coding
      ),
    by = c("transcript" = "ensembl_transcript_id")
  ) %>%
  rename(
    NMD_region_start = can_region_start,
    NMD_region_end = can_region_end,
    cds_ptc_loc = ptc_pos,
  ) %>%
  dplyr::select(-any_of(c("genomic_loc", "PTC_pos"))) # remove the original genomic_pos column since we have renamed it to PTC_genomic_pos

gnomad_fs_all <- gnomad_fs_filtered2 %>%
  left_join(
    gnomad_fs_dis %>%
      dplyr::select(Variant_Key, ptc_pos, cds_mutation_loc),
    by = c("key" = "Variant_Key")
  ) %>%
  left_join(
    fs_nmdesc_df %>%
      dplyr::select(
        ensembl_transcript_id,
        median_can_region_start,
        median_can_region_end,
        fs_nmdesc_cds,
        coding
      ),
    by = c("transcript" = "ensembl_transcript_id")
  ) %>%
  # rename for median_can_region_start and median_can_region_end as NMD_region_start and NMD_region_end
  rename(
    NMD_region_start = median_can_region_start,
    NMD_region_end = median_can_region_end,
    cds_ptc_loc = ptc_pos,
  ) %>%
  dplyr::select(-any_of(c("genomic_loc", "PTC_pos"))) # remove the original genomic_pos column since we have renamed it to PTC_genomic_pos

fs_all <- fs_variants %>%
  left_join(
    fs_dis %>%
      dplyr::select(Variant_Key, ptc_pos, cds_mutation_loc),
    by = c("key" = "Variant_Key")
  ) %>%
  left_join(
    fs_nmdesc_df %>%
      dplyr::select(
        ensembl_transcript_id,
        median_can_region_start,
        median_can_region_end,
        fs_nmdesc_cds,
        coding
      ),
    by = c("transcript" = "ensembl_transcript_id")
  ) %>%
  rename(
    NMD_region_start = median_can_region_start,
    NMD_region_end = median_can_region_end,
    cds_ptc_loc = ptc_pos,
    nmdesc_cds = fs_nmdesc_cds
  ) %>%
  dplyr::select(-any_of(c("genomic_loc", "PTC_pos"))) # remove the original genomic_pos column since we have renamed it to PTC_genomic_pos

gnomad_fs_all2 <- gnomad_fs_all %>%
  dplyr::select(-any_of(c('type','id','chrom','source','mutation_genomic_pos'))) %>%
  rename(nmdesc_cds = fs_nmdesc_cds)
gnomad_snv_all2 <- gnomad_snv_all %>%
  dplyr::select(-any_of(c('type','id','chrom','source','mutation_genomic_pos'))) %>%
  rename(nmdesc_cds = snv_nmdesc_cds)

#add uniprot id to gnomad_snv_all2 and gnomad_fs_all2
gnomad_snv_all2$uniprotswissprot = getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = gnomad_snv_all2$transcript,
  mart = ensembl
)$uniprotswissprot[match(gnomad_snv_all2$transcript, getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = gnomad_snv_all2$transcript,
  mart = ensembl
)$ensembl_transcript_id)]

gnomad_fs_all2$uniprotswissprot = getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = gnomad_fs_all2$transcript,
  mart = ensembl
)$uniprotswissprot[match(gnomad_fs_all2$transcript, getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = gnomad_fs_all2$transcript,
  mart = ensembl
)$ensembl_transcript_id)]

variants_all = rbind(snv_all, fs_all,gnomad_snv_all2, gnomad_fs_all2)
variants_all$source = rep(c("snv","fs","snv_control","fs_control"), c(nrow(snv_all), nrow(fs_all), nrow(gnomad_snv_all2), nrow(gnomad_fs_all2)))
#merge variants_all1 with variant_all, that is, combine their information
variants_all2 = variants_all %>%
  left_join(
    variants_all1 %>%
      dplyr::select(-any_of(c('type','id','chrom','source','mutation_genomic_pos'))),
    by = c("key" = "Variant_Key")
  )
variants_all2$cds_mutation_loc = variants_all2$cds_mutation_loc.x
variants_all2$cds_end = nchar(variants_all2$coding)
variants_all2$dist_to_cds_end = variants_all2$cds_end - variants_all2$ptc_pos

# # ------------------------------------------------------------------------------
# 4.0.3  PPI + Pfam annotation (BioMart, human_1_ interactome)
# ------------------------------------------------------------------------------

# Pfam domain coordinates per transcript, used to build pfam_fin below.
get_pfam_annotations <- function(transcript_ids, ensembl) {
  pfam_raw <- getBM(
    attributes = c("ensembl_transcript_id", "hgnc_symbol", "pfam",
                   "pfam_start", "pfam_end", "uniprotswissprot"),
    filters = "ensembl_transcript_id",
    values  = unique(transcript_ids),
    mart    = ensembl
  )

  pfam_raw %>%
    filter(!is.na(pfam), pfam != "") %>%
    distinct() %>%
    mutate(uniprot = uniprotswissprot)
}

variant_pfam_ppi <- function(variants_all2,
                             human_1_,
                             pfam_fin,
                             ensembl,
                             out_dir = ".") {
  
  .convert_to_c <- function(x) {
    if (is.na(x) || x == "") return(numeric(0))
    x <- gsub("\\[|\\]|\\s", "", x)
    if (x == "") return(numeric(0))
    vals <- suppressWarnings(as.numeric(unlist(strsplit(x, ","))))
    vals[!is.na(vals)]
  }
  
  # CDS base position -> amino acid position
  .cds2aa <- function(base_idx) ((as.integer(base_idx) - 1L) %/% 3L) + 1L
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- uniprot map from pfam_fin ---
  uniprot_map <- pfam_fin %>%
    filter(!is.na(uniprot), uniprot != "") %>%
    distinct(ensembl_transcript_id, uniprot) %>%
    group_by(ensembl_transcript_id) %>%
    slice(1) %>%
    ungroup()
  
  dat <- variants_all2 %>%
    dplyr::select(-any_of(c(
      "uniprot",
      "variant_ppi_overlap",
      "variant_ppi_nearest_interface_bp",
      "variant_ppi_dist_to_nearest_interface_bp",
      "ptc_aa",
      "max_pfam_end",
      "dist_ptc_to_max_pfam_end_aa",
      "ptc_after_max_pfam_end",
      "ptc_before_max_pfam_end"
    ))) %>%
    left_join(uniprot_map, by = "ensembl_transcript_id") %>%
    mutate(
      variant_ppi_overlap                      = 0L,
      variant_ppi_nearest_interface_bp         = NA_real_,
      variant_ppi_dist_to_nearest_interface_bp = NA_real_
    )
  
  # --- PPI: is there an interface residue downstream of the PTC? ---
  for (i in seq_len(nrow(dat))) {
    
    uid     <- dat$uniprot[[i]]
    cds_ptc <- dat$ptc_pos[[i]]
    
    if (length(uid) != 1 || length(cds_ptc) != 1) next
    if (is.na(uid) || uid == "")                   next
    if (is.na(cds_ptc))                            next
    
    re1 <- unlist(lapply(
      human_1_$interface_residues1[human_1_$uniprot1 == uid], .convert_to_c)) * 3
    re2 <- unlist(lapply(
      human_1_$interface_residues2[human_1_$uniprot2 == uid], .convert_to_c)) * 3
    
    all_iface_bp <- na.omit(c(re1, re2))
    if (length(all_iface_bp) == 0) next
    
    downstream <- all_iface_bp[all_iface_bp >= cds_ptc]
    
    if (length(downstream) > 0) {
      dat$variant_ppi_overlap[i]                      <- 1L
      dat$variant_ppi_nearest_interface_bp[i]         <- min(downstream)
      dat$variant_ppi_dist_to_nearest_interface_bp[i] <- min(downstream) - cds_ptc
    }
  }
  
  # --- Pfam: PTC position relative to the last domain end ---
  max_pfam_end_df <- pfam_fin %>%
    filter(!is.na(pfam_end)) %>%
    group_by(ensembl_transcript_id) %>%
    summarise(max_pfam_end = max(pfam_end), .groups = "drop")
  
  dat <- dat %>%
    mutate(ptc_aa = .cds2aa(ptc_pos)) %>%
    left_join(max_pfam_end_df, by = "ensembl_transcript_id") %>%
    mutate(
      dist_ptc_to_max_pfam_end_aa = ptc_aa - max_pfam_end,
      ptc_after_max_pfam_end      = as.integer(dist_ptc_to_max_pfam_end_aa > 0),
      ptc_before_max_pfam_end     = as.integer(dist_ptc_to_max_pfam_end_aa < 0)
    )
  
  # --- descriptive counts per group ---
  flag_summary <- dat %>%
    group_by(group) %>%
    summarise(
      n            = n(),
      n_ppi        = sum(variant_ppi_overlap     == 1L, na.rm = TRUE),
      prop_ppi     = n_ppi / n,
      n_pfam_before = sum(ptc_before_max_pfam_end == 1L, na.rm = TRUE),
      prop_pfam_before = n_pfam_before / n,
      .groups = "drop"
    )
  
  write.csv(dat,          file.path(out_dir, "variants_annotated.csv"), row.names = FALSE)
  write.csv(flag_summary, file.path(out_dir, "ppi_pfam_flag_summary.csv"), row.names = FALSE)
  
  # --- return the annotated variant table ---
  new_cols <- c(
    "uniprot",
    "variant_ppi_overlap",
    "variant_ppi_nearest_interface_bp",
    "variant_ppi_dist_to_nearest_interface_bp",
    "ptc_aa",
    "max_pfam_end",
    "dist_ptc_to_max_pfam_end_aa",
    "ptc_after_max_pfam_end",
    "ptc_before_max_pfam_end"
  )
  
  dat_collapsed <- dat %>%
    dplyr::select(Variant_Key, ensembl_transcript_id, all_of(new_cols)) %>%
    distinct(Variant_Key, ensembl_transcript_id, .keep_all = TRUE)
  
  variants_all2 %>%
    left_join(dat_collapsed, by = c("Variant_Key", "ensembl_transcript_id"))
}


# ------------------------------------------------------------------------------
# 4.0.4  Motif annotation (LCS / SLiM / MoRF / PTM / NLS flags)
# ------------------------------------------------------------------------------
variants_motif <- function(variants_all2,
                           mart,
                           touni_path,
                           motif_path,
                           lcs_path,
                           lcs_sheet = "G") {
  
  # helper: largest integer in a cell string
  get_max_from_cell <- function(x) {
    if (is.na(x) || x == "") return(NA_real_)
    nums <- str_extract_all(x, "\\d+")[[1]]
    if (length(nums) == 0) return(NA_real_)
    max(as.numeric(nums), na.rm = TRUE)
  }
  
  # step 1: map transcript IDs -> uniprotswissprot
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
  
  # step 2: load supplement files
  touni     <- read_csv(touni_path,  show_col_types = FALSE)
  motif_doc <- read_csv(motif_path,  show_col_types = FALSE)
  LCS_doc   <- read_excel(lcs_path,  sheet = lcs_sheet)
  
  # step 3: reduce cells to max
  motif_max <- motif_doc %>%
    mutate(across(-1, ~ vapply(as.character(.x), get_max_from_cell, numeric(1))))
  
  LCS_max <- LCS_doc %>%
    mutate(across(-1, ~ vapply(as.character(.x), get_max_from_cell, numeric(1))))
  
  # step 4: attach uniprot IDs to motif/LCS tables
  motif_max$uniprot <- touni$`Uniprot ID`[match(motif_max$Protein, touni$Protein)]
  LCS_max$uniprot   <- touni$`Uniprot ID`[match(LCS_max$Protein,   touni$Protein)]
  
  # step 5: add uniprot alias to variants
  variants_all2$uniprot <- variants_all2$uniprotswissprot
  
  # step 6: merge and compute motif flags
  variants_all2$.vm_row <- seq_len(nrow(variants_all2))
  motif_max3 <- merge(motif_max, variants_all2, by = "uniprot")
  
  motif_max3$variant_protein_feature_flag <- motif_max3$`Protein Features` * 3 >= motif_max3$cds_mutation_loc
  motif_max3$variant_domains_flag         <- motif_max3$`Domains`           * 3 >= motif_max3$cds_mutation_loc
  motif_max3$variant_slim_flag            <- motif_max3$`SLiMs`             * 3 >= motif_max3$cds_mutation_loc
  motif_max3$variant_morf_flag            <- motif_max3$MORFs               * 3 >= motif_max3$cds_mutation_loc
  motif_max3$variant_ptm_flag             <- motif_max3$`PTMs`              * 3 >= motif_max3$cds_mutation_loc
  motif_max3$variant_nls_flag             <- motif_max3$`NLSs`              * 3 >= motif_max3$cds_mutation_loc
  
  # step 7: merge and compute LCS flag
  LCS_max3 <- merge(LCS_max, variants_all2, by = "uniprot")
  LCS_max3$variant_LCS_flag <- LCS_max3$`LCSs` * 3 >= LCS_max3$cds_mutation_loc
  
  # step 8: map flags back by VARIANT-level key (was: by uniprot)
  im <- match(variants_all2$.vm_row, motif_max3$.vm_row)
  variants_all2$variant_protein_flag <- motif_max3$variant_protein_feature_flag[im]
  variants_all2$variant_domain_flag  <- motif_max3$variant_domains_flag[im]
  variants_all2$variant_slim_flag    <- motif_max3$variant_slim_flag[im]
  variants_all2$variant_morf_flag    <- motif_max3$variant_morf_flag[im]
  variants_all2$variant_ptm_flag     <- motif_max3$variant_ptm_flag[im]
  variants_all2$variant_nls_flag     <- motif_max3$variant_nls_flag[im]
  il <- match(variants_all2$.vm_row, LCS_max3$.vm_row)
  variants_all2$variant_LCS_flag     <- LCS_max3$variant_LCS_flag[il]
  
  variants_all2$.vm_row <- NULL
  variants_all2
}

annotate_motifs <- function(variants_annotated, ensembl) {
  variants_motif(
    variants_all2 = variants_annotated,
    mart          = ensembl,
    touni_path = file.path(OTHERS_DIR, "NIHMS1818854-supplement-2(A).csv"),
    motif_path = file.path(OTHERS_DIR, "NIHMS1818854-supplement-2(B).csv"),
    lcs_path   = file.path(OTHERS_DIR, "Copy of NIHMS1818854-supplement-2.xls")
  )
}

# ------------------------------------------------------------------------------
# 4.0.5  Final variant table with derived features
# ------------------------------------------------------------------------------

prepare_final_variant_table <- function(variants_all4) {
  variants_all4 %>%
    filter(!duplicated(Variant_Key)) %>%
    mutate(
      is_disease          = as.integer(grepl("disease", group)),
      gene_set            = case_when(
        group %in% c("snv_disease", "snv_control") ~ "SNV",
        group %in% c("fs_disease",  "fs_control")  ~ "FS"
      ),
      dist_to_cds_end_log = log(abs(dist_to_cds_end) + 1)
    )
}

# ------------------------------------------------------------------------------
# 4.0.6  Generic modeling helpers
# ------------------------------------------------------------------------------

# Extract exponentiated OR + 95% CI + p-value for one coefficient of a glm fit
extract_or <- function(fit, coef_name) {
  ci <- tryCatch(confint(fit)[coef_name, ], error = function(e) c(NA, NA))
  data.frame(
    flag      = coef_name,
    OR        = exp(coef(fit)[coef_name]),
    OR_low    = exp(ci[1]),
    OR_high   = exp(ci[2]),
    p_value   = summary(fit)$coefficients[coef_name, "Pr(>|z|)"],
    row.names = NULL
  )
}

# Flags with < 2 unique non-NA values within a gene_set can't be modeled -> skip
flags_without_variation <- function(data, gene_set_label, flag_cols = FLAG_COLS) {
  data %>%
    filter(gene_set == gene_set_label) %>%
    dplyr::select(all_of(flag_cols)) %>%
    summarise(across(everything(), ~ n_distinct(na.omit(.)))) %>%
    pivot_longer(everything(), names_to = "flag", values_to = "n_unique") %>%
    filter(n_unique < 2) %>%
    pull(flag)
}

# Single binary/continuous logistic regression: is_disease ~ flag [+ covariates]
fit_flag_glm <- function(data, flag, covariates = character(0)) {
  df <- data %>%
    dplyr::select(is_disease, all_of(flag), all_of(covariates)) %>%
    filter(if_all(everything(), ~ !is.na(.))) %>%
    mutate(across(all_of(flag), as.integer))
  
  if (nrow(df) == 0 ||
      n_distinct(df[[flag]]) < 2 ||
      n_distinct(df$is_disease) < 2) {
    message(sprintf("Skipping %s: insufficient variation", flag))
    return(NULL)
  }
  
  rhs <- paste(c(flag, covariates), collapse = " + ")
  fit <- tryCatch(
    glm(as.formula(paste("is_disease ~", rhs)), data = df, family = binomial),
    error = function(e) { message(sprintf("Model failed for %s: %s", flag, e$message)); NULL }
  )
  if (is.null(fit)) return(NULL)
  extract_or(fit, flag)
}



# BH-adjust p-values within strata and assign sig stars.
# Pass group_vars = character(0) to treat the whole table as one family.
add_fdr_sig <- function(results, group_vars = c("model")) {
  out <- results
  if (length(group_vars) > 0) out <- out %>% group_by(across(all_of(group_vars)))
  out %>%
    mutate(
      p_adj = p.adjust(p_value, method = FDR_METHOD),
      n_in_family = n(),
      sig   = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    ) %>%
    ungroup()
}

# Shared forest-plot template for OR / ratio results
plot_or_forest <- function(results, facet_var = "model", shape_var = NULL,
                           title = "", subtitle = "") {
  p <- ggplot(results, aes(x = OR, y = reorder(flag, OR), color = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.2,
                   position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c("***" = "red", "**" = "orange",
                                  "*" = "gold", "ns" = "grey60")) +
    scale_x_log10() +
    labs(title = title, subtitle = subtitle,
         x = "OR / ratio (log scale)", y = NULL, color = "Adj. p") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  if (!is.null(shape_var)) {
    p <- p + aes(shape = .data[[shape_var]]) +
      geom_point(size = 3, position = position_dodge(width = 0.5)) +
      scale_shape_manual(values = c("SNV" = 16, "FS" = 17))
  } else {
    p <- p + geom_point(size = 3)
  }
  
  if (!is.null(facet_var)) p <- p + facet_wrap(as.formula(paste("~", facet_var)))
  p
}


# ==============================================================================
# 4.1  UNMATCHED ANALYSIS
# ------------------------------------------------------------------------------
# Comparisons of disease vs. control variants: PTC-distance sanity
# checks, PPI-overlap by interactome source, and the full flag panel
# (unadjusted vs. covariate-adjusted).
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1.1  Simple unadjusted PTC-distance models (sanity check)
# ------------------------------------------------------------------------------

run_dist_sanity_models <- function(variant_data2) {
  variant_data2 <- variant_data2 %>%
    mutate(is_nmdesc = if_else(group %in% c("fs_control", "snv_control"), 0L, 1L))
  
  snv_df <- filter(variant_data2, group %in% c("snv_disease", "snv_control"))
  fs_df  <- filter(variant_data2, group %in% c("fs_disease",  "fs_control"))
  
  list(
    snv_dist        = glm(is_nmdesc ~ dist_to_cds_end,           data = snv_df, family = binomial),
    fs_dist         = glm(is_nmdesc ~ dist_to_cds_end,           data = fs_df,  family = binomial),
    snv_dist_cdsend = glm(is_nmdesc ~ dist_to_cds_end + cds_end, data = snv_df, family = binomial),
    fs_dist_cdsend  = glm(is_nmdesc ~ dist_to_cds_end + cds_end, data = fs_df,  family = binomial)
  )
}

# Tidies the sanity models into one table and BH-adjusts
# across all non-intercept coefficients from all four models.
tidy_dist_models <- function(models, labels = FLAG_LABELS) {
  purrr::imap_dfr(models, function(fit, nm) {
    broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term != "(Intercept)") %>%
      mutate(model = nm)
  }) %>%
    mutate(
      p_adj = p.adjust(p.value, method = FDR_METHOD),
      sig   = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      )
    )
}

# 4.2  MIXED-EFFECT MODEL  
SEPARATION_MIN_CELL <- 5

# 4.2.1  Fit one flag within one gene set
fit_one_glmm <- function(dat, flag, gs) {
  
  is_cont <- flag %in% CONT_PREDICTORS
  df <- dat %>%
    dplyr::select(is_disease, all_of(flag), ensembl_transcript_id) %>%
    filter(!is.na(.data[[flag]]))
  if (!is_cont) df <- df %>% mutate(across(all_of(flag), as.integer))
  
  base <- data.frame(
    flag = flag, gene_set = gs,
    n_variants = nrow(df),
    n_transcripts = n_distinct(df$ensembl_transcript_id),
    OR = NA_real_, OR_low = NA_real_, OR_high = NA_real_, p_value = NA_real_,
    singular = NA, method = NA_character_, row.names = NULL
  )
  
  fit <- tryCatch(
    lme4::glmer(as.formula(paste("is_disease ~", flag, "+ (1 | ensembl_transcript_id)")),
                data = df, family = binomial,
                control = lme4::glmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 2e5))),
    error = function(e) { message(sprintf("[%s] %s failed: %s", gs, flag, e$message)); NULL }
  )
  if (is.null(fit)) { base$method <- "failed"; return(base) }
  
  base$singular <- lme4::isSingular(fit)
  
  tid <- tryCatch(
    broom.mixed::tidy(fit, effects = "fixed", exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == flag),
    error = function(e) NULL
  )
  if (is.null(tid) || nrow(tid) == 0) { base$method <- "failed"; return(base) }
  
  base$OR      <- tid$estimate[1]
  base$OR_low  <- tid$conf.low[1]
  base$OR_high <- tid$conf.high[1]
  base$p_value <- tid$p.value[1]
  base$method  <- "GLMM"
  base
}

# 4.2.2  Driver
run_mixed_effect_flag_analysis <- function(variants_all5,
                                           flag_cols = GENE_MATCHED_FLAGS,
                                           gene_sets = c("SNV", "FS")) {
  
  map_dfr(gene_sets, function(gs) {
    dat <- filter(variants_all5, gene_set == gs)
    map_dfr(intersect(flag_cols, names(dat)), ~ fit_one_glmm(dat, .x, gs))
  })
}

# 4.2.3  Tidy and adjust
tidy_mixed_results <- function(mixed_results, labels = FLAG_LABELS) {
  mixed_results %>%
    mutate(model = "Mixed-effect") %>%
    group_by(gene_set, model) %>%
    mutate(
      p_adj       = ifelse(is.na(p_value), NA_real_,
                           p.adjust(replace(p_value, is.na(p_value), 1),
                                    method = FDR_METHOD)),
      n_in_family = sum(!is.na(p_value))
    ) %>%
    ungroup() %>%
    mutate(
      sig = case_when(
        is.na(p_adj)  ~ "not estimable",
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ "ns"
      ),
      flag_lab = ifelse(flag %in% names(labels), labels[flag], flag)
    )
}

# 4.2.4  Plot
plot_mixed_effect_flags <- function(mixed_results, or_limits = c(0.01, 1000)) {
  d <- tidy_mixed_results(mixed_results) %>%
    filter(!is.na(OR)) %>%
    mutate(
      trustworthy = (method == "GLMM") & !singular,
      OR_low      = pmax(OR_low,  or_limits[1]),
      OR_high     = pmin(OR_high, or_limits[2])
    )
  
  ggplot(d, aes(x = OR, y = reorder(flag_lab, OR), color = sig, shape = gene_set)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.2,
                   position = position_dodge(width = 0.6)) +
    geom_point(aes(alpha = trustworthy), size = 3,
               position = position_dodge(width = 0.6)) +
    scale_color_manual(values = c("***" = "red", "**" = "orange", "*" = "gold",
                                  "ns" = "grey60", "not estimable" = "grey85")) +
    scale_shape_manual(values = c("SNV" = 16, "FS" = 17)) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35),
                       labels = c(`TRUE` = "GLMM", `FALSE` = "singular"),
                       name = "Fit") +
    scale_x_log10(limits = or_limits) +
    labs(
      title    = "Mixed-effect model: disease association per feature",
      subtitle = "Logistic GLMM by gene set, random intercept per transcript (BH-adjusted within gene set)",
      x = "Odds ratio (log scale)", y = NULL, color = "Adj. p", shape = "Gene set"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

# 4.3  HIERARCHICAL BAYESIAN MODEL
# Partial pooling of flag effects across genes within each gene set (SNV/FS).
run_bayesian_by_geneset <- function(variants_all5, flag_cols = GENE_MATCHED_FLAGS) {
  library(brms)
  gene_sets <- c("SNV", "FS")
  
  purrr::map(gene_sets, function(gs) {
    dat_gs <- filter(variants_all5, gene_set == gs)   # keeps *_disease and *_control for this set
    
    purrr::map(flag_cols, function(flag) {
      df <- dat_gs %>%
        dplyr::select(is_disease, all_of(flag), ensembl_transcript_id) %>%
        filter(!is.na(.data[[flag]]))
      if (!(flag %in% CONT_PREDICTORS))
        df <- df %>% mutate(across(all_of(flag), as.integer))
      
      tryCatch(
        brm(as.formula(paste("is_disease ~", flag, "+ (1 | ensembl_transcript_id)")),
            data = df, family = bernoulli(), chains = 4, cores = 4, refresh = 0),
        error = function(e) { message(sprintf("[%s] %s failed: %s", gs, flag, e$message)); NULL }
      )
    }) %>% setNames(flag_cols)
  }) %>% setNames(gene_sets)
}

tidy_bayesian_by_geneset <- function(fit_lists, labels = FLAG_LABELS) {
  purrr::imap_dfr(fit_lists, function(flag_fits, gs) {
    purrr::imap_dfr(flag_fits, function(fit, flag) {
      if (!inherits(fit, "brmsfit")) return(NULL)
      fx  <- brms::fixef(fit, probs = c(0.025, 0.975))
      row <- if (flag %in% rownames(fx)) flag else setdiff(rownames(fx), "Intercept")[1]
      est <- fx[row, ]
      data.frame(
        flag     = flag,
        gene_set = gs,
        OR       = exp(est[["Estimate"]]),
        OR_low   = exp(est[["Q2.5"]]),
        OR_high  = exp(est[["Q97.5"]]),
        row.names = NULL
      )
    })
  }) %>%
    mutate(
      sig  = ifelse(OR_low > 1 | OR_high < 1, "*", "ns"),
      flag = ifelse(flag %in% names(labels), labels[flag], flag)
    )
}

plot_bayesian_by_geneset <- function(fit_lists) {
  plot_or_forest(
    tidy_bayesian_by_geneset(fit_lists),
    facet_var = NULL, shape_var = "gene_set",   # SNV/FS shown as different shapes
    title    = "Bayesian disease-vs-control association by gene set",
    subtitle = "Posterior OR with 95% credible interval; multiplicity handled by partial pooling"
  ) + ggplot2::labs(color = "95% CrI", shape = "Gene set")
}

##run analysis
# --- 4.0: data + annotation ---------------------------------------------------
pfam_fin     <- get_pfam_annotations(variants_all2$ensembl_transcript_id, ensembl)
human_1_     <- read_delim(data_file("human (1).txt"),
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE)
variants_all2 <- variants_all2 %>% filter(!is.na(group))
variants_all2$Variant_Key = variants_all2$key
variants_all3 <- variant_pfam_ppi(
  variants_all2 = variants_all2,
  human_1_      = human_1_,
  pfam_fin      = pfam_fin,
  ensembl       = ensembl
)
variants_all4 <- annotate_motifs(variants_all3, ensembl)
variants_all5 <- prepare_final_variant_table(variants_all4)
#remove cds_mutation_loc.x and cds_mutation_loc.y columns
variants_all5 <- variants_all5 %>% dplyr::select(-cds_mutation_loc.x, -cds_mutation_loc.y)
write.csv(variants_all5, "variants_all0805.csv", row.names = FALSE)

# --- 4.1: unmatched analysis ---------------------------------------------------
dist_models  <- run_dist_sanity_models(variants_all5)
dist_results <- tidy_dist_models(dist_models)          # BH-adjusted table
print(dist_results)
write.csv(dist_results, file.path(OUT_DIR, "dist_sanity_models_fdr.csv"), row.names = FALSE)

# --- 4.2: mixed-effect model ---------------------------------------------------
mixed_results <- run_mixed_effect_flag_analysis(variants_all5)
mixed_tidy    <- tidy_mixed_results(mixed_results)     # adjusted table
write.csv(mixed_tidy, file.path(OUT_DIR, "mixed_effect_fdr.csv"), row.names = FALSE)
plot_mixed_effect_flags(mixed_results)

# --- 4.3: hierarchical Bayesian model -------------------------------------------
bayes_by_gs <- run_bayesian_by_geneset(variants_all5, GENE_MATCHED_FLAGS)
plot_bayesian_by_geneset(bayes_by_gs)

