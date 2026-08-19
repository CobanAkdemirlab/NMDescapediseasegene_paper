# ==============================================================================
# 4. VARIANT-LEVEL COMPARISON
# ------------------------------------------------------------------------------
#   4.0  Config, data loading, annotation pipeline, and shared helper functions
#   4.1  Unmatched analysis
#   4.2  Mixed-effect model
#   4.3  Hierarchical Bayesian model
#   4.4  Bootstrap and gene-matched analysis
#
# MULTIPLE TESTING CORRECTION SUMMARY
#   4.0.3  .pval_df()            -> BH across the Fisher bar-chart comparisons
#   4.0.3  stat_compare_means()  -> BH across the Wilcoxon violin comparisons
#   4.1    tidy_dist_models()    -> BH across the PTC-distance sanity models
#   4.2    add_fdr_sig()         -> BH within model (already present)
#   4.3    none by design: partial pooling handles multiplicity
#   4.4    add_fdr_sig()         -> BH within (gene_set, model) (already present)
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

DATA_DIR    <- "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data"
CLINVAR_DIR <- file.path(DATA_DIR, "clinvar")
OTHERS_DIR  <- file.path(DATA_DIR, "others")
OUT_DIR     <- "plots/pfam_ppi_analysis"

FDR_METHOD <- "BH"   # single place to change the adjustment method

FLAG_COLS <- c(
  "variant_ppi_overlap", "ptc_after_max_pfam_end", "ptc_before_max_pfam_end",
  "variant_protein_flag", "variant_domain_flag", "variant_slim_flag",
  "variant_morf_flag", "variant_ptm_flag", "variant_nls_flag", "variant_LCS_flag"
)
CONT_COL         <- "dist_to_cds_end_log"
CONFOUNDERS      <- c("cds_length")                                    # base adjustment set
CONFOUNDERS_FULL <- c("cds_length", "NMDesc_region_length", "GC_Content")

# Flags for the gene-matched panel: one Pfam flag, one PPI flag, plus motifs
GENE_MATCHED_FLAGS <- c(
  "variant_ppi_overlap",       # PPI: the single overlap flag
  "ptc_before_max_pfam_end",   # Pfam: keep the "domain disrupted" direction
  "variant_protein_flag", "variant_domain_flag", "variant_slim_flag",
  "variant_morf_flag", "variant_ptm_flag", "variant_nls_flag", "variant_LCS_flag"
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
source('/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/variant level/new_create_fasta_functions.R')
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
snv_variants = read.csv('snv_variants20260201_plp_dbh_clinvar.csv')
snv_dis <- create_fasta(snv_variants, output_dir = "snv_disease_fasta_output_test2")

snv_control_variants = read.csv('gnomad_snv_filtered_wald.csv')
gnomad_snv_filtered <- snv_control_variants
gnomad_snv_filtered$key <- gnomad_snv_filtered$id
gnomad_snv_variants <- gnomad_snv_filtered[,  c("transcript", "key")]
gnomad_snv_dis <- create_fasta(gnomad_snv_variants, output_dir = "snv_control_fasta_output_test2")

fs_variants = read.csv('fs_variants20260201_plp_acat_clinvar.csv')
fs_dis <- create_fasta(fs_variants, output_dir = "fs_disease_fasta_output_test2")

fs_control_variants = read.csv('gnomad_fs_filtered_wald.csv')
gnomad_fs_filtered  <- fs_control_variants
gnomad_fs_filtered$key  <- gnomad_fs_filtered$id
#remove inframe indel, get gnomad_fs_filtered2
gnomad_fs_filtered$ref = sub(".*\\|(.+)\\|.*", "\\1", gnomad_fs_filtered$id)
gnomad_fs_filtered$alt = sub(".*\\|.*\\|(.*)", "\\1", gnomad_fs_filtered$id)
gnomad_fs_filtered$len_diff = abs(nchar(gnomad_fs_filtered$ref) - nchar(gnomad_fs_filtered$alt))
gnomad_fs_filtered2 = gnomad_fs_filtered[gnomad_fs_filtered$len_diff %% 3 != 0, ]
gnomad_fs_filtered2$ref = NULL
gnomad_fs_filtered2$alt = NULL
gnomad_fs_filtered2$len_diff = NULL

gnomad_fs_variants  <- gnomad_fs_filtered[,  c("transcript", "key")]
gnomad_fs_dis <- create_fasta(gnomad_fs_variants, output_dir = "fs_control_fasta_output_test2")

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

# ------------------------------------------------------------------------------
# 4.0.3  PPI + Pfam annotation (BioMart, human_1_ interactome)
# ------------------------------------------------------------------------------

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


variant_pfam_ppi <- function(
    variants_all2,
    human_1_,
    pfam_fin,
    ensembl,
    out_dir      = ".",
    group_levels = c("fs_disease", "fs_control", "snv_disease", "snv_control"),
    group_colors = c(
      "fs_disease"  = "#1f77b4",
      "fs_control"  = "#aec7e8",
      "snv_disease" = "#2ca02c",
      "snv_control" = "#98df8a"
    )
) {
  
  .convert_to_c <- function(x) {
    if (is.na(x) || x == "") return(numeric(0))
    x <- gsub("\\[|\\]|\\s", "", x)
    if (x == "") return(numeric(0))
    vals <- suppressWarnings(as.numeric(unlist(strsplit(x, ","))))
    vals[!is.na(vals)]
  }
  
  .convert_to_ints <- function(x) {
    if (is.null(x) || is.na(x) || x == "") return(integer(0))
    s <- trimws(gsub("\\[|\\]", "", x))
    if (s == "") return(integer(0))
    as.integer(trimws(unlist(strsplit(s, ","))))
  }
  
  # CDS base position -> amino acid position
  .cds2aa <- function(base_idx) ((as.integer(base_idx) - 1L) %/% 3L) + 1L
  
  # Fisher's exact test -> p-value
  .fisher_p <- function(n_event_a, n_total_a, n_event_b, n_total_b) {
    mat <- matrix(
      c(n_event_a,   n_total_a - n_event_a,
        n_event_b,   n_total_b - n_event_b),
      nrow = 2, byrow = TRUE
    )
    fisher.test(mat)$p.value
  }
  
  # Build p-value bracket data frame for stat_pvalue_manual()
  # CHANGED: raw Fisher p-values are BH-adjusted across the comparisons in this
  # panel; the bracket label now reports the adjusted value.
  .pval_df <- function(summary_df, prop_col, event_col, total_col,
                       pairs, y_mult = c(1.08, 1.20)) {
    ymax <- max(summary_df[[prop_col]], na.rm = TRUE)
    rows <- lapply(seq_along(pairs), function(i) {
      g1 <- pairs[[i]][1]; g2 <- pairs[[i]][2]
      if (!all(c(g1, g2) %in% summary_df$group)) return(NULL)
      p <- .fisher_p(
        summary_df[[event_col]][summary_df$group == g1],
        summary_df[[total_col]][summary_df$group == g1],
        summary_df[[event_col]][summary_df$group == g2],
        summary_df[[total_col]][summary_df$group == g2]
      )
      mult <- if (i <= length(y_mult)) y_mult[i] else y_mult[length(y_mult)] + 0.12 * (i - length(y_mult))
      data.frame(group1 = g1, group2 = g2,
                 y.position = ymax * mult,
                 p_raw      = p,
                 stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, Filter(Negate(is.null), rows))
    if (is.null(out)) return(NULL)
    out$p_adj <- p.adjust(out$p_raw, method = FDR_METHOD)
    out$label <- paste0("p.adj = ", signif(out$p_adj, 3))
    out
  }
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  #get pfam info
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
  
  #get ppi info
  for (i in seq_len(nrow(dat))) {
    
    uid     <- dat$uniprot[[i]]
    cds_ptc <- dat$ptc_pos[[i]]
    
    if (length(uid) != 1 || length(cds_ptc) != 1) next
    if (is.na(uid)  || uid == "")                  next
    if (is.na(cds_ptc))                            next
    
    re1 <- unlist(lapply(
      human_1_$interface_residues1[human_1_$uniprot1 == uid],
      .convert_to_c
    )) * 3
    
    re2 <- unlist(lapply(
      human_1_$interface_residues2[human_1_$uniprot2 == uid],
      .convert_to_c
    )) * 3
    
    all_iface_bp <- na.omit(c(re1, re2))
    if (length(all_iface_bp) == 0) next
    
    downstream <- all_iface_bp[all_iface_bp >= cds_ptc] #if ppi residue is downstream of PTC, calculate distance
    
    if (length(downstream) > 0) {
      dat$variant_ppi_overlap[i]                      <- 1L
      dat$variant_ppi_nearest_interface_bp[i]         <- min(downstream)
      dat$variant_ppi_dist_to_nearest_interface_bp[i] <- min(downstream) - cds_ptc
    }
  }
  
  #clean data
  present_levels <- unique(as.character(dat$group))
  ordered_levels <- c(
    intersect(group_levels, present_levels),
    setdiff(present_levels, group_levels)
  )
  dat$group <- factor(dat$group, levels = ordered_levels)
  
  # Restrict group_colors to levels actually present
  group_colors <- group_colors[names(group_colors) %in% present_levels]
  
  # Pair adjacent levels for comparisons: (1 vs 2), (3 vs 4), etc.
  comparisons <- lapply(
    seq(1, length(ordered_levels) - 1, by = 2),
    function(k) ordered_levels[k:(k + 1)]
  )
  
  #plot ppi
  ppi_summary <- dat %>%
    group_by(group) %>%
    summarise(
      n       = n(),
      n_match = sum(variant_ppi_overlap == 1L, na.rm = TRUE),
      prop    = n_match / n,
      .groups = "drop"
    )
  
  pv_ppi   <- .pval_df(ppi_summary, "prop", "n_match", "n", comparisons)
  ymax_ppi <- max(ppi_summary$prop, na.rm = TRUE)
  
  p_ppi_prop <- ggplot(ppi_summary, aes(x = group, y = prop, fill = group)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = group_colors, guide = "none") +
    geom_text(aes(label = paste0(n_match, "/", n)), vjust = -0.3, size = 4) +
    stat_pvalue_manual(pv_ppi, label = "label", xmin = "group1", xmax = "group2",
                       y.position = "y.position", tip.length = 0.01) +
    expand_limits(y = ymax_ppi * 1.35) +
    labs(title = "Proportion of variants with downstream PPI interface residue",
         x = NULL, y = "Proportion") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  # CHANGED: Wilcoxon comparisons now BH-adjusted (p.adjust.method + label)
  p_ppi_dist <- ggplot(
    dat %>% filter(!is.na(variant_ppi_dist_to_nearest_interface_bp)),
    aes(x = group, y = variant_ppi_dist_to_nearest_interface_bp, fill = group)
  ) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
    scale_fill_manual(values = group_colors, guide = "none") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       p.adjust.method = FDR_METHOD, label = "p.adj") +
    labs(title = "Distance from PTC to nearest downstream PPI interface residue",
         x = NULL, y = "Distance (bp)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  # CHANGED: Wilcoxon comparisons now BH-adjusted
  p_ppi_log <- ggplot(
    dat %>% filter(!is.na(variant_ppi_dist_to_nearest_interface_bp),
                   variant_ppi_dist_to_nearest_interface_bp > 0),
    aes(x = group, y = variant_ppi_dist_to_nearest_interface_bp, fill = group)
  ) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
    scale_y_log10() +
    scale_fill_manual(values = group_colors, guide = "none") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       p.adjust.method = FDR_METHOD, label = "p.adj") +
    labs(title = "Log10 distance from PTC to nearest downstream PPI interface residue",
         x = NULL, y = "Distance (bp, log10)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  #plot pfam
  pfam_bm <- getBM(
    attributes = c("ensembl_transcript_id", "pfam_end"),
    filters    = "ensembl_transcript_id",
    values     = unique(dat$ensembl_transcript_id),
    mart       = ensembl
  )
  
  max_pfam_end_df <- pfam_bm %>%
    filter(!is.na(pfam_end)) %>%
    group_by(ensembl_transcript_id) %>%
    summarise(max_pfam_end = max(pfam_end), .groups = "drop")
  
  dat <- dat %>%
    mutate(ptc_aa = .cds2aa(ptc_pos)) %>%
    dplyr::select(-any_of(c("max_pfam_end",
                            "dist_ptc_to_max_pfam_end_aa",
                            "ptc_after_max_pfam_end",
                            "ptc_before_max_pfam_end"))) %>%
    left_join(max_pfam_end_df, by = "ensembl_transcript_id") %>%
    mutate(
      dist_ptc_to_max_pfam_end_aa = ptc_aa - max_pfam_end,
      ptc_after_max_pfam_end      = as.integer(dist_ptc_to_max_pfam_end_aa >  0),
      ptc_before_max_pfam_end     = as.integer(dist_ptc_to_max_pfam_end_aa <  0)
    )
  
  dat <- dat %>%
    left_join(
      pfam_fin %>%
        dplyr::select(ensembl_transcript_id, pfam, pfam_start, pfam_end, hgnc_symbol),
      by = "ensembl_transcript_id"
    ) %>%
    mutate(
      in_pfam = !is.na(pfam_start) & !is.na(pfam_end) &
        !is.na(ptc_aa)     & ptc_aa <= pfam_end
    )
  
  pfam_summary <- dat %>%
    group_by(group) %>%
    summarise(
      n           = n(),
      n_after     = sum(ptc_after_max_pfam_end  == 1L, na.rm = TRUE),
      n_before    = sum(ptc_before_max_pfam_end == 1L, na.rm = TRUE),
      prop_after  = n_after  / n,
      prop_before = n_before / n,
      .groups = "drop"
    )
  
  pv_pfam   <- .pval_df(pfam_summary, "prop_after", "n_after", "n", comparisons)
  ymax_pfam <- max(pfam_summary$prop_after, na.rm = TRUE)
  
  # CHANGED: Wilcoxon comparisons now BH-adjusted
  p_pfam_dist <- ggplot(
    dat %>% filter(!is.na(dist_ptc_to_max_pfam_end_aa)),
    aes(x = group, y = dist_ptc_to_max_pfam_end_aa, fill = group)
  ) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
    scale_fill_manual(values = group_colors, guide = "none") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       p.adjust.method = FDR_METHOD, label = "p.adj") +
    labs(title = "Distance from PTC to the largest PFAM end",
         x = NULL, y = "PTC aa position - largest PFAM end (aa)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  # CHANGED: Wilcoxon comparisons now BH-adjusted
  p_pfam_abs <- ggplot(
    dat %>%
      filter(!is.na(dist_ptc_to_max_pfam_end_aa)) %>%
      mutate(abs_dist = abs(dist_ptc_to_max_pfam_end_aa)),
    aes(x = group, y = abs_dist, fill = group)
  ) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
    scale_fill_manual(values = group_colors, guide = "none") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       p.adjust.method = FDR_METHOD, label = "p.adj") +
    labs(title = "Absolute distance from PTC to the largest PFAM end",
         x = NULL, y = "Absolute distance (aa)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  p_pfam_prop <- ggplot(pfam_summary,
                        aes(x = group, y = prop_before, fill = group)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = group_colors, guide = "none") +
    geom_text(aes(label = paste0(n_before, "/", n)), vjust = -0.3, size = 4) +
    stat_pvalue_manual(pv_pfam, label = "label", xmin = "group1", xmax = "group2",
                       y.position = "y.position", tip.length = 0.01) +
    expand_limits(y = ymax_pfam * 1.35) +
    labs(title = "Proportion of variants with PTC influencing PFAM domain",
         x = NULL, y = "Proportion") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  pfam_pct_df <- dat %>%
    count(group, in_pfam) %>%
    group_by(group) %>%
    tidyr::complete(in_pfam = c(FALSE, TRUE), fill = list(n = 0)) %>%
    mutate(prop = n / sum(n),
           pct  = scales::percent(prop, accuracy = 0.1)) %>%
    ungroup()
  
  p_pfam_stack <- ggplot(pfam_pct_df,
                         aes(x = group, y = prop, fill = as.character(in_pfam))) +
    geom_col(color = "grey30") +
    geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "#2ca02c"),
                      labels = c("Outside PFAM", "Influences PFAM"), name = NULL) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(title = "% of variants influencing PFAM domains",
         x = NULL, y = "Fraction of variants") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  #output
  write.csv(dat,          file.path(out_dir, "variants_annotated.csv"),              row.names = FALSE)
  write.csv(ppi_summary,  file.path(out_dir, "ppi_flag_summary.csv"),                row.names = FALSE)
  write.csv(pfam_summary, file.path(out_dir, "pfam_flag_summary.csv"),               row.names = FALSE)
  write.csv(pfam_pct_df,  file.path(out_dir, "pfam_influenced_variant_counts.csv"),  row.names = FALSE)
  # CHANGED: export raw and adjusted p-values for the bar-chart comparisons
  if (!is.null(pv_ppi))  write.csv(pv_ppi,  file.path(out_dir, "ppi_prop_pvalues.csv"),  row.names = FALSE)
  if (!is.null(pv_pfam)) write.csv(pv_pfam, file.path(out_dir, "pfam_prop_pvalues.csv"), row.names = FALSE)
  
  plot_list <- list(
    ppi_proportion    = p_ppi_prop,
    ppi_distance      = p_ppi_dist,
    ppi_distance_log  = p_ppi_log,
    pfam_distance     = p_pfam_dist,
    pfam_abs_distance = p_pfam_abs,
    pfam_proportion   = p_pfam_prop,
    pfam_stacked      = p_pfam_stack
  )
  
  for (nm in names(plot_list)) {
    ggsave(file.path(out_dir, paste0(nm, ".pdf")), plot_list[[nm]], width = 9, height = 5)
    ggsave(file.path(out_dir, paste0(nm, ".png")), plot_list[[nm]], width = 9, height = 5, dpi = 300)
  }
  
  #return new variant_all dataframe
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
  
  # dat may have been expanded by the pfam_fin join (multiple domains per
  # transcript); collapse back to one row per Variant_Key before joining
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
  
  # step 8: map flags back to variants_all2
  variants_all2$variant_protein_flag <- motif_max3$variant_protein_feature_flag[match(variants_all2$uniprot, motif_max3$uniprot)]
  variants_all2$variant_domain_flag  <- motif_max3$variant_domains_flag[match(variants_all2$uniprot, motif_max3$uniprot)]
  variants_all2$variant_slim_flag    <- motif_max3$variant_slim_flag[match(variants_all2$uniprot, motif_max3$uniprot)]
  variants_all2$variant_morf_flag    <- motif_max3$variant_morf_flag[match(variants_all2$uniprot, motif_max3$uniprot)]
  variants_all2$variant_ptm_flag     <- motif_max3$variant_ptm_flag[match(variants_all2$uniprot, motif_max3$uniprot)]
  variants_all2$variant_nls_flag     <- motif_max3$variant_nls_flag[match(variants_all2$uniprot, motif_max3$uniprot)]
  variants_all2$variant_LCS_flag     <- LCS_max3$variant_LCS_flag[match(variants_all2$uniprot, LCS_max3$uniprot)]
  
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

# Standardize group -> is_disease / gene_set encoding, drop duplicated .x/.y cols
encode_groups <- function(data) {
  data %>%
    dplyr::select(-ends_with(".x")) %>%
    rename_with(~ gsub("\\.y$", "", .), ends_with(".y")) %>%
    mutate(
      is_disease = as.integer(group %in% c("snv_disease", "fs_disease")),
      gene_set   = case_when(
        group %in% c("snv_disease", "snv_control") ~ "SNV",
        group %in% c("fs_disease",  "fs_control")  ~ "FS"
      )
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

# Draw one random variant per gene per group, return a wide gene-level pair
sample_one_per_gene <- function(data, disease_grp, control_grp, value_col) {
  disease_samp <- data %>%
    filter(group == disease_grp, !is.na(.data[[value_col]])) %>%
    group_by(ensembl_transcript_id) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    rename(disease = all_of(value_col)) %>%
    dplyr::select(ensembl_transcript_id, disease)
  
  control_samp <- data %>%
    filter(group == control_grp, !is.na(.data[[value_col]])) %>%
    group_by(ensembl_transcript_id) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    rename(control = all_of(value_col)) %>%
    dplyr::select(ensembl_transcript_id, control)
  
  inner_join(disease_samp, control_samp, by = "ensembl_transcript_id")
}

# Gene-matched paired-Wilcoxon test + bootstrap CI for one variable
# (summary_fn = mean for binary flags -> gene-level proportion;
#  summary_fn = median for continuous distance -> gene-level median)
run_gene_matched_one <- function(data, gene_set_label, value_col,
                                 summary_fn = mean, n_boot = 1000, seed = 42) {
  
  disease_grp <- if (gene_set_label == "SNV") "snv_disease" else "fs_disease"
  control_grp <- if (gene_set_label == "SNV") "snv_control" else "fs_control"
  
  gene_summary <- data %>%
    filter(!is.na(.data[[value_col]])) %>%
    group_by(ensembl_transcript_id, group) %>%
    summarise(val = summary_fn(.data[[value_col]]), .groups = "drop")
  
  paired <- inner_join(
    gene_summary %>% filter(group == disease_grp) %>%
      rename(disease = val) %>% dplyr::select(ensembl_transcript_id, disease),
    gene_summary %>% filter(group == control_grp) %>%
      rename(control = val) %>% dplyr::select(ensembl_transcript_id, control),
    by = "ensembl_transcript_id"
  )
  
  message(sprintf("%s / %s: %d paired genes", gene_set_label, value_col, nrow(paired)))
  if (nrow(paired) < 2) return(NULL)
  
  if (all(paired$disease == paired$control)) {
    return(data.frame(flag = value_col, gene_set = gene_set_label,
                      OR = 1, OR_low = NA, OR_high = NA, p_value = 1))
  }
  
  w_test <- wilcox.test(paired$disease, paired$control, paired = TRUE, exact = FALSE)
  ratio  <- summary_fn(paired$disease) / summary_fn(paired$control)
  
  set.seed(seed)
  boot_ratios <- replicate(n_boot, {
    paired_b <- sample_one_per_gene(data, disease_grp, control_grp, value_col)
    if (nrow(paired_b) < 2) return(NA)
    summary_fn(paired_b$disease) / summary_fn(paired_b$control)
  })
  
  data.frame(
    flag      = value_col,
    gene_set  = gene_set_label,
    OR        = ratio,
    OR_low    = quantile(boot_ratios, 0.025, na.rm = TRUE),
    OR_high   = quantile(boot_ratios, 0.975, na.rm = TRUE),
    p_value   = w_test$p.value,
    row.names = NULL
  )
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

# CHANGED: new helper. Tidies the sanity models into one table and BH-adjusts
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

# ==============================================================================
# 4.2  MIXED-EFFECT MODEL   [replaces the previous 4.2]
# ------------------------------------------------------------------------------
# Logistic GLMM with a random intercept per transcript, fitted SEPARATELY within
# each gene set (SNV, FS). The previous version pooled both gene sets into one
# model, which is inconsistent with sections 4.3 and 4.4 and confounds the
# feature effect with systematic differences in PTC position between SNVs and
# frameshift variants.
#
# WHAT THE RANDOM INTERCEPT DOES AND DOES NOT DO
#   It accounts for clustering of variants within transcripts, so standard
#   errors are not inflated by treating correlated variants as independent.
#   It does NOT restrict the comparison to within-gene contrasts: the fixed
#   effect is still estimated using transcripts that appear in only one group,
#   so it mixes between- and within-gene information. For the purely within-gene
#   estimate see the conditional logistic regression in the matched analysis.
#
# MULTIPLICITY
#   BH within each (gene_set, model) stratum, matching 4.3 and 4.4.
#
# DIAGNOSTICS REPORTED PER MODEL
#   n_variants     variants entering the model
#   n_transcripts  transcripts entering the model
#   singular       TRUE if the random-intercept variance collapsed to zero, in
#                  which case the model has degenerated to ordinary logistic
#                  regression and the clustering claim does not hold
#   separation     TRUE if the flag almost perfectly predicts disease status,
#                  in which case the OR is unstable and Firth is used instead
# ==============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

FDR_METHOD <- "BH"

# One flag per Pfam direction: ptc_after and ptc_before are complementary sides
# of the same continuous variable and must not both enter the same family.
FLAG_COLS <- c(
  "variant_ppi_overlap",
  "ptc_before_max_pfam_end",
  "variant_protein_flag",
  "variant_domain_flag",
  "variant_slim_flag",
  "variant_morf_flag",
  "variant_ptm_flag",
  "variant_nls_flag",
  "variant_LCS_flag"
)

# A cell count at or below this in the 2x2 of flag by is_disease triggers the
# Firth-penalised fallback.
SEPARATION_MIN_CELL <- 5


# ------------------------------------------------------------------------------
# 4.2.1  Separation check
# ------------------------------------------------------------------------------

check_separation_2x2 <- function(df, flag, min_cell = SEPARATION_MIN_CELL) {
  tab <- table(df[[flag]], df$is_disease)
  if (any(dim(tab) < 2)) return(list(separated = TRUE, min_cell = 0, tab = tab))
  list(separated = min(tab) <= min_cell, min_cell = min(tab), tab = tab)
}


# ------------------------------------------------------------------------------
# 4.2.2  Fit one flag within one gene set
# ------------------------------------------------------------------------------

fit_one_glmm <- function(dat, flag, gs) {
  
  df <- dat %>%
    dplyr::select(is_disease, all_of(flag), ensembl_transcript_id) %>%
    filter(!is.na(.data[[flag]])) %>%
    mutate(across(all_of(flag), as.integer))
  
  base <- data.frame(
    flag = flag, gene_set = gs,
    n_variants = nrow(df),
    n_transcripts = n_distinct(df$ensembl_transcript_id),
    OR = NA_real_, OR_low = NA_real_, OR_high = NA_real_, p_value = NA_real_,
    min_cell = NA_integer_, separation = NA, singular = NA,
    method = NA_character_, row.names = NULL
  )
  
  if (nrow(df) == 0 ||
      n_distinct(df[[flag]]) < 2 ||
      n_distinct(df$is_disease) < 2) {
    message(sprintf("[%s] %s: insufficient variation, skipped", gs, flag))
    base$method <- "skipped"
    return(base)
  }
  
  sep <- check_separation_2x2(df, flag)
  base$min_cell   <- sep$min_cell
  base$separation <- sep$separated
  
  # --- separation: Firth-penalised logistic, no random effect ----------------
  # logistf cannot carry a random intercept, so the clustering correction is
  # lost. This is flagged in `method` so it is not read as a GLMM result.
  if (sep$separated) {
    message(sprintf("[%s] %s: min cell = %d, using Firth penalisation",
                    gs, flag, sep$min_cell))
    fit <- tryCatch(
      logistf::logistf(as.formula(paste("is_disease ~", flag)), data = df),
      error = function(e) { message(sprintf("  Firth failed: %s", e$message)); NULL }
    )
    if (is.null(fit)) { base$method <- "failed"; return(base) }
    
    ci <- confint(fit)
    base$OR      <- unname(exp(coef(fit)[flag]))
    base$OR_low  <- unname(exp(ci[flag, 1]))
    base$OR_high <- unname(exp(ci[flag, 2]))
    base$p_value <- unname(fit$prob[flag])
    base$method  <- "Firth (no random effect)"
    return(base)
  }
  
  # --- standard path: GLMM ---------------------------------------------------
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


# ------------------------------------------------------------------------------
# 4.2.3  Driver
# ------------------------------------------------------------------------------

run_mixed_effect_flag_analysis <- function(variants_all5,
                                           flag_cols = FLAG_COLS,
                                           gene_sets = c("SNV", "FS")) {
  
  map_dfr(gene_sets, function(gs) {
    dat <- filter(variants_all5, gene_set == gs)
    map_dfr(intersect(flag_cols, names(dat)), ~ fit_one_glmm(dat, .x, gs))
  })
}


# ------------------------------------------------------------------------------
# 4.2.4  Tidy and adjust
# ------------------------------------------------------------------------------
# BH within each gene set. Rows that never produced an estimate are excluded
# from the family so they do not inflate m.

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


# ------------------------------------------------------------------------------
# 4.2.5  Plot
# ------------------------------------------------------------------------------
# Gene sets are separate shapes. Estimates from a singular fit or from the Firth
# fallback are drawn hollow, since in neither case did a random intercept
# actually contribute.

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
                       labels = c(`TRUE` = "GLMM", `FALSE` = "singular / Firth"),
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


# ==============================================================================
# 4.3  HIERARCHICAL BAYESIAN MODEL
# ------------------------------------------------------------------------------
# Partial pooling of flag effects across genes within each gene set (SNV/FS).
#
# NOTE ON MULTIPLICITY: no post-hoc p-value adjustment is applied here by
# design. Multiplicity is handled by partial pooling: the hierarchical prior
# shrinks transcript-level effects toward the group mean, which plays the same
# role as an FDR correction. Applying BH on top of shrinkage would double-count
# the adjustment. Significance is read off the 95% credible interval.
# ==============================================================================

run_bayesian_by_geneset <- function(variants_all5, flag_cols = GENE_MATCHED_FLAGS) {
  library(brms)
  gene_sets <- c("SNV", "FS")
  
  purrr::map(gene_sets, function(gs) {
    dat_gs <- filter(variants_all5, gene_set == gs)   # keeps *_disease and *_control for this set
    
    purrr::map(flag_cols, function(flag) {
      df <- dat_gs %>%
        dplyr::select(is_disease, all_of(flag), ensembl_transcript_id) %>%
        filter(!is.na(.data[[flag]])) %>%
        mutate(across(all_of(flag), as.integer))
      
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

# ==============================================================================
# 4.4  BOOTSTRAP AND GENE-MATCHED ANALYSIS
# ------------------------------------------------------------------------------
# For each gene set (SNV / FS), match disease and control variants by gene and
# compare via paired Wilcoxon on gene-level proportions/medians, with
# bootstrap CIs built by resampling one variant per gene per group.
# Multiplicity: BH within each (gene_set, model) stratum.
# ==============================================================================

run_gene_matched_flag_analysis <- function(variants_all5,
                                           flag_cols = GENE_MATCHED_FLAGS,
                                           cont_col  = CONT_COL,
                                           fdr_group_vars = c("gene_set", "model")) {
  
  gene_sets <- c("SNV", "FS")
  
  results <- map_dfr(gene_sets, function(gs) {
    data          <- filter(variants_all5, gene_set == gs)
    skip_flags    <- flags_without_variation(variants_all5, gs, flag_cols)
    active_flags  <- setdiff(flag_cols, skip_flags)
    
    message(sprintf("Skipping in %s: %s", gs,
                    if (length(skip_flags)) paste(skip_flags, collapse = ", ") else "none"))
    
    unadj <- map_dfr(active_flags, ~ fit_flag_glm(data, .x)) %>%
      mutate(gene_set = gs, model = "Unadjusted")
    unadj_cont <- fit_flag_glm(data, cont_col) %>%
      mutate(gene_set = gs, model = "Unadjusted")
    
    matched <- map_dfr(active_flags, ~ run_gene_matched_one(data, gs, .x, summary_fn = mean)) %>%
      mutate(model = "Gene-matched")
    matched_cont <- run_gene_matched_one(data, gs, cont_col, summary_fn = median) %>%
      mutate(model = "Gene-matched")
    
    bind_rows(unadj, unadj_cont, matched, matched_cont)
  })
  
  results %>%
    add_fdr_sig(group_vars = fdr_group_vars) %>%
    mutate(model = factor(model, levels = c("Unadjusted", "Gene-matched")))
}

plot_gene_matched_flags <- function(combined, labels = FLAG_LABELS) {
  combined <- combined %>%
    mutate(flag = ifelse(flag %in% names(labels), labels[flag], flag))
  plot_or_forest(
    combined, facet_var = "model", shape_var = "gene_set",
    title    = "Disease vs. control variant features",
    subtitle = "Unadjusted logistic regression vs. gene-matched paired Wilcoxon (BH-adjusted p)"
  )
}


# ==============================================================================
# PIPELINE ENTRY POINT
# ==============================================================================

# --- 4.0: data + annotation ---------------------------------------------------
pfam_fin     <- get_pfam_annotations(variants_all2$ensembl_transcript_id, ensembl)
human_1_     <- read_delim(file.path(OTHERS_DIR, "human (1).txt"),
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
dist_results <- tidy_dist_models(dist_models)          # CHANGED: BH-adjusted table
print(dist_results)
write.csv(dist_results, file.path(OUT_DIR, "dist_sanity_models_fdr.csv"), row.names = FALSE)

# --- 4.2: mixed-effect model ---------------------------------------------------
mixed_results <- run_mixed_effect_flag_analysis(variants_all5)
mixed_tidy    <- tidy_mixed_results(mixed_results)     # CHANGED: keep the adjusted table
write.csv(mixed_tidy, file.path(OUT_DIR, "mixed_effect_fdr.csv"), row.names = FALSE)
plot_mixed_effect_flags(mixed_results)

# --- 4.3: hierarchical Bayesian model -------------------------------------------
bayes_by_gs <- run_bayesian_by_geneset(variants_all5, GENE_MATCHED_FLAGS)
plot_bayesian_by_geneset(bayes_by_gs)

# --- 4.4: bootstrap / gene-matched analysis -------------------------------------
gene_matched_results <- run_gene_matched_flag_analysis(variants_all5)
write.csv(gene_matched_results, file.path(OUT_DIR, "gene_matched_fdr.csv"), row.names = FALSE)
plot_gene_matched_flags(gene_matched_results)

# boot_matched <- run_boot_matched_analysis(variants_all5)
# print(as.data.frame(boot_matched), row.names = FALSE)
# write.csv(boot_matched, file.path(OUT_DIR, "boot_matched_4_4.csv"), row.names = FALSE)
# plot_boot_matched(boot_matched)