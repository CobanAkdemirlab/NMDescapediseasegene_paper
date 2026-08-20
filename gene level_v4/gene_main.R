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
# --- 路径解析层（自动插入）------------------------------------------------
# 原来这里是旧机器 /Users/jxu14/ 的绝对路径，换机器必断。改走 paths.R：
#   data_file("x.csv")  按文件名定位，找不到会明确报错
#   out_file("y.csv")   输出到 NMDESC_OUT（默认 ~/Desktop/NMDesc_out）
#   data_root("clinvar") 需要目录而非文件时用
.p <- c("gene level_v3/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../lib/paths.R", "../../../../lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p)) stop("找不到 paths.R —— 请从仓库根目录运行 R")
source(.p[1]); rm(.p)
# --------------------------------------------------------------------------

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
  path_touni          = data_file("NIHMS1818854-supplement-2(A).csv"),
  path_motif          = data_file("NIHMS1818854-supplement-2(B).csv"),
  path_lcs            = data_file("Copy of NIHMS1818854-supplement-2.xls"),
  ppi_file            = "~/Downloads/human (1) (2).txt",
  gtex_path           = data_file("GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct",
                                must = FALSE),   # 可选：缺则 tau=NA
  lof_metrics_path    = data_file("gnomad.v2.1.1.lof_metrics.by_gene.txt"),
  
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
# === 函数库加载（替代原来内联的 10 个函数定义）=========================
# 这些函数现在只有一份权威实现，在 features/functions/ 下。
# 改那里就生效 —— 不再需要同步两处。
.fn_dir <- c("gene level_v3/features/functions", "features/functions", "functions")
.fn_dir <- .fn_dir[dir.exists(.fn_dir)]
if (!length(.fn_dir)) stop("找不到 features/functions/ —— 请从仓库根目录运行 R")
for (.f in list.files(.fn_dir[1], pattern = "\\.R$", full.names = TRUE))
  source(.f)
rm(.f, .fn_dir)
# ======================================================================



# fn: build_gene_all -----------------------------------------------------


# fn: calculate_ppi_degree_centrality ------------------------------------


# fn: plot_gc_content ----------------------------------------------------


# fn: plot_repeat_content ------------------------------------------------


# fn: annotate_motif_flags -----------------------------------------------


# fn: run_pfam_overlap_analysis ------------------------------------------


# fn: run_ppi_overlap_analysis -------------------------------------------


# fn: run_tau_analysis ---------------------------------------------------


# fn: plot_gene_level_features -------------------------------------------


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
