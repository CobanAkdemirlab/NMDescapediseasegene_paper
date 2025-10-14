library(readr)
library(dplyr)
library(ggplot2)
library(tidyr) 

#1. input vep txt file
plus1_vep_cadd_dbnsfp_tss_prot_nearestexon <- read_table("~/Downloads/vep0813/plus1_vep_cadd_dbnsfp_tss_prot_nearestexon.txt", 
                                                              skip = 142)
plus1_control_vep_cadd_dbnsfp_tss_prot_nearestexon <- read_table("~/Downloads/vep0813/plus1_control_vep_cadd_dbnsfp_tss_prot_nearestexon.txt", 
                                                                      skip = 142)
minus1_vep_cadd_dbnsfp_tss_prot_nearestexon <- read_table("~/Downloads/vep0813/minus1_vep_cadd_dbnsfp_tss_prot_nearestexon.txt",
                                                               skip = 142)
minus1_control_vep_cadd_dbnsfp_tss_prot_nearestexon <- read_table("~/Downloads/vep0813/minus1_control_vep_cadd_dbnsfp_tss_prot_nearestexon.txt",
                                                                     skip = 142)
snv_vep_cadd_dbnsfp_tss_prot_nearestexon <- read_delim("~/Downloads/vep0813/snv_vep_cadd_dbnsfp_tss_prot_nearestexon.txt", 
                                                            delim = "\t", escape_double = FALSE, 
                                                            trim_ws = TRUE, skip = 142)
snv_control_vep_cadd_dbnsfp_tss_prot_nearestexon <- read_delim("~/Downloads/vep0813/snv_control_vep_cadd_dbnsfp_tss_prot_nearestexon.txt", 
                                                                 delim = "\t", escape_double = FALSE, 
                                                                 trim_ws = TRUE, skip = 142)

#2. only keep the rows where variant and gene both match minus1_variants
snv_variants <- read_csv("~/Downloads/0406 3/snv_variants0406.csv")
snv_control_variants <- read_csv("~/Downloads/0406 3/snv_gnomAD_variants0406.csv")
minus1_variants <- read_csv("~/Downloads/0406 3/minus1_variants0406.csv")
minus1_control_variants <- read_csv("~/Downloads/0406 3/minus1_control_gnomAD_variants0406.csv")
plus1_variants <- read_csv("~/Downloads/0406 3/plus1_variants0406.csv")
plus1_control_variants <- read_csv("~/Downloads/0406 3/plus1_control_gnomAD_variants0406.csv")
minus1_vep = minus1_vep_cadd_dbnsfp_tss_prot_nearestexon[(which(minus1_vep_cadd_dbnsfp_tss_prot_nearestexon$Feature %in% minus1_variants$transcript)),]
minus1_control_vep = minus1_control_vep_cadd_dbnsfp_tss_prot_nearestexon[(which(minus1_control_vep_cadd_dbnsfp_tss_prot_nearestexon$Feature %in% minus1_control_variants$transcript)),]
plus1_vep = plus1_vep_cadd_dbnsfp_tss_prot_nearestexon[(which(plus1_vep_cadd_dbnsfp_tss_prot_nearestexon$Feature %in% plus1_variants$transcript)),]
plus1_control_vep = plus1_control_vep_cadd_dbnsfp_tss_prot_nearestexon[(which(plus1_control_vep_cadd_dbnsfp_tss_prot_nearestexon$Feature %in% plus1_control_variants$transcript)),]
snv_vep = snv_vep_cadd_dbnsfp_tss_prot_nearestexon[(which(snv_vep_cadd_dbnsfp_tss_prot_nearestexon$Feature %in% snv_variants$transcript)),]
snv_control_vep = snv_control_vep_cadd_dbnsfp_tss_prot_nearestexon[(which(snv_control_vep_cadd_dbnsfp_tss_prot_nearestexon$Feature %in% snv_control_variants$transcript)),]

#only keep the AD variants
AD_gene_ID = getBM(attributes=c("ensembl_gene_id",'hgnc_symbol'),
                   filters = "hgnc_symbol", 
                   values = pli_AD_genes$gene,
                   mart = ensembl)

minus1_vep_AD = minus1_vep[which(minus1_vep$Gene %in% AD_gene_ID$ensembl_gene_id),]
minus1_control_vep_AD = minus1_control_vep[which(minus1_control_vep$Gene %in% AD_gene_ID$ensembl_gene_id),]
plus1_vep_AD = plus1_vep[which(plus1_vep$Gene %in% AD_gene_ID$ensembl_gene_id),]
plus1_control_vep_AD = plus1_control_vep[which(plus1_control_vep$Gene %in% AD_gene_ID$ensembl_gene_id),]
snv_vep_AD = snv_vep[which(snv_vep$Gene %in% AD_gene_ID$ensembl_gene_id),]
snv_control_vep_AD = snv_control_vep[which(snv_control_vep$Gene %in% AD_gene_ID$ensembl_gene_id),]
invalid_vals <- c("-", ".", ".,.", ".,.,.", "invalid_field", "", "NA", "NaN", "nan")

to_num <- function(x){
  x <- suppressWarnings(as.character(x))
  x <- ifelse(is.na(x) | x %in% invalid_vals, NA, x)
  suppressWarnings(as.numeric(x))
}

parse_dbnsfp_numeric <- function(x){
  x <- suppressWarnings(as.character(x))
  x[is.na(x) | x %in% invalid_vals] <- NA_character_
  parts <- strsplit(x, "[,;|]")
  sapply(parts, function(v){
    if (!length(v)) return(NA_real_)
    v <- trimws(v); v <- v[nzchar(v) & !(v %in% invalid_vals)]
    nums <- suppressWarnings(as.numeric(v))
    if (!length(nums) || all(is.na(nums))) NA_real_ else max(nums, na.rm = TRUE)
  })
}

parse_nearest_exon_dist <- function(x){
  x <- suppressWarnings(as.character(x))
  suppressWarnings(as.numeric(vapply(strsplit(x, "\\|"), function(v){
    if (length(v) >= 4) v[[4]] else NA_character_
  }, character(1))))
}

len_change_from_aa <- function(s){
  s <- suppressWarnings(as.character(s))
  s[is.na(s) | s %in% c("", "-")] <- NA_character_
  ref <- sub("/.*$", "", s); alt <- sub("^.*/", "", s)
  ref_len <- nchar(gsub("-", "", ref), allowNA = TRUE)
  alt_len <- nchar(gsub("-", "", alt), allowNA = TRUE)
  alt_len - ref_len
}

coalesce_sift <- function(df){
  out <- rep(NA_real_, nrow(df))
  for (cl in c("SIFT_score","SIFT4G_score")) if (cl %in% names(df)) {
    val <- parse_dbnsfp_numeric(df[[cl]])
    idx <- is.na(out) & !is.na(val); out[idx] <- val[idx]
  }
  out
}

NUMERIC_DBNSFP <- c(
  "REVEL_score","PrimateAI_score","DANN_score","ClinPred_score",
  "BayesDel_addAF_score","BayesDel_noAF_score","LIST_S2_score","M_CAP_score",
  "MetaLR_score","MetaSVM_score","MutationTaster_score","MutationAssessor_score",
  "FATHMM_score","PROVEAN_score","LRT_score",
  "Polyphen2_HDIV_score","Polyphen2_HVAR_score","SIFT4G_score",
  "GERP++_RS","GERP++_NR","phyloP100way_vertebrate","phyloP30way_mammalian",
  "phyloP17way_primate","phastCons100way_vertebrate","phastCons30way_mammalian",
  "phastCons17way_primate","SiPhy_29way_logOdds",
  "dbscSNV_ADA_SCORE","dbscSNV_RF_SCORE"
)

clean_vep <- function(df){
  out <- df %>%
    mutate(
      # distances
      tss_dist      = to_num(TSSDistance),
      tss_dist_abs  = abs(tss_dist),
      nearest_exon_dist     = parse_nearest_exon_dist(NearestExonJB),
      nearest_exon_dist_abs = abs(nearest_exon_dist),
      # CADD + SIFT
      cadd_phred = to_num(CADD_PHRED),
      cadd_raw   = to_num(CADD_RAW),
      sift_score = coalesce_sift(df),
      # protein length change
      prot_len_change = len_change_from_aa(Amino_acids)
    )
  for (nm in NUMERIC_DBNSFP) if (nm %in% names(out)) out[[nm]] <- parse_dbnsfp_numeric(out[[nm]])
  out
}

# Build from your filtered objects in #2
groups_list <- list(
  minus1         = minus1_vep,
  minus1_control = minus1_control_vep,
  plus1          = plus1_vep,
  plus1_control  = plus1_control_vep,
  SNV            = snv_vep,
  SNV_control    = snv_control_vep
)
AD_groups_list <- list(
  minus1         = minus1_vep_AD,
  minus1_control = minus1_control_vep_AD,
  plus1          = plus1_vep_AD,
  plus1_control  = plus1_control_vep_AD,
  SNV            = snv_vep_AD,
  SNV_control    = snv_control_vep_AD
)

cleaned <- lapply(groups_list, clean_vep)
AD_cleaned <- lapply(AD_groups_list, clean_vep)

# Combine; drop Location to avoid hms/character clashes
all_clean <- dplyr::bind_rows(lapply(names(cleaned), function(nm){
  x <- cleaned[[nm]]; x$Group <- nm; dplyr::select(x, -any_of("Location"))
}))
AD_clean <- dplyr::bind_rows(lapply(names(AD_cleaned), function(nm){
  x <- AD_cleaned[[nm]]; x$Group <- nm; dplyr::select(x, -any_of("Location"))
}))
# Group order + colors
all_clean$Group <- factor(all_clean$Group,
                          levels = c("minus1","minus1_control","plus1","plus1_control","SNV","SNV_control")
)
AD_clean$Group <- factor(AD_clean$Group,
                         levels = c("minus1","minus1_control","plus1","plus1_control","SNV","SNV_control")
)
group_colors <- c(
  "minus1"="#1f77b4","minus1_control"="#aec7e8",
  "plus1" ="#ff7f0e","plus1_control" ="#ffbb78",
  "SNV"   ="#2ca02c","SNV_control"   ="#98df8a"
)

theme_compact <- theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Case–control pairs for p-values
pairs_base <- list(
  c("minus1","minus1_control"),
  c("plus1","plus1_control"),
  c("SNV","SNV_control")
)
valid_comparisons <- function(dat){
  present <- unique(dat$Group)
  Filter(function(cp) all(cp %in% present), pairs_base)
}

# Small helper to draw violin+box with optional log10 and p-values
vbox_with_p <- function(df, feature, ylab, log10_y = FALSE, title = feature){
  stopifnot(feature %in% names(df))
  dat <- df %>% transmute(Group, value = .data[[feature]]) %>% filter(!is.na(value))
  if (log10_y) dat <- dat %>% filter(value > 0)
  comps <- valid_comparisons(dat)
  
  p <- ggplot(dat, aes(Group, value, fill = Group)) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.14, outlier.size = 0.4, color = "black") +
    scale_fill_manual(values = group_colors, guide = "none") +
    labs(x = NULL, y = ylab, title = title) +
    theme_compact
  if (log10_y) p <- p + scale_y_continuous(trans = "log10")
  
  if (length(comps)){
    p <- p +
      ggpubr::stat_compare_means(
        comparisons = comps, method = "wilcox.test",
        label = "p.format", hide.ns = TRUE, label.y.npc = 0.98, size = 3
      ) +
      coord_cartesian(clip = "off") +
      scale_y_continuous(expand = expansion(mult = if (log10_y) c(0.05, 0.35) else c(0.05, 0.35)))
  }
  p
}

dir.create("feature_plots", showWarnings = FALSE)



#3. plot tss diatance
p_tss <- vbox_with_p(all_clean, "tss_dist_abs",
                     ylab = "Distance to TSS (bp, log10)", log10_y = TRUE,
                     title = "TSS distance")
p_tss_AD = vbox_with_p(AD_clean, "tss_dist_abs",
                         ylab = "Distance to TSS (bp, log10)", log10_y = TRUE,
                         title = "TSS distance in AD genes")
print(p_tss)
print(p_tss_AD)
ggsave("feature_plots/tss_distance_violin_pvals.pdf", p_tss, width = 7.2, height = 7.8)
ggsave("feature_plots/tss_distance_AD_genes_violin_pvals.pdf", p_tss_AD, width = 7.2, height = 7.8)

#4. plot dbNSFP related featureslibrary(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(scales)

# ---- data prep ----
feat_linear <- intersect(unique(c("cadd_phred","cadd_raw","sift_score", NUMERIC_DBNSFP)),
                         names(all_clean))

long_lin <- all_clean %>%
  select(Group, all_of(feat_linear)) %>%
  pivot_longer(-Group, names_to = "feature", values_to = "value") %>%
  filter(!is.na(value))
AD_long_lin <- AD_clean %>%
  select(Group, all_of(feat_linear)) %>%
  pivot_longer(-Group, names_to = "feature", values_to = "value") %>%
  filter(!is.na(value))

have_any <- long_lin %>%
  group_by(feature) %>%
  summarise(n = sum(!is.na(value)), .groups = "drop") %>%
  filter(n > 0) %>% pull(feature)
AD_have_any <- AD_long_lin %>%
  group_by(feature) %>%
  summarise(n = sum(!is.na(value)), .groups = "drop") %>%
  filter(n > 0) %>% pull(feature)

long_lin <- long_lin %>% filter(feature %in% have_any)
AD_long_lin <- AD_long_lin %>% filter(feature %in% AD_have_any)

# ---- comparisons ----
pairs_base <- list(c("SNV","SNV_Control"))

# ---- theme ----
title_sz      <- 20
axis_title_sz <- 16
axis_text_sz  <- 14
strip_sz      <- 14
pval_sz       <- 4.2

theme_compact_big <- theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1, size = axis_text_sz),
    axis.text.y  = element_text(size = axis_text_sz),
    axis.title.x = element_text(size = axis_title_sz),
    axis.title.y = element_text(size = axis_title_sz),
    strip.text   = element_text(size = strip_sz, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold", size = title_sz)
  )

# ---- pre-compute Wilcoxon p-values per feature ----
pvals <- long_lin %>%
  group_by(feature) %>%
  wilcox_test(value ~ Group, comparisons = pairs_base) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(p_label = rstatix::p_format(p))  # numeric formatted p
AD_pvals = AD_long_lin %>%
  group_by(feature) %>%
  wilcox_test(value ~ Group, comparisons = pairs_base) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(p_label = rstatix::p_format(p))  # numeric formatted p

# assign y positions above the data range in each facet
facet_tops <- long_lin %>%
  group_by(feature) %>%
  summarise(ymax = max(value, na.rm = TRUE), .groups = "drop")
AD_facet_tops = AD_long_lin %>%
  group_by(feature) %>%
  summarise(ymax = max(value, na.rm = TRUE), .groups = "drop")

pvals <- pvals %>%
  left_join(facet_tops, by = "feature") %>%
  group_by(feature) %>%
  mutate(y.position = ymax * 1.05) %>%
  ungroup()
AD_pvals = AD_pvals %>%
  left_join(AD_facet_tops, by = "feature") %>%
  group_by(feature) %>%
  mutate(y.position = ymax * 1.05) %>%
  ungroup()

# ---- plot ----
p_dbnsfp <- ggplot(long_lin, aes(Group, value, fill = Group)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.14, outlier.size = 0.6, color = "black") +
  facet_wrap(~ feature, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.45))) +
  labs(title = "dbNSFP features across variant groups (linear scale)",
       x = NULL, y = "Value") +
  theme_compact_big +
  coord_cartesian(clip = "off") +
  stat_pvalue_manual(
    data        = pvals,
    label       = "p_label",
    xmin        = "group1",
    xmax        = "group2",
    y.position  = "y.position",
    tip.length  = 0.01,
    size        = pval_sz,
    inherit.aes = FALSE
  )
AD_p_dbnsfp = ggplot(AD_long_lin, aes(Group, value, fill = Group)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.14, outlier.size = 0.6, color = "black") +
  facet_wrap(~ feature, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.45))) +
  labs(title = "dbNSFP features across variant groups in AD genes (linear scale)",
       x = NULL, y = "Value") +
  theme_compact_big +
  coord_cartesian(clip = "off") +
  stat_pvalue_manual(
    data        = AD_pvals,
    label       = "p_label",
    xmin        = "group1",
    xmax        = "group2",
    y.position  = "y.position",
    tip.length  = 0.01,
    size        = pval_sz,
    inherit.aes = FALSE
  )

print(p_dbnsfp)
print(AD_p_dbnsfp)

ggsave("feature_plots/dbnsfp_cadd_linear_facets_pvals.pdf", p_dbnsfp, width = 16, height = 12)
ggsave("feature_plots/dbnsfp_cadd_linear_AD_genes_facets_pvals.pdf", AD_p_dbnsfp, width = 16, height = 12)
#5. plot protein length change
all_clean_plc <- all_clean %>%
  mutate(abs_plc = abs(prot_len_change)) %>%
  filter(!is.na(abs_plc), abs_plc > 0)
AD_all_clean_plc = AD_clean %>%
  mutate(abs_plc = abs(prot_len_change)) %>%
  filter(!is.na(abs_plc), abs_plc > 0)

p_plc <- vbox_with_p(all_clean_plc, "abs_plc",
                     ylab = "|Δ protein length| (AA, log10)",
                     log10_y = TRUE, title = "Protein length change")
plc_AD = vbox_with_p(AD_all_clean_plc, "abs_plc",
                        ylab = "|Δ protein length| (AA, log10)",
                        log10_y = TRUE, title = "Protein length change in AD genes")
print(p_plc)
print(plc_AD)
ggsave("feature_plots/protein_length_change_log_violin_pvals.pdf", p_plc, width = 9.5, height = 6.5)
ggsave("feature_plots/protein_length_change_AD_genes_log_violin_pvals.pdf", plc_AD, width = 9.5, height = 6.5)
#6. plot NearestExonJB_distance 
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

# --- parser: pick the entry with the smallest |distance| ---
parse_nearest_exonjb_dist <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  parts <- str_split(x, "&", simplify = FALSE)[[1]]
  df <- map_dfr(parts, function(p) {
    toks <- str_split(p, "\\|", simplify = TRUE)
    if (ncol(toks) < 4) {
      toks <- cbind(toks, matrix(NA_character_, nrow = nrow(toks), ncol = 4 - ncol(toks)))
    }
    tibble(
      exon_id  = toks[,1],
      distance = suppressWarnings(as.numeric(toks[,2])),
      boundary = toks[,3],
      length   = suppressWarnings(as.numeric(toks[,4]))
    )
  }) %>% mutate(dist_abs = abs(distance)) %>% arrange(dist_abs)
  
  df$dist_abs[1] %>% as.numeric()
}

# vectorized wrapper that returns a tibble
vec_nearest_dist <- function(v) tibble(
  nearest_exon_dist_abs_calc = map_dbl(v, parse_nearest_exonjb_dist)
)

# ---- Compute column in both data frames (only fill if missing) ----
if (!("nearest_exon_dist_abs" %in% names(all_clean))) {
  all_clean <- all_clean %>%
    bind_cols(vec_nearest_dist(.$NearestExonJB)) %>%
    mutate(nearest_exon_dist_abs = coalesce(nearest_exon_dist_abs, nearest_exon_dist_abs_calc)) %>%
    select(-nearest_exon_dist_abs_calc)
} else if (!("NearestExonJB" %in% names(all_clean))) {
  warning("all_clean lacks NearestExonJB; using existing nearest_exon_dist_abs as-is.")
}

if (!("nearest_exon_dist_abs" %in% names(AD_clean))) {
  AD_clean <- AD_clean %>%
    bind_cols(vec_nearest_dist(.$NearestExonJB)) %>%
    mutate(nearest_exon_dist_abs = coalesce(nearest_exon_dist_abs, nearest_exon_dist_abs_calc)) %>%
    select(-nearest_exon_dist_abs_calc)
} else if (!("NearestExonJB" %in% names(AD_clean))) {
  warning("AD_clean lacks NearestExonJB; using existing nearest_exon_dist_abs as-is.")
}

# ---- Plot: overall + AD genes ----
p_nearest <- vbox_with_p(
  all_clean, "nearest_exon_dist_abs",
  ylab = "Distance to nearest exon (bp, log10)",
  log10_y = TRUE,
  title = "Nearest exon distance"
)

AD_p_nearest <- vbox_with_p(
  AD_clean, "nearest_exon_dist_abs",
  ylab = "Distance to nearest exon (bp, log10)",
  log10_y = TRUE,
  title = "Nearest exon distance in AD genes"
)

# Optionally print/show
print(p_nearest)
print(AD_p_nearest)
#also plot the distance to exon end, if boundary type is start, use length - distance, if end, use distance
# helper to parse one NearestExonJB string; returns a list(dist_abs, to_end_abs)
parse_nearest_exonjb <- function(x) {
  if (is.na(x) || x == "") return(list(dist_abs = NA_real_, to_end_abs = NA_real_))
  
  # split multiple entries (e.g., "ex1|12|start|150&ex2|7|end|120")
  parts <- str_split(x, "&", simplify = FALSE)[[1]]
  
  # parse each part into tibble rows
  df <- map_dfr(parts, function(p) {
    toks <- str_split(p, "\\|", simplify = TRUE)
    # pad if fewer than 4 tokens
    if (ncol(toks) < 4) {
      toks <- cbind(toks, matrix(NA_character_, nrow = nrow(toks), ncol = 4 - ncol(toks)))
    }
    tibble(
      exon_id  = toks[,1],
      distance = suppressWarnings(as.numeric(toks[,2])),
      boundary = toks[,3],
      length   = suppressWarnings(as.numeric(toks[,4]))
    )
  })
  
  # prefer row with smallest |distance|
  df <- df %>%
    mutate(dist_abs = abs(distance)) %>%
    arrange(dist_abs)
  
  best <- df[1, , drop = FALSE]
  
  # distance to exon end:
  # - if boundary == "start": length - |distance|
  # - if boundary == "end"  : |distance|
  to_end <- dplyr::case_when(
    !is.na(best$boundary) & tolower(best$boundary) == "start" ~ best$length - best$dist_abs,
    !is.na(best$boundary) & tolower(best$boundary) == "end"   ~ best$dist_abs,
    TRUE ~ NA_real_
  )
  
  # clamp negatives to 0 (just in case distance slightly exceeds length)
  to_end <- ifelse(is.na(to_end), NA_real_, pmax(to_end, 0))
  
  list(dist_abs = best$dist_abs %>% as.numeric(), to_end_abs = as.numeric(to_end))
}

# vectorized wrapper
vec_parse_nearest <- function(v) {
  res <- map(v, parse_nearest_exonjb)
  tibble(
    nearest_exon_dist_abs_calc = map_dbl(res, "dist_abs"),
    nearest_exon_to_end_abs    = map_dbl(res, "to_end_abs")
  )
}

# ---- Add columns to your data ----
# Assumes your raw annotation column is named "NearestExonJB"
# If you already have nearest_exon_dist_abs, we keep it and also compute a 'calc' to compare.
all_clean <- all_clean %>%
  bind_cols(vec_parse_nearest(.$NearestExonJB)) %>%
  mutate(
    # prefer existing column if present; otherwise use calculated
    nearest_exon_dist_abs = dplyr::coalesce(nearest_exon_dist_abs, nearest_exon_dist_abs_calc)
  ) %>%
  select(-nearest_exon_dist_abs_calc)

AD_clean <- AD_clean %>%
  bind_cols(vec_parse_nearest(.$NearestExonJB)) %>%
  mutate(
    nearest_exon_dist_abs = dplyr::coalesce(nearest_exon_dist_abs, nearest_exon_dist_abs_calc)
  ) %>%
  select(-nearest_exon_dist_abs_calc)

# ---- Your plots ----
# existing
p_nearest <- vbox_with_p(
  all_clean, "nearest_exon_dist_abs",
  ylab = "Distance to nearest exon (bp, log10)",
  log10_y = TRUE, title = "Nearest exon distance"
)

AD_p_nearest <- vbox_with_p(
  AD_clean, "nearest_exon_dist_abs",
  ylab = "Distance to nearest exon (bp, log10)",
  log10_y = TRUE, title = "Nearest exon distance in AD genes"
)

# new: distance to exon END
p_nearest_end <- vbox_with_p(
  all_clean, "nearest_exon_to_end_abs",
  ylab = "Distance to exon end (bp, log10)",
  log10_y = TRUE, title = "Distance to exon end"
)

AD_p_nearest_end <- vbox_with_p(
  AD_clean, "nearest_exon_to_end_abs",
  ylab = "Distance to exon end (bp, log10)",
  log10_y = TRUE, title = "Distance to exon end in AD genes"
)


print(p_nearest)
print(AD_p_nearest)
#in AD_clean, add a flag to indicate if the mutation is the last exon or not
##add cds location of the start of last exon by gene id using biomart(getBM), compare this with CDS_position
##separate plots for flag=T/F
suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

# ---- 1) Pull exon-level coding segments via biomaRt ----
# If you need a specific Ensembl release, set 'version=' (e.g., 105) in useEnsembl().
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

gene_ids <- unique(AD_clean$Gene)

exon_in = getBM(attributes = c('ensembl_gene_id','ensembl_transcript_id'),
                 filters = 'ensembl_gene_id', values = gene_ids, mart = mart)
exon_in_can = getBM(attributes =c('ensembl_transcript_id','transcript_is_canonical'),
                    filters = 'ensembl_transcript_id',
                    values = exon_in$ensembl_transcript_id,
                    mart = mart
                    )
exon_in_can = exon_in_can[which(exon_in_can$transcript_is_canonical==1),]
exon_in_can_cds = getBM(attributes = c('ensembl_transcript_id','cds_start','cds_end','rank'),
                        filters = 'ensembl_transcript_id',
                        values = exon_in_can$ensembl_transcript_id,
                        mart = mart
                        )
  
exon_last = exon_in_can_cds  %>%
  group_by(ensembl_transcript_id) %>%
  #delete rows with NA
  filter(!is.na(cds_start) & !is.na(cds_end)) %>%
  filter(rank == max(rank)) %>%
  mutate(last_exon_cds_start = cds_start, last_exon_cds_end = cds_end)
exon_last$gene_id = exon_in$ensembl_gene_id[match(exon_last$ensembl_transcript_id, exon_in$ensembl_transcript_id)]

AD_clean$last_exon_start = exon_last$last_exon_cds_start[match(AD_clean$Gene, exon_last$gene_id)]
AD_clean$in_last_exon = ifelse(!is.na(AD_clean$last_exon_start) & !is.na(AD_clean$CDS_position) & AD_clean$CDS_position >= AD_clean$last_exon_start, TRUE, FALSE)

AD_last_true  <- AD_clean %>% filter(in_last_exon %in% TRUE)
AD_last_false <- AD_clean %>% filter(in_last_exon %in% FALSE)
if (!("nearest_exon_to_end_abs" %in% names(AD_clean))) {
  parse_nearest_exonjb <- function(x) {
    if (is.na(x) || x == "") return(list(dist_abs = NA_real_, to_end_abs = NA_real_))
    parts <- str_split(x, "&", simplify = FALSE)[[1]]
    df <- map_dfr(parts, function(p) {
      toks <- str_split(p, "\\|", simplify = TRUE)
      if (ncol(toks) < 4) toks <- cbind(toks, matrix(NA_character_, nrow = nrow(toks), ncol = 4 - ncol(toks)))
      tibble(
        exon_id  = toks[,1],
        distance = suppressWarnings(as.numeric(toks[,2])),
        boundary = toks[,3],
        length   = suppressWarnings(as.numeric(toks[,4]))
      )
    }) %>% mutate(dist_abs = abs(distance)) %>% arrange(dist_abs)
    
    best <- df[1, , drop = FALSE]
    to_end <- dplyr::case_when(
      !is.na(best$boundary) & tolower(best$boundary) == "start" ~ best$length - best$dist_abs,
      !is.na(best$boundary) & tolower(best$boundary) == "end"   ~ best$dist_abs,
      TRUE ~ NA_real_
    )
    to_end <- ifelse(is.na(to_end), NA_real_, pmax(to_end, 0))
    list(dist_abs = as.numeric(best$dist_abs), to_end_abs = as.numeric(to_end))
  }
  
  res <- map(AD_clean$NearestExonJB, parse_nearest_exonjb)
  AD_clean$nearest_exon_to_end_abs <- map_dbl(res, "to_end_abs")
}

yvar_to_plot <- "nearest_exon_to_end_abs"

AD_p_last_true_end <- vbox_with_p(
  AD_last_true, yvar_to_plot,
  ylab = "Distance to exon end (bp, log10)",
  log10_y = TRUE,
  title = "AD genes: distance to exon end (in_last_exon = TRUE)"
)

AD_p_last_false_end <- vbox_with_p(
  AD_last_false, yvar_to_plot,
  ylab = "Distance to exon end (bp, log10)",
  log10_y = TRUE,
  title = "AD genes: distance to exon end (in_last_exon = FALSE)"
)

print(AD_p_last_true_end)
print(AD_p_last_false_end)

yvar_to_plot <- "nearest_exon_dist_abs"

# Plot for variants IN last exon
AD_p_last_true_near <- vbox_with_p(
  AD_last_true, yvar_to_plot,
  ylab = "Nearest exon distance (bp, log10)",
  log10_y = TRUE,
  title = "AD genes: Nearest exon distance (in_last_exon = TRUE)"
)

# Plot for variants NOT in last exon
AD_p_last_false_near <- vbox_with_p(
  AD_last_false, yvar_to_plot,
  ylab = "Nearest exon distance (bp, log10)",
  log10_y = TRUE,
  title = "AD genes: Nearest exon distance (in_last_exon = FALSE)"
)

print(AD_p_last_true_near)
print(AD_p_last_false_near)

print(AD_p_last_true)
print(AD_p_last_false)



ggsave("feature_plots/nearest_exon_distance_violin_pvals.pdf", p_nearest, width = 7.2, height = 7.8)
ggsave("feature_plots/nearest_exon_distance_AD_genes_violin_pvals.pdf", AD_p_nearest, width = 7.2, height = 7.8)
#7. plot ptc to exonjunciton distance

# --- #7a. PTC–variant distance for all 6 groups (with consistent group names) ----

# Helper: process any *_dis data frame into (Group, ptc_variant_dist)
process_dis <- function(dis_df, group_name, variant_keys = NULL) {
  if (is.null(dis_df)) return(NULL)
  df <- dis_df
  if (!is.null(variant_keys)) {
    if (!("Variant_Key" %in% names(df))) return(NULL)
    df <- df[df$Variant_Key %in% variant_keys, , drop = FALSE]
  }
  if (!all(c("Variant_Key","Distance") %in% names(df))) return(NULL)
  df |>
    dplyr::transmute(Group = group_name, ptc_variant_dist = Distance) |>
    dplyr::filter(is.finite(ptc_variant_dist), ptc_variant_dist > 0)
}

# Safe getter
get_or_null <- function(name) if (exists(name, inherits = TRUE)) get(name) else NULL

# Collect available *_dis data frames
minus1_dis         <- get_or_null("minus1_dis")
minus1_control_dis <- get_or_null("minus1_control_dis")
plus1_dis          <- get_or_null("plus1_dis")
plus1_control_dis  <- get_or_null("plus1_control_dis")
snv_dis            <- get_or_null("snv_dis")
snv_control_dis    <- get_or_null("snv_control_dis")

# Build per-group tables (use your variant key vectors if present)
ptc_var_list <- list(
  process_dis(minus1_dis,         "minus1",         if (exists("minus1_variants")) minus1_variants$key else NULL),
  process_dis(minus1_control_dis, "minus1_Control", if (exists("minus1_control_variants")) minus1_control_variants$key else NULL),
  process_dis(plus1_dis,          "plus1",          if (exists("plus1_variants")) plus1_variants$key else NULL),
  process_dis(plus1_control_dis,  "plus1_Control",  if (exists("plus1_control_variants")) plus1_control_variants$key else NULL),
  process_dis(snv_dis,            "SNV",            if (exists("snv_variants")) snv_variants$key else NULL),
  process_dis(snv_control_dis,    "SNV_Control",    if (exists("snv_control_variants")) snv_control_variants$key else NULL)
)

ptc_variant_all <- dplyr::bind_rows(Filter(Negate(is.null), ptc_var_list))

# Define colors & comparisons (names MUST match Group exactly)
group_colors <- c(
  "minus1"         = "#1f77b4",
  "minus1_Control" = "#aec7e8",
  "plus1"          = "#ff7f0e",
  "plus1_Control"  = "#ffbb78",
  "SNV"            = "#2ca02c",
  "SNV_Control"    = "#98df8a"
)

pairs_base <- list(
  c("minus1","minus1_Control"),
  c("plus1","plus1_Control"),
  c("SNV","SNV_Control")
)

valid_comparisons <- function(dat){
  present <- unique(dat$Group)
  Filter(function(cp) all(cp %in% present), pairs_base)
}

# Reuse your violin+box helper
p_ptc_variant_all <- vbox_with_p(
  df      = ptc_variant_all,
  feature = "ptc_variant_dist",
  ylab    = "PTC–variant distance (bp, log10)",
  log10_y = TRUE,
  title   = "PTC–variant distance (all groups)"
)

print(p_ptc_variant_all)
dir.create("feature_plots", showWarnings = FALSE)
ggsave("feature_plots/ptc_variant_distance_all_groups_log_violin_pvals.pdf",
       p_ptc_variant_all, width = 9.5, height = 6.5)

# Save table for reproducibility
readr::write_csv(ptc_variant_all, "feature_plots/ptc_variant_distance_all_groups_table.csv")


