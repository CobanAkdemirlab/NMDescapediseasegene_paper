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
minus1_vep = minus1_vep_cadd_dbnsfp_tss_prot_nearestexon[(which(minus1_vep_cadd_dbnsfp_tss_prot_nearestexon$Feature %in% minus1_variants$transcript)),]
minus1_control_vep = minus1_control_vep_cadd_dbnsfp_tss_prot_nearestexon[(which(minus1_control_vep_cadd_dbnsfp_tss_prot_nearestexon$Feature %in% minus1_control_variants$transcript)),]
plus1_vep = plus1_vep_cadd_dbnsfp_tss_prot_nearestexon[(which(plus1_vep_cadd_dbnsfp_tss_prot_nearestexon$Feature %in% plus1_variants$transcript)),]
plus1_control_vep = plus1_control_vep_cadd_dbnsfp_tss_prot_nearestexon[(which(plus1_control_vep_cadd_dbnsfp_tss_prot_nearestexon$Feature %in% plus1_control_variants$transcript)),]
snv_vep = snv_vep_cadd_dbnsfp_tss_prot_nearestexon[(which(snv_vep_cadd_dbnsfp_tss_prot_nearestexon$Feature %in% snv_variants$transcript)),]
snv_control_vep = snv_control_vep_cadd_dbnsfp_tss_prot_nearestexon[(which(snv_control_vep_cadd_dbnsfp_tss_prot_nearestexon$Feature %in% snv_control_variants$transcript)),]

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

cleaned <- lapply(groups_list, clean_vep)

# Combine; drop Location to avoid hms/character clashes
all_clean <- dplyr::bind_rows(lapply(names(cleaned), function(nm){
  x <- cleaned[[nm]]; x$Group <- nm; dplyr::select(x, -any_of("Location"))
}))

# Group order + colors
all_clean$Group <- factor(all_clean$Group,
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
print(p_tss)
ggsave("feature_plots/tss_distance_violin_pvals.pdf", p_tss, width = 7.2, height = 7.8)


#4. plot dbNSFP related features
feat_linear <- intersect(unique(c("cadd_phred","cadd_raw","sift_score", NUMERIC_DBNSFP)),
                         names(all_clean))

long_lin <- all_clean %>%
  select(Group, all_of(feat_linear)) %>%
  pivot_longer(-Group, names_to = "feature", values_to = "value") %>%
  filter(!is.na(value))

# Keep facets where BOTH groups exist for each case–control (otherwise p-values can’t be computed).
# We’ll draw p-values for any pair that exists; facets remain even if one pair lacks data.
have_any <- long_lin %>%
  group_by(feature) %>%
  summarise(n = sum(!is.na(value)), .groups = "drop") %>%
  filter(n > 0) %>% pull(feature)

long_lin <- long_lin %>% filter(feature %in% have_any)

p_dbnsfp <- ggplot(long_lin, aes(Group, value, fill = Group)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.14, outlier.size = 0.3, color = "black") +
  facet_wrap(~ feature, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, guide = "none") +
  # extra top space so p-value labels don’t get hidden
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(title = "CADD + dbNSFP (numeric) features across groups (linear scale)",
       x = NULL, y = "Value") +
  theme_compact +
  coord_cartesian(clip = "off") +
  ggpubr::stat_compare_means(
    comparisons = pairs_base,
    method = "wilcox.test",
    label = "p.format",
    hide.ns = TRUE,
    label.y.npc = 0.98,
    size = 2.7
  )
print(p_dbnsfp)
ggsave("feature_plots/dbnsfp_cadd_linear_facets_pvals.pdf", p_dbnsfp, width = 16, height = 12)

#5. plot protein length change
all_clean_plc <- all_clean %>%
  mutate(abs_plc = abs(prot_len_change)) %>%
  filter(!is.na(abs_plc), abs_plc > 0)

p_plc <- vbox_with_p(all_clean_plc, "abs_plc",
                     ylab = "|Δ protein length| (AA, log10)",
                     log10_y = TRUE, title = "Protein length change")
print(p_plc)
ggsave("feature_plots/protein_length_change_log_violin_pvals.pdf", p_plc, width = 9.5, height = 6.5)

#6. plot NearestExonJB_distance 
p_nearest <- vbox_with_p(all_clean, "nearest_exon_dist_abs",
                         ylab = "Distance to nearest exon (bp, log10)",
                         log10_y = TRUE, title = "Nearest exon distance")
print(p_nearest)
ggsave("feature_plots/nearest_exon_distance_violin_pvals.pdf", p_nearest, width = 7.2, height = 7.8)

#7. plot ptc to exonjunciton distance
##a. ptc-variant distance


##b. variant-cds end distance
p_nearest[["data"]][["value"]]

