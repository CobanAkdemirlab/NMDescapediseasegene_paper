############################
## 0. Packages
############################
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(broom)
library(broom.mixed)
library(car)

############################
## 1. Prepare data
############################
variants_all2 <- variants_all %>%
  mutate(
    source = factor(source, levels = c("fs", "fs_control", "snv", "snv_control")),
    NMDesc_region_length = NMD_region_end - NMD_region_start + 1
  ) %>%
  filter(
    !is.na(source),
    !is.na(dist_to_cds_end),
    dist_to_cds_end >= 0,
    !is.na(cds_end),
    cds_end > 0,
    !is.na(cds_ptc_loc),
    cds_ptc_loc > 0,
    !is.na(NMD_region_start),
    !is.na(NMD_region_end),
    !is.na(NMDesc_region_length),
    NMDesc_region_length > 0
  ) %>%
  mutate(
    log_dist_to_cds_end = log10(dist_to_cds_end + 1),
    log_cds_end = log10(cds_end),
    log_NMDesc_region_length = log10(NMDesc_region_length),
    dist_to_cds_end_rel = dist_to_cds_end / cds_end,
    ptc_rel_nmdesc = (cds_ptc_loc - NMD_region_start) / NMDesc_region_length
  ) %>%
  filter(
    is.finite(log_dist_to_cds_end),
    is.finite(log_cds_end),
    is.finite(log_NMDesc_region_length),
    is.finite(dist_to_cds_end_rel),
    is.finite(ptc_rel_nmdesc)
  )

############################
## 2. Define shared transcripts for matched analysis
############################
fs_tr <- variants_all2 %>%
  filter(source == "fs") %>%
  pull(transcript) %>%
  unique()

fs_ctrl_tr <- variants_all2 %>%
  filter(source == "fs_control") %>%
  pull(transcript) %>%
  unique()

snv_tr <- variants_all2 %>%
  filter(source == "snv") %>%
  pull(transcript) %>%
  unique()

snv_ctrl_tr <- variants_all2 %>%
  filter(source == "snv_control") %>%
  pull(transcript) %>%
  unique()

shared_fs <- intersect(fs_tr, fs_ctrl_tr)
shared_snv <- intersect(snv_tr, snv_ctrl_tr)

## matched datasets for separate analyses
variants_fs_matched <- variants_all2 %>%
  filter(
    transcript %in% shared_fs,
    source %in% c("fs", "fs_control")
  ) %>%
  droplevels()

variants_snv_matched <- variants_all2 %>%
  filter(
    transcript %in% shared_snv,
    source %in% c("snv", "snv_control")
  ) %>%
  droplevels()

## matched dataset with all four groups shown together in one plot
variants_matched4 <- bind_rows(
  variants_all2 %>%
    filter(
      transcript %in% shared_fs,
      source %in% c("fs", "fs_control")
    ),
  variants_all2 %>%
    filter(
      transcript %in% shared_snv,
      source %in% c("snv", "snv_control")
    )
) %>%
  mutate(source = factor(source, levels = c("fs", "fs_control", "snv", "snv_control")))

############################
## 3. Colors and comparisons
############################
group_colors <- c(
  "fs" = "#4C78A8",
  "fs_control" = "#9FB7D3",
  "snv" = "#59A14F",
  "snv_control" = "#8CD17D"
)

pair_comparisons <- list(
  c("fs", "fs_control"),
  c("snv", "snv_control")
)

############################
## 4. Plot functions
############################
plot_feature_4groups <- function(dat, yvar, ylab, title_text, log_scale = TRUE) {
  p <- ggplot(dat, aes(x = source, y = .data[[yvar]], fill = source)) +
    geom_violin(trim = TRUE, alpha = 0.85, color = "black") +
    geom_boxplot(width = 0.15, outlier.size = 0.8, color = "black") +
    scale_fill_manual(values = group_colors, guide = "none") +
    labs(
      title = title_text,
      x = NULL,
      y = ylab
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.title.y = element_text(face = "bold")
    ) +
    stat_compare_means(
      comparisons = pair_comparisons,
      method = "wilcox.test",
      label = "p.format",
      step.increase = 0.12
    )
  
  if (log_scale) {
    p <- p + scale_y_log10()
  }
  
  return(p)
}

plot_cds_vs_nmdesc <- function(dat, title_text) {
  ggplot(dat, aes(x = cds_end, y = NMDesc_region_length, color = source)) +
    geom_point(alpha = 0.25, size = 1.2) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = group_colors) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    labs(
      title = title_text,
      x = "CDS length (bp, log10 scale)",
      y = "NMDesc region length (bp, log10 scale)",
      color = "Group"
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20)
    )
}

############################
## 5. Plots: unfiltered/full dataset (all four groups)
############################
p_dist_full <- plot_feature_4groups(
  dat = variants_all2,
  yvar = "dist_to_cds_end",
  ylab = "Distance to CDS end (bp, log10 scale)",
  title_text = "Distance from PTC to CDS end across variant groups",
  log_scale = TRUE
)

p_nmdesc_full <- plot_feature_4groups(
  dat = variants_all2,
  yvar = "NMDesc_region_length",
  ylab = "NMDesc region length (bp, log10 scale)",
  title_text = "NMDesc region length across variant groups",
  log_scale = TRUE
)

p_cds_full <- plot_feature_4groups(
  dat = variants_all2,
  yvar = "cds_end",
  ylab = "CDS length (bp, log10 scale)",
  title_text = "CDS length across variant groups",
  log_scale = TRUE
)

p_relation_full <- plot_cds_vs_nmdesc(
  dat = variants_all2,
  title_text = "CDS length vs NMDesc region length (full dataset)"
)

print(p_dist_full)
print(p_nmdesc_full)
print(p_cds_full)
print(p_relation_full)

############################
## 6. Plots: filtered/matched dataset (still all four groups)
############################
p_dist_matched4 <- plot_feature_4groups(
  dat = variants_matched4,
  yvar = "dist_to_cds_end",
  ylab = "Distance to CDS end (bp, log10 scale)",
  title_text = "Matched analysis: distance from PTC to CDS end",
  log_scale = TRUE
)

p_nmdesc_matched4 <- plot_feature_4groups(
  dat = variants_matched4,
  yvar = "NMDesc_region_length",
  ylab = "NMDesc region length (bp, log10 scale)",
  title_text = "Matched analysis: NMDesc region length",
  log_scale = TRUE
)

p_cds_matched4 <- plot_feature_4groups(
  dat = variants_matched4,
  yvar = "cds_end",
  ylab = "CDS length (bp, log10 scale)",
  title_text = "Matched analysis: CDS length",
  log_scale = TRUE
)

p_relation_matched4 <- plot_cds_vs_nmdesc(
  dat = variants_matched4,
  title_text = "CDS length vs NMDesc region length (matched dataset)"
)

print(p_dist_matched4)
print(p_nmdesc_matched4)
print(p_cds_matched4)
print(p_relation_matched4)

############################
## 7. Correlation and linear relation: full + matched
############################
cor_full <- cor.test(
  variants_all2$cds_end,
  variants_all2$NMDesc_region_length,
  method = "spearman"
)
print(cor_full)

fit_cds_nmdesc_full <- lm(
  log_NMDesc_region_length ~ log_cds_end,
  data = variants_all2
)
print(summary(fit_cds_nmdesc_full))

cor_matched <- cor.test(
  variants_matched4$cds_end,
  variants_matched4$NMDesc_region_length,
  method = "spearman"
)
print(cor_matched)

fit_cds_nmdesc_matched <- lm(
  log_NMDesc_region_length ~ log_cds_end,
  data = variants_matched4
)
print(summary(fit_cds_nmdesc_matched))

## optional VIF check
vif_model_full <- lm(
  log_dist_to_cds_end ~ log_cds_end + log_NMDesc_region_length,
  data = variants_all2
)
print(vif(vif_model_full))

############################
## 8. Modeling function
## feature model:
## feature ~ source + log_dist_to_cds_end + log_cds_end + log_NMDesc_region_length
## mixed model:
## feature ~ source + log_dist_to_cds_end + log_cds_end + log_NMDesc_region_length + (1 | transcript)
############################
run_feature_models <- function(dat, feature, feature_type = c("continuous", "binary")) {
  feature_type <- match.arg(feature_type)
  
  stopifnot(feature %in% colnames(dat))
  
  dat2 <- dat %>%
    filter(
      !is.na(.data[[feature]]),
      !is.na(source),
      !is.na(log_dist_to_cds_end),
      !is.na(log_cds_end),
      !is.na(log_NMDesc_region_length),
      !is.na(transcript)
    )
  
  if (feature_type == "binary") {
    dat2 <- dat2 %>%
      filter(.data[[feature]] %in% c(0, 1))
  }
  
  fixed_formula <- as.formula(
    paste0("`", feature, "` ~ source + log_dist_to_cds_end + log_cds_end + log_NMDesc_region_length")
  )
  
  mixed_formula <- as.formula(
    paste0("`", feature, "` ~ source + log_dist_to_cds_end + log_cds_end + log_NMDesc_region_length + (1 | transcript)")
  )
  
  if (feature_type == "continuous") {
    fit_fixed <- lm(fixed_formula, data = dat2)
    fit_mixed <- lmer(mixed_formula, data = dat2, REML = TRUE)
  } else {
    fit_fixed <- glm(fixed_formula, data = dat2, family = binomial())
    fit_mixed <- glmer(mixed_formula, data = dat2, family = binomial())
  }
  
  list(
    data = dat2,
    fit_fixed = fit_fixed,
    fit_mixed = fit_mixed,
    tidy_fixed = tryCatch(broom::tidy(fit_fixed), error = function(e) NULL),
    tidy_mixed = tryCatch(broom.mixed::tidy(fit_mixed, effects = "fixed"), error = function(e) NULL)
  )
}

############################
## 9. Dataset list for modeling
## full/unfiltered and filtered/matched, separately for fs and snv
############################
dat_full_fs <- variants_all2 %>%
  filter(source %in% c("fs", "fs_control")) %>%
  droplevels()

dat_full_snv <- variants_all2 %>%
  filter(source %in% c("snv", "snv_control")) %>%
  droplevels()

dat_matched_fs <- variants_fs_matched %>%
  droplevels()

dat_matched_snv <- variants_snv_matched %>%
  droplevels()

############################
## 10. Put your feature names here
############################
continuous_features <- c(
  # "your_continuous_feature"
)

binary_features <- c(
  # "ppi_overlap",
  # "pfam_overlap",
  # "idr_overlap"
)

############################
## 11. Run all continuous-feature models
############################
continuous_results <- list()

for (feat in continuous_features) {
  continuous_results[[paste0(feat, "_full_fs")]] <- run_feature_models(
    dat = dat_full_fs,
    feature = feat,
    feature_type = "continuous"
  )
  
  continuous_results[[paste0(feat, "_matched_fs")]] <- run_feature_models(
    dat = dat_matched_fs,
    feature = feat,
    feature_type = "continuous"
  )
  
  continuous_results[[paste0(feat, "_full_snv")]] <- run_feature_models(
    dat = dat_full_snv,
    feature = feat,
    feature_type = "continuous"
  )
  
  continuous_results[[paste0(feat, "_matched_snv")]] <- run_feature_models(
    dat = dat_matched_snv,
    feature = feat,
    feature_type = "continuous"
  )
}

############################
## 12. Run all binary-feature models
############################
binary_results <- list()

for (feat in binary_features) {
  binary_results[[paste0(feat, "_full_fs")]] <- run_feature_models(
    dat = dat_full_fs,
    feature = feat,
    feature_type = "binary"
  )
  
  binary_results[[paste0(feat, "_matched_fs")]] <- run_feature_models(
    dat = dat_matched_fs,
    feature = feat,
    feature_type = "binary"
  )
  
  binary_results[[paste0(feat, "_full_snv")]] <- run_feature_models(
    dat = dat_full_snv,
    feature = feat,
    feature_type = "binary"
  )
  
  binary_results[[paste0(feat, "_matched_snv")]] <- run_feature_models(
    dat = dat_matched_snv,
    feature = feat,
    feature_type = "binary"
  )
}

############################
## 13. Print summaries
############################
for (nm in names(continuous_results)) {
  cat("\n==============================\n")
  cat("CONTINUOUS MODEL:", nm, "\n")
  cat("==============================\n")
  
  cat("\n--- Fixed model summary ---\n")
  print(summary(continuous_results[[nm]]$fit_fixed))
  
  cat("\n--- Mixed model summary ---\n")
  print(summary(continuous_results[[nm]]$fit_mixed))
  
  cat("\n--- Fixed model tidy ---\n")
  print(continuous_results[[nm]]$tidy_fixed)
  
  cat("\n--- Mixed model tidy ---\n")
  print(continuous_results[[nm]]$tidy_mixed)
}

for (nm in names(binary_results)) {
  cat("\n==============================\n")
  cat("BINARY MODEL:", nm, "\n")
  cat("==============================\n")
  
  cat("\n--- Fixed model summary ---\n")
  print(summary(binary_results[[nm]]$fit_fixed))
  
  cat("\n--- Mixed model summary ---\n")
  print(summary(binary_results[[nm]]$fit_mixed))
  
  cat("\n--- Fixed model tidy ---\n")
  print(binary_results[[nm]]$tidy_fixed)
  
  cat("\n--- Mixed model tidy ---\n")
  print(binary_results[[nm]]$tidy_mixed)
}

############################
## 14. Optional: export coefficient tables
############################
for (nm in names(continuous_results)) {
  if (!is.null(continuous_results[[nm]]$tidy_fixed)) {
    write.csv(
      continuous_results[[nm]]$tidy_fixed,
      paste0(nm, "_fixed_model.csv"),
      row.names = FALSE
    )
  }
  if (!is.null(continuous_results[[nm]]$tidy_mixed)) {
    write.csv(
      continuous_results[[nm]]$tidy_mixed,
      paste0(nm, "_mixed_model.csv"),
      row.names = FALSE
    )
  }
}

for (nm in names(binary_results)) {
  if (!is.null(binary_results[[nm]]$tidy_fixed)) {
    write.csv(
      binary_results[[nm]]$tidy_fixed,
      paste0(nm, "_fixed_model.csv"),
      row.names = FALSE
    )
  }
  if (!is.null(binary_results[[nm]]$tidy_mixed)) {
    write.csv(
      binary_results[[nm]]$tidy_mixed,
      paste0(nm, "_mixed_model.csv"),
      row.names = FALSE
    )
  }
}

############################
## 15. Optional: save plots
############################
ggsave("dist_to_cds_end_full_4groups.png", p_dist_full, width = 9, height = 6, dpi = 300)
ggsave("nmdesc_region_length_full_4groups.png", p_nmdesc_full, width = 9, height = 6, dpi = 300)
ggsave("cds_length_full_4groups.png", p_cds_full, width = 9, height = 6, dpi = 300)
ggsave("cds_vs_nmdesc_full.png", p_relation_full, width = 8, height = 6, dpi = 300)

ggsave("dist_to_cds_end_matched_4groups.png", p_dist_matched4, width = 9, height = 6, dpi = 300)
ggsave("nmdesc_region_length_matched_4groups.png", p_nmdesc_matched4, width = 9, height = 6, dpi = 300)
ggsave("cds_length_matched_4groups.png", p_cds_matched4, width = 9, height = 6, dpi = 300)
ggsave("cds_vs_nmdesc_matched.png", p_relation_matched4, width = 8, height = 6, dpi = 300)