#This R script is to plot ptc-cds end distance in a matched way
variants_all5_1_0621$group

group_colors2 <- c(
  fs_disease  = "#1f77b4",
  fs_control  = "#aec7e8",
  snv_disease = "#2ca02c",
  snv_control = "#98df8a"
)
ggplot(variants_all5_1_0621, 
       aes(x = group, y = dist_to_cds_end, fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  scale_y_log10() + 
  geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
  scale_fill_manual(values = group_colors2, guide = "none") +
  labs(
    title = "Distance from PTC to CDS end across variant groups",
    x = NULL,
    y = "Distance to CDS end (bp)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  ) + 
  stat_compare_means(
  comparisons = list(
    c("snv_disease", "snv_control"),
    c("fs_disease", "fs_control")
  ),
  method = "wilcox.test",
  label = "p.format"
)

# Make sure groups match your color keys
pfam_all_pervariant_cds$group <- factor(
  pfam_all_pervariant_cds$group,
  levels = c("fs","fs_control",
             "Plus1","Plus1_control",
             "snv","snv_control")
)

diff_ptc = variants_all %>% 
  group_by(transcript) %>% 
  summarise(
    control_ptc_cds = mean(dist_to_cds_end[source %in% c("fs_control", "snv_control")], na.rm = TRUE),
    di_ptc_cds      = mean(dist_to_cds_end[source %in% c("fs", "snv")], na.rm = TRUE),
    diff_ptc_cds    = di_ptc_cds - control_ptc_cds,
    .groups = "drop"
  )
summary(diff_ptc$diff_ptc_cds)

hist(diff_ptc$diff_ptc_cds, breaks = 50, main = "Distribution of PTC-CDS end differences", xlab = "Difference in PTC-CDS end (di - control)")


library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(grid)

# =========================
# 1. Colors
# =========================
group_colors2 <- c(
  fs_control  = "#aec7e8",
  fs_disease  = "#1f77b4",
  snv_control = "#98df8a",
  snv_disease = "#2ca02c"
)

# =========================
# 2. Prepare data
# =========================
plot_dat <- variants_all5_1_0621 %>%
  filter(group %in% c(
    "fs_control", "fs_disease",
    "snv_control", "snv_disease"
  )) %>%
  filter(
    !is.na(uniprot),
    !is.na(dist_to_cds_end),
    dist_to_cds_end > 0
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "fs_control", "fs_disease",
        "snv_control", "snv_disease"
      )
    ),
    variant_type = factor(
      ifelse(grepl("^fs", group), "FS", "SNV"),
      levels = c("FS", "SNV")
    ),
    disease_status = factor(
      ifelse(grepl("disease", group), "disease", "control"),
      levels = c("control", "disease")
    ),
    log_dist = log10(dist_to_cds_end)
  )

# =========================
# 3. Mixed-effect models
# =========================
fit_fs <- lmer(
  log_dist ~ disease_status + (1 | uniprot),
  data = plot_dat %>% filter(variant_type == "FS")
)

fit_snv <- lmer(
  log_dist ~ disease_status + (1 | uniprot),
  data = plot_dat %>% filter(variant_type == "SNV")
)

p_fs <- summary(fit_fs)$coefficients[
  "disease_statusdisease",
  "Pr(>|t|)"
]

p_snv <- summary(fit_snv)$coefficients[
  "disease_statusdisease",
  "Pr(>|t|)"
]

p_fs_label <- format.pval(p_fs, digits = 2, eps = 1e-300)
p_snv_label <- format.pval(p_snv, digits = 2, eps = 1e-300)

# =========================
# 4. Plot limit
# =========================
y_min <- min(plot_dat$dist_to_cds_end, na.rm = TRUE)
y_max <- max(plot_dat$dist_to_cds_end, na.rm = TRUE)

bracket_y <- y_max * 1.05
label_y   <- y_max * 1.20
plot_top  <- y_max * 1.35

# =========================
# 5. Plot
# =========================
ggplot(
  plot_dat,
  aes(x = group, y = dist_to_cds_end, fill = group)
) +
  geom_violin(
    trim = TRUE,
    alpha = 0.7,
    color = "gray30"
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.size = 0.6,
    color = "black",
    alpha = 0.9
  ) +
  annotate(
    "segment",
    x = 1, xend = 2,
    y = bracket_y, yend = bracket_y,
    linewidth = 0.6
  ) +
  annotate(
    "segment",
    x = 1, xend = 1,
    y = bracket_y, yend = bracket_y / 1.08,
    linewidth = 0.6
  ) +
  annotate(
    "segment",
    x = 2, xend = 2,
    y = bracket_y, yend = bracket_y / 1.08,
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 1.5,
    y = label_y,
    label = p_fs_label,
    size = 5
  ) +
  annotate(
    "segment",
    x = 3, xend = 4,
    y = bracket_y, yend = bracket_y,
    linewidth = 0.6
  ) +
  annotate(
    "segment",
    x = 3, xend = 3,
    y = bracket_y, yend = bracket_y / 1.08,
    linewidth = 0.6
  ) +
  annotate(
    "segment",
    x = 4, xend = 4,
    y = bracket_y, yend = bracket_y / 1.08,
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 3.5,
    y = label_y,
    label = p_snv_label,
    size = 5
  ) +
  scale_y_log10(
    limits = c(y_min, plot_top),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  scale_fill_manual(values = group_colors2, guide = "none") +
  labs(
    title = "Distance from PTC to CDS end across variant groups",
    subtitle = "P-values from random-intercept mixed effect models",
    x = NULL,
    y = "Distance to CDS end (bp)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )



library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

set.seed(42)

group_colors2 <- c(
  fs_control  = "#aec7e8",
  fs_disease  = "#1f77b4",
  snv_control = "#98df8a",
  snv_disease = "#2ca02c"
)

boot_dat <- variants_all5_1_0621 %>%
  filter(group %in% c(
    "fs_control", "fs_disease",
    "snv_control", "snv_disease"
  )) %>%
  filter(
    !is.na(uniprot),
    !is.na(dist_to_cds_end),
    dist_to_cds_end > 0
  ) %>%
  mutate(
    variant_type = case_when(
      group %in% c("fs_control", "fs_disease") ~ "FS",
      group %in% c("snv_control", "snv_disease") ~ "SNV"
    ),
    disease_status = case_when(
      group %in% c("fs_disease", "snv_disease") ~ "disease",
      group %in% c("fs_control", "snv_control") ~ "control"
    ),
    group = factor(
      group,
      levels = c(
        "fs_control", "fs_disease",
        "snv_control", "snv_disease"
      )
    ),
    log_dist = log10(dist_to_cds_end)
  )

one_boot_paired <- function(data, variant_type_keep) {
  
  sampled <- data %>%
    filter(variant_type == variant_type_keep) %>%
    group_by(uniprot, group) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  paired_wide <- sampled %>%
    select(uniprot, disease_status, dist_to_cds_end, log_dist) %>%
    pivot_wider(
      names_from = disease_status,
      values_from = c(dist_to_cds_end, log_dist)
    ) %>%
    filter(
      !is.na(dist_to_cds_end_control),
      !is.na(dist_to_cds_end_disease)
    )
  
  if (nrow(paired_wide) < 2) {
    return(tibble(
      variant_type = variant_type_keep,
      n_pair = nrow(paired_wide),
      median_diff_log10 = NA_real_,
      p_value = NA_real_
    ))
  }
  
  test <- wilcox.test(
    paired_wide$log_dist_disease,
    paired_wide$log_dist_control,
    paired = TRUE,
    exact = FALSE
  )
  
  tibble(
    variant_type = variant_type_keep,
    n_pair = nrow(paired_wide),
    median_diff_log10 = median(
      paired_wide$log_dist_disease -
        paired_wide$log_dist_control,
      na.rm = TRUE
    ),
    p_value = test$p.value
  )
}

n_boot <- 1000

boot_results <- map_dfr(
  1:n_boot,
  function(i) {
    bind_rows(
      one_boot_paired(boot_dat, "FS"),
      one_boot_paired(boot_dat, "SNV")
    ) %>%
      mutate(iteration = i)
  }
)

boot_summary <- boot_results %>%
  group_by(variant_type) %>%
  summarise(
    n_boot = n(),
    median_n_pair = median(n_pair, na.rm = TRUE),
    median_diff_log10 = median(median_diff_log10, na.rm = TRUE),
    diff_log10_025 = quantile(median_diff_log10, 0.025, na.rm = TRUE),
    diff_log10_975 = quantile(median_diff_log10, 0.975, na.rm = TRUE),
    median_p = median(p_value, na.rm = TRUE),
    p_025 = quantile(p_value, 0.025, na.rm = TRUE),
    p_975 = quantile(p_value, 0.975, na.rm = TRUE),
    prop_p_less_0.05 = mean(p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

print(boot_summary)

plot_sampled <- boot_dat %>%
  group_by(uniprot, group) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  group_by(uniprot, variant_type) %>%
  filter(n_distinct(disease_status) == 2) %>%
  ungroup()

ggplot(
  plot_sampled,
  aes(x = group, y = dist_to_cds_end, fill = group)
) +
  geom_violin(
    trim = TRUE,
    alpha = 0.7,
    color = "gray30"
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.size = 0.6,
    color = "black",
    alpha = 0.9
  ) +
  scale_y_log10() +
  scale_fill_manual(values = group_colors2, guide = "none") +
  labs(
    title = "Distance from PTC to CDS end across variant groups",
    subtitle = "Bootstrap paired sampling: one variant per gene per group",
    x = NULL,
    y = "Distance to CDS end (bp)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
