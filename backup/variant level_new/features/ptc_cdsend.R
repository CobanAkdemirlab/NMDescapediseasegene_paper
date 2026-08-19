#This R script is to plot ptc-cds end distance in a matched way

ggplot(variants_all, 
       aes(x = source, y = dist_to_cds_end, fill = source)) +
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
    c("snv", "snv_control"),
    c("fs", "fs_control")
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
