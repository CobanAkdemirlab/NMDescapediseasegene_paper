#This R script is to plot NMDesc region length, and NMDesc region length/cds length 


ggplot(data = gene_all, aes(x = group, y = NMDesc_region_length, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  scale_y_log10() +
  
  stat_compare_means(
    comparisons = list(
      c("snv", "snv_control"),
      c("fs", "fs_control")
    ),
    method = "wilcox.test",
    label = "p.format"
  ) +
  
  theme_minimal() +
  labs(
    title = "NMDesc Region Length by Group (log10 scale)",
    x = "Group",
    y = "NMDesc Region Length (log10)"
  )


#plot nmdesc/cds
gene_all$NMDesc_cds_ratio = gene_all$NMDesc_region_length / gene_all$cds_length

ggplot(data = gene_all, aes(x = group, y = NMDesc_cds_ratio, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  
  stat_compare_means(
    comparisons = list(
      c("snv", "snv_control"),
      c("fs", "fs_control")
    ),
    method = "wilcox.test",
    label = "p.format"
  ) +
  
  theme_minimal() +
  labs(
    title = "NMDesc Region Length Ratio by Group",
    x = "Group",
    y = "NMDesc Region Length (log10)"
  )