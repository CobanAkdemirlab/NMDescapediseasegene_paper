library(rstatix)

#compare cds length between disease genes and control genes

#cut gene_all in five vategory by cds length
gene_all4$cds_length 
gene_all <- gene_all %>%
  mutate(
    cds_length = nchar(coding),
    cds_group = cut(
      cds_length,
      breaks = c(0, 1000, 3000, 5000, 10000, Inf),
      labels = c("1-1000", "1000-3000", "3000-5000", "5000-10000", ">10000"),
      right = TRUE,
      include.lowest = TRUE
    )
  )

#only keep group "100-300", "300-500", "500-1000",
df_bin = gene_all %>%
  group_by(group,cds_group) %>%
  summarise(
    mean_cds = mean(cds_length, na.rm = TRUE),
    median_cds = median(cds_length, na.rm = TRUE),
    n = sum(!is.na(NumberSubmitters)),
    min_ns= min(NumberSubmitters, na.rm = TRUE),
    median_ns = median(NumberSubmitters, na.rm = TRUE),
    mean = mean(NumberSubmitters, na.rm = TRUE),
    sd = sd(NumberSubmitters, na.rm = TRUE)
  )

#significant in 2 and 3
df1 <- gene_all %>%
  filter(cds_group %in% c("1-1000"))
df2 <- gene_all %>%
  filter(cds_group %in% c("1000-3000"))
df3 <- gene_all %>%
  filter(cds_group %in% c("3000-5000"))
df4 = gene_all %>%
  filter(cds_group %in% c("5000-10000"))
df5 = gene_all %>%
  filter(cds_group %in% c(">10000"))

df1 %>%
  wilcox_test(
    NumberSubmitters ~ group,
    comparisons = list(
      c("snv", "snv_control"),
      c("fs", "fs_control")
    ),
    p.adjust.method = "BH"
  )
df2 %>%
  wilcox_test(
    NumberSubmitters ~ group,
    comparisons = list(
      c("snv", "snv_control"),
      c("fs", "fs_control")
    ),
    p.adjust.method = "BH"
  )
df3 %>%
  wilcox_test(
    NumberSubmitters ~ group,
    comparisons = list(
      c("snv", "snv_control"),
      c("fs", "fs_control")
    ),
    p.adjust.method = "BH"
  )
df4 %>%
  wilcox_test(
    NumberSubmitters ~ group,
    comparisons = list(
      c("snv", "snv_control"),
      c("fs", "fs_control")
    ),
    p.adjust.method = "BH"
  )
df5 %>%
  wilcox_test(
    NumberSubmitters ~ group,
    comparisons = list(
      c("snv", "snv_control"),
      c("fs", "fs_control")
    ),
    p.adjust.method = "BH"
  )

library(ggplot2)
library(ggpubr)

ggplot(df, aes(x = group, y = NumberSubmitters, fill = group)) +
  geom_boxplot(width = 0.6, outlier.size = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  
  stat_compare_means(
    comparisons = list(
      c("snv", "snv_control"),
      c("fs", "fs_control")
    ),
    method = "wilcox.test",
    label = "p.format"
  ) +
  
  theme_classic() +
  labs(
    title = "Number of Submitters Across Groups",
    x = NULL,
    y = "Number of Submitters"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )


#plot cds_length by cds_group and group
ggplot(gene_all, aes(x = cds_group, y = cds_length, fill = group)) +
  geom_boxplot(width = 0.6, outlier.size = 0.8) +
  
  stat_compare_means(
    comparisons = list(
      c("snv", "snv_control"),
      c("fs", "fs_control")
    ),
    method = "wilcox.test",
    label = "p.format"
  ) +
  
  theme_classic() +
  labs(
    title = "CDS Length Across Groups and CDS Length Categories",
    x = "CDS Length Category",
    y = "CDS Length (bp)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )