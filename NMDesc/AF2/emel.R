library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(patchwork)
WT_vars_NMD_region <- read_csv("~/Downloads/WT-vars_NMD_region.csv")
WT_vars_full_length <- read_csv("~/Downloads/WT-vars_full-length.csv")

#1. subgroup into 6 categories, match by id
process_group <- function(data, filter_keys) {
  data %>%
    separate(id, into = c("chr", "pos", "ref", "alt"), sep = "_", remove = FALSE) %>%
    mutate(pos = as.character(as.integer(pos))) %>%
    mutate(id2 = str_c(chr, ":", pos, "|", ref, "|", alt)) %>%
    filter(id2 %in% filter_keys)
}

# Apply to all 6 groups
WT_NMD_snv             <- process_group(WT_vars_NMD_region, snv_variants$key)
WT_NMD_snv_control     <- process_group(WT_vars_NMD_region, snv_control_variants$key)
WT_NMD_plus1           <- process_group(WT_vars_NMD_region, plus1_variants$key)
WT_NMD_plus1_control   <- process_group(WT_vars_NMD_region, plus1_control_variants$key)
WT_NMD_minus1          <- process_group(WT_vars_NMD_region, minus1_variants$key)
WT_NMD_minus1_control  <- process_group(WT_vars_NMD_region, minus1_control_variants$key)
WT_full_snv             <- process_group(WT_vars_full_length, snv_variants$key)
WT_full_snv_control     <- process_group(WT_vars_full_length, snv_control_variants$key)
WT_full_plus1           <- process_group(WT_vars_full_length, plus1_variants$key)
WT_full_plus1_control   <- process_group(WT_vars_full_length, plus1_control_variants$key)
WT_full_minus1          <- process_group(WT_vars_full_length, minus1_variants$key)
WT_full_minus1_control  <- process_group(WT_vars_full_length, minus1_control_variants$key)

#2. analyze aa content change
##2.1 Define amino acid list
aa_list <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
             "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

##2.2 Define a function to extract AA % in WT and var
extract_aa_comp <- function(df, group_label, suffix = "_WT_FL") {
  # Identify all amino acid composition columns based on suffix
  aa_cols <- grep(paste0("^aa_.*", suffix, "$"), colnames(df), value = TRUE)
  
  if (length(aa_cols) == 0) {
    stop(paste0("No AA composition columns with suffix ", suffix, " found."))
  }
  
  df %>%
    select(any_of(aa_cols)) %>%
    mutate(group = group_label, row = row_number()) %>%
    pivot_longer(cols = -c(group, row), names_to = "aa", values_to = "percent") %>%
    mutate(aa = gsub("aa_|_WT_nmd|_vars_nmd|_WT_FL|_var_FL", "", aa))
}



##2.3 Create combined aa dataframe for all 6 groups
# WT NMD
aa_WT_NMD <- bind_rows(
  extract_aa_comp(WT_NMD_snv, "SNV", "_WT_nmd"),
  extract_aa_comp(WT_NMD_snv_control, "SNV_Control", "_WT_nmd"),
  extract_aa_comp(WT_NMD_plus1, "Plus1", "_WT_nmd"),
  extract_aa_comp(WT_NMD_plus1_control, "Plus1_Control", "_WT_nmd"),
  extract_aa_comp(WT_NMD_minus1, "Minus1", "_WT_nmd"),
  extract_aa_comp(WT_NMD_minus1_control, "Minus1_Control", "_WT_nmd")
)

# Variant NMD
aa_var_NMD <- bind_rows(
  extract_aa_comp(WT_NMD_snv, "SNV", "_vars_nmd"),
  extract_aa_comp(WT_NMD_snv_control, "SNV_Control", "_vars_nmd"),
  extract_aa_comp(WT_NMD_plus1, "Plus1", "_vars_nmd"),
  extract_aa_comp(WT_NMD_plus1_control, "Plus1_Control", "_vars_nmd"),
  extract_aa_comp(WT_NMD_minus1, "Minus1", "_vars_nmd"),
  extract_aa_comp(WT_NMD_minus1_control, "Minus1_Control", "_vars_nmd")
)

# WT Full-Length
aa_WT_FL <- bind_rows(
  extract_aa_comp(WT_full_snv, "SNV", "_WT_FL"),
  extract_aa_comp(WT_full_snv_control, "SNV_Control", "_WT_FL"),
  extract_aa_comp(WT_full_plus1, "Plus1", "_WT_FL"),
  extract_aa_comp(WT_full_plus1_control, "Plus1_Control", "_WT_FL"),
  extract_aa_comp(WT_full_minus1, "Minus1", "_WT_FL"),
  extract_aa_comp(WT_full_minus1_control, "Minus1_Control", "_WT_FL")
)

# Variant Full-Length
aa_var_FL <- bind_rows(
  extract_aa_comp(WT_full_snv, "SNV", "_var_FL"),
  extract_aa_comp(WT_full_snv_control, "SNV_Control", "_var_FL"),
  extract_aa_comp(WT_full_plus1, "Plus1", "_var_FL"),
  extract_aa_comp(WT_full_plus1_control, "Plus1_Control", "_var_FL"),
  extract_aa_comp(WT_full_minus1, "Minus1", "_var_FL"),
  extract_aa_comp(WT_full_minus1_control, "Minus1_Control", "_var_FL")
)

##2.4 Assign group pairs for comparison
# Add 'region' column to each dataset
aa_WT_NMD$region <- "WT_NMD"
aa_var_NMD$region <- "Var_NMD"
aa_WT_FL$region <- "WT_FL"
aa_var_FL$region <- "Var_FL"

# Combine
aa_comp_all <- bind_rows(aa_WT_NMD, aa_var_NMD, aa_WT_FL, aa_var_FL)

aa_comp_all <- aa_comp_all %>%
  mutate(group_pair = case_when(
    group %in% c("SNV", "SNV_Control") ~ "SNV_vs_Control",
    group %in% c("Plus1", "Plus1_Control") ~ "Plus1_vs_Control",
    group %in% c("Minus1", "Minus1_Control") ~ "Minus1_vs_Control"
  ))

##2.5 Perform Wilcoxon tests per amino acid per group_pair
stat_results <- aa_comp_all %>%
  group_by(region, group_pair, aa) %>%
  summarise(
    p = wilcox.test(percent ~ group, data = cur_data(), exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    sig = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01 ~ "**",
      p_adj <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

##2.6 Plot for each group_pairlibrary(ggplot2)
plot_group_region_pair <- function(df, region_name, group_name, stat_df) {
  df_filtered <- df %>% filter(region == region_name, group_pair == group_name)
  stat_filtered <- stat_df %>% filter(region == region_name, group_pair == group_name)
  
  ggplot(df_filtered, aes(x = aa, y = percent, fill = group)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8)) +
    theme_minimal() +
    labs(title = paste(region_name, "-", group_name),
         x = "Amino Acid", y = "% Composition") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c(
      "SNV" = "#2ca02c", "SNV_Control" = "#98df8a",
      "Plus1" = "#ff7f0e", "Plus1_Control" = "#ffbb78",
      "Minus1" = "#1f77b4", "Minus1_Control" = "#aec7e8"
    )) +
    stat_compare_means(aes(group = group), method = "wilcox.test", label = "p.signif",
                       label.y = max(df_filtered$percent, na.rm = TRUE) + 2,
                       comparisons = list(unique(df_filtered$group)))
}

##2.7 Draw plots for all comparisons
# Regions and group pairs to iterate
regions <- unique(aa_comp_all$region)
group_pairs <- unique(aa_comp_all$group_pair)

# Create all plots
plot_list <- list()
for (r in regions) {
  for (gp in group_pairs) {
    plot_name <- paste("plot", r, gp, sep = "_")
    plot_list[[plot_name]] <- plot_group_region_pair(aa_comp_all, r, gp, stat_results)
  }
}

#print all 12 plots in one pdf
pdf("~/Downloads/aa_comp_all_12plots.pdf", width = 10, height = 6)  # adjust size as needed
for (plot_name in names(plot_list)) {
  print(plot_list[[plot_name]])
}
dev.off()

print(plot_list)
ggsave("plot_WT_NMD_SNV_vs_Control.pdf", plot = plot_list[["plot_WT_NMD_SNV_vs_Control"]])
write.csv(aa_comp_all, "~/Downloads/aa_comp_all_regions.csv", row.names = FALSE)
write.csv(stat_results, "~/Downloads/aa_comp_stats_all_regions.csv", row.names = FALSE)



#3. compare net charge
summary(WT_NMD_snv$net_charge_WT_nmd)
summary(WT_NMD_snv_control$net_charge_WT_nmd)
summary(WT_NMD_plus1$net_charge_WT_nmd)
summary(WT_NMD_plus1_control$net_charge_WT_nmd)
summary(WT_NMD_minus1$net_charge_WT_nmd)
summary(WT_NMD_minus1_control$net_charge_WT_nmd)

##3.1: Combine net charge groups into WT and vars df
net_charge_WT_nmd <- bind_rows(
  WT_NMD_snv             %>% select(net_charge_WT_nmd) %>% mutate(group = "SNV"),
  WT_NMD_snv_control     %>% select(net_charge_WT_nmd) %>% mutate(group = "SNV_Control"),
  WT_NMD_plus1           %>% select(net_charge_WT_nmd) %>% mutate(group = "Plus1"),
  WT_NMD_plus1_control   %>% select(net_charge_WT_nmd) %>% mutate(group = "Plus1_Control"),
  WT_NMD_minus1          %>% select(net_charge_WT_nmd) %>% mutate(group = "Minus1"),
  WT_NMD_minus1_control  %>% select(net_charge_WT_nmd) %>% mutate(group = "Minus1_Control")
)
net_charge_var_nmd = bind_rows(
  WT_NMD_snv             %>% select(net_charge_vars_nmd) %>% mutate(group = "SNV"),
  WT_NMD_snv_control     %>% select(net_charge_vars_nmd) %>% mutate(group = "SNV_Control"),
  WT_NMD_plus1           %>% select(net_charge_vars_nmd) %>% mutate(group = "Plus1"),
  WT_NMD_plus1_control   %>% select(net_charge_vars_nmd) %>% mutate(group = "Plus1_Control"),
  WT_NMD_minus1          %>% select(net_charge_vars_nmd) %>% mutate(group = "Minus1"),
  WT_NMD_minus1_control  %>% select(net_charge_vars_nmd) %>% mutate(group = "Minus1_Control")
)
net_charge_WT_FL <- bind_rows(
  WT_full_snv            %>% select(net_charge_WT_FL) %>% mutate(group = "SNV"),
  WT_full_snv_control    %>% select(net_charge_WT_FL) %>% mutate(group = "SNV_Control"),
  WT_full_plus1          %>% select(net_charge_WT_FL) %>% mutate(group = "Plus1"),
  WT_full_plus1_control  %>% select(net_charge_WT_FL) %>% mutate(group = "Plus1_Control"),
  WT_full_minus1         %>% select(net_charge_WT_FL) %>% mutate(group = "Minus1"),
  WT_full_minus1_control %>% select(net_charge_WT_FL) %>% mutate(group = "Minus1_Control")
)
net_charge_var_FL = bind_rows(
  WT_full_snv            %>% select(net_charge_var_FL) %>% mutate(group = "SNV"),
  WT_full_snv_control    %>% select(net_charge_var_FL) %>% mutate(group = "SNV_Control"),
  WT_full_plus1          %>% select(net_charge_var_FL) %>% mutate(group = "Plus1"),
  WT_full_plus1_control  %>% select(net_charge_var_FL) %>% mutate(group = "Plus1_Control"),
  WT_full_minus1         %>% select(net_charge_var_FL) %>% mutate(group = "Minus1"),
  WT_full_minus1_control %>% select(net_charge_var_FL) %>% mutate(group = "Minus1_Control")
)
write.csv(net_charge_WT_nmd, "~/Downloads/net_charge_WT_nmd.csv", row.names = FALSE)
write.csv(net_charge_var_nmd, "~/Downloads/net_charge_var_nmd.csv", row.names = FALSE)
write.csv(net_charge_WT_FL, "~/Downloads/net_charge_WT_FL.csv", row.names = FALSE)
write.csv(net_charge_var_FL, "~/Downloads/net_charge_var_FL.csv", row.names = FALSE)

my_comparisons <- list(
  c("Minus1", "Minus1_Control"),
  c("Plus1", "Plus1_Control"),
  c("SNV", "SNV_Control")
)
##3.2 Plot WT and var net charge for nmd and full
p1 <- ggplot(net_charge_WT_nmd, aes(x = group, y = net_charge_WT_nmd, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  theme_minimal() +
  labs(title = "Net Charge of WT NMD Regions", x = "Group", y = "Net Charge (WT NMD)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "Minus1" = "#1f77b4", "Minus1_Control" = "#aec7e8",
    "Plus1" = "#ff7f0e", "Plus1_Control" = "#ffbb78",
    "SNV" = "#2ca02c", "SNV_Control" = "#98df8a"
  )) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", tip.length = 0.01)

p2 <- ggplot(net_charge_var_nmd, aes(x = group, y = net_charge_vars_nmd, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  theme_minimal() +
  labs(title = "Net Charge of Variants in NMD", x = "Group", y = "Net Charge (Variant NMD)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "Minus1" = "#1f77b4", "Minus1_Control" = "#aec7e8",
    "Plus1" = "#ff7f0e", "Plus1_Control" = "#ffbb78",
    "SNV" = "#2ca02c", "SNV_Control" = "#98df8a"
  )) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", tip.length = 0.01)

p3 <- ggplot(net_charge_WT_FL, aes(x = group, y = net_charge_WT_FL, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  theme_minimal() +
  labs(title = "Net Charge of WT Full-Length Proteins", x = "Group", y = "Net Charge (WT Full-Length)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "Minus1" = "#1f77b4", "Minus1_Control" = "#aec7e8",
    "Plus1" = "#ff7f0e", "Plus1_Control" = "#ffbb78",
    "SNV" = "#2ca02c", "SNV_Control" = "#98df8a"
  )) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", tip.length = 0.01)

p4 <- ggplot(net_charge_var_FL, aes(x = group, y = net_charge_var_FL, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  theme_minimal() +
  labs(title = "Net Charge of Variants in Full-Length", x = "Group", y = "Net Charge (Variant Full-Length)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "Minus1" = "#1f77b4", "Minus1_Control" = "#aec7e8",
    "Plus1" = "#ff7f0e", "Plus1_Control" = "#ffbb78",
    "SNV" = "#2ca02c", "SNV_Control" = "#98df8a"
  )) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", tip.length = 0.01)

#plot these 4 ggplot into one plot
combined_plot <- (p1 | p2) / (p3 | p4)

# Show the plot
print(combined_plot)

# Optionally save it
ggsave("net_charge_comparison_all.png", combined_plot, width = 14, height = 10, dpi = 300)

#4. compare sasa to see if level of exposure changed
sasa_WT_nmd <- bind_rows(
  WT_NMD_snv             %>% select(rel_sasa_WT_nmd) %>% mutate(group = "SNV"),
  WT_NMD_snv_control     %>% select(rel_sasa_WT_nmd) %>% mutate(group = "SNV_Control"),
  WT_NMD_plus1           %>% select(rel_sasa_WT_nmd) %>% mutate(group = "Plus1"),
  WT_NMD_plus1_control   %>% select(rel_sasa_WT_nmd) %>% mutate(group = "Plus1_Control"),
  WT_NMD_minus1          %>% select(rel_sasa_WT_nmd) %>% mutate(group = "Minus1"),
  WT_NMD_minus1_control  %>% select(rel_sasa_WT_nmd) %>% mutate(group = "Minus1_Control")
)
sasa_var_nmd = bind_rows(
  WT_NMD_snv             %>% select(rel_sasa_var_nmd) %>% mutate(group = "SNV"),
  WT_NMD_snv_control     %>% select(rel_sasa_var_nmd) %>% mutate(group = "SNV_Control"),
  WT_NMD_plus1           %>% select(rel_sasa_var_nmd) %>% mutate(group = "Plus1"),
  WT_NMD_plus1_control   %>% select(rel_sasa_var_nmd) %>% mutate(group = "Plus1_Control"),
  WT_NMD_minus1          %>% select(rel_sasa_var_nmd) %>% mutate(group = "Minus1"),
  WT_NMD_minus1_control  %>% select(rel_sasa_var_nmd) %>% mutate(group = "Minus1_Control")
)
abs_sasa_WT_FL = bind_rows(
  WT_full_snv            %>% select(abs_sasa_WT_FL) %>% mutate(group = "SNV"),
  WT_full_snv_control    %>% select(abs_sasa_WT_FL) %>% mutate(group = "SNV_Control"),
  WT_full_plus1          %>% select(abs_sasa_WT_FL) %>% mutate(group = "Plus1"),
  WT_full_plus1_control  %>% select(abs_sasa_WT_FL) %>% mutate(group = "Plus1_Control"),
  WT_full_minus1         %>% select(abs_sasa_WT_FL) %>% mutate(group = "Minus1"),
  WT_full_minus1_control %>% select(abs_sasa_WT_FL) %>% mutate(group = "Minus1_Control")
)
abs_sasa_mean_WT_FL = bind_rows(
  WT_full_snv            %>% select(abs_sasa_mean_WT_FL) %>% mutate(group = "SNV"),
  WT_full_snv_control    %>% select(abs_sasa_mean_WT_FL) %>% mutate(group = "SNV_Control"),
  WT_full_plus1          %>% select(abs_sasa_mean_WT_FL) %>% mutate(group = "Plus1"),
  WT_full_plus1_control  %>% select(abs_sasa_mean_WT_FL) %>% mutate(group = "Plus1_Control"),
  WT_full_minus1         %>% select(abs_sasa_mean_WT_FL) %>% mutate(group = "Minus1"),
  WT_full_minus1_control %>% select(abs_sasa_mean_WT_FL) %>% mutate(group = "Minus1_Control")
)
abs_sasa_var_FL = bind_rows(
  WT_full_snv            %>% select(abs_sasa_var_FL) %>% mutate(group = "SNV"),
  WT_full_snv_control    %>% select(abs_sasa_var_FL) %>% mutate(group = "SNV_Control"),
  WT_full_plus1          %>% select(abs_sasa_var_FL) %>% mutate(group = "Plus1"),
  WT_full_plus1_control  %>% select(abs_sasa_var_FL) %>% mutate(group = "Plus1_Control"),
  WT_full_minus1         %>% select(abs_sasa_var_FL) %>% mutate(group = "Minus1"),
  WT_full_minus1_control %>% select(abs_sasa_var_FL) %>% mutate(group = "Minus1_Control")
)
abs_sasa_mean_var_FL = bind_rows(
  WT_full_snv            %>% select(abs_sasa_mean_var_FL) %>% mutate(group = "SNV"),
  WT_full_snv_control    %>% select(abs_sasa_mean_var_FL) %>% mutate(group = "SNV_Control"),
  WT_full_plus1          %>% select(abs_sasa_mean_var_FL) %>% mutate(group = "Plus1"),
  WT_full_plus1_control  %>% select(abs_sasa_mean_var_FL) %>% mutate(group = "Plus1_Control"),
  WT_full_minus1         %>% select(abs_sasa_mean_var_FL) %>% mutate(group = "Minus1"),
  WT_full_minus1_control %>% select(abs_sasa_mean_var_FL) %>% mutate(group = "Minus1_Control")
)

ggplot(sasa_WT_nmd, aes(x = group, y = rel_sasa_WT_nmd, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  #geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", tip.length = 0.01) +
  labs(title = "Relative SASA (WT NMD Region)", x = "Group", y = "Relative SASA") +
  scale_fill_manual(values = c(
    "Minus1" = "#1f77b4",
    "Minus1_Control" = "#aec7e8",
    "Plus1" = "#ff7f0e",
    "Plus1_Control" = "#ffbb78",
    "SNV" = "#2ca02c",
    "SNV_Control" = "#98df8a"
  )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(sasa_var_nmd, aes(x = group, y = rel_sasa_var_nmd, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  #geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", tip.length = 0.01) +
  labs(title = "Relative SASA (Variant NMD Region)", x = "Group", y = "Relative SASA") +
  scale_fill_manual(values = c(
    "Minus1" = "#1f77b4",
    "Minus1_Control" = "#aec7e8",
    "Plus1" = "#ff7f0e",
    "Plus1_Control" = "#ffbb78",
    "SNV" = "#2ca02c",
    "SNV_Control" = "#98df8a"
  )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(abs_sasa_WT_FL, aes(x = group, y = abs_sasa_WT_FL, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  #geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", tip.length = 0.01) +
  labs(title = "Absolute SASA (Full-Length WT)", x = "Group", y = "Absolute SASA") +
  scale_fill_manual(values = c(
    "Minus1" = "#1f77b4",
    "Minus1_Control" = "#aec7e8",
    "Plus1" = "#ff7f0e",
    "Plus1_Control" = "#ffbb78",
    "SNV" = "#2ca02c",
    "SNV_Control" = "#98df8a"
  )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(abs_sasa_var_FL, aes(x = group, y = abs_sasa_var_FL, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  #geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", tip.length = 0.01) +
  labs(title = "Absolute SASA (Variant Full-Length)", x = "Group", y = "Absolute SASA") +
  scale_fill_manual(values = c(
    "Minus1" = "#1f77b4",
    "Minus1_Control" = "#aec7e8",
    "Plus1" = "#ff7f0e",
    "Plus1_Control" = "#ffbb78",
    "SNV" = "#2ca02c",
    "SNV_Control" = "#98df8a"
  )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("~/Downloads/sasa_var_FL.png", plot = last_plot(), width = 8, height = 6)

ggplot(abs_sasa_mean_WT_FL, aes(x = group, y = abs_sasa_mean_WT_FL, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  #geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", tip.length = 0.01) +
  labs(title = "Mean Absolute SASA (Full-Length WT)", x = "Group", y = "Mean Absolute SASA") +
  scale_fill_manual(values = c(
    "Minus1" = "#1f77b4",
    "Minus1_Control" = "#aec7e8",
    "Plus1" = "#ff7f0e",
    "Plus1_Control" = "#ffbb78",
    "SNV" = "#2ca02c",
    "SNV_Control" = "#98df8a"
  )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(abs_sasa_mean_var_FL, aes(x = group, y = abs_sasa_mean_var_FL, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  #geom_jitter(width = 0.2, alpha = 0.4, size = 0.8) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", tip.length = 0.01) +
  labs(title = "Mean Absolute SASA (Variant Full-Length)", x = "Group", y = "Mean Absolute SASA") +
  scale_fill_manual(values = c(
    "Minus1" = "#1f77b4",
    "Minus1_Control" = "#aec7e8",
    "Plus1" = "#ff7f0e",
    "Plus1_Control" = "#ffbb78",
    "SNV" = "#2ca02c",
    "SNV_Control" = "#98df8a"
  )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

write.csv(abs_sasa_mean_var_FL, "~/Downloads/abs_sasa_mean_var_FL.csv", row.names = FALSE)
write.csv(abs_sasa_mean_WT_FL, "~/Downloads/abs_sasa_mean_WT_FL.csv", row.names = FALSE)

# Colors (matches your scheme)
group_colors <- c(
  "Minus1" = "#1f77b4", "Minus1_Control" = "#aec7e8",
  "Plus1"  = "#ff7f0e", "Plus1_Control"  = "#ffbb78",
  "SNV"    = "#2ca02c", "SNV_Control"    = "#98df8a"
)

# For stat_compare_means: return the right pair for the given group_pair
get_comparisons_for_pair <- function(group_pair) {
  switch(group_pair,
         "SNV_vs_Control"    = list(c("SNV", "SNV_Control")),
         "Plus1_vs_Control"  = list(c("Plus1", "Plus1_Control")),
         "Minus1_vs_Control" = list(c("Minus1", "Minus1_Control")),
         NULL
  )
}

# Safer Wilcoxon (skips if a group is missing or has all NA)
wilcox_safe <- function(df, y, g) {
  yv <- df[[y]]
  gv <- factor(df[[g]])
  if (length(setdiff(levels(gv), unique(gv))) > 0) return(NA_real_)
  if (nlevels(gv) != 2) return(NA_real_)
  if (sum(!is.na(yv) & !is.na(gv)) < 2) return(NA_real_)
  tryCatch(wilcox.test(yv ~ gv, exact = FALSE)$p.value, error = function(e) NA_real_)
}

# --- (A) Make AA boxplots annotate correctly --------------------------------
# Rebuild plotting function so comparisons map to the correct pair for each group_pair
plot_group_region_pair <- function(df, region_name, group_pair, ypad = 2) {
  df_filtered <- df %>% filter(region == region_name, group_pair == group_pair)
  if (nrow(df_filtered) == 0) return(NULL)
  
  # y max per amino acid for clean label placement
  y_max <- df_filtered %>%
    group_by(aa) %>%
    summarise(y = max(percent, na.rm = TRUE) + ypad, .groups = "drop")
  
  p <- ggplot(df_filtered, aes(x = aa, y = percent, fill = group)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), width = 0.7, color = "black") +
    scale_fill_manual(values = group_colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste(region_name, "-", group_pair), x = "Amino Acid", y = "% Composition")
  
  comps <- get_comparisons_for_pair(group_pair)
  if (!is.null(comps)) {
    # one y-position for each x level
    p <- p + ggpubr::stat_compare_means(
      comparisons = comps,
      method = "wilcox.test",
      label = "p.signif",
      tip.length = 0.01,
      hide.ns = TRUE,
      label.y = y_max$y
    )
  }
  p
}

# Rebuild the list (keeps your 12 plots idea)
regions <- unique(aa_comp_all$region)
group_pairs <- unique(aa_comp_all$group_pair)
plot_list <- list()
for (r in regions) {
  for (gp in group_pairs) {
    nm <- paste("plot", r, gp, sep = "_")
    plot_list[[nm]] <- plot_group_region_pair(aa_comp_all, r, gp)
  }
}

#5. plot mean of plddt
plddt_WT_FL <- bind_rows(
  WT_full_snv            %>% select(plddt_mean_WT_FL) %>% mutate(group = "SNV",            region = "WT_FL",   metric = "pLDDT"),
  WT_full_snv_control    %>% select(plddt_mean_WT_FL) %>% mutate(group = "SNV_Control",    region = "WT_FL",   metric = "pLDDT"),
  WT_full_plus1          %>% select(plddt_mean_WT_FL) %>% mutate(group = "Plus1",          region = "WT_FL",   metric = "pLDDT"),
  WT_full_plus1_control  %>% select(plddt_mean_WT_FL) %>% mutate(group = "Plus1_Control",  region = "WT_FL",   metric = "pLDDT"),
  WT_full_minus1         %>% select(plddt_mean_WT_FL) %>% mutate(group = "Minus1",         region = "WT_FL",   metric = "pLDDT"),
  WT_full_minus1_control %>% select(plddt_mean_WT_FL) %>% mutate(group = "Minus1_Control", region = "WT_FL",   metric = "pLDDT")
) %>% rename(value = plddt_mean_WT_FL)

plddt_var_FL <- bind_rows(
  WT_full_snv            %>% select(plddt_mean_var_FL) %>% mutate(group = "SNV",            region = "Var_FL",  metric = "pLDDT"),
  WT_full_snv_control    %>% select(plddt_mean_var_FL) %>% mutate(group = "SNV_Control",    region = "Var_FL",  metric = "pLDDT"),
  WT_full_plus1          %>% select(plddt_mean_var_FL) %>% mutate(group = "Plus1",          region = "Var_FL",  metric = "pLDDT"),
  WT_full_plus1_control  %>% select(plddt_mean_var_FL) %>% mutate(group = "Plus1_Control",  region = "Var_FL",  metric = "pLDDT"),
  WT_full_minus1         %>% select(plddt_mean_var_FL) %>% mutate(group = "Minus1",         region = "Var_FL",  metric = "pLDDT"),
  WT_full_minus1_control %>% select(plddt_mean_var_FL) %>% mutate(group = "Minus1_Control", region = "Var_FL",  metric = "pLDDT")
) %>% rename(value = plddt_mean_var_FL)

plddt_WT_nmd <- bind_rows(
  WT_NMD_snv            %>% select(plddt_mean_WT_nmd)  %>% mutate(group = "SNV",            region = "WT_NMD"),
  WT_NMD_snv_control    %>% select(plddt_mean_WT_nmd)  %>% mutate(group = "SNV_Control",    region = "WT_NMD"),
  WT_NMD_plus1          %>% select(plddt_mean_WT_nmd)  %>% mutate(group = "Plus1",          region = "WT_NMD"),
  WT_NMD_plus1_control  %>% select(plddt_mean_WT_nmd)  %>% mutate(group = "Plus1_Control",  region = "WT_NMD"),
  WT_NMD_minus1         %>% select(plddt_mean_WT_nmd)  %>% mutate(group = "Minus1",         region = "WT_NMD"),
  WT_NMD_minus1_control %>% select(plddt_mean_WT_nmd)  %>% mutate(group = "Minus1_Control", region = "WT_NMD")
) %>% rename(value = plddt_mean_WT_nmd)

plddt_var_nmd <- bind_rows(
  WT_NMD_snv            %>% select(plddt_mean_vars_nmd) %>% mutate(group = "SNV",            region = "Var_NMD"),
  WT_NMD_snv_control    %>% select(plddt_mean_vars_nmd) %>% mutate(group = "SNV_Control",    region = "Var_NMD"),
  WT_NMD_plus1          %>% select(plddt_mean_vars_nmd) %>% mutate(group = "Plus1",          region = "Var_NMD"),
  WT_NMD_plus1_control  %>% select(plddt_mean_vars_nmd) %>% mutate(group = "Plus1_Control",  region = "Var_NMD"),
  WT_NMD_minus1         %>% select(plddt_mean_vars_nmd) %>% mutate(group = "Minus1",         region = "Var_NMD"),
  WT_NMD_minus1_control %>% select(plddt_mean_vars_nmd) %>% mutate(group = "Minus1_Control", region = "Var_NMD")
) %>% rename(value = plddt_mean_vars_nmd)
plddt_WT_FL  <- plddt_WT_FL  %>% mutate(type = "WT")
plddt_var_FL <- plddt_var_FL %>% mutate(type = "Var")

plddt_WT_nmd  <- plddt_WT_nmd  %>% mutate(type = "WT",  metric = "pLDDT")
plddt_var_nmd <- plddt_var_nmd %>% mutate(type = "Var", metric = "pLDDT")

# --- FL-only plot (two panels: WT_FL and Var_FL) ---
plddt_all_FL <- bind_rows(plddt_WT_FL, plddt_var_FL)

p_plddt_FL <- ggplot(plddt_all_FL, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  facet_wrap(~ region + type, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  labs(title = "Mean pLDDT (Full-Length)", x = "Group", y = "Mean pLDDT") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggpubr::stat_compare_means(
    comparisons = my_comparisons, method = "wilcox.test",
    label = "p.signif", tip.length = 0.01
  )
print(p_plddt_FL)
ggsave("~/Downloads/plddt_mean_FL.png", p_plddt_FL, width = 12, height = 4.5, dpi = 300)

# --- NMD-only plot (two panels: WT_NMD and Var_NMD) ---
plddt_all_NMD <- bind_rows(plddt_WT_nmd, plddt_var_nmd)

p_plddt_NMD <- ggplot(plddt_all_NMD, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  facet_wrap(~ region + type, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  labs(title = "Mean pLDDT (NMD Regions)", x = "Group", y = "Mean pLDDT") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggpubr::stat_compare_means(
    comparisons = my_comparisons, method = "wilcox.test",
    label = "p.signif", tip.length = 0.01
  ).
print(p_plddt_NMD)
ggsave("~/Downloads/plddt_mean_NMD.png", p_plddt_NMD, width = 12, height = 4.5, dpi = 300)

write.csv(plddt_all_NMD, "~/Downloads/plddt_all_NMD.csv", row.names = FALSE)
write.csv(plddt_all_FL, "~/Downloads/plddt_all_FL.csv", row.names = FALSE)

#6. plot nmd length change
##use nmd_length_loss_in_WT in WT_NMD_snv
nmd_length_loss <- dplyr::bind_rows(
  WT_NMD_minus1         %>% dplyr::select(nmd_length_loss_in_WT) %>% dplyr::mutate(group = "Minus1"),
  WT_NMD_minus1_control %>% dplyr::select(nmd_length_loss_in_WT) %>% dplyr::mutate(group = "Minus1_Control"),
  WT_NMD_plus1          %>% dplyr::select(nmd_length_loss_in_WT) %>% dplyr::mutate(group = "Plus1"),
  WT_NMD_plus1_control  %>% dplyr::select(nmd_length_loss_in_WT) %>% dplyr::mutate(group = "Plus1_Control"),
  WT_NMD_snv            %>% dplyr::select(nmd_length_loss_in_WT) %>% dplyr::mutate(group = "SNV"),
  WT_NMD_snv_control    %>% dplyr::select(nmd_length_loss_in_WT) %>% dplyr::mutate(group = "SNV_Control")
)

p_nmd_length_loss <- ggplot(nmd_length_loss,
                            aes(x = group, y = nmd_length_loss_in_WT, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  scale_fill_manual(values = group_colors) +
  labs(title = "NMD Length Loss in WT",
       x = "Group", y = "NMD Length Loss (aa)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggpubr::stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.01
  )

print(p_nmd_length_loss)
ggsave("~/Downloads/nmd_length_loss_in_WT.png",
       p_nmd_length_loss, width = 8, height = 5, dpi = 300)


#7. plot pI_var_FL# --- pI plots like pLDDT: WT vs Var × (FL & NMD) ----------------------------

# Safety: colors/comparisons
if (!exists("group_colors")) {
  group_colors <- c(
    "Minus1" = "#1f77b4", "Minus1_Control" = "#aec7e8",
    "Plus1"  = "#ff7f0e", "Plus1_Control"  = "#ffbb78",
    "SNV"    = "#2ca02c", "SNV_Control"    = "#98df8a"
  )
}
if (!exists("my_comparisons")) {
  my_comparisons <- list(
    c("Minus1","Minus1_Control"),
    c("Plus1","Plus1_Control"),
    c("SNV","SNV_Control")
  )
}

# ---------- FL: build WT & Var frames ----------
pI_WT_FL <- dplyr::bind_rows(
  WT_full_minus1         %>% dplyr::select(pI_WT_FL)  %>% dplyr::mutate(group="Minus1",         region="WT_FL",  type="WT"),
  WT_full_minus1_control %>% dplyr::select(pI_WT_FL)  %>% dplyr::mutate(group="Minus1_Control", region="WT_FL",  type="WT"),
  WT_full_plus1          %>% dplyr::select(pI_WT_FL)  %>% dplyr::mutate(group="Plus1",          region="WT_FL",  type="WT"),
  WT_full_plus1_control  %>% dplyr::select(pI_WT_FL)  %>% dplyr::mutate(group="Plus1_Control",  region="WT_FL",  type="WT"),
  WT_full_snv            %>% dplyr::select(pI_WT_FL)  %>% dplyr::mutate(group="SNV",            region="WT_FL",  type="WT"),
  WT_full_snv_control    %>% dplyr::select(pI_WT_FL)  %>% dplyr::mutate(group="SNV_Control",    region="WT_FL",  type="WT")
) %>% dplyr::rename(value = pI_WT_FL)

pI_var_FL <- dplyr::bind_rows(
  WT_full_minus1         %>% dplyr::select(pI_var_FL) %>% dplyr::mutate(group="Minus1",         region="Var_FL", type="Var"),
  WT_full_minus1_control %>% dplyr::select(pI_var_FL) %>% dplyr::mutate(group="Minus1_Control", region="Var_FL", type="Var"),
  WT_full_plus1          %>% dplyr::select(pI_var_FL) %>% dplyr::mutate(group="Plus1",          region="Var_FL", type="Var"),
  WT_full_plus1_control  %>% dplyr::select(pI_var_FL) %>% dplyr::mutate(group="Plus1_Control",  region="Var_FL", type="Var"),
  WT_full_snv            %>% dplyr::select(pI_var_FL) %>% dplyr::mutate(group="SNV",            region="Var_FL", type="Var"),
  WT_full_snv_control    %>% dplyr::select(pI_var_FL) %>% dplyr::mutate(group="SNV_Control",    region="Var_FL", type="Var")
) %>% dplyr::rename(value = pI_var_FL)

pI_all_FL <- dplyr::bind_rows(pI_WT_FL, pI_var_FL)

p_pI_FL <- ggplot(pI_all_FL, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  facet_wrap(~ region + type, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  labs(title = "Isoelectric Point (Full-Length)", x = "Group", y = "pI") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggpubr::stat_compare_means(
    comparisons = my_comparisons, method = "wilcox.test",
    label = "p.signif", tip.length = 0.01
  )
print(p_pI_FL)
ggsave("~/Downloads/pI_FL_faceted.png", p_pI_FL, width = 12, height = 4.5, dpi = 300)

# ---------- NMD: build WT & Var frames (if columns exist) ----------
has_pI_NMD <- all(c("pI_WT_nmd","pI_vars_nmd") %in% names(WT_NMD_snv))

if (has_pI_NMD) {
  pI_WT_nmd <- dplyr::bind_rows(
    WT_NMD_minus1         %>% dplyr::select(pI_WT_nmd)  %>% dplyr::mutate(group="Minus1",         region="WT_NMD", type="WT"),
    WT_NMD_minus1_control %>% dplyr::select(pI_WT_nmd)  %>% dplyr::mutate(group="Minus1_Control", region="WT_NMD", type="WT"),
    WT_NMD_plus1          %>% dplyr::select(pI_WT_nmd)  %>% dplyr::mutate(group="Plus1",          region="WT_NMD", type="WT"),
    WT_NMD_plus1_control  %>% dplyr::select(pI_WT_nmd)  %>% dplyr::mutate(group="Plus1_Control",  region="WT_NMD", type="WT"),
    WT_NMD_snv            %>% dplyr::select(pI_WT_nmd)  %>% dplyr::mutate(group="SNV",            region="WT_NMD", type="WT"),
    WT_NMD_snv_control    %>% dplyr::select(pI_WT_nmd)  %>% dplyr::mutate(group="SNV_Control",    region="WT_NMD", type="WT")
  ) %>% dplyr::rename(value = pI_WT_nmd)
  
  pI_var_nmd <- dplyr::bind_rows(
    WT_NMD_minus1         %>% dplyr::select(pI_vars_nmd) %>% dplyr::mutate(group="Minus1",         region="Var_NMD", type="Var"),
    WT_NMD_minus1_control %>% dplyr::select(pI_vars_nmd) %>% dplyr::mutate(group="Minus1_Control", region="Var_NMD", type="Var"),
    WT_NMD_plus1          %>% dplyr::select(pI_vars_nmd) %>% dplyr::mutate(group="Plus1",          region="Var_NMD", type="Var"),
    WT_NMD_plus1_control  %>% dplyr::select(pI_vars_nmd) %>% dplyr::mutate(group="Plus1_Control",  region="Var_NMD", type="Var"),
    WT_NMD_snv            %>% dplyr::select(pI_vars_nmd) %>% dplyr::mutate(group="SNV",            region="Var_NMD", type="Var"),
    WT_NMD_snv_control    %>% dplyr::select(pI_vars_nmd) %>% dplyr::mutate(group="SNV_Control",    region="Var_NMD", type="Var")
  ) %>% dplyr::rename(value = pI_vars_nmd)
  
  pI_all_NMD <- dplyr::bind_rows(pI_WT_nmd, pI_var_nmd)
  
  p_pI_NMD <- ggplot(pI_all_NMD, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
    facet_wrap(~ region + type, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = group_colors) +
    labs(title = "Isoelectric Point (NMD Regions)", x = "Group", y = "pI") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons, method = "wilcox.test",
      label = "p.signif", tip.length = 0.01
    )
  print(p_pI_NMD)
  ggsave("~/Downloads/pI_NMD_faceted.png", p_pI_NMD, width = 12, height = 4.5, dpi = 300)
  

write.csv(pI_all_NMD, "~/Downloads/pI_all_NMD.csv", row.names = FALSE)
write.csv(pI_all_FL, "~/Downloads/pI_all_FL.csv", row.names = FALSE)
#8. plot aromaticity# --- Aromaticity plots like pLDDT/pI: WT vs Var × (FL & NMD) ----------------
  
  # Safety: colors/comparisons
  if (!exists("group_colors")) {
    group_colors <- c(
      "Minus1" = "#1f77b4", "Minus1_Control" = "#aec7e8",
      "Plus1"  = "#ff7f0e", "Plus1_Control"  = "#ffbb78",
      "SNV"    = "#2ca02c", "SNV_Control"    = "#98df8a"
    )
  }
  if (!exists("my_comparisons")) {
    my_comparisons <- list(
      c("Minus1","Minus1_Control"),
      c("Plus1","Plus1_Control"),
      c("SNV","SNV_Control")
    )
  }
  
  # ---------- FL: build WT & Var frames ----------
  arom_WT_FL <- dplyr::bind_rows(
    WT_full_minus1         %>% dplyr::select(aromaticity_WT_FL) %>% dplyr::mutate(group="Minus1",         region="WT_FL",  type="WT"),
    WT_full_minus1_control %>% dplyr::select(aromaticity_WT_FL) %>% dplyr::mutate(group="Minus1_Control", region="WT_FL",  type="WT"),
    WT_full_plus1          %>% dplyr::select(aromaticity_WT_FL) %>% dplyr::mutate(group="Plus1",          region="WT_FL",  type="WT"),
    WT_full_plus1_control  %>% dplyr::select(aromaticity_WT_FL) %>% dplyr::mutate(group="Plus1_Control",  region="WT_FL",  type="WT"),
    WT_full_snv            %>% dplyr::select(aromaticity_WT_FL) %>% dplyr::mutate(group="SNV",            region="WT_FL",  type="WT"),
    WT_full_snv_control    %>% dplyr::select(aromaticity_WT_FL) %>% dplyr::mutate(group="SNV_Control",    region="WT_FL",  type="WT")
  ) %>% dplyr::rename(value = aromaticity_WT_FL)
  
  arom_var_FL <- dplyr::bind_rows(
    WT_full_minus1         %>% dplyr::select(aromaticity_var_FL) %>% dplyr::mutate(group="Minus1",         region="Var_FL", type="Var"),
    WT_full_minus1_control %>% dplyr::select(aromaticity_var_FL) %>% dplyr::mutate(group="Minus1_Control", region="Var_FL", type="Var"),
    WT_full_plus1          %>% dplyr::select(aromaticity_var_FL) %>% dplyr::mutate(group="Plus1",          region="Var_FL", type="Var"),
    WT_full_plus1_control  %>% dplyr::select(aromaticity_var_FL) %>% dplyr::mutate(group="Plus1_Control",  region="Var_FL", type="Var"),
    WT_full_snv            %>% dplyr::select(aromaticity_var_FL) %>% dplyr::mutate(group="SNV",            region="Var_FL", type="Var"),
    WT_full_snv_control    %>% dplyr::select(aromaticity_var_FL) %>% dplyr::mutate(group="SNV_Control",    region="Var_FL", type="Var")
  ) %>% dplyr::rename(value = aromaticity_var_FL)
  
  arom_all_FL <- dplyr::bind_rows(arom_WT_FL, arom_var_FL)
  
  p_arom_FL <- ggplot(arom_all_FL, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
    facet_wrap(~ region + type, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = group_colors) +
    labs(title = "Aromaticity (Full-Length)", x = "Group", y = "Aromaticity") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggpubr::stat_compare_means(
      comparisons = my_comparisons, method = "wilcox.test",
      label = "p.signif", tip.length = 0.01
    )
  print(p_arom_FL)
  ggsave("~/Downloads/aromaticity_FL_faceted.png", p_arom_FL, width = 12, height = 4.5, dpi = 300)
  
  # ---------- NMD: build WT & Var frames (auto-skip if missing) ----------
  has_arom_NMD <- all(c("aromaticity_WT_nmd","aromaticity_vars_nmd") %in% names(WT_NMD_snv))
  
  if (has_arom_NMD) {
    arom_WT_nmd <- dplyr::bind_rows(
      WT_NMD_minus1         %>% dplyr::select(aromaticity_WT_nmd)  %>% dplyr::mutate(group="Minus1",         region="WT_NMD", type="WT"),
      WT_NMD_minus1_control %>% dplyr::select(aromaticity_WT_nmd)  %>% dplyr::mutate(group="Minus1_Control", region="WT_NMD", type="WT"),
      WT_NMD_plus1          %>% dplyr::select(aromaticity_WT_nmd)  %>% dplyr::mutate(group="Plus1",          region="WT_NMD", type="WT"),
      WT_NMD_plus1_control  %>% dplyr::select(aromaticity_WT_nmd)  %>% dplyr::mutate(group="Plus1_Control",  region="WT_NMD", type="WT"),
      WT_NMD_snv            %>% dplyr::select(aromaticity_WT_nmd)  %>% dplyr::mutate(group="SNV",            region="WT_NMD", type="WT"),
      WT_NMD_snv_control    %>% dplyr::select(aromaticity_WT_nmd)  %>% dplyr::mutate(group="SNV_Control",    region="WT_NMD", type="WT")
    ) %>% dplyr::rename(value = aromaticity_WT_nmd)
    
    arom_var_nmd <- dplyr::bind_rows(
      WT_NMD_minus1         %>% dplyr::select(aromaticity_vars_nmd) %>% dplyr::mutate(group="Minus1",         region="Var_NMD", type="Var"),
      WT_NMD_minus1_control %>% dplyr::select(aromaticity_vars_nmd) %>% dplyr::mutate(group="Minus1_Control", region="Var_NMD", type="Var"),
      WT_NMD_plus1          %>% dplyr::select(aromaticity_vars_nmd) %>% dplyr::mutate(group="Plus1",          region="Var_NMD", type="Var"),
      WT_NMD_plus1_control  %>% dplyr::select(aromaticity_vars_nmd) %>% dplyr::mutate(group="Plus1_Control",  region="Var_NMD", type="Var"),
      WT_NMD_snv            %>% dplyr::select(aromaticity_vars_nmd) %>% dplyr::mutate(group="SNV",            region="Var_NMD", type="Var"),
      WT_NMD_snv_control    %>% dplyr::select(aromaticity_vars_nmd) %>% dplyr::mutate(group="SNV_Control",    region="Var_NMD", type="Var")
    ) %>% dplyr::rename(value = aromaticity_vars_nmd)
    
    arom_all_NMD <- dplyr::bind_rows(arom_WT_nmd, arom_var_nmd)
    
    p_arom_NMD <- ggplot(arom_all_NMD, aes(x = group, y = value, fill = group)) +
      geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
      facet_wrap(~ region + type, nrow = 1, scales = "free_y") +
      scale_fill_manual(values = group_colors) +
      labs(title = "Aromaticity (NMD Regions)", x = "Group", y = "Aromaticity") +
      theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggpubr::stat_compare_means(
        comparisons = my_comparisons, method = "wilcox.test",
        label = "p.signif", tip.length = 0.01
      )
    print(p_arom_NMD)
    ggsave("~/Downloads/aromaticity_NMD_faceted.png", p_arom_NMD, width = 12, height = 4.5, dpi = 300)

    write.csv(arom_all_NMD, "~/Downloads/aromaticity_all_NMD.csv", row.names = FALSE)   
    write.csv(arom_all_FL, "~/Downloads/aromaticity_all_FL.csv", row.names = FALSE)


    
    