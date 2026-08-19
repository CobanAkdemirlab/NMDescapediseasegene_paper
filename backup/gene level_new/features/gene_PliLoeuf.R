#This Rscript is to plot pli/loeuf, and exon count etc info on gene level

library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)

# Load gnomAD constraint metrics
Lof_metrics <- read.delim("/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/others/gnomad.v2.1.1.lof_metrics.by_gene.txt")
# Merge with constraint metrics
pli_all = merge(gene_all, Lof_metrics,by.x = 'hgnc_symbol', by.y = 'gene') %>%
  dplyr::select(hgnc_symbol, pLI, oe_lof_upper,group) 


#2. define ggplot paras
# Define comparisons
comparisons <- list(
  c("fs", "fs_control"),
  c("snv", "snv_control")
)

# Custom colors
group_colors <- c(
  "fs" = "#ff7f0e",
  "fs_control" = "#ffbb78",  # lighter orange
  "snv" = "#2ca02c",
  "snv_control" = "#98df8a"     # lighter green
)

#3. plot pli and loeuf

pli_all$pli_cat = cut(pli_all$pLI, breaks = c(0,0.35,0.66,1), labels = c('Low','Medium','High'))
pli_all$loefu_cat = cut(pli_all$oe_lof_upper, breaks = c(0,0.2,0.6,2), labels = c('Low','Medium','High')) #level of intolerance
write_csv(pli_all, "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/others/pli_all_cate.csv")
# Plot pLI
# Create the plot object
pli_plot <- ggplot(pli_all, aes(x = group, y = pLI, fill = group)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  labs(title = "pLI distribution by gene group",
       y = "pLI",
       x = "group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")
ggplot(pli_AD, aes(x = group, y = pLI, fill = group)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  labs(title = "pLI distribution in AD genes",
       y = "pLI",
       x = "group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")
ggsave("pli_distribution_by_group.pdf", plot = pli_plot, width = 8, height = 5)

loeuf_plot <- ggplot(pli_all, aes(x = group, y = oe_lof_upper, fill = group)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  labs(title = "LOEUF distribution by gene group",
       y = "LOEUF",
       x = "group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")
loeuf_plot <- ggplot(pli_AD, aes(x = group, y = oe_lof_upper, fill = group)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  labs(title = "LOEUF distribution in AD genes",
       y = "LOEUF",
       x = "group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")

ggsave("loeuf_distribution.pdf", plot = loeuf_plot,width = 8, height = 5)

#4. plot cds_length
cds_plot = ggplot(gene_all, aes(x = group, y = cds_length, fill = group)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  labs(title = "CDS Length Distribution in AD genes",
       y = "CDS Length (bp)",
       x = "group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_y_log10() +
  scale_fill_manual(values = group_colors) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")
ggsave("cds_length_distribution.pdf", plot = cds_plot, width = 8, height = 5)

#5. plot exon count
exons <- getBM(
  attributes = c("ensembl_transcript_id", "cds_start", "cds_end", "rank", "strand"),
  filters = "ensembl_transcript_id",
  values = gene_all$ensembl_transcript_id,
  mart = ensembl
)
# add gene name to exon info
exons_joined <- exons %>%
  inner_join(gene_all, by = "ensembl_transcript_id")
exon_counts <- exons_joined %>%
  group_by(ensembl_transcript_id, hgnc_symbol, strand) %>%
  arrange(rank) %>%
  summarise(
    exon_num = max(rank),
    .groups = "drop"
  ) %>%
  dplyr::select(hgnc_symbol,exon_num)
#get group info from pli_all
exon_counts <- exon_counts %>%
  left_join(gene_all %>% dplyr::select(hgnc_symbol, group), by = "hgnc_symbol")

exon_plot = ggplot(exon_counts, aes(x = group, y = exon_num, fill = group)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  scale_y_log10() +
  labs(title = "Exon count in AD Genes",
       y = "number",
       x = "group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")


# 1) Rename groups to Title Case (same mapper you used before)
rename_levels <- function(x) {
  recoded <- dplyr::recode(tolower(x),
                           "minus1" = "Minus1",
                           "minus1_control" = "Minus1_control",
                           "fs" = "fs",
                           "fs_control" = "fs_control",
                           "snv" = "Nonsense",
                           "snv_control" = "Nonsense_control",
                           .default = x
  )
  factor(recoded, levels = c("Minus1","Minus1_control",
                             "fs","fs_control",
                             "Nonsense","Nonsense_control"))
}

exon_counts2 <- exon_counts %>%
  mutate(group = rename_levels(group)) %>%
  filter(!is.na(group)) %>%
  droplevels()

# 2) Rebuild comparisons to match new labels and available groups
# If you already have `comparisons` in old labels, convert them:
comparisons_new <- lapply(comparisons, function(p) as.character(rename_levels(p)))

present_groups <- levels(exon_counts2$group)
comparisons_use <- Filter(function(p) length(p) == 2 && all(p %in% present_groups), comparisons_new)

# 3) Plot (with bold/centered title and bold x-label & ticks)
group_colors2 <- c(
  Minus1="#1f77b4", Minus1_control="#aec7e8",
  fs ="#ff7f0e", fs_control ="#ffbb78",
  Nonsense="#2ca02c", Nonsense_control="#98df8a"
)

p_exon <- ggplot(exon_counts2, aes(x = group, y = exon_num, fill = group)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  scale_y_log10() +
  scale_fill_manual(values = group_colors2, guide = "none") +
  labs(title = "Exon count by Gene Group",
       y = "Number",
       x = "Gene group") +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  ggpubr::stat_compare_means(comparisons = comparisons_use, method = "wilcox.test")

print(p_exon)

ggsave("exon_count_distribution.pdf", plot = exon_plot, width = 8, height = 5)




