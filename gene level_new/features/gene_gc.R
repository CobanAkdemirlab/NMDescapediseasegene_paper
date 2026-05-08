#This R script is to get gene-level GC content

get_gc_content <- function(sequence) {
  bases <- strsplit(toupper(sequence), "")[[1]]
  gc <- sum(bases %in% c("G", "C"))
  return(round(gc / length(bases) * 100, 2))
}

##ovarall gene
gene_all$gc_content = sapply(gene_all$coding, get_gc_content)
##NMDesc region
gene_all$nmdesc_gc_content = sapply(gene_all$nmdesc_cds, get_gc_content)

# boxplot
ggplot(gene_all, aes(x = group, y = gc_content, fill = group)) +
  geom_boxplot(width = 0.7, outlier.shape = 16, outlier.size = 1.5) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format"
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    x = "",
    y = "GC content (%)",
    title = "GC content in NMDesc regions"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

