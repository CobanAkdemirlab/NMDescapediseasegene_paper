#This R script is to calculate ppi degree centrality from gene_all, which contains all disease and control genes
#from clinvar and gnomad, further more, plot the result

library(STRINGdb)
library(igraph)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Step 1: Initialize STRING
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400, input_directory = "")

# Step 2: Map gene_all$hgnc_symbol to STRING
gene_all_wald = data.frame(hgnc_symbol = c(snv_gene,fs_gene), group = c(rep("snv", length(snv_gene)), rep("fs", length(fs_gene))))
gene_all_control = gene_all %>%
  filter(group %in% c("snv_control", "fs_control")) %>%
  select(hgnc_symbol, group) %>%
  distinct()
gene_all2 = rbind(gene_all_wald, gene_all_control)
mapped_genes <- string_db$map(gene_all2, "hgnc_symbol", removeUnmappedRows = TRUE)

# Merge back category
mapped_genes <- left_join(mapped_genes, gene_all2, by = c("hgnc_symbol" = "hgnc_symbol"))

# Step 3: Get interaction edges
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

# Step 4: Build PPI network
edges <- data.frame(from = interactions$from, to = interactions$to)
ppi_network <- graph_from_data_frame(edges, directed = FALSE)
ppi_network <- simplify(ppi_network)

# Step 5: Compute degree centrality
degree_scores <- degree(ppi_network, normalized = TRUE)
degree_df <- data.frame(STRING_id = names(degree_scores), Degree = degree_scores)

# Step 6: Merge centrality scores with category
mapped_genes_degree <- mapped_genes %>%
  select(STRING_id, hgnc_symbol, group.x) %>%
  distinct()
mapped_genes_degree$group <- mapped_genes_degree$group.x
centrality_data <- left_join(degree_df, mapped_genes_degree, by = "STRING_id") %>%
  filter(!is.na(group))  # Remove unmapped categories

# Step 8: Plot
ggplot(centrality_data, aes(x = group, y = Degree, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(
      c("fs", "fs_control"),
      c("snv", "snv_control")
    ),
    label = "p.signif"
  ) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = c(
    "snv" = "#1f77b4",
    "snv_control" = "#aec7e8",
    "fs" = "#ff7f0e",
    "fs_control" = "#ffbb78"
  )) +
  labs(
    title = "Degree Centrality of Genes in STRING PPI Network",
    y = "Degree Centrality (log10)",
    x = "Gene Category"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )
write.csv(centrality_data, "centrality_data.csv", row.names = FALSE)
