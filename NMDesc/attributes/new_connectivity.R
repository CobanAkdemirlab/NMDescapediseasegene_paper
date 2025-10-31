library(STRINGdb)
library(igraph)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Step 1: Initialize STRING
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400, input_directory = "")

# Step 2: Map pli_all$gene to STRING
mapped_genes <- string_db$map(pli_all, "gene", removeUnmappedRows = TRUE)

# Merge back category
mapped_genes <- left_join(mapped_genes, pli_all, by = c("gene" = "gene"))

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
  select(STRING_id, gene, category.x) %>%
  distinct()
mapped_genes_degree$category <- mapped_genes_degree$category.x
centrality_data <- left_join(degree_df, mapped_genes_degree, by = "STRING_id") %>%
  filter(!is.na(category))  # Remove unmapped categories

# Step 7: Order category
centrality_data$category <- factor(
  centrality_data$category,
  levels = c("minus1", "minus1_Control", "plus1", "plus1_Control", "SNV", "SNV_Control")
)

# Step 8: Plot
ggplot(centrality_data, aes(x = category, y = Degree, fill = category)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  stat_compare_means(
    method = "t.test",
    comparisons = list(
      c("plus1", "plus1_Control"),
      c("minus1", "minus1_Control"),
      c("SNV", "SNV_Control")
    ),
    label = "p.signif"
  ) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = c(
    "minus1" = "#1f77b4",
    "minus1_Control" = "#aec7e8",
    "plus1" = "#ff7f0e",
    "plus1_Control" = "#ffbb78",
    "SNV" = "#2ca02c",
    "SNV_Control" = "#98df8a"
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

AD_centrality_data = centrality_data[which(centrality_data$gene %in% pli_AD$gene),]
g1 = ggplot(AD_centrality_data, aes(x = category, y = Degree, fill = category)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  ggpubr::stat_compare_means(
    method = "t.test",
    size = 6,
    comparisons = list(
      c("plus1", "plus1_Control"),
      c("minus1", "minus1_Control"),
      c("SNV", "SNV_Control")
    ),
    label = "p.signif"
  ) +
  scale_y_continuous(trans = "log10") +
  scale_fill_manual(values = c(
    "minus1" = "#1f77b4",
    "minus1_Control" = "#aec7e8",
    "plus1" = "#ff7f0e",
    "plus1_Control" = "#ffbb78",
    "SNV" = "#2ca02c",
    "SNV_Control" = "#98df8a"
  )) +
  scale_x_discrete(labels = c(
    "minus1" = "Minus1",
    "minus1_Control" = "Minus1_Control",
    "plus1" = "Plus1",
    "plus1_Control" = "Plus1_Control",
    "SNV" = "Nonsense",
    "SNV_Control" = "Nonsense_Control"
  )) +
  labs(
    title = "Degree Centrality of Genes in STRING PPI Network",
    y = "Degree Centrality (log10)",
    x = "Gene Category"
  ) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid       = element_blank(),
    axis.text.x  = element_text(angle = 45,size=12, hjust = 1, face = "bold"),
    axis.text.y  = element_text(size=12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none",
    plot.title = element_text(size=16, hjust = 0.5, face = "bold")
  ) +
  theme(panel.border = element_rect(colour = "black", fill = NA))
ggsave(g1, file = "AD_centrality_data.pdf")

