# Required packages (checked on source; missing packages fail immediately)
for (.pkg in c(
  "STRINGdb",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "igraph"
)) if (!requireNamespace(.pkg, quietly = TRUE))
  stop(sprintf("Missing package %s -- install.packages(\"%s\") or BiocManager::install(\"%s\")", .pkg, .pkg, .pkg))
for (.pkg in c(
  "STRINGdb",
  "dplyr",
  "ggplot2",
  "ggpubr",
  "igraph"
)) suppressPackageStartupMessages(library(.pkg, character.only = TRUE))
rm(.pkg)

calculate_ppi_degree_centrality <- function(
    gene_input,
    output_csv      = "centrality_data.csv",
    output_fig      = "centrality_plot.png",
    string_version  = "11.5",
    species         = 9606,
    score_threshold = 400
) {
  # STRING setup
  string_db <- STRINGdb$new(
    version         = string_version,
    species         = species,
    score_threshold = score_threshold,
    input_directory = ""
  )
  
  # Map genes to STRING IDs. string_db$map() keeps the input columns
  # (including `group`) and adds STRING_id, so no second join is needed.
  mapped_genes <- string_db$map(gene_input, "hgnc_symbol", removeUnmappedRows = TRUE)
  
  # Interactions -> undirected network
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  ppi_network <- graph_from_data_frame(
    data.frame(from = interactions$from, to = interactions$to),
    directed = FALSE
  )
  ppi_network <- simplify(ppi_network)
  
  # Degree centrality
  degree_scores <- degree(ppi_network, normalized = TRUE)
  degree_df     <- data.frame(STRING_id = names(degree_scores), Degree = degree_scores)
  
  # Merge with group labels
  mapped_genes_degree <- mapped_genes %>%
    dplyr::select(STRING_id, hgnc_symbol, group) %>%
    dplyr::distinct()
  
  centrality_data <- left_join(degree_df, mapped_genes_degree, by = "STRING_id") %>%
    dplyr::filter(!is.na(group))
  
  # Plot
  p <- ggplot(centrality_data, aes(x = group, y = Degree, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    stat_compare_means(
      method      = "t.test",
      comparisons = list(c("fs", "fs_control"), c("snv", "snv_control")),
      label       = "p.signif"
    ) +
    scale_y_continuous(trans = "log10") +
    scale_fill_manual(values = CONFIG$group_colors) +
    labs(title = "Degree Centrality of Genes in STRING PPI Network",
         y = "Degree Centrality (log10)", x = "Gene Category") +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.title.y    = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )
  
  write.csv(centrality_data, output_csv, row.names = FALSE)
  message("CSV saved to: ", output_csv)
  ggsave(output_fig, plot = p, width = 7, height = 6, dpi = 300)
  message("Figure saved to: ", output_fig)
  
  invisible(list(plot = p, data = centrality_data))
}
