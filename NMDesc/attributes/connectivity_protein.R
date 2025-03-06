#find connectivity in protein-protein interaction network
library(STRINGdb)
library(igraph)
library(ggplot2)
library(dplyr)


string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400, input_directory = "")
mapped_plus1_can_genes <- string_db$map(data.frame(gene = plus1_can_gene0115$V1), "gene", removeUnmappedRows = TRUE)
mapped_plus1_css_genes <- string_db$map(data.frame(gene = plus1_css_gene0115$V1), "gene", removeUnmappedRows = TRUE)
mapped_plus1_long_genes <- string_db$map(data.frame(gene = plus1_long_gene0115$V1), "gene", removeUnmappedRows = TRUE)
mapped_plus2_css_genes <- string_db$map(data.frame(gene = plus2_css_gene0115$V1), "gene", removeUnmappedRows = TRUE)
mapped_plus2_can_genes <- string_db$map(data.frame(gene = plus2_can_gene0115$V1), "gene", removeUnmappedRows = TRUE)
mapped_plus2_long_genes <- string_db$map(data.frame(gene = plus2_long_gene0115$V1), "gene", removeUnmappedRows = TRUE)
mapped_plus1_trigger_genes <- string_db$map(data.frame(gene = plus1_trigger_gene0115$V1), "gene", removeUnmappedRows = TRUE)
mapped_plus2_trigger_genes <- string_db$map(data.frame(gene = plus2_trigger_gene0115$V1), "gene", removeUnmappedRows = TRUE)
mapped_snv_can_genes <- string_db$map(data.frame(gene = snv_plp_ptc_p1120_NMDenriched2_all$X1), "gene", removeUnmappedRows = TRUE)
mapped_snv_css_genes <- string_db$map(data.frame(gene = snv_plp_ptc_p1120_NMDenriched2_css$X1), "gene", removeUnmappedRows = TRUE)
mapped_snv_long_genes <- string_db$map(data.frame(gene = snv_plp_ptc_p1120_NMDenriched2_long$X1), "gene", removeUnmappedRows = TRUE)
mapped_snv_trig_genes <- string_db$map(data.frame(gene = snv_plp_ptc_p1120_NMDenriched2_trig$X1), "gene", removeUnmappedRows = TRUE)
mapped_control_genes <- string_db$map(data.frame(gene = control$gene[-1]), "gene", removeUnmappedRows = TRUE)
# Get interaction data for the mapped genes

# Define the function
calculate_degree_centrality <- function(mapped_genes) {
  # Get interactions from STRING database
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  
  # Check if interactions are found
  if (nrow(interactions) == 0) {
    stop("No interactions found for the provided STRING IDs.")
  }
  
  # Create edges data frame
  edges <- data.frame(from = interactions$from, to = interactions$to)
  
  # Build PPI network
  ppi_network <- graph_from_data_frame(edges, directed = FALSE)
  ppi_network <- simplify(ppi_network)
  
  # Calculate degree centrality
  degree_centrality <- degree(ppi_network,normalized = T)

  # Convert to a data frame for output
  centrality_df <- data.frame(
    Node = names(degree_centrality),
    Degree = degree_centrality
  )
  return(centrality_df)
}
plus1_can_centrality_df <- calculate_degree_centrality(mapped_plus1_can_genes)
plus1_css_centrality_df <- calculate_degree_centrality(mapped_plus1_css_genes)
plus1_long_centrality_df <- calculate_degree_centrality(mapped_plus1_long_genes)
plus2_css_centrality_df <- calculate_degree_centrality(mapped_plus2_css_genes)
plus2_can_centrality_df <- calculate_degree_centrality(mapped_plus2_can_genes)
plus2_long_centrality_df <- calculate_degree_centrality(mapped_plus2_long_genes)
plus1_trigger_centrality_df <- calculate_degree_centrality(mapped_plus1_trigger_genes)
plus2_trigger_centrality_df <- calculate_degree_centrality(mapped_plus2_trigger_genes)
snv_can_centrality_df <- calculate_degree_centrality(mapped_snv_can_genes)
snv_css_centrality_df <- calculate_degree_centrality(mapped_snv_css_genes)
snv_long_centrality_df <- calculate_degree_centrality(mapped_snv_long_genes)
snv_trig_centrality_df <- calculate_degree_centrality(mapped_snv_trig_genes)
control_centrality_df <- calculate_degree_centrality(mapped_control_genes)


# Combine all centrality data frames into one
all_centrality_df <- bind_rows(
  plus1_can_centrality_df %>% mutate(Group = "plus1_can"),
  plus1_css_centrality_df %>% mutate(Group = "plus1_css"),
  plus1_long_centrality_df %>% mutate(Group = "plus1_long"),
  plus2_css_centrality_df %>% mutate(Group = "plus2_css"),
  plus2_can_centrality_df %>% mutate(Group = "plus2_can"),
  plus2_long_centrality_df %>% mutate(Group = "plus2_long"),
  plus1_trigger_centrality_df %>% mutate(Group = "plus1_trigger"),
  plus2_trigger_centrality_df %>% mutate(Group = "plus2_trigger"),
  snv_can_centrality_df %>% mutate(Group = "snv_can"),
  snv_css_centrality_df %>% mutate(Group = "snv_css"),
  snv_long_centrality_df %>% mutate(Group = "snv_long"),
  snv_trig_centrality_df %>% mutate(Group = "snv_trig"),
  control_centrality_df %>% mutate(Group = "control")
)
all_centrality_df <- all_centrality_df %>%
  mutate(
    group = case_when(
      grepl("can", Group) ~ "can",
      grepl("css", Group) ~ "css",
      grepl("long", Group) ~ "long",
      grepl("trig", Group) ~ "trigger",
      Group == "control_df" ~ "control"
    )
  )
all_centrality_df$Group <- factor(
  all_centrality_df$Group,
  levels = all_centrality_df %>%
    arrange(group, Group) %>%
    pull(Group) %>% 
    unique()
)

# Create a boxplot to compare Degree Centrality across groups
ggplot(all_centrality_df, aes(x = Group, y = Degree, fill = group)) +
  geom_violin() +
  #geom_jitter(width = 0.06, color = "darkblue", alpha = 0.7) +
  labs(
    title = "Comparison of Degree Centrality Across Groups",
    x = "Group",
    y = "Degree Centrality"
  ) +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3", name = "Group") +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")








int_net <- rba_string_interactions_network(ids = mapped_genes,
                                           species = 9606,
                                           required_score = 500)

# View degree centrality for all genes
print(degree_centrality)



graph = ppi_network
dc2 = centr_degree(ppi_network, mode = 'in')
V(graph)$indegree <- centr_degree(graph, mode = "in")$res

nodes <- get.data.frame(graph, what="vertices")
nodes <- data.frame(id = nodes$name, title = nodes$name, group = nodes$indegree, indegree = nodes$indegree)
setnames(nodes, "indegree", "in-degree centrality")
nodes <- nodes[order(nodes$id, decreasing = F),]

edges <- get.data.frame(graph, what="edges")[1:2]

visNetwork(nodes, edges, height = "500px", width = "100%") %>%
  visOptions(selectedBy = "in-degree centrality", highlightNearest = TRUE, nodesIdSelection = TRUE)%>%
  visPhysics(stabilization = FALSE)%>%
  visEdges(arrows = "to")