library(ggplot2)
library(dplyr)
library(readr)

# Step 1: Load GTEx median TPM expression matrix
gtex <- read_tsv("/Users/jxu14/Downloads/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct", skip = 2)

# Step 2: Clean Ensembl IDs and aggregate by gene symbol
gtex$Name <- sub("\\..*", "", gtex$Name)

gene_expr <- gtex %>%
  select(-Name) %>%
  group_by(Description) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  rename(Gene = Description)

# Step 3: Save matrix as CSV (optional)
write.csv(gene_expr, "gene.matrix.csv", row.names = FALSE)

# Step 4: Prepare gene matrix for tau calculation
gene.df <- read.csv("gene.matrix.csv", stringsAsFactors = FALSE)
gene.matrix <- as.matrix(gene.df[, -1])  # drop Gene column

rownames(gene.matrix) <- gene.df$Gene   # set gene names as rownames

# Step 5: Define tau function
tau <- function(x) {
  if (any(is.na(x))) stop("NA values found.")
  if (any(x < 0)) stop("Negative values found. Maybe data is log-transformed?")
  sum(1 - x / max(x)) / (length(x) - 1)
}

# Step 6: Compute tau score for each gene
tau_scores <- apply(gene.matrix, 1, tau)

# Step 7: Combine pli data
pli_all <- bind_rows(minus1_pli, plus1_pli, snv_pli,
                     minus1_control_pli, plus1_control_pli, snv_control_pli)

# Step 8: Filter and add tau scores
tau_all <- pli_all %>%
  filter(gene %in% names(tau_scores)) %>%
  mutate(tau = tau_scores[gene])

# Step 9: Set category order for plotting
tau_all$category <- factor(tau_all$category,
                           levels = c("minus1", "minus1_Control",
                                      "plus1", "plus1_Control",
                                      "SNV", "SNV_Control"))

# Step 10: Violin + boxplot
library(ggplot2)
library(ggpubr)

# Define custom colors
group_colors <- c(
  "minus1" = "#1f77b4",
  "minus1_Control" = "#aec7e8",
  "plus1" = "#ff7f0e",
  "plus1_Control" = "#ffbb78",
  "SNV" = "#2ca02c",
  "SNV_Control" = "#98df8a"
)

# Define comparisons
comparisons <- list(
  c("plus1", "plus1_Control"),
  c("minus1", "minus1_Control"),
  c("SNV", "SNV_Control")
)

# Plot with comparisons
ggplot(tau_all, aes(x = category, y = tau, fill = category)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = group_colors) +
  scale_x_discrete(labels = c(
    "plus1" = "Plus1",
    "plus1_Control" = "Plus1_Control",
    "minus1" = "Minus1",
    "minus1_Control" = "Minus1_Control",
    "SNV" = "Nonsense",
    "SNV_Control" = "Nonsense_Control"
  )) +
  theme_minimal() +
  labs(y = "Tissue Specificity (Tau)", x = "Category") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )
