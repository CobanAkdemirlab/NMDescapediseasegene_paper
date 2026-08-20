library(ggplot2)
library(ggpubr)
library(dplyr)
library(readr)
# --- 路径解析层（自动插入）------------------------------------------------
# 原来这里是旧机器 /Users/jxu14/ 的绝对路径，换机器必断。改走 paths.R：
#   data_file("x.csv")  按文件名定位，找不到会明确报错
#   out_file("y.csv")   输出到 NMDESC_OUT（默认 ~/Desktop/NMDesc_out）
#   data_root("clinvar") 需要目录而非文件时用
.p <- c("gene level_v3/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../lib/paths.R", "../../../../lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p)) stop("找不到 paths.R —— 请从仓库根目录运行 R")
source(.p[1]); rm(.p)
# --------------------------------------------------------------------------


# Step 1: Load GTEx median TPM expression matrix
gtex <- read_tsv(missing_file("GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct", was = "/Users/jxu14/Downloads/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct"), skip = 2)

# Step 2: Clean Ensembl IDs and aggregate by gene symbol
gtex$Name <- sub("\\..*", "", gtex$Name)

# No need to re-load dplyr; just namespace the calls
gene_expr <- gtex %>%
  dplyr::select(-tidyselect::any_of("Name")) %>%          # drop Name if it exists
  dplyr::group_by(Description) %>%
  dplyr::summarise(
    dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::rename(Gene = Description)


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

# Step 7: Filter and add tau scores

tau_all <- gene_all %>%
  filter(hgnc_symbol %in% names(tau_scores)) %>%
  mutate(tau = tau_scores[hgnc_symbol])

# Step 9: Set group order for plotting
tau_all$group <- factor(tau_all$group,
                           levels = c("fs", "fs_control",
                                      "snv", "snv_control"))

# Step 10: Violin + boxplot

# Define custom colors
group_colors <- c(
  "fs" = "#1f77b4",
  "fs_control" = "#aec7e8",
  "snv" = "#2ca02c",
  "snv_control" = "#98df8a"
)

# Define comparisons
comparisons <- list(
  c("fs", "fs_control"),
  c("snv", "snv_control")
)

# Plot with comparisons
ggplot(tau_all, aes(x = group, y = tau, fill = group)) +
geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  ggpubr::stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = group_colors) +
  scale_x_discrete(labels = c(
    "fs" = "fs",
    "fs_control" = "fs_Control",
    "snv" = "Nonsense",
    "snv_control" = "Nonsense_Control"
  )) +
  theme_minimal() +
  labs(
    title = "Tissue specificity (Tau) across gene categories",
    y = "Tissue Specificity (Tau)",
    x = "Gene group"
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
theme(panel.border = element_rect(colour = "black", fill = NA))
write.csv(tau_all, "tau_all.csv", row.names = FALSE)
