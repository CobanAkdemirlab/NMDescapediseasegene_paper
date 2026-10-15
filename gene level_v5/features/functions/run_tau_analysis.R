# Required packages, checked immediately on source
for (.pkg in c(
  "dplyr",
  "ggplot2",
  "ggpubr",
  "readr"
)) if (!requireNamespace(.pkg, quietly = TRUE))
  stop(sprintf("Missing package %s -- install.packages(\"%s\") or BiocManager::install(\"%s\")", .pkg, .pkg, .pkg))
for (.pkg in c(
  "dplyr",
  "ggplot2",
  "ggpubr",
  "readr"
)) suppressPackageStartupMessages(library(.pkg, character.only = TRUE))
rm(.pkg)

run_tau_analysis <- function(gene_all, gtex_path, output_prefix = "tau") {

  # GTEx input is optional: returns tau as NA if missing.
  # Matches annotate_tau() behavior in gene_main_dbh.R.
  # gene_all_matched.csv / gene_all_random.csv / tau_all.csv
  # already contain computed tau; only needed to recompute with new GTEx.
  if (is.null(gtex_path) || length(gtex_path) != 1 || is.na(gtex_path) ||
      !file.exists(gtex_path)) {
    message("  Skipping tissue-specific tau -- GTEx file unavailable (tau set to NA)")
    # Return structure matches normal path: list(tau_all, tau_scores, plot)
    # otherwise caller's $tau_all would be NULL
    return(invisible(list(
      tau_all    = gene_all %>% dplyr::mutate(tau = NA_real_),
      tau_scores = NULL,
      plot       = NULL)))
  }

  # 1. GTEx median TPM
  gtex <- read_tsv(gtex_path, skip = 2)
  gtex$Name <- sub("\\..*", "", gtex$Name)
  
  gene_expr <- gtex %>%
    dplyr::select(-tidyselect::any_of("Name")) %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    dplyr::rename(Gene = Description)
  
  write.csv(gene_expr, paste0(output_prefix, "_gene_matrix.csv"), row.names = FALSE)
  
  # 2. Build expression matrix
  gene.df     <- as.data.frame(gene_expr)
  gene.matrix <- as.matrix(gene.df[, -1])
  rownames(gene.matrix) <- gene.df$Gene
  
  # 3. Tau (tissue specificity)
  tau <- function(x) {
    if (any(is.na(x))) stop("NA values found.")
    if (any(x < 0))    stop("Negative values found. Maybe data is log-transformed?")
    sum(1 - x / max(x)) / (length(x) - 1)
  }
  tau_scores <- apply(gene.matrix, 1, tau)
  
  # 4. Attach tau to gene_all
  tau_all <- gene_all %>%
    filter(hgnc_symbol %in% names(tau_scores)) %>%
    mutate(tau = tau_scores[hgnc_symbol])
  
  # 5. Group order / colors
  tau_all$group <- factor(tau_all$group, levels = c("fs", "fs_control", "snv", "snv_control"))
  group_colors <- c("fs" = "#1f77b4", "fs_control" = "#aec7e8",
                    "snv" = "#2ca02c", "snv_control" = "#98df8a")
  comparisons <- list(c("fs", "fs_control"), c("snv", "snv_control"))
  
  # 6. Violin plot
  p <- ggplot(tau_all, aes(x = group, y = tau, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") +
    scale_fill_manual(values = group_colors) +
    scale_x_discrete(labels = c("fs" = "fs", "fs_control" = "fs_Control",
                                "snv" = "Nonsense", "snv_control" = "Nonsense_Control")) +
    theme_minimal() +
    labs(title = "Tissue specificity (Tau) across gene categories",
         y = "Tissue Specificity (Tau)", x = "Gene group") +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          legend.position = "none",
          plot.title   = element_text(hjust = 0.5, face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA))
  print(p)
  
  # 7. Save
  write.csv(tau_all, paste0(output_prefix, "_all.csv"), row.names = FALSE)
  ggsave(paste0(output_prefix, "_violin_plot.png"), p, width = 8, height = 6, dpi = 300)
  
  invisible(list(tau_all = tau_all, tau_scores = tau_scores, plot = p))
}
