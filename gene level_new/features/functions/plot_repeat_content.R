plot_repeat_content <- function(
    gene_all,
    comparisons = list(
      c("fs", "fs_control"),
      c("snv", "snv_control")
    ),
    output_csv = "repeat_content.csv",
    output_fig = "repeat_content.png"
) {
  library(ggplot2)
  library(ggpubr)
  library(stringr)
  
  # ── Helpers ───────────────────────────────────────────────────────────────
  detect_repeats <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(0)
    matches <- str_extract_all(sequence, "([ATGC]{1,6})\\1+")[[1]]
    if (length(matches) == 0) return(0)
    sum(nchar(matches))
  }
  
  repeat_fraction <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(NA_real_)
    detect_repeats(sequence) / nchar(sequence)
  }
  
  detect_homopolymer <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(0)
    matches <- str_extract_all(sequence, "([ATGC])\\1{3,}")[[1]]
    if (length(matches) == 0) return(0)
    sum(nchar(matches))
  }
  
  homopolymer_fraction <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(NA_real_)
    detect_homopolymer(sequence) / nchar(sequence)
  }
  
  # ── Compute features ──────────────────────────────────────────────────────
  gene_all$repeat_fraction          <- sapply(gene_all$coding,      repeat_fraction)
  gene_all$nmdesc_repeat_fraction   <- sapply(gene_all$nmdesc_cds,  repeat_fraction)
  gene_all$homopolymer_fraction     <- sapply(gene_all$coding,      homopolymer_fraction)
  gene_all$nmdesc_homopolymer_fraction <- sapply(gene_all$nmdesc_cds, homopolymer_fraction)
  
  # ── Shared plot aesthetics ────────────────────────────────────────────────
  group_colours <- c(
    "snv"         = "#1f77b4",
    "snv_control" = "#aec7e8",
    "fs"          = "#ff7f0e",
    "fs_control"  = "#ffbb78"
  )
  
  make_plot <- function(y_var, y_label, title) {
    ggplot(gene_all, aes(x = group, y = .data[[y_var]], fill = group)) +
      geom_boxplot(width = 0.7, outlier.shape = 16, outlier.size = 1.5) +
      stat_compare_means(
        comparisons = comparisons,
        method      = "wilcox.test",
        label       = "p.format"
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      scale_fill_manual(values = group_colours) +
      labs(x = "", y = y_label, title = title) +
      theme_bw() +
      theme(
        plot.title      = element_text(hjust = 0.5),
        legend.position = "none"
      )
  }
  
  # ── Four panels ───────────────────────────────────────────────────────────
  p1 <- make_plot("repeat_fraction",             "Repeat fraction",      "Overall Repeat Content")
  p2 <- make_plot("nmdesc_repeat_fraction",      "Repeat fraction",      "NMD-escape Repeat Content")
  p3 <- make_plot("homopolymer_fraction",        "Homopolymer fraction", "Overall Homopolymer Content")
  p4 <- make_plot("nmdesc_homopolymer_fraction", "Homopolymer fraction", "NMD-escape Homopolymer Content")
  
  combined <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
  
  # ── Save outputs ──────────────────────────────────────────────────────────
  out_cols <- c("hgnc_symbol", "ensembl_transcript_id", "group",
                "repeat_fraction", "nmdesc_repeat_fraction",
                "homopolymer_fraction", "nmdesc_homopolymer_fraction")
  write.csv(gene_all[, out_cols], output_csv, row.names = FALSE)
  message("CSV saved to: ", output_csv)
  
  ggsave(output_fig, plot = combined, width = 12, height = 10, dpi = 300)
  message("Figure saved to: ", output_fig)
  
  invisible(list(plot = combined, data = gene_all))
}


# ── Usage ─────────────────────────────────────────────────────────────────
result <- plot_repeat_content(
  gene_all   = gene_all,
  output_csv = "repeat_content.csv",
  output_fig = "repeat_content.png"
)