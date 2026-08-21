plot_gc_content <- function(
    gene_all,
    comparisons = list(
      c("fs", "fs_control"),
      c("snv", "snv_control")
    ),
    output_csv = "gc_content.csv",
    output_fig = "gc_content.png"
) {
  library(ggplot2)
  library(ggpubr)
  
  # ── GC content helper ─────────────────────────────────────────────────────
  get_gc_content <- function(sequence) {
    if (is.na(sequence) || sequence == "") return(NA_real_)
    bases <- strsplit(toupper(sequence), "")[[1]]
    round(sum(bases %in% c("G", "C")) / length(bases) * 100, 2)
  }
  
  # ── Compute GC content ────────────────────────────────────────────────────
  gene_all$gc_content       <- sapply(gene_all$coding,      get_gc_content)
  gene_all$nmdesc_gc_content <- sapply(gene_all$nmdesc_cds, get_gc_content)
  
  # ── Plot ──────────────────────────────────────────────────────────────────
  make_plot <- function(y_var, y_label, title) {
    ggplot(gene_all, aes(x = group, y = .data[[y_var]], fill = group)) +
      geom_boxplot(width = 0.7, outlier.shape = 16, outlier.size = 1.5) +
      stat_compare_means(
        comparisons = comparisons,
        method      = "wilcox.test",
        label       = "p.format"
      ) +
      scale_y_continuous(labels = function(x) paste0(x, "%")) +
      scale_fill_manual(values = c(
        "snv"         = "#1f77b4",
        "snv_control" = "#aec7e8",
        "fs"          = "#ff7f0e",
        "fs_control"  = "#ffbb78"
      )) +
      labs(x = "", y = y_label, title = title) +
      theme_bw() +
      theme(
        plot.title      = element_text(hjust = 0.5),
        legend.position = "none"
      )
  }
  
  p1 <- make_plot("gc_content",        "GC content (%)", "Overall Gene GC Content")
  p2 <- make_plot("nmdesc_gc_content", "GC content (%)", "GC Content in NMD-escape Regions")
  
  combined <- ggarrange(p1, p2, ncol = 2, nrow = 1)
  
  # ── Save outputs ──────────────────────────────────────────────────────────
  write.csv(
    gene_all[, c("hgnc_symbol", "ensembl_transcript_id", "group",
                 "gc_content", "nmdesc_gc_content")],
    output_csv, row.names = FALSE
  )
  message("CSV saved to: ", output_csv)
  
  ggsave(output_fig, plot = combined, width = 10, height = 6, dpi = 300)
  message("Figure saved to: ", output_fig)
  
  invisible(list(plot = combined, data = gene_all))
}

