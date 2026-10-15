#!/usr/bin/env Rscript
# Figure S3 - paralog / isoform / complex, three panels on one canvas, tagged A-C.
# The shipped merge_paralog_isoform_complex.R looks for a ggplot named `p` and for a
# file `complex_comparison_unmatched.R`; the three scripts actually name their final
# plot `pA` and the third script is `complex_comparison.R`. Everything else is the
# shipped logic: source each script with ggsave stubbed, strip its titles, patchwork.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
SUB <- Sys.getenv("S3_SUBTITLE", "")
get_plot <- function(script, obj = "pA") {
  e <- new.env(parent = globalenv())
  e$ggsave <- function(...) invisible(NULL)
  sys.source(script, envir = e)
  if (!exists(obj, envir = e)) stop(sprintf("'%s' not found after sourcing %s", obj, script))
  get(obj, envir = e)
}
strip <- function(g) g + labs(title = NULL, subtitle = NULL, caption = NULL, tag = NULL) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank())
pA <- strip(get_plot("paralog_comparison.R"))
pB <- strip(get_plot("isoform_comparison.R"))
pC <- strip(get_plot("complex_comparison.R"))
merged <- (pA | pB | pC) +
  plot_annotation(tag_levels = "A",
    title = "Figure S3 - Paralog, isoform and protein-complex counts of NMD-escape disease genes and their length-matched controls",
    subtitle = SUB,
    caption = paste("A, B: paired Wilcoxon signed-rank, unadjusted P shown, sqrt y-axis.",
                    "C: Wilcoxon rank-sum (unpaired) on binary CORUM membership, Wilson 95% CI.",
                    "BH-adjusted P in the three *_comparison_stats.csv files."),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(colour = "grey35", size = 10),
                  plot.caption = element_text(colour = "grey35", size = 8.5))) &
  theme(plot.tag = element_text(face = "bold", size = 18))
ggsave("FigureS3_paralog_isoform_complex.png", merged, width = 18, height = 6, dpi = 200, bg = "white")
ggsave("FigureS3_paralog_isoform_complex.pdf", merged, width = 18, height = 6, bg = "white")
cat("saved FigureS3_paralog_isoform_complex.png / .pdf\n")
