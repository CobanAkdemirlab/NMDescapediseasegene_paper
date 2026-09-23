#!/usr/bin/env Rscript
# =============================================================================
# merge_paralog_isoform_complex.R
# Merge the paralog, isoform, and complex comparisons into ONE figure (A/B/C),
# with NO panel titles/subtitles and NO caption notes (clean panels only).
#
# This sources the three original scripts to reuse their data-loading and
# per-panel plot builders, then strips titles/subtitles/captions and combines
# them with patchwork, tagged A, B, C.
#
# Requirements (same as the originals, run locally):
#   - matched_pairs_full.csv                       (length-matched pairs)
#   - paralog_counts.csv / isoform_counts.csv      (biomaRt caches; auto-built)
#   - CORUM coreComplexes table (see complex_comparison.R: CORUM_PATH)
#   - R packages: ggplot2, patchwork, biomaRt (for the cache builds)
#
# Output: merged_paralog_isoform_complex.pdf / .png
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

# ---- 1. Build each panel by sourcing the originals in a child environment ----
# Each original script assigns its final ggplot to an object (p_paralog / p_iso /
# p_cplx below). We source with local=TRUE so their side effects stay contained,
# then pull out the ggplot object. Adjust the object names if yours differ.
get_plot <- function(script, plot_obj) {
  e <- new.env(parent = globalenv())
  sys.source(script, envir = e)
  if (!exists(plot_obj, envir = e))
    stop(sprintf("'%s' not found after sourcing %s \u2014 set plot_obj to the ggplot object name in that script.", plot_obj, script))
  get(plot_obj, envir = e)
}

# NOTE: set these to the ggplot variable names used at the end of each script.
p_paralog <- get_plot("paralog_comparison.R", "p")     # e.g. the final ggplot
p_iso     <- get_plot("isoform_comparison.R", "p")
p_cplx    <- get_plot("complex_comparison_unmatched.R", "p")

# ---- 2. Strip titles / subtitles / captions from every panel -----------------
strip <- function(g) g +
  labs(title = NULL, subtitle = NULL, caption = NULL, tag = NULL) +
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank())

p_paralog <- strip(p_paralog)
p_iso     <- strip(p_iso)
p_cplx    <- strip(p_cplx)

# ---- 3. Merge A | B | C ------------------------------------------------------
merged <- (p_paralog | p_iso | p_cplx) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 20))

ggsave("merged_paralog_isoform_complex.pdf", merged, width = 18, height = 6, device = cairo_pdf)
ggsave("merged_paralog_isoform_complex.png", merged, width = 18, height = 6, dpi = 200)
cat("Wrote merged_paralog_isoform_complex.pdf / .png\n")
