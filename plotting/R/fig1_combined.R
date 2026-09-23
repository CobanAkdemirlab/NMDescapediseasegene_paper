#!/usr/bin/env Rscript
# =============================================================================
# fig1_combined.R  —  assemble the full Figure (panels A-K) from the three
# sub-figure scripts, so all panels share one canvas and one A..K tag sequence.
#
#   A B C  length-distribution violins (CDS, NMD region, NMD/CDS ratio)
#   D E    covariate-balance Love plots (Frameshift, Stopgain)
#   F G    PPI degree centrality, Tissue specificity (Tau)
#   H I    pLI category distribution (Frameshift, Stopgain)
#   J K    LOEUF category distribution (Frameshift, Stopgain)
#
# It re-runs the three validated scripts (saving disabled) and reuses their
# panel objects verbatim, so the merged figure stays identical to the pieces.
#
# Requires, in the working directory:
#   scripts : fig_matching_combined_fixed.R, updated_builder_ppitau.R,
#             updated_fig_pli_loeuf.R
#   data    : matched_pairs_full.csv, gene_all_AD.csv,
#             ppi_tau_4groups.csv (PPI/Tau cache),
#             pli_loeuf_genes.csv (pLI/LOEUF cache)
# Outputs : Figure1_combined.png / .pdf
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
  library(scales); library(readr); library(forcats); library(grid)
})

SCRIPT_DIR <- Sys.getenv("SCRIPT_DIR", ".")   # where the three sub-scripts live

# ---- source a script in an isolated env, with ggsave/cat/message silenced ----
run_env <- function(file) {
  e <- new.env(parent = globalenv())
  e$ggsave  <- function(...) invisible(NULL)          # don't write sub-figures
  e$message <- function(...) invisible(NULL)
  e$cat     <- function(...) invisible(NULL)
  sys.source(file.path(SCRIPT_DIR, file), envir = e)
  e
}

message("Building panels A-E (matching) ...")
mAE <- run_env("fig_matching_combined_fixed.R")   # pA pB pC pD pE
message("Building panels F-G (PPI/Tau) ...")
mFG <- run_env("updated_builder_ppitau.R")        # pA(->F) pB(->G)
message("Building panels H-K (pLI/LOEUF) ...")
mHK <- run_env("fig_pli_loeuf.R")         # pA pB pC pD (->H I J K)

# ---- collect + re-tag A..K (strip any tag the sub-scripts set) --------------
retag <- function(p, tag) p + labs(tag = tag) +
  theme(plot.tag = element_text(face = "bold", size = 24),
        plot.tag.position = c(0, 1))

A <- retag(mAE$pA,"A"); B <- retag(mAE$pB,"B"); C <- retag(mAE$pC,"C")
D <- retag(mAE$pD,"D"); E <- retag(mAE$pE,"E")
F_ <- retag(mFG$pA,"F"); G <- retag(mFG$pB,"G")
H <- retag(mHK$pA,"H"); I <- retag(mHK$pB,"I")
J <- retag(mHK$pC,"J"); K <- retag(mHK$pD,"K")

# keep just one shared legend (the constraint panels carry it); drop it elsewhere
A<-A+theme(legend.position="none"); B<-B+theme(legend.position="none"); C<-C+theme(legend.position="none")
D<-D+theme(legend.position="none"); E<-E+theme(legend.position="none")
F_<-F_+theme(legend.position="none"); G<-G+theme(legend.position="none")
H<-H+theme(legend.position="none"); I<-I+theme(legend.position="none")
J<-J+theme(legend.position="none")   # K keeps the legend -> collected to bottom

# ---- section headers --------------------------------------------------------
hdr <- function(txt) wrap_elements(full = textGrob(
  txt, x = 0.5, hjust = 0.5, gp = gpar(fontface = "bold", fontsize = 21)))

final <-
  hdr("Length-matched controls reproduce the size distribution of PTC variants") /
  (A | B | C) /
  (D | E) /
  hdr("Protein-network and tissue-expression properties") /
  (F_ | G) /
  hdr("Constraint (pLI / LOEUF) category distributions") /
  (H | I) /
  (J | K) +
  plot_layout(heights = c(0.10, 1, 0.85, 0.10, 1, 0.10, 1, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.text = element_text(size = 13))

ggsave("Figure1_combined.png", final, width = 16, height = 23, dpi = 300,
       bg = "white", limitsize = FALSE)
ggsave("Figure1_combined.pdf", final, width = 16, height = 23,
       bg = "white", limitsize = FALSE)
cat("saved Figure1_combined.png / .pdf\n")
