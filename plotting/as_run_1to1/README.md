# Gene-level manuscript figures, as run

Gene set 4 (top 20% rare synonymous: 94 frameshift, 64 stopgain) at 1:1 matching.

Run order: `01_...` first, which writes the three canonical tables every later step reads,
then the per-figure scripts in any order. Each script reads its inputs from the working
directory. `03_fig3_combine_AK.R` and `10_figS3_merge_ABC.R` take `SCRIPT_DIR` pointing at
the folder holding their sub-scripts; the latter also reads the three
`10_figS3_*_comparison.R` from the working directory.

| Figure | Script | Inputs |
|---|---|---|
| all gene-level | `01_match_1to1_and_build_tables.py` | `union_all_gene_withfeatures_fixed.csv`, `gene_all_AD_top20.csv`, the case list |
| 2, frameshift | `02_fig2_GO_frameshift_always.R` | `Figure2_GO/` |
| 2, stopgain | `02_fig2_GO_stopgain.R` | `Figure2_GO/` |
| 3 A-E | `03_fig3_panels_AE_matching.R` | `Figure3/` |
| 3 F-G | `03_fig3_build_string_tau_cache.py` then `03_fig3_panels_FG_ppi_tau.R` | `Figure3/` plus the STRING and GTEx files below |
| 3 H-K | `03_fig3_panels_HK_pli_loeuf.R` | `Figure3/pli_loeuf_genes.csv` |
| 3 A-K composite | `03_fig3_combine_AK.R` | re-runs the three above |
| 4, masked and unmasked | `04_fig4_gene_flags.py` | `Figure4/` |
| 5 and 5 position-adjusted | `05_fig5_variant_flags_lme.py` | `Figure5/` |
| 8 | `08_fig8_ndd_enrichment.R` | `Figure8/` |
| S1 | `06_figS1_go_robustness_permutation.R` | `FigureS1_GOperm/` |
| S2 and S9 | `07_figS2_S9_repeat_homopolymer.R` | `FigureS2_S9/` |
| S3 A-C | `10_figS3_fetch_paralog_isoform_counts.py`, the three `10_figS3_*_comparison.R`, then `10_figS3_merge_ABC.R` | `FigureS3/` |
| S5 | `11_figS5_build_rows_and_lme.py` then `11_figS5_plot.py` | `FigureS5/` |

Figures 6, 6E, 7, S6 and S7 are variant-level and carry no matching ratio; their scripts are
in `../variant_level/`.

## Inputs held outside this folder by size

- STRING v11.5 and GTEx tau, read by `03_fig3_build_string_tau_cache.py`:
  `9606.protein.info.v11.5.txt.gz`, `9606.protein.aliases.v11.5.txt.gz`,
  `9606.protein.links.v11.5.txt.gz`, `Tau_gene_V8_derived.csv`. The script writes a cache
  that `03_fig3_panels_FG_ppi_tau.R` reads.
- `gnomad.v4.0.constraint_metrics.tsv` (86 MB) is the upstream source of
  `Figure3/pli_loeuf_genes.csv`; the H-K script reads the cache.

## How the scripts behave

- `02_fig2_GO_frameshift_always.R` writes its terms table and draws a figure whether or not
  any term reaches FDR < 5%. On this gene set the frameshift stratum reaches none, and the
  table is written with a header only.
- `03_fig3_build_string_tau_cache.py` builds the STRING degree and GTEx tau cache that
  panels F-G read; `degree_centrality` in the gene table is NA when the STRING download is
  unavailable.
- Figure S3 panels A and B drop a pair when either gene has no count. One control gene
  (TCF7L2) has no paralogue count, so panel A rests on 93 / 64 pairs against 94 / 64 in
  panels B and C. Panel C uses an unpaired rank-sum on binary CORUM membership.
- `11_figS5_build_rows_and_lme.py` reads the IDR catalogue built over all four arms, so a
  protein with no Region >= 20 aa is a true zero.
- Gene-level Pfam follows `max_pfam_end * 3 >= NMD_region_start` with NA preserved,
  computed in `01_match_1to1_and_build_tables.py`.
