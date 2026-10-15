# Manuscript figures

Every figure in the manuscript, with the script that drew it.

All computed figures use gene set 4 (`4_raresyn_top20`, top 20% rare synonymous: 94
frameshift and 64 stopgain case genes) at 1:1 CDS-length matching. Gene-level panels state
94 / 64 pairs; variant-level panels carry no matching ratio.

`MANIFEST.csv` lists each file with its byte size and a sha256 prefix.

## Figure to script

| figure | what | script |
|---|---|---|
| 1 | pipeline schematic (panels A, B) | drawn, not computed |
| 2 | GO:BP enrichment above the AD background, fold-ranked | `../../plotting/as_run_1to1/02_fig2_GO_frameshift_always.R`, `02_fig2_GO_stopgain.R` |
| 3 | length-matched control design and gene-level properties | `../../plotting/as_run_1to1/03_fig3_panels_AE_matching.R`, `03_fig3_panels_FG_ppi_tau.R` (with `03_fig3_build_string_tau_cache.py`), `03_fig3_panels_HK_pli_loeuf.R`, `03_fig3_combine_AK.R` |
| 4 | gene-level flag effects | `../../plotting/as_run_1to1/04_fig4_gene_flags.py` |
| 5 | variant-level flag effects, random intercept on gene | `../../plotting/as_run_1to1/05_fig5_variant_flags_lme.py` |
| 6 | post-frameshift tail composition (A-D) | `../../plotting/variant_level/composition_lmm_run.py` then `plot_composition_run.py` |
| 7 | whole-protein composition, stopgain | the same two scripts as Figure 6, one run writes both |
| 8 | NDD risk-gene enrichment | `../../plotting/as_run_1to1/08_fig8_ndd_enrichment.R` |
| 9 | RARB condensate microscopy | experimental |
| S1 | GO permutation robustness | `../../plotting/as_run_1to1/06_figS1_go_robustness_permutation.R` |
| S2 | repeat, homopolymer and GC content, paired | `../../plotting/as_run_1to1/07_figS2_S9_repeat_homopolymer.R` |
| S3 | paralog, isoform and protein-complex counts (A-C) | `../../plotting/as_run_1to1/10_figS3_fetch_paralog_isoform_counts.py`, the three `10_figS3_*_comparison.R`, then `10_figS3_merge_ABC.R` |
| S5 | IDR overlap of the NMD-escape region | `../../plotting/as_run_1to1/11_figS5_build_rows_and_lme.py` then `11_figS5_plot.py` |
| S6 | per-gene frameshift ranking | `../../plotting/variant_level/plot_rank.py <comp_variants.csv> fs <PS_proteins.txt> <out>` |
| S7 | per-gene stopgain ranking | the same script with `sg` |

The manuscript numbers run 1-9 and S1, S2, S3, S5, S6, S7. The repository's S7, S8 and S9
are the manuscript's S6, S7 and S2.

## Figure 1 and Figure 9

Figure 1 is a drawn schematic, shipped as its two panels: `Figure1_panelA.png` (the
canonical 50-bp rule) and `Figure1_panelB.png` (the pipeline, gene-level enrichment at
BH FDR < 0.2). Figure 9 is experimental microscopy. Neither is produced by this pipeline.

## Parameters the figure scripts take

- `composition_lmm_run.py` takes `MIN_TAIL = 20`; Figures 6, 7, S6 and S7 use that value.
- `04_fig4_gene_flags.py` takes the escape coordinates from the design row,
  `NMD_region_start` and `NMD_region_end` in the gene table.
- `03_fig3_combine_AK.R` and `10_figS3_merge_ABC.R` take `SCRIPT_DIR` pointing at the
  folder holding their sub-scripts.
- Gene-level Pfam follows `max_pfam_end * 3 >= NMD_region_start` with NA preserved,
  computed in `01_match_1to1_and_build_tables.py`.
