# NMDesc: NMD Escape Annotation & Feature Extraction Pipeline

![badge](https://img.shields.io/badge/status-active-brightgreen)
![badge](https://img.shields.io/badge/R->=4.2-blue)
![badge](https://img.shields.io/badge/data-ClinVar-orange)
![badge](https://img.shields.io/badge/purpose-NMD%20annotation-purple)

The **NMDesc pipeline** annotates and analyzes **premature termination codon (PTC) variants from ClinVar and gnomAD**, classifying them by whether they **escape Nonsense-Mediated Decay (NMD)** under canonical **Exon Junction Complex (EJC) rules**.  
It additionally extracts **gene-, variant-, and protein-level features** from multiple genomic and structural databases.

---

## Features

- Canonical NMD escape determination using EJC rules  
- Automated extraction of:
  - Gene-level features (pLI, LOEUF, enrichment analysis, tau, etc.)
  - Variant-level features (VEP annotation, CDS position of the PTC, variant distance to CDS end)
  - Protein-level features (IDRs, Pfam, AlphaFold2)
- FASTA and VCF generation from key 
- Modular script design for flexible expansion  

---

## Directory Structure

This project includes gene level(NMDesc disease genes and control disease genes), variant level(NMDesc variants from clinvar and gnomad) and protein level analysis.

```text
NMDescapediseasegene_paper-main/
│
├── main.R
│
├── scripts/
│   │
│   ├── preprocessing/
│   │   ├── preprocess_clinvar.R
│   │   ├── preprocess_gnomad.R
│   │   ├── filter_ptc_variants.R
│   │   ├── build_txdb.R
│   │   ├── generate_nmdesc_regions.R
│   │   ├── annotate_transcripts.R
│   │   ├── extract_canonical_transcripts.R
│   │   └── prepare_synonymous_controls.R
│   │
│   ├── variant_level/
│   │   ├── calculate_NMD_escape.R
│   │   ├── calculate_frameshift_PTCs.R
│   │   ├── calculate_ptc_distance.R
│   │   ├── calculate_nmdesc_region_length.R
│   │   ├── compute_gc_content.R
│   │   ├── compute_repeat_content.R
│   │   ├── motif_overlap_analysis.R
│   │   ├── LCS_analysis.R
│   │   ├── transcript_matched_analysis.R
│   │   ├── repeated_sampling_1000x.R
│   │   ├── paired_variant_comparison.R
│   │   ├── variant_feature_matrix.R
│   │   ├── synonymous_normalization.R
│   │   ├── plus1_plus2_comparison.R
│   │   ├── variant_QC.R
│   │   └── variant_filtering_pipeline.R
│   │
│   ├── gene_level/
│   │   ├── disease_gene_enrichment.R
│   │   ├── calculate_pLI_LOEUF.R
│   │   ├── transcript_matching.R
│   │   ├── tau_expression_analysis.R
│   │   ├── gene_feature_matrix.R
│   │   ├── gene_level_statistics.R
│   │   ├── gene_motif_analysis.R
│   │   ├── gene_LCS_analysis.R
│   │   ├── gene_gc_content_analysis.R
│   │   ├── OMIM_AD_filtering.R
│   │   ├── ClinVar_gene_summary.R
│   │   ├── gnomAD_gene_summary.R
│   │   └── gene_QC.R
│   │
│   ├── protein_level/
│   │   ├── PFAM_overlap_analysis.R
│   │   ├── PFAM_distance_analysis.R
│   │   ├── PPI_overlap_analysis.R
│   │   ├── STRINGdb_degree_centrality.R
│   │   ├── AlphaFold_feature_extraction.R
│   │   ├── SASA_analysis.R
│   │   ├── IDR_analysis.R
│   │   ├── phase_separation_analysis.R
│   │   ├── PICNIC_score_analysis.R
│   │   ├── hydrophobic_cluster_analysis.R
│   │   ├── sticker_feature_analysis.R
│   │   ├── prion_like_domain_analysis.R
│   │   ├── protein_charge_analysis.R
│   │   ├── amino_acid_composition_analysis.R
│   │   ├── interface_residue_overlap.R
│   │   ├── protein_structure_mapping.R
│   │   └── protein_QC.R
│   │
│   ├── statistics/
│   │   ├── paired_wilcoxon_tests.R
│   │   ├── mcnemar_tests.R
│   │   ├── exact_binomial_tests.R
│   │   ├── Wald_logOR_analysis.R
│   │   ├── FDR_correction.R
│   │   ├── mixed_effect_models.R
│   │   ├── bootstrap_analysis.R
│   │   ├── permutation_tests.R
│   │   ├── regression_models.R
│   │   ├── feature_correlation_analysis.R
│   │   ├── enrichment_statistics.R
│   │   ├── sensitivity_analysis.R
│   │   └── model_comparison_analysis.R
│   │
│   ├── plotting/
│   │   ├── plot_dist_to_cds_end.R
│   │   ├── plot_nmdesc_region_length.R
│   │   ├── plot_cds_length.R
│   │   ├── plot_gc_content.R
│   │   ├── plot_repeat_content.R
│   │   ├── plot_pfams.R
│   │   ├── plot_ppi_overlap.R
│   │   ├── plot_tau_violin.R
│   │   ├── plot_LOEUF_pLI.R
│   │   ├── plot_phase_separation_features.R
│   │   ├── plot_resampling_results.R
│   │   ├── plot_matched_analysis.R
│   │   ├── plot_foldchange_histograms.R
│   │   ├── plot_centrality_results.R
│   │   ├── plot_variant_density.R
│   │   ├── plot_protein_features.R
│   │   ├── plot_feature_heatmaps.R
│   │   └── generate_manuscript_figures.R
│   │
│   │
│   └── QC/
│       ├── clinvar_QC.R
│       ├── gnomad_QC.R
│       ├── transcript_QC.R
│       ├── FASTA_QC.R
│       ├── PFAM_QC.R
│       ├── PPI_QC.R
│       ├── matching_QC.R
│       └── statistical_QC.R
│
└── 
```

---

## Installation

### 1. Install R (≥ 4.2)

Download R from: <https://www.r-project.org/>

### 2. Install required R packages

```{r install-packages, eval=FALSE}
install.packages(c(
  "tidyverse", "data.table", "biomaRt",
  "stringr", "jsonlite", "readr", "ggplot2",
  "scales", "ggpubr" 
))
```

### 3. Optional external tools

| Tool | Purpose |
|------|---------|
| **VEP (Variant Effect Predictor)** | Variant functional annotation |
| **AlphaFold2 models** | Protein structural feature extraction |
| **MetaPredict** | Intrinsic disorder prediction |

---

## Quick Start

```{r step1, eval=FALSE}
source("main.R")
```

This script generate all core variant objects used throughout the NMDesc pipeline.

---

## Variant Objects & Usage

### Example: `snv_variants`

| Output | Generated From   | Used For |
|--------|------------------|----------|
| FASTA  | `snv_variants` | IDR analysis, AlphaFold2 inputs |
| VCF    | `snv_variants` | VEP functional annotation |

#### Example: FASTA generation

```{r fasta-example, eval=FALSE}
# library(stringr)
 snv_dis = create_fasta(snv_variants, output_dir = "snv_test_fasta_output")
```

#### Example: VCF generation

```{r vcf-example, eval=FALSE}
 vcf_df <- snv_control_variants %>%
  extract(key, into = c("CHROM", "POS", "REF", "ALT"), regex = "([^:]+):([0-9]+)\\|#([^|]+)\\|(.+)", remove = FALSE) %>%
  mutate(
    ID = ".",
    QUAL = ".",
    FILTER = "PASS",
    INFO = paste0("TRANSCRIPT=", transcript, ";UNIPROT=", uniprotswissprot)
  ) %>%
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
```

---

## Workflow Diagram

```text
ClinVar pipeline
────────────────

ClinVar
  │
  ├─ Select germline variants
  │      ResetID_clinvar_Clnsig.txt
  │
  ├─ NMD annotation
  │      Clinvar_1120.rds
  │
  ├─ Select snv/fs, plp/vus, ptc, nmdesc
  │      Snv_plp_ptc_res1120.rds
  │
  ├─ Add p value
  │      Snv_plp_ptc_p1122.rds
  │
  └─ Get_NMD_enrichment (modify output txt name)
         plus1_can_gene0217.txt


R helper scripts / metadata
───────────────────────────

Gene-level
  ├─ BM.info4(cds, exon_chrom, transcript_id, rank)
  └─ Snv_tx (canonical transcript names)

Variant-level
  ├─ Get_snv_variant_new.R
  │      (remove repeated steps for creating res file)
  └─ Snv_variants(snv_variants0406.csv,
                  includes uniprot id, transcript, key)

Key mapping
  └─ Snv_key_to_transcript   (key is not unique)


FASTA / VCF branches
────────────────────

From Snv_variants:

  ├─ Create_fasta
  │      → FASTA files (e.g. Minus1_dis)
  │           → IDR analysis / AF2 analysis
  │
  └─ csv2vep
         → Snv.vcf
              (includes key, transcript & uniprot;
               variants identified by key)
         → VEP
              Snv_NMD_result3_vep.txt

```

---

## Downstream Analyses

### 1. IDR Prediction
- MetaPredict


### 2. AlphaFold2 Structural Feature Extraction
- pLDDT
- secondary structure
- SASA

### 3. VEP Functional Annotation
- Consequence terms
- Nearest exon junction boundary
- dbNSFP features: Condel_score, GERP scores etc

### 4. Protein domain analysis
- PPI
- PFAM
- SLM, NLS



---

## Output Summary

| Folder | Description |
|--------|-------------|
| `gene_results/` | gene level features |
| `variant_results/` | variant level features |
| `fasta/` | FASTA files for protein-based analyses |
| `vcf/` | VCF files for VEP input |
| `idr/` | Intrinsic disorder predictions & plots |
| `af2/` | AlphaFold2 structural features |
| `vep/` | VEP annotations and processed tables |

---


## Contact

**Maintainer:** Jiaoyang Xu (JXU)  
Email: [jiaoyang.xu@uth.tmc.edu]
