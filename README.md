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
NMDescapediseasegene_paper/
├── gene level_v4/
│   ├── QC/
│   │   ├── check.R
│   │   ├── combine_gene.R
│   │   ├── compare.R
│   │   ├── compare_submitters.R
│   │   ├── cross_check.R
│   │   ├── download_gnomad_syn_vcf.R
│   │   ├── download_syn2.R
│   │   ├── negative_control.R
│   │   ├── parse_gnomad_syn_vcf.R
│   │   ├── process_syn2.R
│   │   ├── select_AD.R
│   │   └── select_AD_variants.R
│   ├── control genes/
│   │   ├── frameshift/
│   │   │   └── get_fs_control_gene.R
│   │   └── snv/
│   │       ├── get_gnomAD_control.R
│   │       ├── get_gnomAD_control_step2.R
│   │       ├── get_snv_control_gene.R
│   │       └── get_snv_gnomAD_control_step2.R
│   ├── disease genes/
│   │   ├── framesift/
│   │   │   ├── bind_result_dbh.R
│   │   │   ├── compare+1-1.R
│   │   │   ├── frameshift_code.R
│   │   │   ├── fs_transcript_level.R
│   │   │   ├── process_syn.R
│   │   │   └── variant_level.R
│   │   └── snv/
│   │       ├── main.R        <- entry point
│   │       ├── clinar_step1_NMD.R
│   │       ├── ClinVar_NMD.R
│   │       ├── ClinVar_step2_NMD.R
│   │       ├── extract_enriched.R
│   │       ├── get_NMD_enrich2.R
│   │       ├── get_NMD_enrichment_DBH.R
│   │       ├── get_pvalue.R
│   │       ├── NMD_annotate.R
│   │       └── process_syn.R
│   ├── features/
│   │   ├── functions/
│   │   │   ├── annotate_motif_flags.R
│   │   │   ├── build_gene_all.R
│   │   │   ├── calculate_ppi_degree_centrality.R
│   │   │   ├── plot_gc_content.R
│   │   │   ├── plot_gene_level_features.R
│   │   │   ├── plot_repeat_content.R
│   │   │   ├── run_pfam_overlap_analysis.R
│   │   │   ├── run_ppi_overlap_analysis.R
│   │   │   └── run_tau_analysis.R
│   │   ├── big_pIc.R
│   │   ├── gene_connectivity.R
│   │   ├── gene_gc.R
│   │   ├── gene_motif.R
│   │   ├── gene_NmdescRegion.R
│   │   ├── gene_overview.R
│   │   ├── gene_pfam.R
│   │   ├── gene_PliLoeuf.R
│   │   ├── gene_plot.R
│   │   ├── gene_plot_match.R
│   │   ├── gene_plot_match2.R
│   │   ├── gene_ppi.R
│   │   ├── gene_repeat.R
│   │   ├── gene_tau.R
│   │   ├── GO_enrich.R
│   │   ├── inheritance.R
│   │   ├── new_loc.R
│   │   └── plot_supplemental_fig3_gene_level.R
│   ├── lib/
│   │   ├── get_statistics.R
│   │   └── paths.R
│   └── gene_main_dbh.R        <- entry point
├── protein level_v4/
│   ├── AF2/
│   │   └── AF2_draw.R
│   └── fasta/
│       ├── create_fs_control.R
│       └── new_create_fasta.R
├── variant level_v4/
│   ├── QC/
│   │   ├── adjust_variant.R
│   │   ├── check_pogz.R
│   │   ├── check_variant_dis.R
│   │   ├── clean_variant_AD.R
│   │   ├── combine_variant.R
│   │   ├── compare_D.R
│   │   ├── plot_negative_control.R
│   │   ├── regression.R
│   │   └── resample.R
│   ├── clinvar/
│   │   ├── frameshift/
│   │   │   └── get_fs_variant_new.R
│   │   └── snv/
│   │       ├── get_snv_control_variant_new.R
│   │       └── get_snv_variant_new.R
│   ├── emel/
│   │   ├── sequence/
│   │   │   ├── sequence-parsev2_main.R
│   │   │   ├── sequence-peptides_main.R
│   │   │   ├── sequence-peptides_read.R
│   │   │   ├── sequence_cider_main.R
│   │   │   ├── sequence_cider_match.R
│   │   │   ├── sequence_cider_read.R
│   │   │   ├── sequence_cider_table1.R
│   │   │   ├── sequence_parsev2_hier.R
│   │   │   ├── sequence_parsev2_match.R
│   │   │   ├── sequence_parsev2_read.R
│   │   │   ├── sequence_parsev2_table1.R
│   │   │   ├── sequence_peptide_hier.R
│   │   │   ├── sequence_peptide_match.R
│   │   │   └── sequence_peptides_table1.R
│   │   └── structure/
│   │       ├── structure_hier.R
│   │       ├── structure_match.R
│   │       ├── structure_mixed.R
│   │       ├── structure_PAE.R
│   │       ├── structure_plddt.R
│   │       ├── structure_read.R
│   │       ├── structure_table1.R
│   │       └── structure_table2.R
│   ├── features/
│   │   ├── IDR/
│   │   │   └── idr/
│   │   │       ├── compare_length/
│   │   │       │   ├── compare_idr_length.R
│   │   │       │   ├── compare_idr_snv_control.R
│   │   │       │   ├── idr_difff.R
│   │   │       │   ├── idr_plot.R
│   │   │       │   └── wildtype_idr.R
│   │   │       ├── overlap/
│   │   │       │   ├── get_fs_idr_match.R
│   │   │       │   ├── get_snv_idr_match.R
│   │   │       │   ├── get_snv_idr_match2.R
│   │   │       │   └── idr_main.R
│   │   │       ├── with_in/
│   │   │       │   └── idr_match2.R
│   │   │       ├── get_snv_control_idr.R
│   │   │       ├── idr_merge.R
│   │   │       ├── idr_output.R
│   │   │       └── quality_filter_idr.R
│   │   ├── functions/
│   │   │   ├── variant_pfam_ppi.R
│   │   │   └── variants_motif.R
│   │   ├── vep/
│   │   │   ├── process_vep.R
│   │   │   └── vep_draw.R
│   │   ├── clean_motif.R
│   │   ├── create_fasta.R
│   │   ├── csv2vcf.R
│   │   ├── do_AD.R
│   │   ├── do_motif.R
│   │   ├── het_6fasta.R
│   │   ├── merge_loc_uni.R
│   │   ├── number_submitters.R
│   │   ├── parse.R
│   │   ├── pfam_ppi_variant.R
│   │   ├── plus1_control.R
│   │   ├── transcript_features.R
│   │   ├── variant_motif.R
│   │   ├── variant_plot.R
│   │   └── variant_plot2.R
│   ├── gnomad/
│   │   ├── frameshift/
│   │   │   └── get_fs_control_variant_new.R
│   │   ├── snv/
│   │   │   └── get_gnomAD_control.R
│   │   └── gnomAD_downloaddata.R
│   ├── variant_main_DBH.R        <- entry point
│   └── new_create_fasta_functions.R
└── run_analysis.R        <- entry point
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

#### Example: Calculate PPI Degree centrality

```{r vcf-example, eval=FALSE}
 calculate_ppi_degree_centrality(
  gene_all,
  output_csv = "wald_ppi_degree_centrality_results.csv"
)
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

### 1. Protein domain analysis
- PPI
- PFAM
- SLM, NLS

### 2. Potential confounders
- CDS length
- PTC distance to CDS end
- GC content

### 3. AlphaFold2 Structural Feature Extraction
- pLDDT
- secondary structure
- SASA

### 4. VEP Functional Annotation
- Consequence terms
- Nearest exon junction boundary
- dbNSFP features: Condel_score, GERP scores etc





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


## Execution Order

The three analysis levels are run independently. Each has its own entry point,
and all of them expect `paths.R` to resolve the inputs.

```r
# 1. Gene level
source("gene level_v4/disease genes/snv/main.R")   # ClinVar SNV enrichment
source("gene level_v4/gene_main_dbh.R")            # gene-level feature assembly

# 2. Variant level
source("variant level_v4/variant_main_DBH.R")      # variant feature matrix + models

# 3. Protein level
source("protein level_v4/fasta/new_create_fasta.R")  # FASTA for IDR / AF2 input
source("protein level_v4/AF2/AF2_draw.R")            # AF2 structural figures
```


```text
ClinVar / gnomAD download
        |
        v
aenmd NMD-escape annotation      (ClinVar_NMD.R, NMD_annotate.R)
        |
        v
PTC filtering, canonical transcripts   (clinar_step1_NMD.R, ClinVar_step2_NMD.R)
        |
        +--> gene-level enrichment      (get_NMD_enrich2.R, get_pvalue.R)
        +--> variant feature matrix     (transcript_features.R, variant_motif.R)
        +--> control sets               (control genes/, gnomad/)
        |
        v
feature extraction  (pfam, ppi, motif, gc, repeat, tau, IDR, VEP, AF2)
        |
        v
statistics and figures   (features/*_plot*.R, emel/, QC/)
```

External steps are run outside R and their output read back in:

| Step | Tool | Consumed by |
|------|------|-------------|
| Variant annotation | VEP | `variant level_v4/features/vep/process_vep.R` |
| Disorder prediction | metapredict | `variant level_v4/features/IDR/idr/` |
| Structure prediction | AlphaFold2 | `protein level_v4/AF2/AF2_draw.R` |
| Sequence properties | CIDER / PARSE | `variant level_v4/emel/sequence/` |

---

## Statistical Methods in the Code

Tests and models used, with the number of scripts referencing each.

| Method | Scripts | Typical use |
|--------|---------|-------------|
| Wilcoxon | 37 | Paired and unpaired feature comparisons |
| BH / FDR | 24 | Multiple-testing correction across feature panels |
| Fisher exact | 13 | Gene- and variant-level enrichment tables |
| binomial exact | 7 | Observed vs expected NMD-escape counts |
| GLMM (lme4) | 7 | Mixed-effect models with gene or protein random intercepts |
| Bayesian GLMM (brms) | 4 | Sensitivity analysis for degenerate fits |
| McNemar | 3 | Paired categorical flag comparisons |
| chi-squared | 3 | Categorical contingency tests |
| t-test | 3 | Continuous feature comparisons |
| logistic regression | 3 | Binary outcome modelling on variant features |
| STRINGdb centrality | 3 | PPI degree centrality per gene |
| DiscreteFDR | 2 | FDR correction for discrete test statistics |
| Spearman | 2 | Rank correlation between p-value methods |
| ACAT | 1 | Aggregated Cauchy combination of per-transcript p-values |
| Tarone | 1 | Filtering structurally unreachable discrete tests |
| Kruskal-Wallis | 1 | Multi-group comparisons |
| bootstrap/resampling | 1 | Repeated sampling of one variant per protein |

---

## Contact

**Maintainer:** Jiaoyang Xu (JXU)  
Email: [jiaoyang.xu@uth.tmc.edu]

---
