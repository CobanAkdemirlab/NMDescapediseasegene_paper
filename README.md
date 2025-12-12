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
  - Gene-level features (pLI, LOEUF, eenrichment analysis, tau, etc.)
  - Variant-level features (VEP annotation, CDS position of the PTC, variant distance to CDS end)
  - Protein-level features (IDRs, Pfam, AlphaFold2)
- FASTA and VCF generation from key 
- Modular script design for flexible expansion  

---

## Directory Structure

```text
NMDesc(make it match the folder structure)/
│
├── README.md
├── Clinvar_step1_NMD.R
├── Clinvar_step2_NMD.R
│
├── data(add to git)/
│   ├── clinvar_raw/
│   ├── gnomAD_variants/
│   └── external_annotations/
│
├── results/
│   ├── variant_sets/
│   ├── fasta/
│   ├── vcf/
│   ├── idr/
│   ├── af2/
│   └── vep/
│
└── scripts/
├── idr/
│   ├── predict_disorder.py
│   ├── predict_disorder_can.py
│   ├── get_snv_control_idr.R
│   ├── idr_merge.R
│   ├── idr_output.R
│   ├── quality_filter_idr.R
│   ├── README.md
│
├── af2/
│   ├── AF2/
│   ├── AF2_draw.R
│   ├── af2_feature_extraction.R
│   ├── nmd_alphafold_aa_composition.ipynb
│   ├── nmd_alphafold_predicted_properties.ipynb
│   ├── Density/
│   ├── attributes/
│   └── README.md
│
├── vep/
│   ├── vep_processing.R
│   ├── VEP/
│   ├── clean_variant_AD.R
│   ├── clinar_step1_NMD.R
│   └── README.md
│
├── frameshift/
│   ├── frameshift/
│   ├── frameshift_code.R
│   ├── create_fs_control.R
│   ├── get_fs_control.R
│   ├── get_fs_gnomAD_control.R
│   ├── get_fs_gnomAD_control_step1.R
│   ├── get_fs_gnomAD_control_step2.R
│   ├── get_fs_control_variant_new.R
│   ├── get_fs_variant_new.R
│   └── README.md
│
├── gnomad/
│   ├── gnomAD/
│   ├── get_gnomAD_control.R
│   ├── get_gnomAD_control_step2.R
│   ├── get_snv_gnomAD_control_step2.R
│   ├── get_snv_control_variant_new.R
│   ├── get_snv_variant_new.R
│   ├── get_variants.R
│   └── README.md
│
├── NMD/
│   ├── NMDesc
│   ├── WT_var_NMD_4c_AD.csv
│   ├── ClinVar_NMD.R
│   ├── ClinVar_step2_NMD.R
│   ├── NMD_annotate.R
│   ├── get_NMD_enrich2
│   ├── get_NMD_enrichment.R
│   ├── get_statistics.R
│   ├── get_pvalue.R
│   └── README.md
│
├── variant_processing/
│   ├── gene_level.R
│   ├── variant_level.R
│   ├── transcript_level.R
│   ├── select_AD.R
│   ├── select_AD_variants.R
│   ├── clean_variant_AD.R
│   ├── new_create_fasta.R
│   ├── create_fasta_new/
│   ├── compare.R
│   ├── bind_result.R
│   ├── cross_check.R
│   └── README.md
│
├── plotting/
│   ├── plots/
│   ├── compare_length/
│   ├── overlap/
│   ├── with_in/
│   └── plotting_functions.R
│
└── main.R

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

### Step 1 — Annotate ClinVar PTC variants

```{r step1, eval=FALSE}
source("Clinvar_step1_NMD.R")
```

### Step 2 — Apply canonical NMD rules (EJC model)

```{r step2, eval=FALSE}
source("Clinvar_step2_NMD.R")
```

These two scripts generate all core variant objects used throughout the NMDesc pipeline.

---

## Variant Objects & Usage

### Example: `plus1_variants`

| Output | Generated From   | Used For |
|--------|------------------|----------|
| FASTA  | `plus1_variants` | IDR analysis, AlphaFold2 inputs |
| VCF    | `plus1_variants` | VEP functional annotation |

#### Example: FASTA generation

```{r fasta-example, eval=FALSE}
# library(stringr)
 minus1_dis = create_fasta(minus1_variants, output_dir = "minus1_test_fasta_output")
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



---

## Output Summary

| Folder | Description |
|--------|-------------|
| `variant_sets/` | Final NMD classification results |
| `fasta/` | FASTA files for protein-based analyses |
| `vcf/` | VCF files for VEP input |
| `idr/` | Intrinsic disorder predictions & plots |
| `af2/` | AlphaFold2 structural features |
| `vep/` | VEP annotations and processed tables |

---


## Contact

**Maintainer:** Jiaoyang Xu (JXU)  
Email: [jiaoyang.xu@uth.tmc.edu]
