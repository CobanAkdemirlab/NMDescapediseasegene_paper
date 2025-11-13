# NMDesc: NMD Escape Annotation & Feature Extraction Pipeline

![badge](https://img.shields.io/badge/status-active-brightgreen)
![badge](https://img.shields.io/badge/R->=4.2-blue)
![badge](https://img.shields.io/badge/data-ClinVar-orange)
![badge](https://img.shields.io/badge/purpose-NMD%20annotation-purple)

The **NMDesc pipeline** annotates and analyzes **premature termination codon (PTC) variants from ClinVar**, classifying them by whether they **escape Nonsense-Mediated Decay (NMD)** under canonical **Exon Junction Complex (EJC) rules**.  
It additionally extracts **gene-, variant-, and protein-level features** from multiple genomic and structural databases.

---

## 📌 Features

- Canonical NMD escape determination using EJC rules  
- Automated extraction of:
  - Gene-level features (constraint metrics, LOEUF, etc.)
  - Variant-level features (VEP, CDS position, NMD rules)
  - Protein-level features (IDRs, Pfam, AlphaFold2)
- FASTA and VCF generation  
- Modular script design for flexible expansion  

---

## 📁 Directory Structure

```text
NMDesc/
│
├── README.md
├── Clinvar_step1_NMD.R
├── Clinvar_step2_NMD.R
│
├── data/
│   ├── clinvar_raw/
│   ├── canonical_transcripts/
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
    ├── idr_analysis.R
    ├── af2_feature_extraction.R
    ├── vep_processing.R
    └── plotting/
```

---

## 🚀 Installation

### 1. Install R (≥ 4.2)

Download R from: <https://www.r-project.org/>

### 2. Install required R packages

```{r install-packages, eval=FALSE}
install.packages(c(
  "tidyverse", "data.table", "biomaRt",
  "stringr", "jsonlite", "readr", "ggplot2"
))
```

### 3. Optional external tools

| Tool | Purpose |
|------|---------|
| **VEP (Variant Effect Predictor)** | Variant functional annotation |
| **AlphaFold2 models** | Protein structural feature extraction |
| **MetaPredict / IUPred2A** | Intrinsic disorder prediction |

---

## 🔥 Quick Start

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

## 📌 Variant Objects & Usage

### `plus1_variants`

| Output | Generated From   | Used For |
|--------|------------------|----------|
| FASTA  | `plus1_variants` | IDR analysis, AlphaFold2 inputs |
| VCF    | `plus1_variants` | VEP functional annotation |

#### Example: FASTA generation

```{r fasta-example, eval=FALSE}
# library(seqinr)
# write.fasta(sequences = prot_seqs,
#             names     = prot_ids,
#             file.out  = "results/fasta/plus1.fasta")
```

#### Example: VCF generation

```{r vcf-example, eval=FALSE}
# write_vcf(plus1_variants, "results/vcf/plus1.vcf")
```

---

## 📊 Workflow Diagram

```text
ClinVar PTC Variants
        │
        ├── Step 1: Annotation (Clinvar_step1_NMD.R)
        │
        └── Step 2: Canonical NMD classification (Clinvar_step2_NMD.R)
                  │
         ┌────────┴─────────┐
         │                  │
   plus1_variants      other variant sets
         │
         ├── FASTA → IDR analysis / AF2 analysis
         └── VCF   → VEP annotation
```

---

## 🔬 Downstream Analyses

### 1. IDR Prediction
- MetaPredict
- IUPred2A

### 2. AlphaFold2 Structural Feature Extraction
- pLDDT
- secondary structure
- SASA

### 3. VEP Functional Annotation
- Consequence terms
- Nearest exon junction boundary
- Impact categories

---

## 📦 Output Summary

| Folder | Description |
|--------|-------------|
| `variant_sets/` | Final NMD classification results |
| `fasta/` | FASTA files for protein-based analyses |
| `vcf/` | VCF files for VEP input |
| `idr/` | Intrinsic disorder predictions & plots |
| `af2/` | AlphaFold2 structural features |
| `vep/` | VEP annotations and processed tables |

---


## 📫 Contact

**Maintainer:** Jiaoyang (JXU)  
Email: [jiaoyang.xu@uth.tmc.edu]
