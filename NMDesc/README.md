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

### 1. Install R (≥4.2)
Download: https://www.r-project.org/

### 2. Install required R packages
```r
install.packages(c(
  "tidyverse", "data.table", "biomaRt",
  "stringr", "jsonlite", "readr", "ggplot2"
))

