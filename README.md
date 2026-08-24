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

The layout below is the current v4 tree. Entry points are marked; `backup/`
holds earlier versions and is omitted here.

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


## Contact

**Maintainer:** Jiaoyang Xu (JXU)  
Email: [jiaoyang.xu@uth.tmc.edu]

---

## Configuration and Data Layout

All data paths resolve through `gene level_v4/lib/paths.R`. Scripts source it
near the top and then refer to files by name rather than by absolute path:

```r
source("gene level_v4/lib/paths.R")

d <- data_file("genemap2.txt")        # locate an input by filename
x <- read.csv(data_file("BM_info.csv"))
out <- out_file("results.csv")        # path under the output directory
root <- data_root("clinvar")          # a whole input directory
```

`paths.R` indexes these roots, in priority order:

| Root name | Default location |
|-----------|------------------|
| `clinvar` | `~/Desktop/clinvar` |
| `new_clinvar` | `~/Desktop/new_clinvar` |
| `gnomAD_snv` | `~/Desktop/gnomAD_snv` |
| `idr` | `~/Desktop/idr` |
| `syn` | `~/Desktop/syn` |
| `legacy` | `machine-specific archive directory` |

Two environment variables override the defaults:

| Variable | Effect |
|----------|--------|
| `NMDESC_DATA` | Prepends a directory to the search order, taking priority over all roots |
| `NMDESC_OUT` | Output directory; defaults to `~/Desktop/clinvar/output` |

```sh
export NMDESC_DATA=/path/to/nmdesc_data
export NMDESC_OUT=/path/to/nmdesc_output
```

The file index is cached as `.nmdesc_path_index.rds`. After adding or moving
input files, rebuild it:

```r
refresh_data_index()
check_data()          # report which expected inputs are present
```

---

## Package Dependencies

Packages referenced by two or more scripts in the v4 tree.

**CRAN**

```r
install.packages(c(
  "BiocManager", "DiscreteFDR", "STRINGdb", "bayesplot", "brms", "broom",
  "broom.mixed", "data.table", "dplyr", "enrichR", "flextable", "future",
  "genekitr", "ggplot2", "ggpubr", "glue", "gt", "gtsummary", "here",
  "igraph", "jsonlite", "lme4", "lmerTest", "patchwork", "posterior",
  "purrr", "readr", "readxl", "rstatix", "scales", "stringr", "tibble",
  "tidyr", "tidyselect", "tidyverse"
))
```

**Bioconductor**

```r
BiocManager::install(c(
  "AnnotationDbi", "BSgenome.Hsapiens.UCSC.hg38", "Biostrings",
  "GenomeInfoDb", "GenomicFeatures", "GenomicRanges", "IRanges",
  "S4Vectors", "SummarizedExperiment", "VariantAnnotation", "biomaRt",
  "org.Hs.eg.db", "txdbmaker"
))
```

`aenmd` supplies the NMD-escape rule annotations and is not on CRAN; the
transcript models come from a companion data package. Both install from GitHub:

```r
remotes::install_github("kostkalab/aenmd")
remotes::install_github("kostkalab/aenmd_data", subdir = "aenmd.data.ensdb.v105")
```

`aenmd.data.ensdb.v105` is the Ensembl v105 high-confidence transcript set
(protein-coding, transcript support level 1). A containerized command-line
interface is available at `kostkalab/aenmd_cli`.

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

Within a level the general dependency order is:

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

## Script Index (v4)

Every script in the v4 tree, grouped by analysis level, with a one-line
description of what it does. Entry points are marked **entry**.

### Gene level

**Shared library** — `gene level_v4/lib/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `get_statistics.R` | 142 | Computes summary statistics for gene-level NMD escape data |
| `paths.R` | 155 | Locates and resolves data file paths across the pipeline |

**Disease genes, SNV** — `gene level_v4/disease genes/snv/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `main.R` **entry** | 1066 | Runs full NMD/PPI enrichment pipeline across ClinVar, GTEx, and supplementary data |
| `ClinVar_NMD.R` | 154 | Downloads ClinVar, re-annotates with aenmd, keeps PTC/NMD-escape pathogenic variants |
| `ClinVar_step2_NMD.R` | 418 | Builds transcript/UTR/CDS sequences and merges ClinVar with OMIM data |
| `NMD_annotate.R` | 51 | Overlaps ClinVar variants with aenmd transcript exons for NMD annotation |
| `clinar_step1_NMD.R` | 80 | Annotates BCFtools-processed ClinVar variants with aenmd PTC calls |
| `extract_enriched.R` | 58 | Computes observed vs expected enrichment ratios and flags underpowered tests |
| `get_NMD_enrich2.R` | 117 | Calculates NMD escape enrichment statistics per gene |
| `get_NMD_enrichment_DBH.R` | 85 | Runs binomial and Fisher tests with BH correction for NMD enrichment |
| `get_pvalue.R` | 329 | Extracts p-values and exports significance tiers for OMIM AD genes |
| `process_syn.R` | 10 | Counts synonymous variants within specified genomic regions |

**Disease genes, frameshift** — `gene level_v4/disease genes/framesift/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `bind_result_dbh.R` | 554 | Binds PTC/NMD results with gene info and combines p-values |
| `compare+1-1.R` | 151 | Compares plus1 vs plus2 frameshift results by transcript |
| `frameshift_code.R` | 334 | Classifies frameshift variants and predicts NMD escape from ClinVar |
| `fs_transcript_level.R` | 214 | Computes PTC and NMD-escape indicators at transcript level |
| `process_syn.R` | 63 | Counts synonymous variants falling within specified regions |
| `variant_level.R` | 191 | Determines strand, locates variants in exons, classifies frameshift type |

**Control genes, SNV** — `gene level_v4/control genes/snv/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `get_gnomAD_control.R` | 175 | Builds gnomAD control SNV PTC set excluding snv_can genes |
| `get_gnomAD_control_step2.R` | 188 | Builds gnomAD control SNV PTC set excluding NMD_can_esc genes |
| `get_snv_control_gene.R` | 33 | Filters clinvar SNV PTC genes to AD canonical-transcript control set |
| `get_snv_gnomAD_control_step2.R` | 189 | Builds gnomAD control SNV PTC set excluding NMD_can_esc genes |

**Control genes, frameshift** — `gene level_v4/control genes/frameshift/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `get_fs_control_gene.R` | 79 | Filters control gene list excluding NMD-escape frameshift genes and single-exon genes |

**Gene features** — `gene level_v4/features/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `GO_enrich.R` | 59 | Runs GO enrichment on ClinVar gene sets and writes results |
| `big_pIc.R` | 120 | Builds big-picture summary plots across ClinVar gene categories |
| `gene_NmdescRegion.R` | 47 | Plots NMDesc region length and region-to-CDS length ratio |
| `gene_PliLoeuf.R` | 209 | Plots pLI, LOEUF, exon count and CDS length by gene group |
| `gene_connectivity.R` | 76 | Computes PPI degree centrality for disease and control genes |
| `gene_gc.R` | 283 | Computes gene-level and NMDesc region GC content and plots it |
| `gene_motif.R` | 77 | Counts motifs per gene from external supplementary datasets |
| `gene_overview.R` | 29 | Summarizes variable distributions across ClinVar variant dataframes |
| `gene_pfam.R` | 220 | Computes per-gene Pfam domain overlap and plots results |
| `gene_plot.R` | 312 | Plots gene- and variant-level features across fs/snv groups |
| `gene_plot_match.R` | 327 | Plots matched gene-level enrichment with McNemar significance |
| `gene_plot_match2.R` | 346 | Plots four-panel matched gene-level enrichment figure |
| `gene_ppi.R` | 61 | Plots percentage of genes present in PPI network |
| `gene_repeat.R` | 173 | Computes repeat and homopolymer content for gene sequences |
| `gene_tau.R` | 108 | Computes tau tissue-specificity scores from gene expression matrix |
| `inheritance.R` | 109 | Checks OMIM inheritance pattern for a gene list |
| `new_loc.R` | 102 | Derives new gene feature locations from hg38 sequence data |
| `plot_supplemental_fig3_gene_level.R` | 335 | Builds Supplemental Figure 3 gene-level enrichment barplot |

**Gene feature functions** — `gene level_v4/features/functions/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `annotate_motif_flags.R` | 113 | Annotates genes with motif flags and max scores from UniProt |
| `build_gene_all.R` | 185 | Builds combined gene-level table from PTC and case/control gene files |
| `calculate_ppi_degree_centrality.R` | 86 | Maps genes to STRING IDs and computes PPI degree centrality |
| `plot_gc_content.R` | 55 | Computes and plots GC content across gene regions |
| `plot_gene_level_features.R` | 91 | Merges LOF metrics, bins pLI/LOEUF, and plots gene-level features |
| `plot_repeat_content.R` | 87 | Detects repeat and homopolymer fractions and plots repeat content |
| `run_pfam_overlap_analysis.R` | 142 | Overlaps PFAM domains with NMD-escape regions per record |
| `run_ppi_overlap_analysis.R` | 74 | Flags genes whose PPI interface residues overlap NMD-escape region |
| `run_tau_analysis.R` | 95 | Computes tissue-specificity tau from GTEx expression data |

**QC and negative controls** — `gene level_v4/QC/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `check.R` | 146 | Parses variant positions and computes CDS distance for QC |
| `combine_gene.R` | 305 | Combines gene-level PTC/NMDesc region and enrichment data |
| `compare.R` | 128 | Compares NMD-enriched variant lists against sequence features |
| `compare_submitters.R` | 14 | Compares number of submitters across variant groups |
| `cross_check.R` | 11 | Cross-checks NMD result data |
| `download_gnomad_syn_vcf.R` | 28 | Streams and filters gnomAD VCF for PASS synonymous variants |
| `download_syn2.R` | 38 | Downloads and processes synonymous variant data |
| `negative_control.R` | 627 | Finds plp/benign PTC variants for control gene list analysis |
| `parse_gnomad_syn_vcf.R` | 55 | Parses gnomAD synonymous variants and maps CDS locations |
| `process_syn2.R` | 53 | Processes gnomAD exome synonymous variants against CDS coordinates |
| `select_AD.R` | 46 | Selects AD control genes from candidate gene lists |
| `select_AD_variants.R` | 60 | Selects AD and control gene variants from gnomAD/ClinVar |

**Gene-level entry point** — `gene level_v4/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `gene_main_dbh.R` **entry** | 974 | Matches disease genes to control genes by CDS length and analyzes NMD escape |

### Variant level

**ClinVar SNV variants** — `variant level_v4/clinvar/snv/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `get_snv_control_variant_new.R` | 41 | Builds control variant set from canonical transcripts of control genes |
| `get_snv_variant_new.R` | 59 | Extracts ClinVar SNV PTC/NMD-escape variants for canonical transcripts of gene list |

**ClinVar frameshift variants** — `variant level_v4/clinvar/frameshift/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `get_fs_variant_new.R` | 89 | Extracts frameshift variants from ClinVar, removing in-frame indels |

**gnomAD download** — `variant level_v4/gnomad/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `gnomAD_downloaddata.R` | 83 | Downloads and annotates gnomAD variant data with aenmd NMD escape predictions |

**gnomAD SNV controls** — `variant level_v4/gnomad/snv/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `get_gnomAD_control.R` | 139 | Builds gnomAD-derived control variant sets for NMDesc PTC enrichment comparisons |

**gnomAD frameshift controls** — `variant level_v4/gnomad/frameshift/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `get_fs_control_variant_new.R` | 29 | Builds control frameshift variant set using biomaRt and dplyr |

**Variant features** — `variant level_v4/features/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `clean_motif.R` | 193 | Computes motif counts on variant level from external datasource |
| `create_fasta.R` | 186 | Translates sequences and finds PTC to build FASTA files |
| `csv2vcf.R` | 22 | Converts CSV key column into VCF format |
| `do_AD.R` | 103 | Filters AD genes and prepares data for Venn plot |
| `do_motif.R` | 103 | Indexes FASTA files with positions and batch-plots motif results |
| `het_6fasta.R` | 128 | Retrieves canonical transcript info and writes gene variant CSV |
| `merge_loc_uni.R` | 32 | Merges uni data with loc data into combined CSV |
| `number_submitters.R` | 144 | Compares CDS length between disease and control genes by category |
| `parse.R` | 108 | Standardizes and plots PARSE score data by category |
| `pfam_ppi_variant.R` | 619 | Analyzes and plots Pfam and PPI features on variant level |
| `plus1_control.R` | 119 | Computes NMD escape region from exon and CDS data |
| `transcript_features.R` | 270 | Derives canonical transcript and exon count features |
| `variant_motif.R` | 192 | Computes motif counts on variant level from external datasource |
| `variant_plot.R` | 364 | Plots continuous and binary variant-level motif and LCS features |
| `variant_plot2.R` | 330 | Builds variant-level enrichment barplot with p-value brackets |

**Variant feature functions** — `variant level_v4/features/functions/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `variant_pfam_ppi.R` | 371 | Maps CDS positions to amino acids and tests Pfam/PPI enrichment via Fisher's exact test |
| `variants_motif.R` | 110 | Annotates variants with motif and LCS flags using biomaRt UniProt mapping |

**VEP annotation** — `variant level_v4/features/vep/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `process_vep.R` | 207 | Filters rows where variant and gene match minus1_variants |
| `vep_draw.R` | 752 | Plots VEP/dbNSFP features like CADD, exon distance, protein length change with p-values |

**IDR core** — `variant level_v4/features/IDR/idr/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `get_snv_control_idr.R` | 31 | Joins gene from Uniprot IDs into snv_control_idr for wildtype disordered domain boundaries |
| `idr_merge.R` | 41 | Merges Uniprot_ID-derived gene and Disordered_Domain_Boundaries into snv_control_idr |
| `idr_output.R` | 23 | Converts list-columns to character and writes IDR data to file |
| `quality_filter_idr.R` | 36 | Reads, dedupes, computes nmdesc_length2, and filters IDR text files |

**IDR overlap** — `variant level_v4/features/IDR/idr/overlap/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `get_fs_idr_match.R` | 26 | Matches frameshift variants to IDR regions per gene transcript |
| `get_snv_idr_match.R` | 53 | Matches SNV IDR locations to canonical transcript CDS per gene |
| `get_snv_idr_match2.R` | 55 | Matches SNV IDR locations to CDS using resolved data paths |
| `idr_main.R` | 133 | Calculates IDR disorder changes after NMD for SNVs and frameshifts |

**IDR length comparison** — `variant level_v4/features/IDR/idr/compare_length/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `compare_idr_length.R` | 53 | Parses JSON IDR boundaries and computes IDR length |
| `compare_idr_snv_control.R` | 53 | Parses JSON IDR boundaries and computes control IDR length |
| `idr_difff.R` | 66 | Computes IDR length difference between control and variant by gene |
| `idr_plot.R` | 109 | Combines IDR datasets and plots length differences with p-values |
| `wildtype_idr.R` | 48 | Merges UniProt IDR boundaries into SNV data by gene |

**IDR within-domain** — `variant level_v4/features/IDR/idr/with_in/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `idr_match2.R` | 110 | Flags whether PTC variants fall inside IDR regions |

**Sequence-property models** — `variant level_v4/emel/sequence/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `sequence-parsev2_main.R` | 529 | Fits LMM/GLMM on PARSE features and plots forest plot |
| `sequence-peptides_main.R` | 446 | Fits LMM on peptide properties and plots forest plot |
| `sequence-peptides_read.R` | 162 | Reads peptide sequence data and plots group panels |
| `sequence_cider_main.R` | 779 | Fits LMM on CIDER features by region and plots forest plot |
| `sequence_cider_match.R` | 119 | Computes paired VAR-WT CIDER differences via repeated resampling |
| `sequence_cider_read.R` | 177 | Reads CIDER sequence data and plots group panels |
| `sequence_cider_table1.R` | 229 | Builds summary table of CIDER properties with pairwise p-values |
| `sequence_parsev2_hier.R` | 279 | Fits Bayesian hierarchical model on PARSE data with brms |
| `sequence_parsev2_match.R` | 118 | Computes paired VAR-WT PARSE differences via Wilcoxon test |
| `sequence_parsev2_read.R` | 180 | Reads PARSE v2 data, cleans columns, and plots panels |
| `sequence_parsev2_table1.R` | 296 | Builds summary table of PARSE properties with pairwise p-values |
| `sequence_peptide_hier.R` | 225 | Fits Bayesian hierarchical model on peptide data with brms |
| `sequence_peptide_match.R` | 133 | Computes paired VAR-WT peptide differences via repeated resampling |
| `sequence_peptides_table1.R` | 314 | Builds merged table of VAR-WT peptide property differences |

**Structure-property models** — `variant level_v4/emel/structure/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `structure_PAE.R` | 278 | Plots PAE comparisons with FDR-corrected significance stars |
| `structure_hier.R` | 232 | Fits hierarchical Bayesian models of secondary structure variables |
| `structure_match.R` | 129 | Runs paired Wilcoxon tests on NMDPos VAR-WT structural values |
| `structure_mixed.R` | 314 | Fits linear mixed models and plots forest plot of deltaVARWT |
| `structure_plddt.R` | 274 | Plots pLDDT comparisons with FDR-corrected significance stars |
| `structure_read.R` | 491 | Reads structural feature data and plots pLDDT and SASA panels |
| `structure_table1.R` | 240 | Builds Table 1 summarizing structural variables by cohort |
| `structure_table2.R` | 211 | Builds Table 2 of delta structural variables between cohorts |

**QC, resampling, regression** — `variant level_v4/QC/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `adjust_variant.R` | 547 | Plots CDS length, NMDesc region length, and distance to PTC across groups |
| `check_pogz.R` | 188 | Resolves data file paths and checks POGZ variant data |
| `check_variant_dis.R` | 37 | Filters variant groups to genes and UniProts shared across groups |
| `clean_variant_AD.R` | 41 | Cleans and deduplicates AD variant lists for plus1 strand |
| `combine_variant.R` | 125 | Combines ClinVar PLP and gnomAD benign PTC variants into variant_all |
| `compare_D.R` | 548 | Compares fs_control vs fs_control2 and snv_control vs snv_control2 with plots |
| `plot_negative_control.R` | 372 | Compares control sets on pfam_ppi domains and motifs with plots |
| `regression.R` | 513 | Fits regression models and plots CDS/NMDesc length comparisons across groups |
| `resample.R` | 307 | Resamples matched variants 1000 times, one per transcript per source |

**Variant-level entry point** — `variant level_v4/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `variant_main_DBH.R` **entry** | 852 | Runs variant-level comparison of PTC/NMD-escape features and fits GLMs |
| `new_create_fasta_functions.R` | 353 | Builds FASTA sequences and locates PTCs from variant coordinates |

### Protein level

**FASTA generation** — `protein level_v4/fasta/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `create_fs_control.R` | 142 | Builds frameshift control set excluding NMDesc, non-canonical, single-exon transcripts |
| `new_create_fasta.R` | 460 | Generates FASTA sequences with PTC, GC and repeat content annotations |

**AlphaFold2 figures** — `protein level_v4/AF2/`

| Script | Lines | Purpose |
|--------|-------|---------|
| `AF2_draw.R` | 818 | Plots AA composition, charge, pI, aromaticity across NMD/FL regions |

### Repository root

**Top-level driver** — repository root

| Script | Lines | Purpose |
|--------|-------|---------|
| `run_analysis.R` **entry** | 83 | Locates raw data, checks packages, and runs the analysis pipeline |

