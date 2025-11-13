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

