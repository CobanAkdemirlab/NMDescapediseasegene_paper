
# --- 路径解析层（自动插入）------------------------------------------------
# 原来这里是旧机器 /Users/jxu14/ 的绝对路径，换机器必断。改走 paths.R：
#   data_file("x.csv")  按文件名定位，找不到会明确报错
#   out_file("y.csv")   输出到 NMDESC_OUT（默认 ~/Desktop/NMDesc_out）
#   data_root("clinvar") 需要目录而非文件时用
.p <- c("gene level_v3/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../lib/paths.R", "../../../../lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p)) stop("找不到 paths.R —— 请从仓库根目录运行 R")
source(.p[1]); rm(.p)
# --------------------------------------------------------------------------

#1. input human_1
human_1_ <- read_delim(data_file("human (1).txt"), 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

#2. find POGZ related records

##2.1 find uniprot id for pogz
pogz_uniprot <- "Q7Z3K3"

##2.2 match by uniprot id
pogz_ppi_raw <- human_1_ %>%
  filter(uniprot1 == pogz_uniprot | uniprot2 == pogz_uniprot) %>%
  mutate(
    ppi_partner = case_when(
      uniprot1 == pogz_uniprot ~ uniprot2,
      uniprot2 == pogz_uniprot ~ uniprot1
    ),
    
    # Select the interface-residue column corresponding specifically to POGZ
    pogz_interface_residues = case_when(
      uniprot1 == pogz_uniprot ~ interface_residues1,
      uniprot2 == pogz_uniprot ~ interface_residues2
    )
  ) %>%
  select(
    uniprot1, uniprot2, source,
    ppi_partner,
    pogz_interface_residues
  )

#1389-1402

#3. biomart 1302

# ── 3. Retrieve POGZ Pfam-domain coordinates ──────────────────────────────────
pogz_pfam_raw <- getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "transcript_is_canonical",
    "uniprotswissprot",
    "pfam",
    "pfam_start",
    "pfam_end"
  ),
  filters = "hgnc_symbol",
  values = "POGZ",
  mart = ensembl
)

# Keep only rows with a Pfam annotation
pogz_pfam <- pogz_pfam_raw %>%
  filter(
    !is.na(pfam),
    pfam != "",
    !is.na(pfam_start),
    !is.na(pfam_end)
  ) %>%
  mutate(
    pfam_start = as.integer(pfam_start),
    pfam_end   = as.integer(pfam_end),
    domain_type = "Pfam"
  ) %>%
  arrange(ensembl_transcript_id, pfam_start, pfam_end) %>%
  distinct()

print(pogz_pfam, n = Inf)

# ── 4. Identify the reviewed POGZ UniProt accession ───────────────────────────
# Canonical POGZ UniProt ID
pogz_uniprot <- "Q7Z3K3"

# Check whether BioMart returns the same UniProt ID
pogz_pfam %>%
  distinct(hgnc_symbol, ensembl_transcript_id,
           transcript_is_canonical, uniprotswissprot)

# ── 5. Extract POGZ PPI-interface records ─────────────────────────────────────
pogz_ppi_raw <- human_1_ %>%
  filter(uniprot1 == pogz_uniprot | uniprot2 == pogz_uniprot) %>%
  mutate(
    ppi_partner = case_when(
      uniprot1 == pogz_uniprot ~ uniprot2,
      uniprot2 == pogz_uniprot ~ uniprot1
    ),
    pogz_interface_residues = case_when(
      uniprot1 == pogz_uniprot ~ interface_residues1,
      uniprot2 == pogz_uniprot ~ interface_residues2
    )
  ) %>%
  select(
    uniprot1, uniprot2, source,
    ppi_partner,
    pogz_interface_residues
  )

print(pogz_ppi_raw, n = Inf)

# ── 6. Convert POGZ PPI residues into individual amino-acid positions ─────────
# Example: "[269, 270, 271]" becomes three rows: 269, 270, 271

pogz_ppi_residues <- pogz_ppi_raw %>%
  mutate(
    residue_list = str_extract_all(
      as.character(pogz_interface_residues),
      "\\d+"
    )
  ) %>%
  unnest_longer(residue_list) %>%
  transmute(
    gene = "POGZ",
    uniprot = pogz_uniprot,
    ppi_partner,
    source,
    aa_position = as.integer(residue_list),
    domain_type = "PPI_interface"
  ) %>%
  distinct() %>%
  arrange(aa_position)

print(pogz_ppi_residues, n = Inf)

# ── 7. Summarize PPI-interface locations by interaction partner ───────────────
# NOTE: This gives the range spanning observed interface residues.
# It is NOT a true protein domain annotation.

pogz_ppi_ranges <- pogz_ppi_residues %>%
  group_by(ppi_partner, source) %>%
  summarise(
    ppi_start = min(aa_position),
    ppi_end = max(aa_position),
    n_interface_residues = n_distinct(aa_position),
    interface_positions = paste(sort(unique(aa_position)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(ppi_start, ppi_end)

print(pogz_ppi_ranges, n = Inf)

# ── 8. Create one combined table for plotting/annotation ──────────────────────
pogz_domain_locations <- bind_rows(
  pogz_pfam %>%
    transmute(
      gene = hgnc_symbol,
      uniprot = uniprotswissprot,
      feature_type = "Pfam_domain",
      feature_name = pfam,
      start_aa = pfam_start,
      end_aa = pfam_end,
      ppi_partner = NA_character_,
      source = NA_character_
    ),
  
  pogz_ppi_ranges %>%
    transmute(
      gene = "POGZ",
      uniprot = pogz_uniprot,
      feature_type = "PPI_interface_range",
      feature_name = paste0("PPI with ", ppi_partner),
      start_aa = ppi_start,
      end_aa = ppi_end,
      ppi_partner = ppi_partner,
      source = source
    )
) %>%
  arrange(start_aa, end_aa)

print(pogz_domain_locations, n = Inf)

# ── 9. Save outputs ──────────────────────────────────────────────────────────
write_csv(pogz_pfam, "POGZ_Pfam_domains_biomart.csv")
write_csv(pogz_ppi_residues, "POGZ_PPI_interface_residues.csv")
write_csv(pogz_ppi_ranges, "POGZ_PPI_interface_ranges.csv")
write_csv(pogz_domain_locations, "POGZ_combined_Pfam_PPI_locations.csv")