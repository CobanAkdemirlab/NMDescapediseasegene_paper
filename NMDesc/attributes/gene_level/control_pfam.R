# Query Pfam domain locations from Ensembl
snv_control1_pfam <- getBM(
  attributes = c("ensembl_transcript_id", "pfam", 
                 "pfam_start", "pfam_end"),
  filters = "ensembl_transcript_id",
  values = snv_control1_tr6$ensembl_transcript_id,  # TP53 Gene ID
  mart = ensembl
)

# Merge with snv_control1_tr6 by ensembl_transcript_id, keeping all rows from snv_control1_tr6
snv_control1_pfam2 <- merge(snv_control1_pfam, snv_control1_tr6, by = "ensembl_transcript_id", all.y = TRUE)

# Identify overlaps between Pfam domain and aft.NMD.ind region
snv_control1_pfam2$pfam_match <- ifelse(
  snv_control1_pfam2$pfam_end >= snv_control1_pfam2$aft.NMD.ind / 3 & 
    snv_control1_pfam2$pfam_start <= snv_control1_pfam2$aft.NMD.ind / 3, 
  1, 0
)

# Keep only matched rows
snv_control1_pfam3 <- snv_control1_pfam2[snv_control1_pfam2$pfam_match == 1, ]

# Count unique HGNC symbols
length(unique(snv_control1_pfam3$hgnc_symbol))
