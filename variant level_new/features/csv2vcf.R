# Split `key` into CHROM, POS, REF, ALT using regular expression
vcf_df <- snv_control_variants %>%
  extract(key, into = c("CHROM", "POS", "REF", "ALT"), regex = "([^:]+):([0-9]+)\\|([^|]+)\\|(.+)", remove = FALSE) %>%
  mutate(
    ID = ".",
    QUAL = ".",
    FILTER = "PASS",
    INFO = paste0("TRANSCRIPT=", transcript, ";UNIPROT=", uniprotswissprot)
  ) %>%
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)

# Write VCF header
vcf_header <- c(
  "##fileformat=VCFv4.2",
  '##INFO=<ID=TRANSCRIPT,Number=1,Type=String,Description="Ensembl Transcript ID">',
  '##INFO=<ID=UNIPROT,Number=1,Type=String,Description="UniProt Swiss-Prot ID">',
  paste0("#", paste(colnames(vcf_df), collapse = "\t"))
)

# Write to VCF file
writeLines(c(vcf_header, apply(vcf_df, 1, paste, collapse = "\t")), "snv_control.vcf")
