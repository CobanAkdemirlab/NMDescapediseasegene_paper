res = readRDS('~/Desktop/new_clinvar/raw_data/clinvar_1112.rds')
res2 = unlist(res)
plp = grep('pathogenic',res2@elementMetadata@listData[["clinsig"]],ignore.case = T)

#germline

#PTC
ptc = which(res2@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]]==T)

#snv
snv = which(res2@elementMetadata@listData[["type"]]=='snv')
fs = which(res2@elementMetadata@listData[["type"]] %in% c('ins','del'))

#remove single exon
no_single = which(res2@elementMetadata@listData[["res_aenmd"]]@listData[["is_single"]] == F)

combined_ind = intersect(plp,ptc)
snv_ind = intersect(combined_ind,snv)
snv_ind2 = intersect(snv_ind,no_single)

#get the canonical transcript
transcript_list = unique(res2@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][snv_ind2])

gene_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "transcript_is_canonical", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = transcript_list,
  mart = mart
)
gene_list = gene_mapping[which(gene_mapping$transcript_is_canonical == 1),]

#exclude any gene genes with any NMD_can_esc mutation
snv_res = res2[snv_ind2]
snv_res_can = snv_res[which(snv_res@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% gene_list$ensembl_transcript_id),]
gene2remove.ind = which(snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)
gene2remove = unique(snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][gene2remove.ind])
#remove genes2remove remove from gene_list$ensemble_transcript_id
genes_remain_in_clinvar = gene_list[-which(gene_list$ensembl_transcript_id %in% gene2remove),]
snv_gnomAD_control = snv_NMD_results[which(snv_NMD_results$transcript %in% genes_remain_in_clinvar$ensembl_transcript_id),]
write.csv(snv_gnomAD_control, 'snv_gnomAD_control_variant.csv', row.names = FALSE)
#get hgnc_symbol of snv_gnomAD_control
snv_gnomAD_control_hgnc = getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = snv_gnomAD_control$transcript,
  mart = mart
)
write.csv(unique(snv_gnomAD_control_hgnc$hgnc_symbol), 'snv_gnomAD_control_genes.csv', row.names = FALSE)
#get high frequency variants using annovar



# 1. Extract necessary information from the 'id' column
# Split the id column into CHROM, POS, REF, and ALT
# Correctly extract CHROM, POS, REF, and ALT from the id column
snv_NMD_results3_vcf <- snv_NMD_results3 %>%
  mutate(
    CHROM = sapply(strsplit(id, ":"), `[`, 1),  # Extract chromosome
    POS = sapply(strsplit(sapply(strsplit(id, ":"), `[`, 2), "\\|"), `[`, 1),  # Extract position
    REF = sapply(strsplit(id, "\\|"), `[`, 2),  # Extract REF allele
    ALT = sapply(strsplit(id, "\\|"), `[`, 3)   # Extract ALT allele
  )


# 2. Add VCF mandatory columns
snv_NMD_results3_vcf$ID <- "."  # No variant ID
snv_NMD_results3_vcf$QUAL <- "."  # No quality score
snv_NMD_results3_vcf$FILTER <- "PASS"  # All variants pass filter
snv_NMD_results3_vcf$INFO <- "."  # No additional information

# 3. Select and reorder the columns in the standard VCF format
vcf_output <- snv_NMD_results3_vcf %>%
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)

# 4. Write VCF header and data to a file
vcf_file <- "snv_NMD_results3.vcf"

# Write VCF header
vcf_header <- c(
  "##fileformat=VCFv4.2",
  "##source=snv_NMD_results3",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
)

# Write the header and data to the VCF file
writeLines(vcf_header, con = vcf_file)

# Append VCF data
write.table(
  vcf_output,
  file = vcf_file,
  append = TRUE,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# ✅ 5. Confirm the file was created
cat("VCF file created successfully at:", vcf_file, "\n")



snv_NMD_results4 = 
#select random 200 rows of snv_NMD_results
snv_NMD_result2 = snv_NMD_results3[sample(nrow(snv_NMD_results3), 200),]

# Load required libraries
library(data.table)
library(dplyr)

# 1. Prepare Variant Data (Example)
# Replace this with your actual transcript ID variant data
# Assuming you have variant data extracted from transcript IDs
transcript_variants <- data.frame(
  CHR = c("1", "2", "3"),  # Chromosome
  POS = c(123456, 234567, 345678),  # Position
  REF = c("A", "G", "T"),  # Reference allele
  ALT = c("T", "C", "A"),  # Alternate allele
  stringsAsFactors = FALSE
)

# Create additional fields for ANNOVAR input
transcript_variants$ID <- rep(".", nrow(transcript_variants))  # Placeholder for ID
transcript_variants$ZYG <- rep("Het", nrow(transcript_variants))  # Zygosity
transcript_variants$FILTER <- rep("PASS", nrow(transcript_variants))  # Filter status

# 2. Write VCF File for ANNOVAR
vcf_file <- "variants_for_annovar.vcf"
annovar_input <- cbind(
  transcript_variants$CHR,
  transcript_variants$POS,
  transcript_variants$POS,
  transcript_variants$REF,
  transcript_variants$ALT,
  transcript_variants$ID,
  transcript_variants$ZYG,
  transcript_variants$FILTER
)

# Save the file in VCF format
write.table(
  annovar_input, 
  file = vcf_file, 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = FALSE
)

# 3. Convert VCF to ANNOVAR format using system() call
# Replace '/path/to/annovar/' with the actual ANNOVAR installation directory
system("perl /path/to/annovar/convert2annovar.pl -format vcf4 variants_for_annovar.vcf > variants_for_annovar.avinput")

# 4. Run ANNOVAR to Annotate Variant Frequencies
# Replace '/path/to/annovar/' and '/path/to/annovar/humandb/' with actual paths
system("perl /path/to/annovar/annotate_variation.pl -build hg38 -out variant_frequency_results -dbtype gnomad211_exome variants_for_annovar.avinput /path/to/annovar/humandb/")

# 5. Read ANNOVAR Output into R
# Assuming the output file is named 'variant_frequency_results.gnomad211_exome_dropped'
variant_freq_file <- "variant_frequency_results.gnomad211_exome_dropped"

if (file.exists(variant_freq_file)) {
  variant_freq <- fread(variant_freq_file)
  
  # 6. Merge Frequency Data with Original Variant Data
  # Construct keys for merging
  transcript_variants$key <- paste0(
    transcript_variants$CHR, ":", transcript_variants$POS, "_", transcript_variants$REF, ">", transcript_variants$ALT
  )
  
  variant_freq$key <- paste0(
    variant_freq$V1, ":", variant_freq$V2, "_", variant_freq$V4, ">", variant_freq$V5
  )
  
  # Merge frequency annotation into the original dataset
  merged_data <- merge(transcript_variants, variant_freq, by = "key", all.x = TRUE)
  
  # 7. Save the Annotated Variants with Frequencies to CSV
  write.csv(merged_data, "variant_frequencies_with_transcripts.csv", row.names = FALSE)
  
  # Print the merged data
  print("Variant annotation completed. Frequency data saved to CSV.")
  head(merged_data)
  
} else {
  stop("ANNOVAR output file not found. Please check if the annotation ran successfully.")
}
