# Load required libraries
library(biomaRt)
library(dplyr)
library(readr)

# Get transcript-to-gene mapping from BioMart
transcript_gene_map <- getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  mart = ensembl
)

# Rename columns for merging
colnames(transcript_gene_map) <- c("transcript", "gene")



snv_gnomAD_control_variant <- read_csv("~/Downloads/variant list/snv_gnomAD_control_variant.csv")
snv_control_variant <- snv_gnomAD_control_variant %>%
  left_join(transcript_gene_map, by = "transcript")

snv_control_ADgenes <- read_csv("~/Desktop/autism/data/snv_control_ADgenes.csv")
# Keep only variants where gene is in snv_can_ADgenes$gene
filtered_variants <- snv_control_variant[which(snv_control_variant$gene %in% snv_control_ADgenes$gene),]
nrow(filtered_variants)
#write result
write.csv(filtered_variants,file = "/Users/jxu14/Desktop/autism/data/snv_control_ADvariants.csv", row.names = FALSE)



tr_id = data.frame(transcript = res2@ranges@NAMES,
id = res2@elementMetadata@listData[["key"]])
tr_id_snv = tr_id[which(tr_id$id %in% snv_variants$x),]

snv_variant <- tr_id_snv  %>%
  distinct(id, .keep_all = TRUE)

snv_variant2 <- snv_variant %>%
  left_join(transcript_gene_map, by = "transcript")

filtered_variants <- snv_variant2[which(snv_variant2$gene %in% snv_can_ADgenes$gene),]
nrow(filtered_variants)
#write result
write.csv(filtered_variants,file = "/Users/jxu14/Desktop/autism/data/snv_ADvariants.csv", row.names = FALSE)




