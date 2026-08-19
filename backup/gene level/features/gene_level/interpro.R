

# Query InterPro domain locations
snv_can_interpro <- getBM(
  attributes = c("ensembl_transcript_id", "interpro", "interpro_description", 
                 "interpro_start", "interpro_end"),
  filters = "ensembl_transcript_id",
  values = snv_can_tr6$ensembl_transcript_id,  # TP53 Gene ID
  mart = ensembl
)

#select the rows that can NMDesc region overlap with interpro_start and interpro_end
##merge with snv_can_tr6, by ensembl_transcript_id, keep all rows in snv_can_interpro
snv_can_interpro2 <- merge(snv_can_interpro, snv_can_tr6, by = "ensembl_transcript_id", all.y = TRUE)
#if end > aft.NMD.ind & start < aft.NMD.ind then match
snv_can_interpro2$interpro_match = ifelse(snv_can_interpro2$interpro_end >= snv_can_interpro2$aft.NMD.ind/3 & snv_can_interpro2$interpro_start <= (snv_can_interpro2$aft.NMD.ind/3-50), 1, 0)
#keep matched rows
snv_can_interpro3 <- snv_can_interpro2[snv_can_interpro2$interpro_match == 1,]
length(unique(snv_can_interpro3$hgnc_symbol))
write.csv(snv_can_interpro3, file = "snv_can_interpro3.csv", row.names = FALSE)

#transform protein loc to chrom loc 
# Check the result
library(biomaRt)

# Function to map protein position to genomic coordinate
map_protein_to_genome <- function(transcript_id, protein_position) {
  # Retrieve exon data for the given transcript
  exon_data <- getBM(
    attributes = c(
      'ensembl_transcript_id', 'ensembl_peptide_id',
      'chromosome_name', 'exon_chrom_start', 'exon_chrom_end',
      'cds_start', 'cds_end', 'strand', 'phase'
    ),
    filters = 'ensembl_transcript_id',
    values = transcript_id,
    mart = ensembl
  )
  #remove rows contain NA
  exon_data = exon_data[complete.cases(exon_data),]
  
  if (nrow(exon_data) == 0) {
    stop("No exon data found for the provided transcript ID.")
  }
  
  # Convert protein position to CDS position
  cds_position <- protein_position * 3 + 1
  
  # Identify the exon containing the CDS position
  mapped_exon <- exon_data[
    exon_data$cds_start <= cds_position & exon_data$cds_end >= cds_position, 
  ]
  
  
  if (nrow(mapped_exon) == 0) {
    return("CDS position not found within the exons.")
  }
  
  # Map CDS position to genomic position based on strand orientation
  if (mapped_exon$strand[1] == 1) {  # Positive strand
    genomic_position <- mapped_exon$exon_chrom_start + (cds_position - mapped_exon$cds_start)
  } else {  # Negative strand
    genomic_position <- mapped_exon$exon_chrom_end - (cds_position - mapped_exon$cds_start)
  }
  
  # Return the genomic coordinates
  result <- list(
    chromosome = mapped_exon$chromosome_name[1],
    genomic_position = genomic_position,
    strand = ifelse(mapped_exon$strand[1] == 1, "+", "-"),
    chrom = mapped_exon$chromosome_name[1]
  )
  
  return(result)
}
map_protein_to_genome(snv_can_interpro3$ensembl_transcript_id[1], snv_can_interpro3$C_start_IDR[1])

for(i in 1:10){
  #store the output to snv_can_interpro3
  snv_can_interpro3$genomic_IDR_start[i] = map_protein_to_genome(snv_can_interpro3$ensembl_transcript_id[i], snv_can_interpro3$C_start_IDR[i])$genomic_position
  snv_can_interpro3$genomic_IDR_end[i] = map_protein_to_genome(snv_can_interpro3$ensembl_transcript_id[i], snv_can_interpro3$C_end_IDR[i])$genomic_position
  snv_can_interpro3$strand[i] = map_protein_to_genome(snv_can_interpro3$ensembl_transcript_id[i], snv_can_interpro3$C_start_IDR[i])$strand
  snv_can_interpro3$chrom[i] = map_protein_to_genome(snv_can_interpro3$ensembl_transcript_id[i], snv_can_interpro3$C_start_IDR[i])$chrom
}
snv_can_interpro3 = snv_can_interpro3[,-17]
write.csv(snv_can_interpro3[1:10,], file = "test_paint.csv", row.names = FALSE)
# Example Usage
# transcript_id <- "ENST00000357654"
# protein_position <- 150
# map_protein_to_genome(transcript_id, protein_position)




#regroup the description by gene

