#Filter p/lp, canonical, not ptc, remove transcripts of NMDesc genes

noptc = which(res2@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]]==F)
plp = grep('pathogenic',res2@elementMetadata@listData[["clinsig"]],ignore.case = T)
res_plp_noptc = res2[intersect(noptc,plp)]
del = which(res_plp_noptc@elementMetadata@listData[["type"]] == 'del')
ins = which(res_plp_noptc@elementMetadata@listData[["type"]] == 'ins')
fs_control = res_plp_noptc[union(del,ins)]

transcript_control_set = unique(fs_control@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]])
can.info = getBM(
  attributes = c("ensembl_transcript_id","transcript_is_canonical"), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values =  transcript_control_set,                          # Transcript ID
  mart = ensembl                                         # Database connection
)
can.info2 = can.info[which(can.info$ensembl_transcript_id %in% transcript_control_set),]
#exclude not canonial transcripts
transcript_set2 = can.info2[which(can.info2$transcript_is_canonical==1),"ensembl_transcript_id"]
BM.info <- getBM(
  attributes = c("ensembl_transcript_id","rank",'cds_start','cds_end','exon_chrom_start','exon_chrom_end'), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values =  transcript_control_set,                          # Transcript ID
  mart = ensembl                                         # Database connection
)
#exclude single (coding) exon genes
exon.num = BM.info %>% group_by(ensembl_transcript_id) %>% summarise(max_rank = sum(!is.na(cds_start)))
merge = merge(BM.info,exon.num,by='ensembl_transcript_id')
transcript_set_controlfinal = unique(merge[which(merge$max_rank>1),"ensembl_transcript_id"])

library(GenomicRanges)
library(biomaRt)


# Output storage
distance_list <- data.frame(Variant_Key = character(), Distance = numeric(), ptc_loc = numeric(), cds_mutation_loc = numeric(),stringsAsFactors = FALSE)

# Get all transcript IDs from fs_control2
transcript_ids <- unique(sample(transcript_set_controlfinal,1000))
fs_control2 = fs_control[which(fs_control@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% transcript_ids)]

# 1. Batch query: Get CDS sequences for all transcripts
cds_seq_df <- getBM(
  attributes = c("ensembl_transcript_id", "coding"),
  filters = "ensembl_transcript_id",
  values = transcript_ids,
  mart = ensembl
)
cds_seq_map <- setNames(cds_seq_df$coding, cds_seq_df$ensembl_transcript_id)

# 2. Batch query: Get exon data for all transcripts
exon_data_df <- getBM(
  attributes = c("ensembl_transcript_id", "chromosome_name", "exon_chrom_start", 
                 "cds_start", "cds_end", "exon_chrom_end", "rank", "strand"),
  filters = "ensembl_transcript_id",
  values = transcript_ids,
  mart = ensembl
)
exon_data_df <- exon_data_df[complete.cases(exon_data_df),]  # Remove rows with NA
exon_data_df <- exon_data_df[order(exon_data_df$cds_start),]  # Sort by CDS start
exon_data_map <- split(exon_data_df, exon_data_df$ensembl_transcript_id)  # Split for lookup
output_dir <- "fs_control1"
# 3. Process each variant
for (i in 1:length(fs_control2)) {
  # 1. Get genomic pos from res file
  transcript_id <- fs_control2@ranges@NAMES[i]
  key = fs_control2@elementMetadata@listData[["id"]][i]
  variant_info <- unlist(strsplit(key, "_"))
  chromosome <- variant_info[1]
  position <- as.numeric(variant_info[2])  # Genomic position
  ref_allele <- variant_info[3]
  alt_allele <- variant_info[4]
  
  # 2. Get original CDS sequence
  if (!transcript_id %in% names(cds_seq_map)) next  # Skip if missing
  original_seq <- cds_seq_map[[transcript_id]]
  
  # 3. Get exon data
  if (!transcript_id %in% names(exon_data_map)) next  # Skip if missing
  exon_data <- exon_data_map[[transcript_id]]
  
  # 4. Transform genomic position to CDS position
  get_re_loc <- function(location) {
    cumulative_cds_length <- 0
    cds_position <- NA
    
    for (i in seq_along(exon_data$cds_start)) {
      if (location >= exon_data$exon_chrom_start[i] && location <= exon_data$exon_chrom_end[i]) {
        cds_position <- cumulative_cds_length + (location - exon_data$exon_chrom_start[i] + 1)
        break
      }
      cumulative_cds_length <- cumulative_cds_length + (exon_data$cds_end[i] - exon_data$cds_start[i] + 1)
    }
    return(cds_position)
  }
  cds_mutation_loc <- get_re_loc(position)
  
  # 5. Modify coding sequence
  mutated_seq <- paste0(
    substr(original_seq, 1, cds_mutation_loc - 1),  # Before mutation
    alt_allele,  # Insert mutation
    substr(original_seq, cds_mutation_loc + nchar(ref_allele), nchar(original_seq))  # After mutation
  )
  
  # 6. Generate new protein sequence, cut at first PTC
  find_ptc <- function(sequence) {
    stop_codons <- c("TAA", "TAG", "TGA")  # Stop codons in cDNA
    seq_length <- nchar(sequence)
    
    for (i in seq(1, seq_length - 2, by = 3)) {
      codon <- substr(sequence, i, i + 2)
      if (codon %in% stop_codons) {
        return(i)  # Return position of the stop codon
      }
    }
    return(NA)  # No PTC found
  }
  
  ptc_pos <- find_ptc(mutated_seq)
  
  if (!is.na(ptc_pos)) {
    truncated_seq <- substr(mutated_seq, 1, ptc_pos + 2)  # Include stop codon
  } else {
    truncated_seq <- mutated_seq  # No PTC, keep full sequence
  }
  
  # 7. Store distance between mutation and PTC
  distance <- if (!is.na(ptc_pos)) ptc_pos - cds_mutation_loc else NA
  distance_list <- rbind(distance_list, data.frame(Variant_Key = key, Distance = distance, ptc_loc = ptc_pos,cds_mutation_loc=cds_mutation_loc, stringsAsFactors = FALSE))
  
  # 8. Write to individual FASTA file
  fasta_filename <- file.path(output_dir, paste0(key, ".fasta"))
  sink(fasta_filename)
  cat(">", key, "\n", truncated_seq, "\n", sep="")
  sink()
}

distance_list = distance_list[which(distance_list$Distance>=0),]
mean(distance_list$Distance, na.rm = TRUE)
median(distance_list$Distance, na.rm = TRUE)
hist(distance_list$Distance, breaks = 20, xlab = "Distance", main = "Distribution of Distances between Mutation and PTC")
