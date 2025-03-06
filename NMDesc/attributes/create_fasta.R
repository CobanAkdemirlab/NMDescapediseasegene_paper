create_fasta = function(GR,output_dir = "snv_control2"){
  reverse_complement <- function(sequence) {
    sequence <- rev(strsplit(sequence, NULL)[[1]])  # Reverse sequence
    complement_map <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G")
    sequence <- sapply(sequence, function(base) complement_map[base])
    return(paste(sequence, collapse = ""))
  }
  
  translate_sequence <- function(dna_seq) {
    # Define standard genetic code
    codon_table <- c(
      "TTT" = "F", "TTC" = "F", "TTA" = "L", "TTG" = "L",
      "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
      "ATT" = "I", "ATC" = "I", "ATA" = "I", "ATG" = "M",
      "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
      "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
      "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
      "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
      "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
      "TAT" = "Y", "TAC" = "Y", "TAA" = "", "TAG" = "",
      "CAT" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
      "AAT" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
      "GAT" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
      "TGT" = "C", "TGC" = "C", "TGA" = "", "TGG" = "W",
      "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
      "AGT" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
      "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
    )
    
    # Split sequence into codons
    codons <- substring(dna_seq, seq(1, nchar(dna_seq) - 2, by = 3), seq(3, nchar(dna_seq), by = 3))
    
    # Translate codons to amino acids
    protein_seq <- sapply(codons, function(codon) codon_table[[codon]], USE.NAMES = FALSE)
    
    # Join translated amino acids into a protein sequence
    return(paste(protein_seq, collapse = ""))
  }
  # Initialize dataframe to store results
  distance_list <- data.frame(Variant_Key = character(), strand = numeric(),
                              Distance = numeric(), cds_mutation_loc = numeric(),
                              ptc_loc = numeric(), stringsAsFactors = FALSE)
  #select the ones in gene_list2(clinvar can control)
  snv_NMD_results3 = snv_NMD_results[which(snv_NMD_results$transcript %in% gene_list2$ensembl_transcript_id),]
  #select random 200 rows of snv_NMD_results
  snv_NMD_result2 = snv_NMD_results3[sample(nrow(snv_NMD_results3), 200),]
  
  # Get all transcript IDs from GR
  #transcript_ids <- unique(GR@ranges@NAMES)
  transcript_ids = snv_NMD_result2$transcript
  
  #ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # 1. Get CDS sequences for all transcripts
  cds_seq_df <- getBM(
    attributes = c("ensembl_transcript_id", "coding"),
    filters = "ensembl_transcript_id",
    values = transcript_ids,
    mart = ensembl
  )
  cds_seq_map <- setNames(cds_seq_df$coding, cds_seq_df$ensembl_transcript_id)
  
  # 2. Get exon data for all transcripts
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
  
  # Function to map genomic position to CDS position
  get_re_loc <- function(location, exon_data) {
    cumulative_cds_length <- 0
    cds_position <- NA
    
    for (i in seq_along(exon_data$cds_start)) {
      if (location >= exon_data$exon_chrom_start[i] && location <= exon_data$exon_chrom_end[i]) {
        if(strand == -1) {
          cds_position <- cumulative_cds_length + (exon_data$exon_chrom_end[i] - location + 1)
        } else
          cds_position <- cumulative_cds_length + (location - exon_data$exon_chrom_start[i] + 1)
        break
      }
      cumulative_cds_length <- cumulative_cds_length + (exon_data$cds_end[i] - exon_data$cds_start[i] + 1)
    }
    return(cds_position)
  }
  
  
  # Function to find the first premature termination codon (PTC)
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
  
  # 3. Process each variant
  #for (i in seq_along(GR)) {
  for(i in 1:nrow(snv_NMD_result2)){
    # Extract variant details
    #transcript_id <- GR@ranges@NAMES[i]
    transcript_id = snv_NMD_result2$transcript[i]
    #key <- GR@elementMetadata@listData[["id"]][i]
    key = snv_NMD_result2$id[i]
    #variant_info <- unlist(strsplit(key, "_"))
    # Assuming 'key' contains a string like "3:151436983|G|T"
    variant_info <- unlist(strsplit(key, ":"))  # Split chromosome from the rest
    chromosome <- variant_info[1]  # "3"
    
    # Now split the rest by '|'
    position_alleles <- unlist(strsplit(variant_info[2], "\\|"))  # Need double backslash for regex
    position <- as.numeric(position_alleles[1])  # 151436983
    ref_allele <- position_alleles[2]  # "G"
    alt_allele <- position_alleles[3]  # "T"
    
    #chromosome <- variant_info[1]
    #position <- as.numeric(variant_info[2])  # Genomic position
    #ref_allele <- variant_info[3]
    #alt_allele <- variant_info[4]
    
    # Get original CDS sequence
    if (!transcript_id %in% names(cds_seq_map)) next  # Skip if missing
    original_seq <- cds_seq_map[[transcript_id]]
    
    # Get exon data
    if (!transcript_id %in% names(exon_data_map)) next  # Skip if missing
    exon_data <- exon_data_map[[transcript_id]]
    strand <- exon_data$strand[1]
    
    # Adjust for negative strand
    if (strand == -1) {
      original_seq <- reverse_complement(original_seq)
    }
    
    # Map genomic position to CDS position
    cds_mutation_loc <- get_re_loc(position, exon_data)
    if (strand == -1) {
      cds_mutation_loc <- nchar(original_seq) - cds_mutation_loc + 1
    }
    # Modify coding sequence
    mutated_seq <- paste0(
      substr(original_seq, 1, cds_mutation_loc - 1),  # Before mutation
      alt_allele,  # Insert mutation
      substr(original_seq, cds_mutation_loc + nchar(ref_allele), nchar(original_seq))  # After mutation
    )
    
    # Reverse complement back if on negative strand
    if (strand == -1) {
      mutated_seq <- reverse_complement(mutated_seq)
    }
    
    # Identify premature termination codon (PTC)
    ptc_pos <- find_ptc(mutated_seq)
    if (strand == -1) {
      ptc_pos <- nchar(mutated_seq) - ptc_pos + 1
    }
    
    # Truncate sequence at PTC if found
    truncated_seq <- if (!is.na(ptc_pos)) substr(mutated_seq, 1, ptc_pos + 2) else mutated_seq
    
    # Compute distance from mutation to PTC
    distance <- if (!is.na(ptc_pos)) ptc_pos - cds_mutation_loc else NA
    
    # Store results
    distance_list <- rbind(distance_list, data.frame(
      Variant_Key = key, strand = strand, ptc_pos = ptc_pos, 
      cds_mutation_loc = cds_mutation_loc, Distance = distance, stringsAsFactors = FALSE
    ))
    protein_seq <- translate_sequence(truncated_seq)
    # Write to FASTA file
    fasta_filename <- file.path(output_dir, paste0(key, ".fasta"))
    sink(fasta_filename)
    cat(">", key, "\n", protein_seq, "\n", sep = "")
    sink()
  }
  return(distance_list)
}