library(stringr)
library(biomaRt)

# --- helpers ---

get_re_loc <- function(location, exon_data, strand) {
  # Always iterate in CDS (transcript) order: start of CDS -> end of CDS
  exon_data <- exon_data[order(exon_data$cds_start, decreasing = FALSE), ]
  
  cumulative_cds_length <- 0
  for (i in 1:nrow(exon_data)) {
    # Skip non-coding rows if present
    if (is.na(exon_data$cds_start[i]) || is.na(exon_data$cds_end[i])) next
    
    # Is the genomic 'location' within this exon?
    if (location >= exon_data$exon_chrom_start[i] && location <= exon_data$exon_chrom_end[i]) {
      if (strand == -1) {
        # On negative strand, transcript 5'->3' runs from exon_end down to exon_start
        offset_within_exon <- exon_data$exon_chrom_end[i] - location + 1
      } else {
        # On positive strand, transcript 5'->3' runs from exon_start up to exon_end
        offset_within_exon <- location - exon_data$exon_chrom_start[i] + 1
      }
      return(cumulative_cds_length + offset_within_exon)
    } else {
      # Add only the coding portion length in this exon
      cumulative_cds_length <- cumulative_cds_length + (exon_data$cds_end[i] - exon_data$cds_start[i] + 1)
    }
  }
  return(NA)  # if not found
}

reverse_complement <- function(sequence) {
  chars <- strsplit(sequence, NULL)[[1]]
  chars <- rev(chars)
  complement_map <- c("A"="T","T"="A","G"="C","C"="G",
                      "a"="t","t"="a","g"="c","c"="g")
  # fallback to 'N' for any unexpected base
  comp <- vapply(chars, function(b) if (!is.na(complement_map[b])) complement_map[b] else "N", character(1))
  paste(comp, collapse = "")
}

translate_sequence <- function(dna_seq) {
  codon_table <- c(
    "TTT"="F","TTC"="F","TTA"="L","TTG"="L",
    "CTT"="L","CTC"="L","CTA"="L","CTG"="L",
    "ATT"="I","ATC"="I","ATA"="I","ATG"="M",
    "GTT"="V","GTC"="V","GTA"="V","GTG"="V",
    "TCT"="S","TCC"="S","TCA"="S","TCG"="S",
    "CCT"="P","CCC"="P","CCA"="P","CCG"="P",
    "ACT"="T","ACC"="T","ACA"="T","ACG"="T",
    "GCT"="A","GCC"="A","GCA"="A","GCG"="A",
    "TAT"="Y","TAC"="Y","TAA"="","TAG"="",
    "CAT"="H","CAC"="H","CAA"="Q","CAG"="Q",
    "AAT"="N","AAC"="N","AAA"="K","AAG"="K",
    "GAT"="D","GAC"="D","GAA"="E","GAG"="E",
    "TGT"="C","TGC"="C","TGA"="","TGG"="W",
    "CGT"="R","CGC"="R","CGA"="R","CGG"="R",
    "AGT"="S","AGC"="S","AGA"="R","AGG"="R",
    "GGT"="G","GGC"="G","GGA"="G","GGG"="G"
  )
  stop_at <- nchar(dna_seq) - (nchar(dna_seq) %% 3)
  if (stop_at < 3) return("")
  codons <- substring(dna_seq, seq(1, stop_at - 2, by = 3), seq(3, stop_at, by = 3))
  aas <- vapply(codons, function(c) if (!is.null(codon_table[[c]])) codon_table[[c]] else "X", character(1))
  paste(aas, collapse = "")
}

find_ptc <- function(sequence, cut) {
  stop_codons <- c("TAA","TAG","TGA")
  for (i in seq(1, nchar(sequence) - 2, by = 3)) {
    codon <- substr(sequence, i, i + 2)
    # First in-frame stop codon whose END is at/after the cut position
    if (codon %in% stop_codons && (i + 2) >= cut) return(i)  # return codon START index
  }
  return(NA)
}

# --- main ---

create_fasta <- function(variants,output_dir = "plus1_fasta_output") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  distance_list <- data.frame(
    Variant_Key = character(), strand = numeric(),
    cds_mutation_loc = numeric(), ptc_pos = numeric(),
    Distance = numeric(), stringsAsFactors = FALSE
  )
  
  transcript_ids <- variants$transcript
  
  # 1) CDS sequences
  cds_seq_df <- getBM(
    attributes = c("ensembl_transcript_id", "coding"),
    filters    = "ensembl_transcript_id",
    values     = transcript_ids,
    mart       = ensembl
  )
  cds_seq_df <- cds_seq_df[!is.na(cds_seq_df$coding) & cds_seq_df$coding != "", ]
  cds_seq_map <- setNames(cds_seq_df$coding, cds_seq_df$ensembl_transcript_id)
  
  # 2) Exon data
  exon_data_df <- getBM(
    attributes = c("ensembl_transcript_id", "chromosome_name", "exon_chrom_start",
                   "cds_start", "cds_end", "exon_chrom_end", "rank", "strand"),
    filters    = "ensembl_transcript_id",
    values     = transcript_ids,
    mart       = ensembl
  )
  exon_data_df <- exon_data_df[complete.cases(exon_data_df[, c("ensembl_transcript_id","exon_chrom_start","exon_chrom_end","cds_start","cds_end","strand")]), ]
  exon_data_map <- split(exon_data_df, exon_data_df$ensembl_transcript_id)
  
  for (i in 1:nrow(variants)) {
    tryCatch({
      transcript_id <- variants$transcript[i]
      key           <- variants$key[i]
      
      # Parse "chr:pos|ref|alt"
      variant_info     <- strsplit(key, ":", fixed = TRUE)[[1]]
      position_alleles <- strsplit(variant_info[2], "\\|")[[1]]
      position         <- as.numeric(position_alleles[1])
      ref_allele       <- position_alleles[2]
      alt_allele       <- position_alleles[3]
      
      if (!transcript_id %in% names(cds_seq_map)) next
      original_seq <- cds_seq_map[[transcript_id]]
      
      if (!transcript_id %in% names(exon_data_map)) next
      exon_data <- exon_data_map[[transcript_id]]
      strand    <- exon_data$strand[1]
      
      # Map genomic position to CDS position
      cds_mutation_loc <- get_re_loc(position, exon_data, strand)
      if (is.na(cds_mutation_loc)) stop("Could not map genomic position to CDS.")
      
      # Reverse complement alleles for -1 strand
      if (strand == -1) {
        ref_allele <- reverse_complement(ref_allele)
        alt_allele <- reverse_complement(alt_allele)
      }
      
      # Edit CDS
      mutated_seq <- paste0(
        substr(original_seq, 1, cds_mutation_loc - 1),
        alt_allele,
        substr(original_seq, cds_mutation_loc + nchar(ref_allele), nchar(original_seq))
      )
      
      # Find first in-frame stop
      ptc_pos <- find_ptc(mutated_seq, cut = cds_mutation_loc)
      
      # Distance in codons
      if (!is.na(ptc_pos)) {
        mut_codon_start <- cds_mutation_loc - ((cds_mutation_loc - 1) %% 3)
        distance <- max(0, ptc_pos - mut_codon_start)
      } else {
        distance <- NA_real_
      }
      
      # --- FASTA writing for frameshifts ---
      if (nchar(ref_allele) != nchar(alt_allele)) {
        truncated_seq <- if (!is.na(ptc_pos)) substr(mutated_seq, 1, ptc_pos + 2) else mutated_seq
        protein_seq <- translate_sequence(truncated_seq)
        
        fasta_filename <- file.path(output_dir, paste0(gsub("[:|]", "_", key), ".fasta"))
        sink(fasta_filename)
        cat(">", key, "\n", protein_seq, "\n", sep = "")
        sink()
      }
      
      distance_list <- rbind(
        distance_list,
        data.frame(
          Variant_Key      = key,
          strand           = strand,
          cds_mutation_loc = cds_mutation_loc,
          ptc_pos          = ptc_pos,
          Distance         = distance,
          stringsAsFactors = FALSE
        )
      )
      
    }, error = function(e) {
      cat("  Error at index", i, "-", variants$key[i], ":", conditionMessage(e), "\n")
    })
  }
  
  return(distance_list)
}

esembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
minus1_control_variants$key = minus1_control_variants$id
plus1_control_variants$key = plus1_control_variants$id
snv_control_variants$key = snv_control_variants$id

plus1_dis = create_fasta(plus1_variants)
plus1_control_dis = create_fasta(plus1_control_variants, output_dir = "plus1_control_fasta_output")
minus1_dis = create_fasta(minus1_variants, output_dir = "minus1_test_fasta_output")
minus1_control_dis = create_fasta(minus1_control_variants, output_dir = "minus1_control_fasta_output")
snv_dis <- create_fasta(snv_variants,output_dir = "snv_fasta_output")
snv_control_dis = create_fasta(snv_control_variants, output_dir = "snv_control_fasta_output")