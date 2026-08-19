library(stringr)
library(biomaRt)

get_re_loc <- function(location, exon_data, strand) {
  exon_data <- exon_data[order(exon_data$cds_start, decreasing = FALSE), ]
  
  cumulative_cds_length <- 0
  for (i in 1:nrow(exon_data)) {
    if (is.na(exon_data$cds_start[i]) || is.na(exon_data$cds_end[i])) next
    
    if (location >= exon_data$exon_chrom_start[i] && location <= exon_data$exon_chrom_end[i]) {
      if (strand == -1) {
        offset_within_exon <- exon_data$exon_chrom_end[i] - location + 1
      } else {
        offset_within_exon <- location - exon_data$exon_chrom_start[i] + 1
      }
      return(cumulative_cds_length + offset_within_exon)
    } else {
      cumulative_cds_length <- cumulative_cds_length + (exon_data$cds_end[i] - exon_data$cds_start[i] + 1)
    }
  }
  return(NA)
}

reverse_complement <- function(sequence) {
  chars <- strsplit(sequence, NULL)[[1]]
  chars <- rev(chars)
  complement_map <- c("A"="T","T"="A","G"="C","C"="G",
                      "a"="t","t"="a","g"="c","c"="g")
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
    if (codon %in% stop_codons && (i + 2) >= cut) return(i) #we are forcing the code to only find PTC after the mutation location
  }
  return(NA)
}

find_ptc2 <- function(sequence, cut) {
  stop_codons <- c("TAA","TAG","TGA")
  for (i in seq(1, nchar(sequence) - 2, by = 3)) {
    codon <- substr(sequence, i, i + 2)
    if (codon %in% stop_codons){
      if((i + 2) < cut){
        print(i)
      }  #see if the first PTC found is before the mutation location or not 
      return(i)
    }  
  }
  return(NA)
}

#should we look at overall GC content or just the GC in NMDesc region?
#both compare difference
get_gc_content <- function(sequence) {
  bases <- strsplit(toupper(sequence), "")[[1]]
  gc <- sum(bases %in% c("G", "C"))
  return(round(gc / length(bases) * 100, 2))
}

get_repeat_content <- function(sequence) {
  sequence <- toupper(sequence)
  seq_len  <- nchar(sequence)
  
  # ---------------------------------------------------------------
  # Core detection: mirrors Python regex '([ATGC]{1,6}?)\1+'
  # with finditer(..., overlapped=True)
  # In R, we use gregexpr with perl=TRUE for each unit length 1-6,
  # scanning from every position to capture overlapping matches
  # ---------------------------------------------------------------
  
  motif_hits <- list()
  
  for (unit_len in 1:6) {
    # Build pattern equivalent to ([ATGC]{unit_len})\1+
    pattern <- paste0("([ATGC]{", unit_len, "})\\1+")
    
    # gregexpr with perl=TRUE finds non-overlapping matches natively;
    # to mirror overlapped=TRUE, we shift the start position by 1
    # and re-search, collecting all unique hit coordinates
    
    for (offset in 0:(unit_len - 1)) {
      # Prepend 'offset' dummy chars so matches shift by offset bases
      padded_seq <- paste0(strrep("N", offset), sequence)
      
      m <- gregexpr(pattern, padded_seq, perl = TRUE)[[1]]
      if (m[1] == -1) next
      
      starts     <- as.integer(m) - offset          # correct for padding
      lengths    <- attr(m, "match.length")
      
      for (k in seq_along(starts)) {
        s <- starts[k]
        e <- s + lengths[k] - 1
        if (s < 1 || e > seq_len) next              # discard out-of-bounds
        
        tract     <- substr(sequence, s, e)
        unit_motif <- substr(sequence, s, s + unit_len - 1)
        if (!grepl("^[ATGC]+$", unit_motif)) next   # skip N-containing motifs
        
        motif_hits[[length(motif_hits) + 1]] <- list(
          motif        = unit_motif,
          unit_len     = unit_len,
          repeat_count = lengths[k] %/% unit_len,
          start        = s,
          end          = e,
          tract_len    = lengths[k]
        )
      }
    }
  }
  
  # ---------------------------------------------------------------
  # Empty result guard
  # ---------------------------------------------------------------
  if (length(motif_hits) == 0) {
    return(list(
      total_repeat_bp       = 0,
      repeat_fraction       = 0,
      n_repeat_tracts       = 0,
      longest_tract         = 0,
      homopolymer_fraction  = 0,
      dinucleotide_fraction = 0,
      dominant_motif        = NA_character_,
      repeat_length_dist    = integer(0)
    ))
  }
  
  hits_df <- do.call(rbind, lapply(motif_hits, as.data.frame))
  hits_df <- hits_df[!duplicated(hits_df[, c("start","end","unit_len")]), ]  # remove offset duplicates
  
  # ---------------------------------------------------------------
  # Per-category base counts (before interval merging, for fractions)
  # ---------------------------------------------------------------
  hp_bp <- sum(hits_df$tract_len[hits_df$unit_len == 1], na.rm = TRUE)
  di_bp <- sum(hits_df$tract_len[hits_df$unit_len == 2], na.rm = TRUE)
  
  # ---------------------------------------------------------------
  # Merge overlapping intervals ??? unique repeat-covered bases
  # (same logic as the paper's DRL: each base counted once)
  # ---------------------------------------------------------------
  hits_sorted    <- hits_df[order(hits_df$start), ]
  merged_end     <- -1L
  total_repeat_bp <- 0L
  
  for (j in seq_len(nrow(hits_sorted))) {
    s <- hits_sorted$start[j]
    e <- hits_sorted$end[j]
    if (s > merged_end) {
      total_repeat_bp <- total_repeat_bp + (e - s + 1L)
      merged_end <- e
    } else if (e > merged_end) {
      total_repeat_bp <- total_repeat_bp + (e - merged_end)
      merged_end <- e
    }
  }
  
  # ---------------------------------------------------------------
  # Repeat length distribution (DRL) ??? tract lengths as a vector,
  # mirroring what the paper measures per genome segment.
  # Bootstrap CIs can be computed on this vector downstream.
  # ---------------------------------------------------------------
  repeat_length_dist <- sort(hits_df$tract_len)
  
  # ---------------------------------------------------------------
  # Dominant motif
  # ---------------------------------------------------------------
  motif_totals   <- tapply(hits_df$tract_len, hits_df$motif, sum)
  dominant_motif <- names(which.max(motif_totals))
  
  return(list(
    total_repeat_bp       = total_repeat_bp,
    repeat_fraction       = round(total_repeat_bp / seq_len * 100, 2),
    n_repeat_tracts       = nrow(hits_df),
    longest_tract         = max(hits_df$tract_len),
    homopolymer_fraction  = round(hp_bp / seq_len * 100, 2),
    dinucleotide_fraction = round(di_bp / seq_len * 100, 2),
    dominant_motif        = dominant_motif,
    repeat_length_dist    = repeat_length_dist   # DRL vector for bootstrapping
  ))
}

bootstrap_repeat_ci <- function(repeat_length_dist, n_boot = 1000, ci = 0.95) {
  # repeat_length_dist: the $repeat_length_dist vector from get_repeat_content()
  # Mirrors: randomly sample tract lengths with replacement, repeat 1000x,
  # then derive 95% CI by removing top/bottom 25 counts per length bin
  
  if (length(repeat_length_dist) == 0) return(NULL)
  
  max_len    <- max(repeat_length_dist)
  boot_dists <- matrix(0L, nrow = n_boot, ncol = max_len)
  
  for (b in 1:n_boot) {
    sampled <- sample(repeat_length_dist, replace = TRUE)
    counts  <- tabulate(sampled, nbins = max_len)
    boot_dists[b, ] <- counts
  }
  
  # Remove top and bottom 25 counts per bin, then take min/max ??? 95% CI
  trim    <- 25
  ci_low  <- apply(boot_dists, 2, function(x) min(sort(x)[(trim + 1):(n_boot - trim)]))
  ci_high <- apply(boot_dists, 2, function(x) max(sort(x)[(trim + 1):(n_boot - trim)]))
  
  data.frame(
    tract_length = 1:max_len,
    ci_low       = ci_low,
    ci_high      = ci_high
  )
}

create_fasta <- function(variants, output_dir = "fasta_output") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  distance_list <- data.frame(
    Variant_Key           = character(),
    strand                = numeric(),
    cds_mutation_loc      = numeric(),
    ptc_pos               = numeric(),
    Distance              = numeric(),
    GC_Content            = numeric(),
    repeat_fraction       = numeric(),
    homopolymer_fraction  = numeric(),
    dinucleotide_fraction = numeric(),
    n_repeat_tracts       = numeric(),
    longest_tract         = numeric(),
    dominant_motif        = character(),
    stringsAsFactors      = FALSE
  )
  
  transcript_ids <- variants$transcript
  
  #get CDS sequences
  cds_seq_df <- getBM(
    attributes = c("ensembl_transcript_id", "coding"),
    filters    = "ensembl_transcript_id",
    values     = transcript_ids,
    mart       = ensembl
  )
  cds_seq_df  <- cds_seq_df[!is.na(cds_seq_df$coding) & cds_seq_df$coding != "", ]
  cds_seq_map <- setNames(cds_seq_df$coding, cds_seq_df$ensembl_transcript_id)
  
  #get exon data
  exon_data_df <- getBM(
    attributes = c("ensembl_transcript_id", "chromosome_name", "exon_chrom_start",
                   "cds_start", "cds_end", "exon_chrom_end", "rank", "strand"),
    filters    = "ensembl_transcript_id",
    values     = transcript_ids,
    mart       = ensembl
  )
  exon_data_df <- exon_data_df[complete.cases(exon_data_df[, c("ensembl_transcript_id",
                                                               "exon_chrom_start","exon_chrom_end",
                                                               "cds_start","cds_end","strand")]), ]
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
      
      # Write FASTA for all variant types (SNV + frameshift)
      truncated_seq <- if (!is.na(ptc_pos)) substr(mutated_seq, 1, ptc_pos + 2) else mutated_seq
      
      #get GC content
      gc_content <- get_gc_content(truncated_seq)
      
      #get NMDesc region length
      ##for snv
      NMDesc_length = nchar(original_seq) - ptc_pos + 1
      ##for fs
      NMDesc_length = PTC_info[which(PTC_info$transcript == transcript_id),'can_region']
      
      #get repeat / homopolymer content 
      repeat_info <- get_repeat_content(truncated_seq)
      
      protein_seq   <- translate_sequence(truncated_seq)
      
      fasta_filename <- file.path(output_dir, paste0(gsub("[:|]", "_", key), ".fasta"))
      sink(fasta_filename)
      cat(">", key, "\n", protein_seq, "\n", sep = "")
      sink()
      
      distance_list <- rbind(
        distance_list,
        data.frame(
          Variant_Key           = key,
          strand                = strand,
          cds_mutation_loc      = cds_mutation_loc,
          ptc_pos               = ptc_pos,
          Distance              = distance,
          GC_Content            = gc_content,
          repeat_fraction       = repeat_info$repeat_fraction,
          homopolymer_fraction  = repeat_info$homopolymer_fraction,
          dinucleotide_fraction = repeat_info$dinucleotide_fraction,
          n_repeat_tracts       = repeat_info$n_repeat_tracts,
          longest_tract         = repeat_info$longest_tract,
          dominant_motif        = repeat_info$dominant_motif,
          stringsAsFactors      = FALSE
        )
      )
      
    }, error = function(e) {
      cat("  Error at index", i, "-", variants$key[i], ":", conditionMessage(e), "\n")
    })
  }
  
  return(distance_list)
}

