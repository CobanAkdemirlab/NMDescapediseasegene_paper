#This R script is to get repeat / homopolymer content from gene_all and plot it
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
  # Merge overlapping intervals → unique repeat-covered bases
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
  # Repeat length distribution (DRL) — tract lengths as a vector,
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
  
  # Remove top and bottom 25 counts per bin, then take min/max → 95% CI
  trim    <- 25
  ci_low  <- apply(boot_dists, 2, function(x) min(sort(x)[(trim + 1):(n_boot - trim)]))
  ci_high <- apply(boot_dists, 2, function(x) max(sort(x)[(trim + 1):(n_boot - trim)]))
  
  data.frame(
    tract_length = 1:max_len,
    ci_low       = ci_low,
    ci_high      = ci_high
  )
}

detect_repeats <- function(sequence) {
  if (is.na(sequence) || sequence == "") return(0)
  
  pattern <- "([ATGC]{1,6})\\1+"
  
  matches <- str_extract_all(sequence, pattern)[[1]]
  
  if (length(matches) == 0) return(0)
  
  total_repeat_length <- sum(nchar(matches))
  
  return(total_repeat_length)
}

repeat_fraction <- function(sequence) {
  if (is.na(sequence) || sequence == "") return(NA)
  
  rep_len <- detect_repeats(sequence)
  
  return(rep_len / nchar(sequence))
}

detect_homopolymer <- function(sequence) {
  if (is.na(sequence) || sequence == "") return(0)
  
  pattern <- "([ATGC])\\1{3,}"
  
  matches <- str_extract_all(sequence, pattern)[[1]]
  
  if (length(matches) == 0) return(0)
  
  return(sum(nchar(matches)))
}

homopolymer_fraction <- function(sequence) {
  if (is.na(sequence) || sequence == "") return(NA)
  
  hp_len <- detect_homopolymer(sequence)
  
  return(hp_len / nchar(sequence))
}

gene_all$nmdesc_repeat_fraction <- sapply(gene_all$nmdesc_cds, repeat_fraction)
gene_all$repeat_fraction <- sapply(gene_all$coding, repeat_fraction)
gene_all$nmdesc_homopolymer_fraction <- sapply(gene_all$nmdesc_cds, homopolymer_fraction)
gene_all$homopolymer_fraction <- sapply(gene_all$coding, homopolymer_fraction)

#plot by group