
parse_position <- function(key_col) {
  as.integer(gsub("^[^:]+:(0*)(\\d+)\\|.*", "\\2", key_col))
}

snv_variants <- snv_variants %>%
  mutate(genomic_loc = parse_position(key))

gnomad_snv_filtered <- gnomad_snv_filtered %>%
  mutate(genomic_loc = parse_position(key))

gnomad_fs_filtered <- gnomad_fs_filtered %>%
  mutate(genomic_loc = parse_position(id))

fs_variants2 <- fs_variants2 %>%
  mutate(genomic_loc = parse_position(key))

all_transcripts <- unique(c(
  snv_variants$transcript,
  gnomad_snv_filtered$transcript,
  gnomad_fs_filtered$transcript,
  fs_variants2$transcript
))

exon_info <- getBM(
  attributes = c(
    "ensembl_transcript_id",
    "strand",
    "exon_chrom_start",
    "exon_chrom_end",
    "genomic_coding_start",
    "genomic_coding_end",
    "cds_start",        # CDS内的累计起始位置
    "cds_end",          # CDS内的累计终止位置
    "rank"              # 外显子顺序
  ),
  filters = "ensembl_transcript_id",
  values  = unique(all_transcripts),
  mart    = mart
) %>%
  # 只保留有CDS注释的外显子（非UTR）
  filter(!is.na(cds_start), !is.na(cds_end)) %>%
  arrange(ensembl_transcript_id, rank)

cds_lengths <- exon_info %>%
  group_by(ensembl_transcript_id) %>%
  summarise(
    cds_total_length = max(cds_end),
    .groups = "drop"
  )

get_re_loc <- function(location, transcript_match) {
  df <- exon_info[exon_info$ensembl_transcript_id == transcript_match, ]
  if (nrow(df) == 0) return(NA_real_)
  
  is_minus <- unique(df$strand) == -1
  
  for (i in seq_len(nrow(df))) {
    start <- df$genomic_coding_start[i]
    end   <- df$genomic_coding_end[i]
    
    if (location >= min(start, end) && location <= max(start, end)) {
      if (is_minus) {
        offset <- max(start, end) - location
      } else {
        offset <- location - min(start, end)
      }
      return(df$cds_start[i] + offset)
    }
  }
  return(NA_real_)
}

add_cds_dist <- function(df, label) {
  
  df$cds_loc <- vapply(
    seq_len(nrow(df)),
    function(i) {
      loc <- df$genomic_loc[i]
      tx  <- df$transcript[i]
      
      if (is.na(loc) || is.na(tx) || tx == "") return(NA_real_)
      
      tryCatch(
        get_re_loc(loc, tx),
        error = function(e) NA_real_
      )
    },
    FUN.VALUE = numeric(1)
  )
  
  df <- df %>%
    left_join(cds_lengths,
              by = c("transcript" = "ensembl_transcript_id")) %>%
    mutate(
      dist_to_cds_end = cds_total_length - cds_loc,
      source          = label
    )
  
  return(df)
}

snv_variants_out     <- add_cds_dist(snv_variants,        "snv_variants")
gnomad_snv_out       <- add_cds_dist(gnomad_snv_filtered, "gnomad_snv")
gnomad_fs_out        <- add_cds_dist(gnomad_fs_filtered,  "gnomad_fs")
fs_variants2_out     <- add_cds_dist(fs_variants2,        "fs_variants2")

all_results <- bind_rows(
  snv_variants_out,
  gnomad_snv_out,
  gnomad_fs_out,
  fs_variants2_out
)

all_results %>%
  group_by(source) %>%
  summarise(
    n           = n(),
    mean_dist   = mean(dist_to_cds_end, na.rm = TRUE),
    median_dist = median(dist_to_cds_end, na.rm = TRUE),
    n_NA        = sum(is.na(dist_to_cds_end)),
    n_negative  = sum(dist_to_cds_end < 0, na.rm = TRUE)
  ) 


unique_tx <- c("ENST00000295754", "ENST00000231790", 
               "ENST00000423572", "ENST00000639785")

strand_info <- getBM(
  attributes = c("ensembl_transcript_id", "strand", 
                 "hgnc_symbol", "chromosome_name"),
  filters    = "ensembl_transcript_id",
  values     = unique_tx,
  mart       = mart
)



getSequence(
  id      = "ENST00000263253",
  type    = "ensembl_transcript_id",
  seqType = "peptide",
  mart    = mart
)

