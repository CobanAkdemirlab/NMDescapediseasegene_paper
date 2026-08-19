library(biomaRt)
library(dplyr)
get_frameshift_type2 <- function(key) {
  parts <- str_split(key, "\\|")[[1]]
  ref <- parts[2]
  alt <- parts[3]
  
  if (nchar(ref) > nchar(alt)) {
    type <- 'del'
  } else if (nchar(alt) > nchar(ref)) {
    type <- 'ins'
  } else {
    return('error')  # no size difference, not indel
  }
  
  yu <- (abs(nchar(ref) - nchar(alt))) %% 3
  
  if ((type == 'del' & yu == 1) | (type == 'ins' & yu == 2)) {
    return('minus1')
  } else if ((type == 'del' & yu == 2) | (type == 'ins' & yu == 1)) {
    return('plus1')
  } else {
    return('error')
  }
}
fs_gnomAD_variants$fs_type <- mapply(
  get_frameshift_type2,
  fs_gnomAD_variants$id
)

# Subset by plus1 and minus1
plus1_gnomAD_variants <- fs_gnomAD_variants %>%
  filter(fs_type == 'plus1')

minus1_gnomAD_variants <- fs_gnomAD_variants %>%
  filter(fs_type == 'minus1')
# PLUS1 CONTROL
plus1_control_gene <- plus1_control_pli$gene

plus1_control_tr <- getBM(
  attributes = c('ensembl_transcript_id', 'hgnc_symbol', 'uniprotswissprot', 'transcript_is_canonical'),
  filters = 'hgnc_symbol',
  values = plus1_control_gene,
  mart = ensembl
) %>%
  filter(transcript_is_canonical == 1)

plus1_control_tr_ids <- plus1_control_tr$ensembl_transcript_id

plus1_control_variants <- plus1_gnomAD_variants %>%
  filter(transcript %in% plus1_control_tr_ids) %>%
  left_join(
    plus1_control_tr %>% select(ensembl_transcript_id, uniprotswissprot),
    by = c("transcript" = "ensembl_transcript_id")
  )

write.csv(plus1_control_variants, 'plus1_control_gnomAD_variants0406.csv', row.names = FALSE)


# MINUS1 CONTROL
minus1_control_gene <- minus1_control_pli$gene

minus1_control_tr <- getBM(
  attributes = c('ensembl_transcript_id', 'hgnc_symbol', 'uniprotswissprot', 'transcript_is_canonical'),
  filters = 'hgnc_symbol',
  values = minus1_control_gene,
  mart = ensembl
) %>%
  filter(transcript_is_canonical == 1)

minus1_control_tr_ids <- minus1_control_tr$ensembl_transcript_id

minus1_control_variants <- minus1_gnomAD_variants %>%
  filter(transcript %in% minus1_control_tr_ids) %>%
  left_join(
    minus1_control_tr %>% select(ensembl_transcript_id, uniprotswissprot),
    by = c("transcript" = "ensembl_transcript_id")
  )

write.csv(minus1_control_variants, 'minus1_control_gnomAD_variants0406.csv', row.names = FALSE)
