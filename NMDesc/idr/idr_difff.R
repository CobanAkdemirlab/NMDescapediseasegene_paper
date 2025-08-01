library(jsonlite)

#create minus1_control_idr
minus1_control_uniprot_ids_out <- read_delim("~/Downloads/idr0512/minus1_control_uniprot_ids_out.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)
minus1_control_uniprot_ids_out$gene = getBM(attributes = c('hgnc_symbol','uniprotswissprot'),
                                   filters = 'uniprotswissprot',
                                   values = minus1_control_uniprot_ids_out$UniProt_ID,
                                   mart = ensembl)$hgnc_symbol
#from minus1_control_idr$Transcript_ID get gene
tr2gene = getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = unique(minus1_control_idr$Transcript_ID),
  mart = ensembl
)
minus1_control_idr <- merge(minus1_control_idr, tr2gene, by.x = "Transcript_ID", by.y = "ensembl_transcript_id", all.x = TRUE)
minus1_control_idr$gene = minus1_control_idr$hgnc_symbol


minus1_control_idr$wild_Disordered_Domain_Boundaries = minus1_control_uniprot_ids_out$Disordered_Domain_Boundaries[match(minus1_control_idr$gene, minus1_control_uniprot_ids_out$gene)]



# Length computation
compute_idr_length <- function(bounds) {
  if (is.null(bounds) || length(bounds) == 0) return(0)
  
  tryCatch({
    if (is.matrix(bounds) && ncol(bounds) == 2) {
      # multiple [start, end] pairs
      return(sum(bounds[,2] - bounds[,1] + 1))
    } else if (is.numeric(bounds) && length(bounds) == 2) {
      # single [start, end] pair returned as a vector
      return(bounds[2] - bounds[1] + 1)
    } else if (is.list(bounds)) {
      # maybe a list of vectors (rare case)
      return(sum(sapply(bounds, function(pair) {
        if (is.numeric(pair) && length(pair) == 2) {
          return(pair[2] - pair[1] + 1)
        } else {
          return(0)
        }
      })))
    } else {
      return(0)
    }
  }, error = function(e) {
    return(0)
  })
}

# Wild boundaries
minus1_control_idr$wild_boundaries_list <- lapply(minus1_control_idr$wild_Disordered_Domain_Boundaries, parse_json_boundaries)
minus1_control_idr$wild_IDR_total_length <- sapply(minus1_control_idr$wild_boundaries_list, compute_idr_length)

# Disorder boundaries
minus1_control_idr$disorder_boundaries_list <- lapply(minus1_control_idr$Disorder_Domains, parse_json_boundaries)
minus1_control_idr$mut_disorder_IDR_total_length <- sapply(minus1_control_idr$disorder_boundaries_list, compute_idr_length)

#compare mut with wild idr total length
minus1_control_idr$IDR_length_difference <- minus1_control_idr$mut_disorder_IDR_total_length - minus1_control_idr$wild_IDR_total_length
minus1_control_mean_difference <- mean(minus1_control_idr$IDR_length_difference, na.rm = TRUE)
minus1_control_median_difference <- median(minus1_control_idr$IDR_length_difference, na.rm = TRUE)