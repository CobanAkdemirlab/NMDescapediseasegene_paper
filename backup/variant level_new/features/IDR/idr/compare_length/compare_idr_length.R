library(jsonlite)

# Safe parsing
parse_json_boundaries <- function(x) {
  if (is.na(x) || x == "No Data" || trimws(x) == "") return(NULL)
  tryCatch({
    fromJSON(x)
  }, error = function(e) {
    return(NULL)
  })
}

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
snv_idr$wild_boundaries_list <- lapply(snv_idr$wild_Disordered_Domain_Boundaries, parse_json_boundaries)
snv_idr$wild_IDR_total_length <- sapply(snv_idr$wild_boundaries_list, compute_idr_length)

# Disorder boundaries
snv_idr$disorder_boundaries_list <- lapply(snv_idr$Disorder_Domains, parse_json_boundaries)
snv_idr$mut_disorder_IDR_total_length <- sapply(snv_idr$disorder_boundaries_list, compute_idr_length)

#compare mut with wild idr total length
snv_idr$IDR_length_difference <- snv_idr$mut_disorder_IDR_total_length - snv_idr$wild_IDR_total_length
snv_mean_difference <- mean(snv_idr$IDR_length_difference, na.rm = TRUE)
snv_median_difference <- median(snv_idr$IDR_length_difference, na.rm = TRUE)
