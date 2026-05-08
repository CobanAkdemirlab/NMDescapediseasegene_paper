#do motif
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
})

# base->aa index helper
aa_index <- function(base_idx) ((base_idx - 1L) %/% 3L) + 1L

# Build an index of FASTA files and attach aa_pos from the distance_list
index_fasta_with_positions <- function(mut_dir, dist_tbl) {
  files <- list.files(mut_dir, pattern = "\\.fasta$", full.names = TRUE)
  if (!length(files)) return(tibble(path = character(), Variant_Key = character(), aa_pos = integer()))
  
  tibble(path = files) %>%
    mutate(Variant_Key = str_remove(basename(path), "\\.fasta$")) %>%
    left_join(dist_tbl %>% select(Variant_Key, cds_mutation_loc), by = "Variant_Key") %>%
    mutate(aa_pos = ifelse(!is.na(cds_mutation_loc), aa_index(cds_mutation_loc), NA_integer_))
}

# Batch plot for a create_fasta() output directory
plot_batch_from_create_fasta <- function(mut_dir, dist_tbl, window = 60, out_dir = file.path("plots", basename(mut_dir))) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  idx <- index_fasta_with_positions(mut_dir, dist_tbl)
  
  if (!nrow(idx)) {
    message("No FASTA files found in: ", mut_dir)
    return(invisible(list()))
  }
  
  res <- vector("list", nrow(idx))
  names(res) <- idx$Variant_Key
  
  for (i in seq_len(nrow(idx))) {
    mut_fa <- idx$path[i]
    aa_pos <- idx$aa_pos[i]
    key    <- idx$Variant_Key[i]
    
    out_base <- file.path(out_dir, key)
    p <- plot_variant_level_motifs(
      mut_fasta = mut_fa,
      aa_pos    = if (!is.na(aa_pos)) aa_pos else NULL,
      window    = window,
      out_file  = out_base
    )
    res[[i]] <- p
  }
  invisible(res)
}

# ==== Run for each group (using your create_fasta outputs) ====

# plus1
plots_plus1 <- plot_batch_from_create_fasta(
  mut_dir  = "plus1_fasta_output",
  dist_tbl = plus1_dis,
  window   = 60,
  out_dir  = "plots/plus1"
)

# plus1_control
plots_plus1_ctl <- plot_batch_from_create_fasta(
  mut_dir  = "plus1_control_fasta_output",
  dist_tbl = plus1_control_dis,
  window   = 60,
  out_dir  = "plots/plus1_control"
)

# minus1
plots_minus1 <- plot_batch_from_create_fasta(
  mut_dir  = "minus1_test_fasta_output",
  dist_tbl = minus1_dis,
  window   = 60,
  out_dir  = "plots/minus1"
)

# minus1_control
plots_minus1_ctl <- plot_batch_from_create_fasta(
  mut_dir  = "minus1_control_fasta_output",
  dist_tbl = minus1_control_dis,
  window   = 60,
  out_dir  = "plots/minus1_control"
)

# snv
plots_snv <- plot_batch_from_create_fasta(
  mut_dir  = "snv_fasta_output",
  dist_tbl = snv_dis,
  window   = 60,
  out_dir  = "plots/snv"
)

# snv_control
plots_snv_ctl <- plot_batch_from_create_fasta(
  mut_dir  = "snv_control_fasta_output",
  dist_tbl = snv_control_dis,
  window   = 60,
  out_dir  = "plots/snv_control"
)

message("All plots saved under ./plots/<group>/")
