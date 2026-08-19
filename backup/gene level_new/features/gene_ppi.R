#This R script is to plot percentage of genes that appear in ppi network
library(data.table)
human_1_ <- fread("~/Downloads/human (1) (2).txt", data.table = FALSE)

#this function is to convert string like "[1, 2, 3]" to numeric vector c(1, 2, 3)
convert_to_c <- function(x) {
  x <- gsub("\\[|\\]", "", x) # Remove square brackets
  if (x == "") return(c())     # Handle empty brackets
  as.numeric(trimws(unlist(strsplit(x, ",")))) # Split by comma, trim whitespace, and convert to numeric
}

gene_all$matched_uniprot <- 0

for (i in seq_len(nrow(gene_all))) {
  
  uid <- gene_all$uniprot[i]
  aft_NMD_ind <- gene_all$NMD_region_start[i]
  
  if (is.na(uid) || uid == "") next
  if (is.na(aft_NMD_ind)) next
  
  re_1 <- unlist(lapply(
    human_1_$interface_residues1[human_1_$uniprot1 == uid],
    convert_to_c
  )) * 3 - 2
  
  re_2 <- unlist(lapply(
    human_1_$interface_residues2[human_1_$uniprot2 == uid],
    convert_to_c
  )) * 3 - 2
  
  if (length(re_1) == 0 && length(re_2) == 0) next
  
  if (any(re_1 >= aft_NMD_ind, na.rm = TRUE) ||
      any(re_2 >= aft_NMD_ind, na.rm = TRUE)) {
    gene_all$matched_uniprot[i] <- 1
  }
}

gene_all %>% group_by(group) %>%
  summarise(
    total_genes = n(),
    matched_genes = sum(matched_uniprot),
    percentage = matched_genes / total_genes * 100
  )

#rename match_uniprot to ppi_overlap
gene_all <- gene_all %>%
  rename(ppi_overlap = matched_uniprot)
