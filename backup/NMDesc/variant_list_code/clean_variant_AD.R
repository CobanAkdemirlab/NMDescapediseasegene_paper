library(readr)

#clean variant list

#emel clean
WT_vars_full_length_corrected_10_ <- read_csv("~/Downloads/WT-vars_full-length_corrected[10].csv")
WT_vars_NMD_region_corrected_65_ <- read_csv("~/Downloads/WT-vars_NMD_region_corrected[65].csv")

#+1 strand
WT_vars_full_length_plus1_strand_16_ <- read_csv("~/Downloads/WT_vars_full_length_plus1_strand[16].csv")
WT_vars_NMD_region_plus1_strand_53_ <- read_csv("~/Downloads/WT_vars_NMD_region_plus1_strand[53].csv")

WT_var_NMD_2c_key = intersect(WT_vars_NMD_region_corrected_65_$id, WT_vars_NMD_region_plus1_strand_53_$id)
WT_var_full_2c_key = intersect(WT_vars_full_length_corrected_10_$id, WT_vars_full_length_plus1_strand_16_$id)

WT_var_NMD_2c       <- WT_vars_NMD_region_corrected_65_[
  WT_vars_NMD_region_corrected_65_$id %in% WT_var_NMD_2c_key, ]
WT_var_full_2c       <- WT_vars_full_length_corrected_10_[
  WT_vars_full_length_corrected_10_$id %in% WT_var_full_2c_key, ]

#AD variants
WT_var_NMD_2c <- WT_var_NMD_2c %>%
  tidyr::separate(id, into = c("chr","pos","ref","alt"), sep = "_", remove = FALSE) %>%
  mutate(pos = as.integer(pos),
         key = paste0(chr, ":", sprintf("%09d", pos), "|", ref, "|", alt))
WT_var_NMD_3c <- WT_var_NMD_2c %>%
  filter(key %in% ppi_all_AD$Variant_Key)

WT_var_full_2c <- WT_var_full_2c %>%
  tidyr::separate(id, into = c("chr","pos","ref","alt"), sep = "_", remove = FALSE) %>%
  mutate(pos = as.integer(pos),
         key = paste0(chr, ":", sprintf("%09d", pos), "|", ref, "|", alt))
WT_var_full_3c <- WT_var_full_2c %>%
  filter(key %in% ppi_all_AD$Variant_Key)

#unique
WT_var_NMD_4c <- WT_var_NMD_3c[!duplicated(WT_var_NMD_3c$key), ]
WT_var_full_4c <- WT_var_full_3c[!duplicated(WT_var_full_3c$key), ]

write_csv(WT_var_NMD_4c, file = "~/Downloads/WT_var_NMD_4c_AD.csv")
write_csv(WT_var_full_4c, file = "~/Downloads/WT_var_full_4c_AD.csv")