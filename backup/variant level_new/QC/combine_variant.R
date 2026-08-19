#variants are selected from same genes, but plp ptc from clinvar and (benign) ptc from gnomad, This R script is to combine them into variant_all
library(dplyr)
snv_all <- snv_variants2 %>%
  left_join(
    snv_dis %>%
      select(Variant_Key, ptc_pos, cds_mutation_loc),
    by = c("key" = "Variant_Key")
  ) %>%
  left_join(
    snv_nmdesc_df %>%
      select(
        ensembl_transcript_id,
        can_region_start,
        can_region_end,
        snv_nmdesc_cds,
        coding
      ),
    by = c("transcript" = "ensembl_transcript_id")
  ) %>%
  #can_region_start and can_region_end as NMD_region_start and NMD_region_end
  rename(
    NMD_region_start = can_region_start,
    NMD_region_end = can_region_end,
    cds_ptc_loc = ptc_pos,
    nmdesc_cds = snv_nmdesc_cds,
  ) %>%
  select(-any_of(c("genomic_loc", "PTC_pos"))) # remove the original genomic_pos column since we have renamed it to PTC_genomic_pos
write.csv(snv_all, 'snv_all.csv', row.names = FALSE)

gnomad_snv_all <- gnomad_snv_filtered2 %>%
  left_join(
    gnomad_snv_dis %>%
      select(Variant_Key, ptc_pos, cds_mutation_loc),
    by = c("key" = "Variant_Key")
  ) %>%
  left_join(
    snv_nmdesc_df %>%
      select(
        ensembl_transcript_id,
        can_region_start,
        can_region_end,
        snv_nmdesc_cds,
        coding
      ),
    by = c("transcript" = "ensembl_transcript_id")
  ) %>%
  rename(
    NMD_region_start = can_region_start,
    NMD_region_end = can_region_end,
    mutation_genomic_pos = genomic_loc,
    cds_ptc_loc = ptc_pos,
  ) %>%
  select(-any_of(c("genomic_loc", "PTC_pos"))) # remove the original genomic_pos column since we have renamed it to PTC_genomic_pos
write.csv(gnomad_snv_all, 'gnomad_snv_all.csv', row.names = FALSE)

gnomad_fs_all <- gnomad_fs_filtered2 %>%
  left_join(
    gnomad_fs_dis %>%
      select(Variant_Key, ptc_pos, cds_mutation_loc),
    by = c("key" = "Variant_Key")
  ) %>%
  left_join(
    fs_nmdesc_df %>%
      select(
        ensembl_transcript_id,
        median_can_region_start,
        median_can_region_end,
        fs_nmdesc_cds,
        coding
      ),
    by = c("transcript" = "ensembl_transcript_id")
  ) %>%
  # rename for median_can_region_start and median_can_region_end as NMD_region_start and NMD_region_end
  rename(
    NMD_region_start = median_can_region_start,
    NMD_region_end = median_can_region_end,
    mutation_genomic_pos = genomic_loc,
    cds_ptc_loc = ptc_pos,
  ) %>%
  select(-any_of(c("genomic_loc", "PTC_pos"))) # remove the original genomic_pos column since we have renamed it to PTC_genomic_pos

write.csv(gnomad_fs_all, 'gnomad_fs_all.csv', row.names = FALSE)



fs_all <- fs_variants2 %>%
  left_join(
    fs_dis %>%
      select(Variant_Key, ptc_pos, cds_mutation_loc),
    by = c("key" = "Variant_Key")
  ) %>%
  left_join(
    fs_nmdesc_df %>%
      select(
        ensembl_transcript_id,
        median_can_region_start,
        median_can_region_end,
        fs_nmdesc_cds,
        coding
      ),
    by = c("transcript" = "ensembl_transcript_id")
  ) %>%
  rename(
    NMD_region_start = median_can_region_start,
    NMD_region_end = median_can_region_end,
    cds_ptc_loc = ptc_pos,
    nmdesc_cds = fs_nmdesc_cds
  ) %>%
  select(-any_of(c("genomic_loc", "PTC_pos"))) # remove the original genomic_pos column since we have renamed it to PTC_genomic_pos
write.csv(fs_all, 'fs_all.csv', row.names = FALSE)

gnomad_fs_all2 <- gnomad_fs_all %>%
select(-any_of(c('type','id','chrom','source','mutation_genomic_pos'))) %>%
  rename(nmdesc_cds = fs_nmdesc_cds)
gnomad_snv_all2 <- gnomad_snv_all %>%
  select(-any_of(c('type','id','chrom','source','mutation_genomic_pos'))) %>%
  rename(nmdesc_cds = snv_nmdesc_cds)

variants_all = rbind(snv_all, fs_all,gnomad_snv_all2, gnomad_fs_all2)

variants_all$source = rep(c("snv","fs","snv_control","fs_control"), c(nrow(snv_all), nrow(fs_all), nrow(gnomad_snv_all2), nrow(gnomad_fs_all2)))
variants_all2$cds_end = nchar(variants_all2$coding)
variants_all2$dist_to_cds_end = variants_all2$cds_end - variants_all2$ptc_pos
write.csv(variants_all, 'variants_all.csv', row.names = FALSE)
