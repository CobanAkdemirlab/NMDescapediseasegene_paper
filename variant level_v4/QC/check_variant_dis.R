motif_max_AD$group =  variants_all5_1_$group[match(motif_max_AD$key, variants_all5_1_$Variant_Key)]
table(motif_max_AD$group)


variants_all5_1_ %>%
  group_by(group) %>%
  summarise(
    n_protein = n_distinct(uniprot)
  )

#for variants_all5_1_, only keep the gene in gene_motif_flags
table(gene_motif_flags$group)
snv_part = variants_all5_1_ %>%
  #in each group, only keep the uniprot that also appeard in this group
  filter(uniprot %in% gene_motif_flags$uniprot[which(gene_motif_flags$group %in% c('snv_control','snv'))]) 
fs_part = variants_all5_1_ %>%
  #in each group, only keep the uniprot that also appeard in this group
  filter(uniprot %in% gene_motif_flags$uniprot[which(gene_motif_flags$group %in% c('fs_control','fs'))]) 
snv_part2 = snv_part[which(snv_part$group %in% c('snv_control','snv_disease')),]
fs_part2 = fs_part[which(fs_part$group %in% c('fs_control','fs_disease')),]
variants_all5_1_0621 = rbind(snv_part2, fs_part2)
table(variants_all5_1_0621$group)
write.csv(variants_all5_1_0621, "variants_all5_1_0621.csv", row.names = FALSE)

u1 <- variants_all5_1_ %>%
  filter(group == "snv_disease") %>%
  pull(uniprot) %>%
  unique()

u2 <- gene_motif_flags %>%
  filter(group == "snv") %>%
  pull(uniprot) %>%
  unique()

length(intersect(u1, u2))  # ??????
length(setdiff(u1, u2))    # ??????variants_all5_1_
length(setdiff(u2, u1))    # ??????gene_motif_flags