#This file shows the variable distribution in different dataframes

#uniprot id
snv_variants2 = read.csv('snv_variants20260201_plp_clinvar.csv')
snv_variants2$key
snv_variants2$uniprotswissprot
snv_variants2$transcript
fs_variants2 = read.csv('fs_variants20260201_plp_clinvar.csv')
gnomad_snv_filtered2
gnomad_fs_filtered2

#ptc loc
fs_dis
snv_dis$Variant_Key
snv_dis$ptc_pos
snv_dis$cds_mutation_loc
gnomad_fs_dis
gnomad_snv_dis

#nmdesc loc: the genes are the same for clinvar and gnomad since we use the same genes for analysis
snv_nmdesc_df$can_region_start
snv_nmdesc_df$can_region_end
snv_nmdesc_df$snv_nmdesc_cds
snv_nmdesc_df$coding
snv_nmdesc_df$ensembl_transcript_id
fs_nmdesc_df


