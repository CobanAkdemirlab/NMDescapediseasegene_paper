##get snv variants from clinvar

#input gene list
#snv_gene = read.csv('snv_plp_ptc_nmdesc_can_p_f_syn_20260201_NMDesc_enriched_can_AD_p_0.8.csv', header = T)$V1
snv_gene = wald_AD_genes
#get canonical transcript of snv_gene$hgnc
snv_tr = getBM(attributes=c('ensembl_transcript_id','hgnc_symbol'),
               filters = 'hgnc_symbol', values = snv_gene, mart = ensembl)
snv_can_tr = getBM(attributes=c('ensembl_transcript_id','transcript_is_canonical'),
                   filters = 'ensembl_transcript_id', values = snv_tr$ensembl_transcript_id, mart = ensembl)
snv_can_tr = snv_can_tr[which(snv_can_tr$transcript_is_canonical == 1),1]

#in clinvar snv ptc, nmdesc variants, select those corresponding to the canonical transcript of gene list
snv_res_can = readRDS('snv_plp_ptc_nmdesc_can_filtered20260201.rds')
snv_vus_res_can = readRDS('snv_vus_ptc_nmdesc_can_filtered20260201.rds')
snv_benign_res_can = readRDS('snv_benign_ptc_nmdesc_can_filtered20260201.rds')
#get NMDesc variants 
snv_variants = data.frame(transcript = snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][which(snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% snv_can_tr)],
                          key = snv_res_can@elementMetadata@listData[["key"]][which(snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% snv_can_tr)])
snv_vus_variants = data.frame(transcript = snv_vus_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][which(snv_vus_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% snv_can_tr)],
                          key = snv_vus_res_can@elementMetadata@listData[["key"]][which(snv_vus_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% snv_can_tr)])
snv_benign_variants = data.frame(transcript = snv_benign_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][which(snv_benign_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% snv_can_tr)],
                          key = snv_benign_res_can@elementMetadata@listData[["key"]][which(snv_benign_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% snv_can_tr)])
length(unique(snv_can_tr))
length(unique(snv_variants$transcript))
#get uniprot id using biomart
unique_transcripts <- unique(snv_variants$transcript)
unique_benign_transcripts <- unique(snv_benign_variants$transcript)
uniprot_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = unique_transcripts,
  mart = ensembl
)
benign_uniprot_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = unique_benign_transcripts,
  mart = ensembl
)
snv_variants<- merge(
  snv_variants,
  uniprot_mapping,
  by.x = "transcript",
  by.y = "ensembl_transcript_id",
  all.x = TRUE
)
snv_benign_variants<- merge(
  snv_benign_variants,
  benign_uniprot_mapping,
  by.x = "transcript",
  by.y = "ensembl_transcript_id",
  all.x = TRUE
)
#do the same for benign variants


write.csv(snv_variants, 'snv_variants20260201_plp_wald_clinvar.csv', row.names = FALSE)
write.csv(snv_benign_variants, 'snv_variants20260201_benign_wald_clinvar.csv', row.names = FALSE)