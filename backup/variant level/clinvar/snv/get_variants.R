#get snv variants from clinvar
snv_gene = snv_plp_ptc_p1120_NMDenriched2_can$X1
#get canonical transcript of snv_gene$hgnc
snv_tr = getBM(attributes=c('ensembl_transcript_id','hgnc_symbol'),
                     filters = 'hgnc_symbol', values = snv_gene, mart = mart)
snv_can_tr = getBM(attributes=c('ensembl_transcript_id','transcript_is_canonical'),
                     filters = 'ensembl_transcript_id', values = snv_tr$ensembl_transcript_id, mart = mart)
snv_can_tr = snv_can_tr[which(snv_can_tr$transcript_is_canonical == 1),1]
snv_variants = snv_res_can[which(snv_res_can@ranges@NAMES %in% snv_can_tr)]
write.csv(snv_variants@elementMetadata@listData[["key"]], 'snv_variants.csv', row.names = FALSE)


