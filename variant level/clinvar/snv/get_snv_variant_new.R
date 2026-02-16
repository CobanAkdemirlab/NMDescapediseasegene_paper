##get snv variants from clinvar

#input gene list
snv_gene = snv_pli$gene
#get canonical transcript of snv_gene$hgnc
snv_tr = getBM(attributes=c('ensembl_transcript_id','hgnc_symbol'),
               filters = 'hgnc_symbol', values = snv_gene, mart = ensembl)
snv_can_tr = getBM(attributes=c('ensembl_transcript_id','transcript_is_canonical'),
                   filters = 'ensembl_transcript_id', values = snv_tr$ensembl_transcript_id, mart = ensembl)
snv_can_tr = snv_can_tr[which(snv_can_tr$transcript_is_canonical == 1),1]

#in clinvar snv ptc, nmdesc variants, select those corresponding to the canonical transcript of gene list
snv = which(res2@elementMetadata@listData[["type"]]=='snv')
#remove single exon
no_single = which(res2@elementMetadata@listData[["res_aenmd"]]@listData[["is_single"]] == F)
plp = grep('pathogenic',res2@elementMetadata@listData[["clinsig"]],ignore.case = T)
ptc = which(res2@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]]==T)
combined_ind = intersect(plp,ptc)
snv_ind = intersect(combined_ind,snv)
snv_ind2 = intersect(snv_ind,no_single)
snv_res = res2[snv_ind2]
#get NMDesc variants by canonical rule
snv_res_can = snv_res[which(snv_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |snv_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)]
snv_variants = data.frame(transcript = snv_res_can@ranges@NAMES[which(snv_res_can@ranges@NAMES %in% snv_can_tr)],
                          key = snv_res_can@elementMetadata@listData[["key"]][which(snv_res_can@ranges@NAMES %in% snv_can_tr)])
#get uniprot id using biomart
unique_transcripts <- unique(snv_variants$transcript)
uniprot_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = unique_transcripts,
  mart = ensembl
)
snv_variants<- merge(
  snv_variants,
  uniprot_mapping,
  by.x = "transcript",
  by.y = "ensembl_transcript_id",
  all.x = TRUE
)


write.csv(snv_variants, 'snv_variants0406.csv', row.names = FALSE)
