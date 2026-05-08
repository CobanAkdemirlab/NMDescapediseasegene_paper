##exclude any gene genes with any frameshift NMD_can_esc mutation
snv_gene2remove.ind = which(snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)
snv_gene2remove = unique(snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][gene2remove.ind])
#get hgnc_symbol of gene2remove using getBM
snv_gene2remove = getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = snv_gene2remove,
  mart = mart
)
#remove genes2remove from snv genes with plp ptc variants in clinvar
snv_gene_list = unique(res_snv_plp_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]])

snv_genes_remain_in_clinvar = snv_gene_list[-which(snv_gene_list %in% snv_gene2remove$ensembl_transcript_id)]

#only keep the genes that passed the variant_summary filter
v_ch = read.csv('v_ch20260201.csv')
snv_genes_remain_in_clinvar2 = snv_genes_remain_in_clinvar[which(snv_genes_remain_in_clinvar %in% v_ch$tx_id2)]

length(unique(snv_genes_remain_in_clinvar2)) 

#filter for AD genes
omim_AD_symbols_tx = getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol",'transcript_is_canonical'),
  filters = "hgnc_symbol",
  values = omim_AD_symbols,
  mart = mart
)
#keep canonical transcript only
omim_AD_symbols_tx = omim_AD_symbols_tx[which(omim_AD_symbols_tx$transcript_is_canonical == 1),]
snv_control_gene = unique(snv_genes_remain_in_clinvar2)
snv_control_gene_AD = snv_control_gene[which(snv_control_gene %in% omim_AD_symbols_tx$ensembl_transcript_id)]
write.csv(unique(snv_control_gene_AD), 'snv_control_genes_AD.csv', row.names = FALSE)
write.csv(unique(snv_genes_remain_in_clinvar2), 'snv_control_genes.csv', row.names = FALSE)