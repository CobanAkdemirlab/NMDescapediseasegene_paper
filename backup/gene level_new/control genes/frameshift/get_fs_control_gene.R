library(stringr)

##exclude any gene with any plp frameshift NMD_can_esc mutation
fs_res = read_rds('fs.rds'). #this file is already plp and ptc
fs_res_can = fs_res[which(fs_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |fs_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)]

fs_gene2remove.ind = which(fs_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |fs_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)
fs_gene2remove = unique(fs_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][gene2remove.ind])
#get hgnc_symbol of gene2remove using getBM
fs_gene2remove$hgnc = getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = fs_gene2remove,
  mart = mart
)
#remove genes2remove from gene_list
fs_gene_list = unique(fs@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]])
fs_genes_remain_in_clinvar = fs_gene_list[-which(fs_gene_list %in% fs_gene2remove)]
length(unique(fs_genes_remain_in_clinvar)) 

#keep canonical transcripts only using getBM
fs_can.info = getBM(
  attributes = c("ensembl_transcript_id","transcript_is_canonical"), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values =  fs_genes_remain_in_clinvar,                          # Transcript ID
  mart = ensembl                                         # Database connection
)

fs_canonical_transcripts <- fs_can.info$ensembl_transcript_id[
  fs_can.info$transcript_is_canonical == 1
]

fs_genes_remain_in_clinvar <- fs_genes_remain_in_clinvar[
  fs_genes_remain_in_clinvar %in% fs_canonical_transcripts
]

length(unique(fs_genes_remain_in_clinvar))


#remove single exon genes using getBM
fs_BM.infoo <- getBM(
  attributes = c("ensembl_transcript_id","rank",'cds_start','cds_end','exon_chrom_start','exon_chrom_end'), # Attributes to retrieve
  filters = "ensembl_transcript_id",                  
  values =  fs_genes_remain_in_clinvar,                
  mart = ensembl                                         
)

fs_exon_counts <- aggregate(rank ~ ensembl_transcript_id, data = fs_BM.infoo, FUN = max)
colnames(fs_exon_counts) <- c("ensembl_transcript_id", "max_exon_rank")

# Keep only multi-exon transcripts (max rank > 1)
fs_multi_exon_transcripts <- fs_exon_counts$ensembl_transcript_id[
  fs_exon_counts$max_exon_rank > 1
]

fs_genes_remain_in_clinvar <- fs_genes_remain_in_clinvar[
  fs_genes_remain_in_clinvar %in% fs_multi_exon_transcripts
]

length(unique(fs_genes_remain_in_clinvar))


#get their hgnc_symbol
fs_genes_remain_in_clinvar_hgnc = getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = fs_genes_remain_in_clinvar,
  mart = mart
)
length(unique(fs_genes_remain_in_clinvar_hgnc$hgnc_symbol))

write.csv(unique(fs_genes_remain_in_clinvar_hgnc$hgnc_symbol),'fs_control_gene.csv')

#filter for AD genes
fs_control_gene = unique(fs_genes_remain_in_clinvar_hgnc$ensembl_transcript_id)
fs_control_gene_AD = fs_control_gene[which(fs_control_gene %in% omim_AD_symbols_tx$ensembl_transcript_id)]
write.csv(unique(fs_control_gene_AD), 'fs_control_genes_AD.csv', row.names = FALSE)