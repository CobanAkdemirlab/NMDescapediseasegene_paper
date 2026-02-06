library(biomaRt)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
variant_count = 0
gene_count = 0
names = paste('chr',seq(1,22),'.rds',sep='')
fs_NMD_df <- data.frame(transcript = character(), id = character(), stringsAsFactors = FALSE)
for(outfilename in names){
  l2 <- readRDS(outfilename)
  #filter only fs variant in l
  #l2 = unlist(l)
  ins.ind = which(l2@elementMetadata@listData[["type"]] == 'ins')
  del.ind = which(l2@elementMetadata@listData[["type"]] == 'del')
  fs.ind = union(ins.ind,del.ind)
  #get ptc
  ptc.ind = which(l2@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]] == T)
  #get canonical transcript
  tr_set = unique(l2@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]])
  can.info = getBM(attributes=c('ensembl_transcript_id','transcript_is_canonical'),
                   filters = 'ensembl_transcript_id',
                   values = tr_set,
                   mart = ensembl)
  tr_set2 = can.info[which(can.info$transcript_is_canonical == 1),'ensembl_transcript_id']
  #remove single exon variant
  BM.info = getBM(
    attributes=c('ensembl_transcript_id','cds_start'),
    filters = 'ensembl_transcript_id',
    values = tr_set2,
    mart = ensembl)
  exon.num = BM.info %>% group_by(ensembl_transcript_id) %>% summarise(max_rank = sum(!is.na(cds_start)))
  merge = merge(BM.info,exon.num,by='ensembl_transcript_id')
  tr_set3 = unique(merge[which(merge$max_rank > 1),'ensembl_transcript_id'])
  
  fs = l2[intersect(fs.ind,ptc.ind)]
  fs = fs[fs@elementMetadata@listData[['res_aenmd']]@listData[['transcript']] %in% tr_set3]
  #get NMDesc variant by canonical rule
  fs_NMD = fs[fs@elementMetadata@listData[['res_aenmd']]@listData[['is_last']] | fs@elementMetadata@listData[['res_aenmd']]@listData[['is_penultimate']]]
  #give count of variants
  variant_count = variant_count + length(fs_NMD)
  #give count of genes
  gene_count = gene_count + length(unique(l2@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]]))
  
  temp_df <- data.frame(
    transcript = fs_NMD@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]],
    id = fs_NMD@elementMetadata@listData[["key"]],
    stringsAsFactors = FALSE
  )
  
  # Append to main dataframe
  fs_NMD_df <- rbind(fs_NMD_df, temp_df)
}
write.csv(fs_NMD_df, '/home/jxu/fs_gnomAD_variants.csv', row.names = FALSE)