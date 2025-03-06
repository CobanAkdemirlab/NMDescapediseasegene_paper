variant_count = 0
gene_count = 0
names = paste('chr',seq(1,22),'.rds',sep='')
snv_NMD_df <- data.frame(transcript = character(), id = character(), stringsAsFactors = FALSE)
for(outfilename in names){
  l2 <- readRDS(outfilename)
  #filter only snv variant in l
  #l2 = unlist(l)
  snv.ind = which(l2@elementMetadata@listData[["type"]] == 'snv')
  #get ptc
  ptc.ind = which(l2@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]] == T)
  #get benign variant
  #benign.ind = grep('benign',l2@elementMetadata@listData[["clinsig"]],ignore.case=T)
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
  
  
  snv = l2[intersect(snv.ind,ptc.ind)]
  snv = snv[snv@elementMetadata@listData[['res_aenmd']]@listData[['transcript']] %in% tr_set3]
  #get NMDesc variant by canonical rule
  snv_NMD = snv[snv@elementMetadata@listData[['res_aenmd']]@listData[['is_last']] | snv@elementMetadata@listData[['res_aenmd']]@listData[['is_penultimate']]]
  #give count of variants
  variant_count = variant_count + length(snv_NMD)
  #give count of genes
  gene_count = gene_count + length(unique(l2@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]]))
  
  temp_df <- data.frame(
    transcript = snv_NMD@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]],
    id = snv_NMD@elementMetadata@listData[["key"]],
    stringsAsFactors = FALSE
  )
  
  # Append to main dataframe
  snv_NMD_df <- rbind(snv_NMD_df, temp_df)
}

perl convert2annovar.pl -format vcf4 snv_NMD_results3.vcf > snv_NMD.avinput
sed -i 's/^1/chr1/; s/^2/chr2/; s/^3/chr3/; s/^4/chr4/; s/^5/chr5/; s/^6/chr6/; s/^7/chr7/; s/^8/chr8/; s/^9/chr9/; s/^10/chr10/; s/^11/chr11/; s/^12/chr12/; s/^13/chr13/; s/^14/chr14/; s/^15/chr15/; s/^16/chr16/; s/^17/chr17/; s/^18/chr18/; s/^19/chr19/; s/^20/chr20/; s/^21/chr21/; s/^22/chr22/; s/^X/chrX/; s/^Y/chrY/' ~/gnomAD/annovar/snv_NMD.avinput
perl ~/gnomAD/annovar/annotate_variation.pl   -buildver hg38   -out ~/gnomAD/snv_NMD_gnomAD41   -filter   -dbtype gnomad41_exome   -maf 0.01   ~/gnomAD/annovar/snv_NMD.avinput   ~/gnomAD/annovar/humandb/
