#This script is to read in data download from clinvar, annotate with aenmd, and select fs/snv variants, separate by plp/vus, and canonical/non-canonical transcripts
#It's also include codes to get genes that are enriched with NMDesc variants

library(aenmd.data.ensdb.v105)
library(aenmd)
library(GenomicRanges)
library(tidyr)
library(GenomicFeatures)
library(VariantAnnotation)
library(AnnotationDbi)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biomaRt)
library(enrichR)
library(aenmd)

#filter for germline variants
resetID_clinvar_20241015_Clnsig_Origin1 = resetID_clinvar_20241015_Clnsig_Origin[which(resetID_clinvar_20241015_Clnsig_Origin$ORIGIN==1),]
resetID_clinvar_20241015_Clnsig_filtered = resetID_clinvar_20241015_Clnsig_Origin1[,-7]
write.table(resetID_clinvar_20241015_Clnsig_filtered,'~/Desktop/new_clinvar/resetID_clinvar_20241015_Clnsig_filtered.txt',quote = F,row.names = F,sep = '\t')
NMD_annotate('/Users/jxu14/Desktop/new_clinvar/resetID_clinvar_20241015_Clnsig_filtered.txt','~/Desktop/new_clinvar/clinvar_1120.rds')
res = readRDS('~/Desktop/new_clinvar/clinvar_1120.rds')
res2 = unlist(res)

#separate for snv and framshift
snv = which(res2@elementMetadata@listData[["type"]] == 'snv')
plp = grep('pathogenic',res2@elementMetadata@listData[["clinsig"]],ignore.case = T) #there are pathogenic and Pathogenic
ptc = which(res2@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]]==T)
ptc_nmd = res2[ptc]
plp_ptc = intersect(plp,ptc)
plp_ptc_snv = res2[intersect(plp_ptc,snv)]
plp_ptc_nmdesc = plp_ptc_snv[which(plp_ptc_snv@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]]==T | plp_ptc_snv@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]] == T | plp_ptc_snv@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]] == T | plp_ptc_snv@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] == T)]
temp = unique(plp_ptc_nmdesc@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]])
temp2 = getBM(attributes = c('ensembl_transcript_id','transcript_is_canonical','hgnc_symbol'),
      filters = 'ensembl_transcript_id',
      values = temp,
      mart = ensembl) #get canonical transcript info
temp3 = temp2[which(temp2$transcript_is_canonical == 1),]
snv_plp_ptc_canonical = plp_ptc_nmdesc[which(plp_ptc_nmdesc@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% temp3$ensembl_transcript_id),]


vus = which(res2@elementMetadata@listData[["clinsig"]] == 'Uncertain_significance')
vus_ptc = intersect(vus,ptc)
vus_ptc_nmd = res2[vus_ptc]
length(which(vus_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]]==T | vus_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]] == T | vus_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]] == T | vus_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] == T))

snv_plp_ptc = intersect(snv,plp_ptc)
can = which(plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]]==T | plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]] == T)
length(which(plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]]==T | plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]] == T))/length(plp_ptc_nmd)
length(which(plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]]==T | plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] == T))/length(plp_ptc_nmd)
length(which(plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]]==F & plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]] == F & plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]]==F & plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] == F))/length(plp_ptc_nmd)

snv_vus_ptc = intersect(snv,vus_ptc)
plp_nmd = res2[plp]
snv_plp_ptc_res = res2[snv_plp_ptc]
snv_vus_ptc_res = res2[snv_vus_ptc]
snv_nmd = res2[snv]
saveRDS(snv_nmd,'~/Desktop/new_clinvar/snv_nmd1120.rds')
saveRDS(temp3,'~/Desktop/new_clinvar/raw_data/snv_plp_ptc_res1120.rds') #temp3 is snv_plp_ptc_canonical

saveRDS(snv_vus_ptc_res,'~/Desktop/new_clinvar/raw_data/snv_vus_ptc_res1120.rds')
snv_vus_ptc_tr_set = snv_vus_ptc_res@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]]
snv_plp_ptc_tr_set = snv_plp_ptc_res@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]]
snv_vus_ptc_tr_can = getBM(attributes = c('ensembl_transcript_id','transcript_is_canonical'),filters = 'ensembl_transcript_id',values = snv_vus_ptc_tr_set,mart = ensembl)
snv_plp_ptc_tr_can = getBM(attributes = c('ensembl_transcript_id','transcript_is_canonical'),filters = 'ensembl_transcript_id',values = snv_plp_ptc_tr_set,mart = ensembl)
snv_plp_ptc_can = snv_plp_ptc_res[which(snv_plp_ptc_tr_can$transcript_is_canonical == 1)]
snv_vus_ptc_can = snv_vus_ptc_res[which(snv_vus_ptc_tr_can$transcript_is_canonical == 1)]
saveRDS(snv_plp_ptc_canonical,'~/Desktop/new_clinvar/raw_data/snv_plp_ptc_can1120.rds')
saveRDS(snv_vus_ptc_can,'~/Desktop/new_clinvar/raw_data/snv_vus_ptc_can1120.rds') 

#select include canonical frameshift variants
res_plp_ptc = res2[intersect(ptc,plp)]
del = which(res_plp_ptc@elementMetadata@listData[["type"]] == 'del')
ins = which(res_plp_ptc@elementMetadata@listData[["type"]] == 'ins')
fs = res_plp_ptc[union(del,ins)]
transcript_set = unique(fs@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]])
can.info = getBM(
  attributes = c("ensembl_transcript_id","transcript_is_canonical"), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values =  transcript_set,                          # Transcript ID
  mart = ensembl                                         # Database connection
)
#exclude not canonial transcripts
transcript_set2 = can.info[which(can.info$transcript_is_canonical==1),"ensembl_transcript_id"]
BM.info <- getBM(
  attributes = c("ensembl_transcript_id","rank",'cds_start','cds_end','exon_chrom_start','exon_chrom_end'), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values =  transcript_set2,                          # Transcript ID
  mart = ensembl                                         # Database connection
)
#exclude single (coding) exon genes
exon.num = BM.info %>% group_by(ensembl_transcript_id) %>% summarise(max_rank = sum(!is.na(cds_start)))
merge = merge(BM.info,exon.num,by='ensembl_transcript_id')
transcript_set3 = unique(merge[which(merge$max_rank>1),"ensembl_transcript_id"])
cds.info = getBM(
  attributes = c("ensembl_transcript_id","coding"), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values =  transcript_set3,                          # Transcript ID
  mart = ensembl                                         # Database connection
)
BM.info2 = BM.info[which(BM.info$ensembl_transcript_id %in% transcript_set3),]
#calculate cds length, keep the coding exons
BM.info2$exon_length = BM.info2$cds_end - BM.info2$cds_start + 1
BM.info3 = BM.info2[which(BM.info2$exon_length>0),]
BM.info4 = BM.info3 %>% group_by(ensembl_transcript_id) %>% summarise(cds_length = sum(exon_length))
BM.info4 = merge(BM.info3,BM.info4,by='ensembl_transcript_id')
#according to the new transcript list, filter input canonical variants
fs2 = fs[which(fs@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% transcript_set3),]
saveRDS(fs2,'~/Desktop/new_clinvar/fs2.rds') 

fs_plus1 = fs2[which(fs_NMD_result$type == 'plus1'),]
fs_plus2 = fs2[which(fs_NMD_result$type == 'plus2'),]
saveRDS(fs_plus1,'~/Desktop/new_clinvar/fs_plus1.rds')
saveRDS(fs_plus2,'~/Desktop/new_clinvar/fs_plus2.rds')
get_pvalue('~/Desktop/new_clinvar/fs_plus1.rds','~/Desktop/new_clinvar/fs_p_plus1.rds')
get_pvalue('~/Desktop/new_clinvar/fs_plus2.rds','~/Desktop/new_clinvar/fs_p_plus2.rds')
get_pvalue('~/Desktop/new_clinvar/raw_data/snv_plp_ptc_can1120.rds','~/Desktop/new_clinvar/snv_nmd_p0411.rds')
get_pvalue('~/Desktop/new_clinvar/raw_data/snv_plp_ptc_res1120.rds','~/Desktop/new_clinvar/raw_data/snv_plp_ptc_p1122.rds')

snv_plp_ptc_p1122 = readRDS('~/Desktop/new_clinvar/raw_data/snv_plp_ptc_p1122.rds')
snv_plp_ptc_res1120 = readRDS('~/Desktop/new_clinvar/raw_data/snv_plp_ptc_res1120.rds')

#some test code to find specific genes
#use getBM, find the transcript name for NOVA2
gene_list = getBM(attributes = c('hgnc_symbol','ensembl_transcript_id',"transcript_is_canonical"), filters = 'hgnc_symbol', values = 'NOVA2', mart = mart)
for (i in 1:6000) {
  tryCatch({
    if (snv_plp_ptc_p1122[[i]][["hgnc_symbol"]] == "MPZ") {
      print(i)
    }
  }, error = function(e) {
    # Do nothing, simply continue to the next iteration
  })
}


temp = readRDS('~/Desktop/new_clinvar/raw_data/snv_plp_ptc_p1122.rds')
#get geen enrichment
for(i in 1:4){
  get_NMD_enrichment('~/Desktop/new_clinvar/snv_nmd_p0411.rds',p_cutoff = 0.9,filter_type = typelist[2])
  get_NMD_enrichment('~/Desktop/new_clinvar/fs_p_plus1.rds',p_cutoff = 0.9,filter_type = typelist[i])
  get_NMD_enrichment('~/Desktop/new_clinvar/fs_p_plus2.rds',p_cutoff = 0.9,filter_type = typelist[i])
}

