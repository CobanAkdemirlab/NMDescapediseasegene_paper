##Library aenmd into R
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

typelist = c('all','can','css','long','trig')

resetID_clinvar_20241015_Clnsig_Origin1 = resetID_clinvar_20241015_Clnsig_Origin[which(resetID_clinvar_20241015_Clnsig_Origin$ORIGIN==1),]
resetID_clinvar_20241015_Clnsig_filtered = resetID_clinvar_20241015_Clnsig_Origin1[,-7]
write.table(resetID_clinvar_20241015_Clnsig_filtered,'~/Desktop/new_clinvar/resetID_clinvar_20241015_Clnsig_filtered.txt',quote = F,row.names = F,sep = '\t')
NMD_annotate('/Users/jxu14/Desktop/new_clinvar/resetID_clinvar_20241015_Clnsig_filtered.txt','~/Desktop/new_clinvar/clinvar_1120.rds')
res = readRDS('~/Desktop/new_clinvar/clinvar_1120.rds')
res2 = unlist(res)


snv = which(res2@elementMetadata@listData[["type"]] == 'snv')
plp = grep('pathogenic',res2@elementMetadata@listData[["clinsig"]],ignore.case = T)
ptc = which(res2@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]]==T)
ptc_nmd = res2[ptc]
plp_ptc = intersect(plp,ptc)
plp_ptc_nmd = res2[plp_ptc]
plp_ptc_nmdesc = plp_ptc_nmd[which(plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]]==T | plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]] == T | plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]] == T | plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] == T)]
temp = unique(plp_ptc_nmdesc@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]])
temp2 = getBM(attributes = c('ensembl_transcript_id','transcript_is_canonical','hgnc_symbol'),
      filters = 'ensembl_transcript_id',
      values = temp,
      mart = ensembl)
temp3 = temp2[which(temp2$transcript_is_canonical == 1),]
vus = which(res2@elementMetadata@listData[["clinsig"]] == 'Uncertain_significance')
vus_ptc = intersect(vus,ptc)
vus_ptc_nmd = res2[vus_ptc]
length(which(vus_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]]==T | vus_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]] == T | vus_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]] == T | vus_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] == T))

snv_plp_ptc = intersect(snv,plp_ptc)
can = which(plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]]==T | plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]] == T)
non_can = which(plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]]==T | plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] == T)
long = which(plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] == T)
length(intersect(can,non_can))/204795
length(intersect(can,long))/204795
length(which(plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]]==T | plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]] == T))/length(plp_ptc_nmd)
length(which(plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]]==T | plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] == T))/length(plp_ptc_nmd)
length(which(plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]]==F & plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]] == F & plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]]==F & plp_ptc_nmd@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] == F))/length(plp_ptc_nmd)

snv_vus_ptc = intersect(snv,vus_ptc)
plp_nmd = res2[plp]
snv_plp_ptc_res = res2[snv_plp_ptc]
snv_vus_ptc_res = res2[snv_vus_ptc]
snv_nmd = res2[snv]
saveRDS(snv_nmd,'~/Desktop/new_clinvar/snv_nmd1120.rds')
saveRDS(snv_plp_ptc_res,'~/Desktop/new_clinvar/raw_data/snv_plp_ptc_res1120.rds')
saveRDS(snv_vus_ptc_res,'~/Desktop/new_clinvar/raw_data/snv_vus_ptc_res1120.rds')
snv_vus_ptc_tr_set = snv_vus_ptc_res@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]]
snv_plp_ptc_tr_set = snv_plp_ptc_res@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]]
snv_vus_ptc_tr_can = getBM(attributes = c('ensembl_transcript_id','transcript_is_canonical'),filters = 'ensembl_transcript_id',values = snv_vus_ptc_tr_set,mart = ensembl)
snv_plp_ptc_tr_can = getBM(attributes = c('ensembl_transcript_id','transcript_is_canonical'),filters = 'ensembl_transcript_id',values = snv_plp_ptc_tr_set,mart = ensembl)
snv_plp_ptc_can = snv_plp_ptc_res[which(snv_plp_ptc_tr_can$transcript_is_canonical == 1)]
snv_vus_ptc_can = snv_vus_ptc_res[which(snv_vus_ptc_tr_can$transcript_is_canonical == 1)]
saveRDS(snv_plp_ptc_can,'~/Desktop/new_clinvar/raw_data/snv_plp_ptc_can1120.rds')
saveRDS(snv_vus_ptc_can,'~/Desktop/new_clinvar/raw_data/snv_vus_ptc_can1120.rds')
saveRDS(fs2,'~/Desktop/new_clinvar/fs2.rds')
fs_plus1 = fs2[which(fs_NMD_result$type == 'plus1'),]
fs_plus2 = fs2[which(fs_NMD_result$type == 'plus2'),]
saveRDS(fs_plus1,'~/Desktop/new_clinvar/fs_plus1.rds')
saveRDS(fs_plus2,'~/Desktop/new_clinvar/fs_plus2.rds')
get_pvalue('~/Desktop/new_clinvar/fs_plus1.rds','~/Desktop/new_clinvar/fs_p_plus1.rds')
get_pvalue('~/Desktop/new_clinvar/fs_plus2.rds','~/Desktop/new_clinvar/fs_p_plus2.rds')
get_pvalue('~/Desktop/new_clinvar/snv_nmd1120.rds','~/Desktop/new_clinvar/snv_nmd_p1122.rds')
get_pvalue('~/Desktop/new_clinvar/raw_data/snv_plp_ptc_res1120.rds','~/Desktop/new_clinvar/raw_data/snv_plp_ptc_p1122.rds')

snv_plp_ptc_p1122 = readRDS('~/Desktop/new_clinvar/raw_data/snv_plp_ptc_p1122.rds')
snv_plp_ptc_res1120 = readRDS('~/Desktop/new_clinvar/raw_data/snv_plp_ptc_res1120.rds')
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

for(i in 1:4){
  #get_NMD_enrichment('~/Desktop/new_clinvar/raw_data/snv_plp_ptc_p1120.rds',p_cutoff = 0.9,filter_type = typelist[i])
  get_NMD_enrichment('~/Desktop/new_clinvar/fs_p_plus1.rds',p_cutoff = 0.9,filter_type = typelist[i])
  get_NMD_enrichment('~/Desktop/new_clinvar/fs_p_plus2.rds',p_cutoff = 0.9,filter_type = typelist[i])
}



add_OMIM(gene_list)
cut=0.9
clinvar_genes = get_NMD_enrichment_mul_fliter('/Users/jxu14/Library/Mobile Documents/com~apple~CloudDocs/Desktop/桌面 - q的MacBook Pro/clinvar/rds_with_pvalue/ClinVar_p_NMD0906.rds',p_cutoff = cut,filter_type = 'all')
gnomAD_genes = get_gnomAD_de_mult(p_cut = 1-cut, filter_type = 'all')
Overlap_per(clinvar_genes,gnomAD_genes)

max_per = 0
max_min_var = 0
max_min_length = 0
max_min_MAF = 0
max_cut = 0
max_filter_type = 'all'
for(min_var in 1:4){
  for(min_length in 1:4){
    for(min_MAF in seq(from=0.01, to = 0.1, by =0.01)){
      for(cut in seq(from = 0.8,to=0.95,by=0.05)){
        for(filter_type in c('all','can','css','long')){
          clinvar_genes = get_NMD_enrichment_mul_fliter('/Users/jxu14/Library/Mobile Documents/com~apple~CloudDocs/Desktop/桌面 - q的MacBook Pro/clinvar/rds_with_pvalue/ClinVar_p_NMD0906.rds',p_cutoff = cut,filter_type = filter_type)
          gnomAD_genes = get_gnomAD_de_mult(p_cut = 1-cut, filter_type = filter_type)
          new_per = Overlap_per(clinvar_genes,gnomAD_genes)
          if(new_per > max_per){
            max_per = new_per
            max_min_var = min_var
            max_min_length = min_length
            max_min_MAF = min_MAF
            max_cut = cut
            max_filter_type = filter_type
          }
        }
      }
    }
  }
}

#how to reduce the time for runnning the gnomAD data? maybe get a dataframe with only all required properties, and then filter out genes 
#by certain properties. To do that, we muat first be sure about which properties are required, like, which is the quality control property?
#but one gene may have multiple variants, each variant has a MAF, thus it might not be workable

enrich_path = "/Users/jxu14/Downloads/pcut_90%/p/Clinvar_p_NMD0906_NMDenriched_can.txt"
enrich1(enrich_path)
