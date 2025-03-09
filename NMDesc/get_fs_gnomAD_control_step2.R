library(stringr)
res = readRDS('~/Desktop/new_clinvar/raw_data/clinvar_1112.rds')
res2 = unlist(res)
plp = grep('pathogenic',res2@elementMetadata@listData[["clinsig"]],ignore.case = T)

#germline

#PTC
ptc = which(res2@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]]==T)

#fs
fs = which(res2@elementMetadata@listData[["type"]] %in% c('ins','del'))

#remove single exon
no_single = which(res2@elementMetadata@listData[["res_aenmd"]]@listData[["is_single"]] == F)

combined_ind = intersect(plp,ptc)
fs_ind = intersect(combined_ind,fs)
fs_ind2 = intersect(fs_ind,no_single)

#get the canonical transcript
transcript_list = unique(res2@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][fs_ind2])

gene_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "transcript_is_canonical", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = transcript_list,
  mart = mart
)
gene_list = gene_mapping[which(gene_mapping$transcript_is_canonical == 1),]

#get all fs variants
fs_res = res2[fs_ind2]
fs_res_can = fs_res[which(fs_res@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% gene_list$ensembl_transcript_id),]
#select plus1 variants
get_frameshift_type = function(res,x)
{
  #retrive ref and alt info
  temp = res[x]
  key = temp@elementMetadata@listData[["id"]]
  ref = str_split(key, "_")[[1]][3]
  alt = str_split(key, "_")[[1]][4]
  type = temp@elementMetadata@listData[["type"]]
  yu = (abs(nchar(ref)-nchar(alt)))%%3
  #if delete 1(4,7) or insert 2 bp, then plus1; 
  if(type=='del' & yu==1 | type=='ins' & yu==2){
    plus_type = 'plus1'
  }
  #if delete 2 or insert 1 bp, then plus2
  else if(type=='del' & yu==2 | type=='ins' & yu==1){
    plus_type = 'plus2'
  }else{
    plus_type = 'error'
  }
  return(plus_type)
}
#run get_frameshift_type on fs_res_can using loop
fs_type = rep(NA,length(fs_res_can))
for(i in 1:length(fs_res_can)){
  fs_type[i] = get_frameshift_type(fs_res_can,i)
}
plus1_res_can = fs_res_can[which(fs_type=='plus1')]
##exclude any gene genes with any plus1 NMD_can_esc mutation
gene2remove.ind = which(plus1_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |plus1_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)
gene2remove = unique(plus1_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][gene2remove.ind])
#remove genes2remove from gene_list$ensemble_transcript_id
genes_remain_in_clinvar = gene_list[-which(gene_list$ensembl_transcript_id %in% gene2remove),]

##filter for gnomad control
fs_gnomAD_variants <- read_csv("~/Desktop/fs_gnomAD_variants.csv")
#get tyoe of fs_gnomAD_variants based on id
fs_gnomAD_variants_type = rep(NA,length(fs_gnomAD_variants$id))
get_frameshift_type2 = function(df,x)
{
  #retrive ref and alt info
  key = df[i,2]$id
  parts <- strsplit(key, "\\|")[[1]]
  ref <- parts[2] 
  alt = parts[3]

  if(nchar(ref)>nchar(alt)){
    type = 'del'
    }else if(nchar(ref)<nchar(alt)){
      type = 'ins'
    }
  yu = (abs(nchar(ref)-nchar(alt)))%%3
  #if delete 1(4,7) or insert 2 bp, then plus1; 
  if(type=='del' & yu==1 | type=='ins' & yu==2){
    plus_type = 'plus1'
  }
  #if delete 2 or insert 1 bp, then plus2
  else if(type=='del' & yu==2 | type=='ins' & yu==1){
    plus_type = 'plus2'
  }else{
    plus_type = 'error'
  }
  return(plus_type)
}
for(i in 1:length(fs_gnomAD_variants$id)){
  fs_gnomAD_variants_type[i] = get_frameshift_type2(fs_gnomAD_variants,i)
}
plus1_NMD_results = fs_gnomAD_variants[which(fs_gnomAD_variants_type=='plus1'),] #this is from gnomAD

plus1_gnomAD_control = plus1_NMD_results[which(plus1_NMD_results$transcript %in% genes_remain_in_clinvar$ensembl_transcript_id),]
write.csv(plus1_gnomAD_control, 'plus1_gnomAD_control_variant.csv', row.names = FALSE)
#get hgnc_symbol of plus1_gnomAD_control
plus1_gnomAD_control_hgnc = getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = unique(plus1_gnomAD_control$transcript),
  mart = mart
)
write.csv(unique(plus1_gnomAD_control_hgnc$hgnc_symbol), 'plus1_gnomAD_control_genes.csv', row.names = FALSE)
