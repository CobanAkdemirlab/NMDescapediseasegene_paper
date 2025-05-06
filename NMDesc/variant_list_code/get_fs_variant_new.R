##get fs variants from clinvar
library(stringr)
#input gene list
plus1_gene = plus1_pli$gene
minus1_gene = minus1_pli$gene
#get canonical transcript of fs_gene$hgnc
plus1_tr = getBM(attributes=c('ensembl_transcript_id','hgnc_symbol'),
               filters = 'hgnc_symbol', values = plus1_gene, mart = ensembl)
minus1_tr = getBM(attributes=c('ensembl_transcript_id','hgnc_symbol'),
                filters = 'hgnc_symbol', values = minus1_gene, mart = ensembl)
plus1_can_tr = getBM(attributes=c('ensembl_transcript_id','transcript_is_canonical'),
                   filters = 'ensembl_transcript_id', values = plus1_tr$ensembl_transcript_id, mart = ensembl)
minus1_can_tr = getBM(attributes=c('ensembl_transcript_id','transcript_is_canonical'),
                        filters = 'ensembl_transcript_id', values = minus1_tr$ensembl_transcript_id, mart = ensembl)

plus1_can_tr = plus1_can_tr[which(plus1_can_tr$transcript_is_canonical == 1),1]
minus1_can_tr = minus1_can_tr[which(minus1_can_tr$transcript_is_canonical == 1),1]

#in clinvar fs ptc, nmdesc variants, select those corresponding to the canonical transcript of gene list
fs = which(res2@elementMetadata@listData[["type"]] %in% c('ins','del'))
#remove single exon
no_single = which(res2@elementMetadata@listData[["res_aenmd"]]@listData[["is_single"]] == F)
combined_ind = intersect(plp,ptc)
fs_ind = intersect(combined_ind,fs)
fs_ind2 = intersect(fs_ind,no_single)
fs_res = res2[fs_ind2]
#use fs_res to get plus1_res and minus1_res
##1. categorize by ref and alt in key
get_frameshift_type = function(key,type)
{
  ref = str_split(key, "\\|")[[1]][2]
  alt = str_split(key, "\\|")[[1]][3]
  yu = (abs(nchar(ref)-nchar(alt)))%%3
  #if delete 1(4,7) or insert 2 bp, then plus1; 
  if(type=='del' & yu==1 | type=='ins' & yu==2){
    plus_type = 'minus1'
  }
  #if delete 2 or insert 1 bp, then plus2
  else if(type=='del' & yu==2 | type=='ins' & yu==1){
    plus_type = 'plus1'
  }else{
    plus_type = 'error'
  }
  return(plus_type)
}
type <- mapply(get_frameshift_type,
                  fs_res@elementMetadata@listData[["key"]],
                  fs_res@elementMetadata@listData[["type"]])
plus1_res = fs_res[which(type == 'plus1')]
minus1_res = fs_res[which(type == 'minus1')]
#select NMDesc variants by canonical rule
plus1_res_can = plus1_res[which(plus1_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |plus1_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)]
minus1_res_can = minus1_res[which(minus1_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |minus1_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)]
plus1_variants = data.frame(transcript = plus1_res_can@ranges@NAMES[which(plus1_res_can@ranges@NAMES %in% plus1_can_tr)],
                          key = plus1_res_can@elementMetadata@listData[["key"]][which(plus1_res_can@ranges@NAMES %in% plus1_can_tr)])
minus1_variants = data.frame(transcript = minus1_res_can@ranges@NAMES[which(minus1_res_can@ranges@NAMES %in% minus1_can_tr)],
                            key = minus1_res_can@elementMetadata@listData[["key"]][which(minus1_res_can@ranges@NAMES %in% minus1_can_tr)])
#get uniprot id using biomart
plus1_unique_transcripts <- unique(plus1_variants$transcript)
minus1_unique_transcripts <- unique(minus1_variants$transcript)
plus1_uniprot_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = plus1_unique_transcripts,
  mart = ensembl
)
minus1_uniprot_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = minus1_unique_transcripts,
  mart = ensembl
)
plus1_variants<- merge(
  plus1_variants,
  plus1_uniprot_mapping,
  by.x = "transcript",
  by.y = "ensembl_transcript_id",
  all.x = TRUE
)
minus1_variants<- merge(
  minus1_variants,
  minus1_uniprot_mapping,
  by.x = "transcript",
  by.y = "ensembl_transcript_id",
  all.x = TRUE
)


write.csv(plus1_variants, 'plus1_variants0406.csv', row.names = FALSE)
write.csv(minus1_variants, 'minus1_variants0406.csv', row.names = FALSE)
