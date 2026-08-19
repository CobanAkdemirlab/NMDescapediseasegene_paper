library(dplyr)
library(stringr)
library(biomaRt)

bp=50
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

res = readRDS('~/Desktop/new_clinvar/clinvar_1120.rds')
res2 = unlist(res)
ptc = which(res2@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]]==T)
plp = grep('pathogenic',res2@elementMetadata@listData[["clinsig"]],ignore.case = T)
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
get_PTC_plus_ind <- function(input_seqs,type='plus1')
{
  if(type=='plus1'){
    para = 2
  }else if(type=='plus2'){
    para = 0
  }
  a <- unlist(input_seqs)
  ind <- unlist(gregexpr('TAG|TAA|TGA',as.character(a)))
  plus.ind <- ind[which(ind%%3==para)]
  collapsed.plus=paste(plus.ind,collapse=',')
  return(collapsed.plus)
}
get_NMDesc_PTC = function(input_seqs,bp,variant_transcript,type='plus1')
{
  if(type=='plus1'){
    para = 2
  }else if(type=='plus2'){
    para = 0
  }
  a <- unlist(input_seqs)
  ind <- unlist(gregexpr('TAG|TAA|TGA',as.character(a)))
  plus1.ind <- ind[which(ind%%3==para)]
  df <- BM.info4[which(BM.info4$ensembl_transcript_id==variant_transcript),]
  #remove NA rows in df
  df = df[!is.na(df$cds_start),]
  ### get the length of each exon
  exon.length <- (df$cds_end-df$cds_start)+1
  
  ##rule 1 The PTC located in the last coding exon and rule 2 The PTC located within d_pen bp upstream of the penultimate exon boundary (penultimate exon rule; default: d_pen = 50)
  #check pen exon
  pen.length = exon.length[length(exon.length)-1]
  if(pen.length<bp){
    aft.NMD.ind = df$cds_start[length(exon.length)-1]
  }else{
    aft.NMD.ind <- sum(exon.length[1:length(exon.length)-1])-bp 
  }
  NMD_status_can = plus1.ind>=aft.NMD.ind
  #rule 3 The PTC located within d_css bp downstream of the coding start site (css rule default: d_css = 150)
  NMD_status_css = plus1.ind<=150
  #rule 4 The PTC located within an exon spanning more than 407bp 
  ##find which exon>407
  exon.long <- which(exon.length>407)
  #exclude the last exon(remove length(exon.length) from exon.long)
  exon.long <- exon.long[exon.long != length(exon.length)]
  #exclude the first exon
  exon.long <- exon.long[exon.long != 1]
  ##find if the PTC is in that long exon
  long.start = df$cds_start[exon.long]
  long.end = df$cds_end[exon.long]
  if(length(exon.long)==0){
    long.start = 0
    long.end = 0
  }
  NMD_status_long = rep(FALSE,length(plus1.ind))
  for(k in 1:length(plus1.ind)){
    NMD_status_long[k] = plus1.ind[k]>=long.start & plus1.ind[k]<=long.end
  }
  
  
  #return all the rules
  #NMD_status = NMD_status_can | NMD_status_css | NMD_status_long
  NMD_status = c(NMD_status_can,NMD_status_css,NMD_status_long)
  
  return(paste(NMD_status,collapse=','))
}
get_NMDesc_PTC_region = function(PTC_ind,variant_transcript,PTC_status)
{
  PTC_status_vector = as.logical(strsplit(PTC_status, ",")[[1]])
  can_PTC_status_vector = PTC_status_vector[1:(length(PTC_status_vector)/3)]
  css_PTC_status_vector = PTC_status_vector[(length(PTC_status_vector)/3+1):(length(PTC_status_vector)/3*2)]
  long_PTC_status_vector = PTC_status_vector[(length(PTC_status_vector)/3*2+1):length(PTC_status_vector)]
  PTC_ind_vector = as.numeric(unlist(strsplit(PTC_ind, ",")))
  can_change_points = which(diff(can_PTC_status_vector) != 0)
  can_change_points_loc = PTC_ind_vector[can_change_points]
  css_change_points = which(diff(css_PTC_status_vector) != 0)
  css_change_points_loc = PTC_ind_vector[css_change_points]
  long_change_points = which(diff(long_PTC_status_vector) != 0)
  long_change_points_loc = PTC_ind_vector[long_change_points]
  #find all the F->T change point
  long_FT_change_points <- long_change_points_loc[long_PTC_status_vector[long_change_points] == FALSE]
  #find all the T->F change point
  long_TF_change_points <- long_change_points_loc[long_PTC_status_vector[long_change_points] == TRUE]
  can_region = max((PTC_ind_vector[length(PTC_ind_vector)] - (can_change_points_loc+1)),0)
  css_region = max((css_change_points_loc - 1),0)
  long_region = max(long_TF_change_points - (long_FT_change_points+1),0)
  #if the penultimate exon is long, -150bp from the long.region
  df <- BM.info4[which(BM.info4$ensembl_transcript_id==variant_transcript),]
  #remove NA rows in df
  df = df[!is.na(df$cds_start),]
  ### get the length of each exon
  exon.length <- (df$cds_end-df$cds_start)+1
  pen.length = exon.length[length(exon.length)-1]
  if(pen.length>407){
    long_region = long_region - 150
  }
  
  if(length(long_region) == 0 | is.na(long_region)){
    long_region = 0
  }
  if(length(css_region) == 0 | is.na(css_region)){
    css_region = 0
  }
  if(length(can_region) == 0 | is.na(can_region)){
    can_region = 0
  }
  list(
    can_region = can_region,
    css_region = css_region,
    long_region = long_region
  )
}

PTC_info = data.frame(PTC_loc = rep(NA,2*length(transcript_set3)),
                      type = rep(NA,2*length(transcript_set3)),
                      transcript = rep(NA,2*length(transcript_set3)),
                      PTC_NMDesc = rep(NA,2*length(transcript_set3)),
                      can_region = rep(NA,2*length(transcript_set3)),
                      css_region = rep(NA,2*length(transcript_set3)),
                      long_region = rep(NA,2*length(transcript_set3)))
counter = 1
for (i in 1:length(transcript_set3)) {
  for (j in c('plus1', 'plus2')) {
      transcript <- transcript_set3[i]
      input_seqs <- cds.info[which(cds.info$ensembl_transcript_id == transcript), 'coding']
      PTC_ind <- get_PTC_plus_ind(input_seqs, type = j)
      PTC_status <- get_NMDesc_PTC(input_seqs, bp, transcript, type = j)
      region_set <- get_NMDesc_PTC_region(PTC_ind,transcript, PTC_status)
      # Update PTC_info
      PTC_info[counter, 'PTC_loc'] <- PTC_ind
      PTC_info[counter, 'type'] <- j
      PTC_info[counter, 'transcript'] <- transcript
      PTC_info[counter, 'PTC_NMDesc'] <- PTC_status
      PTC_info[counter, "can_region"] <- region_set$can_region
      PTC_info[counter, "css_region"] <- region_set$css_region
      PTC_info[counter, "long_region"] <- region_set$long_region
      counter <- counter + 1
  }
}

write.csv(PTC_info, file = 'PTC_info0217.csv', row.names = FALSE)
