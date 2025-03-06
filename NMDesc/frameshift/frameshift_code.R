library(dplyr)
library(stringr)
library(biomaRt)

bp=50
#ensembl2 <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")

res = readRDS('~/Desktop/new_clinvar/clinvar_1120.rds')
res2 = unlist(res)
ptc = which(res2@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]]==T)
plp = grep('pathogenic',res2@elementMetadata@listData[["clinsig"]],ignore.case = T)
#res_plp = res2[plp]
res_plp = res2[plp]
del = which(res2@elementMetadata@listData[["type"]] == 'del')
ins = which(res2@elementMetadata@listData[["type"]] == 'ins')
fs = res2[union(del,ins)]
del_plp = which(res_plp@elementMetadata@listData[["type"]] == 'del')
ins_plp = which(res_plp@elementMetadata@listData[["type"]] == 'ins')
fs_plp = res_plp[union(del_plp,ins_plp)]

transcript_set = unique(fs@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]])
transcript_set_plp = unique(fs_plp@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]])
can.info = getBM(
  attributes = c("ensembl_transcript_id","transcript_is_canonical"), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values =  transcript_set,                          # Transcript ID
  mart = ensembl                                         # Database connection
)
can.info2 = can.info[which(can.info$ensembl_transcript_id %in% transcript_set_plp),]
#exclude not canonial transcripts
transcript_set2 = can.info2[which(can.info2$transcript_is_canonical==1),"ensembl_transcript_id"]
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
transcript_set3_plp = transcript_set3[which(transcript_set3 %in% transcript_set_plp)]
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
BM.info5 = BM.info4[which(BM.info4$ensembl_transcript_id %in% transcript_set_plp),]

#sort BM.info5 by rank
BM.info5 = BM.info5[order(BM.info5$ensembl_transcript_id, BM.info5$rank), ]

#according to the new transcript list, filter input variants
#fs1.5 = fs[which(fs@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% transcript_set3_plp_plp),]
fs2 = fs[which(fs@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% transcript_set3_plp_plp),]

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
  df <- BM.info5[which(BM.info5$ensembl_transcript_id==variant_transcript),]
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
  df <- BM.info5[which(BM.info5$ensembl_transcript_id==variant_transcript),]
  #remove NA rows in df
  df = df[!is.na(df$cds_start),]
  ### get the length of each exon
  exon.length <- (df$cds_end-df$cds_start)+1
  pen.length = exon.length[length(exon.length)-1]
  if(pen.length>407){
    long_region = long_region - 150
  }
  long_region = max(long_region,0)
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

PTC_info = data.frame(PTC_loc = rep(NA,2*length(transcript_set3_plp_plp)),
                      type = rep(NA,2*length(transcript_set3_plp_plp)),
                      transcript = rep(NA,2*length(transcript_set3_plp_plp)),
                      PTC_NMDesc = rep(NA,2*length(transcript_set3_plp_plp)),
                      can_region = rep(NA,2*length(transcript_set3_plp_plp)),
                      css_region = rep(NA,2*length(transcript_set3_plp_plp)),
                      long_region = rep(NA,2*length(transcript_set3_plp_plp)))
counter = 1
for (i in 1:length(transcript_set3_plp_plp)) {
  for (j in c('plus1', 'plus2')) {
    transcript <- transcript_set3_plp_plp[i]
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

get_re_loc = function(location,transcript_match){
  cumulative_cds_length <- 0
  cds_position <- NA
  cds_start = BM.info5[which(BM.info5$ensembl_transcript_id == transcript_match),'cds_start']
  cds_end = BM.info5[which(BM.info5$ensembl_transcript_id == transcript_match),'cds_end']
  exon_start = BM.info5[which(BM.info5$ensembl_transcript_id == transcript_match),'exon_chrom_start']
  exon_end = BM.info5[which(BM.info5$ensembl_transcript_id == transcript_match),'exon_chrom_end']
  ex.length = BM.info5[which(BM.info5$ensembl_transcript_id == transcript_match),'exon_length']
  for (i in seq_along(cds_start)) {
    # Check if the location falls within the current exon
    if (location >= exon_start[i] && location <= exon_end[i]) {
      # Calculate the position within the CDS
      cds_position <- cumulative_cds_length + min((location - exon_start[i] + 1),ex.length[i]) #if this exon has nocoding region
      break
    }
    
    # Update cumulative CDS length
    cumulative_cds_length <- cumulative_cds_length + (cds_end[i] - cds_start[i] + 1)
  }
  return(cds_position)
}
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
get_NMDesc_variant = function(variant_loc_rel,PTC_loc,PTC_staus)
{
  #match variant to the closest PTC
  PTC_loc2 = as.numeric(unlist(strsplit(PTC_loc, ",")))
  #transform variant_loc to make it comparable with PTC_loc
  dist = PTC_loc2-variant_loc_rel
  ##PTC should be after the variant
  dist2 = dist[dist>0]
  ##if all PTC is before the variant, then return all FALSE, which means no NMDesc
  if(length(dist2)==0){
    return(c(FALSE,FALSE,FALSE))
  }
  min.dist = min(dist2)
  min.ind = which(dist==min.dist)
  PTC_status_vector = as.logical(strsplit(PTC_status, ",")[[1]])
  #NMDesc_status = PTC_status_vector[min.ind]
  jump = length(PTC_status_vector)/3
  NMDesc_can = PTC_status_vector[min.ind]
  NMDesc_css = PTC_status_vector[min.ind+jump]
  NMDesc_long = PTC_status_vector[min.ind+jump*2]
  return(c(NMDesc_can,NMDesc_css,NMDesc_long))
}

key_set = fs2@elementMetadata@listData[["id"]]
loc_set = rep(NA,length(key_set))
for(i in 1:length(key_set)){
  loc_set[i] = str_split(key_set[i], "_")[[1]][2]
}

fs_NMD_result = data.frame(NMDesc_can = rep(NA,length(fs2)),
                           NMDesc_css = rep(NA,length(fs2)),
                           NMDesc_long = rep(NA,length(fs2)),
                           variant_loc = rep(NA,length(fs2)),
                           variant_re_loc = rep(NA,length(fs2)),
                           type = rep(NA,length(fs2)),
                           transcript_id = rep(NA,length(fs2)))
for(i in 1:length(fs2)){
  #print(i)
  fs_type = get_frameshift_type(fs2,i)
  variant_loc = loc_set[i]
  variant_transcript = fs2@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][i]
  variant_re_loc = get_re_loc(as.numeric(variant_loc),variant_transcript)
  #get PTC_ind and PTC_status from PTC_info
  PTC_ind = PTC_info[which(PTC_info$transcript == variant_transcript & PTC_info$type == fs_type),'PTC_loc']
  PTC_status = PTC_info[which(PTC_info$transcript == variant_transcript & PTC_info$type == fs_type),'PTC_NMDesc']
  variant_status = get_NMDesc_variant(variant_re_loc,PTC_ind,PTC_status)
  fs_NMD_result[i,'NMDesc_can'] = variant_status[1]
  fs_NMD_result[i,'NMDesc_css'] = variant_status[2]
  fs_NMD_result[i,'NMDesc_long'] = variant_status[3]
  fs_NMD_result[i,'type'] = fs_type
  fs_NMD_result[i,'transcript_id'] = variant_transcript
  fs_NMD_result[i,'variant_loc'] = variant_loc
  fs_NMD_result[i,'variant_re_loc'] = variant_re_loc
}

transcript_object <- list()

for (transcript in transcript_set3_plp) {
  # For each transcript, organize the data by plus
  PTC_plus1_idx <- which(PTC_info$transcript == transcript & PTC_info$type == "plus1")
  PTC_plus2_idx <- which(PTC_info$transcript == transcript & PTC_info$type == "plus2") 
  va_plus1_idx <- which(fs_NMD_result$transcript_id == transcript & fs_NMD_result$type == "plus1")
  va_plus2_idx <- which(fs_NMD_result$transcript_id == transcript & fs_NMD_result$type == "plus2")
  
  # Add data for each plus
  transcript_object[[transcript]] <- list(
    plus1 = list(
      PTC_loc = PTC_info$PTC_loc[PTC_plus1_idx],
      PTC_NMDesc = PTC_info$NMDesc[PTC_plus1_idx],
      can_region = PTC_info$can_region[PTC_plus1_idx],
      css_region = PTC_info$css_region[PTC_plus1_idx],
      long_region = PTC_info$long_region[PTC_plus1_idx],
      NMDesc_can = fs_NMD_result$NMDesc_can[va_plus1_idx],
      NMDesc_css = fs_NMD_result$NMDesc_css[va_plus1_idx],
      NMDesc_long = fs_NMD_result$NMDesc_long[va_plus1_idx]
    ),
    plus2 = list(
      PTC_loc = PTC_info$PTC_loc[PTC_plus2_idx],
      PTC_NMDesc = PTC_info$NMDesc[PTC_plus2_idx],
      can_region = PTC_info$can_region[PTC_plus2_idx],
      css_region = PTC_info$css_region[PTC_plus2_idx],
      long_region = PTC_info$long_region[PTC_plus2_idx],
      NMDesc_can = fs_NMD_result$NMDesc_can[va_plus2_idx],
      NMDesc_css = fs_NMD_result$NMDesc_css[va_plus2_idx],
      NMDesc_long = fs_NMD_result$NMDesc_long[va_plus2_idx]
    )
  )
}

