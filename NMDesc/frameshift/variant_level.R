get_re_loc = function(location,transcript_match){
  cumulative_cds_length <- 0
  cds_position <- NA
  cds_start = BM.info4[which(BM.info4$ensembl_transcript_id == transcript_match),'cds_start']
  cds_end = BM.info4[which(BM.info4$ensembl_transcript_id == transcript_match),'cds_end']
  exon_start = BM.info4[which(BM.info4$ensembl_transcript_id == transcript_match),'exon_chrom_start']
  exon_end = BM.info4[which(BM.info4$ensembl_transcript_id == transcript_match),'exon_chrom_end']
  
  for (i in seq_along(cds_start)) {
    # Check if the location falls within the current exon
    if (location >= exon_start[i] && location <= exon_end[i]) {
      # Calculate the position within the CDS
      cds_position <- cumulative_cds_length + (location - exon_start[i] + 1)
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
  print(i)
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
write.csv(fs_NMD_result,'fs_NMD_result0217.csv',row.names=FALSE)