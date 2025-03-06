library(dplyr)
library(stringr)
library(readr)

get_second_largest_number <- function(boundaries) {
  if (boundaries == "No Data" || is.na(boundaries)) {
    return(NA)  # Handle missing or "No Data" cases
  }
  
  # Extract numbers using regex
  numbers <- as.numeric(unlist(str_extract_all(boundaries, "\\d+")))
  
  # Return the maximum value
  #return(max(numbers, na.rm = TRUE))
  
  #get the second largest
  return(sort(numbers, decreasing = TRUE)[2])
}
get_largest_number <- function(boundaries) {
  if (boundaries == "No Data" || is.na(boundaries)) {
    return(NA)  # Handle missing or "No Data" cases
  }
  
  # Extract numbers using regex
  numbers <- as.numeric(unlist(str_extract_all(boundaries, "\\d+")))
  
  # Return the maximum value
  return(max(numbers, na.rm = TRUE))
  
}
calculate_snv_aft_NMD <- function(df, d_pen=50) {
  df <- df[order(df$cds_start),]  # Ensure sorting by cds_start
  exon.length <- (df$cds_end - df$cds_start) + 1
  
  if (length(exon.length) < 2) {
    # If there is only one exon, use cds_start as aft.NMD.ind
    return(data.frame(ensembl_transcript_id = unique(df$ensembl_transcript_id), aft.NMD.ind = df$cds_start[1]))
  }
  
  pen.length <- exon.length[length(exon.length) - 1]
  
  if (pen.length < d_pen) {
    aft.NMD.ind <- df$cds_start[length(exon.length) - 1]
  } else {
    aft.NMD.ind <- sum(exon.length[1:(length(exon.length) - 1)]) - d_pen
  }
  
  cds_length = max(df$cds_end)
  return(data.frame(ensembl_transcript_id = unique(df$ensembl_transcript_id), aft.NMD.ind = aft.NMD.ind, cds_length = cds_length))
}
calculate_snv_aft_NMD2 <- function(df, d_pen=50) {
  df <- df[order(df$cds_start),]  # Ensure sorting by cds_start
  exon.length <- (df$cds_end - df$cds_start) + 1
  
  if (length(exon.length) < 2) {
    # If there is only one exon, use cds_start as aft.NMD.ind
    return(data.frame(ensembl_transcript_id = unique(df$ensembl_transcript_id), aft.NMD.ind = df$cds_start[1]))
  }
  
  pen.length <- exon.length[length(exon.length) - 1]
  
  if (pen.length < d_pen) {
    aft.NMD.ind <- df$cds_start[length(exon.length) - 1]
  } else {
    aft.NMD.ind <- sum(exon.length[1:(length(exon.length) - 1)]) - d_pen
  }
  
  cds_length = max(df$cds_end)
  return(data.frame(ensembl_transcript_id = unique(df$ensembl_transcript_id), aft.NMD.ind = aft.NMD.ind, cds_length = cds_length))
}
calculate_fs_aft_NMD <- function(fs_transcript_id, type='plus2') {
  df2 = PTC_info[which(PTC_info$transcript == fs_transcript_id & PTC_info$type==type),]
  PTC_ind <- df2$PTC_loc
  PTC_status <- df2$PTC_NMDesc
  PTC_status_vector = as.logical(strsplit(PTC_status, ",")[[1]])
  can_PTC_status_vector = PTC_status_vector[1:(length(PTC_status_vector)/3)]
  PTC_ind_vector = as.numeric(unlist(strsplit(PTC_ind, ",")))
  can_change_points = which(diff(can_PTC_status_vector) != 0)
  can_change_points_loc = PTC_ind_vector[can_change_points]
  
  NMDesc.start = (can_change_points_loc+1)
  NMDesc.end = PTC_ind_vector[length(PTC_ind_vector)] 
  #cds_length = max(df$cds_end)
  return(data.frame(Transcript_ID = fs_transcript_id, 
                    NMDesc.start = NMDesc.start, 
                    NMDesc.end = NMDesc.end
                    #cds_length = cds_length
                    ))
}
key_to_transcript <- data.frame(
  Variant_Key = fs2@elementMetadata@listData[["id"]],
  Transcript_ID = fs2@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]],
  stringsAsFactors = FALSE
)

source("~/Desktop/idr/clean_code/get_snv_idr_match.R")
source("~/Desktop/idr/clean_code/get_fs_idr_match.R")
#generate new fasta file
snv.ind = which(res_plp_ptc@elementMetadata@listData[["type"]] == 'snv')
snv = res_plp_ptc[snv.ind]
plus1 = fs2[which(fs_NMD_result$type=='plus1')]
plus2 = fs2[which(fs_NMD_result$type=='plus2')]
#select NMDesc variants
snv_NMD = snv[snv@elementMetadata@listData[['res_aenmd']]@listData[['is_last']] | snv@elementMetadata@listData[['res_aenmd']]@listData[['is_penultimate']]]
snv_nonNMD = snv[!snv@elementMetadata@listData[['res_aenmd']]@listData[['is_last']] & !snv@elementMetadata@listData[['res_aenmd']]@listData[['is_penultimate']]]
plus1_NMD = plus1[plus1@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] | plus1@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]]
plus1_nonNMD = plus1[!plus1@elementMetadata@listData[['res_aenmd']]@listData[['is_last']] & !plus1@elementMetadata@listData[['res_aenmd']]@listData[['is_penultimate']]]
plus2_NMD = plus2[plus2@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] | plus2@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]]
#select variants corresponding to the gene list
snv_NMD_list = snv_NMD[match(snv_NMD@ranges@NAMES, snv_can_list, nomatch = 0) > 0]
snv_nonNMD_list = snv_nonNMD[match(snv_nonNMD@ranges@NAMES, snv_can_list, nomatch = 0) > 0]
plus1_NMD_list = plus1_NMD[match(plus1_NMD@ranges@NAMES, plus1_can_list, nomatch = 0) > 0]
plus1_nonNMD_list = plus1_nonNMD[match(plus1_nonNMD@ranges@NAMES, plus1_can_list, nomatch = 0) > 0]
plus2_NMD_list = plus2_NMD[match(plus2_NMD@ranges@NAMES, plus2_can_list, nomatch = 0) > 0]
snv_nf = create_fasta(snv_NMD_list,output_dir = "/Users/jxu14/Desktop/idr/snv_can")
plus1_nf = create_fasta(plus1_NMD_list,output_dir = "/Users/jxu14/Desktop/idr/plus1_can")
plus2_nf = create_fasta(plus2_NMD_list,output_dir = "/Users/jxu14/Desktop/idr/plus2_can")
#test if the overlap length >= 20*3bp
get_snv_idr_match("~/Downloads/long_idr_region/snv_control2_out.txt")
plus1_pre = get_fs_idr_match("/Users/jxu14/Desktop/idr/output/plus1_disorder_predictions.csv",type='plus1')
plus2_pre = get_fs_idr_match("/Users/jxu14/Desktop/idr/output/plus2_disorder_predictions.csv",type='plus2')