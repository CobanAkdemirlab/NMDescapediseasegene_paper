#input idr loc file path and gene list type
get_snv_idr_match = function(input_path="~/Downloads/long_idr_region/snv_control2_out.txt"){
  idr_loc = read_delim(input_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # from idr loc, get c-start and c-end loc
  idr_loc = idr_loc %>%
    filter(!is.na(Disordered_Domain_Boundaries) & Disordered_Domain_Boundaries != "[]") %>%  # Remove empty disorder domains
    mutate(C_end_IDR = sapply(Disordered_Domain_Boundaries, get_largest_number),
           C_start_IDR = sapply(Disordered_Domain_Boundaries, get_second_largest_number))
  
  
  #for each uniprot id, find gene id
  idr_loc$hgnc_symbol = getBM(attributes = c('hgnc_symbol','uniprotswissprot'),
                                   filters = 'uniprotswissprot',
                                   values = idr_loc$UniProt_ID,
                                   mart = ensembl)$hgnc_symbol
  
  #for each gene, find it's canonical transcript
  gene2tr = getBM(attributes = c('hgnc_symbol','ensembl_transcript_id','transcript_is_canonical'),
                                      filters = 'hgnc_symbol',
                                      values = idr_loc$hgnc_symbol,
                                      mart = ensembl)
  gene2tr = gene2tr[which(gene2tr$transcript_is_canonical== 1),]
  #merge idr_loc with gene2tr by hgnc_symbol, keep one row for each gene
  idr_loc2 = merge(gene2tr,idr_loc,by = 'hgnc_symbol')
  idr_loc2 = idr_loc2[!duplicated(idr_loc2$hgnc_symbol),]
  #from transcript id, get cds information
  tr2cds = getBM(attributes = c('ensembl_transcript_id','cds_start','cds_end'),
                          filters = 'ensembl_transcript_id',
                          values = idr_loc2$ensembl_transcript_id,
                          mart = ensembl)
  tr2cds <- tr2cds[!is.na(tr2cds$cds_start),]
  
  aft_NMD_results <- tr2cds %>%
    group_by(ensembl_transcript_id) %>%
    group_split() %>%
    lapply(calculate_snv_aft_NMD, d_pen) %>%
    bind_rows() %>%
    left_join(tr2cds, by = "ensembl_transcript_id")  %>%
    select(ensembl_transcript_id, aft.NMD.ind, cds_length) %>%
    distinct(ensembl_transcript_id, .keep_all = TRUE) 

  #merge with idr_loc
  idr_NMD = merge(aft_NMD_results,idr_loc2,by = 'ensembl_transcript_id')
  #match NMD and idr
  idr_NMD$match_start = (idr_NMD$aft.NMD.ind - idr_NMD$C_start_IDR*3) >= 20*3
  idr_NMD$match_end = idr_NMD$C_end_IDR*3 > idr_NMD$aft.NMD.ind
  idr_NMD$match = idr_NMD$match_start & idr_NMD$match_end
  #return the count of matched idr and total gene count
  return(data.frame(matched_idr = sum(idr_NMD$match,na.rm=T), total_gene = nrow(idr_NMD)))
  
}
