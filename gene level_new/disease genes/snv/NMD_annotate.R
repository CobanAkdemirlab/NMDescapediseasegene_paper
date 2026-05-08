NMD_annotate = function(raw_name, rds_name){
  clinvarClinsig <- read.delim(raw_name)
  clinvarClinsig = resetID_clinvar_20241015_Clnsig_filtered
  #Get total variants in Clinvar:
  Variants_in_the_original_dataset <- nrow(clinvarClinsig)
  #remove the NA in clinvarClinsig$CHROM
  clinvarClinsig <- clinvarClinsig[!is.na(clinvarClinsig$CHROM),]
  
  #Make a GRange object of all the clinvar variants without filtering.
  clinvar_gr_total  <- GenomicRanges::GRanges(clinvarClinsig$CHROM,
                                              IRanges::IRanges(clinvarClinsig$POS,
                                                               clinvarClinsig$POS))

  clinvar_gr_total$ref     <- clinvarClinsig$REF |> Biostrings::DNAStringSet()
  clinvar_gr_total$alt     <- clinvarClinsig$ALT |> Biostrings::DNAStringSet()
  clinvar_gr_total$id      <- clinvarClinsig$ID
  clinvar_gr_total$clinsig      <- clinvarClinsig$CLNSIG
  
  #Retrieve the base AENMD transcript set information: 
  AENMD_base_exonset <- future::value(aenmd:::._EA_exn_grl)
  
  #unlist and sort
  AENMD_base_exonset_unlist_srtd <- sort(GenomeInfoDb::sortSeqlevels(unlist(AENMD_base_exonset)))
  
  #find overlap between clinvar variants and the exons of the transcript set:
  ov <- GenomicRanges::findOverlaps(clinvar_gr_total, AENMD_base_exonset_unlist_srtd)
  
  #find the number of unique hits, such that a variant can only be counted once:
  Variants_coding_our_tx_set <- length(unique(S4Vectors::queryHits(ov)))
  
  #make a data.frame:
  Total_var_stats <- data.frame(functional_class = c("Variants_in_the_original_dataset" , "Variants_coding_our_tx_set"),
                                value = c(Variants_in_the_original_dataset, Variants_coding_our_tx_set), percent = c(NA, Variants_coding_our_tx_set/Variants_in_the_original_dataset*100)
  )
  NA_ind <- which(clinvarClinsig$ALT != "." & !stringr::str_detect(clinvarClinsig$ALT, "R|Y|S|W|K|M|B|D|H|V|N"))
  
  #now put together the GRanges without NA to put into AENMD
  #we do it this way to include the clinsig annotation
  clinvar_gr  <- GenomicRanges::GRanges(clinvarClinsig$CHROM[NA_ind],
                                        IRanges::IRanges(clinvarClinsig$POS[NA_ind],
                                                         clinvarClinsig$POS[NA_ind]))
  clinvar_gr$ref     <- clinvarClinsig$REF[NA_ind] |> Biostrings::DNAStringSet()
  clinvar_gr$alt     <- clinvarClinsig$ALT[NA_ind] |> Biostrings::DNAStringSet()
  clinvar_gr$id      <- clinvarClinsig$ID[NA_ind]
  clinvar_gr$clinsig      <- clinvarClinsig$CLNSIG[NA_ind]
  
  #Process variants
  clinvar_proc <- process_variants(clinvar_gr)
  res_clinvar_default <- aenmd::annotate_nmd(clinvar_proc)
  saveRDS(res_clinvar_default, file = rds_name)
}