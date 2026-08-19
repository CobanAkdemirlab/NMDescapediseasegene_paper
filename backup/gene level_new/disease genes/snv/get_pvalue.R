get_pvalue = function(rds_name, outfilename){
  #saveDb(txdb, 'txdb.Ensembl105.sqlite')
  ###get DNAstring of UTR seqs
  res_p = readRDS(rds_name)
  #rds_name = 'snv_plp_ptc_nmdesc_can_filtered20260201.rds'
  txnames <- unique(res_p@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]])
  txnames.list <- list()
  
  for (i in 1: length(txnames)){
    print(i)
    #i=471
    txname=txnames[i]
    txnames.list[[i]] <- list ()
    txnames.list[[i]]$txname <- txname
    BM.ind <- which(BM.info$ensembl_transcript_id==txname)
    hgnc_symbol <- BM.info[BM.ind,'hgnc_symbol']
    txnames.list[[i]]$hgnc_symbol <- hgnc_symbol
    txnames.list[[i]]$is_canonical <- BM.info[BM.ind,'transcript_is_canonical']
    sub.ind <- which(res_p@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]]==txname) 
    res_p.ind <- res_p[sub.ind,]
    ptc.gen <- which(res_p.ind$res_aenmd$is_ptc==FALSE | res_p.ind$res_aenmd$is_single==TRUE)
    #delete the non-PTC variants and the single exon variants from ptc.gen
    if(length(ptc.gen) >0){
      res_p.ind <- res_p.ind[-ptc.gen,]
    }
    can.PTC.ind <- which(res_p.ind$res_aenmd$is_penultimate==TRUE | res_p.ind$res_aenmd$is_last==TRUE)
    if(length(can.PTC.ind)>0){
      can.PTC <- length (unique(sapply(1:length(names(res_p.ind[can.PTC.ind])), function(x)
      {
        res_p.ind$key[can.PTC.ind[x]]
      }))) } else {
        can.PTC <-0 
      }
  
    ens.ind <- which(all.df.sub$TXNAME==txname)
    tx.cor <- all.df.sub[ens.ind,]
    
    ### remove the non-coding exons
    tx.cor <- tx.cor[!is.na(tx.cor$CDSSTART),]
    #if this is not single exon
    if (nrow(tx.cor)>1) { 
      txnames.list[[i]]$txcor <- tx.cor
      txnames.list[[i]]$can.PTC <- can.PTC
      
      ex.length <- tx.cor$CDSEND-tx.cor$CDSSTART

      if(ex.length[length(ex.length)-1]<50) {
        
        NMD.lastexon <- ex.length[length(ex.length)] + ex.length[length(ex.length)-1]
      } else {
        
        NMD.lastexon <- ex.length[length(ex.length)] + 50
      }
      
      cds.start = 1  
      cds.end = sum(ex.length)
      NMDesc.start = cds.end - NMD.lastexon
      NMDesc.end = cds.end
      chrom = as.numeric(tx.cor$CDSCHROM[1])
      can.PTC.density <- can.PTC/get_syn_count(chrom,NMDesc.start,NMDesc.end)
      all.gr <-   GRanges(seqnames = c(tx.cor$CDSCHROM), ranges = IRanges(tx.cor$CDSSTART,tx.cor$CDSEND),
                          strand = tx.cor$CDSSTRAND[1])
          
      all.PTC <- length(unique(as.character(res_p.ind$key)))
      all.PTC.density <- all.PTC/get_syn_count(chrom,cds.start,cds.end)
      
      n_can <- get_syn_count(chrom, NMDesc.start, NMDesc.end)
      n_all <- get_syn_count(chrom, cds.start, cds.end)
      
      p0 <- all.PTC / n_all
      
      if (is.na(n_can) || n_can <= 0 || can.PTC > n_can || is.na(n_all) || n_all <= 0 || is.na(p0) || p0 < 0 || p0 > 1) {
        can.pvalue <- NA_real_
      } else {
        can.pvalue <- binom.test(can.PTC, n_can, p0, alternative = "greater")$p.value
      }
      
      txnames.list[[i]]$can.PTC.density <- can.PTC.density
      txnames.list[[i]]$all.PTC.density <- all.PTC.density
      txnames.list[[i]]$can.pvalue <- can.pvalue

    }
  }
  saveRDS(txnames.list,outfilename)
}
