get_pvalue = function(rds_name, outfilename){
  #saveDb(txdb, 'txdb.Ensembl105.sqlite')
  ###get DNAstring of UTR seqs
  chr1 = readRDS(rds_name)
  
  txnames <- unique(names(chr1))
  txnames.list <- list()
  
  for (i in 1: length(txnames)){
    #print(i)
    txname=txnames[i]
    txnames.list[[i]] <- list ()
    txnames.list[[i]]$txname <- txname
    BM.ind <- which(BM.info$ensembl_transcript_id==txname)
    hgnc_symbol <- BM.info[BM.ind,'hgnc_symbol']
    txnames.list[[i]]$hgnc_symbol <- hgnc_symbol
    txnames.list[[i]]$is_canonical <- BM.info[BM.ind,'transcript_is_canonical']
    sub.ind <- which(names(chr1)==txname)
    chr1.ind <- chr1[sub.ind,]
    ptc.gen <- which(chr1.ind$res_aenmd$is_ptc==FALSE | chr1.ind$res_aenmd$is_single==TRUE | chr1.ind$filter!='PASS')
    if(length(ptc.gen) >0){
      chr1.ind <- chr1.ind[-ptc.gen,]
    }
    can.PTC.ind <- which(chr1.ind$res_aenmd$is_penultimate==TRUE | chr1.ind$res_aenmd$is_last==TRUE)
    if(length(can.PTC.ind)>0){
      can.PTC <- length (unique(sapply(1:length(names(chr1.ind[can.PTC.ind])), function(x)
      {
        chr1.ind$key[can.PTC.ind[x]]
      }))) } else {
        can.PTC <-0 
      }
    css.PTC.ind <- which(chr1.ind$res_aenmd$is_cssProximal==TRUE)
    if(length(css.PTC.ind)>0){
      css.PTC <- length (unique(sapply(1:length(names(chr1.ind[css.PTC.ind])), function(x)
      {
        chr1.ind$key[css.PTC.ind[x]]
        
      }))) } else {
        css.PTC <-0 
      }
    long.PTC.ind <- which(chr1.ind$res_aenmd$is_407plus==TRUE)
    if(length(long.PTC.ind>0)){
      long.PTC <- length (unique(sapply(1:length(names(chr1.ind[long.PTC.ind])), function(x)
      {
        chr1.ind$key[long.PTC.ind[x]]
        
      }))) } else {
        long.PTC <-0 
      }
    ens.ind <- which(all.df.sub$TXNAME==txname)
    tx.cor <- all.df.sub[ens.ind,]
    
    ### remove the non-coding exons
    tx.cor <- tx.cor[!is.na(tx.cor$CDSSTART),]
    if (nrow(tx.cor)>1) { #if this is not single exon
      txnames.list[[i]]$txcor <- tx.cor
      txnames.list[[i]]$can.PTC <- can.PTC
      txnames.list[[i]]$css.PTC <- css.PTC
      txnames.list[[i]]$long.PTC <- long.PTC
      
      long.ex.ind <- which(tx.cor$CDSEND-tx.cor$CDSSTART>=407)
      ex.length <- tx.cor$CDSEND-tx.cor$CDSSTART
      NMD.longexon <- sum(ex.length[long.ex.ind])
      NMD.CSS <- 150
      if(ex.length[length(ex.length)-1]<50) {
        
        NMD.lastexon <- ex.length[length(ex.length)] + ex.length[length(ex.length)-1]
      } else {
        
        NMD.lastexon <- ex.length[length(ex.length)] + 50
      }
      can.PTC.density <- can.PTC/NMD.lastexon
      CSS.PTC.density <- css.PTC/NMD.CSS
      if(NMD.longexon>0){
        long.PTC.density <- long.PTC/NMD.longexon
      } else {
        
        long.PTC.density <-NA
        
      }
      txnames.list[[i]]$can.PTC.density <- can.PTC.density
      txnames.list[[i]]$CSS.PTC.density <- CSS.PTC.density
      txnames.list[[i]]$long.PTC.density <- long.PTC.density
     
      ## get the remaining region
      
      ### construct a grange list object for each transcript's CSS, long-exon, last and final exon region
      
      ### for each transcript's CSS
      
      if (length(which(cumsum(ex.length)>=150)) >0){
        css.ex.ind <- min(which(cumsum(ex.length)>=150))
        if(css.ex.ind==1 & tx.cor$CDSSTRAND[1]=='-') {
          css.gr <- tx.cor[which(tx.cor$EXONRANK==1),]
          css.gr$CDSSTART <- css.gr$CDSEND-150+1
        } else if (css.ex.ind==1 & tx.cor$CDSSTRAND[1]=='+') {
          css.gr <- tx.cor[which(tx.cor$EXONRANK==1),]
          css.gr$CDSEND <- css.gr$CDSSTART+150-1
        } else if (css.ex.ind>1 & tx.cor$CDSSTRAND[1]=='-') {
          css.gr <- tx.cor[1:css.ex.ind,]
          dif <- 150-cumsum(ex.length)[css.ex.ind-1]
          css.gr$CDSSTART[css.ex.ind] <- css.gr$CDSEND[css.ex.ind]-dif+1
        } else if (css.ex.ind>1 & tx.cor$CDSSTRAND[1]=='+') {
          css.gr <- tx.cor[1:css.ex.ind,]
          dif <- 150-cumsum(ex.length)[css.ex.ind-1]
          css.gr$CDSEND[css.ex.ind] <- css.gr$CDSSTART[css.ex.ind]+dif-1
        } 
      } else {
        css.ex.ind <- which(cumsum(ex.length)>=150)
        css.gr <-tx.cor[css.ex.ind,]
      }
      css.long.ind <- which(ex.length>=407)
      long.gr <- tx.cor[css.long.ind,]
      
      if(ex.length[length(ex.length)-1]<50) {
        last.gr <- tx.cor[(length(ex.length)-1):length(ex.length),] 
      } else if (ex.length[length(ex.length)-1]>=50 & tx.cor$CDSSTRAND[1]=='-') {
        last.gr <- tx.cor[(length(ex.length)-1):length(ex.length),] 
        last.gr$CDSSTART[1] <- last.gr$CDSEND[1]-50+1
      } else if (ex.length[length(ex.length)-1]>=50 & tx.cor$CDSSTRAND[1]=='+') {
        last.gr <- tx.cor[(length(ex.length)-1):length(ex.length),] 
        last.gr$CDSEND[1] <- last.gr$CDSSTART[1]+50-1
      }
      
      txnames.list[[i]]$css.gr <- css.gr
      txnames.list[[i]]$long.gr <- long.gr
      txnames.list[[i]]$last.gr <- last.gr
      
      NMD.esc.gr <-
        GRanges(seqnames = c(css.gr$CDSCHROM,long.gr$CDSCHROM,last.gr$CDSCHROM), ranges = c(IRanges(css.gr$CDSSTART,css.gr$CDSEND),IRanges(long.gr$CDSSTART,long.gr$CDSEND),IRanges(last.gr$CDSSTART,last.gr$CDSEND)),
                strand = tx.cor$CDSSTRAND[1], subtype=c(rep('css',nrow(css.gr)),rep('long',nrow(long.gr)),rep('last',nrow(last.gr))))
      all.gr <-   GRanges(seqnames = c(tx.cor$CDSCHROM), ranges = IRanges(tx.cor$CDSSTART,tx.cor$CDSEND),
                          strand = tx.cor$CDSSTRAND[1])
      NMD.trig.gr <- setdiff(all.gr,NMD.esc.gr)
      txnames.list[[i]]$NMD.trig.gr <- NMD.trig.gr
      trig.PTC.ind <- which(chr1.ind$res_aenmd$is_penultimate==FALSE & chr1.ind$res_aenmd$is_last==FALSE & chr1.ind$res_aenmd$is_cssProximal==FALSE & chr1.ind$res_aenmd$is_407plus==FALSE )
      if(length(trig.PTC.ind)>0){
        trig.PTC <- length (unique(sapply(1:length(names(chr1.ind[trig.PTC.ind])), function(x)
        {
          #strsplit(names(chr1.ind[trig.PTC.ind])[x],'\\|')[[1]][2]
          chr1.ind[trig.PTC.ind[x],]$key
          
        }))) } else {
          trig.PTC <-0 
        }
            txnames.list[[i]]$trig.PTC <- trig.PTC
      trig.PTC.density <- trig.PTC/ sum(width(NMD.trig.gr))
      
      all.PTC <- length (unique(sapply(1:length(names(chr1.ind)), function(x)
      {
        chr1.ind[x,]$key
        
      }))) 
      
      all.PTC.density <- length(all.PTC)/sum(ex.length)
      txnames.list[[i]]$trig.PTC.density <- trig.PTC.density
      if (nrow(css.gr) >0) {
        css.pvalue <- binom.test(css.PTC,sum(css.gr$CDSEND-css.gr$CDSSTART),all.PTC/sum(ex.length),alternative='less')$p.value
      } else {
        css.pvalue='NA'
      }
      can.pvalue <- binom.test(can.PTC,sum(last.gr$CDSEND-last.gr$CDSSTART),all.PTC/sum(ex.length),alternative='less')$p.value
      if (nrow(long.gr) >0) {
        long.pvalue <- binom.test(long.PTC,sum(long.gr$CDSEND-long.gr$CDSSTART),all.PTC/sum(ex.length),alternative='less')$p.value
        
      } else {
        long.pvalue='NA'
        
      }
      txnames.list[[i]]$all.PTC.density <- all.PTC.density
      txnames.list[[i]]$css.pvalue <- css.pvalue
      txnames.list[[i]]$can.pvalue <- can.pvalue
      txnames.list[[i]]$long.pvalue <- long.pvalue
    }
  }
  saveRDS(txnames.list,outfilename)
}
