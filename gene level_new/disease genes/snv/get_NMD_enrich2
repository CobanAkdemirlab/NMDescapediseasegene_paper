get_NMD_enrichment = function(rds_name,ADfilter = 'off',p_cutoff=0.1,filter_type){
  l <- readRDS(rds_name)
  if(ADfilter == 'on'){
    genes.mat <- NULL
    genes.mat$gene = data.frame('gene' = unlist(sapply(l, function(x) x$hgnc_symbol)))
    genes.mat.1 <- merge(genes.mat,Lof_metrics,by='gene')
    colnames(mim2gene)[which(names(mim2gene)=='X..MIM.Number')] <- 'MIM.Number'
    merge_OMIM_mim2gene <- merge(mim2gene,OMIM,by='MIM.Number')
    colnames(merge_OMIM_mim2gene)[which(names(merge_OMIM_mim2gene)=='Approved.Gene.Symbol..HGNC.')] <- 'Raw_gene_name'
    genes.mat.1$Raw_gene_name <- genes.mat.1$gene
    
    data.new.merged <- merge(genes.mat.1,merge_OMIM_mim2gene,by='Raw_gene_name',all.x=TRUE)
    colnames(data.new.merged)[which(names(data.new.merged)=='Phenotypes')] <- 'OMIM_disorder'
    
    domi.ind <- grep('dominant',data.new.merged$OMIM_disorder)
    rec.ind <- grep('recessive',data.new.merged$OMIM_disorder)
    
    length(intersect(domi.ind,rec.ind))
    domi.number <- length(domi.ind)-length(intersect(domi.ind,rec.ind))
    rec.number <- length(rec.ind)-length(intersect(domi.ind,rec.ind))
    data.new.merged$inheritance <- rep('Unknown',nrow(data.new.merged))
    data.new.merged$inheritance[setdiff(domi.ind,intersect(domi.ind,rec.ind))] <- 'AD'
    data.new.merged$inheritance[setdiff(rec.ind,intersect(rec.ind,domi.ind))] <- 'AR'
    data.new.merged$inheritance[intersect(domi.ind,rec.ind)] <- 'AD/AR'
    AD_gene = data.new.merged$gene[which(data.new.merged$inheritance == 'AD')]
    l2 = l[vapply(l, function(x) !is.null(x$hgnc_symbol) && !is.na(x$hgnc_symbol) && x$hgnc_symbol %in% AD_gene, logical(1))]
  }

  p_indexes <- list()
  genes <- vector()

  can.pvalues <-  sapply( 1:length(l), function(p)
  {
    as.numeric(l[[p]]$can.pvalue)
  })
  can_cut = quantile(unlist(can.pvalues),p_cutoff,na.rm=TRUE)
  cat('can_cut:', can_cut)
  
  css.pvalues <-  sapply( 1:length(l), function(p)
  {
    as.numeric(l[[p]]$css.pvalue)
  })
  css_cut = quantile(unlist(css.pvalues),p_cutoff,na.rm=TRUE)
  cat('\n css_cut:', css_cut)
  
  long.pvalues <-  sapply( 1:length(l), function(p)
  {
    as.numeric(l[[p]]$long.pvalue)
  })
  long_cut = quantile(unlist(long.pvalues),p_cutoff,na.rm=TRUE)
  cat('\n long_cut:', long_cut)
  
  trig.pvalues <-  sapply( 1:length(l), function(p)
  {
    as.numeric(l[[p]]$trig.pvalue)
  })
  trig_cut = quantile(unlist(trig.pvalues),p_cutoff,na.rm=TRUE)
  cat('\n trig_cut:', trig_cut)
  
  filter_type = 'can'
  for (p in 1:length(l)) {
    if(filter_type == 'all'){
      p_indexes[[p]] <- (as.numeric(l[[p]]$css.pvalue)>= css_cut | as.numeric(l[[p]]$can.pvalue) <= 0.05 | as.numeric(l[[p]]$long.pvalue)>=long_cut) & l[[p]]$is_canonical
      #p_indexes[[p]] <- (as.numeric(l[[p]]$css.pvalue)>= 0.95 | as.numeric(l[[p]]$can.pvalue)>= 0.95 | as.numeric(l[[p]]$long.pvalue)>=0.95) & l[[p]]$is_canonical
    } else if(filter_type == 'can'){
      p_indexes[[p]] <- (as.numeric(l[[p]]$can.pvalue) <=  0.05) 
      #p_indexes[[p]] <- (as.numeric(l[[p]]$can.pvalue) >= 0.95) & l[[p]]$is_canonical
    } else if(filter_type == 'css'){
      p_indexes[[p]] <- (as.numeric(l[[p]]$css.pvalue) >= css_cut) & l[[p]]$is_canonical
      #p_indexes[[p]] <- (as.numeric(l[[p]]$css.pvalue) >= 0.95) & l[[p]]$is_canonical
    } else if(filter_type == 'long'){
      p_indexes[[p]] <- (as.numeric(l[[p]]$long.pvalue) >= long_cut) & l[[p]]$is_canonical
      #p_indexes[[p]] <- (as.numeric(l[[p]]$long.pvalue) >= 0.95) & l[[p]]$is_canonical
      }else if(filter_type == 'trig'){
        p_indexes[[p]] <- (as.numeric(l[[p]]$trig.pvalue) >= trig_cut) & l[[p]]$is_canonical
      }
  }
 
  ind.1 <- which(p_indexes==TRUE)
  new.trs <- unlist(sapply( ind.1, function(x) {
    l[[x]]$txname
  }))
  #use getBM, transform txname to hgnc_symbol
  genes = getBM(
    attributes = c("external_gene_name", "ensembl_transcript_id", "transcript_is_canonical"),
    filters = "ensembl_transcript_id",
    values = unique(new.trs),
    mart = ensembl
  )
  genes$gene = genes$external_gene_name
  genes <- c(new.genes,genes)
  #from genes select AD ones
  genes.mat.1 <- merge(genes,Lof_metrics,by='gene')
  colnames(mim2gene)[which(names(mim2gene)=='X..MIM.Number')] <- 'MIM.Number'
  merge_OMIM_mim2gene <- merge(mim2gene,OMIM,by='MIM.Number')
  colnames(merge_OMIM_mim2gene)[which(names(merge_OMIM_mim2gene)=='Approved.Gene.Symbol..HGNC.')] <- 'Raw_gene_name'
  genes.mat.1$Raw_gene_name <- genes.mat.1$gene
  
  data.new.merged <- merge(genes.mat.1,merge_OMIM_mim2gene,by='Raw_gene_name',all.x=TRUE)
  colnames(data.new.merged)[which(names(data.new.merged)=='Phenotypes')] <- 'OMIM_disorder'
  
  domi.ind <- grep('dominant',data.new.merged$OMIM_disorder)
  rec.ind <- grep('recessive',data.new.merged$OMIM_disorder)
  
  length(intersect(domi.ind,rec.ind))
  domi.number <- length(domi.ind)-length(intersect(domi.ind,rec.ind))
  rec.number <- length(rec.ind)-length(intersect(domi.ind,rec.ind))
  data.new.merged$inheritance <- rep('Unknown',nrow(data.new.merged))
  data.new.merged$inheritance[setdiff(domi.ind,intersect(domi.ind,rec.ind))] <- 'AD'
  data.new.merged$inheritance[setdiff(rec.ind,intersect(rec.ind,domi.ind))] <- 'AR'
  data.new.merged$inheritance[intersect(domi.ind,rec.ind)] <- 'AD/AR'
  AD_gene = data.new.merged$gene[which(data.new.merged$inheritance == 'AD')]
  
  out_name = paste(gsub('.rds','_',rds_name),'NMDenriched_snv_',filter_type,'.txt',sep = '')
  write.csv(AD_gene,out_name)
}
