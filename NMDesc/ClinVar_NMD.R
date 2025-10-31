setwd("/Users/qkelly/Desktop/autism/data")
res_clinvar_default_1 = readRDS("clinvar0716_1.rds")
res_clinvar_default_2 = readRDS("clinvar0716_2.rds")
res_clinvar_default_3 = readRDS("clinvar0716_3.rds")
res_clinvar_default = c(res_clinvar_default_1, res_clinvar_default_2, res_clinvar_default_3)
rm(res_clinvar_default_1, res_clinvar_default_2, res_clinvar_default_3)
res_clinvar_default = unlist(res_clinvar_default)
### get the ptc-generating variants
### look at the pathogenic ones

table(res_clinvar_default@elementMetadata@listData[["clinsig"]])

pat.ind <- grep('Pathogenic',res_clinvar_default@elementMetadata@listData[["clinsig"]])
ptc.gen <- which(res_clinvar_default@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]]==TRUE & res_clinvar_default@elementMetadata@listData[["res_aenmd"]]@listData[["is_single"]]==FALSE)

pat_ptc <- intersect(pat.ind,ptc.gen)
clinvar_pat_ptc = res_clinvar_default[pat_ptc]

### escape ones
length(which(clinvar_pat_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] | clinvar_pat_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] | clinvar_pat_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]] | clinvar_pat_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]))

### for each transcript

#txnames <- unique(unlist(sapply(1:length(names(res_clinvar_default)),function(x){
  #strsplit(names(res_clinvar_default)[x],'\\|')[[1]][1]
  
#})))

txnames <- unique(names(clinvar_pat_ptc))

txnames.list <- list()
source('/Users/qkelly/Desktop/clinvar/get_statistics.R')
txnames.list = get_statistics(txnames, clinvar_pat_ptc)

outfilename <- paste('~/Downloads/ClinVar_NMD.rds',sep='')
saveRDS(txnames.list,outfilename)

genes <- vector()
allgenes <- vector()

  
  outfilename = '/Users/qkelly/Desktop/clinvar/Clinvar_0820NMD.rds'
  l <- readRDS(outfilename)
  p_indexes <- list()
  for ( p in 1:length(l))
  {
    p_indexes[[p]] <- (as.numeric(l[[p]]$css.pvalue)>=0.95 | as.numeric(l[[p]]$can.pvalue)>=0.95 | as.numeric(l[[p]]$long.pvalue)>=0.87) & l[[p]]$is_canonical
    }
  
  can.pvalues <-  sapply( 1:length(l), function(p)
  {
    as.numeric(l[[p]]$can.pvalue)
    #p_indexes[[p]] <- grep('ADGRG1',l[[p]]$hgnc_symbol)
    #p_indexes[[p]] <- grep('ENST00000378888)',l[[p]]$txname)
    
  })
  quantile(unlist(can.pvalues),0.90,na.rm=TRUE)
  
  css.pvalues <-  sapply( 1:length(l), function(p)
  {
    as.numeric(l[[p]]$css.pvalue)
    #p_indexes[[p]] <- grep('',l[[p]]$hgnc_symbol)
    #p_indexes[[p]] <- grep('ENST00000378888)',l[[p]]$txname)
    
  })
  quantile(unlist(css.pvalues),0.90,na.rm=TRUE)
  
  long.pvalues <-  sapply( 1:length(l), function(p)
  {
    as.numeric(l[[p]]$long.pvalue)
    #p_indexes[[p]] <- grep('ADGRG1',l[[p]]$hgnc_symbol)
    #p_indexes[[p]] <- grep('ENST00000378888)',l[[p]]$txname)
    
  })
  quantile(unlist(long.pvalues),0.90,na.rm=TRUE)
  
  all.genes <- unlist(sapply( 1:length(l), function(x) {
    
    l[[x]]$hgnc_symbol
  }))
  
  #ind <- which(sapply( p_indexes, function(x) length(x) > 0 )==TRUE)
  ind.1 <- which(p_indexes==TRUE)
  
  new.genes <- unlist(sapply( ind.1, function(x) {
    
    l[[x]]$hgnc_symbol
  }))
  
  genes <- c(new.genes,genes)
  #allgenes <- c(all.genes, allgenes)
  writeLines(genes,'~/Desktop/clinvar/ClinVar_pathogenic_NMDescapeneriched.genes.txt')
  writeLines(all.genes,'~/Desktop/clinvar/ClinVar_pathogenic_allgenes.txt')
  
  
  ### merge with OMIM and PLI and LOEUFF scores
  
 Lof_metrics <- read.table('~/Downloads/gnomad.v2.1.1.lof_metrics.by_gene.txt',header=TRUE,sep='\t')
 genes.mat <- NULL
 genes.mat$gene <- unique(all.genes)
 genes.mat.1 <- merge(genes.mat,Lof_metrics,by='gene')
 genes.mat.1$NMDescape.enr <- sapply(1:nrow(genes.mat.1),function(x){
   genes.mat.1$gene[x]%in%unique(genes)
   
   
 })
 
 ind.enr <- which(genes.mat.1$NMDescape.enr=='TRUE')
 boxplot(genes.mat.1[ind.enr,'pLI'])
 
 ind <- which(genes.mat.1$pLI<0.9 | genes.mat.1$oe_lof_upper>=0.6 )
 genes.mat.1$tolerant <- sapply(1:nrow(genes.mat.1),function(x){
   x%in%ind
   
   
 })
 
 gnomAD_NMDescape_enr <-unique(readLines('~/Desktop/gnomAD_NMD/gnomAD_NMDescape_depleted_genes.txt'))
 
 genes.mat.1$gnomAD.NMDescape.enr <-  sapply(1:nrow(genes.mat.1),function(x){
   genes.mat.1$gene[x]%in%unique(gnomAD_NMDescape_enr)
   
   
 })
 ## THEN ADD THE OMIM INFORMATION
 
 OMIM <- read.csv('~/Downloads/OMIM/genemap2.txt',sep='\t',stringsAsFactors=FALSE,header=TRUE,skip=3)
 mim2gene <- read.csv('~/Downloads/OMIM/mim2gene.txt',sep='\t',stringsAsFactors=FALSE,header=TRUE,skip=4)
 
 #mim2gene$MIM.Number <- mim2gene$X..MIM.Number
 colnames(mim2gene)[which(names(mim2gene)=='X..MIM.Number')] <- 'MIM.Number'
 merge_OMIM_mim2gene <- merge(mim2gene,OMIM,by='MIM.Number')
 #merge_OMIM_mim2gene$Raw_gene_name <- merge_OMIM_mim2gene$Approved.Gene.Symbol..HGNC.
 colnames(merge_OMIM_mim2gene)[which(names(merge_OMIM_mim2gene)=='Approved.Gene.Symbol..HGNC.')] <- 'Raw_gene_name'
 genes.mat.1$Raw_gene_name <- genes.mat.1$gene
 
 data.new.merged <- merge(genes.mat.1,merge_OMIM_mim2gene,by='Raw_gene_name',all.x=TRUE)
 #data.new.merged$OMIM_disorder <- data.new.merged$Phenotypes
 colnames(data.new.merged)[which(names(data.new.merged)=='Phenotypes')] <- 'OMIM_disorder'
   
 domi.ind <- grep('dominant',data.new.merged$OMIM_disorder)
 rec.ind <- grep('recessive',data.new.merged$OMIM_disorder)
 length(intersect(domi.ind,rec.ind))
 domi.number <- length(domi.ind)-length(intersect(domi.ind,rec.ind))
 rec.number <- length(rec.ind)-length(intersect(domi.ind,rec.ind))
 data.new.merged$inheritance <- rep('Unknow',nrow(data.new.merged))
 data.new.merged$inheritance[setdiff(domi.ind,intersect(domi.ind,rec.ind))] <- 'AD'
 data.new.merged$inheritance[setdiff(rec.ind,intersect(rec.ind,domi.ind))] <- 'AR'
 data.new.merged$inheritance[intersect(domi.ind,rec.ind)] <- 'AD/AR'
 
 
 domi.mat <- data.new.merged[domi.ind,]
 
  
 write.table( data.new.merged,'~/Desktop/CLinvar_genes_NMDdata.txt',sep='\t')
