library(GenomicFeatures)
library(VariantAnnotation)
library(AnnotationDbi)
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

txdb <- makeTxDbFromEnsembl("Homo sapiens", release=105)

#saveDb(txdb, 'txdb.Ensembl105.sqlite')

ensgene <- txdb

seqlevels(ensgene) <- paste('chr',seqlevels(ensgene),sep='')

chr.list <- c(paste('chr',1:22,sep=''))

ensgene.sub <- keepSeqlevels(ensgene,chr.list)
ensgene <- ensgene.sub

cds_seqs <- GenomicFeatures::extractTranscriptSeqs(genome, cdsBy(ensgene, by="tx",use.names=TRUE))
threeutr.grange <-threeUTRsByTranscript(ensgene, use.names=TRUE)
fiveutr.grange<-fiveUTRsByTranscript(ensgene, use.names=TRUE)
introns.grange<- intronsByTranscript(ensgene, use.names=TRUE)

###get DNAstring of UTR seqs
five_utr_seqs <- extractTranscriptSeqs(Hsapiens, fiveutr.grange) #what is hsapiens?
three_utr_seqs <- extractTranscriptSeqs(Hsapiens, threeutr.grange)

keys <- names(cds_seqs)
cols <- columns(ensgene)
all.df <- select(ensgene, keys = keys, columns = cols, keytype="TXNAME")
## this has length of 304803 so most contain NAs or duplicates?
##cds_seqs length is 28856
all.df.sub <- all.df[which(all.df$TXNAME!='NA'),]
### do it for each chromosome
library(biomaRt)
mart <- useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
BM.info <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","hgnc_symbol","transcript_is_canonical"),mart=mart)



setwd("/Users/qkelly/Desktop/clinvar/0819")
res_clinvar_default_1 = readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/桌面 - q的MacBook Pro/clinvar/0819/clinvar0819_1.rds")
res_clinvar_default_2 = readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/桌面 - q的MacBook Pro/clinvar/0819/clinvar0819_2.rds")
res_clinvar_default_3 = readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/桌面 - q的MacBook Pro/clinvar/0819/clinvar0819_3.rds")
res_clinvar_default = c(res_clinvar_default_1, res_clinvar_default_2, res_clinvar_default_3)
rm(res_clinvar_default_1, res_clinvar_default_2, res_clinvar_default_3)

#res_clinvar_default = readRDS("~/Desktop/clinvar/aenmdDefaultSettings_Clinvar20230616_ClnsigClinDN_ensmbl105_results.rds")

res_clinvar_default = unlist(res_clinvar_default)

### get the ptc-generating variants
### look at the pathogenic ones

table(res_clinvar_default@elementMetadata@listData[["clinsig"]])

pat.ind <- grep('Pathogenic',res_clinvar_default@elementMetadata@listData[["clinsig"]])
vus.ind <- grep('Uncertain_significance',res_clinvar_default@elementMetadata@listData[["clinsig"]])

ptc.gen <- which(res_clinvar_default@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]]==TRUE & res_clinvar_default@elementMetadata@listData[["res_aenmd"]]@listData[["is_single"]]==FALSE)
table(res_clinvar_default[ptc.gen]@elementMetadata@listData[["clinsig"]])

pat_ptc <- intersect(pat.ind,ptc.gen)
vus_ptc <- intersect(vus.ind,ptc.gen)
clinvar_pat_ptc = res_clinvar_default[pat_ptc]
clinvar_vus_ptc = res_clinvar_default[vus_ptc]
clinvar_pat = res_clinvar_default[pat.ind]
clinvar_vus = res_clinvar_default[vus.ind]

### escape ones in Pathogenic
#can: last+ pen
p_can = length(which(clinvar_pat_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] | clinvar_pat_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]))
p_can_per = p_can/length(clin)
#all rules
p_all = length(which(clinvar_pat_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] | clinvar_pat_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] | clinvar_pat_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]] | clinvar_pat_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]))
p_all_per = p_all/length(pat.ind)
### escape ones in VUS
#can: last+ pen
v_can = length(which(clinvar_vus_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] | clinvar_vus_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]))
v_can_per = v_can/length(vus.ind)
#all rules
#v_all = length(which(clinvar_vus@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] | clinvar_vus@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] | clinvar_vus@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]] | clinvar_vus@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]))
v_all = length(which(clinvar_vus_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]] | clinvar_vus_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] | clinvar_vus_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]] | clinvar_vus_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]))
v_all_per = v_all/length(clinvar_vus_ptc)
v_all_per = v_all/length(vus.ind)

#txnames <- unique(unlist(sapply(1:length(names(res_clinvar_default)),function(x){
  #strsplit(names(res_clinvar_default)[x],'\\|')[[1]][1]
  
#})))
chr1 = clinvar_vus_ptc
txnames <- unique(names(chr1))
txnames.list <- list()

for (i in 1: length(txnames)){
  print(i)
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
  if (nrow(tx.cor)>1) {
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

outfilename <- 'ClinVar_vus_NMD0906.rds'
saveRDS(txnames.list,outfilename)

genes <- vector()
allgenes <- vector()

  
  outfilename = 'Clinvar_p_NMD0906.rds'
  l <- readRDS(outfilename)
  p_indexes <- list()
  for (p in 1:length(l)) {
     p_indexes[[p]] <- (as.numeric(l[[p]]$css.pvalue)>=0.95 | as.numeric(l[[p]]$can.pvalue)>=0.95 | as.numeric(l[[p]]$long.pvalue)>=0.87) & l[[p]]$is_canonical
    
    # only can
    #p_indexes[[p]] <- (as.numeric(l[[p]]$can.pvalue) >= 0.95) & l[[p]]$is_canonical
    
    # only css rule  
     #p_indexes[[p]] <- (as.numeric(l[[p]]$css.pvalue) >= 0.95) & l[[p]]$is_canonical
    
    # only long rule
     #p_indexes[[p]] <- (as.numeric(l[[p]]$long.pvalue) >= 0.87) & l[[p]]$is_canonical
    
    # p_indexes[[p]] <- grep('NOVA2', l[[p]]$hgnc_symbol)
    # p_indexes[[p]] <- grep('ENST00000378888)', l[[p]]$txname)
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
  writeLines(genes,'~/Desktop/clinvar/ClinVar_all_vus_NMDescapeneriched.genes.txt')
  writeLines(all.genes,'~/Desktop/clinvar/ClinVar_all_vus_allgenes.txt')
  
  
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

 gene_column <- CLinvar_genes_NMDdata$gene
 
 write.table(gene_column, file = "clinvar_gene_list.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
 
