#txdb = makeTxDbFromGFF('~/Desktop/NMDrevisioncode/gencode.v26.primary_assembly.annotation.gtf.gz')
library(AnnotationDbi)
#saveDb(txdb, '~/Desktop/NMDrevisioncode/txdb.gencode26.sqlite')
ensgene <- txdb
chr.list <- paste('chr',1:22,sep='')
ensgene.sub <- keepSeqlevels(ensgene,chr.list)
ensgene <- ensgene.sub
cds_seqs <- extractTranscriptSeqs(Hsapiens, cdsBy(ensgene, by="tx",use.names=TRUE))
utr.grange <-threeUTRsByTranscript(ensgene, use.names=TRUE)
fiveutr.grange<-fiveUTRsByTranscript(ensgene, use.names=TRUE)
introns.grange<- intronsByTranscript(ensgene, use.names=TRUE)
exons.grange<-exonsBy(ensgene, use.names=TRUE)
cds.grange<- cdsBy(ensgene, by="tx",use.names=TRUE)

### get the nuclear localization signal in the first 200 bp
cds.resized=resize(cds.grange,width=200,fix='start')
nucloc <- read.csv('~/Desktop/unipLocSignal.bed.csv',header=FALSE,sep='\t',stringsAsFactors=FALSE)
#nucloc <- nucloc[which(nucloc$V22=='short sequence motif'),]
gr2 <- GRanges(seqnames = nucloc$V1, ranges = IRanges(start = nucloc$V2, end=nucloc$V3))
hits <- unique(queryHits(findOverlaps(cds.resized,gr2)))
cds.nucloc <- names(cds.grange)[hits]
merged.2$cds.nucloc <- rep('NA',nrow(merged.2))
merged.2[which(merged.2$txnames%in%cds.nucloc),'cds.nucloc'] <- 'There is a nuclear localization signal in the first 200 bp'
save(merged.2, file='~/Desktop/NMDrevisioncode/hg38_seqfeatures.RData')
variants.features.fr$cds.nucloc <- rep('NA',nrow(variants.features.fr))
variants.features.fr[which(variants.features.fr$txnames%in%cds.nucloc),'cds.nucloc'] <- 'There is a nuclear localization signal in the first 200 bp'
