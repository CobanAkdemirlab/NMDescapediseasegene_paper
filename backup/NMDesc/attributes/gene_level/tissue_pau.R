#### tissue specificity
gene.df = read.table('~/Downloads/gene.matrix.csv',sep=',',header=TRUE,stringsAsFactors=FALSE) # load table
gene.matrix = as.matrix(gene.df[,2:17]) # cast the FPKM columns to a matrix
rownames(gene.matrix) = gene.df$X # reassign the gene symbols as the row names for the matrix
barplot(gene.matrix['PRNP',],cex.names=.6,ylab='FPKMs',main='PRNP',col='yellow',las=3) # cex.names shrinks tissue labels, las=3 makes them perpendicular

#### for tau measure

tau<-function(x){
  if(any(is.na(x))) stop('NA\'s need to be 0.')
  if(any(x<0)) stop('Negative input values not permitted. Maybe data is log transformed?')
  d<-0
  for(i in 1:length(x)) d<-d+(1-x[i]/max(x))
  t<-d/(length(x)-1)
  t<-sum(1-x/max(x))/(length(x)-1)
}


tau <- apply(gene.matrix, 1, tau)

gene_ts_common <- tau[which(names(tau)%in%common_top5genes)]
gene_ts_exac <- tau[which(names(tau)%in%exac_top5genes)]
gene_ts_aric <- tau[which(names(tau)%in%aric_top5genes)]
genes_tolerant_ts <- tau[names(tau)%in%lof_tolerant]

allgenes_ts <- tau[which(names(tau)%in% syn.var$gene)]
genes_intolerant_ts <- tau[names(tau)%in%genes_intolerant]

boxplot(gene_ts_common,gene_ts_exac,gene_ts_aric,genes_intolerant_ts, genes_tolerant_ts, allgenes_ts)


density= data.frame(density=c(gene_ts_common,gene_ts_exac,gene_ts_aric,genes_tolerant_ts,genes_intolerant_ts,
                              allgenes_ts),category=c(rep('Depleted in ARIC and ExAC',length(gene_ts_common )),rep('Depleted only in ExAC',length(gene_ts_exac)),rep('Depleted only in ARIC',length(gene_ts_aric)),rep('LOF-tolerant genes',length(genes_tolerant_ts)), rep('Intolerant genes',length(genes_intolerant_ts)),rep('All genes',length(allgenes_ts)) ))

ggplot(density) + aes(density, fill=category)+ geom_density(alpha = 0.2)+ xlim(0,1)


p <- ggplot(density) + aes(x=category, y=density) +
  geom_violin(trim=FALSE) +  geom_boxplot(width=0.1)
p +theme(legend.position='none',axis.text.x = element_text(size = 10,angle=90, face='bold'),axis.title.y = element_text(size = 10,angle=90, face='bold')) +
  ylab('Tissue specificity')
