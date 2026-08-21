#big picture, do the four clinvar ones first
library(enrichR)
library(biomaRt)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(genekitr)

omim <- readLines("~/Desktop/autism/data/genemap2.txt")
omim_hg38_data <- read.delim("~/Desktop/autism/data/genemap2.txt", skip = 3, header = TRUE, nrows=length(omim) - 80)
final_omim_hg38_data <- dplyr::select(omim_hg38_data, X..Chromosome, Genomic.Position.Start, Genomic.Position.End, MIM.Number, Approved.Gene.Symbol, Gene.Locus.And.Other.Related.Symbols, Gene.Name, Entrez.Gene.ID, Ensembl.Gene.ID, Phenotypes, Comments)
#select pathogenic/lp
final_omim_hg38_data = final_omim_hg38_data[which(nchar(final_omim_hg38_data[,'Phenotypes'])>1),]
colnames(final_omim_hg38_data) <- c("chr", "start", "end", "Gene_stable_ID", "approved_gene_symbol", "other_gene_symbol", "gene_name", "entrez_gene_id", "Ensembl_gene_ID", "Phenotypes", "Comments")
final_omim_hg38_data$Inheritance_pattern <- ""
final_omim_hg38_data$Inheritance_pattern[grep("Autosomal dominant", final_omim_hg38_data$Phenotypes)] <- "Autosomal dominant"
final_omim_hg38_data$Inheritance_pattern[grep("Autosomal recessive", final_omim_hg38_data$Phenotypes)] <- paste0(final_omim_hg38_data$Inheritance_pattern[grep("Autosomal recessive", final_omim_hg38_data$Phenotypes)], ", ", "Autosomal recessive")
final_omim_hg38_data$Inheritance_pattern[grep("Digenic recessive", final_omim_hg38_data$Phenotypes)] <- "Digenic recessive"
final_omim_hg38_data$Inheritance_pattern[grep("X-linked recessive", final_omim_hg38_data$Phenotypes)] <- "X-linked recessive"
final_omim_hg38_data$Inheritance_pattern[grep("X-linked dominant", final_omim_hg38_data$Phenotypes)] <- "X-linked dominant"
final_omim_hg38_data$Inheritance_pattern[grep("Y-linked", final_omim_hg38_data$Phenotypes)] <- "Y-linked"
final_omim_hg38_data$Inheritance_pattern[which(final_omim_hg38_data$Inheritance_pattern  == ", Autosomal recessive")] = "Autosomal recessive"
final_omim_hg38_data$Inheritance_pattern[which(final_omim_hg38_data$Inheritance_pattern  == "")] = "Other"
table(final_omim_hg38_data$Inheritance_pattern)
omim2 = final_omim_hg38_data %>% select(Ensembl_gene_ID, Inheritance_pattern)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","phenotype_description"), filters = "ensembl_gene_id", values = omim2$Ensembl_gene_ID, mart = ensembl)
omim2 = merge(omim2, genes, by.x = "Ensembl_gene_ID", by.y = "ensembl_gene_id")
omim3 = unique(omim2$hgnc_symbol)
write.csv(omim3, "~/Desktop/new_clinvar/snv_plp_nmd_p_NMDenriched_omim.txt", row.names = FALSE,col.names = F, quote = F)

#object_list = c('all','css','long','can','control','trig')
#object_list = c('all','css','long','can')
object_list = c('css','long','can','trigger')
getpath1 = function(input_name){
  #output_name = paste0("~/Desktop/new_clinvar/snv_list/list4/snv_plp_ptc_p1120_NMDenriched2_",input_name,".txt")
  output_name = paste0("~/Desktop/frameshift/plus1_",input_name,"_gene0115.txt")
  return(output_name)
}
getpath2 = function(input_name){
  #output_name = paste0("~/Desktop/new_clinvar/snv_list/list4/snv_plp_ptc_p1120_NMDenriched2_",input_name,".txt")
  output_name = paste0("~/Desktop/frameshift/plus2_",input_name,"_gene0115.txt")
  return(output_name)
}
path_list1 = sapply(object_list, getpath1)
path_list2 = sapply(object_list, getpath2)
path_list = c(path_list1,path_list2)
 

#1. length of coding sequence
gnomad.v2.1.1.lof_metrics.by_gene <- read.delim("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/桌面 - q的MacBook Pro/autism/data/gnomad.v2.1.1.lof_metrics.by_gene.txt")
i=1
g <- list()
source("~/Desktop/cds_length/code/boxplot_clinvar.R")
for(path in path_list2){
  g[[i]] = boxplot_clinvar(path)
  i = i+1
}

#2. AR/AD
Lof_metrics <- read.delim("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/桌面 - q的MacBook Pro/autism/data/gnomad.v2.1.1.lof_metrics.by_gene.txt")
OMIM <- read.csv('/Users/jxu14/Desktop/autism/data/genemap2.txt',sep='\t',stringsAsFactors=FALSE,header=TRUE,skip=3)
mim2gene <- read.csv('/Users/jxu14/Desktop/autism/data//mim2gene.txt',sep='\t',stringsAsFactors=FALSE,header=TRUE,skip=4)
source("~/Desktop/omim/omim_check.R")
i=1
for(path in path_list){
   g[[i]] = omim_check(path)
   i = i+1
}
grid.arrange(g[[1]],g[[2]],g[[3]],g[[4]],
             g[[5]],g[[6]],g[[7]],g[[8]],
             ncol = 4, nrow = 2)
#3. GO
source("~/Desktop/enrich/enrich2.R")
g[[9]] = enrich2(path_list[1])
g[[10]] = enrich2(path_list[2])
g[[11]] = enrich2(path_list[3])
g[[12]] = enrich2(path_list[4])
g[[13]] = enrich2(path_list[5])
g[[14]] = enrich2(path_list[6])
g[[15]] = enrich2(path_list[7])
g[[16]] = enrich2(path_list[8])
grid.arrange(g[[9]],g[[10]],
             g[[13]],g[[14]],
             ncol = 2, nrow = 2)
grid.arrange(g[[11]],g[[12]],
             g[[15]],g[[16]],
             ncol = 2, nrow = 2)

for(i in 1:5){
  writeLines(clinvar_all$Genes[i],paste0('/Users/jxu14/Desktop/enrich/clinvar_plp_all_',clinvar_all$Term[i],'.txt'))
  writeLines(clinvar_can$Genes[i],paste0('/Users/jxu14/Desktop/enrich/clinvar_plp_can_',clinvar_can$Term[i],'.txt'))
  writeLines(clinvar_css$Genes[i],paste0('/Users/jxu14/Desktop/enrich/clinvar_plp_css_',clinvar_css$Term[i],'.txt'))
  writeLines(clinvar_long$Genes[i],paste0('/Users/jxu14/Desktop/enrich/clinvar_plp_long_',clinvar_long$Term[i],'.txt'))
  writeLines(omim$Genes[i],paste0('/Users/jxu14/Desktop/enrich/omim_plp_',omim$Term[i],'.txt'))
}



print(g[[7]])
print(g[[8]])
print(g[[9]])
print(g[[10]])
print(g[[11]])
print(g[[12]])
