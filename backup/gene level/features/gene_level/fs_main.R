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

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#object_list = c('all','css','long','can','control','trig')
object_list = c('all','css','long','can')
getpath2 = function(input_name){
  #output_name = paste0("~/Desktop/new_clinvar/snv_list/list4/snv_plp_ptc_p1120_NMDenriched2_",input_name,".txt")
  output_name = paste0("~/Desktop/new_clinvar/snv_list/list4/snv_plp_nmd_p_NMDenriched_",input_name,".txt")
  return(output_name)
}
path_list2 = sapply(object_list, getpath2)

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
i=7
for(path in path_list2){
  g[[i]] = omim_check("~/Desktop/frameshift/plus1_trigger_gene0108.txt")
  i = i+1
}
grid.arrange(g[[7]],g[[8]],g[[9]],g[[10]],g[[11]],g[[12]],
             ncol = 6, nrow = 1)
#3. GO
source("~/Desktop/enrich/enrich2.R")
g[[13]] = enrich2(path_list2[1])
g[[14]] = enrich2(path_list2[2])
g[[15]] = enrich2(path_list2[3])
g[[16]] = enrich2(path_list2[4])
g[[17]] = enrich2(path_list2[5])
g[[18]] = enrich2(path_list2[6])
grid.arrange(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]],g[[6]],
             g[[7]],g[[8]],g[[9]],g[[10]],g[[11]],g[[12]],
             ncol = 6, nrow = 2)
g[[13]],g[[14]],g[[15]],g[[16]],g[[17]],g[[18]],
ncol = 6, nrow = 3)

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
