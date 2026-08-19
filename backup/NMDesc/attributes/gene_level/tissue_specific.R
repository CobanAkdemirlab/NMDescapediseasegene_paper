#tissue specific analysis
library(readr)
library(ggplot2)
gtex_data <- read_tsv("/Users/jxu14/Desktop/big_picture/pfam/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", skip = 2)
#mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
calculate_tau_list <- function(gene_list){
  tau_data = data.frame(Gene = character(), Tau = numeric())
  gene_id_list <- getBM(
    attributes = c("external_gene_name", "ensembl_gene_id"),
    filters = "external_gene_name",
    values = gene_list,
    mart = mart
  )
  for(gene_id in gene_id_list$ensembl_gene_id){
    #1. Get expression data for each tissue of this gene
    expression_data = gtex_data[which(substr(gtex_data$Name,1,15)==gene_id),]
    expression_value = as.numeric(expression_data[1,-c(1,2)]) #row 1 and row 2 are gene id and gene name
    #2. Calculate tau value for this gene by expression data across different tissues
    expression_value = log2(expression_value)
    max_exp <- max(expression_value)
    tau <- sum(1 - (expression_value / max_exp)) / (length(expression_value) - 1)
    new_tau_data <- data.frame(Gene = gene_id, Tau = tau)
    #remove where expression value < 0
    new_tau_data = new_tau_data[expression_value >= 0,]
    #add new_tau_data to tau_data
    tau_data = rbind(tau_data,new_tau_data)
  }
  #remove the rows where tau is NA
  tau_data <- tau_data[!is.na(tau_data$Tau),]
  return(tau_data)
}

object_list = c('css','long','can')
getpath_snv= function(input_name){
  output_name = paste0("~/Desktop/new_clinvar/snv_list/list4/snv_plp_ptc_p1120_NMDenriched2_",input_name,".txt")
  return(output_name)
}
getpath_plus1 = function(input_name){
  output_name = paste0("~/Desktop/frameshift/plus1_",input_name,"_gene0115.txt")
  return(output_name)
}
getpath_plus2 = function(input_name){
  output_name = paste0("~/Desktop/frameshift/plus2_",input_name,"_gene0115.txt")
  return(output_name)
}
path_list_snv = sapply(object_list, getpath_snv)
path_list_plus1 = sapply(object_list, getpath_plus1)
path_list_plus2 = sapply(object_list, getpath_plus2)
tau_data = data.frame(Gene = character(), Tau = numeric(), object = character(),group = character())
for(i in 1:length(object_list)){
  object = object_list[i]
  path = path_list_snv[i]
  gene_list = read.table(path,quote="\"", comment.char="",col.names = 'gene')
  tau_df = calculate_tau_list(gene_list$gene)
  #add a col, where the value is the object_list[i]
  tau_df = cbind(tau_df, object = paste0(object,"_snv"),group = object)
  tau_data = rbind(tau_data, tau_df)
  
  path = path_list_plus1[i]
  gene_list = read.table(path,quote="\"", comment.char="",col.names = 'gene')
  tau_df = calculate_tau_list(gene_list$gene)
  #add a col, where the value is the object_list[i]
  tau_df = cbind(tau_df, object=paste0(object,"_plus1"),group = object)
  tau_data = rbind(tau_data, tau_df)
  
  path = path_list_plus2[i]
  gene_list = read.table(path,quote="\"", comment.char="",col.names = 'gene')
  tau_df = calculate_tau_list(gene_list$gene)
  #add a col, where the value is the object_list[i]
  tau_df = cbind(tau_df, object = paste0(object,"_plus2"),group = object)
  tau_data = rbind(tau_data, tau_df)
}                  

control_tau = calculate_tau_list(control$gene[-1])
tau_df = cbind(control_tau, object = "control",group = "control")
tau_data = rbind(tau_data, tau_df)


#3. Plot tau value distribution by violin plot
ggplot(tau_data, aes(x = object, y = Tau,fill = group)) +
  #set y rang to 0,1
  ylim(0,1.25) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9) +
  labs(title = "Tau Value Distribution (expression log2 transform and fliter by >=0)", y = "Tau Value", x = "") +
  theme_classic() +
  scale_fill_brewer(palette = "Set2", name = "") +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")





