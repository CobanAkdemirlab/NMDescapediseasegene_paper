enrich2 = function(input_path){
  #import clinvar data
  ClinVar0820_all_pathogenic_NMDescapeneriched.genes <- read.table(input_path, quote="\"", comment.char="")
  search_type = function(input_string) {
    if (grepl("all", input_string)) {
      return("clinvar_all")
    } else if (grepl("css", input_string)) {
      return("clinvar_css")
    } else if (grepl("can", input_string)) {
      return("clinvar_can")
    } else if (grepl("long", input_string)) {
      return("clinvar_long")
    } else if (grepl("gen", input_string)) {
      return("genome")
    } else if (grepl("omim", input_string)) {
      return("omim")
    } else {
      return(NA)  # Return NA if none of the strings are found
    }
  }
  type = search_type(input_path)
  gene_list <- ClinVar0820_all_pathogenic_NMDescapeneriched.genes$V1
  enriched <- enrichr(gene_list, databases = c("GO_Biological_Process_2023"))
  #print(enriched)
  
  # Assuming 'enriched' contains the enrichment results for a specific database
  enriched_result <- enriched[["GO_Biological_Process_2023"]]  # or the respective result you want to visualize
  enriched_result = enriched_result[which(enriched_result$P.value<0.05),]
  p_values <- enriched_result$P.value
  terms <- enriched_result$Term
  #delete what's inside () in enriched_result$Term 
  term2 = gsub("\\(.*\\)", "", terms)
  enriched_result$Term = term2
  
  #wirte enriched_result to local
  #write.table(enriched_result, file=paste0('/Users/jxu14/Desktop/enrich/',type,'.txt'), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
  
  #pdf(paste0('/Users/jxu14/Desktop/enrich/clinvar_',type,'.pdf'))
  # Adjust margins to increase space for long terms on the left
  #par(mar = c(5, 15, 4, 2))  # Set bottom, left, top, and right margins respectively
  
  # Create the bar plot with modified margins
  #barplot(height = -log10(p_values), names.arg = terms, horiz = TRUE, las = 2,
  #        col = "skyblue", main = paste("Gene Enrichment Analysis of Clinvar",type), 
  #        xlab = "-log10(P-value)", ylab = "")
  hg_gs <- geneset::getGO(org = "human",ont = "bp")
  entrez_id <- mapIds(org.Hs.eg.db, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  ego <- genORA(entrez_id, geneset = hg_gs, p_cutoff = 0.1, q_cutoff = 0.1)
  ego <- ego[1:5, ]
  p = genekitr:: plotEnrich(ego,
                            plot_type = "dot",
                            scale_ratio = 0.4,
                            stats_metric = "pvalue",
                            term_metric = "RichFactor")
  
  
  return(p)
  #dev.off()
  
}
