#omim check
omim_check <- function(path){
  search1 = function(path) {
    if (grepl("all", path)) {
      return("clinvar_all")
    } else if (grepl("css", path)) {
      return("clinvar_css")
    } else if (grepl("can", path)) {
      return("clinvar_can")
    } else if (grepl("long", path)) {
      return("clinvar_long")
    } else if (grepl("gen", path)) {
      return("genome")
    } else if (grepl("control", path)) {
      return("control")
    } else {
      return('trig')  # Return NA if none of the strings are found
    }
  }
  obj = search1(path)
  genes.mat <- NULL
  #read in gene list according to path
  input_gene <- read.table(path,quote="\"", comment.char="",col.names = 'gene')
  genes.mat$gene <- input_gene$gene
  genes.mat.1 <- merge(genes.mat,Lof_metrics,by='gene')
  colnames(mim2gene)[which(names(mim2gene)=='X..MIM.Number')] <- 'MIM.Number'
  merge_OMIM_mim2gene <- merge(mim2gene,OMIM,by='MIM.Number')
  colnames(merge_OMIM_mim2gene)[which(names(merge_OMIM_mim2gene)=='Approved.Gene.Symbol..HGNC.')] <- 'Raw_gene_name'
  genes.mat.1$Raw_gene_name <- genes.mat.1$gene

  data.new.merged <- merge(genes.mat.1,merge_OMIM_mim2gene,by='Raw_gene_name',all.x=TRUE)
  colnames(data.new.merged)[which(names(data.new.merged)=='Phenotypes')] <- 'OMIM_disorder'

  final_omim_hg38_data$Inheritance_pattern <- ""
  final_omim_hg38_data$Inheritance_pattern[grep("Autosomal dominant", final_omim_hg38_data$Phenotypes)] <- "Autosomal dominant"
  final_omim_hg38_data$Inheritance_pattern[grep("Autosomal recessive", final_omim_hg38_data$Phenotypes)] <- paste0(final_omim_hg38_data$Inheritance_pattern[grep("Autosomal recessive", final_omim_hg38_data$Phenotypes)], ", ", "Autosomal recessive")
  final_omim_hg38_data$Inheritance_pattern[grep("Digenic recessive", final_omim_hg38_data$Phenotypes)] <- "Digenic recessive"
  final_omim_hg38_data$Inheritance_pattern[grep("X-linked recessive", final_omim_hg38_data$Phenotypes)] <- "X-linked recessive"
  final_omim_hg38_data$Inheritance_pattern[grep("X-linked dominant", final_omim_hg38_data$Phenotypes)] <- "X-linked dominant"
  final_omim_hg38_data$Inheritance_pattern[grep("Y-linked", final_omim_hg38_data$Phenotypes)] <- "Y-linked"
  
  domi.ind <- grep('Autosomal dominant',data.new.merged$OMIM_disorder,ignore.case = T)
  rec.ind <- grep('Autosomal recessive',data.new.merged$OMIM_disorder,ignore.case = T)
  xr.ind <- grep('X-linked recessive',data.new.merged$OMIM_disorder,ignore.case = T)
  xd.ind <- grep('X-linked dominant',data.new.merged$OMIM_disorder,ignore.case = T)
  yl.ind <- grep('Y-linked',data.new.merged$OMIM_disorder,ignore.case = T)
  dr.ind <- grep('Digenic recessive',data.new.merged$OMIM_disorder,ignore.case = T)
  
  data.new.merged$inheritance <- rep('Unknown',nrow(data.new.merged))
  data.new.merged$inheritance[setdiff(domi.ind,intersect(domi.ind,rec.ind))] <- 'AD'
  data.new.merged$inheritance[setdiff(rec.ind,intersect(rec.ind,domi.ind))] <- 'AR'
  data.new.merged$inheritance[intersect(domi.ind,rec.ind)] <- 'AD/AR'
  data.new.merged$inheritance[xr.ind] <- 'XR'
  data.new.merged$inheritance[xd.ind] <- 'XD'
  data.new.merged$inheritance[yl.ind] <- 'Y-linked'
  data.new.merged$inheritance[dr.ind] <- 'Digenic recessive'

  col_set = ifelse(obj == 'clinvar_all','Blues',
                 ifelse(obj == 'clinvar_css','Reds',
                        ifelse(obj == 'clinvar_long','Greens',
                               ifelse(obj == 'clinvar_can','Purples',
                                      ifelse(obj == 'control', 'Oranges',
                                             'Greys'))))) #change it, organe are similar to red

  #plot it by pie chart
  data.new.merged$inheritance <- factor(data.new.merged$inheritance,levels=c('AD','AR','AD/AR','Unknown','XR','XD','Y-linked','Digenic recessive'))
  data.new.merged = data.new.merged[!is.na(data.new.merged$inheritance),]
  
  all_invs_category_2_omim_inh_pattern <- data.new.merged %>%
    group_by(inheritance) %>%
    summarise(n = n()) %>%
    mutate(value = round((n / sum(n)) * 100, 1))
  df_cat2 <- all_invs_category_2_omim_inh_pattern %>% 
    mutate(csum = rev(cumsum(rev(value))), 
           pos = value/2 + lead(csum, 1),
           pos = if_else(is.na(pos), value/2, pos))
  max_value <- max(df_cat2$value)
  g=ggplot(all_invs_category_2_omim_inh_pattern, aes(x = "" , y = value, fill = inheritance)) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = col_set, direction = -1) +
    geom_label_repel(data = df_cat2,
                     aes(y = pos, label = paste0(value, "%")),
                     size = 3, nudge_x = 1, show.legend = F,color = "goldenrod1") +
    guides(fill = guide_legend(title = "Inheritance")) +
    theme_void()
  g2 = ggplot(data.new.merged,aes(x='',fill=inheritance)) + geom_bar(width = 1) + 
   coord_polar("y") + theme_minimal() + 
   #change the color to gradient, like blue, light blue etc
   scale_fill_brewer(palette = col_set, direction = -1) +
   #don't show the legend
   #theme(legend.position = "none") +
   #show the percentage 
   geom_text(aes(label = scales::percent(..count.. / sum(..count..))), stat = "count", position = position_stack(vjust = 0.6),size=1) +
   labs(title = obj,x='') +
   theme(
      axis.text = element_blank(), 
      axis.ticks = element_blank(), 
      axis.title = element_blank(), 
      axis.line = element_blank(), 
      panel.grid = element_blank(),   # Remove grid lines
      panel.border = element_blank()
    )



  return(g)
}
