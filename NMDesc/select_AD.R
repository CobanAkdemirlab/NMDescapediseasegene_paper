#draw figures for genes

##0. filter for automsomal dominat genes 
plus2_gnomAD_control_genes <- read_csv("~/Downloads/gene list/plus2_gnomAD_control_genes.csv", 
     col_names =  'hgnc_symbol')


# Rename hgnc_symbol to gene in plus2_gnomAD_control_genes
colnames(plus2_gnomAD_control_genes)[colnames(plus2_gnomAD_control_genes) == "hgnc_symbol"] <- "gene"

# Ensure both columns are character type
plus2_gnomAD_control_genes$gene <- as.character(plus2_gnomAD_control_genes$gene)
merge_OMIM_mim2gene$gene <- as.character(merge_OMIM_mim2gene$gene)

# Merge the datasets
ad_genes <- merge(plus2_gnomAD_control_genes, merge_OMIM_mim2gene, by = "gene")

# Filter only Autosomal Dominant genes
ad_genes <- ad_genes %>% filter(inheritance == "AD")

# View result
write.csv(ad_genes, file = "/Users/jxu14/Desktop/autism/data/plus2_control_ADgenes.csv", row.names = FALSE)

##1.input 6 list


##2. get attributes
att = get_inher(input_file)

##3. draw figures
#layout: 2*3

ggplot()