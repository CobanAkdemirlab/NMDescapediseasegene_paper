
# --- Path resolution layer (auto-inserted) ------------------------------
# Resolve paths via paths.R for portability across machines:
#   data_file("x.csv")  locate by filename, errors clearly if not found
#   out_file("y.csv")   output to NMDESC_OUT (default ~/Desktop/NMDesc_out)
#   data_root("clinvar") use when a directory is needed instead of a file
.p <- c("gene level_v3/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../lib/paths.R", "../../../../lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p)) stop("paths.R not found -- run R from the repository root")
source(.p[1]); rm(.p)
# --------------------------------------------------------------------------

#draw figures for genes

##0. filter for automsomal dominat genes 
plus2_gnomAD_control_genes <- read_csv("~/Downloads/gene list/plus2_gnomAD_control_genes.csv", 
     col_names =  'hgnc_symbol')


# Rename hgnc_symbol to gene
colnames(plus2_gnomAD_control_genes)[colnames(plus2_gnomAD_control_genes) == "hgnc_symbol"] <- "gene"

# Ensure both columns are character type
plus2_gnomAD_control_genes$gene <- as.character(plus2_gnomAD_control_genes$gene)
merge_OMIM_mim2gene$gene <- as.character(merge_OMIM_mim2gene$gene)

# Merge the datasets
ad_genes <- merge(plus2_gnomAD_control_genes, merge_OMIM_mim2gene, by = "gene")

# Filter only Autosomal Dominant genes
ad_genes <- ad_genes %>% filter(inheritance == "AD")

# View result
write.csv(ad_genes, file = data_file("plus2_control_ADgenes.csv", must = FALSE), row.names = FALSE)

##1.input 6 list


##2. get attributes
att = get_inher(input_file)

##3. draw figures
#layout: 2*3

ggplot()