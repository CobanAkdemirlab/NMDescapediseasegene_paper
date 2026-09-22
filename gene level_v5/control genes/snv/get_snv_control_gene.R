# paths.R supplies data_file() and out_file(), which resolve every filename
# below against the data roots rather than the working directory.
.p <- c("lib/paths.R", "gene level_v5/lib/paths.R", "../lib/paths.R",
        "../../lib/paths.R", "../../../lib/paths.R",
        "../../../gene level_v5/lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p) && !exists("data_file"))
  stop("paths.R not found -- run R from the repository root", call. = FALSE)
if (length(.p)) source(.p[1])
rm(.p)

library(biomaRt)
library(dplyr)
library(readr)

# Two distinct sets, both built by gene_get_main.R (step 1):
#   res_snv_plp_ptc  the control pool -- snv_ind & plp_ind & ptc_ind, every
#                    transcript carrying a P/LP PTC SNV, with no canonical
#                    restriction; canonical filtering happens further down
#                    via getBM
#   snv_res_can      the transcripts to exclude -- the NMD-escaping, canonical
#                    subset of that pool
# The frameshift side has the same shape, with fs.rds as the pool and
# fs_res_can as its escaping subset. The two must stay distinct: one object
# serving as both would make the pool and the exclusion set identical and leave
# no controls behind.
# Each is read from the file step 1 writes when it is absent from the session.
if (!exists("res_snv_plp_ptc"))
  res_snv_plp_ptc <- readRDS(data_file("snv_plp_ptc20260201.rds"))
if (!exists("snv_res_can"))
  snv_res_can <- readRDS(data_file("snv_plp_ptc_nmdesc_can20260201.rds"))
if (!exists("mart"))
  mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
if (!exists("omim_AD_symbols"))
  omim_AD_symbols <- read.csv(data_file("omim_AD_symbols.csv"),
                              header = TRUE)$hgnc_symbol

##exclude any gene genes with any frameshift NMD_can_esc mutation
snv_gene2remove.ind = which(snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)
snv_gene2remove = unique(snv_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][snv_gene2remove.ind])
#get hgnc_symbol of gene2remove using getBM
snv_gene2remove = getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = snv_gene2remove,
  mart = mart
)
#remove genes2remove from snv genes with plp ptc variants in clinvar
snv_gene_list = unique(res_snv_plp_ptc@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]])
# setdiff() rather than snv_gene_list[-which(...)]: with nothing to remove,
# which() gives integer(0) and x[-integer(0)] returns an empty vector, which
# would drop the entire pool.
snv_genes_remain_in_clinvar = setdiff(snv_gene_list, snv_gene2remove$ensembl_transcript_id)

#only keep the genes that passed the variant_summary filter
v_ch = read.csv(data_file("v_ch20260201.csv"))
snv_genes_remain_in_clinvar2 = snv_genes_remain_in_clinvar[which(snv_genes_remain_in_clinvar %in% v_ch$tx_id2)]

length(unique(snv_genes_remain_in_clinvar2)) 

#filter for AD genes
omim_AD_symbols_tx = getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol",'transcript_is_canonical'),
  filters = "hgnc_symbol",
  values = omim_AD_symbols,
  mart = mart
)
#keep canonical transcript only
omim_AD_symbols_tx = omim_AD_symbols_tx[which(omim_AD_symbols_tx$transcript_is_canonical == 1),]
snv_control_gene = unique(snv_genes_remain_in_clinvar2)
snv_control_gene_AD = snv_control_gene[which(snv_control_gene %in% omim_AD_symbols_tx$ensembl_transcript_id)]
write.csv(unique(snv_control_gene_AD), out_file("snv_control_genes_AD.csv"), row.names = FALSE)
write.csv(unique(snv_genes_remain_in_clinvar2), out_file("snv_control_genes.csv"), row.names = FALSE)