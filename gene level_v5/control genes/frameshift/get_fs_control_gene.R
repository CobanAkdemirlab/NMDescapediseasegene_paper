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

library(stringr)
library(readr)
library(biomaRt)
library(dplyr)

# `ensembl`, `mart`, `v_ch` and `omim_AD_symbols_tx` come from step 1 and from
# get_snv_control_gene.R when those ran in this session; each is rebuilt here
# otherwise, so the two control-gene scripts can run in either order.
if (!exists("ensembl"))
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
if (!exists("mart")) mart <- ensembl
if (!exists("v_ch")) v_ch <- read.csv(data_file("v_ch20260201.csv"))
if (!exists("omim_AD_symbols"))
  omim_AD_symbols <- read.csv(data_file("omim_AD_symbols.csv"),
                              header = TRUE)$hgnc_symbol
if (!exists("omim_AD_symbols_tx")) {
  omim_AD_symbols_tx <- getBM(
    attributes = c("ensembl_transcript_id", "hgnc_symbol", "transcript_is_canonical"),
    filters    = "hgnc_symbol",
    values     = omim_AD_symbols,
    mart       = mart)
  omim_AD_symbols_tx <- omim_AD_symbols_tx[which(omim_AD_symbols_tx$transcript_is_canonical == 1), ]
}

##exclude any gene with any plp frameshift NMD_can_esc mutation
fs_res = read_rds(data_file("fs.rds"))  # already restricted to plp and ptc
fs_res_can = fs_res[which(fs_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |fs_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)]

fs_gene2remove.ind = which(fs_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |fs_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)
fs_gene2remove = unique(fs_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][fs_gene2remove.ind])
# fs_gene2remove is an atomic vector of transcript ids, matched directly
# against the pool below; no symbol lookup is needed.
#remove genes2remove from gene_list
fs_gene_list = unique(fs_res@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]])
# setdiff() rather than fs_gene_list[-which(...)]: with nothing to remove,
# which() gives integer(0) and x[-integer(0)] returns an empty vector.
fs_genes_remain_in_clinvar = setdiff(fs_gene_list, fs_gene2remove)
length(unique(fs_genes_remain_in_clinvar)) 

#keep canonical transcripts only using getBM
fs_can.info = getBM(
  attributes = c("ensembl_transcript_id","transcript_is_canonical"), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values =  fs_genes_remain_in_clinvar,                          # Transcript ID
  mart = ensembl                                         # Database connection
)

fs_canonical_transcripts <- fs_can.info$ensembl_transcript_id[
  fs_can.info$transcript_is_canonical == 1
]

fs_genes_remain_in_clinvar <- fs_genes_remain_in_clinvar[
  fs_genes_remain_in_clinvar %in% fs_canonical_transcripts
]

length(unique(fs_genes_remain_in_clinvar))


#remove single exon genes using getBM
fs_BM.infoo <- getBM(
  attributes = c("ensembl_transcript_id","rank",'cds_start','cds_end','exon_chrom_start','exon_chrom_end'), # Attributes to retrieve
  filters = "ensembl_transcript_id",                  
  values =  fs_genes_remain_in_clinvar,                
  mart = ensembl                                         
)

fs_exon_counts <- aggregate(rank ~ ensembl_transcript_id, data = fs_BM.infoo, FUN = max)
colnames(fs_exon_counts) <- c("ensembl_transcript_id", "max_exon_rank")

# Keep only multi-exon transcripts (max rank > 1)
fs_multi_exon_transcripts <- fs_exon_counts$ensembl_transcript_id[
  fs_exon_counts$max_exon_rank > 1
]

fs_genes_remain_in_clinvar <- fs_genes_remain_in_clinvar[
  fs_genes_remain_in_clinvar %in% fs_multi_exon_transcripts
]

length(unique(fs_genes_remain_in_clinvar))


#get their hgnc_symbol
fs_genes_remain_in_clinvar_hgnc = getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = fs_genes_remain_in_clinvar,
  mart = mart
)
length(unique(fs_genes_remain_in_clinvar_hgnc$hgnc_symbol))

write.csv(unique(fs_genes_remain_in_clinvar_hgnc$hgnc_symbol), out_file("fs_control_gene.csv"))

fs_genes_remain_in_clinvar2 = fs_genes_remain_in_clinvar[which(fs_genes_remain_in_clinvar %in% v_ch$tx_id2)]

#filter for AD genes
fs_control_gene = unique(fs_genes_remain_in_clinvar2)
fs_control_gene_AD = fs_control_gene[which(fs_control_gene %in% omim_AD_symbols_tx$ensembl_transcript_id)]
write.csv(unique(fs_control_gene_AD), out_file("fs_control_genes_AD.csv"), row.names = FALSE)