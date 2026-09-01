##get fs variants from clinvar
library(stringr)
library(readr)
# --- Path resolution layer (auto-inserted) ------------------------------------------------
#   data_file("x.csv")  locate by filename, error clearly if not found
#   out_file("y.csv")   output to NMDESC_OUT (default ~/Desktop/NMDesc_out)
#   data_root("clinvar") use when a directory is needed instead of a file
.p <- c("gene level_v3/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../lib/paths.R", "../../../../lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p)) stop("paths.R not found -- run R from the repository root")
source(.p[1]); rm(.p)
# --------------------------------------------------------------------------

#input gene list
#fs_gene = read.csv('fs_can_syn_AD_gene_filtered_greater_20260201.txt', header = TRUE, stringsAsFactors = FALSE)
#fs_gene = read_csv(data_file("fs_can_AD_FDR0.05_gene.csv", must = FALSE))
fs_gene = read.csv('fs_can_AD_acat_FDR0.20_all.txt')$hgnc_symbol
#get canonical transcript of fs_gene$hgnc
fs_tr = getBM(attributes=c('ensembl_transcript_id','hgnc_symbol'),
               filters = 'hgnc_symbol', values = fs_gene, mart = ensembl)
fs_can_tr = getBM(attributes=c('ensembl_transcript_id','transcript_is_canonical'),
                   filters = 'ensembl_transcript_id', values = fs_tr$ensembl_transcript_id, mart = ensembl)

fs_can_tr = fs_can_tr[which(fs_can_tr$transcript_is_canonical == 1),1]

#in clinvar fs ptc, nmdesc variants, select those corresponding to the canonical transcript of gene list
fs_res = read_rds('fs.rds') #fs = res_plp_ptc[union(del_ind,ins_ind)]
fs_vus_res = read_rds('fs_vus.rds')
fs_benign_res = read_rds('fs_benign.rds')
#select NMDesc variants by canonical rule
fs_res_can = fs_res[which(fs_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |fs_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)]
fs_vus_res_can = fs_vus_res[which(fs_vus_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |fs_vus_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)]
fs_benign_res_can = fs_benign_res[which(fs_benign_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T |fs_benign_res@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)]
fs_variants = data.frame(transcript = fs_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][which(fs_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% fs_can_tr)],
                          key = fs_res_can@elementMetadata@listData[["key"]][which(fs_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% fs_can_tr)])
fs_vus_variants = data.frame(transcript = fs_vus_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][which(fs_vus_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% fs_can_tr)],
                          key = fs_vus_res_can@elementMetadata@listData[["key"]][which(fs_vus_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% fs_can_tr)])
fs_benign_variants = data.frame(transcript = fs_benign_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]][which(fs_benign_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% fs_can_tr)],
                          key = fs_benign_res_can@elementMetadata@listData[["key"]][which(fs_benign_res_can@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% fs_can_tr)])
length(unique(fs_variants$transcript))
#get uniprot id using biomart
fs_unique_transcripts <- unique(fs_variants$transcript)
fs_benign_unique_transcripts <- unique(fs_benign_variants$transcript)
fs_uniprot_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = fs_unique_transcripts,
  mart = ensembl
)
fs_benign_uniprot_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = fs_benign_unique_transcripts,
  mart = ensembl
)
fs_variants<- merge(
  fs_variants,
  fs_uniprot_mapping,
  by.x = "transcript",
  by.y = "ensembl_transcript_id",
  all.x = TRUE
)
fs_benign_variants<- merge(
  fs_benign_variants,
  fs_benign_uniprot_mapping,
  by.x = "transcript",
  by.y = "ensembl_transcript_id",
  all.x = TRUE
)


write.csv(fs_variants, 'fs_variants20260201_plp_acat_clinvar.csv', row.names = FALSE)
write.csv(fs_benign_variants, 'fs_variants20260201_benign_acat_clinvar.csv', row.names = FALSE)

remove_inframe <- function(df) {
  ref <- sub(".*\\|(.+)\\|.*", "\\1", df$key)
  alt <- sub(".*\\|.*\\|(.*)", "\\1", df$key)
  len_diff <- abs(nchar(ref) - nchar(alt))
  df[len_diff %% 3 != 0, ]   # keep length diff not a multiple of 3 (true frameshift)
}

# filter all three data frames
fs_variants2        <- remove_inframe(fs_variants)
write.csv(fs_variants2, 'fs_variants20260201_plp_clinvar_acat_0.2_filtered.csv', row.names = FALSE)
fs_vus_variants2    <- remove_inframe(fs_vus_variants)
write.csv(fs_vus_variants2, 'fs_variants20260201_vus_clinvar_acat_filtered.csv', row.names = FALSE)
fs_benign_variants <- remove_inframe(fs_benign_variants)
