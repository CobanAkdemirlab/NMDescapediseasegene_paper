library(readxl)
nmd_esc_strict_list <- read_excel("~/Downloads/List of NMD escape disease genes/nmd_esc_strict_list.xlsx")
NMD_esc_extensive_list <- read_excel("~/Desktop/NMD_esc_extensive_list.xlsx")

plus1_css_gene0115 = read.table("~/Desktop/frameshift/plus1_css_gene0115.txt", header = F, sep = "\t")
plus1_can_gene0115 = read.table("~/Desktop/frameshift/plus1_can_gene0115.txt", header = F, sep = "\t")
plus1_long_gene0115 = read.table("~/Desktop/frameshift/plus1_long_gene0115.txt", header = F, sep = "\t")
plus2_css_gene0115 = read.table("~/Desktop/frameshift/plus2_css_gene0115.txt", header = F, sep = "\t")
plus2_can_gene0115 = read.table("~/Desktop/frameshift/plus2_can_gene0115.txt", header = F, sep = "\t")
plus2_long_gene0115 = read.table("~/Desktop/frameshift/plus2_long_gene0115.txt", header = F, sep = "\t")
plus1_trigger_gene0115 = read.table("~/Desktop/frameshift/plus1_trigger_gene0115.txt", header = F, sep = "\t")
plus2_trigger_gene0115 = read.table("~/Desktop/frameshift/plus2_trigger_gene0115.txt", header = F, sep = "\t")
snv_plp_ptc_p1120_NMDenriched2_all <- read_csv("~/Desktop/new_clinvar/snv_list/list4/snv_plp_ptc_p1120_NMDenriched2_all.txt", 
                                               col_names = FALSE)
snv_plp_ptc_p1120_NMDenriched2_can <- read_csv("~/Desktop/new_clinvar/snv_list/list4/snv_plp_ptc_p1120_NMDenriched2_can.txt", 
                                               col_names = FALSE)
snv_plp_ptc_p1120_NMDenriched2_css <- read_csv("~/Desktop/new_clinvar/snv_list/list4/snv_plp_ptc_p1120_NMDenriched2_css.txt", 
                                               col_names = FALSE)
snv_plp_ptc_p1120_NMDenriched2_long <- read_csv("~/Desktop/new_clinvar/snv_list/list4/snv_plp_ptc_p1120_NMDenriched2_long.txt", 
                                               col_names = FALSE)
snv_plp_ptc_p1120_NMDenriched2_trig <- read_csv("~/Desktop/new_clinvar/snv_list/list4/snv_plp_ptc_p1120_NMDenriched2_trig.txt", 
                                               col_names = FALSE)
fs_p_plus1_NMDenriched2_can = read.table("~/Desktop/new_clinvar/fs_p_plus1_NMDenriched2_can.txt", header = F, sep = "\t")
fs_p_plus1_NMDenriched2_css = read.table("~/Desktop/new_clinvar/fs_p_plus1_NMDenriched2_css.txt", header = F, sep = "\t")
fs_p_plus1_NMDenriched2_long = read.table("~/Desktop/new_clinvar/fs_p_plus1_NMDenriched2_long.txt", header = F, sep = "\t")
fs_p_plus2_NMDenriched2_can = read.table("~/Desktop/new_clinvar/fs_p_plus2_NMDenriched2_can.txt", header = F, sep = "\t")
fs_p_plus2_NMDenriched2_css = read.table("~/Desktop/new_clinvar/fs_p_plus2_NMDenriched2_css.txt", header = F, sep = "\t")
fs_p_plus2_NMDenriched2_long = read.table("~/Desktop/new_clinvar/fs_p_plus2_NMDenriched2_long.txt", header = F, sep = "\t")

length(which(plus1_can_gene0115$V1 %in% fs_p_plus1_NMDenriched2_can$V1))
length(which(plus1_css_gene0115$V1 %in% fs_p_plus1_NMDenriched2_css$V1))
length(which(plus1_long_gene0115$V1 %in% fs_p_plus1_NMDenriched2_long$V1))
length(which(plus2_can_gene0115$V1 %in% fs_p_plus2_NMDenriched2_can$V1))
length(which(plus2_css_gene0115$V1 %in% fs_p_plus2_NMDenriched2_css$V1))
length(which(plus2_long_gene0115$V1 %in% fs_p_plus2_NMDenriched2_long$V1))
load("/Users/jxu15/Downloads/hg38_seqfeatures.RData")
gene_name <- getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = sub("\\.\\d+$", "", merged.2$txnames),
  mart = ensembl
)
tr_id = getBM(
  attributes = c("hgnc_symbol","ensembl_transcript_id","transcript_is_canonical"),
  filters = "hgnc_symbol",
  values = snv_plp_ptc_p1120_NMDenriched2_can$X1,
  mart = ensembl
)
snv_can_can_id = getBM(
  attributes = c("ensembl_transcript_id","transcript_is_canonical"),
  filters = "ensembl_transcript_id",
  values = tr_id$ensembl_transcript_id,
  mart = ensembl
)
merged.2$gene = gene_name
merged.2$tr_sim = sub("\\.\\d+$", "", merged.2$txnames)
can.tr = getBM(
  attributes = c("transcript_is_canonical"),
  filters = "ensembl_transcript_id",
  values = merged.2$tr_sim,
  mart = ensembl
)
merged.2.can =merged.2[,]
plus1_can_df = merged.2[which(merged.2$tr_sim %in% plus1_can_list),]
plus1_css_df = merged.2[which(merged.2$tr_sim %in% plus1_css_list),]
plus1_long_df = merged.2[which(merged.2$tr_sim %in% plus1_long_list),]


sum(plus1_css_gene0115$V1 %in% NMD_esc_extensive_list$Gene)
sum(plus1_css_gene0108$V1 %in% NMD_esc_extensive_list$Gene)
sum(plus1_can_gene0115$V1 %in% NMD_esc_extensive_list$Gene)
sum(plus1_long_gene0115$V1 %in% NMD_esc_extensive_list$Gene)
sum(plus2_css_gene0115$V1 %in% NMD_esc_extensive_list$Gene)
sum(plus2_can_gene0115$V1 %in% NMD_esc_extensive_list$Gene)
sum(plus2_long_gene0115$V1 %in% NMD_esc_extensive_list$Gene)
sum(plus1_trigger_gene0115$V1 %in% NMD_esc_extensive_list$Gene)
sum(plus2_trigger_gene0115$V1 %in% NMD_esc_extensive_list$Gene)
sum(snv_plp_ptc_p1120_NMDenriched2_all$X1 %in% NMD_esc_extensive_list$Gene)
sum(snv_plp_ptc_p1120_NMDenriched2_can$X1 %in% NMD_esc_extensive_list$Gene)
sum(snv_plp_ptc_p1120_NMDenriched2_css$X1 %in% NMD_esc_extensive_list$Gene)
sum(snv_plp_ptc_p1120_NMDenriched2_long$X1 %in% NMD_esc_extensive_list$Gene)
sum(snv_plp_ptc_p1120_NMDenriched2_trig$X1 %in% NMD_esc_extensive_list$Gene)

getBM(
  attributes = c("ensembl_transcript_id",'cds_start','cds_end'), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values = 'ENST00000255476',                          # Transcript ID
  mart = ensembl                                         # Database connection
)

sum(plus1_css_gene0115$V1 %in% nmd_esc_strict_list$Gene)
sum(plus1_css_gene0108$V1 %in% nmd_esc_strict_list$Gene)
sum(plus1_can_gene0115$V1 %in% nmd_esc_strict_list$Gene)
sum(plus1_can_gene0108$V1 %in% nmd_esc_strict_list$Gene)
sum(plus1_long_gene0115$V1 %in% nmd_esc_strict_list$Gene)
sum(plus2_css_gene0115$V1 %in% nmd_esc_strict_list$Gene)
sum(plus2_can_gene0115$V1 %in% nmd_esc_strict_list$Gene)
sum(plus2_can_gene0115$V1 %in% nmd_esc_strict_list$Gene)
sum(plus2_long_gene0115$V1 %in% nmd_esc_strict_list$Gene)
sum(plus1_trigger_gene0115$V1 %in% nmd_esc_strict_list$Gene)
sum(plus2_trigger_gene0115$V1 %in% nmd_esc_strict_list$Gene)
sum(snv_plp_ptc_p1120_NMDenriched2_all$X1 %in% nmd_esc_strict_list$Gene)
sum(snv_plp_ptc_p1120_NMDenriched2_can$X1 %in% nmd_esc_strict_list$Gene)
sum(snv_plp_ptc_p1120_NMDenriched2_css$X1 %in% nmd_esc_strict_list$Gene)
sum(snv_plp_ptc_p1120_NMDenriched2_long$X1 %in% nmd_esc_strict_list$Gene)
sum(snv_plp_ptc_p1120_NMDenriched2_trig$X1 %in% nmd_esc_strict_list$Gene)

gene_set = getBM(
  attributes = c("ensembl_transcript_id","hgnc_symbol"), # Attributes to retrieve
  filters = "ensembl_transcript_id",                    # Filter to query
  values = transcript_set3,                          # Transcript ID
  mart = ensembl                                         # Database connection
)
gene_set2 = gene_set[which(gene_set$ensembl_transcript_id %in% transcript_set3_plp),]
sum(gene_set2$hgnc_symbol %in% nmd_esc_strict_list$Gene)
sum(gene_set2$hgnc_symbol %in% NMD_esc_extensive_list$Gene)