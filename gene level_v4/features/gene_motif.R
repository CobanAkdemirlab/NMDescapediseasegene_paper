#This R script is to get motif count based on outside datasource
library(stringr)
library(dplyr)
library(readxl)
library(readr)
# --- 路径解析层（自动插入）------------------------------------------------
# 原来这里是旧机器 /Users/jxu14/ 的绝对路径，换机器必断。改走 paths.R：
#   data_file("x.csv")  按文件名定位，找不到会明确报错
#   out_file("y.csv")   输出到 NMDESC_OUT（默认 ~/Desktop/NMDesc_out）
#   data_root("clinvar") 需要目录而非文件时用
.p <- c("gene level_v3/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../lib/paths.R", "../../../../lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p)) stop("找不到 paths.R —— 请从仓库根目录运行 R")
source(.p[1]); rm(.p)
# --------------------------------------------------------------------------

touni <- read_csv(data_file("NIHMS1818854-supplement-2(A).csv"))
motif_doc<- read_csv(data_file("NIHMS1818854-supplement-2(B).csv"))
LCS_doc <- read_excel(data_file("Copy of NIHMS1818854-supplement-2.xls"), 
                      sheet = "G")
#select largest value in each cell
get_max_from_cell <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  nums <- str_extract_all(x, "\\d+")[[1]]
  if (length(nums) == 0) return(NA_real_)
  max(as.numeric(nums), na.rm = TRUE)
}

motif_max <- motif_doc %>%
  mutate(across(-1, ~ vapply(as.character(.x), get_max_from_cell, numeric(1))))
LCS_max <- LCS_doc %>%
  mutate(across(-1, ~ vapply(as.character(.x), get_max_from_cell, numeric(1))))
motif_max$uniprot = touni$`Uniprot ID`[match(motif_max$Protein, touni$Protein)]
LCS_max$uniprot = touni$`Uniprot ID`[match(LCS_max$Protein, touni$Protein)]

#combine motif value with variant loc and nmdesc loc, match by uniprot id
#this is protein level, which is gene level, so for the variant list, we compare with variant loc
#motif_A$uniprot_ID

uni_map <- getBM(
  attributes = c("hgnc_symbol", "uniprotswissprot"),
  filters    = "hgnc_symbol",
  values     = unique(gene_all$hgnc_symbol),
  mart       = ensembl
)

gene_all$uniprot <- uni_map$uniprotswissprot[
  match(gene_all$hgnc_symbol, uni_map$hgnc_symbol)
]
merge(motif_max,gene_all, by = "uniprot") -> motif_max3


motif_max3$gene_protein_flag = motif_max3$`Protein Features`*3 >= motif_max3$NMD_region_start
motif_max3$gene_domains_flag = motif_max3$`Domains`*3 >= motif_max3$NMD_region_start
motif_max3$gene_slim_flag = motif_max3$`SLiMs`*3 >= motif_max3$NMD_region_start
motif_max3$gene_morf_flag = motif_max3$MORFs*3 >= motif_max3$NMD_region_start
motif_max3$gene_ptm_flag = motif_max3$`PTMs`*3 >= motif_max3$NMD_region_start
motif_max3$gene_nls_flag = motif_max3$`NLSs`*3 >= motif_max3$NMD_region_start
LCS_max3 = merge(LCS_max, gene_all, by = "uniprot")
LCS_max3$gene_LCS_flag = LCS_max3$`LCSs`*3 >= LCS_max3$NMD_region_start

#check if all gene is omim AD
sum(gene_all$hgnc_symbol %in% omim_AD_symbols)

write.csv(LCS_max3, "gene_LCS_flags.csv", row.names = FALSE)
LCS_max3 = read.csv("gene_LCS_flags.csv")
write.csv(motif_max3, "gene_motif_flags.csv", row.names = FALSE)

#merge the flags back to gene_all by uniprot
gene_all$gene_protein_flag = motif_max3$gene_protein_flag[match(gene_all$uniprot, motif_max3$uniprot)]
gene_all$gene_domains_flag = motif_max3$gene_domains_flag[match(gene_all$uniprot, motif_max3$uniprot)]
gene_all$gene_slim_flag = motif_max3$gene_slim_flag[match(gene_all$uniprot, motif_max3$uniprot)]
gene_all$gene_morf_flag = motif_max3$gene_morf_flag[match(gene_all$uniprot, motif_max3$uniprot)]
gene_all$gene_ptm_flag = motif_max3$gene_ptm_flag[match(gene_all$uniprot, motif_max3$uniprot)]
gene_all$gene_nls_flag = motif_max3$gene_nls_flag[match(gene_all$uniprot, motif_max3$uniprot)]
gene_all$gene_LCS_flag = LCS_max3$gene_LCS_flag[match(gene_all$uniprot, LCS_max3$uniprot)]
