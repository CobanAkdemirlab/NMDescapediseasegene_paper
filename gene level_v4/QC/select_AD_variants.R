# Load required libraries
library(biomaRt)
library(dplyr)
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


# Get transcript-to-gene mapping from BioMart
transcript_gene_map <- getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  mart = ensembl
)

# Rename columns for merging
colnames(transcript_gene_map) <- c("transcript", "gene")



snv_gnomAD_control_variant <- read_csv(data_file("snv_gnomAD_control_variant.csv", must = FALSE))
snv_control_variant <- snv_gnomAD_control_variant %>%
  left_join(transcript_gene_map, by = "transcript")

snv_control_ADgenes <- read_csv("~/Desktop/autism/data/snv_control_ADgenes.csv")
# Keep only variants where gene is in snv_can_ADgenes$gene
filtered_variants <- snv_control_variant[which(snv_control_variant$gene %in% snv_control_ADgenes$gene),]
nrow(filtered_variants)
#write result
write.csv(filtered_variants,file = data_file("snv_control_ADvariants.csv", must = FALSE), row.names = FALSE)



tr_id = data.frame(transcript = res2@ranges@NAMES,
id = res2@elementMetadata@listData[["key"]])
tr_id_snv = tr_id[which(tr_id$id %in% snv_variants$x),]

snv_variant <- tr_id_snv  %>%
  distinct(id, .keep_all = TRUE)

snv_variant2 <- snv_variant %>%
  left_join(transcript_gene_map, by = "transcript")

filtered_variants <- snv_variant2[which(snv_variant2$gene %in% snv_can_ADgenes$gene),]
nrow(filtered_variants)
#write result
write.csv(filtered_variants,file = data_file("snv_ADvariants.csv", must = FALSE), row.names = FALSE)




