library(data.table)
library(stringr)
library(biomaRt)
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

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
for(i in 16:22){
  chr.name = paste0(data_file("gnomad.exomes.v4.1.sites.chr", must = FALSE),i,".cut.syn.vcf")
  chr.syn <- fread(
         cmd = paste("grep -v '^#'", shQuote(chr.name)),
         sep = "\t",
         header = FALSE,
         fill = TRUE,
         showProgress = TRUE
     )
  chr.syn <- chr.syn[, 1:8]
  setnames(chr.syn, c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"))
  syn = chr.syn
  rm(chr.syn)
  # key: chr:000963993|C|T  (POS 固定 9 位补0)
  syn[, key := sprintf("%s:%09d|%s|%s", CHROM, POS, REF, ALT)]
  syn[, is_syn := str_detect(INFO, "synonymous_variant")]
  #check for synonymous variants
  #table(syn$is_syn, useNA = "ifany")
  syn = syn[is_syn == TRUE]
  #get transcript id, only keep enst
  syn[, transcript_id := str_extract(INFO, "ENST[0-9]+")]
  #delete info
  syn = syn[, -'INFO']
  #filter for canonical transcript
  tx_set = unique(na.omit(syn$transcript_id))
  tx_can = getBM(attributes = c("ensembl_transcript_id", "transcript_is_canonical"),
                 filters = "ensembl_transcript_id",
                 values = tx_set,
                 mart = ensembl)
  tx_can_set = tx_can[which(tx_can$transcript_is_canonical == 1), 'ensembl_transcript_id']
  syn_can = syn[transcript_id %in% tx_can_set]
  rm(syn)
  write.csv(syn_can, paste0(data_file("gnomad.exomes.v4.1.sites.chr", must = FALSE),i,".cut.syn.can.csv"), row.names = FALSE)
  rm(syn_can)
  rm(tx_can)
  rm(tx_can_set)
}
