#combine uni with loc
library(data.table)
library(dplyr)
library(stringr)
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

plus1_control_variants0406 <- read_csv(data_file("plus1_control_gnomAD_variants0406.csv", must = FALSE))
plus1_control_data <- fread("~/Downloads/pro_loc/plus1_control_data.txt")
plus1_control_variants0406 = plus1_control_variants0406[-3621,]
#in plus1_control_data, remove rows where protein position is -
plus1_control_data <- plus1_control_data %>%
  filter(Protein_position != "-")


plus1_control_merged_uni_loc = data.frame(
  variant_key = plus1_control_variants0406$id,
  up_va = plus1_control_data$Uploaded_variation,
  transcript_id = plus1_control_variants0406$transcript,
  location = plus1_control_data$Protein_position
)
write.csv(plus1_control_merged_uni_loc, "plus1_control_merged_uni_loc.csv", row.names = FALSE)

