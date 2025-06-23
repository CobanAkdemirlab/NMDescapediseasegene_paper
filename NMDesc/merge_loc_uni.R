#combine uni with loc
library(data.table)
library(dplyr)
library(stringr)
plus1_control_variants0406 <- read_csv("~/Downloads/0406 2/plus1_control_gnomAD_variants0406.csv")
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

