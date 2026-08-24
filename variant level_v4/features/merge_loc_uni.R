#combine uni with loc
library(data.table)
library(dplyr)
library(stringr)
# --- Path resolution layer (auto-inserted) ------------------------------------------------
#   data_file("x.csv") locates by filename, errors clearly if missing
#   out_file("y.csv")   writes to NMDESC_OUT (default ~/Desktop/NMDesc_out)
#   data_root("clinvar") use when a directory is needed, not a file
.p <- c("gene level_v3/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../lib/paths.R", "../../../../lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p)) stop("paths.R not found -- run R from the repository root")
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

