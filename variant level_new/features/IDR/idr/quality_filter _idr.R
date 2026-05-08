minus1_control_idr <- read_delim("minus1_control_idr.txt", 
                                 +     delim = "\t", escape_double = FALSE, 
                                 +     trim_ws = TRUE)

library(readr)

# Define input .txt files and output .txt paths
groups_txt <- list(
  minus1_control = "minus1_control_idr.txt",
  minus1 = "minus1_idr.txt",
  plus1_control = "plus1_control_idr.txt",
  plus1 = "plus1_idr.txt",
  snv_control = "snv_control_idr.txt",
  snv = "snv_idr.txt"
)

output_dir <- "~/Downloads/idr_txt/"

for (group_name in names(groups_txt)) {
  # Read the .txt file
  file_path <- groups_txt[[group_name]]
  idr_data <- read_delim(file_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # Remove duplicated rows
  idr_data_uni <- idr_data[!duplicated(idr_data), ]
  
  # Compute nmdesc_length2 and filter
  idr_data_uni$nmdesc_length2 <- idr_data_uni$NMDesc.end - idr_data_uni$NMDesc.start + 1
  idr_data_uni <- idr_data_uni[idr_data_uni$nmdesc_length2 == idr_data_uni$nmdesc_length, ]
  
  # Write filtered data to .txt file
  out_path <- file.path(output_dir, paste0(group_name, "_idr0603_uni.txt"))
  write_delim(idr_data_uni, out_path, delim = "\t")
}
