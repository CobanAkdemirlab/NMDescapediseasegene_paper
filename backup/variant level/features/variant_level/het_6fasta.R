# Load biomaRt
library(biomaRt)


# Define your genes
genes <- c("EZH2", "FBN1", "KAT6B", "FGF14", "RARB", "SIX3", "KCNQ2")
 genes %in% omim_AD_symbols
# Retrieve canonical transcript information
canonical_tx <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "ensembl_transcript_id", 
                 "transcript_is_canonical"),
  filters = "hgnc_symbol",
  values = genes,
  mart = mart
)

# Filter only canonical (transcript_is_canonical = 1)
canonical_tx <- canonical_tx[which(canonical_tx$transcript_is_canonical == 1), ]

# Show final canonical list
canonical_tx

tx6 = all_variants[which(all_variants$transcript %in% canonical_tx$ensembl_transcript_id),]

#add hgncsymbol
tx6$hgnc_symbol = canonical_tx$hgnc_symbol[match(tx6$transcript, canonical_tx$ensembl_transcript_id)]

write_csv(tx6, "/Users/jxu14/Downloads/genes7_canonical_variants.csv")

library(dplyr)
library(stringr)
library(readr)
library(fs)

# ----------------------------
# 1) Parse tx6$key -> chr, pos, ref, alt
# ----------------------------
# Expected examples:
#   "5:041739447|CTCAA|C"
#   "chr5:41739447|CTCAA|C"
# Adjust the regex if your keys are different.
parse_key <- function(x) {
  # remove "chr" if present; capture CHR:POS|REF|ALT
  x2 <- str_replace(toupper(x), "^CHR", "")
  m  <- str_match(x2, "^([0-9XYMT]+):0*([0-9]+)\\|([ACGTN]+)\\|([ACGTN]+)$")
  tibble(
    key = x,
    chr = m[,2],
    pos = suppressWarnings(as.integer(m[,3])),
    ref = m[,4],
    alt = m[,5]
  )
}

tx6_parsed <- parse_key(tx6$key)

# ----------------------------
# 2) Parse FASTA filenames in the directory
# ----------------------------
fasta_dir <- "/Users/jxu14/Downloads/Archive 4 (1)/plus1_fasta_output"

# List all .fasta (non-recursive; set recursive=TRUE if there are subfolders)
fasta_files <- dir(fasta_dir, pattern = "[.]fasta$", full.names = TRUE, recursive = FALSE)

# Filenames are assumed like "<CHR>_<zero*POS>_<REF>_<ALT>.fasta"
# Be robust to zero padding on position.
parse_fname <- function(path) {
  fname <- path_file(path)
  # Capture groups before ".fasta"
  # Examples this matches: "17_00770305_TGA_T.fasta", "3_01014662_CAA_C.fasta"
  m <- str_match(fname, "^([0-9XYMT]+)_0*([0-9]+)_([A-Z]+)_([A-Z]+)[.]FASTA$" )
  tibble(
    file_path = path,
    file_name = fname,
    chr = m[,2],
    pos = suppressWarnings(as.integer(m[,3])),
    ref = m[,4],
    alt = m[,5]
  )
}

files_df <- bind_rows(lapply(fasta_files, parse_fname)) %>%
  filter(!is.na(chr), !is.na(pos), !is.na(ref), !is.na(alt))

# ----------------------------
# 3) Join tx6 rows to FASTA paths
# ----------------------------
tx6_with_fasta <- tx6 %>%
  # add parsed columns
  bind_cols(tx6_parsed %>% select(chr, pos, ref, alt)) %>%
  # normalize REF/ALT to upper just in case
  mutate(ref = toupper(ref), alt = toupper(alt)) %>%
  left_join(files_df, by = c("chr","pos","ref","alt"))

# Inspect unmatched (if any)
unmatched <- tx6_with_fasta %>% filter(is.na(file_path))
if (nrow(unmatched) > 0) {
  message("WARNING: ", nrow(unmatched), " variants did not find a matching FASTA file. ",
          "Check chr/pos/ref/alt or filename pattern/zero padding.")
  print(unmatched %>% select(key, chr, pos, ref, alt) %>% head(20))
}

# Show matched rows with their fasta paths
tx6_with_fasta_selected <- tx6_with_fasta %>%
  select(key, transcript, hgnc_symbol, chr, pos, ref, alt, file_path)

print(tx6_with_fasta_selected, n = 50)

# ----------------------------
# 4) (Optional) Copy matched FASTAs to a new output folder
# ----------------------------
# out_dir <- "/Users/jxu14/Downloads/plus1_fasta_matched"
# dir_create(out_dir)
# to_copy <- tx6_with_fasta_selected %>% filter(!is.na(file_path)) %>% pull(file_path) %>% unique()
# file_copy(to_copy, out_dir, overwrite = FALSE)
# message("Copied ", length(to_copy), " FASTA files to: ", out_dir)

