#This script is to read in synonmous variants and get it's count in certain regions
library(data.table)
library(stringr)
library(biomaRt)
# --- Path resolution layer -------------------------------------------------
# Resolves paths via paths.R instead of absolute paths
#   data_file("x.csv") locates by filename, errors if not found
#   out_file("y.csv") writes to NMDESC_OUT (default ~/Desktop/NMDesc_out)
#   data_root("clinvar") use when a directory is needed instead of a file
.p <- c("gene level_v5/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../lib/paths.R", "../../../../lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p)) stop("Could not find paths.R -- run R from the repository root")
source(.p[1]); rm(.p)
# --------------------------------------------------------------------------


#read in all variants in syn folder
syn_dir <- data_root("clinvar")
files <- list.files(
  syn_dir,
  pattern = "\\.(csv|tsv|txt)$",
  full.names = TRUE
)
read_one <- function(f) {
  is_csv <- grepl("\\.csv$", f, ignore.case = TRUE)
  dt <- fread(f, sep = if (is_csv) "," else "\t", showProgress = FALSE)
  dt[, source_file := basename(f)]
  dt
}
#combine them into syn_all
syn_all <- rbindlist(lapply(files, read_one), fill = TRUE, use.names = TRUE)
#filter for canonical transcript using getBM
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
syn_all_tx_set = unique(na.omit(syn_all$transcript_id))
syn_all_can = getBM(attributes = c("ensembl_transcript_id", "transcript_is_canonical"),
                    filters = "ensembl_transcript_id",
                    values = syn_all_tx_set,
                    mart = ensembl)
syn_all_can_set = syn_all_can[which(syn_all_can$transcript_is_canonical == 1), 'ensembl_transcript_id']
syn_all_can_only = syn_all[transcript_id %in% syn_all_can_set]

#make chr1 -> 1 and as.numeric
syn_all_can_only$CHROM = as.numeric(gsub("chr", "", syn_all_can_only$CHROM))
#get cds.loc, add the whole back
syn_all_can_only[, hgvsc_value := sub(".*:(c\\.[^ ]+)", "\\1", HGVSc)]
syn_all_can_only[, cds_pos := as.integer(str_extract(hgvsc_value, "(?<=c\\.)-?[0-9]+"))]
write.csv(syn_all_can_only, "syn_all_can_only.csv", row.names = FALSE)
rm(syn_all, syn_all_can, syn_all_can_set)
rm(syn_all_tx_set)
syn_all = read.csv("syn_all_can_only.csv")
#get syn variants in certain region
##input: chrom, region.start region.end, output: syn.count
.syn_index <- new.env(parent = emptyenv())
get_syn_count = function(chrom, region.start, region.end, transcript) {
  if (missing(transcript) || is.null(transcript) || is.na(transcript))
    stop("get_syn_count() needs the transcript ID so it counts only that transcript's variants")
  if (!exists("syn_all")) stop("syn_all is not loaded -- run process_syn.R first")
  # Split positions by transcript once per syn_all, then reuse
  key <- nrow(syn_all)
  if (!identical(.syn_index$key, key)) {
    ids <- sub("\\.[0-9]+$", "", as.character(syn_all$transcript_id))
    pos <- suppressWarnings(as.integer(syn_all$cds_pos))
    ok  <- !is.na(ids) & !is.na(pos)
    .syn_index$by_tx <- split(pos[ok], ids[ok])
    .syn_index$key   <- key
  }
  p <- .syn_index$by_tx[[sub("\\.[0-9]+$", "", transcript)]]
  if (is.null(p)) return(0L)
  sum(p >= region.start & p <= region.end)
}