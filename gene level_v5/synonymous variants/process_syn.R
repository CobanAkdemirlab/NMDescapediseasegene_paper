#This script is to read in synonymous variants and get their count in certain regions
#Input: gnomad.exomes.v4.1.syn.mane.AFlt1pct.csv made by parse_gnomad_syn_vcf.R
library(data.table)
library(stringr)
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

#read in the single synonymous file (replaces reading every csv/tsv/txt in a folder)
SYN_CSV    <- "gnomad.exomes.v4.1.syn.mane.AFlt1pct.csv"
GNOMAD_DIR <- "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/gnomad"
syn_path <- tryCatch(data_file(SYN_CSV),
                     error = function(e) file.path(GNOMAD_DIR, SYN_CSV))
if (!file.exists(syn_path)) stop("Cannot find ", SYN_CSV, " -- run parse_gnomad_syn_vcf.R first")
syn_all <- fread(syn_path, showProgress = FALSE)
syn_all[, source_file := basename(syn_path)]

#canonical filter: no longer needed.
#The Hail step already kept only MANE Select transcripts, which are the Ensembl
#canonical transcripts for protein-coding genes, so the biomaRt/getBM step is dropped.
syn_all_can_only <- syn_all

#make chr1 -> 1 and as.numeric
syn_all_can_only$CHROM = as.numeric(gsub("chr", "", syn_all_can_only$CHROM))
#get cds.loc, add the whole back
syn_all_can_only[, hgvsc_value := sub(".*:(c\\.[^ ]+)", "\\1", HGVSc)]
syn_all_can_only[, cds_pos := as.integer(str_extract(hgvsc_value, "(?<=c\\.)-?[0-9]+"))]
cat("Variants:", nrow(syn_all_can_only),
    " | missing cds_pos:", sum(is.na(syn_all_can_only$cds_pos)), "\n")
write.csv(syn_all_can_only, "syn_all_can_only_0928.csv", row.names = FALSE)
rm(syn_all, syn_all_can_only)
syn_all = read.csv("syn_all_can_only_0928.csv")
#get syn variants in certain region
##input: chrom, region.start region.end,transcript output: syn.count
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
