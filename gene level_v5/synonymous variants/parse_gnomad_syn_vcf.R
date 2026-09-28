#This script reads in rare synonymous variants from gnomAD v4.1 exomes.
#Input is the Hail output of gnomad_synonymous_rare_hail_v2.py, which already applied:
#  PASS only, AC > 0, AF < 1%, MANE Select transcripts with only "synonymous_variant"
#So the old steps (grep VCF per chromosome, detect synonymous in INFO, pull the first ENST,
#parse HGVSc out of the vep= block) are no longer needed; the columns come straight from Hail.
library(data.table)
library(stringr)

GNOMAD_DIR <- "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/gnomad"
IN_FILE    <- file.path(GNOMAD_DIR, "gnomad_exomes_v4.1_synonymous_AFlt1pct.tsv.bgz")
OUT_FILE   <- file.path(GNOMAD_DIR, "gnomad.exomes.v4.1.syn.mane.AFlt1pct.csv")
CHROMS     <- paste0("chr", 1:22)   # same chromosomes as the old 1:22 loop

if (!file.exists(IN_FILE)) stop("Cannot find ", IN_FILE)

#quote = "" because the alleles column looks like ["C","T"]
hail <- fread(
  cmd = paste("gzip -dc", shQuote(IN_FILE)),
  sep = "\t",
  header = TRUE,
  quote = "",
  colClasses = "character",
  na.strings = c("NA", ""),
  showProgress = TRUE
)

al <- gsub('\\[|\\]|"', "", hail$alleles)
syn <- data.table(
  CHROM         = sub(":.*$", "", hail$locus),
  POS           = as.integer(sub("^.*:", "", hail$locus)),
  ID            = NA_character_,   # not exported by Hail
  REF           = sub(",.*$", "", al),
  ALT           = sub("^[^,]*,", "", al),
  QUAL          = NA_character_,   # not exported by Hail
  FILTER        = hail$FILTER,
  transcript_id = sub("\\.[0-9]+$", "", hail$transcript),
  HGVSc         = hail$HGVSc,      # e.g. ENST00000269305.9:c.1101G>A
  HGVSp         = hail$HGVSp,
  gene          = hail$gene,
  gene_id       = hail$gene_id,
  AC            = as.integer(hail$AC),
  AN            = as.integer(hail$AN),
  AF            = as.numeric(hail$AF),
  nhomalt       = as.integer(hail$nhomalt)
)
rm(hail, al)

#keep autosomes only, as before
syn <- syn[CHROM %in% CHROMS]

# key: chr1:000963993|C|T (POS zero-padded to 9 digits), same as before
syn[, key := sprintf("%s:%09d|%s|%s", CHROM, POS, REF, ALT)]
syn[, is_syn := TRUE]

#same column order as the old per-chromosome csv files, extra columns at the end
setcolorder(syn, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                   "key", "is_syn", "transcript_id", "HGVSc"))
syn <- syn[order(as.integer(sub("^chr", "", CHROM)), POS)]

#quick checks
cat("Rows:", nrow(syn), "\n")
cat("Unique variants:", uniqueN(syn$key), "\n")
cat("Missing HGVSc:", sum(is.na(syn$HGVSc)), "\n")
print(table(factor(syn$CHROM, levels = CHROMS)))

fwrite(syn, OUT_FILE)
cat("Wrote", OUT_FILE, "\n")
rm(syn)
