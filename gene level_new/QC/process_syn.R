#This script is to read in synonmous variants and get it's count in certain regions
library(data.table)
library(stringr)
library(biomaRt)

#read in all variants in syn folder
syn_dir <- "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/gnomad/syn"
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
syn_all = read.csv("/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/clinvar/syn_all_can_only.csv")
#get syn variants in certain region
##input: chrom, region.start region.end, output: syn.count
get_syn_count = function(chrom, region.start, region.end) {
  syn.sub = syn_all[
    syn_all$CHROM == chrom &
      syn_all$cds_pos >= region.start &
      syn_all$cds_pos < region.end,
  ]
  syn.count = nrow(syn.sub)
  return(syn.count)
}
