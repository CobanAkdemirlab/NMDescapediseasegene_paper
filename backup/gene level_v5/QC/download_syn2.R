library(data.table)
library(stringr)

for(i in 1:22){
  
  chr.name <- paste0("gnomad.exomes.v4.1.sites.chr", i, ".cut.syn.vcf")
  
  chr.syn <- fread(
    cmd = paste("grep -v '^#'", shQuote(chr.name)),
    sep = "\t",
    header = FALSE,
    fill = TRUE,
    showProgress = TRUE
  )
  
  chr.syn <- chr.syn[, 1:8]
  setnames(chr.syn, c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"))
  
  syn <- chr.syn
  rm(chr.syn)
  
  syn[, key := sprintf("%s:%09d|%s|%s", CHROM, POS, REF, ALT)]
  syn[, is_syn := str_detect(INFO, "synonymous_variant")]
  syn <- syn[is_syn == TRUE]
  
  syn[, HGVSc := str_extract(INFO, "ENST[0-9]+\\.[0-9]+:c\\.[^|;]+")]
  syn[, transcript_id := str_extract(HGVSc, "ENST[0-9]+")]
  
  syn[, INFO := NULL]
  
  write.csv(
    syn,
    paste0("gnomad.exomes.v4.1.sites.chr", i, ".cut.syn.csv"),
    row.names = FALSE
  )
  
  rm(syn)
}