#This script is to read in synonmous variants 
library(data.table)
library(stringr)
library(biomaRt)
for(i in 1:22){
  chr.name = paste0("gnomad.exomes.v4.1.sites.chr",i,".cut.syn.vcf")
  chr.syn <- fread(
    cmd = paste("grep -v '^#'", shQuote(chr.name)),
    sep = "\t",
    header = FALSE,
    fill = TRUE,
    showProgress = TRUE
  )
  chr.syn <- chr.syn[, 1:8]
  setnames(chr.syn, c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"))
  syn = chr.syn
  rm(chr.syn)
  # key: chr:000963993|C|T  (POS 固定 9 位补0)
  syn[, key := sprintf("%s:%09d|%s|%s", CHROM, POS, REF, ALT)]
  syn[, is_syn := str_detect(INFO, "synonymous_variant")]
  #check for synonymous variants
  #table(syn$is_syn, useNA = "ifany")
  syn = syn[is_syn == TRUE]
  #get transcript id, only keep enst
  syn[, transcript_id := str_extract(INFO, "ENST[0-9]+")]
  #get cds loc
  syn[, HGVSc := {
    enst <- transcript_id
    
    # get vep=... block (until ';')
    vep_blob <- sub(".*\\bvep=([^;]+).*", "\\1", INFO)
    vep_blob[vep_blob == INFO] <- NA_character_  # no vep= -> NA
    
    mapply(function(vb, e) {
      if (is.na(vb) || is.na(e)) return(NA_character_)
      items <- strsplit(vb, ",", fixed = TRUE)[[1]]
      hit <- items[grepl(e, items, fixed = TRUE)]
      if (length(hit) == 0) return(NA_character_)
      
      # prefer ENST*.ver:c....
      mm <- str_match(hit[1], "(ENST[0-9]+\\.[0-9]+:c\\.[^|]+)")
      if (!is.na(mm[1,2])) return(mm[1,2])
      
      # fallback: sometimes NM_...:c....
      mm2 <- str_match(hit[1], "((ENST|NM_)[^|]*:c\\.[^|]+)")
      if (!is.na(mm2[1,2])) return(mm2[1,2])
      
      NA_character_
    }, vep_blob, enst, USE.NAMES = FALSE)
  }]
  #delete info
  syn = syn[, -'INFO']
  write.csv(syn, paste0("gnomad.exomes.v4.1.sites.chr",i,".cut.syn.csv"), row.names = FALSE)
  rm(syn)
}