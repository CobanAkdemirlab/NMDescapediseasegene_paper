### go to R
library(aenmd)
library(GenomicRanges)
Run_Annotate_and_Save <- function(vcf_rng_proc, outputname) {
  vcf_rng_an  <- aenmd::annotate_nmd(vcf_rng_proc,rettype="gr")
  saveRDS(vcf_rng_an, file= paste0(outputname, '.rds'))
  #return(vcf_rng_an)
}
chrs <- c(1:22)
for (chr in chrs) {
print(chr)
chr.no <- paste('chr',chr,sep='')
link=paste0('https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.',chr.no,'.vcf.bgz')
system(paste0('wget ',link,sep=' '))

### get the header of a vcf file

infile <- paste0("gnomad.exomes.v4.1.sites.",chr.no,'.vcf.bgz')
vcffile <- paste0("gnomad.exomes.v4.1.sites.",chr.no,'.vcf')
cutfile <- paste0("gnomad.exomes.v4.1.sites.",chr.no,'.cut.vcf')
truncfile <- paste0("gnomad.exomes.v4.1.sites.",chr.no,'.cut.trunc.vcf')
trunc1file <- paste0("gnomad.exomes.v4.1.sites.",chr.no,'.cut.trunc.1.vcf')


system(paste0("gunzip -c ", infile, " | grep '#' > header.vcf"))
system(paste0("gunzip -c ", infile, " > " , vcffile))
system(paste0("grep -v '#' ", vcffile, " | grep 'PASS' > " , cutfile))
system(paste0("grep -e 'stop_gained' -e 'frameshift' ", cutfile, " > " , truncfile))
system(paste0("cat header.vcf ", truncfile, " > " , trunc1file))








vcf_file <- trunc1file
vcf      <- aenmd:::parse_vcf_VariantAnnotation(vcf_file,param = VariantAnnotation::ScanVcfParam(geno=NA))
outputname <- chr.no

vcf_rng <- vcf$vcf_rng
vcf_rng$ref     <- vcf_rng$ref |> Biostrings::DNAStringSet()
vcf_rng$alt     <- vcf_rng$alt |> Biostrings::DNAStringSet()

seqlevels(vcf_rng) <- gsub('chr','',seqlevels(vcf_rng))
seqnames(vcf_rng) <- gsub('chr','',seqnames(vcf_rng))

seqlevelsStyle(vcf_rng) <- 'NCBI' #gives you 1.....X
#seqlevelsStyle(vcf_rng) <- 'UCSC' #gives you chr1. ....X

#Funciton to make unique IDs
make_keys <- function(vcf_rng){
  starts <- GenomicRanges::start(vcf_rng) |> stringr::str_pad(9L, pad="0")
  keys <- paste0(GenomicRanges::seqnames(vcf_rng), ":", starts,"|" ,vcf_rng$ref, "|", vcf_rng$alt)
  return(keys)
}
#Make keys, which are unique IDs
vcf_rng$key <- aenmd:::make_keys(vcf_rng)
#vcf_ifo <- VariantAnnotation::info(vcf)

genome(vcf_rng) <- gsub('gnomAD_','',genome(vcf_rng))

vcf_rng_proc <- aenmd::process_variants(vcf_rng)


Run_Annotate_and_Save(vcf_rng_proc,outputname)
rm(vcffile)

#- split ranges and info
#vcf_rng <- SummarizedExperiment::rowRanges(vcf)
#colnames( vcf_rng |> S4Vectors::mcols() ) <- vcf_rng |> S4Vectors::mcols() |> 
  #colnames() |> janitor::make_clean_names()

#Some variants have the alt of "<DEL>" for some weird reason, just make those alt alleles empty:
#vcf_rng$alt <- ifelse(vcf_rng$alt == "<DEL>", "", vcf_rng$alt)
#make sure that ref and alt are class biostring:




}

