library(GenomicRanges)
chrs <- 1:22
for (chr in chrs) {
  print(chr)
  chr.no <- paste('chr', chr, sep = '')
  link <- paste0('https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.', chr.no, '.vcf.bgz')
  system(paste0('wget ', link, sep=' '))
  
  infile <- paste0("gnomad.exomes.v4.1.sites.", chr.no, '.vcf.bgz')
  synfile <- paste0("gnomad.exomes.v4.1.sites.", chr.no, '.cut.syn.vcf')      
  
  # Extract header directly from compressed file
  system(paste0("gunzip -c ", infile, " | grep '#' > header.vcf"))
  
  # Stream filter: PASS + synonymous_variant (NO large .vcf created)
  system(paste0(
    "gunzip -c ", infile,
    " | grep -v '#' | grep 'PASS' | grep 'synonymous_variant' > body.vcf"
  ))
  
  # Combine header + filtered variants
  system(paste0("cat header.vcf body.vcf > ", synfile))
  
  # Clean temporary body file
  system("rm body.vcf")
}

