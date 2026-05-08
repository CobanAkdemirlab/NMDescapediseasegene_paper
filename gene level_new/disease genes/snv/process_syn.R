#This script is to read in synonmous variants and get it's count in certain regions
get_syn_count = function(chrom, region.start, region.end) {
  syn.sub = syn_all[
    syn_all$CHROM == chrom &
      syn_all$cds_pos >= region.start &
      syn_all$cds_pos < region.end,
  ]
  syn.count = nrow(syn.sub)
  return(syn.count)
}