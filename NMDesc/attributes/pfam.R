#change how we define long exon.

#double check single exon gene

#use clinvar_can as example
#A. prepare the data
library(dplyr)
library(PFAM.db)
library(biomaRt)
library(UpSetR)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
pfam_final = data.frame(gene = character(),pfam = character(),ana_type = character(), pfam_type = character())
map_cds_to_chromosome <- function(cds_position, exon_data, strand) {
  cumulative_length <- 0
  
  for (i in seq_len(nrow(exon_data))) {
    exon <- exon_data[i, ] #exon contain the current checking exon info
    exon_length <- exon$cds_end - exon$cds_start + 1 #coding sequence length of this exon
    
    # until the pfam_cds_position falls within the current exon
    if (cumulative_length + exon_length >= cds_position) {
      # Position within the exon
      position_in_exon <- cds_position - cumulative_length #loc of cds - sum of cds from the passed exon = number of bp from the start of the coding start side exon
      
      # Determine chromosomal position based on strand
      if (strand == 1) {  # positive strand
        chromosomal_position <- exon$exon_chrom_start + position_in_exon - 1 #this assumes that the start of the exon is the start of the coding sequence, which is wrong
        if(exon$rank == 1){
          chromosomal_position = chromosomal_position + exon$cdna_coding_start 
        }
      } else {  # negative strand
        chromosomal_position <- exon$exon_chrom_end - position_in_exon + 1
        #if the pfam_cds falls in the first exon, correct for the UTR at the first exon
        if(exon$rank == 1){
          chromosomal_position = chromosomal_position - exon$cdna_coding_start 
        }
      }
      
      return(chromosomal_position) #stop the loop
    }
    cumulative_length <- cumulative_length + exon_length #sum of the coding sequence length of exons passed
  }
  
  return(NA)  # If not found within exons
  
}
calculate_overall_overlap_percentage <- function(pfam_region, NMDesc_region) {
  total_overlap_length <- 0
  total_shorter_length <- 0
  
  for (i in 1:nrow(pfam_region)) {
    pfam_start <- pfam_region$start[i]
    pfam_end <- pfam_region$end[i]
    
    for (j in 1:nrow(NMDesc_region)) {
      NMDesc_start <- NMDesc_region$start[j]
      NMDesc_end <- NMDesc_region$end[j]
      
      # Calculate overlap
      overlap_start <- max(pfam_start, NMDesc_start)
      overlap_end <- min(pfam_end, NMDesc_end)
      
      if (overlap_start < overlap_end) {  # There is an overlap
        overlap_length <- overlap_end - overlap_start + 1
        pfam_length <- pfam_end - pfam_start + 1
        NMDesc_length <- NMDesc_end - NMDesc_start + 1
        shorter_length <- min(pfam_length, NMDesc_length)
        
        # Add to totals
        total_overlap_length <- total_overlap_length + overlap_length
        total_shorter_length <- total_shorter_length + shorter_length
      }
    }
  }
  
  # Calculate overall overlap percentage
  overall_overlap_percentage <- (total_overlap_length / total_shorter_length) * 100
  return(overall_overlap_percentage)
}

get_overlap = function(gene_name='NOVA2',type = 'can',cut_off = 0.5){
  #1. get pfam region
  ##1.1 get pfam loc of amino acid and transcript info
  overlap_flag = 0
  pfam_data <- getBM(
    attributes = c('hgnc_symbol',"ensembl_gene_id","pfam", "pfam_start", "pfam_end","transcript_is_canonical","ensembl_transcript_id"),
    filters = "hgnc_symbol",
    values = gene_name,
    mart = ensembl
  )
  pfam_data = pfam_data[which(pfam_data$transcript_is_canonical==1),]
  if(sum(!is.na(pfam_data$pfam)) == 0){ #if this gene contain no pfam domain,return 0
    return(overlap_flag)
  }
  pfam_type = sapply(pfam_data$pfam, function(x) PFAMDE[[x]])
  pfam_output = data.frame(gene = rep(gene_name,length(pfam_type)),pfam = pfam_data$pfam, ana_type = rep(type,length(pfam_type)), pfam_type)
  pfam_output =  pfam_output[which(pfam_output$pfam_type != 'NA'),]
  pfam_final <<- rbind(pfam_final,pfam_output)
  transcript_id <- unique(pfam_data$ensembl_transcript_id) #this is the canonical transcript id for the target gene
  exon_data <- getBM(
    attributes = c("cds_start", "cds_end","exon_chrom_start", "exon_chrom_end", 
                   "ensembl_exon_id",'cdna_coding_start','cdna_coding_end', "strand", "rank"),
    filters = "ensembl_transcript_id",
    values = transcript_id,
    mart = ensembl
  )#cds is transcript level,represents the start of the coding sequence within that specific exon. Exon_chrom_start and exon_chrom_end is chromosome level (absolute positions on the chromosome).
  #only keep the coding exon info(at least contain part of the coding sequence)
  exon_data <- exon_data %>%
    filter(!is.na(cds_start) & !is.na(cds_end) & !is.na(exon_chrom_start) & !is.na(exon_chrom_end)) 
  exon_data$exon_length = exon_data$exon_chrom_end - exon_data$exon_chrom_start + 1
  exon_data$cds_length = exon_data$cds_end - exon_data$cds_start + 1 #for each exon, the length on chrom might be different from the length of cds
  exon_data$cDNA_length = exon_data$cdna_coding_end - exon_data$cdna_coding_start + 1
  ##1.2 map portein level loc of pfam to transcript level loc of pfam
  pfam_data$pfam_cds_start <- pfam_data$pfam_start * 3 
  pfam_data$pfam_cds_end <- pfam_data$pfam_end * 3
  
  #1.3 map /transcript level loc of pfam/ to /chrom level loc of pfam/
  
  # Apply the function to each row in pfam_data for cds_start and cds_end
  pfam_data <- pfam_data %>%
    rowwise() %>%
    mutate(
      pfam_chromosomal_start = map_cds_to_chromosome(pfam_cds_start, exon_data, exon_data$strand[1]),
      pfam_chromosomal_end = map_cds_to_chromosome(pfam_cds_end, exon_data, exon_data$strand[1])
    ) %>%
    ungroup()
  
  #make sure start<end, if not swap
  for(j in 1:nrow(pfam_data)){
    if(pfam_data$pfam_chromosomal_start[j] > pfam_data$pfam_chromosomal_end[j]){
      temp = pfam_data$pfam_chromosomal_start[j]
      pfam_data$pfam_chromosomal_start[j] = pfam_data$pfam_chromosomal_end[j]
      pfam_data$pfam_chromosomal_end[j] = temp
    }
  }
  
  #2.get NMDesc region for snv
  if(type == 'can'){
    ##2.1 last coding exon
    exon_data = exon_data[which(!is.na(exon_data$cds_start)),]
    cds_length = exon_data$cds_end[which(exon_data$rank == max(exon_data$rank))] - exon_data$cds_start[which(exon_data$rank == max(exon_data$rank))]
    if(exon_data$strand[1] == -1){
      last_begin = exon_data$exon_chrom_end[which(exon_data$rank == max(exon_data$rank))]
      last_end = last_begin - cds_length
    }else{
      last_begin = exon_data$exon_chrom_start[which(exon_data$rank == max(exon_data$rank))]
      last_end = last_begin + cds_length
    } 
    
    ##2.2 penultimate exon, if -1 strand, start and end should swap
    if(exon_data$strand[1] == -1){
      pen_end = exon_data$exon_chrom_start[which(exon_data$rank == max(exon_data$rank)-1)]
      pen_begin = pen_end + 50
    }else{
      pen_end = exon_data$exon_chrom_end[which(exon_data$rank == max(exon_data$rank)-1)]
      pen_begin = pen_end - 50
    }
    NMDesc_region = data.frame(start = c(last_begin,pen_begin), end = c(last_end,pen_end))
  
    } else if(type == 'css'){
    ##2.3 CSS proximal exon: d_css bp downstream of the coding start site 
    if(exon_data$strand[1] == -1){
      css_begin = exon_data$exon_chrom_end[which(exon_data$rank == 1)] - exon_data$cdna_coding_start[which(exon_data$rank == 1)]
      css_end = css_begin - 150
    }else{
      css_begin = exon_data$exon_chrom_start[which(exon_data$rank == 1)] + exon_data$cdna_coding_start[which(exon_data$rank == 1)]
        #already filter only coding exon, so if this is NA, means NMDesc is no-coding
      css_end = css_begin + 150
    }
    NMDesc_region = data.frame(start = css_begin, end = css_end)
  
    } else if(type == 'long'){
    ##2.4 407plus exon, any exon longer than 407bp, change it to only coding part
    exon_data$coding_exon_length = exon_data$cds_end - exon_data$cds_start + 1
    long_rank = exon_data$rank[which(exon_data$coding_exon_length > 407)]
    #exclud any long_rank that is 1 or max(rank)
    long_rank = long_rank[which(long_rank != max(exon_data$rank) & long_rank != 1)]
    long_rank = unique(long_rank)
    long_cds_start = exon_data$cds_start[which(exon_data$rank == long_rank)]
    long_cds_start = unique(long_cds_start)
    long_cds_end = exon_data$cds_end[which(exon_data$rank == long_rank)]
    long_cds_end = unique(long_cds_end)
    long_chrom_start = 0
    long_chrom_end = 0
    if(length(long_cds_start) != 0){
      long_chrom_start = map_cds_to_chromosome(long_cds_start, exon_data, exon_data$strand[1])
      long_chrom_end = map_cds_to_chromosome(long_cds_end, exon_data, exon_data$strand[1])
      exon_long_start = exon_data$exon_chrom_start[which(exon_data$exon_length > 407)]
      }
    NMDesc_region = data.frame(start = long_chrom_start, end = long_chrom_end)
    }else if(type == 'all'){
      NMDesc_region = data.frame(start = c(last_begin,pen_begin,css_begin,long_chrom_start), end = c(last_end,pen_end,css_end,long_chrom_end))
  }
  
  #if NMDesc_region is empty, return 0
  if(nrow(NMDesc_region) == 0){
    return(0)
  }
  
  ##2.5 make start < end always
  for(j in 1:nrow(NMDesc_region)){
    if(NMDesc_region$start[j] > NMDesc_region$end[j]){
      temp = NMDesc_region$start[j]
      NMDesc_region$start[j] = NMDesc_region$end[j]
      NMDesc_region$end[j] = temp
    }
  }
  
  #3. find the overlap percentage of NMDesc_region and pfam_region
  pfam_region = data.frame(start = pfam_data$pfam_chromosomal_start, end = pfam_data$pfam_chromosomal_end)
  
  overall_overlap_percentage <- calculate_overall_overlap_percentage(pfam_region, NMDesc_region)
  if (!is.na(overall_overlap_percentage) && overall_overlap_percentage > cut_off*100) {
    overlap_flag <- 1
  }
  return(overlap_flag)
}
get_mean_overlap = function(genelist, list_type = 'can'){
  overlap_flag = rep(0,length(genelist))
  
  for(i in 1:length(genelist)){
    #print(i)
    overlap_flag[i] <- tryCatch(
      {
        get_overlap(genelist[i],type =list_type)  # Try to get the overlap
      },
      error = function(e) {
        message("Error in get_overlap() for gene ", genelist[i], ": ", e$message)
        NA  # Return NA if there is an error
      }
    )
  }
  #write.table(overlap_flag, file = paste0('overlap_flag',object,'.txt'), quote = FALSE, row.names = FALSE, col.names = FALSE)
  mean_flag = mean(overlap_flag, na.rm = TRUE)
  return(mean_flag)
}

object_list = c('all','css','long','can','control')
getpath2 = function(input_name){
  output_name = paste0("~/Desktop/new_clinvar/snv_list/list2/snv_plp_nmd_p_NMDenriched_",input_name,".txt")
  return(output_name)
}
path_list2 = sapply(object_list, getpath2)
all_list = read.table(path_list2[1],quote="\"", comment.char="",col.names = 'gene')
css_list = read.table(path_list2[2],quote="\"", comment.char="",col.names = 'gene')
long_list = read.table(path_list2[3],quote="\"", comment.char="",col.names = 'gene')
can_list = read.table(path_list2[4],quote="\"", comment.char="",col.names = 'gene')
control_list = read.table(path_list2[5],quote="\"", comment.char="",col.names = 'gene')
#make an upsetR plot
gene_lists <- list(
  all_list = all_list$gene,
  css_list = css_list$gene,
  long_list = long_list$gene,
  can_list = can_list$gene
)

# Convert to binary matrix
gene_matrix <- fromList(gene_lists)

# Generate the UpSet plot
upset(gene_matrix, 
      sets = c('all_list', 'css_list', 'long_list', 'can_list'),
      keep.order = TRUE, 
      main.bar.color = "darkblue", 
      sets.bar.color = "lightblue", 
      order.by = "freq")

#generate omim list that does not overlap with clinvar list
all_list = read.table(path_list2[1],quote="\"", comment.char="",col.names = 'gene')
css_list = read.table(path_list2[2],quote="\"", comment.char="",col.names = 'gene')
long_list = read.table(path_list2[3],quote="\"", comment.char="",col.names = 'gene')
can_list = read.table(path_list2[4],quote="\"", comment.char="",col.names = 'gene')
omim_list = read.table(path_list2[5],quote="\"", comment.char="",col.names = 'gene')
clinvar_list = unique(c(all_list$gene,css_list$gene,long_list$gene,can_list$gene))
omim_list2 = omim_list[which(!omim_list$gene %in% clinvar_list),]


#B. run the code
for(i in 1:4){
  pfam_final = data.frame(gene = character(),pfam = character(),ana_type = character(), pfam_type = character())
  object = object_list[i]
  path = path_list2[i]
  gene_list = read.table(path,quote="\"", comment.char="",col.names = 'gene')
  print(nrow(gene_list))
  op = get_mean_overlap(gene_list$gene[211:311],list_type = object_list[i])
  cat(object,':',op,'\n')
  write.csv(pfam_final, file = paste0("~/Desktop/pfam/pfam1118_clinvar_",object,".csv"), row.names = FALSE)
}      
#for omim
path = path_list2[5]
gene_list = read.table(path,quote="\"", comment.char="",col.names = 'gene')
gene_list = gene_list[-1,]
gene_list = data.frame(gene = gene_list)
#start timing
start_time = Sys.time()
for(k in 1:4){
  pfam_final = data.frame(gene = character(),pfam = character(),ana_type = character(), pfam_type = character())
  op = get_mean_overlap(gene_list,list_type = object_list[k])
  cat(object_list[k],':',op,'\n')
  object = object_list[k]
  write.csv(pfam_final, file = paste0("~/Desktop/pfam/pfam_control_",object,".csv"), row.names = FALSE)
}
#end timing
end_time = Sys.time()
#print the time
end_time - start_time


#all : 0.2076229 ;css : 0.110331 ;long : 0.106319 ;can : 0.1404213 
#all : 0.2414487 ;css : 0.1076459 ;long : 0.1529175;can : 0.1730382 
#all : 0.2347041 ;css : 0.1233701 ;long : 0.1444333 ;can : 0.1494483 
#all : 0.2744702 ;css : 0.1440081 ;long : 0.1782477 can : 0.1993958
#all : 0.3489168 ;css : 0.1334045 ;long : 0.2292111 ;can : 0.2665245 

all_per = (0.2076229*1000+0.2414487*1000+0.2347041*1000+0.2744702*1000+0.3489168*949)/4949
css_per = (0.110331*1000+0.1076459*1000+0.1233701*1000+0.1440081*1000+0.1334045*949)/4949
long_per = (0.106319*1000+0.1529175*1000+0.1444333*1000+0.1782477*1000+0.2292111*949)/4949
can_per = (0.1404213*1000+0.1730382*1000+0.1494483*1000+0.1993958*1000+0.2665245*949)/4949


clinvar_overlap = data.frame(type = c('clinvar_all','clinvar_css','clinvar_long','clinvar_can'),percentage = c(0.23,0.14,0.30,0.15))
omim_overlap = data.frame(type = c('omim_all','omim_css','omim_long','omim_can'),percentage = c(0.26,0.12,0.16,0.19))
#use ggplot to give bar plot comparing pfam_overlap with omim_overlap
clinvar_overlap$dataset <- 'ClinVar'
omim_overlap$dataset <- 'OMIM'
combined_data <- bind_rows(clinvar_overlap, omim_overlap) %>%
  mutate(group = case_when(
    type %in% c('clinvar_all', 'omim_all') ~ 'all',
    type %in% c('clinvar_css', 'omim_css') ~ 'css',
    type %in% c('clinvar_long', 'omim_long') ~ 'long',
    type %in% c('clinvar_can', 'omim_can') ~ 'can'
  ))

# Plot with custom spacing between groups
ggplot(combined_data, aes(x = type, y = percentage, fill = dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_x_discrete(limits = c('clinvar_all', 'omim_all', '', 
                              'clinvar_css', 'omim_css', '', 
                              'clinvar_long', 'omim_long', '', 
                              'clinvar_can', 'omim_can')) +
  labs(title = "Comparison of ClinVar and OMIM Overlap Percentages with PFAM",
       x = "Type",
       y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()) 





library(ggplot2)
ggplot(pfam_overlap, aes(x = type, y = percentage)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  #geom_text(aes(label = paste0(round(percentage, 2), "%")), vjust = -0.5) +
  theme_minimal() +
  #make it horizontal
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "",
       x = "", y = "") +
  theme_bw()


#upsetR plot