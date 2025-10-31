plus1_can_uni = getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = plus1_can$gene,
  mart = ensembl
)
plus1_css_uni = getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = plus1_css$gene,
  mart = ensembl
)
plus1_long_uni = getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = plus1_long$gene,
  mart = ensembl
)
plus1_can_uniprot = unique(plus1_can_uni$uniprotswissprot)
plus1_css_uniprot = unique(plus1_css_uni$uniprotswissprot)
plus1_long_uniprot = unique(plus1_long_uni$uniprotswissprot)


convert_to_c <- function(x) {
  x <- gsub("\\[|\\]", "", x) # Remove square brackets
  if (x == "") return(c())     # Handle empty brackets
  as.numeric(trimws(unlist(strsplit(x, ",")))) # Split by comma, trim whitespace, and convert to numeric
}
plus1.can.tr.set = getBM(
  attributes = c("transcript_is_canonical",'ensembl_transcript_id','hgnc_symbol'),
  filters = "hgnc_symbol",
  values = plus1_can_uni$hgnc_symbol,
  mart = ensembl
)
plus1.css.tr.set = getBM(
  attributes = c("transcript_is_canonical",'ensembl_transcript_id','hgnc_symbol'),
  filters = "hgnc_symbol",
  values = plus1_css_uni$hgnc_symbol,
  mart = ensembl
)
plus1.long.tr.set = getBM(
  attributes = c("transcript_is_canonical",'ensembl_transcript_id','hgnc_symbol'),
  filters = "hgnc_symbol",
  values = plus1_long_uni$hgnc_symbol,
  mart = ensembl
)
plus1.can.tr.set = plus1.can.tr.set[which(plus1.can.tr.set$transcript_is_canonical == 1),]
plus1.css.tr.set = plus1.css.tr.set[which(plus1.css.tr.set$transcript_is_canonical == 1),]
plus1.long.tr.set = plus1.long.tr.set[which(plus1.long.tr.set$transcript_is_canonical == 1),]
plus1.can.exon.set = getBM(
  attributes = c("ensembl_transcript_id","cds_start","cds_end","uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = plus1.can.tr.set$ensembl_transcript_id,
  mart = ensembl
)
plus1.css.exon.set = getBM(
  attributes = c("ensembl_transcript_id","cds_start","cds_end","uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = plus1.css.tr.set$ensembl_transcript_id,
  mart = ensembl
)
plus1.long.exon.set = getBM(
  attributes = c("ensembl_transcript_id","cds_start","cds_end","uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = plus1.long.tr.set$ensembl_transcript_id,
  mart = ensembl
)
plus1_can_re_1_flag = rep(0,length(plus1_can_uniprot))
plus1_can_re_2_flag = rep(0,length(plus1_can_uniprot))
plus1_css_re_1_flag = rep(0,length(plus1_css_uniprot))
plus1_css_re_2_flag = rep(0,length(plus1_css_uniprot))
plus1_long_re_1_flag = rep(0,length(plus1_long_uniprot))
plus1_long_re_2_flag = rep(0,length(plus1_long_uniprot))
for (i in 2:length(plus1_can_uniprot)) {
  tryCatch({
    # Get residues location
    re_1_ind = which(human_1_$uniprot1 == plus1_can_uniprot[i])
    re_2_ind = which(human_1_$uniprot2 == plus1_can_uniprot[i])
    re_1_set = human_1_$interface_residues1[re_1_ind]
    re_2_set = human_1_$interface_residues2[re_2_ind]
    re_1 = unlist(lapply(re_1_set, convert_to_c))*3
    re_2 = unlist(lapply(re_2_set, convert_to_c))*3
    # Get exon information
    hgnc = plus1_can_uni$hgnc_symbol[which(plus1_can_uni$uniprotswissprot == plus1_can_uniprot[i])]
    ## Get canonical transcript ID, then find exon start and end
    plus1.can.tr = plus1.can.tr.set[which(plus1.can.tr.set$hgnc_symbol == hgnc), 2][1]
    plus1.can.exon = plus1.can.exon.set[which(plus1.can.exon.set$ensembl_transcript_id == plus1.can.tr), ]
    
    # Remove rows with NA
    plus1.can.exon = plus1.can.exon[complete.cases(plus1.can.exon), ]
    
    # Use exon information to check if the residues are in the NMD escape region
    exon.length <- (plus1.can.exon$cds_end - plus1.can.exon$cds_start) + 1
    
    ## Rule 1: The PTC located in the last coding exon
    ## Rule 2: The PTC located within d_pen bp upstream of the penultimate exon boundary
    pen.length = exon.length[length(exon.length) - 1]
    if (pen.length < bp) {
      aft.NMD.ind = plus1.can.exon$cds_start[length(exon.length) - 1]
    } else {
      aft.NMD.ind <- sum(exon.length[1:(length(exon.length) - 1)]) - bp
    }
    
    NMD_status_can.1 = re_1 >= aft.NMD.ind
    NMD_status_can.2 = re_2 >= aft.NMD.ind
    if (sum(NMD_status_can.1) > 0) {
      plus1_can_re_1_flag[i] = 1
    }
    if (sum(NMD_status_can.2) > 0) {
      plus1_can_re_2_flag[i] = 1
    }
  }, error = function(e) {
    # Handle error gracefully (e.g., log it, skip iteration)
    message(sprintf("Error at index %d: %s", i, e$message))
  })
}

for (i in 2:length(plus1_css_uniprot)) {
  tryCatch({
    # Get residues location
    re_1_ind = which(human_1_$uniprot1 == plus1_css_uniprot[i])
    re_2_ind = which(human_1_$uniprot2 == plus1_css_uniprot[i])
    re_1_set = human_1_$interface_residues1[re_1_ind]
    re_2_set = human_1_$interface_residues2[re_2_ind]
    re_1 = unlist(lapply(re_1_set, convert_to_c))*3
    re_2 = unlist(lapply(re_2_set, convert_to_c))*3
    #rule 3 The PTC located within d_css bp downstream of the coding start site (css rule default: d_css = 150)
    NMD_status_css.1 = re_1<=150
    NMD_status_css.2 = re_2<=150
    if (sum(NMD_status_css.1) > 0) {
      plus1_css_re_1_flag[i] = 1
    }
    if (sum(NMD_status_css.2) > 0) {
      plus1_css_re_2_flag[i] = 1
    }
  }, error = function(e) {
    # Handle error gracefully (e.g., log it, skip iteration)
    message(sprintf("Error at index %d: %s", i, e$message))
  })
}

for (i in 2:length(plus1_long_uniprot)) {
  tryCatch({
    # Get residues location
    re_1_ind = which(human_1_$uniprot1 == plus1_long_uniprot[i])
    re_2_ind = which(human_1_$uniprot2 == plus1_long_uniprot[i])
    re_1_set = human_1_$interface_residues1[re_1_ind]
    re_2_set = human_1_$interface_residues2[re_2_ind]
    re_1 = unlist(lapply(re_1_set, convert_to_c))*3
    re_2 = unlist(lapply(re_2_set, convert_to_c))*3
    # Get exon information
    hgnc = plus1_long_uni$hgnc_symbol[which(plus1_long_uni$uniprotswissprot == plus1_long_uniprot[i])]
    ## Get longonical transcript ID, then find exon start and end
    plus1.long.tr = plus1.long.tr.set[which(plus1.long.tr.set$hgnc_symbol == hgnc), 2][1]
    plus1.long.exon = plus1.long.exon.set[which(plus1.long.exon.set$ensembl_transcript_id == plus1.long.tr), ]
    
    # Remove rows with NA
    plus1.long.exon = plus1.long.exon[complete.cases(plus1.long.exon), ]
    
    # Use exon information to check if the residues are in the NMD escape region
    exon.length <- (plus1.long.exon$cds_end - plus1.long.exon$cds_start) + 1
    
    #rule 4 The PTC located within an exon spanning more than 407bp 
    ##find which exon>407
    exon.long <- which(exon.length>407)
    #exclude the last exon(remove length(exon.length) from exon.long)
    exon.long <- exon.long[exon.long != length(exon.length)]
    #exclude the first exon
    exon.long <- exon.long[exon.long != 1]
    ##find if the PTC is in that long exon
    long.start = plus1.can.exon$cds_start[exon.long]
    long.end = plus1.can.exon$cds_end[exon.long]
    if(length(exon.long)==0){
      long.start = 0
      long.end = 0
    }
    NMD_status_long.1 = rep(FALSE,length(re_1))
    NMD_status_long.2 = rep(FALSE,length(re_2))
    for(k in 1:length(re_1)){
      NMD_status_long.1[k] = re_1[k]>=long.start & re_1[k]<=long.end
    }
    for(k in 1:length(re_2)){
      NMD_status_long.2[k] = re_2[k]>=long.start & re_2[k]<=long.end
    }
    if (sum(NMD_status_long.1) > 0) {
      plus1_long_re_1_flag[i] = 1
    }
    if (sum(NMD_status_long.2) > 0) {
      plus1_long_re_2_flag[i] = 1
    }
  }, error = function(e) {
    # Handle error gracefully (e.g., log it, skip iteration)
    message(sprintf("Error at index %d: %s", i, e$message))
  })
}
plus1_can_matched = plus1_can_uniprot[which(plus1_can_re_1_flag > 0 | plus1_can_re_2_flag > 0)]
length(plus1_can_matched)
plus1_css_matched = plus1_css_uniprot[which(plus1_css_re_1_flag > 0 | plus1_css_re_2_flag > 0)]
length(plus1_css_matched)
plus1_long_matched = plus1_long_uniprot[which(plus1_long_re_1_flag > 0 | plus1_long_re_2_flag > 0)]
length(plus1_long_matched)
