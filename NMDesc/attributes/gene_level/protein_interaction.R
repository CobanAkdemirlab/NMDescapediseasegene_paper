plus1_can_uniprot = getBM(
       attributes = c("hgnc_symbol", "uniprotswissprot"),
       filters = "hgnc_symbol",
       values = plus1_can_gene$hgnc_symbol,
       mart = ensembl
   )
plus1_can_uniprot = unique(plus1_can_uniprot$uniprotswissprot)
length(plus1_can_uniprot)
length(which(human_1_$uniprot1 %in% plus1_can_uniprot & human_1_$uniprot2 %in% plus1_can_uniprot))


# are there any  interaction interface residues in the nMD escape region of those genes?
#for snv, just treat the residues like PTC
##get residues location
snv_can_uni = getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = snv_plp_ptc_p1120_NMDenriched2_can$X1,
  mart = ensembl
)
snv_css_uni = getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = snv_plp_ptc_p1120_NMDenriched2_css$X1,
  mart = ensembl
)
snv_long_uni = getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = snv_plp_ptc_p1120_NMDenriched2_long$X1,
  mart = ensembl
)
snv_control_uni = getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "hgnc_symbol",
  values = control$gene[-1],
  mart = ensembl
)
snv_can_uniprot = unique(snv_can_uni$uniprotswissprot)
snv_css_uniprot = unique(snv_css_uni$uniprotswissprot)
snv_long_uniprot = unique(snv_long_uni$uniprotswissprot)
snv_control_uniprot = unique(snv_control_uni$uniprotswissprot)

convert_to_c <- function(x) {
  x <- gsub("\\[|\\]", "", x) # Remove square brackets
  if (x == "") return(c())     # Handle empty brackets
  as.numeric(trimws(unlist(strsplit(x, ",")))) # Split by comma, trim whitespace, and convert to numeric
}
snv.can.tr.set = getBM(
  attributes = c("transcript_is_canonical",'ensembl_transcript_id','hgnc_symbol'),
  filters = "hgnc_symbol",
  values = snv_can_uni$hgnc_symbol,
  mart = ensembl
)
snv.css.tr.set = getBM(
  attributes = c("transcript_is_canonical",'ensembl_transcript_id','hgnc_symbol'),
  filters = "hgnc_symbol",
  values = snv_css_uni$hgnc_symbol,
  mart = ensembl
)
snv.long.tr.set = getBM(
  attributes = c("transcript_is_canonical",'ensembl_transcript_id','hgnc_symbol'),
  filters = "hgnc_symbol",
  values = snv_long_uni$hgnc_symbol,
  mart = ensembl
)
snv.control.tr.set = getBM(
  attributes = c("transcript_is_canonical",'ensembl_transcript_id','hgnc_symbol'),
  filters = "hgnc_symbol",
  values = snv_control_uni$hgnc_symbol,
  mart = ensembl
)
snv.can.tr.set = snv.can.tr.set[which(snv.can.tr.set$transcript_is_canonical == 1),]
snv.css.tr.set = snv.css.tr.set[which(snv.css.tr.set$transcript_is_canonical == 1),]
snv.long.tr.set = snv.long.tr.set[which(snv.long.tr.set$transcript_is_canonical == 1),]
snv.control.tr.set = snv.control.tr.set[which(snv.control.tr.set$transcript_is_canonical == 1),]
snv.can.exon.set = getBM(
  attributes = c("ensembl_transcript_id","cds_start","cds_end","uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = snv.can.tr.set$ensembl_transcript_id,
  mart = ensembl
)
snv.css.exon.set = getBM(
  attributes = c("ensembl_transcript_id","cds_start","cds_end","uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = snv.css.tr.set$ensembl_transcript_id,
  mart = ensembl
)
snv.long.exon.set = getBM(
  attributes = c("ensembl_transcript_id","cds_start","cds_end","uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = snv.long.tr.set$ensembl_transcript_id,
  mart = ensembl
)
snv.control.exon.set = getBM(
  attributes = c("ensembl_transcript_id","cds_start","cds_end","uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = snv.control.tr.set$ensembl_transcript_id,
  mart = ensembl
)
snv_can_re_1_flag = rep(0,length(snv_can_uniprot))
snv_can_re_2_flag = rep(0,length(snv_can_uniprot))
snv_css_re_1_flag = rep(0,length(snv_css_uniprot))
snv_css_re_2_flag = rep(0,length(snv_css_uniprot))
snv_long_re_1_flag = rep(0,length(snv_long_uniprot))
snv_long_re_2_flag = rep(0,length(snv_long_uniprot))
snv_control_re_1_flag = rep(0,length(snv_control_uniprot))
snv_control_re_2_flag = rep(0,length(snv_control_uniprot))
for (i in 2:length(snv_control_uniprot)) {
  tryCatch({
    # Get residues location
    re_1_ind = which(human_1_$uniprot1 == snv_control_uniprot[i])
    re_2_ind = which(human_1_$uniprot2 == snv_control_uniprot[i])
    re_1_set = human_1_$interface_residues1[re_1_ind]
    re_2_set = human_1_$interface_residues2[re_2_ind]
    re_1 = unlist(lapply(re_1_set, convert_to_c))*3
    re_2 = unlist(lapply(re_2_set, convert_to_c))*3
    # Get exon information
    hgnc = snv_control_uni$hgnc_symbol[which(snv_control_uni$uniprotswissprot == snv_control_uniprot[i])]
    ## Get controlonical transcript ID, then find exon start and end
    snv.control.tr = snv.control.tr.set[which(snv.control.tr.set$hgnc_symbol == hgnc), 2][1]
    snv.control.exon = snv.control.exon.set[which(snv.control.exon.set$ensembl_transcript_id == snv.control.tr), ]
    
    # Remove rows with NA
    snv.control.exon = snv.control.exon[complete.cases(snv.control.exon), ]
    
    # Use exon information to check if the residues are in the NMD escape region
    exon.length <- (snv.control.exon$cds_end - snv.control.exon$cds_start) + 1
    
    ## Rule 1: The PTC located in the last coding exon
    ## Rule 2: The PTC located within d_pen bp upstream of the penultimate exon boundary
    pen.length = exon.length[length(exon.length) - 1]
    if (pen.length < bp) {
      aft.NMD.ind = snv.control.exon$cds_start[length(exon.length) - 1]
    } else {
      aft.NMD.ind <- sum(exon.length[1:(length(exon.length) - 1)]) - bp
    }
    
    NMD_status_control.1 = re_1 >= aft.NMD.ind
    NMD_status_control.2 = re_2 >= aft.NMD.ind
    if (sum(NMD_status_control.1) > 0) {
      snv_control_re_1_flag[i] = 1
    }
    if (sum(NMD_status_control.2) > 0) {
      snv_control_re_2_flag[i] = 1
    }
  }, error = function(e) {
    # Handle error gracefully (e.g., log it, skip iteration)
    message(sprintf("Error at index %d: %s", i, e$message))
  })
}

for (i in 2:length(snv_css_uniprot)) {
  tryCatch({
    # Get residues location
    re_1_ind = which(human_1_$uniprot1 == snv_css_uniprot[i])
    re_2_ind = which(human_1_$uniprot2 == snv_css_uniprot[i])
    re_1_set = human_1_$interface_residues1[re_1_ind]
    re_2_set = human_1_$interface_residues2[re_2_ind]
    re_1 = unlist(lapply(re_1_set, convert_to_c))*3
    re_2 = unlist(lapply(re_2_set, convert_to_c))*3
    #rule 3 The PTC located within d_css bp downstream of the coding start site (css rule default: d_css = 150)
    NMD_status_css.1 = re_1<=150
    NMD_status_css.2 = re_2<=150
    if (sum(NMD_status_css.1) > 0) {
      snv_css_re_1_flag[i] = 1
    }
    if (sum(NMD_status_css.2) > 0) {
      snv_css_re_2_flag[i] = 1
    }
  }, error = function(e) {
    # Handle error gracefully (e.g., log it, skip iteration)
    message(sprintf("Error at index %d: %s", i, e$message))
  })
}

for (i in 2:length(snv_long_uniprot)) {
  tryCatch({
    # Get residues location
    re_1_ind = which(human_1_$uniprot1 == snv_long_uniprot[i])
    re_2_ind = which(human_1_$uniprot2 == snv_long_uniprot[i])
    re_1_set = human_1_$interface_residues1[re_1_ind]
    re_2_set = human_1_$interface_residues2[re_2_ind]
    re_1 = unlist(lapply(re_1_set, convert_to_c))*3
    re_2 = unlist(lapply(re_2_set, convert_to_c))*3
    # Get exon information
    hgnc = snv_long_uni$hgnc_symbol[which(snv_long_uni$uniprotswissprot == snv_long_uniprot[i])]
    ## Get longonical transcript ID, then find exon start and end
    snv.long.tr = snv.long.tr.set[which(snv.long.tr.set$hgnc_symbol == hgnc), 2][1]
    snv.long.exon = snv.long.exon.set[which(snv.long.exon.set$ensembl_transcript_id == snv.long.tr), ]
    
    # Remove rows with NA
    snv.long.exon = snv.long.exon[complete.cases(snv.long.exon), ]
    
    # Use exon information to check if the residues are in the NMD escape region
    exon.length <- (snv.long.exon$cds_end - snv.long.exon$cds_start) + 1
    
    #rule 4 The PTC located within an exon spanning more than 407bp 
    ##find which exon>407
    exon.long <- which(exon.length>407)
    #exclude the last exon(remove length(exon.length) from exon.long)
    exon.long <- exon.long[exon.long != length(exon.length)]
    #exclude the first exon
    exon.long <- exon.long[exon.long != 1]
    ##find if the PTC is in that long exon
    long.start = snv.control.exon$cds_start[exon.long]
    long.end = snv.control.exon$cds_end[exon.long]
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
      snv_long_re_1_flag[i] = 1
    }
    if (sum(NMD_status_long.2) > 0) {
      snv_long_re_2_flag[i] = 1
    }
  }, error = function(e) {
    # Handle error gracefully (e.g., log it, skip iteration)
    message(sprintf("Error at index %d: %s", i, e$message))
  })
}
snv_control_matched = snv_control_uniprot[which(snv_control_re_1_flag > 0 | snv_control_re_2_flag > 0)]
length(snv_control_matched)
snv_css_matched = snv_css_uniprot[which(snv_css_re_1_flag > 0 | snv_css_re_2_flag > 0)]
length(snv_css_matched)
snv_long_matched = snv_long_uniprot[which(snv_long_re_1_flag > 0 | snv_long_re_2_flag > 0)]
length(snv_long_matched)


##for framshift, use PTCs as boundary


