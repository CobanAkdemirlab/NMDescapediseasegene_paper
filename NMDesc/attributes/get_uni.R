#wild type uniport id

# Write unique UniProt IDs to a text file
writeLines(unique(na.omit(minus1_can_uni$uniprotswissprot)), "minus1_uniprot_ids.txt")
writeLines(unique(na.omit(minus1_control_uni$uniprotswissprot)), "minus1_control_uniprot_ids.txt")
writeLines(unique(na.omit(plus1_can_uni$uniprotswissprot)), "plus1_uniprot_ids.txt")
writeLines(unique(na.omit(plus1_control_uni$uniprotswissprot)), "plus1_control_uniprot_ids.txt")
writeLines(unique(na.omit(snv_can_uni$uniprotswissprot)), "snv_uniprot_ids.txt")
writeLines(unique(na.omit(snv_control_uni$uniprotswissprot)), "snv_control_uniprot_ids.txt")
