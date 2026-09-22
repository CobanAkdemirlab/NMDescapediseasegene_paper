
#number of submitters on variant level
submit_info$NumberSubmitters
submit_info$tx_id2
submit_info$key

#variant list
variants_all5_1_0621$Variant_Key
variants_all5_1_0621$NumberSubmitters = submit_info$NumberSubmitters[match(variants_all5_1_0621$Variant_Key, submit_info$key)]
variants_all5_1_0621_clinvar  = variants_all5_1_0621[which(variants_all5_1_0621$group %in% c('snv_disease','fs_disease')),]
sum(is.na(variants_all5_1_0621_clinvar$NumberSubmitters))
sum(variants_all5_1_0621$NumberSubmitters > 1, na.rm=T)
#see if number of submitters is different in groups
length(unique(submit_info$tx_id2))