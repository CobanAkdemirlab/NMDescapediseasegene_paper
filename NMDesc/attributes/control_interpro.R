

# Query InterPro domain locations
snv_control1_interpro <- getBM(
  attributes = c("ensembl_transcript_id", "interpro", "interpro_description", 
                 "interpro_start", "interpro_end"),
  filters = "ensembl_transcript_id",
  values = snv_control1_tr6$ensembl_transcript_id,  # TP53 Gene ID
  mart = ensembl
)

#select the rows that control1 NMDesc region overlap with interpro_start and interpro_end
##merge with snv_control1_tr6, by ensembl_transcript_id, keep all rows in snv_control1_interpro
snv_control1_interpro2 <- merge(snv_control1_interpro, snv_control1_tr6, by = "ensembl_transcript_id", all.y = TRUE)
#if end > aft.NMD.ind & start < aft.NMD.ind then match
snv_control1_interpro2$interpro_match = ifelse(snv_control1_interpro2$interpro_end >= snv_control1_interpro2$aft.NMD.ind/3 & snv_control1_interpro2$interpro_start <= (snv_control1_interpro2$aft.NMD.ind/3-50), 1, 0)
#keep matched rows
snv_control1_interpro3 <- snv_control1_interpro2[snv_control1_interpro2$interpro_match == 1,]
length(unique(snv_control1_interpro3$hgnc_symbol))
write.csv(snv_control1_interpro3, file = "snv_control1_interpro3.csv", row.names = FALSE)