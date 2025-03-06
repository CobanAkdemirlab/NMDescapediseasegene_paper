#merge PTC_info with BM.info4
PTC_info2 = merge(PTC_info,BM.info5,by.x = 'transcript',by.y='ensembl_transcript_id')
PTC_info2$tr_region = PTC_info2$cds_length - PTC_info2$can_region - PTC_info2$css_region - PTC_info2$long_region
#if PTC_info2$transcript %in% first_long$ensembl_transcript_id, then plus 150 to PTC_info2$tr_region 
#PTC_info2$tr_region[PTC_info2$transcript %in% first_long$ensembl_transcript_id] = PTC_info2$tr_region[PTC_info2$transcript %in% first_long$ensembl_transcript_id] + 150

