#input data
plus1_can_gene0217 = read.table("~/Desktop/frameshift/plus1_can_gene0217.txt", header = T, sep = "\t")
plus2_can_gene0217 = read.table("~/Desktop/frameshift/plus2_can_gene0217.txt", header = T, sep = "\t")
snv_plp_ptc_p1120_NMDenriched2_can <- read_csv("~/Desktop/new_clinvar/snv_list/list4/snv_plp_ptc_p1120_NMDenriched2_can.txt", 
                                               col_names = FALSE)


#get canonical transcript
plus1_can_list = getBM(
  attributes = c("ensembl_transcript_id", 'transcript_is_canonical'),
  filters = "hgnc_symbol",
  values = plus1_can_gene0217$hgnc_symbol,
  mart = ensembl
)
plus1_can_list = plus1_can_list$ensembl_transcript_id[which(plus1_can_list$transcript_is_canonical == 1)]
plus2_can_list = getBM(
  attributes = c("ensembl_transcript_id", 'transcript_is_canonical'),
  filters = "hgnc_symbol",
  values = plus2_can_gene0217$hgnc_symbol,
  mart = ensembl
)
plus2_can_list = plus2_can_list$ensembl_transcript_id[which(plus2_can_list$transcript_is_canonical == 1)]
snv_can_list = getBM(
  attributes = c("ensembl_transcript_id", 'transcript_is_canonical'),
  filters = "hgnc_symbol",
  values = snv_plp_ptc_p1120_NMDenriched2_can$X1,
  mart = ensembl
)
snv_can_list = snv_can_list$ensembl_transcript_id[which(snv_can_list$transcript_is_canonical == 1)]

plus1_can_df = merged.2[which(merged.2$tr_sim %in% plus1_can_list),]
plus2_can_df = merged.2[which(merged.2$tr_sim %in% plus2_can_list),]
snv_can_df = merged.2[which(merged.2$tr_sim %in% snv_can_list),]

#exon.count
plus1_can_df$type <- 'plus1_can_df'
plus2_can_df$type <- 'plus2_can_df'
snv_can_df$type <- 'snv_can_df'
combined_df <- bind_rows(
  plus1_can_df, 
  plus2_can_df, 
  snv_can_df
)
combined_df <- combined_df %>%
  mutate(
    group = case_when(
      grepl("can", type) ~ "can",
      type == "control_df" ~ "control"
    )
  )
combined_df$type <- factor(
  combined_df$type,
  levels = combined_df %>%
    arrange(group, type) %>%
    pull(type) %>% 
    unique()
)
ggplot(combined_df, aes(x=type,y = exon_count, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Comparison of Exon Count",
    x = "Type",
    y = "exon_count"
  ) +
  theme_minimal() +
  #make y log scale
  scale_y_log10() +
  scale_fill_brewer(palette = "Set2", name = "") +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")


#threeUTR_AU_content
ggplot(combined_df, aes(x=type, y = threeUTR_AU_content, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Comparison of threeUTR_AU_content",
    x = "Exon Count",
    y = "Frequency"
  ) +
  theme_minimal() +
  #make x log scale
  #scale_x_log10() +
  scale_fill_brewer(palette = "Set2", name = "") +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")

#threeUTR_GC_content
ggplot(combined_df, aes(x=type, y = threeUTR_GC_content, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Comparison of threeUTR_GC_content",
    x = "Exon Count",
    y = "Frequency"
  ) +
  theme_minimal() +
  #make x log scale
  #scale_x_log10() +
  scale_fill_brewer(palette = "Set2", name = "") +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")

#threeUTR_UC_content 
ggplot(combined_df, aes(x=type, y = threeUTR_UC_content, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Comparison of threeUTR_UC_content",
    x = "Exon Count",
    y = "Frequency"
  ) +
  theme_minimal() +
  #make x log scale
  #scale_x_log10() +
  scale_fill_brewer(palette = "Set2", name = "") +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")

#threeUTR_length 
ggplot(combined_df, aes(x=type, y = threeUTR_length, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Comparison of threeUTR_length",
    x = "Exon Count",
    y = "Frequency"
  ) +
  theme_minimal() +
  #make y log scale
  scale_y_log10() +
  scale_fill_brewer(palette = "Set2", name = "") +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")

#ThreeUTR_AUcontentlast200
ggplot(combined_df, aes(x=type, y = ThreeUTR_AUcontentlast200, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Comparison of ThreeUTR_AUcontentlast200",
    x = "Exon Count",
    y = "Frequency"
  ) +
  theme_minimal() +
  #make x log scale
  #scale_x_log10() +
  scale_fill_brewer(palette = "Set2", name = "") +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")


#ThreeUTR_AUcontentfirst200
ggplot(combined_df, aes(x=type, y = ThreeUTR_AUcontentfirst200, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Comparison of ThreeUTR_AUcontentfirst200",
    x = "Exon Count",
    y = "Frequency"
  ) +
  theme_minimal() +
  #make x log scale
  #scale_x_log10() +
  scale_fill_brewer(palette = "Set2", name = "") +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")

#ThreeUTR_GCcontentlast200
ggplot(combined_df, aes(x=type, y = ThreeUTR_GCcontentlast200, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Comparison of ThreeUTR_GCcontentlast200",
    x = "Exon Count",
    y = "Frequency"
  ) +
  theme_minimal() +
  #make x log scale
  #scale_x_log10() +
  scale_fill_brewer(palette = "Set2", name = "") +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")


#ThreeUTR_GCcontentfirst200
ggplot(combined_df, aes(x=type, y = ThreeUTR_GCcontentfirst200, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Comparison of ThreeUTR_GCcontentfirst200",
    x = "Exon Count",
    y = "Frequency"
  ) +
  theme_minimal() +
  #make x log scale
  #scale_x_log10() +
  scale_fill_brewer(palette = "Set2", name = "") +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")

#ThreeUTR_UCcontentlast200
ggplot(combined_df, aes(x=type, y = ThreeUTR_UCcontentlast200, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Comparison of ThreeUTR_UCcontentlast200",
    x = "Exon Count",
    y = "Frequency"
  ) +
  theme_minimal() +
  #make x log scale
  #scale_x_log10() +
  scale_fill_brewer(palette = "Set2", name = "") +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")

#ThreeUTR_UCcontentfirst200
ggplot(combined_df, aes(x=type, y = ThreeUTR_UCcontentfirst200, fill = group)) +
  geom_boxplot() +
  labs(
    title = "Comparison of ThreeUTR_UCcontentfirst200",
    x = "Exon Count",
    y = "Frequency"
  ) +
  theme_minimal() +
  #make x log scale
  #scale_x_log10() +
  scale_fill_brewer(palette = "Set2", name = "") +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")

#fiveutrseqs.uORF
binom.test(table(plus1_can_df$fiveutrseqs.uORF)[2],nrow(plus1_can_df),table(control_df$fiveutrseqs.uORF)[2]/nrow(control_df))
binom.test(table(plus1_css_df$fiveutrseqs.uORF)[2],nrow(plus1_css_df),table(control_df$fiveutrseqs.uORF)[2]/nrow(control_df))
binom.test(table(plus1_long_df$fiveutrseqs.uORF)[2],nrow(plus1_long_df),table(control_df$fiveutrseqs.uORF)[2]/nrow(control_df))

#threeUTR.introns 
#found sig p
binom.test(table(plus1_can_df$threeUTR.introns)[2],nrow(plus1_can_df),table(control_df$threeUTR.introns)[2]/nrow(control_df))
binom.test(table(plus1_css_df$threeUTR.introns)[2],nrow(plus1_css_df),table(control_df$threeUTR.introns)[2]/nrow(control_df))
binom.test(table(plus1_long_df$threeUTR.introns)[2],nrow(plus1_long_df),table(control_df$threeUTR.introns)[2]/nrow(control_df))
binom.test(table(plus2_can_df$threeUTR.introns)[2],nrow(plus2_can_df),table(control_df$threeUTR.introns)[2]/nrow(control_df))
binom.test(table(plus2_css_df$threeUTR.introns)[2],nrow(plus2_css_df),table(control_df$threeUTR.introns)[2]/nrow(control_df))
binom.test(table(plus2_long_df$threeUTR.introns)[2],nrow(plus2_long_df),table(control_df$threeUTR.introns)[2]/nrow(control_df))

#threeUTR.PTBP1
binom.test(table(plus1_can_df$threeUTR.PTBP1)[2],nrow(plus1_can_df),table(control_df$threeUTR.PTBP1)[2]/nrow(control_df))
binom.test(table(plus1_css_df$threeUTR.PTBP1)[2],nrow(plus1_css_df),table(control_df$threeUTR.PTBP1)[2]/nrow(control_df))
binom.test(table(plus1_long_df$threeUTR.PTBP1)[2],nrow(plus1_long_df),table(control_df$threeUTR.PTBP1)[2]/nrow(control_df))

#cds.nucloc 
binom.test(sum(plus1_can_df$cds.nucloc != 'NA'),nrow(plus1_can_df),sum(control_df$cds.nucloc != 'NA')/nrow(control_df))
binom.test(sum(plus1_css_df$cds.nucloc != 'NA'),nrow(plus1_css_df),sum(control_df$cds.nucloc != 'NA')/nrow(control_df))
binom.test(sum(plus1_long_df$cds.nucloc != 'NA'),nrow(plus1_long_df),sum(control_df$cds.nucloc != 'NA')/nrow(control_df))
