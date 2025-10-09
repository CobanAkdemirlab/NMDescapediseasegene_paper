#ptc-cds end distance

#get cds_end of pfam_all_pervariant$transcript from biomart
cds_end <- getBM(
  attributes = c("ensembl_transcript_id", "cds_end"),
  filters    = "ensembl_transcript_id",
  values     = unique(pfam_all_pervariant$transcript),
  mart       = ensembl
)
#group by transcript_id and get max cds_end
cds_end <- cds_end %>%
  group_by(ensembl_transcript_id) %>%
  summarise(cds_end = max(cds_end, na.rm = TRUE))
#add cds_end to pfam_all_pervariant
pfam_all_pervariant_cds = merge(pfam_all_pervariant,cds_end,by.x='transcript',by.y='ensembl_transcript_id',all.x=TRUE)

#plot distance to cds end
pfam_all_pervariant_cds$distance_to_cds_end = pfam_all_pervariant_cds$cds_end - pfam_all_pervariant_cds$ptc_pos
#only keep those with distance_to_cds_end > 0
pfam_all_pervariant_cds = pfam_all_pervariant_cds[pfam_all_pervariant_cds$distance_to_cds_end > 0,]

group_colors2 <- c(
  "Minus1" = "#1f77b4",
  "Minus1_Control" = "#aec7e8",
  "Plus1" = "#ff7f0e",
  "Plus1_Control" = "#ffbb78",
  "Nonsense" = "#2ca02c",
  "Nonsense_Control" = "#98df8a"
)
# violin + boxplot overlay
ggplot(pfam_all_pervariant_cds, 
                 aes(x = group, y = distance_to_cds_end, fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
  scale_fill_manual(values = group_colors2, guide = "none") +
  labs(
    title = "Distance from PTC to CDS end across variant groups",
    x = NULL,
    y = "Distance to CDS end (bp)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  )


library(ggplot2)
library(ggpubr)

# Make sure groups match your color keys
pfam_all_pervariant_cds$group <- factor(
  pfam_all_pervariant_cds$group,
  levels = c("Minus1","Minus1_Control",
             "Plus1","Plus1_Control",
             "Nonsense","Nonsense_Control")
)

# Case-control comparisons
pairs_base <- list(
  c("Minus1", "Minus1_Control"),
  c("Plus1", "Plus1_Control"),
  c("Nonsense", "Nonsense_Control")
)

# Plot with p-values
ggplot(pfam_all_pervariant_cds, 
       aes(x = group, y = distance_to_cds_end, fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  scale_y_continuous(limits = c(0, 1000), expand = expansion(mult = c(0.02, 0.08)))+
  geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
  scale_fill_manual(values = group_colors2, guide = "none") +
  labs(
    title = "Distance from PTC to CDS end across variant groups",
    x = NULL,
    y = "Distance to CDS end (bp)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title  = element_text(hjust = 0.5, face = "bold")
  ) +
  #scale_y_log10() +  # Log scale for y-axis
  ggpubr::stat_compare_means(
    comparisons = pairs_base,
    method      = "wilcox.test",
    label       = "p.format",   # change to "p.signif" for stars (*, **, ***)
    hide.ns     = TRUE,
    label.y     = c(
      max(pfam_all_pervariant_cds$distance_to_cds_end, na.rm = TRUE) * 1.05,
      max(pfam_all_pervariant_cds$distance_to_cds_end, na.rm = TRUE) * 1.15,
      max(pfam_all_pervariant_cds$distance_to_cds_end, na.rm = TRUE) * 1.25
    )
  )
