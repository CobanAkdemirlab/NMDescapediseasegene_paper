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

#in pfam_all_pervariant_cds, filter for the emel's clean list
pfam_2ac = pfam_all_pervariant_cds[which(pfam_all_pervariant_cds$Variant_Key %in% WT_var_full_2ac$key),]
pfam_ad = pfam_all_pervariant_cds[which(pfam_all_pervariant_cds$Variant_Key %in% ppi_all_AD$Variant_Key),]
# Plot with p-values
ggplot(pfam_2ac, 
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



plot_df <- pfam_ad %>% filter(is.finite(distance_to_cds_end))

# 1) Pick the shown upper limit with a little headroom
y_max_data  <- max(plot_df$distance_to_cds_end, na.rm = TRUE)
y_upper     <- y_max_data * 1.15   # 15% headroom so text/brackets never touch the top

# 2) Put labels comfortably inside (e.g., 82–94% of the panel height)
n_comp      <- length(pairs_base)
label_y     <- seq(0.82, 0.94, length.out = n_comp) * y_upper

ggplot(plot_df, aes(x = group, y = distance_to_cds_end, fill = group)) +
  geom_violin(trim = TRUE, alpha = 0.7, na.rm = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black", na.rm = TRUE) +
  scale_fill_manual(values = group_colors2, guide = "none") +
  coord_cartesian(ylim = c(0, y_upper), clip = "on") +   # keep inside the border
  theme_classic(base_size = 14) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA),
    axis.text.x       = element_text(angle = 30, hjust = 1),
    plot.title        = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(title = "Distance from PTC to CDS end across variant groups",
       x = NULL, y = "Distance to CDS end (bp)") +
  ggpubr::stat_compare_means(
    comparisons      = pairs_base,
    method           = "wilcox.test",
    label            = "p.format",
    hide.ns          = FALSE,                 # show everything
    y.position       = label_y,               # bracket height
    label.y          = label_y - 0.01*y_upper,# nudge text below bracket
    step.increase    = 0,                     # don't auto-bump past the limit
    bracket.nudge.y  = 0,
    tip.length       = 0.01,
    vjust            = 0                      # push text downward a bit
  )

#do AD part, filter for ppi_all_AD$Variant_Key
#get canonical transcript for AD genes
AD_tr = getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_is_canonical"),
  filters    = "ensembl_gene_id",
  values     = unique(ppi_all_AD$Gene_ID),
  mart       = ensembl
)
AD_plot_df = plot_df %>% filter(Variant_Key %in% ppi_all_AD$Variant_Key)
ggplot(pfam_ad, 
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

