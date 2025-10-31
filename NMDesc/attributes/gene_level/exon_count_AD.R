
canonical_tx <- getBM(
  attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_is_canonical"),
  filters    = "hgnc_symbol",
  values     = unique(pli_AD$gene),
  mart       = ensembl
) %>%
  filter(transcript_is_canonical == 1) %>%
  distinct(hgnc_symbol, ensembl_transcript_id, .keep_all = TRUE)

cds_info <- getBM(
  attributes = c("ensembl_transcript_id", "cds_length"),
  filters    = "ensembl_transcript_id",
  values     = canonical_tx$ensembl_transcript_id,
  mart       = ensembl
)

cds_lengths <- canonical_tx %>%
  inner_join(cds_info, by = "ensembl_transcript_id") %>%
  transmute(gene = hgnc_symbol, cds_length = as.numeric(cds_length)) %>%
  group_by(gene) %>%
  slice_max(order_by = cds_length, n = 1, with_ties = FALSE) %>%
  ungroup()

pli_AD_cds <- pli_AD %>%
  left_join(cds_lengths, by = "gene")

# -------------------------
# 2) Exon counts for canonical transcripts
# -------------------------
exons <- getBM(
  attributes = c("ensembl_transcript_id", "cds_start", "cds_end", "rank", "strand"),
  filters    = "ensembl_transcript_id",
  values     = canonical_tx$ensembl_transcript_id,
  mart       = ensembl
)

exons_joined <- exons %>%
  inner_join(canonical_tx, by = "ensembl_transcript_id")

exon_counts <- exons_joined %>%
  group_by(ensembl_transcript_id, hgnc_symbol, strand) %>%
  summarise(exon_num = max(rank), .groups = "drop") %>%
  transmute(gene = hgnc_symbol, exon_num = as.numeric(exon_num)) %>%
  left_join(pli_AD %>% select(gene, category), by = "gene")

# -------------------------
# 3) Collapse categories to Case vs Control
# -------------------------
collapse_to_two_groups <- function(x) {
  lx <- tolower(as.character(x))
  grp <- dplyr::case_when(
    str_detect(lx, "control") ~ "Control",
    str_detect(lx, "minus1|plus1|snv|nonsense") ~ "Case",
    TRUE ~ NA_character_
  )
  factor(grp, levels = c("Case", "Control"))
}

# CDS: Case vs Control
pli_AD_cds2 <- pli_AD_cds %>%
  mutate(group2 = collapse_to_two_groups(category)) %>%
  filter(!is.na(group2), is.finite(cds_length), cds_length > 0)

# Exons: Case vs Control
exon_counts2 <- exon_counts %>%
  mutate(group2 = collapse_to_two_groups(category)) %>%
  filter(!is.na(group2), is.finite(exon_num), exon_num > 0)

# -------------------------
# 4) Plots with Wilcoxon p-values
# -------------------------
two_colors <- c(Case = "#4C78A8", Control = "#9EC3E6")
comparisons_two <- list(c("Case", "Control"))

# CDS length plot
cds_plot2 <- ggplot(pli_AD_cds2, aes(x = group2, y = cds_length, fill = group2)) +
  geom_boxplot(width = 0.22, color = "black", outlier.size = 0.9) +
  scale_y_log10() +
  scale_fill_manual(values = two_colors, guide = "none") +
  labs(
    title = "CDS Length in AD Genes",
    x = NULL, y = "CDS Length (bp)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  stat_compare_means(comparisons = comparisons_two, method = "wilcox.test", label = "p.format")

print(cds_plot2)
ggsave("cds_length_case_control.pdf", plot = cds_plot2, width = 6, height = 4)

# Exon count plot
exon_plot2 <- ggplot(exon_counts2, aes(x = group2, y = exon_num, fill = group2)) +
  geom_boxplot(width = 0.22, color = "black", outlier.size = 0.9) +
  scale_y_log10() +
  scale_fill_manual(values = two_colors, guide = "none") +
  labs(
    title = "Exon Count in AD Genes",
    x = NULL, y = "Number of Exons"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  stat_compare_means(comparisons = comparisons_two, method = "wilcox.test", label = "p.format")

print(exon_plot2)
ggsave("exon_count_case_control.pdf", plot = exon_plot2, width = 6, height = 4)







# --- Collapse to Case vs Control ---
collapse_to_two_groups <- function(x) {
  lx <- tolower(as.character(x))
  grp <- dplyr::case_when(
    str_detect(lx, "control") ~ "Control",
    str_detect(lx, "minus1|plus1|snv|nonsense") ~ "Case",
    TRUE ~ NA_character_
  )
  factor(grp, levels = c("Case", "Control"))
}

AD_region2 <- AD_region %>%
  mutate(group2 = collapse_to_two_groups(category)) %>%
  filter(!is.na(group2), is.finite(nmdesc), nmdesc > 0)

# --- Define colors & comparisons ---
two_colors <- c(Case = "#4C78A8", Control = "#9EC3E6")
comparisons_two <- list(c("Case", "Control"))

# --- Build plot ---
region_plot2 <- ggplot(AD_region2, aes(x = group2, y = nmdesc, fill = group2)) +
  geom_boxplot(width = 0.22, color = "black", outlier.size = 0.8) +
  theme_minimal() +
  scale_y_log10() +
  labs(
    title = "Canonical NMD Escape Region Length in AD Genes ",
    y = "NMDesc Length (bp)",
    x = NULL
  ) +
  theme(
    axis.text.x  = element_text(face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = two_colors, guide = "none") +
  stat_compare_means(comparisons = comparisons_two, method = "binomial.test", label = "p.format")

print(region_plot2)

# --- Save to file (no Cairo device) ---
ggsave("nmdesc_length_case_control.pdf", plot = region_plot2, width = 6, height = 4)

# --- Optional: print p-value to console ---
if (length(unique(AD_region2$group2)) == 2) {
  cat("\nWilcoxon (NMD escape length, Case vs Control):\n")
  print(wilcox.test(nmdesc ~ group2, data = AD_region2, exact = FALSE))
}


AD_unique_Tugce = AD_Tugce %>%
  distinct(Gene_name, .keep_all = TRUE)
write.csv(AD_unique_Tugce, "AD_unique_Tugce.csv", row.names = FALSE)


