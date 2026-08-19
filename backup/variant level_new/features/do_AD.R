#filter for AD genes
pli_all$gene

omim_AD_symbols <- omim2 %>%
  filter(str_detect(Inheritance_pattern, fixed("Autosomal dominant", ignore_case = TRUE))) %>%
  mutate(hgnc_symbol = str_trim(hgnc_symbol)) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct(hgnc_symbol) %>%
  pull(hgnc_symbol)

pli_AD <- pli_all %>%
  filter(gene %in% omim_AD_symbols)

write.csv(pli_AD, file = "pli_AD_genes.csv", row.names = FALSE)

#prepare data for venn plot
table(pli_AD$category)
cat_map <- c(
  "SNV" = "Nonsense", "snv" = "Nonsense", "Nonsense" = "Nonsense",
  "plus1" = "Plus1", "Plus1" = "Plus1",
  "minus1" = "Minus1", "Minus1" = "Minus1"
)

pli_AD2 <- pli_AD %>%
  mutate(category_norm = recode(category, !!!cat_map)) %>%
  filter(!is.na(category_norm)) %>%
  select(gene, category_norm) %>%
  distinct()  # unique (gene, category) pairs

mem <- pli_AD2 %>%
  mutate(val = 1) %>%
  pivot_wider(
    id_cols = gene,
    names_from = category_norm,
    values_from = val,
    values_fill = 0
  ) %>%
  # Ensure all three columns exist
  mutate(
    Nonsense = if (!"Nonsense" %in% names(.)) 0 else Nonsense,
    Plus1    = if (!"Plus1"    %in% names(.)) 0 else Plus1,
    Minus1   = if (!"Minus1"   %in% names(.)) 0 else Minus1
  ) %>%
  mutate(across(c(Nonsense, Plus1, Minus1), ~ . > 0))

N <- mem$Nonsense
P <- mem$Plus1
M <- mem$Minus1

venn_counts <- c(
  `Nonsense only`                 = sum( N & !P & !M ),
  `Plus1 only`                    = sum(!N &  P & !M ),
  `Minus1 only`                   = sum(!N & !P &  M ),
  `Nonsense ∩ Plus1`              = sum( N &  P & !M ),
  `Nonsense ∩ Minus1`             = sum( N & !P &  M ),
  `Plus1 ∩ Minus1`                = sum(!N &  P &  M ),
  `Nonsense ∩ Plus1 ∩ Minus1`     = sum( N &  P &  M )
)

total_genes <- sum(venn_counts)
venn_df <- tibble(
  region = names(venn_counts),
  count  = as.integer(venn_counts),
  pct    = round(100 * count / ifelse(total_genes == 0, 1, total_genes), 1)
)

#create new pfam and ppi data
pfam_all_AD = pfam_all_pervariant[which(pfam_all_pervariant$hgnc_symbol %in% pli_AD$gene),]
pfam_data_AD <- pfam_all_pervariant %>%
  count(group, in_pfam) %>%
  group_by(group) %>%
  mutate(prop = n / sum(n),
         pct  = percent(prop, accuracy = 0.1)) %>%
  ungroup()
pfam_data_AD$group_label = pfam_data_AD$group
pfam_data_AD$matched = pfam_data_AD$in_pfam
pfam_data_AD$total = pfam_data_AD$n
write.csv(pfam_data_AD, file = "pfam_data.csv", row.names = FALSE)

#create ppi_all$hgnc_symbol by  pfam_all_pervariant$hgnc_symbol and pfam_all_pervariant$Variant_Key
ppi_all$hgnc_symbol = pfam_all_pervariant$hgnc_symbol[match(ppi_all$Variant_Key, pfam_all_pervariant$Variant_Key)]

ppi_all_AD = ppi_all[which(ppi_all$hgnc_symbol %in% pli_AD$gene),]
ppi_df_AD <- ppi_all_AD %>%
  count(group, in_interface) %>%
  group_by(group) %>%
  tidyr::complete(in_interface = c(FALSE, TRUE), fill = list(n = 0)) %>%
  mutate(prop = n / sum(n),
         pct  = scales::percent(prop, accuracy = 0.1)) %>%
  ungroup()
#bulid these : group_label, overlap_count, total_count, percent
ppi_df_AD$group_label = ppi_df_AD$group
ppi_df_AD$matched = ppi_df_AD$in_interface
ppi_df_AD$total_count = ppi_df_AD$n
ppi_df_AD$percent = ppi_df_AD$pct
ppi_df_AD$overlap_count = ppi_df_AD$n
ppi_df_AD$condition = 

write.csv(ppi_df_AD, file = "ppi_data.csv", row.names = FALSE)

#idr

