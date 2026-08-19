library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)

#1. input gene lists
minus1_can_gene0412 <- read_csv("/Users/jxu14/Desktop/frameshift/plus1_can_gene0412.txt")
plus1_can_gene0412 <- read_csv("/Users/jxu14/Desktop/frameshift/plus2_can_gene0412.txt")
snv_nmd_p0411_NMDenriched_snv_can <- read_csv("~/Desktop/new_clinvar/snv_nmd_p0411_NMDenriched_snv_can.txt")
snv_nmd_p0411_NMDenriched_snv_can$gene = snv_nmd_p0411_NMDenriched_snv_can$x
minus1_gnomAD_control_genes <- read_csv("~/Downloads/plus1_gnomAD_control_genes.csv",col_names = 'gene',skip=1)
plus1_gnomAD_control_genes <- read_csv("~/Downloads/plus2_gnomAD_control_genes.csv",col_names = 'gene',skip=1)
snv_gnomAD_control_genes <- read_csv("~/Downloads/snv_gnomAD_control_genes.csv",col_names = 'gene',skip=1)

#create can gene set
can_all = unique(c(plus1_can_gene0412$gene, minus1_can_gene0412$gene,
                   snv_nmd_p0411_NMDenriched_snv_can$x))
#remove can_all from plus1_gnomAD_control_genes
control_minus1 = data.frame(gene = unique(minus1_gnomAD_control_genes$gene[!minus1_gnomAD_control_genes$gene %in% can_all]))
control_plus1 = data.frame(gene = unique(plus1_gnomAD_control_genes$gene[!plus1_gnomAD_control_genes$gene %in% can_all]))
control_snv = data.frame(gene = unique(snv_gnomAD_control_genes$gene[!snv_gnomAD_control_genes$gene %in% can_all]))


# Load gnomAD constraint metrics
Lof_metrics <- read.delim("/Users/jxu14/Desktop/autism/data/gnomad.v2.1.1.lof_metrics.by_gene.txt")
# Merge with constraint metrics
minus1_pli <- merge(minus1_can_gene0412, Lof_metrics, by = 'gene') %>%
  dplyr::select(gene, pLI, oe_lof_upper) %>%
  mutate(category = "minus1")

plus1_pli <- merge(plus1_can_gene0412, Lof_metrics, by = 'gene') %>%
  dplyr::select(gene, pLI, oe_lof_upper) %>%
  mutate(category = "plus1")

snv_pli <- merge(snv_nmd_p0411_NMDenriched_snv_can, Lof_metrics, by = 'gene') %>%
  dplyr::select(gene, pLI, oe_lof_upper) %>%
  mutate(category = "SNV")

minus1_control_pli = merge(control_minus1, Lof_metrics, by = 'gene') %>%
  dplyr::select(gene, pLI, oe_lof_upper) %>%
  mutate(category = "minus1_Control")

plus1_control_pli = merge(control_plus1, Lof_metrics, by = 'gene') %>%
  dplyr::select(gene, pLI, oe_lof_upper) %>%
  mutate(category = "plus1_Control")

snv_control_pli = merge(control_snv, Lof_metrics, by = 'gene') %>%
  dplyr::select(gene, pLI, oe_lof_upper) %>%
  mutate(category = "SNV_Control")

# Combine all lists
pli_all <- bind_rows(
  minus1_pli,
  plus1_pli,
  snv_pli,
  minus1_control_pli,
  plus1_control_pli,
  snv_control_pli
)
#use pli_AD instead

#2. define ggplot paras
# Define comparisons
comparisons <- list(
  c("minus1", "minus1_Control"),
  c("plus1", "plus1_Control"),
  c("SNV", "SNV_Control")
)

# Custom colors
group_colors <- c(
  "minus1" = "#1f77b4",
  "minus1_Control" = "#aec7e8",  # lighter blue
  "plus1" = "#ff7f0e",
  "plus1_Control" = "#ffbb78",  # lighter orange
  "SNV" = "#2ca02c",
  "SNV_Control" = "#98df8a"     # lighter green
)

#3. plot pli and loeuf
# Plot pLI
# Create the plot object
pli_plot <- ggplot(pli_all, aes(x = category, y = pLI, fill = category)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  labs(title = "pLI distribution by gene group",
       y = "pLI",
       x = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")
ggplot(pli_AD, aes(x = category, y = pLI, fill = category)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  labs(title = "pLI distribution in AD genes",
       y = "pLI",
       x = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")
ggsave("pli_distribution_by_group.pdf", plot = pli_plot, width = 8, height = 5)

loeuf_plot <- ggplot(pli_all, aes(x = category, y = oe_lof_upper, fill = category)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  labs(title = "LOEUF distribution by gene group",
       y = "LOEUF",
       x = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")
loeuf_plot <- ggplot(pli_AD, aes(x = category, y = oe_lof_upper, fill = category)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  labs(title = "LOEUF distribution in AD genes",
       y = "LOEUF",
       x = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")

ggsave("loeuf_distribution.pdf", plot = loeuf_plot,width = 8, height = 5)

#4. plot cds_length
#use canonical transcript to get cds length
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
canonical_tx <- getBM(
  attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_is_canonical"),
  filters = "hgnc_symbol",
  values = unique(pli_all$gene),
  mart = ensembl
) %>%
  filter(transcript_is_canonical == 1)
cds_info <- getBM(
  attributes = c("ensembl_transcript_id", "cds_length"),
  filters = "ensembl_transcript_id",
  values = canonical_tx$ensembl_transcript_id,
  mart = ensembl
)
cds_lengths <- canonical_tx %>%
  inner_join(cds_info, by = "ensembl_transcript_id") %>%
  dplyr::select(gene = hgnc_symbol, cds_length)

pli_all_cds <- pli_all %>%
  left_join(cds_lengths, by = "gene")
pli_AD_cds <- pli_AD %>%
  left_join(cds_lengths, by = "gene")
cds_plot = ggplot(pli_AD_cds, aes(x = category, y = cds_length, fill = category)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  labs(title = "CDS Length Distribution in AD genes",
       y = "CDS Length (bp)",
       x = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_y_log10() +
  scale_fill_manual(values = group_colors) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")
ggsave("cds_length_distribution.pdf", plot = cds_plot, width = 8, height = 5)

#5. plot exon count
exons <- getBM(
  attributes = c("ensembl_transcript_id", "cds_start", "cds_end", "rank", "strand"),
  filters = "ensembl_transcript_id",
  values = canonical_tx$ensembl_transcript_id,
  mart = ensembl
)
# add gene name to exon info
exons_joined <- exons %>%
  inner_join(canonical_tx, by = "ensembl_transcript_id")
exon_counts <- exons_joined %>%
  group_by(ensembl_transcript_id, hgnc_symbol, strand) %>%
  arrange(rank) %>%
  summarise(
    exon_num = max(rank),
    .groups = "drop"
  ) %>%
  mutate(gene  = hgnc_symbol) %>%
  dplyr::select(gene,exon_num)
#get category info from pli_all
exon_counts <- exon_counts %>%
  left_join(pli_AD %>% dplyr::select(gene, category), by = "gene")

exon_plot = ggplot(exon_counts, aes(x = category, y = exon_num, fill = category)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  scale_y_log10() +
  labs(title = "Exon count in AD Genes",
       y = "number",
       x = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")
library(dplyr)
library(ggplot2)
library(ggpubr)

# 1) Rename groups to Title Case (same mapper you used before)
rename_levels <- function(x) {
  recoded <- dplyr::recode(tolower(x),
                           "minus1" = "Minus1",
                           "minus1_control" = "Minus1_Control",
                           "plus1" = "Plus1",
                           "plus1_control" = "Plus1_Control",
                           "snv" = "Nonsense",
                           "snv_control" = "Nonsense_Control",
                           .default = x
  )
  factor(recoded, levels = c("Minus1","Minus1_Control",
                             "Plus1","Plus1_Control",
                             "Nonsense","Nonsense_Control"))
}

exon_counts2 <- exon_counts %>%
  mutate(category = rename_levels(category)) %>%
  filter(!is.na(category)) %>%
  droplevels()

# 2) Rebuild comparisons to match new labels and available groups
# If you already have `comparisons` in old labels, convert them:
comparisons_new <- lapply(comparisons, function(p) as.character(rename_levels(p)))

present_groups <- levels(exon_counts2$category)
comparisons_use <- Filter(function(p) length(p) == 2 && all(p %in% present_groups), comparisons_new)

# 3) Plot (with bold/centered title and bold x-label & ticks)
group_colors2 <- c(
  Minus1="#1f77b4", Minus1_Control="#aec7e8",
  Plus1 ="#ff7f0e", Plus1_Control ="#ffbb78",
  Nonsense="#2ca02c", Nonsense_Control="#98df8a"
)

p_exon <- ggplot(exon_counts2, aes(x = category, y = exon_num, fill = category)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  scale_y_log10() +
  scale_fill_manual(values = group_colors2, guide = "none") +
  labs(title = "Exon count by Gene Group",
       y = "Number",
       x = "Gene category") +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  ggpubr::stat_compare_means(comparisons = comparisons_use, method = "wilcox.test")

print(p_exon)

ggsave("exon_count_distribution.pdf", plot = exon_plot, width = 8, height = 5)


#6. plot NMDesc region length
#get snv transcript names
snv_tx <- getBM(
  attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_is_canonical"),
  filters = "hgnc_symbol",
  values = c(snv_pli$gene, snv_control_pli$gene),
  mart = ensembl
) %>%
  filter(transcript_is_canonical == 1)
snv_exon_info = exons_joined %>%
  filter(ensembl_transcript_id %in% snv_tx$ensembl_transcript_id) %>%
  dplyr::select(ensembl_transcript_id, cds_start, cds_end, rank, strand, hgnc_symbol)
snv_nmdesc_lengths <- snv_exon_info %>%
  group_by(ensembl_transcript_id, hgnc_symbol, strand) %>%
  arrange(rank) %>%
  summarise(
    last_exon_length = abs(cds_end[which.max(rank)] - cds_start[which.max(rank)]) + 1,
    penult_exon_end = cds_end[rank == max(rank) - 1],
    penult_exon_start = cds_start[rank == max(rank) - 1],
    exon_num = max(rank),
    .groups = "drop"
  ) %>%
  mutate(
    penult_exon_length = abs(penult_exon_end - penult_exon_start) + 1,
    penult_50bp = pmin(50, penult_exon_length),
    nmdesc_length = last_exon_length + penult_50bp
  ) %>%
  dplyr::select(gene = hgnc_symbol, nmdesc_length,exon_num,ensembl_transcript_id)
#add category info
snv_nmdesc_lengths <- snv_nmdesc_lengths %>%
  left_join(pli_all %>% dplyr::select(gene, category), by = "gene") %>%
  #select only category is SNV or SNV_Control
  filter(category %in% c("SNV", "SNV_Control"))

#get fs NMDesc length from PTC_info
minus1_tx <- getBM(
  attributes = c("external_gene_name", "ensembl_transcript_id", "transcript_is_canonical"),
  filters = "external_gene_name",
  values = unique(minus1_pli$gene),
  mart = ensembl
) %>%
  filter(transcript_is_canonical == 1)

minus1.cds.info <- getBM(
  attributes = c("ensembl_transcript_id", "coding"),
  filters = "ensembl_transcript_id",
  values = minus1_tx$ensembl_transcript_id,
  mart = ensembl
)

can_region_minus1 <- data.frame(
  gene = character(),
  ensembl_transcript_id = character(),
  nmdesc_minus1 = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(minus1_tx))) {
  gene <- minus1_tx$external_gene_name[i]
  ttranscript <- minus1_tx$ensembl_transcript_id[i]
  
  region_length <- PTC_info %>%
    filter(transcript == ttranscript, type == "minus1") %>%
    pull(can_region)
  
  if (length(region_length) > 0) {
    can_region_minus1 <- rbind(
      can_region_minus1,
      data.frame(
        gene = gene,
        ensembl_transcript_id = ttranscript,
        nmdesc_minus1 = region_length
      )
    )
  }
}

minus1_pli_with_region <- minus1_pli %>%
  left_join(can_region_minus1, by = "gene")
# Use biomart to get canonical transcript IDs for plus1 genes
plus1_tx <- getBM(
  attributes = c("external_gene_name", "ensembl_transcript_id", "transcript_is_canonical"),
  filters = "external_gene_name",
  values = unique(plus1_pli$gene),
  mart = ensembl
) %>%
  filter(transcript_is_canonical == 1)

can_region_plus1 <- data.frame(
  gene = character(),
  ensembl_transcript_id = character(),
  nmdesc_plus1 = numeric(),
  stringsAsFactors = FALSE
)

plus1.cds.info <- getBM(
  attributes = c("ensembl_transcript_id", "coding"),
  filters = "ensembl_transcript_id",
  values = plus1_tx$ensembl_transcript_id,
  mart = ensembl
)

for (i in seq_len(nrow(plus1_tx))) {
  gene <- plus1_tx$external_gene_name[i]
  ttranscript <- plus1_tx$ensembl_transcript_id[i]
  
  region_length <- PTC_info %>%
    filter(transcript == ttranscript, type == "plus1") %>%
    pull(can_region)
  
  if (length(region_length) > 0) {
    can_region_plus1 <- rbind(
      can_region_plus1,
      data.frame(
        gene = gene,
        ensembl_transcript_id = ttranscript,
        nmdesc_plus1 = region_length
      )
    )
  }
}

plus1_pli_with_region <- plus1_pli %>%
  left_join(can_region_plus1, by = "gene")
fs_nmdesc <- bind_rows(
  can_region_minus1 %>% mutate(category = "minus1") %>% rename(nmdesc = nmdesc_minus1),
  can_region_plus1 %>% mutate(category = "plus1") %>% rename(nmdesc = nmdesc_plus1)
)
#use self-define function to get NMDesc length info for fs control
control_genes <- unique(plus1_control_pli$gene)
# Get canonical transcripts
control_tx <- getBM(
  attributes = c("external_gene_name", "ensembl_transcript_id", "transcript_is_canonical"),
  filters = "external_gene_name",
  values = control_genes,
  mart = ensembl
) %>% filter(transcript_is_canonical == 1)

# Get CDS sequence for canonical transcripts
cds_info <- getBM(
  attributes = c("ensembl_transcript_id", "coding"),
  filters = "ensembl_transcript_id",
  values = control_tx$ensembl_transcript_id,
  mart = ensembl
)

# Get cds info
exon_info <- getBM(
  attributes = c("ensembl_transcript_id", "rank", "cds_start", "cds_end"),
  filters = "ensembl_transcript_id",
  values = control_tx$ensembl_transcript_id,
  mart = ensembl
)

# Join with gene symbols
control_tx_map <- control_tx %>%
  dplyr::select(gene = external_gene_name, ensembl_transcript_id)

# Prepare final results
can_region_control <- data.frame(
  gene = character(),
  ensembl_transcript_id = character(),
  nmdesc = numeric(),
  category = character()
)


# Loop through transcripts and compute can_region
for (i in seq_len(nrow(control_tx))) {
  gene <- control_tx$external_gene_name[i]
  transcript <- control_tx$ensembl_transcript_id[i]
  
  # Get CDS sequence
  seq_entry <- cds_info %>%
    filter(ensembl_transcript_id == transcript)
  
  # Get coding exon
  df <- exon_info %>%
    filter(ensembl_transcript_id == transcript) %>%
    filter(!is.na(cds_start) & !is.na(cds_end)) %>%
    arrange(rank)
  
  if (nrow(df) < 2) next  # Skip if fewer than 2 coding exons
  
  input_seq <- seq_entry$coding
  
  # Use tryCatch to skip errors gracefully
  tryCatch({
    PTC_ind <- get_PTC_plus_ind(input_seq, type = "plus1")
    PTC_status <- get_NMDesc_PTC(input_seq, bp = 50, transcript, type = "plus1")
    region_set <- get_NMDesc_PTC_region(PTC_ind, transcript, PTC_status)
    
    can_region_control <- rbind(
      can_region_control,
      data.frame(
        gene = gene,
        ensembl_transcript_id = transcript,
        nmdesc = region_set$can_region,
        category = "plus1_Control"
      )
    )
  }, error = function(e) {
    message(sprintf("Skipped %s due to error: %s", transcript, e$message))
    # Just skip on error
  })
}

# Prepare results
can_region_minus1_control <- data.frame(
  gene = character(),
  ensembl_transcript_id = character(),
  nmdesc = numeric(),
  category = character()
)

for (i in seq_len(nrow(control_tx))) {
  gene <- control_tx$external_gene_name[i]
  transcript <- control_tx$ensembl_transcript_id[i]
  
  seq_entry <- cds_info %>%
    filter(ensembl_transcript_id == transcript)
  
  df <- exon_info %>%
    filter(ensembl_transcript_id == transcript) %>%
    filter(!is.na(cds_start) & !is.na(cds_end)) %>%
    arrange(rank)
  
  if (nrow(df) < 2) next
  
  input_seq <- seq_entry$coding
  
  tryCatch({
    PTC_ind <- get_PTC_plus_ind(input_seq, type = "minus1")
    PTC_status <- get_NMDesc_PTC(input_seq, bp = 50, transcript, type = "minus1")
    region_set <- get_NMDesc_PTC_region(PTC_ind, transcript, PTC_status)
    
    can_region_minus1_control <- rbind(
      can_region_minus1_control,
      data.frame(
        gene = gene,
        ensembl_transcript_id = transcript,
        nmdesc = region_set$can_region,
        category = "minus1_Control"
      )
    )
  }, error = function(e) {
    message(sprintf("Skipped %s due to error: %s", transcript, e$message))
  })
}

fs_control_nmdesc <- bind_rows(
  can_region_minus1_control %>% mutate(category = "minus1_Control"),
  can_region_control %>% mutate(category = "plus1_Control") 
)
# Combine all NMDesc regions
region_combined <- bind_rows(
  snv_nmdesc_lengths %>% rename(nmdesc = nmdesc_length),
  fs_nmdesc,
  fs_control_nmdesc
)
AD_region = region_combined %>%
  filter(gene %in% pli_AD$gene)
#plot result
region_plot = ggplot(AD_region, aes(x = category, y = nmdesc, fill = category)) +
  geom_boxplot(width = 0.1, color = "black") +
  theme_minimal() +
  scale_y_log10() +
  labs(title = "Canonical NMD Escape Region Length in AD genes",
       y = "NMDesc Length (bp)",
       x = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = group_colors) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")
ggsave("nmdesc_length_distribution.pdf", plot = region_plot, width = 8, height = 5)

#7. plot NMDesc region length/cds length ratio
length_ratio_df <- pli_all %>% #this is gene level
  left_join(region_combined, by = "gene") %>%
  left_join(cds_lengths, by = "gene") %>%
  filter(!is.na(nmdesc), !is.na(cds_length), cds_length > 0) %>%
  mutate(
    ratio = nmdesc / cds_length
  )

length_ratio_df <- pli_all_with_tx %>%
  left_join(region_combined, by = c("gene", "ensembl_transcript_id")) %>%
  left_join(cds_lengths, by = c("gene", "ensembl_transcript_id")) %>%
  filter(!is.na(nmdesc), !is.na(cds_length), cds_length > 0) %>%
  mutate(
    ratio = nmdesc / cds_length
  )


# Plot the ratio
length_ratio_df$nc = length_ratio_df$nmdesc/length_ratio_df$cds_length
#keep one row for each gene
length_ratio_df <- length_ratio_df %>%
  arrange(gene) %>%            
  distinct(gene, .keep_all = TRUE)

ratio_plot = ggplot(length_ratio_df, aes(x = category.x, y = nc, fill = category.x)) +
  geom_boxplot(width = 0.4, color = "black", outlier.shape = NA) +
  #geom_jitter(width = 0.2, size = 1.2, alpha = 0.6, color = "black") +
  scale_y_log10() +
  theme_minimal() +
  labs(
    title = "Ratio of Canonical NMDesc Length to CDS Length",
    y = "NMDesc / CDS Length",
    x = "Gene Category"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")
# libs
library(dplyr)
library(ggplot2)
library(ggpubr)

# Helper: rename + order labels
rename_levels <- function(x) {
  recoded <- dplyr::recode(tolower(x),
                           "minus1" = "Minus1",
                           "minus1_control" = "Minus1_Control",
                           "plus1" = "Plus1",
                           "plus1_control" = "Plus1_Control",
                           "snv" = "Nonsense",
                           "snv_control" = "Nonsense_Control",
                           .default = x
  )
  factor(recoded, levels = c("Minus1","Minus1_Control","Plus1","Plus1_Control","Nonsense","Nonsense_Control"))
}

# Apply renaming to category.x
length_ratio_df2 <- length_ratio_df %>%
  mutate(category.x = rename_levels(category.x)) %>%
  filter(!is.na(category.x)) %>%
  droplevels()

# Update comparisons to new labels and keep only valid pairs
comparisons_new <- lapply(comparisons, function(p) as.character(rename_levels(p)))
present_groups  <- levels(length_ratio_df2$category.x)
comparisons_use <- Filter(function(p) length(p)==2 && all(p %in% present_groups), comparisons_new)

# Colors mapped to new labels (edit if you already have a palette)
group_colors2 <- c(
  Minus1="#1f77b4", Minus1_Control="#aec7e8",
  Plus1 ="#ff7f0e", Plus1_Control ="#ffbb78",
  Nonsense="#2ca02c", Nonsense_Control="#98df8a"
)

# Plot
AD_ratio = length_ratio_df2 %>%
  filter(gene %in% pli_AD$gene)
p_lenratio <- ggplot(AD_ratio, aes(x = category.x, y = nc, fill = category.x)) +
  geom_boxplot(width = 0.4, color = "black", outlier.shape = NA) +
  scale_y_log10() +
  theme_minimal() +
  labs(
    title = "Ratio of Canonical NMDesc Length to CDS Length in AD genes",
    y = "NMDesc / CDS Length",
    x = "Gene category"                  # <- x-axis label text
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),  # <- bold ticks
    axis.title.x = element_text(size = 12, face = "bold"),              # <- bold x label
    plot.title   = element_text(hjust = 0.5, face = "bold"),            # <- bold centered title
    legend.position = "none"
  ) +
  scale_fill_manual(values = group_colors2) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  ggpubr::stat_compare_means(comparisons = comparisons_use, method = "wilcox.test")

print(p_lenratio)

ggplot(length_ratio_df, aes(x = category.x, y = nc, fill = category.x)) +
  geom_boxplot(width = 0.4, color = "black", outlier.shape = NA) +
  # geom_jitter(width = 0.2, size = 1.2, alpha = 0.6, color = "black") +
  scale_y_log10() +
  theme_minimal() +
  labs(
    title = "Ratio of Canonical NMDesc Length to CDS Length",
    y = "NMDesc / CDS Length",
    x = "Gene Category"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_manual(values = group_colors) +
  scale_x_discrete(labels = c(
    "SNV" = "Nonsense",
    "SNV_Control" = "Nonsense_Control",
    "plus1" = "Plus1",
    "plus1_Control" = "Plus1_Control",
    "minus1" = "Minus1",
    "minus1_Control" = "Minus1_Control"
  )) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")

ggsave("nmdesc_cds_length_ratio.pdf", plot = ratio_plot, width = 8, height = 5)

#7. plot pfam domain overlap
pfam_domains <- getBM(
  attributes = c("ensembl_transcript_id","pfam", "pfam_start", "pfam_end"),
  filters = "ensembl_transcript_id",
  values = canonical_tx$ensembl_transcript_id,
  mart = ensembl
)
#merge pfam_domains with exons, only keep rows in pfam_domains
pfam_merged = merge(
  pfam_domains,
  exons_joined,
  by = "ensembl_transcript_id",
  all.x = TRUE,all.y = FALSE
)
pfam_fin = pfam_merged %>%
  arrange(ensembl_transcript_id) %>%            
  distinct(ensembl_transcript_id, .keep_all = TRUE)

pfam_fin %>% 
  group_by(ensembl_transcript_id) %>%
  mutate(
    exon_length = abs(cds_end - cds_start) + 1,
    cucds_start = cumsum(lag(exon_length, default = 0)) + 1,
    cucds_end = cds_start + exon_length - 1,
    num_exon = row_number(),
  ) %>%
  ungroup()
# get NMDesc region
nmdesc_regions <- exons_joined %>%
  group_by(ensembl_transcript_id) %>%
  arrange(rank) %>%
  summarise(
    gene = first(hgnc_symbol),
    last_exon_start = cds_start[which.max(rank)],
    last_exon_end = cds_end[which.max(rank)],
    penult_exon_start = cds_start[rank == max(rank) - 1],
    penult_exon_end = cds_end[rank == max(rank) - 1],
    .groups = "drop"
  ) %>%
  mutate(
    penult_50bp_start = pmax(penult_exon_end - 49, penult_exon_start), #if pen exon < 50bp, use it all
    nmdesc_start = penult_50bp_start,
    nmdesc_end = last_exon_end
  ) %>%
  dplyr::select(ensembl_transcript_id, gene, nmdesc_start, nmdesc_end)


# compare pfam with NMDesc region
overlap_df <- inner_join(nmdesc_regions, pfam_domains, by = "ensembl_transcript_id") %>%
  mutate(
    overlap_start = pmax(nmdesc_start, pfam_start*3),
    overlap_end = pmin(nmdesc_end, pfam_end*3),
    overlap_length = pmax(overlap_end - overlap_start + 1, 0),
    nmdesc_length = nmdesc_end - nmdesc_start + 1,
    overlap_percent = ifelse(overlap_length > 0, overlap_length / nmdesc_length, 0)
  )

# sort by gene
overlap_summary <- overlap_df %>%
  group_by(ensembl_transcript_id) %>%
  summarise(overlap_flag = ifelse(any(overlap_length > 20), 1, 0), .groups = "drop") %>%
  inner_join(canonical_tx, by = "ensembl_transcript_id") %>%
  mutate(gene = hgnc_symbol) %>%
  distinct(gene, .keep_all = TRUE)


pli_overlap <- pli_all %>%
  left_join(overlap_summary, by = "gene") %>%
  #remove rows where overlap_flag is NA
  filter(!is.na(overlap_flag)) %>%
  mutate(
    overlap_flag = as.factor(overlap_flag),
    category = factor(category, levels = c("minus1", "plus1", "SNV", "minus1_Control", "plus1_Control", "SNV_Control"))
  )


pli_overlap %>% group_by(category) %>%
  summarise(
    overlap_count = sum(overlap_flag == 1),
    total_count = n(),
    overlap_fraction = overlap_count / total_count,
    .groups = "drop"
  )

pfam_plot = ggplot(pli_overlap, aes(x = category, fill = as.factor(overlap_flag))) +
  geom_bar(position = "fill", color = "black") +
  scale_fill_manual(
    values = c("0" = "lightgray", "1" = "steelblue"),
    labels = c("0" = "", "1" = "Yes"),
    breaks = "1"  # only show "1" in legend
  ) + 
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Fraction of Genes with Pfam Domain Overlap in NMDesc Region",
    x = "Gene Category",
    y = "Percent of Genes",
    fill = "Overlap"
  ) +
  scale_fill_manual(values = c("0" = "lightgray", "1" = "steelblue"), labels = c("No", "Yes")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("pfam_overlap_distribution.pdf", plot = pfam_plot, width = 8, height = 5)

pli_gene_AD = pli_overlap[pli_overlap$gene %in% pli_AD$gene,]
write.csv(pli_gene_AD, "pfam_overlap_AD_genes.csv", row.names = FALSE)