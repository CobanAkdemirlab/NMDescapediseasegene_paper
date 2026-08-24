
# --- path resolution layer (auto-inserted) ------------------------------------------------
# Resolves paths via paths.R:
#   data_file("x.csv")  locate by filename, errors clearly if missing
#   out_file("y.csv")   output to NMDESC_OUT (default ~/Desktop/NMDesc_out)
#   data_root("clinvar") use when a directory is needed instead of a file
.p <- c("gene level_v3/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../lib/paths.R", "../../../../lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p)) stop("paths.R not found -- run R from the repository root")
source(.p[1]); rm(.p)
# --------------------------------------------------------------------------

#This Rscript is to combine gene list from disease and control, then add necessary variables including NMDesc_region_length
snv_gene           = read_csv(data_file("snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201_NMDesc_wald_enriched_can.txt"), 
                              col_names = FALSE)
fs_gene            = read_csv("fs_can_AD_FDR0.05_wald_gene.csv")
snv_control_gene_AD = read_csv("snv_control_genes_AD.csv")
fs_control_gene_AD  = read_csv("fs_control_genes_AD.csv")
PTC_info = read.csv('PTC_info20260201_region.csv')

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

fs_control_gene = getBM(
  attributes = c("hgnc_symbol",'ensembl_transcript_id'),
  filters = "ensembl_transcript_id",
  values = unique(fs_control_gene_AD),
  mart = ensembl
)
snv_control_gene = getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = unique(snv_control_gene_AD),
  mart = ensembl
)

#get cds sequence
snv_gene_cds = getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_transcript_id",
    "transcript_is_canonical",
    "coding"
  ),
  filters = "hgnc_symbol",
  values = snv_gene,
  mart = ensembl
)
fs_gene_cds <- getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_transcript_id",
    "transcript_is_canonical",
    "coding"
  ),
  filters = "hgnc_symbol",
  values = fs_gene,
  mart = ensembl
)
fs_control_gene_cds = getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_transcript_id",
    "transcript_is_canonical",
    "coding"
  ),
  filters = "hgnc_symbol",
  values = fs_control_gene$hgnc_symbol,
  mart = ensembl
)
snv_control_gene_cds = getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_transcript_id",
    "transcript_is_canonical",
    "coding"
  ),
  filters = "ensembl_transcript_id",
  values = snv_control_gene$ensembl_transcript_id,
  mart = ensembl
)
#keep canonical transcript only
snv_gene_cds_can <- snv_gene_cds %>%
  filter(transcript_is_canonical == 1)
fs_gene_cds_can <- fs_gene_cds %>%
  filter(transcript_is_canonical == 1)
fs_control_gene_cds_can <- fs_control_gene_cds %>%
  filter(transcript_is_canonical == 1)
snv_control_gene_cds_can <- snv_control_gene_cds %>%
  filter(transcript_is_canonical == 1)
#remove empty CDS
snv_gene_cds_can <- snv_gene_cds_can %>%
  filter(!is.na(coding), coding != "")
fs_gene_cds_can <- fs_gene_cds_can %>%
  filter(!is.na(coding), coding != "")
fs_control_gene_cds_can <- fs_control_gene_cds_can %>%
  filter(!is.na(coding), coding != "")
snv_control_gene_cds_can <- snv_control_gene_cds_can %>%
  filter(!is.na(coding), coding != "")


##get snv NMDesc region: last exon + 50bp upstream of penutimate exon boundary
#get the cds location of exons
snv_exon_info = getBM(
  attributes = c(
    "ensembl_transcript_id",
    "rank",
    "cds_start",
    "cds_end"
  ),
  filters = "hgnc_symbol",
  values = unique(snv_gene),
  mart = ensembl
)
snv_control_exon_info = getBM(
  attributes = c(
    "ensembl_transcript_id",
    "rank",
    "cds_start",
    "cds_end"
  ),
  filters = "ensembl_transcript_id",
  values = unique(snv_control_gene$ensembl_transcript_id),
  mart = ensembl
)

snv_tx_info = getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_transcript_id",
    "transcript_is_canonical"
  ),
  filters = "hgnc_symbol",
  values = unique(snv_gene),
  mart = ensembl
)
snv_control_tx_info = getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_transcript_id",
    "transcript_is_canonical"
  ),
  filters = "ensembl_transcript_id",
  values = unique(snv_control_gene$ensembl_transcript_id),
  mart = ensembl
)
snv_tx_info = snv_tx_info %>%
  filter(transcript_is_canonical == 1) %>%
  distinct(hgnc_symbol, ensembl_transcript_id, .keep_all = TRUE)
snv_control_tx_info = snv_control_tx_info %>%
  filter(transcript_is_canonical == 1) %>%
  distinct(hgnc_symbol, ensembl_transcript_id, .keep_all = TRUE)
snv_exon_info = snv_exon_info %>%
  inner_join(snv_tx_info, by = "ensembl_transcript_id") %>%
  filter(!is.na(rank), !is.na(cds_start), !is.na(cds_end))
snv_control_exon_info = snv_control_exon_info %>%
  inner_join(snv_control_tx_info, by = "ensembl_transcript_id") %>%
  filter(!is.na(rank), !is.na(cds_start), !is.na(cds_end))

snv_nmdesc_info = 
  snv_exon_info %>% 
  group_by(ensembl_transcript_id) %>%
  summarise(nmdesc_end = max(cds_end), #cds end 
            last_exon_length = cds_end[which.max(rank)] - cds_start[which.max(rank)] + 1,
            nmdesc_start = max(cds_end) - 50 - last_exon_length) #cds_end - (length of last exon + 50bp)
snv_control_nmdesc_info = 
  snv_control_exon_info %>% 
  group_by(ensembl_transcript_id) %>%
  summarise(nmdesc_end = max(cds_end), #cds end 
            last_exon_length = cds_end[which.max(rank)] - cds_start[which.max(rank)] + 1,
            nmdesc_start = max(cds_end) - 50 - last_exon_length) #cds_end - (length of last exon + 50bp)

snv_nmdesc_df = snv_gene_cds_can %>%
  inner_join(snv_nmdesc_info, by = "ensembl_transcript_id")
snv_control_nmdesc_df = snv_control_gene_cds_can %>%
  inner_join(snv_control_nmdesc_info, by = "ensembl_transcript_id")

snv_nmdesc_df$can_region_start = snv_nmdesc_df$nmdesc_start
snv_control_nmdesc_df$can_region_start = snv_control_nmdesc_df$nmdesc_start
snv_nmdesc_df$can_region_end = snv_nmdesc_df$nmdesc_end
snv_control_nmdesc_df$can_region_end = snv_control_nmdesc_df$nmdesc_end

snv_nmdesc_df$snv_nmdesc_cds = mapply(
  substr,
  snv_nmdesc_df$coding,
  snv_nmdesc_df$can_region_start,
  snv_nmdesc_df$can_region_end
)
snv_control_nmdesc_df$snv_control_nmdesc_cds = mapply(
  substr,
  snv_control_nmdesc_df$coding,
  snv_control_nmdesc_df$can_region_start,
  snv_control_nmdesc_df$can_region_end
)

snv_nmdesc_df$snv_nmdesc_gc_content <- sapply(snv_nmdesc_df$snv_nmdesc_cds, get_gc_content)
snv_control_nmdesc_df$snv_control_nmdesc_gc_content <- sapply(snv_control_nmdesc_df$snv_control_nmdesc_cds, get_gc_content)


## get fs NMDesc region: derive from PTC_info$can_region_start/can_region_end, using median for plus1 and plus2
PTC_combined = PTC_info %>%
  dplyr::group_by(transcript) %>%
  dplyr::summarise(
    median_can_region_start = round(median(can_region_start, na.rm = TRUE)),
    median_can_region_end = round(median(can_region_end, na.rm = TRUE)),
    .groups = "drop"
  )

fs_nmdesc_df = fs_gene_cds_can %>%
  dplyr::inner_join(PTC_combined, by = c("ensembl_transcript_id" = "transcript"))
fs_control_nmdesc_df = fs_control_gene_cds_can %>%
  dplyr::inner_join(PTC_combined, by = c("ensembl_transcript_id" = "transcript"))

#get cds sequence in the NMDesc region from fs_gene_cds_can$coding
fs_nmdesc_df$fs_nmdesc_cds = mapply(
  substr,
  fs_nmdesc_df$coding,
  fs_nmdesc_df$median_can_region_start,
  fs_nmdesc_df$median_can_region_end
)
fs_control_nmdesc_df$fs_control_nmdesc_cds = mapply(
  substr,
  fs_control_nmdesc_df$coding,
  fs_control_nmdesc_df$median_can_region_start,
  fs_control_nmdesc_df$median_can_region_end
)

fs_nmdesc_df$fs_nmdesc_gc_content <- sapply(fs_nmdesc_df$fs_nmdesc_cds, get_gc_content)
fs_control_nmdesc_df$fs_control_nmdesc_gc_content <- sapply(fs_control_nmdesc_df$fs_control_nmdesc_cds, get_gc_content)


fs_nmdesc_df2 = fs_nmdesc_df %>%
  rename(NMD_region_start = median_can_region_start,
         NMD_region_end = median_can_region_end,
         nmdesc_cds = fs_nmdesc_cds,
         nmdesc_gc_content = fs_nmdesc_gc_content) 
fs_control_nmdesc_df2 = fs_control_nmdesc_df %>%
  rename(NMD_region_start = median_can_region_start,
         NMD_region_end = median_can_region_end,
         nmdesc_cds = fs_control_nmdesc_cds,
         nmdesc_gc_content = fs_control_nmdesc_gc_content) 
snv_nmdesc_df2 = snv_nmdesc_df %>%
  rename(NMD_region_start = can_region_start,
         NMD_region_end = can_region_end,
         nmdesc_cds = snv_nmdesc_cds,
         nmdesc_gc_content = snv_nmdesc_gc_content) %>%
  dplyr::select(-any_of(c("last_exon_length", "nmdesc_start",'nmdesc_end')))
snv_control_nmdesc_df2 = snv_control_nmdesc_df %>%
  rename(NMD_region_start = can_region_start,
         NMD_region_end = can_region_end,
         nmdesc_cds = snv_control_nmdesc_cds,
         nmdesc_gc_content = snv_control_nmdesc_gc_content) %>%
  dplyr::select(-any_of(c("last_exon_length", "nmdesc_start",'nmdesc_end')))
fs_nmdesc_df2 <- fs_nmdesc_df2 %>%
  mutate(transcript_is_canonical = as.character(transcript_is_canonical))

fs_control_nmdesc_df2 <- fs_control_nmdesc_df2 %>%
  mutate(transcript_is_canonical = as.character(transcript_is_canonical))

snv_nmdesc_df2 <- snv_nmdesc_df2 %>%
  mutate(transcript_is_canonical = as.character(transcript_is_canonical))

snv_control_nmdesc_df2 <- snv_control_nmdesc_df2 %>%
  mutate(transcript_is_canonical = as.character(transcript_is_canonical))

gene_all <- bind_rows(
  fs_nmdesc_df2 %>% mutate(group = "fs"),
  fs_control_nmdesc_df2 %>% mutate(group = "fs_control"),
  snv_nmdesc_df2 %>% mutate(group = "snv"),
  snv_control_nmdesc_df2 %>% mutate(group = "snv_control")
)

gene_all$NMDesc_region_length = gene_all$NMD_region_end - gene_all$NMD_region_start + 1
#use canonical transcript to get cds length
cds_info <- getBM(
  attributes = c("ensembl_transcript_id", "cds_length"),
  filters = "ensembl_transcript_id",
  values = gene_all$ensembl_transcript_id,
  mart = ensembl
)
gene_all$cds_length = cds_info$cds_length[match(gene_all$ensembl_transcript_id,cds_info$ensembl_transcript_id)]
gene_all$row_id <- seq_len(nrow(gene_all))

# extract transcript ids
tx_ids <- unique(gene_all$ensembl_transcript_id)

# query transcript -> UniProt Swiss-Prot
tx_uniprot <- getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters    = "ensembl_transcript_id",
  values     = tx_ids,
  mart       = ensembl
)


# merge back into gene_all
gene_all <- gene_all %>%
  left_join(tx_uniprot, by = "ensembl_transcript_id")

#rename uniprotswissprot to uniprot for easier use
gene_all <- gene_all %>%
  rename(uniprot = uniprotswissprot)


write.csv(gene_all, "gene_all0407.csv", row.names = FALSE)