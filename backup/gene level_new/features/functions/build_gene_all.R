# ── Helper: fetch all required info for a gene list ────────────────────────
get_gene_info <- function(
    hgnc_symbols = NULL,        # provide either hgnc_symbols OR ensembl_tx_ids
    ensembl_tx_ids = NULL,
    mart,
    PTC_combined = NULL,        # required only for fs / fs_control groups
    group_type = c("snv", "fs", "snv_control", "fs_control")
) {
  group_type <- match.arg(group_type)
  filter_by_symbol <- !is.null(hgnc_symbols)
  filter_val  <- if (filter_by_symbol) unique(hgnc_symbols) else unique(ensembl_tx_ids)
  filter_name <- if (filter_by_symbol) "hgnc_symbol" else "ensembl_transcript_id"
  
  # ── 1. CDS sequence (canonical only, non-empty) ───────────────────────────
  cds_df <- getBM(
    attributes = c("hgnc_symbol", "ensembl_transcript_id",
                   "transcript_is_canonical", "coding"),
    filters    = filter_name,
    values     = filter_val,
    mart       = mart
  ) %>%
    filter(transcript_is_canonical == 1,
           !is.na(coding), coding != "")
  
  # ── 2. NMD-escape region ──────────────────────────────────────────────────
  if (group_type %in% c("snv", "snv_control")) {
    
    # Exon-level CDS coordinates
    exon_info <- getBM(
      attributes = c("ensembl_transcript_id", "rank", "cds_start", "cds_end"),
      filters    = filter_name,
      values     = filter_val,
      mart       = mart
    )
    
    # Canonical transcript map
    tx_info <- getBM(
      attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_is_canonical"),
      filters    = filter_name,
      values     = filter_val,
      mart       = mart
    ) %>%
      filter(transcript_is_canonical == 1) %>%
      distinct(hgnc_symbol, ensembl_transcript_id, .keep_all = TRUE)
    
    exon_info <- exon_info %>%
      inner_join(tx_info, by = "ensembl_transcript_id") %>%
      filter(!is.na(rank), !is.na(cds_start), !is.na(cds_end))
    
    # NMD-escape region = last exon + 50 bp upstream of penultimate exon boundary
    nmdesc_info <- exon_info %>%
      group_by(ensembl_transcript_id) %>%
      summarise(
        nmdesc_end        = max(cds_end),
        last_exon_length  = cds_end[which.max(rank)] - cds_start[which.max(rank)] + 1,
        nmdesc_start      = max(cds_end) - 50 - last_exon_length,
        .groups = "drop"
      )
    
    nmdesc_df <- cds_df %>%
      inner_join(nmdesc_info, by = "ensembl_transcript_id") %>%
      mutate(
        NMD_region_start = nmdesc_start,
        NMD_region_end   = nmdesc_end,
        nmdesc_cds       = mapply(substr, coding, nmdesc_start, nmdesc_end)
      ) %>%
      select(-any_of(c("last_exon_length", "nmdesc_start", "nmdesc_end")))
    
  } else {
    # fs / fs_control: use median PTC region across plus1 / plus2 alleles
    stopifnot(!is.null(PTC_combined),
              all(c("transcript", "median_can_region_start",
                    "median_can_region_end") %in% names(PTC_combined)))
    
    nmdesc_df <- cds_df %>%
      inner_join(PTC_combined, by = c("ensembl_transcript_id" = "transcript")) %>%
      mutate(
        NMD_region_start  = median_can_region_start,
        NMD_region_end    = median_can_region_end,
        nmdesc_cds        = mapply(substr, coding,
                                   median_can_region_start, median_can_region_end)
      ) %>%
      select(-median_can_region_start, -median_can_region_end)
  }
  
  nmdesc_df %>%
    mutate(
      transcript_is_canonical = as.character(transcript_is_canonical),
      group = group_type
    )
}


# ── Build gene_all ──────────────────────────────────────────────────────────
build_gene_all <- function(
    snv_gene,
    fs_gene,
    snv_control_gene_AD,   # vector of transcript IDs
    fs_control_gene_AD,    # vector of transcript IDs
    PTC_info,
    mart,
    output_csv = "gene_all.csv"
) {
  # ── Pre-compute PTC median region (shared by fs and fs_control) ───────────
  PTC_combined <- PTC_info %>%
    group_by(transcript) %>%
    summarise(
      median_can_region_start = round(median(can_region_start, na.rm = TRUE)),
      median_can_region_end   = round(median(can_region_end,   na.rm = TRUE)),
      .groups = "drop"
    )
  
  # ── Fetch info for each group ─────────────────────────────────────────────
  snv_df         <- get_gene_info(hgnc_symbols   = snv_gene,              mart = mart, group_type = "snv")
  fs_df          <- get_gene_info(hgnc_symbols   = fs_gene,               mart = mart, group_type = "fs",          PTC_combined = PTC_combined)
  snv_control_df <- get_gene_info(ensembl_tx_ids = snv_control_gene_AD,   mart = mart, group_type = "snv_control")
  fs_control_df  <- get_gene_info(ensembl_tx_ids = fs_control_gene_AD,    mart = mart, group_type = "fs_control",  PTC_combined = PTC_combined)
  
  # ── Combine ───────────────────────────────────────────────────────────────
  gene_all <- bind_rows(snv_df, fs_df, snv_control_df, fs_control_df) %>%
    mutate(NMDesc_region_length = NMD_region_end - NMD_region_start + 1,
           row_id = row_number())
  
  # ── Append CDS length ─────────────────────────────────────────────────────
  cds_len <- getBM(
    attributes = c("ensembl_transcript_id", "cds_length"),
    filters    = "ensembl_transcript_id",
    values     = unique(gene_all$ensembl_transcript_id),
    mart       = mart
  )
  gene_all <- gene_all %>%
    left_join(cds_len, by = "ensembl_transcript_id")
  
  # ── Append UniProt Swiss-Prot ID ──────────────────────────────────────────
  uniprot <- getBM(
    attributes = c("ensembl_transcript_id", "uniprotswissprot"),
    filters    = "ensembl_transcript_id",
    values     = unique(gene_all$ensembl_transcript_id),
    mart       = mart
  )
  gene_all <- gene_all %>%
    left_join(uniprot, by = "ensembl_transcript_id") %>%
    rename(uniprot = uniprotswissprot)
  #filter for AD genes
  omim_AD_symbols = read.csv('omim_AD_symbols.csv',header=F)$V1
  gene_all <- gene_all %>%
    filter(hgnc_symbol %in% omim_AD_symbols)
  # ── Save & return ─────────────────────────────────────────────────────────
  write.csv(gene_all, output_csv, row.names = FALSE)
  message("Saved: ", output_csv)
  invisible(gene_all)
}


# ── Usage ───────────────────────────────────────────────────────────────────
PTC_info <-  read.csv('PTC_info20260201_region.csv')
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 105)
gene_all <- build_gene_all(
  snv_gene           = read_csv("~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/lists/FDR0.05/snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201_NMDesc_wald_enriched_can.txt", 
                                 col_names = FALSE),
  fs_gene            = read_csv("fs_can_AD_FDR0.05_wald_gene.csv"),
  snv_control_gene_AD = read_csv("snv_control_genes_AD.csv"),
  fs_control_gene_AD  = read_csv("fs_control_genes_AD.csv"),
  PTC_info           = PTC_info,
  mart               = ensembl,
  output_csv         = "gene_all.csv"
)