library(dplyr)
# --- Path resolution layer (auto-inserted) --------------------------------
# Path helpers from paths.R:
#   data_file("x.csv") locates by filename, errors clearly if missing
#   out_file("y.csv")  writes to NMDESC_OUT (default ~/Desktop/NMDesc_out)
#   data_root("clinvar") for directories instead of files
.p <- c("gene level_v3/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../lib/paths.R", "../../../../lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p)) stop("paths.R not found -- run R from the repository root")
source(.p[1]); rm(.p)
# --------------------------------------------------------------------------


###1. get snv variants in gnomAD using ptc_can_df
#input snv disease gene list
#snv_gene = read_csv("snv_plp_ptc_nmdesc_can_p_f_syn_20260201_NMDesc_enriched_can_AD_p_0.8.csv")
snv_gene = read_csv(data_file("snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201_NMDesc_binom_enriched_can.txt"), 
                   col_names = FALSE)
ptc_can_NMD_df2 =  read_csv("~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/gnomad/snv_fs/ptc_can_NMD_df.csv")
#filter for gnomad variants corresponding to the gene list
gnomad_snv = ptc_can_NMD_df2[which(ptc_can_NMD_df2$transcript %in% snv_can_tr & ptc_can_NMD_df2$type == 'snv'),]
#remove plp and vus variants in clinvar
clinvar_snv_plp = snv_variants$key
clinvar_snv_vus = snv_vus_variants$key
gnomad_snv_filtered = gnomad_snv[!gnomad_snv$id %in% clinvar_snv_plp & !gnomad_snv$id %in% clinvar_snv_vus,]

###2. get frameshift variants in gnomAD using ptc_can_df
#input fs disease gene list
#fs_gene = read.csv('fs_can_syn_AD_gene_filtered_greater_20260201.txt', header = TRUE, stringsAsFactors = FALSE)
fs_gene = read_csv(data_file("fs_can_AD_FDR0.05_gene.csv", must = FALSE))
#get canonical transcripts
fs_can_tr <- getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "transcript_is_canonical"
  ),
  filters = "hgnc_symbol",
  values = fs_gene$hgnc_symbol,
  mart = ensembl
)

fs_can_tr <- fs_can_tr %>%
  filter(transcript_is_canonical == 1)

gnomad_fs = ptc_can_NMD_df2[which(ptc_can_NMD_df2$transcript %in% fs_can_tr$ensembl_transcript_id & ptc_can_NMD_df2$type != 'snv'),]
length(unique(gnomad_fs$transcript))
#remove inframe frameshift variants based on the key
# Extract ref and alt from id, compute length diff
gnomad_fs_filtered = read.csv('gnomad_fs_filtered.csv')
gnomad_fs_filtered$ref = sub(".*\\|(.+)\\|.*", "\\1", gnomad_fs_filtered$id)
gnomad_fs_filtered$alt = sub(".*\\|.*\\|(.*)", "\\1", gnomad_fs_filtered$id)
gnomad_fs_filtered$len_diff = abs(nchar(gnomad_fs_filtered$ref) - nchar(gnomad_fs_filtered$alt))

# Keep variants where length diff is not a multiple of 3 (true frameshift)
gnomad_fs_filtered = gnomad_fs_filtered[gnomad_fs_filtered$len_diff %% 3 != 0, ]

# Remove temporary columns
gnomad_fs_filtered$ref = NULL
gnomad_fs_filtered$alt = NULL
gnomad_fs_filtered$len_diff = NULL

#add uniprot id using getBM
gnomad_fs_unique = unique(gnomad_fs_filtered$transcript)
gnomad_snv_unique = unique(gnomad_snv_filtered$transcript)
gnomad_fs_uniprot_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = gnomad_fs_unique,
  mart = ensembl
)
gnomad_snv_uniprot_mapping <- getBM(
  attributes = c("ensembl_transcript_id", "uniprotswissprot"),
  filters = "ensembl_transcript_id",
  values = gnomad_snv_unique,
  mart = ensembl
)
gnomad_fs_filtered2 = merge(
  gnomad_fs_filtered,
  gnomad_fs_uniprot_mapping,
  by.x = "transcript",
  by.y = "ensembl_transcript_id",
  all.x = TRUE
)
gnomad_snv_filtered2 = merge(
  gnomad_snv_filtered,
  gnomad_snv_uniprot_mapping,
  by.x = "transcript",
  by.y = "ensembl_transcript_id",
  all.x = TRUE
)

write.csv(gnomad_snv_filtered, 'gnomad_snv_filtered_wald.csv', row.names = FALSE)
write.csv(gnomad_fs_filtered, 'gnomad_fs_filtered_wald.csv', row.names = FALSE)
length(unique(gnomad_fs_filtered_wald$transcript))

#for each tr, do a table to show how many clinvar plp and gnomad benign variants
snv_summary = snv_variants %>%
  left_join(snv_tr, by = c("transcript" = "ensembl_transcript_id")) %>%
  filter(hgnc_symbol %in% snv_gene$V1) %>%
  group_by(hgnc_symbol) %>%
  summarise(
    n_clinvar_snv = n_distinct(key),
    .groups = "drop"
  ) %>%
  left_join(
    gnomad_snv_filtered %>%
      left_join(snv_tr, by = c("transcript" = "ensembl_transcript_id")) %>%
      filter(hgnc_symbol %in% snv_gene$V1) %>%
      group_by(hgnc_symbol) %>%
      summarise(n_gnomad_snv = n_distinct(id), .groups = "drop"),
    by = "hgnc_symbol"
  ) 

fs_summary = fs_variants2 %>%
  left_join(fs_tr, by = c("transcript" = "ensembl_transcript_id")) %>%
  filter(hgnc_symbol %in% fs_gene$hgnc_symbol) %>%
  group_by(hgnc_symbol) %>%
  summarise(
    n_clinvar_fs = n_distinct(key),
    .groups = "drop"
  ) %>%
  left_join(
    gnomad_fs_filtered %>%
      left_join(fs_tr, by = c("transcript" = "ensembl_transcript_id")) %>%
      filter(hgnc_symbol %in% fs_gene$hgnc_symbol) %>%
      group_by(hgnc_symbol) %>%
      summarise(n_gnomad_fs = n_distinct(id), .groups = "drop"),
    by = "hgnc_symbol"
  ) 

print(snv_summary)
print(fs_summary)

write.csv(snv_summary, 'snv_gene_summary_wald.csv', row.names = FALSE)
write.csv(fs_summary, 'fs_gene_summary_wald.csv', row.names = FALSE)
#remove extra long ones