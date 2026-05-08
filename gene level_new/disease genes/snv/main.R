##Library aenmd into R
library(aenmd.data.ensdb.v105)
library(aenmd)
library(GenomicRanges)
library(tidyr)
library(GenomicFeatures)
library(VariantAnnotation)
library(AnnotationDbi)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biomaRt)
library(enrichR)
library(dplyr)
library(readr)
library(GenomicFeatures)
library(VariantAnnotation)
library(AnnotationDbi)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)

#RARB, DAGLA, MN1

#download gnomad again

genome <- BSgenome.Hsapiens.UCSC.hg38

txdb <- makeTxDbFromEnsembl("Homo sapiens", release=105)

#saveDb(txdb, 'txdb.Ensembl105.sqlite')

ensgene <- txdb

seqlevels(ensgene) <- paste('chr',seqlevels(ensgene),sep='')

chr.list <- c(paste('chr',1:22,sep=''))

ensgene.sub <- keepSeqlevels(ensgene,chr.list)
ensgene <- ensgene.sub

cds_seqs <- GenomicFeatures::extractTranscriptSeqs(genome, cdsBy(ensgene, by="tx",use.names=TRUE))
threeutr.grange <-threeUTRsByTranscript(ensgene, use.names=TRUE)
fiveutr.grange<-fiveUTRsByTranscript(ensgene, use.names=TRUE)
introns.grange<- intronsByTranscript(ensgene, use.names=TRUE)

###get DNAstring of UTR seqs
five_utr_seqs <- extractTranscriptSeqs(Hsapiens, fiveutr.grange) #what is hsapiens?
three_utr_seqs <- extractTranscriptSeqs(Hsapiens, threeutr.grange)

keys <- names(cds_seqs)
cols <- columns(ensgene)
all.df <- AnnotationDbi::select(ensgene, keys = keys, columns = cols, keytype="TXNAME")
## this has length of 304803 so most contain NAs or duplicates?
##cds_seqs length is 28856
all.df.sub <- all.df[which(all.df$TXNAME!='NA'),]
rm(all.df)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#tx_vec <- unique(as.character(res_p@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]]))
#tx_vec <- tx_vec[!is.na(tx_vec)]
BM.info <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","hgnc_symbol","transcript_is_canonical"),mart=ensembl)
#BM.info <- getBM(
#  attributes = c("ensembl_transcript_id", 'transcript_is_canonical','hgnc_symbol'),
#  filters    = "ensembl_transcript_id",
#  values     = tx_vec,
#  mart       = ensembl
#)
write.csv(BM.info,'BM_info.csv',row.names = F)

vcf_file = "clinvar_20260201.vcf.gz"
vcf = aenmd:::parse_vcf_VariantAnnotation(vcf_file)
vcf_rng = vcf$vcf_rng
#rm(vcf)
#add clinical significance
clnsig_list <- info(vcf$vcf_obj)$CLNSIG 
mcols(vcf_rng)$CLNSIG <- clnsig_list
vcf_rng_fil = process_variants(vcf_rng)
#- filter out variants with ill-defined alternative allele
ind_out =  Biostrings::vcountPattern("N", vcf_rng_fil$alt) > 0
vcf_rng_fil = vcf_rng_fil[!ind_out]
#- back to the original workflow
res = annotate_nmd(vcf_rng_fil, rettype="gr")

snv_ind = which(res@elementMetadata@listData[["type"]] == 'snv')
clnsig_str = sapply(res$CLNSIG, function(x) paste(as.character(x), collapse="|"))
plp_ind = which(grepl("pathogenic", clnsig_str, ignore.case = TRUE))
length(unique(res@elementMetadata@listData[["key"]][plp_ind]))

benign_ind = which(grepl("benign", clnsig_str, ignore.case = TRUE))
ptc_ind = which(res@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]]==T)
length(unique(res@elementMetadata@listData[["key"]][intersect(plp_ind,ptc_ind)]))

nmdesc_ind = which(res@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T | res@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)

#filter for canonical transcript using getBM
tx_vec <- res@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]]
tx_n <- data.frame(transcript = tx_vec) %>%
  filter(!is.na(transcript)) %>%
  group_by(transcript) %>%
  summarise(n = n(), .groups = "drop")

tx_can <- getBM(
  attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_is_canonical"),
  filters    = "ensembl_transcript_id",
  values     = unique(tx_n$transcript),
  mart       = ensembl
) %>%
  filter(transcript_is_canonical == 1) %>%
  distinct(hgnc_symbol, ensembl_transcript_id)
tx_n$canonical = ifelse(tx_n$transcript %in% tx_can$ensembl_transcript_id, "canonical", "non_canonical")
tx_n$hgnc_symbol = tx_can$hgnc_symbol[match(tx_n$transcript, tx_can$ensembl_transcript_id)]
can_ind = which(res@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% tx_can$ensembl_transcript_id)
#combine the filters
plp_ptc_ind = intersect(plp_ind,ptc_ind)
benign_ptc_ind = intersect(benign_ind,ptc_ind)
snv_plp_ptc_ind = intersect(snv_ind,plp_ptc_ind)
res_snv_plp_ptc = res[snv_plp_ptc_ind]
snv_benign_ptc_ind = intersect(snv_ind,benign_ptc_ind)
snv_plp_ptc_nmdesc_ind = intersect(snv_plp_ptc_ind,nmdesc_ind)
snv_benign_ptc_nmdesc_ind = intersect(snv_benign_ptc_ind,nmdesc_ind)
snv_plp_ptc_nmdesc_can_ind = intersect(snv_plp_ptc_nmdesc_ind,can_ind)
snv_benign_ptc_nmdesc_can_ind = intersect(snv_benign_ptc_nmdesc_ind,can_ind)

snv_plp_ptc_nmdesc_can = res[snv_plp_ptc_nmdesc_can_ind]
snv_benign_ptc_nmdesc_can = res[snv_benign_ptc_nmdesc_can_ind]
length(unique(snv_plp_ptc_nmdesc_can@ranges@NAMES))
length(unique(snv_plp_ptc_nmdesc_can@elementMetadata@listData[["key"]]))

#Quality Control, match by key
variant_summary <- read_delim("variant_summary.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

#get the total number of submitters for each gene
gene_all$NumberSubmitters = v_ch$NumberSubmitters[match(gene_all$k, v_ch$key)]
gene_all %>% group_by(group) %>% summarise(mean_submitters = mean(NumberSubmitters, na.rm = TRUE), median_submitters = median(NumberSubmitters, na.rm = TRUE), low_submitters = sum(NumberSubmitters <= 1, na.rm = TRUE))
#remove the low submitters, which may be more likely to be false positives
gene_all %>% filter(NumberSubmitters > 1) %>% group_by(group) %>% summarise(mean_submitters = mean(NumberSubmitters, na.rm = TRUE), median_submitters = median(NumberSubmitters, na.rm = TRUE), low_submitters = sum(NumberSubmitters <= 1, na.rm = TRUE))
gene_all4 = gene_all %>% filter(NumberSubmitters >= 3)
table(gene_all4$group)
#plot the number of submitters
ggplot(gene_all4, aes(x = group, y = NumberSubmitters,color = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Group",
    y = "Number of Submitters",
    title = "Number of Submitters by Group"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )
plot(x=gene_all4$NumberSubmitters, y=gene_all4$cds_length)


summary(gene_can_AD_plp_ptc_uni$NumberSubmitters)
#plot NumberSubmittters against cds_length
ggplot(gene_can_AD_plp_ptc_uni, aes(x = cds_length, y = NumberSubmitters)) +
  geom_point(color = "blue", alpha = 0.6) +
  theme_classic() +
  labs(
    x = "CDS Length",
    y = "Number of Submitters",
    title = "CDS Length vs Number of Submitters"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )
#add line of best fit and 95% CI
ggplot(gene_can_AD_plp_ptc_uni, aes(x = cds_length, y = NumberSubmitters)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", fill = "lightpink", se = TRUE) +
  theme_classic() +
  labs(
    x = "CDS Length",
    y = "Number of Submitters",
    title = "CDS Length vs Number of Submitters"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )



gene_can_AD_plp_ptc_nmdesc_uni$NumberSubmitters = variant_summary$NumberSubmitters[match(gene_can_AD_plp_ptc_nmdesc_uni$key, variant_summary$key)]
summary(gene_can_AD_plp_ptc_nmdesc_uni$NumberSubmitters)

# do other filters in variant_summary

#ReviewStatus does not contain either of “no assertion” or “no interpretation”
#variant_summary$ReviewStatus
v_re = variant_summary %>%
  filter(
    !str_detect(tolower(ReviewStatus), "no assertion"),
    !str_detect(tolower(ReviewStatus), "no interpretation")
  )
v_re[which(v_re$`#AlleleID` == '2831360'),]
#v_re2 = variant_summary_match %>%
#  filter(
#    !str_detect(tolower(ReviewStatus), "no assertion"),
#    !str_detect(tolower(ReviewStatus), "no interpretation")
#  )
#nrow(variant_summary_match) - nrow(v_re2)

#ClinicalSignificance does not contain any of “not provided”, “drug
#response”, “other”, “risk”, “low penetrance”, “conflicting”, “affects”, “association”, “protective”, “confers
#sensitivity”
#variant_summary$ClinicalSignificance
exclude_terms <- c(
  "not provided", "drug response", "other", "risk", "low penetrance",
  "conflicting", "affects", "association", "protective", "confers sensitivity"
)
pattern_exclude <- paste(exclude_terms, collapse = "|")
v_cs <- v_re %>%
  filter(!str_detect(ClinicalSignificance, regex(pattern_exclude, ignore_case = TRUE)))
v_cs[which(v_cs$`#AlleleID` == '2831360'),]
#v_cs2 = v_re2 %>%
#  filter(
#    !str_detect(tolower(ClinicalSignificance), pattern_exclude)
#  )
#nrow(v_re2) - nrow(v_cs2)
rm(v_re)
#Assembly == GRCh38
v_gr = v_cs %>%
  filter(str_trim(tolower(Assembly)) == "grch38")
v_gr[which(v_gr$`#AlleleID` == '2831360'),]
#v_gr2 = v_cs2 %>%
#  filter(str_trim(tolower(Assembly)) == "grch38")
#nrow(v_cs2) - nrow(v_gr2)
#autosomal contigs (chr1-22) only
v_ch = v_gr %>%
  mutate(
    Chromosome_clean = str_remove(tolower(str_trim(Chromosome)), "^chr")
  ) %>%
  filter(Chromosome_clean %in% as.character(1:22)) %>%
  dplyr::select(-Chromosome_clean)
v_ch[which(v_ch$`#AlleleID` == '2831360'),]
#v_ch2 = v_gr2 %>%
#  mutate(
#    Chromosome_clean = str_remove(tolower(str_trim(Chromosome)), "^chr")
#  ) %>%
#  filter(Chromosome_clean %in% as.character(1:22)) %>%
#  select(-Chromosome_clean)
#nrow(v_gr2) - nrow(v_ch2)
#add canonical and non_canonical transcript tag use Name
#substract from Name(eg. NM_017547.4 from NM_017547.4(FOXRED1):c.694C>T (p.Gln232Ter))
v_ch$tx_id1 = sub("\\(.*", "", v_ch$Name)
#v_ch2$tx_id1 = sub("\\(.*", "", v_ch2$Name)
v_ch = v_ch %>%
  mutate(
    tx_id1_noversion = sub("\\..*$", "", tx_id1)
  )
#v_ch2 = v_ch2 %>%
#  mutate(
#    tx_id1_noversion = sub("\\..*$", "", tx_id1)
#  )
#presence of a valid RefSeq ID
v_ch = v_ch %>%
  filter(str_detect(tx_id1_noversion, "^NM_"))
v_ch[which(v_ch$`#AlleleID` == '2831360'),]
#v_ch2 = v_ch2 %>%
#  filter(str_detect(tx_id1_noversion, "^NM_"))

#transfer to transcript id in ensembl, also get whether the transcript is canonical
#this is not a 1-1 relationship, (1-many?)
map_no_version = getBM(
  attributes = c("refseq_mrna", "ensembl_transcript_id", "transcript_is_canonical", "hgnc_symbol"),
  filters    = "refseq_mrna",
  values     = unique(v_ch$tx_id1_noversion),
  mart       = ensembl
)
#map_no_version2 = getBM(
#  attributes = c("refseq_mrna", "ensembl_transcript_id", "transcript_is_canonical", "hgnc_symbol"),
#  filters    = "refseq_mrna",
#  values     = unique(v_ch2$tx_id1_noversion),
#  mart       = ensembl
#)
v_ch$tx_id2 = map_no_version$ensembl_transcript_id[match(v_ch$tx_id1_noversion, map_no_version$refseq_mrna)]
#v_ch2$tx_id2 = map_no_version2$ensembl_transcript_id[match(v_ch2$tx_id1_noversion, map_no_version2$refseq_mrna)]
v_ch$tx_canonical = ifelse(map_no_version$transcript_is_canonical[match(v_ch$tx_id1_noversion, map_no_version$refseq_mrna)] == 1, "canonical", "non_canonical")
#v_ch2$tx_canonical = ifelse(map_no_version2$transcript_is_canonical[match(v_ch2$tx_id1_noversion, map_no_version2$refseq_mrna)] == 1, "canonical", "non_canonical")
#rename the NA into 0
v_ch$tx_canonical[is.na(v_ch$tx_canonical)] = "non_canonical"
v_ch[which(v_ch$`#AlleleID` == '2831360'),]
#v_ch2$tx_canonical[is.na(v_ch2$tx_canonical)] = "non_canonical"
#table(v_ch2$tx_canonical)
#create key in variant_summary file
v_ch$key <- paste0(
  v_ch$Chromosome, ":",
  sprintf("%09d", v_ch$PositionVCF), "|",
  v_ch$ReferenceAlleleVCF, "|",
  v_ch$AlternateAlleleVCF
)

variant_summary$key = paste0(
  variant_summary$Chromosome, ":",
  sprintf("%09d", variant_summary$PositionVCF), "|",
  variant_summary$ReferenceAlleleVCF, "|",
  variant_summary$AlternateAlleleVCF
)

#remove v_cs,v_gr, v_re to save memory
rm(v_cs, v_gr, v_re)
#output v_ch
write.csv(v_ch,'v_ch20260201.csv',row.names = F)
v_ch = read.csv('v_ch20260201.csv')
#return this information to res, build res3
#the problem is, not all tr in res appeared in this variant_summary
#res3 = res[which(res@ranges@NAMES %in% v_ch$Submitted_variation_ID)]
snv_plp_ptc_nmdesc_can$NumberSubmitters = v_ch$NumberSubmitters[match(snv_plp_ptc_nmdesc_can@elementMetadata@listData[["key"]], v_ch$key)]
snv_plp_ptc_nmdesc_can$NumberSubmitters = variant_summary$NumberSubmitters[match(snv_plp_ptc_nmdesc_can@elementMetadata@listData[["key"]], variant_summary$key)]
variant_summary_match = variant_summary[match(snv_plp_ptc_nmdesc_can@elementMetadata@listData[["key"]], variant_summary$key),]


table(snv_plp_ptc_nmdesc_can$NumberSubmitters,useNA = "ifany")

tx_n$canonical_by_vs = v_ch$tx_canonical[match(tx_n$gene, v_ch$key)]

saveRDS(snv_plp_ptc_nmdesc_can,'snv_plp_ptc_nmdesc_can20260201.rds')
saveRDS(snv_benign_ptc_nmdesc_can,'snv_benign_ptc_nmdesc_can20260201.rds')
#filter for the variants in v_ch2
snv_plp_ptc_nmdesc_can_filtered = snv_plp_ptc_nmdesc_can[which(snv_plp_ptc_nmdesc_can@elementMetadata@listData[["key"]] %in% v_ch2$key),]
snv_benign_ptc_nmdesc_can_filtered = snv_benign_ptc_nmdesc_can[which(snv_benign_ptc_nmdesc_can@elementMetadata@listData[["key"]] %in% v_ch$key),]
saveRDS(snv_plp_ptc_nmdesc_can_filtered,'snv_plp_ptc_nmdesc_can_filtered20260201.rds')
saveRDS(snv_benign_ptc_nmdesc_can_filtered,'snv_benign_ptc_nmdesc_can_filtered20260201.rds')
get_pvalue('snv_plp_ptc_nmdesc_can_filtered20260201.rds','snv_plp_ptc_nmdesc_can_p_f_syn_20260201.rds')
get_pvalue_wald('snv_plp_ptc_nmdesc_can_filtered20260201.rds','snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201.rds')
res_wald_p = readRDS('snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201.rds')
res_p1 = readRDS('snv_plp_ptc_nmdesc_can_filtered20260201.rds')
res_p_syn = readRDS('snv_plp_ptc_nmdesc_can_p_f_syn_20260201.rds')
p_set = NULL
for(i in 1:790){
  p_set = rbind(p_set,(res_p_syn[[i]][["can.pvalue"]]))
}
write.csv(p_set,'p_less.csv',row.names = F)

#change p cutoff
get_NMD_enrichment_wald('snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201.rds',FDR = 0.05,filter_type = 'can')
ppi_AD_genes <- read_csv("~/Desktop/frameshift/ppi_AD_genes.csv")
old_snv_AD = ppi_AD_genes %>% filter(group =='Nonsense') %>% dplyr::select(hgnc_symbol) 
length(unique(old_snv_AD$hgnc_symbol))
snv_can_gene = read.csv("snv_plp_ptc_nmdesc_can_p_f_syn_20260201_NMDesc_enriched_can.txt",header=F)
#filter for AD genes
omim_AD_symbols = read.csv('omim_AD_symbols.csv',header=F)$V1
wald_genes = read.csv('snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201_NMDesc_wald_enriched_can.txt',header=F)$V1
binom_genes = read.csv('snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201_NMDesc_binom_enriched_can.txt',header=F)$V1
wald_AD_genes = wald_genes[wald_genes %in% omim_AD_symbols$x]
binom_AD_genes = binom_genes[binom_genes %in% omim_AD_symbols$x]
snv_gene = gene_all$hgnc_symbol[gene_all$group == 'snv']
fs_gene = gene_all$hgnc_symbol[gene_all$group == 'fs']
length(wald_AD_genes)
length(binom_AD_genes)
length(snv_gene)
length(intersect(wald_AD_genes, binom_AD_genes))
length(intersect(wald_AD_genes, snv_gene))
length(intersect(binom_AD_genes, snv_gene))

snv_AD_can_gene = snv_can_gene %>% filter(V1 %in% omim_AD_symbols)
write.csv(snv_AD_can_gene,'snv_plp_ptc_nmdesc_can_p_f_syn_20260201_NMDesc_enriched_can_AD_p_0.8.csv',row.names = F)
length(intersect(old_snv_AD$hgnc_symbol, snv_can_gene$V1))


get_snv_variant_new.R
get_fs_variant_new.R
get_gnomad_control.R

build_gene_all.R
-------------------------
  #gene_level comparision
gene_all = read.csv('gene_all.csv')

calculate_ppi_degree_centrality(
  gene_all,
  output_csv = "wald_ppi_degree_centrality_results.csv"
)

plot_gc_content(
  gene_all    = gene_all,
  output_csv  = "gc_content.csv",
  output_fig  = "gc_content.png"
)

plot_repeat_content(
  gene_all   = gene_all,
  output_csv = "repeat_content.csv",
  output_fig = "repeat_content.png"
)

annotate_motif_flags(
  gene_all         = gene_all,
  path_touni       = "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/others/NIHMS1818854-supplement-2(A).csv",
  path_motif       = "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/others/NIHMS1818854-supplement-2(B).csv",
  path_LCS         = "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/others/Copy of NIHMS1818854-supplement-2.xls",
  mart             = ensembl,
  output_motif_csv = "gene_motif_flags.csv",
  output_lcs_csv   = "gene_LCS_flags.csv"
)

run_pfam_overlap_analysis(
  gene_all      = gene_all,
  ensembl       = ensembl,
  output_prefix = "pfam_overlap"   # 输出文件前缀，可自定义
)

run_ppi_overlap_analysis(
  gene_all      = gene_all,
  ppi_file_path = "~/Downloads/human (1) (2).txt",
  output_prefix = "ppi_overlap"
)

run_tau_analysis(
  gene_all      = gene_all,
  gtex_path     = "/Users/jxu14/Downloads/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct",
  output_prefix = "tau"
)

plot_gene_level_features(
  gene_all = gene_all,
  lof_metrics_path = "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/others/gnomad.v2.1.1.lof_metrics.by_gene.txt",
  ensembl = NULL,
  out_dir = ".",
  prefix = "gene_level"
)

#do a prelim regression to see the relationship
wald_ppi_degree_centrality <- read_csv("wald_ppi_degree_centrality_results.csv")
gc_content <- read_csv("gc_content.csv")
repeat_content <- read_csv("repeat_content.csv")
gene_motif_flags <- read_csv("gene_motif_flags.csv")
gene_LCS_flags   <- read_csv("gene_LCS_flags.csv")
pfam_overlap   <- read_csv("pfam_overlap_gene_all.csv")
ppi_overlap  <- read_csv("ppi_overlap_gene_all.csv")
tau_results  <- read_csv("tau_gene_matrix.csv")
gene_level <- read_csv("gene_level_pli_loeuf_category.csv")

#merge them to one dataframe
gene_all_merged <- gene_all %>%
  left_join(wald_ppi_degree_centrality %>% select(hgnc_symbol, group, Degree),
            by = c("hgnc_symbol", "group")) %>%
  left_join(gc_content %>% select(ensembl_transcript_id, group, gc_content, nmdesc_gc_content),
            by = c("ensembl_transcript_id", "group")) %>%
  left_join(repeat_content %>% select(ensembl_transcript_id, group, repeat_fraction, nmdesc_repeat_fraction,
                                            homopolymer_fraction, nmdesc_homopolymer_fraction),
            by = c("ensembl_transcript_id", "group")) %>%
  left_join(gene_motif_flags %>% select(ensembl_transcript_id, group, gene_protein_flag, gene_domains_flag,
                                              gene_slim_flag, gene_morf_flag, gene_ptm_flag, gene_nls_flag),
            by = c("ensembl_transcript_id", "group")) %>%
  left_join(gene_LCS_flags %>% select(ensembl_transcript_id, group, gene_LCS_flag),
            by = c("ensembl_transcript_id", "group")) %>%
  left_join(pfam_overlap %>% select(ensembl_transcript_id, group, pfam_overlap_length,
                                          pfam_overlap_flag, pfam_overlap_fraction, n_overlapping_pfam),
            by = c("ensembl_transcript_id", "group")) %>%
  left_join(ppi_overlap %>% select(ensembl_transcript_id, group, ppi_overlap),
            by = c("ensembl_transcript_id", "group")) %>%
  left_join(gene_level %>% select(hgnc_symbol, group, pLI, oe_lof_upper, pli_cat, loeuf_cat),
            by = c("hgnc_symbol", "group"))

#do regression
model_data <- gene_all_merged %>%
  mutate(is_nmdesc = if_else(group %in%  c("fs_control","snv_control"), 0L, 1L)) 

write.csv(model_data,'gene_model_data.csv',row.names = F)
#ppi degree centrality
model1_snv <- glm(is_nmdesc ~ Degree, data = model_data[which(model_data$group %in% c('snv','snv_control')),], family = binomial)
model2_snv <- glm(is_nmdesc ~ cds_length + NMDesc_region_length + Degree, data = model_data[which(model_data$group %in% c('snv','snv_control')),], family = binomial)
model1_fs <- glm(is_nmdesc ~ Degree, data = model_data[which(model_data$group %in% c('fs','fs_control')),], family = binomial)
model2_fs <- glm(is_nmdesc ~ cds_length + NMDesc_region_length + Degree, data = model_data[which(model_data$group %in% c('fs','fs_control')),], family = binomial)
# Results
tidy(model1_snv, exponentiate = TRUE, conf.int = TRUE) %>%
  arrange(p.value)
tidy(model2_snv, exponentiate = TRUE, conf.int = TRUE) %>%
  arrange(p.value)
tidy(model1_fs, exponentiate = TRUE, conf.int = TRUE) %>%
  arrange(p.value)
tidy(model2_fs, exponentiate = TRUE, conf.int = TRUE) %>%
  arrange(p.value)

-----------------
  #variant_level comparison
  snv_variants = read.csv('snv_variants20260201_plp_wald_clinvar.csv')
  fs_variants = read.csv('fs_variants20260201_plp_wald_clinvar.csv')
  snv_control_variants = read.csv('gnomad_snv_filtered_wald.csv')
  fs_control_variants = read.csv('gnomad_fs_filtered_wald.csv')
  
  variant_data = bind_rows(
    snv_variants %>% mutate(group = 'snv'),
    fs_variants %>% mutate(group = 'fs'),
    snv_control_variants %>% mutate(group = 'snv_control'),
    fs_control_variants %>% mutate(group = 'fs_control')
  )
  
  #add attributes
  new_create_fasta.R #add cds_mutation_loc etc
  combine_variant.R #combine variant level info with gene level info
  
  #do a regression on dist to cds end using variants_all2
  variant_data2 = variants_all2
  variant_data2 = variant_data2 %>%
    mutate(is_nmdesc = if_else(group %in%  c("fs_control","snv_control"), 0L, 1L))
  model_dist_snv <- glm(is_nmdesc ~ dist_to_cds_end, data = variant_data2[which(variant_data2$group %in% c('snv_disease','snv_control')),], family = binomial)
  model_dist_fs <- glm(is_nmdesc ~ dist_to_cds_end, data = variant_data2[which(variant_data2$group %in% c('fs_disease','fs_control')),], family = binomial)
  tidy(model_dist_snv, exponentiate = TRUE, conf.int = TRUE) %>%
    arrange(p.value)
  tidy(model_dist_fs, exponentiate = TRUE, conf.int = TRUE) %>%
    arrange(p.value)
  #for model2 add cds_end
  model_dist2_snv <- glm(is_nmdesc ~ dist_to_cds_end + cds_end, data = variant_data2[which(variant_data2$group %in% c('snv_disease','snv_control')),], family = binomial)
  model_dist2_fs <- glm(is_nmdesc ~ dist_to_cds_end + cds_end, data = variant_data2[which(variant_data2$group %in% c('fs_disease','fs_control')),], family = binomial)
  tidy(model_dist2_snv, exponentiate = TRUE, conf.int = TRUE) %>%
    arrange(p.value)
  tidy(model_dist2_fs, exponentiate = TRUE, conf.int = TRUE) %>%
    arrange(p.value)
  
  #add variant level pfam and ppi flags
  pfam_raw <- getBM(
    attributes = c(
      "ensembl_transcript_id",
      "hgnc_symbol",
      "pfam",
      "pfam_start",
      "pfam_end",
      "uniprotswissprot"
    ),
    filters = "ensembl_transcript_id",
    values  = unique(variants_all2$ensembl_transcript_id),
    mart    = ensembl
  )
  
  pfam_fin <- pfam_raw %>%
    filter(!is.na(pfam), pfam != "") %>%
    distinct()
  pfam_fin$uniprot = pfam_fin$uniprotswissprot
  
  human_1_ <- read_delim("~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/others/human (1).txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
  variant_pfam_ppi(
       variants_all2 = variants_all2,
       human_1_      = human_1_,
       pfam_fin      = pfam_fin,
       ensembl       = ensembl,
       out_dir       = "plots/pfam_ppi_analysis"
   )
  
  
  ------------
    #control for covariates
    
    #control for cds length
    
    #control for nmdesc region 
    
    #control for GC content
    
    #control for transcript(variant level)
    
    
  
