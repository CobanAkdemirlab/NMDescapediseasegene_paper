
.p <- c("gene level_v3/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../gene level_v3/lib/paths.R", "lib/paths.R")
.p <- .p[file.exists(.p)]
source(.p[1])
# ------------------------------------------------------------

# --- load functions-----------------------------------------------
.fn_dir <- c("gene level_v3/features/functions", "../../features/functions",
             "../features/functions", "features/functions")
.fn_dir <- .fn_dir[dir.exists(.fn_dir)]
for (.f in list.files(.fn_dir[1], pattern = "\\.R$", full.names = TRUE)) source(.f)
rm(.f, .fn_dir)

for (.dep in c("get_pvalue.R", "extract_enriched.R", "process_syn.R")) {
  .cand <- c(file.path("gene level_v3/disease genes/snv", .dep), .dep,
             file.path("../snv", .dep))
  .cand <- .cand[file.exists(.cand)]
  if (length(.cand)) source(.cand[1]) else message("  can't find ", .dep)
}
rm(.dep, .cand)
# --------------------------------------------------------------

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
library(ggplot2)

#1. NMDesc annotation
# setwd("/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/clinvar")
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl = mart
genome <- BSgenome.Hsapiens.UCSC.hg38
txdb <- txdbmaker::makeTxDbFromEnsembl("Homo sapiens", release=105)
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
all.df.sub <- all.df[which(all.df$TXNAME!='NA'),]
rm(all.df)
BM.info <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","hgnc_symbol","transcript_is_canonical"),mart=ensembl)
#write.csv(BM.info,'BM_info.csv',row.names = F)

vcf_file = data_file("clinvar_20260201.vcf.gz")
vcf = aenmd:::parse_vcf_VariantAnnotation(vcf_file)
vcf_rng = vcf$vcf_rng
#add clinical significance
clnsig_list <- info(vcf$vcf_obj)$CLNSIG 
rm(vcf)
mcols(vcf_rng)$CLNSIG <- clnsig_list
vcf_rng_fil = process_variants(vcf_rng)
#- filter out variants with ill-defined alternative allele
ind_out =  Biostrings::vcountPattern("N", vcf_rng_fil$alt) > 0
vcf_rng_fil = vcf_rng_fil[!ind_out]
#- back to the original workflow
res = annotate_nmd(vcf_rng_fil, rettype="gr")


#2. get nmdesc enriched genes
snv_ind = which(res@elementMetadata@listData[["type"]] == 'snv')
fs_ind = which(res@elementMetadata@listData[["type"]] %in% c('del','ins'))
clnsig_str = sapply(res$CLNSIG, function(x) paste(as.character(x), collapse="|"))
plp_ind = which(grepl("pathogenic", clnsig_str, ignore.case = TRUE))
length(unique(res@elementMetadata@listData[["key"]][plp_ind]))

benign_ind = which(grepl("benign", clnsig_str, ignore.case = TRUE))
ptc_ind = which(res@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]]==T)
length(unique(res@elementMetadata@listData[["key"]][intersect(plp_ind,ptc_ind)]))

nmdesc_ind = which(res@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]] == T | res@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]]==T)
length(unique(res@elementMetadata@listData[["key"]][nmdesc_ind]))
#filter for canonical transcript
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
can_ind = which(res@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] %in% tx_can$ensembl_transcript_id)
#combine the filters
length(unique(res@elementMetadata@listData[["key"]][plp_ind]))

plp_ptc_ind = intersect(plp_ind,ptc_ind)
length(unique(res@elementMetadata@listData[["key"]][plp_ptc_ind]))
plp_ptc_can_ind = intersect(plp_ptc_ind,can_ind)
length(unique(res@elementMetadata@listData[["key"]][plp_ptc_can_ind]))
plp_ptc_nmdesc_can_ind = intersect(plp_ptc_can_ind,nmdesc_ind)
benign_ptc_ind = intersect(benign_ind,ptc_ind)
snv_plp_ptc_ind = intersect(snv_ind,plp_ptc_ind)
snv_plp_ptc = res[snv_plp_ptc_ind]
snv_plp_ptc_can_ind = intersect(snv_plp_ptc_ind,can_ind)
snv_plp_ptc_can = res[snv_plp_ptc_can_ind]
snv_benign_ptc_ind = intersect(snv_ind,benign_ptc_ind)
snv_plp_ptc_nmdesc_ind = intersect(snv_plp_ptc_ind,nmdesc_ind)
snv_benign_ptc_nmdesc_ind = intersect(snv_benign_ptc_ind,nmdesc_ind)
snv_plp_ptc_nmdesc_can_ind = intersect(snv_plp_ptc_nmdesc_ind,can_ind)
fs_plp_ptc_nmdesc_can_ind = intersect(fs_ind,plp_ptc_nmdesc_can_ind)
snv_benign_ptc_nmdesc_can_ind = intersect(snv_benign_ptc_nmdesc_ind,can_ind)

snv_plp_ptc_nmdesc_can = res[snv_plp_ptc_nmdesc_can_ind]
fs_plp_ptc_nmdesc_can = res[fs_plp_ptc_nmdesc_can_ind]
snv_benign_ptc_nmdesc_can = res[snv_benign_ptc_nmdesc_can_ind]
plp_ptc_nmdesc_can = res[plp_ptc_nmdesc_can_ind]
length(unique(plp_ptc_nmdesc_can@elementMetadata@listData[["key"]]))
length(unique(snv_plp_ptc_nmdesc_can@ranges@NAMES))
length(unique(snv_plp_ptc_nmdesc_can@elementMetadata@listData[["key"]]))

#Quality Control, match by key
variant_summary <- read_delim("variant_summary.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
#ReviewStatus does not contain either of “no assertion” or “no interpretation”
#variant_summary$ReviewStatus
v_re = variant_summary %>%
  filter(
    !str_detect(tolower(ReviewStatus), "no assertion"),
    !str_detect(tolower(ReviewStatus), "no interpretation")
  )
v_re[which(v_re$`#AlleleID` == '2831360'),]

#ClinicalSignificance does not contain any of “not provided”, “drug
#response”, “other”, “risk”, “low penetrance”, “conflicting”, “affects”, “association”, “protective”, “confers
#sensitivity”
exclude_terms <- c(
  "not provided", "drug response", "other", "risk", "low penetrance",
  "conflicting", "affects", "association", "protective", "confers sensitivity"
)
pattern_exclude <- paste(exclude_terms, collapse = "|")
v_cs <- v_re %>%
  filter(!str_detect(ClinicalSignificance, regex(pattern_exclude, ignore_case = TRUE)))
v_cs[which(v_cs$`#AlleleID` == '2831360'),]

#Assembly == GRCh38
v_gr = v_cs %>%
  filter(str_trim(tolower(Assembly)) == "grch38")
v_gr[which(v_gr$`#AlleleID` == '2831360'),]

#autosomal contigs (chr1-22) only
v_ch = v_gr %>%
  mutate(
    Chromosome_clean = str_remove(tolower(str_trim(Chromosome)), "^chr")
  ) %>%
  filter(Chromosome_clean %in% as.character(1:22)) %>%
  dplyr::select(-Chromosome_clean)
v_ch[which(v_ch$`#AlleleID` == '2831360'),]

#add canonical and non_canonical transcript tag use Name
#substract from Name(eg. NM_017547.4 from NM_017547.4(FOXRED1):c.694C>T (p.Gln232Ter))
v_ch$tx_id1 = sub("\\(.*", "", v_ch$Name)
#v_ch2$tx_id1 = sub("\\(.*", "", v_ch2$Name)
v_ch = v_ch %>%
  mutate(
    tx_id1_noversion = sub("\\..*$", "", tx_id1)
  )

#presence of a valid RefSeq ID
v_ch = v_ch %>%
  filter(str_detect(tx_id1_noversion, "^NM_"))
v_ch[which(v_ch$`#AlleleID` == '2831360'),]

#transfer to transcript id in ensembl, also get whether the transcript is canonical
#this is not a 1-1 relationship, (1-many?)
map_no_version = getBM(
  attributes = c("refseq_mrna", "ensembl_transcript_id", "transcript_is_canonical", "hgnc_symbol"),
  filters    = "refseq_mrna",
  values     = unique(v_ch$tx_id1_noversion),
  mart       = ensembl
)

v_ch$tx_id2 = map_no_version$ensembl_transcript_id[match(v_ch$tx_id1_noversion, map_no_version$refseq_mrna)]
v_ch$tx_canonical = ifelse(map_no_version$transcript_is_canonical[match(v_ch$tx_id1_noversion, map_no_version$refseq_mrna)] == 1, "canonical", "non_canonical")
#rename the NA into 0
v_ch$tx_canonical[is.na(v_ch$tx_canonical)] = "non_canonical"
v_ch[which(v_ch$`#AlleleID` == '2831360'),]

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
res2 = res[which(res@elementMetadata@listData[["key"]] %in% v_ch$key)]
snv_plp_ptc_nmdesc_can$NumberSubmitters = v_ch$NumberSubmitters[match(snv_plp_ptc_nmdesc_can@elementMetadata@listData[["key"]], v_ch$key)]
snv_plp_ptc_nmdesc_can$NumberSubmitters = variant_summary$NumberSubmitters[match(snv_plp_ptc_nmdesc_can@elementMetadata@listData[["key"]], variant_summary$key)]
variant_summary_match = variant_summary[match(snv_plp_ptc_nmdesc_can@elementMetadata@listData[["key"]], variant_summary$key),]

table(snv_plp_ptc_nmdesc_can$NumberSubmitters,useNA = "ifany")

tx_n$canonical_by_vs = v_ch$tx_canonical[match(tx_n$gene, v_ch$key)]

saveRDS(snv_plp_ptc_nmdesc_can,'snv_plp_ptc_nmdesc_can20260201.rds')
saveRDS(snv_benign_ptc_nmdesc_can,'snv_benign_ptc_nmdesc_can20260201.rds')
#filter for the variants in v_ch2
snv_plp_ptc_nmdesc_can_filtered = snv_plp_ptc_nmdesc_can[which(snv_plp_ptc_nmdesc_can@elementMetadata@listData[["key"]] %in% v_ch$key),]
snv_plp_ptc_can_filtered = snv_plp_ptc_can[which(snv_plp_ptc_can@elementMetadata@listData[["key"]] %in% v_ch$key),]
saveRDS(snv_plp_ptc_can_filtered,'snv_plp_ptc_can_filtered20260201.rds')
fs_plp_ptc_nmdesc_can_filtered = fs_plp_ptc_nmdesc_can[which(fs_plp_ptc_nmdesc_can@elementMetadata@listData[["key"]] %in% v_ch$key),]
snv_benign_ptc_nmdesc_can_filtered = snv_benign_ptc_nmdesc_can[which(snv_benign_ptc_nmdesc_can@elementMetadata@listData[["key"]] %in% v_ch$key),]
saveRDS(snv_plp_ptc_nmdesc_can_filtered,'snv_plp_ptc_nmdesc_can_filtered20260201_check.rds')
temp1 = readRDS('snv_plp_ptc_nmdesc_can_filtered20260201.rds')
temp2 = readRDS('snv_plp_ptc_nmdesc_can20260201.rds')
saveRDS(snv_benign_ptc_nmdesc_can_filtered,'snv_benign_ptc_nmdesc_can_filtered20260201.rds')
get_pvalue('snv_plp_ptc_can_filtered20260201.rds',
                          'snv_plp_ptc_nmdesc_can_filtered20260201.rds',
                           'snv_plp_ptc_nmdesc_can_p_f_syn_20260201_AD_BH_FDR020.rds',
                            restrict_symbols = omim_AD_symbols)
                  
txnames.list <- readRDS('snv_plp_ptc_nmdesc_can_p_f_syn_20260201_AD_FDR020.rds')
rest.all <- sapply(txnames.list, function(x) if(is.null(x$rest.PTC)) NA else x$rest.PTC)
summary(rest.all)
quantile(rest.all, c(0.1, 0.25, 0.5, 0.75), na.rm = TRUE)
res3 <- extract_enriched(txnames.list, fdr.method = "dbh")
res3$genes                    # 富集基因名向量
head(res3$enriched, 20)       # 按 p 值排序的转录本明细
write.csv(res3$enriched, "enriched_genes_dbh.csv", row.names = FALSE)

# 三种方法的基因集对比
g.dbh <- extract_enriched(txnames.list, "dbh")$genes
g.dby <- extract_enriched(txnames.list, "dby")$genes
g.bh  <- extract_enriched(txnames.list, "bh")$genes
length(g.dbh); length(g.dby); length(g.bh)
length(intersect(g.dbh, g.bh))

re1 = readRDS('snv_plp_ptc_nmdesc_can_p_f_syn_20260201_Jul30.rds')
get_pvalue_wald('snv_plp_ptc_nmdesc_can_filtered20260201_check.rds','snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201_Jul30.rds')
res_wald_p = readRDS('snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201.rds')
res_p1 = readRDS('snv_plp_ptc_nmdesc_can_filtered20260201.rds')
res_p_syn = readRDS('snv_plp_ptc_nmdesc_can_p_f_syn_20260201.rds')
p_set = NULL
for(i in 1:790){
  p_set = rbind(p_set,(res_p_syn[[i]][["can.pvalue"]]))
}
write.csv(p_set,'p_less.csv',row.names = F)

#get enriched genes
# get_NMD_enrichment_wald('snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201.rds',FDR = 0.05,filter_type = 'can')

.e <- c("gene level_v3/disease genes/snv/get_NMD_enrichment_DBH.R",
        "get_NMD_enrichment_DBH.R", "../snv/get_NMD_enrichment_DBH.R")
.e <- .e[file.exists(.e)]
if (!length(.e)) stop("找不到 get_NMD_enrichment_DBH.R")
source(.e[1]); rm(.e)
get_NMD_enrichment_DBH()
snv_can_gene = read.csv("snv_plp_ptc_nmdesc_can_p_f_syn_20260201_NMDesc_enriched_can.txt",header=F)

#filter for AD genes
omim_AD_symbols = read.csv('omim_AD_symbols.csv',header=F)$V1
dbh_genes = read.csv('snv_plp_ptc_nmdesc_can_p_f_syn_20260201_Jul30_NMDesc_dbh_enriched_can.txt',header=F)$V1
wald_genes = read.csv('snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201_NMDesc_wald_enriched_can.txt',header=F)$V1
binom_genes = read.csv('snv_plp_ptc_nmdesc_can_wald_p_f_syn_20260201_NMDesc_binom_enriched_can.txt',header=F)$V1
wald_AD_genes = wald_genes[wald_genes %in% omim_AD_symbols[-1]]
dbh_AD_genes = dbh_genes[dbh_genes %in% omim_AD_symbols[-1]]
write.csv(wald_AD_genes,'snv_wald_AD.csv',row.names = F)
write.csv(dbh_AD_genes,'snv_dbh_AD.csv',row.names = F)
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


#get the total number of submitters for each gene
gene_all$NumberSubmitters = v_ch$NumberSubmitters[match(gene_all$k, v_ch$key)]
submit_info = v_ch[,c('NumberSubmitters','key','tx_id2','tx_canonical')]
write.csv(submit_info,'submit_info.csv',row.names = F)
gene_all %>% group_by(group) %>% summarise(mean_submitters = mean(NumberSubmitters, na.rm = TRUE), median_submitters = median(NumberSubmitters, na.rm = TRUE), low_submitters = sum(NumberSubmitters <= 1, na.rm = TRUE))
#remove the low submitters, which may be more likely to be false positives
gene_all %>% filter(NumberSubmitters > 1) %>% group_by(group) %>% summarise(mean_submitters = mean(NumberSubmitters, na.rm = TRUE), median_submitters = median(NumberSubmitters, na.rm = TRUE), low_submitters = sum(NumberSubmitters <= 1, na.rm = TRUE))
table(gene_all$group)

gene_can_AD_plp_ptc_nmdesc_uni$NumberSubmitters = variant_summary$NumberSubmitters[match(gene_can_AD_plp_ptc_nmdesc_uni$key, variant_summary$key)]
summary(gene_can_AD_plp_ptc_nmdesc_uni$NumberSubmitters)


get_snv_variant_new.R
get_fs_variant_new.R
get_gnomad_control.R

build_gene_all.R #(adapted from combine_gene.R)
-------------------------
#3. gene_level comparision
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
  path_touni       = data_file("NIHMS1818854-supplement-2(A).csv"),
  path_motif       = data_file("NIHMS1818854-supplement-2(B).csv"),
  path_LCS         = data_file("Copy of NIHMS1818854-supplement-2.xls"),
  mart             = ensembl,
  output_motif_csv = "gene_motif_flags.csv",
  output_lcs_csv   = "gene_LCS_flags.csv"
)

run_pfam_overlap_analysis(
  gene_all      = gene_all,
  ensembl       = ensembl,
  output_prefix = "pfam_overlap"   
)

run_ppi_overlap_analysis(
  gene_all      = gene_all,
  ppi_file_path = data_file("human (1).txt"),
  output_prefix = "ppi_overlap"
)

run_tau_analysis(
  gene_all      = gene_all,
  gtex_path     = data_file("GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct",
                                must = FALSE),   
  output_prefix = "tau"
)

plot_gene_level_features(
  gene_all = gene_all,
  lof_metrics_path = data_file("gnomad.v2.1.1.lof_metrics.by_gene.txt"),
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

#do multiple correction

#do regression
model_data <- gene_all_merged %>%
  mutate(is_nmdesc = if_else(group %in%  c("fs_control","snv_control"), 0L, 1L)) 

write.csv(model_data,'gene_model_data.csv',row.names = F)
model_data <- read_csv('gene_model_data.csv')
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


#clean data
model_data <- model_data %>%
  mutate(
    is_disease = as.integer(group %in% c("snv", "fs")),
    gene_set   = case_when(
      group %in% c("snv", "snv_control") ~ "SNV",
      group %in% c("fs",  "fs_control")  ~ "FS"
    )
  )

flag_cols <- c("gene_protein_flag", "gene_domains_flag", "gene_slim_flag",
               "gene_morf_flag", "gene_ptm_flag", "gene_nls_flag",
               "gene_LCS_flag", "pfam_overlap_flag", "ppi_overlap")

#confounders <- "cds_length + NMDesc_region_length"
confounders <- "cds_length + NMDesc_region_length + gc_content"

#check if in each feature they have enough data points
check_flags <- function(data, gene_set_label) {
  data %>%
    filter(gene_set == gene_set_label) %>%
    select(all_of(flag_cols)) %>%
    mutate(across(everything(), as.integer)) %>%
    summarise(across(everything(), ~ n_distinct(na.omit(.)))) %>%
    pivot_longer(everything(), names_to = "flag", values_to = "n_unique") %>%
    filter(n_unique < 2) %>%
    pull(flag)
}

skip_snv <- check_flags(model_data, "SNV")
skip_fs  <- check_flags(model_data, "FS")

message("Skipping in SNV: ", if (length(skip_snv)) paste(skip_snv, collapse = ", ") else "none")
message("Skipping in FS: ",  if (length(skip_fs))  paste(skip_fs,  collapse = ", ") else "none")

#get OR
extract_or <- function(fit, coef_name) {
  ci <- tryCatch(confint(fit)[coef_name, ], error = function(e) c(NA, NA))
  data.frame(
    flag      = coef_name,
    OR        = exp(coef(fit)[coef_name]),
    OR_low    = exp(ci[1]),
    OR_high   = exp(ci[2]),
    p_value   = summary(fit)$coefficients[coef_name, "Pr(>|z|)"],
    row.names = NULL
  )
}

run_unadjusted <- function(data, gene_set_label, skip_flags = character(0)) {
  
  active_flags <- setdiff(flag_cols, skip_flags)
  
  lapply(active_flags, function(col) {
    df <- data %>%
      select(is_disease, all_of(col)) %>%
      filter(!is.na(.data[[col]])) %>%
      mutate(across(all_of(col), as.integer))
    
    if (nrow(df) == 0 ||
        length(unique(df[[col]])) < 2 ||
        length(unique(df$is_disease)) < 2) {
      message(sprintf("Skipping %s in %s: no variation", col, gene_set_label))
      return(NULL)
    }
    
    fit <- tryCatch(
      glm(as.formula(paste("is_disease ~", col)), data = df, family = binomial),
      error = function(e) {
        message(sprintf("Model failed for %s in %s: %s", col, gene_set_label, e$message))
        NULL
      }
    )
    if (is.null(fit)) return(NULL)
    extract_or(fit, col)
  }) %>%
    bind_rows() %>%
    mutate(gene_set = gene_set_label, model = "Unadjusted")
}

run_adjusted <- function(data, gene_set_label, skip_flags = character(0)) {
  
  active_flags <- setdiff(flag_cols, skip_flags)
  
  lapply(active_flags, function(col) {
    df <- data %>%
      select(is_disease, all_of(col), cds_length, NMDesc_region_length, gc_content) %>%
      filter(!is.na(.data[[col]]),
             !is.na(cds_length),
             !is.na(NMDesc_region_length),
             !is.na(gc_content)) %>%
      mutate(across(all_of(col), as.integer))
    
    if (nrow(df) == 0 ||
        length(unique(df[[col]])) < 2 ||
        length(unique(df$is_disease)) < 2) {
      message(sprintf("Skipping %s in %s: no variation", col, gene_set_label))
      return(NULL)
    }
    
    fit <- tryCatch(
      glm(as.formula(paste("is_disease ~", col, "+", confounders)),
          data = df, family = binomial),
      error = function(e) {
        message(sprintf("Model failed for %s in %s: %s", col, gene_set_label, e$message))
        NULL
      }
    )
    if (is.null(fit)) return(NULL)
    extract_or(fit, col)
  }) %>%
    bind_rows() %>%
    mutate(gene_set = gene_set_label, model = "Adjusted")
}

unadj_snv <- run_unadjusted(filter(model_data, gene_set == "SNV"), "SNV", skip_snv)
unadj_fs  <- run_unadjusted(filter(model_data, gene_set == "FS"),  "FS",  skip_fs)
adj_snv   <- run_adjusted(filter(model_data, gene_set == "SNV"),   "SNV", skip_snv)
adj_fs    <- run_adjusted(filter(model_data, gene_set == "FS"),    "FS",  skip_fs)

combined <- bind_rows(unadj_snv, unadj_fs, adj_snv, adj_fs) %>%
  group_by(gene_set, model) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    sig   = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  ) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = c("Unadjusted", "Adjusted")))


ggplot(combined, aes(x = OR, y = reorder(flag, OR), color = sig, shape = gene_set)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.2,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  facet_wrap(~ model) +
  scale_color_manual(values = c("***" = "red", "**" = "orange",
                                "*"   = "gold", "ns" = "grey60")) +
  scale_shape_manual(values = c("SNV" = 16, "FS" = 17)) +
  scale_x_log10() +
  labs(
    title    = "Gene-level feature association with disease (model_data)",
    subtitle = "Adjusted for cds_length + NMDesc_region_length + GC content",
    x        = "Odds Ratio (log scale)", y = NULL,
    color    = "Adj. p", shape = "Gene set"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))


-----------------
#4. variant_level comparison
##the variant level comparision has 4 parts: unmatched, mixed-effect model, hierarchical baysian model,bootstrap and matched analysis
  # setwd('/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/clinvar/')
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
##4.1 unmatched analysis
  
##4.2 mixed-effect model
  
##4.3 hierarchial mode
  
##4.4 bootstrap and matched analysis
  
  
  
  #do a regression on dist to cds end using variants_all
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
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
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
    values  = unique(variants_all5_1_0621$ensembl_transcript_id),
    mart    = ensembl
  )
  
  pfam_fin <- pfam_raw %>%
    filter(!is.na(pfam), pfam != "") %>%
    distinct()
  pfam_fin$uniprot = pfam_fin$uniprotswissprot
  
  human_1_ <- read_delim(data_file("human (1).txt"), 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
  variants_all3_PDB = variant_pfam_ppi(
       variants_all2 = variants_all5,
       human_1_      = human_1_[which(human_1_$source == 'PDB'),],
       pfam_fin      = pfam_fin,
       ensembl       = ensembl,
       out_dir       = "plots/pfam_ppi_analysis_PDB"
   )
  variants_all3_PIONEER = variant_pfam_ppi(
    variants_all2 = variants_all5,
    human_1_      = human_1_[which(human_1_$source == 'PIONEER'),],
    pfam_fin      = pfam_fin,
    ensembl       = ensembl,
    out_dir       = "plots/pfam_ppi_analysis_PIONEER"
  )
  variants_all3_HM = variant_pfam_ppi(
    variants_all2 = variants_all5,
    human_1_      = human_1_[which(human_1_$source == 'HM'),],
    pfam_fin      = pfam_fin,
    ensembl       = ensembl,
    out_dir       = "plots/pfam_ppi_analysis_HM"
  )
  
  str(variants_all3_PDB)
  
#define functions
  
  # clean duplicated rows
  clean_xy_cols <- function(data) {
    data %>%
      select(-ends_with(".x")) %>%
      rename_with(~ gsub("\\.y$", "", .), ends_with(".y"))
  }
  
  # code is_disease and gene_set
  encode_groups <- function(data) {
    data %>%
      mutate(
        is_disease = as.integer(group %in% c("snv_disease", "fs_disease")),
        gene_set   = case_when(
          group %in% c("snv_disease", "snv_control") ~ "SNV",
          group %in% c("fs_disease",  "fs_control")  ~ "FS"
        )
      )
  }
  
  # Unadjusted logistic regression
  run_unadjusted_ppi <- function(data, gene_set_label) {
    df <- data %>%
      filter(!is.na(variant_ppi_overlap))
    
    if (nrow(df) == 0 ||
        length(unique(df$variant_ppi_overlap)) < 2 ||
        length(unique(df$is_disease)) < 2) return(NULL)
    
    fit <- tryCatch(
      glm(is_disease ~ variant_ppi_overlap, data = df, family = binomial),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    
    ci <- tryCatch(confint(fit)["variant_ppi_overlap", ], error = function(e) c(NA, NA))
    data.frame(
      gene_set = gene_set_label,
      model    = "Unadjusted",
      OR       = exp(coef(fit)["variant_ppi_overlap"]),
      OR_low   = exp(ci[1]),
      OR_high  = exp(ci[2]),
      p_value  = summary(fit)$coefficients["variant_ppi_overlap", "Pr(>|z|)"]
    )
  }
  
  #run models
  run_all <- function(data, source_label) {
    data <- clean_xy_cols(data) %>% encode_groups()
    
    message(sprintf("\n══ %s ══", source_label))
    print(table(data$group))
    
    bind_rows(
      run_unadjusted_ppi(filter(data, gene_set == "SNV"), "SNV"),
      run_unadjusted_ppi(filter(data, gene_set == "FS"),  "FS")
    ) %>%
      mutate(source = source_label)
  }
  
  # ══════════════════════════════════════════════════════════════════════════════
  # run on PDB, PIONEER, and HM datasets
  # ══════════════════════════════════════════════════════════════════════════════
  results_pdb     <- run_all(variants_all3_PDB,     "PDB")
  results_pioneer <- run_all(variants_all3_PIONEER,  "PIONEER")
  results_hm      <- run_all(variants_all3_HM,       "HM")
  
  # combine results
  results_all <- bind_rows(results_pdb, results_pioneer, results_hm) %>%
    mutate(
      model  = factor(model,  levels = c("Unadjusted", "Gene-matched")),
      source = factor(source, levels = c("PDB", "HM", "PIONEER")),
      sig    = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  
  print(results_all)
  
  # ══════════════════════════════════════════════════════════════════════════════
  # PLOT
  # ══════════════════════════════════════════════════════════════════════════════
  ggplot(results_all, aes(x = OR, y = source, color = sig, shape = gene_set)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.2,
                   position = position_dodge(width = 0.5)) +
    geom_point(size = 4, position = position_dodge(width = 0.5)) +
    facet_wrap(~ model) +
    scale_color_manual(values = c("***" = "red", "**" = "orange",
                                  "*"   = "gold", "ns" = "grey60")) +
    scale_shape_manual(values = c("SNV" = 16, "FS" = 17)) +
    scale_x_log10() +
    labs(
      title    = "PPI overlap across PPI sources: Unadjusted vs Gene-matched",
      subtitle = "Unadjusted: logistic regression | Gene-matched: paired Wilcoxon\nPDB = 实验结构 | HM = 同源建模 | PIONEER = 计算预测",
      x        = "OR / Mean proportion ratio (log scale)",
      y        = "PPI source", color = "p-value", shape = "Gene set"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  #add motif flags
  variants_all4_PDB <- variants_motif(
    variants_all2 = variants_all3,
    mart          = ensembl,
    touni_path    = data_file("NIHMS1818854-supplement-2(A).csv"),
    motif_path    = data_file("NIHMS1818854-supplement-2(B).csv"),
    lcs_path      = data_file("Copy of NIHMS1818854-supplement-2.xls")
  )
  variants_all4_PIONEER <- variants_motif(
    variants_all2 = variants_all3_PIONEER,
    mart          = ensembl,
    touni_path    = data_file("NIHMS1818854-supplement-2(A).csv"),
    motif_path    = data_file("NIHMS1818854-supplement-2(B).csv"),
    lcs_path      = data_file("Copy of NIHMS1818854-supplement-2.xls")
  )
  variants_all4_HM <- variants_motif(
    variants_all2 = variants_all3_HM,
    mart          = ensembl,
    touni_path    = data_file("NIHMS1818854-supplement-2(A).csv"),
    motif_path    = data_file("NIHMS1818854-supplement-2(B).csv"),
    lcs_path      = data_file("Copy of NIHMS1818854-supplement-2.xls")
  )
  
  variants_all5 <- variants_all4[!duplicated(variants_all4$Variant_Key), ]
  
  write.csv(variants_all5,'variants_all5.csv',row.names = F)
 variants_all5 = read.csv('variants_all5.csv')
    
    #control for covariates
  
  variants_all5$is_disease          <- as.integer(grepl("disease", variants_all5$group))
  variants_all5$dist_to_cds_end_log <- log(abs(variants_all5$dist_to_cds_end) + 1)
  
  flag_cols   <- c("variant_ppi_overlap", "ptc_after_max_pfam_end", "ptc_before_max_pfam_end",
                   "variant_protein_flag", "variant_domain_flag", "variant_slim_flag",
                   "variant_morf_flag", "variant_ptm_flag", "variant_nls_flag", "variant_LCS_flag")
  confounders <- "cds_length"
  cont_col    <- "dist_to_cds_end_log"
  
  extract_or <- function(fit, coef_name) {
    ci <- confint(fit)[coef_name, ]
    data.frame(
      flag    = coef_name,
      OR      = exp(coef(fit)[coef_name]),
      OR_low  = exp(ci[1]),
      OR_high = exp(ci[2]),
      p_value = summary(fit)$coefficients[coef_name, "Pr(>|z|)"],
      row.names = NULL
    )
  }
  
  # ── unadjusted flag models ─────────────────────────────────────────────────────
  unadj_stats <- lapply(flag_cols, function(col) {
    df  <- variants_all5 %>%
      select(is_disease, all_of(col)) %>%
      filter(!is.na(.data[[col]])) %>%
      mutate(across(all_of(col), as.integer))
    fit <- glm(as.formula(paste("is_disease ~", col)), data = df, family = binomial)
    extract_or(fit, col)
  }) %>% bind_rows() %>% mutate(model = "Unadjusted")
  
  # ── adjusted flag models ───────────────────────────────────────────────────────
  adj_stats <- lapply(flag_cols, function(col) {
    df  <- variants_all5 %>%
      select(is_disease, all_of(col), cds_length, NMDesc_region_length, GC_Content) %>%
      filter(!is.na(.data[[col]])) %>%
      mutate(across(all_of(col), as.integer))
    fit <- glm(as.formula(paste("is_disease ~", col, "+", confounders)), data = df, family = binomial)
    extract_or(fit, col)
  }) %>% bind_rows() %>% mutate(model = "Adjusted")
  
  # ── cont: filter out NA/Inf before fitting ─────────────────────────────────────
  df_cont <- variants_all5 %>%
    filter(!is.na(dist_to_cds_end_log), is.finite(dist_to_cds_end_log),
           !is.na(is_disease), !is.na(cds_length),
           !is.na(NMDesc_region_length), !is.na(GC_Content))
  
  # ── unadjusted dist_to_cds_end_log ────────────────────────────────────────────
  fit_u      <- glm(is_disease ~ dist_to_cds_end_log, data = df_cont, family = binomial)
  cont_unadj <- extract_or(fit_u, cont_col) %>% mutate(model = "Unadjusted")
  
  # ── adjusted dist_to_cds_end_log ──────────────────────────────────────────────
  fit_a      <- glm(as.formula(paste("is_disease ~", cont_col, "+", confounders)),
                    data = df_cont, family = binomial)
  cont_adj   <- extract_or(fit_a, cont_col) %>% mutate(model = "Adjusted")
  
  # ── combine, adjust p-values ───────────────────────────────────────────────────
  combined <- bind_rows(unadj_stats, adj_stats, cont_unadj, cont_adj) %>%
    group_by(model) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"),
           sig   = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                             p_adj < 0.05  ~ "*",   TRUE          ~ "ns")) %>%
    ungroup()
  
   ggplot(combined, aes(x = OR, y = reorder(flag, OR), color = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.2) +
    geom_point(size = 3) +
    facet_wrap(~ model) +
    scale_color_manual(values = c("***" = "red", "**" = "orange",
                                  "*"   = "gold", "ns" = "grey60")) +
    scale_x_log10() +
    labs(title    = "Unadjusted vs Adjusted OR for disease association",
         subtitle  = "Adjusted for cds_length",
         x = "Odds Ratio (log scale)", y = NULL, color = "Adj. p") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

   
   
   ----------------------
#match by gene
   
   # ── Encode ────────────────────────────────────────────────────────────────────
   variants_all5 <- variants_all5 %>%
     mutate(
       is_disease          = as.integer(group %in% c("snv_disease", "fs_disease")),
       gene_set            = case_when(
         group %in% c("snv_disease", "snv_control") ~ "SNV",
         group %in% c("fs_disease",  "fs_control")  ~ "FS"
       ),
       dist_to_cds_end_log = log(abs(dist_to_cds_end) + 1)
     )
   
   flag_cols <- c("variant_ppi_overlap", "ptc_after_max_pfam_end", "ptc_before_max_pfam_end",
                  "variant_protein_flag", "variant_domain_flag", "variant_slim_flag",
                  "variant_morf_flag", "variant_ptm_flag", "variant_nls_flag", "variant_LCS_flag")
   cont_col  <- "dist_to_cds_end_log"
   
   # ── Check which flags are usable per gene set ─────────────────────────────────
   check_flags <- function(data, gene_set_label) {
     data %>%
       filter(gene_set == gene_set_label) %>%
       select(all_of(flag_cols)) %>%
       summarise(across(everything(), ~ n_distinct(na.omit(.)))) %>%
       pivot_longer(everything(), names_to = "flag", values_to = "n_unique") %>%
       filter(n_unique < 2) %>%
       pull(flag)
   } #if one flag has less than 2 unique values (after removing NAs) in a gene set, it can't be used in the model for that set
   
   skip_snv <- check_flags(variants_all5, "SNV")
   skip_fs  <- check_flags(variants_all5, "FS")
   
   message("Skipping in SNV: ", if (length(skip_snv)) paste(skip_snv, collapse = ", ") else "none")
   message("Skipping in FS: ",  if (length(skip_fs))  paste(skip_fs,  collapse = ", ") else "none")
   
   # ── OR extractor ──────────────────────────────────────────────────────────────
   extract_or <- function(fit, coef_name) {
     ci <- confint(fit)[coef_name, ]
     data.frame(
       flag      = coef_name,
       OR        = exp(coef(fit)[coef_name]),
       OR_low    = exp(ci[1]),
       OR_high   = exp(ci[2]),
       p_value   = summary(fit)$coefficients[coef_name, "Pr(>|z|)"],
       row.names = NULL
     )
   }
   
   # ══════════════════════════════════════════════════════════════════════════════
   # PART 1 — UNADJUSTED (logistic regression)
   # ══════════════════════════════════════════════════════════════════════════════
   run_unadjusted <- function(data, gene_set_label, skip_flags = character(0)) {
     
     active_flags <- setdiff(flag_cols, skip_flags)
     
     # binary flags
     flag_res <- lapply(active_flags, function(col) {
       df <- data %>%
         select(is_disease, all_of(col)) %>%
         filter(!is.na(.data[[col]])) %>%
         mutate(across(all_of(col), as.integer))
       
       if (nrow(df) == 0 ||
           length(unique(df[[col]])) < 2 ||
           length(unique(df$is_disease)) < 2) {
         message(sprintf("Skipping %s in %s: no variation", col, gene_set_label))
         return(NULL)
       }
       
       fit <- tryCatch(
         glm(as.formula(paste("is_disease ~", col)), data = df, family = binomial),
         error = function(e) {
           message(sprintf("Model failed for %s in %s: %s", col, gene_set_label, e$message))
           NULL
         }
       )
       if (is.null(fit)) return(NULL)
       extract_or(fit, col)
     }) %>% bind_rows()
     
     # continuous
     df_cont <- data %>%
       filter(!is.na(dist_to_cds_end_log), is.finite(dist_to_cds_end_log))
     
     cont_res <- tryCatch({
       extract_or(
         glm(is_disease ~ dist_to_cds_end_log, data = df_cont, family = binomial),
         cont_col
       )
     }, error = function(e) {
       message(sprintf("Continuous model failed in %s: %s", gene_set_label, e$message))
       NULL
     })
     
     bind_rows(flag_res, cont_res) %>%
       mutate(gene_set = gene_set_label, model = "Unadjusted")
   }
   
   # ══════════════════════════════════════════════════════════════════════════════
   # RUN MODELS
   # ══════════════════════════════════════════════════════════════════════════════
   unadj_snv   <- run_unadjusted(filter(variants_all5, gene_set == "SNV"), "SNV", skip_snv)
   unadj_fs    <- run_unadjusted(filter(variants_all5, gene_set == "FS"),  "FS",  skip_fs)
   
   # ══════════════════════════════════════════════════════════════════════════════
   # COMBINE & BH-ADJUST
   # ══════════════════════════════════════════════════════════════════════════════
   combined <- bind_rows(unadj_snv, unadj_fs) %>%
     group_by(gene_set, model) %>%
     mutate(
       p_adj = p.adjust(p_value, method = "BH"),
       sig   = case_when(
         p_adj < 0.001 ~ "***",
         p_adj < 0.01  ~ "**",
         p_adj < 0.05  ~ "*",
         TRUE          ~ "ns"
       )
     ) %>%
     ungroup() %>%
   mutate(model = factor(model, levels = c("Unadjusted")))
   
   # ══════════════════════════════════════════════════════════════════════════════
   # PLOT
   # ══════════════════════════════════════════════════════════════════════════════
   ggplot(combined, aes(x = OR, y = reorder(flag, OR), color = sig, shape = gene_set)) +
     geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
     geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.2,
                    position = position_dodge(width = 0.5)) +
     geom_point(size = 3, position = position_dodge(width = 0.5)) +
     facet_wrap(~ model) +
     scale_color_manual(values = c("***" = "red", "**" = "orange",
                                   "*"   = "gold", "ns" = "grey60")) +
     scale_shape_manual(values = c("SNV" = 16, "FS" = 17)) +
     scale_x_log10() +
     labs(
       title    = "Unadjusted vs Gene-matched analysis",
       subtitle = "Binary flags: paired Wilcoxon on gene-level proportions | Continuous: paired Wilcoxon on gene-level means",
       x        = "OR (unadjusted) / Mean proportion ratio (gene-matched, log scale)",
       y        = NULL,
       color    = "Adj. p",
       shape    = "Gene set"
     ) +
     theme_minimal(base_size = 12) +
     theme(plot.title = element_text(face = "bold"))
  
