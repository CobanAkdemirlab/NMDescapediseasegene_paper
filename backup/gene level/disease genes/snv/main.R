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
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
tx_vec <- unique(as.character(res@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]]))
tx_vec <- tx_vec[!is.na(tx_vec)]
BM.info <- getBM(
  attributes = c("ensembl_transcript_id", "cds_start", "cds_end", "chromosome_name", "strand"),
  filters    = "ensembl_transcript_id",
  values     = tx_vec,
  mart       = ensembl
)


typelist = c('all','can','css','long','trig')

vcf_file = "clinvar_20260201.vcf.gz"
vcf = aenmd:::parse_vcf_VariantAnnotation(vcf_file)
vcf_rng = vcf$vcf_rng
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
plp_ind = grep('pathogenic',res@elementMetadata@listData[["CLNSIG"]]@unlistData,ignore.case = T)
ptc_ind = which(res@elementMetadata@listData[["res_aenmd"]]@listData[["is_ptc"]]==T)
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

plp_ptc_ind = intersect(plp_ind,ptc_ind)
snv_plp_ptc_ind = intersect(snv_ind,plp_ptc_ind)
snv_plp_ptc_nmdesc_ind = intersect(snv_plp_ptc_ind,nmdesc_ind)
snv_plp_ptc_nmdesc_can_ind = intersect(snv_plp_ptc_nmdesc_ind,can_ind)

snv_plp_ptc_nmdesc_can = res[snv_plp_ptc_nmdesc_can_ind]
length(unique(snv_plp_ptc_nmdesc_can@ranges@NAMES))
length(unique(snv_plp_ptc_nmdesc_can@elementMetadata@listData[["key"]]))

#Quality Control, should we match by key or gene id?
variant_summary <- read_delim("variant_summary.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
gene_snv_plp_ptc_nmdesc_can$NumberSubmitters = variant_summary$NumberSubmitters[match(snv_plp_ptc_nmdesc_can$hgnc_symbol, variant_summary$GeneSymbol)]
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



gene_can_AD_plp_ptc_nmdesc_uni$NumberSubmitters = variant_summary$NumberSubmitters[match(gene_can_AD_plp_ptc_nmdesc_uni$hgnc_symbol, variant_summary$GeneSymbol)]
summary(gene_can_AD_plp_ptc_nmdesc_uni$NumberSubmitters)

# do other filters in variant_summary

#ReviewStatus does not contain either of “no assertion” or “no interpretation”
variant_summary$ReviewStatus
v_re = variant_summary %>%
  filter(
    !str_detect(tolower(ReviewStatus), "no assertion"),
    !str_detect(tolower(ReviewStatus), "no interpretation")
  )

#ClinicalSignificance does not contain any of “not provided”, “drug
#response”, “other”, “risk”, “low penetrance”, “conflicting”, “affects”, “association”, “protective”, “confers
#sensitivity”
variant_summary$ClinicalSignificance
exclude_terms <- c(
  "not provided", "drug response", "other", "risk", "low penetrance",
  "conflicting", "affects", "association", "protective", "confers sensitivity"
)
pattern_exclude <- paste(exclude_terms, collapse = "|")
v_cs = v_re %>%
  filter(
    !str_detect(tolower(ClinicalSignificance), pattern_exclude)
  )

#Assembly == GRCh38
variant_summary$Assembly
v_gr = v_cs %>%
  filter(str_trim(tolower(Assembly)) == "grch38")

#autosomal contigs (chr1-22) only
variant_summary$Chromosome
v_ch = v_gr %>%
  mutate(
    Chromosome_clean = str_remove(tolower(str_trim(Chromosome)), "^chr")
  ) %>%
  filter(Chromosome_clean %in% as.character(1:22)) %>%
  select(-Chromosome_clean)

#add canonical and non_canonical transcript tag use Name
variant_summary$Name
#substract from Name(eg. NM_017547.4 from NM_017547.4(FOXRED1):c.694C>T (p.Gln232Ter))
v_ch$tx_id1 = sub("\\(.*", "", v_ch$Name)
v_ch = v_ch %>%
  mutate(
    tx_id1_noversion = sub("\\..*$", "", tx_id1)
  )
#presence of a valid RefSeq ID
v_ch = v_ch %>%
  filter(str_detect(tx_id1_noversion, "^NM_"))
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

#return this information to res, build res3
#the problem is, not all tr in res appeared in this variant_summary
#res3 = res[which(res@ranges@NAMES %in% v_ch$Submitted_variation_ID)]

tx_n$canonical_by_vs = v_ch$tx_canonical[match(tx_n$gene, v_ch$tx_id2)]

saveRDS(snv_plp_ptc_nmdesc_can,'snv_plp_ptc_nmdesc_can20260201.rds')

get_pvalue('snv_plp_ptc_nmdesc_can20260201.rds','snv_plp_ptc_nmdesc_can_p_20260201.rds')
res_p = readRDS('snv_plp_ptc_nmdesc_can_p_20260201.rds')

get_NMD_enrichment('snv_plp_ptc_nmdesc_can_p_20260201.rds',p_cutoff = 0.9,ADfilter = 'off',filter_type = typelist[2])

