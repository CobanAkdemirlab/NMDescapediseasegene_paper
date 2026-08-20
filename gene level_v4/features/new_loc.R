library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(AnnotationDbi)
# --- 路径解析层（自动插入）------------------------------------------------
# 原来这里是旧机器 /Users/jxu14/ 的绝对路径，换机器必断。改走 paths.R：
#   data_file("x.csv")  按文件名定位，找不到会明确报错
#   out_file("y.csv")   输出到 NMDESC_OUT（默认 ~/Desktop/NMDesc_out）
#   data_root("clinvar") 需要目录而非文件时用
.p <- c("gene level_v3/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../lib/paths.R", "../../../../lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p)) stop("找不到 paths.R —— 请从仓库根目录运行 R")
source(.p[1]); rm(.p)
# --------------------------------------------------------------------------

ensgene <- txdb
chr.list <- paste('chr',1:22,sep='')
ensgene.sub <- keepSeqlevels(ensgene,chr.list)
ensgene <- ensgene.sub
cds_seqs <- extractTranscriptSeqs(Hsapiens, cdsBy(ensgene, by="tx",use.names=TRUE))
utr.grange <-threeUTRsByTranscript(ensgene, use.names=TRUE)
fiveutr.grange<-fiveUTRsByTranscript(ensgene, use.names=TRUE)
introns.grange<- intronsByTranscript(ensgene, use.names=TRUE)
exons.grange<-exonsBy(ensgene, use.names=TRUE)
cds.grange<- cdsBy(ensgene, by="tx",use.names=TRUE)

cds.resized=resize(cds.grange,width=200,fix='start')
nucloc <- read.csv('~/Desktop/unipLocSignal.bed.csv',header=FALSE,sep='\t',stringsAsFactors=FALSE)
#nucloc <- nucloc[which(nucloc$V22=='short sequence motif'),]
gr2 <- GRanges(seqnames = nucloc$V1, ranges = IRanges(start = nucloc$V2, end=nucloc$V3))
hits <- unique(queryHits(findOverlaps(cds.resized,gr2)))
cds.nucloc <- names(cds.grange)[hits]
merged.2$cds.nucloc <- rep('NA',nrow(merged.2))
merged.2[which(merged.2$txnames%in%cds.nucloc),'cds.nucloc'] <- 'There is a nuclear localization signal in the first 200 bp'
save(merged.2, file='~/Desktop/NMDrevisioncode/hg38_seqfeatures.RData')
variants.features.fr$cds.nucloc <- rep('NA',nrow(variants.features.fr))
variants.features.fr[which(variants.features.fr$txnames%in%cds.nucloc),'cds.nucloc'] <- 'There is a nuclear localization signal in the first 200 bp'

load(data_file("hg38_seqfeatures.RData", must = FALSE))
tx2gene <- getBM(
  attributes = c("ensembl_transcript_id", "hgnc_symbol"),
  filters = "ensembl_transcript_id",
  values = unique(substr(merged.2$txnames,1,15)),
  mart = ensembl
)
merged.2$txnames <- substr(merged.2$txnames,1,15)

merged.2 <- merge(merged.2, tx2gene, by.x = "txnames", by.y = "ensembl_transcript_id", all.x = TRUE)
merged_nls <- merged.2 %>%
  filter(cds.nucloc == "There is a nuclear localization signal in the first 200 bp") %>%
  dplyr::select(gene = hgnc_symbol) %>%
  distinct()
pli_all$cds_nls <- ifelse(pli_all$gene %in% merged_nls$gene, "Yes", "No")
library(ggplot2)
library(dplyr)

summary_tbl <- pli_all %>%
  group_by(category, cds_nls) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(category) %>%
  mutate(percentage = n / sum(n))

ggplot(summary_tbl, aes(x = category, y = percentage, fill = cds_nls)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("Yes" = "#377EB8", "No" = "#E41A1C")) +
  labs(
    title = "Nuclear Localization Signal in First 200 bp of CDS",
    x = "Gene Category",
    y = "Percentage of Genes",
    fill = "NLS Present"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

summary_yes <- summary_tbl %>%
  filter(cds_nls == "Yes")

# Create wide contingency table for Fisher’s test
fisher_matrix <- pli_all %>%
  group_by(category, cds_nls) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = cds_nls, values_from = count, values_fill = 0)

# Create pairwise comparisons
categories <- unique(fisher_matrix$category)
comparisons <- combn(categories, 2, simplify = FALSE)

# Run Fisher’s exact test
fisher_results <- lapply(comparisons, function(pair) {
  group1 <- fisher_matrix %>% filter(category == pair[1]) %>% select(Yes, No)
  group2 <- fisher_matrix %>% filter(category == pair[2]) %>% select(Yes, No)
  
  mat <- rbind(as.numeric(group1), as.numeric(group2))
  pval <- fisher.test(mat)$p.value
  
  data.frame(group1 = pair[1], group2 = pair[2], p.value = pval)
})

fisher_df <- do.call(rbind, fisher_results)

