# ============================================================================
# match_and_plot_fromlists.R
# Same matching + plotting as match_and_plot.R, but the CASE genes are taken
# from the two txt lists (so the new gene lists drive the analysis directly),
# rather than from the `group` column of gene_all_AD.csv.
#
# gene_all_AD.csv still supplies the COVARIATES for every gene (cases + the
# candidate control pool) and must contain:
#   hgnc_symbol, ensembl_transcript_id, cds_length,
#   NMD_region_start, NMD_region_end, NMDesc_region_length
#   and a `group` column identifying the control pools
#   (fs_control / snv_control).  If those pool labels are absent, any AD gene
#   not in either case list is treated as a candidate control.
#
# Genes present in BOTH lists (n = 50 here) are matched in BOTH classes.
# Case genes are always removed from the control pools.
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(scales); library(ggpubr); library(clue)
})

input_csv <- "gene_all_AD.csv"
fs_list   <- "fs_can_AD_acat_FDR0_20_all.txt"
snv_list  <- "snv_can_ADrestricted_bh_FDR0_20_all.txt"
out_dir   <- "."
set.seed(42)

read_syms <- function(p){ x <- trimws(readLines(p)); x <- x[nzchar(x)]
  if (tolower(x[1]) %in% c("hgnc_symbol","protein","gene","symbol")) x <- x[-1]; unique(x) }
fs_syms  <- read_syms(fs_list)
snv_syms <- read_syms(snv_list)

df <- read.csv(input_csv, stringsAsFactors = FALSE) %>%
  mutate(nmd_cds_ratio = NMDesc_region_length / cds_length)
if (!"group" %in% names(df)) df$group <- NA_character_
# keep one row per gene (longest CDS) so a symbol maps to one covariate set
df <- df %>% arrange(hgnc_symbol, desc(cds_length)) %>% distinct(hgnc_symbol, .keep_all = TRUE)

case_rows <- function(syms){
  m <- df[df$hgnc_symbol %in% syms, , drop = FALSE]
  missing <- setdiff(syms, df$hgnc_symbol)
  if (length(missing)) message(sprintf("  %d case genes not in gene_all_AD.csv: %s",
                                        length(missing), paste(missing, collapse=", ")))
  m
}
pool_rows <- function(ctrl_label, case_syms){
  p <- if (any(df$group == ctrl_label, na.rm = TRUE))
         df[df$group %in% ctrl_label, , drop = FALSE]
       else df[!(df$hgnc_symbol %in% union(fs_syms, snv_syms)), , drop = FALSE]
  p[!(p$hgnc_symbol %in% case_syms), , drop = FALSE]         # never match a case to a case
}

match_optimal <- function(cases, pool, label){
  stopifnot(nrow(cases) > 0, nrow(pool) >= nrow(cases))
  fc <- data.frame(cds = log10(cases$cds_length), ratio = cases$nmd_cds_ratio)
  fp <- data.frame(cds = log10(pool$cds_length),  ratio = pool$nmd_cds_ratio)
  combined <- rbind(fc, fp); mu <- colMeans(combined)
  sdv <- apply(combined, 2, sd); sdv[sdv == 0] <- 1
  fcz <- scale(fc, center = mu, scale = sdv); fpz <- scale(fp, center = mu, scale = sdv)
  cost <- as.matrix(dist(rbind(fcz, fpz)))[seq_len(nrow(fcz)),
                                           nrow(fcz) + seq_len(nrow(fpz))]^2
  a <- as.integer(clue::solve_LSAP(cost))
  data.frame(case_group = label,
    case_hgnc = cases$hgnc_symbol, case_transcript = cases$ensembl_transcript_id,
    case_cds_length = cases$cds_length,
    case_NMD_region_start = cases$NMD_region_start, case_NMD_region_end = cases$NMD_region_end,
    case_NMDesc_region_length = cases$NMDesc_region_length, case_nmd_cds_ratio = cases$nmd_cds_ratio,
    control_group = pool$group[a], control_hgnc = pool$hgnc_symbol[a],
    control_transcript = pool$ensembl_transcript_id[a], control_cds_length = pool$cds_length[a],
    control_NMD_region_start = pool$NMD_region_start[a], control_NMD_region_end = pool$NMD_region_end[a],
    control_NMDesc_region_length = pool$NMDesc_region_length[a], control_nmd_cds_ratio = pool$nmd_cds_ratio[a],
    cds_log10_diff = abs(log10(cases$cds_length) - log10(pool$cds_length[a])),
    ratio_abs_diff = abs(cases$nmd_cds_ratio - pool$nmd_cds_ratio[a]),
    stringsAsFactors = FALSE)
}

message("Matching SNV ..."); snv_pairs <- match_optimal(case_rows(snv_syms), pool_rows("snv_control", snv_syms), "snv")
message("Matching FS  ..."); fs_pairs  <- match_optimal(case_rows(fs_syms),  pool_rows("fs_control",  fs_syms),  "fs")
all_pairs <- bind_rows(snv_pairs, fs_pairs)

write.csv(all_pairs, file.path(out_dir,"matched_pairs_full.csv"), row.names = FALSE)
write.csv(all_pairs %>% transmute(matched_to_case_group = case_group,
    hgnc_symbol = control_hgnc, ensembl_transcript_id = control_transcript,
    cds_length = control_cds_length, NMDesc_region_length = control_NMDesc_region_length,
    nmd_cds_ratio = control_nmd_cds_ratio),
  file.path(out_dir,"control_gene_list.csv"), row.names = FALSE)

cat("\n--- Match quality ---\n")
for (lab in c("snv","fs")){ m <- filter(all_pairs, case_group == lab)
  cat(sprintf("  %s (n=%d): MWU p_CDS=%.3f p_ratio=%.3f median|log10 CDS diff|=%.3f median|ratio diff|=%.3f\n",
      lab, nrow(m), wilcox.test(m$case_cds_length,m$control_cds_length)$p.value,
      wilcox.test(m$case_nmd_cds_ratio,m$control_nmd_cds_ratio)$p.value,
      median(m$cds_log10_diff), median(m$ratio_abs_diff))) }

# ---- plotting is identical to match_and_plot.R (sourced from section 3 on) ----
# The four plot groups are built from the txt-defined cases + matched controls:
fs_case  <- df %>% filter(hgnc_symbol %in% fs_syms)  %>% transmute(hgnc_symbol, ensembl_transcript_id, cds_length, NMDesc_region_length, nmd_cds_ratio, plot_group="Frameshift")
snv_case <- df %>% filter(hgnc_symbol %in% snv_syms) %>% transmute(hgnc_symbol, ensembl_transcript_id, cds_length, NMDesc_region_length, nmd_cds_ratio, plot_group="Stopgain")
fs_ctrl  <- fs_pairs  %>% transmute(hgnc_symbol=control_hgnc, ensembl_transcript_id=control_transcript, cds_length=control_cds_length, NMDesc_region_length=control_NMDesc_region_length, nmd_cds_ratio=control_nmd_cds_ratio, plot_group="FS Control")
snv_ctrl <- snv_pairs %>% transmute(hgnc_symbol=control_hgnc, ensembl_transcript_id=control_transcript, cds_length=control_cds_length, NMDesc_region_length=control_NMDesc_region_length, nmd_cds_ratio=control_nmd_cds_ratio, plot_group="Stopgain Control")
plot_df <- bind_rows(fs_case, fs_ctrl, snv_case, snv_ctrl) %>%
  mutate(plot_group = factor(plot_group, levels=c("Frameshift","FS Control","Stopgain","Stopgain Control")))
n_by <- plot_df %>% count(plot_group)
subtitle_str <- paste(paste(sprintf("%s (n=%d)", n_by$plot_group, n_by$n), collapse=" \u00b7 "), "\u00b7 Wilcoxon p shown")

group_cols <- c("Frameshift"="#3A6FB0","FS Control"="#A9C8E4","Stopgain"="#C03B3B","Stopgain Control"="#F2B89A")
theme_nmd <- theme_minimal(base_size=16) + theme(plot.title=element_text(face="bold",hjust=0.5,size=18),
  plot.subtitle=element_text(hjust=0.5,color="grey40",size=13), legend.position="bottom",
  legend.text=element_text(size=14), panel.grid.minor=element_blank(),
  panel.grid.major=element_line(linewidth=0.4), axis.title=element_text(size=15),
  axis.text=element_text(size=13), axis.line=element_line(linewidth=0.7), axis.ticks=element_line(linewidth=0.7))
comparisons <- list(c("Frameshift","FS Control"), c("Stopgain","Stopgain Control"))
violin_layer <- function(yvar, ylab, title, sub, log_y=TRUE){
  p <- ggplot(plot_df, aes(x=plot_group, y=.data[[yvar]], fill=plot_group, colour=plot_group)) +
    geom_violin(alpha=0.55, trim=FALSE, scale="width", linewidth=0.9) +
    geom_jitter(width=0.13, size=2.0, alpha=0.55, shape=16) +
    geom_boxplot(width=0.16, fill="white", colour="grey20", outlier.shape=NA, linewidth=0.7, alpha=0.9) +
    scale_fill_manual(values=group_cols, guide="none") + scale_colour_manual(values=group_cols, guide="none") +
    stat_compare_means(comparisons=comparisons, method="wilcox.test", label="p.format", size=4.8, tip.length=0.01) +
    labs(title=title, subtitle=sub, x=NULL, y=ylab) + theme_nmd
  if (log_y) p <- p + scale_y_log10(labels=comma); p }
p_cds   <- violin_layer("cds_length","CDS length (bp)","CDS Length","log10 scale",TRUE)
p_nmd   <- violin_layer("NMDesc_region_length","NMD escape region length (bp)","NMD Escape Region Length","log10 scale",TRUE)
p_ratio <- violin_layer("nmd_cds_ratio","NMD escape / CDS length","NMD / CDS Length Ratio","proportion of CDS in NMD escape region",FALSE)
plot_df <- plot_df %>% mutate(is_case = plot_group %in% c("Frameshift","Stopgain"))
p_scatter <- ggplot(plot_df, aes(x=cds_length, y=NMDesc_region_length, colour=plot_group, fill=plot_group, shape=plot_group)) +
  geom_abline(data=data.frame(frac=c(.10,.25,.50)), aes(slope=frac, intercept=0), linetype="dashed", colour="grey55", linewidth=0.7) +
  geom_point(size=3.6, stroke=1.0, alpha=0.85) + scale_x_log10(labels=comma) + scale_y_log10(labels=comma) +
  scale_colour_manual(values=group_cols) +
  scale_fill_manual(values=c("Frameshift"="#3A6FB0","FS Control"=NA,"Stopgain"="#C03B3B","Stopgain Control"=NA)) +
  scale_shape_manual(values=c("Frameshift"=21,"FS Control"=21,"Stopgain"=24,"Stopgain Control"=24)) +
  labs(title="CDS Length vs NMD Escape Region Length",
       subtitle="Dashed lines = 10%, 25%, 50% of CDS \u00b7 log10 axes \u00b7 filled = study, open = control",
       x="CDS length (bp, log10)", y="NMD escape region length (bp, log10)", colour=NULL, fill=NULL, shape=NULL) +
  theme_nmd + guides(fill="none", colour=guide_legend(override.aes=list(shape=c(21,21,24,24), size=4.5, fill=c("#3A6FB0",NA,"#C03B3B",NA))))
ecdf_theme <- theme_nmd + theme(legend.title=element_blank())
p_ecdf_len <- ggplot(plot_df, aes(x=NMDesc_region_length, colour=plot_group, linetype=plot_group)) +
  stat_ecdf(geom="step", linewidth=1.4) + scale_x_log10(labels=comma) + scale_y_continuous(labels=percent) +
  scale_colour_manual(values=group_cols) +
  scale_linetype_manual(values=c("Frameshift"="solid","FS Control"="dashed","Stopgain"="solid","Stopgain Control"="dashed")) +
  labs(title="ECDF \u2014 NMD Escape Length", x="NMD escape length (bp, log10)", y="Cumulative proportion") + ecdf_theme
p_ecdf_ratio <- ggplot(plot_df, aes(x=nmd_cds_ratio, colour=plot_group, linetype=plot_group)) +
  stat_ecdf(geom="step", linewidth=1.4) + scale_y_continuous(labels=percent) + scale_colour_manual(values=group_cols) +
  scale_linetype_manual(values=c("Frameshift"="solid","FS Control"="dashed","Stopgain"="solid","Stopgain Control"="dashed")) +
  labs(title="ECDF \u2014 NMD / CDS Fraction", x="NMD escape / CDS ratio", y="Cumulative proportion") + ecdf_theme
top_row <- (p_cds | p_nmd | p_ratio)
bottom_row <- (p_ecdf_len | p_ecdf_ratio) + plot_layout(guides="collect") & theme(legend.position="bottom")
final <- (top_row / p_scatter / bottom_row) + plot_layout(heights=c(1,1.2,0.9)) +
  plot_annotation(title="NMD Escape Region vs CDS Length \u2014 4 Gene Groups", subtitle=subtitle_str,
    theme=theme(plot.title=element_text(face="bold",hjust=0.5,size=22), plot.subtitle=element_text(hjust=0.5,color="grey40",size=14)))
ggsave(file.path(out_dir,"nmd_4group_plot.pdf"), final, width=18, height=20, device=cairo_pdf)
ggsave(file.path(out_dir,"nmd_4group_plot.png"), final, width=18, height=20, dpi=200)
cat("\nWrote: matched_pairs_full.csv, control_gene_list.csv, nmd_4group_plot.pdf/png\n")
