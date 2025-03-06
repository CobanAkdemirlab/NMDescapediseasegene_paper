library(ggplot2)
library(ggpubr)


plus1_css = read.table("~/Desktop/frameshift/plus1_css_gene0115.txt", quote="\"", comment.char="",col.names = 'gene')
plus1_can = read.table("~/Desktop/frameshift/plus1_can_gene0115.txt", quote="\"", comment.char="",col.names = 'gene')
plus1_long = read.table("~/Desktop/frameshift/plus1_long_gene0115.txt", quote="\"", comment.char="",col.names = 'gene')
plus1_trigger = read.table("~/Desktop/frameshift/plus1_trigger_gene0115.txt", quote="\"", comment.char="",col.names = 'gene')
plus2_css = read.table("~/Desktop/frameshift/plus2_css_gene0115.txt", quote="\"", comment.char="",col.names = 'gene')
plus2_can = read.table("~/Desktop/frameshift/plus2_can_gene0115.txt", quote="\"", comment.char="",col.names = 'gene')
plus2_long = read.table("~/Desktop/frameshift/plus2_long_gene0115.txt", quote="\"", comment.char="",col.names = 'gene')
plus2_trigger = read.table("~/Desktop/frameshift/plus2_trigger_gene0115.txt", quote="\"", comment.char="",col.names = 'gene')
control = read.table("~/Desktop/new_clinvar/snv_list/list4/snv_plp_ptc_p1120_NMDenriched2_control.txt", quote="\"", comment.char="",col.names = 'gene')


gnomad.v2.1.1.lof_metrics.by_gene <- read.delim("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/桌面 - q的MacBook Pro/autism/data/gnomad.v2.1.1.lof_metrics.by_gene.txt")
plus1_css_combined = merge(plus1_css,gnomad.v2.1.1.lof_metrics.by_gene[,c('gene','cds_length')], by="gene")
plus1_can_combined = merge(plus1_can,gnomad.v2.1.1.lof_metrics.by_gene[,c('gene','cds_length')], by="gene")
plus1_long_combined = merge(plus1_long,gnomad.v2.1.1.lof_metrics.by_gene[,c('gene','cds_length')], by="gene")
plus1_trigger_combined = merge(plus1_trigger,gnomad.v2.1.1.lof_metrics.by_gene[,c('gene','cds_length')], by="gene")
plus2_css_combined = merge(plus2_css,gnomad.v2.1.1.lof_metrics.by_gene[,c('gene','cds_length')], by="gene")
plus2_can_combined = merge(plus2_can,gnomad.v2.1.1.lof_metrics.by_gene[,c('gene','cds_length')], by="gene")
plus2_long_combined = merge(plus2_long,gnomad.v2.1.1.lof_metrics.by_gene[,c('gene','cds_length')], by="gene")
plus2_trigger_combined = merge(plus2_trigger,gnomad.v2.1.1.lof_metrics.by_gene[,c('gene','cds_length')], by="gene")
control_combined = merge(control,gnomad.v2.1.1.lof_metrics.by_gene[,c('gene','cds_length')], by="gene")

plus1_css_combined$obj = 'plus1_css'
plus1_can_combined$obj = 'plus1_can'
plus1_long_combined$obj = 'plus1_long'
plus1_trigger_combined$obj = 'plus1_trigger'
plus2_css_combined$obj = 'plus2_css'
plus2_can_combined$obj = 'plus2_can'
plus2_long_combined$obj = 'plus2_long'
plus2_trigger_combined$obj = 'plus2_trigger'
control_combined$obj = 'control'

data_combined = rbind(plus1_css_combined,plus1_can_combined,plus1_long_combined,plus1_trigger_combined,plus2_css_combined,plus2_can_combined,plus2_long_combined,plus2_trigger_combined,control_combined)
  
g = ggplot(data_combined, aes(x = obj, y = cds_length, color = obj)) +
  geom_boxplot() +
  #set y lim 0-100000
  #coord_cartesian(ylim = c(200, 20000)) +
  #scale_color_manual(values = c('cornflowerblue', 'brown1', 'chartreuse', 'darkorchid1', 'goldenrod1','grey')) +
  #add p value 
  scale_y_log10() +
  labs(title = '', y = "CDS Length", x = "") +
  theme_minimal() +
  stat_compare_means(aes(label = ..p.format..), 
                     method = "t.test", 
                     ref.group = "control", 
                     label = "p.signif")

#why pli low in trig is high
#vcf-avision-somatic information key-2
#intre disorder region
#VEP NMD plugin


ggplot(data_combined, aes(x = obj, y = cds_length, color = obj)) +
  geom_boxplot() +
  scale_color_manual(values = c('cornflowerblue', 'brown1', 'chartreuse', 'darkorchid1', 'goldenrod1')) +
  scale_y_log10() +
  labs(title = '', y = "CDS Length", x = "") +
  theme_minimal() +
  stat_compare_means(method = "t.test", label = "p.signif") 


pairwise_comparisons <- compare_means(
  cds_length ~ obj, 
  data = data_combined, 
  method = "t.test"
)

# Determine y positions for each comparison
pairwise_comparisons$y.position <- max(data_combined$cds_length) * 1.1 # adjust multiplier if needed

# Plot with ggboxplot and pairwise comparison annotations
g <- ggboxplot(
  data_combined, 
  x = "obj", 
  y = "cds_length", 
  color = "obj",
  add = "jitter"
) +
  scale_y_log10() +
  labs(title = '', y = "CDS Length", x = "") +
  scale_color_manual(values = c('cornflowerblue', 'brown1', 'chartreuse', 'darkorchid1', 'goldenrod1')) +
  stat_compare_means(method = "t.test", label = "p.signif", paired = FALSE) +
  theme_minimal()


comparisons <- list(
  c("all", "control"),
  c("css", "control"),
  c("long", "control"),
  c("can", "control")
)

data_noout = data_combined[data_combined$cds_length < 100000,]
# Plot with pairwise comparisons and hiding non-significant labels
g <- ggboxplot(
  data_noout, 
  x = "obj", 
  y = "cds_length", 
  color = "obj",
  add = "jitter"
) +
  scale_y_log10() +
  #scale_y_log10(limits = c(1, max(data_combined$cds_length) * 1.5)) +
  labs(title = '', y = "CDS Length", x = "") +
  scale_color_manual(values = c('cornflowerblue', 'brown1', 'chartreuse', 'darkorchid1', 'goldenrod1')) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif", hide.ns = TRUE) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15, face = "bold"), # make x-axis label bigger and bold
    axis.text.y = element_text(size = 15, face = "bold")  # make y-axis label bigger and bold
  )

t.test(data_noout[data_noout$obj == 'all',]$cds_length, data_noout[data_noout$obj == 'css',]$cds_length)
