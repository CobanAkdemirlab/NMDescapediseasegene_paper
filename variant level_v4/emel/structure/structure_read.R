library(tidyverse)
library(ggpubr)
library(patchwork)
library(dplyr)
library(ggplot2)
library(scales)

VAR_WT_structural_results <- read_csv("~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/structure/VAR-WT_structural_results.csv")
VARs_full_NMD_structural <- read_csv("~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/structure/VARs_full_NMD_structural.csv")
WT_full_NMD_structural <- read_csv("~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/structure/WT_full_NMD_structural.csv")

str(VAR_WT_structural_results)
str(VARs_full_NMD_structural)
str(WT_full_NMD_structural)

# ── 1. Define pLDDT features ──────────────────────────────────────────────────

plddt_features <- tribble(
  ~col,                          ~label,
  "VAR_plddt_full_mean",         "VAR pLDDT – Full (mean)",
  "VAR_plddt_NMDPos_mean",       "VAR pLDDT – NMD (mean)",
  "VAR_plddt_DivergentPos_mean", "VAR pLDDT – Divergent (mean)",
  "WT_plddt_full_mean",          "WT pLDDT – Full (mean)",
  "WT_plddt_NMDPos_mean",        "WT pLDDT – NMD (mean)",
  "WT_plddt_DivergentPos_mean",  "WT pLDDT – Divergent (mean)"
)

# ── 2. Prepare data ───────────────────────────────────────────────────────────
## 2.1 rename the folders
plddt_data <- VAR_WT_structural_results %>%
  mutate(
    group = factor(source_folder,
                   levels = c("snv_disease", "snv_control",
                              "fs_disease",  "fs_control"),
                   labels = c("SNV\nDisease", "SNV\nControl",
                              "FS\nDisease",  "FS\nControl"))
  )

## 2.2 select disease variants present in new variant list
snv_genes <- read_csv("snv_can_ADrestricted_bh_FDR0.20_all.txt")
fs_genes <- read_csv("fs_can_AD_acat_FDR0.20_all.txt")
all_genes <- bind_rows(snv_genes, fs_genes)
snv_variants = read_csv("snv_variants20260201_plp_bh0.2_clinvar.csv")
fs_variants = read_csv("fs_variants20260201_plp_clinvar_acat_0.2_filtered.csv")
all_variants <- bind_rows(snv_variants, fs_variants)
plddt_data_disease <- plddt_data %>%
  filter(key %in% all_variants$key)
plddt_data_control <- plddt_data %>%
  filter(source_folder %in% c("snv_control", "fs_control"))
plddt_data <- bind_rows(plddt_data_disease, plddt_data_control)

# ── 3. Color palette ──────────────────────────────────────────────────────────

panel_colors <- c(
  "SNV\nDisease" = "#D6604D",
  "SNV\nControl" = "#4393C3",
  "FS\nDisease"  = "#F4A582",
  "FS\nControl"  = "#92C5DE"
)

# ── 4. Plot function ──────────────────────────────────────────────────────────

plot_one_feature <- function(df, col_name, feat_label) {
  
  plot_df <- df %>%
    select(group, value = all_of(col_name)) %>%
    filter(!is.na(value))
  
  get_pval <- function(data, grp1, grp2) {
    v1 <- data %>% filter(group == grp1) %>% pull(value)
    v2 <- data %>% filter(group == grp2) %>% pull(value)
    if (length(v1) < 3 || length(v2) < 3) return(NA_real_)
    wilcox.test(v1, v2)$p.value
  }
  
  p_snv <- get_pval(plot_df, "SNV\nDisease", "SNV\nControl")
  p_fs  <- get_pval(plot_df, "FS\nDisease",  "FS\nControl")
  
  fmt_p <- function(p) {
    if (is.na(p))  return("n.s.")
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    return("n.s.")
  }
  
  y_sig <- quantile(plot_df$value, 0.98, na.rm = TRUE) * 1.05
  
  sig_df <- tibble(
    x     = c(1.5, 3.5),
    y     = y_sig,
    label = c(fmt_p(p_snv), fmt_p(p_fs))
  )
  
  ggplot(plot_df, aes(x = group, y = value, fill = group)) +
    geom_violin(alpha = 0.55, trim = TRUE, scale = "width", color = NA) +
    geom_boxplot(width = 0.18, outlier.size = 0.4,
                 outlier.alpha = 0.3, alpha = 0.85, color = "grey30") +
    # reference lines for pLDDT confidence thresholds
    geom_hline(yintercept = 90, linetype = "dotted",
               color = "#2166AC", linewidth = 0.4, alpha = 0.7) +
    geom_hline(yintercept = 70, linetype = "dotted",
               color = "#4DAC26", linewidth = 0.4, alpha = 0.7) +
    geom_hline(yintercept = 50, linetype = "dotted",
               color = "#D01C8B", linewidth = 0.4, alpha = 0.7) +
    geom_vline(xintercept = 2.5, linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_text(data = sig_df,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              size = 3.5, fontface = "bold", color = "grey20") +
    scale_fill_manual(values = panel_colors) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(NA, NA),
                       breaks = c(25, 50, 70, 90)) +
    labs(title = feat_label, x = NULL, y = "pLDDT score") +
    theme_bw(base_size = 10) +
    theme(
      legend.position    = "none",
      plot.title         = element_text(size = 9, face = "bold", hjust = 0.5),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x        = element_text(size = 8),
      plot.margin        = margin(4, 6, 4, 6)
    )
}

# ── 5. Build all plots ────────────────────────────────────────────────────────

plddt_plots <- map2(plddt_features$col, plddt_features$label,
                    ~ plot_one_feature(plddt_data, .x, .y))

# ── 6. Arrange: VAR (top row) | WT (bottom row) ───────────────────────────────
plddt_grid <- wrap_plots(plddt_plots, ncol = 3) +
  plot_annotation(
    title    = "pLDDT Comparison: SNV Disease vs Control  |  FS Disease vs Control",
    subtitle = paste0(
      "pLDDT confidence thresholds: ",
      "\u2022 >90 (blue dotted) = very high  ",
      "\u2022 70-90 (green dotted) = confident  ",
      "\u2022 50-70 = low  ",
      "\u2022 <50 (pink dotted) = very low\n",
      "* p<0.05  ** p<0.01  *** p<0.001  (Wilcoxon rank-sum test)"
    ),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 8.5, color = "grey40")
    )
  )

# ── 7. Save ───────────────────────────────────────────────────────────────────

ggsave("pLDDT_all_features_4panels.pdf", plddt_grid,
       width = 15, height = 10)

ggsave("pLDDT_all_features_4panels.png", plddt_grid,
       width = 15, height = 10, dpi = 150)






# =====================================================
# Calculate ΔSASA
# =====================================================

sasa_delta_df <- VAR_WT_structural_results %>%
  mutate(
    
    delta_total_SASA_full =
      VAR_SASA_total_full -
      WT_SASA_total_full,
    
    delta_total_SASA_NMD =
      VAR_SASA_total_NMDPos -
      WT_SASA_total_NMDPos,
    
    delta_total_SASA_Divergent =
      VAR_SASA_total_DivergentPos -
      WT_SASA_total_DivergentPos
    
  ) %>%
  
  mutate(
    
    group = factor(
      source_folder,
      
      levels = c(
        "snv_disease",
        "snv_control",
        "fs_disease",
        "fs_control"
      ),
      
      labels = c(
        "SNV\nDisease",
        "SNV\nControl",
        "FS\nDisease",
        "FS\nControl"
      )
      
    )
    
  )

sasa_data_disease <- sasa_delta_df %>%
  filter(key %in% all_variants$key)
sasa_data_control <- sasa_delta_df %>%
  filter(source_folder %in% c("snv_control", "fs_control"))
sasa_delta_df <- bind_rows(sasa_data_disease, sasa_data_control)

# =====================================================
# Colors
# =====================================================
my_cols <- c(
  
  "SNV\nDisease"="#E87A72",
  "SNV\nControl"="#6BAED6",
  
  "FS\nDisease"="#F4A582",
  "FS\nControl"="#92C5DE"
  
)

# =====================================================
# Plot function
# =====================================================

plot_delta <- function(
    data,
    yvar,
    title
){
  
  ggplot(
    
    data,
    
    aes(
      
      x=group,
      y=.data[[yvar]],
      fill=group
      
    )
    
  )+
    
    
    geom_violin(
      
      trim=FALSE,
      alpha=.75,
      color=NA
      
    )+
    
    
    geom_boxplot(
      
      width=.15,
      alpha=.7,
      outlier.shape=NA,
      color="gray30"
      
    )+
    
    
    geom_hline(
      
      yintercept=0,
      linetype="dashed"
      
    )+
    
    
    stat_compare_means(
      
      comparisons=list(
        
        c(
          "SNV\nDisease",
          "SNV\nControl"
        ),
        
        c(
          "FS\nDisease",
          "FS\nControl"
        )
        
      ),
      
      method="wilcox.test",
      label="p.signif"
      
    )+
    
    
    scale_fill_manual(
      
      values=my_cols
      
    )+
    
    
    # -------- pseudo log scale --------
  
  scale_y_continuous(
    
    trans =
      pseudo_log_trans(
        
        base=10
        
      )
    
  )+
    
    
    theme_bw(
      
      base_size=12
      
    )+
    
    
    theme(
      
      plot.title =
        element_text(
          
          hjust=.5,
          face="bold"
          
        ),
      
      legend.position="none",
      
      panel.grid.minor=
        element_blank(),
      
      panel.grid.major.x=
        element_blank(),
      
      axis.title.x=
        element_blank(),
      
      axis.title.y=
        element_text(
          
          face="bold"
          
        )
      
    )+
    
    
    labs(
      
      title=title,
      
      y=
        expression(
          Delta*" Total SASA"
        )
      
    )
  
}



# =====================================================
# Generate panels
# =====================================================


p1 <-
  
  plot_delta(
    
    sasa_delta_df,
    
    "delta_total_SASA_full",
    
    expression(
      Delta*" Total SASA - Full ("*ring(A)^2*")"
    )
    
  )



p2 <-
  
  plot_delta(
    
    sasa_delta_df,
    
    "delta_total_SASA_NMD",
    
    expression(
      Delta*" Total SASA - NMD ("*ring(A)^2*")"
    )
    
  )



p3 <-
  
  plot_delta(
    
    sasa_delta_df,
    
    "delta_total_SASA_Divergent",
    
    expression(
      Delta*" Total SASA - Divergent ("*ring(A)^2*")"
    )
    
  )



# =====================================================
# Combine
# =====================================================

p_all <-
  
  p1 +
  p2 +
  p3 +
  
  plot_annotation(
    
    title=
      
      "ΔSASA Comparison:
SNV Disease vs Control |
FS Disease vs Control",
    
    theme=
      
      theme(
        
        plot.title=
          element_text(
            
            hjust=.5,
            
            face="bold",
            
            size=16
            
          )
        
      )
    
  )



print(
  
  p_all
  
)



# =====================================================
# Save
# =====================================================

ggsave(
  
  "delta_SASA_logscale.png",
  
  p_all,
  
  width=15,
  
  height=5,
  
  dpi=300
  
)