library(tidyverse)
library(lme4)
library(lmerTest)
library(gtsummary)
library(gt)
library(ggpubr)
library(patchwork)

# ══════════════════════════════════════════════════════════════
# 0. Parameters
# ══════════════════════════════════════════════════════════════

MIN_PEPTIDE_LEN <- 10
FDR_CUTOFF <- 0.05

base_dir <- "~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/sequence_analysis/CIDER"

# This script assumes all_variants$key has already been created/loaded,
# as in the PARSE v2 workflow.
if (!exists("all_variants") || !"key" %in% names(all_variants)) {
  stop("Please load/create all_variants$key before running this script.")
}


# ══════════════════════════════════════════════════════════════
# 1. Read CIDER results
# ══════════════════════════════════════════════════════════════

cider_div <- read_csv(
  file.path(base_dir, "NMD_region_DivergentPos_CIDER.csv"),
  show_col_types = FALSE
)

cider_nmd <- read_csv(
  file.path(base_dir, "NMD_region_NMDPos_CIDER.csv"),
  show_col_types = FALSE
)


# ══════════════════════════════════════════════════════════════
# 2. Variant key, cohort filtering, and region filtering
# ══════════════════════════════════════════════════════════════

make_key <- function(fn) {
  p <- str_split_fixed(fn, "_", 4)
  paste0(p[, 1], ":", as.numeric(p[, 2]), "|", p[, 3], "|", p[, 4])
}

prep_cider <- function(d, region_name,
                       apply_len_filter = TRUE,
                       drop_snv = FALSE) {
  
  d <- d %>%
    mutate(
      key = make_key(file_name),
      group = case_when(
        source_folder == "snv_disease" ~ "SNV Disease",
        source_folder == "snv_control" ~ "SNV Control",
        source_folder == "fs_disease"  ~ "FS Disease",
        source_folder == "fs_control"  ~ "FS Control",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group))
  
  n0 <- nrow(d)
  
  # Disease variants must be in all_variants;
  # all control variants are retained.
  d <- d %>%
    filter(
      source_folder %in% c("snv_control", "fs_control") |
        key %in% all_variants$key
    )
  
  n1 <- nrow(d)
  
  d <- d %>% filter(!is.na(var_length))
  
  if (apply_len_filter) {
    d <- d %>% filter(var_length >= MIN_PEPTIDE_LEN)
  }
  
  n2 <- nrow(d)
  
  # Divergent region: nonsense SNVs do not generate a novel shifted peptide.
  if (drop_snv) {
    d <- d %>% filter(str_starts(source_folder, "fs_"))
  }
  
  n3 <- nrow(d)
  
  cat(sprintf(
    "[%-10s] raw %d -> cohort %d -> length %d -> after SNV exclusion %d\n",
    region_name, n0, n1, n2, n3
  ))
  
  d %>%
    mutate(
      region = region_name,
      group = factor(
        group,
        levels = c("SNV Disease", "SNV Control",
                   "FS Disease", "FS Control")
      ) %>%
        droplevels()
    )
}

cat("=== Cleaning record ===\n")

df_div <- prep_cider(
  cider_div,
  "Divergent",
  apply_len_filter = TRUE,
  drop_snv = TRUE
)

df_nmd <- prep_cider(
  cider_nmd,
  "NMD",
  apply_len_filter = TRUE,
  drop_snv = FALSE
)


# ══════════════════════════════════════════════════════════════
# 3. Diagnostics
# ══════════════════════════════════════════════════════════════

cat("\n=== Disease-variant key matching ===\n")

key_diag <- imap_dfr(
  list(Divergent = cider_div, NMD = cider_nmd),
  function(d, nm) {
    d %>%
      filter(str_detect(source_folder, "disease")) %>%
      mutate(matched = make_key(file_name) %in% all_variants$key) %>%
      count(source_folder, matched) %>%
      pivot_wider(
        names_from = matched,
        values_from = n,
        names_prefix = "matched_",
        values_fill = 0
      ) %>%
      mutate(file = nm, .before = 1)
  }
)

print(as.data.frame(key_diag))

cat("\nUnmatched disease file_name examples:\n")

bind_rows(
  cider_div %>% mutate(region = "Divergent"),
  cider_nmd %>% mutate(region = "NMD")
) %>%
  filter(
    str_detect(source_folder, "disease"),
    !make_key(file_name) %in% all_variants$key
  ) %>%
  select(region, source_folder, file_name) %>%
  slice_head(n = 10) %>%
  print(n = 10)

cat("\n=== Sample size by region/group ===\n")

bind_rows(df_div, df_nmd) %>%
  count(region, group) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
  as.data.frame() %>%
  print()


# ══════════════════════════════════════════════════════════════
# 4. CIDER variable definition
# ══════════════════════════════════════════════════════════════
#
# † = absolute / strongly length-dependent measure
# ‡ = composition, fraction, or sequence-pattern measure
#
# For absolute measures, the main LMM is accompanied by a
# length-adjusted sensitivity model.
#
# CIDER already provides normalized composition measures such as
# FCR, charge fractions, disorder-promoting fraction, hydropathy,
# kappa, Omega, delta, and deltaMax, so no extra "/length" variables
# need to be created.

var_info <- tribble(
  ~var,                           ~label,                                      ~type,
  "length",                       "Sequence Length (aa)",                       "absolute",
  "fcr_7",                        "Fraction Charged Residues (FCR)",            "normalized",
  "isoelectric_point",            "Isoelectric Point",                         "normalized",
  "molecular_weight",             "Molecular Weight",                          "absolute",
  "count_neg",                    "Negative Residue Count",                     "absolute",
  "count_pos",                    "Positive Residue Count",                     "absolute",
  "count_neut",                   "Neutral Residue Count",                      "absolute",
  "fraction_negative",            "Fraction Negative Residues",                 "normalized",
  "fraction_positive",            "Fraction Positive Residues",                 "normalized",
  "fraction_expanding_7",         "Fraction Expanding Residues",                "normalized",
  "fraction_disorder_promoting",  "Fraction Disorder-Promoting Residues",       "normalized",
  "mean_hydropathy",              "Mean Hydropathy",                            "normalized",
  "uversky_hydropathy",           "Uversky Hydropathy",                        "normalized",
  "kappa",                        "Kappa",                                      "normalized",
  "Omega",                        "Omega",                                      "normalized",
  "delta",                        "Delta",                                      "normalized",
  "deltaMax",                     "DeltaMax",                                   "normalized"
) %>%
  mutate(
    label_disp = paste0(
      label,
      if_else(type == "absolute", "", "")
    )
  )

var_col <- function(v) paste0("var_", v)
wt_col  <- function(v) paste0("wt_", v)


# ══════════════════════════════════════════════════════════════
# 5. Descriptive summaries: VAR, WT, and Δ(VAR − WT)
# ══════════════════════════════════════════════════════════════

q1 <- function(x) quantile(x, 0.25, na.rm = TRUE, names = FALSE)
q3 <- function(x) quantile(x, 0.75, na.rm = TRUE, names = FALSE)

make_desc <- function(d, region_name) {
  
  map_dfr(var_info$var, function(v) {
    
    vc <- var_col(v)
    wc <- wt_col(v)
    
    d %>%
      transmute(
        group,
        var_value = .data[[vc]],
        wt_value  = .data[[wc]],
        delta     = .data[[vc]] - .data[[wc]]
      ) %>%
      group_by(group) %>%
      summarise(
        n_var = sum(!is.na(var_value)),
        var_median = median(var_value, na.rm = TRUE),
        var_q1 = q1(var_value),
        var_q3 = q3(var_value),
        
        n_wt = sum(!is.na(wt_value)),
        wt_median = median(wt_value, na.rm = TRUE),
        wt_q1 = q1(wt_value),
        wt_q3 = q3(wt_value),
        
        n_delta = sum(!is.na(delta)),
        delta_median = median(delta, na.rm = TRUE),
        delta_q1 = q1(delta),
        delta_q3 = q3(delta),
        .groups = "drop"
      ) %>%
      mutate(
        region = region_name,
        var = v,
        .before = 1
      )
  }) %>%
    left_join(
      var_info %>% select(var, label, label_disp, type),
      by = "var"
    )
}

desc_all <- bind_rows(
  make_desc(df_div, "Divergent"),
  make_desc(df_nmd, "NMD")
)

write_csv(desc_all, "CIDER_descriptive_summary.csv")

cat("\n=== CIDER descriptive summary ===\n")
desc_all %>%
  select(region, group, label, var_median, var_q1, var_q3,
         wt_median, wt_q1, wt_q3,
         delta_median, delta_q1, delta_q3) %>%
  print(n = Inf)


# ══════════════════════════════════════════════════════════════
# 6. Linear mixed-effects models
# ══════════════════════════════════════════════════════════════
#
# Main:
#   value ~ group + (1 | uniprot_id)
#
# Length-adjusted sensitivity analysis for absolute measures:
#   value ~ group + var_length + (1 | uniprot_id)
#
# As in the PARSE v2 workflow, only VAR metrics are modeled.
# WT is retained for descriptive VAR/WT/Δ summaries.
#
# Estimate = Disease − Control.

fit_lmm <- function(d, v, g_dis, g_con, fam,
                    region_name, adjust_len = FALSE) {
  
  col <- var_col(v)
  
  dd <- d %>%
    select(
      uniprot_id,
      group,
      value = all_of(col),
      var_length
    ) %>%
    filter(
      group %in% c(g_dis, g_con),
      !is.na(value),
      is.finite(value),
      !is.na(uniprot_id)
    ) %>%
    mutate(
      group = factor(group, levels = c(g_con, g_dis))
    )
  
  if (
    nrow(dd) < 20 ||
    n_distinct(dd$uniprot_id) < 3 ||
    n_distinct(dd$group) < 2
  ) {
    return(NULL)
  }
  
  if (adjust_len) {
    dd <- dd %>% filter(!is.na(var_length), is.finite(var_length))
    form <- value ~ group + var_length + (1 | uniprot_id)
  } else {
    form <- value ~ group + (1 | uniprot_id)
  }
  
  m <- tryCatch(
    lmerTest::lmer(form, data = dd),
    error = function(e) NULL
  )
  
  if (is.null(m)) return(NULL)
  
  cf <- coef(summary(m))
  r <- grep("^group", rownames(cf))[1]
  
  if (is.na(r)) return(NULL)
  
  est <- cf[r, "Estimate"]
  se  <- cf[r, "Std. Error"]
  
  tibble(
    region = region_name,
    family = fam,
    var = v,
    model = if (adjust_len) "LMM + length" else "LMM",
    estimate = est,
    se = se,
    ci_lo = est - 1.96 * se,
    ci_hi = est + 1.96 * se,
    p_raw = cf[r, "Pr(>|t|)"],
    sd_value = sd(dd$value, na.rm = TRUE),
    n_obs = nrow(dd),
    n_genes = n_distinct(dd$uniprot_id)
  )
}


run_region <- function(d, region_name) {
  
  pairs <- list(
    SNV = c("SNV Disease", "SNV Control"),
    FS  = c("FS Disease",  "FS Control")
  )
  
  pairs <- pairs[
    map_lgl(pairs, ~ all(.x %in% levels(d$group)))
  ]
  
  absolute_vars <- var_info$var[var_info$type == "absolute"]
  
  bind_rows(
    # Main LMM: all metrics
    imap_dfr(
      pairs,
      function(gp, fam) {
        map_dfr(
          var_info$var,
          ~ fit_lmm(
            d = d,
            v = .x,
            g_dis = gp[1],
            g_con = gp[2],
            fam = fam,
            region_name = region_name,
            adjust_len = FALSE
          )
        )
      }
    ),
    
    # Sensitivity analysis:
    # length-dependent outcomes adjusted for peptide length.
    # Length itself is not adjusted for itself.
    imap_dfr(
      pairs,
      function(gp, fam) {
        map_dfr(
          setdiff(absolute_vars, "length"),
          ~ fit_lmm(
            d = d,
            v = .x,
            g_dis = gp[1],
            g_con = gp[2],
            fam = fam,
            region_name = region_name,
            adjust_len = TRUE
          )
        )
      }
    )
  )
}


res_all <- bind_rows(
  run_region(df_div, "Divergent"),
  run_region(df_nmd, "NMD")
) %>%
  group_by(family, model) %>%
  mutate(
    q_fdr = p.adjust(p_raw, method = "BH")
  ) %>%
  ungroup() %>%
  left_join(
    var_info %>% select(var, label, label_disp, type),
    by = "var"
  ) %>%
  mutate(
    family = factor(family, levels = c("SNV", "FS")),
    region = factor(region, levels = c("Divergent", "NMD")),
    
    std_est = if_else(
      is.finite(sd_value) & sd_value > 0,
      estimate / sd_value,
      NA_real_
    ),
    std_lo = if_else(
      is.finite(sd_value) & sd_value > 0,
      ci_lo / sd_value,
      NA_real_
    ),
    std_hi = if_else(
      is.finite(sd_value) & sd_value > 0,
      ci_hi / sd_value,
      NA_real_
    ),
    
    sig = case_when(
      is.na(q_fdr) ~ "",
      q_fdr < 0.001 ~ "***",
      q_fdr < 0.01  ~ "**",
      q_fdr < FDR_CUTOFF ~ "*",
      TRUE ~ ""
    )
  )

write_csv(res_all, "CIDER_model_results.csv")


# ══════════════════════════════════════════════════════════════
# 7. Key comparison:
#    absolute measures vs normalized measures vs length adjustment
# ══════════════════════════════════════════════════════════════

cat("\n=== Key comparison: length effect vs composition/pattern effect ===\n")

res_all %>%
  filter(
    model %in% c("LMM", "LMM + length")
  ) %>%
  select(
    region, family, type, var, model,
    estimate, q_fdr, sig
  ) %>%
  mutate(
    across(
      c(estimate, q_fdr),
      ~ signif(.x, 3)
    )
  ) %>%
  arrange(region, family, type, var, model) %>%
  print(n = Inf)

cat("\n=== Significant main-model results ===\n")

res_all %>%
  filter(
    model == "LMM",
    q_fdr < FDR_CUTOFF
  ) %>%
  select(
    region, family, label, type,
    estimate, ci_lo, ci_hi, p_raw, q_fdr,
    n_obs, n_genes
  ) %>%
  mutate(
    across(
      c(estimate, ci_lo, ci_hi, p_raw, q_fdr),
      ~ signif(.x, 3)
    )
  ) %>%
  arrange(q_fdr) %>%
  print(n = Inf)


# ══════════════════════════════════════════════════════════════
# 8. GT result tables
# ══════════════════════════════════════════════════════════════

fmt_p <- function(x) {
  ifelse(
    is.na(x),
    "—",
    ifelse(
      x < 0.001,
      formatC(x, format = "e", digits = 1),
      formatC(x, format = "f", digits = 3)
    )
  )
}

make_tbl <- function(res, model_name, title_text, subtitle_text) {
  
  d <- res %>% filter(model == model_name)
  
  if (nrow(d) == 0) return(NULL)
  
  d %>%
    mutate(
      ci = sprintf(
        "(%s, %s)",
        signif(ci_lo, 3),
        signif(ci_hi, 3)
      )
    ) %>%
    select(
      region,
      label_disp,
      family,
      estimate,
      ci,
      p_raw,
      q_fdr,
      sig
    ) %>%
    pivot_wider(
      names_from = family,
      values_from = c(estimate, ci, p_raw, q_fdr, sig),
      names_sep = "_"
    ) %>%
    arrange(region, label_disp) %>%
    select(
      region,
      label_disp,
      any_of(c(
        "estimate_SNV", "ci_SNV", "p_raw_SNV", "q_fdr_SNV", "sig_SNV",
        "estimate_FS",  "ci_FS",  "p_raw_FS",  "q_fdr_FS",  "sig_FS"
      ))
    ) %>%
    gt(groupname_col = "region") %>%
    fmt_number(
      columns = any_of(c("estimate_SNV", "estimate_FS")),
      n_sigfig = 3
    ) %>%
    fmt(
      columns = any_of(c(
        "p_raw_SNV", "q_fdr_SNV",
        "p_raw_FS", "q_fdr_FS"
      )),
      fns = fmt_p
    ) %>%
    sub_missing(missing_text = "—") %>%
    cols_label(
      label_disp = "Variable",
      estimate_SNV = "Estimate",
      ci_SNV = "95% CI",
      p_raw_SNV = "p",
      q_fdr_SNV = "q (FDR)",
      sig_SNV = "",
      estimate_FS = "Estimate",
      ci_FS = "95% CI",
      p_raw_FS = "p",
      q_fdr_FS = "q (FDR)",
      sig_FS = ""
    ) %>%
    tab_spanner(
      label = md("**SNV**"),
      columns = any_of(c(
        "estimate_SNV", "ci_SNV", "p_raw_SNV",
        "q_fdr_SNV", "sig_SNV"
      ))
    ) %>%
    tab_spanner(
      label = md("**FS**"),
      columns = any_of(c(
        "estimate_FS", "ci_FS", "p_raw_FS",
        "q_fdr_FS", "sig_FS"
      ))
    ) %>%
    tab_header(
      title = title_text,
      subtitle = subtitle_text
    ) %>%
    tab_source_note(
      md(
        paste0(
          "Random intercept for gene (UniProt ID) accounts for multiple ",
          "variants per gene. q = Benjamini-Hochberg FDR within each ",
          "comparison and model family. ",
          "\\* q<0.05, \\*\\* q<0.01, \\*\\*\\* q<0.001. ",
          "† absolute/length-dependent measure. ",
          "‡ composition/fraction/pattern measure. ",
          "Divergent region: SNV excluded. Peptides < ",
          MIN_PEPTIDE_LEN,
          " aa excluded."
        )
      )
    )
}

tbl_lmm <- make_tbl(
  res_all,
  "LMM",
  "CIDER — linear mixed-effects models",
  "Estimate = Disease − Control, original units"
)

tbl_adj <- make_tbl(
  res_all,
  "LMM + length",
  "CIDER — length-adjusted sensitivity analysis",
  "value ~ group + length + (1 | gene); Estimate = Disease − Control"
)

tbl_lmm
tbl_adj

if (!is.null(tbl_lmm)) {
  gtsave(tbl_lmm, "CIDER_table_LMM.html")
  try(gtsave(tbl_lmm, "CIDER_table_LMM.docx"), silent = TRUE)
}

if (!is.null(tbl_adj)) {
  gtsave(tbl_adj, "CIDER_table_LMM_length_adjusted.html")
  try(gtsave(tbl_adj, "CIDER_table_LMM_length_adjusted.docx"), silent = TRUE)
}


# ══════════════════════════════════════════════════════════════
# 9. Forest plot
# ══════════════════════════════════════════════════════════════
res_all3 = res_all %>%
  filter(region != "NMD") 
p_forest <- res_all %>%
  filter(
    model == "LMM",
    !is.na(std_est),
    is.finite(std_est)
  ) %>%
  mutate(
    label_disp = factor(
      label_disp,
      levels = rev(var_info$label_disp)
    )
  ) %>%
  ggplot(
    aes(
      x = std_est,
      y = label_disp,
      color = family
    )
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray55"
  ) +
  geom_errorbarh(
    aes(
      xmin = std_lo,
      xmax = std_hi
    ),
    height = 0,
    linewidth = 0.6,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    aes(shape = q_fdr < FDR_CUTOFF),
    size = 2.4,
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    aes(label = sig),
    size = 3.2,
    hjust = -0.5,
    vjust = 0.8,
    position = position_dodge(width = 0.6),
    show.legend = FALSE
  ) +
  facet_wrap(~ region, nrow = 1) +
  scale_color_manual(
    values = c(
      SNV = "#2C7FB8",
      FS = "#C0392B"
    ),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(
      `TRUE` = 16,
      `FALSE` = 1
    ),
    labels = c(
      `TRUE` = paste0("q < ", FDR_CUTOFF),
      `FALSE` = "n.s."
    ),
    name = NULL
  ) +
  labs(
    title = "CIDER sequence properties: Disease vs Control",
    subtitle = "Mixed-effects estimates standardized by outcome SD; bars = 95% CI",
    x = "Standardized estimate (Disease − Control)",
    y = NULL,
    caption = paste0(
      "Linear mixed-effects model: value ~ group + (1 | gene). Filled points ",
      "BH-adjusted p < 0.05"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(
      fill = "gray95",
      color = NA
    ),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(
      hjust = 0,
      size = 8,
      color = "gray40"
    )
  )

print(p_forest)

ggsave(
  "CIDER_forest.pdf",
  p_forest,
  width = 11,
  height = 9,
  dpi = 300,
  device = cairo_pdf
)

ggsave(
  "CIDER_forest.png",
  p_forest,
  width = 11,
  height = 9,
  dpi = 300
)

