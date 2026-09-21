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
base_dir <- "~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/sequence_analysis/parse_v2"


# ══════════════════════════════════════════════════════════════
# 1. Read in and standardize column names
# ══════════════════════════════════════════════════════════════
# ...4/...5 = VAR, ...9/...10 = WT (by column order)

read_parse <- function(fname, len_var, len_wt) {
  read_csv(file.path(base_dir, fname), show_col_types = FALSE) %>%
    rename(
      var_ps_dist = `Σ classifier distance P...4`,
      var_ps_idr  = `length of longest PS IDR...5`,
      wt_ps_dist  = `Σ classifier distance P...9`,
      wt_ps_idr   = `length of longest PS IDR...10`,
      uniprot_id  = uniprot
    ) %>%
    rename(var_len = all_of(len_var), wt_len = all_of(len_wt))
}

parse_full <- read_parse("full_length_PARSE_v2.csv",
                         "var_length", "wt_length")
parse_div  <- read_parse("nmd_length_DivergentPos_PARSE_v2.csv",
                         "var_NMD_length", "wt_NMD_length")
parse_nmd  <- read_parse("nmd_length_NMDPos_PARSE_v2.csv",
                         "var_NMD_length_old", "wt_NMD_length_old")


# ══════════════════════════════════════════════════════════════
# 2. ID splitting, cohort filtering, density-normalized metrics
# ══════════════════════════════════════════════════════════════

make_key <- function(fn) {
  p <- str_split_fixed(fn, "_", 4)
  paste0(p[, 1], ":", as.numeric(p[, 2]), "|", p[, 3], "|", p[, 4])
}

prep_parse <- function(d, region_name, apply_len_filter = TRUE,
                       drop_snv = FALSE) {
  
  d <- d %>%
    separate(id, into = c("source_folder", "file_name"),
             sep = "-", extra = "merge", remove = FALSE) %>%
    mutate(
      key   = make_key(file_name),
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
  d <- bind_rows(
    d %>% filter(key %in% all_variants$key),
    d %>% filter(source_folder %in% c("snv_control", "fs_control"))
  )
  n1 <- nrow(d)
  
  d <- d %>% filter(!is.na(var_len))
  if (apply_len_filter) d <- d %>% filter(var_len >= MIN_PEPTIDE_LEN)
  n2 <- nrow(d)
  
  # Divergent region: nonsense SNVs yield no novel peptide (median length 0), excluded
  if (drop_snv) {
    d <- d %>% filter(str_starts(source_folder, "fs_"))
  }
  n3 <- nrow(d)
  
  cat(sprintf("[%-10s] raw %d -> cohort %d -> length %d -> after SNV exclusion %d\n",
              region_name, n0, n1, n2, n3))
  
  d %>%
    mutate(
      region = region_name,
      # ── Density normalization: distinguish shorter sequence from lower intrinsic PS propensity ──
      ps_dist_per100 = var_ps_dist / var_len * 100,
      idr_fraction   = var_ps_idr  / var_len,
      wt_ps_dist_per100 = wt_ps_dist / wt_len * 100,
      wt_idr_fraction   = wt_ps_idr  / wt_len,
      group = factor(group, levels = c("SNV Disease", "SNV Control",
                                       "FS Disease", "FS Control")) %>%
        droplevels()
    )
}

cat("=== Filtering Log ===\n")
df_full <- prep_parse(parse_full, "Full",      apply_len_filter = FALSE)
df_div  <- prep_parse(parse_div,  "Divergent", apply_len_filter = TRUE,
                      drop_snv = TRUE)
df_nmd  <- prep_parse(parse_nmd,  "NMD",       apply_len_filter = TRUE)


# ══════════════════════════════════════════════════════════════
# 3. Diagnostics: key match rate (disease group n is low, check for missed matches)
# ══════════════════════════════════════════════════════════════

cat("\n=== Key Match Diagnostics (Disease Group) ===\n")
key_diag <- imap_dfr(
  list(Full = parse_full, Divergent = parse_div, NMD = parse_nmd),
  function(d, nm) {
    d %>%
      separate(id, c("sf", "fn"), sep = "-", extra = "merge") %>%
      filter(str_detect(sf, "disease")) %>%
      mutate(matched = make_key(fn) %in% all_variants$key) %>%
      count(sf, matched) %>%
      pivot_wider(names_from = matched, values_from = n,
                  names_prefix = "matched_", values_fill = 0) %>%
      mutate(file = nm, .before = 1)
  }
)
print(as.data.frame(key_diag))

cat("\nUnmatched file_name examples (to check indel notation):\n")
parse_full %>%
  separate(id, c("sf", "fn"), sep = "-", extra = "merge") %>%
  filter(str_detect(sf, "disease"),
         !make_key(fn) %in% all_variants$key) %>%
  slice_head(n = 8) %>% pull(fn) %>% print()

cat("\n=== Sample Size by Region and Group ===\n")
bind_rows(df_full, df_div, df_nmd) %>%
  count(region, group) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
  as.data.frame() %>% print()

cat("\n=== Zero-Value Proportions ===\n")
bind_rows(df_full, df_div, df_nmd) %>%
  group_by(region, group) %>%
  summarise(n = n(),
            pct_zero_idr  = round(100 * mean(var_ps_idr  == 0), 1),
            pct_zero_dist = round(100 * mean(var_ps_dist == 0), 1),
            .groups = "drop") %>%
  as.data.frame() %>% print()
cat("\n")


# ══════════════════════════════════════════════════════════════
# 4. Variable definitions
# ══════════════════════════════════════════════════════════════
# Absolute measures (†) scale with sequence length; density metrics (‡) are normalized,
# key to whether phase-separation propensity itself has changed.

var_info <- tribble(
  ~var,              ~label,                             ~type,
  "len",             "Sequence Length (aa)",              "absolute",
  "ps_dist",         "Σ Classifier Distance (PS)",        "absolute",
  "ps_idr",          "Longest PS IDR (aa)",               "absolute",
  "ps_dist_per100",  "Σ Classifier Distance / 100 aa",    "density",
  "idr_fraction",    "Longest PS IDR / length",           "density"
) %>%
  mutate(label_disp = paste0(label, if_else(type == "absolute", "", "")))

# Column names used for modeling (density metrics lack var_ prefix)
col_of <- function(v) if_else(v %in% c("ps_dist_per100", "idr_fraction"),
                              v, paste0("var_", v))


# ══════════════════════════════════════════════════════════════
# 5. Descriptive tables
# ══════════════════════════════════════════════════════════════

make_desc <- function(d, region_name) {
  
  d2 <- d %>%
    mutate(
      delta_len          = var_len     - wt_len,
      delta_ps_dist      = var_ps_dist - wt_ps_dist,
      delta_ps_idr       = var_ps_idr  - wt_ps_idr,
      delta_dist_per100  = ps_dist_per100 - wt_ps_dist_per100,
      delta_idr_fraction = idr_fraction   - wt_idr_fraction,
      var_has_idr = factor(if_else(var_ps_idr > 0, "Yes", "No"),
                           levels = c("No", "Yes")),
      wt_has_idr  = factor(if_else(wt_ps_idr  > 0, "Yes", "No"),
                           levels = c("No", "Yes"))
    ) %>%
    select(group, var_len, wt_len, delta_len,
           var_ps_dist, wt_ps_dist, delta_ps_dist,
           var_ps_idr, wt_ps_idr, delta_ps_idr,
           ps_dist_per100, wt_ps_dist_per100, delta_dist_per100,
           idr_fraction, wt_idr_fraction, delta_idr_fraction,
           var_has_idr, wt_has_idr)
  
  grp <- c(
    var_len = "Length (†)", wt_len = "Length (†)", delta_len = "Length (†)",
    var_ps_dist = "Σ Classifier Distance (†)",
    wt_ps_dist = "Σ Classifier Distance (†)",
    delta_ps_dist = "Σ Classifier Distance (†)",
    var_ps_idr = "Longest PS IDR (†)", wt_ps_idr = "Longest PS IDR (†)",
    delta_ps_idr = "Longest PS IDR (†)",
    ps_dist_per100 = "Density: distance / 100 aa (‡)",
    wt_ps_dist_per100 = "Density: distance / 100 aa (‡)",
    delta_dist_per100 = "Density: distance / 100 aa (‡)",
    idr_fraction = "Density: IDR fraction (‡)",
    wt_idr_fraction = "Density: IDR fraction (‡)",
    delta_idr_fraction = "Density: IDR fraction (‡)",
    var_has_idr = "PS IDR present", wt_has_idr = "PS IDR present"
  )
  
  d2 %>%
    tbl_summary(
      by = group,
      statistic = list(all_continuous()  ~ "{median} ({p25}, {p75})",
                       all_categorical() ~ "{n} ({p}%)"),
      digits  = list(all_continuous() ~ 2,
                     c(idr_fraction, wt_idr_fraction,
                       delta_idr_fraction) ~ 3),
      missing = "no",
      label = list(
        var_len ~ "VAR Length (aa)", wt_len ~ "WT Length (aa)",
        delta_len ~ "Δ Length (aa)",
        var_ps_dist ~ "VAR Σ Classifier Distance",
        wt_ps_dist ~ "WT Σ Classifier Distance",
        delta_ps_dist ~ "Δ Σ Classifier Distance",
        var_ps_idr ~ "VAR Longest PS IDR (aa)",
        wt_ps_idr ~ "WT Longest PS IDR (aa)",
        delta_ps_idr ~ "Δ Longest PS IDR (aa)",
        ps_dist_per100 ~ "VAR Distance / 100 aa",
        wt_ps_dist_per100 ~ "WT Distance / 100 aa",
        delta_dist_per100 ~ "Δ Distance / 100 aa",
        idr_fraction ~ "VAR IDR fraction",
        wt_idr_fraction ~ "WT IDR fraction",
        delta_idr_fraction ~ "Δ IDR fraction",
        var_has_idr ~ "VAR has PS IDR", wt_has_idr ~ "WT has PS IDR"
      )
    ) %>%
    add_overall() %>%
    modify_table_body(~ .x %>% mutate(groupname_col = grp[variable])) %>%
    modify_header(label ~ "**Variable**") %>%
    modify_caption(paste0("**PARSE v2 — ", region_name,
                          " region (descriptive)**")) %>%
    bold_labels()
}

desc_full <- make_desc(df_full, "Full")
desc_div  <- make_desc(df_div,  "Divergent (FS only)")
desc_nmd  <- make_desc(df_nmd,  "NMD")

desc_full; desc_div; desc_nmd

desc_full %>% as_gt() %>% gtsave("PARSE_descriptive_Full.html")
desc_div  %>% as_gt() %>% gtsave("PARSE_descriptive_Divergent.html")
desc_nmd  %>% as_gt() %>% gtsave("PARSE_descriptive_NMD.html")


# ══════════════════════════════════════════════════════════════
# 6. Models
# ══════════════════════════════════════════════════════════════
# LMM          : value ~ group + (1 | gene)                 primary analysis
# LMM+length   : value ~ group + var_len + (1 | gene)       sensitivity (absolute only)
# GLMM         : (value > 0) ~ group + (1 | gene)           zero-inflation auxiliary
#
# Model VAR only: WT is constant within gene, absorbed by random intercept, so Delta model and VAR
# model give mathematically identical group effects.

fit_lmm <- function(d, v, g_dis, g_con, fam, region_name, adjust_len = FALSE) {
  
  col <- col_of(v)
  dd <- d %>%
    select(uniprot_id, group, value = all_of(col), var_len) %>%
    filter(group %in% c(g_dis, g_con), !is.na(value), is.finite(value)) %>%
    mutate(group = factor(group, levels = c(g_con, g_dis)))
  
  if (nrow(dd) < 20 || n_distinct(dd$uniprot_id) < 3 ||
      n_distinct(dd$group) < 2) return(NULL)
  
  form <- if (adjust_len) value ~ group + var_len + (1 | uniprot_id)
  else            value ~ group + (1 | uniprot_id)
  
  m <- tryCatch(lmerTest::lmer(form, data = dd), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  
  cf <- coef(summary(m)); r <- grep("^group", rownames(cf))[1]
  est <- cf[r, "Estimate"]; se <- cf[r, "Std. Error"]
  
  tibble(region = region_name, family = fam, var = v,
         model = if (adjust_len) "LMM + length" else "LMM",
         estimate = est, se = se,
         ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
         p_raw = cf[r, "Pr(>|t|)"], sd_value = sd(dd$value),
         n_obs = nrow(dd), n_genes = n_distinct(dd$uniprot_id))
}

fit_glmm <- function(d, v, g_dis, g_con, fam, region_name) {
  
  col <- col_of(v)
  dd <- d %>%
    select(uniprot_id, group, value = all_of(col)) %>%
    filter(group %in% c(g_dis, g_con), !is.na(value)) %>%
    mutate(group = factor(group, levels = c(g_con, g_dis)),
           y = as.integer(value > 0))
  
  if (nrow(dd) < 20 || n_distinct(dd$y) < 2 ||
      n_distinct(dd$uniprot_id) < 3) return(NULL)
  
  m <- tryCatch(
    glmer(y ~ group + (1 | uniprot_id), data = dd, family = binomial,
          control = glmerControl(optimizer = "bobyqa")),
    error = function(e) NULL, warning = function(w) NULL)
  if (is.null(m)) return(NULL)
  
  cf <- coef(summary(m)); r <- grep("^group", rownames(cf))[1]
  est <- cf[r, "Estimate"]; se <- cf[r, "Std. Error"]
  
  tibble(region = region_name, family = fam, var = v,
         model = "GLMM (presence)", estimate = est, se = se,
         ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
         p_raw = cf[r, "Pr(>|z|)"], sd_value = NA_real_,
         n_obs = nrow(dd), n_genes = n_distinct(dd$uniprot_id))
}

run_region <- function(d, region_name) {
  
  pairs <- list(SNV = c("SNV Disease", "SNV Control"),
                FS  = c("FS Disease",  "FS Control"))
  pairs <- pairs[map_lgl(pairs, ~ all(.x %in% levels(d$group)))]
  
  abs_vars <- var_info$var[var_info$type == "absolute"]
  den_vars <- var_info$var[var_info$type == "density"]
  
  bind_rows(
    # Primary analysis: all variables
    imap_dfr(pairs, function(gp, fam)
      map_dfr(var_info$var, ~ fit_lmm(d, .x, gp[1], gp[2], fam, region_name))),
    # Sensitivity: absolute measures plus length covariate (length itself excluded)
    imap_dfr(pairs, function(gp, fam)
      map_dfr(setdiff(abs_vars, "len"),
              ~ fit_lmm(d, .x, gp[1], gp[2], fam, region_name,
                        adjust_len = TRUE))),
    # Zero-inflation auxiliary
    imap_dfr(pairs, function(gp, fam)
      map_dfr(c("ps_idr", "ps_dist"),
              ~ fit_glmm(d, .x, gp[1], gp[2], fam, region_name)))
  )
}

res_all <- bind_rows(
  run_region(df_full, "Full"),
  run_region(df_div,  "Divergent"),
  run_region(df_nmd,  "NMD")
) %>%
  group_by(family, model) %>%            # independent correction within each comparison and model family
  mutate(q_fdr = p.adjust(p_raw, method = "BH")) %>%
  ungroup() %>%
  left_join(var_info %>% select(var, label, label_disp, type), by = "var") %>%
  mutate(
    family  = factor(family, levels = c("SNV", "FS")),
    region  = factor(region, levels = c("Full", "Divergent", "NMD")),
    std_est = if_else(model == "GLMM (presence)", estimate, estimate / sd_value),
    std_lo  = if_else(model == "GLMM (presence)", ci_lo,    ci_lo    / sd_value),
    std_hi  = if_else(model == "GLMM (presence)", ci_hi,    ci_hi    / sd_value),
    sig = case_when(
      is.na(q_fdr)  ~ "",
      q_fdr < 0.001 ~ "***",
      q_fdr < 0.01  ~ "**",
      q_fdr < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

write_csv(res_all, "PARSE_model_results.csv")


# ══════════════════════════════════════════════════════════════
# 7. Key comparison: absolute vs density vs length-adjusted
# ══════════════════════════════════════════════════════════════
# Density metric or length-adjusted still significant -> intrinsic PS propensity changed
# Only absolute measure significant -> pure truncation effect

cat("=== Key Comparison: Truncation Effect vs Intrinsic Propensity ===\n")
res_all %>%
  filter(model %in% c("LMM", "LMM + length"),
         var %in% c("ps_dist", "ps_idr", "ps_dist_per100", "idr_fraction")) %>%
  select(region, family, var, model, estimate, q_fdr, sig) %>%
  mutate(across(c(estimate, q_fdr), ~ signif(.x, 3))) %>%
  arrange(region, family, var, model) %>%
  print(n = Inf)

cat("\n=== All Model Results ===\n")
res_all %>%
  select(region, family, model, label, estimate, p_raw, q_fdr, sig,
         n_obs, n_genes) %>%
  mutate(across(c(estimate, p_raw, q_fdr), ~ signif(.x, 3))) %>%
  arrange(model, region, family, q_fdr) %>%
  print(n = Inf)


# ══════════════════════════════════════════════════════════════
# 8. Tables
# ══════════════════════════════════════════════════════════════

fmt_p <- function(x) {
  ifelse(is.na(x), "—",
         ifelse(x < 0.001, formatC(x, format = "e", digits = 1),
                formatC(x, format = "f", digits = 3)))
}

make_tbl <- function(res, model_name, title_text, subtitle_text) {
  
  d <- res %>% filter(model == model_name)
  if (nrow(d) == 0) return(NULL)
  
  d %>%
    mutate(ci = sprintf("(%s, %s)", signif(ci_lo, 3), signif(ci_hi, 3))) %>%
    select(region, label_disp, family, estimate, ci, p_raw, q_fdr, sig) %>%
    pivot_wider(names_from = family,
                values_from = c(estimate, ci, p_raw, q_fdr, sig),
                names_sep = "_") %>%
    arrange(region) %>%
    select(region, label_disp,
           any_of(c("estimate_SNV","ci_SNV","p_raw_SNV","q_fdr_SNV","sig_SNV",
                    "estimate_FS","ci_FS","p_raw_FS","q_fdr_FS","sig_FS"))) %>%
    gt(groupname_col = "region") %>%
    fmt_number(columns = any_of(c("estimate_SNV","estimate_FS")),
               n_sigfig = 3) %>%
    fmt(columns = any_of(c("p_raw_SNV","q_fdr_SNV","p_raw_FS","q_fdr_FS")),
        fns = fmt_p) %>%
    sub_missing(missing_text = "—") %>%
    cols_label(
      label_disp = "Variable",
      estimate_SNV = "Estimate", ci_SNV = "95% CI",
      p_raw_SNV = "p", q_fdr_SNV = "q (FDR)", sig_SNV = "",
      estimate_FS = "Estimate", ci_FS = "95% CI",
      p_raw_FS = "p", q_fdr_FS = "q (FDR)", sig_FS = ""
    ) %>%
    tab_spanner(label = md("**SNV**"),
                columns = any_of(c("estimate_SNV","ci_SNV","p_raw_SNV",
                                   "q_fdr_SNV","sig_SNV"))) %>%
    tab_spanner(label = md("**FS**"),
                columns = any_of(c("estimate_FS","ci_FS","p_raw_FS",
                                   "q_fdr_FS","sig_FS"))) %>%
    tab_header(title = title_text, subtitle = subtitle_text) %>%
    tab_source_note(md(paste0(
      "Random intercept for gene (UniProt ID) accounts for multiple variants ",
      "per gene. q: Benjamini-Hochberg FDR within each comparison and model ",
      "family. \\* q<0.05, \\*\\* q<0.01, \\*\\*\\* q<0.001. ",
      "† absolute measure — scales with sequence length. ",
      "‡ density measure — normalized to length; these test whether ",
      "phase-separation propensity per residue differs, independent of how ",
      "much sequence was lost. ",
      "Divergent region: SNV excluded (nonsense variants produce no divergent ",
      "peptide; median length 0). Region peptides < ", MIN_PEPTIDE_LEN,
      " aa excluded. WT metrics are constant within gene and were not modeled."
    )))
}

tbl_lmm  <- make_tbl(res_all, "LMM",
                     "PARSE v2 — linear mixed-effects models",
                     "Estimate = Disease − Control, original units")
tbl_adj  <- make_tbl(res_all, "LMM + length",
                     "PARSE v2 — length-adjusted sensitivity analysis",
                     "value ~ group + length + (1 | gene); Estimate = Disease − Control")
tbl_glmm <- make_tbl(res_all, "GLMM (presence)",
                     "PARSE v2 — presence of a non-zero score (logistic GLMM)",
                     "Estimate = log-odds ratio, Disease vs Control")

tbl_lmm; tbl_adj; tbl_glmm

gtsave(tbl_lmm,  "PARSE_table_LMM.html")
gtsave(tbl_adj,  "PARSE_table_LMM_length_adjusted.html")
gtsave(tbl_glmm, "PARSE_table_GLMM.html")
try(gtsave(tbl_lmm,  "PARSE_table_LMM.docx"), silent = TRUE)
try(gtsave(tbl_adj,  "PARSE_table_LMM_length_adjusted.docx"), silent = TRUE)
try(gtsave(tbl_glmm, "PARSE_table_GLMM.docx"), silent = TRUE)


# ══════════════════════════════════════════════════════════════
# 9. Forest plot: absolute and density metrics side by side
# ══════════════════════════════════════════════════════════════
res_all2 <- res_all %>%
  filter(region != "NMD") 
p_forest <- res_all %>%
  filter(model == "LMM", !is.na(estimate)) %>%
  mutate(label_disp = factor(label_disp, levels = rev(var_info$label_disp))) %>%
  ggplot(aes(x = std_est, y = label_disp, color = family)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55") +
  geom_errorbarh(aes(xmin = std_lo, xmax = std_hi), height = 0,
                 linewidth = 0.6, position = position_dodge(width = 0.6)) +
  geom_point(aes(shape = q_fdr < 0.05), size = 2.4,
             position = position_dodge(width = 0.6)) +
  geom_text(aes(label = sig), size = 3.2, hjust = -0.5, vjust = 0.8,
            position = position_dodge(width = 0.6), show.legend = FALSE) +
  facet_wrap(~ region, nrow = 1) +
  scale_color_manual(values = c(SNV = "#2C7FB8", FS = "#C0392B"), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."),
                     name = NULL) +
  labs(
    title    = "PARSE v2 phase-separation metrics: Disease vs Control",
    subtitle = "Mixed-effects estimates standardized by outcome SD; bars = 95% CI",
    x = "Standardized estimate (Disease − Control)", y = NULL,
    caption = paste0(
      "Linear mixed-effects model: value ~ group + (1 | gene). Filled points ",
      "BH-adjusted p < 0.05"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background   = element_rect(fill = "gray95", color = NA),
    strip.text         = element_text(face = "bold"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position    = "top",
    plot.title         = element_text(face = "bold"),
    plot.caption       = element_text(hjust = 0, size = 8, color = "gray40")
  )

print(p_forest)
ggsave("PARSE_forest.pdf", p_forest, width = 11, height = 8,
       dpi = 300, device = cairo_pdf)
ggsave("PARSE_forest.png", p_forest, width = 11, height = 8, dpi = 300)

