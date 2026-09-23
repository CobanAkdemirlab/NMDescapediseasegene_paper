# =============================================================================
# Repeat & Homopolymer Content — Paired Analysis (4 Gene Groups)
#
# Compares within matched pairs:
#   * CDS repeat fraction / NMD-escape-region repeat fraction
#   * CDS homopolymer fraction / NMD-escape-region homopolymer fraction
#   * Within-gene NMD-vs-CDS difference (escape - whole CDS), per metric
#
# Inputs:
#   matching_log_v2.csv     <- study_gene, ctrl_gene, group ("FS"/"Stopgain"),
#                              study_transcript, ctrl_transcript
#   repeat_content.xlsx     <- ensembl_transcript_id, group, repeat_fraction,
#                              nmdesc_repeat_fraction, homopolymer_fraction,
#                              nmdesc_homopolymer_fraction
#
# Outputs:
#   repeat_homopolymer_paired.pdf / .png    <- main 6-panel figure
#   repeat_homopolymer_scatter.pdf / .png   <- CDS-vs-NMD scatter, faceted
#   repeat_homopolymer_stats.csv            <- paired Wilcoxon stats
#
# Test: paired Wilcoxon signed-rank, BH-adjusted within metric.
# Style: unified box + jittered points + paired connecting lines; pale
#        per-class shaded backgrounds; dashed zero line on the Δ panels.
# =============================================================================

suppressPackageStartupMessages({
  cran <- c("ggplot2","dplyr","readr","tidyr","forcats","patchwork","scales")
  for (p in cran)
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(ggplot2); library(dplyr); library(readr); library(tidyr)
  library(forcats); library(patchwork); library(scales)
})

# Robust .xlsx reader: try readxl, fall back to openxlsx. Avoids a hard
# dependency on one Excel backend (readxl can fail to load on some systems).
read_xlsx_robust <- function(path) {
  has_readxl <- isTRUE(tryCatch(requireNamespace("readxl", quietly = TRUE),
                                error = function(e) FALSE))
  if (has_readxl) {
    out <- tryCatch(readxl::read_excel(path), error = function(e) NULL)
    if (!is.null(out)) return(as.data.frame(out))
  }
  if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
  openxlsx::read.xlsx(path)
}

# ── 1. Load ─────────────────────────────────────────────────────────────────
cat("Reading inputs...\n")
# Single-file input: gene_flags_per_gene.csv already contains the matched pairs
# (design / pair_id / stratum / role) AND the per-gene repeat & homopolymer
# metrics, so no separate pairs file or repeat_content.xlsx is needed.
IN <- "gene_flags_per_gene.csv"
metric_cols <- c("repeat_fraction", "nmdesc_repeat_fraction",
                 "homopolymer_fraction", "nmdesc_homopolymer_fraction")

raw <- read_csv(IN, show_col_types = FALSE) %>% filter(design == "matched")
stopifnot(all(c("pair_id","stratum","role", metric_cols) %in% names(raw)))

mkside <- function(df, role_val, prefix)
  df %>% filter(role == role_val) %>%
    dplyr::select(pair_id, stratum, dplyr::all_of(metric_cols)) %>%
    dplyr::rename_with(~ paste0(prefix, "_", .x), dplyr::all_of(metric_cols))

pp <- mkside(raw, "case", "study") %>%
  inner_join(mkside(raw, "control", "ctrl"), by = c("pair_id", "stratum")) %>%
  mutate(
    group = dplyr::recode(stratum, FS = "FS", SNV = "Stopgain"),
    study_rep_diff = study_nmdesc_repeat_fraction      - study_repeat_fraction,
    ctrl_rep_diff  = ctrl_nmdesc_repeat_fraction       - ctrl_repeat_fraction,
    study_hp_diff  = study_nmdesc_homopolymer_fraction - study_homopolymer_fraction,
    ctrl_hp_diff   = ctrl_nmdesc_homopolymer_fraction  - ctrl_homopolymer_fraction)
pp$pair_id <- seq_len(nrow(pp))
cat(sprintf("Pairs: FS=%d, Stopgain=%d\n",
            sum(pp$group == "FS"), sum(pp$group == "Stopgain")))

# ── 2. Plotting metadata ────────────────────────────────────────────────────
grp_levels <- c("Frameshift","FS Control","Stopgain","Stopgain Control")
grp_cols   <- c("Frameshift"       = "#2166AC",
                "FS Control"       = "#92C5DE",
                "Stopgain"         = "#B2182B",
                "Stopgain Control" = "#F4A582")
# pale per-class backgrounds (left pair = blue, right pair = red)
bg_fs  <- "#2166AC"; bg_sv <- "#B2182B"; bg_alpha <- 0.06

build_long2 <- function(pp, study_col, ctrl_col, metric_label) {
  fs  <- pp %>% filter(group == "FS")
  snv <- pp %>% filter(group == "Stopgain")
  bind_rows(
    tibble(pair_id = fs$pair_id,  Group = "Frameshift",
           value = fs[[study_col]],  metric = metric_label),
    tibble(pair_id = fs$pair_id,  Group = "FS Control",
           value = fs[[ctrl_col]],   metric = metric_label),
    tibble(pair_id = snv$pair_id, Group = "Stopgain",
           value = snv[[study_col]], metric = metric_label),
    tibble(pair_id = snv$pair_id, Group = "Stopgain Control",
           value = snv[[ctrl_col]],  metric = metric_label)
  ) %>%
    mutate(Group = factor(Group, levels = grp_levels))
}

# ── 3. Paired Wilcoxon helper + stats table ─────────────────────────────────
paired_w <- function(s, c) {
  ok <- !is.na(s) & !is.na(c); s <- s[ok]; c <- c[ok]
  if (length(s) < 3 || all(s == c))
    return(list(p = NA_real_, n = length(s),
                med_s = if (length(s)) median(s) else NA,
                med_c = if (length(c)) median(c) else NA))
  pv <- suppressWarnings(wilcox.test(s, c, paired = TRUE, exact = FALSE)$p.value)
  list(p = pv, n = length(s), med_s = median(s), med_c = median(c))
}

stats_rows <- list(); k <- 0
for (mt in list(
  c("CDS repeat",            "study_repeat_fraction",            "ctrl_repeat_fraction"),
  c("NMD-escape repeat",     "study_nmdesc_repeat_fraction",     "ctrl_nmdesc_repeat_fraction"),
  c("Repeat NMD-CDS \u0394", "study_rep_diff",                   "ctrl_rep_diff"),
  c("CDS homopolymer",       "study_homopolymer_fraction",       "ctrl_homopolymer_fraction"),
  c("NMD-escape homopolymer","study_nmdesc_homopolymer_fraction","ctrl_nmdesc_homopolymer_fraction"),
  c("Homopolymer NMD-CDS \u0394","study_hp_diff",                "ctrl_hp_diff"))) {
  for (cmp in list(c("Frameshift", "FS"), c("Stopgain", "Stopgain"))) {
    sub <- pp %>% filter(group == cmp[2])
    res <- paired_w(sub[[mt[2]]], sub[[mt[3]]]); k <- k + 1
    stats_rows[[k]] <- tibble(metric = mt[1], comparison = sprintf("%s vs control", cmp[1]),
                              n_pairs = res$n, med_case = res$med_s, med_ctrl = res$med_c, p_raw = res$p)
  }
}
stats_df <- bind_rows(stats_rows) %>%
  group_by(metric) %>% mutate(p_BH = p.adjust(p_raw, method = "BH")) %>% ungroup() %>%
  mutate(signif = case_when(is.na(p_BH) ~ "NA", p_BH < 0.001 ~ "***",
                            p_BH < 0.01 ~ "**", p_BH < 0.05 ~ "*", TRUE ~ "ns"))
write_csv(stats_df, "repeat_homopolymer_stats.csv")
cat("\n=== Paired Wilcoxon (BH within metric) ===\n"); print(as.data.frame(stats_df))

# ── 4. Theme ────────────────────────────────────────────────────────────────
theme_big <- theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
        plot.subtitle = element_text(size = 11.5, hjust = 0.5, color = "grey30"),
        plot.tag = element_text(face = "bold", size = 19),
        plot.tag.position = c(0.04, 0.98),
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text = element_text(size = 13, color = "black"),
        axis.line = element_line(linewidth = 0.9),
        axis.ticks = element_line(linewidth = 0.9),
        panel.border = element_rect(linewidth = 0.9, color = "grey30"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")

# Device-robust saver. On macOS use the native quartz device (renders Unicode,
# needs no XQuartz/cairo); elsewhere use cairo. Avoids "failed to load cairo DLL".
save_fig <- function(path, plot, width, height, dpi = 200) {
  is_pdf <- grepl("\\.pdf$", path, ignore.case = TRUE)
  if (isTRUE(capabilities("aqua"))) {                 # macOS, no XQuartz needed
    if (is_pdf) grDevices::quartz(file = path, type = "pdf", width = width, height = height)
    else        grDevices::quartz(file = path, type = "png", width = width, height = height, dpi = dpi)
    print(plot); grDevices::dev.off()
  } else if (is_pdf) {                                 # Linux / Windows
    ggsave(path, plot, width = width, height = height, device = cairo_pdf)
  } else {
    ggsave(path, plot, width = width, height = height, dpi = dpi, type = "cairo")
  }
  invisible(path)
}

# ── 5. ONE unified panel: shaded bg + paired lines + box + jitter + brackets ─
#      (zero_line = TRUE adds the dashed 0 reference used by the Δ panels)
p_panel <- function(study_col, ctrl_col, metric_label, ylab,
                    title_str, sub_str, tag_letter,
                    percent_y = TRUE, zero_line = FALSE) {
  d <- build_long2(pp, study_col, ctrl_col, metric_label)

  # pair-connecting segments (study <-> matched control)
  fs_lines <- pp %>% filter(group == "FS") %>%
    transmute(pair_id, x_start = "Frameshift", x_end = "FS Control",
              y_start = .data[[study_col]], y_end = .data[[ctrl_col]])
  sv_lines <- pp %>% filter(group == "Stopgain") %>%
    transmute(pair_id, x_start = "Stopgain", x_end = "Stopgain Control",
              y_start = .data[[study_col]], y_end = .data[[ctrl_col]])
  lines <- bind_rows(fs_lines, sv_lines) %>%
    mutate(xs = match(x_start, grp_levels), xe = match(x_end, grp_levels))

  fs_st <- paired_w(pp[[study_col]][pp$group == "FS"],       pp[[ctrl_col]][pp$group == "FS"])
  sv_st <- paired_w(pp[[study_col]][pp$group == "Stopgain"], pp[[ctrl_col]][pp$group == "Stopgain"])
  fmt_p <- function(p) if (is.na(p)) "p = n/a" else if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)

  ymax <- max(d$value, na.rm = TRUE); ymin <- min(d$value, na.rm = TRUE)
  yspan <- ymax - ymin; yb <- ymax + yspan * 0.10
  brackets <- tibble(x1 = c(1, 3), x2 = c(2, 4), y = yb,
                     lab = c(fmt_p(fs_st$p), fmt_p(sv_st$p)),
                     sig = c(!is.na(fs_st$p) && fs_st$p < 0.05,
                             !is.na(sv_st$p) && sv_st$p < 0.05))

  d$xn <- as.numeric(d$Group)   # numeric positions -> no discrete/continuous mixing
  p <- ggplot(d, aes(x = xn, y = value, fill = Group, colour = Group, group = Group)) +
    # pale per-class background panels
    annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
             fill = bg_fs, alpha = bg_alpha) +
    annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
             fill = bg_sv, alpha = bg_alpha)
  if (zero_line)
    p <- p + geom_hline(yintercept = 0, colour = "grey50",
                        linetype = "dashed", linewidth = 0.6)
  p <- p +
    geom_segment(data = lines, inherit.aes = FALSE,
                 aes(x = xs, xend = xe, y = y_start, yend = y_end,
                     group = pair_id),
                 colour = "grey60", linewidth = 0.4, alpha = 0.45) +
    geom_boxplot(aes(group = Group), width = 0.42, outlier.shape = NA,
                 linewidth = 1.0, alpha = 0.55, fatten = 1.4) +
    geom_jitter(width = 0.08, size = 1.6, alpha = 0.55, shape = 16) +
    geom_segment(data = brackets, inherit.aes = FALSE,
                 aes(x = x1, xend = x2, y = y, yend = y),
                 colour = "grey20", linewidth = 0.7) +
    geom_segment(data = brackets, inherit.aes = FALSE,
                 aes(x = x1, xend = x1, y = y, yend = y - yspan*0.015),
                 colour = "grey20", linewidth = 0.7) +
    geom_segment(data = brackets, inherit.aes = FALSE,
                 aes(x = x2, xend = x2, y = y, yend = y - yspan*0.015),
                 colour = "grey20", linewidth = 0.7) +
    geom_text(data = brackets, inherit.aes = FALSE,
              aes(x = (x1 + x2)/2, y = y + yspan*0.025, label = lab,
                  fontface = ifelse(sig, "bold", "plain")),
              size = 4.4, colour = "grey10") +
    scale_fill_manual(values = grp_cols) +
    scale_colour_manual(values = grp_cols) +
    scale_x_continuous(breaks = 1:4, labels = grp_levels, limits = c(0.5, 4.5)) +
    labs(title = title_str, subtitle = sub_str, x = NULL, y = ylab, tag = tag_letter) +
    theme_big
  if (percent_y) p <- p + scale_y_continuous(labels = label_percent(accuracy = 0.1))
  p
}

# ── 6. Build the six panels (C & F = Δ panels with zero line) ────────────────
MINUS <- "\u2212"   # true minus sign, to match "NMD-Escape − CDS"
p1 <- p_panel("study_repeat_fraction", "ctrl_repeat_fraction",
              "CDS repeat", "Repeat fraction",
              "CDS Repeat Fraction",
              "Fraction of full CDS in repeat elements", "A")
p2 <- p_panel("study_nmdesc_repeat_fraction", "ctrl_nmdesc_repeat_fraction",
              "NMD-esc repeat", "Repeat fraction",
              "NMD-Escape Repeat Fraction",
              "Fraction of NMD-escape region in repeat elements", "B")
p3 <- p_panel("study_rep_diff", "ctrl_rep_diff",
              "Repeat \u0394", paste0("NMD-escape ", MINUS, " CDS repeat fraction"),
              paste0("Repeat \u0394 (NMD-Escape ", MINUS, " CDS)"),
              "Positive = NMD-escape region repeat-richer than full CDS", "C",
              zero_line = TRUE)
p4 <- p_panel("study_homopolymer_fraction", "ctrl_homopolymer_fraction",
              "CDS homopolymer", "Homopolymer fraction",
              "CDS Homopolymer Fraction",
              "Fraction of full CDS in homopolymer runs", "D")
p5 <- p_panel("study_nmdesc_homopolymer_fraction", "ctrl_nmdesc_homopolymer_fraction",
              "NMD-esc homopolymer", "Homopolymer fraction",
              "NMD-Escape Homopolymer Fraction",
              "Fraction of NMD-escape region in homopolymer runs", "E")
p6 <- p_panel("study_hp_diff", "ctrl_hp_diff",
              "Homopolymer \u0394", paste0("NMD-escape ", MINUS, " CDS homopolymer fraction"),
              paste0("Homopolymer \u0394 (NMD-Escape ", MINUS, " CDS)"),
              "Positive = NMD-escape homopolymer-richer than full CDS", "F",
              zero_line = TRUE)

n_fs  <- sum(pp$group == "FS"); n_snv <- sum(pp$group == "Stopgain")
fig <- (p1 | p2 | p3) / (p4 | p5 | p6) +
  plot_annotation(
    title = "Repeat & Homopolymer Content \u2014 Paired Analysis",
    subtitle = sprintf(
      "Frameshift vs FS Control (n=%d pairs)  \u00b7  Stopgain vs Stopgain Control (n=%d pairs)  \u00b7  Paired Wilcoxon signed-rank, BH-adjusted within metric",
      n_fs, n_snv),
    theme = theme(plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
                  plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey30")))

save_fig("repeat_homopolymer_paired.pdf", fig, 22, 13)
save_fig("repeat_homopolymer_paired.png", fig, 22, 13, dpi = 200)

# ── 7. CDS-vs-NMD scatter, faceted by group (unchanged) ─────────────────────
make_scatter_df <- function(pp, study_x, ctrl_x, study_y, ctrl_y, metric) {
  fs  <- pp %>% filter(group == "FS"); snv <- pp %>% filter(group == "Stopgain")
  bind_rows(
    tibble(Group = "Frameshift",       x = fs[[study_x]],  y = fs[[study_y]],  metric = metric),
    tibble(Group = "FS Control",       x = fs[[ctrl_x]],   y = fs[[ctrl_y]],   metric = metric),
    tibble(Group = "Stopgain",         x = snv[[study_x]], y = snv[[study_y]], metric = metric),
    tibble(Group = "Stopgain Control", x = snv[[ctrl_x]],  y = snv[[ctrl_y]],  metric = metric)
  ) %>% mutate(Group = factor(Group, levels = grp_levels))
}
scat_rep <- make_scatter_df(pp, "study_repeat_fraction","ctrl_repeat_fraction",
                                "study_nmdesc_repeat_fraction","ctrl_nmdesc_repeat_fraction","Repeat")
scat_hp  <- make_scatter_df(pp, "study_homopolymer_fraction","ctrl_homopolymer_fraction",
                                "study_nmdesc_homopolymer_fraction","ctrl_nmdesc_homopolymer_fraction","Homopolymer")
mk_scatter <- function(d, label, title_str) {
  ggplot(d, aes(x = x, y = y, colour = Group, fill = Group, shape = Group)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "dashed", linewidth = 0.6) +
    geom_point(size = 2.4, stroke = 0.9, alpha = 0.78) +
    facet_wrap(~ Group, nrow = 1) +
    scale_colour_manual(values = grp_cols) +
    scale_fill_manual(values = c("Frameshift"="#2166AC","FS Control"=NA,"Stopgain"="#B2182B","Stopgain Control"=NA)) +
    scale_shape_manual(values = c("Frameshift"=21,"FS Control"=21,"Stopgain"=24,"Stopgain Control"=24)) +
    scale_x_continuous(labels = label_percent(accuracy = 1)) +
    scale_y_continuous(labels = label_percent(accuracy = 1)) +
    labs(title = title_str,
         subtitle = "Dashed line = NMD-escape region matches full-CDS composition (y = x)",
         x = sprintf("%s fraction in CDS", label),
         y = sprintf("%s fraction in NMD-escape region", label)) +
    theme_big +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "grey95", linewidth = 0.5),
          strip.text = element_text(face = "bold", size = 13))
}
scat_fig <- mk_scatter(scat_rep, "Repeat", "CDS vs NMD-Escape: Repeat Fraction") /
            mk_scatter(scat_hp,  "Homopolymer", "CDS vs NMD-Escape: Homopolymer Fraction") +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "Within-Gene Repeat & Homopolymer Composition (NMD-Escape vs Full CDS)",
    theme = theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5)))
save_fig("repeat_homopolymer_scatter.pdf", scat_fig, 20, 11)
save_fig("repeat_homopolymer_scatter.png", scat_fig, 20, 11, dpi = 200)

cat("\n=== Done ===\n")
cat("  repeat_homopolymer_paired.pdf / .png   <- 6-panel main figure\n")
cat("  repeat_homopolymer_scatter.pdf / .png  <- CDS-vs-NMD scatter\n")
cat("  repeat_homopolymer_stats.csv           <- paired Wilcoxon stats\n")
