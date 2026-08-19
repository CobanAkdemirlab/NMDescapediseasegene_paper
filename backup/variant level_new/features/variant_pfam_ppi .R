# ============================================================
# analyze_pfam_ppi()
#
# Analyzes PFAM domain overlap and PPI interface overlap for
# variants, produces summary plots, and returns annotated data.
#
# Actual variants_all2 column mapping:
#   Variant_Key           → unique variant identifier
#   ensembl_transcript_id → transcript ID
#   ptc_pos               → PTC position in CDS (bp)
#   group                 → sample group (e.g. "fs_disease")
#   (uniprot is joined in from pfam_fin via ensembl_transcript_id)
#
# Required inputs:
#   variants_all2  — data.frame (see str() above for columns)
#   human_1_       — data.frame with columns:
#                     uniprot1, uniprot2,
#                     interface_residues1, interface_residues2
#   pfam_fin       — data.frame with columns:
#                     ensembl_transcript_id, pfam, pfam_start,
#                     pfam_end, hgnc_symbol, uniprot
#   ensembl        — a biomaRt Mart object (from useMart())
#
# Optional:
#   out_dir        — directory to write CSVs and plots (default: ".")
#   group_levels   — preferred factor order for the `group` column;
#                    levels absent from the data are silently dropped
#   group_colors   — named color vector keyed to group values
#
# Returns a named list:
#   $variants     — annotated data frame
#   $ppi_summary  — per-group PPI flag counts/proportions
#   $pfam_summary — per-group PFAM flag counts/proportions
#   $plots        — named list of ggplot objects
# ============================================================

library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(biomaRt)

variant_pfam_ppi <- function(
    variants_all2,
    human_1_,
    pfam_fin,
    ensembl,
    out_dir      = ".",
    group_levels = c("fs_disease", "fs_control", "snv_disease", "snv_control"),
    group_colors = c(
      "fs_disease"  = "#1f77b4",
      "fs_control"  = "#aec7e8",
      "snv_disease" = "#2ca02c",
      "snv_control" = "#98df8a"
    )
) {
  
  # ----------------------------------------------------------
  # 0. Internal helpers
  # ----------------------------------------------------------
  
  # "[12,34,56]" -> numeric vector (used for bp conversion)
  .convert_to_c <- function(x) {
    if (is.na(x) || x == "") return(numeric(0))
    x <- gsub("\\[|\\]|\\s", "", x)
    if (x == "") return(numeric(0))
    vals <- suppressWarnings(as.numeric(unlist(strsplit(x, ","))))
    vals[!is.na(vals)]
  }
  
  # "[12,34,56]" -> integer vector
  .convert_to_ints <- function(x) {
    if (is.null(x) || is.na(x) || x == "") return(integer(0))
    s <- trimws(gsub("\\[|\\]", "", x))
    if (s == "") return(integer(0))
    as.integer(trimws(unlist(strsplit(s, ","))))
  }
  
  # CDS base position -> amino acid position
  .cds2aa <- function(base_idx) ((as.integer(base_idx) - 1L) %/% 3L) + 1L
  
  # Fisher's exact test -> p-value
  .fisher_p <- function(n_event_a, n_total_a, n_event_b, n_total_b) {
    mat <- matrix(
      c(n_event_a,   n_total_a - n_event_a,
        n_event_b,   n_total_b - n_event_b),
      nrow = 2, byrow = TRUE
    )
    fisher.test(mat)$p.value
  }
  
  # Build p-value bracket data frame for stat_pvalue_manual()
  .pval_df <- function(summary_df, prop_col, event_col, total_col,
                       pairs, y_mult = c(1.08, 1.20)) {
    ymax <- max(summary_df[[prop_col]], na.rm = TRUE)
    rows <- lapply(seq_along(pairs), function(i) {
      g1 <- pairs[[i]][1]; g2 <- pairs[[i]][2]
      if (!all(c(g1, g2) %in% summary_df$group)) return(NULL)
      p <- .fisher_p(
        summary_df[[event_col]][summary_df$group == g1],
        summary_df[[total_col]][summary_df$group == g1],
        summary_df[[event_col]][summary_df$group == g2],
        summary_df[[total_col]][summary_df$group == g2]
      )
      mult <- if (i <= length(y_mult)) y_mult[i] else y_mult[length(y_mult)] + 0.12 * (i - length(y_mult))
      data.frame(group1 = g1, group2 = g2,
                 y.position = ymax * mult,
                 label      = paste0("p = ", signif(p, 3)),
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, Filter(Negate(is.null), rows))
  }
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ----------------------------------------------------------
  # 1. Validate required columns
  # ----------------------------------------------------------
  required_cols <- c("Variant_Key", "ensembl_transcript_id", "ptc_pos", "group")
  missing_cols  <- setdiff(required_cols, names(variants_all2))
  if (length(missing_cols) > 0)
    stop("variants_all2 is missing columns: ", paste(missing_cols, collapse = ", "))
  
  if (!"uniprot" %in% names(pfam_fin))
    stop("pfam_fin must contain a 'uniprot' column to map transcripts to Uniprot IDs.")
  
  # ----------------------------------------------------------
  # 2. Working copy + join uniprot from pfam_fin
  # ----------------------------------------------------------
  uniprot_map <- pfam_fin %>%
    filter(!is.na(uniprot), uniprot != "") %>%
    distinct(ensembl_transcript_id, uniprot) %>%
    group_by(ensembl_transcript_id) %>%
    slice(1) %>%
    ungroup()
  
  dat <- variants_all2 %>%
    left_join(uniprot_map, by = "ensembl_transcript_id") %>%
    mutate(
      variant_ppi_overlap                      = 0L,
      variant_ppi_nearest_interface_bp         = NA_real_,
      variant_ppi_dist_to_nearest_interface_bp = NA_real_
    )
  
  # ----------------------------------------------------------
  # 3. PPI: flag variants with a downstream interface residue
  #    interface residues (aa) x 3 -> bp; compare with ptc_pos (bp)
  # ----------------------------------------------------------
  for (i in seq_len(nrow(dat))) {
    
    uid     <- dat$uniprot[[i]]
    cds_ptc <- dat$ptc_pos[[i]]
    
    if (length(uid) != 1 || length(cds_ptc) != 1) next
    if (is.na(uid)  || uid == "")                  next
    if (is.na(cds_ptc))                            next
    
    re1 <- unlist(lapply(
      human_1_$interface_residues1[human_1_$uniprot1 == uid],
      .convert_to_c
    )) * 3
    
    re2 <- unlist(lapply(
      human_1_$interface_residues2[human_1_$uniprot2 == uid],
      .convert_to_c
    )) * 3
    
    all_iface_bp <- na.omit(c(re1, re2))
    if (length(all_iface_bp) == 0) next
    
    downstream <- all_iface_bp[all_iface_bp >= cds_ptc]
    
    if (length(downstream) > 0) {
      dat$variant_ppi_overlap[i]                      <- 1L
      dat$variant_ppi_nearest_interface_bp[i]         <- min(downstream)
      dat$variant_ppi_dist_to_nearest_interface_bp[i] <- min(downstream) - cds_ptc
    }
  }
  
  # ----------------------------------------------------------
  # 4. Factor-ise group column (keep only levels present in data)
  # ----------------------------------------------------------
  present_levels <- unique(as.character(dat$group))
  ordered_levels <- c(
    intersect(group_levels, present_levels),
    setdiff(present_levels, group_levels)
  )
  dat$group <- factor(dat$group, levels = ordered_levels)
  
  # Restrict group_colors to levels actually present
  group_colors <- group_colors[names(group_colors) %in% present_levels]
  
  # Pair adjacent levels for comparisons: (1 vs 2), (3 vs 4), etc.
  comparisons <- lapply(
    seq(1, length(ordered_levels) - 1, by = 2),
    function(k) ordered_levels[k:(k + 1)]
  )
  
  # ----------------------------------------------------------
  # 5. PPI summary & plots
  # ----------------------------------------------------------
  ppi_summary <- dat %>%
    group_by(group) %>%
    summarise(
      n       = n(),
      n_match = sum(variant_ppi_overlap == 1L, na.rm = TRUE),
      prop    = n_match / n,
      .groups = "drop"
    )
  
  pv_ppi   <- .pval_df(ppi_summary, "prop", "n_match", "n", comparisons)
  ymax_ppi <- max(ppi_summary$prop, na.rm = TRUE)
  
  p_ppi_prop <- ggplot(ppi_summary, aes(x = group, y = prop, fill = group)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = group_colors, guide = "none") +
    geom_text(aes(label = paste0(n_match, "/", n)), vjust = -0.3, size = 4) +
    stat_pvalue_manual(pv_ppi, label = "label", xmin = "group1", xmax = "group2",
                       y.position = "y.position", tip.length = 0.01) +
    expand_limits(y = ymax_ppi * 1.35) +
    labs(title = "Proportion of variants with downstream PPI interface residue",
         x = NULL, y = "Proportion") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  p_ppi_dist <- ggplot(
    dat %>% filter(!is.na(variant_ppi_dist_to_nearest_interface_bp)),
    aes(x = group, y = variant_ppi_dist_to_nearest_interface_bp, fill = group)
  ) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
    scale_fill_manual(values = group_colors, guide = "none") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       label = "p.format") +
    labs(title = "Distance from PTC to nearest downstream PPI interface residue",
         x = NULL, y = "Distance (bp)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  p_ppi_log <- ggplot(
    dat %>% filter(!is.na(variant_ppi_dist_to_nearest_interface_bp),
                   variant_ppi_dist_to_nearest_interface_bp > 0),
    aes(x = group, y = variant_ppi_dist_to_nearest_interface_bp, fill = group)
  ) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
    scale_y_log10() +
    scale_fill_manual(values = group_colors, guide = "none") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       label = "p.format") +
    labs(title = "Log10 distance from PTC to nearest downstream PPI interface residue",
         x = NULL, y = "Distance (bp, log10)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  # ----------------------------------------------------------
  # 6. PFAM: join max pfam_end per transcript via biomaRt
  # ----------------------------------------------------------
  pfam_bm <- getBM(
    attributes = c("ensembl_transcript_id", "pfam_end"),
    filters    = "ensembl_transcript_id",
    values     = unique(dat$ensembl_transcript_id),
    mart       = ensembl
  )
  
  max_pfam_end_df <- pfam_bm %>%
    filter(!is.na(pfam_end)) %>%
    group_by(ensembl_transcript_id) %>%
    summarise(max_pfam_end = max(pfam_end), .groups = "drop")
  
  dat <- dat %>%
    mutate(ptc_aa = .cds2aa(ptc_pos)) %>%
    left_join(max_pfam_end_df, by = "ensembl_transcript_id") %>%
    mutate(
      dist_ptc_to_max_pfam_end_aa = ptc_aa - max_pfam_end,
      ptc_after_max_pfam_end      = as.integer(dist_ptc_to_max_pfam_end_aa >  0),
      ptc_before_max_pfam_end     = as.integer(dist_ptc_to_max_pfam_end_aa <  0)
    )
  
  # Join per-domain pfam info (expands rows if multiple domains per transcript)
  dat <- dat %>%
    left_join(
      pfam_fin %>%
        select(ensembl_transcript_id, pfam, pfam_start, pfam_end, hgnc_symbol),
      by = "ensembl_transcript_id"
    ) %>%
    mutate(
      in_pfam = !is.na(pfam_start) & !is.na(pfam_end) &
        !is.na(ptc_aa)     & ptc_aa <= pfam_end
    )
  
  # ----------------------------------------------------------
  # 7. PFAM summary & plots
  # ----------------------------------------------------------
  pfam_summary <- dat %>%
    group_by(group) %>%
    summarise(
      n           = n(),
      n_after     = sum(ptc_after_max_pfam_end  == 1L, na.rm = TRUE),
      n_before    = sum(ptc_before_max_pfam_end == 1L, na.rm = TRUE),
      prop_after  = n_after  / n,
      prop_before = n_before / n,
      .groups = "drop"
    )
  
  pv_pfam   <- .pval_df(pfam_summary, "prop_after", "n_after", "n", comparisons)
  ymax_pfam <- max(pfam_summary$prop_after, na.rm = TRUE)
  
  p_pfam_dist <- ggplot(
    dat %>% filter(!is.na(dist_ptc_to_max_pfam_end_aa)),
    aes(x = group, y = dist_ptc_to_max_pfam_end_aa, fill = group)
  ) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
    scale_fill_manual(values = group_colors, guide = "none") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       label = "p.format") +
    labs(title = "Distance from PTC to the largest PFAM end",
         x = NULL, y = "PTC aa position - largest PFAM end (aa)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  p_pfam_abs <- ggplot(
    dat %>%
      filter(!is.na(dist_ptc_to_max_pfam_end_aa)) %>%
      mutate(abs_dist = abs(dist_ptc_to_max_pfam_end_aa)),
    aes(x = group, y = abs_dist, fill = group)
  ) +
    geom_violin(trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
    scale_fill_manual(values = group_colors, guide = "none") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                       label = "p.format") +
    labs(title = "Absolute distance from PTC to the largest PFAM end",
         x = NULL, y = "Absolute distance (aa)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  p_pfam_prop <- ggplot(pfam_summary,
                        aes(x = group, y = prop_before, fill = group)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = group_colors, guide = "none") +
    geom_text(aes(label = paste0(n_before, "/", n)), vjust = -0.3, size = 4) +
    stat_pvalue_manual(pv_pfam, label = "label", xmin = "group1", xmax = "group2",
                       y.position = "y.position", tip.length = 0.01) +
    expand_limits(y = ymax_pfam * 1.35) +
    labs(title = "Proportion of variants with PTC influencing PFAM domain",
         x = NULL, y = "Proportion") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  pfam_pct_df <- dat %>%
    count(group, in_pfam) %>%
    group_by(group) %>%
    tidyr::complete(in_pfam = c(FALSE, TRUE), fill = list(n = 0)) %>%
    mutate(prop = n / sum(n),
           pct  = scales::percent(prop, accuracy = 0.1)) %>%
    ungroup()
  
  p_pfam_stack <- ggplot(pfam_pct_df,
                         aes(x = group, y = prop, fill = as.character(in_pfam))) +
    geom_col(color = "grey30") +
    geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "#2ca02c"),
                      labels = c("Outside PFAM", "Influences PFAM"), name = NULL) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(title = "% of variants influencing PFAM domains",
         x = NULL, y = "Fraction of variants") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title  = element_text(hjust = 0.5, face = "bold"))
  
  # ----------------------------------------------------------
  # 8. Save outputs
  # ----------------------------------------------------------
  write.csv(dat,          file.path(out_dir, "variants_annotated.csv"),              row.names = FALSE)
  write.csv(ppi_summary,  file.path(out_dir, "ppi_flag_summary.csv"),                row.names = FALSE)
  write.csv(pfam_summary, file.path(out_dir, "pfam_flag_summary.csv"),               row.names = FALSE)
  write.csv(pfam_pct_df,  file.path(out_dir, "pfam_influenced_variant_counts.csv"),  row.names = FALSE)
  
  plot_list <- list(
    ppi_proportion    = p_ppi_prop,
    ppi_distance      = p_ppi_dist,
    ppi_distance_log  = p_ppi_log,
    pfam_distance     = p_pfam_dist,
    pfam_abs_distance = p_pfam_abs,
    pfam_proportion   = p_pfam_prop,
    pfam_stacked      = p_pfam_stack
  )
  
  for (nm in names(plot_list)) {
    ggsave(file.path(out_dir, paste0(nm, ".pdf")), plot_list[[nm]], width = 9, height = 5)
    ggsave(file.path(out_dir, paste0(nm, ".png")), plot_list[[nm]], width = 9, height = 5, dpi = 300)
  }
  
  invisible(list(
    variants     = dat,
    ppi_summary  = ppi_summary,
    pfam_summary = pfam_summary,
    plots        = plot_list
  ))
}


# ----------------------------------------------------------
# Example usage
# ----------------------------------------------------------
# Check actual group values first:
#   unique(variants_all2$group)
#
# results <- variant_pfam_ppi(
#   variants_all2 = variants_all2,
#   human_1_      = human_1_,
#   pfam_fin      = pfam_fin,
#   ensembl       = ensembl,
#   out_dir       = "plots/pfam_ppi_analysis",
#   group_levels  = c("fs_disease", "fs_control"),  # match your actual group names
#   group_colors  = c("fs_disease" = "#1f77b4", "fs_control" = "#aec7e8")
# )
#
# print(results$plots$ppi_proportion)
# head(results$variants)