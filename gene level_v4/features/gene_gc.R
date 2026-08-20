#This R script is to get gene-level GC content
matching_log_v2 <- read_csv("~/Downloads/matching_log_v2.csv")
str(matching_log_v2)
get_gc_content <- function(sequence) {
  bases <- strsplit(toupper(sequence), "")[[1]]
  gc <- sum(bases %in% c("G", "C"))
  return(round(gc / length(bases) * 100, 2))
}

gene_all = gene_all_AD_in
##ovarall gene
gene_all$gc_content = sapply(gene_all$coding, get_gc_content)
##NMDesc region
gene_all$nmdesc_gc_content = sapply(gene_all$nmdesc_cds, get_gc_content)

# boxplot
ggplot(gene_all, aes(x = group, y = gc_content, fill = group)) +
  geom_boxplot(width = 0.7, outlier.shape = 16, outlier.size = 1.5) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format"
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    x = "",
    y = "GC content (%)",
    title = "GC content in NMDesc regions"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )


# ── Libraries ────────────────────────────────────────────────────────────────
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── Read matched transcript pairs ────────────────────────────────────────────
matching_log_v2 <- read_csv("~/Downloads/matching_log_v2.csv")

# ── Function: calculate GC percentage safely ─────────────────────────────────
get_gc_content <- function(sequence) {
  
  if (is.na(sequence) || !nzchar(sequence)) {
    return(NA_real_)
  }
  
  bases <- strsplit(toupper(as.character(sequence)), "")[[1]]
  
  # Keep only standard DNA bases
  valid_bases <- bases[bases %in% c("A", "T", "G", "C")]
  
  if (length(valid_bases) == 0) {
    return(NA_real_)
  }
  
  gc <- sum(valid_bases %in% c("G", "C"))
  
  round(gc / length(valid_bases) * 100, 2)
}

# ── Calculate GC content for all transcripts ─────────────────────────────────
gene_gc <- gene_all_AD_in %>%
  mutate(
    group_key = tolower(as.character(group)),
    
    # GC content for the full coding sequence
    gc_content = vapply(coding, get_gc_content, numeric(1)),
    
    # GC content for the NMDesc retained CDS region
    nmdesc_gc_content = vapply(nmdesc_cds, get_gc_content, numeric(1))
  )

# Choose the GC variable to analyze:
# "nmdesc_gc_content" = GC content in NMDesc region
# "gc_content"        = GC content across entire CDS
gc_variable <- "nmdesc_gc_content"

# ── Keep one GC value per transcript and group ───────────────────────────────
transcript_gc <- gene_gc %>%
  filter(!is.na(ensembl_transcript_id)) %>%
  transmute(
    ensembl_transcript_id,
    group_key,
    gc_value = .data[[gc_variable]]
  ) %>%
  group_by(ensembl_transcript_id, group_key) %>%
  summarise(
    gc_value = if (any(!is.na(gc_value))) {
      first(gc_value[!is.na(gc_value)])
    } else {
      NA_real_
    },
    .groups = "drop"
  )

# ── Define disease/control groups according to matching file ─────────────────
match_pairs <- matching_log_v2 %>%
  mutate(
    match_class = case_when(
      tolower(group) %in% c("stopgain", "nonsense", "snv") ~ "Nonsense",
      tolower(group) %in% c("frameshift", "fs") ~ "Frameshift",
      TRUE ~ as.character(group)
    ),
    
    study_group = case_when(
      match_class == "Nonsense" ~ "snv",
      match_class == "Frameshift" ~ "fs",
      TRUE ~ NA_character_
    ),
    
    control_group = case_when(
      match_class == "Nonsense" ~ "snv_control",
      match_class == "Frameshift" ~ "fs_control",
      TRUE ~ NA_character_
    ),
    
    pair_id = row_number()
  ) %>%
  filter(!is.na(study_group), !is.na(control_group))

# ── Attach study and matched-control GC values ───────────────────────────────
study_gc <- transcript_gc %>%
  rename(
    study_transcript = ensembl_transcript_id,
    study_group = group_key,
    study_gc = gc_value
  )

control_gc <- transcript_gc %>%
  rename(
    ctrl_transcript = ensembl_transcript_id,
    control_group = group_key,
    control_gc = gc_value
  )

matched_gc <- match_pairs %>%
  left_join(
    study_gc,
    by = c("study_transcript", "study_group")
  ) %>%
  left_join(
    control_gc,
    by = c("ctrl_transcript", "control_group")
  ) %>%
  filter(
    !is.na(study_gc),
    !is.na(control_gc)
  )

# Number of usable matched pairs
matched_gc %>%
  count(match_class, name = "n_matched_pairs") %>%
  print()

# ── Paired Wilcoxon tests: raw p-values, no BH adjustment ────────────────────
paired_results <- matched_gc %>%
  group_by(match_class) %>%
  group_modify(~ {
    
    test_result <- wilcox.test(
      .x$study_gc,
      .x$control_gc,
      paired = TRUE,
      exact = FALSE
    )
    
    tibble(
      n_pairs = nrow(.x),
      median_disease = median(.x$study_gc, na.rm = TRUE),
      median_control = median(.x$control_gc, na.rm = TRUE),
      median_difference = median(.x$study_gc - .x$control_gc, na.rm = TRUE),
      p_value = test_result$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_label = paste0(
      "Paired Wilcoxon\np = ",
      format.pval(p_value, digits = 3, eps = 0.001)
    )
  )

print(paired_results)

# ── Convert matched data to long format for plotting ─────────────────────────
matched_gc_long <- matched_gc %>%
  select(pair_id, match_class, study_gc, control_gc) %>%
  pivot_longer(
    cols = c(study_gc, control_gc),
    names_to = "status",
    values_to = "gc_value"
  ) %>%
  mutate(
    status = recode(
      status,
      study_gc = "Disease",
      control_gc = "Matched control"
    ),
    status = factor(
      status,
      levels = c("Disease", "Matched control")
    ),
    match_class = factor(
      match_class,
      levels = c("Frameshift", "Nonsense")
    )
  )

# ── Set p-value annotation positions ─────────────────────────────────────────
y_positions <- matched_gc_long %>%
  group_by(match_class) %>%
  summarise(
    y_position = max(gc_value, na.rm = TRUE) +
      0.10 * diff(range(gc_value, na.rm = TRUE)),
    .groups = "drop"
  )

paired_results_plot <- paired_results %>%
  left_join(y_positions, by = "match_class") %>%
  mutate(
    x_position = 1.5
  )

# ── Matched GC-content figure ────────────────────────────────────────────────
p_matched_gc <- ggplot(
  matched_gc_long,
  aes(x = status, y = gc_value)
) +
  # Each line represents one matched disease-control transcript pair
  geom_line(
    aes(group = pair_id),
    linewidth = 0.35,
    alpha = 0.35,
    color = "grey45"
  ) +
  # No individual dots
  geom_boxplot(
    aes(fill = status),
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_text(
    data = paired_results_plot,
    aes(
      x = x_position,
      y = y_position,
      label = p_label
    ),
    inherit.aes = FALSE,
    size = 3.5,
    vjust = 0
  ) +
  facet_wrap(
    ~ match_class,
    scales = "free_y"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.05, 0.22))
  ) +
  labs(
    x = NULL,
    y = "GC content (%)",
    title = "Matched comparison of GC content in NMDesc regions"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold")
  )

print(p_matched_gc)

write.csv(matched_gc, file = "matched_gc.csv")