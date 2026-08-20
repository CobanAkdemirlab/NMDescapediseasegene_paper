library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Match plus1 and plus2 by transcript
pair_df <- PTC_info %>%
  filter(
    can_region_start != 1,
    type %in% c("plus1", "plus2")
  ) %>%
  select(transcript, type, can_region_start) %>%
  distinct() %>%
  group_by(transcript) %>%
  filter(n_distinct(type) == 2) %>%
  ungroup() %>%
  pivot_wider(
    names_from = type,
    values_from = can_region_start
  ) %>%
  filter(!is.na(plus1), !is.na(plus2)) %>%
  mutate(
    pair_id = row_number(),
    diff = plus1 - plus2,
    abs_diff = abs(diff)
  )

# 2. Difference magnitude statistics
cor_val <- cor(pair_df$plus1, pair_df$plus2)

median_abs_diff <- median(pair_df$abs_diff)
mean_abs_diff <- mean(pair_df$abs_diff)

pct_same <- mean(pair_df$diff == 0) * 100
pct_5aa <- mean(pair_df$abs_diff <= 5) * 100
pct_30aa <- mean(pair_df$abs_diff <= 30) * 100

summary_stats <- pair_df %>%
  summarise(
    n_transcripts = n(),
    correlation = cor(plus1, plus2),
    median_abs_diff = median(abs_diff),
    mean_abs_diff = mean(abs_diff),
    pct_identical = mean(diff == 0) * 100,
    pct_within_5aa = mean(abs_diff <= 5) * 100,
    pct_within_50aa = mean(abs_diff <= 50) * 100
  )

print(summary_stats)

# 3. Long format for plotting
plot_df <- pair_df %>%
  pivot_longer(
    cols = c(plus1, plus2),
    names_to = "Frame",
    values_to = "can_region_start"
  ) %>%
  mutate(
    Frame = factor(
      Frame,
      levels = c("plus1", "plus2"),
      labels = c("Plus1", "Minus1")
    )
  )

# 4. Y-axis limit
y_upper <- quantile(plot_df$can_region_start, 0.95, na.rm = TRUE)

# 5. Paired boxplot
p <- ggplot(plot_df, aes(x = Frame, y = can_region_start)) +
  
  geom_line(
    aes(group = pair_id),
    color = "grey45",
    alpha = 0.12,
    linewidth = 0.35
  ) +
  
  geom_boxplot(
    aes(fill = Frame),
    width = 0.35,
    alpha = 0.85,
    outlier.shape = NA,
    linewidth = 1.1
  ) +
  
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 3.5,
    fill = "white"
  ) +
  
  scale_fill_manual(
    values = c(
      "Plus1" = "#2166AC",
      "Minus1" = "#67A9CF"
    )
  ) +
  
  coord_cartesian(
    ylim = c(0, y_upper * 1.12)
  ) +
  
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = y_upper * 1.03,
    yend = y_upper * 1.03,
    linewidth = 0.8
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = y_upper * 1.08,
    label = paste0(
      "Median |Δ| = ",
      round(median_abs_diff, 1),
      " aa"
    ),
    size = 5.5,
    fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = 1.5,
    y = y_upper * 0.97,
    label = paste0(
      "n = ", nrow(pair_df),
      " matched transcripts\n"
    ),
    size = 4.5
  ) +
  
  labs(
    x = NULL,
    y = "Canonical NMDesc region start",
    title = "Canonical NMDesc Region is Similar for +1 and −1 Frames"
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

print(p)