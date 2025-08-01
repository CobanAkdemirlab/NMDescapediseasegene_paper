library(ggplot2)
library(dplyr)

snv_idr0603 <- readRDS("snv_idr0603.rds")
snv_control_idr0603 <- readRDS("snv_control_idr0603.rds")
plus1_idr0603 <- readRDS("plus1_idr0603.rds")
plus1_control_idr0603 <- readRDS("plus1_control_idr0603.rds")
minus1_idr0603 <- readRDS("minus1_idr0603.rds")
minus1_control_idr0603 <- readRDS("minus1_control_idr0603.rds")

# Step 1: Add group and Is_Control indicator
snv_idr0603$group <- "SNV"; snv_idr0603$Is_Control <- FALSE
snv_control_idr0603$group <- "SNV"; snv_control_idr0603$Is_Control <- TRUE
plus1_idr0603$group <- "plus1"; plus1_idr0603$Is_Control <- FALSE
plus1_control_idr0603$group <- "plus1"; plus1_control_idr0603$Is_Control <- TRUE
minus1_idr0603$group <- "minus1"; minus1_idr0603$Is_Control <- FALSE
minus1_control_idr0603$group <- "minus1"; minus1_control_idr0603$Is_Control <- TRUE

# Step 2: Combine and add bar_group + control label
combined_df <- bind_rows(
  snv_idr0603,
  snv_control_idr0603,
  plus1_idr0603,
  plus1_control_idr0603,
  minus1_idr0603,
  minus1_control_idr0603
) %>%
  mutate(
    bar_group = ifelse(Is_Control, paste0(group, "_Control"), group),
    Is_Control_label = ifelse(Is_Control, "Benign", "P/LP")
  )
combined_df$group <- recode(combined_df$group,
                            "minus1" = "Minus1",
                            "plus1" = "Plus1",
                            "SNV" = "Nonsense")

combined_df$group <- factor(combined_df$group,
                            levels = c("Minus1", "Plus1", "Nonsense"))

#remove rows where IDR_length_difference is NA or 0
combined_df <- combined_df %>%
  filter(!is.na(IDR_length_difference) & IDR_length_difference != 0)



# Step 3: Compute p-values per group
get_pval <- function(subdf) {
  tryCatch({
    wilcox.test(IDR_length_difference ~ Is_Control, data = subdf)$p.value
  }, error = function(e) NA)
}

pval_df <- combined_df %>%
  group_by(group) %>%
  summarise(pval = get_pval(cur_data()), .groups = "drop") %>%
  mutate(label = paste0("p = ", format(pval, digits = 2, scientific = TRUE)))

# Step 4: Define your 6-color group scheme
group_colors <- c(
  "minus1" = "#1f77b4",
  "minus1_Control" = "#aec7e8",
  "plus1" = "#ff7f0e",
  "plus1_Control" = "#ffbb78",
  "SNV" = "#2ca02c",
  "SNV_Control" = "#98df8a"
)

# Step 1: Create transformed y for plotting
combined_df_trans <- combined_df %>%
  filter(!is.na(IDR_length_difference) & IDR_length_difference != 0) %>%
  mutate(nmdesc_length = NMDesc.end - NMDesc.start) %>%
  mutate(IDR_log_trans = log10(abs(IDR_length_difference))) %>%
  mutate(IDR_log_norm_trans = log10(abs(IDR_length_difference/nmdesc_length)))

#recalculate p for IDR_log_norm_trans
get_pval_norm <- function(subdf) {
  tryCatch({
    wilcox.test(IDR_log_norm_trans ~ Is_Control, data = subdf)$p.value
  }, error = function(e) NA)
}

pval_df <- combined_df_trans %>%
  group_by(group) %>%
  summarise(pval = get_pval_norm(cur_data()), .groups = "drop") %>%
  mutate(label = paste0("p = ", format(pval, digits = 2, scientific = TRUE)))


  

# Step 2: Plot using the transformed column
ggplot(combined_df_trans, aes(x = Is_Control_label, y = IDR_log_trans, fill = bar_group)) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  facet_wrap(~ group, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  geom_text(data = pval_df, aes(x = 1.5, y = Inf, label = label),
            vjust = 1.5, inherit.aes = FALSE, size = 4.5) +
  labs(
    x = NULL,
    y = "Signed log10(IDR Length Difference)",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "none"
  )
