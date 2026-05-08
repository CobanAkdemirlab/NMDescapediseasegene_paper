snv_idr2 = merge(snv_idr0603,snv_merged_uni_loc,by.x = 'Variant_Key',by.y = 'variant_key')
#compare variant loc with wild_IDR, create flag. flag=1 if location within wild_Disordered_Domain_Boundaries

# Function to check if location is within any boundary
is_in_idr <- function(loc, boundary_list) {
  if (is.null(boundary_list)) return(0)
  any(apply(boundary_list, 1, function(b) loc >= b[1] && loc <= b[2])) * 1
}

# Apply the function row-wise
snv_idr2$within_wild_IDR <- mapply(is_in_idr, snv_idr2$location, snv_idr2$wild_boundaries_list)


snv_control_idr2 = merge(snv_control_idr0603,snv_control_merged_uni_loc,by.x = 'Variant_Key',by.y = 'variant_key')
snv_control_idr2$within_wild_IDR <- mapply(is_in_idr, snv_control_idr2$location, snv_control_idr2$wild_boundaries_list)


plus1_control_merged_uni_loc <- read_csv("~/Downloads/uni_loc/plus1_control_merged_uni_loc.csv")
plus1_control_idr2 = merge(plus1_control_idr0603,plus1_control_merged_uni_loc,by.x = 'Variant_Key',by.y = 'variant_key')
plus1_control_idr2$within_wild_IDR <- mapply(is_in_idr, plus1_control_idr2$location, plus1_control_idr2$wild_boundaries_list)


plus1_merged_uni_loc <- read_csv("~/Downloads/uni_loc/plus1_merged_uni_loc.csv")
plus1_idr2 = merge(plus1_idr0603,plus1_merged_uni_loc,by.x = 'Variant_Key',by.y = 'variant_key')
plus1_idr2$within_wild_IDR <- mapply(is_in_idr, plus1_idr2$location, plus1_idr2$wild_boundaries_list)


minus1_control_merged_uni_loc <- read_csv("~/Downloads/uni_loc/minus1_control_merged_uni_loc.csv")
minus1_control_idr2 = merge(minus1_control_idr0603,minus1_control_merged_uni_loc,by.x = 'Variant_Key',by.y = 'variant_key')
minus1_control_idr2$within_wild_IDR <- mapply(is_in_idr, minus1_control_idr2$location, minus1_control_idr2$wild_boundaries_list)


minus1_merged_uni_loc <- read_csv("~/Downloads/uni_loc/minus1_merged_uni_loc.csv")
minus1_idr2 = merge(minus1_idr0603,minus1_merged_uni_loc,by.x = 'Variant_Key',by.y = 'variant_key')
minus1_idr2$within_wild_IDR <- mapply(is_in_idr, minus1_idr2$location, minus1_idr2$wild_boundaries_list)


table(minus1_idr2$within_wild_IDR)
table(minus1_control_idr2$within_wild_IDR)
table(plus1_idr2$within_wild_IDR)
table(plus1_control_idr2$within_wild_IDR)
table(snv_idr2$within_wild_IDR)
table(snv_control_idr2$within_wild_IDR)

library(ggplot2)
library(dplyr)

# Create a summary data frame manually
idr_summary <- data.frame(
  group = c("minus1", "minus1_control", "plus1", "plus1_control", "snv", "snv_control"),
  in_idr = c(219, 1596, 173, 1277, 72, 392),
  total = c(293+219, 2621+1596, 362+173, 2367+1277, 666+72, 3363+392)
)

# Add percentage column
idr_summary <- idr_summary %>%
  mutate(percent_in_idr = 100 * in_idr / total)

library(ggplot2)
library(dplyr)

# Create the summary data frame
idr_summary <- data.frame(
  group = c("minus1", "minus1_Control", "plus1", "plus1_Control", "SNV", "SNV_Control"),
  in_idr = c(219, 1596, 173, 1277, 72, 392),
  total = c(293+219, 2621+1596, 362+173, 2367+1277, 666+72, 3363+392)
)

# Calculate percent in IDR
idr_summary <- idr_summary %>%
  mutate(percent_in_idr = 100 * in_idr / total)

# Define group colors
group_colors <- c(
  "minus1" = "#1f77b4",
  "minus1_Control" = "#aec7e8",
  "plus1" = "#ff7f0e",
  "plus1_Control" = "#ffbb78",
  "SNV" = "#2ca02c",
  "SNV_Control" = "#98df8a"
)

# Plot
ggplot(idr_summary, aes(x = group, y = percent_in_idr, fill = group)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Percentage of Variants within Wild-Type IDR",
       x = "Group",
       y = "Percent within IDR (%)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = group_colors)

