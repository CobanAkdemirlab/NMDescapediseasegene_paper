library(readr)
Kat6b_PARSEV2_SCORES <- read_csv("~/Desktop/Kat6b_PARSEV2_SCORES.csv")
parse_score = Kat6b_PARSEV2_SCORES[,1:3]
parse_score[52,3] = 'WT'

library(ggplot2)
library(dplyr)

# --- Prepare data ---
df <- df %>%
  rename(ParseV2 = `Parse V2 scores`,
         Category = Category)

# Standardize category names
df$Category <- tolower(df$Category)
df$Category[df$Category == "snv"] <- "nonsense"

# Capitalize
df$Category <- recode(df$Category,
                      "minus1" = "Minus1",
                      "plus1" = "Plus1",
                      "nonsense" = "Nonsense",
                      "wt" = "WT")

# Factor order
df$Category <- factor(df$Category,
                      levels = c("Minus1", "Plus1", "Nonsense", "WT"))

# Define colors
group_colors <- c("Minus1" = "#1f77b4",
                  "Plus1" = "#ff7f0e",
                  "Nonsense" = "#2ca02c",
                  "WT" = "#7f7f7f")

# --- Get WT value (ensure not NA) ---
wt_value <- df %>%
  filter(Category == "WT") %>%
  pull(ParseV2) %>%
  mean(na.rm = TRUE)

cat("WT value used:", wt_value, "\n")  # sanity check

# --- Plot ---library(ggplot2)
library(dplyr)

# --- Prepare data ---
df <- df %>%
  rename(ParseV2 = `Parse V2 scores`,
         Category = Category) %>%
  filter(!is.na(Category), !is.na(ParseV2))   # remove NA rows

# Standardize category names
df$Category <- tolower(df$Category)
df$Category[df$Category == "snv"] <- "nonsense"

# Capitalize display names
df$Category <- recode(df$Category,
                      "minus1" = "Minus1",
                      "plus1"  = "Plus1",
                      "nonsense" = "Nonsense",
                      "wt" = "WT")

# Set factor order
df$Category <- factor(df$Category,
                      levels = c("Minus1", "Plus1", "Nonsense", "WT"))

# Define colors
group_colors <- c("Minus1" = "#1f77b4",
                  "Plus1"  = "#ff7f0e",
                  "Nonsense" = "#2ca02c",
                  "WT" = "#7f7f7f")

# Manual WT reference line
wt_value <- 667.4

# --- Plot ---
ggplot(df, aes(x = Category, y = ParseV2, fill = Category)) +
  geom_hline(yintercept = wt_value,
             color = group_colors["WT"],
             linetype = "dashed",
             linewidth = 1.3) +
  geom_boxplot(alpha = 0.3, width = 0.45, outlier.shape = NA) +
  geom_jitter(aes(color = Category),
              width = 0.15, size = 3, alpha = 0.9) +
  scale_fill_manual(values = group_colors, na.translate = FALSE) +
  scale_color_manual(values = group_colors, na.translate = FALSE) +
  labs(
    title = "Parse V2 Scores across Variant Categories",
    subtitle = "Dashed line = WT reference (667.4)",
    x = "Variant Category",
    y = "Parse V2 Score"
  ) +
  coord_cartesian(ylim = c(0, 700)) +
  theme_classic(base_size = 16) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
    axis.text.x   = element_text(face = "bold", size=16),
    axis.text.y   = element_text(face = "bold",size=16),
    axis.title.x  = element_text(face = "bold"),
    axis.title.y  = element_text(face = "bold")
  )

write.csv(omim_AD_symbols, file = "omim_AD_symbols.csv")

Tugce_genelist_fortheppaer <- read_excel("~/Desktop/Tugce_genelist_fortheppaer.xlsx")
AD_Tugce = Tugce_genelist_fortheppaer[which(Tugce_genelist_fortheppaer$Gene_name %in% omim_AD_symbols),]
write_csv(AD_Tugce, file = "AD_Tugce.csv")