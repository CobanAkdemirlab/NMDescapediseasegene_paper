library(tidyverse)

# assume uniprot and gene id are 1-1
# therefore use uniprot to match with the PARSE data
# d = var - wt

#1. prepare source_folder and create VAR - WT values

nmd_length_NMDPos_PARSE_v2 = nmd_length_NMDPos_PARSE_v2 %>%
  mutate(
    source_folder = str_extract(
      id,
      "^(fs_disease|fs_control|snv_disease|snv_control)"
    ),
    
    d_NMD_length =
      var_NMD_length_old - wt_NMD_length_old,
    
    d_classifier_distance_P =
      `Σ classifier distance P...4` -
      `Σ classifier distance P...9`,
    
    d_longest_PS_IDR =
      `length of longest PS IDR...5` -
      `length of longest PS IDR...10`
  )

d_vars = c(
  "d_NMD_length",
  "d_classifier_distance_P",
  "d_longest_PS_IDR"
)


#2. paired Wilcoxon function

get_paired_p = function(data, d_var, group1, group2) {
  
  tmp = data %>%
    filter(source_folder %in% c(group1, group2)) %>%
    select(uniprot, source_folder, value = all_of(d_var)) %>%
    filter(!is.na(value)) %>%
    pivot_wider(
      names_from = source_folder,
      values_from = value
    ) %>%
    drop_na(all_of(c(group1, group2)))
  
  p = wilcox.test(
    tmp[[group1]],
    tmp[[group2]],
    paired = TRUE
  )$p.value
  
  return(p)
}


#3. randomly select one variant for each uniprot per group, repeat 1000 times

set.seed(123)

all_p_results_PARSE = map_dfr(1:1000, function(i) {
  
  sampled_data = nmd_length_NMDPos_PARSE_v2 %>%
    group_by(uniprot, source_folder) %>%
    sample_n(1) %>%
    ungroup()
  
  p_results = tibble(
    d_var = d_vars
  ) %>%
    rowwise() %>%
    mutate(
      fs_p = get_paired_p(
        sampled_data,
        d_var,
        "fs_disease",
        "fs_control"
      ),
      
      snv_p = get_paired_p(
        sampled_data,
        d_var,
        "snv_disease",
        "snv_control"
      )
    ) %>%
    ungroup() %>%
    mutate(
      fs_FDR = p.adjust(fs_p, method = "BH"),
      snv_FDR = p.adjust(snv_p, method = "BH"),
      
      fs_significant_FDR_0.05 = fs_FDR <= 0.05,
      snv_significant_FDR_0.05 = snv_FDR <= 0.05,
      
      iteration = i
    )
  
  p_results
})

all_p_results_PARSE


#4. output summary

fdr_summary_PARSE = all_p_results_PARSE %>%
  group_by(d_var) %>%
  summarise(
    median_fs_p = median(fs_p, na.rm = TRUE),
    median_snv_p = median(snv_p, na.rm = TRUE),
    median_fs_FDR = median(fs_FDR, na.rm = TRUE),
    median_snv_FDR = median(snv_FDR, na.rm = TRUE),
    .groups = "drop"
  )

fdr_summary_PARSE