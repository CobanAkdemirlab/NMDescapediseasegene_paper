library(tidyverse)

# ── Structural data: NMDPos VAR - WT analysis only ─────────────
# assume uniprot_id and gene id are 1-1

#1. create NMDPos VAR - WT values

props = c(
  "plddt_NMDPos_mean",
  "rel_asa_NMDPos_mean",
  "abs_sasa_dssp_NMDPos_mean",
  "sec_struc_NMDPos_helix_percent",
  "sec_struc_NMDPos_beta_percent",
  "sec_struc_NMDPos_coil_percent",
  "SASA_total_NMDPos",
  "pae_NMDPos_mean"
)

for (p in props) {
  VAR_WT_structural_results[[paste0("d_", p)]] =
    VAR_WT_structural_results[[paste0("VAR_", p)]] -
    VAR_WT_structural_results[[paste0("WT_", p)]]
}

d_vars = paste0("d_", props)


#2. paired Wilcoxon function

get_paired_p = function(data, d_var, group1, group2) {
  
  tmp = data %>%
    filter(source_folder %in% c(group1, group2)) %>%
    select(uniprot_id, source_folder, value = all_of(d_var)) %>%
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


#3. repeat random sampling 1000 times

set.seed(123)

all_p_results_structural_NMD = map_dfr(1:1000, function(i) {
  
  sampled_data = VAR_WT_structural_results %>%
    group_by(uniprot_id, source_folder) %>%
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

all_p_results_structural_NMD


#4. output summary

structural_summary_NMD = all_p_results_structural_NMD %>%
  group_by(d_var) %>%
  summarise(
    median_fs_p = median(fs_p, na.rm = TRUE),
    median_snv_p = median(snv_p, na.rm = TRUE),
    
    median_fs_FDR = median(fs_FDR, na.rm = TRUE),
    median_snv_FDR = median(snv_FDR, na.rm = TRUE),
    
    prop_fs_FDR_lt_0.05 = mean(fs_FDR < 0.05, na.rm = TRUE),
    prop_snv_FDR_lt_0.05 = mean(snv_FDR < 0.05, na.rm = TRUE),
    
    .groups = "drop"
  )

structural_summary_NMD


#5. save output

write_csv(
  all_p_results_structural_NMD,
  "all_p_results_structural_NMD_1000_repeats.csv"
)

write_csv(
  structural_summary_NMD,
  "structural_summary_NMD_1000_repeats.csv"
)