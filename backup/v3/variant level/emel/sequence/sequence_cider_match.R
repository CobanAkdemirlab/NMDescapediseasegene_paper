library(tidyverse)

# assume uniprot id and gene id are 1-1
# therefore use uniprot id to match with the CIDER data
# d = var - wt

#1. create the var-wt value

props = c(
  "length",
  "fcr_7",
  "isoelectric_point",
  "molecular_weight",
  "count_neg",
  "count_pos",
  "count_neut",
  "fraction_negative",
  "fraction_positive",
  "fraction_expanding_7",
  "fraction_disorder_promoting",
  "mean_hydropathy",
  "uversky_hydropathy",
  "kappa",
  "Omega",
  "delta",
  "deltaMax"
)

for (p in props) {
  NMD_region_NMDPos_CIDER[[paste0("d_", p)]] =
    NMD_region_NMDPos_CIDER[[paste0("var_", p)]] -
    NMD_region_NMDPos_CIDER[[paste0("wt_", p)]]
}

d_vars = paste0("d_", props)


#2. randomly select one variant for each uniprot id per group, repeat this process 1000 times

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


set.seed(123)

all_p_results_CIDER = map_dfr(1:1000, function(i) {
  
  sampled_data = NMD_region_NMDPos_CIDER %>%
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

all_p_results_CIDER


#3. write output: d_var, median_fs_p, median_snv_p, median_fs_FDR, median_snv_FDR

fdr_summary_CIDER = all_p_results_CIDER %>%
  group_by(d_var) %>%
  summarise(
    median_fs_p = median(fs_p, na.rm = TRUE),
    median_snv_p = median(snv_p, na.rm = TRUE),
    median_fs_FDR = median(fs_FDR, na.rm = TRUE),
    median_snv_FDR = median(snv_FDR, na.rm = TRUE),
    .groups = "drop"
  )

fdr_summary_CIDER