str(NMD_region_NMDPos_peptide_props)

#assume uniprot id and gene id are 1-1, therefore use uniprot id to match with the peptide data
#each uniprot id, randomly select one variant, repeat this process 1000 times
#compare the var-wt value of properties

#1. create the var-wt value
props = c(
  "aliphatic_index",
  "boman",
  "charge",
  "entropy",
  "instability_index",
  "isoelectric_point",
  "molecular_weight",
  "mz",
  "longest_run",
  "max_frequency",
  "energy_cost",
  "nutrient_cost",
  "mass_shift",
  "hydrophobicity",
  "hydrophobic_moment"
)

for (p in props) {
  NMD_region_NMDPos_peptide_props[[paste0("d_", p)]] =
    NMD_region_NMDPos_peptide_props[[paste0("var_", p)]] -
    NMD_region_NMDPos_peptide_props[[paste0("wt_", p)]]
}

d_vars = paste0("d_", props)

#look at the distribution of mean -> normal distribution
mean_dis = NMD_region_NMDPos_peptide_props %>%
  group_by(uniprot_id, source_folder) %>%
  summarise(across(all_of(d_vars), mean, na.rm = TRUE), .groups = "drop") 
fs_control_boman = mean_dis %>%
  filter(source_folder == "fs_control") %>%
  pull(d_boman)
hist(fs_control_boman, breaks = 30, main = "Distribution of mean d_boman for fs_control", xlab = "mean d_boman")
shapiro.test(fs_control_boman)  # if p > 0.05, normality holds
#look at the distribution within each uniprot id
NMD_region_NMDPos_peptide_props %>%
  group_by(uniprot_id, source_folder) %>% #
  

#2. randomly select one variant for each uniprot id per group, repeat this process 100 times
set.seed(123)
sampled_data = NMD_region_NMDPos_peptide_props %>%
  group_by(uniprot_id, source_folder) %>%
  sample_n(1) %>%
  ungroup()

#3. use paired Wilcoxon to compare the d value(d_aliphatic_index, d_boman) from fs vs fs_control and snv vs snv_control, match by uniprot id, get p value, adjust multiple test by BH
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

all_p_results = map_dfr(1:1000, function(i) {
  
  sampled_data = NMD_region_NMDPos_peptide_props %>%
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

all_p_results

#4. write output: d_var, median_fs_FDR, median_snv_FDR
fdr_summary = all_p_results %>%
  group_by(d_var) %>%
  summarise(
    median_fs_p = median(fs_p, na.rm = TRUE),
    median_snv_p = median(snv_p, na.rm = TRUE),
    median_fs_FDR = median(fs_FDR, na.rm = TRUE),
    median_snv_FDR = median(snv_FDR, na.rm = TRUE),
    .groups = "drop"
  )

fdr_summary
