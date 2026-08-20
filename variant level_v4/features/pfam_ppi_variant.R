#This R script is to analyze and plot pfam and ppi on variant level
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
variants_all$key
variants_all$transcript
variants_all$cds_ptc_loc

 pfam_fin$pfam_start
 pfam_fin$pfam_end
 pfam_fin$ensembl_transcript_id
 
 ###--------------------------------------------------
 ### 0. helper function: convert interface residue string to numeric vector
 ### example input: "[12,34,56]"
 ###--------------------------------------------------
 convert_to_c <- function(x) {
   if (is.na(x) || x == "") return(numeric(0))
   
   x <- gsub("\\[|\\]", "", x)
   x <- gsub(" ", "", x)
   
   if (x == "") return(numeric(0))
   
   vals <- unlist(strsplit(x, ","))
   vals <- suppressWarnings(as.numeric(vals))
   vals <- vals[!is.na(vals)]
   
   return(vals)
 }
 
 ###--------------------------------------------------
 ### 1. make sure your input data has these columns:
 ### variants_all$uniprot
 ### variants_all$cds_ptc_loc
 ### variants_all$source
 ### and human_1_ has:
 ### human_1_$uniprot1, human_1_$uniprot2
 ### human_1_$interface_residues1, human_1_$interface_residues2
 ###--------------------------------------------------
 
 variants_all <- variants_all %>%
   mutate(
     matched_uniprot = 0,
     nearest_interface_bp = NA_real_,
     dist_to_nearest_interface_bp = NA_real_
   )
 
 ###--------------------------------------------------
 ### 2. compare interface residues with cds_ptc_loc
 ### logic kept the same as before:
 ### if any interface residue position (converted to bp) >= cds_ptc_loc,
 ### then matched_uniprot = 1
 ###--------------------------------------------------
 for (i in seq_len(nrow(variants_all))) {
   
   uid <- variants_all$uniprot[i]
   cds_ptc_loc <- variants_all$cds_ptc_loc[i]
   
   if (is.na(uid) || uid == "") next
   if (is.na(cds_ptc_loc)) next
   
   re_1 <- unlist(lapply(
     human_1_$interface_residues1[human_1_$uniprot1 == uid],
     convert_to_c
   )) * 3
   
   re_2 <- unlist(lapply(
     human_1_$interface_residues2[human_1_$uniprot2 == uid],
     convert_to_c
   )) * 3
   
   all_interface_bp <- c(re_1, re_2)
   all_interface_bp <- all_interface_bp[!is.na(all_interface_bp)]
   
   if (length(all_interface_bp) == 0) next
   
   downstream_interface_bp <- all_interface_bp[all_interface_bp >= cds_ptc_loc]
   
   if (length(downstream_interface_bp) > 0) {
     variants_all$matched_uniprot[i] <- 1
     variants_all$nearest_interface_bp[i] <- min(downstream_interface_bp, na.rm = TRUE)
     variants_all$dist_to_nearest_interface_bp[i] <- 
       min(downstream_interface_bp, na.rm = TRUE) - cds_ptc_loc
   }
 }
 
 ###--------------------------------------------------
 ### 3. basic check
 ###--------------------------------------------------
 table(variants_all$matched_uniprot, useNA = "ifany")
 summary(variants_all$dist_to_nearest_interface_bp)
 
 #rename matched_uniprot to variant_ppi_overlap
 variants_all <- variants_all %>%
   rename(variant_ppi_overlap = matched_uniprot,
          variant_ppi_nearest_interface_bp = nearest_interface_bp,
          variant_ppi_dist_to_nearest_interface_bp = dist_to_nearest_interface_bp)
 
 ###--------------------------------------------------
 ### 4. set group order and colors
 ###--------------------------------------------------
 variants_all$source <- factor(
   variants_all$source,
   levels = c("snv", "snv_control", "fs", "fs_control")
 )
 
 group_colors2 <- c(
   "snv" = "#2ca02c",
   "snv_control" = "#98df8a",
   "fs" = "#1f77b4",
   "fs_control" = "#aec7e8"
 )
 
 comparisons <- list(
   c("snv", "snv_control"),
   c("fs", "fs_control")
 )
 
 ###--------------------------------------------------
 ### 5. plot 1: proportion with downstream interface residue
 ###--------------------------------------------------
 ppi_flag_summary <- variants_all %>%
   group_by(source) %>%
   summarise(
     n = n(),
     n_match = sum(matched_uniprot == 1, na.rm = TRUE),
     prop_match = n_match / n,
     .groups = "drop"
   )
 
 tab_snv <- matrix(c(
   ppi_flag_summary$n_match[ppi_flag_summary$source == "snv"],
   ppi_flag_summary$n[ppi_flag_summary$source == "snv"] -
     ppi_flag_summary$n_match[ppi_flag_summary$source == "snv"],
   ppi_flag_summary$n_match[ppi_flag_summary$source == "snv_control"],
   ppi_flag_summary$n[ppi_flag_summary$source == "snv_control"] -
     ppi_flag_summary$n_match[ppi_flag_summary$source == "snv_control"]
 ), nrow = 2, byrow = TRUE)
 
 tab_fs <- matrix(c(
   ppi_flag_summary$n_match[ppi_flag_summary$source == "fs"],
   ppi_flag_summary$n[ppi_flag_summary$source == "fs"] -
     ppi_flag_summary$n_match[ppi_flag_summary$source == "fs"],
   ppi_flag_summary$n_match[ppi_flag_summary$source == "fs_control"],
   ppi_flag_summary$n[ppi_flag_summary$source == "fs_control"] -
     ppi_flag_summary$n_match[ppi_flag_summary$source == "fs_control"]
 ), nrow = 2, byrow = TRUE)
 
 p_snv <- fisher.test(tab_snv)$p.value
 p_fs  <- fisher.test(tab_fs)$p.value
 
 ymax <- max(ppi_flag_summary$prop_match, na.rm = TRUE)
 
 pval_df <- data.frame(
   group1 = c("snv", "fs"),
   group2 = c("snv_control", "fs_control"),
   y.position = c(ymax * 1.08, ymax * 1.20),
   label = c(
     paste0("p = ", signif(p_snv, 3)),
     paste0("p = ", signif(p_fs, 3))
   )
 )
 
 p1 <- ggplot(ppi_flag_summary, aes(x = source, y = prop_match, fill = source)) +
   geom_col(width = 0.7) +
   scale_fill_manual(values = group_colors2, guide = "none") +
   labs(
     title = "Proportion of variants with downstream PPI interface residue",
     x = NULL,
     y = "Proportion"
   ) +
   theme_minimal(base_size = 14) +
   theme(
     axis.text.x = element_text(angle = 30, hjust = 1),
     plot.title = element_text(hjust = 0.5, face = "bold")
   ) +
   geom_text(
     aes(label = paste0(n_match, "/", n)),
     vjust = -0.3,
     size = 4
   ) +
   stat_pvalue_manual(
     pval_df,
     label = "label",
     xmin = "group1",
     xmax = "group2",
     y.position = "y.position",
     tip.length = 0.01
   ) +
   expand_limits(y = ymax * 1.30)
 
 print(p1)
 
 ###--------------------------------------------------
 ### 6. plot 2: distance from cds_ptc_loc to nearest downstream interface residue
 ### only for matched variants
 ###--------------------------------------------------
 p2 <- ggplot(
   variants_all %>% filter(!is.na(dist_to_nearest_interface_bp)),
   aes(x = source, y = dist_to_nearest_interface_bp, fill = source)
 ) +
   geom_violin(trim = TRUE, alpha = 0.7) +
   geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
   scale_fill_manual(values = group_colors2, guide = "none") +
   labs(
     title = "Distance from PTC to nearest downstream PPI interface residue",
     x = NULL,
     y = "Distance (bp)"
   ) +
   theme_minimal(base_size = 14) +
   theme(
     axis.text.x = element_text(angle = 30, hjust = 1),
     plot.title = element_text(hjust = 0.5, face = "bold")
   ) +
   stat_compare_means(
     comparisons = comparisons,
     method = "wilcox.test",
     label = "p.format"
   )
 
 print(p2)
 
 ###--------------------------------------------------
 ### 7. plot 3: log10 distance
 ### only for matched variants with distance > 0
 ###--------------------------------------------------
 p3 <- ggplot(
   variants_all %>% filter(!is.na(dist_to_nearest_interface_bp), dist_to_nearest_interface_bp > 0),
   aes(x = source, y = dist_to_nearest_interface_bp, fill = source)
 ) +
   geom_violin(trim = TRUE, alpha = 0.7) +
   geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
   scale_y_log10() +
   scale_fill_manual(values = group_colors2, guide = "none") +
   labs(
     title = "Log10 distance from PTC to nearest downstream PPI interface residue",
     x = NULL,
     y = "Distance (bp, log10)"
   ) +
   theme_minimal(base_size = 14) +
   theme(
     axis.text.x = element_text(angle = 30, hjust = 1),
     plot.title = element_text(hjust = 0.5, face = "bold")
   ) +
   stat_compare_means(
     comparisons = comparisons,
     method = "wilcox.test",
     label = "p.format"
   )
 
 print(p3)
 
 ###--------------------------------------------------
 ### 8. save results
 ###--------------------------------------------------
 write.csv(variants_all, "variants_all0407", row.names = FALSE)
 
 
 pfam_domains2 <- getBM(
   attributes = c("ensembl_transcript_id", "pfam_end"),
   filters    = "ensembl_transcript_id",
   values     = unique(variants_all$transcript),
   mart       = ensembl
 )
 
 max_pfam_end_df <- pfam_domains2 %>%
   filter(!is.na(pfam_end)) %>%
   group_by(ensembl_transcript_id) %>%
   summarise(max_pfam_end = max(pfam_end), .groups = "drop")
 
 variants_all <- variants_all %>%
   mutate(ptc_aa = ceiling(cds_ptc_loc / 3)) %>%
   left_join(max_pfam_end_df, by = c("transcript" = "ensembl_transcript_id")) %>%
   mutate(
     dist_ptc_to_max_pfam_end_aa = ptc_aa - max_pfam_end,
     ptc_after_max_pfam_end = ifelse(dist_ptc_to_max_pfam_end_aa > 0, 1, 0),
     ptc_before_max_pfam_end = ifelse(dist_ptc_to_max_pfam_end_aa < 0, 1, 0)
   )
  variants_all$ptc_before_max_pfam_end = variants_all %>% mutate(ptc_before_max_pfam_end = ifelse(dist_ptc_to_max_pfam_end_aa < 0, 1, 0)) %>% pull(ptc_before_max_pfam_end)
 
  write.csv(variants_all, "variants_all_with_pfam_distance.csv", row.names = FALSE)
 library(dplyr)
 library(ggplot2)
 library(ggpubr)
 
 # set group order
 variants_all$source <- factor(
   variants_all$source,
   levels = c("snv", "snv_control", "fs", "fs_control")
 )
 
 group_colors2 <- c(
   "snv" = "#2ca02c",
   "snv_control" = "#98df8a",
   "fs" = "#1f77b4",
   "fs_control" = "#aec7e8"
 )
 
 comparisons <- list(
   c("snv", "snv_control"),
   c("fs", "fs_control")
 )
 
 ### 1. violin + boxplot for distance to largest PFAM end
 p1 <- ggplot(
   variants_all %>% filter(!is.na(dist_ptc_to_max_pfam_end_aa)),
   aes(x = source, y = dist_ptc_to_max_pfam_end_aa, fill = source)
 ) +
   geom_violin(trim = TRUE, alpha = 0.7) +
   geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
   scale_fill_manual(values = group_colors2, guide = "none") +
   labs(
     title = "Distance from PTC to the largest PFAM end",
     x = NULL,
     y = "PTC aa position - largest PFAM end (aa)"
   ) +
   theme_minimal(base_size = 14) +
   theme(
     axis.text.x = element_text(angle = 30, hjust = 1),
     plot.title = element_text(hjust = 0.5, face = "bold")
   ) +
   stat_compare_means(
     comparisons = comparisons,
     method = "wilcox.test",
     label = "p.format"
   )
 
 print(p1)
 
 ### 2. violin + boxplot for absolute distance
 p2 <- ggplot(
   variants_all %>%
     filter(!is.na(dist_ptc_to_max_pfam_end_aa)) %>%
     mutate(abs_dist_ptc_to_max_pfam_end_aa = abs(dist_ptc_to_max_pfam_end_aa)),
   aes(x = source, y = abs_dist_ptc_to_max_pfam_end_aa, fill = source)
 ) +
   geom_violin(trim = TRUE, alpha = 0.7) +
   geom_boxplot(width = 0.12, outlier.size = 0.6, color = "black") +
   scale_fill_manual(values = group_colors2, guide = "none") +
   labs(
     title = "Absolute distance from PTC to the largest PFAM end",
     x = NULL,
     y = "Absolute distance (aa)"
   ) +
   theme_minimal(base_size = 14) +
   theme(
     axis.text.x = element_text(angle = 30, hjust = 1),
     plot.title = element_text(hjust = 0.5, face = "bold")
   ) +
   stat_compare_means(
     comparisons = comparisons,
     method = "wilcox.test",
     label = "p.format"
   )
 
 print(p2)
 
 ### 3. proportion barplot for PTC after largest PFAM end
 pfam_flag_summary <- variants_all %>%
   group_by(source) %>%
   summarise(
     n = n(),
     n_after = sum(ptc_after_max_pfam_end == 1, na.rm = TRUE),
     n_before = sum(ptc_before_max_pfam_end == 1, na.rm = TRUE),
     prop_after = n_after / n,
     prop_before = n_before / n,
     .groups = "drop"
   )
 
 # fisher tests
 tab_snv <- matrix(c(
   pfam_flag_summary$n_after[pfam_flag_summary$source == "snv"],
   pfam_flag_summary$n[pfam_flag_summary$source == "snv"] -
     pfam_flag_summary$n_after[pfam_flag_summary$source == "snv"],
   pfam_flag_summary$n_after[pfam_flag_summary$source == "snv_control"],
   pfam_flag_summary$n[pfam_flag_summary$source == "snv_control"] -
     pfam_flag_summary$n_after[pfam_flag_summary$source == "snv_control"]
 ), nrow = 2, byrow = TRUE)
 
 tab_fs <- matrix(c(
   pfam_flag_summary$n_after[pfam_flag_summary$source == "fs"],
   pfam_flag_summary$n[pfam_flag_summary$source == "fs"] -
     pfam_flag_summary$n_after[pfam_flag_summary$source == "fs"],
   pfam_flag_summary$n_after[pfam_flag_summary$source == "fs_control"],
   pfam_flag_summary$n[pfam_flag_summary$source == "fs_control"] -
     pfam_flag_summary$n_after[pfam_flag_summary$source == "fs_control"]
 ), nrow = 2, byrow = TRUE)
 
 p_snv <- fisher.test(tab_snv)$p.value
 p_fs  <- fisher.test(tab_fs)$p.value
 
 ymax <- max(pfam_flag_summary$prop_after, na.rm = TRUE)
 
 pval_df <- data.frame(
   group1 = c("snv", "fs"),
   group2 = c("snv_control", "fs_control"),
   y.position = c(ymax * 1.08, ymax * 1.20),
   label = c(
     paste0("p = ", signif(p_snv, 3)),
     paste0("p = ", signif(p_fs, 3))
   )
 )
 
 p3 <- ggplot(pfam_flag_summary, aes(x = source, y = prop_before, fill = source)) +
   geom_col(width = 0.7) +
   scale_fill_manual(values = group_colors2, guide = "none") +
   labs(
     title = "Proportion of variants with PTC influencing PFAM domain",
     x = NULL,
     y = "Proportion"
   ) +
   theme_minimal(base_size = 14) +
   theme(
     axis.text.x = element_text(angle = 30, hjust = 1),
     plot.title = element_text(hjust = 0.5, face = "bold")
   ) +
   geom_text(
     aes(label = paste0(n_before, "/", n)),
     vjust = -0.3,
     size = 4
   ) +
   stat_pvalue_manual(
     pval_df,
     label = "label",
     xmin = "group1",
     xmax = "group2",
     y.position = "y.position",
     tip.length = 0.01
   ) +
   expand_limits(y = ymax * 1.30)
 
 print(p3)
 
 
 #This function transform cds loc to aa loc, which is corresponding to pfam loc
 cds2aa = function(base_idx) ((base_idx - 1L) %/% 3L) + 1L
 
 #This function combine snv_dis and snv_variants and pfam_fin to compare pfam loc with variant loc
 compare_pfam_mut = function(dis_df = snv_dis, var_df = snv_variants){
 
   ##2. join variants_all and pfam_fin by transcript id, so that pfam domain information is added
   variants_all3 <- variants_all %>%
     dplyr::left_join(
       pfam_fin %>% dplyr::select(ensembl_transcript_id, pfam, pfam_start, pfam_end, hgnc_symbol),
       by = c("transcript" = "ensembl_transcript_id")
     )

   ##3. compare mut_aa_pos with pfam_start/end -> new column "in_pfam"
   new_df3 <- new_df2 %>%
     dplyr::mutate(
       in_pfam = dplyr::if_else(!is.na(pfam_start) & !is.na(pfam_end) &
                                  !is.na(mut_aa_pos) &
                                  #mut_aa_pos >= pfam_start & ,
                                  mut_aa_pos <= pfam_end,
                                TRUE, FALSE)
     )
   return(new_df3)
 }
 
 #compare pfam end with var loc
 pfam_all_pervariant <- bind_rows(
   compare_pfam_mut(minus1_dis,         minus1_variants)         %>% mutate(group = "Minus1"),
   compare_pfam_mut(minus1_control_dis, minus1_control_variants) %>% mutate(group = "Minus1_Control"),
   compare_pfam_mut(plus1_dis,          plus1_variants)          %>% mutate(group = "Plus1"),
   compare_pfam_mut(plus1_control_dis,  plus1_control_variants)  %>% mutate(group = "Plus1_Control"),
   compare_pfam_mut(snv_dis,            snv_variants)            %>% mutate(group = "Nonsense"),
   compare_pfam_mut(snv_control_dis,    snv_control_variants)    %>% mutate(group = "Nonsense_Control")
 )

 ##4. ggplot the comparison result
 # overall counts (% inside PFAM vs outside)
 pfam_palette <- c(`TRUE` = "#2ca02c", `FALSE` = "#d62728")

 pfam_percent_plot <- function(pfam_all_pervariant,
                               inside_col = "#1b9e77",   # teal (CB-friendly)
                               outside_col = "lightblue")  # orange (CB-friendly)
 {
   pfam_df <- pfam_all_pervariant %>%
     count(group, in_pfam) %>%
     group_by(group) %>%
     mutate(prop = n / sum(n),
            pct  = percent(prop, accuracy = 0.1)) %>%
     ungroup()
   
   pfam_palette <- c(`FALSE` = outside_col, `TRUE` = inside_col)
   write.csv(prop_df, "plots/pfam_compare/pfam_influenced_variant_counts.csv", row.names = FALSE)
   ggplot(prop_df, aes(x = group, y = prop, fill = as.character(in_pfam))) +
     geom_col(color = "grey30") +
     geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 3) +
     scale_fill_manual(values = pfam_palette,
                       labels = c("Outside PFAM","Influence PFAM"), name = NULL) +
     scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
     labs(x = NULL, y = "Fraction of variants",
          title = "% of variants which PFAM domains were influenced") +
     theme_minimal(base_size = 12) +
     theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
     theme(plot.title = element_text(hjust = 0.5, face = "bold"))
 }

 p <- pfam_percent_plot(pfam_all_pervariant)
 print(p)
 
 write.csv(prop_df, "plots/pfam_compare/pfam_variant_counts.csv", row.names = FALSE)
 
 
 #ppi

 #This function clean the format of ppi interface residues
 convert_to_ints = function(x) {
   if (is.null(x) || is.na(x) || x == "") return(integer(0))
   s <- gsub("\\[|\\]", "", x)
   s <- trimws(s)
   if (s == "") return(integer(0))
   as.integer(trimws(unlist(strsplit(s, ","))))
 }
 
 #This function reformed human_1 to a list of uniprot ids and their interface residues
 clean_human_1 = function(human_df) {
   stopifnot(all(c("uniprot1","uniprot2","interface_residues1","interface_residues2") %in% names(human_df)))
   
   d1 <- human_df %>%
     transmute(uniprot = uniprot1,
               residues = lapply(interface_residues1, convert_to_ints))
   d2 <- human_df %>%
     transmute(uniprot = uniprot2,
               residues = lapply(interface_residues2, convert_to_ints))
   bind_rows(d1, d2) %>%
     group_by(uniprot) %>%
     summarise(residues = list(sort(unique(unlist(residues)))), .groups = "drop") %>%
     { setNames(.$residues, .$uniprot) }
 }
 
 #This function is to join the distance df with variant df, and compare mut_aa_pos with max(ppi_residue)
 ppi_overlap_for_group <- function(group_name, variants_tbl, distance_tbl, iface_index) {
   v <- distance_tbl %>%
     dplyr::left_join(variants_tbl, by = c("Variant_Key" = "key")) %>%
     dplyr::mutate(aa_pos = cds2aa(as.integer(cds_mutation_loc))) 
   
   # flag inside interface by mut_aa_pos <= max(ppi_residue)
   in_iface_vec <- mapply(function(u, pos) {
     if (is.na(u) || !(u %in% names(iface_index)) || is.na(pos)) return(FALSE)
     iface_max <- max(iface_index[[u]], na.rm = TRUE)
     pos <= iface_max
   }, v$uniprot, v$aa_pos)
   
   v %>%
     transmute(Variant_Key, group = group_name, in_interface = as.logical(in_iface_vec)) %>%
     distinct()
 }
 
 
 
 
 

 #Interface index from human_1_
 iface_index <- clean_human_1(human_1_)
 #Run for each group
 v1 = ppi_overlap_for_group("Minus1",           minus1_variants,         minus1_dis,         iface_index)
 v2 = ppi_overlap_for_group("Minus1_Control",   minus1_control_variants, minus1_control_dis, iface_index)
 v3 = ppi_overlap_for_group("Plus1",            plus1_variants,          plus1_dis,          iface_index)
 v4 = ppi_overlap_for_group("Plus1_Control",    plus1_control_variants,  plus1_control_dis,  iface_index)
 v5 = ppi_overlap_for_group("Nonsense",         snv_variants,            snv_dis,            iface_index)
 v6 = ppi_overlap_for_group("Nonsense_Control", snv_control_variants,    snv_control_dis,    iface_index)
 ppi_all <- bind_rows(
   v1, v2, v3, v4, v5, v6
 )
 
 ppi_all$group <- factor(ppi_all$group,
                         levels = c("Minus1","Minus1_Control","Plus1","Plus1_Control","Nonsense","Nonsense_Control"))
 
 #select AD part
 ppi_all_AD = ppi_all[which(ppi_all$Variant_Key %in% ppi_all_AD$Variant_Key), ]
 #plot the ppi result
 pfam_like_palette <- c(`FALSE` = "lightblue", `TRUE` = "#2ca02c")  # outside red, inside green
 
 ppi_df <- ppi_all %>%
   count(group, in_interface) %>%
   group_by(group) %>%
   tidyr::complete(in_interface = c(FALSE, TRUE), fill = list(n = 0)) %>%
   mutate(prop = n / sum(n),
          pct  = scales::percent(prop, accuracy = 0.1)) %>%
   ungroup()
 ppi_df_AD <- ppi_all_AD %>%
   count(group, in_interface) %>%
   group_by(group) %>%
   tidyr::complete(in_interface = c(FALSE, TRUE), fill = list(n = 0)) %>%
   mutate(prop = n / sum(n),
          pct  = scales::percent(prop, accuracy = 0.1)) %>%
   ungroup()
 write.csv(ppi_df, "ppi_interface_variant_counts.csv", row.names = FALSE)
 write.csv(ppi_df_AD, "ppi_interface_variant_counts_AD.csv", row.names = FALSE)
 ppi_all$hgnc_symbol = pfam_all_pervariant$hgnc_symbol[match(ppi_all$Variant_Key, pfam_all_pervariant$Variant_Key)]
 ppi_all_AD$hgnc_symbol = pfam_all_pervariant$hgnc_symbol[match(ppi_all_AD$Variant_Key, pfam_all_pervariant$Variant_Key)]
 write.csv(ppi_all_AD, "ppi_all_AD_with_gene_symbols.csv", row.names = FALSE)
 p_ppi_percent <- ggplot(ppi_df_AD,
                         aes(x = group, y = prop, fill = as.character(in_interface))) +
   geom_col(color = "grey30") +
   geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 3) +
   scale_fill_manual(values = pfam_like_palette,
                     labels = c("Outside PPI interface","Inside PPI interface"),
                     name = NULL) +
   scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
   labs(x = NULL, y = "Fraction of variants",
        title = "% of variants inside PPI interfaces vs outside") +
   theme_minimal(base_size = 12) +
   theme(axis.text.x = element_text(angle = 30, hjust = 1),
         plot.title  = element_text(hjust = 0.5))  # centered title
 p_ppi_percent +
   theme(plot.title = element_text(hjust = 0.5, face = "bold"))
 print(p_ppi_percent)
 
 # (Optional) save
 dir.create("plots/ppi_compare", recursive = TRUE, showWarnings = FALSE)
 ggsave("plots/ppi_compare/ppi_inside_vs_outside_stacked.pdf", p_ppi_percent, width = 9, height = 4.5)
 ggsave("plots/ppi_compare/ppi_inside_vs_outside_stacked.png", p_ppi_percent, width = 9, height = 4.5, dpi = 300)
