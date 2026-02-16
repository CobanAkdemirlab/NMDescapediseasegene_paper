
 library(dplyr)
#This R script is to analyze and plot pfam and ppi on variant level

library(stringr)
library(ggplot2)
library(scales)

 pfam_fin$pfam_start
 pfam_fin$pfam_end
 pfam_fin$ensembl_transcript_id
 
 snv_variants$transcript
 snv_variants$key
 
 snv_dis$Variant_Key
 snv_dis$cds_mutation_loc
 
 #This function transform cds loc to aa loc, which is corresponding to pfam loc
 cds2aa = function(base_idx) ((base_idx - 1L) %/% 3L) + 1L
 
 #This function combine snv_dis and snv_variants and pfam_fin to compare pfam loc with variant loc
 compare_pfam_mut = function(dis_df = snv_dis, var_df = snv_variants){
  ##1. join snv_dis and snv_variants by key -> snv_new, so that each cds loc is matched with a key
  new_df <- dis_df %>%
     dplyr::left_join(
       var_df %>% dplyr::select(transcript, key, uniprotswissprot),
       by = c("Variant_Key" = "key")
     ) %>%
     dplyr::mutate(mut_aa_pos = cds2aa(as.integer(cds_mutation_loc)))
   
   ##2. join snv_new and pfam_fin by transcript id -> snv_new2, so that pfam domain information is added
   new_df2 <- new_df %>%
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
