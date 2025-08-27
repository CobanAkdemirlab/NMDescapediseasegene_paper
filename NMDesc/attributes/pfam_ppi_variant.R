library(dplyr)
library(stringr)
library(ggplot2)
library(scales)

#pfam and ppi on variant level
 pfam_fin$pfam_start
 pfam_fin$pfam_end
 pfam_fin$ensembl_transcript_id
 
 snv_variants$transcript
 snv_variants$key
 
 snv_dis$Variant_Key
 snv_dis$cds_mutation_loc
 
 #1.join snv_dis and snv_variant by key to snv_new
 #2.join snv_new and pfam_fin by transcript id to snv_new2
 #3.compare cds_mutation_loc with pfam_start and pfam_end to get a new column "in_pfam" = 1/0 in snv_new2
 #4. ggplot the comparison result
 
 
 # base→AA helper
 aa_index <- function(base_idx) ((base_idx - 1L) %/% 3L) + 1L
 
 ##1. join snv_dis and snv_variants by key -> snv_new
 snv_new <- snv_dis %>%
   dplyr::left_join(
     snv_variants %>% dplyr::select(transcript, key, uniprotswissprot),
     by = c("Variant_Key" = "key")
   ) %>%
   dplyr::mutate(aa_pos = aa_index(as.integer(cds_mutation_loc)))
 
 ##2. join snv_new and pfam_fin by transcript id -> snv_new2
 snv_new2 <- snv_new %>%
   dplyr::left_join(
     pfam_fin %>% dplyr::select(ensembl_transcript_id, pfam, pfam_start, pfam_end, hgnc_symbol),
     by = c("transcript" = "ensembl_transcript_id")
   )
 
 ##3. compare cds_mutation_loc (as aa_pos) with pfam_start/end -> new column "in_pfam"
 snv_new2 <- snv_new2 %>%
   dplyr::mutate(
     in_pfam = dplyr::if_else(!is.na(pfam_start) & !is.na(pfam_end) &
                                !is.na(aa_pos) &
                                aa_pos >= pfam_start & aa_pos <= pfam_end,
                              TRUE, FALSE)
   )
 
 # Optional: collapse to per-variant summary to avoid PFAM-duplicate rows
 snv_by_variant <- snv_new2 %>%
   dplyr::group_by(Variant_Key, transcript, aa_pos, uniprotswissprot, hgnc_symbol) %>%
   dplyr::summarise(
     in_pfam = any(in_pfam, na.rm = TRUE),
     pfam_hits = paste(unique(na.omit(ifelse(in_pfam, pfam, NA_character_))), collapse = ";"),
     .groups = "drop"
   )
 
 ##4. ggplot the comparison result
 # overall counts (% inside PFAM vs outside)
 pfam_palette <- c(`TRUE` = "#2ca02c", `FALSE` = "#d62728")
 
 count_df <- snv_by_variant %>%
   dplyr::count(in_pfam) %>%
   dplyr::mutate(prop = n / sum(n),
                 label = ifelse(in_pfam, "Inside PFAM", "Outside PFAM"))
 
 p_overall <- ggplot(count_df, aes(x = label, y = prop, fill = as.character(in_pfam))) +
   geom_col(color = "grey30") +
   geom_text(aes(label = percent(prop, accuracy = 0.1)), vjust = -0.2, size = 3) +
   scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
   scale_fill_manual(values = pfam_palette, guide = "none") +
   labs(x = NULL, y = "Fraction of variants",
        title = "Variants inside PFAM domains (AA-based)") +
   theme_minimal(base_size = 12)
 
 print(p_overall)
 library(dplyr); library(ggplot2); library(scales)
 
 pfam_percent_plot <- function(pfam_all_pervariant,
                               inside_col = "#1b9e77",   # teal (CB-friendly)
                               outside_col = "#d95f02")  # orange (CB-friendly)
 {
   prop_df <- pfam_all_pervariant %>%
     count(group, in_pfam) %>%
     group_by(group) %>%
     mutate(prop = n / sum(n),
            pct  = percent(prop, accuracy = 0.1)) %>%
     ungroup()
   
   pfam_palette <- c(`FALSE` = outside_col, `TRUE` = inside_col)
   
   ggplot(prop_df, aes(x = group, y = prop, fill = as.character(in_pfam))) +
     geom_col(color = "grey30") +
     geom_text(aes(label = pct), position = position_stack(vjust = 0.5), size = 3) +
     scale_fill_manual(values = pfam_palette,
                       labels = c("Outside PFAM","Inside PFAM"), name = NULL) +
     scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
     labs(x = NULL, y = "Fraction of variants",
          title = "% of variants inside PFAM vs outside (by group)") +
     theme_minimal(base_size = 12) +
     theme(axis.text.x = element_text(angle = 30, hjust = 1))
 }
 # teal/orange (default)
 p <- pfam_percent_plot(pfam_all_pervariant)
 
 # your previous green/red scheme
 p <- pfam_percent_plot(pfam_all_pervariant,
                        inside_col = "#2ca02c",
                        outside_col = "#d62728")
 
 # or any custom pair
 p <- pfam_percent_plot(pfam_all_pervariant,
                        inside_col = "#4CAF50",
                        outside_col = "#6B7280")
 print(p)
 
 
 library(dplyr)
 library(tidyr)
 library(ggplot2)
 library(scales)
 
 # 1) Rename groups (case-insensitive) and set order
 pfam_all_pervariant2 <- pfam_all_pervariant %>%
   mutate(group = case_when(
     tolower(group) == "minus1"            ~ "Minus1",
     tolower(group) == "minus1_control"    ~ "Minus1_Control",
     tolower(group) == "plus1"             ~ "Plus1",
     tolower(group) == "plus1_control"     ~ "Plus1_Control",
     tolower(group) %in% c("snv","nonsense")                 ~ "Nonsense",
     tolower(group) %in% c("snv_control","nonsense_control") ~ "Nonsense_Control",
     TRUE ~ as.character(group)
   )) %>%
   mutate(group = factor(group,
                         levels = c("Minus1","Minus1_Control",
                                    "Plus1","Plus1_Control",
                                    "Nonsense","Nonsense_Control")))
 
 # 2) Compute % inside vs outside (ensure zeros show up)
 prop_df <- pfam_all_pervariant2 %>%
   count(group, in_pfam) %>%
   complete(group,
            in_pfam = c(FALSE, TRUE),
            fill = list(n = 0)) %>%
   group_by(group) %>%
   mutate(prop = n / sum(n),
          pct  = percent(prop, accuracy = 0.1)) %>%
   ungroup()
 
 # 3) Plot (% inside PFAM vs outside), green=inside, red=outside
 pfam_palette <- c(`FALSE` = "lightblue", `TRUE` = "#2ca02c")
 
 p_pfam_percent <- ggplot(prop_df,
                          aes(x = group, y = prop, fill = as.character(in_pfam))) +
   geom_col(color = "grey30") +
   geom_text(aes(label = pct),
             position = position_stack(vjust = 0.5), size = 3) +
   scale_fill_manual(values = pfam_palette,
                     labels = c("Outside PFAM","Inside PFAM"), name = NULL) +
   scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
   labs(x = NULL, y = "Fraction of variants",
        title = "% of variants inside PFAM vs outside") +
   theme_minimal(base_size = 12) +
   theme(axis.text.x = element_text(angle = 30, hjust = 1))
 
 p_pfam_percent +
     theme(plot.title = element_text(hjust = 0.5, face = "bold"))
 
 
 #ppi
 suppressPackageStartupMessages({
   library(dplyr); library(tidyr); library(purrr)
   library(stringr); library(ggplot2); library(scales)
 })
 
 # =========================
 # PPI interface helper bits
 # =========================
 
 # parse "[12, 34, 56]" -> integer vector c(12,34,56); robust to "", NA
 convert_to_ints <- function(x) {
   if (is.null(x) || is.na(x) || x == "") return(integer(0))
   s <- gsub("\\[|\\]", "", x)
   s <- trimws(s)
   if (s == "") return(integer(0))
   as.integer(trimws(unlist(strsplit(s, ","))))
 }
 
 # Build UniProt -> set of interface AA positions (union over all interactions)
 make_interface_index <- function(human_df) {
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
 
 # base-index -> amino-acid index (1-based)
 aa_index <- function(base_idx) ((base_idx - 1L) %/% 3L) + 1L
 
 # Map variants to UniProt (prefer variants_tbl$uniprotswissprot; fallback via transcript mapping if provided)
 ensure_uniprot <- function(df, tr2uni = NULL) {
   if ("uniprotswissprot" %in% names(df) && any(!is.na(df$uniprotswissprot) & df$uniprotswissprot != "")) {
     df %>% mutate(uniprot = uniprotswissprot)
   } else if (!is.null(tr2uni) && all(c("transcript","uniprotswissprot") %in% names(tr2uni))) {
     df %>%
       left_join(tr2uni %>% select(transcript, uniprotswissprot), by = "transcript") %>%
       mutate(uniprot = uniprotswissprot)
   } else {
     warning("No UniProt mapping found; in_interface will be FALSE for all.")
     df %>% mutate(uniprot = NA_character_)
   }
 }
 
 # One-group overlap (Inside PPI interface vs Outside)
 ppi_overlap_for_group <- function(group_name, variants_tbl, distance_tbl, iface_index, tr2uni = NULL) {
   stopifnot(all(c("Variant_Key","cds_mutation_loc") %in% names(distance_tbl)))
   stopifnot(all(c("transcript","key") %in% names(variants_tbl)))
   
   v <- distance_tbl %>%
     left_join(select(variants_tbl, transcript, key, everything()),
               by = c("Variant_Key" = "key")) %>%
     mutate(aa_pos = aa_index(as.integer(cds_mutation_loc))) %>%
     ensure_uniprot(tr2uni = tr2uni)
   
   # flag inside interface (treat unknown UniProt as outside = FALSE)
   in_iface_vec <- mapply(function(u, pos) {
     if (is.na(u) || !(u %in% names(iface_index)) || is.na(pos)) return(FALSE)
     pos %in% iface_index[[u]]
   }, v$uniprot, v$aa_pos)
   
   v %>%
     transmute(Variant_Key, group = group_name, in_interface = as.logical(in_iface_vec)) %>%
     distinct()
 }
 
 # =========================
 # Build the per-variant table for all 6 groups
 # =========================
 
 # 1) Interface index from human_1_
 iface_index <- make_interface_index(human_1_)
 
 # 2) (Optional) transcript -> UniProt mapping if some variants tables lack uniprotswissprot
 #    Uncomment and fill if needed:
 # tr2uni <- getBM(attributes = c("ensembl_transcript_id","uniprotswissprot"),
 #                 filters    = "ensembl_transcript_id",
 #                 values     = unique(c(plus1_variants$transcript,
 #                                       plus1_control_variants$transcript,
 #                                       minus1_variants$transcript,
 #                                       minus1_control_variants$transcript,
 #                                       snv_variants$transcript,
 #                                       snv_control_variants$transcript)),
 #                 mart = ensembl) %>%
 #           rename(transcript = ensembl_transcript_id)
 
 tr2uni <- NULL  # set to the mapping above if needed
 
 # 3) Run for each group
 ppi_all <- bind_rows(
   ppi_overlap_for_group("Minus1",           minus1_variants,         minus1_dis,         iface_index, tr2uni),
   ppi_overlap_for_group("Minus1_Control",   minus1_control_variants, minus1_control_dis, iface_index, tr2uni),
   ppi_overlap_for_group("Plus1",            plus1_variants,          plus1_dis,          iface_index, tr2uni),
   ppi_overlap_for_group("Plus1_Control",    plus1_control_variants,  plus1_control_dis,  iface_index, tr2uni),
   ppi_overlap_for_group("Nonsense",         snv_variants,            snv_dis,            iface_index, tr2uni),
   ppi_overlap_for_group("Nonsense_Control", snv_control_variants,    snv_control_dis,    iface_index, tr2uni)
 )
 
 ppi_all$group <- factor(ppi_all$group,
                         levels = c("Minus1","Minus1_Control","Plus1","Plus1_Control","Nonsense","Nonsense_Control"))
 
 # =========================
 # Plot: % inside PPI interface vs outside (stacked)
 # =========================
 
 pfam_like_palette <- c(`FALSE` = "lightblue", `TRUE` = "#2ca02c")  # outside red, inside green
 
 prop_df <- ppi_all %>%
   count(group, in_interface) %>%
   group_by(group) %>%
   complete(in_interface = c(FALSE, TRUE), fill = list(n = 0)) %>%
   mutate(prop = n / sum(n),
          pct  = percent(prop, accuracy = 0.1)) %>%
   ungroup()
 
 p_ppi_percent <- ggplot(prop_df,
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
 
 

 