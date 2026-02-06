#input data
library(stringr)
library(dplyr)
library(readxl)
touni <- read_csv("~/Downloads/NIHMS1818854-supplement-2(A).csv")
motif_doc<- read_csv("~/Downloads/NIHMS1818854-supplement-2(B).csv")
LCS_doc <- read_excel("~/Desktop/Copy of NIHMS1818854-supplement-2.xls", 
                                                     sheet = "G")
#select largest value in each cell
get_max_from_cell <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  nums <- str_extract_all(x, "\\d+")[[1]]
  if (length(nums) == 0) return(NA_real_)
  max(as.numeric(nums), na.rm = TRUE)
}

motif_max <- motif_doc %>%
  mutate(across(-1, ~ vapply(as.character(.x), get_max_from_cell, numeric(1))))
LCS_max <- LCS_doc %>%
  mutate(across(-1, ~ vapply(as.character(.x), get_max_from_cell, numeric(1))))
motif_max$uniprot = touni$`Uniprot ID`[match(motif_max$Protein, touni$Protein)]
LCS_max$uniprot = touni$`Uniprot ID`[match(LCS_max$Protein, touni$Protein)]

#combine motif value with variant loc and nmdesc loc, match by uniprot id
#this is protein level, which is gene level, so for the variant list, we compare with variant loc
#motif_A$uniprot_ID
motif_max$uniprot

merge(motif_max,pfam_all_pervariant, by.x = "uniprot", by.y = "uniprotswissprot") -> motif_max3
nmdesc_loc  = rbind(snv_idr[,c('Variant_Key','NMDesc.start')],
                    snv_control_idr[,c('Variant_Key','NMDesc.start')],
                    plus1_idr[,c('Variant_Key','NMDesc.start')],
                    plus1_control_idr[,c('Variant_Key','NMDesc.start')],
                    minus1_idr[,c('Variant_Key','NMDesc.start')],
                    minus1_control_idr[,c('Variant_Key','NMDesc.start')])

motif_max3$nmdesc_start <- nmdesc_loc$NMDesc.start[match(motif_max3$Variant_Key, nmdesc_loc$Variant_Key)]
motif_max3$protein_feature_flag = motif_max3$`Protein Features`*3 >= motif_max3$cds_mutation_loc
motif_max3$protein_flag2 = motif_max3$`Protein Features`*3 >= motif_max3$nmdesc_start
motif_max3$domains_flag = motif_max3$`Domains`*3 >= motif_max3$cds_mutation_loc
motif_max3$domains_flag2 = motif_max3$`Domains`*3 >= motif_max3$nmdesc_start
motif_max3$slim_flag = motif_max3$`SLiMs`*3 >= motif_max3$cds_mutation_loc #variant location on cds scale
motif_max3$slim_flag2 = motif_max3$`SLiMs`*3 >= motif_max3$nmdesc_start
motif_max3$morf_flag = motif_max3$MORFs*3 >= motif_max3$cds_mutation_loc
motif_max3$morf_flag2 = motif_max3$MORFs*3 >= motif_max3$nmdesc_start
motif_max3$ptm_flag = motif_max3$`PTMs`*3 >= motif_max3$cds_mutation_loc
motif_max3$ptm_flag2 = motif_max3$`PTMs`*3 >= motif_max3$nmdesc_start
motif_max3$nls_flag = motif_max3$`NLSs`*3 >= motif_max3$cds_mutation_loc
motif_max3$nls_flag2 = motif_max3$`NLSs`*3 >= motif_max3$nmdesc_start
LCS_max3 = merge(LCS_max, pfam_all_pervariant, by.x = "uniprot", by.y = "uniprotswissprot")
LCS_max3$nmdesc_start <- nmdesc_loc$NMDesc.start[match(LCS_max3$Variant_Key, nmdesc_loc$Variant_Key)]
LCS_max3$LCS_flag = LCS_max3$`LCSs`*3 >= LCS_max3$cds_mutation_loc
LCS_max3$LCS_flag2 = LCS_max3$`LCSs`*3 >= LCS_max3$nmdesc_start
#filter for AD variants
LCS_max_AD = LCS_max3[which(LCS_max3$Variant_Key %in% ppi_all_AD$Variant_Key),]
motif_max3_AD = motif_max3[which(motif_max3$Variant_Key %in% ppi_all_AD$Variant_Key),]
write.csv(LCS_max_AD, "LCS_max_AD.csv", row.names = FALSE)
write.csv(motif_max3_AD, "motif_max_AD.csv", row.names = FALSE)
write.csv(LCS_max3, "LCS_max3.csv", row.names = FALSE)
write.csv(motif_max3, "motif_max3.csv", row.names = FALSE)
#select flags
#motif_flags = motif_max3 %>% select(key, uniprot, ends_with("_flag"),group)
#plot the flags by group

all_dis$cds_mutation_loc #use this one
all_dis$Variant_Key
key_to_transcript$Transcript_ID

tx_ids <- key_to_transcript$Transcript_ID

# query Ensembl for UniProt mappings
bm_raw <- getBM(
  attributes = c("ensembl_transcript_id", 
                 "uniprotswissprot"),
  filters    = "ensembl_transcript_id",
  values     = tx_ids,
  mart       = mart
)
key_to_transcript$uniprot <- bm_raw$uniprotswissprot[match(key_to_transcript$Transcript_ID, bm_raw$ensembl_transcript_id)]
to_pipe_style <- function(x, pad = 9) {
  x <- as.character(x)
  # If already in pipe style, return as-is
  if (str_detect(x, ":" ) && str_detect(x, "\\|")) return(x)
  
  parts <- str_split(x, "_", n = 4, simplify = TRUE)
  if (ncol(parts) < 4) return(NA_character_)
  
  chrom <- str_remove(parts[1], "^chr")           # drop optional "chr"
  pos   <- suppressWarnings(as.integer(parts[2]))
  ref   <- parts[3]
  alt   <- parts[4]
  
  if (is.na(pos)) return(NA_character_)
  pos_padded <- sprintf(paste0("%0", pad, "d"), pos)
  
  paste0(chrom, ":", pos_padded, "|", ref, "|", alt)
}
key_to_transcript$key2 = sapply(key_to_transcript$Variant_Key, to_pipe_style, character(1))
#dis_2 = merge(all_dis, key_to_transcript, by.x = "Variant_Key", by.y = "key2", all.x = TRUE)

#WT_var_NMD_2ac$key2 = sapply(WT_var_NMD_2ac$id, to_pipe_style, character(1))
#dis_2 = merge(all_dis, WT_var_NMD_2ac, by.x = "Variant_Key", by.y = "key2", all.x = TRUE)

dis_2 = merge(pfam_all_pervariant, motif_max3, by.x = "uniprotswissprot", by.y = "uniprot")
w

dis_3 = dis_2
dis_3 = dis_3 %>%select(key, uniprotswissprot, cds_mutation_loc, group.x, ends_with("_flag"))

#dis_3 = merge(dis_2, motif_max, by = "uniprot", all.x = TRUE)



-------
#or
all_vep = rbind(minus1_vep, minus1_control_vep, plus1_vep, plus1_control_vep) #but vep do not have snv data

all_vep$Protein_position
all_vep$`#Uploaded_variation`
all_vep$Gene
all_vep$UNIPROT_ISOFORM
motif_doc$Protein
-------------
  make_flag_plot <- function(flag_name) {
    df_flag <- mf_long %>% dplyr::filter(flag == flag_name)
    
    # TRUE rate per group
    rate_df <- df_flag %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        n_non_na  = sum(!is.na(value)),
        true_rate = ifelse(n_non_na == 0, NA_real_, mean(value, na.rm = TRUE)),
        .groups = "drop"
      )
    
    # comparisons + p-values
    comps <- tibble::tibble(
      group1 = c("Minus1","Plus1","Nonsense"),
      group2 = c("Minus1_Control","Plus1_Control","Nonsense_Control")
    ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        p = pairwise_p(df_flag, group1, group2),
        y.position = {
          y1 <- rate_df$true_rate[rate_df$group == group1]
          y2 <- rate_df$true_rate[rate_df$group == group2]
          y  <- suppressWarnings(max(c(y1, y2), na.rm = TRUE))
          if (is.infinite(y)) NA_real_ else min(1, y + 0.05)
        },
        p.signif = p_to_signif(p)
      ) %>%
      dplyr::ungroup() %>%
      # keep only valid rows
      dplyr::filter(!is.na(y.position), !is.na(p)) %>%
      dplyr::select(group1, group2, y.position, p, p.signif)
    
    ggplot(rate_df, aes(x = group, y = true_rate, fill = group)) +
      geom_col() +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      scale_fill_manual(values = group_palette) +
      scale_x_discrete(labels = x_labels) +
      labs(x = NULL, y = "Proportion TRUE", title = flag_name) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      ggpubr::stat_pvalue_manual(
        data = comps,
        label = "p.signif",   # use "p" for numeric p-values
        tip.length = 0.01,
        bracket.size = 0.5,
        step.increase = 0.08
      )
  }

#save the plots
pdf("flag_plots_all.pdf", width = 6, height = 4)
for (nm in names(plots)) print(plots[[nm]])
dev.off()


