############____________________________________#########
############____________________________________#########
############____________________________________#########
#Code for annotating all clinvar variants and analyzing the data
#
#
setwd("/Users/qkelly/Desktop/clinvar/0819")
##Library aenmd into R
library(aenmd.data.ensdb.v105)
library(aenmd)
library(GenomicRanges)
library(tidyr)

##Read in the downloaded and BCFtools processed Clinvar data
clinvarClinsig <- read.delim("resetID_clinvar_20240819_Clnsig.txt")

#_____________________________________________________________________________
###- 1. First we are going to retrieve how many total variants there are, then  
### how many total variants overlap with the coding regions of AENMD's default 
### transcript set.
#_____________________________________________________________________________
#Get total variants in Clinvar:
Variants_in_the_original_dataset <- nrow(clinvarClinsig)

#Make a GRange object of all the clinvar variants without filtering.
clinvar_gr_total  <- GenomicRanges::GRanges(clinvarClinsig$CHR,
                                            IRanges::IRanges(clinvarClinsig$POS,
                                                             clinvarClinsig$POS))
clinvar_gr_total$ref     <- clinvarClinsig$REF |> Biostrings::DNAStringSet()
clinvar_gr_total$alt     <- clinvarClinsig$ALT |> Biostrings::DNAStringSet()
clinvar_gr_total$id      <- clinvarClinsig$ID
clinvar_gr_total$clinsig      <- clinvarClinsig$CLNSIG

#Retrieve the base AENMD transcript set information: 
AENMD_base_exonset <- future::value(aenmd:::._EA_exn_grl)

#unlist and sort
AENMD_base_exonset_unlist_srtd <- sort(GenomeInfoDb::sortSeqlevels(unlist(AENMD_base_exonset)))

#find overlap between clinvar variants and the exons of the transcript set:
ov <- GenomicRanges::findOverlaps(clinvar_gr_total, AENMD_base_exonset_unlist_srtd)

#find the number of unique hits, such that a variant can only be counted once:
Variants_coding_our_tx_set <- length(unique(S4Vectors::queryHits(ov)))

#make a data.frame:
Total_var_stats <- data.frame(functional_class = c("Variants_in_the_original_dataset" , "Variants_coding_our_tx_set"),
                              value = c(Variants_in_the_original_dataset, Variants_coding_our_tx_set), percent = c(NA, Variants_coding_our_tx_set/Variants_in_the_original_dataset*100)
)
#_____________________________________________________________________________
#_____________________________________________________________________________

#_____________________________________________________________________________
###- 2. Next, we are going to run AENMD on the Clinvar data
#_____________________________________________________________________________
#get NA indexes and indexes where the alt is not a standard base (A, G, T, C)
NA_ind <- which(clinvarClinsig$ALT != "." & !stringr::str_detect(clinvarClinsig$ALT, "R|Y|S|W|K|M|B|D|H|V|N"))

#now put together the GRanges without NA to put into AENMD
#we do it this way to include the clinsig annotation
clinvar_gr  <- GenomicRanges::GRanges(clinvarClinsig$CHR[NA_ind],
                                      IRanges::IRanges(clinvarClinsig$POS[NA_ind],
                                                       clinvarClinsig$POS[NA_ind]))
clinvar_gr$ref     <- clinvarClinsig$REF[NA_ind] |> Biostrings::DNAStringSet()
clinvar_gr$alt     <- clinvarClinsig$ALT[NA_ind] |> Biostrings::DNAStringSet()
clinvar_gr$id      <- clinvarClinsig$ID[NA_ind]
clinvar_gr$clinsig      <- clinvarClinsig$CLNSIG[NA_ind]

#Process variants
clinvar_proc <- process_variants(clinvar_gr)

##- Annotate variants
#Our default setting annotation
res_clinvar_default_1 <- aenmd::annotate_nmd(clinvar_proc[1:floor(length(clinvar_proc)/3)])
saveRDS(res_clinvar_default_1, file = "clinvar0819_1.rds")
res_clinvar_default_2 <- aenmd::annotate_nmd(clinvar_proc[63240:126478])
saveRDS(res_clinvar_default_2, file = "clinvar0819_2.rds")
res_clinvar_default_3 <- aenmd::annotate_nmd(clinvar_proc[(2*floor(length(clinvar_proc)/3)+1):length(clinvar_proc)])
saveRDS(res_clinvar_default_3, file = "clinvar0819_3.rds")
