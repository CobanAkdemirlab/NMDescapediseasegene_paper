These files are used to **get P/LP and benign variants** corresponding to the gene list from Clinvar and GnomAD, then **create fasta files** from PTC variants based on their key and transcript id.

get_XX_variant_new.R should be the first R script to run. Note the difference between **+1 and -1 strand directions**.

If you want to filter for **AD variants**, run clean_variant_AD.R.
