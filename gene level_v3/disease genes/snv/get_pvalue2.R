get_pvalue2 <- function(
    rds_name,
    outfilename,
    fdr_threshold = 0.05
) {
  
  res_p <- readRDS(rds_name)
  
  txnames <- unique(
    res_p@elementMetadata@listData[["res_aenmd"]]@
      listData[["transcript"]]
  )
  
  txnames.list <- vector("list", length(txnames))
  
  for (i in seq_along(txnames)) {
    
    print(i)
    
    txname <- txnames[i]
    
    txnames.list[[i]] <- list()
    txnames.list[[i]]$txname <- txname
    
    # ------------------------------------------------------------
    # Transcript and gene information
    # ------------------------------------------------------------
    
    BM.ind <- which(
      BM.info$ensembl_transcript_id == txname
    )
    
    if (length(BM.ind) > 0) {
      
      txnames.list[[i]]$hgnc_symbol <-
        BM.info[BM.ind[1], "hgnc_symbol"]
      
      txnames.list[[i]]$is_canonical <-
        BM.info[BM.ind[1], "transcript_is_canonical"]
      
    } else {
      
      txnames.list[[i]]$hgnc_symbol <- NA_character_
      txnames.list[[i]]$is_canonical <- NA
    }
    
    # ------------------------------------------------------------
    # Select variants belonging to the transcript
    # ------------------------------------------------------------
    
    sub.ind <- which(
      res_p@elementMetadata@listData[["res_aenmd"]]@
        listData[["transcript"]] == txname
    )
    
    res_p.ind <- res_p[sub.ind, ]
    
    # Remove non-PTC variants and variants from single-exon genes
    ptc.gen <- which(
      res_p.ind$res_aenmd$is_ptc == FALSE |
        res_p.ind$res_aenmd$is_single == TRUE
    )
    
    if (length(ptc.gen) > 0) {
      res_p.ind <- res_p.ind[-ptc.gen, ]
    }
    
    # ------------------------------------------------------------
    # Count PTCs in the canonical NMDesc region
    # ------------------------------------------------------------
    
    can.PTC.ind <- which(
      res_p.ind$res_aenmd$is_penultimate == TRUE |
        res_p.ind$res_aenmd$is_last == TRUE
    )
    
    if (length(can.PTC.ind) > 0) {
      
      can.PTC <- length(
        unique(
          as.character(
            res_p.ind$key[can.PTC.ind]
          )
        )
      )
      
    } else {
      
      can.PTC <- 0
    }
    
    # Total unique PTC variants in the transcript
    all.PTC <- length(unique(as.character(res_p.ind$key)))
    
    txnames.list[[i]]$can.PTC <- can.PTC
    txnames.list[[i]]$all.PTC <- all.PTC
    
    # ------------------------------------------------------------
    # Obtain transcript CDS coordinates
    # ------------------------------------------------------------
    
    ens.ind <- which(
      all.df.sub$TXNAME == txname
    )
    
    tx.cor <- all.df.sub[ens.ind, ]
    
    # Remove noncoding exons
    tx.cor <- tx.cor[
      !is.na(tx.cor$CDSSTART),
    ]
    
    # Only analyze multi-exon transcripts
    if (nrow(tx.cor) > 1) {
      
      txnames.list[[i]]$txcor <- tx.cor
      
      # ----------------------------------------------------------
      # Calculate CDS and canonical NMDesc-region length
      # ----------------------------------------------------------
      
      ex.length <- tx.cor$CDSEND - tx.cor$CDSSTART
      
      penultimate.length <- ex.length[
        length(ex.length) - 1
      ]
      
      last.exon.length <- ex.length[
        length(ex.length)
      ]
      
      if (penultimate.length < 50) {
        
        NMD.lastexon <-
          last.exon.length +
          penultimate.length
        
      } else {
        
        NMD.lastexon <-
          last.exon.length + 50
      }
      
      cds.start <- 1
      cds.end <- sum(ex.length)
      
      NMDesc.start <- cds.end - NMD.lastexon
      NMDesc.end <- cds.end
      
      chrom <- as.numeric(
        tx.cor$CDSCHROM[1]
      )
      
      # ----------------------------------------------------------
      # Count synonymous variants
      # ----------------------------------------------------------
      
      n_can <- get_syn_count(
        chrom,
        NMDesc.start,
        NMDesc.end
      )
      
      n_all <- get_syn_count(
        chrom,
        cds.start,
        cds.end
      )
      
      txnames.list[[i]]$chrom <- chrom
      txnames.list[[i]]$cds.start <- cds.start
      txnames.list[[i]]$cds.end <- cds.end
      txnames.list[[i]]$NMDesc.start <- NMDesc.start
      txnames.list[[i]]$NMDesc.end <- NMDesc.end
      txnames.list[[i]]$n_can <- n_can
      txnames.list[[i]]$n_all <- n_all
      
      # ----------------------------------------------------------
      # Variant densities
      # ----------------------------------------------------------
      
      if (!is.na(n_can) && n_can > 0) {
        
        can.PTC.density <- can.PTC / n_can
        
      } else {
        
        can.PTC.density <- NA_real_
      }
      
      if (!is.na(n_all) && n_all > 0) {
        
        all.PTC.density <- all.PTC / n_all
        
      } else {
        
        all.PTC.density <- NA_real_
      }
      
      # ----------------------------------------------------------
      # Exact binomial enrichment test
      #
      # x = PTCs observed in the NMDesc region
      # n = total PTCs in the transcript
      # p = expected NMDesc proportion based on synonymous variants
      # ----------------------------------------------------------
      
      if (
        is.na(n_can) ||
        is.na(n_all) ||
        n_all <= 0 ||
        n_can < 0 ||
        n_can > n_all ||
        all.PTC <= 0 ||
        can.PTC < 0 ||
        can.PTC > all.PTC
      ) {
        
        expected.proportion <- NA_real_
        expected.PTC <- NA_real_
        fold.enrichment <- NA_real_
        can.pvalue <- NA_real_
        
      } else {
        
        expected.proportion <- n_can / n_all
        expected.PTC <- all.PTC * expected.proportion
        
        if (expected.PTC > 0) {
          
          fold.enrichment <- can.PTC / expected.PTC
          
        } else {
          
          fold.enrichment <- NA_real_
        }
        
        # Handle boundary probabilities explicitly
        if (expected.proportion == 0) {
          
          can.pvalue <- ifelse(
            can.PTC > 0,
            0,
            1
          )
          
        } else if (expected.proportion == 1) {
          
          can.pvalue <- 1
          
        } else {
          
          can.pvalue <- binom.test(
            x = can.PTC,
            n = all.PTC,
            p = expected.proportion,
            alternative = "greater"
          )$p.value
        }
      }
      
      # ----------------------------------------------------------
      # Save transcript-level results
      # ----------------------------------------------------------
      
      txnames.list[[i]]$can.PTC.density <-
        can.PTC.density
      
      txnames.list[[i]]$all.PTC.density <-
        all.PTC.density
      
      txnames.list[[i]]$expected.proportion <-
        expected.proportion
      
      txnames.list[[i]]$expected.PTC <-
        expected.PTC
      
      txnames.list[[i]]$fold.enrichment <-
        fold.enrichment
      
      txnames.list[[i]]$can.pvalue <-
        can.pvalue
      
    } else {
      
      # Transcript is missing or has only one coding exon
      txnames.list[[i]]$can.pvalue <- NA_real_
    }
  }
  
  # ============================================================
  # Benjamini-Hochberg FDR correction
  # ============================================================
  
  raw.pvalues <- vapply(
    txnames.list,
    function(x) {
      
      if (
        is.null(x$can.pvalue) ||
        length(x$can.pvalue) == 0
      ) {
        
        return(NA_real_)
        
      } else {
        
        return(as.numeric(x$can.pvalue[1]))
      }
    },
    numeric(1)
  )
  
  valid.tests <- which(
    !is.na(raw.pvalues) &
      is.finite(raw.pvalues)
  )
  
  fdr.values <- rep(
    NA_real_,
    length(raw.pvalues)
  )
  
  if (length(valid.tests) > 0) {
    
    fdr.values[valid.tests] <- p.adjust(
      raw.pvalues[valid.tests],
      method = "BH"
    )
  }
  
  # Add FDR results to each transcript
  for (i in seq_along(txnames.list)) {
    
    txnames.list[[i]]$can.fdr <-
      fdr.values[i]
    
    txnames.list[[i]]$can.fdr.significant <-
      if (
        is.na(fdr.values[i])
      ) {
        NA
      } else {
        fdr.values[i] < fdr_threshold
      }
  }
  
  # Store information about the multiple-testing procedure
  attr(txnames.list, "multiple_testing_method") <-
    "Benjamini-Hochberg"
  
  attr(txnames.list, "fdr_threshold") <-
    fdr_threshold
  
  attr(txnames.list, "number_of_tests") <-
    length(valid.tests)
  
  # Save results
  saveRDS(
    txnames.list,
    outfilename
  )
  
  invisible(txnames.list)
}