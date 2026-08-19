# =============================================================================
# get_pvalue (v3): restrict the test universe UPFRONT to canonical transcripts
# of a pre-specified gene set (e.g. OMIM autosomal-dominant genes).
#
# The restriction happens right after txnames is extracted, so:
#   * the loop only runs over canonical AD transcripts (faster),
#   * every downstream correction (BH / binom-BH / DBH / DBY) is automatically
#     computed within this universe — no post-hoc subsetting needed,
#   * the discrete DBH/DBY supports are, by construction, built from exactly
#     the corrected test set.
#
# `restrict_symbols = NULL` reproduces the old genome-wide behavior.
#
# usage:
#   omim_AD_symbols <- read.csv('omim_AD_symbols.csv', header = TRUE)$hgnc_symbol
#   res <- get_pvalue('snv_plp_ptc_can_filtered20260201.rds',
#                     'snv_plp_ptc_nmdesc_can_filtered20260201.rds',
#                     'snv_plp_ptc_nmdesc_can_p_f_syn_20260201_ADrestricted.rds',
#                     restrict_symbols = omim_AD_symbols)
# =============================================================================

get_pvalue3 = function(rds_name, rds_name2, outfilename,
                      restrict_symbols = NULL, alpha = 0.05){
  res_p = readRDS(rds_name)
  res_p2 = readRDS(rds_name2)
  txnames <- unique(res_p2@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]]) #only circle the variants with NMDesc PTC
  
  # ---- NEW: restrict txnames to canonical transcripts of restrict_symbols ----
  if (!is.null(restrict_symbols)) {
    restrict_symbols <- unique(trimws(as.character(restrict_symbols)))
    restrict_symbols <- restrict_symbols[!is.na(restrict_symbols) &
                                           nchar(restrict_symbols) > 0]
    
    keep.tx <- BM.info$ensembl_transcript_id[
      BM.info$hgnc_symbol %in% restrict_symbols &
        !is.na(BM.info$transcript_is_canonical) &
        BM.info$transcript_is_canonical == 1
    ]
    
    n.before <- length(txnames)
    txnames  <- txnames[txnames %in% keep.tx]
    message(sprintf("Restriction: %d NMDesc transcripts -> %d canonical transcripts of %d supplied symbols",
                    n.before, length(txnames), length(restrict_symbols)))
    if (length(txnames) == 0)
      stop("No transcripts left after restriction — check that BM.info symbols match restrict_symbols.")
  }
  # ---------------------------------------------------------------------------
  
  txnames.list <- list()
  
  for (i in 1:length(txnames)){
    print(i)
    txname = txnames[i]
    txnames.list[[i]] <- list()
    txnames.list[[i]]$txname <- txname
    BM.ind <- which(BM.info$ensembl_transcript_id == txname)
    hgnc_symbol <- BM.info[BM.ind, 'hgnc_symbol']
    txnames.list[[i]]$hgnc_symbol <- hgnc_symbol
    txnames.list[[i]]$is_canonical <- BM.info[BM.ind, 'transcript_is_canonical']
    sub.ind <- which(res_p@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]] == txname)
    res_p.ind <- res_p[sub.ind,]
    ptc.gen <- which(res_p.ind$res_aenmd$is_ptc == FALSE | res_p.ind$res_aenmd$is_single == TRUE)
    if(length(ptc.gen) > 0){
      res_p.ind <- res_p.ind[-ptc.gen,]
    }
    can.PTC.ind <- which(res_p.ind$res_aenmd$is_penultimate == TRUE | res_p.ind$res_aenmd$is_last == TRUE)
    if(length(can.PTC.ind) > 0){
      can.PTC <- length(unique(sapply(1:length(names(res_p.ind[can.PTC.ind])), function(x)
      {
        res_p.ind$key[can.PTC.ind[x]]
      }))) } else {
        can.PTC <- 0
      }
    
    ens.ind <- which(all.df.sub$TXNAME == txname)
    tx.cor <- all.df.sub[ens.ind,]
    tx.cor <- tx.cor[!is.na(tx.cor$CDSSTART),]
    
    if (nrow(tx.cor) > 1) {
      txnames.list[[i]]$txcor <- tx.cor
      txnames.list[[i]]$can.PTC <- can.PTC
      
      # FIX 1: exon length is a closed interval, add 1
      ex.length <- tx.cor$CDSEND - tx.cor$CDSSTART + 1
      if(ex.length[length(ex.length)-1] < 50) {
        NMD.lastexon <- ex.length[length(ex.length)] + ex.length[length(ex.length)-1]
      } else {
        NMD.lastexon <- ex.length[length(ex.length)] + 50
      }
      
      cds.start = 1
      cds.end = sum(ex.length)
      NMDesc.start = cds.end - NMD.lastexon
      NMDesc.end = cds.end
      chrom = as.numeric(tx.cor$CDSCHROM[1])
      can.PTC.density <- can.PTC / get_syn_count(chrom, NMDesc.start, NMDesc.end)
      all.gr <- GRanges(seqnames = c(tx.cor$CDSCHROM),
                        ranges = IRanges(tx.cor$CDSSTART, tx.cor$CDSEND),
                        strand = tx.cor$CDSSTRAND[1])
      
      all.PTC <- length(unique(as.character(res_p.ind$key)))
      all.PTC.density <- all.PTC / get_syn_count(chrom, cds.start, cds.end)
      
      n_can <- get_syn_count(chrom, NMDesc.start, NMDesc.end)
      n_all <- get_syn_count(chrom, cds.start, cds.end)
      
      # FIX 2: baseline must EXCLUDE the region being tested, otherwise
      # can.PTC inflates its own null (conservative, and length-dependent)
      rest.PTC <- all.PTC - can.PTC
      n_rest   <- n_all - n_can
      p0 <- rest.PTC / n_rest
      
      if (is.na(n_can) || n_can <= 0 || can.PTC > n_can ||
          is.na(n_rest) || n_rest <= 0 || rest.PTC > n_rest ||
          is.na(p0) || p0 < 0 || p0 > 1) {
        can.pvalue        <- NA_real_
        can.fisher.pvalue <- NA_real_
      } else {
        can.pvalue <- binom.test(can.PTC, n_can, p0, alternative = "greater")$p.value
        
        # FIX 3: Fisher's exact test — valid with zero cells and small counts,
        # and treats both counts as random rather than fixing p0
        can.fisher.pvalue <- fisher.test(
          matrix(c(can.PTC, n_can - can.PTC,
                   rest.PTC, n_rest - rest.PTC), nrow = 2),
          alternative = "greater")$p.value
      }
      
      txnames.list[[i]]$can.PTC.density   <- can.PTC.density
      txnames.list[[i]]$all.PTC.density   <- all.PTC.density
      txnames.list[[i]]$can.pvalue        <- can.pvalue
      txnames.list[[i]]$can.fisher.pvalue <- can.fisher.pvalue
      txnames.list[[i]]$all.PTC  <- all.PTC
      txnames.list[[i]]$rest.PTC <- rest.PTC
      txnames.list[[i]]$n_can    <- n_can
      txnames.list[[i]]$n_rest   <- n_rest
    }
  }
  
  # ---------- multiple testing correction ----------
  # (universe is already restricted, so all tests below are within it)
  extract <- function(field) {
    sapply(txnames.list, function(x) {
      v <- x[[field]]
      if (is.null(v) || length(v) == 0) NA_real_ else as.numeric(v)
    })
  }
  
  p.fisher <- extract("can.fisher.pvalue")
  p.binom  <- extract("can.pvalue")
  ok <- !is.na(p.fisher)
  n.test <- sum(ok)
  
  # standard BH — valid but conservative under discreteness; kept as reference
  fdr.bh <- rep(NA_real_, length(p.fisher))
  fdr.bh[ok] <- p.adjust(p.fisher[ok], method = "BH")
  
  fdr.binom <- rep(NA_real_, length(p.binom))
  ok.b <- !is.na(p.binom)
  fdr.binom[ok.b] <- p.adjust(p.binom[ok.b], method = "BH")
  
  # discrete FDR using the test-specific p-value supports
  sig.dbh  <- rep(FALSE,    length(p.fisher))
  sig.dby  <- rep(FALSE,    length(p.fisher))
  padj.dbh <- rep(NA_real_, length(p.fisher))
  padj.dby <- rep(NA_real_, length(p.fisher))
  min.p    <- rep(NA_real_, length(p.fisher))
  
  if (requireNamespace("DiscreteFDR", quietly = TRUE) && n.test > 0) {
    counts <- t(sapply(which(ok), function(i) {
      x <- txnames.list[[i]]
      c(x$can.PTC, x$n_can - x$can.PTC,
        x$rest.PTC, x$n_rest - x$rest.PTC)
    }))
    
    fp <- DiscreteFDR::fisher.pvalues.support(counts, alternative = "greater")
    
    # minimum attainable p-value = smallest value in each test's support
    min.p[ok] <- sapply(fp$support, min)
    
    # ret.crit.consts = FALSE so DiscreteFDR >= 2.0 returns $Adjusted
    res.dbh <- DiscreteFDR::DBH(fp$raw, fp$support, alpha = alpha,
                                ret.crit.consts = FALSE)
    res.dby <- DiscreteFDR::DBY(fp$raw, fp$support, alpha = alpha,
                                ret.crit.consts = FALSE)
    
    sig.dbh[which(ok)[res.dbh$Indices]] <- TRUE
    sig.dby[which(ok)[res.dby$Indices]] <- TRUE
    
    # adjusted p-values, in the original input order of fp$raw
    grab_adj <- function(res) {
      for (slot in c("Adjusted", "adjusted", "Adjusted.pvalues")) {
        if (!is.null(res[[slot]])) return(as.numeric(res[[slot]]))
      }
      NULL
    }
    
    adj.dbh <- grab_adj(res.dbh)
    adj.dby <- grab_adj(res.dby)
    
    if (!is.null(adj.dbh) && length(adj.dbh) == n.test) {
      padj.dbh[ok] <- adj.dbh
      if (sum(padj.dbh < alpha, na.rm = TRUE) != sum(sig.dbh))
        warning("DBH: adjusted p-values and $Indices disagree — check package version / ordering.")
    } else {
      warning("DBH: no $Adjusted slot returned; update DiscreteFDR ",
              "(install.packages('DiscreteFDR')).")
    }
    if (!is.null(adj.dby) && length(adj.dby) == n.test) {
      padj.dby[ok] <- adj.dby
      if (sum(padj.dby < alpha, na.rm = TRUE) != sum(sig.dby))
        warning("DBY: adjusted p-values and $Indices disagree — check package version / ordering.")
    } else {
      warning("DBY: no $Adjusted slot returned; update DiscreteFDR.")
    }
    
  } else if (n.test > 0) {
    warning("DiscreteFDR not installed; reporting standard BH only. install.packages('DiscreteFDR')")
    for (i in which(ok)) {
      x <- txnames.list[[i]]
      tot.PTC <- x$can.PTC + x$rest.PTC
      k <- min(tot.PTC, x$n_can)
      min.p[i] <- fisher.test(
        matrix(c(k, x$n_can - k,
                 tot.PTC - k, x$n_rest - (tot.PTC - k)), nrow = 2),
        alternative = "greater")$p.value
    }
  }
  
  for (i in seq_along(txnames.list)) {
    if (length(txnames.list[[i]]) == 0) next
    txnames.list[[i]]$sig.dbh          <- sig.dbh[i]    # PRIMARY: discrete BH, FDR < alpha
    txnames.list[[i]]$sig.dby          <- sig.dby[i]    # discrete BY, any dependence
    txnames.list[[i]]$padj.dbh         <- padj.dbh[i]   # DBH adjusted p-value
    txnames.list[[i]]$padj.dby         <- padj.dby[i]   # DBY adjusted p-value
    txnames.list[[i]]$fdr.bh           <- fdr.bh[i]     # standard BH, reference
    txnames.list[[i]]$fdr.binom        <- fdr.binom[i]  # binomial test, comparison
    txnames.list[[i]]$min.attainable.p <- min.p[i]
  }
  
  message(sprintf("Tests performed: %d", n.test))
  message(sprintf("Cannot reach p < 0.05 at any effect size: %d (%.1f%%)",
                  sum(min.p > 0.05, na.rm = TRUE),
                  100 * mean(min.p[ok] > 0.05, na.rm = TRUE)))
  message(sprintf("--- Significant at FDR < %.2f ---", alpha))
  message(sprintf("  DBH  (discrete)            : %d", sum(sig.dbh)))
  message(sprintf("  DBY  (discrete, any dep.)  : %d", sum(sig.dby)))
  message(sprintf("  BH   (standard, reference) : %d", sum(fdr.bh    < alpha, na.rm = TRUE)))
  message(sprintf("  BH   (binomial test)       : %d", sum(fdr.binom < alpha, na.rm = TRUE)))
  
  saveRDS(txnames.list, outfilename)
  invisible(txnames.list)
}