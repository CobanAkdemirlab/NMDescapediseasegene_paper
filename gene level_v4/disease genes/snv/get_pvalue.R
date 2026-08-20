# --- 路径解析层（自动插入）---------------------------------
# 数据文件用 data_file("文件名") 定位；输出用 out_file("文件名")
# 换数据位置只改 gene level_v3/lib/paths.R 的 DATA_ROOTS
.p <- c("gene level_v3/lib/paths.R", "../lib/paths.R", "../../lib/paths.R",
        "../../../gene level_v3/lib/paths.R", "lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p)) stop("找不到 paths.R —— 请从仓库根目录运行 R")
source(.p[1])
# ------------------------------------------------------------

# =============================================================================
# get_pvalue (v5): AD-restricted universe; PRIMARY = Fisher exact + standard
# Benjamini-Hochberg, two-tier reporting at FDR 0.05 / 0.20.
#
#   Tier 1  significant : fdr.bh < 0.05
#   Tier 2  suggestive  : 0.05 <= fdr.bh < 0.20
#   Main list           : fdr.bh < 0.20  (union of the two tiers)
#
# Rationale for the 0.20 suggestive threshold (k/0 mechanism; for Methods):
#   Dominant NMD-escape genes deposit essentially all P/LP PTCs inside the
#   escape region (escape -> stable protein -> pathogenic -> in ClinVar;
#   elsewhere -> degraded -> absent), so their characteristic pattern is
#   can.PTC = k, rest.PTC = 0 ("k/0"). Fisher's test conditions on the PTC
#   total, so a k/0 gene's minimum attainable p is ~ f^k (f = escape-region
#   share of synonymous sites, ~0.1-0.15): with k <= 2 such genes can NEVER
#   pass FDR 0.05 at ~300 tests, regardless of effect size — i.e. FDR 0.05
#   structurally excludes the archetype most characteristic of the mechanism.
#   FDR 0.20 makes k = 2 reachable while Fisher retains exact calibration
#   (a binomial alternative degenerates at k/0 since p0 = 0). Empirically,
#   the FDR < 0.20 list was enriched against all three reference sets after
#   background correction (tested-universe hypergeometric: extensive 1.25x
#   p = 0.018; curated NMD-escape 1.59x p = 0.023; NDD 1.43x p = 0.008),
#   with the 0.05-0.20 band matching the significant tier's precision.
#
# Why standard BH (not discrete DBH) for the primary tiers:
#   At alpha = 0.20 the discrete correction admitted an additional tail
#   (standard-BH-FDR 0.20-0.35) whose reference-set precision was at
#   background level; standard BH places the 0.20 boundary at the empirical
#   signal edge. Discrete DBH/DBY are retained as comparison tracks; BH's
#   conservatism under discreteness makes the primary lists strictly valid.
#
# Comparison tracks kept per transcript: DBH (step-up flag + step-down
# adjusted p), DBY (flag + adjusted p), binomial BH.
#
# usage:
#   omim_AD_symbols <- read.csv('omim_AD_symbols.csv', header = TRUE)$hgnc_symbol
#   res <- get_pvalue('snv_plp_ptc_can_filtered20260201.rds',
#                     'snv_plp_ptc_nmdesc_can_filtered20260201.rds',
#                     'snv_plp_ptc_nmdesc_can_p_f_syn_20260201_AD_BH020.rds',
#                     restrict_symbols = omim_AD_symbols)
#   tab <- export_tiers(res)   # writes tier gene lists + annotated CSV
# =============================================================================

get_pvalue = function(rds_name, rds_name2, outfilename,
                      restrict_symbols = NULL,
                      alpha_main = 0.20, alpha_tier = 0.05){
  res_p = readRDS(rds_name)
  res_p2 = readRDS(rds_name2)
  txnames <- unique(res_p2@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]]) #only circle the variants with NMDesc PTC
  
  # ---- restrict txnames to canonical transcripts of restrict_symbols ----
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
  # -----------------------------------------------------------------------
  
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
      
      # ---- k/0 archetype annotation (evidence trail for the 0.20 tier) --
      txnames.list[[i]]$f_region  <- n_can / n_all
      txnames.list[[i]]$is_k0     <- (rest.PTC == 0 && can.PTC > 0)
      txnames.list[[i]]$archetype <- sprintf("%d/%d", can.PTC, rest.PTC)
    }
  }
  
  # ---------- multiple testing correction ----------
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
  
  # ---- PRIMARY: standard BH on Fisher p, two tiers -----------------------
  fdr.bh <- rep(NA_real_, length(p.fisher))
  fdr.bh[ok] <- p.adjust(p.fisher[ok], method = "BH")
  
  tier <- rep(NA_character_, length(p.fisher))
  tier[ok] <- ifelse(fdr.bh[ok] < alpha_tier, "significant",
                     ifelse(fdr.bh[ok] < alpha_main, "suggestive", "ns"))
  
  # binomial BH (comparison)
  fdr.binom <- rep(NA_real_, length(p.binom))
  ok.b <- !is.na(p.binom)
  fdr.binom[ok.b] <- p.adjust(p.binom[ok.b], method = "BH")
  
  # ---- comparison tracks: discrete FDR (DBH / DBY) -----------------------
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
    
    # DBH: step-up flag at alpha_tier; step-down for adjusted p-values
    res.su <- DiscreteFDR::DBH(fp$raw, fp$support, alpha = alpha_tier,
                               ret.crit.consts = FALSE)
    sig.dbh[which(ok)[res.su$Indices]] <- TRUE
    
    res.sd <- DiscreteFDR::DBH(fp$raw, fp$support, alpha = alpha_tier,
                               direction = "sd", ret.crit.consts = FALSE)
    if (!is.null(res.sd$Adjusted) && length(res.sd$Adjusted) == n.test)
      padj.dbh[ok] <- as.numeric(res.sd$Adjusted)
    
    res.dby <- DiscreteFDR::DBY(fp$raw, fp$support, alpha = alpha_tier,
                                ret.crit.consts = FALSE)
    sig.dby[which(ok)[res.dby$Indices]] <- TRUE
    if (!is.null(res.dby$Adjusted) && length(res.dby$Adjusted) == n.test)
      padj.dby[ok] <- as.numeric(res.dby$Adjusted)
    
  } else if (n.test > 0) {
    warning("DiscreteFDR not installed; discrete comparison tracks skipped. install.packages('DiscreteFDR')")
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
  
  # structural reachability at alpha_tier (k/0 evidence trail):
  # conservative BH effective threshold ~ alpha_tier * n_sig / n_test
  n.sig05 <- sum(tier == "significant", na.rm = TRUE)
  thr05   <- alpha_tier * max(n.sig05, 1) / max(n.test, 1)
  
  for (i in seq_along(txnames.list)) {
    if (length(txnames.list[[i]]) == 0) next
    txnames.list[[i]]$tier             <- tier[i]         # PRIMARY: significant / suggestive / ns
    txnames.list[[i]]$fdr.bh           <- fdr.bh[i]       # PRIMARY: standard BH (tiering basis)
    txnames.list[[i]]$sig.dbh          <- sig.dbh[i]      # DBH flag at alpha_tier (comparison)
    txnames.list[[i]]$padj.dbh         <- padj.dbh[i]     # DBH step-down adjusted p
    txnames.list[[i]]$sig.dby          <- sig.dby[i]      # DBY flag
    txnames.list[[i]]$padj.dby         <- padj.dby[i]     # DBY adjusted p
    txnames.list[[i]]$fdr.binom        <- fdr.binom[i]    # binomial BH (comparison)
    txnames.list[[i]]$min.attainable.p <- min.p[i]
    txnames.list[[i]]$unreachable.at.005 <-
      !is.na(min.p[i]) && min.p[i] > thr05
  }
  
  # tier composition summary
  is.k0.v <- sapply(txnames.list, function(x)
    if (length(x) == 0 || is.null(x$is_k0)) NA else x$is_k0)
  unr.v   <- sapply(txnames.list, function(x) isTRUE(x$unreachable.at.005))
  n.sug     <- sum(tier == "suggestive", na.rm = TRUE)
  n.sug.k0  <- sum(tier == "suggestive" & is.k0.v, na.rm = TRUE)
  n.sug.unr <- sum(tier == "suggestive" & unr.v,   na.rm = TRUE)
  
  message(sprintf("Tests performed: %d", n.test))
  message(sprintf("Cannot reach p < 0.05 at any effect size: %d (%.1f%%)",
                  sum(min.p > 0.05, na.rm = TRUE),
                  100 * mean(min.p[ok] > 0.05, na.rm = TRUE)))
  message("--- PRIMARY: Fisher + standard BH, two tiers ---")
  message(sprintf("  Tier 1 significant (FDR < %.2f)         : %d", alpha_tier, n.sig05))
  message(sprintf("  Tier 2 suggestive  (%.2f <= FDR < %.2f)  : %d", alpha_tier, alpha_main, n.sug))
  message(sprintf("    of which k/0 archetype                : %d", n.sug.k0))
  message(sprintf("    structurally unreachable at %.2f      : %d", alpha_tier, n.sug.unr))
  message(sprintf("  Main list (FDR < %.2f) total            : %d", alpha_main, n.sig05 + n.sug))
  message(sprintf("--- Comparison tracks at FDR < %.2f ---", alpha_tier))
  message(sprintf("  DBH  (discrete)            : %d", sum(sig.dbh)))
  message(sprintf("  DBY  (discrete, any dep.)  : %d", sum(sig.dby)))
  message(sprintf("  BH   (binomial test)       : %d", sum(fdr.binom < alpha_tier, na.rm = TRUE)))
  
  saveRDS(txnames.list, outfilename)
  invisible(txnames.list)
}

# =============================================================================
# Post-processing: tiered gene lists + full annotated table
# =============================================================================
export_tiers <- function(txnames.list, prefix = "snv_can_ADrestricted_bh") {
  tab <- do.call(rbind, lapply(txnames.list, function(x) {
    if (is.null(x$can.fisher.pvalue) || is.na(x$can.fisher.pvalue)) return(NULL)
    data.frame(
      txname = x$txname,
      hgnc_symbol = paste(unique(x$hgnc_symbol), collapse = ";"),
      archetype = x$archetype, is_k0 = x$is_k0,
      can.PTC = x$can.PTC, rest.PTC = x$rest.PTC,
      fisher_p = x$can.fisher.pvalue,
      fdr.bh = x$fdr.bh,
      padj.dbh = x$padj.dbh, padj.dby = x$padj.dby,
      fdr.binom = x$fdr.binom,
      min.attainable.p = x$min.attainable.p,
      unreachable.at.005 = x$unreachable.at.005,
      tier = x$tier, stringsAsFactors = FALSE)
  }))
  tab <- tab[order(tab$fdr.bh, tab$fisher_p), ]
  tab$rank <- seq_len(nrow(tab))
  write.csv(tab, paste0(prefix, "_full_annotated.csv"), row.names = FALSE)
  for (tr in c("significant", "suggestive")) {
    g <- unique(tab$hgnc_symbol[tab$tier == tr])
    writeLines(c("hgnc_symbol", g), sprintf("%s_%s.txt", prefix, tr))
    cat(sprintf("%-11s: %3d genes\n", tr, length(g)))
  }
  g.all <- unique(tab$hgnc_symbol[tab$tier %in% c("significant","suggestive")])
  writeLines(c("hgnc_symbol", g.all), sprintf("%s_FDR0.20_all.txt", prefix))
  cat(sprintf("%-11s: %3d genes\n", "combined", length(g.all)))
  invisible(tab)
}