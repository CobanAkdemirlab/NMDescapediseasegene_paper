get_pvalue_wald = function(rds_name, outfilename){
  res_p = readRDS(rds_name)
  txnames <- unique(res_p@elementMetadata@listData[["res_aenmd"]]@listData[["transcript"]])
  txnames.list <- list()
  get_syn_count = function(chrom, region.start, region.end) {
    syn.sub = syn_all[
      syn_all$CHROM == chrom &
        syn_all$cds_pos >= region.start &
        syn_all$cds_pos < region.end,
    ]
    syn.count = nrow(syn.sub)
    return(syn.count)
  }
  syn_all = read.csv("/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/clinvar/syn_all_can_only.csv")
  BM.info = read.csv("/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/clinvar/BM_info.csv")
  all.df.sub = read.csv("/Users/jxu14/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data/clinvar/all_df_sub.csv")
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
      can.PTC <- length(unique(sapply(1:length(names(res_p.ind[can.PTC.ind])), function(x){
        res_p.ind$key[can.PTC.ind[x]]
      })))
    } else {
      can.PTC <- 0
    }
    
    ens.ind <- which(all.df.sub$TXNAME == txname)
    tx.cor <- all.df.sub[ens.ind,]
    tx.cor <- tx.cor[!is.na(tx.cor$CDSSTART),]
    
    if (nrow(tx.cor) > 1) {
      txnames.list[[i]]$txcor <- tx.cor
      txnames.list[[i]]$can.PTC <- can.PTC
      
      ex.length <- tx.cor$CDSEND - tx.cor$CDSSTART
      if(ex.length[length(ex.length)-1] < 50){
        NMD.lastexon <- ex.length[length(ex.length)] + ex.length[length(ex.length)-1]
      } else {
        NMD.lastexon <- ex.length[length(ex.length)] + 50
      }
      
      cds.start    = 1
      cds.end      = sum(ex.length)
      NMDesc.start = cds.end - NMD.lastexon
      NMDesc.end   = cds.end
      chrom        = as.numeric(tx.cor$CDSCHROM[1])
      
      can.PTC.density <- can.PTC / get_syn_count(chrom, NMDesc.start, NMDesc.end)
      all.gr <- GRanges(seqnames = c(tx.cor$CDSCHROM),
                        ranges   = IRanges(tx.cor$CDSSTART, tx.cor$CDSEND),
                        strand   = tx.cor$CDSSTRAND[1])
      
      all.PTC         <- length(unique(as.character(res_p.ind$key)))
      all.PTC.density <- all.PTC / get_syn_count(chrom, cds.start, cds.end)
      
      n_can <- get_syn_count(chrom, NMDesc.start, NMDesc.end)
      n_all <- get_syn_count(chrom, cds.start, cds.end)
      p0    <- all.PTC / n_all
      
      # ── 1. 原有：二项检验 ─────────────────────────────────────────────────
      if (is.na(n_can) || n_can <= 0 || can.PTC > n_can ||
          is.na(n_all) || n_all <= 0 || is.na(p0) || p0 < 0 || p0 > 1) {
        can.pvalue <- NA_real_
      } else {
        can.pvalue <- binom.test(can.PTC, n_can, p0, alternative = "greater")$p.value
      }
      
      # ── 2. 新增：Wald log OR 检验 ─────────────────────────────────────────
      # 构建 2x2 列联表：
      #                PTC数量        同义位点数（分母）
      # NMD逃逸区域    can.PTC        n_can
      # CDS其余区域    rest.PTC       n_rest
      
      rest.PTC <- all.PTC - can.PTC    # CDS其余区域的PTC数
      n_rest   <- n_all - n_can        # CDS其余区域的同义位点数
      
      if (is.na(n_can) || is.na(n_rest) ||
          n_can <= 0   || n_rest <= 0   ||
          can.PTC > n_can || rest.PTC > n_rest) {
        
        can.wald.logOR  <- NA_real_
        can.wald.se     <- NA_real_
        can.wald.pvalue <- NA_real_
        can.wald.CI.low <- NA_real_
        can.wald.CI.up  <- NA_real_
        
      } else {
        
        # 加0.5连续性校正（Haldane-Anscombe），避免log(0)
        a <- can.PTC + 0.5              # NMD逃逸区 PTC
        b <- n_can - can.PTC + 0.5     # NMD逃逸区 非PTC
        c <- rest.PTC + 0.5            # CDS其余区 PTC
        d <- n_rest - rest.PTC + 0.5   # CDS其余区 非PTC
        
        logOR    <- log(a * d / (b * c))
        se_logOR <- sqrt(1/a + 1/b + 1/c + 1/d)  # Woolf公式
        z_stat   <- logOR / se_logOR
        
        can.wald.logOR  <- logOR
        can.wald.se     <- se_logOR
        can.wald.pvalue <- pnorm(z_stat, lower.tail = FALSE)  # 单侧p值：enrichment
        can.wald.CI.low <- logOR - 1.96 * se_logOR      # 95%置信区间下限
        can.wald.CI.up  <- logOR + 1.96 * se_logOR      # 95%置信区间上限
      }
      
      # ── 保存所有输出 ──────────────────────────────────────────────────────
      txnames.list[[i]]$can.PTC.density  <- can.PTC.density
      txnames.list[[i]]$all.PTC.density  <- all.PTC.density
      txnames.list[[i]]$can.pvalue       <- can.pvalue       # 二项检验p值
      txnames.list[[i]]$can.wald.logOR   <- can.wald.logOR   # log OR值
      txnames.list[[i]]$can.wald.se      <- can.wald.se      # log OR的标准误
      txnames.list[[i]]$can.wald.pvalue  <- can.wald.pvalue  # Wald检验p值
      txnames.list[[i]]$can.wald.CI.low  <- can.wald.CI.low  # 95%置信区间下限
      txnames.list[[i]]$can.wald.CI.up   <- can.wald.CI.up   # 95%置信区间上限
    }
  }
  saveRDS(txnames.list, outfilename)
}