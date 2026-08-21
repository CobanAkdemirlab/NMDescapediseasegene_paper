get_NMD_enrichment_wald = function(rds_name, FDR = 0.05, filter_type){
  l <- readRDS(rds_name)

  
  # ── 辅助函数：从 p_indexes 提取基因名并写出 ──────────────────────────────
  extract_and_write <- function(p_indexes, l, out_name) {
    ind.1 <- which(unlist(p_indexes) == TRUE)
    if (length(ind.1) == 0) {
      cat('⚠️ 没有找到符合条件的基因，输出空文件:', out_name, '\n')
      writeLines(character(0), out_name)
      return(invisible(character(0)))
    }
    new.genes <- unlist(sapply(ind.1, function(x) {
      sym <- l[[x]]$hgnc_symbol
      if (is.null(sym) || length(sym) == 0) return(NULL)
      as.character(sym)
    }))
    new.genes <- new.genes[!is.na(new.genes) & nchar(new.genes) > 0]
    cat('✅ 找到基因数:', length(new.genes), '→', out_name, '\n')
    writeLines(new.genes, out_name)
    return(invisible(new.genes))
  }
  
  # ── 1. 基于 can.pvalue（二项检验）做 BH 校正 ─────────────────────────────
  can.pvalues <- sapply(1:length(l), function(p) as.numeric(l[[p]]$can.pvalue))
  can.pvalues.adj <- p.adjust(can.pvalues, method = "BH")
  for (p in 1:length(l)) l[[p]]$can.pvalue.adj <- can.pvalues.adj[p]
  cat('二项检验 FDR < cut 的转录本数:', sum(can.pvalues.adj < FDR, na.rm=TRUE), '\n')
  
  # ── 2. 基于 can.wald.pvalue（Wald log OR）做 BH 校正 ─────────────────────
  can.wald.pvalues <- sapply(1:length(l), function(p) as.numeric(l[[p]]$can.wald.pvalue))
  can.wald.pvalues.adj <- p.adjust(can.wald.pvalues, method = "BH")
  for (p in 1:length(l)) l[[p]]$can.wald.pvalue.adj <- can.wald.pvalues.adj[p]
  cat('Wald检验 FDR < cut 的转录本数:', sum(can.wald.pvalues.adj < FDR, na.rm=TRUE), '\n')
  
  # ── 3. 构建两组 p_indexes ─────────────────────────────────────────────────
  p_indexes_binom <- list()   # 二项检验
  p_indexes_wald  <- list()   # Wald log OR
  
  for (p in 1:length(l)) {
    is_can <- isTRUE(as.logical(l[[p]]$is_canonical))
    
    if(filter_type == 'can'){
      p_indexes_binom[[p]] <- isTRUE(as.numeric(l[[p]]$can.pvalue.adj)      < FDR) & is_can
      p_indexes_wald[[p]]  <- isTRUE(as.numeric(l[[p]]$can.wald.pvalue.adj) < FDR) & is_can
      
    } else if(filter_type == 'all'){
      p_indexes_binom[[p]] <- (as.numeric(l[[p]]$css.pvalue)  >= css_cut |
                                 as.numeric(l[[p]]$can.pvalue.adj) <= can_cut |
                                 as.numeric(l[[p]]$long.pvalue) >= long_cut) & is_can
      p_indexes_wald[[p]]  <- (as.numeric(l[[p]]$css.pvalue)  >= css_cut |
                                 as.numeric(l[[p]]$can.wald.pvalue.adj) <= can_wald_cut |
                                 as.numeric(l[[p]]$long.pvalue) >= long_cut) & is_can
      
    } else if(filter_type == 'css'){
      p_indexes_binom[[p]] <- (as.numeric(l[[p]]$css.pvalue) >= css_cut) & is_can
      p_indexes_wald[[p]]  <- p_indexes_binom[[p]]   # css 不涉及 wald，两组相同
      
    } else if(filter_type == 'long'){
      p_indexes_binom[[p]] <- (as.numeric(l[[p]]$long.pvalue) >= long_cut) & is_can
      p_indexes_wald[[p]]  <- p_indexes_binom[[p]]
      
    } else if(filter_type == 'trig'){
      p_indexes_binom[[p]] <- (as.numeric(l[[p]]$trig.pvalue) >= trig_cut) & is_can
      p_indexes_wald[[p]]  <- p_indexes_binom[[p]]
    }
  }
  
  # ── 4. 写出两个基因列表 ───────────────────────────────────────────────────
  prefix <- gsub('.rds', '_', rds_name)
  
  # 二项检验基因列表
  out_binom <- paste0(prefix, 'NMDesc_binom_enriched_', filter_type, '.txt')
  extract_and_write(p_indexes_binom, l, out_binom)
  
  # Wald log OR 基因列表
  out_wald <- paste0(prefix, 'NMDesc_wald_enriched_', filter_type, '.txt')
  extract_and_write(p_indexes_wald, l, out_wald)
}