get_NMD_enrichment_DBH = function(rds_name = 'snv_plp_ptc_nmdesc_can_p_f_syn_20260201_Jul30.rds',
                                   FDR = 0.05, filter_type = 'can'){
  l <- readRDS(rds_name)
  
  # ── 安全提取：字段缺失时返回 NA，避免 sapply 退化成 list ────────────────
  getnum <- function(field) sapply(l, function(x){
    v <- x[[field]]; if (is.null(v) || length(v) == 0) NA_real_ else as.numeric(v[1])
  })
  getlog <- function(field) sapply(l, function(x){
    v <- x[[field]]; if (is.null(v) || length(v) == 0) FALSE else isTRUE(as.logical(v[1]))
  })
  
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
    new.genes <- unique(new.genes)
    cat('✅ 找到基因数:', length(new.genes), '→', out_name, '\n')
    writeLines(new.genes, out_name)
    return(invisible(new.genes))
  }
  
  # ── 1. 二项检验 p 值 → BH ────────────────────────────────────────────────
  can.pvalues     <- getnum('can.pvalue')
  can.pvalues.adj <- p.adjust(can.pvalues, method = "BH")
  for (p in 1:length(l)) l[[p]]$can.pvalue.adj <- can.pvalues.adj[p]
  cat('二项检验 FDR < cut 的转录本数:', sum(can.pvalues.adj < FDR, na.rm=TRUE), '\n')
  
  # ── 2. Fisher 精确检验 p 值 → BH（替代原来的 Wald，本文件无 wald 字段）──
  can.fisher.pvalues     <- getnum('can.fisher.pvalue')
  can.fisher.pvalues.adj <- p.adjust(can.fisher.pvalues, method = "BH")
  for (p in 1:length(l)) l[[p]]$can.fisher.pvalue.adj <- can.fisher.pvalues.adj[p]
  cat('Fisher检验 FDR < cut 的转录本数:', sum(can.fisher.pvalues.adj < FDR, na.rm=TRUE), '\n')
  
  # ── 3. 离散 FDR：直接用 get_pvalue 里存好的标记 ─────────────────────────
  sig.dbh <- getlog('sig.dbh')
  sig.dby <- getlog('sig.dby')
  cat('DBH 显著转录本数:', sum(sig.dbh), '  DBY 显著转录本数:', sum(sig.dby), '\n')
  
  # ── 4. 构建各方法的 p_indexes ────────────────────────────────────────────
  p_indexes_binom  <- list()
  p_indexes_fisher <- list()
  p_indexes_dbh    <- list()
  p_indexes_dby    <- list()
  
  for (p in 1:length(l)) {
    is_can <- isTRUE(as.logical(l[[p]]$is_canonical))
    
    if(filter_type == 'can'){
      p_indexes_binom[[p]]  <- isTRUE(can.pvalues.adj[p]        < FDR) & is_can
      p_indexes_fisher[[p]] <- isTRUE(can.fisher.pvalues.adj[p] < FDR) & is_can
      p_indexes_dbh[[p]]    <- sig.dbh[p] & is_can
      p_indexes_dby[[p]]    <- sig.dby[p] & is_can
      
    } else if(filter_type == 'nofilter'){   # 不限 canonical
      p_indexes_binom[[p]]  <- isTRUE(can.pvalues.adj[p]        < FDR)
      p_indexes_fisher[[p]] <- isTRUE(can.fisher.pvalues.adj[p] < FDR)
      p_indexes_dbh[[p]]    <- sig.dbh[p]
      p_indexes_dby[[p]]    <- sig.dby[p]
      
    } else {
      stop("本文件只含 can.* 字段，filter_type 仅支持 'can' 或 'nofilter'")
    }
  }
  
  # ── 5. 写出四个基因列表 ──────────────────────────────────────────────────
  prefix <- gsub('.rds', '_', rds_name)
  
  extract_and_write(p_indexes_binom,  l, paste0(prefix, 'NMDesc_binom_enriched_',  filter_type, '.txt'))
  extract_and_write(p_indexes_fisher, l, paste0(prefix, 'NMDesc_fisher_enriched_', filter_type, '.txt'))
  extract_and_write(p_indexes_dbh,    l, paste0(prefix, 'NMDesc_dbh_enriched_',    filter_type, '.txt'))
  extract_and_write(p_indexes_dby,    l, paste0(prefix, 'NMDesc_dby_enriched_',    filter_type, '.txt'))
  
  invisible(l)
}