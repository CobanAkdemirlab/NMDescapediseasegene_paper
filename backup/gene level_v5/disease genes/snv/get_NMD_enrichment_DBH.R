get_NMD_enrichment_DBH = function(rds_name = 'snv_plp_ptc_nmdesc_can_p_f_syn_20260201_Jul30.rds',
                                   FDR = 0.05, filter_type = 'can'){
  l <- readRDS(rds_name)
  
  # ── Safe extraction: returns NA for missing fields ────────────────
  getnum <- function(field) sapply(l, function(x){
    v <- x[[field]]; if (is.null(v) || length(v) == 0) NA_real_ else as.numeric(v[1])
  })
  getlog <- function(field) sapply(l, function(x){
    v <- x[[field]]; if (is.null(v) || length(v) == 0) FALSE else isTRUE(as.logical(v[1]))
  })
  
  # ── Helper: extract gene names from p_indexes and write out ──────────────────────────────
  extract_and_write <- function(p_indexes, l, out_name) {
    ind.1 <- which(unlist(p_indexes) == TRUE)
    if (length(ind.1) == 0) {
      cat('No genes found matching criteria, writing empty file:', out_name, '\n')
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
    cat('Genes found:', length(new.genes), '->', out_name, '\n')
    writeLines(new.genes, out_name)
    return(invisible(new.genes))
  }
  
  # ── 1. Binomial test p-values -> BH ────────────────────────────────────────────────
  can.pvalues     <- getnum('can.pvalue')
  can.pvalues.adj <- p.adjust(can.pvalues, method = "BH")
  for (p in 1:length(l)) l[[p]]$can.pvalue.adj <- can.pvalues.adj[p]
  cat('Transcripts with binomial FDR < cut:', sum(can.pvalues.adj < FDR, na.rm=TRUE), '\n')
  
  # ── 2. Fisher exact test p-values -> BH (no wald field in this file) ──
  can.fisher.pvalues     <- getnum('can.fisher.pvalue')
  can.fisher.pvalues.adj <- p.adjust(can.fisher.pvalues, method = "BH")
  for (p in 1:length(l)) l[[p]]$can.fisher.pvalue.adj <- can.fisher.pvalues.adj[p]
  cat('Transcripts with Fisher FDR < cut:', sum(can.fisher.pvalues.adj < FDR, na.rm=TRUE), '\n')
  
  # ── 3. Discrete FDR: uses flags stored by get_pvalue ─────────────────────────
  sig.dbh <- getlog('sig.dbh')
  sig.dby <- getlog('sig.dby')
  cat('DBH significant transcripts:', sum(sig.dbh), '  DBY significant transcripts:', sum(sig.dby), '\n')
  
  # ── 4. Build p_indexes for each method ────────────────────────────────────────────
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
      
    } else if(filter_type == 'nofilter'){   # not restricted to canonical
      p_indexes_binom[[p]]  <- isTRUE(can.pvalues.adj[p]        < FDR)
      p_indexes_fisher[[p]] <- isTRUE(can.fisher.pvalues.adj[p] < FDR)
      p_indexes_dbh[[p]]    <- sig.dbh[p]
      p_indexes_dby[[p]]    <- sig.dby[p]
      
    } else {
      stop("This file only contains can.* fields; filter_type supports only 'can' or 'nofilter'")
    }
  }
  
  # ── 5. Write out four gene lists ──────────────────────────────────────────────────
  prefix <- gsub('.rds', '_', rds_name)
  
  extract_and_write(p_indexes_binom,  l, paste0(prefix, 'NMDesc_binom_enriched_',  filter_type, '.txt'))
  extract_and_write(p_indexes_fisher, l, paste0(prefix, 'NMDesc_fisher_enriched_', filter_type, '.txt'))
  extract_and_write(p_indexes_dbh,    l, paste0(prefix, 'NMDesc_dbh_enriched_',    filter_type, '.txt'))
  extract_and_write(p_indexes_dby,    l, paste0(prefix, 'NMDesc_dby_enriched_',    filter_type, '.txt'))
  
  invisible(l)
}