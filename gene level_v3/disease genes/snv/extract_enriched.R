extract_enriched <- function(txnames.list, fdr.method = c("dbh", "dby", "bh"),
                             canonical.only = TRUE) {
  fdr.method <- match.arg(fdr.method)
  
  get <- function(x, field, default = NA) {
    v <- x[[field]]
    if (is.null(v) || length(v) == 0) default else v[1]
  }
  
  df <- do.call(rbind, lapply(seq_along(txnames.list), function(i) {
    x <- txnames.list[[i]]
    if (length(x) == 0 || is.null(x$can.fisher.pvalue)) return(NULL)
    data.frame(
      txname       = get(x, "txname", NA_character_),
      hgnc_symbol  = get(x, "hgnc_symbol", NA_character_),
      is_canonical = get(x, "is_canonical", NA),
      can.PTC      = get(x, "can.PTC", NA_real_),
      rest.PTC     = get(x, "rest.PTC", NA_real_),
      n_can        = get(x, "n_can", NA_real_),
      n_rest       = get(x, "n_rest", NA_real_),
      p.fisher     = get(x, "can.fisher.pvalue", NA_real_),
      fdr.bh       = get(x, "fdr.bh", NA_real_),
      sig.dbh      = get(x, "sig.dbh", FALSE),
      sig.dby      = get(x, "sig.dby", FALSE),
      min.p        = get(x, "min.attainable.p", NA_real_),
      stringsAsFactors = FALSE
    )
  }))
  
  # observed vs expected enrichment ratio
  df$OR <- with(df, (can.PTC / n_can) / (rest.PTC / n_rest))
  # flag tests that cannot reach significance regardless of effect size
  df$unreachable <- df$min.p > 0.05
  
  if (canonical.only) {
    df <- df[which(df$is_canonical == 1 | df$is_canonical == TRUE), ]
  }
  
  sig <- switch(fdr.method,
                dbh = df$sig.dbh,
                dby = df$sig.dby,
                bh  = df$fdr.bh < 0.05)
  sig[is.na(sig)] <- FALSE
  
  enriched <- df[sig, ]
  enriched <- enriched[order(enriched$p.fisher), ]
  
  message(sprintf("Transcripts tested: %d", nrow(df)))
  message(sprintf("Unreachable (min attainable p > 0.05): %d (%.1f%%)",
                  sum(df$unreachable, na.rm = TRUE),
                  100 * mean(df$unreachable, na.rm = TRUE)))
  message(sprintf("Enriched by %s: %d transcripts, %d unique genes",
                  toupper(fdr.method), nrow(enriched),
                  length(unique(na.omit(enriched$hgnc_symbol)))))
  
  list(all = df, enriched = enriched,
       genes = sort(unique(na.omit(enriched$hgnc_symbol))))
}