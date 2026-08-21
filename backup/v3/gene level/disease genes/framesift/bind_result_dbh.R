###############################################################################
# Frameshift NMDesc enrichment pipeline  (v3)
#
# 与 snv 管线（get_pvalue v5）对称的最终版本：
#   * 检验全集【预先】限定为 OMIM 常染色体显性基因的经典转录本
#     （不再是全基因组检验后再交集）
#   * 主方法：ACAT 合并 (plus1, plus2 的 Fisher p) + Tarone 可检验性过滤 + BH
#       Tier 1  significant : FDR < 0.05
#       Tier 2  suggestive  : 0.05 <= FDR < 0.20
#       主列表              : FDR < 0.20
#   * 对照路径：Simes(+/-Tarone)、pooled Fisher + DBH/DBY/BH、pooled binomial+BH
#
# 假设不变：只要 plus1 或 plus2 任一类移码在 NMDesc 区域富集，即认为该转录本
# 富集（并-交检验）。ACAT 与 Simes 同为 min 型合并，符合该假设；ACAT 在任意
# 相关结构下渐近有效，本数据不含 p=1（can.PTC=0 已跳过），无需上尾截断。
#
# 0.20 提示层阈值的理由（与 snv 管线一致，写入 Methods）：
#   显性 NMD 逃逸基因的致病截短集中在逃逸区（k/0 原型），Fisher 以计数总和为
#   条件，最小可达 p ≈ f^k；k<=2 的原型基因在 FDR 0.05 下结构性不可达。
#   0.20 使 k=2 可达，同时 Fisher/ACAT 保持校准。0.05–0.20 提示层单独标注。
#
# 环境中需已存在：
#   get_syn_count(chrom, start, end)  -> 区间内同义变异数
#   ensembl                           -> biomaRt mart 对象
#   omim_AD_symbols.csv               -> snv 管线第 9 步已生成；若无则见 0b
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(biomaRt)
})

ALPHA_MAIN <- 0.20    # 主列表阈值（提示层上界）
ALPHA_TIER <- 0.05    # 显著层阈值
HAS_DFDR   <- requireNamespace("DiscreteFDR", quietly = TRUE)


###############################################################################
# 0. 输入
###############################################################################

PTC_info      <- read.csv('PTC_info20260201_region.csv')
fs_NMD_result <- read.csv('fs_NMD_result20260201.csv')

# ---- 0b. OMIM-AD 基因全集（与 snv 管线共用同一文件，保证两条管线一致）----
if (!file.exists("omim_AD_symbols.csv"))
  stop("omim_AD_symbols.csv not found — run the OMIM step (snv pipeline section 9) first, ",
       "so both pipelines use the identical AD gene universe.")
omim_AD_symbols <- read.csv("omim_AD_symbols.csv", header = TRUE)$hgnc_symbol
omim_AD_symbols <- unique(trimws(omim_AD_symbols))
omim_AD_symbols <- omim_AD_symbols[!is.na(omim_AD_symbols) & nchar(omim_AD_symbols) > 0]


###############################################################################
# 1. 每转录本 / 每移码类型的变异计数
###############################################################################

count_fs <- fs_NMD_result %>%
  filter(NMDesc_can == TRUE) %>%
  group_by(transcript_id, type) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = type, values_from = count, values_fill = 0)

if ("plus1" %in% names(count_fs)) count_fs <- dplyr::rename(count_fs, plus1_count = plus1)
if ("plus2" %in% names(count_fs)) count_fs <- dplyr::rename(count_fs, plus2_count = plus2)


###############################################################################
# 2. 构建每转录本对象 —— 全集预先限定为 AD 基因的经典转录本
###############################################################################

transcript_set3 <- unique(PTC_info$transcript)

# ---- NEW: 限定到 OMIM-AD 基因的经典转录本（与 snv v5 的限定逻辑一致）----
tx_map <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol",
                               "transcript_is_canonical"),
                filters = "ensembl_transcript_id",
                values  = transcript_set3, mart = ensembl)

keep.tx <- tx_map$ensembl_transcript_id[
  tx_map$hgnc_symbol %in% omim_AD_symbols &
    !is.na(tx_map$transcript_is_canonical) &
    tx_map$transcript_is_canonical == 1
]

n.before        <- length(transcript_set3)
transcript_set3 <- transcript_set3[transcript_set3 %in% keep.tx]
cat(sprintf("Restriction: %d NMDesc transcripts -> %d canonical transcripts of %d AD symbols\n",
            n.before, length(transcript_set3), length(omim_AD_symbols)))
if (length(transcript_set3) == 0)
  stop("No transcripts left after AD restriction — check symbol matching.")

transcript_object <- list()

for (transcript in transcript_set3) {
  
  PTC_plus1_idx <- which(PTC_info$transcript == transcript & PTC_info$type == "plus1")
  PTC_plus2_idx <- which(PTC_info$transcript == transcript & PTC_info$type == "plus2")
  va_plus1_idx  <- which(fs_NMD_result$transcript_id == transcript & fs_NMD_result$type == "plus1")
  va_plus2_idx  <- which(fs_NMD_result$transcript_id == transcript & fs_NMD_result$type == "plus2")
  
  transcript_object[[transcript]] <- list(
    plus1 = list(
      PTC_loc     = PTC_info$PTC_loc[PTC_plus1_idx],
      PTC_NMDesc  = PTC_info$NMDesc[PTC_plus1_idx],
      can_region  = PTC_info$can_region[PTC_plus1_idx],
      css_region  = PTC_info$css_region[PTC_plus1_idx],
      long_region = PTC_info$long_region[PTC_plus1_idx],
      NMDesc_can  = fs_NMD_result$NMDesc_can[va_plus1_idx],
      NMDesc_css  = fs_NMD_result$NMDesc_css[va_plus1_idx],
      NMDesc_long = fs_NMD_result$NMDesc_long[va_plus1_idx]
    ),
    plus2 = list(
      PTC_loc     = PTC_info$PTC_loc[PTC_plus2_idx],
      PTC_NMDesc  = PTC_info$NMDesc[PTC_plus2_idx],
      can_region  = PTC_info$can_region[PTC_plus2_idx],
      css_region  = PTC_info$css_region[PTC_plus2_idx],
      long_region = PTC_info$long_region[PTC_plus2_idx],
      NMDesc_can  = fs_NMD_result$NMDesc_can[va_plus2_idx],
      NMDesc_css  = fs_NMD_result$NMDesc_css[va_plus2_idx],
      NMDesc_long = fs_NMD_result$NMDesc_long[va_plus2_idx]
    )
  )
}

plus1_NMD_result <- lapply(transcript_object, function(x) x$plus1)
plus2_NMD_result <- lapply(transcript_object, function(x) x$plus2)


###############################################################################
# 3. 单表检验 + 逐转录本循环（不变）
###############################################################################

nmd_test_one <- function(can.PTC, rest.PTC, n_can, n_rest) {
  
  out <- list(binom_p = NA_real_, fisher_p = NA_real_, odds_ratio = NA_real_)
  
  if (any(is.na(c(can.PTC, rest.PTC, n_can, n_rest)))) return(out)
  if (n_can <= 0 || n_rest <= 0)                       return(out)
  if (can.PTC < 0 || rest.PTC < 0)                     return(out)
  if (can.PTC > n_can || rest.PTC > n_rest)            return(out)
  
  p0 <- rest.PTC / n_rest                    # 基线排除被检验区域
  if (is.na(p0) || p0 < 0 || p0 > 1)                   return(out)
  
  out$binom_p <- binom.test(can.PTC, n_can, p0, alternative = "greater")$p.value
  
  ft <- fisher.test(matrix(c(can.PTC,  n_can  - can.PTC,
                             rest.PTC, n_rest - rest.PTC), nrow = 2),
                    alternative = "greater")
  out$fisher_p   <- ft$p.value
  out$odds_ratio <- unname(ft$estimate)
  
  out
}


run_fs_pvalue <- function(NMD_result, PTC_info, type_label, cds.info, verbose = TRUE) {
  
  n  <- length(NMD_result)
  tx <- names(NMD_result)
  
  df <- data.frame(
    transcript = tx,
    NMDesc_can_count = NA_integer_, all_variant = NA_integer_, can_region = NA_integer_,
    can.PTC = NA_integer_, rest.PTC = NA_integer_,
    n_can = NA_real_, n_rest = NA_real_, n_all = NA_real_,
    binom_p = NA_real_, fisher_p = NA_real_, odds_ratio = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(n)) {
    
    if (verbose && i %% 200 == 0) cat(type_label, i, "/", n, "\n")
    
    transcript <- tx[i]
    
    fs_NMD_can           <- NMD_result[[i]][["NMDesc_can"]]
    NMDesc_can_count     <- sum(fs_NMD_can == TRUE,  na.rm = TRUE)
    non_NMDesc_can_count <- sum(fs_NMD_can == FALSE, na.rm = TRUE)
    all.variant          <- NMDesc_can_count + non_NMDesc_can_count
    
    cr <- NMD_result[[i]][["can_region"]]
    can_region_val <- if (is.null(cr) || length(cr) == 0) NA_integer_ else cr[1]
    
    df$NMDesc_can_count[i] <- NMDesc_can_count
    df$all_variant[i]      <- all.variant
    df$can_region[i]       <- can_region_val
    
    cds.rows <- which(cds.info$ensembl_transcript_id == transcript)
    if (length(cds.rows) == 0) next
    
    ptc.rows <- which(PTC_info$transcript == transcript & PTC_info$type == type_label)
    if (length(ptc.rows) == 0) next
    
    cds.start <- 1
    cds.end   <- max(cds.info$cds_end[cds.rows], na.rm = TRUE)
    chrom     <- cds.info$chromosome_name[cds.rows][1]
    
    NMDesc.start <- suppressWarnings(min(PTC_info$can_region_start[ptc.rows], na.rm = TRUE))
    NMDesc.end   <- suppressWarnings(max(PTC_info$can_region_end[ptc.rows],   na.rm = TRUE))
    if (!is.finite(NMDesc.start) || !is.finite(NMDesc.end)) next
    
    n_can  <- get_syn_count(chrom, NMDesc.start, NMDesc.end)
    n_all  <- get_syn_count(chrom, cds.start,    cds.end)
    n_rest <- n_all - n_can
    
    df$can.PTC[i]  <- NMDesc_can_count
    df$rest.PTC[i] <- all.variant - NMDesc_can_count
    df$n_can[i]    <- n_can
    df$n_rest[i]   <- n_rest
    df$n_all[i]    <- n_all
    
    if (all.variant == 0)                                next
    if (NMDesc_can_count <= 0)                           next
    if (is.na(can_region_val) || can_region_val <= 0)    next
    
    tt <- nmd_test_one(df$can.PTC[i], df$rest.PTC[i], n_can, n_rest)
    df$binom_p[i]    <- tt$binom_p
    df$fisher_p[i]   <- tt$fisher_p
    df$odds_ratio[i] <- tt$odds_ratio
  }
  
  # 每个检验的最小可达 p 值（Tarone 过滤用）
  df$min_p <- NA_real_
  ok <- !is.na(df$fisher_p)
  if (HAS_DFDR && any(ok)) {
    counts <- as.matrix(data.frame(
      a = df$can.PTC[ok],  b = df$n_can[ok]  - df$can.PTC[ok],
      c = df$rest.PTC[ok], d = df$n_rest[ok] - df$rest.PTC[ok]))
    fp <- DiscreteFDR::fisher.pvalues.support(counts, alternative = "greater")
    df$min_p[ok] <- sapply(fp$support, min)
  }
  
  df
}


###############################################################################
# 4. 分别运行 plus1 / plus2
###############################################################################

get_cds_info <- function(tx_ids) {
  getBM(attributes = c("ensembl_transcript_id", "chromosome_name", "cds_start", "cds_end"),
        filters = "ensembl_transcript_id", values = tx_ids, mart = ensembl)
}

cds.info1 <- get_cds_info(names(plus1_NMD_result))
cds.info2 <- get_cds_info(names(plus2_NMD_result))

plus1_results_df <- run_fs_pvalue(plus1_NMD_result, PTC_info, "plus1", cds.info1)
plus2_results_df <- run_fs_pvalue(plus2_NMD_result, PTC_info, "plus2", cds.info2)

saveRDS(plus1_results_df, "plus1_fs_results20260201_AD.rds")
saveRDS(plus2_results_df, "plus2_fs_results20260201_AD.rds")


###############################################################################
# 5. 合并函数（不变）
###############################################################################

# ---- Simes：k=2 时即 min(2*p_(1), p_(2))，对 p=1 天然免疫 ----
simes_combine <- function(p) {
  p <- sort(p[!is.na(p) & is.finite(p)])
  if (length(p) == 0) return(NA_real_)
  min(length(p) * p / seq_along(p))
}

# ---- ACAT：柯西合并（主方法）----
# 注：本数据集不含 p=1（can.PTC=0 的转录本在 run_fs_pvalue 中已被 next 跳过），
#     故无需上尾截断。
acat_combine <- function(p, w = NULL) {
  keep <- !is.na(p) & is.finite(p)
  p <- p[keep]
  k <- length(p)
  if (k == 0) return(NA_real_)
  if (is.null(w)) w <- rep(1/k, k) else { w <- w[keep]; w <- w / sum(w) }
  
  p <- pmax(p, 1e-300)
  
  x    <- ifelse(p < 1e-15, 1 / (p * pi), tan((0.5 - p) * pi))
  stat <- sum(w * x)
  if (stat > 1e15) 1 / (stat * pi) else 0.5 - atan(stat) / pi
}


###############################################################################
# 6. 合并 plus1 + plus2
###############################################################################

keep_cols <- c("transcript", "can.PTC", "rest.PTC", "n_can", "n_rest", "n_all",
               "fisher_p", "binom_p", "odds_ratio", "min_p")

comb <- merge(plus1_results_df[, keep_cols],
              plus2_results_df[, keep_cols],
              by = "transcript", all = TRUE, suffixes = c("_p1", "_p2"))

pick <- function(a, b) ifelse(!is.na(a), a, b)

comb$n_can  <- pick(comb$n_can_p1,  comb$n_can_p2)
comb$n_rest <- pick(comb$n_rest_p1, comb$n_rest_p2)
comb$n_all  <- pick(comb$n_all_p1,  comb$n_all_p2)

# 有几个移码类型真正参与了检验
comb$n_classes <- (!is.na(comb$fisher_p_p1)) + (!is.na(comb$fisher_p_p2))
comb$one_only  <- comb$n_classes == 1

# ---- 6a. 主方法：ACAT（min 型，任意相关下渐近有效）----
comb$acat_p <- apply(comb[, c("fisher_p_p1", "fisher_p_p2")], 1,
                     function(r) acat_combine(as.numeric(r)))

# ACAT 对每个输入单调递减 → 在最小可达输入处取值即为最小可达输出
comb$min_acat <- mapply(function(m1, m2) {
  m <- c(m1, m2); m <- m[!is.na(m)]
  if (length(m) == 0) return(NA_real_)
  acat_combine(m)
}, comb$min_p_p1, comb$min_p_p2)

# ---- 6b. 敏感性：Simes ----
comb$simes_p <- apply(comb[, c("fisher_p_p1", "fisher_p_p2")], 1,
                      function(r) simes_combine(as.numeric(r)))

comb$min_simes <- mapply(function(m1, m2) {
  m <- c(m1, m2); m <- m[!is.na(m)]
  if (length(m) == 0) return(NA_real_)
  if (length(m) == 1) return(m)
  min(2 * min(m), max(m), 1)
}, comb$min_p_p1, comb$min_p_p2)

# ---- 6c. 对照：pooled 计数单次检验（检验"整体富集"，假设不同）----
comb$can.PTC  <- rowSums(cbind(comb$can.PTC_p1,  comb$can.PTC_p2),  na.rm = TRUE)
comb$rest.PTC <- rowSums(cbind(comb$rest.PTC_p1, comb$rest.PTC_p2), na.rm = TRUE)

no_data <- is.na(comb$n_can_p1) & is.na(comb$n_can_p2)
comb$can.PTC[no_data]  <- NA_integer_
comb$rest.PTC[no_data] <- NA_integer_

pooled <- mapply(nmd_test_one, comb$can.PTC, comb$rest.PTC, comb$n_can, comb$n_rest,
                 SIMPLIFY = FALSE)
comb$pooled_fisher_p <- vapply(pooled, function(x) x$fisher_p,   numeric(1))
comb$pooled_binom_p  <- vapply(pooled, function(x) x$binom_p,    numeric(1))
comb$pooled_OR       <- vapply(pooled, function(x) x$odds_ratio, numeric(1))

# ---- 6d. 驱动来源 + k/0 原型标注 ----
comb$driver <- with(comb, ifelse(
  n_classes == 0, NA_character_,
  ifelse(n_classes == 1, "single_class",
         ifelse(fisher_p_p1 < fisher_p_p2, "plus1",
                ifelse(fisher_p_p2 < fisher_p_p1, "plus2", "tie")))))

comb$archetype <- sprintf("%d/%d", comb$can.PTC, comb$rest.PTC)
comb$is_k0     <- !is.na(comb$rest.PTC) & comb$rest.PTC == 0 &
  !is.na(comb$can.PTC)  & comb$can.PTC  > 0


###############################################################################
# 7. 多重检验校正 —— 主方法 ACAT + Tarone-BH，双层
###############################################################################

# ---- 7a. Tarone/Gilbert 过滤 + BH ----
tarone_bh <- function(p, min_p, alpha = ALPHA_MAIN) {
  out      <- rep(NA_real_, length(p))
  testable <- !is.na(p) & !is.na(min_p) & min_p <= alpha
  if (any(testable)) out[testable] <- p.adjust(p[testable], method = "BH")
  attr(out, "n_testable") <- sum(testable)
  out
}

# 主方法
comb$fdr_acat_filt  <- tarone_bh(comb$acat_p,  comb$min_acat)

# 对照
comb$fdr_acat       <- p.adjust(comb$acat_p,  method = "BH")
comb$fdr_simes_filt <- tarone_bh(comb$simes_p, comb$min_simes)
comb$fdr_simes      <- p.adjust(comb$simes_p, method = "BH")

# ---- 主方法双层标签 ----
comb$tier <- ifelse(is.na(comb$fdr_acat_filt), "ns",
                    ifelse(comb$fdr_acat_filt < ALPHA_TIER, "significant",
                           ifelse(comb$fdr_acat_filt < ALPHA_MAIN, "suggestive", "ns")))

# 提示层里 0.05 下结构性不可达的转录本（k/0 论证的证据链）
n.sig05 <- sum(comb$tier == "significant", na.rm = TRUE)
n.testable <- attr(comb$fdr_acat_filt, "n_testable")
thr05   <- ALPHA_TIER * max(n.sig05, 1) / max(n.testable, 1)
comb$unreachable.at.005 <- !is.na(comb$min_acat) & comb$min_acat > thr05

# ---- 7b. pooled 对照路径：DBH / DBY / BH ----
comb$sig.dbh          <- FALSE
comb$sig.dby          <- FALSE
comb$min_pooled       <- NA_real_
comb$fdr_pooled       <- p.adjust(comb$pooled_fisher_p, method = "BH")
comb$fdr_pooled_binom <- p.adjust(comb$pooled_binom_p,  method = "BH")

ok     <- !is.na(comb$pooled_fisher_p)
n.test <- sum(ok)

if (HAS_DFDR && n.test > 0) {
  counts <- as.matrix(data.frame(
    a = comb$can.PTC[ok],  b = comb$n_can[ok]  - comb$can.PTC[ok],
    c = comb$rest.PTC[ok], d = comb$n_rest[ok] - comb$rest.PTC[ok]))
  fp <- DiscreteFDR::fisher.pvalues.support(counts, alternative = "greater")
  comb$min_pooled[ok] <- sapply(fp$support, min)
  comb$sig.dbh[which(ok)[DiscreteFDR::DBH(fp$raw, fp$support, alpha = ALPHA_MAIN)$Indices]] <- TRUE
  comb$sig.dby[which(ok)[DiscreteFDR::DBY(fp$raw, fp$support, alpha = ALPHA_MAIN)$Indices]] <- TRUE
}

cat('\n================ Frameshift NMDesc enrichment (AD-restricted) ================\n')
cat(sprintf('检验的转录本数（AD 全集）      : %d\n', sum(comb$n_classes > 0, na.rm = TRUE)))
cat(sprintf('  仅单一移码类型有数据        : %d\n', sum(comb$one_only, na.rm = TRUE)))
cat(sprintf('  ACAT 可检验（Tarone）       : %d\n', attr(comb$fdr_acat_filt, "n_testable")))
cat(sprintf('--- 主方法 ACAT + Tarone-BH，双层 ---\n'))
cat(sprintf('  Tier 1 significant (FDR < %.2f)        : %d\n', ALPHA_TIER, n.sig05))
cat(sprintf('  Tier 2 suggestive  (%.2f <= FDR < %.2f) : %d\n', ALPHA_TIER, ALPHA_MAIN,
            sum(comb$tier == "suggestive", na.rm = TRUE)))
cat(sprintf('    其中 k/0 原型                        : %d\n',
            sum(comb$tier == "suggestive" & comb$is_k0, na.rm = TRUE)))
cat(sprintf('    0.05 下结构性不可达                  : %d\n',
            sum(comb$tier == "suggestive" & comb$unreachable.at.005, na.rm = TRUE)))
cat(sprintf('  主列表 (FDR < %.2f) 合计               : %d\n', ALPHA_MAIN,
            sum(comb$tier %in% c("significant","suggestive"))))
cat(sprintf('--- 对照路径，FDR < %.2f 显著数 ---\n', ALPHA_MAIN))
cat(sprintf('  ACAT  + BH（未过滤）        : %d\n', sum(comb$fdr_acat        < ALPHA_MAIN, na.rm = TRUE)))
cat(sprintf('  Simes + Tarone-BH           : %d\n', sum(comb$fdr_simes_filt  < ALPHA_MAIN, na.rm = TRUE)))
cat(sprintf('  Simes + BH                  : %d\n', sum(comb$fdr_simes       < ALPHA_MAIN, na.rm = TRUE)))
cat(sprintf('  pooled + DBH                : %d\n', sum(comb$sig.dbh)))
cat(sprintf('  pooled + DBY                : %d\n', sum(comb$sig.dby)))
cat(sprintf('  pooled + BH                 : %d\n', sum(comb$fdr_pooled      < ALPHA_MAIN, na.rm = TRUE)))
cat(sprintf('  pooled binomial + BH        : %d\n', sum(comb$fdr_pooled_binom< ALPHA_MAIN, na.rm = TRUE)))
cat('==============================================================================\n\n')

saveRDS(comb,   "fs_combined_results20260201_AD_v3.rds")
write.csv(comb, "fs_combined_results20260201_AD_v3.csv", row.names = FALSE)


###############################################################################
# 8. 转录本 -> 基因：主方法双层输出 + 对照方法输出
#    （全集已是 AD，无需再做 OMIM 交集）
###############################################################################

tx2gene <- function(tx_ids) {
  if (length(tx_ids) == 0) return(character(0))
  bm <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"),
              filters = "ensembl_transcript_id", values = tx_ids, mart = ensembl)
  g <- unique(bm$hgnc_symbol)
  g[!is.na(g) & nchar(g) > 0]
}

# ---- 主方法：双层 ----
tier_sets <- list(
  significant = comb$transcript[comb$tier == "significant"],
  suggestive  = comb$transcript[comb$tier == "suggestive"],
  all         = comb$transcript[comb$tier %in% c("significant","suggestive")]
)
tier_genes <- lapply(tier_sets, tx2gene)

writeLines(c("hgnc_symbol", tier_genes$significant),
           sprintf("fs_can_AD_acat_significant_FDR%.2f.txt", ALPHA_TIER))
writeLines(c("hgnc_symbol", tier_genes$suggestive),
           sprintf("fs_can_AD_acat_suggestive_%.2f-%.2f.txt", ALPHA_TIER, ALPHA_MAIN))
writeLines(c("hgnc_symbol", tier_genes$all),
           sprintf("fs_can_AD_acat_FDR%.2f_all.txt", ALPHA_MAIN))
cat(sprintf('主方法 significant : %4d genes\n', length(tier_genes$significant)))
cat(sprintf('主方法 suggestive  : %4d genes\n', length(tier_genes$suggestive)))
cat(sprintf('主方法 main (0.20) : %4d genes\n', length(tier_genes$all)))

# ---- 对照方法（均在 ALPHA_MAIN 下）----
sig_sets <- list(
  acat       = comb$transcript[!is.na(comb$fdr_acat)         & comb$fdr_acat         < ALPHA_MAIN],
  simes_filt = comb$transcript[!is.na(comb$fdr_simes_filt)   & comb$fdr_simes_filt   < ALPHA_MAIN],
  simes      = comb$transcript[!is.na(comb$fdr_simes)        & comb$fdr_simes        < ALPHA_MAIN],
  dbh        = comb$transcript[comb$sig.dbh],
  dby        = comb$transcript[comb$sig.dby],
  pooled     = comb$transcript[!is.na(comb$fdr_pooled)       & comb$fdr_pooled       < ALPHA_MAIN],
  binom      = comb$transcript[!is.na(comb$fdr_pooled_binom) & comb$fdr_pooled_binom < ALPHA_MAIN]
)
gene_sets <- lapply(sig_sets, tx2gene)

for (m in names(gene_sets)) {
  fn <- sprintf("fs_can_AD_%s_FDR%.2f.txt", m, ALPHA_MAIN)
  writeLines(c("hgnc_symbol", gene_sets[[m]]), fn)
  cat(sprintf('%-11s : %4d genes -> %s\n', m, length(gene_sets[[m]]), fn))
}

saveRDS(c(tier_genes, gene_sets), "AD_gene_sets_v3.rds")


###############################################################################
# 9. 全注释表（补充表底稿；与 snv 管线格式对齐）
###############################################################################

bm_all <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol"),
                filters = "ensembl_transcript_id",
                values = comb$transcript, mart = ensembl)
ann <- merge(comb, bm_all, by.x = "transcript", by.y = "ensembl_transcript_id",
             all.x = TRUE)
ann <- ann[order(ann$fdr_acat_filt, ann$acat_p),
           c("hgnc_symbol","transcript","archetype","is_k0","driver",
             "can.PTC","rest.PTC",
             "fisher_p_p1","fisher_p_p2","acat_p","fdr_acat_filt",
             "simes_p","fdr_simes_filt","pooled_fisher_p","sig.dbh","sig.dby",
             "min_acat","unreachable.at.005","tier")]
ann$rank <- seq_len(nrow(ann))
write.csv(ann, "fs_AD_acat_full_annotated.csv", row.names = FALSE)


###############################################################################
# 10. 诊断（对齐主方法为 ACAT）
###############################################################################

cat('\n--- 诊断 ---\n')

n_p1 <- sum(comb$fisher_p_p1 == 1 | comb$fisher_p_p2 == 1, na.rm = TRUE)
cat(sprintf('至少一类 p=1 的转录本 : %d / %d (%.1f%%)\n',
            n_p1, sum(comb$n_classes == 2, na.rm = TRUE),
            100 * n_p1 / max(sum(comb$n_classes == 2, na.rm = TRUE), 1)))

has1 <- with(comb, !is.na(fisher_p_p1) & !is.na(fisher_p_p2) &
               (fisher_p_p1 == 1 | fisher_p_p2 == 1))
cat(sprintf('ACAT vs Simes  Spearman (全部)      : %.3f\n',
            cor(comb$acat_p, comb$simes_p, use = "complete.obs", method = "spearman")))
cat(sprintf('ACAT vs Simes  Spearman (剔除 p=1)  : %.3f\n',
            cor(comb$acat_p[!has1], comb$simes_p[!has1],
                use = "complete.obs", method = "spearman")))
cat(sprintf('ACAT vs pooled Spearman             : %.3f\n',
            cor(comb$acat_p, comb$pooled_fisher_p, use = "complete.obs", method = "spearman")))

sig_tx <- comb$transcript[comb$tier %in% c("significant","suggestive")]
cat('\n主结果（ACAT+Tarone，FDR<0.20）显著转录本的驱动来源:\n')
print(table(comb$driver[comb$transcript %in% sig_tx], useNA = "ifany"))

cat('\nn_all 分布（若按染色体近似恒定，说明 cds.start=1 与基因组坐标不匹配）:\n')
print(summary(comb$n_all))
cat('n_can 分布:\n')
print(summary(comb$n_can))

methods <- names(gene_sets)
overlap <- matrix(0L, length(methods), length(methods), dimnames = list(methods, methods))
for (a in methods) for (b in methods) overlap[a, b] <- length(intersect(gene_sets[[a]], gene_sets[[b]]))
cat('\n--- 方法间基因集重叠（对角线为集合大小）---\n')
print(overlap)