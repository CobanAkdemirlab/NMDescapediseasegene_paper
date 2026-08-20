# =============================================================================
# nmdesc_paired_analysis.R  (v2)
# -----------------------------------------------------------------------------
# NMDesc 疾病基因 vs 对照基因的成对特征比较，四组结构：
#   snv / snv_control / fs / fs_control
#
# 主分析（matched）
#   1. 按 CDS 长度为每个疾病基因匹配一个对照基因（1:1，无放回）
#   2. 用匹配结果构建 gene_all，每行带 pair_id 与 group
#   3. 通过 getBM 获取特征所需信息
#   4. 在配对内部比较特征（二分类用 McNemar/精确二项，连续用 Wilcoxon 符号秩）
#   5. BH 校正，报告显著特征
#
# 敏感性分析（random）
#   同样 1:1，但对照随机抽取、不按 CDS 匹配；配对是任意的，故改用非配对检验
#   （Fisher / Wilcoxon 秩和），同样做 BH 校正。
#
# 两者的对比本身就是结果：主分析控制了 CDS 长度这个混杂，敏感性分析没有。
# 只在敏感性分析里显著的特征，很可能是 CDS 长度的代理而非真实信号。
#
# v2 修正：
#   - fetch_canonical_cds 拆成两次 getBM。cds_length 属于 sequences 属性页，
#     hgnc_symbol / transcript_is_canonical 属于 feature_page，同一次请求不能跨页
#   - annotate_pfam 同样做了跨页防护
#   - 四组的描述性统计、分组输出、组间重叠诊断
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(readxl)
  library(stringr); library(ggplot2); library(patchwork)
  library(data.table); library(biomaRt)
  library(GenomicRanges); library(IRanges)
})

select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
mutate <- dplyr::mutate; summarise <- dplyr::summarise
setdiff <- dplyr::setdiff; union <- dplyr::union; intersect <- dplyr::intersect

# --- 路径解析层（必须在 CONFIG 之前，CONFIG 里用到 data_file()）-----------
# 数据文件用 data_file("文件名") 定位；换数据位置只改 gene level_v3/lib/paths.R
.p <- c("gene level_v3/lib/paths.R", "lib/paths.R", "../lib/paths.R",
        "../../gene level_v3/lib/paths.R")
.p <- .p[file.exists(.p)]
if (!length(.p)) stop("找不到 paths.R —— 请从仓库根目录运行 R")
source(.p[1]); rm(.p)
# --------------------------------------------------------------------------


# =============================================================================
# 0. 配置
# =============================================================================

CONFIG <- list(
  ensembl_version  = NULL,          # NULL = 当前版本
  ensembl_mirror   = NULL,          # "useast" / "uswest" / "asia"
  chunk_sequence   = 25,            # 含 coding 属性（完整 CDS 序列）
  chunk_light      = 300,           # 其余属性
  
  ptc_info         = "PTC_info20260201_region.csv",
  omim_ad_symbols  = "omim_AD_symbols.csv",
  snv_gene_list    = "snv_dbh_AD.csv",
  fs_gene_list     = "fs_can_AD_gene20260201_simes_FDR0.05.txt",
  
  # 外部特征数据（缺失则跳过对应特征）
  # 原来这些是旧机器 /Users/jxu14/ 的绝对路径，换机器必断。
  # 现在走 data_file()：按文件名在 DATA_ROOTS 里查找，找不到返回 NA，
  # 下游的 file.exists() 检查会自动跳过对应特征。
  # 注意 ppi_file 原为 "human (1) (2).txt" —— 那是 "human (1).txt" 的重复下载副本。
  ppi_file         = data_file("human (1).txt"),
  # GTEx 是【可选】输入：annotate_tau() 里有 file.exists() 检查，缺失就把 tau 填 NA
  # 继续跑。所以这里必须用 data_file(must=FALSE)（返回 NA），不能用 missing_file()
  # —— 后者会 stop()，构造 CONFIG 时就崩，下游的跳过逻辑永远走不到。
  gtex_path        = data_file("GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct",
                               must = FALSE),
  lof_metrics_path = data_file("gnomad.v2.1.1.lof_metrics.by_gene.txt"),
  path_touni       = data_file("NIHMS1818854-supplement-2(A).csv"),
  path_motif       = data_file("NIHMS1818854-supplement-2(B).csv"),
  path_lcs         = data_file("Copy of NIHMS1818854-supplement-2.xls"),
  
  pool_cache       = "omim_AD_cds_length.csv",
  
  caliper_frac     = 0.5,           # |log2(对照/病例)| 上限，0.5 ≈ 1.41 倍
  alpha            = 0.05,
  seed             = 1
)

GROUP_LEVELS <- c("snv", "snv_control", "fs", "fs_control")

GROUP_COLORS <- c("snv" = "#2ca02c", "snv_control" = "#98df8a",
                  "fs"  = "#D07A3A", "fs_control"  = "#E2C7B5")

CONTINUOUS_FEATURES <- c(
  "gc_content", "nmdesc_gc_content",
  "repeat_fraction", "nmdesc_repeat_fraction",
  "homopolymer_fraction", "nmdesc_homopolymer_fraction",
  "pfam_overlap_fraction", "pfam_overlap_length",
  "exon_num", "degree_centrality", "tau", "pLI", "oe_lof_upper"
)

BINARY_FEATURES <- c(
  "pfam_overlap_flag", "ppi_overlap",
  "gene_protein_flag", "gene_domains_flag", "gene_slim_flag",
  "gene_morf_flag", "gene_ptm_flag", "gene_nls_flag", "gene_LCS_flag"
)


# =============================================================================
# 1. biomaRt 辅助
# =============================================================================

connect_ensembl <- function(version = CONFIG$ensembl_version,
                            mirror  = CONFIG$ensembl_mirror) {
  args <- list(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  if (!is.null(version)) args$version <- version
  if (!is.null(mirror))  args$mirror  <- mirror
  mart <- do.call(useEnsembl, args)
  
  probe <- tryCatch(
    getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id"),
          filters = "hgnc_symbol", values = "MYH9", mart = mart),
    error = function(e) stop("Ensembl 连接失败: ", e$message))
  message(sprintf("Ensembl 已连接（MYH9 探针返回 %d 个转录本）", nrow(probe)))
  mart
}


## getBM 用 POST 提交；coding 属性返回完整 CDS 序列，几百个基因就是几 MB，
## 服务器会返回 HTTP 405。分批 + 重试。
## 跨属性页的错误无法靠重试解决，直接抛出并给出可操作的提示。
getBM_chunked <- function(attributes, filters, values, mart,
                          chunk = CONFIG$chunk_light, pause = 0.2,
                          max_retry = 3, verbose = TRUE) {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & nchar(values) > 0]
  if (length(values) == 0) return(data.frame())
  
  idx <- split(seq_along(values), ceiling(seq_along(values) / chunk))
  out <- lapply(seq_along(idx), function(k) {
    if (verbose && length(idx) > 1)
      cat(sprintf("    chunk %d/%d (%d)\n", k, length(idx), length(idx[[k]])))
    res <- NULL
    for (attempt in seq_len(max_retry)) {
      res <- tryCatch(
        getBM(attributes = attributes, filters = filters,
              values = values[idx[[k]]], mart = mart),
        error = function(e) {
          if (grepl("multiple attribute pages", e$message))
            stop("属性跨页: ", paste(attributes, collapse = ", "),
                 "\n需拆成多次 getBM 后按 ID 合并。",
                 "\n用 listAttributes(mart)[, c('name','page')] 查看归属页。",
                 call. = FALSE)
          message("      retry ", attempt, ": ", e$message); NULL
        })
      if (!is.null(res)) break
      Sys.sleep(2 * attempt)
    }
    if (is.null(res)) warning("chunk ", k, " 失败，", length(idx[[k]]), " 个值被丢弃")
    Sys.sleep(pause)
    res
  })
  dplyr::bind_rows(out)
}


# =============================================================================
# 2. 基因列表读取
# =============================================================================

read_gene_list <- function(path, label = NULL) {
  if (!file.exists(path)) stop("文件不存在: ", path)
  x <- if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    if (ncol(df) == 0) character(0)
    else if ("hgnc_symbol" %in% names(df)) df$hgnc_symbol else df[[1]]
  } else readLines(path, warn = FALSE)
  
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nchar(x) > 0]
  x <- x[!x %in% c("hgnc_symbol", "x", "V1", "transcript", "ensembl_transcript_id")]
  x <- unique(x)
  if (!is.null(label)) message(sprintf("%-14s %4d 个  (%s)", label, length(x), basename(path)))
  x
}


# =============================================================================
# 3. 对照构建
# =============================================================================

## cds_length 在 sequences 属性页，hgnc_symbol / transcript_is_canonical 在
## feature_page，同一次 getBM 不能跨页 —— 拆两步，第二步只查 canonical 转录本
fetch_canonical_cds <- function(symbols, mart, chunk = CONFIG$chunk_light) {
  symbols <- unique(trimws(as.character(symbols)))
  symbols <- symbols[!is.na(symbols) & nchar(symbols) > 0]
  if (length(symbols) == 0) return(data.frame())
  
  message("  [1/2] canonical 转录本 ...")
  tx <- getBM_chunked(
    attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_is_canonical"),
    filters = "hgnc_symbol", values = symbols, mart = mart, chunk = chunk) %>%
    filter(transcript_is_canonical == 1) %>%
    distinct(hgnc_symbol, .keep_all = TRUE) %>%
    select(hgnc_symbol, ensembl_transcript_id)
  
  if (nrow(tx) == 0) return(data.frame())
  
  message("  [2/2] CDS 长度 ...")
  len <- getBM_chunked(
    attributes = c("ensembl_transcript_id", "cds_length"),
    filters = "ensembl_transcript_id", values = tx$ensembl_transcript_id,
    mart = mart, chunk = chunk) %>%
    filter(!is.na(cds_length), cds_length > 0) %>%
    distinct(ensembl_transcript_id, .keep_all = TRUE)
  
  tx %>% inner_join(len, by = "ensembl_transcript_id") %>%
    select(hgnc_symbol, ensembl_transcript_id, cds_length)
}


## 候选池缓存：OMIM AD 基因的 CDS 长度在多次运行间不变
get_pool_cds <- function(omim_AD_symbols, mart, cache = CONFIG$pool_cache) {
  if (file.exists(cache)) {
    message("读取候选池缓存: ", cache)
    return(read.csv(cache, stringsAsFactors = FALSE))
  }
  message("首次构建候选池 CDS 长度表（此后走缓存）...")
  pool <- fetch_canonical_cds(omim_AD_symbols, mart)
  write.csv(pool, cache, row.names = FALSE)
  pool
}


## 主分析对照：按 CDS 长度最近邻，1:1 无放回
##
## 匹配尺度用 log2 比值而非绝对 bp 差 —— CDS 长度跨数量级，300 vs 400 bp 是
## 巨大差异，10000 vs 10100 bp 几乎相同。
## 处理顺序"极端优先"：贪心匹配对顺序敏感，若先匹配中等长度基因，最长和最短
## 的只能捡剩下的，所以按距池子中位数的距离降序处理。
build_controls_matched <- function(case_df, pool_df,
                                   caliper_frac = CONFIG$caliper_frac,
                                   seed = CONFIG$seed, label = "matched") {
  set.seed(seed)
  if (nrow(case_df) == 0 || nrow(pool_df) == 0) {
    warning(label, ": 病例或候选池为空"); return(data.frame())
  }
  
  pool_median <- stats::median(pool_df$cds_length)
  case_df <- case_df[order(-abs(log2(case_df$cds_length / pool_median))), , drop = FALSE]
  
  available <- rep(TRUE, nrow(pool_df))
  pairs <- vector("list", nrow(case_df))
  
  for (i in seq_len(nrow(case_df))) {
    d <- abs(log2(pool_df$cds_length / case_df$cds_length[i]))
    d[!available] <- Inf
    j <- which.min(d)
    if (length(j) == 0 || !is.finite(d[j])) next
    if (!is.null(caliper_frac) && d[j] > caliper_frac) next
    available[j] <- FALSE
    
    pairs[[i]] <- data.frame(
      disease_gene = case_df$hgnc_symbol[i],
      disease_tx   = case_df$ensembl_transcript_id[i],
      disease_cds  = case_df$cds_length[i],
      control_gene = pool_df$hgnc_symbol[j],
      control_tx   = pool_df$ensembl_transcript_id[j],
      control_cds  = pool_df$cds_length[j],
      log2_ratio   = log2(pool_df$cds_length[j] / case_df$cds_length[i]),
      stringsAsFactors = FALSE)
  }
  bind_rows(pairs)
}


## 敏感性分析对照：随机抽取，1:1 无放回，不看 CDS 长度
build_controls_random <- function(case_df, pool_df,
                                  seed = CONFIG$seed, label = "random") {
  set.seed(seed + 1000)
  if (nrow(case_df) == 0 || nrow(pool_df) == 0) {
    warning(label, ": 病例或候选池为空"); return(data.frame())
  }
  n <- min(nrow(case_df), nrow(pool_df))
  if (n < nrow(case_df))
    warning(label, ": 候选池不足，仅抽取 ", n, " / ", nrow(case_df))
  
  pick <- sample(seq_len(nrow(pool_df)), n, replace = FALSE)
  case_use <- case_df[seq_len(n), , drop = FALSE]
  
  data.frame(
    disease_gene = case_use$hgnc_symbol,
    disease_tx   = case_use$ensembl_transcript_id,
    disease_cds  = case_use$cds_length,
    control_gene = pool_df$hgnc_symbol[pick],
    control_tx   = pool_df$ensembl_transcript_id[pick],
    control_cds  = pool_df$cds_length[pick],
    log2_ratio   = log2(pool_df$cds_length[pick] / case_use$cds_length),
    stringsAsFactors = FALSE)
}


report_match_quality <- function(pairs, n_case, label) {
  message("\n--- ", label, " 匹配质量 ---")
  message(sprintf("  配对数: %d / %d", nrow(pairs), n_case))
  if (nrow(pairs) == 0) return(invisible(NULL))
  message(sprintf("  |log2(对照/病例)| 中位数 %.3f, 最大 %.3f (最差 %.2f 倍)",
                  stats::median(abs(pairs$log2_ratio)), max(abs(pairs$log2_ratio)),
                  2^max(abs(pairs$log2_ratio))))
  ks <- suppressWarnings(stats::ks.test(pairs$disease_cds, pairs$control_cds))
  message(sprintf("  CDS 长度 KS 检验: p = %.4f  （匹配组期望不显著）", ks$p.value))
  invisible(ks)
}


## 一个 stratum（SNV 或 FS）的两套对照
build_both_controls <- function(disease_genes, pool_all, mart,
                                exclude_genes = character(0),
                                require_transcripts = NULL,
                                label = "SNV") {
  
  message("\n=== 构建对照: ", label, " ===")
  
  case_df <- fetch_canonical_cds(disease_genes, mart)
  message(sprintf("  病例中有 canonical CDS 的: %d / %d", nrow(case_df), length(disease_genes)))
  
  pool_df <- pool_all %>%
    filter(!hgnc_symbol %in% union(disease_genes, exclude_genes))
  message(sprintf("  候选池（OMIM AD 去除两个疾病列表后）: %d", nrow(pool_df)))
  
  # FS 层：NMDesc 区域通过 PTC_info 推导，不在 PTC_info 里的对照转录本会被
  # inner_join 静默丢弃，所以候选池必须先做限制
  if (!is.null(require_transcripts)) {
    n0 <- nrow(pool_df)
    pool_df <- pool_df %>% filter(ensembl_transcript_id %in% require_transcripts)
    message(sprintf("  限制到 PTC_info 转录本: %d -> %d", n0, nrow(pool_df)))
    if (nrow(pool_df) < nrow(case_df))
      warning(label, ": 候选池(", nrow(pool_df), ")小于病例数(", nrow(case_df),
              ")，部分病例无法匹配")
  }
  
  matched <- build_controls_matched(case_df, pool_df, label = paste(label, "matched"))
  random  <- build_controls_random(case_df, pool_df, label = paste(label, "random"))
  
  report_match_quality(matched, nrow(case_df), paste(label, "matched"))
  report_match_quality(random,  nrow(case_df), paste(label, "random"))
  
  list(matched = matched, random = random, case_df = case_df, pool_df = pool_df)
}


# =============================================================================
# 4. gene_all 组装（四组结构 + pair_id）
# =============================================================================

## 把配对表摊平成长表：每个 pair 两行（case + control），group 为四组之一
pairs_to_long <- function(pairs, stratum) {
  if (nrow(pairs) == 0) return(data.frame())
  pairs$pair_id <- sprintf("%s_%04d", stratum, seq_len(nrow(pairs)))
  
  bind_rows(
    data.frame(pair_id = pairs$pair_id, stratum = stratum, role = "case",
               hgnc_symbol = pairs$disease_gene,
               ensembl_transcript_id = pairs$disease_tx,
               cds_length = pairs$disease_cds, stringsAsFactors = FALSE),
    data.frame(pair_id = pairs$pair_id, stratum = stratum, role = "control",
               hgnc_symbol = pairs$control_gene,
               ensembl_transcript_id = pairs$control_tx,
               cds_length = pairs$control_cds, stringsAsFactors = FALSE)
  ) %>%
    mutate(group = factor(paste0(tolower(stratum), ifelse(role == "control", "_control", "")),
                          levels = GROUP_LEVELS))
}


## 补齐 CDS 序列、NMDesc 区域、uniprot
assemble_gene_all <- function(long_df, mart, PTC_combined, label = "") {
  
  message("\n=== 组装 gene_all", if (nchar(label)) paste0(" [", label, "]") else "", " ===")
  if (nrow(long_df) == 0) return(data.frame())
  
  tx <- unique(long_df$ensembl_transcript_id)
  
  message("  拉取 CDS 序列 ...")
  cds_df <- getBM_chunked(
    attributes = c("ensembl_transcript_id", "coding"),
    filters = "ensembl_transcript_id", values = tx, mart = mart,
    chunk = CONFIG$chunk_sequence) %>%
    filter(!is.na(coding), coding != "", coding != "Sequence unavailable") %>%
    distinct(ensembl_transcript_id, .keep_all = TRUE)
  
  message("  拉取外显子结构 ...")
  exon_df <- getBM_chunked(
    attributes = c("ensembl_transcript_id", "rank", "cds_start", "cds_end"),
    filters = "ensembl_transcript_id", values = tx, mart = mart,
    chunk = CONFIG$chunk_light) %>%
    filter(!is.na(rank), !is.na(cds_start), !is.na(cds_end))
  
  nmdesc_exon <- exon_df %>%
    group_by(ensembl_transcript_id) %>%
    summarise(exon_num         = max(rank, na.rm = TRUE),
              nmdesc_end       = max(cds_end),
              last_exon_length = cds_end[which.max(rank)] - cds_start[which.max(rank)] + 1,
              .groups = "drop") %>%
    mutate(nmdesc_start_exon = pmax(1, nmdesc_end - 50 - last_exon_length))
  
  message("  拉取 UniProt ID ...")
  uni_df <- getBM_chunked(
    attributes = c("ensembl_transcript_id", "uniprotswissprot"),
    filters = "ensembl_transcript_id", values = tx, mart = mart,
    chunk = CONFIG$chunk_light) %>%
    filter(uniprotswissprot != "") %>%
    distinct(ensembl_transcript_id, .keep_all = TRUE)
  
  out <- long_df %>%
    left_join(cds_df,      by = "ensembl_transcript_id") %>%
    left_join(nmdesc_exon, by = "ensembl_transcript_id") %>%
    left_join(uni_df,      by = "ensembl_transcript_id") %>%
    rename(uniprot = uniprotswissprot) %>%
    left_join(PTC_combined, by = c("ensembl_transcript_id" = "transcript"))
  
  # NMDesc 区域：FS 层用 PTC_info，SNV 层用外显子结构
  out <- out %>%
    mutate(
      NMD_region_start = ifelse(stratum == "FS" & !is.na(median_can_region_start),
                                median_can_region_start, nmdesc_start_exon),
      NMD_region_end   = ifelse(stratum == "FS" & !is.na(median_can_region_end),
                                median_can_region_end, nmdesc_end),
      NMDesc_region_length = NMD_region_end - NMD_region_start + 1,
      nmdesc_cds = ifelse(is.na(coding) | is.na(NMD_region_start), NA_character_,
                          substr(coding, NMD_region_start, NMD_region_end)),
      row_id = row_number()) %>%
    select(-any_of(c("median_can_region_start", "median_can_region_end",
                     "nmdesc_start_exon", "last_exon_length", "nmdesc_end")))
  
  # 只保留两侧信息都齐全的配对 —— 配对检验要求成对存在
  ok_pairs <- out %>%
    filter(!is.na(coding), !is.na(NMD_region_start), !is.na(NMDesc_region_length)) %>%
    count(pair_id) %>% filter(n == 2) %>% pull(pair_id)
  
  n_before <- length(unique(out$pair_id))
  out <- out %>% filter(pair_id %in% ok_pairs) %>%
    mutate(group = factor(as.character(group), levels = GROUP_LEVELS))
  
  message(sprintf("  完整配对: %d / %d", length(ok_pairs), n_before))
  print(table(out$group))
  
  out
}


# =============================================================================
# 5. 特征标注
# =============================================================================

annotate_sequence_features <- function(gene_all) {
  gc <- function(s) {
    if (is.na(s) || s == "") return(NA_real_)
    b <- strsplit(toupper(s), "")[[1]]
    round(sum(b %in% c("G", "C")) / length(b) * 100, 2)
  }
  rep_frac <- function(s) {
    if (is.na(s) || s == "") return(NA_real_)
    m <- str_extract_all(s, "([ATGC]{1,6})\\1+")[[1]]
    if (length(m) == 0) 0 else sum(nchar(m)) / nchar(s)
  }
  homo_frac <- function(s) {
    if (is.na(s) || s == "") return(NA_real_)
    m <- str_extract_all(s, "([ATGC])\\1{3,}")[[1]]
    if (length(m) == 0) 0 else sum(nchar(m)) / nchar(s)
  }
  
  gene_all %>%
    mutate(gc_content                  = vapply(coding,     gc,        numeric(1)),
           nmdesc_gc_content           = vapply(nmdesc_cds, gc,        numeric(1)),
           repeat_fraction             = vapply(coding,     rep_frac,  numeric(1)),
           nmdesc_repeat_fraction      = vapply(nmdesc_cds, rep_frac,  numeric(1)),
           homopolymer_fraction        = vapply(coding,     homo_frac, numeric(1)),
           nmdesc_homopolymer_fraction = vapply(nmdesc_cds, homo_frac, numeric(1)))
}


annotate_pfam <- function(gene_all, mart) {
  message("  PFAM 结构域 ...")
  tx <- unique(gene_all$ensembl_transcript_id)
  
  # 若 pfam / pfam_start / pfam_end 分属不同页，退回到两次请求后合并
  pfam <- tryCatch(
    getBM_chunked(attributes = c("ensembl_transcript_id", "pfam", "pfam_start", "pfam_end"),
                  filters = "ensembl_transcript_id", values = tx, mart = mart,
                  chunk = CONFIG$chunk_light),
    error = function(e) {
      if (!grepl("属性跨页", e$message)) stop(e)
      message("    pfam 属性跨页，改为分两次请求")
      a <- getBM_chunked(c("ensembl_transcript_id", "pfam"),
                         "ensembl_transcript_id", tx, mart, chunk = CONFIG$chunk_light)
      b <- getBM_chunked(c("ensembl_transcript_id", "pfam_start", "pfam_end"),
                         "ensembl_transcript_id", tx, mart, chunk = CONFIG$chunk_light)
      dplyr::inner_join(a, b, by = "ensembl_transcript_id")
    })
  
  if (nrow(pfam) == 0)
    return(gene_all %>% mutate(pfam_overlap_length = NA_real_,
                               pfam_overlap_flag = NA_integer_,
                               pfam_overlap_fraction = NA_real_))
  
  gp <- gene_all %>%
    select(row_id, ensembl_transcript_id, NMD_region_start, NMD_region_end,
           NMDesc_region_length) %>%
    left_join(pfam, by = "ensembl_transcript_id") %>%
    mutate(pf_start_bp = pfam_start * 3 - 2, pf_end_bp = pfam_end * 3,
           ov_start = pmax(NMD_region_start, pf_start_bp),
           ov_end   = pmin(NMD_region_end,   pf_end_bp),
           ov_valid = !is.na(ov_start) & !is.na(ov_end) & ov_start <= ov_end)
  
  summ <- bind_rows(lapply(split(seq_len(nrow(gp)), gp$row_id), function(idx) {
    d  <- gp[idx, , drop = FALSE]
    nm <- d$NMDesc_region_length[1]
    dv <- d[d$ov_valid %in% TRUE, , drop = FALSE]
    if (nrow(dv) == 0 || is.na(nm) || nm <= 0)
      return(data.frame(row_id = d$row_id[1], pfam_overlap_length = 0,
                        pfam_overlap_flag = 0L, pfam_overlap_fraction = 0))
    gr  <- GenomicRanges::GRanges("x", IRanges::IRanges(dv$ov_start, dv$ov_end))
    tot <- sum(IRanges::width(GenomicRanges::reduce(gr)))
    data.frame(row_id = d$row_id[1], pfam_overlap_length = tot,
               pfam_overlap_flag = as.integer(tot > 20),
               pfam_overlap_fraction = tot / nm)
  }))
  
  gene_all %>% left_join(summ, by = "row_id")
}


annotate_ppi_interface <- function(gene_all, ppi_file) {
  if (!file.exists(path.expand(ppi_file))) {
    message("  跳过 PPI 界面残基（文件不存在）")
    return(gene_all %>% mutate(ppi_overlap = NA_integer_))
  }
  message("  PPI 界面残基 ...")
  hp <- fread(path.expand(ppi_file), data.table = FALSE)
  
  to_num <- function(x) {
    x <- gsub("\\[|\\]", "", x)
    if (is.na(x) || x == "") return(numeric(0))
    as.numeric(trimws(unlist(strsplit(x, ","))))
  }
  
  gene_all$ppi_overlap <- 0L
  for (i in seq_len(nrow(gene_all))) {
    uid <- gene_all$uniprot[i]; aft <- gene_all$NMD_region_start[i]
    if (is.na(uid) || uid == "" || is.na(aft)) { gene_all$ppi_overlap[i] <- NA_integer_; next }
    r1 <- unlist(lapply(hp$interface_residues1[hp$uniprot1 == uid], to_num)) * 3 - 2
    r2 <- unlist(lapply(hp$interface_residues2[hp$uniprot2 == uid], to_num)) * 3 - 2
    if (length(r1) == 0 && length(r2) == 0) next
    if (any(r1 >= aft, na.rm = TRUE) || any(r2 >= aft, na.rm = TRUE))
      gene_all$ppi_overlap[i] <- 1L
  }
  gene_all
}


annotate_string_degree <- function(gene_all, score_threshold = 400) {
  if (!requireNamespace("STRINGdb", quietly = TRUE) ||
      !requireNamespace("igraph", quietly = TRUE)) {
    message("  跳过 STRING 度中心性（缺少 STRINGdb / igraph）")
    return(gene_all %>% mutate(degree_centrality = NA_real_))
  }
  message("  STRING 度中心性 ...")
  res <- tryCatch({
    sdb <- STRINGdb::STRINGdb$new(version = "11.5", species = 9606,
                                  score_threshold = score_threshold, input_directory = "")
    mp  <- sdb$map(as.data.frame(gene_all), "hgnc_symbol", removeUnmappedRows = TRUE)
    itx <- sdb$get_interactions(mp$STRING_id)
    g   <- igraph::simplify(igraph::graph_from_data_frame(
      data.frame(from = itx$from, to = itx$to), directed = FALSE))
    deg <- igraph::degree(g, normalized = TRUE)
    mp %>% select(hgnc_symbol, STRING_id) %>% distinct() %>%
      mutate(degree_centrality = unname(deg[STRING_id])) %>%
      group_by(hgnc_symbol) %>%
      summarise(degree_centrality = suppressWarnings(max(degree_centrality, na.rm = TRUE)),
                .groups = "drop") %>%
      mutate(degree_centrality = ifelse(is.finite(degree_centrality),
                                        degree_centrality, NA_real_))
  }, error = function(e) { message("  STRING 失败: ", e$message); NULL })
  
  if (is.null(res)) return(gene_all %>% mutate(degree_centrality = NA_real_))
  gene_all %>% left_join(res, by = "hgnc_symbol")
}


annotate_tau <- function(gene_all, gtex_path) {
  if (!file.exists(gtex_path)) {
    message("  跳过组织特异性 tau（GTEx 文件不存在）")
    return(gene_all %>% mutate(tau = NA_real_))
  }
  message("  GTEx tau ...")
  gtex <- read_tsv(gtex_path, skip = 2, show_col_types = FALSE)
  expr <- gtex %>%
    select(-tidyselect::any_of("Name")) %>%
    group_by(Description) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  m <- as.matrix(expr[, -1]); rownames(m) <- expr$Description
  m <- m[stats::complete.cases(m) & apply(m, 1, max) > 0, , drop = FALSE]
  
  tau_fun <- function(x) sum(1 - x / max(x)) / (length(x) - 1)
  tau_scores <- apply(m, 1, tau_fun)
  
  gene_all %>% mutate(tau = unname(tau_scores[hgnc_symbol]))
}


annotate_constraint <- function(gene_all, lof_path) {
  if (!file.exists(lof_path)) {
    message("  跳过 gnomAD 约束指标（文件不存在）")
    return(gene_all %>% mutate(pLI = NA_real_, oe_lof_upper = NA_real_))
  }
  message("  gnomAD pLI / LOEUF ...")
  lof <- read.delim(lof_path) %>%
    select(gene, pLI, oe_lof_upper) %>%
    distinct(gene, .keep_all = TRUE)
  gene_all %>% left_join(lof, by = c("hgnc_symbol" = "gene"))
}


# 注意：本函数与 features/functions/annotate_motif_flags.R 里的同名函数**不兼容** ——
# 那个是 7 参数版（多 mart / output_motif_csv / output_lcs_csv，并会写出两个 csv），
# 本函数是 4 参数版，只返回 gene_all，且带"补充文件缺失就填 NA 跳过"的保护。
# 为避免加载函数库时被覆盖，这里改名为 annotate_motif_flags_dbh。
annotate_motif_flags_dbh <- function(gene_all, path_touni, path_motif, path_lcs) {
  flags <- c("gene_protein_flag","gene_domains_flag","gene_slim_flag",
             "gene_morf_flag","gene_ptm_flag","gene_nls_flag","gene_LCS_flag")
  
  if (!all(file.exists(c(path_touni, path_motif, path_lcs)))) {
    message("  跳过 motif / LCS flags（补充文件不存在）")
    for (f in flags) gene_all[[f]] <- NA_integer_
    return(gene_all)
  }
  message("  motif / LCS flags ...")
  
  get_max <- function(x) {
    if (is.na(x) || x == "") return(NA_real_)
    n <- str_extract_all(x, "\\d+")[[1]]
    if (length(n) == 0) NA_real_ else max(as.numeric(n))
  }
  parse_max <- function(df) df %>% mutate(across(-1, ~ vapply(as.character(.x), get_max, numeric(1))))
  
  touni <- read_csv(path_touni, show_col_types = FALSE)
  mot   <- read_csv(path_motif, show_col_types = FALSE) %>% parse_max() %>%
    mutate(uniprot = touni$`Uniprot ID`[match(Protein, touni$Protein)])
  lcs   <- read_excel(path_lcs, sheet = "G") %>% parse_max() %>%
    mutate(uniprot = touni$`Uniprot ID`[match(Protein, touni$Protein)])
  
  add_flag <- function(ga, src, src_col, flag_name) {
    if (!src_col %in% names(src)) { ga[[flag_name]] <- NA_integer_; return(ga) }
    v <- src[[src_col]][match(ga$uniprot, src$uniprot)]
    ga[[flag_name]] <- as.integer(v * 3 >= ga$NMD_region_start)
    ga
  }
  
  gene_all %>%
    add_flag(mot, "Protein Features", "gene_protein_flag") %>%
    add_flag(mot, "Domains",          "gene_domains_flag") %>%
    add_flag(mot, "SLiMs",            "gene_slim_flag")    %>%
    add_flag(mot, "MORFs",            "gene_morf_flag")    %>%
    add_flag(mot, "PTMs",             "gene_ptm_flag")     %>%
    add_flag(mot, "NLSs/NESs",        "gene_nls_flag")     %>%
    add_flag(lcs, "LCSs",             "gene_LCS_flag")
}


annotate_all <- function(gene_all, mart) {
  message("\n=== 特征标注 ===")
  gene_all %>%
    annotate_sequence_features() %>%
    annotate_pfam(mart) %>%
    annotate_ppi_interface(CONFIG$ppi_file) %>%
    annotate_string_degree() %>%
    annotate_tau(CONFIG$gtex_path) %>%
    annotate_constraint(CONFIG$lof_metrics_path) %>%
    annotate_motif_flags_dbh(CONFIG$path_touni, CONFIG$path_motif, CONFIG$path_lcs)
}


# =============================================================================
# 6. 四组的描述性统计与分组输出
# =============================================================================

summarise_by_group <- function(gene_all, label) {
  feats <- intersect(c(CONTINUOUS_FEATURES, BINARY_FEATURES), names(gene_all))
  feats <- feats[vapply(feats, function(f) any(!is.na(gene_all[[f]])), logical(1))]
  
  desc <- gene_all %>%
    group_by(group) %>%
    summarise(n_genes       = n_distinct(hgnc_symbol),
              n_rows        = n(),
              cds_median    = stats::median(cds_length, na.rm = TRUE),
              cds_IQR       = stats::IQR(cds_length, na.rm = TRUE),
              nmdesc_median = stats::median(NMDesc_region_length, na.rm = TRUE),
              .groups = "drop")
  
  cat("\n=== 四组基本情况 [", label, "] ===\n", sep = "")
  print(as.data.frame(desc), row.names = FALSE, digits = 4)
  
  feat_tab <- gene_all %>%
    select(group, all_of(feats)) %>%
    pivot_longer(-group, names_to = "feature", values_to = "value") %>%
    filter(!is.na(value)) %>%
    group_by(feature, group) %>%
    summarise(n     = n(),
              value = if (dplyr::first(feature) %in% BINARY_FEATURES)
                mean(value >= 1) else stats::median(value),
              .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = c(n, value))
  
  cat("\n=== 各特征按组汇总（二分类=阳性比例，连续=中位数）[", label, "] ===\n", sep = "")
  print(as.data.frame(feat_tab), row.names = FALSE, digits = 4)
  
  write_csv(desc,     sprintf("group_summary_%s.csv", label))
  write_csv(feat_tab, sprintf("group_features_%s.csv", label))
  
  # 两个对照组来自同一候选池，可能有重叠 —— 不影响各自的配对检验有效性，
  # 但若之后要把两层合并分析，重叠基因会被计两次
  sc <- unique(gene_all$hgnc_symbol[gene_all$group == "snv_control"])
  fc <- unique(gene_all$hgnc_symbol[gene_all$group == "fs_control"])
  cat(sprintf("\nsnv_control 与 fs_control 重叠基因数: %d\n", length(intersect(sc, fc))))
  
  sd_ <- unique(gene_all$hgnc_symbol[gene_all$group == "snv"])
  fd_ <- unique(gene_all$hgnc_symbol[gene_all$group == "fs"])
  cat(sprintf("snv 与 fs 疾病基因重叠数: %d\n", length(intersect(sd_, fd_))))
  
  invisible(list(desc = desc, features = feat_tab))
}


write_by_group <- function(gene_all, label) {
  cat("\n--- 分组输出 [", label, "] ---\n", sep = "")
  for (g in GROUP_LEVELS) {
    sub <- gene_all %>% filter(group == g)
    if (nrow(sub) == 0) { message(sprintf("  %-12s 空", g)); next }
    fn <- sprintf("gene_all_%s_%s.csv", label, g)
    write_csv(sub, fn)
    cat(sprintf("  %-12s %3d 基因 -> %s\n", g, n_distinct(sub$hgnc_symbol), fn))
  }
}


## 显著性标签：p 值来自对应分析的检验结果（matched 用配对，random 用非配对）
fmt_p_label <- function(p, alpha = CONFIG$alpha) {
  if (length(p) == 0 || is.na(p)) return("p = NA")
  star <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < alpha) "*" else "ns"
  if (p < 0.001) paste0("p<0.001 ", star) else sprintf("p=%.3f %s", p, star)
}

## Brackets: snv vs snv_control (x=1,2) and fs vs fs_control (x=3,4)
add_group_brackets <- function(g, ymax, yrange, res_feat, alpha = CONFIG$alpha) {
  y_line <- ymax + 0.08 * yrange
  y_tick <- ymax + 0.05 * yrange
  y_text <- ymax + 0.14 * yrange
  
  get_p <- function(s) {
    r <- res_feat[res_feat$stratum == s, "p_adj_BH", drop = TRUE]
    if (length(r) == 0) NA_real_ else r[1]
  }
  
  for (bk in list(list(x1 = 1, x2 = 2, s = "SNV"),
                  list(x1 = 3, x2 = 4, s = "FS"))) {
    g <- g +
      annotate("segment", x = bk$x1, xend = bk$x2, y = y_line, yend = y_line, linewidth = 0.5) +
      annotate("segment", x = bk$x1, xend = bk$x1, y = y_tick, yend = y_line, linewidth = 0.5) +
      annotate("segment", x = bk$x2, xend = bk$x2, y = y_tick, yend = y_line, linewidth = 0.5) +
      annotate("text", x = (bk$x1 + bk$x2) / 2, y = y_text,
               label = fmt_p_label(get_p(bk$s), alpha),
               size = 2.9, fontface = "bold")
  }
  g + coord_cartesian(ylim = c(NA, ymax + 0.22 * yrange))
}


plot_four_groups <- function(gene_all, features, label, res = NULL,
                             ncol = 3, alpha = CONFIG$alpha) {
  features <- intersect(features, names(gene_all))
  features <- features[vapply(features, function(f) any(!is.na(gene_all[[f]])), logical(1))]
  if (length(features) == 0) return(invisible(NULL))
  
  plots <- lapply(features, function(f) {
    d <- gene_all %>% filter(!is.na(.data[[f]]))
    if (nrow(d) == 0) return(NULL)
    
    res_feat <- if (is.null(res)) data.frame() else res[res$feature == f, , drop = FALSE]
    
    if (f %in% BINARY_FEATURES) {
      s <- d %>% group_by(group) %>%
        summarise(pct = mean(.data[[f]] >= 1) * 100, n = n(),
                  pos = sum(.data[[f]] >= 1), .groups = "drop")
      g <- ggplot(s, aes(group, pct, fill = group)) +
        geom_col(width = 0.7, colour = "grey40") +
        geom_text(aes(label = paste0(pos, "/", n)), vjust = -0.4, size = 2.8) +
        scale_fill_manual(values = GROUP_COLORS, guide = "none") +
        labs(title = f, x = NULL, y = "Positive (%)") +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1),
              plot.title = element_text(face = "bold", size = 11))
      if (nrow(res_feat) > 0) g <- add_group_brackets(g, max(s$pct, na.rm = TRUE), 100, res_feat, alpha)
      g
      
    } else {
      v   <- d[[f]]
      # 95th percentile as the ceiling: long-tailed features would otherwise
      # push the brackets off-canvas and flatten the boxes
      yt  <- stats::quantile(v, 0.95, na.rm = TRUE)
      yb  <- stats::quantile(v, 0.05, na.rm = TRUE)
      rng <- max(yt - yb, .Machine$double.eps)
      
      g <- ggplot(d, aes(group, .data[[f]], fill = group)) +
        geom_boxplot(width = 0.6, outlier.size = 0.6) +
        scale_fill_manual(values = GROUP_COLORS, guide = "none") +
        labs(title = f, x = NULL, y = NULL) +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1),
              plot.title = element_text(face = "bold", size = 11))
      if (nrow(res_feat) > 0) g <- add_group_brackets(g, yt, rng, res_feat, alpha)
      g
    }
  })
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) return(invisible(NULL))
  
  test_note <- if (is.null(res)) "" else
    sprintf("; %s, BH-adjusted (%d tests per stratum)",
            paste(unique(res$test[!is.na(res$test)]), collapse = " / "),
            max(res$n_in_stratum, na.rm = TRUE))
  
  fig <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      title = sprintf("Feature distribution across four groups [%s]", label),
      subtitle = paste0("Brackets compare snv vs snv_control and fs vs fs_control", test_note),
      theme = theme(plot.title = element_text(face = "bold", size = 16),
                    plot.subtitle = element_text(size = 11)))
  
  fn <- sprintf("four_groups_%s.pdf", label)
  ggsave(fn, fig, width = 4.2 * ncol,
         height = 3.6 * ceiling(length(plots) / ncol), limitsize = FALSE)
  message("Saved: ", fn)
  invisible(fig)
}


# =============================================================================
# 7. 检验框架
# =============================================================================

make_pair_table <- function(gene_all, feature) {
  gene_all %>%
    select(pair_id, stratum, role, all_of(feature)) %>%
    pivot_wider(names_from = role, values_from = all_of(feature)) %>%
    filter(!is.na(case), !is.na(control))
}


## 配对检验（主分析）
##   二分类 -> 不一致配对上的精确二项检验（小样本）或 McNemar
##   连续   -> Wilcoxon 符号秩检验
test_paired <- function(gene_all, feature, stratum_now, exact_cutoff = 25) {
  
  d <- make_pair_table(gene_all, feature) %>% filter(stratum == stratum_now)
  n <- nrow(d)
  if (n < 3) return(NULL)
  
  if (feature %in% BINARY_FEATURES) {
    cs <- as.integer(d$case >= 1); ct <- as.integer(d$control >= 1)
    tab <- table(factor(cs, levels = c(0,1)), factor(ct, levels = c(0,1)))
    b <- unname(tab["0","1"]); cc <- unname(tab["1","0"])
    disc <- b + cc
    
    if (disc == 0)
      return(data.frame(feature = feature, stratum = stratum_now, n_pairs = n,
                        test = "none (无不一致配对)", statistic = NA_real_,
                        effect = NA_real_, p_value = NA_real_,
                        case_summary = mean(cs), control_summary = mean(ct),
                        stringsAsFactors = FALSE))
    
    if (disc < exact_cutoff) {
      p <- binom.test(cc, disc, 0.5)$p.value; tst <- "exact binomial (paired)"
    } else {
      p <- mcnemar.test(tab, correct = TRUE)$p.value; tst <- "McNemar"
    }
    data.frame(feature = feature, stratum = stratum_now, n_pairs = n,
               test = tst, statistic = disc,
               effect = if (b > 0) cc / b else NA_real_,
               p_value = p,
               case_summary = mean(cs), control_summary = mean(ct),
               stringsAsFactors = FALSE)
    
  } else {
    diff <- d$case - d$control
    if (all(diff == 0, na.rm = TRUE)) return(NULL)
    wt <- suppressWarnings(wilcox.test(d$case, d$control, paired = TRUE))
    data.frame(feature = feature, stratum = stratum_now, n_pairs = n,
               test = "Wilcoxon signed-rank", statistic = unname(wt$statistic),
               effect = stats::median(diff, na.rm = TRUE),
               p_value = wt$p.value,
               case_summary = stats::median(d$case, na.rm = TRUE),
               control_summary = stats::median(d$control, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }
}


## 非配对检验（敏感性分析）
##   二分类 -> Fisher 精确检验
##   连续   -> Wilcoxon 秩和检验
test_unpaired <- function(gene_all, feature, stratum_now) {
  
  d <- gene_all %>% filter(stratum == stratum_now, !is.na(.data[[feature]]))
  x <- d[[feature]][d$role == "case"]
  y <- d[[feature]][d$role == "control"]
  if (length(x) < 3 || length(y) < 3) return(NULL)
  
  if (feature %in% BINARY_FEATURES) {
    tab <- matrix(c(sum(x >= 1), sum(x < 1), sum(y >= 1), sum(y < 1)), nrow = 2)
    if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) return(NULL)
    ft <- fisher.test(tab)
    data.frame(feature = feature, stratum = stratum_now, n_pairs = NA_integer_,
               test = "Fisher exact", statistic = NA_real_,
               effect = unname(ft$estimate), p_value = ft$p.value,
               case_summary = mean(x >= 1), control_summary = mean(y >= 1),
               stringsAsFactors = FALSE)
  } else {
    wt <- suppressWarnings(wilcox.test(x, y))
    data.frame(feature = feature, stratum = stratum_now, n_pairs = NA_integer_,
               test = "Wilcoxon rank-sum", statistic = unname(wt$statistic),
               effect = stats::median(x, na.rm = TRUE) - stats::median(y, na.rm = TRUE),
               p_value = wt$p.value,
               case_summary = stats::median(x, na.rm = TRUE),
               control_summary = stats::median(y, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }
}


run_all_tests <- function(gene_all, paired = TRUE) {
  feats <- intersect(c(CONTINUOUS_FEATURES, BINARY_FEATURES), names(gene_all))
  feats <- feats[vapply(feats, function(f) any(!is.na(gene_all[[f]])), logical(1))]
  
  grid <- expand.grid(feature = feats, stratum = unique(gene_all$stratum),
                      stringsAsFactors = FALSE)
  
  bind_rows(lapply(seq_len(nrow(grid)), function(k) {
    if (paired) test_paired(gene_all,   grid$feature[k], grid$stratum[k])
    else        test_unpaired(gene_all, grid$feature[k], grid$stratum[k])
  }))
}


# =============================================================================
# 8. 多重检验校正
#
# 家族按 stratum 划分：SNV 与 FS 是两个独立的科学问题，各自校正；
# 同时给出跨 stratum 的全局校正供参考。
# 特征之间高度相关（GC 与 repeat、pLI 与 LOEUF、PFAM 三个视角），BH 在
# 正相关（PRDS）下有效；需要更保守时看 p_adj_holm。
# =============================================================================

correct_and_report <- function(res, label, alpha = CONFIG$alpha, output_csv = NULL) {
  if (nrow(res) == 0) { message(label, ": 无检验结果"); return(res) }
  
  res <- res %>%
    group_by(stratum) %>%
    mutate(n_in_stratum = n(),
           p_adj_BH     = p.adjust(p_value, method = "BH"),
           p_adj_holm   = p.adjust(p_value, method = "holm")) %>%
    ungroup() %>%
    mutate(p_adj_global = p.adjust(p_value, method = "BH"),
           sig = case_when(is.na(p_adj_BH)   ~ "NA",
                           p_adj_BH < 0.001  ~ "***",
                           p_adj_BH < 0.01   ~ "**",
                           p_adj_BH < alpha  ~ "*",
                           TRUE              ~ "ns")) %>%
    arrange(stratum, p_value)
  
  cat("\n==================== ", label, " ====================\n", sep = "")
  cat(sprintf("检验总数: %d\n", nrow(res)))
  print(res %>% group_by(stratum) %>%
          summarise(n_tests = n(),
                    raw_sig = sum(p_value  < alpha, na.rm = TRUE),
                    BH_sig  = sum(p_adj_BH < alpha, na.rm = TRUE),
                    .groups = "drop"))
  
  sig_rows <- res %>% filter(!is.na(p_adj_BH), p_adj_BH < alpha)
  if (nrow(sig_rows) > 0) {
    cat("\n--- BH FDR <", alpha, "的显著特征 ---\n")
    print(sig_rows %>%
            select(stratum, feature, test, case_summary, control_summary,
                   effect, p_value, p_adj_BH, sig) %>%
            as.data.frame(), row.names = FALSE, digits = 4)
  } else {
    cat("\n--- 无特征通过 BH FDR <", alpha, "---\n")
  }
  cat("=========================================================\n")
  
  if (!is.null(output_csv)) { write_csv(res, output_csv); message("已保存: ", output_csv) }
  res
}


# =============================================================================
# 9. 主流程
# =============================================================================

PTC_info <- read.csv(CONFIG$ptc_info)

PTC_combined <- PTC_info %>%
  group_by(transcript) %>%
  summarise(median_can_region_start = round(median(can_region_start, na.rm = TRUE)),
            median_can_region_end   = round(median(can_region_end,   na.rm = TRUE)),
            .groups = "drop")

snv_gene        <- read_gene_list(CONFIG$snv_gene_list,   "snv (DBH)")
fs_gene         <- read_gene_list(CONFIG$fs_gene_list,    "fs (Simes)")
omim_AD_symbols <- read_gene_list(CONFIG$omim_ad_symbols, "OMIM AD")

cat(sprintf("\nsnv n = %d | fs n = %d | 重叠 = %d\n",
            length(snv_gene), length(fs_gene), length(intersect(snv_gene, fs_gene))))

pool_all <- get_pool_cds(omim_AD_symbols, ensembl)

ctrl_snv <- build_both_controls(snv_gene, pool_all, ensembl,
                                exclude_genes = fs_gene, label = "SNV")
ctrl_fs  <- build_both_controls(fs_gene, pool_all, ensembl,
                                exclude_genes = snv_gene,
                                require_transcripts = unique(PTC_info$transcript),
                                label = "FS")

write_csv(ctrl_snv$matched, "pairs_snv_matched.csv")
write_csv(ctrl_fs$matched,  "pairs_fs_matched.csv")
write_csv(ctrl_snv$random,  "pairs_snv_random.csv")
write_csv(ctrl_fs$random,   "pairs_fs_random.csv")


# ---- 9a. 主分析：CDS 匹配对照 + 配对检验 -----------------------------------

long_matched <- bind_rows(pairs_to_long(ctrl_snv$matched, "SNV"),
                          pairs_to_long(ctrl_fs$matched,  "FS"))

gene_all_matched <- assemble_gene_all(long_matched, ensembl, PTC_combined, "matched") %>%
  annotate_all(ensembl)

write_csv(gene_all_matched, "gene_all_matched.csv")
write_by_group(gene_all_matched, "matched")
summarise_by_group(gene_all_matched, "matched")
res_matched <- run_all_tests(gene_all_matched, paired = TRUE) %>%
  correct_and_report("主分析：CDS 匹配对照，配对检验",
                     output_csv = "results_matched_paired.csv")

plot_four_groups(gene_all_matched, c(CONTINUOUS_FEATURES, BINARY_FEATURES),
                 "matched", res = res_matched)


# ---- 9b. 敏感性分析：随机对照 + 非配对检验 ---------------------------------

long_random <- bind_rows(pairs_to_long(ctrl_snv$random, "SNV"),
                         pairs_to_long(ctrl_fs$random,  "FS"))

gene_all_random <- assemble_gene_all(long_random, ensembl, PTC_combined, "random") %>%
  annotate_all(ensembl)

write_csv(gene_all_random, "gene_all_random.csv")
write_by_group(gene_all_random, "random")
summarise_by_group(gene_all_random, "random")
res_random <- run_all_tests(gene_all_random, paired = FALSE) %>%
  correct_and_report("敏感性分析：随机对照，非配对检验",
                     output_csv = "results_random_unpaired.csv")

plot_four_groups(gene_all_random, c(CONTINUOUS_FEATURES, BINARY_FEATURES),
                 "random", res = res_random)


# =============================================================================
# 10. 两套分析的对比
#
# 只在随机对照里显著的特征，很可能是 CDS 长度的代理而非真实信号；
# 两者都显著的最可信；只在匹配分析里显著的，说明 CDS 长度原本在掩盖信号
# （负混杂），也值得关注。
# =============================================================================

cmp <- full_join(
  res_matched %>% select(stratum, feature, p_matched = p_value,
                         padj_matched = p_adj_BH, eff_matched = effect),
  res_random  %>% select(stratum, feature, p_random = p_value,
                         padj_random = p_adj_BH, eff_random = effect),
  by = c("stratum", "feature")) %>%
  mutate(verdict = case_when(
    !is.na(padj_matched) & padj_matched < CONFIG$alpha &
      !is.na(padj_random) & padj_random < CONFIG$alpha ~ "两者均显著（最可信）",
    !is.na(padj_matched) & padj_matched < CONFIG$alpha ~ "仅匹配分析显著（CDS 长度曾掩盖信号）",
    !is.na(padj_random)  & padj_random  < CONFIG$alpha ~ "仅随机对照显著（可能是 CDS 长度的代理）",
    TRUE ~ "均不显著")) %>%
  arrange(stratum, padj_matched)

cat("\n\n==================== 主分析 vs 敏感性分析 ====================\n")
print(as.data.frame(cmp), row.names = FALSE, digits = 4)
cat("==============================================================\n")
print(table(cmp$verdict))
write_csv(cmp, "results_comparison.csv")


# =============================================================================
# 11. 配对分布图（仅主分析的显著特征）
# =============================================================================

sig_feats <- res_matched %>% filter(!is.na(p_adj_BH), p_adj_BH < CONFIG$alpha)

if (nrow(sig_feats) > 0) {
  plots <- lapply(seq_len(nrow(sig_feats)), function(k) {
    f <- sig_feats$feature[k]; s <- sig_feats$stratum[k]
    d <- make_pair_table(gene_all_matched, f) %>% filter(stratum == s)
    if (nrow(d) == 0) return(NULL)
    
    dl <- d %>% pivot_longer(c(case, control), names_to = "role", values_to = "value")
    
    ggplot(dl, aes(x = role, y = value, group = pair_id)) +
      geom_line(alpha = 0.25, colour = "grey50") +
      geom_point(aes(colour = role), size = 1.8, alpha = 0.8) +
      scale_colour_manual(values = c(case = "#D07A3A", control = "#8ABA67"), guide = "none") +
      labs(title = sprintf("%s [%s]", f, s),
           subtitle = sprintf("%s, BH p = %.3g", sig_feats$test[k], sig_feats$p_adj_BH[k]),
           x = NULL, y = f) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 12))
  })
  plots <- Filter(Negate(is.null), plots)
  
  if (length(plots) > 0) {
    fig <- patchwork::wrap_plots(plots, ncol = min(3, length(plots))) +
      patchwork::plot_annotation(
        title = "Paired distribution of significant features (matched analysis)",
        subtitle = "Each line is one CDS-length-matched case-control pair",
        theme = theme(plot.title = element_text(face = "bold", size = 16)))
    ggsave("paired_significant_features.pdf", fig,
           width = 4.2 * min(3, length(plots)),
           height = 3.6 * ceiling(length(plots) / 3), limitsize = FALSE)
    message("已保存: paired_significant_features.pdf")
  }
} else {
  message("主分析无显著特征，跳过绘图")
}