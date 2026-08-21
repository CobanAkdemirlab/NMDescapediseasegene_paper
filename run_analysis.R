# =============================================================================


# ---- 1) 定位仓库根目录 ------------------------------------------------------
.find_repo <- function() {
  p <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NA)
  if (!is.na(p) && dir.exists(file.path(p, "gene level_v3"))) return(p)
  for (cand in c("~/Desktop/NMDesc/repo_cleaned", getwd(),
                 dirname(getwd()))) {
    cand <- path.expand(cand)
    if (dir.exists(file.path(cand, "gene level_v3"))) return(cand)
  }
  stop("找不到仓库目录。请手动设置：REPO <- \"~/Desktop/NMDesc/repo_cleaned\"",
       call. = FALSE)
}
REPO <- if (exists("REPO")) path.expand(REPO) else .find_repo()
cat("仓库:", REPO, "\n")

# ---- 2) 检查包 --------------------------------------------------------------
need_main  <- c("dplyr", "tidyr", "readr", "ggplot2", "purrr", "patchwork", "lme4")
need_bayes <- c("brms")     # 只有贝叶斯敏感性分析需要；主分析不受影响

miss_main  <- need_main [!vapply(need_main,  requireNamespace, logical(1), quietly = TRUE)]
miss_bayes <- need_bayes[!vapply(need_bayes, requireNamespace, logical(1), quietly = TRUE)]

if (length(miss_main)) {
  cat("\n主分析缺这些包，先在 Console 跑：\n")
  cat(sprintf('  install.packages(c(%s))\n',
              paste0('"', miss_main, '"', collapse = ", ")))
  stop("缺包，装完再运行", call. = FALSE)
}
RUN_BAYES <- length(miss_bayes) == 0
if (!RUN_BAYES) {
  cat("\n注意：缺", paste(miss_bayes, collapse = "、"),
      "—— 贝叶斯敏感性分析会跳过。\n")
  cat("  想跑的话先装（brms 要编译，约 10-20 分钟）：\n")
  cat(sprintf('  install.packages(c(%s))\n',
              paste0('"', miss_bayes, '"', collapse = ", ")))
  cat("  主分析（mixed effect）不受影响，继续。\n")
}

# ---- 3) 定位数据 ------------------------------------------------------------
# gene level 用你已配好的对；variant level 用整合表 + IDR 指标
DATA_DIRS <- path.expand(c("~/Desktop/clinvar", "~/Desktop/clinvar/derived"))
find_data <- function(fn) {
  for (d in DATA_DIRS) {
    p <- file.path(d, fn)
    if (file.exists(p)) return(p)
  }
  NA_character_
}
F_MATCHED <- find_data("gene_all_matched.csv")
F_RANDOM  <- find_data("gene_all_random.csv")
F_TABLE   <- find_data("analysis_table.csv")
F_IDR     <- find_data("idr_metrics.csv")

cat("\n=== 数据文件 ===\n")
for (nm in c(gene_matched = F_MATCHED, gene_random = F_RANDOM,
             analysis_table = F_TABLE, idr_metrics = F_IDR)) NULL
for (k in c("F_MATCHED", "F_RANDOM", "F_TABLE", "F_IDR")) {
  v <- get(k)
  cat(sprintf("  %-14s %s\n", k,
              if (is.na(v)) "缺失" else sub(path.expand("~"), "~", v)))
}
if (is.na(F_MATCHED)) stop("找不到 gene_all_matched.csv", call. = FALSE)

# ---- 4) 输出目录 ------------------------------------------------------------
OUTDIR <- file.path(REPO, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
cat("\n输出目录:", sub(path.expand("~"), "~", OUTDIR), "\n")

# ---- 5) 跑分析 --------------------------------------------------------------
run_one <- function(script, args, label) {
  path <- file.path(REPO, "scripts", script)
  if (!file.exists(path)) { cat("\n跳过", label, "—— 找不到", script, "\n"); return(invisible()) }
  cat("\n", strrep("=", 70), "\n", label, "\n", strrep("=", 70), "\n", sep = "")
  # 用子进程跑，避免脚本之间的变量互相污染
  old <- commandArgs
  res <- tryCatch({
    assign("commandArgs", function(trailingOnly = FALSE) args, envir = globalenv())
    sys.source(path, envir = new.env(parent = globalenv()))
    "OK"
  }, error = function(e) paste("失败:", conditionMessage(e)))
  assign("commandArgs", old, envir = globalenv())
  cat("\n[", label, "]", res, "\n")
}

run_one("fig_gene_level.R",
        c(F_MATCHED, if (is.na(F_RANDOM)) "NA" else F_RANDOM, OUTDIR),
        "Gene level：基因层面特征，CDS 长度配对")

if (!is.na(F_TABLE) && !is.na(F_IDR)) {
  run_one("fig_variant_level.R",
          c(F_TABLE, F_IDR, OUTDIR, if (RUN_BAYES) "bayes" else "nobayes"),
          "Variant level：mixed effect（主）+ Bayesian（敏感性）")
} else {
  cat("\n跳过 variant level —— 缺 analysis_table.csv 或 idr_metrics.csv\n")
}

# ---- 6) 汇总 ----------------------------------------------------------------
cat("\n", strrep("=", 70), "\n生成的文件\n", strrep("=", 70), "\n", sep = "")
fs <- list.files(OUTDIR, pattern = "[.](pdf|png|csv)$", full.names = TRUE)
if (length(fs)) {
  for (f in fs[order(file.mtime(fs), decreasing = TRUE)])
    cat(sprintf("  %8.0f KB  %s\n", file.size(f) / 1024, basename(f)))
  cat("\n在 RStudio 里打开：\n")
  cat(sprintf('  system("open \\"%s\\"")\n', OUTDIR))
} else cat("  没有生成文件 —— 上面应该有报错信息\n")
