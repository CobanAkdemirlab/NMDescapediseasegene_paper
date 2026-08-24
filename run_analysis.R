# ---- 1) look for where is the raw data ------------------------------------------------------
.find_repo <- function() {
  p <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NA)
  if (!is.na(p) && dir.exists(file.path(p, "gene level_v3"))) return(p)
  for (cand in c("~/Desktop/NMDesc/repo_cleaned", getwd(),
                 dirname(getwd()))) {
    cand <- path.expand(cand)
    if (dir.exists(file.path(cand, "gene level_v3"))) return(cand)
  }
  stop("can't locate",
       call. = FALSE)
}
REPO <- if (exists("REPO")) path.expand(REPO) else .find_repo()

# ---- 2) check necessary packages --------------------------------------------------------------
need_main  <- c("dplyr", "tidyr", "readr", "ggplot2", "purrr", "patchwork", "lme4")
need_bayes <- c("brms")    

miss_main  <- need_main [!vapply(need_main,  requireNamespace, logical(1), quietly = TRUE)]
miss_bayes <- need_bayes[!vapply(need_bayes, requireNamespace, logical(1), quietly = TRUE)]

# ---- 3) find raw file ------------------------------------------------------------
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

cat("\n=== data file ===\n")
for (nm in c(gene_matched = F_MATCHED, gene_random = F_RANDOM,
             analysis_table = F_TABLE, idr_metrics = F_IDR)) NULL
for (k in c("F_MATCHED", "F_RANDOM", "F_TABLE", "F_IDR")) {
  v <- get(k)
  cat(sprintf("  %-14s %s\n", k,
              if (is.na(v)) "miss" else sub(path.expand("~"), "~", v)))
}
if (is.na(F_MATCHED)) stop("can't find gene_all_matched.csv", call. = FALSE)

# ---- 4) output dir ------------------------------------------------------------
OUTDIR <- file.path(REPO, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
cat("\n output loc:", sub(path.expand("~"), "~", OUTDIR), "\n")

# ---- 5) run analysis --------------------------------------------------------------
run_one <- function(script, args, label) {
  path <- file.path(REPO, "scripts", script)
  if (!file.exists(path)) { cat("\n跳过", label, "—— 找不到", script, "\n"); return(invisible()) }
  cat("\n", strrep("=", 70), "\n", label, "\n", strrep("=", 70), "\n", sep = "")
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
        "Gene level：matched by CDS length")


  run_one("fig_variant_level.R",
          c(F_TABLE, F_IDR, OUTDIR, if (RUN_BAYES) "bayes" else "nobayes"),
          "Variant level：mixed effect（main）+ Bayesian（sensitivity）")

# ---- 6) main ----------------------------------------------------------------
cat("\n", strrep("=", 70), "\ngenerated files\n", strrep("=", 70), "\n", sep = "")
fs <- list.files(OUTDIR, pattern = "[.](pdf|png|csv)$", full.names = TRUE)
if (length(fs)) {
  for (f in fs[order(file.mtime(fs), decreasing = TRUE)])
    cat(sprintf("  %8.0f KB  %s\n", file.size(f) / 1024, basename(f)))
  cat(sprintf('  system("open \\"%s\\"")\n', OUTDIR))
} else cat("  missing files \n")
