# ---- 1) look for where is the raw data ------------------------------------------------------
.find_repo <- function() {
  p <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NA)
  if (!is.na(p) && dir.exists(file.path(p, "gene level_v4"))) return(p)
  for (cand in c("~/Desktop/NMDesc/repo_cleaned", getwd(),
                 dirname(getwd()))) {
    cand <- path.expand(cand)
    if (dir.exists(file.path(cand, "gene level_v4"))) return(cand)
  }
  stop("can't locate",
       call. = FALSE)
}
REPO <- if (exists("REPO")) path.expand(REPO) else .find_repo()
# Absolute, because run_one() moves the working directory.
REPO <- normalizePath(REPO, mustWork = TRUE)

# The analysis scripts look for paths.R and helpers under "* level_v3" names and
# under a top-level "lib". These links point both at the v4 directories, so
# data_file() is available instead of each script's bare-relative-path fallback.
.link_compat <- function(repo) {
  pairs <- lapply(c("gene", "variant", "protein"),
                  function(v) c(paste0(v, " level_v4"), paste0(v, " level_v3")))
  pairs <- c(pairs, list(c("gene level_v4/lib", "lib")))
  for (p in pairs) {
    src <- file.path(repo, p[1]); dst <- file.path(repo, p[2])
    if (dir.exists(src) && !file.exists(dst))
      tryCatch(file.symlink(src, dst),
               warning = function(w) cat("  no link for", p[2], "\n"))
  }
}
.link_compat(REPO)

# paths.R supplies data_file(), data_root() and out_dir() to every script.
source(file.path(REPO, "gene level_v4/lib/paths.R"))

# ---- 2) check necessary packages --------------------------------------------------------------
need_main  <- c("dplyr", "tidyr", "readr", "ggplot2", "purrr", "patchwork", "lme4")
need_bayes <- c("brms")    

miss_main  <- need_main [!vapply(need_main,  requireNamespace, logical(1), quietly = TRUE)]
miss_bayes <- need_bayes[!vapply(need_bayes, requireNamespace, logical(1), quietly = TRUE)]

if (length(miss_main))  cat("missing packages:", paste(miss_main, collapse = ", "), "\n")
RUN_BAYES <- length(miss_bayes) == 0
cat("Bayesian sensitivity model:", if (RUN_BAYES) "on" else "off (brms absent)", "\n")

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
if (is.na(F_MATCHED))
  cat("  gene_all_matched.csv absent; scripts that need it will report it\n")

# ---- 4) output dir ------------------------------------------------------------
OUTDIR <- file.path(REPO, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
cat("\n output loc:", sub(path.expand("~"), "~", OUTDIR), "\n")

# ---- 5) run analysis --------------------------------------------------------------
# script is a path relative to the repository root. Each script resolves its own
# helpers through paths relative to its own directory, so the working directory
# moves there for the duration of the call.
# Stages are given by file name, located inside the repository, so a script that
# moves between directories still runs. Several names may be supplied; the first
# one present wins, which covers a rename that has not reached every copy yet.
.locate <- function(names) {
  for (n in names) {
    if (grepl("/", n) && file.exists(file.path(REPO, n))) return(file.path(REPO, n))
    hits <- list.files(REPO, pattern = paste0("^", basename(n), "$"),
                       recursive = TRUE, full.names = TRUE)
    hits <- hits[!grepl("/backup/", hits)]
    if (length(hits)) return(hits[1])
  }
  NA_character_
}

run_one <- function(script, label, args = NULL) {
  path <- .locate(script)
  if (is.na(path)) { cat("\nskipping", label, "-- not found:",
                         paste(script, collapse = " / "), "\n")
                     return(invisible()) }
  cat("\n", strrep("=", 70), "\n", label, "\n", strrep("=", 70), "\n", sep = "")
  wd <- getwd(); old <- commandArgs
  res <- tryCatch({
    setwd(dirname(path))
    if (!is.null(args))
      assign("commandArgs", function(trailingOnly = FALSE) args, envir = globalenv())
    # Global environment, because the stage scripts call source() themselves and
    # that lands in globalenv: a separate environment splits their definitions
    # from the objects those definitions read.
    sys.source(basename(path), envir = globalenv())
    "OK"
  }, error = function(e) paste("failed:", conditionMessage(e)))
  setwd(wd)
  assign("commandArgs", old, envir = globalenv())
  cat("\n[", label, "]", res, "\n")
}

# Acquisition first, then comparison. The second name on each line is the
# pre-rename file, used when a checkout still carries it.
STAGES <- list(
  list(c("gene_get_main.R", "main.R"),
       "1/4  Disease gene lists: NMDesc enrichment, AD-restricted candidates"),
  list(c("variant_get_main.R"),
       "2/4  Variant sets: ClinVar P/LP and gnomAD controls on the same transcripts"),
  list(c("gene_compare_main.R", "gene_main_dbh.R"),
       "3/4  Gene level: matched by CDS length, feature comparison within pairs"),
  list(c("variant_compare_main.R", "variant_main_DBH.R"),
       "4/4  Variant level: mixed-effect model, Bayesian sensitivity model"))

for (s in STAGES) run_one(s[[1]], s[[2]])

# ---- 6) list outputs --------------------------------------------------------
# Scripts write to REPO/figures, to out_dir() from paths.R, to an "out" directory
# beside themselves, or into their own directory when the path is relative.
cat("\n", strrep("=", 70), "\ngenerated files\n", strrep("=", 70), "\n", sep = "")
beside <- list.dirs(REPO, recursive = TRUE, full.names = TRUE)
beside <- beside[basename(beside) == "out" & !grepl("/backup/", beside)]
beside <- normalizePath(beside, mustWork = FALSE)
stage_dirs <- unlist(lapply(STAGES, function(s) {
  p <- .locate(s[[1]])
  # normalized, so a v3 compatibility link and its v4 target count once
  if (is.na(p)) NULL else normalizePath(dirname(p), mustWork = FALSE)
}))
dirs <- unique(c(OUTDIR, tryCatch(out_dir(), error = function(e) NULL),
                 beside, stage_dirs))
dirs <- dirs[dir.exists(dirs)]
for (d in dirs) {
  fs <- list.files(d, pattern = "[.](pdf|png|csv)$", full.names = TRUE)
  cat("\n ", sub(path.expand("~"), "~", d), "\n")
  if (length(fs)) {
    for (f in fs[order(file.mtime(fs), decreasing = TRUE)])
      cat(sprintf("  %8.0f KB  %s\n", file.size(f) / 1024, basename(f)))
  } else cat("   no files\n")
}
