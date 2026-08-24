# =============================================================================
# paths.R -- data file locator
# -----------------------------------------------------------------------------
# Why this file exists:
#   Edit only this file to change data locations.
#
# Usage (at top of script):
#   source("gene level_v3/lib/paths.R")
#   d <- data_file("genemap2.txt")        # automatically locates it
#   x <- read.csv(data_file("BM_info.csv"))
#
# Reports missing files clearly; never fails silently.
# =============================================================================

# Data root directories, searched in priority order.
DATA_ROOTS <- c(
  clinvar    = path.expand("~/Desktop/clinvar"),
  new_clinvar= path.expand("~/Desktop/new_clinvar"),
  gnomAD_snv = path.expand("~/Desktop/gnomAD_snv"),
  idr        = path.expand("~/Desktop/idr"),
  syn        = path.expand("~/Desktop/syn"),
  legacy     = path.expand("~/Desktop/桌面 - q的MacBook Pro")
)
# Overridable via env var: export NMDESC_DATA=/other/location
if (nzchar(Sys.getenv("NMDESC_DATA"))) {
  DATA_ROOTS <- c(main = Sys.getenv("NMDESC_DATA"), DATA_ROOTS)
}

.pathcache <- new.env(parent = emptyenv())

# Index cache file (avoids rescanning 20,000+ files)
# Tries data dir -> home dir -> temp dir; uses first writable
.cache_path <- function() {
  cands <- c(file.path(DATA_ROOTS[[1]], ".nmdesc_path_index.rds"),
             file.path(path.expand("~"), ".nmdesc_path_index.rds"),
             file.path(tempdir(), "nmdesc_path_index.rds"))
  for (p in cands) {
    d <- dirname(p)
    if (dir.exists(d) && file.access(d, 2) == 0) return(p)
  }
  cands[length(cands)]
}
.CACHE_TTL <- 7 * 24 * 3600     # rebuilds automatically after 7 days

# Build index: read disk cache, rescan if expired or missing
.build_index <- function(verbose = FALSE, force = FALSE) {
  if (!force && !is.null(.pathcache$idx)) return(invisible(.pathcache$idx))

  CF <- .cache_path()
  if (!force && file.exists(CF)) {
    age <- as.numeric(difftime(Sys.time(), file.mtime(CF), units = "secs"))
    if (age < .CACHE_TTL) {
      cc <- tryCatch(readRDS(CF), error = function(e) NULL)
      if (!is.null(cc) && identical(cc$roots, DATA_ROOTS)) {
        .pathcache$idx <- cc$idx
        if (verbose) message(sprintf("  Using cached index (%d files, built %.0f hours ago)",
                                     length(cc$idx), age / 3600))
        return(invisible(cc$idx))
      }
    }
  }

  if (verbose) message("  Scanning data directories to build index (~30s first time)...")
  idx <- list()
  for (nm in names(DATA_ROOTS)) {
    root <- DATA_ROOTS[[nm]]
    if (!dir.exists(root)) next
    fs <- list.files(root, recursive = TRUE, full.names = TRUE,
                     all.files = FALSE, no.. = TRUE)
    fs <- fs[!dir.exists(fs)]
    for (f in fs) {
      b <- basename(f)
      if (is.null(idx[[b]])) idx[[b]] <- f    # first match wins: DATA_ROOTS order sets priority
    }
    if (verbose) message(sprintf("  %s: %d files", nm, length(fs)))
  }
  .pathcache$idx <- idx
  tryCatch({ saveRDS(list(idx = idx, roots = DATA_ROOTS), .cache_path())
             if (verbose) message("  Index cached: ", .cache_path())
           }, error = function(e) if (verbose) message("  Cache write failed (does not affect usage)"))
  invisible(idx)
}

#' Rebuild index after data files move
refresh_data_index <- function() invisible(.build_index(verbose = TRUE, force = TRUE))

#' Locate data file by filename
#' @param name File name (basename only; subdirectory does not matter)
#' @param must Whether to error if not found (FALSE returns NA)
data_file <- function(name, must = TRUE) {
  idx <- .build_index()
  b <- basename(name)
  hit <- idx[[b]]
  if (!is.null(hit) && file.exists(hit)) return(hit)
  if (must) {
    stop(sprintf(paste0(
      "Data file not found: %s\n",
      "  Searched: %s\n",
      "  Fix: place the file in ~/Desktop/clinvar/, ",
      "or export NMDESC_DATA=<directory>"),
      b, paste(DATA_ROOTS[dir.exists(DATA_ROOTS)], collapse = ", ")),
      call. = FALSE)
  }
  NA_character_
}

#' Output directory (results written here)
out_dir <- function(sub = "") {
  base <- if (nzchar(Sys.getenv("NMDESC_OUT"))) Sys.getenv("NMDESC_OUT") else
          file.path(path.expand("~/Desktop/clinvar"), "output")
  d <- if (nzchar(sub)) file.path(base, sub) else base
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
out_file <- function(name, sub = "") file.path(out_dir(sub), basename(name))

#' Health check: report which files exist, which are missing
#' Get actual path of a named data root directory
#' For cases needing a directory rather than a file (e.g. setwd or file.path)
#' @param name Name in DATA_ROOTS, e.g. "clinvar"; defaults to first existing root
data_root <- function(name = NULL) {
  if (is.null(name)) {
    ex <- DATA_ROOTS[dir.exists(DATA_ROOTS)]
    if (!length(ex)) stop("No directory in DATA_ROOTS exists -- check paths.R", call. = FALSE)
    return(unname(ex[1]))
  }
  if (!name %in% names(DATA_ROOTS))
    stop(sprintf("No root directory named '%s' in DATA_ROOTS. Available: %s",
                 name, paste(names(DATA_ROOTS), collapse = ", ")), call. = FALSE)
  p <- unname(DATA_ROOTS[[name]])
  if (!dir.exists(p))
    warning(sprintf("Data root directory '%s' does not exist: %s", name, p), call. = FALSE)
  p
}

check_data <- function(names) {
  res <- data.frame(file = names, found = FALSE, path = NA_character_,
                    stringsAsFactors = FALSE)
  for (i in seq_along(names)) {
    p <- data_file(names[i], must = FALSE)
    if (!is.na(p)) { res$found[i] <- TRUE; res$path[i] <- p }
  }
  res
}

#' Placeholder for a file known missing locally
#' Valid syntax; errors only when called, stating what is missing
missing_file <- function(name, was = NULL) {
  stop(sprintf(paste0(
    "This file is not present locally: %s\n",
    if (!is.null(was)) sprintf("  Original path: %s\n", was) else "",
    "  place it in ~/Desktop/clinvar/ once found (then refresh_data_index())"),
    name), call. = FALSE)
}
