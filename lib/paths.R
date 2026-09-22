# =============================================================================
# paths.R -- data file locator
# -----------------------------------------------------------------------------
# Why this file exists:
#   Edit only this file to change data locations.
#
# Usage (at top of script):
#   source("gene level_v5/lib/paths.R")
#   d <- data_file("genemap2.txt")        # automatically locates it
#   x <- read.csv(data_file("BM_info.csv"))
#
# Reports missing files clearly; never fails silently.
# =============================================================================

# Single data root: every input is resolved from one directory inside the repo,
# <repo>/data level_v5. This file lives at <repo>/gene level_v5/lib/paths.R, so the
# repo root is three levels up from it; the location is found by scanning the
# source() call stack for this file's own path, which works however deep the
# sourcing chain is. REPO (set by run_analysis.R) and the Desktop layout are
# fallbacks for the case where the stack carries no file name.
.find_data_root <- function() {
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of) && nzchar(of) && grepl("paths\\.R$", of)) {
      return(file.path(dirname(dirname(dirname(normalizePath(of)))), "data level_v5"))
    }
  }
  if (exists("REPO", envir = globalenv())) {
    return(file.path(get("REPO", envir = globalenv()), "data level_v5"))
  }
  path.expand("~/Desktop/NMDesc/repo_v5/data level_v5")
}
DATA_ROOTS <- c(main = .find_data_root())
rm(.find_data_root)
# Overridable via env var: export NMDESC_DATA=/other/location
if (nzchar(Sys.getenv("NMDESC_DATA"))) {
  DATA_ROOTS <- c(main = Sys.getenv("NMDESC_DATA"))
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
  # A miss can mean the index predates the file, so rescan once before failing.
  # Only misses pay for this; a hit returns from the cached index.
  if (is.null(hit) || !file.exists(hit)) {
    idx <- .build_index(force = TRUE)
    hit <- idx[[b]]
  }
  if (!is.null(hit) && file.exists(hit)) return(hit)
  if (must) {
    stop(sprintf(paste0(
      "Data file not found: %s\n",
      "  Searched: %s\n",
      "  Fix: place the file in the data root, ",
      "or export NMDESC_DATA=<directory>"),
      b, paste(DATA_ROOTS[dir.exists(DATA_ROOTS)], collapse = ", ")),
      call. = FALSE)
  }
  NA_character_
}

#' Output directory (results written here)
out_dir <- function(sub = "") {
  # Defaults to <data root>/output, which keeps it inside the search path: the
  # files one step writes are found by the next step's data_file().
  base <- if (nzchar(Sys.getenv("NMDESC_OUT"))) Sys.getenv("NMDESC_OUT") else
          file.path(unname(DATA_ROOTS[[1]]), "output")
  d <- if (nzchar(sub)) file.path(base, sub) else base
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
out_file <- function(name, sub = "") file.path(out_dir(sub), basename(name))

#' Alias for out_dir(), usable as a default argument
#' A function whose own parameter is named `out_dir` cannot default to
#' `out_dir()`: the default expression would resolve to that parameter's own
#' promise. Such formals use `results_dir()`.
results_dir <- function(sub = "") out_dir(sub)

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
  # There is one root now. Older call sites ask for it by its former name
  # ("clinvar", "idr", ...); they all resolve to the single root.
  p <- if (name %in% names(DATA_ROOTS)) unname(DATA_ROOTS[[name]]) else unname(DATA_ROOTS[[1]])
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

#' Join a freshly computed result table against a stored one and report drift
#' Both comparison scripts source this file, so the check lives here.
#' @param new   the table this run produced
#' @param stored_name  file name of a previous run's table, located by data_file()
#' @param keys  columns identifying a row in both tables
#' @param sig_col  column holding the significance mark, compared for flips
#' Writes compare_<stored_name> to out_dir() and returns it, or NULL when the
#' stored table is not on a data root.
compare_to_stored <- function(new, stored_name, keys, sig_col = "sig") {
  p <- data_file(stored_name, must = FALSE)
  if (is.na(p)) {
    message("  no stored ", stored_name, " on a data root, comparison skipped")
    return(invisible(NULL))
  }
  old <- utils::read.csv(p)
  if (!all(keys %in% names(old)) || !all(keys %in% names(new))) {
    message("  ", stored_name, " does not carry ", paste(keys, collapse = " + "),
            ", comparison skipped")
    return(invisible(NULL))
  }
  shared <- setdiff(intersect(names(old), names(new)), keys)
  num <- shared[vapply(shared, function(c) is.numeric(old[[c]]) && is.numeric(new[[c]]),
                       logical(1))]
  m <- merge(old[c(keys, shared)], new[c(keys, shared)], by = keys,
             suffixes = c("_stored", "_now"), all = TRUE)

  cat(sprintf("\n--- %s vs this run: %d rows matched, %d stored-only, %d new-only\n",
              stored_name, sum(stats::complete.cases(m[paste0(num[1], c("_stored", "_now"))])),
              sum(is.na(m[[paste0(num[1], "_now")]])),
              sum(is.na(m[[paste0(num[1], "_stored")]]))))
  for (c in num) {
    a <- m[[paste0(c, "_stored")]]; b <- m[[paste0(c, "_now")]]
    ok <- !is.na(a) & !is.na(b)
    if (!any(ok)) next
    cat(sprintf("  %-18s max |diff| %-12.4g identical %d/%d\n", c,
                max(abs(a[ok] - b[ok])), sum(a[ok] == b[ok]), sum(ok)))
  }
  if (sig_col %in% shared) {
    s1 <- m[[paste0(sig_col, "_stored")]]; s2 <- m[[paste0(sig_col, "_now")]]
    flips <- which(!is.na(s1) & !is.na(s2) & s1 != s2)
    cat(sprintf("  significance flips: %d\n", length(flips)))
    if (length(flips))
      print(m[flips, c(keys, paste0(sig_col, c("_stored", "_now")))], row.names = FALSE)
  }
  utils::write.csv(m, out_file(paste0("compare_", stored_name)), row.names = FALSE)
  invisible(m)
}

#' Placeholder for a file known missing locally
#' Valid syntax; errors only when called, stating what is missing
missing_file <- function(name, was = NULL) {
  stop(sprintf(paste0(
    "This file is not present locally: %s\n",
    if (!is.null(was)) sprintf("  Original path: %s\n", was) else "",
    "  place it in the data root once found (then refresh_data_index())"),
    name), call. = FALSE)
}
