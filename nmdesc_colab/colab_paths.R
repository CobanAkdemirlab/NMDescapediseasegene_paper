# colab_paths.R -- path mapping layer for running this repo on Google Colab.
# The v4 scripts read data through absolute paths on the original author's
# machine. DATA_DIR points at one folder holding those inputs; run_script()
# rewrites the prefixes on a working copy and sources that copy.

DATA_DIR <- Sys.getenv("NMDESC_DATA", "/content/drive/MyDrive/NMDesc_data")
REPO_DIR <- Sys.getenv("NMDESC_REPO", "/content/NMDescapediseasegene_paper")
WORK_DIR <- Sys.getenv("NMDESC_WORK", "/content/work")
dir.create(WORK_DIR, showWarnings = FALSE, recursive = TRUE)

# Longest prefixes first, so nested ones are not shadowed.
PREFIXES <- c(
  "~/Desktop/NMDescapediseasegene_paper-main/new_NMDesc/data",
  "/Users/jxu14/Desktop/NMDescapediseasegene_paper-main",
  "~/Desktop/new_clinvar/raw_data",
  "~/Desktop/new_clinvar/snv_list/list4",
  "~/Desktop/new_clinvar",
  "/Users/jxu14/Desktop/enrich",
  "/Users/jxu14/Desktop/autism",
  "/Users/qkelly/Desktop/clinvar",
  "~/Desktop/clinvar",
  "~/Downloads",
  "~/Desktop",
  "/Users/jxu14/Desktop",
  "/Users/qkelly/Desktop"
)

map_path <- function(x) {
  for (p in PREFIXES) {
    if (startsWith(x, p)) return(file.path(DATA_DIR, sub("^/+", "", substring(x, nchar(p) + 1))))
  }
  x
}

# Every data reference a script makes, resolved against DATA_DIR.
script_inputs <- function(script) {
  txt <- readLines(file.path(REPO_DIR, script), warn = FALSE)
  m <- regmatches(txt, gregexpr('"[^"]*\\.(csv|txt|tsv|rds|fasta|fa|xlsx)"', txt, ignore.case = TRUE))
  refs <- unique(gsub('"', "", unlist(m)))
  if (!length(refs)) return(data.frame(ref = character(), resolved = character(), exists = logical()))
  res <- vapply(refs, function(r)
    if (startsWith(r, "/") || startsWith(r, "~")) map_path(r) else file.path(DATA_DIR, basename(r)),
    character(1))
  data.frame(ref = refs, resolved = unname(res), exists = file.exists(unname(res)),
             row.names = NULL, stringsAsFactors = FALSE)
}

# Report what a script needs before spending time on it.
check_script <- function(script) {
  d <- script_inputs(script)
  cat(sprintf("%s: %d data refs, %d present, %d missing\n",
              script, nrow(d), sum(d$exists), sum(!d$exists)))
  if (any(!d$exists)) print(d[!d$exists, c("ref", "resolved")], row.names = FALSE)
  invisible(d)
}

# Rewrite absolute prefixes on a copy, then source the copy.
run_script <- function(script, echo = TRUE) {
  src <- file.path(REPO_DIR, script)
  txt <- readLines(src, warn = FALSE)
  for (p in PREFIXES) txt <- gsub(p, DATA_DIR, txt, fixed = TRUE)
  txt <- gsub("setwd\\(", "# setwd(", txt)
  dst <- file.path(WORK_DIR, gsub("[/ ]", "_", script))
  writeLines(txt, dst)
  source(dst, echo = echo, max.deparse.length = 200)
  invisible(dst)
}

# Sources a script where it lives, with the working directory set to its own
# directory, so relative helper lookups inside the script resolve. Use this for
# entry points; run_script() suits scripts that carry hardcoded data paths.
run_in_place <- function(script, echo = FALSE) {
  path <- file.path(REPO_DIR, script)
  wd <- getwd(); on.exit(setwd(wd), add = TRUE)
  setwd(dirname(path))
  source(basename(path), echo = echo, max.deparse.length = 200)
  invisible(path)
}

cat("colab_paths.R loaded\n  DATA_DIR:", DATA_DIR, "\n  REPO_DIR:", REPO_DIR, "\n")
