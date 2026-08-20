# =============================================================================
# paths.R —— 数据文件定位层（唯一需要改的地方）
# -----------------------------------------------------------------------------
# 为什么要这个文件：
#   原来脚本里有 270 处硬编码路径，同一个 genemap2.txt 出现过 5 种不同写法
#   （~/Downloads/OMIM/、~/Desktop/autism/data/、/Users/jxu14/... 等）。
#   换机器、挪文件夹，就得逐个改。现在只改这里。
#
# 用法（在脚本开头）：
#   source("gene level_v3/lib/paths.R")
#   d <- data_file("genemap2.txt")        # 自动找到它在哪
#   x <- read.csv(data_file("BM_info.csv"))
#
# 找不到时会明确报错并告诉你去哪放，不会静默返回坏路径。
# =============================================================================

# 数据根目录：按优先级搜索。改这里就能整体搬家。
DATA_ROOTS <- c(
  clinvar    = path.expand("~/Desktop/clinvar"),
  new_clinvar= path.expand("~/Desktop/new_clinvar"),
  gnomAD_snv = path.expand("~/Desktop/gnomAD_snv"),
  idr        = path.expand("~/Desktop/idr"),
  syn        = path.expand("~/Desktop/syn"),
  legacy     = path.expand("~/Desktop/桌面 - q的MacBook Pro")
)
# 环境变量可覆盖：export NMDESC_DATA=/别的位置
if (nzchar(Sys.getenv("NMDESC_DATA"))) {
  DATA_ROOTS <- c(main = Sys.getenv("NMDESC_DATA"), DATA_ROOTS)
}

.pathcache <- new.env(parent = emptyenv())

# 索引缓存文件（避免每次都扫 2 万多个文件）
# 依次尝试：数据目录 -> 家目录 -> 临时目录，取第一个可写的
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
.CACHE_TTL <- 7 * 24 * 3600     # 7 天后自动重建

# 建索引：先读磁盘缓存，过期或不存在才重扫
.build_index <- function(verbose = FALSE, force = FALSE) {
  if (!force && !is.null(.pathcache$idx)) return(invisible(.pathcache$idx))

  CF <- .cache_path()
  if (!force && file.exists(CF)) {
    age <- as.numeric(difftime(Sys.time(), file.mtime(CF), units = "secs"))
    if (age < .CACHE_TTL) {
      cc <- tryCatch(readRDS(CF), error = function(e) NULL)
      if (!is.null(cc) && identical(cc$roots, DATA_ROOTS)) {
        .pathcache$idx <- cc$idx
        if (verbose) message(sprintf("  用缓存索引 (%d 个文件名, %.0f 小时前建)",
                                     length(cc$idx), age / 3600))
        return(invisible(cc$idx))
      }
    }
  }

  if (verbose) message("  扫描数据目录建索引（首次约 30 秒）...")
  idx <- list()
  for (nm in names(DATA_ROOTS)) {
    root <- DATA_ROOTS[[nm]]
    if (!dir.exists(root)) next
    fs <- list.files(root, recursive = TRUE, full.names = TRUE,
                     all.files = FALSE, no.. = TRUE)
    fs <- fs[!dir.exists(fs)]
    for (f in fs) {
      b <- basename(f)
      if (is.null(idx[[b]])) idx[[b]] <- f    # 先到先得：DATA_ROOTS 顺序即优先级
    }
    if (verbose) message(sprintf("  %s: %d 个文件", nm, length(fs)))
  }
  .pathcache$idx <- idx
  tryCatch({ saveRDS(list(idx = idx, roots = DATA_ROOTS), .cache_path())
             if (verbose) message("  索引已缓存: ", .cache_path())
           }, error = function(e) if (verbose) message("  缓存写入失败（不影响使用）"))
  invisible(idx)
}

#' 数据文件挪动过就调一次，重建索引
refresh_data_index <- function() invisible(.build_index(verbose = TRUE, force = TRUE))

#' 按文件名定位数据文件
#' @param name 文件名（只要 basename，不用管在哪个子目录）
#' @param must 找不到时是否报错（FALSE 则返回 NA）
data_file <- function(name, must = TRUE) {
  idx <- .build_index()
  b <- basename(name)
  hit <- idx[[b]]
  if (!is.null(hit) && file.exists(hit)) return(hit)
  if (must) {
    stop(sprintf(paste0(
      "找不到数据文件: %s\n",
      "  已搜索: %s\n",
      "  解决办法：把该文件放进 ~/Desktop/clinvar/，",
      "或 export NMDESC_DATA=<存放目录>"),
      b, paste(DATA_ROOTS[dir.exists(DATA_ROOTS)], collapse = ", ")),
      call. = FALSE)
  }
  NA_character_
}

#' 输出目录（结果写这里，不再散落）
out_dir <- function(sub = "") {
  base <- if (nzchar(Sys.getenv("NMDESC_OUT"))) Sys.getenv("NMDESC_OUT") else
          file.path(path.expand("~/Desktop/clinvar"), "output")
  d <- if (nzchar(sub)) file.path(base, sub) else base
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
out_file <- function(name, sub = "") file.path(out_dir(sub), basename(name))

#' 体检：报告哪些文件在、哪些缺
#' 取某个已命名数据根目录的实际路径
#' 用于少数确实需要"目录"而非"文件"的场景（例如原代码里的 setwd 或 file.path 拼接）。
#' 新写的代码建议直接用 data_file("文件名")，不依赖目录结构。
#' @param name DATA_ROOTS 里的名字，如 "clinvar"；缺省取第一个存在的根目录
data_root <- function(name = NULL) {
  if (is.null(name)) {
    ex <- DATA_ROOTS[dir.exists(DATA_ROOTS)]
    if (!length(ex)) stop("DATA_ROOTS 里没有一个目录存在 —— 请检查 paths.R", call. = FALSE)
    return(unname(ex[1]))
  }
  if (!name %in% names(DATA_ROOTS))
    stop(sprintf("DATA_ROOTS 里没有名为 '%s' 的根目录。可用：%s",
                 name, paste(names(DATA_ROOTS), collapse = ", ")), call. = FALSE)
  p <- unname(DATA_ROOTS[[name]])
  if (!dir.exists(p))
    warning(sprintf("数据根目录 '%s' 不存在：%s", name, p), call. = FALSE)
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

#' 已知本机缺失的文件占位
#' 语法合法、可被解析；真正执行到才报错，并说清楚缺什么、原来在哪
missing_file <- function(name, was = NULL) {
  stop(sprintf(paste0(
    "此文件本机没有: %s\n",
    if (!is.null(was)) sprintf("  原路径（旧机器）: %s\n", was) else "",
    "  找到后放进 ~/Desktop/clinvar/ 即可（然后 refresh_data_index()）"),
    name), call. = FALSE)
}
