#!/usr/bin/env Rscript
# Redraw figures from the cached ggplot objects, without re-running any analysis.
#
# Every ggplot figure saved through save_plot() is also written to
# cache/figures/<figure file name>.rds. This script loads those objects and
# re-saves the same plots at whatever size the registry in
# pipeline/R/figure_sizes.R currently specifies.
#
# Usage (from the project root):
#   Rscript scripts/redraw_figures.R                    # every cached figure
#   Rscript scripts/redraw_figures.R s2_go_composite    # one key or several
#   Rscript scripts/redraw_figures.R GSEA_GO_BP_composite_Treatment_vs_Control
#   Rscript scripts/redraw_figures.R --list             # keys, sizes, paths
#
# Figures drawn with pheatmap or base graphics have no cached object. Change
# their size in the registry and rerun only the section that owns them, e.g.
#   Rscript pipeline/run_analysis.R --sections 1 --no-report

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(".", winslash = "/")

source(file.path(project_root, "pipeline", "R", "figure_sizes.R"))

# A ggraph or patchwork object needs its package loaded before it can be drawn.
for (pkg in c("ggplot2", "scales", "ggraph", "patchwork", "ggrepel")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

cache_dir <- Sys.getenv("CLAWSGO_FIGURE_CACHE", "cache/figures")
if (!grepl("^([A-Za-z]:)?/", cache_dir)) cache_dir <- file.path(project_root, cache_dir)
files <- if (dir.exists(cache_dir)) list.files(cache_dir, "\\.rds$", full.names = TRUE) else character(0)
if (!length(files)) {
  stop("No cached figures in ", cache_dir,
       ". Run the analysis (or one section) once with the current code first.",
       call. = FALSE)
}

# index.csv (written by cache_figure()) lists key / file / size arguments, so
# --list and single-figure redraws do not have to deserialize every plot.
parse_args_txt <- function(txt) {
  txt <- if (is.na(txt)) "" else txt
  if (!nzchar(txt)) return(list())
  parts <- strsplit(txt, ";", fixed = TRUE)[[1]]
  out <- lapply(parts, function(p) {
    kv <- strsplit(p, "=", fixed = TRUE)[[1]]
    suppressWarnings(as.numeric(strsplit(kv[2], "|", fixed = TRUE)[[1]]))
  })
  names(out) <- vapply(strsplit(parts, "=", fixed = TRUE), `[`, character(1), 1)
  out
}

index_file <- file.path(cache_dir, "index.csv")
index <- if (file.exists(index_file)) {
  try(utils::read.csv(index_file, stringsAsFactors = FALSE), silent = TRUE)
} else NULL
if (inherits(index, "try-error")) index <- NULL

# Fall back to reading the plot files when no index exists (e.g. a cache
# produced by an older version).
load_records <- function(rds_files) {
  records <- lapply(rds_files, readRDS)
  records[vapply(records, function(r) is.list(r) && !is.null(r$plot), logical(1))]
}

if ("--list" %in% args) {
  cat(sprintf("%-22s %-58s %s\n", "KEY", "FIGURE", "CURRENT SIZE (in)"))
  if (!is.null(index) && nrow(index)) {
    for (i in seq_len(nrow(index))) {
      sz <- do.call(fig_size, c(list(index$key[i]), parse_args_txt(index$size_args[i])))
      cat(sprintf("%-22s %-58s %.2f x %.2f\n", index$key[i],
                  tools::file_path_sans_ext(index$rds[i]), sz[1], sz[2]))
    }
  } else {
    records <- load_records(files)
    for (rec in records) {
      sz <- do.call(fig_size, c(list(rec$key), rec$size_args))
      cat(sprintf("%-22s %-58s %.2f x %.2f\n", rec$key,
                  tools::file_path_sans_ext(basename(rec$path)), sz[1], sz[2]))
    }
  }
  quit(save = "no")
}

wanted <- setdiff(args, "--list")
# Resolve the requested keys to cache files first, then load only those.
target_files <- if (!is.null(index) && nrow(index)) {
  stems <- tools::file_path_sans_ext(index$rds)
  hit <- if (length(wanted)) (index$key %in% wanted) | (stems %in% wanted) else rep(TRUE, nrow(index))
  if (!any(hit)) {
    stop("No cached figure matched: ", paste(wanted, collapse = ", "),
         ". Use --list to see the available keys.", call. = FALSE)
  }
  file.path(cache_dir, index$rds[hit])
} else {
  files
}
records <- load_records(target_files)
keys <- vapply(records, function(r) r$key, character(1))
stems <- tools::file_path_sans_ext(basename(vapply(records, function(r) r$path, character(1))))
keep <- if (length(wanted)) (keys %in% wanted) | (stems %in% wanted) else rep(TRUE, length(records))
if (!any(keep)) {
  stop("No cached figure matched: ", paste(wanted, collapse = ", "),
       ". Use --list to see the available keys.", call. = FALSE)
}

message("Redrawing ", sum(keep), " figure(s) from ", cache_dir)
for (i in which(keep)) {
  rec <- records[[i]]
  sz <- do.call(fig_size, c(list(rec$key), rec$size_args))
  ok <- tryCatch({
    ggplot2::ggsave(rec$path, plot = rec$plot, width = sz[1], height = sz[2],
                    units = "in", dpi = 320, bg = "white")
    TRUE
  }, error = function(e) {
    warning("Could not redraw ", basename(rec$path), ": ", conditionMessage(e))
    FALSE
  })
  if (ok) {
    message(sprintf("  %-58s %.2f x %.2f in", basename(rec$path), sz[1], sz[2]))
  }
}
message("Done. Figures were written in place; rebuild the atlas with ",
        "`Rscript scripts/build_figure_atlas.R` if you need the PDF updated.")
