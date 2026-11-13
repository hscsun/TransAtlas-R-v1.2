#!/usr/bin/env Rscript
# Regenerate the fixed-template results report without re-running any analysis.
#
# Reads the tables under results/<data>/ and provenance/config_used.yml, then
# fills the same SUMMARY_TEMPLATE used inside the pipeline.
#
# Usage:
#   Rscript scripts/build_summary_report.R
#   Rscript scripts/build_summary_report.R --results results/RNAseqSampleData
#   Rscript scripts/build_summary_report.R --out results/RNAseqSampleData/report \
#       --file results_report_review.docx

suppressPackageStartupMessages({
  ok <- requireNamespace("yaml", quietly = TRUE)
})
if (!ok) stop("Package 'yaml' is required: install.packages('yaml')", call. = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1L]]
}

root <- get_flag("root", NULL)
if (is.null(root)) {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  root <- if (length(script_arg)) {
    normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), ".."),
                  mustWork = FALSE)
  } else normalizePath(".", mustWork = FALSE)
}

results_dir <- get_flag("results", file.path(root, "results", "RNAseqSampleData"))
results_dir <- normalizePath(results_dir, mustWork = TRUE)
config_file <- get_flag("config",
                        file.path(results_dir, "provenance", "config_used.yml"))
if (!file.exists(config_file)) {
  stop("Configuration snapshot not found: ", config_file,
       "\nRun the pipeline once, or pass --config explicitly.", call. = FALSE)
}
out_dir <- get_flag("out", file.path(results_dir, "report"))
file_name <- get_flag("file", "results_report.docx")

source(file.path(root, "pipeline", "R", "helpers.R"))
source(file.path(root, "pipeline", "R", "summary_report.R"))

cfg <- yaml::read_yaml(config_file)
if (is.null(cfg$report)) cfg$report <- list()
cfg$report$summary_report <- TRUE

render_summary_report(cfg, results_dir, out_dir, file_name = file_name)
