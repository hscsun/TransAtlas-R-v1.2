# Command-line entry point.
#
#   Rscript pipeline/run_analysis.R --counts input/NIH3T3.xlsx \
#           --samples input/NIH3T3_metadata.csv --species mouse --contrast FBS,Control
#
#   Rscript pipeline/run_analysis.R --help     (or -h) prints this block
#
#   Rscript pipeline/run_analysis.R --config my_settings.yml --counts data.xlsx \
#           --species human --sections 1,2,3
#
# Every flag overrides one argument of run_all(); anything left out keeps the
# built-in default. No configuration file is read (--config remains available
# as an optional YAML overlay).
#
# Input / design   --counts --samples|--metadata --species mouse|human|rat
#                  --group-column --contrast A,B --contrast-name LABEL --gene-id-column
#                  --gene-symbol-column --annotation-column
#                  --counts-sheet --count-columns --expression-matrix
#                  --gene-length-column --gene-annotation (GTF/GFF fallback)
# DE               --log2fc --min-fpkm --alpha --min-count --min-samples
# TFs              --top-tfs --confidence A,B --tf-min-targets
#                  --tf-ora-min-overlap --tf-coherence-min-targets --tf
#                  --tf-dorothea-dir
# Display          --max-genes --top-n-labels --up-cutoff --volcano-thresholds
#                  --expression-log10-threshold (legacy alias: --expression-log2-threshold)
# Enrichment       --low-cutoff (GSEA only; 0 = no lower bound)
#                  --go-pvalue --go-qvalue   (--gene-set-lfc-threshold: deprecated, use --log2fc)
#                  --go-bp-pvalue --go-bp-qvalue --go-simplify --gsea-min-size
#                  --gsea-max-size --gsea-pvalue --gsea-top-sets --no-msigdb
#                  --msigdb-category --no-pathview --no-pathview-download
#                  --pathview-max-pathways
# Regulatory       --network-layouts --network-max-labels --cluster-min-jaccard
#                  --cluster-max --activity-methods
#                  --activity-min-size --activity-alpha --no-gsva
#                  --gsva-collection --gsva-min-size --gsva-max-size
# QC               --no-qc --qc-plots pca,correlation,...
# Output           --out --sections 1,2,3 --no-report --report-title
#                  --report-author --report-language --pdf --no-summary-report
#                  --kegg-cache-dir
#                  --kegg-timeout --no-package-check --no-install-packages
#                  --stop-on-error --keep-existing --seed --config
#
# TF-target regulons come from the DoRothEA data of the `dorothea` package, so
# only the count matrix and the sample metadata have to be supplied.
#
# The interactive path (file dialog) is pipeline/run_all.R used from an R
# session; this script takes its inputs from the command line.

args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(name) {
  i <- match(paste0("--", name), args)
  if (!is.na(i) && i < length(args) && !grepl("^--", args[i + 1L])) return(args[i + 1L])
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) return(sub(paste0("^--", name, "="), "", hit[1L]))
  NULL
}
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else "pipeline/run_analysis.R"
project_root <- normalizePath(file.path(dirname(normalizePath(script_path, mustWork = FALSE)), ".."),
                              winslash = "/", mustWork = TRUE)

if (any(args %in% c("--help", "-h"))) {
  header <- readLines(normalizePath(script_path, mustWork = FALSE), warn = FALSE)
  header <- header[seq_len(which(!grepl("^#", header))[1L] - 1L)]
  cat(paste(sub("^# ?", "", header), collapse = "\n"), "\n", sep = "")
  quit(save = "no", status = 0)
}

source(file.path(project_root, "pipeline", "R", "figure_sizes.R"))
source(file.path(project_root, "pipeline", "R", "helpers.R"))
source(file.path(project_root, "pipeline", "R", "context.R"))
source(file.path(project_root, "pipeline", "R", "kegg_module.R"))
source(file.path(project_root, "pipeline", "R", "deseq2_module.R"))
source(file.path(project_root, "pipeline", "R", "provenance.R"))
source(file.path(project_root, "pipeline", "R", "qc_module.R"))
source(file.path(project_root, "pipeline", "R", "run_pipeline.R"))

has_flag <- function(name) any(args == paste0("--", name))
list_flag <- function(name) {
  v <- get_flag(name)
  if (is.null(v)) NULL else as.character(trimws(strsplit(v, ",")[[1L]]))
}
num_flag <- function(name) {
  v <- get_flag(name)
  if (is.null(v)) NULL else as.numeric(v)
}
int_flag <- function(name) {
  v <- get_flag(name)
  if (is.null(v)) NULL else as.integer(v)
}

# Only flags that were actually given become arguments, so build_run_config()
# applies its built-in defaults to everything else.
overrides <- local({
  out <- list()
  add <- function(key, value) if (!is.null(value)) out[[key]] <<- value

  add("counts_file", get_flag("counts"))
  add("samples", get_flag("samples") %||% get_flag("metadata"))
  add("species", get_flag("species"))
  add("group_column", get_flag("group-column"))
  add("contrast", list_flag("contrast"))
  add("contrast_name", get_flag("contrast-name"))
  add("log2fc", num_flag("log2fc"))
  add("min_fpkm", num_flag("min-fpkm"))
  add("alpha", num_flag("alpha"))
  add("min_count", int_flag("min-count"))
  add("min_samples", int_flag("min-samples"))
  add("gene_id_column", get_flag("gene-id-column"))
  add("gene_symbol_column", get_flag("gene-symbol-column"))
  add("annotation_column", get_flag("annotation-column"))
  add("counts_sheet", get_flag("counts-sheet"))
  add("count_columns", list_flag("count-columns"))
  add("expression_matrix_file", get_flag("expression-matrix"))
  add("expression_matrix_unit", get_flag("expression-matrix-unit"))
  add("gene_annotation", get_flag("gene-annotation"))
  add("gene_length_column", get_flag("gene-length-column"))
  add("top_tfs", int_flag("top-tfs"))
  add("confidence", if (is.null(get_flag("confidence"))) NULL else
    toupper(list_flag("confidence")))

  add("tf_min_targets", int_flag("tf-min-targets"))
  add("tf_ora_min_overlap", int_flag("tf-ora-min-overlap"))
  add("tf_coherence_min_targets", int_flag("tf-coherence-min-targets"))
  add("tf_file", get_flag("tf"))
  add("tf_dorothea_dir", get_flag("tf-dorothea-dir"))
  add("max_genes", int_flag("max-genes"))
  add("top_n_labels", int_flag("top-n-labels"))
  add("volcano_label_expression_n", int_flag("volcano-label-expression-n"))
  add("volcano_label_effect_n", int_flag("volcano-label-effect-n"))
  add("up_cutoff", num_flag("up-cutoff"))
  add("low_cutoff", num_flag("low-cutoff"))
  add("log2fc_cap", num_flag("log2fc-cap"))
  add("volcano_lfc_thresholds", if (is.null(get_flag("volcano-thresholds"))) NULL else
    suppressWarnings(as.numeric(list_flag("volcano-thresholds"))))
  add("expression_log10_threshold", num_flag("expression-log10-threshold") %||%
    num_flag("expression-log2-threshold"))

  # Two cut-offs of older versions are now one each: the DEG cut (--log2fc) and
  # the expression floor (--min-fpkm). The old flags are accepted so archived
  # command lines still start, but they no longer change the analysis.
  if (!is.null(get_flag("gene-set-lfc-threshold"))) {
    warning("--gene-set-lfc-threshold is deprecated and ignored; set --log2fc instead.",
            call. = FALSE)
  }
  if (!is.null(get_flag("tf-min-fpkm"))) {
    warning("--tf-min-fpkm is deprecated and ignored; set --min-fpkm instead.", call. = FALSE)
  }

  add("go_pvalue_cutoff", num_flag("go-pvalue"))
  add("go_qvalue_cutoff", num_flag("go-qvalue"))
  add("go_bp_pvalue_cutoff", num_flag("go-bp-pvalue"))
  add("go_bp_qvalue_cutoff", num_flag("go-bp-qvalue"))
  add("go_simplify_cutoff", num_flag("go-simplify"))
  add("gsea_min_size", int_flag("gsea-min-size"))
  add("gsea_max_size", int_flag("gsea-max-size"))
  add("gsea_pvalue_cutoff", num_flag("gsea-pvalue"))
  add("gsea_top_sets", int_flag("gsea-top-sets"))
  if (has_flag("no-msigdb")) add("run_msigdb", FALSE)
  add("msigdb_category", get_flag("msigdb-category"))
  if (has_flag("no-pathview")) add("run_pathview", FALSE)
  if (has_flag("no-pathview-download")) add("pathview_download", FALSE)
  add("pathview_max_pathways", int_flag("pathview-max-pathways"))
  add("network_layouts", list_flag("network-layouts"))
  add("network_max_labels", int_flag("network-max-labels"))
  add("cluster_min_jaccard", num_flag("cluster-min-jaccard"))
  add("cluster_max", int_flag("cluster-max"))
  add("activity_methods", list_flag("activity-methods"))
  add("activity_min_size", int_flag("activity-min-size"))
  add("activity_alpha", num_flag("activity-alpha"))
  if (has_flag("no-gsva")) add("run_gsva", FALSE)
  add("gsva_collection", get_flag("gsva-collection"))
  add("gsva_min_size", int_flag("gsva-min-size"))
  add("gsva_max_size", int_flag("gsva-max-size"))
  if (has_flag("no-ppi")) add("run_ppi", FALSE)
  add("ppi_score_threshold", num_flag("ppi-score"))
  add("ppi_max_nodes", int_flag("ppi-max-nodes"))
  add("ppi_label_n", int_flag("ppi-label-n"))
  add("ppi_max_modules", int_flag("ppi-max-modules"))
  add("ppi_module_annotation_alpha", num_flag("ppi-module-annotation-alpha"))
  add("ppi_cache_dir", get_flag("ppi-cache-dir"))
  add("ppi_timeout", num_flag("ppi-timeout"))
  if (has_flag("no-qc")) add("qc", FALSE)
  add("qc_plots", list_flag("qc-plots"))
  add("sections", list_flag("sections"))
  if (has_flag("no-report")) add("report", FALSE)
  if (has_flag("no-summary-report")) add("summary_report", FALSE)
  add("report_title", get_flag("report-title"))
  add("report_author", get_flag("report-author"))
  add("report_language", get_flag("report-language"))
  if (has_flag("pdf")) add("generate_pdf", TRUE)
  add("kegg_cache_dir", get_flag("kegg-cache-dir"))
  add("kegg_download_timeout", num_flag("kegg-timeout"))
  if (has_flag("no-package-check")) add("check_packages", FALSE)
  if (has_flag("no-install-packages")) add("install_missing_packages", FALSE)
  if (has_flag("stop-on-error")) add("continue_on_error", FALSE)
  add("output_dir", get_flag("out"))
  if (has_flag("keep-existing")) add("overwrite", FALSE)
  add("seed", int_flag("seed"))
  out
})

cfg <- do.call(build_run_config, c(
  overrides,
  list(project_root = project_root, config_file = get_flag("config"))))

message("  output   : ", cfg$project$output_dir)
result <- run_pipeline(cfg, project_root)
if (!identical(result$status, "completed")) quit(save = "no", status = 1L)
