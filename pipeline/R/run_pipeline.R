# Per-run warning records. Stage scopes restore themselves when a step fails.

with_run_stage <- function(stage, out_dir, code) {
  state <- getOption("transatlas.warning_state")
  if (is.null(state)) return(force(code))
  previous <- list(stage = state$stage, out_dir = state$out_dir)
  on.exit({ state$stage <- previous$stage; state$out_dir <- previous$out_dir }, add = TRUE)
  state$stage <- stage
  state$out_dir <- safe_dir(out_dir)
  path <- file.path(state$out_dir, "warnings.csv")
  if (!path %in% state$initialized_paths) {
    state$initialized_paths <- c(state$initialized_paths, path)
    utils::write.csv(state$rows[FALSE, , drop = FALSE], path, row.names = FALSE)
  }
  force(code)
}

record_pipeline_warning <- function(w) {
  state <- getOption("transatlas.warning_state")
  if (is.null(state)) return(invisible(NULL))
  row <- data.frame(
    run_id = state$run_id, sequence = nrow(state$rows) + 1L,
    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stage = state$stage, warning_class = paste(class(w), collapse = "/"),
    message = conditionMessage(w),
    call = paste(deparse(conditionCall(w)), collapse = " "),
    stringsAsFactors = FALSE
  )
  state$rows <- rbind(state$rows, row)
  paths <- unique(c(state$path, file.path(state$out_dir, "warnings.csv")))
  for (path in paths) {
    utils::write.table(row, path, sep = ",", row.names = FALSE,
                       col.names = !file.exists(path), append = file.exists(path),
                       qmethod = "double", na = "")
  }
  invisible(row)
}

finalize_warning_records <- function(state) {
  rows <- state$rows
  if (nrow(rows)) {
    key <- interaction(rows$stage, rows$warning_class, rows$message, drop = TRUE)
    summary <- do.call(rbind, lapply(split(rows, key), function(x) {
      data.frame(stage = x$stage[1L], warning_class = x$warning_class[1L],
                 message = x$message[1L], occurrences = nrow(x),
                 first_sequence = min(x$sequence), stringsAsFactors = FALSE)
    }))
    summary <- summary[order(summary$first_sequence), , drop = FALSE]
  } else {
    summary <- data.frame(stage = character(), warning_class = character(),
                          message = character(), occurrences = integer(),
                          first_sequence = integer())
  }
  utils::write.csv(summary, file.path(dirname(state$path), "warning_summary.csv"), row.names = FALSE)
  jsonlite::write_json(list(run_id = state$run_id, warnings = nrow(rows),
                            distinct_stage_messages = nrow(summary),
                            records = "warnings.csv", summary = "warning_summary.csv"),
                       file.path(dirname(state$path), "status_warnings.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  if (nrow(rows)) message("Recorded ", nrow(rows), " warning(s): ", state$path)
  invisible(summary)
}

run_pipeline <- function(cfg, project_root = attr(cfg, "project_root")) {
  output_dir <- safe_dir(cfg$project$output_dir)
  if (!isTRUE(cfg$project$overwrite) && file.exists(file.path(output_dir, "run_manifest.json"))) {
    stop("Output directory contains a previous run. Set overwrite = TRUE or choose a new output_dir.",
         call. = FALSE)
  }
  state <- new.env(parent = emptyenv())
  state$run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%OS6", tz = "UTC"), "-", Sys.getpid())
  state$stage <- "setup"
  state$out_dir <- output_dir
  state$path <- file.path(output_dir, "warnings.csv")
  state$initialized_paths <- state$path
  state$rows <- data.frame(run_id = character(), sequence = integer(), timestamp_utc = character(),
                           stage = character(), warning_class = character(), message = character(),
                           call = character(), stringsAsFactors = FALSE)
  utils::write.csv(state$rows, state$path, row.names = FALSE)
  old <- options(transatlas.warning_state = state)
  on.exit(options(old), add = TRUE)
  on.exit(finalize_warning_records(state), add = TRUE)
  withCallingHandlers(run_pipeline_impl(cfg, project_root), warning = record_pipeline_warning)
}

# Shared execution path: builds a run configuration from a handful of arguments
# and runs the ported analysis sections.
#
# Two entry points call it:
#   pipeline/run_all.R      - run_all(counts_file = ..., species = ..., ...)
#   pipeline/run_analysis.R - command-line wrapper (--counts, --species, ...)
#
# Everything the user normally edits lives in the run_all() call. The defaults
# for every field are built into default_config() in helpers.R, so no
# configuration file is read for a normal run.

# Sections run in the order a paper reports them: differential expression, then
# functional enrichment, then transcriptional regulation, then the protein-level
# extensions. The keys are the numbers used on the command line and in the
# output folder names.
SECTION_REGISTRY <- list(
  `1`  = list(file = "s01_screening.R",           fn = "run_screening",  slug = "screening"),
  `2`  = list(file = "s02_go_ora.R",              fn = "run_go_ora",     slug = "go_ora"),
  `3`  = list(file = "s03_gsea_go.R",             fn = "run_gsea_go",    slug = "gsea_go"),
  `4`  = list(file = "s04_kegg_pathview.R",       fn = "run_kegg",       slug = "kegg_pathview"),
  `5`  = list(file = "s05_tf_regulon.R",          fn = "run_tf_regulon", slug = "tf_regulon"),
  `6`  = list(file = "s06_regulatory_advanced.R", fn = "run_regulatory_activity", slug = "regulatory_activity"),
  `7`  = list(file = "s07_tf_network.R",          fn = "run_tf_network", slug = "tf_network"),
  `8`  = list(file = "s08_tf_clusters.R",         fn = "run_tf_modules", slug = "tf_modules"),
  `9`  = list(file = "s09_geneset_heatmaps.R",    fn = "run_geneset_heatmaps", slug = "geneset_heatmaps"),
  `10` = list(file = "s10_ppi.R",                 fn = "run_ppi",        slug = "ppi_network"),
  `11` = list(file = "s11_subcellular.R",         fn = "run_subcellular", slug = "subcellular_localization")
)

# ---------------------------------------------------------------------------
# Input selection helpers
# ---------------------------------------------------------------------------

# Interactive selection when a file is not named. Rscript runs are not
# interactive, so the message has to say what to pass instead.
choose_input_file <- function(what = "input file") {
  if (!interactive()) {
    stop("No ", what, " was supplied and this R session is not interactive.\n",
         "  Pass its path explicitly (e.g. run_all(counts_file = \"path/to/counts.xlsx\")),\n",
         "  or start the workflow from an interactive R / RStudio session to pick it in a dialog.",
         call. = FALSE)
  }
  path <- file.choose()
  if (!nzchar(path)) stop("No ", what, " selected.", call. = FALSE)
  path
}

# Column names without reading the whole table, used to match a metadata file
# against the sample columns of the count matrix.
peek_columns <- function(path, sheet = NULL) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) return(character())
    sheets <- readxl::excel_sheets(path)
    sh <- if (!is.null(sheet) && sheet %in% sheets) sheet else sheets[[1L]]
    dat <- readxl::read_excel(path, sheet = sh, n_max = 0, .name_repair = "unique")
  } else if (ext == "csv") {
    dat <- utils::read.csv(path, nrows = 0, check.names = FALSE)
  } else if (ext == "txt") {
    dat <- utils::read.delim(path, nrows = 0, check.names = FALSE)
  } else {
    return(character())
  }
  names(dat)
}

# Metadata files are looked for next to the data first, then in the project
# folders. A candidate is accepted only when its sample_id column matches the
# count-matrix columns exactly, so an unrelated file in the same folder is
# skipped instead of being silently applied to the wrong samples.
find_metadata_file <- function(counts_path, sample_columns, project_root, explicit = NULL,
                               annotation_columns = character()) {
  if (!is.null(explicit)) {
    if (is.data.frame(explicit)) return(explicit)
    return(resolve_project_path(explicit, project_root))
  }
  dirs <- unique(c(dirname(counts_path), file.path(project_root, "config"),
                   file.path(project_root, "input"), project_root))
  stem <- tools::file_path_sans_ext(basename(counts_path))
  named <- unlist(lapply(c("csv", "xlsx", "xls", "tsv", "txt"), function(ext) {
    file.path(dirname(counts_path), paste0(c(paste0(stem, "_metadata"),
                                           paste0(stem, "_samples"), "sample_metadata"), ".", ext))
  }))
  expected <- setdiff(sample_columns, unique(c(annotation_columns, "Exonic.gene.sizes",
                                                "gene_length", "GeneLength")))
  others <- unlist(lapply(dirs, function(d) {
    sort(list.files(d, pattern = "(metadata|samples|design)[^/]*\\.(csv|xlsx|xls|tsv|txt)$",
                    full.names = TRUE, ignore.case = TRUE))
  }))
  for (path in unique(c(named, others))) {
    if (!file.exists(path)) next
    meta <- tryCatch(read_metadata_table(path), error = record_optional_error)
    if (is.null(meta) || !"sample_id" %in% names(meta)) next
    ids <- trimws(as.character(meta$sample_id))
    if (length(ids) >= 2L && !anyNA(ids) && !anyDuplicated(ids) &&
        setequal(ids, expected)) return(path)
  }
  NULL
}

# results/<name of the data file>/ - one folder per input workbook, so two
# experiments never overwrite each other.
default_output_dir <- function(counts_path, project_root) {
  stem <- tools::file_path_sans_ext(basename(counts_path))
  stem <- gsub("[^A-Za-z0-9._-]+", "_", stem)
  file.path(project_root, "results", stem)
}

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Turn the run_all() arguments into the cfg list every section consumes.
# The defaults are built in (default_config() in helpers.R); no configuration
# file is read unless one is passed explicitly as config_file.
build_run_config <- function(counts_file = NULL,
                             samples = NULL,
                             species = "mouse",
                             group_column = "condition",
                             contrast = NULL,
                             contrast_name = NULL,
                             log2fc = 0.75,
                             min_fpkm = 5,
                             alpha = 0.05,
                             min_count = 10,
                             min_samples = 2,
                             gene_id_column = "gene_id",
                             gene_symbol_column = "GeneSymbol",
                             annotation_column = "Description",
                             counts_sheet = NULL,
                             count_columns = NULL,
                             expression_matrix_file = NULL,
                             expression_matrix_unit = "FPKM",
                             gene_annotation = NULL,
                             gene_length_column = "Exonic.gene.sizes",
                             top_tfs = 30,
                             confidence = c("A", "B"),
                             tf_min_fpkm = NULL,
                             tf_min_targets = 5,
                             tf_ora_min_overlap = 3,
                             tf_coherence_min_targets = 5,
                             tf_file = NULL,
                             tf_dorothea_dir = NULL,
                             max_genes = 30,
                             top_n_labels = 20,
                             volcano_label_expression_n = 5,
                             volcano_label_effect_n = 5,
                             up_cutoff = 2.5,
                             low_cutoff = 0.15,
                             log2fc_cap = NULL,
                             volcano_lfc_thresholds = c(0.15, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5),
                             expression_log10_threshold = NULL,
                             gene_set_lfc_threshold = NULL,
                             go_pvalue_cutoff = 0.05,
                             go_qvalue_cutoff = 0.05,
                             go_bp_pvalue_cutoff = 0.01,
                             go_bp_qvalue_cutoff = 0.01,
                             go_simplify_cutoff = 0.5,
                             gsea_min_size = 15,
                             gsea_max_size = 500,
                             gsea_pvalue_cutoff = 0.05,
                             gsea_top_sets = 9,
                             run_msigdb = TRUE,
                             msigdb_category = "H",
                             pathview_download = TRUE,
                             run_pathview = NULL,
                             pathview_max_pathways = 12,
                             network_layouts = c("fr", "kk", "stress", "circle"),
                             network_max_labels = 30,
                             cluster_tfs = 30,
                             cluster_min_jaccard = 0.1,
                             cluster_max = 8,
                             activity_methods = c("ulm", "aucell", "wmean"),
                             activity_min_size = 5,
                             activity_alpha = NULL,
                             run_gsva = TRUE,
                             gsva_collection = "H",
                             gsva_min_size = 10,
                             gsva_max_size = 500,
                             run_ppi = TRUE,
                             ppi_score_threshold = 0.4,
                             ppi_max_nodes = 300,
                             ppi_label_n = 4,
                             ppi_max_modules = 10,
                             ppi_module_annotation_alpha = 0.2,
                             ppi_cache_dir = NULL,
                             ppi_timeout = 120,
                             qc = TRUE,
                             qc_plots = NULL,
                             sections = NULL,
                             report = TRUE,
                             report_title = NULL,
                             report_author = NULL,
                             report_language = "en",
                             generate_pdf = FALSE,
                             summary_report = TRUE,
                             kegg_cache_dir = NULL,
                             kegg_download_timeout = 240,
                             check_packages = TRUE,
                             install_missing_packages = TRUE,
                             continue_on_error = TRUE,
                             output_dir = NULL,
                             overwrite = TRUE,
                             seed = 20261013,
                             project_root = normalizePath(".", mustWork = FALSE),
                             config_file = NULL) {
  cfg <- default_config()
  if (!is.null(config_file)) {
    # Optional legacy overlay: a YAML file merged over the built-in defaults.
    # Nothing in a normal run reads it.
    cfg <- modifyList(cfg, read_config(resolve_project_path(config_file, project_root)))
  }
  preset <- species_preset(species)

  # --- input files ---------------------------------------------------------
  counts_path <- counts_file %||% choose_input_file("count matrix file")
  counts_path <- resolve_project_path(counts_path, project_root)
  if (!file.exists(counts_path)) stop("Count matrix does not exist: ", counts_path, call. = FALSE)
  cfg$input$counts_file <- counts_path
  cfg$input$counts_sheet <- counts_sheet
  cfg$input$gene_id_column <- gene_id_column
  cfg$input$gene_symbol_column <- gene_symbol_column
  cfg$input$annotation_column <- annotation_column
  cfg$input$count_columns <- count_columns
  cfg$input$expression_matrix_file <- expression_matrix_file
  expression_matrix_unit <- match.arg(expression_matrix_unit, c("FPKM", "TPM", "normalized_counts"))
  cfg$input$expression_matrix_unit <- expression_matrix_unit
  # Optional path to the GTF/GFF that defines the exon-union gene lengths used
  # for FPKM. Without one the annotation is found through an option, an
  # environment variable or a *.gtf beside the count matrix, and otherwise from
  # the species annotation package.
  cfg$input$gene_annotation <- if (is.null(gene_annotation)) NULL else
    resolve_project_path(gene_annotation, project_root)
  # Column of the count workbook that carries the exon-union gene lengths. It is
  # the primary source for the FPKM display matrix; NULL or "none" switches the
  # column off and leaves the computed matrix to the annotation route.
  cfg$input$gene_length_column <- gene_length_column

  sample_columns <- count_columns %||% peek_columns(counts_path, counts_sheet)
  meta <- find_metadata_file(counts_path, sample_columns, project_root, explicit = samples,
                             annotation_columns = c(gene_id_column, gene_symbol_column,
                                                    annotation_column, gene_length_column))
  if (is.data.frame(meta)) {
    cfg$input$sample_metadata_file <- meta            # written after the output dir exists
  } else if (is.null(meta)) {
    meta <- if (interactive()) choose_input_file("sample metadata file (sample_id, condition)") else NULL
    if (is.null(meta)) {
      stop("No sample metadata found for ", basename(counts_path), ".\n",
           "  Expected CSV/Excel metadata with 'sample_id' matching the count columns, next to the data file\n",
           "  (e.g. ", tools::file_path_sans_ext(basename(counts_path)), "_metadata.csv),\n",
           "  or pass one: run_all(samples = \"path/to/metadata.csv\").", call. = FALSE)
    }
    cfg$input$sample_metadata_file <- meta
  } else {
    cfg$input$sample_metadata_file <- meta
  }

  # --- TF resource: DoRothEA package data --------------------------------
  # Two collections cover every supported species: dorothea_hs for human and
  # dorothea_mm for mouse and rat. A custom table passed as `tf_file` still
  # wins, and is the only case where an extra file is read.
  cfg$input$tf_source <- if (is.null(tf_file)) "dorothea" else "file"
  cfg$input$tf_database_file <- if (is.null(tf_file)) NULL else resolve_project_path(tf_file, project_root)
  cfg$input$tf_species <- preset$tf_dataset
  cfg$input$tf_dorothea_dir <- tf_dorothea_dir %||% cfg$input$tf_dorothea_dir %||% "data/dorothea"
  cfg$regulatory$confidence_keep <- toupper(as.character(confidence))

  # --- output folder -------------------------------------------------------
  out <- output_dir %||% default_output_dir(counts_path, project_root)
  out <- resolve_project_path(out, project_root)
  if (!isTRUE(overwrite) && dir.exists(out) && length(list.files(out)) > 0L) {
    base <- out
    i <- 2L
    while (dir.exists(out <- paste0(base, "_run", i)) && i < 100L) i <- i + 1L
  }
  cfg$project$output_dir <- out
  cfg$project$overwrite <- isTRUE(overwrite)
  cfg$project$seed <- as.integer(seed)

  # --- species -------------------------------------------------------------
  cfg$species <- modifyList(cfg$species, preset)
  cfg$analysis$msigdb_db_species <- preset$msigdb_db_species

  # --- design and thresholds ----------------------------------------------
  nonnegative <- list(min_count = min_count, top_n_labels = top_n_labels,
                      volcano_label_expression_n = volcano_label_expression_n,
                      volcano_label_effect_n = volcano_label_effect_n)
  positive <- list(min_samples = min_samples, max_genes = max_genes, top_tfs = top_tfs,
                   tf_min_targets = tf_min_targets, tf_ora_min_overlap = tf_ora_min_overlap)
  for (nm in names(c(nonnegative, positive))) {
    value <- c(nonnegative, positive)[[nm]]
    if (is.null(value) && nm == "top_tfs") next
    minimum <- if (nm %in% names(nonnegative)) 0 else 1
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < minimum || value != floor(value))
      stop(nm, " must be a finite integer >= ", minimum, ".", call. = FALSE)
  }
  if (!is.null(min_fpkm) && (!is.numeric(min_fpkm) || length(min_fpkm) != 1L || !is.finite(min_fpkm) || min_fpkm < 0))
    stop("min_fpkm must be NULL or a finite non-negative number.", call. = FALSE)
  for (nm in c("activity_alpha", "ppi_module_annotation_alpha")) {
    value <- get(nm)
    if (!is.null(value) && (length(value) != 1L || !is.finite(value) || value <= 0 || value > 1))
      stop(nm, " must be between 0 and 1.", call. = FALSE)
  }
  meta_preview <- read_metadata_table(cfg$input$sample_metadata_file)
  stop_if_missing(meta_preview$sample_id, "sample_id metadata")
  if (!group_column %in% names(meta_preview)) {
    stop("The metadata has no column '", group_column, "'. Available: ",
         paste(setdiff(names(meta_preview), "sample_id"), collapse = ", "), call. = FALSE)
  }
  groups <- sort(unique(as.character(meta_preview[[group_column]])))
  if (!is.null(contrast)) {
    if (length(contrast) != 2L) stop("contrast must be c(numerator, denominator).", call. = FALSE)
    numerator <- as.character(contrast[1L])
    denominator <- as.character(contrast[2L])
    miss <- setdiff(c(numerator, denominator), groups)
    if (length(miss)) {
      stop("Contrast level(s) absent from the metadata: ", paste(miss, collapse = ", "),
           ". Levels are: ", paste(groups, collapse = ", "), call. = FALSE)
    }
  } else {
    if (length(groups) != 2L) {
      stop("The metadata has ", length(groups), " groups (", paste(groups, collapse = ", "),
           "); this workflow compares exactly two. Pass contrast = c(\"groupA\", \"groupB\").",
           call. = FALSE)
    }
    # Default: the later level (alphabetically) over the first, so
    # Control/Treatment becomes Treatment vs Control.
    denominator <- groups[1L]
    numerator <- groups[2L]
  }
  cfg$analysis$group_column <- group_column
  cfg$analysis$design <- paste0("~ ", group_column)
  cfg$analysis$reference_condition <- denominator
  contrast_label <- trimws(as.character(contrast_name %||% ""))
  cfg$analysis$contrasts <- list(list(name = if (nzchar(contrast_label)) contrast_label else
                                        paste0(numerator, "_vs_", denominator),
                                      variable = group_column,
                                      numerator = numerator, denominator = denominator))
  cfg$analysis$lfc_threshold <- as.numeric(log2fc)
  cfg$analysis$alpha <- as.numeric(alpha)
  cfg$analysis$min_count <- as.integer(min_count)
  cfg$analysis$min_samples <- as.integer(min_samples)
  cfg$analysis$max_display_features <- as.integer(max_genes)
  cfg$analysis$top_n_heatmap <- as.integer(max_genes)
  cfg$analysis$top_n_labels <- as.integer(top_n_labels)
  cfg$analysis$volcano_label_expression_n <- as.integer(volcano_label_expression_n)
  cfg$analysis$volcano_label_effect_n <- as.integer(volcano_label_effect_n)
  if (!is.null(log2fc_cap)) {
    warning("log2fc_cap is a deprecated display-only alias; use up_cutoff.", call. = FALSE)
    up_cutoff <- log2fc_cap
  }
  if (!is.null(up_cutoff) && (length(up_cutoff) != 1L || !is.finite(up_cutoff) || up_cutoff <= 0))
    stop("up_cutoff must be NULL or a finite positive number.", call. = FALSE)
  if (!is.null(low_cutoff) && (length(low_cutoff) != 1L || !is.finite(low_cutoff) || low_cutoff < 0))
    stop("low_cutoff must be NULL or a finite non-negative number.", call. = FALSE)
  cfg$analysis["up_cutoff"] <- list(if (is.null(up_cutoff)) NULL else as.numeric(up_cutoff))
  cfg$analysis["low_cutoff"] <- list(if (is.null(low_cutoff)) NULL else as.numeric(low_cutoff))
  cfg$analysis$log2fc_cap <- NULL
  cfg$analysis$volcano_x_limit <- if (is.null(up_cutoff)) NULL else as.numeric(up_cutoff) + 0.1
  thresholds <- suppressWarnings(as.numeric(unlist(volcano_lfc_thresholds)))
  if (!length(thresholds) || any(!is.finite(thresholds) | thresholds <= 0)) {
    stop("volcano_lfc_thresholds must contain finite positive numbers (CLI: --volcano-thresholds 0.5,0.75,1).", call. = FALSE)
  }
  if (length(log2fc) != 1L || !is.finite(log2fc) || log2fc <= 0) stop("log2fc must be positive.", call. = FALSE)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1.", call. = FALSE)
  cfg$analysis$volcano_lfc_thresholds <- sort(unique(c(thresholds, as.numeric(log2fc))))
  cfg$analysis$expression_log10_threshold <- if (is.null(expression_log10_threshold)) NULL else
    as.numeric(expression_log10_threshold)
  # Primary DEG membership is separate from the GSEA-only lower cutoff and
  # the display-only upper cutoff.
  if (!is.null(gene_set_lfc_threshold)) {
    warning("gene_set_lfc_threshold is deprecated and ignored; set log2fc instead.",
            call. = FALSE)
  }
  cfg$analysis$min_fpkm <- if (is.null(min_fpkm)) NULL else as.numeric(min_fpkm)
  cfg$analysis$go_pvalue_cutoff <- as.numeric(go_pvalue_cutoff)
  cfg$analysis$go_qvalue_cutoff <- as.numeric(go_qvalue_cutoff)
  cfg$analysis$go_bp_pvalue_cutoff <- as.numeric(go_bp_pvalue_cutoff)
  cfg$analysis$go_bp_qvalue_cutoff <- as.numeric(go_bp_qvalue_cutoff)
  cfg$analysis$go_simplify_cutoff <- as.numeric(go_simplify_cutoff)
  cfg$analysis$gsea_min_size <- as.integer(gsea_min_size)
  cfg$analysis$gsea_max_size <- as.integer(gsea_max_size)
  cfg$analysis$gsea_pvalue_cutoff <- as.numeric(gsea_pvalue_cutoff)
  cfg$analysis$gsea_top_sets <- as.integer(gsea_top_sets)
  cfg$analysis$run_msigdb <- isTRUE(run_msigdb)
  cfg$analysis$msigdb_category <- msigdb_category
  # pathview_download is the documented switch; run_pathview is the older name
  # and wins only when a script still supplies it.
  cfg$analysis$run_pathview <- if (!is.null(run_pathview)) {
    isTRUE(run_pathview)
  } else {
    isTRUE(pathview_download)
  }
  cfg$analysis$pathview_download <- cfg$analysis$run_pathview
  cfg$analysis$pathview_max_pathways <- as.integer(pathview_max_pathways)
  if (!is.null(tf_min_fpkm)) {
    warning("tf_min_fpkm is deprecated and ignored; set min_fpkm instead.", call. = FALSE)
  }
  cfg$regulatory$tf_min_fpkm <- cfg$analysis$min_fpkm
  cfg$regulatory$min_targets <- as.integer(tf_min_targets)
  cfg$regulatory$ora_min_overlap <- as.integer(tf_ora_min_overlap)
  cfg$regulatory$coherence_min_targets <- as.integer(tf_coherence_min_targets)
  cfg$regulatory$network_layouts <- as.list(as.character(network_layouts))
  cfg$regulatory$network_max_labels <- as.integer(network_max_labels)
  cfg$regulatory$activity_methods <- as.list(as.character(activity_methods))
  cfg$regulatory$activity_min_size <- as.integer(activity_min_size)
  cfg$regulatory$activity_alpha <- as.numeric(activity_alpha %||% alpha)
  cfg$pathway_activity$run <- isTRUE(run_gsva)
  cfg$pathway_activity$collection <- gsva_collection
  cfg$pathway_activity$min_size <- as.integer(gsva_min_size)
  cfg$pathway_activity$max_size <- as.integer(gsva_max_size)
  cfg$ppi$run <- isTRUE(run_ppi)
  cfg$ppi$score_threshold <- as.numeric(ppi_score_threshold)
  cfg$ppi$max_nodes <- as.integer(ppi_max_nodes)
  cfg$ppi$label_n <- as.integer(ppi_label_n)
  cfg$ppi$max_modules <- as.integer(ppi_max_modules)
  cfg$ppi$module_annotation_alpha <- as.numeric(ppi_module_annotation_alpha)
  cfg$ppi$cache_dir <- ppi_cache_dir
  cfg$ppi$timeout <- as.numeric(ppi_timeout)
  cfg$report$language <- report_language
  cfg$report$generate_html <- isTRUE(report)
  cfg$report$generate_pdf <- isTRUE(generate_pdf)
  cfg$report$summary_report <- isTRUE(summary_report)
  if (!is.null(report_title)) cfg$report$title <- report_title
  if (!is.null(report_author)) cfg$report$author <- report_author
  cfg$runtime$kegg_cache_dir <- kegg_cache_dir
  cfg$runtime$kegg_download_timeout <- as.numeric(kegg_download_timeout)
  cfg$runtime$continue_on_optional_module_error <- isTRUE(continue_on_error)

  # --- quality control -----------------------------------------------------
  qc_flags <- cfg$quality_control
  if (!isTRUE(qc)) {
    for (nm in names(qc_flags)) qc_flags[[nm]] <- FALSE
  } else if (!is.null(qc_plots)) {
    wanted <- tolower(as.character(qc_plots))
    unknown_qc <- setdiff(wanted, names(QC_PLOT_FLAGS))
    if (length(unknown_qc)) {
      stop("Unknown QC plot(s): ", paste(unknown_qc, collapse = ", "),
           ". Use: ", paste(names(QC_PLOT_FLAGS), collapse = ", "), call. = FALSE)
    }
    for (nm in names(QC_PLOT_FLAGS)) qc_flags[[QC_PLOT_FLAGS[[nm]]]] <- nm %in% wanted
  }
  cfg$quality_control <- qc_flags

  # --- transcription factors ----------------------------------------------
  # One number controls both figures: the dot plot shows at most this many
  # regulons, and the network draws the same set.
  if (!is.null(top_tfs)) {
    cfg$regulatory$tf_dotplot_max <- as.integer(top_tfs)
    cfg$regulatory$network_tfs <- as.integer(top_tfs)
    cfg$regulatory$cluster_tfs <- as.integer(top_tfs)
  }
  if (!is.null(cluster_tfs)) cfg$regulatory$cluster_tfs <- as.integer(cluster_tfs)
  if (!is.null(cluster_min_jaccard)) cfg$regulatory$cluster_min_jaccard <- as.numeric(cluster_min_jaccard)
  if (!is.null(cluster_max)) cfg$regulatory$cluster_max <- as.integer(cluster_max)
  cfg$runtime$check_packages <- isTRUE(check_packages)
  cfg$runtime$install_missing_packages <- isTRUE(install_missing_packages)

  # --- sections ------------------------------------------------------------
  requested <- if (is.null(sections)) as.character(unlist(cfg$sections$run %||% 1:11)) else
    as.character(sections)
  unknown <- setdiff(requested, names(SECTION_REGISTRY))
  if (length(unknown)) stop("Unknown section(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  cfg$sections$run <- as.integer(requested)

  # Paths must be absolute before the TF resource is probed for.
  cfg <- resolve_config_paths(cfg, project_root)
  tf_ok <- if (identical(cfg$input$tf_source, "file")) {
    !is.null(cfg$input$tf_database_file) && file.exists(cfg$input$tf_database_file)
  } else {
    dorothea_available(cfg)
  }
  if (!tf_ok && any(c("5", "6", "7", "8", "9") %in% requested)) {
    warning("DoRothEA regulons are unavailable for ", preset$name,
            "; TF-dependent outputs will be marked as skipped. Install dorothea or keep data/dorothea/ in the project.",
            call. = FALSE)
  }
  cfg
}

# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

run_pipeline_impl <- function(cfg, project_root = attr(cfg, "project_root")) {
  project_root <- project_root %||% "."
  old_options <- options(transatlas.project_root = project_root,
                         transatlas.diagnostics_file = getOption("transatlas.diagnostics_file"))
  on.exit(options(old_options), add = TRUE)
  # The registry is refreshed at the start of every run, whichever entry point
  # (run_all, run_analysis.R, a direct run_pipeline call) started it.
  if (exists("reload_figure_sizes", mode = "function")) {
    reload_figure_sizes(project_root)
  }
  if (exists("reset_figure_size_log", mode = "function")) reset_figure_size_log()
  requested <- as.character(unlist(cfg$sections$run))
  output_dir <- safe_dir(cfg$project$output_dir)
  run_issues <- list()
  finished <- FALSE
  options(transatlas.diagnostics_file = file.path(output_dir, "optional_errors.log"))
  write_status(file.path(output_dir, "run_status.json"), "pipeline", "running")
  on.exit({
    if (!finished) write_status(file.path(output_dir, "run_status.json"), "pipeline", "failed",
                                "Run stopped before completion; inspect the console and section status files.",
                                list(issues = run_issues))
  }, add = TRUE)
  # Every package is checked (and installed when missing) before any input is
  # read, so a run either has its dependencies or stops with the exact list.
  if (isTRUE(cfg$runtime$check_packages)) with_run_stage("dependencies", output_dir, check_and_install_packages(cfg))
  if (!isTRUE(cfg$project$overwrite) && file.exists(file.path(output_dir, "run_manifest.json"))) {
    stop("Output directory contains a previous run. Set overwrite = TRUE or choose a new output_dir.",
         call. = FALSE)
  }
  if (!is.null(cfg$project$seed)) set.seed(as.integer(cfg$project$seed))
  # One thread only: several packages here (fgsea, clusterProfiler) fall back to
  # forked workers, and a forked worker that inherits an open SQLite or network
  # handle can block indefinitely on this host.
  old_parallel <- options(mc.cores = 1L)
  on.exit(options(old_parallel), add = TRUE)
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    previous_backend <- BiocParallel::bpparam()
    on.exit(BiocParallel::register(previous_backend, default = TRUE), add = TRUE)
    BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
  }

  message("Reading raw count matrix from ", basename(cfg$input$counts_file), " ...")
  counts_obj <- with_run_stage("input", output_dir, read_counts_workbook(cfg))

  # The metadata may have arrived as a data.frame; persist it so the run is
  # reproducible from its own output folder.
  if (is.data.frame(cfg$input$sample_metadata_file)) {
    meta_path <- file.path(output_dir, "sample_metadata.csv")
    utils::write.csv(cfg$input$sample_metadata_file, meta_path, row.names = FALSE)
    cfg$input$sample_metadata_file <- meta_path
  }
  meta <- with_run_stage("metadata", output_dir, read_sample_metadata(cfg, colnames(counts_obj$counts)))
  group_col <- cfg$analysis$group_column %||% "condition"
  for (contrast in cfg$analysis$contrasts) {
    variable <- contrast$variable %||% group_col
    if (!variable %in% names(meta)) stop("Contrast variable is missing: ", variable, call. = FALSE)
    if (!all(c(contrast$numerator, contrast$denominator) %in% levels(meta[[variable]]))) {
      stop("Contrast levels are absent from metadata for ", contrast$name, call. = FALSE)
    }
  }
  utils::write.csv(meta, file.path(output_dir, "sample_metadata_resolved.csv"), row.names = TRUE)

  # DESeq2 is the shared upstream step; cache it so individual sections can be
  # re-run without repeating the model fit.
  cache_dir <- safe_dir(file.path(output_dir, ".cache"))
  cache_file <- file.path(cache_dir, "deseq_object.rds")
  signature <- deseq_cache_signature(cfg, counts_obj, meta)
  cached <- NULL
  if (file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = record_optional_error)
    if (!is.null(cached) && identical(cached$signature, signature)) {
      message("Reusing cached DESeq2 fit.")
      deseq_obj <- cached$object
      with_run_stage("deseq2", file.path(output_dir, "00_differential_expression"),
        write_cached_deseq2_tables(cfg, counts_obj, deseq_obj,
                                  file.path(output_dir, "00_differential_expression")))
    } else {
      message("DESeq2 cache invalidated: input, metadata, model settings or package versions changed.")
      cached <- NULL
    }
  }
  if (is.null(cached)) {
    deseq_obj <- with_run_stage("deseq2", file.path(output_dir, "00_differential_expression"),
        run_deseq2(cfg, counts_obj, meta, file.path(output_dir, "00_differential_expression")))
    saveRDS(list(object = deseq_obj, signature = signature), cache_file)
  }
  qc_dir <- safe_dir(file.path(output_dir, "00_differential_expression"))
  tryCatch({
    with_run_stage("quality_control", qc_dir, run_qc_module(cfg, deseq_obj, counts_obj, meta, qc_dir))
    write_status(file.path(qc_dir, "status_qc.json"), "quality_control",
                 if (any(vapply(cfg$quality_control, isTRUE, logical(1)))) "completed" else "skipped",
                 if (!any(vapply(cfg$quality_control, isTRUE, logical(1)))) "QC disabled by configuration" else NULL)
  }, error = function(e) {
    write_status(file.path(qc_dir, "status_qc.json"), "quality_control", "failed", conditionMessage(e))
    record_optional_error(e)
    run_issues[["quality_control"]] <<- conditionMessage(e)
  })

  ctx <- with_run_stage("context", qc_dir, build_context(cfg, counts_obj, meta, deseq_obj))
  # The screened gene tables and the ranked identifier vector are built by
  # section 1, but every later section consumes them. When a run starts after
  # section 1, rebuild them here so each section stays independently runnable.
  if (!("1" %in% requested)) {
    ctx <- attach_ranked_context(cfg, ctx, screening_tables(cfg, ctx))
  }
  # Same idea for the sections that consume another section's table: reload it
  # from disk when its producer is not part of this run. Only empty slots are
  # filled, so a full run is unaffected.
  ctx <- hydrate_context_from_disk(cfg, ctx, output_dir)
  if (is.null(ctx$regulons) && any(c("5", "6", "7", "8", "9") %in% requested)) {
    # Without the TF resource sections 5-9 can only write a "skipped"
    # status, leaving what looks like an empty output directory. Say exactly
    # which resource is missing, rather than leaving a one-line status to guess.
    tf_source <- cfg$input$tf_source %||% "dorothea"
    where <- if (identical(tf_source, "file")) {
      paste0("file ", cfg$input$tf_database_file %||% "(unset)")
    } else {
      paste0("DoRothEA dataset ", cfg$input$tf_species, " from the dorothea package or ",
             file.path(cfg$input$tf_dorothea_dir %||% "data/dorothea", paste0("dorothea_", cfg$input$tf_species, ".rda")))
    }
    missing_tf <- paste0(
      "The TF-target resource could not be loaded, but section(s) ",
      paste(intersect(c("5", "6", "7", "8", "9"), requested), collapse = ", "),
      " were requested.\n",
      "  Expected: ", where)
    warning(missing_tf, call. = FALSE)
  }

  run_section <- function(n) {
    entry <- SECTION_REGISTRY[[n]]
    path <- file.path(project_root, "pipeline", "sections", entry$file)
    out_dir <- safe_dir(file.path(output_dir, sprintf("section%02d_%s", as.integer(n), entry$slug)))
    status_path <- file.path(out_dir, "status.json")
    old_diagnostics <- options(transatlas.diagnostics_file = file.path(out_dir, "optional_errors.log"))
    on.exit(options(old_diagnostics), add = TRUE)
    write_status(status_path, paste0("section", n), "running")
    message("=== Section ", n, ": ", entry$slug, " ===")
    start <- Sys.time()
    result <- with_run_stage(paste0("section", n, "_", entry$slug), out_dir, tryCatch({
      if (!file.exists(path)) stop("Section file missing: ", entry$file)
      section_env <- new.env(parent = environment(run_pipeline))
      sys.source(path, envir = section_env)
      fn <- get(entry$fn, envir = section_env, mode = "function", inherits = FALSE)
      value <- fn(cfg, ctx, out_dir)
      if (!is.list(value) || !all(c("de", "meta", "contrast_name") %in% names(value))) {
        stop("Section returned an invalid context; expected a context list with de, meta and contrast_name.")
      }
      st <- jsonlite::read_json(status_path)
      if (identical(st$status, "running")) write_status(status_path, paste0("section", n), "completed")
      if (identical(st$status, "failed")) stop(st$message %||% "Section reported failure.")
      value
    }, error = function(e) {
      write_status(status_path, paste0("section", n), "failed", conditionMessage(e))
      run_issues[[paste0("section", n)]] <<- conditionMessage(e)
      if (n == "1" || !isTRUE(cfg$runtime$continue_on_optional_module_error)) stop(e)
      warning("Section ", n, " (", entry$slug, ") failed: ", conditionMessage(e), call. = FALSE)
      ctx
    }))
    message(sprintf("    section %d finished in %.1f s", as.integer(n),
                    as.numeric(difftime(Sys.time(), start, units = "secs"))))
    result
  }

  for (n in requested) ctx <- run_section(n)

  with_run_stage("provenance", file.path(output_dir, "provenance"),
    write_provenance(cfg, project_root, counts_obj, meta, file.path(output_dir, "provenance"), deseq_obj))
  # Sizes actually used, one row per figure: the registry values read at the
  # start of the run, so a figure's dimensions can be checked after the fact.
  sizes <- figure_size_log()
  if (nrow(sizes)) {
    utils::write.csv(sizes, file.path(output_dir, "provenance", "figure_sizes.csv"),
                     row.names = FALSE)
  }
  report_dir <- file.path(output_dir, "report")
  if (isTRUE(cfg$report$generate_html) &&
      file.exists(file.path(project_root, "pipeline", "R", "report_module.R"))) {
    run_report_step(
      "manuscript_html", report_dir, "status_manuscript.json",
      target = file.path(report_dir, "manuscript_draft.html"),
      code = {
        source(file.path(project_root, "pipeline", "R", "report_module.R"), local = TRUE)
        render_manuscript(cfg, project_root, output_dir, report_dir)
      })
  }
  # Fixed-template result report: fills the same sections on every run.
  summary_path <- file.path(project_root, "pipeline", "R", "summary_report.R")
  if (isTRUE(cfg$report$summary_report) && file.exists(summary_path)) {
    run_report_step(
      "summary_report_docx", report_dir, "status_summary.json",
      target = file.path(report_dir, "results_report.docx"),
      code = {
        source(summary_path, local = TRUE)
        render_summary_report(cfg, output_dir, report_dir)
      })
  } else if (isTRUE(cfg$report$summary_report)) {
    warning("summary_report is TRUE but pipeline/R/summary_report.R is missing; ",
            "the Word report was not generated.", call. = FALSE)
  }
  for (status_file in c(if (isTRUE(cfg$report$generate_html)) "status_manuscript.json",
                        if (isTRUE(cfg$report$summary_report)) "status_summary.json")) {
    path <- file.path(report_dir, status_file)
    st <- if (file.exists(path)) jsonlite::read_json(path) else NULL
    if (is.null(st) || !identical(st$status, "completed")) {
      run_issues[[status_file]] <- if (is.null(st)) "Report status missing" else st$message %||% st$status
    }
  }
  status <- if (length(run_issues)) "completed_with_errors" else "completed"
  write_status(file.path(output_dir, "run_status.json"), "pipeline", status,
               details = list(issues = run_issues, requested_sections = requested))
  finished <- TRUE
  message(if (length(run_issues)) paste0("Analysis finished with ", length(run_issues), " failed step(s). Results: ") else
            "Integrated analysis completed. Results: ", output_dir)
  invisible(list(cfg = cfg, ctx = ctx, output_dir = output_dir, status = status, issues = run_issues))
}
