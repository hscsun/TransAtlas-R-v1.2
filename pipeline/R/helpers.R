# Reusable helpers for the configurable bulk RNA-seq workflow.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

stop_if_missing <- function(x, label) {
  if (length(x) == 0L || any(is.na(x)) || any(!nzchar(as.character(x)))) {
    stop("Missing or invalid ", label, call. = FALSE)
  }
}

safe_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

write_status <- function(path, module, status, message = NULL, details = list()) {
  payload <- c(list(module = module, status = status, message = message,
                    timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), details)
  jsonlite::write_json(payload, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
}

# Run one report-building step without hiding why it failed. Earlier versions
# wrapped these calls in try(..., silent = TRUE); when officer or flextable was
# missing on macOS the DOCX simply never appeared and the run still ended with
# "completed". The step now records its own status, keeps the full error message
# and a short session snapshot in <step>_error.txt, and emits a warning that is
# visible in an interactive session as well as in the console log.
run_report_step <- function(label, report_dir, status_name, target, code) {
  report_dir <- safe_dir(report_dir)
  error_file <- file.path(report_dir, paste0(label, "_error.txt"))
  warning_file <- file.path(report_dir, paste0(label, "_warnings.txt"))
  status_file <- file.path(report_dir, status_name)
  unlink(error_file)
  unlink(warning_file)
  unlink(status_file)
  stage <- if (exists("with_run_stage", mode = "function")) with_run_stage else
    function(stage, out_dir, code) force(code)
  result <- tryCatch(
    withCallingHandlers(
      stage(label, report_dir, force(code)),
      warning = function(w) {
        # Warnings do not abort report generation, but they are preserved so a
        # half-written document can be diagnosed afterwards.
        if (exists("record_pipeline_warning", mode = "function")) record_pipeline_warning(w)
        cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "WARNING:", conditionMessage(w),
            "\n", file = warning_file, append = TRUE)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e)
  if (!inherits(result, "error") && !is.null(target) && !file.exists(target)) {
    result <- simpleError(paste0("Report builder did not produce the requested file: ", target))
  }
  if (inherits(result, "error")) {
    msg <- conditionMessage(result)
    writeLines(c(
      paste0("Report step: ", label),
      paste0("Time: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      paste0("Target: ", target %||% "(not written)"),
      paste0("R: ", R.version.string),
      paste0("Platform: ", R.version$platform),
      paste0("Library paths: ", paste(.libPaths(), collapse = ", ")),
      "",
      "Error message:",
      msg,
      "",
      "Traceback:",
      paste(utils::capture.output(traceback(2)), collapse = "\n")
    ), error_file)
    write_status(status_file, label, "failed", msg,
                 list(target = target, error_file = basename(error_file),
                      r_version = R.version.string, platform = R.version$platform))
    warning(sprintf("%s failed: %s\n  Details written to %s", label, msg, error_file),
            call. = FALSE)
    return(invisible(NULL))
  }
  # The builder may have written a richer status (section and figure counts).
  # Keep it instead of replacing it with the wrapper's minimal record.
  # A builder that failed before writing its own status leaves no file; read
  # only when it exists, so the run does not end with a spurious
  # "cannot open file ... status_manuscript.json" warning.
  existing <- if (file.exists(status_file)) {
    tryCatch(jsonlite::read_json(status_file, simplifyVector = FALSE),
             error = record_optional_error)
  } else NULL
  if (is.null(existing) || !identical(existing$status, "completed")) {
    write_status(status_file, label, "completed", NULL,
                 list(target = if (is.null(target)) NA_character_ else basename(target)))
  }
  result
}

# Cache the actual selected count matrix, complete metadata and model settings.
# Explicit schema changes invalidate older cache formats.
deseq_cache_signature <- function(cfg, counts_obj, meta) {
  key <- list(schema = 2L, counts = counts_obj$counts, annotations = counts_obj$annotations,
              metadata = meta,
              analysis = cfg$analysis[c("design", "group_column", "reference_condition", "contrasts",
                "min_count", "min_samples", "alpha", "independent_filtering", "lfc_shrinkage", "lfc_threshold")],
              seed = cfg$project$seed,
              versions = vapply(c("DESeq2", "apeglm", "SummarizedExperiment"), function(p) {
                if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "unavailable"
              }, character(1)))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(key, path, compress = FALSE, version = 2)
  list(schema = 2L, md5 = unname(tools::md5sum(path)))
}

# Preserve the reason when an optional plot, lookup or method cannot run.
record_optional_error <- function(e) {
  path <- getOption("transatlas.diagnostics_file")
  msg <- paste0(format(Sys.time(), tz = "UTC", usetz = TRUE), " | ",
                paste(deparse(conditionCall(e)), collapse = " "), " | ", conditionMessage(e))
  if (!is.null(path)) cat(msg, "\n", file = path, append = TRUE)
  message("    Optional operation unavailable: ", conditionMessage(e))
  NULL
}

close_plot_device <- function(device) {
  if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
  invisible(NULL)
}

read_config <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Package 'yaml' is required.", call. = FALSE)
  cfg <- yaml::read_yaml(path)
  cfg$.config_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  cfg
}

# Every setting the workflow needs, in code. A run no longer reads a
# configuration file: build_run_config() starts from this list and applies the
# arguments of run_all(). Values are the defaults documented on run_all().
default_config <- function() {
  list(
    project = list(name = "bulk_RNAseq_integrated", output_dir = NULL,
                   overwrite = TRUE, seed = 20261013L),
    input = list(
      counts_file = NULL, counts_sheet = NULL,
      gene_id_column = "gene_id", gene_symbol_column = "GeneSymbol",
      annotation_column = "Description", expression_matrix_file = NULL, expression_matrix_unit = "FPKM",
      gene_annotation = NULL, gene_length_column = "Exonic.gene.sizes",
      count_columns = NULL,
      sample_metadata_file = NULL,
      tf_source = "dorothea", tf_species = "mm",
      tf_dorothea_dir = "data/dorothea", tf_database_file = NULL),
    species = species_preset("mouse"),
    analysis = list(
      group_column = "condition", design = "~ condition",
      reference_condition = "Control", min_replicates_per_group = 2L,
      contrasts = list(list(name = "Treatment_vs_Control", variable = "condition",
                            numerator = "Treatment", denominator = "Control")),
      min_count = 10L, min_samples = 2L, independent_filtering = TRUE,
      lfc_shrinkage = "apeglm", alpha = 0.05, lfc_threshold = 0.75,
      volcano_lfc_thresholds = list(0.15, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5),
      # Legacy AnaDGEsS1 rule (log10(baseMean) <= 2 hidden from the up/down
      # classes). Off by default so the FPKM floor in `min_fpkm` is the only
      # expression threshold; set it to 2 to restore the old behaviour.
      expression_log10_threshold = NULL,
      up_cutoff = 2.5, low_cutoff = 0.15, volcano_x_limit = 2.6, volcano_y_limit = NULL,
      # Parenthesis threshold used for the chromosome landscape (section 1);
      # kept separate from the primary DEG cut so the map can show the broader
      # response while the tables stay at the reporting threshold.
      chromosome_lfc = 0.75,
      # Section 1 draws one landscape per threshold so the effect-size
      # dependence of the chromosomal pattern is visible directly.
      chromosome_lfc_thresholds = c(0.5, 0.75, 1),
      top_n_heatmap = 30L, top_n_labels = 20L, max_display_features = 30L,
      volcano_label_expression_n = 5L, volcano_label_effect_n = 5L,
      heatmap_transform = "vst",
      # One expression floor for every section: the higher of a gene's two
      # group-mean FPKM values must reach this value to enter the analysis.
      min_fpkm = 5,
      go_pvalue_cutoff = 0.05, go_qvalue_cutoff = 0.05,
      go_bp_pvalue_cutoff = 0.01, go_bp_qvalue_cutoff = 0.01,
      go_simplify_cutoff = 0.5,
      gsea_min_size = 15L, gsea_max_size = 500L, gsea_pvalue_cutoff = 0.05,
      gsea_top_sets = 9L, gsea_all_sets_max_pages = 60L,
      run_pathview = TRUE, pathview_download = TRUE, pathview_max_pathways = 12L,
      run_msigdb = TRUE, msigdb_category = "H", msigdb_subcategory = NULL,
      msigdb_db_species = "MM"),
    sections = list(run = as.list(1:11)),
    quality_control = list(
      run_pca = TRUE, run_sample_correlation = TRUE, run_library_size = TRUE,
      run_count_distribution = TRUE, run_ma_plot = TRUE, run_dispersion_plot = TRUE,
      run_independent_filter_plot = TRUE, run_sample_distance_heatmap = TRUE,
      run_rlog_vst_comparison = FALSE),
    regulatory = list(
      confidence_keep = list("A", "B"), min_targets = 5L, ora_min_overlap = 3L,
      tf_dotplot_max = 30L, tf_chromosome_n = 10L,
      gsea_min_size = 10L, gsea_max_size = 500L,
      network_tfs = 15L, network_layouts = list("fr", "kk", "stress", "circle"),
      network_max_labels = 30L,
      cluster_tfs = 30L, cluster_min_jaccard = 0.1, cluster_max = 8L,
      cluster_max_edges = 120L,
      activity_methods = list("ulm", "aucell", "wmean"),
      activity_min_size = 5L, activity_alpha = 0.05,
      coherence_min_targets = 5L),
    pathway_activity = list(run = TRUE, method = "gsva", collection = "H",
                            min_size = 10L, max_size = 500L),
    ppi = list(run = TRUE, score_threshold = 0.4, max_nodes = 300L,
               label_n = 4L, max_modules = 10L, module_go = TRUE, module_annotation_alpha = 0.2,
               cache_dir = NULL, timeout = 120),
    report = list(language = "en", author = "ClawsGO Science Agent",
                  title = "Bulk RNA-seq treatment-versus-control demonstration",
                  generate_html = TRUE, generate_pdf = FALSE,
                  summary_report = TRUE,
                  include_journal_recommendations = TRUE,
                  include_workflow_figure = TRUE),
    runtime = list(use_renv = FALSE, bioconductor_version = "3.23", r_version = "4.6",
                   check_packages = TRUE, install_missing_packages = TRUE,
                   internet_required_for_annotation = TRUE,
                   continue_on_optional_module_error = TRUE,
                   kegg_cache_dir = NULL, kegg_download_timeout = 240)
  )
}

# Names of the quality-control figures accepted by run_all(qc_plots = ...).
QC_PLOT_FLAGS <- c(
  pca = "run_pca", correlation = "run_sample_correlation",
  library_size = "run_library_size", count_distribution = "run_count_distribution",
  ma_plot = "run_ma_plot", dispersion = "run_dispersion_plot",
  independent_filter = "run_independent_filter_plot",
  sample_distance = "run_sample_distance_heatmap",
  rlog_vst = "run_rlog_vst_comparison")

# Chromosome lengths (Mb) of the current reference assemblies, used by the
# chromosome landscape so every track spans the real chromosome rather than the
# span covered by the genes that happen to be on the array. Gene coordinates
# come from the OrgDb CHRLOC map.
reference_chrom_lengths <- function(common_name = "mouse") {
  mouse <- c(`1` = 195.47, `2` = 181.74, `3` = 159.97, `4` = 156.51, `5` = 151.83,
             `6` = 149.74, `7` = 145.44, `8` = 129.40, `9` = 124.60, `10` = 130.69,
             `11` = 122.08, `12` = 120.13, `13` = 120.42, `14` = 124.90, `15` = 104.04,
             `16` = 98.20, `17` = 95.20, `18` = 90.70, `19` = 61.43,
             X = 169.48, Y = 91.47, MT = 0.0163)
  human <- c(setNames(c(248.96, 242.19, 198.30, 190.21, 181.54, 170.81, 159.35,
                        145.14, 138.39, 133.80, 135.09, 133.28, 114.36, 107.04,
                        101.99, 90.34, 83.26, 80.37, 58.62, 64.44, 46.71, 50.82),
                      as.character(1:22)),
             X = 156.04, Y = 57.23, MT = 0.0166)
  rat <- c(setNames(c(282.76, 266.44, 171.74, 187.13, 173.11, 147.30, 143.02,
                      129.04, 114.19, 112.63, 87.69, 46.63, 111.15, 112.19,
                      110.72, 90.24, 90.86, 83.83, 62.25, 56.22),
                    as.character(1:20)),
           X = 159.97, Y = 33.36, MT = 0.0163)
  switch(common_name, human = human, rat = rat, mouse)
}

# Louvain communities can fragment a dense graph into many tiny groups. Merge
# the smallest community into the community it shares most edges with until at
# most `max_n` remain; isolated communities fall back to the largest module.
reduce_communities <- function(g, membership, max_n = 10L) {
  membership <- as.integer(membership)
  max_n <- max(1L, as.integer(max_n))
  repeat {
    tab <- table(membership)
    if (length(tab) <= max_n) break
    small <- as.integer(names(tab)[which.min(tab)])
    idx <- which(membership == small)
    nb <- unique(unlist(igraph::neighborhood(g, order = 1, nodes = idx)))
    nb <- setdiff(nb, idx)
    target <- if (length(nb)) {
      vt <- table(membership[nb])
      as.integer(names(vt)[which.max(vt)])
    } else {
      as.integer(names(tab)[which.max(tab)])
    }
    membership[idx] <- target
  }
  membership
}

canonical_chromosomes <- function(common_name = "mouse") {
  n <- switch(common_name, human = 22L, rat = 20L, 19L)
  c(as.character(seq_len(n)), "X", "Y", "MT")
}

# Map gene symbols to CHRLOC coordinates through the species OrgDb. Only the
# supplied genes are queried, so this stays cheap even on a full transcriptome.
chromosome_positions <- function(orgdb, genes, common_name = "mouse") {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) return(NULL)
  loc <- suppressWarnings(AnnotationDbi::select(
    orgdb, keys = genes, columns = c("CHRLOC", "CHRLOCEND"), keytype = "SYMBOL"))
  loc <- loc[!is.na(loc$CHRLOC) & !is.na(loc$CHRLOCCHR), , drop = FALSE]
  if (!nrow(loc)) return(NULL)
  loc$chr <- sub("^chr", "", as.character(loc$CHRLOCCHR))
  loc$start <- abs(as.numeric(loc$CHRLOC))
  loc$end <- abs(as.numeric(loc$CHRLOCEND))
  loc$strand <- ifelse(as.numeric(loc$CHRLOC) < 0, "-", "+")
  # Several transcripts can share one symbol; keep the most 5' start.
  loc <- loc[order(loc$SYMBOL, loc$start), , drop = FALSE]
  loc <- loc[!duplicated(loc$SYMBOL), , drop = FALSE]
  loc <- loc[loc$chr %in% canonical_chromosomes(common_name), , drop = FALSE]
  if (!nrow(loc)) return(NULL)
  data.frame(gene = loc$SYMBOL, chr = loc$chr, start = loc$start,
             end = loc$end, strand = loc$strand, stringsAsFactors = FALSE)
}

# One horizontal track per chromosome, up-regulated genes drawn as stems above
# the baseline and down-regulated genes below it, stem height proportional to
# |log2FC|. `features` needs a gene and a log2 fold-change column.
chromosome_landscape <- function(features, orgdb, common_name = "mouse",
                                 cap = 2.5, title = NULL, subtitle = NULL,
                                 legend_title = "Direction") {
  features <- features[!is.na(features$gene) & !is.na(features$lfc), , drop = FALSE]
  features <- features[!duplicated(features$gene), , drop = FALSE]
  pos <- chromosome_positions(orgdb, features$gene, common_name)
  if (is.null(pos) || !nrow(pos)) return(NULL)
  pos <- merge(features[, c("gene", "lfc")], pos, by = "gene")
  pos <- pos[is.finite(pos$lfc) & pos$lfc != 0, , drop = FALSE]
  if (!nrow(pos)) return(NULL)

  chr_levels <- canonical_chromosomes(common_name)
  ref_len <- reference_chrom_lengths(common_name)
  chr_levels <- chr_levels[chr_levels %in% unique(pos$chr)]
  observed <- tapply(pos$end, pos$chr, max, na.rm = TRUE) / 1e6
  len <- as.numeric(ref_len[chr_levels])
  len[is.na(len)] <- as.numeric(observed[chr_levels][is.na(len)])

  base <- stats::setNames(-(seq_along(chr_levels) - 1L), chr_levels)
  # The subtitle carries the threshold and the direction rule; at the enlarged
  # font it is wider than the canvas, so break it into two lines.
  if (!is.null(subtitle)) subtitle <- paste(strwrap(subtitle, width = 92), collapse = "\n")
  cap <- if (is.finite(cap) && cap > 0) cap else max(abs(pos$lfc))
  pos$height <- 0.42 * pmin(abs(pos$lfc), cap) / cap
  pos$base <- base[pos$chr]
  pos$tip <- pos$base + ifelse(pos$lfc > 0, pos$height, -pos$height)
  pos$mb <- pos$start / 1e6
  pos$dir <- ifelse(pos$lfc > 0, "Up", "Down")

  tracks <- data.frame(chr = chr_levels, base = unname(base[chr_levels]),
                       length_mb = len, stringsAsFactors = FALSE)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = tracks,
                          ggplot2::aes(x = 0, xend = length_mb, y = base, yend = base),
                          linewidth = 2.6, colour = "grey88", lineend = "round") +
    ggplot2::geom_segment(data = pos,
                          ggplot2::aes(x = mb, xend = mb, y = base, yend = tip, colour = dir),
                          linewidth = 0.35, alpha = 0.85) +
    ggplot2::geom_point(data = pos,
                        ggplot2::aes(x = mb, y = tip, colour = dir), size = 0.7) +
    ggplot2::geom_text(data = tracks,
                       ggplot2::aes(x = length_mb, y = base, label = sprintf("%.0f Mb", length_mb)),
                       hjust = -0.08, size = 4.4, colour = "grey35") +
    ggplot2::scale_colour_manual(values = c(Up = "#B2182B", Down = "#2166AC"),
                                 name = legend_title) +
    ggplot2::scale_y_continuous(breaks = tracks$base, labels = tracks$chr,
                                expand = ggplot2::expansion(add = c(0.8, 0.8))) +
    ggplot2::scale_x_continuous(labels = function(x) x,
                                expand = ggplot2::expansion(mult = c(0.01, 0.06))) +
    ggplot2::labs(x = "Chromosome position (Mb)", y = NULL, title = title, subtitle = subtitle) +
    publication_theme(base_size = 16) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 13.5),
                   axis.text.y = ggplot2::element_text(face = "bold", size = 14),
                   axis.title = ggplot2::element_text(face = "bold", size = 16),
                   plot.title = ggplot2::element_text(face = "bold", size = 18.5),
                   plot.subtitle = ggplot2::element_text(size = 16),
                   legend.text = ggplot2::element_text(size = 13.5),
                   legend.title = ggplot2::element_text(face = "bold", size = 16),
                   legend.position = "top")
  list(plot = p, positions = pos[, c("gene", "chr", "start", "end", "strand", "lfc", "dir")],
       mapped = nrow(pos), unmapped = nrow(features) - nrow(pos),
       chromosomes = length(unique(pos$chr)))
}

resolve_project_path <- function(path, project_root) {
  if (is.data.frame(path)) return(path)
  if (is.null(path) || length(path) == 0L) return(NULL)
  if (!is.character(path) || length(path) != 1L) stop("A path must be a single character string.", call. = FALSE)
  if (is.na(path) || !nzchar(path)) return(NULL)
  path <- path.expand(path)
  if (grepl("^([A-Za-z]:[/\\\\]|/)", path)) return(normalizePath(path, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(project_root, path), winslash = "/", mustWork = FALSE)
}

resolve_config_paths <- function(cfg, project_root) {
  cfg$project$output_dir <- resolve_project_path(cfg$project$output_dir, project_root)
  cfg$input$counts_file <- resolve_project_path(cfg$input$counts_file, project_root)
  cfg$input$sample_metadata_file <- resolve_project_path(cfg$input$sample_metadata_file, project_root)
  cfg$input$tf_database_file <- resolve_project_path(cfg$input$tf_database_file, project_root)
  cfg$input$tf_dorothea_dir <- resolve_project_path(cfg$input$tf_dorothea_dir, project_root)
  cfg$input$expression_matrix_file <- resolve_project_path(cfg$input$expression_matrix_file, project_root)
  cfg$input$gene_annotation <- resolve_project_path(cfg$input$gene_annotation, project_root)
  # Kept on the object so modules that cache static reference data (the KEGG
  # pathway tables) know where the project lives, not where R was started.
  attr(cfg, "project_root") <- project_root
  cfg
}

# One word selects the annotation resources that match an organism and its
# identifier system, so `run_all(species = "mouse")` replaces the whole
# species block the configuration file used to carry.
species_preset <- function(name = "mouse") {
  presets <- list(
    mouse = list(common_name = "mouse",
                 name = "Mus musculus", latin = "mus musculus", orgdb = "org.Mm.eg.db",
                 kegg_code = "mmu", ensembl_prefix = "ENSMUSG", gene_id_type = "ENSEMBL",
                 symbol_type = "SYMBOL", entrez_type = "ENTREZID",
                 pathview_species = "mmu", kegg_key_type = "ncbi-geneid",
                 pathview_id_type = "ENTREZ", msigdb_db_species = "MM",
                 tf_dataset = "mm",
                 string_taxon = 10090L),
    human = list(common_name = "human",
                 name = "Homo sapiens", latin = "homo sapiens", orgdb = "org.Hs.eg.db",
                 kegg_code = "hsa", ensembl_prefix = "ENSG", gene_id_type = "ENSEMBL",
                 symbol_type = "SYMBOL", entrez_type = "ENTREZID",
                 pathview_species = "hsa", kegg_key_type = "ncbi-geneid",
                 pathview_id_type = "ENTREZ", msigdb_db_species = "HS",
                 tf_dataset = "hs",
                 string_taxon = 9606L),
    rat = list(common_name = "rat",
               name = "Rattus norvegicus", latin = "rattus norvegicus", orgdb = "org.Rn.eg.db",
               kegg_code = "rno", ensembl_prefix = "ENSRNOG", gene_id_type = "ENSEMBL",
               symbol_type = "SYMBOL", entrez_type = "ENTREZID",
               pathview_species = "rno", kegg_key_type = "ncbi-geneid",
               pathview_id_type = "ENTREZ", msigdb_db_species = "RN",
               # DoRothEA ships human and mouse collections only; the mouse
               # regulons are the ortholog reference for rat, whose gene symbols
               # follow the same spelling.
               tf_dataset = "mm", string_taxon = 10116L)
  )
  aliases <- c(mouse = "mouse", mm = "mouse", mus = "mouse", "mus musculus" = "mouse",
               human = "human", hs = "human", hsa = "human", homo = "human",
               "homo sapiens" = "human",
               rat = "rat", rn = "rat", rno = "rat", rattus = "rat", "rattus norvegicus" = "rat")
  key <- tolower(trimws(as.character(name)))[1L]
  hit <- aliases[[key]]
  if (is.null(hit)) {
    stop("Unknown species: ", name, ". Use one of: mouse, human, rat.", call. = FALSE)
  }
  presets[[hit]]
}

read_counts_workbook <- function(cfg) {
  path <- cfg$input$counts_file
  stop_if_missing(path, "counts matrix path")
  if (!file.exists(path)) stop("Counts matrix does not exist: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required for Excel input.", call. = FALSE)
    sheets <- readxl::excel_sheets(path)
    sheet <- cfg$input$counts_sheet %||% sheets[[1L]]
    if (!sheet %in% sheets) stop("Sheet not found: ", sheet, call. = FALSE)
    dat <- as.data.frame(readxl::read_excel(path, sheet = sheet, .name_repair = "unique"))
  } else if (ext %in% c("csv", "txt")) {
    sheet <- NA_character_
    dat <- if (ext == "csv") utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) else
      utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported count-matrix format. Use .xlsx, .xls, .csv, or .txt.", call. = FALSE)
  }
  gene_id <- cfg$input$gene_id_column
  symbol <- cfg$input$gene_symbol_column %||% ""
  count_cols <- unlist(cfg$input$count_columns %||% character())
  if (length(count_cols) == 0L) {
    # Prefer sample IDs from the metadata file. This avoids accidentally treating
    # numeric annotation columns (such as gene lengths) as samples.
    metadata_path <- cfg$input$sample_metadata_file
    metadata_ids <- character()
    if (!is.null(metadata_path)) {
      metadata_preview <- read_metadata_table(metadata_path)
      if ("sample_id" %in% names(metadata_preview)) metadata_ids <- as.character(metadata_preview$sample_id)
    }
    count_cols <- intersect(metadata_ids, names(dat))
    if (length(count_cols) < 2L) {
      excluded <- unique(c(gene_id, symbol, cfg$input$annotation_column,
                           "Exonic.gene.sizes", "gene_length", "GeneLength"))
      count_cols <- names(dat)[vapply(dat, is.numeric, logical(1)) & !names(dat) %in% excluded]
    }
    if (length(count_cols) < 2L) stop("Could not infer at least two sample columns; set input.count_columns explicitly.", call. = FALSE)
  }
  missing_cols <- setdiff(c(gene_id, count_cols), names(dat))
  if (length(missing_cols)) stop("Missing input columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  count_df <- dat[, count_cols, drop = FALSE]
  count_mat <- as.matrix(count_df)
  suppressWarnings(storage.mode(count_mat) <- "numeric")
  if (anyNA(count_mat)) stop("Count matrix contains non-numeric or missing values.", call. = FALSE)
  if (any(count_mat < 0) || any(abs(count_mat - round(count_mat)) > 1e-8)) {
    stop("DESeq2 input must be non-negative integer raw counts.", call. = FALSE)
  }
  row_ids <- as.character(dat[[gene_id]])
  if (anyNA(row_ids) || any(!nzchar(row_ids))) stop("Gene identifiers contain missing values.", call. = FALSE)
  if (anyDuplicated(row_ids)) {
    # Sum duplicate rows at the stable ID level; preserve the first annotation.
    sums <- rowsum(count_mat, group = row_ids, reorder = FALSE)
    first <- match(rownames(sums), row_ids)
    keep_annotation <- intersect(c(symbol, cfg$input$annotation_column), names(dat))
    dat2 <- dat[first, c(gene_id, keep_annotation), drop = FALSE]
    dat2[[gene_id]] <- rownames(sums)
    dat <- dat2
    count_mat <- sums
  }
  rownames(count_mat) <- as.character(dat[[gene_id]])
  colnames(count_mat) <- count_cols
  annotation_cols <- intersect(c(gene_id, symbol, cfg$input$annotation_column), names(dat))
  annotations <- dat[, annotation_cols, drop = FALSE]
  list(counts = count_mat, annotations = annotations, raw = dat, sheet = sheet)
}

# Read a TF-target table from Excel, CSV or tab-delimited text by extension.
read_tf_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required for Excel TF tables.", call. = FALSE)
    as.data.frame(readxl::read_excel(path, .name_repair = "unique"))
  } else if (ext == "csv") {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

# DoRothEA regulons ship with the `dorothea` package: dorothea_hs (human) and
# dorothea_mm (mouse; the same collection serves rat, whose gene symbols use the
# same spelling). The package data are used when the package is installed,
# otherwise the identical .rda objects are read from the copy bundled in
# data/dorothea/. Either way no TF table has to come from the input folder.
dorothea_dataset <- function(cfg) {
  ds <- tolower(cfg$input$tf_species %||% "mm")
  if (!ds %in% c("hs", "mm")) {
    stop("Unsupported DoRothEA dataset: ", ds, " (use 'hs' or 'mm').", call. = FALSE)
  }
  ds
}

dorothea_data_path <- function(cfg, ds) {
  dir <- cfg$input$tf_dorothea_dir %||% NULL
  if (is.null(dir) || !nzchar(dir)) return(NULL)
  path <- file.path(dir, paste0("dorothea_", ds, ".rda"))
  if (file.exists(path)) path else NULL
}

dorothea_available <- function(cfg) {
  ds <- tryCatch(dorothea_dataset(cfg), error = function(e) NA_character_)
  if (is.na(ds)) return(FALSE)
  requireNamespace("dorothea", quietly = TRUE) || !is.null(dorothea_data_path(cfg, ds))
}

read_dorothea_table <- function(cfg) {
  ds <- dorothea_dataset(cfg)
  obj <- paste0("dorothea_", ds)
  env <- new.env()
  from_package <- FALSE
  if (requireNamespace("dorothea", quietly = TRUE)) {
    from_package <- tryCatch({
      suppressWarnings(utils::data(list = obj, package = "dorothea", envir = env))
      exists(obj, envir = env, inherits = FALSE)
    }, error = function(e) FALSE)
  }
  if (!from_package) {
    path <- dorothea_data_path(cfg, ds)
    if (is.null(path)) {
      stop("DoRothEA data for '", ds, "' not found.\n",
           "  Install the package (BiocManager::install(\"dorothea\")) or point\n",
           "  input.tf_dorothea_dir at the folder holding dorothea_", ds, ".rda.", call. = FALSE)
    }
    load(path, envir = env)
    if (!exists(obj, envir = env, inherits = FALSE)) {
      stop("Bundled DoRothEA file does not contain ", obj, ": ", path, call. = FALSE)
    }
  }
  df <- as.data.frame(env[[obj]], stringsAsFactors = FALSE)
  attr(df, "dataset") <- obj
  attr(df, "from_package") <- from_package
  df
}

read_metadata_table <- function(path) {
  if (is.data.frame(path)) return(as.data.frame(path, stringsAsFactors = FALSE))
  stop_if_missing(path, "sample metadata path")
  if (!file.exists(path)) stop("Sample metadata file does not exist: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' is required for Excel metadata.", call. = FALSE)
    return(as.data.frame(readxl::read_excel(path, .name_repair = "unique")))
  }
  if (ext == "csv") return(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  if (ext %in% c("tsv", "txt")) return(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE))
  stop("Unsupported metadata format. Use .csv, .xlsx, .xls, .tsv or .txt.", call. = FALSE)
}

read_sample_metadata <- function(cfg, sample_names) {
  meta <- read_metadata_table(cfg$input$sample_metadata_file)
  stop_if_missing(meta$sample_id, "sample_id metadata")
  meta$sample_id <- trimws(as.character(meta$sample_id))
  stop_if_missing(meta$sample_id, "sample_id metadata")
  if (anyDuplicated(meta$sample_id)) stop("sample_id values must be unique.", call. = FALSE)
  if (!setequal(meta$sample_id, sample_names)) {
    stop("Metadata sample_id values do not exactly match count-matrix columns.", call. = FALSE)
  }
  meta <- meta[match(sample_names, meta$sample_id), , drop = FALSE]
  group_col <- cfg$analysis$group_column %||% "condition"
  design_vars <- all.vars(stats::as.formula(cfg$analysis$design))
  missing_design <- setdiff(design_vars, names(meta))
  if (length(missing_design)) stop("Metadata lacks design variables: ", paste(missing_design, collapse = ", "), call. = FALSE)
  for (v in design_vars) meta[[v]] <- as.factor(meta[[v]])
  if (!group_col %in% names(meta)) {
    stop("Metadata lacks the configured group column: ", group_col, call. = FALSE)
  }
  if (nlevels(meta[[group_col]]) != 2L) {
    stop("This two-group workflow requires exactly two levels in '", group_col, "'.", call. = FALSE)
  }
  ref <- cfg$analysis$reference_condition %||% levels(meta[[group_col]])[1L]
  if (!ref %in% levels(meta[[group_col]])) {
    stop("reference_condition is not a level of '", group_col, "': ", ref, call. = FALSE)
  }
  meta[[group_col]] <- stats::relevel(meta[[group_col]], ref = ref)
  if (any(table(meta[[group_col]]) < (cfg$analysis$min_replicates_per_group %||% 2L))) {
    stop("Each group must have at least min_replicates_per_group biological replicates.", call. = FALSE)
  }
  rownames(meta) <- meta$sample_id
  meta
}

make_unique_ranked_vector <- function(symbols, values) {
  ok <- !is.na(symbols) & nzchar(symbols) & is.finite(values)
  x <- data.frame(id = as.character(symbols[ok]), value = as.numeric(values[ok]))
  x <- x[order(abs(x$value), decreasing = TRUE), , drop = FALSE]
  x <- x[!duplicated(x$id), , drop = FALSE]
  out <- x$value
  names(out) <- x$id
  sort(out, decreasing = TRUE)
}

result_to_df <- function(res, annotations, gene_id_col, symbol_col) {
  d <- as.data.frame(res)
  d[[gene_id_col]] <- rownames(d)
  ann <- annotations[, intersect(c(gene_id_col, symbol_col), names(annotations)), drop = FALSE]
  d <- merge(d, ann, by = gene_id_col, all.x = TRUE, sort = FALSE)
  d
}

p_adjusted_or_one <- function(x) {
  y <- as.numeric(x)
  y[is.na(y)] <- 1
  y[y <= 0] <- .Machine$double.xmin
  -log10(pmin(y, 1))
}

# Where redrawable ggplot objects are cached. scripts/redraw_figures.R reads
# this directory and re-saves the same plots at whatever size the registry in
# pipeline/R/figure_sizes.R currently says, without touching the analysis.
figure_cache_dir <- function() {
  dir <- Sys.getenv("CLAWSGO_FIGURE_CACHE", "cache/figures")
  dir <- resolve_project_path(dir, getOption("transatlas.project_root", getwd()))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir
}

# A ggplot object carries the environment it was built in. Inside a section
# that environment holds the full context (DESeq2 object, count matrices,
# regulons), so serialising it would write hundreds of megabytes and dominate
# the run time. The data are already inside the plot object, so the environment
# can point at the global environment, which serialises as a reference.
strip_plot_env <- function(plot) {
  # A layer mapping is a list of quosures whose environments point back at the
  # section function. When those are serialized, every object of the section
  # (the full DESeq2 result, the VST matrix, ...) travels with the figure --
  # a single 6x6 inch panel reached 55 MB. The mappings only need the layer
  # data, so the environments can safely be repointed at the global env.
  strip_quosures <- function(x) {
    if (rlang::is_quosure(x)) return(rlang::quo_set_env(x, globalenv()))
    if (is.list(x) && !is.data.frame(x)) return(lapply(x, strip_quosures))
    x
  }
  # ggplot2 3.5 stores mappings as S7 objects; `m[] <-` keeps that class,
  # whereas rebuilding the list would fail the object validation.
  strip_mapping <- function(m) {
    if (is.null(m)) return(NULL)
    m[] <- lapply(m, strip_quosures)
    m
  }
  if (inherits(plot, "ggplot")) {
    plot$plot_env <- globalenv()
    plot$mapping <- strip_mapping(plot$mapping)
    if (length(plot$layers)) {
      plot$layers <- lapply(plot$layers, function(layer) {
        if (!is.null(layer$mapping)) layer$mapping <- strip_mapping(layer$mapping)
        layer
      })
    }
  }
  if (inherits(plot, "patchwork") && !is.null(plot$patches$plots)) {
    plot$patches$plots <- lapply(plot$patches$plots, strip_plot_env)
  }
  plot
}

cache_figure <- function(key, plot, path, size_args = list()) {
  rds <- file.path(figure_cache_dir(),
                   paste0(tools::file_path_sans_ext(basename(path)), ".rds"))
  ok <- try(saveRDS(list(key = key, plot = strip_plot_env(plot),
                         path = normalizePath(path, mustWork = FALSE),
                         size_args = size_args), rds,
                    # gzip keeps the cache write fast (xz spends ~50 s on the
                    # 20,000-point volcano alone); the files are larger but read
                    # back in a fraction of the time.
                    compress = "gzip"), silent = TRUE)
  if (!inherits(ok, "try-error")) index_cached_figure(key, rds, path, size_args)
  invisible(NULL)
}

# A one-row-per-figure index so `redraw_figures.R --list` and single-key
# redraws never have to deserialize every cached plot (the objects are large
# and reading them all takes minutes).
figure_index_file <- function() file.path(figure_cache_dir(), "index.csv")

index_cached_figure <- function(key, rds, path, size_args = list()) {
  idx_file <- figure_index_file()
  old <- if (file.exists(idx_file)) {
    try(utils::read.csv(idx_file, stringsAsFactors = FALSE), silent = TRUE)
  } else NULL
  if (inherits(old, "try-error")) old <- NULL
  args_txt <- if (length(size_args)) {
    paste(names(size_args), vapply(size_args, function(v) paste(format(v), collapse = "|"),
                                   character(1)), sep = "=", collapse = ";")
  } else ""
  new <- data.frame(key = key %||% "",
                    rds = basename(rds),
                    path = normalizePath(path, mustWork = FALSE),
                    size_args = args_txt,
                    stringsAsFactors = FALSE)
  out <- if (is.null(old) || !nrow(old)) new else {
    old <- old[old$rds != new$rds, , drop = FALSE]
    rbind(old, new)
  }
  out <- out[order(out$rds), , drop = FALSE]
  utils::write.csv(out, idx_file, row.names = FALSE)
  invisible(NULL)
}

parse_size_args <- function(txt) {
  txt <- txt %||% ""
  if (is.na(txt) || !nzchar(txt)) return(list())
  parts <- strsplit(txt, ";", fixed = TRUE)[[1]]
  out <- lapply(parts, function(p) {
    kv <- strsplit(p, "=", fixed = TRUE)[[1]]
    v <- strsplit(kv[2], "|", fixed = TRUE)[[1]]
    suppressWarnings(as.numeric(v))
  })
  names(out) <- vapply(strsplit(parts, "=", fixed = TRUE), `[`, character(1), 1)
  out
}

# The size of every figure comes from FIGURE_SIZES in pipeline/R/figure_sizes.R.
# `key` names the entry; `...` is forwarded to function entries whose height
# depends on the number of rows. Numeric width/height arguments are still
# accepted so a one-off call can override the registry, but section code passes
# a key and nothing else.
save_plot <- function(plot, path, width = NULL, height = NULL, key = NULL,
                      dpi = 320, ...) {
  size_args <- list(...)
  if (is.null(width) || is.null(height)) {
    sz <- do.call(fig_size, c(list(key), size_args))
    if (is.null(width)) width <- sz[[1L]]
    if (is.null(height)) height <- sz[[2L]]
  }
  if (!is.null(key)) record_figure_size(key, c(width, height), path)
  ggplot2::ggsave(path, plot = plot, width = width, height = height,
                  units = "in", dpi = dpi, bg = "white")
  if (!is.null(key)) cache_figure(key, plot, path, size_args)
  invisible(path)
}

# Display clipping only. Calculation tables and ranking vectors keep true effects.
display_lfc_limit <- function(cfg, values = numeric()) {
  cap <- cfg$analysis$up_cutoff
  if (is.null(cap) && "log2fc_cap" %in% names(cfg$analysis)) cap <- cfg$analysis$log2fc_cap
  if (!is.null(cap)) return(as.numeric(cap))
  max(c(abs(values[is.finite(values)]), 0.5))
}

cap_log2fc <- function(x, cap) {
  if (is.null(cap) || length(cap) == 0L) return(x)
  cap <- suppressWarnings(as.numeric(cap))
  if (!is.finite(cap) || cap <= 0) return(x)
  pmax(pmin(x, cap), -cap)
}

# The number of features a heat map or dot plot may draw. Kept in one place so
# every figure in the pipeline respects the same ceiling.
max_display_features <- function(cfg, default = 30L) {
  as.integer(cfg$analysis$max_display_features %||% default)
}

# Pick gene sets for a composite GSEA panel so that both directions are shown.
#
# Taking the N most significant sets by adjusted P alone fills the panel with
# whichever direction happens to dominate - on a strongly proliferative contrast
# every panel came out all-up. Half the slots go to the strongest positively
# enriched sets and half to the strongest negatively enriched ones, so a reader
# sees both the induced and the repressed programme in one figure.
#
# When `gene_sets` is supplied, a set is also skipped when more than
# `max_overlap` of the genes of the smaller set are shared with a set already
# chosen, so the panel shows distinct pathways rather than the same programme
# under several names.
balanced_set_selection <- function(df, n_total = 9L, p_col = "p.adjust", nes_col = "NES",
                                   gene_sets = NULL, max_overlap = 0.4) {
  if (is.null(df) || nrow(df) == 0L) return(integer(0))
  id <- as.character(if ("ID" %in% names(df)) df$ID else rownames(df))
  p <- suppressWarnings(as.numeric(df[[p_col]]))
  nes <- suppressWarnings(as.numeric(df[[nes_col]]))
  if (all(is.na(nes))) return(utils::head(order(p, na.last = TRUE), n_total))
  n_total <- min(as.integer(n_total), nrow(df))
  too_similar <- function(candidate, chosen) {
    if (is.null(gene_sets) || !length(chosen)) return(FALSE)
    a <- gene_sets[[candidate]]
    if (is.null(a) || !length(a)) return(FALSE)
    any(vapply(chosen, function(s) {
      b <- gene_sets[[s]]
      if (is.null(b) || !length(b)) return(FALSE)
      length(intersect(a, b)) / min(length(a), length(b)) > max_overlap
    }, logical(1)))
  }
  n_up <- ceiling(n_total / 2)
  n_down <- n_total - n_up
  take <- function(sign, k, chosen = character(0)) {
    if (k <= 0L) return(character(0))
    ok <- which(!is.na(nes) & sign * nes > 0)
    if (!length(ok)) return(character(0))
    out <- character(0)
    for (i in ok[order(p[ok], na.last = TRUE)]) {
      if (length(out) >= k) break
      if (too_similar(id[i], c(chosen, out))) next
      out <- c(out, id[i])
    }
    out
  }
  up <- take(1, n_up)
  down <- take(-1, n_down, chosen = up)
  # A one-sided result still fills the panel from the other side, so the figure
  # never comes back with fewer sets than were asked for when it can avoid it.
  short <- n_total - length(c(up, down))
  if (short > 0L) {
    chosen <- c(up, down)
    rest <- setdiff(id, chosen)
    extra <- character(0)
    for (i in rest[order(p[match(rest, id)], na.last = TRUE)]) {
      if (length(extra) >= short) break
      if (too_similar(i, c(chosen, extra))) next
      extra <- c(extra, i)
    }
    c(up, down, extra)
  } else {
    c(up, down)
  }
}

# Composite GSEA panel selection. The panel is meant to show distinct
# biological programmes, so a set is skipped when more than `max_overlap` of the
# smaller set's genes are shared with a set already chosen. Within each
# direction, sets are ranked by the strength of the enrichment score (largest
# positive NES first, most negative NES first), and only sets with more than
# `min_size` genes compete. The number of up and down panels follows the
# positive/negative NES counts, so a strongly one-sided result is not forced
# into a misleading half-and-half layout; one slot is still reserved for each
# direction when both are present.
diverse_gsea_selection <- function(df, gene_sets = NULL, n_total = 9L,
                                   min_size = 50L, max_overlap = 0.4,
                                   size_col = "setSize", nes_col = "NES",
                                   id_col = "ID") {
  if (is.null(df) || !nrow(df)) return(character(0))
  id <- as.character(if (id_col %in% names(df)) df[[id_col]] else rownames(df))
  nes <- suppressWarnings(as.numeric(df[[nes_col]]))
  size <- if (size_col %in% names(df)) {
    suppressWarnings(as.numeric(df[[size_col]]))
  } else {
    rep(NA_real_, length(id))
  }
  keep <- !is.na(nes)
  if (any(!is.na(size))) keep <- keep & (is.na(size) | size > min_size)
  if (sum(keep) < 2L) keep <- !is.na(nes)
  id <- id[keep]
  nes <- nes[keep]
  if (!length(id)) return(character(0))
  n_total <- min(as.integer(n_total), length(id))
  n_pos <- sum(nes > 0)
  n_neg <- sum(nes < 0)
  if (n_pos == 0L && n_neg == 0L) return(utils::head(id, n_total))

  n_up <- round(n_total * n_pos / max(1L, n_pos + n_neg))
  if (n_pos > 0L && n_neg > 0L) n_up <- min(max(n_up, 1L), n_total - 1L)
  n_up <- min(n_up, n_pos)
  n_down <- min(n_total - n_up, n_neg)

  too_similar <- function(candidate, chosen) {
    if (is.null(gene_sets) || !length(chosen)) return(FALSE)
    a <- gene_sets[[candidate]]
    if (is.null(a) || !length(a)) return(FALSE)
    any(vapply(chosen, function(s) {
      b <- gene_sets[[s]]
      if (is.null(b) || !length(b)) return(FALSE)
      length(intersect(a, b)) / min(length(a), length(b)) > max_overlap
    }, logical(1)))
  }
  pick <- function(candidates, k, chosen = character(0)) {
    out <- character(0)
    for (g in candidates) {
      if (length(out) >= k) break
      if (too_similar(g, c(chosen, out))) next
      out <- c(out, g)
    }
    out
  }
  up_ids <- id[nes > 0][order(-nes[nes > 0])]
  down_ids <- id[nes < 0][order(nes[nes < 0])]
  up <- pick(up_ids, n_up)
  down <- pick(down_ids, n_down, chosen = up)
  chosen <- c(up, down)
  short <- n_total - length(chosen)
  if (short > 0L) {
    rest <- c(setdiff(up_ids, chosen), setdiff(down_ids, chosen))
    chosen <- c(chosen, pick(rest, short, chosen = chosen))
  }
  chosen
}

# Pick the terms shown on a GO over-representation dot plot. `enrichGO` keeps
# parent and child terms side by side, so the first 20 rows by adjusted P are
# often two or three variants of the same programme. Terms are read in the order
# given (the result is already sorted by adjusted P) and a term is skipped when
# more than `max_overlap` of the genes of the smaller term are shared with a
# term already chosen. Returns the term IDs, in selection order, so the caller
# can pass them to enrichplot::dotplot(showCategory = ...).
diverse_term_selection <- function(ids, gene_sets, n_total = 20L, max_overlap = 0.4) {
  ids <- as.character(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) return(character(0))
  n_total <- min(as.integer(n_total), length(ids))
  out <- character(0)
  for (g in ids) {
    if (length(out) >= n_total) break
    a <- gene_sets[[g]]
    shared <- FALSE
    if (!is.null(a) && length(a) && length(out)) {
      shared <- any(vapply(out, function(s) {
        b <- gene_sets[[s]]
        if (is.null(b) || !length(b)) return(FALSE)
        length(intersect(a, b)) / min(length(a), length(b)) > max_overlap
      }, logical(1)))
    }
    if (!shared) out <- c(out, g)
  }
  # A very redundant ontology can yield fewer than `n_total` distinct terms;
  # the panel is then drawn with those rather than topped up with the
  # near-duplicates the rule just removed.
  out
}
# redundant set from every cluster for the enrichment network. Sets are
# clustered on Jaccard distance of their gene memberships (average linkage,
# at most `max_clusters` groups); inside a group the most significant sets are
# taken first, and a set is skipped when more than `max_overlap` of the smaller
# set is shared with one already picked for that group. Clustering first keeps
# one dominant programme (cell cycle, ribosome biogenesis) from filling the
# whole network; the overlap rule then removes near-duplicates within a group.
cluster_gsea_terms <- function(df, gene_sets, max_clusters = 12L,
                               per_cluster = 10L, max_overlap = 0.6,
                               id_col = "ID") {
  if (is.null(df) || !nrow(df) || is.null(gene_sets) || !length(gene_sets)) return(NULL)
  id <- as.character(if (id_col %in% names(df)) df[[id_col]] else rownames(df))
  n <- length(id)
  members <- lapply(id, function(x) {
    m <- gene_sets[[x]]
    if (is.null(m)) character(0) else unique(as.character(m))
  })
  names(members) <- id
  sim <- matrix(0, n, n, dimnames = list(id, id))
  for (i in seq_len(n)) {
    a <- members[[i]]
    if (!length(a)) next
    for (j in seq_len(i)) {
      b <- members[[j]]
      if (!length(b)) next
      sim[i, j] <- sim[j, i] <- length(intersect(a, b)) / length(union(a, b))
    }
  }
  k <- min(as.integer(max_clusters), n)
  cluster_raw <- if (k >= 2L) {
    hc <- stats::hclust(stats::as.dist(1 - sim), method = "average")
    stats::cutree(hc, k = k)
  } else {
    rep(1L, n)
  }
  p_adj <- suppressWarnings(as.numeric(df$p.adjust))
  p_val <- suppressWarnings(as.numeric(df$pvalue))
  if (!length(p_adj)) p_adj <- rep(NA_real_, n)
  if (!length(p_val)) p_val <- rep(NA_real_, n)
  ord <- order(p_adj, p_val, na.last = TRUE)
  # Number the groups C1.. by the best adjusted P of any member, so the legend
  # follows the strength of the evidence rather than cutree's internal order.
  first_pos <- vapply(split(seq_along(ord), cluster_raw[ord]), min, numeric(1))
  cluster_order <- as.integer(names(first_pos)[order(first_pos)])
  cluster <- paste0("C", match(cluster_raw, cluster_order))

  chosen <- character(0)
  rows <- list()
  for (cl in unique(cluster[ord])) {
    idx <- ord[cluster[ord] == cl]
    picked <- character(0)
    for (i in idx) {
      if (length(picked) >= per_cluster) break
      a <- members[[id[i]]]
      if (length(a) < 1L) next
      ov <- if (length(picked)) {
        max(vapply(picked, function(s) {
          length(intersect(a, members[[s]])) / min(length(a), length(members[[s]]))
        }, numeric(1)))
      } else 0
      if (ov > max_overlap) next
      picked <- c(picked, id[i])
      rows[[length(rows) + 1L]] <- data.frame(
        ID = id[i], Description = as.character(df$Description[i]),
        cluster = cl, rank_in_cluster = length(picked),
        setSize = if ("setSize" %in% names(df)) df$setSize[i] else NA_integer_,
        NES = if ("NES" %in% names(df)) suppressWarnings(as.numeric(df$NES[i])) else NA_real_,
        p.adjust = p_adj[i], pvalue = p_val[i],
        max_overlap_previous = ov, stringsAsFactors = FALSE)
    }
    chosen <- c(chosen, picked)
  }
  if (!length(chosen)) return(NULL)
  selected <- do.call(rbind, rows)
  list(selected = selected, members = members[chosen],
       sim = sim[chosen, chosen, drop = FALSE],
       n_clusters = length(unique(selected$cluster)))
}

# Resolve a sans-serif family that is actually installed on this machine.
# Arial is preferred for journal submissions, but plots must still render on
# hosts where it is absent; Liberation Sans is metric-compatible with Arial.
publication_font <- function() {
  cached <- getOption("clawsgo.publication_font")
  if (!is.null(cached)) return(cached)
  candidates <- c("Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans", "sans")
  available <- tryCatch({
    if (requireNamespace("systemfonts", quietly = TRUE)) unique(systemfonts::system_fonts()$family) else character()
  }, error = function(e) character())
  resolved <- if (length(available)) {
    hit <- candidates[candidates %in% available]
    if (length(hit)) hit[[1L]] else "sans"
  } else "sans"
  options(clawsgo.publication_font = resolved)
  resolved
}

publication_theme <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size, base_family = publication_font()) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

safe_write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

format_number <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

# clusterProfiler's KEGG interface does not accept the OrgDb key type directly:# Entrez identifiers must be requested as "ncbi-geneid", otherwise every KEGG
# query fails at the REST URL level. Pathview uses a separate identifier vocabulary.
kegg_key_type <- function(cfg, key_type) {
  explicit <- cfg$species$kegg_key_type
  if (!is.null(explicit) && nzchar(explicit)) return(explicit)
  if (identical(key_type, "ENTREZID")) "ncbi-geneid" else key_type
}

pathview_id_type <- function(cfg, key_type) {
  explicit <- cfg$species$pathview_id_type
  if (!is.null(explicit) && nzchar(explicit)) return(explicit)
  switch(key_type,
         ENTREZID = "ENTREZ",
         ENSEMBL = "ENSEMBL",
         UNIPROT = "UNIPROT",
         SYMBOL = "SYMBOL",
         key_type)
}

# Map the configured symbol column onto the configured annotation key type,
# returning a unique symbol -> key table so downstream modules can look up
# identifiers without depending on the enrichment module's intermediate state.
symbol_to_key_map <- function(symbols, orgdb, symbol_type, key_type) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  if (length(symbols) == 0L) return(NULL)
  map <- tryCatch(clusterProfiler::bitr(symbols, fromType = symbol_type, toType = key_type,
                                        OrgDb = orgdb, drop = FALSE), error = record_optional_error)
  if (is.null(map) || nrow(map) == 0L) return(NULL)
  map[!duplicated(map[[symbol_type]]), , drop = FALSE]
}

# msigdbr renamed `category`/`subcategory` to `collection`/`subcollection` in
# version 10; support both so the module runs against pinned older versions.
msigdbr_gene_sets <- function(species, collection, subcollection = NULL, db_species = NULL) {
  args <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% args) {
    if (!is.null(db_species) && nzchar(db_species)) {
      msigdbr::msigdbr(db_species = db_species, species = species,
                       collection = collection, subcollection = subcollection)
    } else {
      msigdbr::msigdbr(species = species, collection = collection, subcollection = subcollection)
    }
  } else {
    msigdbr::msigdbr(species = species, category = collection, subcategory = subcollection)
  }
}

# Keep the requested collection exact: H and mouse-native MH are distinct.
resolve_msigdb_gene_sets <- function(species, collection, subcollection = NULL,
                                      db_species = NULL, provenance_file = NULL) {
  modern <- "collection" %in% names(formals(msigdbr::msigdbr))
  requested_db <- db_species %||% "HS"
  candidates <- if (modern) unique(c(requested_db, "HS")) else "HS"
  candidates <- candidates[candidates %in% c("HS", "MM")]
  attempts <- list()
  result <- NULL
  selected_db <- NULL
  for (db in candidates) {
    available <- tryCatch({
      args <- if ("db_species" %in% names(formals(msigdbr::msigdbr_collections))) list(db_species = db) else list()
      do.call(msigdbr::msigdbr_collections, args)
    }, error = function(e) {
      warning("Cannot inspect MSigDB collections for ", db, ": ", conditionMessage(e), call. = FALSE)
      NULL
    })
    if (!is.null(available)) {
      collection_col <- if ("gs_collection" %in% names(available)) "gs_collection" else "gs_cat"
      sub_col <- if ("gs_subcollection" %in% names(available)) "gs_subcollection" else "gs_subcat"
      keep <- available[[collection_col]] == collection
      if (!is.null(subcollection)) keep <- keep & available[[sub_col]] == subcollection
      if (!any(keep)) {
        attempts[[db]] <- "Requested collection/subcollection is unavailable"
        next
      }
    }
    result <- tryCatch(msigdbr_gene_sets(species, collection, subcollection,
                                         if (modern) db else NULL), error = function(e) {
      attempts[[db]] <<- conditionMessage(e)
      warning("MSigDB ", db, " lookup failed: ", conditionMessage(e), call. = FALSE)
      NULL
    })
    if (!is.null(result) && nrow(result)) { selected_db <- db; break }
  }
  if (is.null(result) || !nrow(result)) {
    if (!is.null(provenance_file)) jsonlite::write_json(
      list(status = "unavailable", requested_collection = collection,
           requested_subcollection = subcollection, requested_database = requested_db,
           target_species = species, attempts = attempts), provenance_file,
      auto_unbox = TRUE, pretty = TRUE, null = "null")
    stop("No MSigDB gene sets available for collection ", collection,
         if (!is.null(subcollection)) paste0(" / ", subcollection),
         ". Mouse-native Hallmark uses MH; human Hallmark uses H.", call. = FALSE)
  }
  membership <- unique(result[, c("gs_name", "gene_symbol"), drop = FALSE])
  membership <- membership[order(membership$gs_name, membership$gene_symbol), , drop = FALSE]
  hash_file <- tempfile(fileext = ".csv")
  on.exit(unlink(hash_file), add = TRUE)
  utils::write.csv(membership, hash_file, row.names = FALSE, na = "")
  native_species <- if (selected_db == "MM") "Mus musculus" else "Homo sapiens"
  provenance <- list(
    status = "completed", requested_collection = collection,
    requested_subcollection = subcollection, requested_database = requested_db,
    database = selected_db,
    database_version = if ("db_version" %in% names(result)) unique(result$db_version) else "unavailable",
    collection = collection, subcollection = subcollection,
    database_species = native_species, target_species = species,
    ortholog_mapping = !identical(species, native_species),
    mapping_policy = if (identical(species, native_species)) "native identifiers" else "msigdbr ortholog mapping",
    package_version = as.character(utils::packageVersion("msigdbr")),
    gene_sets = length(unique(result$gs_name)), genes = length(unique(result$gene_symbol)),
    membership_rows = nrow(membership), membership_md5 = unname(tools::md5sum(hash_file)),
    attempts = attempts
  )
  if (!is.null(provenance_file)) jsonlite::write_json(provenance, provenance_file,
    auto_unbox = TRUE, pretty = TRUE, null = "null")
  attr(result, "msigdb_provenance") <- provenance
  result
}

# Preserve descriptive results when a Welch test is undefined.
differential_activity <- function(activity_long, group_of, groups) {
  rows <- lapply(split(activity_long, activity_long$source), function(x) {
    group <- unname(group_of[as.character(x$condition)])
    a <- x$score[!is.na(group) & group == groups[1L]]
    b <- x$score[!is.na(group) & group == groups[2L]]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    ma <- if (length(a)) mean(a) else NA_real_
    mb <- if (length(b)) mean(b) else NA_real_
    row <- data.frame(source = as.character(x$source[1L]), mean_reference = ma,
                       mean_treatment = mb, difference = mb - ma,
                       statistic = NA_real_, p_value = NA_real_, n_reference = length(a),
                       n_treatment = length(b), test_status = "not_testable",
                       test_reason = "", stringsAsFactors = FALSE)
    if (length(a) < 2L || length(b) < 2L) {
      row$test_reason <- "Fewer than two finite observations in one or both groups"
      return(row)
    }
    se <- sqrt(stats::var(a) / length(a) + stats::var(b) / length(b))
    if (!is.finite(se) || se == 0 || se < 10 * .Machine$double.eps * max(abs(ma), abs(mb))) {
      row$test_reason <- "Within-group variation is zero or numerically negligible; Welch test is undefined"
      return(row)
    }
    test <- tryCatch(stats::t.test(b, a), error = function(e) {
      record_optional_error(e)
      e
    })
    if (inherits(test, "error")) {
      row$test_status <- "error"; row$test_reason <- conditionMessage(test)
    } else if (!is.finite(test$p.value)) {
      row$test_reason <- "Welch test returned a non-finite P value"
    } else {
      row$statistic <- unname(test$statistic); row$p_value <- test$p.value
      row$test_status <- "tested"
    }
    row
  })
  if (!length(rows)) return(NULL)
  df <- do.call(rbind, rows)
  df$padj <- NA_real_
  tested <- df$test_status == "tested" & is.finite(df$p_value)
  df$padj[tested] <- stats::p.adjust(df$p_value[tested], method = "BH")
  df[order(df$p_value, na.last = TRUE), , drop = FALSE]
}

# A palette wide enough for the functional classes and stable between figures,
# so the same class keeps the same colour in every panel of a report.
class_palette <- function(n) {
  # Saturated qualitative palette: the previous Tableau-like set had two muted
  # earth tones and a grey, which read as one another in the class legend.
  base <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
            "#00CED1", "#F781BF", "#A65628", "#7F7F00", "#1B9E77",
            "#D37295", "#A0CBE8", "#FFBE7D", "#8CD17D", "#B6992D",
            "#499894", "#FABFD2", "#86BCB6")
  rep(base, length.out = max(1L, n))
}

# Diverging ramp shared by every expression heat map: purple for the low end
# through a light middle to orange for the high end (ColorBrewer PuOr). It
# replaces the green-white-blue ramp, whose two ends both read as "blue-green"
# once pheatmap stretches the row z-scores over the whole cell range.
heatmap_palette <- function(n = 100L) {
  grDevices::colorRampPalette(c("#542788", "#8073AC", "#B2ABD2", "#D8DAEB", "#F7F7F7",
                                "#FEE0B6", "#FDB863", "#E08214", "#B35806"))(n)
}

# Collapse a vector of term names into legend-sized labels.
wrap_label <- function(x, width = 42L) {
  vapply(x, function(s) paste(strwrap(as.character(s), width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

# Group the enriched terms of one ontology into functional classes and draw the
# two network figures on the user's own ggraph idiom.
#
# enrichplot::emapplot gave every term the same visual weight and let the labels
# collide; enrichplot::cnetplot left the gene nodes grey, so the effect direction
# of the genes - the point of the figure - was invisible. Both are replaced here.
# The terms are clustered on Wang semantic similarity into ten classes, each
# class gets one colour, only the class representative is labelled, and every
# gene is coloured by its log2 fold change with the green-white-purple ramp used
# by the user's own NetForTFs.R.
# Identifier translation used by the enrichment tables: Entrez (or whichever
# key type the organism uses) -> gene symbol.
key_to_symbol_map <- function(gene_map, key_type, symbol_type) {
  if (is.null(gene_map) || is.null(key_type) || is.null(symbol_type)) return(NULL)
  if (!all(c(key_type, symbol_type) %in% names(gene_map))) return(NULL)
  key <- as.character(gene_map[[key_type]])
  sym <- as.character(gene_map[[symbol_type]])
  ok <- !is.na(key) & !is.na(sym) & nzchar(key) & nzchar(sym)
  map <- stats::setNames(sym[ok], key[ok])
  map[!duplicated(names(map))]
}

# Convert a "/"-separated id string (GSEA core_enrichment, GO geneID) to symbols.
# Ids without a mapping are kept as they are, so nothing is silently dropped.
ids_to_symbols <- function(ids, map, sep = "/") {
  if (is.null(map) || !length(map)) return(ids)
  vapply(as.character(ids), function(x) {
    if (is.na(x) || !nzchar(x)) return(x)
    parts <- strsplit(x, sep, fixed = TRUE)[[1L]]
    hit <- unname(map[parts])
    parts[!is.na(hit)] <- hit[!is.na(hit)]
    paste(parts, collapse = sep)
  }, character(1), USE.NAMES = FALSE)
}

go_functional_classes <- function(ego, cfg, fc_symbols = NULL, n_terms = 100L, n_classes = 10L,
                                  sim_cut = 0.35, max_genes = 100L, contrast = NULL,
                                  font = publication_font()) {
  if (!requireNamespace("ggraph", quietly = TRUE) ||
      !requireNamespace("tidygraph", quietly = TRUE)) return(NULL)

  df <- as.data.frame(ego)
  df <- df[!is.na(df$p.adjust) & !is.na(df$Description), , drop = FALSE]
  if (!nrow(df)) return(NULL)
  df <- df[order(df$p.adjust), , drop = FALSE]
  df <- utils::head(df[!duplicated(df$Description), , drop = FALSE], n_terms)
  n_classes <- as.integer(min(n_classes, max(2L, nrow(df) %/% 3L)))
  if (nrow(df) < n_classes + 2L) return(NULL)

  ts <- tryCatch(enrichplot::pairwise_termsim(ego), error = record_optional_error)
  sim <- if (!is.null(ts)) ts@termsim else NULL
  if (is.null(sim)) return(NULL)
  keep <- intersect(df$Description, rownames(sim))
  if (length(keep) < n_classes + 2L) return(NULL)
  m <- sim[keep, keep, drop = FALSE]
  m[is.na(m)] <- 0
  df <- df[match(keep, df$Description), , drop = FALSE]

  # K-means on the classical-MDS embedding of the semantic distances. Cutting the
  # dendrogram instead produced one 60-term cell-cycle class plus a tail of
  # singletons, which is not what "a dozen functional classes" should look like.
  emb <- tryCatch(stats::cmdscale(stats::as.dist(1 - m), k = min(8L, nrow(m) - 1L)),
                  error = record_optional_error)
  if (is.null(emb)) return(NULL)
  set.seed(as.integer(cfg$project$seed %||% 1L))
  km <- tryCatch(stats::kmeans(emb, centers = n_classes, nstart = 50L), error = record_optional_error)
  if (is.null(km)) return(NULL)
  cls <- as.integer(km$cluster)

  # Order the classes by size, largest first, so the legend reads sensibly.
  sizes <- vapply(sort(unique(cls)), function(g) sum(cls == g), numeric(1))
  order_ids <- sort(unique(cls))[order(-sizes)]
  cls <- match(cls, order_ids)
  label_of <- vapply(seq_along(order_ids), function(g) {
    i <- which(cls == g)
    df$Description[i][which.min(df$p.adjust[i])]
  }, character(1))
  pal <- class_palette(length(order_ids))
  names(pal) <- as.character(seq_along(order_ids))
  legend_lab <- stats::setNames(wrap_label(label_of), as.character(seq_along(order_ids)))

  nodes <- data.frame(name = df$Description, count = df$Count,
                      cls = factor(cls, levels = seq_along(order_ids)),
                      label = ifelse(df$Description %in% label_of, df$Description, NA_character_),
                      stringsAsFactors = FALSE)
  idx <- which(m >= sim_cut & upper.tri(m), arr.ind = TRUE)
  edges <- data.frame(from = nodes$name[idx[, 1]], to = nodes$name[idx[, 2]],
                      weight = m[idx], stringsAsFactors = FALSE)
  g <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

  # A force-directed layout leaves the terms that no other term is similar to
  # (isolated nodes) scattered across the canvas, and igraph packs disconnected
  # components into their own grid cells. Laying the graph out once and then
  # tucking each isolated term next to the rest of its own class keeps the figure
  # to one tidy blob.
  lay <- igraph::layout_with_fr(tidygraph::as.igraph(g), niter = 2000)
  deg <- igraph::degree(tidygraph::as.igraph(g))
  for (v in which(deg == 0L)) {
    mates <- which(deg > 0L & cls == cls[v])
    base <- if (length(mates)) colMeans(lay[mates, , drop = FALSE]) else c(0, 0)
    ang <- 2 * pi * (v %% 12L) / 12
    lay[v, ] <- base + c(cos(ang), sin(ang)) * 0.45
  }

  # The class names live in the legend only. Drawing the representative term on
  # the canvas as well - which is what the first version did - put a bold label
  # on top of the class it names and duplicated the legend text.
  # The class legend is the only place where the class names appear, so its text
  # and keys are deliberately large: 50% larger than the previous 13/14 pt and
  # 18 pt key. The key glyphs are overridden to 1.5x their natural size so the
  # legend dots stay comparable with the largest nodes in the panel.
  class_breaks <- size_breaks(nodes$count)
  class_key_sizes <- seq(3.3, 10.5, length.out = max(1L, length(class_breaks)))
  cluster_plot <- ggraph::ggraph(g, layout = lay) +
    ggraph::geom_edge_link0(ggplot2::aes(edge_width = weight), colour = "grey55", alpha = 0.7) +
    ggraph::scale_edge_width(range = c(0.15, 1.0), guide = "none") +
    ggraph::geom_node_point(ggplot2::aes(size = count, fill = cls),
                            shape = 21, colour = "skyblue", stroke = 0.9, alpha = 0.95) +
    ggplot2::scale_size_continuous(range = c(2.2, 9), breaks = class_breaks,
                                   name = "Genes in term",
                                   guide = ggplot2::guide_legend(
                                     override.aes = list(size = class_key_sizes, shape = 21,
                                                         fill = "grey70", colour = "skyblue",
                                                         stroke = 0.9))) +
    ggplot2::scale_fill_manual(values = pal, name = "Functional class", labels = legend_lab) +
    # The ten class symbols in the colour legend are drawn two steps larger than
    # the default key glyph, so the colour of every class is easy to match.
    ggplot2::guides(fill = ggplot2::guide_legend(
      override.aes = list(size = 6.5, shape = 21, colour = "skyblue",
                          stroke = 0.9, alpha = 0.95))) +
    ggraph::theme_graph(base_family = font) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 19.5),
                   legend.title = ggplot2::element_text(size = 21, face = "bold"),
                   legend.key.height = ggplot2::unit(27, "pt"),
                   legend.key.width = ggplot2::unit(27, "pt"),
                   legend.spacing.y = ggplot2::unit(4, "pt"),
                   plot.title = ggplot2::element_text(face = "bold", size = 15)) +
    ggplot2::labs(title = paste0("GO functional classes of the enriched terms",
                                 if (!is.null(contrast)) paste0(": ", contrast) else ""))

  # --- per-class term-gene networks -------------------------------------------
  # One figure per functional class. Within a class the term holding the most
  # genes is chosen as the hub, and its genes are ranked by |log2FC| so the 50
  # strongest effects (up and down) are shown. The earlier single panel gave
  # every class an equal gene quota; with a 50-gene budget per class the largest
  # term is the natural representative and each class gets its own figure.
  gene_plots <- list()
  gene_terms <- list()
  for (g in seq_along(order_ids)) {
    rows <- df[cls == g, , drop = FALSE]
    if (!nrow(rows)) next
    # "Genes in term" is the Count column: the term with the largest gene set.
    rows <- rows[order(-rows$Count, rows$p.adjust), , drop = FALSE]
    term_row <- rows[1L, , drop = FALSE]
    term_label <- as.character(term_row$Description[1])
    term_genes <- unique(strsplit(as.character(term_row$geneID[1]), "/")[[1]])
    term_genes <- term_genes[nzchar(term_genes)]
    if (!length(term_genes)) next
    fc_vec <- if (!is.null(fc_symbols)) unname(fc_symbols[toupper(term_genes)]) else
      rep(NA_real_, length(term_genes))
    names(fc_vec) <- term_genes
    # Rank by magnitude so both strongly induced and strongly repressed genes
    # make the cut; genes without a fold change go last.
    ord <- order(-abs(ifelse(is.na(fc_vec), 0, fc_vec)))
    sel <- utils::head(term_genes[ord], as.integer(max_genes))
    fc_sel <- as.numeric(fc_vec[sel])
    # Clockwise from the top: repressed genes first, induced genes last, so the
    # two directions occupy opposite arcs of the star.
    ord2 <- order(ifelse(is.na(fc_sel), 0, fc_sel))
    sel <- sel[ord2]
    fc_sel <- fc_sel[ord2]
    n_g <- length(sel)
    ang <- seq(pi / 2, pi / 2 + 2 * pi, length.out = n_g + 1L)[seq_len(n_g)]
    lay2 <- rbind(c(0, 0), cbind(cos(ang), sin(ang)))
    nodes2 <- data.frame(
      name = c(paste0("T|", term_label), paste0("G|", sel)),
      fc = c(NA_real_, fc_sel),
      is_term = c(TRUE, rep(FALSE, n_g)),
      label = c(NA_character_, sel),
      stringsAsFactors = FALSE)
    edges2 <- data.frame(from = paste0("T|", term_label), to = paste0("G|", sel),
                         stringsAsFactors = FALSE)
    g2 <- tidygraph::tbl_graph(nodes = nodes2, edges = edges2, directed = FALSE)
    fc_limit <- display_lfc_limit(cfg, fc_sel)
    p <- ggraph::ggraph(g2, layout = lay2) +
      ggraph::geom_edge_link(colour = "grey70", alpha = 0.45, width = 0.4,
                             end_cap = ggraph::circle(3.0, "mm")) +
      # The term is a square painted in its class colour (na.value of the fold
      # change ramp); genes are circles on the green-white-purple ramp.
      ggraph::geom_node_point(ggplot2::aes(size = is_term, colour = fc, shape = is_term)) +
      ggplot2::scale_shape_manual(values = c(`TRUE` = 15, `FALSE` = 19), guide = "none") +
      ggplot2::scale_size_manual(values = c(`TRUE` = 8, `FALSE` = 3.2), guide = "none") +
      ggplot2::scale_colour_gradient2(low = "green4", mid = "grey88", high = "purple3",
                                      midpoint = 0, limits = c(-fc_limit, fc_limit), oob = scales::squish,
                                      na.value = pal[[g]],
                                      name = expression("log"[2] * " fold change")) +
      ggraph::geom_node_text(ggplot2::aes(label = label), repel = TRUE,
                             size = 4.3, colour = "black", max.overlaps = Inf, na.rm = TRUE,
                             point.padding = ggplot2::unit(0.12, "lines")) +
      # The hub's function name is carried by the title; a label inside the
      # panel would have 50 spokes running through it, and a filled box would
      # hide part of the network.
      ggraph::theme_graph(base_family = font) +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 12),
                     legend.title = ggplot2::element_text(size = 13, face = "bold"),
                     plot.title = ggplot2::element_text(face = "bold", size = 14),
                     plot.subtitle = ggplot2::element_text(size = 10.5, colour = "grey25"),
                     plot.caption = ggplot2::element_text(size = 9, colour = "grey40")) +
      ggplot2::labs(
        title = paste0("GO term-gene network - class ", g, ": ", term_label),
        subtitle = paste0("Term with the most genes in this class (", term_row$Count[1],
                          " genes); showing the top ", n_g, " of ", length(term_genes),
                          " by |log2FC|"),
        caption = "Square = GO term (class colour); circle = gene, coloured by log2FC.",
        x = NULL, y = NULL)
    gene_plots[[g]] <- p
    gene_terms[[g]] <- data.frame(
      class = g, class_leader = label_of[g], term = term_label,
      term_count = as.integer(term_row$Count[1]),
      genes_in_term = length(term_genes), genes_shown = n_g,
      top_gene = sel[which.max(abs(ifelse(is.na(fc_sel), 0, fc_sel)))],
      stringsAsFactors = FALSE)
  }

  # One row per term, carrying the class it was assigned to, the class
  # representative, and the gene symbols of the term (readable = TRUE already
  # stores symbols in geneID). Rows are grouped by class so the table can be read
  # class by class, matching the figure.
  table <- data.frame(
    Description = df$Description,
    class = cls,
    class_representative = label_of[cls],
    genes = gsub("/", "; ", as.character(df$geneID)),
    Count = df$Count,
    p.adjust = df$p.adjust,
    stringsAsFactors = FALSE)
  table <- table[order(table$class, table$p.adjust), , drop = FALSE]
  rownames(table) <- NULL
  list(cluster_plot = cluster_plot, table = table,
       gene_plots = gene_plots,
       gene_terms = if (length(gene_terms)) do.call(rbind, gene_terms) else NULL,
       class_labels = label_of)
}

# ---------------------------------------------------------------------------
# Legend size breaks
# ---------------------------------------------------------------------------
# A size legend is only honest when every break lies inside the range of the
# plotted values. pretty() routinely returns an endpoint just past the maximum
# (data max 96 -> break at 100), and that key is then drawn larger than any dot
# in the panel - the mismatch a reader notices first. Clipping the breaks to the
# data range keeps each key dot directly comparable with the panel.
size_breaks <- function(x, n = 4L, min_breaks = 3L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(numeric(0))
  rng <- range(x)
  if (rng[1] == rng[2]) return(rng[1])
  br <- pretty(rng, n = n)
  # The data minimum and maximum are always keys: the smallest key then matches
  # the smallest dot and the largest key the largest dot, which is what makes a
  # size legend readable against the panel.
  br <- sort(unique(c(rng[1], br[br >= rng[1] & br <= rng[2]], rng[2])))
  if (length(br) < min_breaks) {
    br <- unique(c(rng[1], rng[1] + diff(rng) / 2, rng[2]))
  }
  br
}

# ---------------------------------------------------------------------------
# Dependency check
# ---------------------------------------------------------------------------
# Every package the workflow can load. `required` is the differential-expression
# core, which has no fallback: if one of those is missing after the install
# attempt, the run stops with the list. Everything else only disables the
# section or figure that uses it, and those sections already write a "skipped"
# status when their package is absent. `dorothea` is deliberately absent: the
# regulons fall back to the bundled data/dorothea/*.rda copies, so the package
# is an optimisation rather than a requirement.
# Annotation package with exon coordinates for the species. Kept here as well
# as in context.R (gene_length_annotation) so that the dependency check does not
# depend on the order in which the two files are sourced.
gene_length_package <- function(preset) {
  switch(tolower(preset$common_name %||% "mouse"),
         mouse = "EnsDb.Mmusculus.v79",
         human = "EnsDb.Hsapiens.v86",
         rat = "TxDb.Rnorvegicus.UCSC.rn6.refGene",
         character())
}

pipeline_packages <- function(species = "mouse") {
  preset <- tryCatch(species_preset(species), error = function(e) species_preset("mouse"))
  list(
    required = c("ggplot2", "scales", "tidyr", "jsonlite", "yaml", "readxl",
                 "pheatmap", "DESeq2", "apeglm", "SummarizedExperiment",
                 "BiocParallel", "rlang",
                 # The fixed-template Word report (summary_report = TRUE) is
                 # part of a default run, so its writers are required too.
                 "officer", "flextable",
                 # Used by the dependency checker itself to install
                 # Bioconductor packages.
                 "BiocManager"),
    optional = c("writexl", "reshape2", "gplots", "systemfonts", "RColorBrewer",
                 "UpSetR", "igraph", "ggraph", "tidygraph", "clusterProfiler",
                 "enrichplot", "DOSE", "fgsea", "msigdbr", "GSVA", "decoupleR",
                 "pathview", "GO.db", "AnnotationDbi", preset$orgdb,
                 "rmarkdown", "knitr", "patchwork", "cluster",
                 # Drawing code that lives behind a package's Suggests list is
                 # not installed with that package: enrichplot::ridgeplot draws
                 # through ggridges, section 4 can use ggupset, section 1 labels
                 # through ggrepel. dorothea supplies the regulons; the bundled
                 # data/dorothea/*.rda files remain the offline fallback.
                 "ggridges", "ggrepel", "ggupset", "dorothea",
                 # Figure atlas (scripts/build_figure_atlas.R): image handling
                 # and PDF concatenation.
                 "magick", "pdftools",
                 # FPKM displays: exon-union gene lengths are read from a GTF or
                 # from the species annotation package when no expression matrix
                 # is supplied. Optional: a supplied matrix needs none of these.
                 "rtracklayer", "GenomicRanges", "S4Vectors", "GenomicFeatures",
                 gene_length_package(preset))
  )
}

# Packages served by Bioconductor. install.packages() cannot see most of them,
# so the two repositories are installed separately; retrying a Bioconductor
# package on CRAN only replaces the real build error with "package not
# available", which is what made earlier failed installs hard to diagnose.
bioc_pipeline_packages <- function() {
  c("DESeq2", "apeglm", "SummarizedExperiment", "BiocParallel",
    "clusterProfiler", "enrichplot", "DOSE", "fgsea", "GSVA", "decoupleR",
    "pathview", "GO.db", "AnnotationDbi", "dorothea",
    "org.Mm.eg.db", "org.Hs.eg.db", "org.Rn.eg.db",
    "rtracklayer", "GenomicRanges", "S4Vectors", "GenomicFeatures")
}

package_table <- function(pkgs, required) {
  data.frame(
    package = pkgs,
    required = pkgs %in% required,
    installed = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE),
    version = vapply(pkgs, function(p) {
      tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
    }, character(1)),
    stringsAsFactors = FALSE)
}

# A system-wide R library is often read-only (the ClawsGO host ships its stack
# in /usr/local/lib/R/site-library), and install.packages() then fails with
# "'lib = ...' is not writable" before it even reaches the repository. Prepend a
# writable user library so every later install and requireNamespace() sees it.
ensure_writable_library <- function(quiet = TRUE) {
  libs <- .libPaths()
  writable <- libs[file.access(libs, mode = 2L) == 0L]
  if (length(writable)) {
    target <- normalizePath(writable[1L], winslash = "/", mustWork = FALSE)
    # install.packages() and BiocManager::install() use .libPaths()[1], so a
    # writable library further down the list has to be moved to the front.
    if (!identical(target, normalizePath(libs[1L], winslash = "/", mustWork = FALSE))) {
      .libPaths(c(target, .libPaths()))
    }
    return(target)
  }
  user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
  user_lib <- strsplit(user_lib, .Platform$path.sep, fixed = TRUE)[[1L]]
  user_lib <- user_lib[nzchar(user_lib)]
  if (!length(user_lib)) {
    user_lib <- file.path(path.expand("~"), "R", paste0("library-", getRversion()))
  }
  target <- user_lib[1L]
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  if (file.access(target, mode = 2L) != 0L) return(NULL)
  .libPaths(c(target, .libPaths()))
  if (!quiet) message("Dependency check: installing into user library ", target)
  normalizePath(target, winslash = "/", mustWork = FALSE)
}

# Install one package and keep the conditions it raised, so a failed build can
# report what actually went wrong instead of a bare "still missing".
install_one_package <- function(pkg, bioc, repos) {
  notes <- character()
  tryCatch(
    withCallingHandlers(
      if (bioc) {
        BiocManager::install(pkg, update = FALSE, ask = FALSE, dependencies = TRUE)
      } else {
        utils::install.packages(pkg, repos = repos, dependencies = TRUE)
      },
      warning = function(w) {
        notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning")
      },
      message = function(m) {
        notes <<- c(notes, conditionMessage(m)); invokeRestart("muffleMessage")
      }
    ),
    error = function(e) notes <<- c(notes, conditionMessage(e))
  )
  list(installed = requireNamespace(pkg, quietly = TRUE), notes = notes)
}

# Release number of the installed Bioconductor stack. BiocVersion pins it
# locally; BiocManager::version() would otherwise fetch config.yaml from the
# network before any mirror has been chosen.
bioc_release_version <- function() {
  if (requireNamespace("BiocVersion", quietly = TRUE)) {
    v <- as.character(utils::packageVersion("BiocVersion"))
    parts <- strsplit(v, ".", fixed = TRUE)[[1L]]
    if (length(parts) >= 2L) return(paste(parts[1:2], collapse = "."))
  }
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = 10)
  tryCatch(as.character(BiocManager::version()), error = function(e) NA_character_)
}

# Pick the first Bioconductor mirror whose package index answers quickly.
# bioconductor.org can take tens of seconds per request from some networks, so
# an install started against it looks hung even when the library is writable;
# the mirrors below carry the same 3.x release.
choose_bioc_mirror <- function(quiet = TRUE) {
  candidates <- unique(c(getOption("BioC_mirror", "https://bioconductor.org"),
                         "https://mirrors.tuna.tsinghua.edu.cn/bioconductor",
                         "https://mirrors.westlake.edu.cn/bioconductor"))
  ver <- bioc_release_version()
  if (is.na(ver) || !nzchar(ver)) return(candidates[1L])
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  for (mirror in candidates) {
    if (!nzchar(mirror)) next
    url <- paste0(sub("/+$", "", mirror), "/packages/", ver, "/bioc/src/contrib/PACKAGES")
    options(timeout = 10)
    con <- NULL
    ok <- tryCatch({
      con <- url(url, open = "rb")
      # Read a slice large enough to measure throughput, not just the time to
      # the first byte: bioconductor.org can answer a 1 KB probe instantly and
      # then stall on the full index.
      length(readBin(con, "raw", n = 200000L)) > 0L
    }, error = function(e) FALSE, warning = function(w) FALSE,
    finally = if (!is.null(con)) try(close(con), silent = TRUE))
    if (isTRUE(ok)) {
      if (!quiet && !identical(mirror, candidates[1L])) {
        message("Dependency check: using Bioconductor mirror ", mirror)
      }
      return(mirror)
    }
  }
  candidates[1L]
}

# Install a set of packages, routing Bioconductor names through BiocManager and
# the rest through install.packages(). A batch install is tried first because
# dependency resolution is cheaper that way; whatever the batch leaves behind is
# retried one package at a time, so one broken build does not take the rest of
# the list down with it.
install_r_packages <- function(missing, quiet = FALSE) {
  missing <- unique(missing)
  if (!length(missing)) {
    return(list(installed = character(), failed = character(), notes = character()))
  }
  install_lib <- ensure_writable_library(quiet = quiet)
  if (is.null(install_lib)) {
    return(list(installed = character(), failed = missing,
                notes = "No writable R library: every entry in .libPaths() is read-only and the user library could not be created."))
  }
  if (!quiet) message("Dependency check: installing into ", install_lib)
  cran <- c(CRAN = "https://cloud.r-project.org")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    if (!quiet) message("Dependency check: installing BiocManager first ...")
    try(utils::install.packages("BiocManager", repos = cran), silent = TRUE)
  }
  have_bioc <- requireNamespace("BiocManager", quietly = TRUE)
  old_mirror <- getOption("BioC_mirror")
  if (have_bioc) {
    chosen_mirror <- choose_bioc_mirror(quiet = quiet)
    if (!is.null(chosen_mirror) && nzchar(chosen_mirror)) {
      options(BioC_mirror = chosen_mirror)
    }
    on.exit(options(BioC_mirror = old_mirror), add = TRUE)
  }
  repos <- if (have_bioc) {
    tryCatch(BiocManager::repositories(), error = function(e) cran)
  } else {
    cran
  }
  if (is.null(repos[["CRAN"]]) || is.na(repos[["CRAN"]])) {
    repos[["CRAN"]] <- unname(cran)
  }
  bioc_set <- bioc_pipeline_packages()
  bioc_missing <- intersect(missing, bioc_set)
  cran_missing <- setdiff(missing, bioc_set)
  notes <- character()

  batch <- function(pkgs, bioc) {
    if (!length(pkgs)) return(invisible(NULL))
    if (!quiet) {
      message("Dependency check: installing ", if (bioc) "Bioconductor" else "CRAN",
              " package(s): ", paste(pkgs, collapse = ", "))
    }
    tryCatch(
      withCallingHandlers(
        if (bioc) {
          if (!have_bioc) stop("BiocManager is not available, cannot install ", paste(pkgs, collapse = ", "))
          BiocManager::install(pkgs, update = FALSE, ask = FALSE, dependencies = TRUE)
        } else {
          utils::install.packages(pkgs, repos = repos, dependencies = TRUE)
        },
        warning = function(w) {
          notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning")
        },
        message = function(m) {
          notes <<- c(notes, conditionMessage(m)); invokeRestart("muffleMessage")
        }
      ),
      error = function(e) notes <<- c(notes, paste0("batch install: ", conditionMessage(e)))
    )
    invisible(NULL)
  }
  batch(bioc_missing, TRUE)
  batch(cran_missing, FALSE)

  installed_now <- function(pkgs) {
    pkgs[vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  }
  # bioconductor.org is slow from some networks (the index alone can exceed the
  # default 60 s), so a failed Bioconductor batch is retried through the
  # mirrors that keep the same 3.x release. CRAN packages get the same treatment
  # when the configured mirror does not carry them.
  still_bioc <- setdiff(bioc_missing, installed_now(bioc_missing))
  if (length(still_bioc)) {
    old_timeout <- getOption("timeout")
    options(timeout = max(300, old_timeout))
    on.exit(options(timeout = old_timeout), add = TRUE)
    for (mirror in c("https://mirrors.tuna.tsinghua.edu.cn/bioconductor",
                     "https://mirrors.westlake.edu.cn/bioconductor")) {
      if (!length(still_bioc)) break
      if (!quiet) message("Dependency check: retrying Bioconductor packages via ", mirror)
      options(BioC_mirror = mirror)
      batch(still_bioc, TRUE)
      still_bioc <- setdiff(still_bioc, installed_now(still_bioc))
    }
  }
  still_cran <- setdiff(cran_missing, installed_now(cran_missing))
  if (length(still_cran)) {
    if (!quiet) message("Dependency check: retrying CRAN packages on https://cloud.r-project.org")
    repos[["CRAN"]] <- "https://cloud.r-project.org"
    batch(still_cran, FALSE)
  }

  still <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  failed <- character()
  for (pkg in still) {
    res <- install_one_package(pkg, pkg %in% bioc_set, repos)
    if (!res$installed) {
      failed <- c(failed, pkg)
      keep <- utils::tail(res$notes[nzchar(res$notes)], 6L)
      if (length(keep)) notes <- c(notes, paste0(pkg, ": ", keep))
    }
  }
  list(installed = setdiff(missing, failed), failed = unique(failed),
       notes = unique(notes))
}

# Check every package, install what is missing, check again, and only then let
# the analysis start. Returns the final status table. A missing required package
# stops the run; a missing optional one is a warning, because the affected
# section already knows how to skip itself.
check_and_install_packages <- function(cfg, install = NULL, quiet = FALSE) {
  install <- if (is.null(install)) isTRUE(cfg$runtime$install_missing_packages) else isTRUE(install)
  plan <- pipeline_packages(cfg$species$common_name %||% "mouse")
  required <- unique(plan$required)
  pkgs <- unique(c(required, plan$optional))
  st <- package_table(pkgs, required)
  missing <- st$package[!st$installed]
  if (!quiet) {
    message("Dependency check: ", sum(st$installed), "/", nrow(st), " packages already installed.")
  }
  if (length(missing)) {
    if (!install) {
      message("Dependency check: ", length(missing), " missing package(s), installation disabled: ",
              paste(missing, collapse = ", "))
    } else {
      message("Dependency check: installing ", length(missing), " missing package(s) ...")
      res <- install_r_packages(missing, quiet = quiet)
      if (length(res$failed) && !quiet) {
        message("Dependency check: ", length(res$failed),
                " package(s) could not be installed: ", paste(res$failed, collapse = ", "))
        if (length(res$notes)) {
          message(paste0("  ", res$notes, collapse = "\n"))
        }
        message("  Bioconductor packages must come from BiocManager; if a build failed, ",
                "check the compiler/system libraries and run\n",
                "  BiocManager::install(c(",
                paste0('"', res$failed, '"', collapse = ", "), "), dependencies = TRUE)")
      }
      st <- package_table(pkgs, required)
    }
  }
  still_missing <- st$package[!st$installed]
  missing_required <- st$package[st$required & !st$installed]
  missing_optional <- st$package[!st$required & !st$installed]
  if (!quiet) {
    message("Dependency check: ", sum(st$installed), "/", nrow(st), " packages available",
            if (length(still_missing)) paste0("; still missing: ", paste(still_missing, collapse = ", ")) else ".")
    writers <- c("officer", "flextable")
    have_writers <- writers[vapply(writers, requireNamespace, logical(1), quietly = TRUE)]
    message("Word report writers: ",
            if (length(have_writers) == length(writers)) "officer and flextable available"
            else paste0("UNAVAILABLE (", paste(setdiff(writers, have_writers), collapse = ", "),
                        "); results_report.docx cannot be generated"))
  }
  # Record the versions actually used next to the run manifest, so a result can
  # be traced back to a package set without rerunning sessionInfo().
  out_dir <- cfg$project$output_dir
  if (!is.null(out_dir) && nzchar(out_dir %||% "")) {
    prov <- safe_dir(file.path(out_dir, "provenance"))
    try(utils::write.csv(st, file.path(prov, "packages.csv"), row.names = FALSE), silent = TRUE)
  }
  if (length(missing_required)) {
    report_writers <- intersect(missing_required, c("officer", "flextable"))
    mac_hint <- if (length(report_writers) &&
                    .Platform$OS.type == "unix" && Sys.info()[["sysname"]] == "Darwin") {
      paste0("\n  macOS: run `xcode-select --install` before installing ",
             paste(report_writers, collapse = " and "),
             "; both need a C++ toolchain (systemfonts/textshaping) and are ",
             "the packages that write results_report.docx.")
    } else ""
    stop("Required R package(s) missing: ", paste(missing_required, collapse = ", "),
         ".\n  Install them with BiocManager::install(c(",
         paste0('"', missing_required, '"', collapse = ", "), ")) and rerun.",
         mac_hint,
         if (length(report_writers))
           "\n  To run the analysis without the Word report, use summary_report = FALSE."
         else "",
         call. = FALSE)
  }
  if (length(missing_optional)) {
    warning("Optional package(s) missing, their sections/figures will be skipped: ",
            paste(missing_optional, collapse = ", "), call. = FALSE)
    if ("dorothea" %in% missing_optional) {
      message("  dorothea is unavailable; the bundled data/dorothea/*.rda regulons will be used instead.")
    }
  }
  invisible(st)
}
