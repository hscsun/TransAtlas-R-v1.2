# run_all() - the workflow's single entry point.
#
# In R / RStudio, from this project's root folder:
#
#   setwd("/Users/changye/Documents/TransAtlas-R-v1.2")
#   source("pipeline/run_all.R")
#   run_all(species = "mouse",log2fc = 0.75)                                   # pick the count matrix in a file dialog
#
# or give more things on one line:
#
#   run_all(counts_file = "input/NIH3T3.xlsx",
#           samples = "input/NIH3T3_metadata.csv", species = "mouse",
#           contrast = c("FBS", "Control"), log2fc = 0.75)
#
# The output goes to results/<name of the data file>/. There is no
# configuration file: every setting is an argument here, and the built-in
# defaults (default_config() in pipeline/R/helpers.R) apply to anything omitted.
#
# The two things a user has to provide are the count matrix and the sample
# metadata. TF-target regulons are read from the DoRothEA data of the
# `dorothea` package (hs for human, mm for mouse and rat).

# Locate the project root from this file rather than from the working
# directory, so `source("path/to/run_all.R")` works from anywhere.
run_all_project_root <- function() {
  candidates <- character()
  for (i in rev(seq_len(sys.nframe()))) {
    ofile <- sys.frame(i)$ofile
    if (is.character(ofile) && length(ofile) == 1L && nzchar(ofile)) {
      candidates <- c(candidates, file.path(dirname(ofile), ".."))
    }
  }
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) candidates <- c(candidates, file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."))
  candidates <- c(candidates, getwd())
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "pipeline", "R", "helpers.R"))) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Cannot locate the project root. Source pipeline/run_all.R by its full path.", call. = FALSE)
}

.clawsgo_root <- run_all_project_root()
for (.f in c("figure_sizes.R", "helpers.R", "context.R", "kegg_module.R", "deseq2_module.R",
             "provenance.R", "qc_module.R", "run_pipeline.R")) {
  source(file.path(.clawsgo_root, "pipeline", "R", .f))
}
rm(.f)

# Two input files are all a run needs: the count matrix and the sample table.
# The TF-target regulons come from the DoRothEA data in the `dorothea` package
# (human = dorothea_hs, mouse and rat = dorothea_mm), so nothing else has to be
# prepared. There is no configuration file: every setting below has a built-in
# default (default_config() in pipeline/R/helpers.R) and can be changed here.
#
# Input and design
# counts_file  Path to the raw-count matrix (.xlsx/.xls/.csv/.txt). NULL opens a
#              file-selection dialog (interactive sessions only).
# samples      Sample grouping table: a path, or a data.frame with sample_id and
#              condition (plus optional replicate). NULL looks for
#              <data>_metadata.csv or sample_metadata.csv next to the data file,
#              then in config/.
# species      "mouse", "human" or "rat" - picks the OrgDb, KEGG code, identifier
#              types, MSigDB database and DoRothEA collection (hs / mm) in one
#              word. Rat uses the mouse regulons.
# group_column Column of the metadata that holds the two groups. Default
#              "condition".
# contrast     c(numerator, denominator), e.g. c("Treatment", "Control"). NULL
#              derives it from the metadata levels: the second level
#              alphabetically is compared against the first.
#
# Differential expression
# log2fc       The one |log2 fold change| cut-off of the workflow. It calls a
#              gene differentially expressed and defines DEG gene sets.
#              GSEA has its own low_cutoff. Default 0.75.
# alpha        Adjusted-P cut-off. Default 0.05.
# min_fpkm     Expression floor shared by every section: the higher of a gene's
#              two group-mean FPKM values must reach this value, otherwise the
#              gene is dropped before screening, enrichment and the regulatory
#              sections. NULL switches the floor off. Default 5.
# min_count    Minimum raw count a gene needs in at least min_samples samples to
#              be tested at all (DESeq2 pre-filter). Default 10.
# min_samples  Number of samples that must reach min_count. Default 2.
#
# Input columns
# gene_id_column, gene_symbol_column, annotation_column
#              Names of the annotation columns in the count workbook. The
#              workbook: when no expression matrix is supplied, FPKM is computed
#              from the exon-union lengths in gene_length_column, and only when
#              that column is absent from a GENCODE/Ensembl GTF (the
#              gene_annotation argument, options(transatlas.gene_annotation=),
#              a *.gtf next to the counts, or the TRANSATLAS_GENE_ANNOTATION
#              variable) or from the species annotation package.
# gene_length_column
#              Workbook column with the exon-union length of each gene.
#              Default "Exonic.gene.sizes"; NULL or "none" disables the column.
# gene_annotation
#              Path to the GTF/GFF used for the exon-union lengths when the
#              workbook carries no length column. NULL uses the search order
#              above.
# counts_sheet Sheet of an Excel workbook. NULL = first sheet.
# count_columns Sample columns of the count matrix. NULL = inferred from the
#              metadata.
# expression_matrix_file
#              Optional expression matrix (FPKM/TPM) used by the display
#              sections. NULL looks for <counts name>_fpkm.csv next to the
#              count workbook, then computes FPKM from the raw counts and the
#              annotation exon lengths, and only then falls back to the DESeq2
#              normalised counts for displays when min_fpkm = NULL. An active
#              FPKM floor requires a usable FPKM matrix.
# expression_matrix_unit
#              Unit of a supplied matrix, default "FPKM". TPM and normalized
#              counts require min_fpkm = NULL.
#
# Transcription factors
# top_tfs      Number of transcription factors analysed and drawn (dot plot and
#              network). Default 30.
# confidence   DoRothEA confidence tiers to keep. Default c("A", "B").
# tf_min_targets, tf_ora_min_overlap, tf_coherence_min_targets
#              Regulon size floor, minimum DEGs inside a regulon for ORA, and
#              minimum targets for the coherence score.
# tf_file      Optional custom TF-target table (TF, target, optional
#              confidence/mor). NULL uses the DoRothEA data for the species.
# tf_dorothea_dir
#              Folder holding the bundled DoRothEA .rda files. Default
#              data/dorothea.
# cluster_min_jaccard, cluster_max
#              Section 8 clusters the screened TFs by the overlap of their DE
#              target sets. TFs sharing at least this Jaccard index are linked,
#              and at most this many modules are cut. Default 0.1 and 8.
#
# Display
# max_genes    Upper bound on genes drawn in any heat map or dot plot. Default 30.
# top_n_labels Maximum labels on the volcano / expression-effect panels. Default 20.
# volcano_label_expression_n, volcano_label_effect_n
#              Candidate labels per direction by expression and effect, default 5 each.
# up_cutoff    Display-only symmetric log2FC limit, default 2.5 (-2.5 to 2.5).
#              NULL uses an automatic display range. Calculations and CSVs stay raw.
# low_cutoff   GSEA-only minimum absolute log2FC, default 0.15. NULL/0 disables
#              this lower bound; gene padj < alpha is still required.
# log2fc_cap   Deprecated display-only alias of up_cutoff.
# volcano_lfc_thresholds
#              Screening ladder, default 0.15-2.5, always including log2fc.
#              Independent of the GSEA low_cutoff parameter.
# expression_log10_threshold
#              Legacy AnaDGEsS1 rule that hid genes with log10(baseMean) <= 2
#              from the up/down classes. Off (NULL) by default, so `min_fpkm`
#              is the only expression threshold; pass 2 to restore it.
#
# Runtime
# check_packages, install_missing_packages
#              On start, check every package the run can use; install what is
#              missing with BiocManager/CRAN, check again, then analyse. A
#              missing required package stops the run with the exact list.
#              Default both TRUE; set check_packages = FALSE to skip the check.
#
# Enrichment
# go_pvalue_cutoff, go_qvalue_cutoff, go_bp_pvalue_cutoff, go_bp_qvalue_cutoff
#              GO over-representation cut-offs.
# go_simplify_cutoff  Semantic-similarity cutoff for GO simplification.
# gsea_min_size, gsea_max_size, gsea_pvalue_cutoff, gsea_top_sets
#              GSEA size window, reporting cutoff and number of sets drawn on
#              the composite panels.
# run_msigdb, msigdb_category  MSigDB collection GSEA.
# pathview_download, pathview_max_pathways
#              Draw pathview maps for the leading KEGG pathways. pathview fetches
#              each pathway's KGML and native map from rest.kegg.jp, so
#              pathview_download = FALSE skips those downloads together with the
#              maps; the section's enrichment, lollipop, GSEA panel and UpSet
#              plot still run. Default TRUE (download and draw) and 12 maps.
# run_pathview Deprecated alias of pathview_download, kept for older scripts;
#              when supplied it takes precedence.
#
# Regulatory and advanced
# network_layouts     Layouts for the TF-target network. Default all four.
# network_max_labels  Number of target labels in the network. Default 30.
# activity_methods    decoupleR methods for regulon activity.
# activity_min_size, activity_alpha  Regulon size floor and significance level.
#              NULL activity_alpha inherits alpha.
# run_gsva, gsva_collection, gsva_min_size, gsva_max_size  GSVA pathway activity.
#
# Protein-protein interactions
# run_ppi      TRUE queries STRING v12 for the DEG-encoded proteins and draws
#              the interaction network and its hub proteins.
# ppi_score_threshold  Minimum STRING combined score. Default 0.4.
# ppi_max_nodes  Maximum number of proteins submitted, ranked by significance
#              and effect size. Default 300.
# ppi_label_n  Maximum hub proteins labelled per module. Default 4; large modules
#              (>= 20 proteins) always receive four labels, smaller modules fewer.
# ppi_max_modules  Maximum number of Louvain modules; smaller communities are
#              merged into their best-connected neighbour beyond this ceiling.
#              Default 10. Each module is annotated with its top GO BP terms.
# ppi_cache_dir, ppi_timeout  Where STRING responses are cached and how long a
#              request may take. Default cache/string and 120 s.
#
# Quality control
# qc           TRUE draws the configured QC figures; FALSE skips them.
# qc_plots     Optional subset: any of "pca", "correlation", "library_size",
#              "count_distribution", "ma_plot", "dispersion",
#              "independent_filter", "sample_distance", "rlog_vst".
#
# Output and runtime
# sections     Sections to run, e.g. c(1, 2, 3). NULL runs the configured set.
# report       TRUE renders the HTML report at the end of the run.
# report_title, report_author, report_language, generate_pdf  Report metadata.
# summary_report  TRUE also writes the fixed-template Word report
#              (report/results_report.docx) filled with this run's results.
# kegg_cache_dir, kegg_download_timeout  Where the KEGG tables are cached and
#              how long the first download may take.
# continue_on_error  TRUE keeps going when an optional module fails.
# output_dir   Output folder. NULL uses results/<name of the data file>/.
# overwrite    TRUE reuses the output folder; FALSE appends _run2, _run3, ...
# seed         Random seed for the reproducible steps.
# config_file  Optional YAML overlay merged over the built-in defaults. Not
#              needed for a normal run.
run_all <- function(counts_file = NULL,
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
                    tf_min_fpkm = NULL,   # deprecated: use min_fpkm
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
                    gene_set_lfc_threshold = NULL,   # deprecated: use log2fc
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
                    project_root = NULL,
                    config_file = NULL) {
  project_root <- project_root %||% .clawsgo_root
  # Figure sizes are read from pipeline/R/figure_sizes.R on every call, not once
  # when run_all.R was sourced: edit a width or height there and the next
  # run_all() uses it without re-sourcing or restarting the session.
  reload_figure_sizes(project_root)
  reset_figure_size_log()
  cfg <- build_run_config(
    counts_file = counts_file, samples = samples, species = species,
    group_column = group_column, contrast = contrast, contrast_name = contrast_name,
    log2fc = log2fc, alpha = alpha, min_count = min_count, min_samples = min_samples,
    min_fpkm = min_fpkm,
    gene_id_column = gene_id_column, gene_symbol_column = gene_symbol_column,
    annotation_column = annotation_column,
    counts_sheet = counts_sheet, count_columns = count_columns,
    expression_matrix_file = expression_matrix_file, expression_matrix_unit = expression_matrix_unit,
    gene_annotation = gene_annotation,
    gene_length_column = gene_length_column,
    top_tfs = top_tfs, confidence = confidence, tf_min_fpkm = tf_min_fpkm,
    tf_min_targets = tf_min_targets, tf_ora_min_overlap = tf_ora_min_overlap,
    tf_coherence_min_targets = tf_coherence_min_targets,
    tf_file = tf_file, tf_dorothea_dir = tf_dorothea_dir,
    max_genes = max_genes, top_n_labels = top_n_labels,
    up_cutoff = up_cutoff, low_cutoff = low_cutoff, log2fc_cap = log2fc_cap,
    volcano_label_expression_n = volcano_label_expression_n,
    volcano_label_effect_n = volcano_label_effect_n,
    volcano_lfc_thresholds = volcano_lfc_thresholds,
    expression_log10_threshold = expression_log10_threshold,
    gene_set_lfc_threshold = gene_set_lfc_threshold,
    go_pvalue_cutoff = go_pvalue_cutoff, go_qvalue_cutoff = go_qvalue_cutoff,
    go_bp_pvalue_cutoff = go_bp_pvalue_cutoff, go_bp_qvalue_cutoff = go_bp_qvalue_cutoff,
    go_simplify_cutoff = go_simplify_cutoff,
    gsea_min_size = gsea_min_size, gsea_max_size = gsea_max_size,
    gsea_pvalue_cutoff = gsea_pvalue_cutoff, gsea_top_sets = gsea_top_sets,
    run_msigdb = run_msigdb, msigdb_category = msigdb_category,
    pathview_download = pathview_download, run_pathview = run_pathview,
    pathview_max_pathways = pathview_max_pathways,
    network_layouts = network_layouts, network_max_labels = network_max_labels,
    cluster_tfs = top_tfs, cluster_min_jaccard = cluster_min_jaccard,
    cluster_max = cluster_max,
    activity_methods = activity_methods, activity_min_size = activity_min_size,
    activity_alpha = activity_alpha,
    run_gsva = run_gsva, gsva_collection = gsva_collection,
    gsva_min_size = gsva_min_size, gsva_max_size = gsva_max_size,
    run_ppi = run_ppi, ppi_score_threshold = ppi_score_threshold,
    ppi_max_nodes = ppi_max_nodes, ppi_label_n = ppi_label_n,
    ppi_max_modules = ppi_max_modules, ppi_module_annotation_alpha = ppi_module_annotation_alpha,
    ppi_cache_dir = ppi_cache_dir, ppi_timeout = ppi_timeout,
    qc = qc, qc_plots = qc_plots,
    sections = sections, report = report,
    report_title = report_title, report_author = report_author,
    report_language = report_language, generate_pdf = generate_pdf,
    summary_report = summary_report,
    kegg_cache_dir = kegg_cache_dir, kegg_download_timeout = kegg_download_timeout,
    check_packages = check_packages, install_missing_packages = install_missing_packages,
    continue_on_error = continue_on_error,
    output_dir = output_dir, overwrite = overwrite, seed = seed,
    project_root = project_root, config_file = config_file)
  message("  counts   : ", cfg$input$counts_file)
  message("  samples  : ", if (is.data.frame(cfg$input$sample_metadata_file)) "in-memory data.frame" else cfg$input$sample_metadata_file)
  message("  species  : ", cfg$species$name, " (", cfg$species$orgdb, ", ",
          cfg$species$kegg_code, ", DoRothEA ", cfg$input$tf_species, ")")
  message("  expression: max group-mean FPKM >= ", cfg$analysis$min_fpkm %||% "off")
  message("  contrast : ", cfg$analysis$contrasts[[1L]]$name)
  message("  thresholds: |log2FC| >= ", cfg$analysis$lfc_threshold,
          ", adjusted P < ", cfg$analysis$alpha,
          ", regulon confidence ", paste(cfg$regulatory$confidence_keep, collapse = "/"))
  message("  GSEA     : padj < ", cfg$analysis$alpha, "; |log2FC| >= ",
          cfg$analysis$low_cutoff %||% 0, "; unbounded effects")
  message("  display  : log2FC limit +/-", cfg$analysis$up_cutoff %||% "auto")
  message("  TFs      : top ", cfg$regulatory$tf_dotplot_max, " of the expressed regulons")
  message("  sections : ", paste(cfg$sections$run, collapse = ", "))
  message("  fig sizes: ", length(FIGURE_SIZES),
          " keys from pipeline/R/figure_sizes.R (reloaded for this run)")
  message("  output   : ", cfg$project$output_dir)
  invisible(run_pipeline(cfg, project_root))
}
