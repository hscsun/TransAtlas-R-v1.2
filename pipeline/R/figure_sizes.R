# Central registry of every figure size written by the pipeline.
#
# All widths and heights are in inches; the dpi stays at the save_plot() default
# of 320. Change a number here, then redraw the affected figures without
# recomputing anything:
#
#   Rscript scripts/redraw_figures.R            # every cached ggplot figure
#   Rscript scripts/redraw_figures.R s3_go_composite
#
# Figures drawn with pheatmap or base graphics cannot be replayed from a cached
# object; change their size here and rerun just the section that owns them:
#
#   Rscript pipeline/run_analysis.R --sections 1 --no-report
#
# Entries may be a c(width, height) vector, or a function of the arguments the
# section passes through save_plot() (used by the panels whose height follows
# the number of rows).

# Standalone sourcing (scripts/redraw_figures.R) does not load helpers.R first.
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
}

FIGURE_SIZES <- list(
  # --- quality control (00_differential_expression) ------------------------
  qc_library_size          = c(6.5, 4.5),
  qc_count_distributions   = c(6.5, 4.5),
  qc_sample_correlation    = c(6, 5),
  qc_sample_distance       = c(6, 5),
  qc_pca                   = c(6.825, 5),
  qc_dispersion            = c(6, 5),
  qc_ma_plot               = c(6, 5),

  # --- section 1: screening ------------------------------------------------
  s1_threshold_counts       = c(6.5, 4.2),
  s1_volcano                = c(6.5, 6.4),
  s1_expression_effect      = c(6.5, 6.0),
  s1_volcano_pdf            = c(6.5, 9),
  s1_top_deg_heatmap        = c(7.5, 10),
  s1_global_deg_heatmap     = c(6, 6),
  s1_chromosome_landscape   = c(14, 8),

  # --- section 2: GO over-representation -----------------------------------
  # The 40% gene-overlap rule can leave well under 20 terms on a panel, so the
  # height follows the number of rows actually drawn instead of leaving blank
  # space below a short list.
  s2_go_dotplot             = function(n_terms = 20) c(9, max(4.5, 0.32 * n_terms + 2)),
  s2_go_classes             = c(19.2, 9.5),
  s2_go_term_gene           = c(12.5, 12),
  s2_richfactor             = c(10, 9),

  # --- section 3: GSEA -----------------------------------------------------
  s3_go_allsets_pdf         = c(8, 5.5),
  # Multi-pathway GSEA panel; 2026-10-18 height raised by 30% (7.5 -> 9.75),
  # 2026-10-20 width raised to 12 so the +2 pt legend labels are not clipped,
  # then narrowed by 25% (12 -> 9) at the user's request; long set names are
  # wrapped onto a second line in the legend so they still fit.
  s3_go_composite           = c(9, 9.75),
  s3_go_enrichment_map      = c(12, 9),
  s3_go_ridgeplot           = c(10, 7.5),
  s3_msigdb_barplot         = c(12, 8),

  # --- section 4: KEGG ------------------------------------------------------
  s4_kegg_lollipop          = c(6.375, 6.5),
  s4_kegg_composite         = c(7.5, 7.5),
  s4_kegg_upset             = c(11.2, 6),
  s4_kegg_upset_full        = c(12.14, 6),

  # --- section 5: TF regulons ----------------------------------------------
  s5_tf_dotplot             = function(n) c(6.4, max(4, 0.272 * n + 1.6)),
  s5_tf_expression          = function(n) c(7.5, max(4.5, 0.32 * n + 2.2)),
  s5_regulon_volcano        = c(8, 6),
  s5_tf_chromosome          = c(14, 8),
  s5_tf_chromosome_summary  = function(n) c(12, max(4, 0.42 * n + 1.8)),

  # --- section 6: regulatory activity --------------------------------------
  s6_activity_volcano      = c(7.5, 5.5),
  s6_activity_heatmap      = c(5.5, 7.5),
  s6_gsva_barplot          = c(8, 6),

  # --- section 7: TF-target network ----------------------------------------
  s7_tf_network             = c(18, 9),

  # --- section 8: TF modules -----------------------------------------------
  s8_tf_modules            = c(12, 9),

  # --- section 9: gene-set heat maps (GO) and dot plots (TF) ---------------
  s9_geneset_heatmap        = c(8.2, 10),
  s9_geneset_dotplot        = c(3.6, 6),

  # --- section 10: protein-protein interaction -----------------------------
  s10_ppi_network           = c(18, 9),
  s10_ppi_hubs              = c(7, 8),
  s10_ppi_functions         = function(n) c(12, max(5, 0.3 * n + 2)),

  # --- section 11: subcellular localization --------------------------------
  s11_subcellular           = c(14, 7.5),
  s11_subcellular_map       = c(9, 7.5),
  s11_subcellular_bars      = c(7, 7.5)
)

# Look a key up and return c(width, height) in inches. Function entries receive
# the extra arguments the section passed through save_plot().
fig_size <- function(key, ...) {
  if (is.null(key) || !key %in% names(FIGURE_SIZES)) {
    stop("Unknown figure-size key: ", key %||% "<NULL>",
         ". Add it to pipeline/R/figure_sizes.R.", call. = FALSE)
  }
  value <- FIGURE_SIZES[[key]]
  size <- if (is.function(value)) value(...) else as.numeric(value)
  if (length(size) != 2L || any(!is.finite(size)) || any(size <= 0)) {
    stop("Figure size for '", key, "' must be two positive numbers (width, height).",
         call. = FALSE)
  }
  size
}

# Re-read the registry from disk. run_all() calls this at the start of every
# invocation, so a width or height edited in this file takes effect on the next
# run without re-sourcing run_all.R (a session keeps the values it loaded first,
# which is the one thing that can make an edited number appear to be ignored).
reload_figure_sizes <- function(project_root = ".", file = NULL, envir = globalenv()) {
  file <- file %||% file.path(project_root, "pipeline", "R", "figure_sizes.R")
  if (!file.exists(file)) {
    warning("Figure-size registry not found at ", file, "; keeping the loaded values.",
            call. = FALSE)
    return(invisible(FIGURE_SIZES))
  }
  sys.source(file, envir = envir, keep.source = FALSE)
  invisible(get("FIGURE_SIZES", envir = envir))
}

# Every size actually used in a run is recorded here and written to
# provenance/figure_sizes.csv, so the fixed sizes the paper relies on can be
# checked against the registry instead of being inferred from the images.
FIGURE_SIZE_LOG <- new.env(parent = emptyenv())

reset_figure_size_log <- function() {
  assign("rows", list(), envir = FIGURE_SIZE_LOG)
  invisible(NULL)
}

record_figure_size <- function(key, size, path = NA_character_) {
  if (is.null(key) || length(size) != 2L) return(invisible(NULL))
  rows <- get0("rows", envir = FIGURE_SIZE_LOG, ifnotfound = list())
  rows[[length(rows) + 1L]] <- data.frame(
    key = as.character(key),
    width_in = as.numeric(size[[1L]]),
    height_in = as.numeric(size[[2L]]),
    file = basename(as.character(path %||% NA_character_)),
    stringsAsFactors = FALSE)
  assign("rows", rows, envir = FIGURE_SIZE_LOG)
  invisible(NULL)
}

figure_size_log <- function() {
  rows <- get0("rows", envir = FIGURE_SIZE_LOG, ifnotfound = list())
  if (!length(rows)) {
    return(data.frame(key = character(0), width_in = numeric(0),
                      height_in = numeric(0), file = character(0),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  out[order(out$key, out$file), , drop = FALSE]
}
