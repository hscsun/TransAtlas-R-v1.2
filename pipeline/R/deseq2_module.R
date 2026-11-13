# Core DESeq2 workflow: raw-count checks, filtering, differential expression.
# Sample-level QC figures live in qc_module.R.

run_deseq2 <- function(cfg, counts_obj, meta, out_dir) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Package 'DESeq2' is required.", call. = FALSE)
  safe_dir(out_dir)
  count_mat <- counts_obj$counts
  min_count <- as.integer(cfg$analysis$min_count %||% 10L)
  min_samples <- as.integer(cfg$analysis$min_samples %||% 2L)
  keep <- rowSums(count_mat >= min_count) >= min_samples
  filtered <- count_mat[keep, , drop = FALSE]
  utils::write.csv(data.frame(gene_id = rownames(count_mat), kept = keep),
                   file.path(out_dir, "gene_filtering.csv"), row.names = FALSE)
  if (nrow(filtered) < 100L) warning("Fewer than 100 genes remain after low-count filtering.")

  design <- stats::as.formula(cfg$analysis$design)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(filtered), colData = meta, design = design)
  dds <- DESeq2::DESeq(dds, quiet = TRUE, parallel = FALSE, minReplicatesForReplace = Inf)
  transform_method <- "vst"
  vsd <- tryCatch(DESeq2::vst(dds, blind = FALSE), error = function(e) {
    message("vst unavailable (", conditionMessage(e), "); using varianceStabilizingTransformation.")
    transform_method <<- "varianceStabilizingTransformation"
    DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  })
  write_status(file.path(out_dir, "status_transform.json"), "expression_transform", "completed",
               details = list(method = transform_method, genes = nrow(dds)))
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  norm_df <- data.frame(gene_id = rownames(norm_counts), norm_counts, check.names = FALSE)
  utils::write.csv(norm_df, file.path(out_dir, "normalized_counts.csv"), row.names = FALSE)

  group_col <- cfg$analysis$group_column %||% "condition"
  if (!group_col %in% names(meta)) stop("Configured group column is missing from metadata: ", group_col, call. = FALSE)
  if (length(levels(meta[[group_col]])) != 2L) stop("This module requires exactly two groups.", call. = FALSE)

  contrast_results <- list()
  contrast_cfg <- cfg$analysis$contrasts %||% list()
  if (length(contrast_cfg) == 0L) stop("At least one contrast is required.", call. = FALSE)
  if (is.list(contrast_cfg) && !is.null(contrast_cfg$name)) contrast_cfg <- list(contrast_cfg)
  for (contrast in contrast_cfg) {
    contrast_name <- contrast$name
    contrast_variable <- contrast$variable %||% group_col
    numerator <- contrast$numerator
    denominator <- contrast$denominator
    if (!contrast_variable %in% names(meta)) stop("Contrast variable is missing from metadata: ", contrast_variable, call. = FALSE)
    if (!is.factor(meta[[contrast_variable]])) meta[[contrast_variable]] <- factor(meta[[contrast_variable]])
    if (!all(c(numerator, denominator) %in% levels(meta[[contrast_variable]]))) {
      warning("Skipping contrast ", contrast_name, ": levels not present in metadata.")
      next
    }
    res <- DESeq2::results(dds, contrast = c(contrast_variable, numerator, denominator),
                           alpha = cfg$analysis$alpha %||% 0.05,
                           independentFiltering = isTRUE(cfg$analysis$independent_filtering))
    res_df <- result_to_df(res, counts_obj$annotations, cfg$input$gene_id_column,
                           cfg$input$gene_symbol_column)
    res_df$contrast <- contrast_name
    res_df$baseMean <- as.numeric(res_df$baseMean)
    res_df$padj_safe <- p_adjusted_or_one(res_df$padj)
    res_df$significant <- !is.na(res_df$padj) & res_df$padj < (cfg$analysis$alpha %||% 0.05) &
      abs(res_df$log2FoldChange) >= (cfg$analysis$lfc_threshold %||% 0.75)
    res_df$group_column <- contrast_variable
    res_df$numerator <- numerator
    res_df$denominator <- denominator
    safe_write_csv(res_df, file.path(out_dir, paste0("DESeq2_", contrast_name, ".csv")))

    # Save shrunken estimates as a separate result for downstream inspection.
    # Screening and the compatibility GSEA ranking use the unshrunken results.
    shrunk <- NULL
    shrinkage_method <- "none"
    shrinkage_requested <- cfg$analysis$lfc_shrinkage %||% "none"
    if (shrinkage_requested == "apeglm" && requireNamespace("apeglm", quietly = TRUE)) {
      coef_name <- grep(paste0("^", contrast_variable, "_", numerator, "_vs_", denominator, "$"),
                        DESeq2::resultsNames(dds), value = TRUE)
      shrunk <- if (length(coef_name) == 1L) {
        tryCatch(DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm"), error = record_optional_error)
      } else {
        tryCatch(DESeq2::lfcShrink(dds, contrast = c(contrast_variable, numerator, denominator),
                                   type = "apeglm"), error = record_optional_error)
      }
      if (!is.null(shrunk)) shrinkage_method <- "apeglm"
    }
    if (is.null(shrunk) && (cfg$analysis$lfc_shrinkage %||% "none") != "none") {
      shrunk <- tryCatch(DESeq2::lfcShrink(dds, contrast = c(contrast_variable, numerator, denominator),
                                           type = "normal"), error = record_optional_error)
      if (!is.null(shrunk)) shrinkage_method <- "normal_fallback"
    }
    if (!is.null(shrunk)) {
      shrink_df <- result_to_df(shrunk, counts_obj$annotations, cfg$input$gene_id_column,
                                cfg$input$gene_symbol_column)
      shrink_df$contrast <- contrast_name
      safe_write_csv(shrink_df, file.path(out_dir, paste0("DESeq2_", contrast_name, "_lfc_shrunk.csv")))
    }
    contrast_results[[contrast_name]] <- list(result = res_df, shrunk = shrunk,
                                              shrinkage_method = shrinkage_method,
                                              numerator = numerator, denominator = denominator)
  }
  list(dds = dds, vsd = vsd, normalized_counts = norm_counts,
       contrast_results = contrast_results, filtered_counts = filtered, keep = keep,
       transform_method = transform_method)
}

# Recreate upstream tables when a valid fit is reused in a new output folder.
write_cached_deseq2_tables <- function(cfg, counts_obj, deseq_obj, out_dir) {
  safe_dir(out_dir)
  safe_write_csv(data.frame(gene_id = rownames(counts_obj$counts), kept = deseq_obj$keep),
                 file.path(out_dir, "gene_filtering.csv"))
  norm_counts <- deseq_obj$normalized_counts
  safe_write_csv(data.frame(gene_id = rownames(norm_counts), norm_counts, check.names = FALSE),
                 file.path(out_dir, "normalized_counts.csv"))
  for (contrast_name in names(deseq_obj$contrast_results)) {
    contrast <- deseq_obj$contrast_results[[contrast_name]]
    safe_write_csv(contrast$result,
                   file.path(out_dir, paste0("DESeq2_", contrast_name, ".csv")))
    if (!is.null(contrast$shrunk)) {
      shrink_df <- result_to_df(contrast$shrunk, counts_obj$annotations,
                               cfg$input$gene_id_column, cfg$input$gene_symbol_column)
      shrink_df$contrast <- contrast_name
      safe_write_csv(shrink_df,
                     file.path(out_dir, paste0("DESeq2_", contrast_name, "_lfc_shrunk.csv")))
    }
  }
  write_status(file.path(out_dir, "status_transform.json"), "expression_transform", "completed",
               details = list(method = deseq_obj$transform_method %||% "cached_transform",
                              source = "cached_fit", genes = nrow(deseq_obj$dds)))
  invisible(out_dir)
}
