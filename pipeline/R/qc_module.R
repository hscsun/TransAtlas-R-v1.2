# Sample-level quality control: a port of the QC diagnostics that the original
# scripts produced implicitly, split out so it can be re-run on its own.
# Base-graphics functions draw onto the current device, so the device is opened
# before the call; otherwise the figure lands in a stray Rplots.pdf.

run_qc_module <- function(cfg, deseq_obj, counts_obj, meta, out_dir) {
  out_dir <- safe_dir(out_dir)
  enabled <- function(flag) isTRUE(cfg$quality_control[[flag]])
  if (!any(vapply(cfg$quality_control, isTRUE, logical(1)))) return(invisible(out_dir))
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required for QC.", call. = FALSE)
  dds <- deseq_obj$dds
  vsd <- deseq_obj$vsd
  count_mat <- counts_obj$counts
  group_col <- cfg$analysis$group_column %||% "condition"

  size_factors <- data.frame(sample_id = names(DESeq2::sizeFactors(dds)),
                             size_factor = as.numeric(DESeq2::sizeFactors(dds)),
                             library_size = colSums(count_mat[, names(DESeq2::sizeFactors(dds)), drop = FALSE]))
  utils::write.csv(size_factors, file.path(out_dir, "sample_size_factors.csv"), row.names = FALSE)

  if (enabled("run_library_size")) {
    lib_plot <- ggplot2::ggplot(size_factors, ggplot2::aes(sample_id, library_size, fill = sample_id)) +
      ggplot2::geom_col(show.legend = FALSE) +
      ggplot2::labs(x = NULL, y = "Raw library size", title = "Sequencing depth") +
      publication_theme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_plot(lib_plot, file.path(out_dir, "qc_library_size.png"), key = "qc_library_size")
  }

  if (enabled("run_count_distribution")) {
    long_counts <- as.data.frame(count_mat)
    long_counts$gene_id <- rownames(count_mat)
    if (requireNamespace("tidyr", quietly = TRUE)) {
      long_counts <- tidyr::pivot_longer(long_counts, -gene_id, names_to = "sample_id", values_to = "count")
    } else {
      long_counts <- reshape(long_counts, varying = colnames(count_mat), v.names = "count",
                             timevar = "sample_id", times = colnames(count_mat), direction = "long")
    }
    count_plot <- ggplot2::ggplot(long_counts, ggplot2::aes(log10(count + 1), colour = sample_id)) +
      ggplot2::geom_density(linewidth = 0.45, show.legend = FALSE) +
      ggplot2::labs(x = "log10(raw count + 1)", y = "Density", title = "Count distributions") +
      publication_theme()
    save_plot(count_plot, file.path(out_dir, "qc_count_distributions.png"), key = "qc_count_distributions")
  }

  vsd_mat <- SummarizedExperiment::assay(vsd)
  sample_cor <- stats::cor(vsd_mat, method = "pearson")
  utils::write.csv(sample_cor, file.path(out_dir, "sample_correlation.csv"))
  dist_mat <- as.matrix(stats::dist(t(vsd_mat)))
  ann <- data.frame(group = meta[colnames(vsd_mat), group_col])
  rownames(ann) <- colnames(vsd_mat)
  if (enabled("run_sample_correlation")) {
    sz <- fig_size("qc_sample_correlation")
    record_figure_size("qc_sample_correlation", sz, "qc_sample_correlation.png")
    pheatmap::pheatmap(sample_cor, filename = file.path(out_dir, "qc_sample_correlation.png"),
                       width = sz[1], height = sz[2], border_color = NA, main = "VST sample correlation",
                       fontsize = 9, family = publication_font(), annotation_col = ann)
  }

  if (enabled("run_sample_distance_heatmap")) {
    sz <- fig_size("qc_sample_distance")
    record_figure_size("qc_sample_distance", sz, "qc_sample_distance.png")
    pheatmap::pheatmap(dist_mat, filename = file.path(out_dir, "qc_sample_distance.png"),
                       width = sz[1], height = sz[2], border_color = NA, main = "VST sample distance",
                       fontsize = 9, family = publication_font(), annotation_col = ann)
  }

  if (enabled("run_pca")) {
    pca_data <- as.data.frame(DESeq2::plotPCA(vsd, intgroup = group_col, returnData = TRUE))
    percent_var <- round(100 * attr(pca_data, "percentVar"), 1)
    pca_data$sample_id <- rownames(pca_data)
    pca_data$group <- meta[pca_data$sample_id, group_col]
    pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(PC1, PC2, colour = group, label = sample_id)) +
      ggplot2::geom_point(size = 3) + ggplot2::geom_text(vjust = -0.8, size = 3.8, show.legend = FALSE) +
      ggplot2::labs(x = paste0("PC1 (", percent_var[1], "%)"), y = paste0("PC2 (", percent_var[2], "%)"),
                    title = "Principal component analysis") + publication_theme(base_size = 13)
    save_plot(pca_plot, file.path(out_dir, "qc_pca.png"), key = "qc_pca")
    utils::write.csv(pca_data, file.path(out_dir, "qc_pca_coordinates.csv"), row.names = FALSE)
  }

  if (enabled("run_dispersion_plot")) {
    sz <- fig_size("qc_dispersion")
    record_figure_size("qc_dispersion", sz, "qc_dispersion.png")
    grDevices::png(file.path(out_dir, "qc_dispersion.png"), width = sz[1] * 300, height = sz[2] * 300, res = 300)
    dispersion_device <- grDevices::dev.cur()
    on.exit(close_plot_device(dispersion_device), add = TRUE)
    DESeq2::plotDispEsts(dds)
    close_plot_device(dispersion_device)
  }

  if (enabled("run_ma_plot")) {
    sz <- fig_size("qc_ma_plot")
    record_figure_size("qc_ma_plot", sz, "qc_ma_plot.png")
    grDevices::png(file.path(out_dir, "qc_ma_plot.png"), width = sz[1] * 300, height = sz[2] * 300, res = 300)
    ma_device <- grDevices::dev.cur()
    on.exit(close_plot_device(ma_device), add = TRUE)
    DESeq2::plotMA(dds, ylim = c(-5, 5), alpha = cfg$analysis$alpha %||% 0.05)
    close_plot_device(ma_device)
  }

  invisible(out_dir)
}
