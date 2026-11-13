# Section 1 - differential-expression screening and volcano plots.
#
# Port of input/AnaDGEsS1.R. The original script chose a DESeq2 result table
# and an FPKM matrix through file.choose(), hard-coded an output directory under
# /Users/changye, and drew its plots with base graphics. Here the inputs come
# from the configuration, the threshold ladder and the five-way gene
# classification are preserved exactly, and every figure is written to the
# section directory instead of a plot window.

run_screening <- function(cfg, ctx, out_dir) {
  out_dir <- safe_dir(out_dir)
  status <- list()
  symbol_col <- ctx$symbols
  gene_id_col <- cfg$input$gene_id_column
  fcap <- display_lfc_limit(cfg, ctx$de$log2FoldChange)

  # The five-way classification and the threshold ladder live in context.R,
  # because the driver rebuilds them when a run starts after this section.
  tab <- screening_tables(cfg, ctx)
  mylist <- tab$mylist
  class_summary <- tab$class_summary
  thresholds <- tab$thresholds
  fold_change <- tab$fold_change
  fold_change_raw <- tab$fold_change_raw
  base_mean <- tab$base_mean
  exp10 <- tab$exp10
  signafi <- tab$signafi
  # Eight genes can reach padj = 0, which makes -log10(padj) infinite. ggplot
  # silently drops those rows and draws the remaining labels' leader lines to
  # the panel edge, so the tallest triangles look cut in half. Cap the plotted
  # value at the largest finite score: the genes stay at the top of the panel,
  # just above every other point, and remain visible.
  finite_sig <- signafi[is.finite(signafi) & signafi < 500]
  sig_cap <- if (length(finite_sig)) max(finite_sig, na.rm = TRUE) else 0
  signafi_plot <- if (is.finite(sig_cap) && sig_cap > 0) pmin(signafi, sig_cap * 1.01) else signafi
  gene_symbol <- tab$gene_symbol
  gene_id <- tab$gene_id
  use_expr_filter <- tab$use_expr_filter
  expression_threshold <- tab$expression_threshold
  sig_cut <- tab$sig_cut
  classify <- tab$classify

  for (i in seq_along(thresholds)) {
    safe_write_csv(mylist[[i]], file.path(out_dir, sprintf("S1_DEG_table_lfc_%.2f.csv", thresholds[i])))
    status[[i]] <- class_summary[i, , drop = FALSE]
  }
  safe_write_csv(class_summary, file.path(out_dir, "S1_class_summary.csv"))

  # Document what the original low-expression rule would have removed.
  low_probe <- exp10 <= 2
  diag <- data.frame(
    rule = "log10(baseMean) <= 2",
    genes_affected = sum(low_probe),
    genes_significant = sum(low_probe & !is.na(signafi) & signafi >= sig_cut),
    genes_significant_abs_lfc_ge_1 = sum(low_probe & !is.na(signafi) & signafi >= sig_cut & abs(fold_change_raw) >= 1),
    filter_active = use_expr_filter,
    note = paste0("AnaDGEsS1.R hid genes with log10(baseMean) <= 2 from the up/down classes. ",
                  "The rule is off by default, so the FPKM floor (min_fpkm) is the only ",
                  "expression threshold; pass run_all(expression_log10_threshold = 2) to restore it.")
  )
  safe_write_csv(diag, file.path(out_dir, "S1_low_expression_rule_effect.csv"))

  # Workbook with one sheet per threshold, matching the original sData.xlsx.
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(mylist, file.path(out_dir, "sData.xlsx"), col_names = TRUE)
  }

  # Threshold ladder overview.
  ladder_long <- rbind(
    data.frame(lfc_limit = class_summary$lfc_limit, direction = "Up", n = class_summary$up),
    data.frame(lfc_limit = class_summary$lfc_limit, direction = "Down", n = class_summary$down)
  )
  ladder_long$lfc_limit <- factor(ladder_long$lfc_limit, levels = thresholds)
  ladder_plot <- ggplot2::ggplot(ladder_long, ggplot2::aes(lfc_limit, n, fill = direction)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::geom_text(ggplot2::aes(label = n), position = ggplot2::position_dodge(width = 0.9),
                       vjust = -0.4, size = 2.6) +
    ggplot2::scale_fill_manual(values = c(Up = "#B2182B", Down = "#2166AC")) +
    ggplot2::labs(x = expression("|log"[2] * " fold change| threshold"), y = "Number of genes",
                  fill = NULL, title = "Differentially expressed genes across screening thresholds") +
    publication_theme()
  save_plot(ladder_plot, file.path(out_dir, "S1_threshold_counts.png"),
            key = "s1_threshold_counts")

  # --- original base-graphics volcano and FC-Exp panels ----------------------
  # cairo_pdf resolves the system sans-serif family; the PostScript pdf()
  # device only knows the base-14 fonts and would drop the family.
  volcano_pdf <- file.path(out_dir, sprintf("S1_volcano_and_FCexp_%s.pdf", ctx$contrast_name))
  sz <- fig_size("s1_volcano_pdf")
  record_figure_size("s1_volcano_pdf", sz, basename(volcano_pdf))
  grDevices::cairo_pdf(volcano_pdf, width = sz[1], height = sz[2], bg = "transparent")
  pdf_device <- grDevices::dev.cur()
  on.exit(close_plot_device(pdf_device), add = TRUE)
  for (i in seq_along(thresholds)) {
    limit <- thresholds[i]
    cls <- classify(limit)
    display <- pmin(pmax(fold_change, -fcap), fcap)
    ymax <- max(signafi[is.finite(signafi) & signafi < 500], na.rm = TRUE)
    # 5% head-room above the tallest point.
    ymax <- max(ymax, sig_cut) * 1.05
    xlim <- range(display, na.rm = TRUE) * 1.05

    graphics::par(oma = c(0, 1, 0, 0), family = publication_font())
    graphics::plot(NA, xlim = xlim, ylim = c(-2.5, ymax),
                   xlab = expression("log"[2] * "(Fold Change)"), ylab = expression("-log"[10] * "(p.adjust)"),
                   main = sprintf("Volcano plot, |log2FC| >= %s", limit),
                   cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, font.axis = 2, cex.main = 1.6)
    for (cat in c("low_expression", "not_significant", "intermediate", "up", "down")) {
      idx <- which(cls == cat)
      if (!length(idx)) next
      if (cat == "low_expression") {
        graphics::points(display[idx], signafi_plot[idx], pch = 16, cex = 0.2, col = "#999999")
      } else if (cat == "not_significant") {
        graphics::points(display[idx], signafi_plot[idx], pch = 19, cex = 0.2, col = "#33CCFF")
      } else if (cat == "intermediate") {
        graphics::points(display[idx], signafi_plot[idx], pch = 19, cex = 0.2, col = "#FF99FF")
      } else if (cat == "up") {
        graphics::points(display[idx], signafi_plot[idx], pch = 2, cex = 0.35, lwd = 1.5, col = "#FF0000")
      } else {
        graphics::points(display[idx], signafi_plot[idx], pch = 6, cex = 0.35, lwd = 1.5, col = "#0000CC")
      }
    }
    graphics::abline(h = sig_cut, lty = 2, col = "grey")
    graphics::abline(v = c(limit, -limit), lty = 2, col = "grey")
    graphics::legend("topleft", legend = c(paste("Up:", sum(cls == "up")), paste("Down:", sum(cls == "down"))),
                     bty = "n", cex = 1.1)

    graphics::par(oma = c(0, 1, 0, 0))
    graphics::plot(NA, xlim = xlim, ylim = c(0, max(exp10, na.rm = TRUE) * 1.05),
                   xlab = expression("log"[2] * "(Fold Change)"), ylab = expression("log"[10] * "(baseMean)"),
                   main = sprintf("FC-Exp plot, |log2FC| >= %s", limit),
                   cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, font.axis = 2, cex.main = 1.6)
    for (cat in c("low_expression", "not_significant", "intermediate", "up", "down")) {
      idx <- which(cls == cat)
      if (!length(idx)) next
      if (cat == "low_expression") {
        graphics::points(display[idx], exp10[idx], pch = 16, cex = 0.2, col = "#999999")
      } else if (cat == "not_significant") {
        graphics::points(display[idx], exp10[idx], pch = 19, cex = 0.2, col = "#33CCFF")
      } else if (cat == "intermediate") {
        graphics::points(display[idx], exp10[idx], pch = 19, cex = 0.2, col = "#FF99FF")
      } else if (cat == "up") {
        graphics::points(display[idx], exp10[idx], pch = 2, cex = 0.35, lwd = 1.5, col = "#FF0000")
      } else {
        graphics::points(display[idx], exp10[idx], pch = 6, cex = 0.35, lwd = 1.5, col = "#0000CC")
      }
    }
    if (use_expr_filter) graphics::abline(h = expression_threshold, lty = 2, col = "grey")
    graphics::abline(v = c(limit, -limit), lty = 2, col = "grey")
    graphics::legend("topleft", legend = c(paste("Up:", sum(cls == "up")), paste("Down:", sum(cls == "down"))),
                     bty = "n", cex = 1.1)
  }
  close_plot_device(pdf_device)
  status$volcano_pdf <- basename(volcano_pdf)

  # --- publication volcano and expression-effect plot -------------------------
  primary <- as.numeric(cfg$analysis$lfc_threshold %||% 0.75)
  # The two screening panels show the whole tested universe, not only the genes
  # that survived the FPKM floor: the genes the floor removed are drawn in grey
  # so the effect of the pre-screen is visible rather than implied by an absent
  # cloud of points.
  all_de <- ctx$de_all %||% ctx$de
  all_ids <- as.character(all_de[[gene_id_col]])
  all_signafi <- -log10(all_de$padj)
  all_lfc_raw <- all_de$log2FoldChange
  all_lfc <- cap_log2fc(all_lfc_raw, cfg$analysis$up_cutoff)
  all_bm <- all_de$baseMean
  all_sym <- as.character(all_de[[symbol_col]])
  min_fpkm <- suppressWarnings(as.numeric(ctx$min_fpkm %||% NA))
  fpkm <- if (!is.null(ctx$expr_fpkm)) as.numeric(ctx$expr_fpkm[all_ids]) else
    rep(NA_real_, length(all_ids))
  below_floor <- is.finite(min_fpkm) & is.finite(fpkm) & fpkm < min_fpkm
  # Classes come from the floor-passing table; every gene the floor removed is
  # labelled low_expression whatever its P value, because it never entered the
  # up/down classes.
  cls_by_id <- stats::setNames(classify(primary), tab$gene_id)
  cls <- unname(cls_by_id[all_ids])
  cls[is.na(cls)] <- "not_significant"
  cls[below_floor] <- "low_expression"
  # Same cap as the volcano draws: -log10(padj) can be infinite.
  finite_all <- all_signafi[is.finite(all_signafi) & all_signafi < 500]
  cap_all <- if (length(finite_all)) max(finite_all, na.rm = TRUE) else 0
  signafi_plot_all <- if (is.finite(cap_all) && cap_all > 0) pmin(all_signafi, cap_all * 1.01) else
    all_signafi
  plot_df <- data.frame(log2FoldChange = all_lfc, raw_lfc = all_lfc_raw,
                        signafi = signafi_plot_all,
                        baseMean = all_bm, fpkm = fpkm, symbol = all_sym,
                        class = cls, stringsAsFactors = FALSE)
  # Labels flag four groups of five genes: the most strongly expressed up- and
  # down-regulated DEGs, and the DEGs with the largest |log2FC| in each
  # direction. Expression ranking uses the mean normalized count; effect ranking
  # uses the uncapped fold change. A gene already picked by expression is not
  # repeated, so each panel carries twenty distinct labels: ten per direction,
  # five of them high-expression and five high-effect.
  n_expr <- as.integer(cfg$analysis$volcano_label_expression_n %||% 5L)
  n_eff <- as.integer(cfg$analysis$volcano_label_effect_n %||% 5L)
  up_i <- which(plot_df$class == "up")
  down_i <- which(plot_df$class == "down")
  top_by <- function(idx, n, score) {
    if (!length(idx) || n <= 0L) return(integer(0))
    utils::head(idx[order(score[idx], decreasing = TRUE)], n)
  }
  up_expr <- top_by(up_i, n_expr, plot_df$baseMean)
  down_expr <- top_by(down_i, n_expr, plot_df$baseMean)
  up_eff <- top_by(setdiff(up_i, up_expr), n_eff, abs(plot_df$raw_lfc))
  down_eff <- top_by(setdiff(down_i, down_expr), n_eff, abs(plot_df$raw_lfc))
  # Alternate directions so a smaller label budget retains both groups.
  interleave <- function(a, b) as.integer(stats::na.omit(as.vector(rbind(
    c(a, rep(NA_integer_, max(0L, length(b) - length(a)))),
    c(b, rep(NA_integer_, max(0L, length(a) - length(b))))))))
  label_idx <- utils::head(c(interleave(up_expr, down_expr), interleave(up_eff, down_eff)),
                          as.integer(cfg$analysis$top_n_labels %||% 20L))
  label_genes <- unique(plot_df$symbol[label_idx])
  label_genes <- label_genes[!is.na(label_genes)]
  plot_df$label <- NA_character_
  label_idx <- label_idx[!duplicated(plot_df$symbol[label_idx])]
  plot_df$label[label_idx] <- plot_df$symbol[label_idx]
  x_limit <- as.numeric(cfg$analysis$volcano_x_limit %||% (max(abs(all_lfc[is.finite(all_lfc)]), 0.5) * 1.05))
  x_cap <- as.numeric(cfg$analysis$up_cutoff %||% NA_real_)
  plot_df$log2FoldChange <- pmin(pmax(plot_df$log2FoldChange, -x_limit), x_limit)
  # Up- and down-regulated genes are told apart by glyph, not only by colour: an
  # upward triangle for up, a downward one for down, as in the original
  # base-graphics volcano (pch 2 / pch 6). The remaining classes stay circular so
  # the two directional classes are the only ones with a distinctive shape.
  class_pal <- c(low_expression = "#BDBDBD", not_significant = "#33CCFF",
                 intermediate = "#FF99FF", up = "#FF0000", down = "#0000CC")
  # Hollow triangles: shapes 2 and 6 are outline-only glyphs, so the interior
  # stays white and only the coloured triangle border carries the direction.
  class_shapes <- c(low_expression = 16, not_significant = 16, intermediate = 16,
                    up = 2, down = 6)
  # Legend order and wording: the floor comes first, then the non-directional
  # classes, then the two calls, whose counts answer "how many genes moved" in
  # the legend itself.
  class_order <- c("low_expression", "not_significant", "intermediate", "up", "down")
  n_up <- sum(plot_df$class == "up")
  n_down <- sum(plot_df$class == "down")
  # Only the two calls carry a count: the floor class is context for the panel,
  # and its size is a property of the filter rather than of the contrast.
  class_labels <- c(low_expression = "Low expression",
                    not_significant = "Not significant",
                    intermediate = "Intermediate",
                    up = sprintf("Up: %d", n_up),
                    down = sprintf("Down: %d", n_down))
  # The genes below the floor are drawn first, smaller and paler: they are
  # context for the calls, not the subject of the panel.
  low_df <- plot_df[plot_df$class == "low_expression", , drop = FALSE]
  main_df <- plot_df[plot_df$class != "low_expression", , drop = FALSE]
  legend_guide <- ggplot2::guide_legend(
    override.aes = list(size = 1.4, stroke = 0.55, alpha = 0.85), order = 1)
  class_scales <- list(
    ggplot2::scale_colour_manual(values = class_pal, name = NULL, breaks = class_order,
                                 labels = class_labels),
    ggplot2::scale_fill_manual(values = class_pal, name = NULL, breaks = class_order,
                               labels = class_labels),
    ggplot2::scale_shape_manual(values = class_shapes, name = NULL, breaks = class_order,
                                labels = class_labels))
  # One short line each: a longer note runs off the right edge of a 6.5-inch
  # panel at 14 pt.
  floor_note <- if (is.finite(min_fpkm))
    paste0("Grey = below FPKM ", format(min_fpkm, trim = TRUE)) else
    "FPKM floor disabled"
  volcano <- ggplot2::ggplot(main_df, ggplot2::aes(log2FoldChange, signafi)) +
    ggplot2::geom_point(data = low_df, ggplot2::aes(log2FoldChange, signafi,
                                                    colour = class, fill = class, shape = class),
                        size = 1.1, stroke = 0.4, alpha = 0.35) +
    ggplot2::geom_point(ggplot2::aes(colour = class, fill = class, shape = class),
                        size = 1.4, stroke = 0.55, alpha = 0.85) +
    class_scales +
    # The points are drawn at size 1.4; without this the merged class legend
    # falls back to ggplot2's much larger default key, so the legend glyphs look
    # heavier than every point in the panel. Colour, fill and shape must share
    # one guide object, or the two identical legends are drawn separately.
    ggplot2::guides(colour = legend_guide,
                    fill = ggplot2::guide_legend(order = 1),
                    shape = ggplot2::guide_legend(order = 1)) +
    ggplot2::geom_vline(xintercept = c(-primary, primary), linetype = 2, colour = "grey40") +
    ggplot2::geom_hline(yintercept = sig_cut, linetype = 2, colour = "grey40") +
    ggplot2::labs(x = if (is.finite(x_cap)) paste0("log2 fold change (capped at ±", x_cap, ")") else
                    expression("log"[2] * " fold change"),
                  y = expression("-log"[10] * " adjusted P"),
                  title = paste0("Volcano plot: ", ctx$contrast_name),
                  subtitle = paste0(floor_note, "; dashed = |log2FC| ", primary,
                                    ", P ", cfg$analysis$alpha %||% 0.05)) +
    # A small expansion keeps labels placed at the left and right edge of the
    # panel (the genes with the largest |log2FC|) inside the plotted area.
    ggplot2::scale_x_continuous(limits = c(-x_limit, x_limit),
                                expand = ggplot2::expansion(mult = 0.06)) +
    # The axis extends 5% of the data range above the top point and 5.75%
    # below it, so the tallest triangles sit close to the top without clipping.
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.0575, 0.05))) +
    publication_theme(base_size = 14)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    # Repel only the labelled rows: handing the layer all 20,000 genes makes
    # ggrepel optimise thousands of empty labels and dominates the run time.
    lab_df <- plot_df[!is.na(plot_df$label), , drop = FALSE]
    volcano <- volcano + ggrepel::geom_text_repel(
      data = lab_df, ggplot2::aes(label = label), size = 3.35, max.overlaps = Inf,
      box.padding = ggplot2::unit(0.35, "lines"),
      point.padding = ggplot2::unit(0.25, "lines"),
      min.segment.length = 0, segment.size = 0.25, segment.colour = "grey55",
      force = 3, max.iter = 20000)
  }
  save_plot(volcano, file.path(out_dir, sprintf("S1_volcano_%s.png", ctx$contrast_name)),
            key = "s1_volcano")

  # The expression axis carries the floor directly, so the cut-off can be drawn
  # where it belongs. FPKM values are floored at 0.01 so a gene with no reads
  # still has a finite position at the left edge of a log axis.
  use_fpkm_axis <- any(is.finite(plot_df$fpkm))
  plot_df$fpkm_axis <- if (use_fpkm_axis) plot_df$fpkm else plot_df$baseMean
  expression_axis_label <- if (use_fpkm_axis) "Maximum group-mean FPKM (log10 scale)" else
    "Mean normalized count (log10 scale)"
  plot_df$fpkm_axis[!is.finite(plot_df$fpkm_axis)] <- NA_real_
  effect_df <- plot_df[is.finite(plot_df$fpkm_axis), , drop = FALSE]
  effect_df$fpkm_axis <- pmax(effect_df$fpkm_axis, 0.01)
  low_df <- effect_df[effect_df$class == "low_expression", , drop = FALSE]
  main_df <- effect_df[effect_df$class != "low_expression", , drop = FALSE]
  effect <- ggplot2::ggplot(main_df, ggplot2::aes(fpkm_axis, log2FoldChange)) +
    ggplot2::geom_point(data = low_df, ggplot2::aes(fpkm_axis, log2FoldChange,
                                                    colour = class, fill = class, shape = class),
                        size = 1.1, stroke = 0.4, alpha = 0.35) +
    ggplot2::geom_point(ggplot2::aes(colour = class, fill = class, shape = class),
                        size = 1.4, stroke = 0.55, alpha = 0.85) +
    # Plain decimal labels: the default log scale switches the whole axis to
    # scientific notation (1e-01, 1e+01), which hides where the floor sits.
    ggplot2::scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000),
                           labels = function(x) format(x, scientific = FALSE,
                                                       drop0trailing = TRUE)) +
    class_scales +
    ggplot2::guides(colour = legend_guide,
                    fill = ggplot2::guide_legend(order = 1),
                    shape = ggplot2::guide_legend(order = 1)) +
    ggplot2::geom_hline(yintercept = c(-primary, primary), linetype = 2, colour = "grey40") +
    ggplot2::labs(x = expression_axis_label,
                  y = if (is.finite(x_cap)) paste0("log2 fold change (capped at ±", x_cap, ")") else
                    expression(log[2] * " fold change"),
                  title = paste0("Expression vs effect size: ", ctx$contrast_name),
                  subtitle = paste0(floor_note, "; dashed = ",
                                    if (is.finite(min_fpkm)) "floor and " else "", "|log2FC| ", primary)) +
    publication_theme(base_size = 14)
  if (is.finite(min_fpkm)) {
    effect <- effect + ggplot2::geom_vline(xintercept = min_fpkm, linetype = 2,
                                           colour = "grey40")
  }
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    # Same twenty genes as the volcano: five high-expression and five
    # high-effect DEGs per direction. Only the labelled rows are passed to
    # ggrepel so the layer never optimises thousands of empty labels.
    lab_df <- effect_df[!is.na(effect_df$label), , drop = FALSE]
    effect <- effect + ggrepel::geom_text_repel(
      data = lab_df,
      ggplot2::aes(x = fpkm_axis, y = log2FoldChange, label = label),
      size = 3.35, max.overlaps = Inf,
      box.padding = ggplot2::unit(0.35, "lines"),
      point.padding = ggplot2::unit(0.25, "lines"),
      min.segment.length = 0, segment.size = 0.25, segment.colour = "grey55",
      force = 3, max.iter = 20000)
  }
  save_plot(effect, file.path(out_dir, sprintf("S1_expression_effect_%s.png", ctx$contrast_name)),
            key = "s1_expression_effect")

  # --- heat maps --------------------------------------------------------------
  # Exactly two, so the same matrix is never drawn twice: one for the top DEGs of
  # the primary threshold (gene names shown) and one global map holding *every*
  # gene passing |log2FC| >= the primary threshold, split into a down-regulated
  # and an up-regulated block. The row-ratio and heatmap.2 variants of the
  # original script showed the same values a second and third time and are gone.
  if (!is.null(ctx$expr_matrix)) {
    ann_col <- data.frame(group = ctx$meta[colnames(ctx$expr_matrix), ctx$group_col])
    rownames(ann_col) <- colnames(ctx$expr_matrix)
    deg_rows <- mylist[[which.min(abs(thresholds - primary))]]

    zscore_matrix <- function(rows) {
      rows <- rows[!duplicated(rows$GenesSymbol), , drop = FALSE]
      symbol_lookup <- match(rows$GenesSymbol, as.character(ctx$annotations[[symbol_col]]))
      gene_idx <- match(as.character(ctx$annotations[[gene_id_col]])[symbol_lookup],
                        rownames(ctx$expr_matrix))
      hm <- ctx$expr_matrix[gene_idx, , drop = FALSE]
      rownames(hm) <- rows$GenesSymbol
      hm <- hm[rowSums(hm, na.rm = TRUE) > 0, , drop = FALSE]
      if (nrow(hm) < 2L) return(NULL)
      z <- t(scale(t(hm)))
      z[!is.finite(z)] <- 0
      z
    }

    # Select genes by their true effect size before applying any display scale.
    top <- deg_rows[order(-abs(deg_rows$GenesFoldChange_raw)), , drop = FALSE]
    top <- utils::head(top, cfg$analysis$top_n_heatmap %||% 40L)
    z_top <- zscore_matrix(top)
    if (!is.null(z_top)) {
      sz_hm <- fig_size("s1_top_deg_heatmap")
      record_figure_size("s1_top_deg_heatmap", sz_hm, "S1_topDEG_heatmap_rowz.png")
      pheatmap::pheatmap(z_top, scale = "none", show_rownames = TRUE, angle_col = 90,
                         color = heatmap_palette(),
                         annotation_col = ann_col, fontsize = 13, fontsize_row = 11,
                         fontsize_col = 13,
                         filename = file.path(out_dir, "S1_topDEG_heatmap_rowz.png"),
                         width = sz_hm[1], height = sz_hm[2], border_color = NA)
    }

    glob <- deg_rows[order(deg_rows$GenesFoldChange_raw), , drop = FALSE]
    z_glob <- zscore_matrix(glob)
    if (!is.null(z_glob) && nrow(z_glob) >= 4L) {
      direction_by_symbol <- stats::setNames(as.character(deg_rows$Direction), deg_rows$GenesSymbol)
      ann_row <- data.frame(
        direction = factor(direction_by_symbol[rownames(z_glob)], levels = c("Down", "Up")))
      rownames(ann_row) <- rownames(z_glob)
      n_down <- sum(ann_row$direction == "Down", na.rm = TRUE)
      gap <- if (n_down > 0L && n_down < nrow(z_glob)) n_down else NULL
      sz_glob <- fig_size("s1_global_deg_heatmap")
      glob_name <- sprintf("S1_DEG_global_heatmap_lfc%.2f.png", primary)
      record_figure_size("s1_global_deg_heatmap", sz_glob, glob_name)
      pheatmap::pheatmap(z_glob, scale = "none", cluster_rows = FALSE, cluster_cols = TRUE,
                         show_rownames = FALSE, angle_col = 90,
                         color = heatmap_palette(),
                         annotation_col = ann_col, annotation_row = ann_row,
                         gaps_row = gap, fontsize = 9, fontsize_col = 9,
                         treeheight_row = 0, border_color = NA,
                         filename = file.path(out_dir, glob_name),
                         width = sz_glob[1], height = sz_glob[2])
      status$global_heatmap_genes <- nrow(z_glob)
    }
  }

  # --- chromosome landscape ---------------------------------------------------
  # The maps use their own, broader thresholds so the genome-wide shape of the
  # response is visible at three effect-size levels, while every table stays at
  # the reporting threshold above. `chromosome_landscape()` draws one track per
  # chromosome with stem height proportional to |log2FC|.
  chrom_status <- tryCatch({
    orgdb_name <- cfg$species$orgdb
    if (is.null(orgdb_name) || !requireNamespace(orgdb_name, quietly = TRUE)) {
      stop("OrgDb ", orgdb_name %||% "<none>", " is not installed")
    }
    orgdb <- getExportedValue(orgdb_name, orgdb_name)
    chrom_lfcs <- as.numeric(unlist(cfg$analysis$chromosome_lfc_thresholds %||%
                                      c(0.5, 0.75, 1)))
    chrom_lfcs <- chrom_lfcs[is.finite(chrom_lfcs)]
    if (!length(chrom_lfcs)) chrom_lfcs <- 0.75
    cap <- display_lfc_limit(cfg, ctx$de$log2FoldChange)
    status_rows <- vector("list", length(chrom_lfcs))
    for (j in seq_along(chrom_lfcs)) {
      chrom_lfc <- chrom_lfcs[j]
      tag <- sprintf("lfc%.2f", chrom_lfc)
      chrom_rows <- mylist[[which.min(abs(thresholds - chrom_lfc))]]
      chrom_rows <- chrom_rows[chrom_rows$Direction %in% c("Up", "Down"), , drop = FALSE]
      chrom_rows <- chrom_rows[!duplicated(chrom_rows$GenesSymbol), , drop = FALSE]
      feats <- data.frame(gene = as.character(chrom_rows$GenesSymbol),
                          lfc = as.numeric(chrom_rows$GenesFoldChange_raw),
                          stringsAsFactors = FALSE)
      land <- chromosome_landscape(
        feats, orgdb, cfg$species$common_name, cap = cap,
        title = paste0("Chromosomal distribution of differentially expressed genes: ",
                       ctx$contrast_name),
        subtitle = paste0("|log2FC| >= ", chrom_lfc, "; stems point up for higher and down for lower ",
                          "expression in ", ctx$numerator, " vs ", ctx$denominator,
                          "; stem height is proportional to |log2FC| (capped at ", cap, ")"),
        legend_title = "DEG direction")
      if (is.null(land)) stop("No chromosome coordinates for |log2FC| >= ", chrom_lfc)
      save_plot(land$plot,
                file.path(out_dir, sprintf("S1_chromosome_landscape_%s_%s.png",
                                           tag, ctx$contrast_name)),
                key = "s1_chromosome_landscape")
      safe_write_csv(land$positions,
                     file.path(out_dir, sprintf("S1_chromosome_positions_%s_%s.csv",
                                                tag, ctx$contrast_name)))
      status_rows[[j]] <- data.frame(
        lfc_cut = chrom_lfc, mapped = land$mapped, chromosomes = land$chromosomes,
        unmapped_degs = land$unmapped, stringsAsFactors = FALSE)
    }
    do.call(rbind, status_rows)
  }, error = function(e) {
    warning("Chromosome landscape skipped: ", conditionMessage(e), call. = FALSE)
    data.frame(lfc_cut = NA_real_, mapped = NA_integer_, chromosomes = NA_integer_,
               unmapped_degs = NA_integer_)
  })
  status$chromosome_landscape <- chrom_status

  ctx$mylist <- mylist
  ctx$class_summary <- class_summary

  # Ranked identifiers consumed by the enrichment and regulatory sections.
  ctx <- attach_ranked_context(cfg, ctx, tab)
  ctx$s1_status <- status

  write_status(file.path(out_dir, "status.json"), "section01_screening", "completed", NULL,
               list(thresholds = thresholds,
                    deg_counts = setNames(class_summary$total_degs, class_summary$lfc_limit),
                    ranked_genes = length(ctx$gene_list),
                    chromosome_degs = sum(chrom_status$mapped, na.rm = TRUE),
                    chromosome_count = if (any(is.finite(chrom_status$chromosomes)))
                      max(chrom_status$chromosomes, na.rm = TRUE) else NA_integer_,
                    chromosome_unmapped = sum(chrom_status$unmapped_degs, na.rm = TRUE),
                    chromosome_lfc = paste(chrom_status$lfc_cut, collapse = ", "),
                    chromosome_thresholds = nrow(chrom_status),
                    chromosome_landscape = chrom_status))
  ctx
}
