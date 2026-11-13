# Section 9 - gene-set display: GO heat maps and TF-regulon dot plots.
#
# Port of input/AnaDGEsS5.R, reorganised so each gene-set family gets the view
# that carries its message: GO terms are shown as row z-scored heat maps across
# samples, TF regulons as expression-deviation dot plots. The section runs after
# the TF analysis (section 5), so the TF sets use the same expressed-TF filter
# as the regulon, activity and network sections.

run_geneset_heatmaps <- function(cfg, ctx, out_dir) {
  out_dir <- safe_dir(out_dir)
  # The display sections prefer a supplied FPKM/TPM matrix. Without one the
  # DESeq2 normalised counts are used instead, with a note in the status file,
  # so the section still renders rather than leaving an empty folder.
  expr_fallback <- FALSE
  if (is.null(ctx$expr_matrix) && !is.null(ctx$normalized)) {
    ctx$expr_matrix <- as.matrix(ctx$normalized)
    ctx$expr_source <- "DESeq2 normalised counts (no expression matrix supplied)"
    expr_fallback <- TRUE
  }
  if (is.null(ctx$expr_matrix)) {
    msg <- "No expression matrix available (set input.expression_matrix_file)"
    message("    section 9 skipped: ", msg)
    write_status(file.path(out_dir, "status.json"), "section09_geneset_heatmaps", "skipped", msg)
    return(ctx)
  }
  symbol_col <- ctx$symbols
  expr <- ctx$expr_matrix
  # The expression matrix and the annotation table are joined on the gene
  # identifier. The fallback matrix (DESeq2 normalised counts) holds only the
  # pre-filtered genes, so its rows are fewer than the annotation table and a
  # positional rowsum would stop with "incorrect length for 'group'".
  expr_by_symbol <- symbol_rowsum(expr, ctx$annotations,
                                  cfg$input$gene_id_column, symbol_col)
  if (is.null(expr_by_symbol) || !nrow(expr_by_symbol)) {
    msg <- paste0("Expression matrix rows could not be matched to symbols in '",
                  symbol_col %||% "(unset)", "' (matrix: ", nrow(expr),
                  " rows; annotation: ", nrow(ctx$annotations), " rows)")
    message("    section 9 skipped: ", msg)
    write_status(file.path(out_dir, "status.json"), "section09_geneset_heatmaps", "skipped", msg)
    return(ctx)
  }
  if (nrow(expr) != nrow(ctx$annotations)) {
    message("    note: expression matrix has ", nrow(expr), " rows and the annotation ",
            nrow(ctx$annotations), "; matched by gene id (",
            nrow(expr_by_symbol), " symbols)")
  }

  gene_sets <- list()
  if (!is.null(ctx$go_ora) && length(ctx$go_ora)) {
    # The ten GO heat maps are chosen across all three ontologies, ranked by
    # adjusted P, so the panel is not filled with one ontology's cell-cycle
    # subtree. A term is skipped when more than half of the smaller gene set is
    # shared with a term already picked.
    go_all <- do.call(rbind, lapply(names(ctx$go_ora), function(ont) {
      tab <- as.data.frame(ctx$go_ora[[ont]])
      if (is.null(tab) || !nrow(tab)) return(NULL)
      tab$ontology <- ont
      tab
    }))
    if (!is.null(go_all) && nrow(go_all)) {
      go_all <- go_all[order(go_all$p.adjust, go_all$pvalue), , drop = FALSE]
      picked_members <- list()
      picked_rows <- list()
      for (i in seq_len(nrow(go_all))) {
        if (length(picked_rows) >= 10L) break
        members <- unique(strsplit(as.character(go_all$geneID[i]), "/")[[1]])
        members <- members[!is.na(members) & nzchar(members)]
        if (length(members) < 3L) next
        if (length(picked_members)) {
          overlap <- vapply(picked_members, function(prev) {
            length(intersect(members, prev)) / min(length(members), length(prev))
          }, numeric(1))
          if (any(overlap > 0.5)) next
        }
        picked_members[[length(picked_members) + 1L]] <- members
        picked_rows[[length(picked_rows) + 1L]] <- go_all[i, , drop = FALSE]
      }
      for (j in seq_along(picked_rows)) {
        row <- picked_rows[[j]]
        set_name <- paste0("GO", row$ontology, "_",
                           gsub("[^A-Za-z0-9]+", "_", substr(row$Description, 1, 40)))
        if (set_name %in% names(gene_sets)) set_name <- paste0(set_name, "_", j)
        gene_sets[[set_name]] <- picked_members[[j]]
      }
    }
  }
  if (!is.null(ctx$regulons)) {
    reg <- ctx$regulons
    reg <- reg[toupper(reg$target) %in% toupper(ctx$gene_df$SYMBOL), , drop = FALSE]
    tab <- sort(table(reg$tf), decreasing = TRUE)
    tf_pick <- names(tab)
    # Apply the same expression floor as section 5: the larger of the two group
    # mean FPKMs must reach min_fpkm. Without this, a dot plot could show a
    # factor the TF analysis itself excludes as unexpressed.
    if (identical(ctx$expr_unit, "FPKM") && !is.null(ctx$expr_matrix)) {
      expr_sym <- symbol_rowsum(ctx$expr_matrix, ctx$annotations,
                                cfg$input$gene_id_column, ctx$symbols, upper = TRUE)
      tf_fpkm <- if (!is.null(expr_sym) && nrow(expr_sym)) {
        t(expr_sym[intersect(toupper(tf_pick), rownames(expr_sym)), , drop = FALSE])
      }
      if (!is.null(tf_fpkm) && nrow(tf_fpkm) && ncol(tf_fpkm)) {
        groups_here <- as.character(ctx$meta[rownames(tf_fpkm), ctx$group_col])
        group_levels <- intersect(ctx$groups, unique(groups_here))
        if (!length(group_levels)) group_levels <- unique(groups_here)
        group_means <- vapply(group_levels, function(g) {
          colMeans(tf_fpkm[groups_here == g, , drop = FALSE], na.rm = TRUE)
        }, numeric(ncol(tf_fpkm)))
        if (is.null(dim(group_means))) group_means <- matrix(group_means, ncol = 1L)
        max_group <- apply(group_means, 1L, max, na.rm = TRUE)
        names(max_group) <- colnames(tf_fpkm)
        min_fpkm <- as.numeric(cfg$analysis$min_fpkm %||% 0)
        tf_pick <- tf_pick[toupper(tf_pick) %in% names(max_group)[max_group >= min_fpkm]]
      }
    }
    for (tf in utils::head(tf_pick, 10L)) {
      gene_sets[[paste0("TF_", tf)]] <- unique(reg$target[reg$tf == tf])
    }
  }
  if (!length(gene_sets)) {
    msg <- "No gene sets available (run sections 2 and 5 first, or supply a TF database)"
    message("    section 9 skipped: ", msg)
    write_status(file.path(out_dir, "status.json"), "section09_geneset_heatmaps", "skipped", msg)
    return(ctx)
  }
  max_genes <- max_display_features(cfg, 40L)

  set_summary <- list()
  skipped_sets <- character(0)
  for (set_name in names(gene_sets)) {
    members <- unique(toupper(gene_sets[[set_name]]))
    idx <- match(members, toupper(rownames(expr_by_symbol)))
    idx <- idx[!is.na(idx)]
    if (length(idx) < 3L) {
      skipped_sets <- c(skipped_sets, paste0(set_name, " (fewer than three members measured)"))
      next
    }
    hm <- expr_by_symbol[idx, , drop = FALSE]
    keep <- rowSums(hm, na.rm = TRUE) > 0
    hm <- hm[keep, , drop = FALSE]
    if (nrow(hm) < 3L) {
      skipped_sets <- c(skipped_sets, paste0(set_name, " (fewer than three expressed genes)"))
      next
    }
    hm <- hm[order(-apply(hm, 1, stats::var)), , drop = FALSE]
    truncated <- nrow(hm) > max_genes
    hm <- utils::head(hm, max_genes)
    ratio <- hm / rowMeans(hm, na.rm = TRUE)
    z <- t(scale(t(hm)))
    z[!is.finite(z)] <- 0
    ann_col <- data.frame(group = ctx$meta[colnames(hm), ctx$group_col])
    rownames(ann_col) <- colnames(hm)
    safe_name <- substr(set_name, 1, 60)

    # One view per gene-set family: GO terms keep the row z-score heat map, TF
    # regulons keep the dot plot. Drawing both for both families duplicated the
    # same genes without adding a distinct reading.
    is_tf <- startsWith(set_name, "TF_")
    draw_error <- NULL
    tryCatch({
      if (!is_tf) {
        sz_hm <- fig_size("s9_geneset_heatmap")
        record_figure_size("s9_geneset_heatmap", sz_hm, paste0("heatmap_rowz_", safe_name, ".png"))
        # pheatmap reserves the annotation-legend column from the default text
        # width, so at a larger fontsize the key runs off the right edge. Draw
        # the returned gtable in a 0.92-wide viewport instead: the heat map
        # keeps its size and the key spills into the spare margin.
        ph_hm <- pheatmap::pheatmap(z, scale = "none", show_rownames = TRUE, angle_col = 90,
                                    color = heatmap_palette(),
                                    annotation_col = ann_col, fontsize = 14,
                                    fontsize_row = 12, fontsize_col = 14,
                                    border_color = NA, silent = TRUE,
                                    main = paste0(set_name, " (row z-score)"))
        hm_path <- file.path(out_dir, paste0("heatmap_rowz_", safe_name, ".png"))
        grDevices::png(hm_path, width = sz_hm[1], height = sz_hm[2], units = "in",
                       res = 320, bg = "white")
                       plot_device <- grDevices::dev.cur()
                       do.call(on.exit, list(substitute(close_plot_device(device), list(device = plot_device)), add = TRUE))
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(x = grid::unit(0.48, "npc"),
                                          width = grid::unit(0.92, "npc")))
        grid::grid.draw(ph_hm$gtable)
        grid::popViewport()
        close_plot_device(plot_device)
      }

      # Melt-based dot plot: point size = expression, fill = deviation from the
      # gene's own mean across samples (original ph1 panel). TF regulons only.
      if (is_tf &&
          requireNamespace("reshape2", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)) {
        exp_long <- reshape2::melt(log10(hm + 1), value.name = "Exp")
        ratio_long <- reshape2::melt(ratio, value.name = "Diff")
        names(exp_long)[1:2] <- c("GeneSymbol", "Sample")
        names(ratio_long)[1:2] <- c("GeneSymbol", "Sample")
        dot_data <- merge(exp_long, ratio_long, by = c("GeneSymbol", "Sample"))
        dot_data$GeneSymbol <- factor(dot_data$GeneSymbol, levels = rev(rownames(hm)))
        dot_data$Sample <- factor(dot_data$Sample, levels = colnames(hm))
        ph <- ggplot2::ggplot(dot_data, ggplot2::aes(Sample, GeneSymbol)) +
          ggplot2::geom_point(ggplot2::aes(fill = Diff, size = Exp), colour = "black", shape = 21) +
          ggplot2::scale_fill_gradient2(low = "#13a075", mid = "white", high = "#e69f03", midpoint = 1,
                                        name = "Fold to row mean") +
          ggplot2::scale_size_continuous(range = c(0.8, 4.2), breaks = size_breaks(dot_data$Exp),
                                         name = if (expr_fallback) "log10(normalised counts + 1)"
                                                else expression(log[10] * "(FPKM+1)")) +
          ggplot2::scale_x_discrete(position = "top") +
          ggplot2::labs(x = NULL, y = NULL,
                        title = paste(strwrap(gsub("_", " ", set_name), width = 42), collapse = "\n")) +
          ggplot2::theme_minimal(base_family = publication_font(), base_size = 12) +
          ggplot2::theme(axis.text.x.top = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5,
                                                                  face = "bold", size = 10),
                         axis.text.y = ggplot2::element_text(size = 9.5),
                         legend.text = ggplot2::element_text(size = 9),
                         legend.title = ggplot2::element_text(size = 10, face = "bold"),
                         legend.key.size = ggplot2::unit(0.85, "lines"),
                         legend.margin = ggplot2::margin(2, 2, 2, 2),
                         plot.title = ggplot2::element_text(size = 11.5, face = "bold",
                                                             margin = ggplot2::margin(b = 22)),
                         plot.margin = ggplot2::margin(4, 6, 4, 4),
                         panel.grid = ggplot2::element_blank())
        # The dot plot is drawn on its own. The original script paired it with a
        # gene dendrogram, but the rows are already ordered by variance and the
        # tree only made the panel narrower without adding a reading the heat map
        # above it does not already give.
        save_plot(ph, file.path(out_dir, paste0("dotplot_", safe_name, ".png")),
                  key = "s9_geneset_dotplot")
      }
    }, error = function(e) draw_error <<- conditionMessage(e))
    if (!is.null(draw_error)) {
      skipped_sets <- c(skipped_sets, paste0(set_name, " (", draw_error, ")"))
      next
    }
    set_summary[[set_name]] <- data.frame(
      set = set_name,
      family = if (is_tf) "TF regulon" else "GO term",
      figure = if (is_tf) "dotplot" else "heatmap_rowz",
      members_found = nrow(hm), truncated = truncated, stringsAsFactors = FALSE)
  }
  if (length(set_summary)) {
    safe_write_csv(do.call(rbind, set_summary), file.path(out_dir, "geneset_heatmap_summary.csv"))
  }

  ctx$geneset_tables <- set_summary
  if (!length(set_summary)) {
    message("    section 9 produced no figures: ", length(skipped_sets),
            " gene set(s) had fewer than three measured members")
  } else {
    message(sprintf("    section 9 rendered %d of %d gene sets", length(set_summary), length(gene_sets)))
  }
  write_status(file.path(out_dir, "status.json"), "section09_geneset_heatmaps", "completed", NULL,
               list(gene_sets = names(gene_sets), rendered = length(set_summary),
                    skipped = skipped_sets, expression_source = ctx$expr_source,
                    normalized_counts_fallback = expr_fallback))
  ctx
}
