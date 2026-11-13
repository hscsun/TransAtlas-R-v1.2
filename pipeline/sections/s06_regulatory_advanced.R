# Section 6 - advanced regulatory analyses.
#
# Two modules that go beyond counting TF targets:
#   A. regulon activity inference (decoupleR: ULM, AUCell, WMean). Target
#      enrichment says a regulon is over-represented; activity inference
#      estimates whether the regulon is switched on in each sample, and the
#      two answers are not the same.
#   B. pathway activity per sample (GSVA), which turns the gene-set analysis
#      into a per-sample score that can be tested between groups.

run_regulatory_activity <- function(cfg, ctx, out_dir) {
  out_dir <- safe_dir(out_dir)
  sub_status <- list()
  activity_alpha <- as.numeric(cfg$regulatory$activity_alpha %||% cfg$analysis$alpha %||% 0.05)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    write_status(file.path(out_dir, "status.json"), "section06_regulatory_activity", "skipped", "ggplot2 unavailable")
    return(ctx)
  }

  # Symbol-indexed expression matrices used by all three modules. The VST rows
  # are filtered gene identifiers, so symbols are looked up by gene ID rather
  # than taken from the annotation table by position.
  gene_id_col <- cfg$input$gene_id_column
  vst_sym <- symbol_rowsum(ctx$vst, ctx$annotations, gene_id_col, ctx$symbols)
  expr_sym <- NULL
  if (!is.null(ctx$expr_matrix)) {
    expr_sym <- symbol_rowsum(ctx$expr_matrix, ctx$annotations, gene_id_col, ctx$symbols)
    if (!is.null(expr_sym)) expr_sym <- log2(expr_sym + 1)
  }
  if (is.null(vst_sym) || !nrow(vst_sym)) {
    msg <- paste0("VST matrix rows could not be matched to symbols in '",
                  ctx$symbols %||% "(unset)", "' of the annotation table")
    message("    section 6 skipped: ", msg)
    write_status(file.path(out_dir, "status.json"), "section06_regulatory_activity", "skipped", msg)
    return(ctx)
  }

  diff_activity <- function(activity_long, group_of, groups, label) {
    result <- differential_activity(activity_long, group_of, groups)
    if (!is.null(result)) {
      sub_status[[paste0(label, "_profiles")]] <<- nrow(result)
      sub_status[[paste0(label, "_tested")]] <<- sum(result$test_status == "tested")
      sub_status[[paste0(label, "_not_testable")]] <<- sum(result$test_status != "tested")
    }
    result
  }

  # ===================== A. TF regulon activity (decoupleR) =================
  if (requireNamespace("decoupleR", quietly = TRUE) && !is.null(ctx$regulons)) {
    if (!"package:decoupleR" %in% search()) {
      suppressPackageStartupMessages(library(decoupleR, quietly = TRUE))
    }
    reg <- ctx$regulons
    network <- data.frame(source = as.character(reg$tf), target = as.character(reg$target),
                          mor = ifelse(is.na(reg$mor), 1, sign(reg$mor)))
    network <- network[network$source %in% rownames(vst_sym) & network$target %in% rownames(vst_sym), , drop = FALSE]
    network <- network[!duplicated(network[, c("source", "target")]), , drop = FALSE]
    if (nrow(network) > 0L) {
      methods <- unlist(cfg$regulatory$activity_methods %||% c("ulm", "aucell", "wmean"))
      min_size <- as.integer(cfg$regulatory$activity_min_size %||% 5L)
      set.seed(as.integer(cfg$project$seed %||% 1L))
      group_of <- stats::setNames(as.character(ctx$meta[[ctx$group_col]]), rownames(ctx$meta))
      for (mth in methods) {
        acts <- tryCatch(switch(mth,
          ulm = decoupleR::run_ulm(mat = vst_sym, network = network, .source = "source",
                                   .target = "target", .mor = "mor", minsize = min_size),
          wmean = decoupleR::run_wmean(mat = vst_sym, network = network, .source = "source",
                                       .target = "target", .mor = "mor", minsize = min_size),
          mlm = decoupleR::run_mlm(mat = vst_sym, network = network, .source = "source",
                                   .target = "target", .mor = "mor", minsize = min_size),
          aucell = decoupleR::run_aucell(mat = vst_sym, network = network, .source = "source",
                                         .target = "target", minsize = min_size, nproc = 1),
          NULL), error = record_optional_error)
        if (is.null(acts) || nrow(acts) == 0L) next
        acts <- as.data.frame(acts)
        score_col <- if ("score" %in% names(acts)) "score" else names(acts)[ncol(acts)]
        long <- data.frame(source = acts$source, condition = acts$condition,
                           score = acts[[score_col]], stringsAsFactors = FALSE)
        # Build the wide matrix by explicit lookup: decoupleR can emit several
        # rows per source/sample pair, which reshape() would silently truncate.
        long <- long[!duplicated(long[, c("source", "condition")]), , drop = FALSE]
        wide <- data.frame(source = unique(long$source), stringsAsFactors = FALSE)
        key <- paste(long$source, long$condition, sep = "\r")
        for (s in colnames(vst_sym)) {
          wide[[s]] <- long$score[match(paste(wide$source, s, sep = "\r"), key)]
        }
        wide <- wide[, c("source", colnames(vst_sym)), drop = FALSE]
        safe_write_csv(wide, file.path(out_dir, paste0("TF_activity_", mth, "_per_sample.csv")))
        safe_write_csv(long, file.path(out_dir, paste0("TF_activity_", mth, "_long.csv")))

        diff <- diff_activity(long, group_of, ctx$groups, mth)
        if (!is.null(diff)) {
          safe_write_csv(diff, file.path(out_dir, paste0("TF_activity_diff_", mth, ".csv")))
          alpha <- activity_alpha
          dd <- diff[is.finite(diff$p_value), , drop = FALSE]
          if (!nrow(dd)) next
          dd$significant <- !is.na(dd$padj) & dd$padj < alpha
          p <- ggplot2::ggplot(dd, ggplot2::aes(difference, -log10(pmax(p_value, .Machine$double.xmin)))) +
            ggplot2::geom_point(ggplot2::aes(colour = significant, size = abs(mean_treatment) + abs(mean_reference)),
                                alpha = 0.8) +
            ggplot2::scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey65"),
                                         name = paste0("padj < ", alpha)) +
            # The dot size used to be drawn without a legend, so the third
            # quantity in the panel could not be read off the figure.
            ggplot2::scale_size_continuous(
              range = c(1, 6),
              breaks = size_breaks(abs(dd$mean_treatment) + abs(dd$mean_reference)),
              labels = function(x) sprintf("%.1f", x),
              name = "Activity magnitude") +
            ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
            ggplot2::labs(x = paste0("Differential regulon activity (", ctx$groups[2L], " - ", ctx$groups[1L], ")"),
                          y = expression("-log"[10] * " P"),
                          title = paste0("TF activity (", toupper(mth), "): ", ctx$contrast_name),
                          subtitle = paste0("Dot size = |activity| in ", ctx$groups[2L],
                                            " + |activity| in ", ctx$groups[1L],
                                            "; colour = differential-activity P")) +
            publication_theme() +
            ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9.5, colour = "grey25"))
          if (requireNamespace("ggrepel", quietly = TRUE) && any(dd$significant)) {
            lab <- utils::head(dd[dd$significant, , drop = FALSE], 15L)
            lab$label <- lab$source
            p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = label), size = 2.6,
                                              max.overlaps = 25, na.rm = TRUE)
          }
          save_plot(p, file.path(out_dir, paste0("TF_activity_volcano_", mth, ".png")),
                    key = "s6_activity_volcano")

          top_tfs <- utils::head(dd$source[order(dd$p_value)], 30L)
          hm <- as.matrix(wide[match(top_tfs, wide$source), colnames(vst_sym), drop = FALSE])
          rownames(hm) <- top_tfs
          hm <- t(scale(t(hm)))
          hm[!is.finite(hm)] <- 0
          ann_col <- data.frame(group = ctx$meta[colnames(hm), ctx$group_col])
          rownames(ann_col) <- colnames(hm)
          sz_hm <- fig_size("s6_activity_heatmap")
          record_figure_size("s6_activity_heatmap", sz_hm, paste0("TF_activity_heatmap_", mth, ".png"))
          pheatmap::pheatmap(hm, cluster_cols = FALSE, cluster_rows = TRUE, annotation_col = ann_col,
                             color = heatmap_palette(),
                             filename = file.path(out_dir, paste0("TF_activity_heatmap_", mth, ".png")),
                             width = sz_hm[1], height = sz_hm[2], border_color = NA, fontsize = 10,
                             fontsize_row = 9.5,
                             main = paste0("Top differential TF activity (", toupper(mth), ")"))
          sub_status[[paste0("activity_", mth)]] <- sum(diff$padj < alpha, na.rm = TRUE)
        }
      }
      ctx$tf_activity_network <- network
    }
  }

  # ===================== B. Pathway activity (GSVA) ========================
  if (isTRUE(cfg$pathway_activity$run) && requireNamespace("GSVA", quietly = TRUE) &&
      requireNamespace("msigdbr", quietly = TRUE) && !is.null(expr_sym)) {
    collection <- cfg$pathway_activity$collection %||% "H"
    sets_db <- tryCatch(resolve_msigdb_gene_sets(cfg$species$name, collection, NULL,
      cfg$analysis$msigdb_db_species,
      file.path(out_dir, paste0("GSVA_", collection, "_source.json"))),
      error = record_optional_error)
    if (!is.null(sets_db) && nrow(sets_db) > 0L) {
      gene_sets <- split(sets_db$gene_symbol, sets_db$gs_name)
      ok <- tryCatch({
        param <- GSVA::gsvaParam(exprData = as.matrix(expr_sym), geneSets = gene_sets,
                                 kcdf = "Gaussian", minSize = as.integer(cfg$pathway_activity$min_size %||% 10L),
                                 maxSize = as.integer(cfg$pathway_activity$max_size %||% 500L))
        es <- GSVA::gsva(param, verbose = FALSE)
        safe_write_csv(data.frame(pathway = rownames(es), es, check.names = FALSE),
                       file.path(out_dir, paste0("GSVA_", collection, "_scores.csv")))
        group_of <- stats::setNames(as.character(ctx$meta[[ctx$group_col]]), rownames(ctx$meta))
        long <- data.frame(pathway = rep(rownames(es), ncol(es)),
                           condition = rep(colnames(es), each = nrow(es)),
                           score = as.vector(es), stringsAsFactors = FALSE)
        names(long)[1] <- "source"
        diff <- diff_activity(long, group_of, ctx$groups, "gsva")
        if (!is.null(diff)) {
          safe_write_csv(diff, file.path(out_dir, paste0("GSVA_", collection, "_diff.csv")))
          dd <- diff[is.finite(diff$p_value), , drop = FALSE]
          dd <- utils::head(dd[order(dd$p_value), , drop = FALSE], 25L)
          if (nrow(dd)) {
          dd$pathway <- factor(dd$source, levels = rev(dd$source))
          p <- ggplot2::ggplot(dd, ggplot2::aes(difference, pathway, fill = padj < activity_alpha)) +
            ggplot2::geom_col(show.legend = FALSE) +
            ggplot2::scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey70")) +
            ggplot2::labs(x = paste0("GSVA score difference (", ctx$groups[2L], " - ", ctx$groups[1L], ")"),
                          y = NULL, title = paste0("Pathway activity (GSVA, ", collection, "): ", ctx$contrast_name)) +
            publication_theme() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6.5))
          save_plot(p, file.path(out_dir, paste0("GSVA_", collection, "_diff_barplot.png")),
                    key = "s6_gsva_barplot")
          }
          sub_status$gsva_significant <- sum(diff$padj < activity_alpha, na.rm = TRUE)
        }
        TRUE
      }, error = function(e) conditionMessage(e))
      if (!isTRUE(ok)) writeLines(as.character(ok), file.path(out_dir, "GSVA_error.txt"))
    }
  }

  details <- c(list(expression_source = ctx$expr_source), sub_status)
  write_status(file.path(out_dir, "status.json"), "section06_regulatory_activity", "completed", NULL, details)
  ctx$advanced_status <- sub_status
  ctx
}
