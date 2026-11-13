# Section 5 - transcription-factor regulon analysis.
#
# Port of input/AnaTFs.R with the statistics upgraded. The original script
# intersected the differentially expressed genes with a TF-target database,
# counted each TF's up- and down-regulated targets, and annotated each regulon
# by GO. This section keeps that table and adds the three tests that make a
# regulon claim defensible:
#   * a hypergeometric over-representation test of DEGs inside each regulon,
#     against the tested-gene universe rather than against nothing;
#   * a regulon-level GSEA (fgsea) over the independently pre-screened list,
#     using low_cutoff and real, unbounded fold changes;
#   * a coherence score, the correlation between the TF's own expression and
#     the mean expression of its regulon across samples.

run_tf_regulon <- function(cfg, ctx, out_dir) {
  out_dir <- safe_dir(out_dir)
  if (is.null(ctx$regulons)) {
    write_status(file.path(out_dir, "status.json"), "section05_tf_regulon", "skipped",
                 paste0("No TF-target resource available (DoRothEA ", cfg$input$tf_species %||% "?", ")"))
    return(ctx)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    write_status(file.path(out_dir, "status.json"), "section05_tf_regulon", "skipped", "ggplot2 unavailable")
    return(ctx)
  }
  alpha <- as.numeric(cfg$analysis$alpha %||% 0.05)
  reg <- ctx$regulons
  reg$tf <- as.character(reg$tf)
  reg$target <- as.character(reg$target)
  reg$tf_key <- toupper(reg$tf)
  reg$target_key <- toupper(reg$target)
  reg <- reg[!duplicated(reg[, c("tf_key", "target_key")]), , drop = FALSE]

  gd <- ctx$gene_df
  gd$SYMBOL_KEY <- toupper(as.character(gd$SYMBOL))
  gd <- gd[!duplicated(gd$SYMBOL_KEY), , drop = FALSE]
  universe <- unique(gd$SYMBOL_KEY)
  lfc_cut <- as.numeric(cfg$analysis$lfc_threshold %||% 0.75)
  deg <- gd[is.finite(gd$padj) & gd$padj < alpha & abs(gd$logFC) >= lfc_cut, , drop = FALSE]
  deg_up <- unique(deg$SYMBOL_KEY[deg$logFC > 0])
  deg_down <- unique(deg$SYMBOL_KEY[deg$logFC < 0])

  # Regulons restricted to genes that were actually tested and measured.
  reg_u <- reg[reg$target_key %in% universe & reg$tf_key %in% universe, , drop = FALSE]
  if (!nrow(reg_u)) {
    write_status(file.path(out_dir, "status.json"), "section05_tf_regulon", "failed",
                 "No TF-target pair mapped onto the tested gene universe")
    return(ctx)
  }
  min_targets <- as.integer(cfg$regulatory$min_targets %||% 5L)
  reg_sizes <- table(reg_u$tf_key)
  keep_tfs <- names(reg_sizes)[reg_sizes >= min_targets]
  reg_u <- reg_u[reg_u$tf_key %in% keep_tfs, , drop = FALSE]

  if (!nrow(reg_u)) {
    write_status(file.path(out_dir, "status.json"), "section05_tf_regulon", "skipped",
                 "No regulon meets tf_min_targets in the analysis universe")
    return(ctx)
  }

  # --- 1. overlap and hypergeometric over-representation ---------------------
  N <- length(universe)
  k_deg <- length(unique(c(deg_up, deg_down)))
  edges <- data.frame(tf = reg_u$tf, tf_key = reg_u$tf_key,
                      target = reg_u$target, target_key = reg_u$target_key,
                      confidence = reg_u$confidence, mor = reg_u$mor,
                      target_log2FC = gd$logFC[match(reg_u$target_key, gd$SYMBOL_KEY)],
                      target_is_deg = reg_u$target_key %in% c(deg_up, deg_down),
                      target_direction = ifelse(reg_u$target_key %in% deg_up, "up",
                                                ifelse(reg_u$target_key %in% deg_down, "down", "ns")),
                      stringsAsFactors = FALSE)
  safe_write_csv(edges, file.path(out_dir, paste0("TF_target_edges_", ctx$contrast_name, ".csv")))
  message("    section 5: TF-target edges written (", nrow(edges), " rows); running regulon over-representation")

  tf_split <- split(edges$target_key, edges$tf_key)
  ora_rows <- lapply(names(tf_split), function(tf) {
    targets <- unique(tf_split[[tf]])
    m <- length(targets)
    overlap <- intersect(targets, c(deg_up, deg_down))
    q <- length(overlap)
    p <- if (q == 0L) 1 else stats::phyper(q - 1L, m, N - m, k_deg, lower.tail = FALSE)
    ov_up <- length(intersect(targets, deg_up))
    ov_down <- length(intersect(targets, deg_down))
    lfc <- gd$logFC[match(targets, gd$SYMBOL_KEY)]
    data.frame(tf = edges$tf[match(tf, edges$tf_key)], tf_key = tf,
               targets_tested = m, targets_deg = q,
               targets_up = ov_up, targets_down = ov_down,
               expected_deg = m * k_deg / N,
               odds_ratio = if (q == 0L) 0 else (q * (N - m - k_deg + q)) / ((m - q) * (k_deg - q)),
               p_value = p,
               mean_target_log2FC = mean(lfc, na.rm = TRUE),
               median_target_log2FC = stats::median(lfc, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  ora <- do.call(rbind, ora_rows)
  ora$padj <- stats::p.adjust(ora$p_value, method = "BH")
  message("    section 5: regulon over-representation done (", nrow(ora), " TFs); running regulon GSEA")

  # --- 2. regulon GSEA over the independently pre-screened ranked list -------------------------
  gsea_rows <- NULL
  if (requireNamespace("fgsea", quietly = TRUE)) {
    # fgsea matches pathway members to the names of the ranked vector, so the
    # members keep the original symbol spelling (rank_symbol is keyed that way).
    # The pathway *names* carry the upper-case analysis key, because they become
    # the join key when the GSEA table is merged into the over-representation
    # table; using the symbol spelling there silently matched nothing.
    rank_symbol <- make_unique_ranked_vector(ctx$gsea_gene_df$SYMBOL, ctx$gsea_gene_df$logFC)
    sets <- lapply(split(edges$target, edges$tf_key), function(x) unique(x))
    fg <- tryCatch(fgsea::fgsea(pathways = sets, stats = rank_symbol,
                                minSize = as.integer(cfg$regulatory$gsea_min_size %||% 10L),
                                maxSize = as.integer(cfg$regulatory$gsea_max_size %||% 500L)),
                   error = record_optional_error)
    if (!is.null(fg) && nrow(fg) > 0L) {
      gsea_rows <- as.data.frame(fg)
      gsea_rows$leadingEdge <- vapply(gsea_rows$leadingEdge, paste, character(1), collapse = "/")
      names(gsea_rows)[names(gsea_rows) == "pathway"] <- "tf_key"
      safe_write_csv(gsea_rows, file.path(out_dir, paste0("TF_regulon_GSEA_", ctx$contrast_name, ".csv")))
    }
  }

  # --- 3. regulon coherence: TF expression vs mean regulon expression --------
  coherence <- data.frame(tf_key = character(), coherence = numeric(), tf_log2FC = numeric())
  if (!is.null(ctx$expr_matrix) && requireNamespace("stats", quietly = TRUE)) {
    expr_sym <- symbol_rowsum(ctx$expr_matrix, ctx$annotations,
                              cfg$input$gene_id_column, ctx$symbols, upper = TRUE)
    if (is.null(expr_sym) || !nrow(expr_sym)) {
      expr_sym <- NULL
    } else {
      expr_sym <- log2(expr_sym + 1)
    }
    coh <- lapply(names(tf_split), function(tf) {
      targets <- intersect(unique(tf_split[[tf]]), rownames(expr_sym))
      if (length(targets) < (cfg$regulatory$coherence_min_targets %||% 5L)) return(NULL)
      if (!(tf %in% rownames(expr_sym))) return(NULL)
      regulon_profile <- colMeans(expr_sym[targets, , drop = FALSE])
      ctf <- stats::cor(expr_sym[tf, ], regulon_profile, method = "pearson")
      data.frame(tf_key = tf, coherence = as.numeric(ctf),
                 tf_log2FC = gd$logFC[match(tf, gd$SYMBOL_KEY)],
                 n_regulon_profiled = length(targets))
    })
    coh <- do.call(rbind, coh)
    if (!is.null(coh)) coherence <- coh
  }

  summary_df <- ora
  if (!is.null(gsea_rows)) {
    # Rename the GSEA columns before merging: `padj` exists on both sides, and
    # relying on merge()'s suffixes would collide with the over-representation
    # padj and overwrite it.
    gsea_merge <- gsea_rows[, c("tf_key", "NES", "ES", "size", "pval", "padj")]
    names(gsea_merge) <- c("tf_key", "NES", "ES", "regulon_size", "gsea_pval", "gsea_padj")
    summary_df <- merge(summary_df, gsea_merge, by = "tf_key", all.x = TRUE)
  }
  if (!"gsea_padj" %in% names(summary_df)) summary_df$gsea_padj <- NA_real_
  if (nrow(coherence)) summary_df <- merge(summary_df, coherence, by = "tf_key", all.x = TRUE)
  summary_df <- summary_df[order(-summary_df$targets_deg), , drop = FALSE]
  safe_write_csv(summary_df, file.path(out_dir, paste0("TF_regulon_summary_", ctx$contrast_name, ".csv")))

  # --- 4. per-TF regulon functional annotation (original AnaTFs step) --------
  go_rows <- list()
  if (requireNamespace("clusterProfiler", quietly = TRUE)) {
    orgdb_name <- cfg$species$orgdb
    if (requireNamespace(orgdb_name, quietly = TRUE)) {
      orgdb <- getExportedValue(orgdb_name, orgdb_name)
      # The overlap test runs on case-folded keys, but the annotation package is
      # keyed by the original symbol spelling, so map back before calling it.
      key2sym <- stats::setNames(gd$SYMBOL, gd$SYMBOL_KEY)
      universe_sym <- unique(key2sym[universe])
      universe_sym <- universe_sym[!is.na(universe_sym)]
      annotate_tfs <- utils::head(summary_df$tf_key[order(-summary_df$targets_deg)],
                                  as.integer(cfg$regulatory$network_tfs %||% 25L))
      for (tf in annotate_tfs) {
        targets <- key2sym[unique(tf_split[[tf]])]
        targets <- unique(targets[!is.na(targets)])
        if (length(targets) < 5L) next
        ego <- tryCatch(clusterProfiler::enrichGO(gene = targets, universe = universe_sym, OrgDb = orgdb,
                                                  keyType = ctx$symbol_type, ont = "BP",
                                                  pAdjustMethod = "BH", pvalueCutoff = alpha,
                                                  qvalueCutoff = alpha, readable = FALSE),
                        error = record_optional_error)
        if (is.null(ego)) next
        df <- as.data.frame(ego)
        if (nrow(df) == 0L) next
        go_rows[[tf]] <- data.frame(tf_key = tf,
                                    top_GO_BP = paste(utils::head(df$Description, 5L), collapse = " / "),
                                    n_terms = nrow(df),
                                    stringsAsFactors = FALSE)
      }
    }
  }
  if (length(go_rows)) {
    go_tab <- do.call(rbind, go_rows)
    safe_write_csv(go_tab, file.path(out_dir, paste0("TF_regulon_GO_BP_", ctx$contrast_name, ".csv")))
    summary_df <- merge(summary_df, go_tab, by = "tf_key", all.x = TRUE)
    safe_write_csv(summary_df, file.path(out_dir, paste0("TF_regulon_summary_", ctx$contrast_name, ".csv")))
  }

  # --- 4b. the TFs' own expression ------------------------------------------
  # A regulon claim is only worth drawing when the factor itself is transcribed,
  # so every selected TF carries its own FPKM, and factors below the configured
  # floor are left out of the figures. The full statistics remain in the summary
  # table with `expressed = FALSE`, so nothing is silently deleted.
  tf_fpkm <- NULL                       # samples x TF, FPKM
  min_fpkm <- as.numeric(cfg$analysis$min_fpkm %||% 0)
  if (identical(ctx$expr_unit, "FPKM") && !is.null(ctx$expr_matrix)) {
    expr_sym <- symbol_rowsum(ctx$expr_matrix, ctx$annotations,
                              cfg$input$gene_id_column, ctx$symbols, upper = TRUE)
    if (!is.null(expr_sym) && nrow(expr_sym)) {
      have <- intersect(unique(summary_df$tf_key), rownames(expr_sym))
      if (length(have)) tf_fpkm <- t(expr_sym[have, , drop = FALSE])
    }
  }
  if (!is.null(tf_fpkm)) {
    # The floor is applied to the larger of the two group means, so a factor
    # that is strongly expressed in only one condition is still drawn; the
    # overall mean can hide exactly the on/off behaviour the comparison is
    # looking for. `tf_fpkm` keeps the overall mean for display, and
    # `tf_fpkm_group_max` records the value the filter used.
    groups_here <- as.character(ctx$meta[rownames(tf_fpkm), ctx$group_col])
    group_levels <- intersect(ctx$groups, unique(groups_here))
    if (!length(group_levels)) group_levels <- unique(groups_here)
    group_means <- vapply(group_levels, function(g) {
      colMeans(tf_fpkm[groups_here == g, , drop = FALSE], na.rm = TRUE)
    }, numeric(ncol(tf_fpkm)))
    if (is.null(dim(group_means))) group_means <- matrix(group_means, ncol = 1L)
    rownames(group_means) <- colnames(tf_fpkm)
    max_group_mean <- apply(group_means, 1L, max, na.rm = TRUE)
    names(max_group_mean) <- colnames(tf_fpkm)
    max_group_mean[!is.finite(max_group_mean)] <- NA_real_
    means <- colMeans(tf_fpkm, na.rm = TRUE)
    summary_df$tf_fpkm <- as.numeric(means[summary_df$tf_key])
    summary_df$tf_fpkm_group_max <- as.numeric(max_group_mean[summary_df$tf_key])
    summary_df$expressed <- !is.na(summary_df$tf_fpkm_group_max) &
      summary_df$tf_fpkm_group_max >= min_fpkm
  } else {
    # No expression matrix: the floor cannot be applied, and the status file
    # records that rather than pretending the filter ran.
    summary_df$tf_fpkm <- NA_real_
    summary_df$tf_fpkm_group_max <- NA_real_
    summary_df$expressed <- TRUE
  }
  n_expressed <- sum(summary_df$expressed, na.rm = TRUE)
  safe_write_csv(summary_df, file.path(out_dir, paste0("TF_regulon_summary_", ctx$contrast_name, ".csv")))
  message("    section 5: ", n_expressed, " of ", nrow(summary_df),
          " TFs pass the group-mean FPKM floor (max group mean >= ", min_fpkm, ")")

  # --- 5. figures -----------------------------------------------------------
  # A dot plot instead of the stacked bar: the same up/down target counts, but the
  # TF's own expression rides along as the fill colour and the regulons that are
  # actually significant are outlined, so one panel answers "how many targets, in
  # which direction, and is the TF itself moving".
  n_show <- as.integer(cfg$regulatory$tf_dotplot_max %||% 30L)
  expressed_df <- summary_df[summary_df$expressed %in% TRUE, , drop = FALSE]
  if (!nrow(expressed_df)) {
    ctx$tf_summary <- summary_df; ctx$tf_edges <- edges; ctx$tf_gsea <- gsea_rows
    write_status(file.path(out_dir, "status.json"), "section05_tf_regulon", "completed",
                 "No TF passes the expression floor; statistics written, figures skipped",
                 list(tfs_tested = nrow(summary_df), tfs_expressed = 0L, tfs_drawn = 0L))
    return(ctx)
  }
  top <- utils::head(expressed_df[order(-expressed_df$targets_deg), , drop = FALSE], n_show)
  top$tf_log2FC <- gd$logFC[match(top$tf_key, gd$SYMBOL_KEY)]
  top$significant <- (top$padj < alpha) |
    (!is.na(top$gsea_padj) & top$gsea_padj < alpha)
  # Labels carry the TF name only; the expression values stay in the tables and
  # in the dedicated expression panel instead of being repeated on every axis.
  tf_label <- top$tf
  dot_df <- rbind(
    data.frame(tf = top$tf, label = tf_label, direction = "Up targets", n = top$targets_up,
               tf_log2FC = top$tf_log2FC, significant = top$significant),
    data.frame(tf = top$tf, label = tf_label, direction = "Down targets", n = top$targets_down,
               tf_log2FC = top$tf_log2FC, significant = top$significant))
  dot_df$tf <- factor(dot_df$tf, levels = rev(top$tf))
  dot_df$label <- factor(dot_df$label, levels = rev(tf_label))
  dot_df$direction <- factor(dot_df$direction, levels = c("Up targets", "Down targets"))
  p1 <- ggplot2::ggplot(dot_df, ggplot2::aes(direction, label)) +
    ggplot2::geom_point(ggplot2::aes(size = n, fill = tf_log2FC, colour = significant),
                        shape = 21, stroke = 0.7) +
    ggplot2::scale_size_continuous(range = c(1.6, 10), breaks = size_breaks(dot_df$n),
                                   name = "DE targets") +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                                  na.value = "grey85", limits = c(-1, 1) * display_lfc_limit(cfg, dot_df$tf_log2FC),
                                  oob = scales::squish, name = expression("TF log"[2] * "FC")) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#333333", `FALSE` = "grey80"),
                                 name = paste0("Regulon q < ", alpha)) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0("TF targets and own expression: ", ctx$contrast_name)) +
    publication_theme(base_size = 12) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 11.5, face = "bold"),
                   axis.text.x = ggplot2::element_text(size = 13, face = "bold"),
                   panel.grid.major.x = ggplot2::element_line(colour = "grey92", linetype = "dotted"))
  save_plot(p1, file.path(out_dir, paste0("TF_dotplot_", ctx$contrast_name, ".png")),
            key = "s5_tf_dotplot", n = nrow(top))

  # Dedicated expression panel: one point per TF and sample, so the abundance the
  # floor was applied to can be read directly instead of being folded into a
  # single mean.
  if (!is.null(tf_fpkm) && nrow(tf_fpkm) > 0L && ncol(tf_fpkm) > 0L) {
    # Highest mean FPKM first, and only the factors the target dot plot shows.
    fpkm_rank <- sort(colMeans(tf_fpkm, na.rm = TRUE), decreasing = TRUE)
    tf_ord <- intersect(names(fpkm_rank), utils::head(top$tf_key, n_show))
    if (length(tf_ord)) {
      mat <- tf_fpkm[, tf_ord, drop = FALSE]
      # Display the same title-case symbols the target dot plot uses; the
      # DoRothEA keys behind tf_ord are upper case.
      tf_show <- top$tf[match(tf_ord, top$tf_key)]
      tf_show[is.na(tf_show)] <- tf_ord[is.na(tf_show)]
      # as.numeric() unwinds a matrix column by column, so each factor's samples
      # come out together: the factor labels repeat (each = samples per factor)
      # while the sample names cycle (times = number of factors). Repeating the
      # whole factor vector instead interleaves the two and lays the points out
      # on a diagonal, one per row.
      tf_expr <- data.frame(
        tf = rep(tf_show, each = nrow(mat)),
        Sample = rep(rownames(mat), times = ncol(mat)),
        FPKM = as.numeric(mat),
        stringsAsFactors = FALSE)
      tf_expr$log_fpkm <- log10(tf_expr$FPKM + 1)
      # The TF's own log2 fold change is constant across its samples; it colours
      # the panel, while the dot area carries the expression level.
      tf_expr$tf_lfc <- rep(top$tf_log2FC[match(tf_ord, top$tf_key)], each = nrow(mat))
      tf_expr$group <- as.character(ctx$meta[tf_expr$Sample, ctx$group_col])
      # The samples keep the metadata order so the two conditions sit side by
      # side under the facet strip.
      tf_expr$tf <- factor(tf_expr$tf, levels = rev(tf_show))
      tf_expr$Sample <- factor(tf_expr$Sample, levels = rownames(mat))
      tf_expr$group <- factor(tf_expr$group, levels = ctx$groups)
      p3 <- ggplot2::ggplot(tf_expr, ggplot2::aes(Sample, tf)) +
        ggplot2::geom_point(ggplot2::aes(size = log_fpkm, fill = tf_lfc),
                            shape = 21, colour = "grey25", stroke = 0.3) +
        ggplot2::scale_size_continuous(range = c(0.8, 8),
                                       breaks = size_breaks(tf_expr$log_fpkm),
                                       labels = function(x) sprintf("%.2f", x),
                                       name = expression("Expression  log"[10] * "(FPKM+1)")) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                      midpoint = 0, limits = c(-1, 1) * display_lfc_limit(cfg, tf_expr$tf_lfc), oob = scales::squish,
                                      name = expression("TF log"[2] * "FC")) +
        ggplot2::facet_grid(. ~ group, scales = "free_x", space = "free_x") +
        ggplot2::labs(x = NULL, y = NULL,
                      title = "Transcription-factor expression",
                      subtitle = paste0("Dot size = log10(FPKM + 1); colour = TF log2 fold change\n",
                                        "Max group-mean FPKM >= ", min_fpkm, " | ", ctx$contrast_name)) +
        publication_theme() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10.5, face = "bold"),
                       axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10.5),
                       strip.text = ggplot2::element_text(face = "bold", size = 11),
                       plot.title = ggplot2::element_text(face = "bold", size = 13),
                       plot.subtitle = ggplot2::element_text(size = 10.5, colour = "grey25"),
                       panel.grid.major.x = ggplot2::element_line(colour = "grey92", linetype = "dotted"))
      save_plot(p3, file.path(out_dir, paste0("TF_expression_dotplot_", ctx$contrast_name, ".png")),
                key = "s5_tf_expression", n = length(tf_ord))
    }
  }

  if (!is.null(summary_df$NES) && any(is.finite(summary_df$NES))) {
    vd <- summary_df[is.finite(summary_df$NES), , drop = FALSE]
    vd$significant <- !is.na(vd$gsea_padj) & vd$gsea_padj < alpha
    vd$mark <- ifelse(vd$expressed %in% TRUE, "expressed", "below FPKM floor")
    p2 <- ggplot2::ggplot(vd, ggplot2::aes(NES, -log10(pmax(p_value, .Machine$double.xmin)))) +
      ggplot2::geom_point(ggplot2::aes(colour = significant, size = targets_tested, shape = mark),
                          alpha = 0.75) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey65"), name = paste0("GSEA padj < ", alpha)) +
      ggplot2::scale_shape_manual(values = c(expressed = 16, `below FPKM floor` = 4),
                                  name = paste0("Max group mean FPKM >= ", min_fpkm)) +
      ggplot2::scale_size_continuous(range = c(1, 6), breaks = size_breaks(vd$targets_tested),
                                     name = "Regulon size") +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
      ggplot2::labs(x = "Regulon NES (GSEA)", y = expression("-log"[10] * " OR-test P"),
                    title = paste0("Regulon enrichment vs regulon activity: ", ctx$contrast_name)) +
      publication_theme()
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      lab <- vd[vd$significant & !is.na(vd$gsea_padj), , drop = FALSE]
      lab <- utils::head(lab[order(-abs(lab$NES)), , drop = FALSE], 15L)
      lab$label <- lab$tf
      p2 <- p2 + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = label), size = 3.2,
                                          max.overlaps = 25, na.rm = TRUE)
    }
    save_plot(p2, file.path(out_dir, paste0("TF_regulon_volcano_", ctx$contrast_name, ".png")),
              key = "s5_regulon_volcano")
  }

  # --- 6. chromosome map of the top regulons ---------------------------------
  # The same stem landscape as the DEG map, but for the genes each of the ten
  # strongest regulons controls. Genes are the DE targets of that TF at the
  # chromosome-map threshold, and the stem height is the target's own |log2FC|,
  # so a TF whose targets cluster on one chromosome is visible at a glance.
  chrom_tf_status <- tryCatch({
    orgdb_name <- cfg$species$orgdb
    if (is.null(orgdb_name) || !requireNamespace(orgdb_name, quietly = TRUE)) {
      stop("OrgDb ", orgdb_name %||% "<none>", " is not installed")
    }
    orgdb <- getExportedValue(orgdb_name, orgdb_name)
    n_tf_chrom <- as.integer(cfg$regulatory$tf_chromosome_n %||% 10L)
    tf_sel <- utils::head(expressed_df$tf_key[order(-expressed_df$targets_deg)], n_tf_chrom)
    key2sym <- stats::setNames(gd$SYMBOL, gd$SYMBOL_KEY)
    cap <- display_lfc_limit(cfg, edges$target_log2FC)
    map_lfc <- lfc_cut
    pos_rows <- list()
    tf_summary_rows <- list()
    for (tf in tf_sel) {
      sub_e <- edges[edges$tf_key == tf & edges$target_key %in% c(deg_up, deg_down), , drop = FALSE]
      sub_e <- sub_e[!duplicated(sub_e$target_key), , drop = FALSE]
      if (!nrow(sub_e)) next
      feats <- data.frame(gene = unname(key2sym[sub_e$target_key]),
                          lfc = as.numeric(sub_e$target_log2FC), stringsAsFactors = FALSE)
      feats <- feats[!is.na(feats$gene), , drop = FALSE]
      tf_name <- summary_df$tf[match(tf, summary_df$tf_key)]
      if (is.na(tf_name) || !nzchar(tf_name)) tf_name <- tf
      land <- chromosome_landscape(
        feats, orgdb, cfg$species$common_name, cap = cap,
        title = paste0("Chromosomal distribution of ", tf_name, " regulon targets: ",
                       ctx$contrast_name),
        subtitle = paste0("DE targets at |log2FC| >= ", map_lfc,
                          "; stem height is proportional to the target's |log2FC| (capped at ",
                          cap, "); top ", n_tf_chrom, " regulons by DE target count"),
        legend_title = "Target direction")
      if (is.null(land)) next
      save_plot(land$plot,
                file.path(out_dir, sprintf("TF_chromosome_%s_%s.png", tf_name, ctx$contrast_name)),
                key = "s5_tf_chromosome")
      pos_rows[[tf]] <- cbind(tf = tf_name, land$positions, stringsAsFactors = FALSE)
      tf_summary_rows[[tf]] <- data.frame(
        tf_key = tf, tf = tf_name, targets_mapped = land$mapped,
        chromosomes = land$chromosomes, stringsAsFactors = FALSE)
    }
    if (length(pos_rows)) {
      pos_all <- do.call(rbind, pos_rows)
      rownames(pos_all) <- NULL
      safe_write_csv(pos_all,
                     file.path(out_dir, paste0("TF_regulon_chromosome_positions_", ctx$contrast_name, ".csv")))
      tf_chrom <- do.call(rbind, tf_summary_rows)

      # Compact companion: how many targets each regulon places on each
      # chromosome, so the ten large maps can be compared in one panel.
      chr_levels <- intersect(canonical_chromosomes(cfg$species$common_name),
                              unique(pos_all$chr))
      mat <- table(factor(pos_all$tf, levels = rev(tf_chrom$tf)),
                   factor(pos_all$chr, levels = chr_levels))
      hm <- as.data.frame(mat, stringsAsFactors = FALSE)
      names(hm) <- c("tf", "chr", "n")
      hm$tf <- factor(hm$tf, levels = rev(tf_chrom$tf))
      hm$chr <- factor(hm$chr, levels = chr_levels)
      tf_order <- tf_chrom$tf
      p_chr <- ggplot2::ggplot(hm, ggplot2::aes(chr, tf, fill = n)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
        ggplot2::geom_text(ggplot2::aes(label = ifelse(n > 0, n, "")), size = 4.5,
                           colour = "grey15") +
        ggplot2::scale_fill_gradient(low = "#F7F7F7", high = "#B2182B",
                                     name = "DE targets") +
        ggplot2::labs(x = "Chromosome", y = NULL,
                      title = paste0("Chromosomal distribution of top regulon targets: ",
                                     ctx$contrast_name),
                      subtitle = paste0("Top ", length(tf_order),
                                        " regulons by DE target count; cell = number of targets on that chromosome")) +
        publication_theme() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10.8),
                       axis.text.y = ggplot2::element_text(size = 10.5, face = "bold"),
                       panel.grid = ggplot2::element_blank())
      save_plot(p_chr, file.path(out_dir, paste0("TF_chromosome_summary_", ctx$contrast_name, ".png")),
                key = "s5_tf_chromosome_summary", n = length(tf_order))
      summary_df <- merge(summary_df, tf_chrom[, c("tf_key", "targets_mapped", "chromosomes")],
                          by = "tf_key", all.x = TRUE)
      safe_write_csv(summary_df, file.path(out_dir, paste0("TF_regulon_summary_", ctx$contrast_name, ".csv")))
      list(tfs = length(tf_order), mapped = sum(tf_chrom$targets_mapped))
    } else {
      list(tfs = 0L, mapped = 0L)
    }
  }, error = function(e) {
    warning("TF chromosome maps skipped: ", conditionMessage(e), call. = FALSE)
    list(tfs = NA_integer_, mapped = NA_integer_)
  })

  ctx$tf_summary <- summary_df
  ctx$tf_edges <- edges
  ctx$tf_gsea <- gsea_rows
  write_status(file.path(out_dir, "status.json"), "section05_tf_regulon", "completed", NULL,
               list(tfs_tested = nrow(summary_df),
                    tfs_ora_significant = sum(summary_df$padj < alpha, na.rm = TRUE),
                    tfs_gsea_significant = if (is.null(gsea_rows)) NA_integer_ else sum(gsea_rows$padj < alpha, na.rm = TRUE),
                    lfc_cut = lfc_cut, universe = N, confidence_tiers = attr(ctx$regulons, "tiers_kept"),
                    min_fpkm = if (is.null(tf_fpkm)) NA_real_ else min_fpkm,
                    tfs_expressed = n_expressed,
                    tfs_drawn = nrow(top),
                    chromosome_tfs = chrom_tf_status$tfs %||% NA_integer_,
                    chromosome_targets = chrom_tf_status$mapped %||% NA_integer_,
                    chromosome_lfc = as.numeric(cfg$analysis$chromosome_lfc %||% 0.75)))
  ctx
}
