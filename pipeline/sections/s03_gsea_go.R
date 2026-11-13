# Section 3 - ranked GSEA over Gene Ontology.
#
# Port of input/AnaDGEsS2.R and input/AnaDGEsS2Mm.R. The original script built
# a ranked Entrez gene list from the threshold-0.25 screening table, ran gseGO,
# wrote one table plus a PDF with one enrichment plot per gene set, and drew a
# multi-set composite. aPEAR's enrichmentNetwork is replaced by enrichplot's
# enrichment map (aPEAR is GitHub-only and not installed here); every other
# output name and layout follows the original.

run_gsea_go <- function(cfg, ctx, out_dir) {
  out_dir <- safe_dir(out_dir)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    write_status(file.path(out_dir, "status.json"), "section03_gsea_go", "skipped", "clusterProfiler unavailable")
    return(ctx)
  }
  orgdb_name <- cfg$species$orgdb
  if (!requireNamespace(orgdb_name, quietly = TRUE)) {
    write_status(file.path(out_dir, "status.json"), "section03_gsea_go", "skipped", paste("Missing", orgdb_name))
    return(ctx)
  }
  orgdb <- getExportedValue(orgdb_name, orgdb_name)
  rank <- ctx$gene_list
  if (is.null(rank) || length(rank) < 100L) {
    write_status(file.path(out_dir, "status.json"), "section03_gsea_go", "skipped", "Ranked gene list unavailable")
    return(ctx)
  }
  min_size <- as.numeric(cfg$analysis$gsea_min_size %||% 15)
  max_size <- as.numeric(cfg$analysis$gsea_max_size %||% 500)
  pcut <- as.numeric(cfg$analysis$gsea_pvalue_cutoff %||% 0.05)
  top_sets <- as.integer(cfg$analysis$gsea_top_sets %||% 9L)
  max_pages <- as.integer(cfg$analysis$gsea_all_sets_max_pages %||% 60L)
  # GSEA reports core_enrichment as the organism's key type (Entrez here); the
  # published tables carry symbols instead, like every other gene list.
  sym_map <- key_to_symbol_map(ctx$gene_map, ctx$key_type, ctx$symbol_type)

  results <- list()
  summary_rows <- list()
  for (ont in c("BP", "MF", "CC")) {
    # fgseaMultilevel estimates small P values by Monte Carlo sampling. Resetting
    # the seed before each ontology makes the significant-set list independent
    # of how much randomness earlier sections consumed, so a full run and a
    # standalone `--sections 3` run produce the same composite panels.
    set.seed(as.integer(cfg$project$seed %||% 1L))
    gse <- tryCatch(clusterProfiler::gseGO(geneList = rank, OrgDb = orgdb,
                                           keyType = ctx$key_type, ont = ont,
                                           minGSSize = min_size, maxGSSize = max_size,
                                           pvalueCutoff = pcut, pAdjustMethod = "BH",
                                           verbose = FALSE, eps = 0),
                    error = record_optional_error)
    if (is.null(gse)) next
    df <- as.data.frame(gse)
    if (nrow(df) == 0L) next
    if ("core_enrichment" %in% names(df)) {
      df$core_enrichment <- ids_to_symbols(df$core_enrichment, sym_map)
    }
    safe_write_csv(df, file.path(out_dir, paste0("GSEA_GO_", ont, "_", ctx$contrast_name, ".csv")))
    results[[ont]] <- gse
    summary_rows[[ont]] <- data.frame(ontology = ont, sets_reported = nrow(df),
                                      sets_significant = sum(df$p.adjust < pcut),
                                      up = sum(df$NES > 0 & df$p.adjust < pcut),
                                      down = sum(df$NES < 0 & df$p.adjust < pcut))
  }
  if (!length(results)) {
    write_status(file.path(out_dir, "status.json"), "section03_gsea_go", "failed", "No GSEA result produced")
    return(ctx)
  }
  safe_write_csv(do.call(rbind, summary_rows), file.path(out_dir, "GSEA_GO_summary.csv"))

  # Per-gene-set enrichment plot, one page per set (original GSEAall.pdf style).
  for (ont in names(results)) {
    gse <- results[[ont]]
    df <- as.data.frame(gse)
    sel <- order(df$p.adjust)[seq_len(min(max_pages, nrow(df)))]
    pdf_path <- file.path(out_dir, sprintf("GSEA_GO_%s_allsets_%s.pdf", ont, ctx$contrast_name))
    sz <- fig_size("s3_go_allsets_pdf")
    record_figure_size("s3_go_allsets_pdf", sz, pdf_path)
    grDevices::cairo_pdf(pdf_path, width = sz[1], height = sz[2], bg = "white")
    plot_device <- grDevices::dev.cur()
    do.call(on.exit, list(substitute(close_plot_device(device), list(device = plot_device)), add = TRUE))
    for (i in sel) {
      title <- sprintf("%s | p.adjust: %.3g", df$Description[i], df$p.adjust[i])
      p <- tryCatch(enrichplot::gseaplot2(gse, geneSetID = i, title = title, base_size = 15),
                    error = record_optional_error)
      if (!is.null(p)) print(p)
    }
    close_plot_device(plot_device)
  }

  # Composite panel of the strongest sets per ontology, half of them induced and
  # half repressed. Taking the N smallest adjusted P values alone produced an
  # all-up panel on this contrast, which hides the down-regulated programme
  # entirely; the direction is also written into the legend label.
  palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999")
  for (ont in names(results)) {
    gse <- results[[ont]]
    df <- as.data.frame(gse)
    # Full gene lists (not just core_enrichment) are needed to reject panels
    # whose pathways overlap heavily. geneInCategory keeps the set IDs as names.
    gene_sets <- tryCatch(clusterProfiler::geneInCategory(gse), error = record_optional_error)
    sel <- diverse_gsea_selection(df, gene_sets, top_sets,
                                  min_size = 50L, max_overlap = 0.4)
    if (!length(sel)) next
    # Record why each set entered the panel: size, NES, and the largest gene
    # overlap with a set already selected (must stay at or below 0.4).
    sel_rows <- df[match(sel, df$ID), c("ID", "Description", "setSize", "NES", "p.adjust"), drop = FALSE]
    sel_rows$direction <- ifelse(sel_rows$NES > 0, "Up", "Down")
    sel_rows$max_overlap_previous <- vapply(seq_along(sel), function(i) {
      if (i == 1L || is.null(gene_sets)) return(0)
      prev <- sel[seq_len(i - 1L)]
      a <- gene_sets[[sel[i]]]
      if (is.null(a) || !length(a)) return(0)
      max(vapply(prev, function(s) {
        b <- gene_sets[[s]]
        if (is.null(b) || !length(b)) return(0)
        length(intersect(a, b)) / min(length(a), length(b))
      }, numeric(1)))
    }, numeric(1))
    safe_write_csv(sel_rows, file.path(out_dir,
                  sprintf("GSEA_GO_%s_composite_sets_%s.csv", ont, ctx$contrast_name)))
    # A copy carries the ▲/▼ legend labels, so the map and ridge plot drawn
    # later from `gse` keep the plain term names.
    gse_lab <- gse
    lab <- gse_lab@result$Description
    lab <- paste0(ifelse(gse_lab@result$NES > 0, "▲ ", "▼ "), lab)
    # The composite canvas is 9 in wide, so long GO names are wrapped onto a
    # second line and the legend is anchored to the top-right corner; the
    # gseaplot2 default (0.8, 0.8) left the longest label running off the edge.
    lab <- wrap_label(lab, width = 30L)
    gse_lab@result$Description <- lab
    p <- tryCatch(enrichplot::gseaplot2(gse_lab, geneSetID = sel, pvalue_table = FALSE, base_size = 16),
                  error = record_optional_error)
    if (is.null(p)) next
    cols <- rep(palette, length.out = length(sel))
    p[[1]] <- p[[1]] + ggplot2::scale_colour_manual(values = cols) +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 13),
                     legend.position = "inside",
                     legend.position.inside = c(0.99, 0.99),
                     legend.justification.inside = c(1, 1),
                     # The legend sits inside the running-score panel; a light
                     # background keeps the curves from crossing the long term
                     # names.
                     legend.background = ggplot2::element_rect(
                       fill = grDevices::adjustcolor("white", alpha.f = 0.88), colour = NA),
                     legend.key = ggplot2::element_rect(
                       fill = grDevices::adjustcolor("white", alpha.f = 0.88), colour = NA))
    p[[2]] <- p[[2]] + ggplot2::scale_colour_manual(values = cols)
    sz <- fig_size("s3_go_composite")
    record_figure_size("s3_go_composite", sz, sprintf("GSEA_GO_%s_composite_%s.png", ont, ctx$contrast_name))
    grDevices::png(file.path(out_dir, sprintf("GSEA_GO_%s_composite_%s.png", ont, ctx$contrast_name)),
                   width = sz[1] * 300, height = sz[2] * 300, res = 300)
                   plot_device <- grDevices::dev.cur()
                   do.call(on.exit, list(substitute(close_plot_device(device), list(device = plot_device)), add = TRUE))
    print(p)
    close_plot_device(plot_device)
  }

  # Enrichment map: the network view of overlapping gene sets. The reported
  # terms are clustered by their gene membership (at most 12 groups), and each
  # group contributes at most ten terms that share no more than 60% of the
  # smaller set with an already-picked term. Without the clustering, one
  # programme (cell cycle or ribosome biogenesis) filled the whole map and the
  # remaining biology was invisible.
  for (ont in names(results)) {
    gse <- results[[ont]]
    df <- as.data.frame(gse)
    if (nrow(df) < 5L) next
    gene_sets <- tryCatch(clusterProfiler::geneInCategory(gse), error = record_optional_error)
    sel <- cluster_gsea_terms(df, gene_sets, max_clusters = 12L,
                              per_cluster = 10L, max_overlap = 0.6)
    if (is.null(sel) || nrow(sel$selected) < 2L) next
    sel_tab <- sel$selected
    sel_tab$cluster_num <- as.integer(sub("^C", "", sel_tab$cluster))
    sel_tab <- sel_tab[order(sel_tab$cluster_num, sel_tab$rank_in_cluster), , drop = FALSE]
    sel_tab$cluster_num <- NULL
    safe_write_csv(sel_tab, file.path(out_dir,
                  sprintf("GSEA_GO_%s_network_sets_%s.csv", ont, ctx$contrast_name)))
    ids <- sel$selected$ID
    sim <- sel$sim
    edges <- list()
    if (length(ids) > 1L) {
      for (i in seq_len(length(ids) - 1L)) {
        js <- sim[i, (i + 1L):length(ids)]
        hit <- which(js >= 0.1)
        if (length(hit)) {
          edges[[length(edges) + 1L]] <- data.frame(
            from = ids[i], to = ids[(i + 1L):length(ids)][hit],
            weight = as.numeric(js[hit]), stringsAsFactors = FALSE)
        }
      }
    }
    edges <- if (length(edges)) do.call(rbind, edges) else NULL
    if (!is.null(edges)) {
      g <- igraph::graph_from_data_frame(edges, vertices = data.frame(name = ids),
                                         directed = FALSE)
    } else {
      g <- igraph::make_empty_graph(n = length(ids), directed = FALSE)
      igraph::V(g)$name <- ids
    }
    lay <- tryCatch(igraph::layout_with_fr(g), error = record_optional_error)
    if (is.null(lay) || any(!is.finite(lay))) {
      ang <- seq_len(length(ids)) * 2 * pi / max(1L, length(ids))
      lay <- cbind(cos(ang), sin(ang))
    }
    node_df <- sel$selected
    cl_levels <- paste0("C", sort(as.integer(unique(sub("^C", "", node_df$cluster)))))
    node_df$cluster <- factor(node_df$cluster, levels = cl_levels)
    node_df$label <- wrap_label(node_df$Description, width = 30)
    idx <- match(igraph::V(g)$name, node_df$ID)
    igraph::V(g)$size <- as.numeric(node_df$setSize[idx])
    igraph::V(g)$cluster <- as.character(node_df$cluster[idx])
    igraph::V(g)$label <- node_df$label[idx]
    p <- ggraph::ggraph(g, layout = lay) +
      ggraph::geom_edge_link(ggplot2::aes(width = weight), colour = "grey65",
                             alpha = 0.6) +
      ggraph::scale_edge_width(range = c(0.3, 1.5), guide = "none") +
      ggraph::geom_node_point(ggplot2::aes(size = size, fill = cluster), shape = 21,
                              colour = "grey20", stroke = 0.45) +
      ggplot2::scale_fill_manual(
        values = stats::setNames(class_palette(nlevels(node_df$cluster)),
                                 levels(node_df$cluster)),
        breaks = levels(node_df$cluster),
        name = "Cluster") +
      ggplot2::scale_size_continuous(range = c(2.2, 9), name = "Genes in set") +
      ggraph::geom_node_text(ggplot2::aes(label = label), repel = TRUE, size = 3.8,
                             family = publication_font(), lineheight = 0.88,
                             max.overlaps = Inf, bg.color = "white", bg.r = 0.11,
                             segment.size = 0.25, min.segment.length = 0,
                             force = 1.6, box.padding = 0.3, point.padding = 0.18) +
      ggplot2::labs(
        title = paste0("GO ", ont, " enrichment network: ", ctx$contrast_name),
        subtitle = paste0(sel$n_clusters, " clusters; up to 10 terms each, ",
                          "sharing <= 60% of the smaller set; edge = shared genes")) +
      ggraph::theme_graph(base_family = publication_font(), base_size = 14) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 17, face = "bold"),
                     plot.subtitle = ggplot2::element_text(size = 13),
                     legend.title = ggplot2::element_text(face = "bold", size = 13),
                     legend.text = ggplot2::element_text(size = 12))
    save_plot(p, file.path(out_dir, sprintf("GSEA_GO_%s_enrichment_map_%s.png", ont, ctx$contrast_name)),
              key = "s3_go_enrichment_map")
  }

  # Alternative network view used by the original script via aPEAR.
  if (requireNamespace("enrichplot", quietly = TRUE)) {
    for (ont in names(results)) {
      gse <- results[[ont]]
      df <- as.data.frame(gse)
      if (nrow(df) < 5L) next
      p <- tryCatch(enrichplot::ridgeplot(gse, showCategory = 20) +
                      ggplot2::labs(title = paste0("GO ", ont, " GSEA ridge plot"), x = "log2 fold change") +
                      publication_theme(base_size = 15) +
                      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10)), error = record_optional_error)
      if (!is.null(p)) save_plot(p, file.path(out_dir, sprintf("GSEA_GO_%s_ridgeplot_%s.png", ont, ctx$contrast_name)),
                 key = "s3_go_ridgeplot")
    }
  }

  # --- extension: MSigDB collection GSEA (fgsea) -----------------------------
  if (isTRUE(cfg$analysis$run_msigdb) && requireNamespace("msigdbr", quietly = TRUE) &&
      requireNamespace("fgsea", quietly = TRUE)) {
    collection <- cfg$analysis$msigdb_category %||% "H"
    sets_db <- tryCatch(resolve_msigdb_gene_sets(cfg$species$name, collection,
      cfg$analysis$msigdb_subcategory, cfg$analysis$msigdb_db_species,
      file.path(out_dir, paste0("GSEA_MSigDB_", collection, "_source.json"))),
      error = record_optional_error)
    if (!is.null(sets_db) && nrow(sets_db) > 0L) {
      # MSigDB collections are keyed by gene symbol, while the ranked list used
      # by the GO/KEGG sections is keyed by Entrez ID; fgsea needs the symbol
      # ranking or it silently returns zero pathways.
      rank_symbol <- make_unique_ranked_vector(ctx$gsea_gene_df$SYMBOL, ctx$gsea_gene_df$logFC)
      sets <- split(sets_db$gene_symbol, sets_db$gs_name)
      fg <- tryCatch(fgsea::fgsea(pathways = sets, stats = rank_symbol, minSize = min_size,
                                  maxSize = max_size), error = record_optional_error)
      if (!is.null(fg) && nrow(fg) > 0L) {
        fg_out <- as.data.frame(fg)
        fg_out$leadingEdge <- vapply(fg_out$leadingEdge, paste, character(1), collapse = "/")
        safe_write_csv(fg_out, file.path(out_dir, paste0("GSEA_MSigDB_", collection, "_", ctx$contrast_name, ".csv")))
        fg_top <- utils::head(fg[order(fg$padj), , drop = FALSE], 20L)
        # The collection name is already in the panel title, so the axis carries
        # the set name alone: HALLMARK_TNFA_SIGNALING_VIA_NFKB -> TNFA_SIGNALING_VIA_NFKB.
        fg_top$pathway <- sub("^HALLMARK_", "", fg_top$pathway)
        fg_top$pathway <- factor(fg_top$pathway, levels = rev(fg_top$pathway))
        fg_top$direction <- ifelse(fg_top$NES > 0, "Activated", "Repressed")
        # Lollipop: the stick carries the sign and magnitude of NES, the head
        # carries the significance, so both can be read in one pass.
        fg_top$significance <- pmin(-log10(pmax(fg_top$padj, .Machine$double.xmin)), 50)
        fg_plot <- ggplot2::ggplot(fg_top, ggplot2::aes(NES, pathway)) +
          ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, y = pathway, yend = pathway,
                                             colour = direction), linewidth = 0.8) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35",
                              linewidth = 0.45) +
          ggplot2::geom_point(ggplot2::aes(size = significance, fill = direction),
                              shape = 21, colour = "grey20", stroke = 0.5) +
          ggplot2::scale_colour_manual(values = c(Activated = "#B2182B", Repressed = "#0072B2"),
                                       name = "Direction") +
          ggplot2::scale_fill_manual(values = c(Activated = "#B2182B", Repressed = "#0072B2"),
                                     guide = "none") +
          ggplot2::scale_size_continuous(range = c(2, 7),
                                         breaks = size_breaks(fg_top$significance),
                                         labels = function(x) ifelse(x >= 10, sprintf("%.0f", x),
                                                                     sprintf("%.1f", x)),
                                         name = expression("-log"[10] * " q")) +
          ggplot2::labs(x = "Normalized enrichment score", y = NULL,
                        title = wrap_label(paste0("MSigDB ", collection, " GSEA: ",
                                                  ctx$contrast_name), width = 30L)) +
          publication_theme(base_size = 13) +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 13),
                         legend.text = ggplot2::element_text(size = 13))
        save_plot(fg_plot, file.path(out_dir, paste0("GSEA_MSigDB_", collection, "_barplot_", ctx$contrast_name, ".png")),
                  key = "s3_msigdb_barplot")
      }
    }
  }

  ctx$gsea_go <- results
  write_status(file.path(out_dir, "status.json"), "section03_gsea_go", "completed", NULL,
               list(ontologies = names(results),
                    significant = lapply(results, function(g) sum(as.data.frame(g)$p.adjust < pcut))))
  ctx
}
