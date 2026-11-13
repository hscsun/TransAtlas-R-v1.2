# Section 2 - GO over-representation analysis.
#
# Port of input/AnaDGEsS3.R / input/AnaDGEsS3Mm.R. The original script took the
# ranked Entrez list, kept genes with |log2 fold change| above a permissive cut
# (0.75 in the mouse variant), ran enrichGO separately for MF, CC and BP with
# the original cutoffs, applied semantic simplification, then drew an
# enrichment map, a rich-factor dot plot grouped by ontology, and a
# concept-network plot. All of those outputs are reproduced here.

run_go_ora <- function(cfg, ctx, out_dir) {
  out_dir <- safe_dir(out_dir)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    write_status(file.path(out_dir, "status.json"), "section02_go_ora", "skipped", "clusterProfiler unavailable")
    return(ctx)
  }
  orgdb_name <- cfg$species$orgdb
  if (!requireNamespace(orgdb_name, quietly = TRUE)) {
    write_status(file.path(out_dir, "status.json"), "section02_go_ora", "skipped", paste("Missing", orgdb_name))
    return(ctx)
  }
  orgdb <- getExportedValue(orgdb_name, orgdb_name)
  rank <- ctx$universe_list
  if (is.null(rank) || length(rank) < 50L) {
    write_status(file.path(out_dir, "status.json"), "section02_go_ora", "skipped", "Eligible tested-gene universe unavailable")
    return(ctx)
  }

  limit <- as.numeric(cfg$analysis$lfc_threshold %||% 0.75)
  gd <- ctx$gene_df
  is_deg <- is.finite(gd$padj) & gd$padj < cfg$analysis$alpha & abs(gd$logFC) >= limit
  genes <- intersect(names(rank), unique(as.character(gd[[ctx$key_type]][is_deg])))
  universe <- names(rank)
  if (length(genes) < 10L) {
    write_status(file.path(out_dir, "status.json"), "section02_go_ora", "skipped",
                 paste0("Only ", length(genes), " genes pass the |log2FC| >= ", limit, " cut"))
    return(ctx)
  }

  # The original cutoffs are preserved: BP is stricter than MF/CC.
  ont_cfg <- list(
    BP = list(p = as.numeric(cfg$analysis$go_bp_pvalue_cutoff %||% 0.01),
              q = as.numeric(cfg$analysis$go_bp_qvalue_cutoff %||% 0.01)),
    MF = list(p = as.numeric(cfg$analysis$go_pvalue_cutoff %||% 0.05),
              q = as.numeric(cfg$analysis$go_qvalue_cutoff %||% 0.05)),
    CC = list(p = as.numeric(cfg$analysis$go_pvalue_cutoff %||% 0.05),
              q = as.numeric(cfg$analysis$go_qvalue_cutoff %||% 0.05))
  )
  simplify_cutoff <- as.numeric(cfg$analysis$go_simplify_cutoff %||% 0.5)

  # Symbol-named fold changes so concept networks can colour nodes by effect,
  # since readable = TRUE stores symbols in the geneID column.
  fc_symbols <- NULL
  if (!is.null(ctx$gene_df)) {
    fc_symbols <- ctx$gene_df$logFC
    names(fc_symbols) <- toupper(as.character(ctx$gene_df$SYMBOL))
    fc_symbols <- fc_symbols[!duplicated(names(fc_symbols))]
  }

  ora <- list()
  simplified <- list()
  summary_rows <- list()
  for (ont in c("BP", "MF", "CC")) {
    ego <- tryCatch(clusterProfiler::enrichGO(gene = genes, universe = universe, OrgDb = orgdb,
                                              keyType = ctx$key_type, ont = ont,
                                              pAdjustMethod = "BH",
                                              pvalueCutoff = ont_cfg[[ont]]$p,
                                              qvalueCutoff = ont_cfg[[ont]]$q,
                                              readable = TRUE),
                    error = record_optional_error)
    if (is.null(ego)) next
    df <- as.data.frame(ego)
    significant_df <- df
    safe_write_csv(significant_df, file.path(out_dir, paste0("ENRICHGO", ont, "significant_", ctx$contrast_name, ".csv")))
    if (nrow(df) == 0L) {
      # No term passes the configured cut-offs. Writing that empty frame left a
      # table with only a header, which reads like a broken run. Re-run the same
      # enrichment without a cut-off so the table still shows every term with its
      # p-value, and say so in the run log; nothing below plots these rows.
      ego_all <- tryCatch(clusterProfiler::enrichGO(gene = genes, universe = universe, OrgDb = orgdb,
                                                     keyType = ctx$key_type, ont = ont,
                                                     pAdjustMethod = "BH",
                                                     pvalueCutoff = 1, qvalueCutoff = 1,
                                                     readable = TRUE),
                          error = record_optional_error)
      if (!is.null(ego_all)) {
        df <- as.data.frame(ego_all)
        message("    GO ", ont, ": no term at the configured cut-off; ",
                nrow(df), " terms written without a cut-off")
      }
    }
    safe_write_csv(df, file.path(out_dir, paste0("ENRICHGO", ont, "alls_", ctx$contrast_name, ".csv")))
    if (!nrow(significant_df)) next
    df <- significant_df
    ora[[ont]] <- ego

    simp <- tryCatch(clusterProfiler::simplify(ego, cutoff = simplify_cutoff, by = "p.adjust",
                                               select_fun = min, measure = "Wang"),
                     error = record_optional_error)
    if (!is.null(simp)) {
      simplified[[ont]] <- simp
      safe_write_csv(as.data.frame(simp), file.path(out_dir, paste0("ENRICHGO", ont, "s_", ctx$contrast_name, ".csv")))
    }
    summary_rows[[ont]] <- data.frame(ontology = ont, terms = nrow(df),
                                      terms_simplified = if (is.null(simp)) NA_integer_ else nrow(as.data.frame(simp)))
  }
  if (!length(ora)) {
    ctx$go_ora <- list()
    write_status(file.path(out_dir, "status.json"), "section02_go_ora", "completed", "No GO term passes the configured cut-offs",
                 list(terms = 0L))
    return(ctx)
  }
  safe_write_csv(do.call(rbind, summary_rows), file.path(out_dir, "GO_ORA_summary.csv"))

  for (ont in names(ora)) {
    ego <- ora[[ont]]
    ego_df <- as.data.frame(ego)
    # enrichGO returns parent and child terms together, so the first 20 rows by
    # adjusted P can be several variants of one programme. Terms whose genes
    # overlap a term already chosen by more than 40% of the smaller set are
    # skipped, and the panel shows the distinct terms that remain.
    gene_sets <- tryCatch(clusterProfiler::geneInCategory(ego), error = record_optional_error)
    sel_ids <- if (!is.null(gene_sets)) {
      diverse_term_selection(ego_df$ID, gene_sets, n_total = 20L, max_overlap = 0.4)
    } else {
      utils::head(as.character(ego_df$ID), 20L)
    }
    if (!length(sel_ids)) next
    # enrichplot selects a character showCategory by Description, not by ID.
    sel_desc <- ego_df$Description[match(sel_ids, as.character(ego_df$ID))]
    sel_desc <- sel_desc[!is.na(sel_desc)]
    if (!length(sel_desc)) next
    # The legend breaks are clipped to the Count range of the terms this panel
    # actually draws (showCategory = 20), not the whole ontology: breaks outside
    # the layer's limits are dropped, so a range taken from all 163 terms left
    # the smallest key at 100 while the panel showed dots down to ~44, which
    # reads as "legend dots are bigger than plot dots".
    count_vec <- ego_df$Count[match(sel_ids, as.character(ego_df$ID))]
    p <- tryCatch(suppressMessages(
                    enrichplot::dotplot(ego, showCategory = sel_desc, font.size = 12) +
                    ggplot2::scale_size_continuous(range = c(2, 7),
                                                   breaks = size_breaks(count_vec),
                                                   name = "Gene count")) +
                    ggplot2::labs(title = paste("GO", ont, "over-representation:", ctx$contrast_name)) +
                    publication_theme() +
                    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 12),
                                   axis.text.x = ggplot2::element_text(size = 11),
                                   legend.text = ggplot2::element_text(size = 11)),
                  error = record_optional_error)
    if (!is.null(p)) save_plot(p, file.path(out_dir, paste0("GO_", ont, "_dotplot_", ctx$contrast_name, ".png")),
                               key = "s2_go_dotplot", n_terms = length(sel_desc))
  }

  # Functional classes of the enriched BP terms, plus the term-gene network with
  # the genes coloured by their fold change. These two replace enrichplot's
  # emapplot and cnetplot, whose labels collided and whose gene nodes were grey.
  classes <- NULL
  if (!is.null(ora$BP)) {
    classes <- tryCatch(go_functional_classes(ora$BP, cfg, fc_symbols = fc_symbols,
                                              contrast = ctx$contrast_name),
                        error = record_optional_error)
    if (!is.null(classes)) {
      if (!is.null(classes$cluster_plot)) {
        save_plot(classes$cluster_plot,
                  file.path(out_dir, paste0("GO_BP_functional_classes_", ctx$contrast_name, ".png")),
                  key = "s2_go_classes")
      }
      if (!is.null(classes$table)) {
        safe_write_csv(classes$table,
                       file.path(out_dir, paste0("GO_BP_functional_classes_", ctx$contrast_name, ".csv")))
      }
      if (!is.null(classes$gene_plots) && length(classes$gene_plots)) {
        # One term-gene network per functional class: the class's largest term
        # as hub, its top 100 genes by |log2FC| around it.
        for (i in seq_along(classes$gene_plots)) {
          save_plot(classes$gene_plots[[i]],
                    file.path(out_dir, sprintf("GO_term_gene_network_class%02d_%s.png",
                                               i, ctx$contrast_name)),
                    key = "s2_go_term_gene")
        }
      }
      if (!is.null(classes$gene_terms)) {
        safe_write_csv(classes$gene_terms,
                       file.path(out_dir, paste0("GO_term_gene_network_terms_", ctx$contrast_name, ".csv")))
      }
    }
  }
  if (is.null(classes)) message("    section 2: functional-class network unavailable (ggraph/pairwise_termsim)")

  # Rich-factor dot plot grouped by ontology (original "go0" panel).
  rf <- do.call(rbind, lapply(names(ora), function(ont) {
    df <- as.data.frame(ora[[ont]])
    if (nrow(df) == 0L) return(NULL)
    df$Group <- ont
    df
  }))
  if (!is.null(rf) && nrow(rf) > 0L) {
    rf$richFactor <- rf$Count / as.numeric(sub("/\\d+", "", rf$BgRatio))
    top_rf <- do.call(rbind, lapply(split(rf, rf$Group), function(x) utils::head(x[order(x$p.adjust), , drop = FALSE], 10L)))
    top_rf$Description <- factor(top_rf$Description, levels = rev(unique(top_rf$Description)))
    p <- ggplot2::ggplot(top_rf, ggplot2::aes(richFactor, Description, colour = Group)) +
      ggplot2::geom_point(ggplot2::aes(size = Count)) +
      ggplot2::scale_size_continuous(range = c(1, 7), breaks = size_breaks(top_rf$Count),
                                     name = "Gene count") +
      ggplot2::labs(x = "Rich factor", y = NULL,
                    title = paste("GO over-representation by rich factor:", ctx$contrast_name)) +
      publication_theme() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 11),
                                            axis.text.x = ggplot2::element_text(size = 11),
                                            legend.text = ggplot2::element_text(size = 13),
                                            legend.title = ggplot2::element_text(size = 13,
                                                                                 face = "bold"),
                                            legend.key.height = ggplot2::unit(13, "pt"))
    save_plot(p, file.path(out_dir, paste0("GO_richfactor_by_ontology_", ctx$contrast_name, ".png")),
              key = "s2_richfactor")
  }

  # The concept network is drawn above from the same BP result, with the gene
  # nodes coloured by fold change; enrichplot::cnetplot is no longer used.

  ctx$go_ora <- ora
  ctx$go_simplified <- simplified
  ctx$go_genes <- genes
  ctx$go_classes <- if (!is.null(classes)) classes$table else NULL
  write_status(file.path(out_dir, "status.json"), "section02_go_ora", "completed", NULL,
               list(lfc_cut = limit, genes_tested = length(genes),
                    ontologies = names(ora),
                    functional_classes = if (is.null(classes)) NA_integer_ else
                      length(unique(classes$table$class)),
                    term_gene_networks = if (is.null(classes)) 0L else length(classes$gene_plots),
                    terms = vapply(ora, function(x) nrow(as.data.frame(x)), integer(1))))
  ctx
}
