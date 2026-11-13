# Section 10 - protein-protein interaction network of the DEG-encoded proteins.
#
# The interaction evidence comes from STRING v12 through its public API, so no
# large reference download is needed: the DEG symbols are sent in batches, the
# returned edge list is cached under cache/string/, and the network is drawn
# with igraph/ggraph. Node colour keeps the expression direction of the gene,
# node size is the STRING degree (number of interaction partners inside the
# network), and Louvain communities are reported as candidate protein modules.

run_ppi <- function(cfg, ctx, out_dir) {
  out_dir <- safe_dir(out_dir)
  if (!isTRUE(cfg$ppi$run %||% TRUE)) {
    write_status(file.path(out_dir, "status.json"), "section10_ppi", "skipped",
                 "run_ppi = FALSE")
    return(ctx)
  }
  needed <- c("igraph", "ggplot2")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    write_status(file.path(out_dir, "status.json"), "section10_ppi", "skipped",
                 paste("Missing package(s):", paste(missing, collapse = ", ")))
    return(ctx)
  }

  # --- DEG list, identical to the screening section ------------------------
  tab <- screening_tables(cfg, ctx)
  thresholds <- tab$thresholds
  primary <- as.numeric(cfg$analysis$lfc_threshold %||% 0.75)
  deg <- tab$mylist[[which.min(abs(thresholds - primary))]]
  deg <- deg[deg$Direction %in% c("Up", "Down"), , drop = FALSE]
  deg <- deg[!duplicated(deg$GenesSymbol), , drop = FALSE]
  if (!nrow(deg)) {
    write_status(file.path(out_dir, "status.json"), "section10_ppi", "skipped",
                 "No differentially expressed genes at the primary threshold")
    return(ctx)
  }

  max_nodes <- as.integer(cfg$ppi$max_nodes %||% 300L)
  deg <- deg[order(deg$GenesQvalue, -abs(deg$GenesFoldChange_raw)), , drop = FALSE]
  if (nrow(deg) > max_nodes) deg <- utils::head(deg, max_nodes)
  symbols <- unique(stats::na.omit(as.character(deg$GenesSymbol)))
  symbols <- symbols[nzchar(symbols)]

  taxon <- cfg$species$string_taxon %||%
    switch(cfg$species$common_name, human = 9606L, rat = 10116L, 10090L)

  # --- STRING API ----------------------------------------------------------
  # The API accepts a carriage-return separated identifier list. Batches of 100
  # keep the URL short enough for the server, and the completed edge table is
  # cached by a hash of the identifier set so a repeated run is offline.
  cache_dir <- cfg$ppi$cache_dir %||% "cache/string"
  cache_dir <- resolve_project_path(cache_dir, attr(cfg, "project_root") %||% getwd())
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  key_file <- tempfile()
  writeLines(sort(symbols), key_file)
  cache_key <- unname(tools::md5sum(key_file))
  unlink(key_file)
  cache_path <- file.path(cache_dir, sprintf("string_%s_%s_%03d.tsv", taxon, cache_key,
                                             round(1000 * (cfg$ppi$score_threshold %||% 0.4))))

  if (file.exists(cache_path)) {
    edges <- utils::read.delim(cache_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    old_timeout <- getOption("timeout")
    options(timeout = max(old_timeout, as.numeric(cfg$ppi$timeout %||% 120)))
    on.exit(options(timeout = old_timeout), add = TRUE)
    score <- round(1000 * as.numeric(cfg$ppi$score_threshold %||% 0.4))
    chunks <- split(symbols, ceiling(seq_along(symbols) / 100L))
    batch_errors <- character()
    parts <- lapply(seq_along(chunks), function(batch) {
      ids <- chunks[[batch]]
      query <- paste(utils::URLencode(ids, reserved = TRUE), collapse = "%0d")
      url <- sprintf(paste0("https://string-db.org/api/tsv/network?identifiers=%s",
                            "&species=%d&required_score=%d&caller_identity=ClawsGO_RNAseq"),
                     query, as.integer(taxon), score)
      out <- tryCatch(
        utils::read.delim(url, check.names = FALSE, stringsAsFactors = FALSE),
        error = function(e) {
          batch_errors <<- c(batch_errors, paste0("Batch ", batch, ": ", conditionMessage(e)))
          NULL
        })
      if (!is.null(out) && nrow(out) && !all(c("preferredName_A", "preferredName_B", "score") %in% names(out))) {
        batch_errors <<- c(batch_errors, paste0("Batch ", batch, ": invalid STRING response: ", paste(names(out), collapse = ", ")))
        return(NULL)
      }
      if (is.null(out) || !nrow(out)) NULL else out
    })
    if (length(batch_errors)) stop(paste(batch_errors, collapse = "; "), call. = FALSE)
    parts <- parts[!vapply(parts, is.null, logical(1))]
    if (!length(parts)) {
      write_status(file.path(out_dir, "status.json"), "section10_ppi", "skipped",
                   "STRING returned no interaction for the DEG list (offline or no edges)")
      return(ctx)
    }
    columns <- unique(unlist(lapply(parts, names)))
    parts <- lapply(parts, function(part) {
      for (column in setdiff(columns, names(part))) part[[column]] <- NA
      part[, columns, drop = FALSE]
    })
    edges <- do.call(rbind, parts)
    edges <- edges[!duplicated(edges[, c("preferredName_A", "preferredName_B")]), , drop = FALSE]
    utils::write.table(edges, cache_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  if (!nrow(edges)) {
    write_status(file.path(out_dir, "status.json"), "section10_ppi", "skipped",
                 "STRING returned an empty edge list")
    return(ctx)
  }
  edges$score <- as.numeric(edges$score)
  edges <- edges[edges$score >= as.numeric(cfg$ppi$score_threshold %||% 0.4), , drop = FALSE]
  if (!nrow(edges)) {
    write_status(file.path(out_dir, "status.json"), "section10_ppi", "skipped",
                 "No STRING edge passed the score threshold")
    return(ctx)
  }

  # --- graph and node annotation ------------------------------------------
  g <- igraph::graph_from_data_frame(
    edges[, c("preferredName_A", "preferredName_B", "score")], directed = FALSE)
  g <- igraph::simplify(g, edge.attr.comb = list(score = "max"))
  deg_lookup <- stats::setNames(deg$GenesFoldChange_raw, toupper(deg$GenesSymbol))
  dir_lookup <- stats::setNames(as.character(deg$Direction), toupper(deg$GenesSymbol))
  node_names <- igraph::V(g)$name
  lfc <- unname(deg_lookup[toupper(node_names)])
  dirn <- unname(dir_lookup[toupper(node_names)])
  dirn[is.na(dirn)] <- "Unclassified"

  set.seed(as.integer(cfg$project$seed %||% 1L))
  comm <- igraph::cluster_louvain(g)
  memb <- as.integer(igraph::membership(comm))
  max_modules <- as.integer(cfg$ppi$max_modules %||% 10L)
  # Louvain's resolution is data-driven, so a dense network can split into a
  # long tail of two-protein communities. Merge the smallest ones into the
  # community they share most edges with until the requested ceiling is met.
  if (length(unique(memb)) > max_modules) {
    memb <- reduce_communities(g, memb, max_modules)
  }
  # Renumber by size: M1 is the largest module.
  sizes <- sort(table(memb), decreasing = TRUE)
  remap <- stats::setNames(seq_along(sizes), names(sizes))
  memb <- as.integer(remap[as.character(memb)])
  nodes <- data.frame(
    gene = node_names,
    log2FC = as.numeric(lfc),
    direction = dirn,
    degree = as.numeric(igraph::degree(g)),
    betweenness = as.numeric(igraph::betweenness(g)),
    module = memb,
    stringsAsFactors = FALSE)
  nodes <- nodes[order(-nodes$degree, -abs(nodes$log2FC)), , drop = FALSE]
  safe_write_csv(nodes, file.path(out_dir, paste0("PPI_nodes_", ctx$contrast_name, ".csv")))

  edge_out <- data.frame(
    from = edges$preferredName_A, to = edges$preferredName_B, score = edges$score,
    stringsAsFactors = FALSE)
  safe_write_csv(edge_out, file.path(out_dir, paste0("PPI_edges_", ctx$contrast_name, ".csv")))

  modules <- do.call(rbind, lapply(sort(unique(nodes$module)), function(m) {
    sub <- nodes[nodes$module == m, , drop = FALSE]
    data.frame(module = m, n_proteins = nrow(sub),
               n_up = sum(sub$direction == "Up"), n_down = sum(sub$direction == "Down"),
               top_hubs = paste(utils::head(sub$gene[order(-sub$degree)], 5L), collapse = ", "),
               stringsAsFactors = FALSE)
  }))
  modules <- modules[order(modules$module), , drop = FALSE]

  # --- functional annotation of every module --------------------------------
  # A community is only a candidate module until it has a function attached.
  # Each module's proteins are tested against the tested-gene universe with GO
  # over-representation in all three ontologies (BP first, then MF and CC), so a
  # module whose biology is enzymatic or structural still gets a description.
  # The best terms travel into the module table, the network labels and a
  # dedicated lollipop panel.
  alpha <- as.numeric(cfg$analysis$alpha %||% 0.05)
  annotation_alpha <- as.numeric(cfg$ppi$module_annotation_alpha %||% 0.2)
  module_go <- NULL
  module_fun <- NULL
  onts <- c("BP", "MF", "CC")
  if (isTRUE(cfg$ppi$module_go %||% TRUE) && requireNamespace("clusterProfiler", quietly = TRUE)) {
    orgdb_name <- cfg$species$orgdb
    if (!is.null(orgdb_name) && requireNamespace(orgdb_name, quietly = TRUE)) {
      orgdb <- getExportedValue(orgdb_name, orgdb_name)
      # When section 10 runs on its own, the ranked context is not attached;
      # fall back to the screened symbols so the GO universe is never empty.
      universe_sym <- if (!is.null(ctx$gene_df) && "SYMBOL" %in% names(ctx$gene_df)) {
        unique(as.character(ctx$gene_df$SYMBOL))
      } else {
        unique(as.character(ctx$sig_symbols))
      }
      universe_sym <- universe_sym[!is.na(universe_sym) & nzchar(universe_sym)]
      key2sym <- stats::setNames(deg$GenesSymbol, toupper(deg$GenesSymbol))
      # GO annotations are built once for the whole universe and reused by every
      # module. enrichGO() remaps the universe on every call, which costs about
      # 90 s per ontology; enricher() against these prebuilt maps takes a few
      # seconds, so annotating ten modules stays practical.
      go_maps <- NULL
      if (requireNamespace("GO.db", quietly = TRUE)) {
        ann <- tryCatch(suppressMessages(AnnotationDbi::select(
          orgdb, keys = universe_sym, columns = c("GOALL", "ONTOLOGYALL"),
          keytype = ctx$symbol_type)), error = record_optional_error)
        if (!is.null(ann) && nrow(ann)) {
          ann <- ann[!is.na(ann$GOALL) & !is.na(ann$ONTOLOGYALL), , drop = FALSE]
          ann <- unique(ann[, c("GOALL", "ONTOLOGYALL", ctx$symbol_type)])
          names(ann)[3L] <- "gene"
          term_names <- tryCatch(suppressMessages(AnnotationDbi::select(
            getExportedValue("GO.db", "GO.db"), keys = unique(ann$GOALL),
            columns = "TERM", keytype = "GOID")), error = record_optional_error)
          if (!is.null(term_names) && nrow(term_names)) {
            names(term_names) <- c("term", "name")
            go_maps <- stats::setNames(lapply(onts, function(ont) {
              t2g <- unique(ann[ann$ONTOLOGYALL == ont, c("GOALL", "gene")])
              names(t2g) <- c("term", "gene")
              list(t2g = t2g, t2n = term_names)
            }), onts)
          }
        }
      }
      go_rows <- list()
      fun_rows <- list()
      for (m in sort(unique(nodes$module))) {
        sub <- nodes[nodes$module == m, , drop = FALSE]
        genes <- unique(unname(key2sym[toupper(sub$gene)]))
        genes <- genes[!is.na(genes) & nzchar(genes)]
        best <- list()
        for (ont in onts) {
          if (length(genes) < 3L || is.null(go_maps[[ont]])) next
          ego <- tryCatch(clusterProfiler::enricher(
            gene = genes, universe = universe_sym,
            TERM2GENE = go_maps[[ont]]$t2g, TERM2NAME = go_maps[[ont]]$t2n,
            pAdjustMethod = "BH", pvalueCutoff = annotation_alpha, qvalueCutoff = annotation_alpha,
            minGSSize = 5, maxGSSize = 500), error = record_optional_error)
          df <- if (is.null(ego)) NULL else as.data.frame(ego)
          if (is.null(df) || !nrow(df)) next
          df$significant <- !is.na(df$p.adjust) & df$p.adjust < alpha
          sig <- df[df$significant, , drop = FALSE]
          if (nrow(sig)) df <- sig
          df <- df[order(df$p.adjust, -df$Count), , drop = FALSE]
          best[[ont]] <- df
          keep <- utils::head(df, 10L)
          go_rows[[paste(m, ont)]] <- data.frame(
            module = m, ontology = ont,
            keep[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue",
                     "p.adjust", "qvalue", "Count", "geneID", "significant")],
            stringsAsFactors = FALSE)
        }
        top_of <- function(ont, n = 5L) {
          if (is.null(best[[ont]])) return(NA_character_)
          paste(utils::head(best[[ont]]$Description, n), collapse = " / ")
        }
        primary <- if (!is.null(best$BP)) "BP" else if (!is.null(best$MF)) "MF" else
          if (!is.null(best$CC)) "CC" else NA_character_
        n_terms <- if (length(best)) sum(vapply(best, nrow, integer(1))) else 0L
        fun_rows[[as.character(m)]] <- data.frame(
          module = m, n_module_genes = length(genes), n_go_terms = as.integer(n_terms),
          top_GO_BP = top_of("BP"), top_GO_MF = top_of("MF"), top_GO_CC = top_of("CC"),
          annotation_ontology = primary,
          top_function = if (is.na(primary)) {
            paste0("isolated group: ", paste(utils::head(genes, 3L), collapse = ", "))
          } else as.character(best[[primary]]$Description[1]),
          best_padj = if (is.na(primary)) NA_real_ else as.numeric(best[[primary]]$p.adjust[1]),
          stringsAsFactors = FALSE)
      }
      if (length(go_rows)) {
        module_go <- do.call(rbind, go_rows)
        rownames(module_go) <- NULL
        safe_write_csv(module_go,
                       file.path(out_dir, paste0("PPI_module_GO_", ctx$contrast_name, ".csv")))
      }
      if (length(fun_rows)) module_fun <- do.call(rbind, fun_rows)
    }
  }
  if (is.null(module_fun)) {
    module_fun <- data.frame(module = modules$module, n_module_genes = NA_integer_,
                             n_go_terms = NA_integer_, top_GO_BP = NA_character_,
                             top_GO_MF = NA_character_, top_GO_CC = NA_character_,
                             annotation_ontology = NA_character_,
                             top_function = NA_character_, best_padj = NA_real_,
                             stringsAsFactors = FALSE)
  }
  modules <- merge(modules, module_fun, by = "module", all.x = TRUE)
  short_fun <- function(x, n = 34L) {
    x <- ifelse(is.na(x), "unannotated", x)
    ifelse(nchar(x) > n, paste0(substr(x, 1L, n - 1L), "…"), x)
  }
  modules$module_label <- paste0("M", modules$module, ": ", short_fun(modules$top_function))
  modules <- modules[order(modules$module), , drop = FALSE]
  safe_write_csv(modules, file.path(out_dir, paste0("PPI_modules_", ctx$contrast_name, ".csv")))

  # --- network figure ------------------------------------------------------
  lay <- igraph::layout_with_fr(g)
  lay_df <- data.frame(name = node_names, x = lay[, 1L], y = lay[, 2L],
                       stringsAsFactors = FALSE)
  lay_df <- merge(lay_df, nodes, by.x = "name", by.y = "gene", all.x = TRUE)
  edge_df <- data.frame(
    x = lay_df$x[match(edges$preferredName_A, lay_df$name)],
    y = lay_df$y[match(edges$preferredName_A, lay_df$name)],
    xend = lay_df$x[match(edges$preferredName_B, lay_df$name)],
    yend = lay_df$y[match(edges$preferredName_B, lay_df$name)],
    score = edges$score)
  edge_df <- edge_df[stats::complete.cases(edge_df), , drop = FALSE]

  # Module membership is encoded as the node border and spelled out in the
  # legend ("M4 locomotion"). Text labels inside the graph collided with the hub
  # names and the leader lines, so the functions live in the legend instead.
  mod_levels <- paste0("M", modules$module)
  mod_labels <- paste0(mod_levels, "  ", short_fun(modules$top_function, 34L))
  mod_values <- stats::setNames(grDevices::hcl.colors(length(mod_levels), "Dark 3"), mod_levels)
  lay_df$module_f <- factor(paste0("M", lay_df$module), levels = mod_levels)

  # Dashed outlines around every module, so the communities read as groups and
  # not only through the node-border colour. A convex hull follows the actual
  # extent of a module and, with only a small outward pad, stays tighter than a
  # covariance ellipse (which grew large enough to overlap the neighbouring
  # modules). Two-node modules get a small capsule-like ellipse because a hull
  # of two points is just a line.
  module_outlines <- do.call(rbind, lapply(split(lay_df, lay_df$module), function(d) {
    xy <- as.matrix(d[, c("x", "y")])
    cen <- colMeans(xy)
    if (nrow(xy) >= 3L) {
      pts <- xy[grDevices::chull(xy), , drop = FALSE]
      span <- max(diff(range(xy[, 1L])), diff(range(xy[, 2L])), 1)
      pad <- 0.03 * span + 0.25
      vec <- pts - matrix(cen, nrow(pts), 2L, byrow = TRUE)
      len <- sqrt(rowSums(vec^2))
      len[len == 0] <- 1
      pts <- pts + vec / len * pad
      pts <- rbind(pts, pts[1L, , drop = FALSE])
    } else {
      v <- if (nrow(xy) == 2L) xy[2L, ] - xy[1L, ] else c(1, 0)
      ang <- atan2(v[2L], v[1L])
      rot <- matrix(c(cos(ang), -sin(ang), sin(ang), cos(ang)), 2L, 2L)
      proj <- (xy - matrix(cen, nrow(xy), 2L, byrow = TRUE)) %*% rot
      rx <- max(abs(proj[, 1L])) + 0.45
      ry <- max(abs(proj[, 2L])) + 0.45
      tt <- seq(0, 2 * pi, length.out = 120L)
      pts <- cbind(rx * cos(tt), ry * sin(tt)) %*% t(rot)
      pts <- sweep(pts, 2L, cen, "+")
    }
    data.frame(x = pts[, 1L], y = pts[, 2L],
               module_f = d$module_f[1L], stringsAsFactors = FALSE)
  }))
  module_layer <- ggplot2::geom_path(
    data = module_outlines,
    ggplot2::aes(x = x, y = y, group = module_f, colour = module_f),
    linetype = "dashed", linewidth = 0.8, alpha = 0.95, show.legend = FALSE)

  # Label at most `ppi$label_n` proteins per module (default 4), always choosing
  # the highest-degree members of that module. Large modules therefore show
  # four hub proteins, while small modules show fewer rather than being drowned
  # out by labels from the dominant communities.
  label_cap <- max(1L, as.integer(cfg$ppi$label_n %||% 4L))
  quota <- function(n) {
    if (n >= 20L) 4L else if (n >= 10L) 3L else if (n >= 5L) 2L else 1L
  }
  lab_df <- do.call(rbind, lapply(split(lay_df, lay_df$module), function(d) {
    k <- min(label_cap, quota(nrow(d)), nrow(d))
    utils::head(d[order(-d$degree), , drop = FALSE], k)
  }))
  rownames(lab_df) <- NULL
  cap <- display_lfc_limit(cfg, lay_df$log2FC)
  net <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edge_df,
                          ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                                       linewidth = score),
                          colour = "grey72", alpha = 0.45) +
    module_layer +
    ggplot2::geom_point(data = lay_df,
                        ggplot2::aes(x = x, y = y, size = degree, fill = log2FC,
                                     colour = module_f),
                        shape = 21, stroke = 1.2) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-cap, cap),
                                  oob = scales::squish,
                                  name = expression("log"[2] * "FC")) +
    ggplot2::scale_size_continuous(range = c(1, 8),
                                   breaks = size_breaks(lay_df$degree),
                                   name = "STRING degree") +
    ggplot2::scale_colour_manual(values = mod_values, labels = mod_labels,
                                 name = "Module (node border)") +
    ggplot2::scale_linewidth_continuous(range = c(0.15, 1.2),
                                        breaks = size_breaks(edge_df$score),
                                        labels = function(x) sprintf("%.2f", x),
                                        name = "STRING score") +
    # Transparent text labels: no white box, no border, only a thin leader line.
    ggrepel::geom_text_repel(data = lab_df,
                             ggplot2::aes(x = x, y = y, label = name),
                             size = 4.2, fontface = "bold", seed = 1,
                             box.padding = 0.55, point.padding = 0.3,
                             min.segment.length = 0, max.overlaps = Inf,
                             segment.size = 0.3, segment.colour = "grey40",
                             force = 10, max.iter = 50000) +
    ggplot2::labs(title = paste0("Protein-protein interaction network of DEG-encoded proteins: ",
                                 ctx$contrast_name),
                  subtitle = paste0("STRING v12, confidence >= ",
                                    cfg$ppi$score_threshold %||% 0.4,
                                    "; node fill = log2FC, node border = module, node size = interaction degree; ",
                                    nrow(nodes), " proteins, ", nrow(edge_out), " interactions, ",
                                    nrow(modules), " modules")) +
    publication_theme() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   plot.caption = ggplot2::element_text(size = 8, hjust = 0, colour = "grey30"),
                   legend.text = ggplot2::element_text(size = 12),
                   legend.title = ggplot2::element_text(size = 14, face = "bold"),
                   legend.position = "right")
  save_plot(net, file.path(out_dir, sprintf("PPI_network_%s.png", ctx$contrast_name)),
            key = "s10_ppi_network")

  # --- hub proteins --------------------------------------------------------
  # Thirty proteins by STRING degree, but not purely by degree: a panel drawn
  # from a strongly induced network can contain almost no down-regulated node.
  # At least 5% of the displayed rows (two of thirty) are reserved for each
  # direction; the lowest-degree rows of the over-represented direction make
  # room for the strongest remaining proteins of the other direction.
  hub_n <- min(30L, nrow(nodes))
  min_each <- ceiling(0.05 * hub_n)
  ranked <- nodes[order(-nodes$degree, -abs(nodes$log2FC)), , drop = FALSE]
  hub <- utils::head(ranked, hub_n)
  for (dir_need in c("Up", "Down")) {
    shortage <- min_each - sum(hub$direction == dir_need, na.rm = TRUE)
    if (shortage <= 0L) next
    extra <- ranked[ranked$direction == dir_need & !(ranked$gene %in% hub$gene), , drop = FALSE]
    extra <- utils::head(extra, shortage)
    if (!nrow(extra)) next
    donor_dir <- setdiff(c("Up", "Down"), dir_need)
    donors <- hub[hub$direction %in% donor_dir, , drop = FALSE]
    donors <- donors[order(donors$degree), , drop = FALSE]
    spare <- sum(hub$direction %in% donor_dir, na.rm = TRUE) - min_each
    n_swap <- min(nrow(extra), max(0L, spare))
    if (n_swap <= 0L) next
    drop_genes <- utils::head(donors$gene, n_swap)
    hub <- hub[!(hub$gene %in% drop_genes), , drop = FALSE]
    hub <- rbind(hub, extra[seq_len(n_swap), , drop = FALSE])
    hub <- hub[order(-hub$degree, -abs(hub$log2FC)), , drop = FALSE]
  }
  hub$gene <- factor(hub$gene, levels = rev(hub$gene))
  hub_plot <- ggplot2::ggplot(hub, ggplot2::aes(degree, gene, fill = direction)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_text(ggplot2::aes(label = degree), hjust = -0.18, size = 3.65) +
    ggplot2::scale_fill_manual(values = c(Up = "#B2182B", Down = "#2166AC"),
                               name = "DEG direction") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(x = "Number of interaction partners (STRING degree)", y = NULL,
                  title = "Hub proteins among the DEG-encoded proteins",
                  subtitle = paste0("Top ", nrow(hub), " hub proteins; at least 5% per direction")) +
    publication_theme(base_size = 13) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 11, face = "bold"))
  save_plot(hub_plot, file.path(out_dir, sprintf("PPI_hub_proteins_%s.png", ctx$contrast_name)),
            key = "s10_ppi_hubs")

  # --- module functions ----------------------------------------------------
  # The two strongest terms per module from its best-covered ontology, one
  # lollipop row each: the panel that turns "community 4" into
  # "ribosome biogenesis".
  if (!is.null(module_go) && nrow(module_go)) {
    top_terms <- do.call(rbind, lapply(split(module_go, module_go$module), function(d) {
      ont <- c("BP", "MF", "CC")[c("BP", "MF", "CC") %in% d$ontology][1]
      d <- d[d$ontology == ont, , drop = FALSE]
      d <- d[order(d$p.adjust, -d$Count), , drop = FALSE]
      utils::head(d, 2L)
    }))
    top_terms$neglog <- -log10(pmax(top_terms$p.adjust, .Machine$double.xmin))
    top_terms$term_label <- paste0("M", top_terms$module, " [", top_terms$ontology, "]  ",
                                   ifelse(nchar(top_terms$Description) > 42,
                                          paste0(substr(top_terms$Description, 1L, 41L), "…"),
                                          top_terms$Description))
    top_terms$module_f <- factor(paste0("M", top_terms$module),
                                 levels = rev(paste0("M", sort(unique(top_terms$module)))))
    top_terms <- top_terms[order(top_terms$module, -top_terms$neglog), , drop = FALSE]
    top_terms$term_label <- factor(top_terms$term_label, levels = rev(top_terms$term_label))
    n_mod <- length(unique(top_terms$module))
    p_fun <- ggplot2::ggplot(top_terms,
                             ggplot2::aes(x = neglog, y = term_label, colour = module_f)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = neglog, yend = term_label),
                            linewidth = 0.6, colour = "grey80") +
      ggplot2::geom_point(ggplot2::aes(size = Count)) +
      ggplot2::scale_colour_hue(name = "Module", h.start = 20, c = 70, l = 55) +
      ggplot2::scale_size_continuous(range = c(2, 7), breaks = size_breaks(top_terms$Count),
                                     name = "Module genes in term") +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
      ggplot2::labs(x = expression("-log"[10] * " adjusted P (GO)"), y = NULL,
                    title = paste0("Functional annotation of PPI modules: ", ctx$contrast_name),
                    subtitle = paste0("Top 2 enriched GO terms per module (best-covered ontology)\n",
                                      n_mod, " annotated modules; colour = module")) +
      # +2 pt throughout. Only the colour guide gets larger keys; the size
      # guide must keep its own point sizes or every key in it looks the same.
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 6.4)),
                      size = ggplot2::guide_legend(override.aes = list(colour = "grey20"))) +
      publication_theme(base_size = 13) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10.6),
                     legend.text = ggplot2::element_text(size = 13.5),
                     legend.title = ggplot2::element_text(size = 15, face = "bold"),
                     legend.key.size = ggplot2::unit(0.8, "cm"),
                     plot.margin = ggplot2::margin(6, 14, 6, 14),
                     panel.grid.major.y = ggplot2::element_blank(),
                     legend.position = "right")
    save_plot(p_fun, file.path(out_dir, sprintf("PPI_module_functions_%s.png", ctx$contrast_name)),
              key = "s10_ppi_functions", n = nrow(top_terms))
  }

  top_hub <- nodes$gene[which.max(nodes$degree)]
  write_status(file.path(out_dir, "status.json"), "section10_ppi", "completed", NULL,
               list(proteins = nrow(nodes), interactions = nrow(edge_out),
                    modules = nrow(modules), max_modules = max_modules,
                    modules_annotated = if (is.null(module_go)) 0L else length(unique(module_go$module)),
                    top_hub = top_hub,
                    score_threshold = as.numeric(cfg$ppi$score_threshold %||% 0.4),
                    source = "STRING v12 API"))
  ctx$ppi_nodes <- nodes
  ctx$ppi_edges <- edge_out
  ctx
}
