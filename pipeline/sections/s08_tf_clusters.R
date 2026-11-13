# Section 8 - TF modules: cluster the screened transcription factors by the
# overlap of their differentially expressed target sets, then draw the module
# network.
#
# The earlier sections answer "which regulons are shifted" one factor at a time.
# This section asks the next question: do the shifted factors work in coherent
# groups? Two TFs whose DE target sets overlap heavily are likely to act on the
# same programme, so their target sets are compared with the Jaccard index and
# the resulting distance matrix is cut into modules (average linkage, number of
# modules chosen by mean silhouette width).
#
# The figure encodes three quantities at once:
#   node size   = number of differentially expressed targets of that TF
#   node colour = module membership
#   edge width  = Jaccard overlap of the two target sets
#
# It is not a physical interaction map: a shared target can be regulated by both
# factors independently, so the modules are hypotheses about co-regulated
# programmes, not evidence of co-binding.

run_tf_modules <- function(cfg, ctx, out_dir) {
  out_dir <- safe_dir(out_dir)
  if (is.null(ctx$tf_summary) || is.null(ctx$tf_edges)) {
    write_status(file.path(out_dir, "status.json"), "section08_tf_modules", "skipped",
                 "Run section 5 first: no TF summary/edge table available")
    return(ctx)
  }
  if (!requireNamespace("ggraph", quietly = TRUE) ||
      !requireNamespace("tidygraph", quietly = TRUE) ||
      !requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("ggrepel", quietly = TRUE) ||
      !requireNamespace("cluster", quietly = TRUE)) {
    write_status(file.path(out_dir, "status.json"), "section08_tf_modules", "skipped",
                 "ggraph/tidygraph/igraph/ggrepel/cluster unavailable")
    return(ctx)
  }

  n_tfs <- as.integer(cfg$regulatory$cluster_tfs %||% cfg$regulatory$tf_dotplot_max %||% 30L)
  min_jaccard <- as.numeric(cfg$regulatory$cluster_min_jaccard %||% 0.1)
  max_modules <- as.integer(cfg$regulatory$cluster_max %||% 8L)
  max_edges <- as.integer(cfg$regulatory$cluster_max_edges %||% 120L)
  min_targets <- max(2L, as.integer(cfg$regulatory$ora_min_overlap %||% 3L))

  sum_df <- ctx$tf_summary
  edges_all <- ctx$tf_edges
  # Only factors that pass the expression floor are clustered, matching the dot
  # plot and the target network; the summary keeps the full statistics.
  if (!is.null(sum_df$expressed)) {
    keep_tf <- as.logical(sum_df$expressed)
    keep_tf[is.na(keep_tf)] <- TRUE
    sum_df <- sum_df[keep_tf, , drop = FALSE]
  }
  ranking <- if (!is.null(sum_df$NES) && any(is.finite(sum_df$NES))) {
    sum_df$tf_key[order(-abs(ifelse(is.finite(sum_df$NES), sum_df$NES, 0)))]
  } else {
    sum_df$tf_key[order(-sum_df$targets_deg)]
  }
  picked <- utils::head(ranking, n_tfs)

  need <- c("tf_key", "target_key", "target_is_deg")
  if (!all(need %in% names(edges_all))) {
    write_status(file.path(out_dir, "status.json"), "section08_tf_modules", "skipped",
                 "TF edge table lacks the columns this section needs")
    return(ctx)
  }
  edges <- edges_all[edges_all$tf_key %in% picked & isTRUE_vec(edges_all$target_is_deg), , drop = FALSE]
  if (!nrow(edges)) {
    write_status(file.path(out_dir, "status.json"), "section08_tf_modules", "skipped",
                 "No differentially expressed TF-target edges")
    return(ctx)
  }
  edges <- edges[!duplicated(edges[, c("tf_key", "target_key")]), , drop = FALSE]
  sets <- split(edges$target_key, edges$tf_key)
  sets <- sets[lengths(sets) >= min_targets]
  if (length(sets) < 3L) {
    write_status(file.path(out_dir, "status.json"), "section08_tf_modules", "skipped",
                 paste0("Fewer than three TFs with >= ", min_targets, " DE targets"))
    return(ctx)
  }
  tfs <- names(sets)
  # Keep the strongest TFs if the target-size floor let more than n_tfs through.
  tfs <- utils::head(tfs[order(-lengths(sets))], n_tfs)
  sets <- sets[tfs]

  # --- Jaccard similarity ---------------------------------------------------
  n <- length(tfs)
  sim <- matrix(1, n, n, dimnames = list(tfs, tfs))
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      a <- sets[[i]]
      b <- sets[[j]]
      inter <- length(intersect(a, b))
      uni <- length(union(a, b))
      sim[i, j] <- sim[j, i] <- if (uni > 0) inter / uni else 0
    }
  }

  # --- module cut -----------------------------------------------------------
  # Average-linkage clustering on 1 - Jaccard, with the number of modules chosen
  # by the mean silhouette width rather than fixed in advance. A single module
  # is the fallback when no split improves on it.
  dist_mat <- stats::as.dist(1 - sim)
  hc <- stats::hclust(dist_mat, method = "average")
  best <- list(k = 1L, label = rep(1L, n), sil = NA_real_)
  if (n >= 4L && max_modules >= 2L) {
    upper <- min(max_modules, n - 1L)
    for (k in 2:upper) {
      lab <- stats::cutree(hc, k = k)
      if (length(unique(lab)) < 2L) next
      sil <- tryCatch(mean(cluster::silhouette(lab, dist_mat)[, 3]), error = function(e) NA_real_)
      if (is.finite(sil) && sil > 0 && (is.na(best$sil) || sil > best$sil)) {
        best <- list(k = as.integer(k), label = lab, sil = sil)
      }
    }
  }
  raw_cluster <- best$label
  # Modules are renamed by size so M1 is always the largest, and the legend
  # reads in the same order as the figure.
  sizes <- sort(table(raw_cluster), decreasing = TRUE)
  remap <- stats::setNames(seq_along(sizes), names(sizes))
  module <- paste0("M", remap[as.character(raw_cluster)])
  module <- factor(module, levels = paste0("M", seq_along(sizes)))
  names(module) <- tfs

  nodes <- data.frame(
    name = tfs,
    module = module[tfs],
    n_targets = as.integer(lengths(sets)),
    n_up = vapply(sets, function(s) sum(edges$target_log2FC[match(s, edges$target_key)] > 0), integer(1)),
    n_down = vapply(sets, function(s) sum(edges$target_log2FC[match(s, edges$target_key)] < 0), integer(1)),
    mean_target_log2FC = vapply(sets, function(s) {
      v <- edges$target_log2FC[match(s, edges$target_key)]
      mean(v, na.rm = TRUE)
    }, numeric(1)),
    nes = if (!is.null(sum_df$NES)) sum_df$NES[match(tfs, sum_df$tf_key)] else NA_real_,
    fpkm = if (!is.null(sum_df$tf_fpkm)) as.numeric(sum_df$tf_fpkm[match(tfs, sum_df$tf_key)]) else NA_real_,
    top_GO_BP = if (!is.null(sum_df$top_GO_BP)) as.character(sum_df$top_GO_BP[match(tfs, sum_df$tf_key)]) else NA_character_,
    stringsAsFactors = FALSE)
  nodes <- nodes[order(nodes$module, -nodes$n_targets), , drop = FALSE]

  # --- edges ----------------------------------------------------------------
  # Every pair above the Jaccard floor is a candidate; the highest overlaps are
  # kept so the figure stays readable, and any TF left without an edge is
  # attached to its nearest neighbour so no module floats unconnected.
  pair_i <- integer(0); pair_j <- integer(0)
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      if (sim[i, j] >= min_jaccard) { pair_i <- c(pair_i, i); pair_j <- c(pair_j, j) }
    }
  }
  links <- data.frame(from = tfs[pair_i], to = tfs[pair_j],
                      jaccard = sim[cbind(pair_i, pair_j)], stringsAsFactors = FALSE)
  links <- links[order(-links$jaccard), , drop = FALSE]
  if (nrow(links) > max_edges) links <- links[seq_len(max_edges), , drop = FALSE]
  linked <- unique(c(links$from, links$to))
  for (tf in setdiff(tfs, linked)) {
    j <- setdiff(tfs, tf)
    k <- j[which.max(sim[tf, j])]
    links <- rbind(links, data.frame(from = tf, to = k, jaccard = sim[tf, k],
                                     stringsAsFactors = FALSE))
  }
  links$jaccard[!is.finite(links$jaccard)] <- 0

  # --- tables ---------------------------------------------------------------
  membership <- nodes[, c("name", "module", "n_targets", "n_up", "n_down",
                          "mean_target_log2FC", "nes", "fpkm", "top_GO_BP")]
  names(membership)[1] <- "tf"
  safe_write_csv(membership, file.path(out_dir, paste0("TF_cluster_membership_", ctx$contrast_name, ".csv")))
  sim_out <- data.frame(tf = rownames(sim), sim, check.names = FALSE)
  safe_write_csv(sim_out, file.path(out_dir, paste0("TF_cluster_similarity_", ctx$contrast_name, ".csv")))
  safe_write_csv(links, file.path(out_dir, paste0("TF_cluster_links_", ctx$contrast_name, ".csv")))

  module_rows <- lapply(levels(module), function(m) {
    idx <- which(nodes$module == m)
    union_targets <- unique(unlist(sets[nodes$name[idx]]))
    rep_tf <- nodes$name[idx][which.max(nodes$n_targets[idx])]
    data.frame(module = m,
               n_tfs = length(idx),
               tfs = paste(nodes$name[idx], collapse = ", "),
               n_targets_union = length(union_targets),
               n_up_union = sum(edges$target_log2FC[match(union_targets, edges$target_key)] > 0, na.rm = TRUE),
               n_down_union = sum(edges$target_log2FC[match(union_targets, edges$target_key)] < 0, na.rm = TRUE),
               representative = rep_tf,
               representative_top_GO_BP = nodes$top_GO_BP[idx][which.max(nodes$n_targets[idx])],
               stringsAsFactors = FALSE)
  })
  module_df <- do.call(rbind, module_rows)
  safe_write_csv(module_df, file.path(out_dir, paste0("TF_cluster_modules_", ctx$contrast_name, ".csv")))

  # --- figure ---------------------------------------------------------------
  g <- tidygraph::tbl_graph(nodes = nodes,
                            edges = links, directed = FALSE)
  pal <- if (nlevels(module) <= 8L) {
    c("#E41A1C", "#1B6FC4", "#4DAF4A", "#984EA3", "#FF7F00",
      "#00A6A6", "#A65628", "#F781BF")[seq_len(nlevels(module))]
  } else {
    grDevices::hcl.colors(nlevels(module), "Dark 3")
  }
  names(pal) <- levels(module)
  # No explicit breaks/limits on the edge-width scale: ggraph draws the edge
  # key glyph inside the node-size keys as soon as that scale is given breaks,
  # which puts a stray line through every circle of the size legend. The default
  # breaks already span the observed range.
  node_breaks <- size_breaks(nodes$n_targets)

  # The layout is computed first so the labels can be placed with ggrepel
  # directly. ggraph 2.2.2 passes repel = TRUE to ggrepel but does not forward
  # the point padding, so labels were printed on top of the large nodes.
  set.seed(as.integer(cfg$project$seed %||% 1L))
  lay <- ggraph::create_layout(g, layout = "fr")

  p <- ggraph::ggraph(lay) +
    # The edge-width legend is suppressed: ggraph 2.2.2 draws its key glyphs
    # inside the node-size keys whenever both legends are present (a line
    # through every size circle). The width still encodes Jaccard on the panel,
    # and the subtitle states the mapping and its range.
    ggraph::geom_edge_link(ggplot2::aes(width = jaccard, alpha = jaccard),
                           colour = "grey60", show.legend = FALSE) +
    ggraph::scale_edge_width_continuous(range = c(0.25, 1.6), guide = "none") +
    ggraph::scale_edge_alpha_continuous(range = c(0.25, 0.8), guide = "none") +
    ggraph::geom_node_point(ggplot2::aes(size = n_targets, fill = module),
                            shape = 21, colour = "grey15", stroke = 0.7) +
    ggplot2::scale_size_continuous(range = c(2.5, 9), breaks = node_breaks,
                                   name = "DE targets") +
    ggplot2::scale_fill_manual(values = pal, name = "TF module",
                               guide = ggplot2::guide_legend(
                                 override.aes = list(size = 5, stroke = 0.7))) +
    # Labels sit outside the bubbles with a short leader line. 6 mm clears the
    # largest node (9 mm diameter) and keeps the text off the dark fills.
    ggrepel::geom_text_repel(data = lay,
                             ggplot2::aes(x, y, label = name),
                             size = 4.25, fontface = "bold", colour = "grey5",
                             max.overlaps = Inf,
                             point.padding = ggplot2::unit(6, "mm"),
                             box.padding = ggplot2::unit(1.3, "mm"),
                             min.segment.length = 0.3,
                             segment.size = 0.3, segment.colour = "grey45",
                             force = 2, max.iter = 20000) +
    ggraph::theme_graph(base_family = publication_font()) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 11.5),
                   legend.title = ggplot2::element_text(size = 12.5, face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold", size = 15),
                   plot.subtitle = ggplot2::element_text(size = 11.5, colour = "grey25")) +
    ggplot2::labs(
      title = paste0("TF modules by shared DE targets: ", ctx$contrast_name),
      subtitle = paste0("Node size = DE targets regulated; colour = module (",
                        nlevels(module), " modules, mean silhouette ",
                        if (is.na(best$sil)) "n/a" else sprintf("%.2f", best$sil),
                        "); edge width = Jaccard overlap (",
                        sprintf("%.2f-%.2f", min(links$jaccard), max(links$jaccard)),
                        ")"))
  save_plot(p, file.path(out_dir, paste0("TF_cluster_network_", ctx$contrast_name, ".png")),
            key = "s8_tf_modules")

  ctx$tf_modules <- list(nodes = nodes, links = links, modules = module_df,
                         similarity = sim, silhouette = best$sil)
  write_status(file.path(out_dir, "status.json"), "section08_tf_modules", "completed", NULL,
               list(tfs = length(tfs), modules = nlevels(module),
                    edges = nrow(links), mean_silhouette = best$sil,
                    min_jaccard = min_jaccard))
  ctx
}

# as.logical() on a column that arrived from CSV can be character ("TRUE").
isTRUE_vec <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  as.logical(x) %in% TRUE
}
