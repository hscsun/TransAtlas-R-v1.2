# Section 7 - TF-target regulatory network.
#
# Port of input/NetForTFs.R. The original script built an edge list around a
# chosen TF, coloured the nodes by log2 fold change and drew the network in
# several ggraph layouts. This section keeps that idiom - ggraph, an edge colour
# for the mode of regulation, a green-white-purple fold-change ramp on the nodes,
# and repelled labels - and feeds it every TF the regulon analysis selected
# instead of one hard-coded factor.

run_tf_network <- function(cfg, ctx, out_dir) {
  out_dir <- safe_dir(out_dir)
  if (is.null(ctx$tf_summary) || is.null(ctx$tf_edges)) {
    write_status(file.path(out_dir, "status.json"), "section07_tf_network", "skipped",
                 "Run section 5 first: no TF summary/edge table available")
    return(ctx)
  }
  if (!requireNamespace("ggraph", quietly = TRUE) || !requireNamespace("tidygraph", quietly = TRUE)) {
    write_status(file.path(out_dir, "status.json"), "section07_tf_network", "skipped",
                 "ggraph/tidygraph unavailable")
    return(ctx)
  }
  n_tfs <- as.integer(cfg$regulatory$network_tfs %||% 25L)
  max_labels <- as.integer(cfg$regulatory$network_max_labels %||% 40L)
  layouts <- unlist(cfg$regulatory$network_layouts %||% c("fr", "kk", "circle", "tree"))

  edges_all <- ctx$tf_edges
  sum_df <- ctx$tf_summary
  # The network draws the same expressed factors the dot plots show: a TF whose
  # own FPKM is below the floor is not a credible regulator to put at the centre
  # of a figure, even when its target overlap is significant. `expressed` arrives
  # as a character column when the summary is reloaded from disk.
  if (!is.null(sum_df$expressed)) {
    keep <- as.logical(sum_df$expressed)
    keep[is.na(keep)] <- TRUE
    sum_df <- sum_df[keep, , drop = FALSE]
  }
  ranking <- if (!is.null(sum_df$NES) && any(is.finite(sum_df$NES))) {
    sum_df$tf_key[order(-abs(ifelse(is.finite(sum_df$NES), sum_df$NES, 0)))]
  } else {
    sum_df$tf_key[order(-sum_df$targets_deg)]
  }
  picked <- utils::head(ranking, n_tfs)
  edges <- edges_all[edges_all$tf_key %in% picked & edges_all$target_is_deg, , drop = FALSE]
  if (nrow(edges) < 5L) {
    write_status(file.path(out_dir, "status.json"), "section07_tf_network", "skipped",
                 "Too few DE TF-target edges to draw a network")
    return(ctx)
  }
  edges <- edges[order(-abs(edges$target_log2FC)), , drop = FALSE]
  edges <- edges[!duplicated(edges[, c("tf_key", "target_key")]), , drop = FALSE]
  # A factor left with only one or two displayed targets cannot show a regulon
  # in a network: its node would be an isolated dot with a couple of spokes.
  # Such factors are removed here (and from the node table), so every TF drawn
  # carries at least three differentially expressed targets.
  target_counts <- table(edges$tf_key)
  keep_tfs <- names(target_counts)[target_counts >= 3L]
  edges <- edges[edges$tf_key %in% keep_tfs, , drop = FALSE]
  picked <- intersect(picked, keep_tfs)
  if (nrow(edges) < 5L) {
    write_status(file.path(out_dir, "status.json"), "section07_tf_network", "skipped",
                 "Too few DE TF-target edges after dropping TFs with fewer than 3 targets")
    return(ctx)
  }

  node_names <- unique(c(edges$tf_key, edges$target_key))
  nd <- data.frame(name = node_names,
                   is_tf = node_names %in% edges$tf_key,
                   lfc = edges$target_log2FC[match(node_names, edges$target_key)],
                   stringsAsFactors = FALSE)
  # Retain regulon activity and the TF's own effect in the exported node table.
  # Display sizes distinguish TFs from targets.
  nd$nes <- if (!is.null(sum_df$NES)) sum_df$NES[match(node_names, sum_df$tf_key)] else NA_real_
  nd$tf_lfc <- if (!is.null(sum_df$tf_log2FC)) sum_df$tf_log2FC[match(node_names, sum_df$tf_key)] else NA_real_
  # The TF's own mean FPKM travels with the node table and is written into the
  # label, so the expression level the figure was filtered on is visible.
  nd$tf_fpkm <- if (!is.null(sum_df$tf_fpkm)) as.numeric(sum_df$tf_fpkm[match(node_names, sum_df$tf_key)]) else NA_real_
  # A TF node that is not a target in this subnetwork still gets its own
  # expression, so no node is left without a value.
  fill <- nd$is_tf & is.na(nd$lfc) & !is.na(nd$tf_lfc)
  nd$lfc[fill] <- nd$tf_lfc[fill]
  # Missing effects remain NA in the export and use the missing-value colour.
  safe_write_csv(edges, file.path(out_dir, paste0("TF_network_edges_", ctx$contrast_name, ".csv")))
  safe_write_csv(nd, file.path(out_dir, paste0("TF_network_nodes_", ctx$contrast_name, ".csv")))

  # Edge direction of regulation comes from the resource's mode-of-regulation
  # column when present; otherwise it follows the target's own DE direction.
  edge_reg <- if (!is.null(edges$mor) && any(!is.na(edges$mor))) {
    ifelse(edges$mor >= 0, "activating", "repressing")
  } else {
    ifelse(edges$target_direction == "up", "activating", "repressing")
  }

  render_layout <- function(layout_name, path) {
    nd$size <- ifelse(nd$is_tf, 7, 2.6)
    # TF labels carry the name only; FPKM stays in the node table.
    nd$tf_label <- ifelse(nd$is_tf, nd$name, NA_character_)
    # Label the 15 most up-regulated and 15 most down-regulated targets (by the
    # effect size, uncapped on the node table). Labelling all of them buries the
    # regulons under hundreds of overlapping names; ranking only by |log2FC|
    # without splitting the sign lets one direction crowd out the other.
    n_per_direction <- max(1L, min(15L, as.integer(max_labels / 2)))
    tg <- which(!nd$is_tf & is.finite(nd$lfc))
    tg_up <- tg[nd$lfc[tg] > 0]
    tg_down <- tg[nd$lfc[tg] < 0]
    lab <- c(utils::head(tg_up[order(-nd$lfc[tg_up])], n_per_direction),
             utils::head(tg_down[order(nd$lfc[tg_down])], n_per_direction))
    nd$target_label <- NA_character_
    nd$target_label[lab] <- nd$name[lab]

    g <- tidygraph::tbl_graph(nodes = nd,
                              edges = data.frame(from = edges$tf_key, to = edges$target_key,
                                                 regulation = edge_reg, stringsAsFactors = FALSE),
                              directed = TRUE)
    # A ring drawn in vertex order puts every TF in one arc - they come first in
    # the node table - which stacks all the labels on top of each other. Spreading
    # the TFs evenly between the targets gives each one its own sector.
    lay <- if (identical(layout_name, "circle")) {
      n <- nrow(nd)
      tf_i <- which(nd$is_tf)
      tg_i <- which(!nd$is_tf)
      ord <- integer(0)
      if (length(tf_i) && length(tg_i)) {
        for (i in seq_along(tf_i)) {
          s <- floor((i - 1) * length(tg_i) / length(tf_i)) + 1L
          e <- floor(i * length(tg_i) / length(tf_i))
          if (e >= s) ord <- c(ord, tg_i[s:e])
          ord <- c(ord, tf_i[i])
        }
        ord <- c(ord, setdiff(seq_len(n), ord))
      } else {
        ord <- seq_len(n)
      }
      ang <- 2 * pi * (seq_len(n) - 1L) / n
      m <- cbind(cos(ang), sin(ang))
      m[ord, ] <- m
      m
    } else {
      layout_name
    }
    # Colour clipping affects this display only; node effects remain true values.
    lim <- display_lfc_limit(cfg, nd$lfc)
    p <- ggraph::ggraph(g, layout = lay) +
      # Saturated red / blue at half opacity: the previous vermillion and dark
      # blue at 0.28 alpha read as two shades of grey once a few hundred edges
      # overlapped, so the mode of regulation - the point of the edge colour -
      # was invisible.
      ggraph::geom_edge_link0(ggplot2::aes(colour = regulation), alpha = 0.55, width = 0.45) +
      ggraph::scale_edge_colour_manual(values = c(activating = "#E41A1C", repressing = "#1B6FC4"),
                                       name = "Regulation",
                                       guide = ggplot2::guide_legend(
                                         override.aes = list(edge_width = 1.8, edge_alpha = 1))) +
      ggraph::geom_node_point(ggplot2::aes(colour = lfc, size = size, shape = is_tf)) +
      ggplot2::scale_shape_manual(values = c(`TRUE` = 15, `FALSE` = 19), name = "Node",
                                  labels = c(`TRUE` = "TF (square)", `FALSE` = "target (circle)")) +
      ggplot2::scale_size_identity() +
      ggplot2::scale_colour_gradient2(low = "green4", mid = "grey88", high = "purple3", midpoint = 0,
                                      limits = c(-lim, lim), oob = scales::squish, na.value = "grey65",
                                      name = expression("log"[2] * " fold change")) +
      ggraph::geom_node_text(ggplot2::aes(label = tf_label), repel = TRUE, size = 4.65,
                             fontface = "bold", colour = "grey5", max.overlaps = Inf,
                             na.rm = TRUE,
                             point.padding = ggplot2::unit(1.6, "mm")) +
      # na.rm = TRUE: only a fraction of the targets is labelled, and the
      # empty labels would otherwise be reported as "Removed N rows" on every
      # layout.
      ggraph::geom_node_text(ggplot2::aes(label = target_label), repel = TRUE, size = 3.85,
                             colour = "grey20", max.overlaps = 40, na.rm = TRUE,
                             point.padding = ggplot2::unit(1.0, "mm")) +
      ggraph::theme_graph(base_family = publication_font()) +
      ggplot2::theme(legend.text = ggplot2::element_text(size = 15),
                     legend.title = ggplot2::element_text(size = 16, face = "bold"),
                     legend.key.size = ggplot2::unit(1.0, "cm"),
                     legend.spacing.y = ggplot2::unit(0.25, "cm"),
                     legend.box.spacing = ggplot2::unit(0.5, "cm"),
                     plot.title = ggplot2::element_text(face = "bold", size = 16),
                     plot.subtitle = ggplot2::element_text(size = 13)) +
      ggplot2::labs(title = paste0("TF-target regulatory network (", layout_name, " layout): ",
                                   ctx$contrast_name),
                    subtitle = paste0("TF squares and target circles have fixed display sizes; the ",
                                      n_per_direction, " most up- and ", n_per_direction,
                                      " most down-regulated targets are labelled"))
    save_plot(p, path, key = "s7_tf_network")
  }

  produced <- character()
  for (lay in layouts) {
    path <- file.path(out_dir, sprintf("TF_network_%s_%s.png", lay, ctx$contrast_name))
    ok <- tryCatch({ render_layout(lay, path); TRUE }, error = function(e) conditionMessage(e))
    if (isTRUE(ok)) produced <- c(produced, basename(path)) else
      writeLines(as.character(ok), file.path(out_dir, paste0("network_error_", lay, ".txt")))
  }

  ctx$tf_network <- list(edges = edges, nodes = nd, layouts = produced)
  # A section that renders nothing must not report "completed": the figure files
  # of a previous run would still be lying in the output directory and look like
  # this run's output.
  if (!length(produced)) {
    errs <- list.files(out_dir, pattern = "^network_error_", full.names = TRUE)
    first <- if (length(errs)) readLines(errs[1L], warn = FALSE)[1L] else "unknown error"
    write_status(file.path(out_dir, "status.json"), "section07_tf_network", "failed",
                 paste0("No layout rendered: ", first),
                 list(tfs = length(picked), edges = nrow(edges), nodes = nrow(nd),
                      layouts_rendered = produced))
    return(ctx)
  }
  unlink(list.files(out_dir, pattern = "^network_error_", full.names = TRUE))
  write_status(file.path(out_dir, "status.json"), "section07_tf_network", "completed", NULL,
               list(tfs = length(picked), edges = nrow(edges), nodes = nrow(nd),
                    layouts_rendered = produced))
  ctx
}
