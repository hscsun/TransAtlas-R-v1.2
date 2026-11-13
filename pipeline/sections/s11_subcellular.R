# Section 11 - subcellular localization of the proteins encoded by the DEGs.
#
# Every differentially expressed gene is mapped, through its Gene Ontology
# cellular-component annotations in the species OrgDb, onto a fixed set of
# high-level compartments (nucleus, cytosol, mitochondrion, ...). The section
# then counts how many up- and down-regulated proteins belong to each
# compartment and compares that with the same count over all tested genes, so a
# compartment that is over-represented among the DEGs stands out from one that
# is merely large.
#
# The figure has two panels: a stylized cell map whose bubbles are placed at the
# compartment they describe, and a diverging bar chart of the same counts. A
# protein can carry annotations for several compartments (a mitochondrial
# ribosomal protein is both), so the compartment counts are annotations, not a
# partition of the DEG list; the table reports the background and the
# enrichment ratio next to the raw counts.


# High-level compartments, each identified by its GO cellular-component term.
# The order is the order of the bars and the legend.

GO_DB_CACHE <- new.env(parent = emptyenv())
SUBCELLULAR_COMPARTMENTS <- c(
  "Extracellular region"   = "GO:0005576",
  "Plasma membrane"        = "GO:0005886",
  "Cytosol"                = "GO:0005829",
  "Cytoplasm"              = "GO:0005737",
  "Nucleus"                = "GO:0005634",
  "Nucleoplasm"            = "GO:0005654",
  "Chromosome"             = "GO:0005694",
  "Mitochondrion"          = "GO:0005739",
  "Endoplasmic reticulum"  = "GO:0005783",
  "Golgi apparatus"        = "GO:0005794",
  "Lysosome"               = "GO:0005764",
  "Endosome"               = "GO:0005768",
  "Cytoskeleton"           = "GO:0005856",
  "Ribosome"               = "GO:0005840",
  "Vesicle"                = "GO:0031982")

# Where each compartment sits on the schematic canvas (x, y in a 0-100 x 0-74
# field); the cell outline and nucleus are drawn around these anchors.
SUBCELLULAR_COORDS <- data.frame(
  compartment = names(SUBCELLULAR_COMPARTMENTS),
  x = c(90, 62, 82, 90, 37, 50, 64, 18, 26, 46, 14, 76, 16, 70, 88),
  y = c(69, 3.5, 24, 40, 45, 31, 46, 26, 14, 15, 48, 53, 9, 60, 12),
  stringsAsFactors = FALSE)

# Terms annotated to a compartment or any of its descendants.
go_terms_in_compartment <- function(terms, compartment) {
  anc <- GO_DB_CACHE$anc
  vapply(terms, function(t) {
    t == compartment || (compartment %in% anc[[t]])
  }, logical(1))
}

# GO.db's ancestor map is large; load it once per session.
load_go_cc_ancestors <- function() {
  if (!is.null(GO_DB_CACHE$anc)) return(GO_DB_CACHE$anc)
  if (!requireNamespace("GO.db", quietly = TRUE)) stop("GO.db is required for subcellular localization.")
  GO_DB_CACHE$anc <- AnnotationDbi::as.list(GO.db::GOCCANCESTOR)
  GO_DB_CACHE$anc
}

run_subcellular <- function(cfg, ctx, out_dir) {
  out_dir <- safe_dir(out_dir)
  orgdb_name <- cfg$species$orgdb
  if (!requireNamespace(orgdb_name, quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("GO.db", quietly = TRUE)) {
    write_status(file.path(out_dir, "status.json"), "section11_subcellular", "skipped",
                 paste0("Localization needs ", orgdb_name, ", AnnotationDbi and GO.db"))
    return(ctx)
  }
  alpha <- as.numeric(cfg$analysis$alpha %||% 0.05)
  lfc <- as.numeric(cfg$analysis$lfc_threshold %||% 0.75)
  de <- ctx$de
  if (is.null(de) || !nrow(de)) {
    write_status(file.path(out_dir, "status.json"), "section11_subcellular", "skipped",
                 "No differential-expression table")
    return(ctx)
  }
  sym_col <- ctx$symbols
  # OrgDb SYMBOL keys are case-sensitive (mouse symbols are title case), so
  # the symbols are used exactly as they appear in the annotation column.
  sym_all <- as.character(de[[sym_col]])
  ok <- !is.na(sym_all) & nzchar(sym_all)
  de <- de[ok, , drop = FALSE]
  sym_all <- sym_all[ok]
  if (!length(sym_all)) {
    write_status(file.path(out_dir, "status.json"), "section11_subcellular", "skipped",
                 "No gene symbols in the differential-expression table")
    return(ctx)
  }
  # The DEG set is the one the screening section reports: adjusted P and effect
  # size, with the AnaDGEsS1 low-expression rule applied when it is switched on.
  # Using every significant gene here would double the list with low-abundance
  # genes the rest of the report does not call differentially expressed.
  sig <- !is.na(de$padj) & de$padj < alpha & abs(de$log2FoldChange) >= lfc
  expr_thr <- cfg$analysis$expression_log10_threshold
  if (!is.null(expr_thr) && is.finite(as.numeric(expr_thr))) {
    sig <- sig & log10(pmax(de$baseMean, .Machine$double.xmin)) > as.numeric(expr_thr)
  }
  deg <- data.frame(symbol = sym_all[sig],
                    direction = ifelse(de$log2FoldChange[sig] > 0, "Up", "Down"),
                    log2FC = de$log2FoldChange[sig],
                    stringsAsFactors = FALSE)
  # One row per symbol: a gene with duplicated annotation rows keeps its
  # strongest effect.
  deg <- deg[order(-abs(deg$log2FC)), , drop = FALSE]
  deg <- deg[!duplicated(deg$symbol), , drop = FALSE]
  if (!nrow(deg)) {
    write_status(file.path(out_dir, "status.json"), "section11_subcellular", "skipped",
                 "No differentially expressed genes at the reporting threshold")
    return(ctx)
  }

  orgdb <- getExportedValue(orgdb_name, orgdb_name)
  ont <- suppressWarnings(tryCatch(
    AnnotationDbi::select(orgdb, keys = unique(sym_all), keytype = "SYMBOL",
                          columns = c("GO", "ONTOLOGY")),
    error = record_optional_error))
  if (is.null(ont) || !nrow(ont)) {
    write_status(file.path(out_dir, "status.json"), "section11_subcellular", "skipped",
                 "OrgDb returned no GO annotations")
    return(ctx)
  }
  ont <- ont[!is.na(ont$GO) & ont$ONTOLOGY == "CC", c("SYMBOL", "GO"), drop = FALSE]
  ont <- unique(ont)
  if (!nrow(ont)) {
    write_status(file.path(out_dir, "status.json"), "section11_subcellular", "skipped",
                 "No cellular-component annotations")
    return(ctx)
  }

  load_go_cc_ancestors()
  terms <- unique(ont$GO)
  term_map <- data.frame(GO = terms, stringsAsFactors = FALSE, check.names = FALSE)
  for (comp in names(SUBCELLULAR_COMPARTMENTS)) {
    term_map[[comp]] <- go_terms_in_compartment(terms, SUBCELLULAR_COMPARTMENTS[[comp]])
  }
  long <- merge(ont, term_map, by = "GO")
  comp_names <- names(SUBCELLULAR_COMPARTMENTS)
  gene_comp <- stats::aggregate(long[, comp_names, drop = FALSE],
                                by = list(SYMBOL = long$SYMBOL), FUN = any)
  # Background: every tested gene that carries at least one CC annotation.
  bg_total <- nrow(gene_comp)
  bg_counts <- vapply(comp_names, function(c) sum(gene_comp[[c]]), integer(1))

  idx <- match(deg$symbol, gene_comp$SYMBOL)
  in_deg <- !is.na(idx)
  n_deg_annotated <- sum(in_deg)
  rows <- list()
  for (c in comp_names) {
    hit <- in_deg & gene_comp[[c]][idx]
    n_up <- sum(hit & deg$direction == "Up")
    n_down <- sum(hit & deg$direction == "Down")
    total <- n_up + n_down
    expected <- if (bg_total > 0) bg_counts[[c]] / bg_total * n_deg_annotated else NA_real_
    rows[[c]] <- data.frame(
      compartment = c, n_up = n_up, n_down = n_down, n_total = total,
      pct_of_degs = 100 * total / nrow(deg),
      background_total = bg_counts[[c]],
      expected = expected,
      enrichment = if (is.finite(expected) && expected > 0) total / expected else NA_real_,
      stringsAsFactors = FALSE)
  }
  loc <- do.call(rbind, rows)
  # One-sided binomial upper-tail test per compartment: the background
  # proportion N_c/N_bg is the null probability that an annotated DEG falls in
  # the compartment, and the number of trials is the annotated DEG count.
  loc$degs_tested <- n_deg_annotated
  loc$p_value <- vapply(seq_len(nrow(loc)), function(k) {
    if (!is.finite(loc$background_total[k]) || loc$background_total[k] <= 0 ||
        !is.finite(n_deg_annotated) || n_deg_annotated <= 0) return(NA_real_)
    stats::pbinom(loc$n_total[k] - 1L, size = n_deg_annotated,
                  prob = loc$background_total[k] / bg_total, lower.tail = FALSE)
  }, numeric(1))
  loc$p_adj <- stats::p.adjust(loc$p_value, method = "BH")
  loc <- loc[order(-loc$n_total), , drop = FALSE]
  safe_write_csv(loc, file.path(out_dir, paste0("subcellular_localization_", ctx$contrast_name, ".csv")))

  # Gene-level annotation: which compartments each DEG protein is annotated to.
  gene_rows <- lapply(seq_len(nrow(deg)), function(i) {
    j <- idx[i]
    if (is.na(j)) return(NULL)
    cs <- comp_names[as.logical(unlist(gene_comp[j, comp_names]))]
    if (!length(cs)) return(NULL)
    data.frame(gene = deg$symbol[i], direction = deg$direction[i], log2FC = deg$log2FC[i],
               compartments = paste(cs, collapse = "; "), stringsAsFactors = FALSE)
  })
  gene_tab <- do.call(rbind, gene_rows)
  if (!is.null(gene_tab)) {
    safe_write_csv(gene_tab, file.path(out_dir, paste0("subcellular_localization_genes_", ctx$contrast_name, ".csv")))
  }

  # --- figure ---------------------------------------------------------------
  map_df <- merge(SUBCELLULAR_COORDS, loc, by = "compartment", all.x = TRUE)
  map_df$n_up[is.na(map_df$n_up)] <- 0L
  map_df$n_down[is.na(map_df$n_down)] <- 0L
  map_df$n_total[is.na(map_df$n_total)] <- 0L
  map_df <- map_df[order(-map_df$n_total), , drop = FALSE]
  frac_up <- map_df$n_up / pmax(map_df$n_up + map_df$n_down, 1)
  map_df$direction <- ifelse(map_df$n_total == 0, "None",
                             ifelse(frac_up >= 0.65, "Up",
                                    ifelse(frac_up <= 0.35, "Down", "Balanced")))
  map_df$direction <- factor(map_df$direction, levels = c("Up", "Down", "Balanced", "None"))
  pal_dir <- c(Up = "#B2182B", Down = "#2166AC", Balanced = "grey62", None = "grey85")
  # Bubble radius in canvas units: max_size = 11 mm over a ~135 mm wide panel at
  # 100 units, so the largest bubble has r ~ 4 units.
  map_df$radius <- 4 * sqrt(map_df$n_total / max(map_df$n_total, 1))
  map_df$lab_y <- ifelse(map_df$y > 55,
                         map_df$y - (map_df$radius + 2.2),
                         map_df$y + (map_df$radius + 2.2))
  map_lab <- map_df[map_df$n_total > 0, , drop = FALSE]

  theta <- seq(0, 2 * pi, length.out = 200)
  cell_poly <- data.frame(x = c(2, 98, 98, 2), y = c(2, 2, 63, 63))
  nucleus_poly <- data.frame(x = 52 + 22 * cos(theta), y = 38 + 15 * sin(theta))

  p_map <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = cell_poly, ggplot2::aes(x, y), fill = "grey98",
                          colour = "grey45", linewidth = 0.9) +
    ggplot2::geom_polygon(data = nucleus_poly, ggplot2::aes(x, y), fill = "grey91",
                          colour = "grey55", linewidth = 0.7) +
    ggplot2::geom_point(data = map_df[map_df$n_total > 0, , drop = FALSE],
                        ggplot2::aes(x, y, size = n_total, fill = direction),
                        shape = 21, colour = "grey15", stroke = 0.6, alpha = 0.88) +
    ggplot2::scale_size_area(max_size = 11,
                             breaks = size_breaks(map_df$n_total[map_df$n_total > 0]),
                             name = "DEG proteins") +
    ggplot2::scale_fill_manual(values = pal_dir, name = "DEG direction") +
    # Labels start just outside the bubble rather than on it: the bubble radius
    # grows with the count, so the offset is computed from the same area rule
    # (radius proportional to sqrt(count)). Points near the top edge are labelled
    # downwards so they stay inside the canvas.
    ggrepel::geom_text_repel(data = map_lab,
                             ggplot2::aes(x, lab_y, label = compartment),
                             size = 3.8, colour = "grey15",
                             max.overlaps = Inf, min.segment.length = Inf,
                             box.padding = ggplot2::unit(0.7, "mm")) +
    ggplot2::coord_equal(xlim = c(0, 100), ylim = c(0, 74), expand = FALSE) +
    ggplot2::labs(title = "Subcellular localization of DEG-encoded proteins",
                  subtitle = paste0("Bubble area = DEGs annotated to the compartment; ",
                                    "a protein can map to several compartments")) +
    ggplot2::theme_void(base_family = publication_font(), base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0),
                   plot.subtitle = ggplot2::element_text(size = 12.5, colour = "grey25", hjust = 0),
                   legend.position = "right",
                   legend.text = ggplot2::element_text(size = 12.5),
                   legend.title = ggplot2::element_text(size = 13.5, face = "bold"),
                   plot.margin = ggplot2::margin(6, 6, 6, 6))

  bar_df <- loc
  bar_df$compartment <- factor(bar_df$compartment, levels = rev(loc$compartment))
  p_bar <- ggplot2::ggplot(bar_df, ggplot2::aes(y = compartment)) +
    ggplot2::geom_col(ggplot2::aes(x = n_up, fill = "Up"), width = 0.72) +
    ggplot2::geom_col(ggplot2::aes(x = -n_down, fill = "Down"), width = 0.72) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(x = ifelse(n_up >= n_down, n_up, -n_down),
                                    label = n_total,
                                    hjust = ifelse(n_up >= n_down, -0.15, 1.15)),
                       size = 3.8, colour = "grey20") +
    # The direction legend is already carried by the cell map next to it.
    ggplot2::scale_fill_manual(values = pal_dir[c("Up", "Down")], guide = "none") +
    ggplot2::scale_x_continuous(labels = abs,
                                expand = ggplot2::expansion(mult = c(0.18, 0.18))) +
    ggplot2::labs(title = "DEGs per compartment", x = "Number of DEG-encoded proteins", y = NULL) +
    publication_theme(base_size = 13) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 12.5),
                   axis.text.x = ggplot2::element_text(size = 12),
                   legend.text = ggplot2::element_text(size = 12.5),
                   legend.title = ggplot2::element_text(size = 13.5, face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold", size = 16),
                   panel.grid.major.y = ggplot2::element_blank())

  fig_path <- file.path(out_dir, paste0("subcellular_localization_", ctx$contrast_name, ".png"))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- p_map + p_bar + patchwork::plot_layout(widths = c(1.25, 1))
    save_plot(combined, fig_path, key = "s11_subcellular")
  } else {
    # Without patchwork the two panels are written separately rather than
    # pretending a combined figure exists.
    save_plot(p_map, file.path(out_dir, paste0("subcellular_map_", ctx$contrast_name, ".png")),
              key = "s11_subcellular_map")
    save_plot(p_bar, file.path(out_dir, paste0("subcellular_bars_", ctx$contrast_name, ".png")),
              key = "s11_subcellular_bars")
  }

  ctx$localization <- loc
  write_status(file.path(out_dir, "status.json"), "section11_subcellular", "completed", NULL,
               list(degs = nrow(deg), annotated = sum(in_deg), background = bg_total,
                    compartments = nrow(loc), figure = if (requireNamespace("patchwork", quietly = TRUE))
                      basename(fig_path) else c("subcellular_map.png", "subcellular_bars.png")))
  ctx
}

