# Section 4 - KEGG pathway analysis and pathway maps.
#
# Port of input/AnaDGEsS4.R / input/AnaDGEsS4Mm.R. The original script ran
# gseKEGG on the ranked Entrez list, wrote a KEGG table, looped pathview over
# every significant pathway, and drew a lollipop plot plus a multi-pathway
# enrichment panel and an UpSet plot. All of that is reproduced; the species
# code, the pathway cap and the identifier vocabulary now come from the config.

run_kegg <- function(cfg, ctx, out_dir) {
  out_dir <- safe_dir(out_dir)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    write_status(file.path(out_dir, "status.json"), "section04_kegg_pathview", "skipped", "clusterProfiler unavailable")
    return(ctx)
  }
  rank <- ctx$gene_list
  if (is.null(rank) || length(rank) < 50L) {
    write_status(file.path(out_dir, "status.json"), "section04_kegg_pathview", "skipped", "Ranked gene list unavailable")
    return(ctx)
  }
  kegg_id_type <- kegg_key_type(cfg, ctx$key_type)
  alpha <- as.numeric(cfg$analysis$gsea_pvalue_cutoff %||% 0.05)
  min_size <- as.numeric(cfg$analysis$gsea_min_size %||% 15)
  max_size <- as.numeric(cfg$analysis$gsea_max_size %||% 500)

  kk <- run_kegg_gsea(cfg, rank, kegg_id_type, alpha, min_size, max_size)
  if (is.null(kk) || nrow(as.data.frame(kk)) == 0L) {
    write_status(file.path(out_dir, "status.json"), "section04_kegg_pathview", "failed",
                 paste0("No KEGG gene sets available: the cached tables under ",
                        kegg_cache_dir(cfg), " are absent and rest.kegg.jp could not ",
                        "deliver the organism's pathway table. Re-run with network ",
                        "access to populate the cache."))
    return(ctx)
  }
  kegg_source <- attr(kk, "kegg_source") %||% "unknown"
  writeLines(kegg_source, file.path(out_dir, paste0("KEGG_source_", ctx$contrast_name, ".txt")))
  message("    KEGG gene sets: ", kegg_source)
  kk_df <- as.data.frame(kk)
  if ("core_enrichment" %in% names(kk_df)) {
    kk_df$core_enrichment <- ids_to_symbols(
      kk_df$core_enrichment,
      key_to_symbol_map(ctx$gene_map, ctx$key_type, ctx$symbol_type))
  }
  safe_write_csv(kk_df, file.path(out_dir, paste0("KEGG_", ctx$contrast_name, ".csv")))
  sig <- kk_df[kk_df$p.adjust < alpha, , drop = FALSE]

  # --- lollipop plot --------------------------------------------------------
  top <- utils::head(kk_df[order(kk_df$p.adjust), , drop = FALSE], 20L)
  top$Description <- factor(top$Description, levels = rev(top$Description))
  lollipop <- ggplot2::ggplot(top, ggplot2::aes(enrichmentScore, Description)) +
    ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = Description), colour = "grey55") +
    ggplot2::geom_point(ggplot2::aes(colour = p.adjust, size = setSize)) +
    ggplot2::scale_colour_gradient(low = "purple", high = "green", trans = "log10", name = "adjusted P") +
    ggplot2::scale_size_continuous(range = c(2, 8), breaks = size_breaks(top$setSize),
                                   name = "Set size") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::labs(x = "Enrichment score", y = NULL,
                  title = wrap_label(paste("KEGG pathway enrichment:", ctx$contrast_name),
                                     width = 26L)) +
    publication_theme() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 11),
                                          axis.text.x = ggplot2::element_text(size = 11),
                                          legend.text = ggplot2::element_text(size = 11))
  save_plot(lollipop, file.path(out_dir, paste0("KEGG_lollipop_", ctx$contrast_name, ".png")),
            key = "s4_kegg_lollipop")

  # --- multi-pathway GSEA panel --------------------------------------------
  # Half the slots go to the induced pathways and half to the repressed ones, so
  # the panel shows both arms of the response rather than the strongest one
  # twice over; a pathway sharing more than 40% of its genes with one already
  # chosen is skipped, and the ▲/▼ prefix in the legend carries the direction.
  n_panel <- min(as.integer(cfg$analysis$gsea_top_sets %||% 9L), nrow(kk_df))
  if (n_panel >= 2L) {
    sel <- balanced_set_selection(kk_df, n_panel,
                                  gene_sets = tryCatch(kk@geneSets, error = record_optional_error),
                                  max_overlap = 0.4)
    kk_lab <- kk
    kk_lab@result$Description <- paste0(ifelse(kk_lab@result$NES > 0, "▲ ", "▼ "),
                                        kk_lab@result$Description)
    palette <- c("#f7ca64", "#43a5bf", "#86c697", "#a670d6", "#ef998a",
                 "#0000FF", "#CD00CD", "#C0FF3E", "#EE0000")
    p <- tryCatch(enrichplot::gseaplot2(kk_lab, geneSetID = sel, pvalue_table = FALSE, base_size = 14),
                  error = record_optional_error)
    if (!is.null(p)) {
      cols <- rep(palette, length.out = length(sel))
      p[[1]] <- p[[1]] + ggplot2::scale_colour_manual(values = cols)
      p[[2]] <- p[[2]] + ggplot2::scale_colour_manual(values = cols)
      sz <- fig_size("s4_kegg_composite")
      record_figure_size("s4_kegg_composite", sz, paste0("KEGG_gsea_panel_", ctx$contrast_name, ".png"))
      grDevices::png(file.path(out_dir, paste0("KEGG_gsea_panel_", ctx$contrast_name, ".png")),
                     width = sz[1] * 300, height = sz[2] * 300, res = 300)
                     plot_device <- grDevices::dev.cur()
                     do.call(on.exit, list(substitute(close_plot_device(device), list(device = plot_device)), add = TRUE))
      print(p)
      close_plot_device(plot_device)
    }
  }

  # --- UpSet plot of pathway gene-set overlap ------------------------------
  # The original AnaDGEsS4.R used enrichplot::upsetplot(), which needs the
  # ggupset package. ggupset cannot be installed on this host, so the same
  # overlap matrix is drawn with UpSetR from the same KEGG gene sets.
  upset_path <- file.path(out_dir, paste0("KEGG_upset_", ctx$contrast_name, ".png"))
  drawn <- FALSE
  if (requireNamespace("enrichplot", quietly = TRUE) &&
      requireNamespace("ggupset", quietly = TRUE)) {
    p <- tryCatch(enrichplot::upsetplot(kk, 20), error = record_optional_error)
    if (!is.null(p)) {
      sz <- fig_size("s4_kegg_upset")
      record_figure_size("s4_kegg_upset", sz, upset_path)
      grDevices::png(upset_path, width = sz[1] * 300, height = sz[2] * 300, res = 300)
      plot_device <- grDevices::dev.cur()
      do.call(on.exit, list(substitute(close_plot_device(device), list(device = plot_device)), add = TRUE))
      print(p)
      close_plot_device(plot_device)
      drawn <- TRUE
    }
  }
  if (!drawn && requireNamespace("UpSetR", quietly = TRUE)) {
    sets <- tryCatch({
      # Twelve pathways keeps the row labels legible; twenty overlap.
      keep <- utils::head(kk@result$ID[order(kk@result$p.adjust)], 12)
      gs <- kk@geneSets[keep]
      gs <- lapply(gs, function(x) unique(as.character(x)))
      gs[!vapply(gs, is.null, logical(1))]
    }, error = record_optional_error)
    if (!is.null(sets) && length(sets) >= 2L) {
      sz <- fig_size("s4_kegg_upset_full")
      record_figure_size("s4_kegg_upset_full", sz, upset_path)
      grDevices::png(upset_path, width = sz[1] * 300, height = sz[2] * 300, res = 300)
      plot_device <- grDevices::dev.cur()
      do.call(on.exit, list(substitute(close_plot_device(device), list(device = plot_device)), add = TRUE))
      graphics::par(family = publication_font())
      # UpSetR draws through its print method, so the object must be printed.
      print(UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), order.by = "freq",
                          mainbar.y.label = "Genes in intersection",
                          sets.x.label = "Genes per KEGG pathway",
                          text.scale = c(1.1, 1.1, 1.0, 1.0, 1.5, 1.1)))
      close_plot_device(plot_device)
      drawn <- TRUE
    }
  }
  if (!drawn) {
    writeLines("Neither enrichplot+ggupset nor UpSetR could draw the overlap matrix.",
               file.path(out_dir, "KEGG_upset_note.txt"))
  }

  # --- pathview maps --------------------------------------------------------
  rendered <- 0L
  errors <- character()
  # pathview needs a real KEGG pathway id (mmu04110) because it fetches that
  # pathway's KGML and native map. The MSigDB fallback has no such ids, so the
  # maps are skipped rather than requested with an id KEGG does not know.
  kegg_ids <- grepl("^[a-z]{3,4}[0-9]{5}$", as.character(kk_df$ID))
  if (!isTRUE(cfg$analysis$run_pathview)) {
    writeLines("pathview maps skipped: pathview_download = FALSE.",
               file.path(out_dir, "pathview_skipped.txt"))
  } else if (!requireNamespace("pathview", quietly = TRUE)) {
    writeLines("pathview maps skipped: the pathview package is not installed.",
               file.path(out_dir, "pathview_skipped.txt"))
  } else if (!any(kegg_ids)) {
    writeLines(paste0("pathview maps skipped: this run's pathway identifiers come from ",
                      kegg_source, " and are not KEGG pathway ids."),
               file.path(out_dir, "pathview_skipped.txt"))
  } else {
    # pathview writes into the working directory and resolves data(bods) only
    # when the package is attached, so it is loaded onto the search path.
    if (!"package:pathview" %in% search()) {
      suppressPackageStartupMessages(library("pathview", character.only = TRUE, quietly = TRUE))
    }
    paths <- utils::head(kk_df[order(kk_df$p.adjust), , drop = FALSE],
                         as.integer(cfg$analysis$pathview_max_pathways %||% 12L))
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(out_dir)
    # pathview fetches the KGML and the native map for every pathway, and those
    # requests are as slow as the KEGG REST ones, so the default 60 s timeout is
    # not enough for a multi-pathway loop.
    old_timeout <- getOption("timeout")
    on.exit(options(timeout = old_timeout), add = TRUE)
    options(timeout = 1800)
    gene_data <- rank[!duplicated(names(rank))]
    for (i in seq_len(nrow(paths))) {
      pid <- paths$ID[i]
      err <- tryCatch({
        pathview::pathview(gene.data = gene_data, pathway.id = pid,
                           species = cfg$species$pathview_species,
                           gene.idtype = cfg$species$pathview_id_type %||% "ENTREZ",
                           limit = list(gene = display_lfc_limit(cfg, gene_data), cpd = 1))
        NULL
      }, error = function(e) conditionMessage(e))
      produced <- list.files(pattern = paste0("^", pid, ".*\\.png$"))
      if (length(produced)) rendered <- rendered + 1L
      else errors <- c(errors, paste0(pid, ": ", err %||% "no output file produced"))
    }
  }
  if (length(errors)) writeLines(errors, file.path(out_dir, paste0("pathview_errors_", ctx$contrast_name, ".txt")))

  ctx$kegg <- kk
  ctx$kegg_df <- kk_df
  write_status(file.path(out_dir, "status.json"), "section04_kegg_pathview", "completed", NULL,
               list(significant_pathways = nrow(sig), rendered_pathways = rendered,
                    pathview_enabled = isTRUE(cfg$analysis$run_pathview),
                    kegg_key_type = kegg_id_type, species = cfg$species$kegg_code,
                    kegg_source = kegg_source))
  ctx
}
