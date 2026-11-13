#!/usr/bin/env Rscript
# Standalone GO term-target gene network.
#
# Draws a bipartite network for an arbitrary list of GO terms: every term is a
# square node, every annotated gene a circle coloured by its log2 fold change,
# and an edge links a term to each of its target genes. The layout is Kamada-
# Kawai (igraph::layout_with_kk), which keeps the terms apart while pulling the
# genes of a shared programme together.
#
# The script reads the differential-expression table the pipeline already wrote,
# so it does not re-run DESeq2. The gene universe is every gene in that table,
# not only the DEGs: a term network is meant to show the complete target set of
# the chosen functions, with the fold change carrying the evidence.
#
# Usage (from the project root):
#   Rscript scripts/go_term_gene_network.R \
#       --terms "GO:0006281,GO:0000278,DNA repair" \
#       --results results/RNAseqSampleData --species mouse \
#       --out results/RNAseqSampleData/go_term_target_network
#
#   --terms        comma/semicolon separated GO ids or exact term names
#   --term-file    text file with one id or name per line (in addition to --terms)
#   --de           explicit DE table; default <results>/00_differential_expression/
#                  unshrunk DESeq2_*.csv, falling back to the shrinkage table
#   --species      mouse (default) | human | rat
#   --out          output directory (default <results>/go_term_target_network)
#   --max-genes    genes drawn per term, strongest |log2FC| first (default 80; 0 = all)
#   --up-cutoff    display-only symmetric log2FC colour limit (default 2.5)
#   --label-top    gene labels drawn, strongest |log2FC| first (default 60)
#   --layout       kk (default) | fr | stress | circle
#   --width/--height  canvas in inches (default 14 x 10)
#   --title        figure title
#   --seed         layout seed (default 1)

args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("-h", "--help"))) {
  cat("GO term-target gene network (Kamada-Kawai layout)\n\n",
      "Usage:\n",
      "  Rscript scripts/go_term_gene_network.R --terms \"GO:0006281,DNA repair\" \\\n",
      "      --results results/RNAseqSampleData --species mouse \\\n",
      "      --out results/RNAseqSampleData/go_term_target_network\n\n",
      "Options:\n",
      "  --terms       comma/semicolon separated GO ids or exact term names\n",
      "  --term-file   text file with one id or name per line\n",
      "  --de          explicit DE table (default: DESeq2 table in --results)\n",
      "  --results     results directory (default results/RNAseqSampleData)\n",
      "  --species     mouse (default) | human | rat\n",
      "  --out         output directory (default <results>/go_term_target_network)\n",
      "  --max-genes   genes drawn per term (default 80; 0 = all)\n",
      "  --up-cutoff   display-only log2FC colour limit (default 2.5)\n",
      "  --label-top   gene labels drawn (default 60)\n",
      "  --layout      kk (default) | fr | stress | circle\n",
      "  --width       canvas width in inches (default 14)\n",
      "  --height      canvas height in inches (default 10)\n",
      "  --title       figure title\n",
      "  --seed        layout seed (default 1)\n", sep = "")
  quit(status = 0)
}

project_root <- normalizePath(".", winslash = "/")
source(file.path(project_root, "pipeline", "R", "helpers.R"))

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  i <- hit[length(hit)] + 1L
  if (i > length(args)) stop(flag, " needs a value", call. = FALSE)
  args[i]
}

has_flag <- function(flag) flag %in% args

terms_arg <- get_arg("--terms", "")
term_file <- get_arg("--term-file", NULL)
results_dir <- get_arg("--results", "results/RNAseqSampleData")
species <- tolower(get_arg("--species", "mouse"))
out_dir <- get_arg("--out", file.path(results_dir, "go_term_target_network"))
max_genes <- as.integer(get_arg("--max-genes", "80"))
label_top <- as.integer(get_arg("--label-top", "60"))
layout_name <- tolower(get_arg("--layout", "kk"))
width <- as.numeric(get_arg("--width", "14"))
height <- as.numeric(get_arg("--height", "10"))
title_arg <- get_arg("--title", NULL)
seed <- as.integer(get_arg("--seed", "1"))
up_cutoff <- as.numeric(get_arg("--up-cutoff", "2.5"))
if (length(up_cutoff) != 1L || !is.finite(up_cutoff) || up_cutoff <= 0)
  stop("--up-cutoff must be positive and finite", call. = FALSE)

if (!nzchar(terms_arg) && is.null(term_file)) {
  stop("Nothing to draw: pass --terms \"GO:0006281,DNA repair\" or --term-file terms.txt",
       call. = FALSE)
}

requested <- character(0)
if (nzchar(terms_arg)) {
  requested <- c(requested, strsplit(terms_arg, "[,;]")[[1]])
}
if (!is.null(term_file)) {
  if (!file.exists(term_file)) stop("Term file not found: ", term_file, call. = FALSE)
  requested <- c(requested, readLines(term_file, warn = FALSE))
}
requested <- trimws(requested)
requested <- requested[nzchar(requested) & !grepl("^#", requested)]
if (!length(requested)) stop("The term list is empty", call. = FALSE)

species_presets <- list(
  mouse = list(orgdb = "org.Mm.eg.db"),
  human = list(orgdb = "org.Hs.eg.db"),
  rat   = list(orgdb = "org.Rn.eg.db"))
if (!species %in% names(species_presets)) {
  stop("Unknown species '", species, "'; use mouse, human or rat", call. = FALSE)
}
orgdb_name <- species_presets[[species]]$orgdb
for (pkg in c(orgdb_name, "AnnotationDbi")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Install it with BiocManager::install('",
         pkg, "').", call. = FALSE)
  }
}

# --- differential-expression table ------------------------------------------
de_file <- get_arg("--de", NULL)
if (is.null(de_file)) {
  de_dir <- file.path(results_dir, "00_differential_expression")
  shrunk <- list.files(de_dir, "^DESeq2_.*_lfc_shrunk\\.csv$", full.names = TRUE)
  plain <- setdiff(list.files(de_dir, "^DESeq2_.*\\.csv$", full.names = TRUE), shrunk)
  de_file <- c(plain, shrunk)[1]
}
if (is.null(de_file) || is.na(de_file) || !file.exists(de_file)) {
  stop("No differential-expression table found. Pass --de <file> or run the pipeline first.",
       call. = FALSE)
}
de <- utils::read.csv(de_file, stringsAsFactors = FALSE, check.names = FALSE)
sym_col <- intersect(c("GeneSymbol", "SYMBOL", "symbol", "gene_symbol"), names(de))[1]
if (is.na(sym_col)) stop("No gene-symbol column in ", de_file, call. = FALSE)
if (!"log2FoldChange" %in% names(de)) stop("No log2FoldChange column in ", de_file, call. = FALSE)
padj_col <- intersect(c("padj", "p.adjust"), names(de))[1]
de$symbol <- as.character(de[[sym_col]])
de <- de[!is.na(de$symbol) & nzchar(de$symbol), , drop = FALSE]
de$log2FC <- suppressWarnings(as.numeric(de$log2FoldChange))
de$padj <- if (!is.na(padj_col)) suppressWarnings(as.numeric(de[[padj_col]])) else NA_real_
# One row per symbol, keeping the strongest effect, so a gene is never drawn twice.
de <- de[order(-abs(ifelse(is.finite(de$log2FC), de$log2FC, 0))), , drop = FALSE]
de <- de[!duplicated(de$symbol), , drop = FALSE]
if (!nrow(de)) stop("The DE table has no usable gene symbols", call. = FALSE)

# --- GO annotation -----------------------------------------------------------
# OrgDb supplies the gene-to-term mapping; the human-readable term name lives in
# GO.db, so a second lookup resolves the ids to their names.
if (!requireNamespace("GO.db", quietly = TRUE)) {
  stop("Package 'GO.db' is required for GO term names. Install it with ",
       "BiocManager::install('GO.db').", call. = FALSE)
}
message("Reading GO annotations from ", orgdb_name, " for ", nrow(de), " genes ...")
orgdb <- getExportedValue(orgdb_name, orgdb_name)
ann <- tryCatch(
  AnnotationDbi::select(orgdb, keys = unique(de$symbol), keytype = "SYMBOL",
                        columns = c("GO", "ONTOLOGY")),
  error = function(e) {
    stop("GO annotation lookup failed in ", orgdb_name, ": ", conditionMessage(e),
         call. = FALSE)
  })
ann <- ann[!is.na(ann$GO), , drop = FALSE]
if (!nrow(ann)) stop(orgdb_name, " returned no GO annotation for these genes", call. = FALSE)
ann <- unique(ann[, c("SYMBOL", "GO", "ONTOLOGY")])

# Term names come from the full GO.db vocabulary, not only from the terms with a
# direct annotation in this gene list: a user copies names straight out of an
# enrichment table, and those are very often parents whose own direct
# annotations are sparse.
all_terms <- tryCatch(
  AnnotationDbi::select(GO.db::GO.db,
                        keys = AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID"),
                        columns = c("TERM", "ONTOLOGY"), keytype = "GOID"),
  error = function(e) NULL)
if (is.null(all_terms) || !nrow(all_terms)) stop("GO.db returned no term names", call. = FALSE)
all_terms <- all_terms[!is.na(all_terms$TERM) & nzchar(all_terms$TERM), , drop = FALSE]
all_terms <- unique(all_terms[, c("GOID", "TERM", "ONTOLOGY")])

resolve_term <- function(x) {
  if (grepl("^GO:[0-9]{7}$", x, ignore.case = TRUE)) {
    id <- toupper(x)
    if (!id %in% all_terms$GOID) stop("Unknown GO id ", id, call. = FALSE)
    return(id)
  }
  hit <- all_terms$GOID[tolower(all_terms$TERM) == tolower(x)]
  if (length(hit)) return(unique(hit)[1])
  near <- unique(all_terms$TERM[tolower(all_terms$TERM) == tolower(x)])
  if (!length(near)) {
    near <- unique(all_terms$TERM[grepl(tolower(x), tolower(all_terms$TERM), fixed = TRUE)])
  }
  if (length(near) == 1L) return(all_terms$GOID[all_terms$TERM == near][1])
  if (length(near) > 1L) {
    stop("Term name '", x, "' is ambiguous; candidates: ",
         paste(utils::head(near, 8L), collapse = " | "), call. = FALSE)
  }
  stop("No GO term matches '", x, "'. Pass a GO id (GO:0006281) or the exact term name.",
       call. = FALSE)
}

term_ids <- vapply(requested, resolve_term, character(1))
term_ids <- unique(term_ids)
term_meta <- all_terms[match(term_ids, all_terms$GOID), , drop = FALSE]
names(term_meta)[1] <- "GO"
message("Drawing ", nrow(term_meta), " GO term(s): ",
        paste(term_meta$TERM, collapse = "; "))

# A term's targets are the genes annotated to the term itself or to any of its
# descendants, which is the same gene set an enrichment test scores.
children_map <- list(
  BP = AnnotationDbi::as.list(GO.db::GOBPCHILDREN),
  MF = AnnotationDbi::as.list(GO.db::GOMFCHILDREN),
  CC = AnnotationDbi::as.list(GO.db::GOCCCHILDREN))
desc_cache <- new.env(parent = emptyenv())
term_with_descendants <- function(id, ont) {
  key <- paste(ont, id)
  if (!is.null(desc_cache[[key]])) return(desc_cache[[key]])
  out <- id
  frontier <- id
  while (length(frontier)) {
    kids <- unique(unlist(children_map[[ont]][frontier], use.names = FALSE))
    kids <- kids[!is.na(kids) & nzchar(kids)]
    kids <- setdiff(kids, out)
    out <- c(out, kids)
    frontier <- kids
  }
  desc_cache[[key]] <- out
  out
}

# --- term-target edges -------------------------------------------------------
fc_of <- stats::setNames(de$log2FC, de$symbol)
padj_of <- stats::setNames(de$padj, de$symbol)
edge_rows <- list()
summary_rows <- list()
for (i in seq_len(nrow(term_meta))) {
  id <- term_meta$GO[i]
  ont <- term_meta$ONTOLOGY[i]
  ids <- term_with_descendants(id, ont)
  genes <- unique(ann$SYMBOL[ann$GO %in% ids])
  genes <- genes[genes %in% de$symbol]
  if (!length(genes)) next
  ord <- order(-abs(ifelse(is.finite(fc_of[genes]), fc_of[genes], 0)))
  keep <- if (is.finite(max_genes) && max_genes > 0L) {
    utils::head(genes[ord], max_genes)
  } else {
    genes[ord]
  }
  edge_rows[[length(edge_rows) + 1L]] <- data.frame(
    term_id = id, term = term_meta$TERM[i], ontology = term_meta$ONTOLOGY[i],
    gene = keep, log2FC = unname(fc_of[keep]), padj = unname(padj_of[keep]),
    stringsAsFactors = FALSE)
  top <- keep[which.max(abs(ifelse(is.finite(fc_of[keep]), fc_of[keep], 0)))]
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    term_id = id, term = term_meta$TERM[i], ontology = term_meta$ONTOLOGY[i],
    genes_in_term = length(genes), genes_drawn = length(keep),
    top_gene = top, top_abs_log2FC = abs(unname(fc_of[top])),
    stringsAsFactors = FALSE)
}
if (!length(edge_rows)) {
  stop("None of the requested terms has an annotated gene in the DE table", call. = FALSE)
}
edges <- do.call(rbind, edge_rows)
summary_tab <- do.call(rbind, summary_rows)

out_dir <- safe_dir(out_dir)

# --- graph and layout --------------------------------------------------------
term_nodes <- data.frame(
  name = paste0("T|", summary_tab$term_id),
  label = summary_tab$term, ontology = summary_tab$ontology,
  count = summary_tab$genes_drawn, is_term = TRUE, log2FC = NA_real_,
  stringsAsFactors = FALSE)
gene_tab <- unique(edges[, c("gene", "log2FC")])
gene_tab <- gene_tab[order(-abs(ifelse(is.finite(gene_tab$log2FC), gene_tab$log2FC, 0))), , drop = FALSE]
gene_nodes <- data.frame(
  name = paste0("G|", gene_tab$gene),
  label = gene_tab$gene, ontology = NA_character_,
  count = NA_integer_, is_term = FALSE, log2FC = gene_tab$log2FC,
  stringsAsFactors = FALSE)
nodes <- rbind(term_nodes, gene_nodes)

el <- data.frame(from = paste0("T|", edges$term_id), to = paste0("G|", edges$gene))
g <- igraph::graph_from_data_frame(el, directed = FALSE, vertices = nodes)

set.seed(seed)
layout_fun <- switch(layout_name,
                     kk = igraph::layout_with_kk,
                     fr = igraph::layout_with_fr,
                     stress = function(g) {
                       if (!requireNamespace("graphlayouts", quietly = TRUE))
                         stop("The stress layout requires graphlayouts.", call. = FALSE)
                       graphlayouts::layout_with_stress(g)
                     },
                     circle = igraph::layout_in_circle,
                     stop("Unknown --layout '", layout_name, "'; use kk, fr, stress or circle",
                          call. = FALSE))
lay <- layout_fun(g)
lay <- apply(lay, 2, function(v) (v - mean(v)) / max(stats::sd(v), 1e-9))
rownames(lay) <- igraph::V(g)$name

# --- figure ------------------------------------------------------------------
node_pos <- data.frame(name = rownames(lay), x = lay[, 1], y = lay[, 2],
                       stringsAsFactors = FALSE)
nodes <- merge(nodes, node_pos, by = "name", sort = FALSE)
ed <- data.frame(from = el$from, to = el$to, stringsAsFactors = FALSE)
ed$x <- nodes$x[match(ed$from, nodes$name)]
ed$y <- nodes$y[match(ed$from, nodes$name)]
ed$xend <- nodes$x[match(ed$to, nodes$name)]
ed$yend <- nodes$y[match(ed$to, nodes$name)]

term_df <- nodes[nodes$is_term, , drop = FALSE]
gene_df <- nodes[!nodes$is_term, , drop = FALSE]
# Keep term colours separate from the green-purple gene fold-change ramp.
term_colours <- c("#FF69B4", "#E69F00", "#00FF7F", "#00F5FF",
                  "#FFE1FF", "#CC3311", "#666666", "#222222")
if (nrow(term_df) > length(term_colours)) {
  n_extra <- nrow(term_df) - length(term_colours)
  extra_hues <- rep(c(35, 245), length.out = n_extra)
  extra_lightness <- rep(seq(40, 70, length.out = ceiling(n_extra / 2)), each = 2)[seq_len(n_extra)]
  term_colours <- c(term_colours, grDevices::hcl(h = extra_hues, c = 55, l = extra_lightness))
}
term_pal <- stats::setNames(term_colours[seq_len(nrow(term_df))], term_df$name)
term_df$label_wrapped <- wrap_label(term_df$label, width = 26L)
# Leave space above the largest square, with the label anchored at its bottom.
term_label_offset <- 0.035 * max(diff(range(nodes$y)), 1)

# Gene labels: the strongest effects only, so a dense term does not turn into a
# wall of text. Every term is always labelled.
label_df <- gene_df[order(-abs(ifelse(is.finite(gene_df$log2FC), gene_df$log2FC, 0))), , drop = FALSE]
if (is.finite(label_top) && label_top > 0L) label_df <- utils::head(label_df, label_top)

fc_cap <- up_cutoff
# A system font discovered by publication_font() may be unavailable to a base
# graphics device. The generic sans family works across R graphics backends.
font_family <- "sans"

p <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = ed,
                        ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = from),
                        linewidth = 0.28, alpha = 0.65) +
  ggplot2::geom_point(data = gene_df,
                      ggplot2::aes(x = x, y = y, fill = log2FC),
                      shape = 21, size = 2.9, colour = "grey25", stroke = 0.35) +
  ggplot2::geom_point(data = term_df,
                      ggplot2::aes(x = x, y = y, size = count, colour = name),
                      shape = 15) +
  ggrepel::geom_text_repel(data = label_df, ggplot2::aes(x = x, y = y, label = label),
                           size = 3.45, family = font_family, fontface = "bold",
                           colour = "grey15", max.overlaps = Inf,
                           min.segment.length = 0, segment.size = 0.25,
                           segment.colour = "grey60", box.padding = 0.18) +
  # Place term names directly above their squares with a transparent background.
  ggplot2::geom_text(data = term_df,
                     ggplot2::aes(x = x, y = y + term_label_offset, label = label_wrapped),
                     size = 4.55, family = font_family, fontface = "bold", colour = "black",
                     hjust = 0.5, vjust = 0) +
  # Term squares and their links share a colour; gene fill encodes fold change.
  ggplot2::scale_fill_gradient2(low = "green4", mid = "grey90", high = "purple3",
                                midpoint = 0, limits = c(-fc_cap, fc_cap), oob = scales::squish,
                                name = expression(log[2] * "(fold change)")) +
  ggplot2::scale_colour_manual(values = term_pal, name = "GO term", drop = FALSE) +
  ggplot2::scale_size_continuous(range = c(1.5, 9),
                                 breaks = size_breaks(term_df$count),
                                 name = "Genes in term") +
  ggplot2::guides(colour = "none", size = "none",
                  fill = ggplot2::guide_colourbar(
                    direction = "vertical", title.position = "left",
                    barwidth = ggplot2::unit(0.34, "cm"),
                    barheight = ggplot2::unit(3.0, "cm"),
                    frame.colour = "grey45", frame.linewidth = 0.3,
                    ticks.colour = "grey35")) +
  ggplot2::theme_void(base_family = font_family, base_size = 13) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 17, hjust = 0),
                 plot.subtitle = ggplot2::element_text(size = 13.5, colour = "grey25", hjust = 0),
                 legend.text = ggplot2::element_text(size = 11.5),
                 legend.title = ggplot2::element_text(size = 12.5, face = "bold", angle = 90,
                                                       hjust = 0.5, vjust = 0.5),
                 # Placed inside the canvas at the lower right, where the layout
                 # leaves free space, so the legend does not take width from the
                 # network.
                 legend.position = "inside",
                 legend.position.inside = c(0.985, 0.04),
                 legend.justification = c(1, 0),
                 legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                 legend.box.background = ggplot2::element_blank(),
                 legend.key = ggplot2::element_rect(fill = "transparent", colour = NA),
                 legend.margin = ggplot2::margin(5, 6, 5, 6),
                 plot.margin = ggplot2::margin(8, 10, 8, 8)) +
  ggplot2::labs(
    title = title_arg %||% "GO term-target gene network",
    subtitle = wrap_label(paste0("Squares mark the ", nrow(term_df),
                      " GO terms, sized by gene count; each term and its links share a colour. Circles mark the ",
                      nrow(gene_df), " target genes (term and descendant annotations). Layout = ",
                      toupper(layout_name), "."), width = 110L),
    x = NULL, y = NULL)

png_path <- file.path(out_dir, "GO_term_target_network.png")
pdf_path <- file.path(out_dir, "GO_term_target_network.pdf")
# Select font-aware devices explicitly instead of relying on local defaults.
if (requireNamespace("ragg", quietly = TRUE)) {
  ggplot2::ggsave(png_path, p, width = width, height = height, units = "in",
                  device = ragg::agg_png, dpi = 320, bg = "white")
} else if (isTRUE(capabilities("cairo"))) {
  ggplot2::ggsave(png_path, p, width = width, height = height, units = "in",
                  device = grDevices::png, type = "cairo", dpi = 320, bg = "white")
} else {
  ggplot2::ggsave(png_path, p, width = width, height = height, units = "in",
                  device = grDevices::png, dpi = 320, bg = "white")
}
# Base PDF supports the generic sans family when Cairo is unavailable.
pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
ggplot2::ggsave(pdf_path, p, width = width, height = height, units = "in",
                device = pdf_device, family = font_family, bg = "white")

# --- tables ------------------------------------------------------------------
node_out <- nodes[, c("name", "label", "is_term", "ontology", "count", "log2FC", "x", "y")]
names(node_out) <- c("node", "label", "is_term", "ontology", "genes_in_term", "log2FC", "x", "y")
safe_write_csv(node_out, file.path(out_dir, "GO_term_target_nodes.csv"))
safe_write_csv(edges, file.path(out_dir, "GO_term_target_edges.csv"))
safe_write_csv(summary_tab, file.path(out_dir, "GO_term_target_summary.csv"))

message("Wrote ", png_path)
message("      ", pdf_path)
message("      GO_term_target_{nodes,edges,summary}.csv")
