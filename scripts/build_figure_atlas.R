#!/usr/bin/env Rscript
# Assemble every figure produced by the pipeline into one browsable PDF atlas
# (pure R implementation).
#
# Usage (from the project root):
#   Rscript scripts/build_figure_atlas.R
#   Rscript scripts/build_figure_atlas.R results/RNAseqSampleData
#   Rscript scripts/build_figure_atlas.R results/RNAseqSampleData out.pdf
#
# Without a result directory the most recent run under results/ is used.
# Figures follow the order of the paper: QC -> differential expression ->
# GO over-representation -> GSEA -> KEGG -> TF regulons -> TF activity ->
# TF network -> TF modules -> gene-set display (GO heat maps, TF dot plots) ->
# protein interactions -> subcellular localization; the original KEGG pathway
# maps form the appendix.

suppressPackageStartupMessages(library(magick))

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(".", winslash = "/")

latest_results_dir <- function(root) {
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[dir.exists(file.path(dirs, "section01_screening"))]
  if (!length(dirs)) stop("No result directory containing section01_screening under ", root)
  dirs[which.max(file.info(file.path(dirs, "section01_screening"))$mtime)]
}

base <- if (length(args) >= 1L) normalizePath(args[1L], winslash = "/", mustWork = TRUE) else
  latest_results_dir(file.path(project_root, "results"))
out <- if (length(args) >= 2L) args[2L] else file.path(base, "RNAseq_results_atlas.pdf")

# (section title, [(caption, relative path), ...]) in paper reading order.
sections <- list(
  list("1. Quality control and sample overview", list(
    c("Sequencing depth (library size)", "00_differential_expression/qc_library_size.png"),
    c("Count distribution", "00_differential_expression/qc_count_distributions.png"),
    c("Principal component analysis (PCA)", "00_differential_expression/qc_pca.png"),
    c("Sample correlation (VST)", "00_differential_expression/qc_sample_correlation.png"),
    c("Sample distance (VST)", "00_differential_expression/qc_sample_distance.png"),
    c("Dispersion estimates", "00_differential_expression/qc_dispersion.png"),
    c("MA plot", "00_differential_expression/qc_ma_plot.png"))),
  list("2. Differential expression (section 1)", list(
    c("Volcano plot (up-regulated triangles point up, down-regulated triangles point down)",
      "section01_screening/S1_volcano_Treatment_vs_Control.png"),
    c("Expression versus effect size", "section01_screening/S1_expression_effect_Treatment_vs_Control.png"),
    c("Differentially expressed genes at each threshold", "section01_screening/S1_threshold_counts.png"),
    c("Heat map of all DEGs (|log2FC| >= 1)", "section01_screening/S1_DEG_global_heatmap_lfc1.png"),
    c("Heat map of the top DEGs (VST, row z-score)", "section01_screening/S1_topDEG_heatmap_rowz.png"))),
  list("3. GO over-representation (section 2)", list(
    c("GO biological-process dot plot", "section02_go_ora/GO_BP_dotplot_Treatment_vs_Control.png"),
    c("GO rich factor by ontology", "section02_go_ora/GO_richfactor_by_ontology_Treatment_vs_Control.png"),
    c("GO functional classes (ten classes, coloured network)",
      "section02_go_ora/GO_BP_functional_classes_Treatment_vs_Control.png"),
    c("GO molecular-function dot plot", "section02_go_ora/GO_MF_dotplot_Treatment_vs_Control.png"),
    c("GO cellular-component dot plot", "section02_go_ora/GO_CC_dotplot_Treatment_vs_Control.png"))),
  list("4. GSEA (section 3)", list(
    c("Composite GSEA panel, GO biological processes",
      "section03_gsea_go/GSEA_GO_BP_composite_Treatment_vs_Control.png"),
    c("Composite GSEA panel, GO molecular functions",
      "section03_gsea_go/GSEA_GO_MF_composite_Treatment_vs_Control.png"),
    c("Composite GSEA panel, GO cellular components",
      "section03_gsea_go/GSEA_GO_CC_composite_Treatment_vs_Control.png"),
    c("GSEA enrichment map, GO biological processes",
      "section03_gsea_go/GSEA_GO_BP_enrichment_map_Treatment_vs_Control.png"),
    c("GO biological-process ridge plot", "section03_gsea_go/GSEA_GO_BP_ridgeplot_Treatment_vs_Control.png"),
    c("MSigDB Hallmark GSEA (NES)", "section03_gsea_go/GSEA_MSigDB_H_barplot_Treatment_vs_Control.png"))),
  list("5. KEGG pathways (section 4)", list(
    c("KEGG enrichment lollipop", "section04_kegg_pathview/KEGG_lollipop_Treatment_vs_Control.png"),
    c("Composite GSEA panel, KEGG pathways",
      "section04_kegg_pathview/KEGG_gsea_panel_Treatment_vs_Control.png"),
    c("UpSet plot of pathway gene overlap", "section04_kegg_pathview/KEGG_upset_Treatment_vs_Control.png"))),
  list("6. Transcription-factor regulons (section 5)", list(
    c("Regulated targets of each TF (dot size = target count, colour = TF log2FC)",
      "section05_tf_regulon/TF_dotplot_Treatment_vs_Control.png"),
    c("TF expression (TFs with maximum group-mean FPKM >= 10)",
      "section05_tf_regulon/TF_expression_dotplot_Treatment_vs_Control.png"),
    c("Regulon enrichment versus regulon activity",
      "section05_tf_regulon/TF_regulon_volcano_Treatment_vs_Control.png"))),
  list("7. Regulatory and pathway activity (section 6)", list(
    c("Differential TF activity volcano plot (ULM)",
      "section06_regulatory_activity/TF_activity_volcano_ulm.png"),
    c("Heat map of differential TF activity (ULM)",
      "section06_regulatory_activity/TF_activity_heatmap_ulm.png"),
    c("Differential pathway activity (GSVA, Hallmark)",
      "section06_regulatory_activity/GSVA_H_diff_barplot.png"))),
  list("8. TF-target regulatory network (section 7)", list(
    c("Force-directed layout (FR)", "section07_tf_network/TF_network_fr_Treatment_vs_Control.png"),
    c("Kamada-Kawai layout", "section07_tf_network/TF_network_kk_Treatment_vs_Control.png"),
    c("Stress layout", "section07_tf_network/TF_network_stress_Treatment_vs_Control.png"),
    c("Circular layout", "section07_tf_network/TF_network_circle_Treatment_vs_Control.png"))),
  list("9. TF module network (section 8)", list(
    c("TF modules obtained by clustering shared differentially expressed targets (size = target count, colour = module, edge = Jaccard overlap)",
      "section08_tf_modules/TF_cluster_network_Treatment_vs_Control.png"))),
  list("10. Gene-set display: GO heat maps (section 9)", list()),
  list("11. Gene-set display: TF-regulon dot plots (section 9)", list()),
  list("12. Protein-protein interaction network (section 10)", list(
    c("STRING interaction network of DEG-encoded proteins (fill = log2FC, size = degree, module functions in the legend)",
      "section10_ppi_network/PPI_network_Treatment_vs_Control.png"),
    c("Top 30 hub proteins by interaction degree (at least 5% per direction)",
      "section10_ppi_network/PPI_hub_proteins_Treatment_vs_Control.png"),
    c("GO functional annotation of each module (two leading terms per module)",
      "section10_ppi_network/PPI_module_functions_Treatment_vs_Control.png"))),
  list("13. Subcellular localization (section 11)", list(
    c("Subcellular localization of DEG-encoded proteins (cell diagram and up/down bar chart)",
      "section11_subcellular_localization/subcellular_localization_Treatment_vs_Control.png")))
)

# Every figure that belongs to a comparison carries the contrast label in its
# name. Read it from the run instead of assuming the sample-data default.
detect_contrast <- function(dir) {
  f <- grep("^S1_volcano_.+\\.png$", list.files(dir), value = TRUE)
  if (!length(f)) return("Treatment_vs_Control")
  sub("^S1_volcano_(.+)\\.png$", "\\1", f[1])
}
contrast <- if (dir.exists(file.path(base, "section01_screening")))
  detect_contrast(file.path(base, "section01_screening")) else "Treatment_vs_Control"
sections <- lapply(sections, function(sec) {
  sec[[2]] <- lapply(sec[[2]], function(item) {
    item[2] <- sub("Treatment_vs_Control", contrast, item[2], fixed = TRUE)
    item
  })
  sec
})

# The global DEG heat map carries the primary |log2FC| cut-off in its name.
glob_hm <- list.files(file.path(base, "section01_screening"),
                      "^S1_DEG_global_heatmap_lfc.*\\.png$", full.names = FALSE)
if (length(glob_hm)) {
  cut <- sub("^S1_DEG_global_heatmap_lfc(.*)\\.png$", "\\1", glob_hm[1])
  sections[[2]][[2]][[5]] <- c(paste0("Heat map of all DEGs (|log2FC| >= ", cut, ")"),
                               file.path("section01_screening", glob_hm[1]))
}

# Section 1 writes one chromosome landscape per effect-size threshold
# (|log2FC| >= 0.5, 0.75, 1.0); insert them after the threshold ladder.
s1_chr <- list()
s1_dir <- file.path(base, "section01_screening")
if (dir.exists(s1_dir)) {
  chr_files <- sort(list.files(s1_dir, "^S1_chromosome_landscape_lfc.*\\.png$", full.names = FALSE))
  s1_chr <- lapply(chr_files, function(f) {
    tag <- sub("^S1_chromosome_landscape_(lfc[0-9.]+)_.*$", "\\1", f)
    tag <- sub("^lfc", "|log2FC| >= ", tag)
    c(paste0("Chromosomal position of DEGs (", tag,
             "; up, above the track; down, below; bar height proportional to |log2FC|)"),
      file.path("section01_screening", f))
  })
}
if (length(s1_chr)) sections[[2]][[2]] <- append(sections[[2]][[2]], s1_chr, after = 3L)

# Section 2 writes one term-gene network per functional class. The class-to-term
# map is in GO_term_gene_network_terms_*.csv, so each caption can name the term
# that was chosen (the one with the most genes in that class).
s2 <- file.path(base, "section02_go_ora")
s2_net <- list()
if (dir.exists(s2)) {
  net_files <- sort(list.files(s2, "^GO_term_gene_network_class[0-9]+_.*\\.png$", full.names = FALSE))
  map_file <- list.files(s2, "^GO_term_gene_network_terms_.*\\.csv$", full.names = TRUE)
  term_map <- if (length(map_file)) tryCatch(utils::read.csv(map_file[1], stringsAsFactors = FALSE),
                                             error = function(e) NULL) else NULL
  s2_net <- lapply(net_files, function(f) {
    cls <- as.integer(sub("^GO_term_gene_network_class([0-9]+)_.*$", "\\1", f))
    term <- if (!is.null(term_map) && any(term_map$class == cls)) {
      term_map$term[term_map$class == cls][1]
    } else paste("class", cls)
    c(paste0("Term-gene network of functional class ", cls, ": ", term,
             " (top 100 genes by |log2FC|)"),
      file.path("section02_go_ora", f))
  })
}
sections[[3]][[2]] <- c(sections[[3]][[2]], s2_net)

# Section 9 draws one figure per gene set: GO terms as row z-scored heat maps,
# TF regulons as expression-deviation dot plots. The atlas separates the two
# families into consecutive sections, and the sets depend on the enrichment
# result, so files are discovered from the output directory.
s9 <- file.path(base, "section09_geneset_heatmaps")
s9_entries <- function(pattern, caption_suffix) {
  if (!dir.exists(s9)) return(list())
  files <- sort(list.files(s9, pattern, full.names = FALSE))
  pretty <- function(stem) gsub("_", " ", stem)
  lapply(files, function(f) {
    stem <- sub("\\.png$", "", f)
    c(paste0(pretty(sub("^(heatmap_rowz_|dotplot_)", "", stem)), " ", caption_suffix),
      file.path("section09_geneset_heatmaps", f))
  })
}
sections[[10]][[2]] <- s9_entries("^heatmap_rowz_GO.*\\.png$", "heat map (row z-score)")
sections[[11]][[2]] <- s9_entries("^dotplot_TF_.*\\.png$", "expression-deviation dot plot")

# Section 6 additionally writes one chromosome map per leading regulon plus a
# TF x chromosome summary; those files depend on which TFs rank highest, so they
# are discovered from the output directory as well.
s5 <- file.path(base, "section05_tf_regulon")
s6_chr <- list()
if (dir.exists(s5)) {
  all_chr <- sort(list.files(s5, "^TF_chromosome_.*\\.png$", full.names = FALSE))
  summ <- all_chr[grepl("^TF_chromosome_summary_", all_chr)]
  maps <- all_chr[!grepl("^TF_chromosome_summary_", all_chr)]
  s6_chr <- c(
    lapply(summ, function(f) c("Chromosomal distribution of targets for the top ten regulons",
                               file.path("section05_tf_regulon", f))),
    lapply(maps, function(f) c(paste0(sub("^TF_chromosome_", "", sub("\\.png$", "", f)),
                                       " target genes on chromosomes"),
                               file.path("section05_tf_regulon", f))))
}
sections[[6]][[2]] <- c(sections[[6]][[2]], s6_chr)

# Each figure is scaled to 2000 px wide (never enlarged), a white band is added
# on top for the section title and caption, and every page is written as a
# one-page PDF before concatenation.
page_width <- 2000L
margin <- 190L
font <- "Noto Sans"
page_files <- character(0)
missing <- character(0)
for (sec in sections) {
  for (item in sec[[2]]) {
    caption <- item[1]
    rel <- item[2]
    path <- file.path(base, rel)
    if (!file.exists(path)) {
      missing <- c(missing, rel)
      next
    }
    img <- magick::image_read(path)
    info <- magick::image_info(img)
    if (info$width > page_width) {
      img <- magick::image_resize(img, paste0(page_width, "x"))
      info <- magick::image_info(img)
    }
    img <- magick::image_extent(
      img, paste0(info$width, "x", as.integer(info$height) + margin),
      gravity = "south", color = "white")
    img <- magick::image_annotate(img, sec[[1]], size = 48, weight = 700,
                                  gravity = "north", location = "+0+16",
                                  font = font, color = "black")
    img <- magick::image_annotate(img, caption, size = 36, weight = 400,
                                  gravity = "north", location = "+0+96",
                                  font = font, color = "grey20")
    page_file <- tempfile(fileext = ".pdf")
    magick::image_write(img, page_file, format = "pdf", density = "200",
                        compression = "jpeg", quality = 92)
    page_files <- c(page_files, page_file)
    rm(img); gc(verbose = FALSE)
  }
}

# The original KEGG pathway maps form the appendix (first six).
kegg <- sort(list.files(file.path(base, "section04_kegg_pathview"),
                        "\\.pathview\\.png$", full.names = TRUE))
for (p in utils::head(kegg, 6L)) {
  img <- magick::image_read(p)
  info <- magick::image_info(img)
  if (info$width > 2200) {
    img <- magick::image_resize(img, "2200x")
    info <- magick::image_info(img)
  }
  img <- magick::image_extent(img, paste0(info$width, "x", as.integer(info$height) + margin),
                              gravity = "south", color = "white")
  img <- magick::image_annotate(img, paste0("Appendix. KEGG pathway map - ", basename(p)),
                                size = 34, weight = 700, gravity = "north",
                                location = "+0+30", font = font, color = "black")
  page_file <- tempfile(fileext = ".pdf")
  magick::image_write(img, page_file, format = "pdf", density = "200",
                      compression = "jpeg", quality = 92)
  page_files <- c(page_files, page_file)
  rm(img); gc(verbose = FALSE)
}

if (!length(page_files)) stop("No images found; run the pipeline before building the atlas.")
# Each page is written as a one-page PDF and the pages are then concatenated,
# which keeps memory flat even for an atlas of 2000-px panels.
n_pages <- length(page_files)
pdftools::pdf_combine(page_files, output = out)
unlink(page_files)

cat("Atlas written:", out, "\n")
cat("Pages:", n_pages, "\n")
if (length(missing)) cat("Missing figures:", paste(missing, collapse = ", "), "\n")
