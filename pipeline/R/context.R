# Shared context for the integrated analysis sections.
#
# Every section in pipeline/sections/ receives (cfg, ctx, out_dir). The context
# carries the objects that the original scripts either produced in an earlier
# section (mylist, geneList, gene_df) or asked the user to pick interactively
# (the FPKM matrix, the TF database).

# --- gene lengths ------------------------------------------------------------

# When the count workbook carries no usable gene lengths, lengths come from
# species annotation: the union of a gene's exons is summed, which is
# the "exon model" length that FPKM is defined on. Reading the exons of a whole
# transcriptome takes about a minute, so the table is cached per annotation
# package under cache/gene_lengths/.
gene_length_cache_dir <- function() {
  dir <- Sys.getenv("CLAWSGO_GENE_LENGTH_CACHE", "cache/gene_lengths")
  dir <- resolve_project_path(dir, getOption("transatlas.project_root", getwd()))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir
}

# EnsDb for mouse and human: it is keyed by Ensembl gene identifiers, which is
# what the count matrix uses. No Ensembl database is distributed for the rat, so
# the UCSC TxDb is used and the identifiers are translated through org.Rn.eg.db.
gene_length_annotation <- function(species) {
  switch(tolower(species$common_name %||% "mouse"),
         mouse = list(package = "EnsDb.Mmusculus.v79", keytype = "ENSEMBL",
                      label = "EnsDb.Mmusculus.v79 (Ensembl 79, GRCm38)"),
         human = list(package = "EnsDb.Hsapiens.v86", keytype = "ENSEMBL",
                      label = "EnsDb.Hsapiens.v86 (Ensembl 86, GRCh38)"),
         rat = list(package = "TxDb.Rnorvegicus.UCSC.rn6.refGene", keytype = "ENTREZID",
                    label = "TxDb.Rnorvegicus.UCSC.rn6.refGene (UCSC rn6)"),
         NULL)
}

# A GENCODE or Ensembl GTF takes precedence when one is available, because it
# carries the version of the annotation that produced the counts. The file is
# looked for through the TRANSATLAS_GENE_ANNOTATION environment variable and
# then next to the count workbook; without one the annotation package is used.
gene_annotation_file <- function(cfg) {
  given <- if (!is.null(cfg)) cfg$input$gene_annotation else NULL
  if (is.character(given) && length(given) && nzchar(given) && file.exists(given)) {
    return(normalizePath(given, winslash = "/", mustWork = FALSE))
  }
  opt <- getOption("transatlas.gene_annotation", "")
  if (is.character(opt) && length(opt) && nzchar(opt) && file.exists(opt)) {
    return(normalizePath(opt, winslash = "/", mustWork = FALSE))
  }
  env <- Sys.getenv("TRANSATLAS_GENE_ANNOTATION", "")
  if (nzchar(env) && file.exists(env)) return(normalizePath(env, winslash = "/", mustWork = FALSE))
  counts_path <- cfg$input$counts_file
  dir <- if (!is.null(counts_path) && nzchar(counts_path %||% "")) dirname(counts_path) else "input"
  hits <- list.files(dir, pattern = "\\.(gtf|gff3?)(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  if (length(hits)) hits[[which.max(file.info(hits)$mtime)]] else NULL
}

gtf_exon_lengths <- function(path) {
  if (!requireNamespace("rtracklayer", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) return(NULL)
  md5 <- tryCatch(unname(tools::md5sum(path)), error = function(e) NA_character_)
  key <- paste0("gtf_", sub("\\.(gtf|gff3?)(\\.gz)?$", "", basename(path), ignore.case = TRUE),
                "_", substr(ifelse(is.na(md5), "unknown", md5), 1, 12))
  cache <- file.path(gene_length_cache_dir(), paste0(key, ".csv"))
  if (file.exists(cache)) {
    tab <- utils::read.csv(cache, stringsAsFactors = FALSE)
    if (nrow(tab) && all(c("gene_id", "length") %in% names(tab))) {
      attr(tab, "label") <- basename(path)
      attr(tab, "keytype") <- "ENSEMBL"
      return(tab)
    }
  }
  tab <- tryCatch({
    format <- if (grepl("\\.gff3?(\\.gz)?$", path, ignore.case = TRUE)) "gff3" else "gtf"
    features <- rtracklayer::import(path, format = format)
    exons <- features[as.character(GenomicRanges::mcols(features)[["type"]]) == "exon"]
    ids <- as.character(GenomicRanges::mcols(exons)[["gene_id"]])
    if (length(ids) != length(exons) || all(is.na(ids))) {
      # GFF3 relates exons to transcripts, and transcripts to genes via Parent.
      attrs <- GenomicRanges::mcols(features)
      feature_ids <- as.character(attrs[["ID"]])
      parents <- attrs[["Parent"]]
      parent_map <- stats::setNames(as.list(parents), feature_ids)
      gene_ids <- feature_ids[as.character(attrs[["type"]]) %in% c("gene", "pseudogene")]
      ancestors <- function(id, seen = character()) {
        if (is.na(id) || !nzchar(id) || id %in% seen) return(character())
        if (id %in% gene_ids) return(id)
        next_ids <- parent_map[[id]]
        unique(unlist(lapply(as.character(next_ids), ancestors, seen = c(seen, id))))
      }
      exon_parents <- as.list(GenomicRanges::mcols(exons)[["Parent"]])
      gene_for_exon <- lapply(exon_parents, function(p) unique(unlist(lapply(as.character(p), ancestors))))
      exons <- exons[rep(seq_along(exons), lengths(gene_for_exon))]
      ids <- as.character(unlist(gene_for_exon))
    }
    ids <- sub("\\.[0-9]+$", "", ids)
    keep <- !is.na(ids) & nzchar(ids)
    exons <- exons[keep]; ids <- ids[keep]
    if (!length(exons)) stop("No exon-to-gene assignments found in the annotation.")
    # S4Vectors::split dispatches to the GRanges method and returns a
    # GRangesList; base::split would return a plain list of GRanges, on which
    # reduce() has no method.
    by_gene <- S4Vectors::split(exons, ids)
    widths <- GenomicRanges::width(GenomicRanges::reduce(by_gene))
    widths <- if (methods::is(widths, "List")) vapply(as.list(widths), sum, numeric(1)) else
      as.numeric(widths)
    data.frame(gene_id = names(widths), length = as.numeric(widths), stringsAsFactors = FALSE)
  }, error = function(e) {
    message("    gene lengths unavailable from ", basename(path), ": ", conditionMessage(e))
    NULL
  })
  if (is.null(tab) || !nrow(tab)) return(NULL)
  tab <- tab[is.finite(tab$length) & tab$length > 0, , drop = FALSE]
  utils::write.csv(tab, cache, row.names = FALSE)
  attr(tab, "label") <- basename(path)
  attr(tab, "keytype") <- "ENSEMBL"
  tab
}

package_exon_lengths <- function(ann) {
  cache <- file.path(gene_length_cache_dir(), paste0(ann$package, ".csv"))
  if (file.exists(cache)) {
    tab <- utils::read.csv(cache, stringsAsFactors = FALSE)
    if (nrow(tab) && all(c("gene_id", "length") %in% names(tab))) {
      attr(tab, "label") <- ann$label
      attr(tab, "keytype") <- ann$keytype
      return(tab)
    }
  }
  if (!requireNamespace("GenomicFeatures", quietly = TRUE) ||
      !requireNamespace("GenomicRanges", quietly = TRUE)) return(NULL)
  db <- tryCatch(getExportedValue(ann$package, ann$package), error = record_optional_error)
  if (is.null(db)) return(NULL)
  tab <- tryCatch({
    exons <- GenomicFeatures::exonsBy(db, by = "gene")
    widths <- GenomicRanges::width(GenomicRanges::reduce(exons))
    widths <- if (methods::is(widths, "List")) vapply(as.list(widths), sum, numeric(1)) else
      as.numeric(widths)
    data.frame(gene_id = names(widths), length = as.numeric(widths),
               stringsAsFactors = FALSE)
  }, error = function(e) {
    message("    gene lengths unavailable from ", ann$package, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(tab) || !nrow(tab)) return(NULL)
  tab <- tab[is.finite(tab$length) & tab$length > 0, , drop = FALSE]
  utils::write.csv(tab, cache, row.names = FALSE)
  attr(tab, "label") <- ann$label
  attr(tab, "keytype") <- ann$keytype
  tab
}

# The length column shipped with the count workbook. It is the first choice for
# the FPKM display matrix because it documents, in the input file itself, the
# annotation version the counts were quantified against. The column is read from
# the unmodified workbook (not from the annotations table, which carries only the
# identifier and symbol columns).
workbook_gene_lengths <- function(cfg, counts_obj) {
  col <- cfg$input$gene_length_column %||% ""
  if (!is.character(col) || length(col) != 1L || !nzchar(col) ||
      tolower(col) %in% c("none", "null", "na", "false")) return(NULL)
  raw <- counts_obj$raw
  if (is.null(raw) || !col %in% names(raw)) return(NULL)
  id_col <- cfg$input$gene_id_column
  ids <- as.character(raw[[id_col]])
  len <- suppressWarnings(as.numeric(raw[[col]]))
  keep <- !is.na(len) & is.finite(len) & len > 0 & !is.na(ids) & nzchar(ids)
  if (!any(keep)) return(NULL)
  tab <- data.frame(gene_id = ids[keep], length = len[keep], stringsAsFactors = FALSE)
  tab <- tab[!duplicated(tab$gene_id), , drop = FALSE]
  attr(tab, "label") <- paste0("workbook column '", col, "'")
  attr(tab, "keytype") <- "ENSEMBL"
  tab
}

annotation_exon_lengths <- function(species, cfg = NULL) {
  gtf <- if (!is.null(cfg)) gene_annotation_file(cfg) else NULL
  if (!is.null(gtf)) {
    tab <- gtf_exon_lengths(gtf)
    if (!is.null(tab)) return(tab)
  }
  ann <- gene_length_annotation(species)
  if (is.null(ann)) return(NULL)
  tab <- package_exon_lengths(ann)
  if (is.null(tab)) message("    No exon-union lengths available from ", ann$package,
                            "; install this annotation package or supply gene_annotation.")
  tab
}

# Lengths are stored under the annotation's key type; when the count matrix uses
# a different identifier the two are joined through the species OrgDb rather
# than by row position.
map_annotation_lengths <- function(lengths, count_ids, species) {
  count_ids <- as.character(count_ids)
  if (!length(count_ids)) return(stats::setNames(numeric(), character()))
  if (is.null(lengths) || !nrow(lengths)) return(stats::setNames(rep(NA_real_, length(count_ids)), count_ids))
  ann_key <- attr(lengths, "keytype") %||% "ENSEMBL"
  direct <- lengths$length[match(count_ids, lengths$gene_id)]
  if (mean(!is.na(direct)) >= 0.5) return(stats::setNames(direct, count_ids))
  # A versioned identifier (ENSMUSG00000000001.4) is matched on its stem.
  bare <- sub("\\..*$", "", count_ids)
  if (any(bare != count_ids)) {
    from_stem <- lengths$length[match(bare, lengths$gene_id)]
    direct[is.na(direct)] <- from_stem[is.na(direct)]
    if (mean(!is.na(direct)) >= 0.5) return(stats::setNames(direct, count_ids))
  }
  from <- species$gene_id_type %||% "ENSEMBL"
  if (identical(from, ann_key) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace(species$orgdb, quietly = TRUE)) {
    return(stats::setNames(direct, count_ids))
  }
  map <- tryCatch(AnnotationDbi::select(getExportedValue(species$orgdb, species$orgdb),
                                        keys = unique(count_ids), keytype = from,
                                        columns = ann_key),
                  error = record_optional_error)
  if (is.null(map) || !nrow(map)) return(stats::setNames(direct, count_ids))
  translated <- map[[ann_key]][match(count_ids, map[[from]])]
  mapped <- lengths$length[match(translated, lengths$gene_id)]
  direct[is.na(direct)] <- mapped[is.na(direct)]
  stats::setNames(direct, count_ids)
}

# FPKM as defined by Mortazavi et al. (2008): fragments per kilobase of exon
# model per million mapped fragments. Genes without an annotated length are
# dropped, because their FPKM would be undefined; the library size is still the
# sum over every gene of the matrix.
compute_fpkm <- function(counts, lengths) {
  counts <- as.matrix(counts)
  if (!nrow(counts) || !ncol(counts)) return(NULL)
  library_size <- colSums(counts)
  if (any(!is.finite(library_size) | library_size <= 0)) stop("FPKM requires positive library sizes in every sample.", call. = FALSE)
  len <- lengths[rownames(counts)]
  keep <- is.finite(len) & len > 0
  if (!any(keep)) return(NULL)
  counts <- counts[keep, , drop = FALSE]
  len <- as.numeric(len[keep])
  fpkm <- sweep(counts, 2, library_size / 1e6, "/")
  fpkm / (len / 1e3)
}

# --- expression matrix -------------------------------------------------------

# Return an expression matrix keyed by gene identifier, ready for the heatmap
# and dot-plot sections. The order of precedence is: an explicitly supplied
# matrix, a matrix named after the count workbook, FPKM computed from the raw
# counts and the exon-union lengths (the workbook's length column first, the
# annotation only when the column is absent), and finally the DESeq2 normalised
# counts (the fallback used by section 9 when no length source is available).
build_expression_matrix <- function(cfg, counts_obj, output_dir = NULL) {
  supplied <- cfg$input$expression_matrix_file
  if (!is.null(supplied) && !file.exists(supplied)) stop("Expression matrix does not exist: ", supplied, call. = FALSE)
  # A matrix named after the count workbook takes precedence over computing
  # one: 3T3.xlsx is paired with 3T3_fpkm.csv, so a plain
  # run_all(counts_file = "input/3T3.xlsx") still displays that matrix.
  if (is.null(supplied) || !nzchar(supplied %||% "")) {
    counts_path <- cfg$input$counts_file
    if (!is.null(counts_path) && nzchar(counts_path %||% "")) {
      guess <- file.path(dirname(counts_path),
                         paste0(tools::file_path_sans_ext(basename(counts_path)), "_fpkm.csv"))
      if (file.exists(guess)) supplied <- guess
    }
  }
  if (!is.null(supplied) && nzchar(supplied %||% "") && file.exists(supplied)) {
    mat <- if (tolower(tools::file_ext(supplied)) %in% c("xlsx", "xls")) {
      as.data.frame(readxl::read_excel(supplied, .name_repair = "unique"))
    } else {
      utils::read.csv(supplied, check.names = FALSE, stringsAsFactors = FALSE)
    }
    gene_id_col <- cfg$input$gene_id_column
    if (!gene_id_col %in% names(mat)) stop("Expression matrix lacks gene ID column: ", gene_id_col, call. = FALSE)
    ids <- as.character(mat[[gene_id_col]])
    if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids))
      stop("Expression matrix gene IDs must be unique and non-missing.", call. = FALSE)
    rownames(mat) <- ids
    mat[[gene_id_col]] <- NULL
    symbol_col <- cfg$input$gene_symbol_column
    if (!is.null(symbol_col) && symbol_col %in% names(mat)) mat[[symbol_col]] <- NULL
    numeric_cols <- names(mat)[vapply(mat, is.numeric, logical(1))]
    mat <- as.matrix(mat[, numeric_cols, drop = FALSE])
    samples <- colnames(counts_obj$counts)
    if (!all(samples %in% colnames(mat))) stop("Expression matrix is missing count-matrix samples.", call. = FALSE)
    mat <- mat[, samples, drop = FALSE]
    if (anyDuplicated(rownames(mat)) || anyNA(rownames(mat))) stop("Expression matrix gene IDs must be unique and non-missing.", call. = FALSE)
    if (any(!is.finite(mat)) || any(mat < 0)) stop("Expression matrix values must be finite and non-negative.", call. = FALSE)
    unit <- cfg$input$expression_matrix_unit %||% "FPKM"
    if (!is.null(output_dir)) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      write_status(file.path(output_dir, "status_fpkm.json"), "expression_matrix",
                   if (identical(unit, "FPKM")) "completed" else "skipped",
                   if (!identical(unit, "FPKM")) paste0("Supplied expression unit: ", unit) else NULL,
                   list(source = basename(supplied), unit = unit, genes_total = nrow(counts_obj$counts),
                        genes_matched = length(intersect(rownames(mat), rownames(counts_obj$counts)))))
    }
    return(list(matrix = mat, source = basename(supplied), unit = unit))
  }

  lengths <- workbook_gene_lengths(cfg, counts_obj)
  if (is.null(lengths)) lengths <- annotation_exon_lengths(cfg$species, cfg)
  if (!is.null(lengths)) {
    len <- map_annotation_lengths(lengths, rownames(counts_obj$counts), cfg$species)
    fpkm <- compute_fpkm(counts_obj$counts, len)
    if (!is.null(fpkm)) {
      source <- paste0("computed FPKM (", attr(lengths, "label"), ", exon-union lengths)")
      if (!is.null(output_dir)) {
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        utils::write.csv(data.frame(gene_id = rownames(fpkm), fpkm, check.names = FALSE),
                         file.path(output_dir, "FPKM_computed.csv"), row.names = FALSE)
        utils::write.csv(data.frame(gene_id = names(len), length = as.numeric(len)),
                         file.path(output_dir, "gene_lengths_used.csv"), row.names = FALSE)
      }
      matched <- sum(is.finite(len) & len > 0)
      message("    FPKM: ", source, "; lengths matched for ", matched, "/", length(len), " genes.")
      if (!is.null(output_dir)) write_status(file.path(output_dir, "status_fpkm.json"), "fpkm", "completed",
        details = list(source = source, genes_total = length(len), genes_matched = matched,
                       genes_without_length = length(len) - matched))
      return(list(matrix = fpkm, source = source, unit = "FPKM"))
    }
  }
  reason <- "No usable exon-union lengths: workbook, GTF/GFF and species annotation routes did not yield matching gene IDs."
  if (!is.null(output_dir)) write_status(file.path(output_dir, "status_fpkm.json"), "fpkm", "unavailable", reason)
  warning(reason, call. = FALSE)
  list(matrix = NULL, source = "unavailable", unit = NA_character_)
}

# Aggregate an expression matrix to gene symbols.
#
# The matrix and the annotation table are joined on the gene identifier, never
# by row position. A supplied expression file, or the DESeq2 normalised counts
# used when no expression matrix is supplied, can carry a different number of
# rows from the annotation table; rowsum() then stops with "incorrect length for
# 'group'". Rows whose identifier is unknown, or whose symbol is empty, are
# dropped, and repeated symbols are averaged.
symbol_rowsum <- function(expr, annotations, id_col, symbol_col, upper = FALSE) {
  expr <- as.matrix(expr)
  ids <- rownames(expr)
  if (is.null(ids) || !nrow(expr) || is.null(annotations) ||
      is.null(id_col) || !id_col %in% names(annotations) ||
      is.null(symbol_col) || !symbol_col %in% names(annotations)) {
    return(NULL)
  }
  ann_ids <- as.character(annotations[[id_col]])
  syms <- as.character(annotations[[symbol_col]])
  sym <- syms[match(ids, ann_ids)]
  if (all(is.na(sym)) && any(toupper(ids) %in% toupper(syms))) {
    # A matrix keyed by symbol rather than by gene identifier.
    sym <- ids
  }
  keep <- !is.na(sym) & nzchar(sym)
  if (!any(keep)) return(NULL)
  out <- rowsum(expr[keep, , drop = FALSE], group = if (isTRUE(upper)) toupper(sym[keep]) else sym[keep],
                reorder = TRUE)
  n_per <- table(if (isTRUE(upper)) toupper(sym[keep]) else sym[keep])
  out / as.numeric(n_per[rownames(out)])
}

# --- identifiers -------------------------------------------------------------

# Symbol -> key map and the reverse annotation frame used by every section.
build_identifier_context <- function(cfg, ctx_annotations) {
  orgdb_name <- cfg$species$orgdb
  key_type <- cfg$species$entrez_type %||% "ENTREZID"
  symbol_type <- cfg$species$symbol_type %||% "SYMBOL"
  symbol_col <- cfg$input$gene_symbol_column
  symbols <- as.character(ctx_annotations[[symbol_col]])
  gene_map <- NULL
  if (requireNamespace("clusterProfiler", quietly = TRUE) &&
      requireNamespace(orgdb_name, quietly = TRUE) && !is.null(symbol_col)) {
    gene_map <- symbol_to_key_map(symbols, getExportedValue(orgdb_name, orgdb_name),
                                  symbol_type, key_type)
  }
  list(gene_map = gene_map, key_type = key_type, symbol_type = symbol_type,
       symbol_col = symbol_col)
}

# --- transcription factors ---------------------------------------------------

# Load the TF-target resource, standardise its column names, and keep the
# configured confidence tiers. Returns NULL (never fabricated edges) when the
# resource is unavailable. By default the DoRothEA regulons are taken from the
# `dorothea` package data (human = dorothea_hs, mouse/rat = dorothea_mm);
# `tf_file` in run_all() is the override for a custom table.
load_tf_regulons <- function(cfg) {
  source <- tolower(cfg$input$tf_source %||% "dorothea")
  db <- tryCatch({
    if (identical(source, "file")) {
      path <- cfg$input$tf_database_file
      if (is.null(path) || !file.exists(path)) return(NULL)
      read_tf_table(path)
    } else {
      read_dorothea_table(cfg)
    }
  }, error = function(e) {
    message("    TF resource unavailable: ", conditionMessage(e))
    NULL
  })
  if (is.null(db) || nrow(db) == 0L) return(NULL)
  names(db) <- tolower(names(db))
  pick <- function(candidates) {
    hit <- intersect(candidates, names(db))
    if (length(hit)) hit[[1L]] else NA_character_
  }
  tf_col <- pick(c("tf", "source", "regulator", "transcription_factor"))
  target_col <- pick(c("target", "target_gene", "gene", "targets"))
  if (is.na(tf_col) || is.na(target_col)) return(NULL)
  if (!identical(tf_col, "tf")) names(db)[names(db) == tf_col] <- "tf"
  if (!identical(target_col, "target")) names(db)[names(db) == target_col] <- "target"
  conf_col <- pick(c("confidence", "confidence_level"))
  mor_col <- pick(c("mor", "mode_of_regulation", "direction"))
  out <- data.frame(tf = as.character(db$tf), target = as.character(db$target),
                    stringsAsFactors = FALSE)
  out$confidence <- if (!is.na(conf_col)) as.character(db[[conf_col]]) else NA_character_
  out$mor <- if (!is.na(mor_col)) suppressWarnings(as.numeric(db[[mor_col]])) else NA_real_
  out <- out[!is.na(out$tf) & !is.na(out$target) & nzchar(out$tf) & nzchar(out$target), , drop = FALSE]
  keep <- cfg$regulatory$confidence_keep %||% c("A", "B")
  if (!is.na(conf_col) && any(!is.na(out$confidence))) {
    out <- out[toupper(out$confidence) %in% toupper(keep), , drop = FALSE]
  }
  out <- out[!duplicated(out[, c("tf", "target")]), , drop = FALSE]
  attr(out, "source_file") <- if (identical(source, "file")) {
    basename(cfg$input$tf_database_file)
  } else {
    paste0(attr(db, "dataset") %||% "dorothea", " (DoRothEA package data)")
  }
  attr(out, "tiers_kept") <- keep
  out
}

# --- context -----------------------------------------------------------------

section_dir <- function(output_dir, n, slug) {
  safe_dir(file.path(output_dir, sprintf("section%02d_%s", as.integer(n), slug)))
}

# One expression floor for the whole workflow: a gene enters the analysis when
# the higher of its two group-mean FPKM values reaches analysis$min_fpkm
# (default 5). Applying it here, on the contrast table, means the screening
# ladder, the ranked lists of the enrichment sections and the regulatory,
# interaction and localisation sections all work on the same gene universe.
# The filter is skipped, with the reason recorded, when the expression matrix
# is not in FPKM (the normalised-count fallback) or when it shares too few
# identifiers with the DE table.
apply_expression_floor <- function(cfg, de, expr_matrix, unit, meta, group_col) {
  min_fpkm <- if (is.null(cfg$analysis$min_fpkm)) NA_real_ else
    suppressWarnings(as.numeric(cfg$analysis$min_fpkm))
  info <- list(min_fpkm = min_fpkm, unit = unit, applied = FALSE,
               tested = nrow(de), kept = nrow(de), dropped = 0L)
  if (!is.finite(min_fpkm) || min_fpkm <= 0) {
    info$reason <- "disabled (min_fpkm is NULL or not positive)"
    return(list(de = de, info = info, max_mean = NULL))
  }
  if (is.null(expr_matrix) || !nrow(expr_matrix)) {
    info$reason <- "no expression matrix was built"
    return(list(de = de, info = info, max_mean = NULL))
  }
  if (!identical(as.character(unit), "FPKM")) {
    info$reason <- paste0("expression matrix is ", unit %||% "unknown", ", not FPKM")
    return(list(de = de, info = info, max_mean = NULL))
  }
  sid <- as.character(meta$sample_id)
  grp <- as.character(meta[[group_col]])
  names(grp) <- sid
  grp <- grp[colnames(expr_matrix)]
  present <- unique(grp[!is.na(grp)])
  names(present) <- present
  if (length(present) < 2L) {
    info$reason <- "fewer than two groups in the expression matrix"
    return(list(de = de, info = info, max_mean = NULL))
  }
  means <- vapply(present, function(g) {
    cols <- which(grp == g)
    if (!length(cols)) return(rep(NA_real_, nrow(expr_matrix)))
    rowMeans(expr_matrix[, cols, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(expr_matrix)))
  means <- matrix(means, ncol = length(present),
                  dimnames = list(rownames(expr_matrix), present))
  max_mean <- means[, 1L]
  if (ncol(means) > 1L) {
    for (j in seq.int(2L, ncol(means))) max_mean <- pmax(max_mean, means[, j], na.rm = TRUE)
  }
  names(max_mean) <- rownames(means)
  # The DE table carries its identifier in a column (its row names are the
  # positional ones of the data.frame), so match on that column.
  id_col <- cfg$input$gene_id_column %||% "gene_id"
  de_ids <- if (id_col %in% names(de)) as.character(de[[id_col]]) else rownames(de)
  common <- intersect(de_ids, names(max_mean))
  info$matched <- length(common)
  if (length(common) < max(10L, floor(0.5 * nrow(de)))) {
    info$reason <- "the expression matrix and the DE table share too few identifiers"
    return(list(de = de, info = info, max_mean = max_mean))
  }
  keep <- common[!is.na(max_mean[common]) & max_mean[common] >= min_fpkm]
  kept_de <- de[de_ids %in% keep, , drop = FALSE]
  info$applied <- TRUE
  info$kept <- nrow(kept_de)
  info$dropped <- nrow(de) - nrow(kept_de)
  info$reason <- NULL
  list(de = kept_de, info = info, max_mean = max_mean)
}

write_expression_floor <- function(cfg, floor, output_dir) {
  if (is.null(floor$max_mean) || is.null(output_dir)) return(invisible(NULL))
  tab <- data.frame(gene_id = names(floor$max_mean),
                    max_group_mean_fpkm = as.numeric(floor$max_mean),
                    kept = as.numeric(floor$max_mean) >= floor$info$min_fpkm,
                    stringsAsFactors = FALSE)
  tab <- tab[order(-tab$max_group_mean_fpkm), , drop = FALSE]
  path <- file.path(output_dir, "expression_floor.csv")
  utils::write.csv(tab, path, row.names = FALSE)
  invisible(path)
}

build_context <- function(cfg, counts_obj, meta, deseq_obj) {
  group_col <- cfg$analysis$group_column %||% "condition"
  groups <- levels(meta[[group_col]])
  contrast_cfg <- cfg$analysis$contrasts[[1L]]
  contrast_name <- contrast_cfg$name
  if (is.null(deseq_obj$contrast_results[[contrast_name]])) {
    stop("The DESeq2 fit holds no result for contrast '", contrast_name,
         "'. It was fitted under a different contrast label; remove the .cache ",
         "folder of the output directory and rerun.", call. = FALSE)
  }
  de <- deseq_obj$contrast_results[[contrast_name]]$result
  # The original scripts work on a per-contrast result table with the columns
  # gene_id / GeneSymbol / baseMean / log2FoldChange / padj.
  de <- de[!is.na(de$padj), , drop = FALSE]
  ids <- build_identifier_context(cfg, counts_obj$annotations)
  expr <- build_expression_matrix(
    cfg, counts_obj,
    output_dir = file.path(cfg$project$output_dir %||% ".", "00_differential_expression"))
  # The shared expression floor is applied once, before any section sees the
  # contrast table. The unfiltered table and the per-gene FPKM values are kept
  # as well: the screening figures draw the genes the floor removed in grey, so
  # a reader can see how much of the tested universe sits below the cut-off.
  de_all <- de
  floor <- apply_expression_floor(cfg, de, expr$matrix, expr$unit, meta,
                                  group_col %||% "condition")
  write_status(file.path(cfg$project$output_dir, "00_differential_expression", "status_expression_floor.json"),
               "expression_floor", if (isTRUE(floor$info$applied)) "completed" else "skipped",
               floor$info$reason, floor$info)
  if (is.finite(floor$info$min_fpkm) && floor$info$min_fpkm > 0 && !isTRUE(floor$info$applied)) {
    stop("Cannot apply min_fpkm: ", floor$info$reason,
         ". Supply matching gene lengths/GTF or an FPKM matrix; set min_fpkm = NULL to explicitly disable the floor.", call. = FALSE)
  }
  de <- floor$de
  write_expression_floor(cfg, floor, file.path(cfg$project$output_dir %||% ".",
                                               "00_differential_expression"))
  ch <- deseq_obj$contrast_results[[contrast_name]]
  ctx <- list(
    cfg = cfg,
    group_col = group_col,
    groups = groups,                 # levels: reference first
    numerator = ch$numerator,
    denominator = ch$denominator,
    contrast_name = contrast_name,
    de = de,                         # DESeq2 table for the contrast
    de_all = de_all,                 # the same table before the expression floor
    expr_fpkm = floor$max_mean,      # max group-mean FPKM per expression-matrix row
    min_fpkm = floor$info$min_fpkm,  # the floor itself (NULL when disabled)
    symbols = ids$symbol_col,
    key_type = ids$key_type,
    symbol_type = ids$symbol_type,
    gene_map = ids$gene_map,
    annotations = counts_obj$annotations,
    count_mat = counts_obj$counts,
    normalized = deseq_obj$normalized_counts,
    vst = SummarizedExperiment::assay(deseq_obj$vsd),
    expr_matrix = expr$matrix,
    expr_source = expr$source,
    expr_unit = expr$unit,
    expr_floor = floor$info,
    meta = meta,
    dds = deseq_obj$dds,
    deseq = deseq_obj,
    regulons = NULL,
    # filled in by section 1 and consumed by sections 2-4 and 9
    mylist = NULL,
    gene_list = NULL,                # ranked Entrez vector
    gene_df = NULL,                  # symbol / logFC / key frame
    sig_symbols = NULL
  )
  if (any(as.integer(unlist(cfg$sections$run)) %in% c(5L, 6L, 7L, 8L, 9L))) {
    ctx$regulons <- load_tf_regulons(cfg)
  }
  ctx
}

# ---------------------------------------------------------------------------
# Multi-threshold screening, shared between the screening section and the driver.
#
# The five-way classification, the threshold ladder and the per-threshold DEG
# tables are a property of the differential-expression result rather than of the
# section that draws them. Section 1 calls this to build its tables and figures;
# the driver calls it when a run starts after section 1, so that sections 2-11
# stay independently runnable. Both callers therefore see identical tables.

screening_tables <- function(cfg, ctx) {
  de <- ctx$de
  symbol_col <- ctx$symbols
  gene_id_col <- cfg$input$gene_id_column

  signafi <- -log10(de$padj)
  fold_change_raw <- de$log2FoldChange
  # All tabulated and computed effects remain unbounded.
  fold_change <- fold_change_raw
  base_mean <- de$baseMean
  exp10 <- log10(pmax(base_mean, .Machine$double.xmin))
  gene_symbol <- as.character(de[[symbol_col]])
  gene_id <- as.character(de[[gene_id_col]])

  thresholds <- sort(unique(c(as.numeric(unlist(cfg$analysis$volcano_lfc_thresholds %||%
                                    c(0.15, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5))),
                              as.numeric(cfg$analysis$lfc_threshold %||% 0.75),
                              as.numeric(unlist(cfg$analysis$chromosome_lfc_thresholds %||%
                                                  c(0.5, 0.75, 1))))))
  # AnaDGEsS1.R hid every gene with log10(baseMean) <= 2 from the up/down
  # classes. That removes genuinely strong low-abundance DEGs (here 941 genes,
  # mostly the most down-regulated ones), so the filter is off by default and
  # can be switched back on by setting analysis.expression_log10_threshold: 2.
  expression_threshold <- cfg$analysis$expression_log10_threshold
  use_expr_filter <- !is.null(expression_threshold) && is.finite(as.numeric(expression_threshold))
  if (use_expr_filter) expression_threshold <- as.numeric(expression_threshold)
  alpha <- as.numeric(cfg$analysis$alpha %||% 0.05)
  sig_cut <- -log10(alpha)

  classify <- function(limit) {
    ok <- !is.na(signafi)
    significant <- ok & de$padj < alpha
    low <- if (use_expr_filter) exp10 <= expression_threshold else rep(FALSE, length(signafi))
    cls <- rep("intermediate", length(signafi))
    cls[low] <- "low_expression"
    cls[!low & !significant] <- "not_significant"
    cls[!low & significant & fold_change_raw >= limit] <- "up"
    cls[!low & significant & fold_change_raw <= -limit] <- "down"
    cls
  }

  mylist <- vector("list", length(thresholds))
  class_summary <- vector("list", length(thresholds))
  for (i in seq_along(thresholds)) {
    limit <- thresholds[i]
    cls <- classify(limit)
    up <- cls == "up"
    down <- cls == "down"
    mylist[[i]] <- data.frame(
      GenesID = c(gene_id[up], gene_id[down]),
      GenesSymbol = c(gene_symbol[up], gene_symbol[down]),
      GenesbaseMean = c(base_mean[up], base_mean[down]),
      GenesFoldChange = c(fold_change[up], fold_change[down]),
      GenesFoldChange_raw = c(fold_change_raw[up], fold_change_raw[down]),
      GenesQvalue = c(de$padj[up], de$padj[down]),
      Direction = c(rep("Up", sum(up)), rep("Down", sum(down))),
      stringsAsFactors = FALSE
    )
    class_summary[[i]] <- data.frame(
      lfc_limit = limit,
      low_expression = sum(cls == "low_expression"),
      not_significant = sum(cls == "not_significant"),
      intermediate = sum(cls == "intermediate"),
      up = sum(up), down = sum(down), total_degs = sum(up) + sum(down)
    )
  }
  names(mylist) <- sprintf("lfc_%.2f", thresholds)
  list(mylist = mylist,
       class_summary = do.call(rbind, class_summary),
       thresholds = thresholds,
       fold_change = fold_change, fold_change_raw = fold_change_raw,
       base_mean = base_mean, exp10 = exp10,
       signafi = signafi, gene_symbol = gene_symbol, gene_id = gene_id,
       use_expr_filter = use_expr_filter, expression_threshold = expression_threshold,
       sig_cut = sig_cut, classify = classify)
}

# Keep the tested universe separate from the GSEA-only pre-screen. No display
# value enters either frame. Duplicate identifiers retain the largest absolute
# finite effect, identically for positive and negative effects.
attach_ranked_context <- function(cfg, ctx, tab) {
  de <- ctx$de
  symbols <- as.character(de[[ctx$symbols]])
  ok <- !is.na(symbols) & nzchar(symbols) & is.finite(de$log2FoldChange)
  gene_df <- data.frame(SYMBOL = symbols[ok], logFC = de$log2FoldChange[ok],
                        padj = de$padj[ok], baseMean = de$baseMean[ok],
                        gene_id = as.character(de[[cfg$input$gene_id_column]][ok]),
                        stringsAsFactors = FALSE)
  if (!is.null(ctx$gene_map)) {
    gene_df <- merge(gene_df, unique(ctx$gene_map), by.x = "SYMBOL",
                     by.y = ctx$symbol_type, all.x = TRUE, sort = FALSE)
  } else {
    gene_df[[ctx$key_type]] <- gene_df$SYMBOL
  }
  expr_thr <- cfg$analysis$expression_log10_threshold
  if (!is.null(expr_thr)) {
    keep <- is.finite(gene_df$baseMean) &
      log10(pmax(gene_df$baseMean, .Machine$double.xmin)) > expr_thr
    gene_df <- gene_df[keep, , drop = FALSE]
  }
  gene_df <- gene_df[order(-abs(gene_df$logFC), gene_df$SYMBOL), , drop = FALSE]
  low <- as.numeric(cfg$analysis$low_cutoff %||% 0)
  eligible <- is.finite(gene_df$padj) & gene_df$padj < cfg$analysis$alpha &
    abs(gene_df$logFC) >= low
  gsea_df <- gene_df[eligible, , drop = FALSE]
  ctx$gene_df <- gene_df
  ctx$gsea_gene_df <- gsea_df
  ctx$gene_list <- make_unique_ranked_vector(gsea_df[[ctx$key_type]], gsea_df$logFC)
  ctx$universe_list <- make_unique_ranked_vector(gene_df[[ctx$key_type]], gene_df$logFC)
  ctx$sig_symbols <- unique(gene_df$SYMBOL[is.finite(gene_df$padj) &
    gene_df$padj < cfg$analysis$alpha & abs(gene_df$logFC) >= cfg$analysis$lfc_threshold])
  provenance_dir <- safe_dir(file.path(cfg$project$output_dir, "provenance"))
  safe_write_csv(gsea_df, file.path(provenance_dir, "gsea_input_genes.csv"))
  safe_write_csv(data.frame(identifier = names(ctx$gene_list), log2FoldChange = unname(ctx$gene_list)),
                 file.path(provenance_dir, "gsea_ranked_vector.csv"))
  write_status(file.path(provenance_dir, "status_gsea_input.json"), "gsea_input", "completed",
               details = list(statistic = "unshrunk log2FoldChange; no upper clipping",
                              gene_p_column = "padj", alpha = cfg$analysis$alpha,
                              low_cutoff = low, input_rows = nrow(gsea_df),
                              ranked_identifiers = length(ctx$gene_list),
                              tested_identifiers = length(ctx$universe_list)))
  ctx
}

# ---------------------------------------------------------------------------
# Sections hand results to each other through the context object, which only
# exists inside one run. Reloading the published tables from the output
# directory keeps a later section runnable on its own (for example section 7,
# which consumes the table section 5 writes). Only empty slots are filled.

hydrate_context_from_disk <- function(cfg, ctx, output_dir) {
  contrast <- ctx$contrast_name
  read_out <- function(...) {
    path <- file.path(output_dir, ...)
    if (file.exists(path)) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) else NULL
  }
  if (is.null(ctx$go_ora) && !2L %in% unlist(cfg$sections$run)) {
    ora <- list()
    for (ont in c("BP", "CC", "MF")) {
      tab <- read_out("section02_go_ora", paste0("ENRICHGO", ont, "significant_", contrast, ".csv"))
      if (is.null(tab)) {
        tab <- read_out("section02_go_ora", paste0("ENRICHGO", ont, "alls_", contrast, ".csv"))
        if (!is.null(tab) && nrow(tab)) {
          p <- if (ont == "BP") cfg$analysis$go_bp_pvalue_cutoff else cfg$analysis$go_pvalue_cutoff
          q <- if (ont == "BP") cfg$analysis$go_bp_qvalue_cutoff else cfg$analysis$go_qvalue_cutoff
          keep <- !is.na(tab$pvalue) & tab$pvalue <= p & !is.na(tab$p.adjust) & tab$p.adjust <= p
          if ("qvalue" %in% names(tab) && any(!is.na(tab$qvalue))) keep <- keep & !is.na(tab$qvalue) & tab$qvalue <= q
          tab <- tab[keep, , drop = FALSE]
        }
      }
      if (!is.null(tab) && nrow(tab)) ora[[ont]] <- tab
    }
    if (length(ora)) ctx$go_ora <- ora
  }
  if (is.null(ctx$tf_summary) && !5L %in% unlist(cfg$sections$run)) {
    ctx$tf_summary <- read_out("section05_tf_regulon", paste0("TF_regulon_summary_", contrast, ".csv"))
  }
  if (is.null(ctx$tf_edges) && !5L %in% unlist(cfg$sections$run)) {
    ctx$tf_edges <- read_out("section05_tf_regulon", paste0("TF_target_edges_", contrast, ".csv"))
  }
  if (is.null(ctx$tf_gsea) && !5L %in% unlist(cfg$sections$run)) {
    ctx$tf_gsea <- read_out("section05_tf_regulon", paste0("TF_regulon_GSEA_", contrast, ".csv"))
  }
  ctx
}
