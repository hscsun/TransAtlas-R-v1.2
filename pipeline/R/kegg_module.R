# KEGG pathway gene sets, cached on disk.
#
# gseKEGG() downloads the organism's <pathway, gene> link table and the
# KEGG <-> Entrez conversion table on every call and keeps them only in a
# session-local environment, so a new R process downloads them again. Both
# tables change only when KEGG publishes a release, and on a slow connection
# the download exceeds R's default 60 s timeout - which is how section 4 came
# to write nothing but a status file. The tables are therefore fetched once
# into cache/kegg/ (with a generous timeout), and every later run reads them
# from disk and never touches the network.
#
# The gene half of a KEGG link pair is "mmu:103988", i.e. the organism prefix
# plus the Entrez GeneID, which is exactly the identifier the ranked list uses,
# so no separate KEGG-to-Entrez conversion is needed.

kegg_rest_url <- "https://rest.kegg.jp"

#' Seconds allowed for one KEGG table download before giving up.
#'
#' rest.kegg.jp serves the full pathway link table slowly, so the default is far
#' above R's 60 s; it is bounded all the same, because a run that blocks forever
#' on an unreachable mirror is worse than a run that falls back and says so.
kegg_download_timeout <- function(cfg) {
  as.numeric(cfg$runtime$kegg_download_timeout %||% 300)
}

#' Marker recording that the KEGG download failed recently.
#'
#' A run that fell back to the MSigDB KEGG collection should not pay the
#' download timeout again on every subsequent run, so the failure is remembered
#' for a day. Pass run_all(kegg_download_timeout = ...) with a larger value, or
#' delete the marker, to retry sooner.
kegg_download_blocked <- function(dir, code, within_hours = 24) {
  marker <- file.path(dir, paste0(code, "_download_failed.txt"))
  if (!file.exists(marker)) return(FALSE)
  age <- as.numeric(difftime(Sys.time(), file.info(marker)$mtime, units = "hours"))
  age < within_hours
}

kegg_mark_download_failed <- function(dir, code, reason) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(paste0("KEGG pathway table download failed at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               paste0("Organism: ", code), paste0("Reason: ", reason),
               "The section fell back to the MSigDB KEGG collection.",
               "Delete this file (or raise runtime.kegg_download_timeout) to retry."),
             file.path(dir, paste0(code, "_download_failed.txt")))
}

#' Directory holding the cached KEGG tables.
kegg_cache_dir <- function(cfg) {
  configured <- cfg$runtime$kegg_cache_dir %||% NULL
  if (!is.null(configured) && nzchar(configured)) return(configured)
  file.path(attr(cfg, "project_root") %||% ".", "cache", "kegg")
}

#' Download a KEGG REST endpoint into the cache unless it is already there.
#'
#' Returns TRUE when the file is present and non-empty afterwards. The timeout
#' is raised for the duration of the call because rest.kegg.jp serves the full
#' pathway table slowly, and a truncated download is discarded rather than
#' cached.
kegg_cached_download <- function(url, dest, timeout = 600) {
  if (file.exists(dest) && file.info(dest)$size > 0) return(TRUE)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  old <- getOption("timeout")
  on.exit(options(timeout = old), add = TRUE)
  options(timeout = timeout)
  tmp <- paste0(dest, ".part")
  unlink(tmp)
  ok <- tryCatch({
    utils::download.file(url, tmp, quiet = TRUE, mode = "wb")
    file.exists(tmp) && file.info(tmp)$size > 0
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (isTRUE(ok)) {
    file.rename(tmp, dest)
  } else {
    unlink(tmp)
  }
  file.exists(dest) && file.info(dest)$size > 0
}

#' Read a cached tab-delimited KEGG table.
kegg_read_table <- function(path) {
  dat <- tryCatch(utils::read.table(path, sep = "\t", quote = "", comment.char = "",
                                    stringsAsFactors = FALSE, fill = TRUE),
                  error = record_optional_error)
  if (is.null(dat) || ncol(dat) < 2L || nrow(dat) == 0L) return(NULL)
  dat <- dat[, 1:2, drop = FALSE]
  dat <- dat[!is.na(dat[[1L]]) & !is.na(dat[[2L]]), , drop = FALSE]
  dat[nzchar(dat[[1L]]) & nzchar(dat[[2L]]), , drop = FALSE]
}

#' Pathway gene sets for the configured organism.
#'
#' Returns list(term2gene, term2name, source) ready for clusterProfiler::GSEA(),
#' or NULL when neither the cache nor the network can supply the tables.
#' `refresh = TRUE` forces a fresh download.
kegg_pathway_sets <- function(cfg, refresh = FALSE, timeout = NULL) {
  code <- cfg$species$kegg_code
  if (is.null(code) || !nzchar(code)) return(NULL)
  if (is.null(timeout)) timeout <- kegg_download_timeout(cfg)
  dir <- kegg_cache_dir(cfg)
  link_file <- file.path(dir, paste0(code, "_pathway_link.txt"))
  list_file <- file.path(dir, paste0(code, "_pathway_list.txt"))
  if (isTRUE(refresh)) unlink(c(link_file, list_file))

  had_cache <- file.exists(link_file) && file.exists(list_file)
  if (!had_cache && kegg_download_blocked(dir, code)) {
    # A recent attempt already timed out; do not pay the timeout again.
    return(NULL)
  }
  # One attempt on the endpoint clusterProfiler itself uses. A second URL form
  # would double the wait for the same table on the same server, and a run that
  # spends ten minutes failing to download is worse than one that falls back in
  # minutes and says which source it used.
  ok_link <- kegg_cached_download(paste0(kegg_rest_url, "/link/", code, "/pathway"),
                                  link_file, timeout = timeout)
  if (!ok_link) {
    kegg_mark_download_failed(dir, code, paste0("no response within ", timeout, " s"))
    return(NULL)
  }
  # The pathway names are small; a failure there leaves the ids as descriptions
  # rather than losing the whole analysis.
  ok_list <- kegg_cached_download(paste0(kegg_rest_url, "/list/pathway/", code),
                                  list_file, timeout = min(timeout, 120))

  link <- kegg_read_table(link_file)
  if (is.null(link)) return(NULL)
  # Drop the "mmu:" / "path:" prefix from both halves, exactly as
  # clusterProfiler's download.KEGG.Path does.
  link[[1L]] <- sub("^[^:]+:", "", link[[1L]])
  link[[2L]] <- sub("^[^:]+:", "", link[[2L]])

  names_df <- if (ok_list) kegg_read_table(list_file) else NULL
  if (!is.null(names_df)) {
    names_df[[1L]] <- sub("^[^:]+:", "", names_df[[1L]])
    # The list endpoint appends " - Mus musculus (house mouse)".
    names_df[[2L]] <- sub("\\s-\\s[^-]+$", "", names_df[[2L]])
    link <- link[link[[1L]] %in% names_df[[1L]], , drop = FALSE]
  }
  if (nrow(link) == 0L) return(NULL)

  term2gene <- data.frame(term = link[[1L]], gene = link[[2L]], stringsAsFactors = FALSE)
  term2gene <- term2gene[!duplicated(term2gene), , drop = FALSE]
  term2name <- if (!is.null(names_df)) {
    data.frame(term = names_df[[1L]], name = names_df[[2L]], stringsAsFactors = FALSE)
  } else {
    data.frame(term = unique(term2gene$term), name = unique(term2gene$term),
               stringsAsFactors = FALSE)
  }
  term2name <- term2name[!duplicated(term2name$term), , drop = FALSE]

  list(term2gene = term2gene, term2name = term2name,
       source = if (had_cache) paste0("KEGG cache (", code, ")")
                else paste0("downloaded from ", basename(kegg_rest_url), " (", code, ")"),
       cache_dir = dir)
}

#' KEGG gene sets from the bundled MSigDB data, as an offline fallback.
#'
#' Used only when rest.kegg.jp cannot deliver the organism's pathway table.
#' MSigDB's KEGG collection covers the same pathways under MSigDB names, and the
#' Entrez identifiers come from the local msigdbr data, so the enrichment itself
#' is unchanged; only the pathway labels differ, and the source is recorded.
msigdb_kegg_sets <- function(cfg) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) return(NULL)
  species <- cfg$species$name %||% "Mus musculus"
  db_species <- cfg$analysis$msigdb_db_species %||% "MM"
  grab <- function(sub) {
    tryCatch(msigdbr::msigdbr(species = species, db_species = db_species,
                              collection = "C2", subcollection = sub),
             error = function(e) tryCatch(msigdbr::msigdbr(species = species, collection = "C2",
                                                           subcollection = sub),
                                          error = function(e2) NULL))
  }
  dat <- grab("CP:KEGG_LEGACY")
  if (is.null(dat) || nrow(dat) == 0L) dat <- grab("CP:KEGG_MEDICUS")
  if (is.null(dat) || nrow(dat) == 0L) return(NULL)
  gene_col <- if ("ncbi_gene" %in% names(dat)) "ncbi_gene" else "gene_symbol"
  term2gene <- data.frame(term = as.character(dat$gs_name),
                          gene = as.character(dat[[gene_col]]),
                          stringsAsFactors = FALSE)
  term2gene <- term2gene[!is.na(term2gene$gene) & nzchar(term2gene$gene), , drop = FALSE]
  term2gene <- term2gene[!duplicated(term2gene), , drop = FALSE]
  if (nrow(term2gene) == 0L) return(NULL)
  # The MSigDB set name is the identifier here ("KEGG_CELL_CYCLE"), which keeps
  # the provenance visible in the KEGG table next to a readable description.
  pretty_name <- function(x) {
    x <- gsub("_", " ", sub("^KEGG_", "", x))
    # Capitalise each word without lowercasing the rest, so DNA and G2M survive.
    paste0(gsub("(^| )([a-z])", "\\1\\U\\2", x, perl = TRUE), " (MSigDB KEGG)")
  }
  term2name <- data.frame(term = unique(term2gene$term), stringsAsFactors = FALSE)
  term2name$name <- pretty_name(term2name$term)
  list(term2gene = term2gene, term2name = term2name,
       source = paste0("MSigDB KEGG collection via msigdbr (offline fallback; ",
                       gene_col, " identifiers)"))
}

#' Run KEGG gene-set enrichment without gseKEGG's per-call download.
#'
#' gseKEGG() is equivalent to GSEA() over the organism's pathway gene sets, so
#' the same analysis runs from the cached tables. Falls back to the MSigDB KEGG
#' collection when neither the cache nor the network can supply the KEGG table;
#' when a fallback is used the source is attached to the result and reported, so
#' a reader can never mistake it for a live KEGG download.
run_kegg_gsea <- function(cfg, rank, kegg_id_type, alpha, min_size, max_size,
                          timeout = NULL) {
  if (is.null(timeout)) timeout <- kegg_download_timeout(cfg)
  # fgsea normalises the enrichment score against random gene sets, so the
  # reported NES depends on the state of the RNG. Re-seeding here makes it
  # depend on the configured seed alone rather than on how many random draws the
  # preceding sections happened to make.
  set.seed(as.integer(cfg$project$seed %||% 1L))

  run_gsea <- function(sets) {
    res <- tryCatch(clusterProfiler::GSEA(
      geneList = rank,
      TERM2GENE = sets$term2gene,
      TERM2NAME = sets$term2name,
      minGSSize = min_size, maxGSSize = max_size,
      pvalueCutoff = alpha, pAdjustMethod = "BH", eps = 0, verbose = FALSE),
      error = record_optional_error)
    if (!is.null(res) && nrow(as.data.frame(res)) > 0L) {
      attr(res, "kegg_source") <- sets$source
      return(res)
    }
    NULL
  }

  sets <- tryCatch(kegg_pathway_sets(cfg, timeout = timeout), error = record_optional_error)
  if (!is.null(sets)) {
    res <- run_gsea(sets)
    if (!is.null(res)) return(res)
  }
  # gseKEGG() downloads the same tables that just failed, so going through it
  # would only repeat the timeout.
  if (!kegg_download_blocked(kegg_cache_dir(cfg), cfg$species$kegg_code)) {
    old <- getOption("timeout")
    on.exit(options(timeout = old), add = TRUE)
    options(timeout = timeout)
    res <- tryCatch(clusterProfiler::gseKEGG(
      geneList = rank, organism = cfg$species$kegg_code, keyType = kegg_id_type,
      minGSSize = min_size, maxGSSize = max_size, pvalueCutoff = alpha,
      pAdjustMethod = "BH", eps = 0, verbose = FALSE),
      error = record_optional_error)
    if (!is.null(res) && nrow(as.data.frame(res)) > 0L) {
      attr(res, "kegg_source") <- "gseKEGG (live download)"
      return(res)
    }
  }
  fallback <- tryCatch(msigdb_kegg_sets(cfg), error = record_optional_error)
  if (!is.null(fallback)) run_gsea(fallback) else NULL
}
