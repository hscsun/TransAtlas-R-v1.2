# Machine-readable and human-readable provenance for each run.

write_provenance <- function(cfg, project_root, counts_obj, meta, out_dir, de_obj = NULL) {
  safe_dir(out_dir)
  input_files <- c(cfg$.config_path, cfg$input$counts_file, cfg$input$sample_metadata_file)
  input_files <- input_files[!is.na(input_files) & file.exists(input_files)]
  checksums <- stats::setNames(as.character(tools::md5sum(input_files)), input_files)
  contrast_summary <- list()
  if (!is.null(de_obj)) {
    contrast_summary <- lapply(de_obj$contrast_results, function(x) {
      d <- x$result
      list(n_genes = nrow(d), n_tested = sum(!is.na(d$padj)),
           n_significant = sum(d$significant, na.rm = TRUE),
           n_up = sum(d$significant & d$log2FoldChange > 0, na.rm = TRUE),
           n_down = sum(d$significant & d$log2FoldChange < 0, na.rm = TRUE),
           shrinkage_method = x$shrinkage_method)
    })
  }
  manifest <- list(
    software_version = trimws(readLines(file.path(project_root, "VERSION"), warn = FALSE)[1L]),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    project_root = normalizePath(project_root, winslash = "/"),
    input = list(counts_file = cfg$input$counts_file, sheet = counts_obj$sheet,
                 n_genes = nrow(counts_obj$counts), n_samples = ncol(counts_obj$counts),
                 sample_names = colnames(counts_obj$counts), metadata = meta),
    species = cfg$species,
    analysis = cfg$analysis,
    input_md5 = checksums,
    contrasts = contrast_summary,
    session_info = capture.output(utils::sessionInfo())
  )
  jsonlite::write_json(manifest, file.path(out_dir, "run_manifest.json"), auto_unbox = TRUE,
                       pretty = TRUE, na = "null", dataframe = "rows")
  writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  invisible(manifest)
}
