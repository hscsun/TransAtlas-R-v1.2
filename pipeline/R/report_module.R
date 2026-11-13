# Render the English manuscript draft after an analysis run.

render_manuscript <- function(cfg, project_root, results_dir, out_dir) {
  safe_dir(out_dir)
  template <- file.path(project_root, "manuscript", "manuscript_draft.Rmd")
  if (!file.exists(template) || !requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("R Markdown template or rmarkdown package unavailable.", call. = FALSE)
  }
  # The report reads the settings it was generated from, not the template: when
  # run_all() supplies the data file, species and thresholds, the template's
  # placeholders are still null and anything that touches them would fail.
  # Writing the resolved configuration keeps the report consistent with the run
  # and makes it self-describing.
  cfg_snapshot <- cfg[setdiff(names(cfg), ".config_path")]
  cfg_path <- file.path(dirname(out_dir), "provenance", "config_used.yml")
  safe_dir(dirname(cfg_path))
  ok_yaml <- tryCatch({ yaml::write_yaml(cfg_snapshot, cfg_path); TRUE }, error = function(e) FALSE)
  if (!ok_yaml) cfg_path <- cfg$.config_path %||% cfg_path
  params <- list(results_dir = normalizePath(results_dir, winslash = "/", mustWork = FALSE),
                config_path = normalizePath(cfg_path, winslash = "/", mustWork = FALSE),
                project_title = cfg$report$title %||% cfg$project$name)
  # knitr rewrites figure paths relative to the report directory, so pandoc
  # needs that directory on its resource path to embed the figures into a
  # single self-contained HTML file.
  report_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  tryCatch({
    rmarkdown::render(template, output_dir = out_dir, params = params, quiet = TRUE,
                      envir = new.env(parent = globalenv()),
                      output_options = list(pandoc_args = c("--resource-path", report_dir)))
    write_status(file.path(out_dir, "status.json"), "report", "completed")
  }, error = function(e) {
    write_status(file.path(out_dir, "status.json"), "report", "failed", conditionMessage(e))
    stop(e)
  })
  invisible(TRUE)
}
