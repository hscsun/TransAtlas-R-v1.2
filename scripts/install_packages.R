# Install every R package the workflow can use.
#
# The list comes from pipeline_packages() in pipeline/R/helpers.R, the same
# function run_all() checks at start-up, so this script cannot drift away from
# what the analysis actually loads. Installations are routed by repository:
# Bioconductor packages (DESeq2, GSVA, dorothea, pathview, the OrgDb packages,
# ...) go through BiocManager, so their Bioconductor dependencies such as
# Rhdf5lib, rhdf5, HDF5Array and SpatialExperiment resolve with the matching
# release; CRAN packages go through install.packages().
#
#   Rscript scripts/install_packages.R
#   Rscript scripts/install_packages.R --species human
#
# The command is safe to repeat: installed packages are skipped.

root <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."),
                         winslash = "/", mustWork = FALSE))
  }
  normalizePath(".", winslash = "/", mustWork = FALSE)
})
source(file.path(root, "pipeline", "R", "helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
species <- "mouse"
i <- match("--species", args)
if (!is.na(i) && length(args) >= i + 1L) species <- args[i + 1L]

plan <- pipeline_packages(species)
pkgs <- unique(c(plan$required, plan$optional))
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

message("Package plan for ", species, ": ", length(pkgs), " packages, ",
        length(missing), " missing.")
if (!length(missing)) {
  message("Nothing to install.")
  quit(save = "no", status = 0L)
}
message("Missing: ", paste(missing, collapse = ", "))

res <- install_r_packages(missing)
if (length(res$failed)) {
  message("\nStill missing: ", paste(res$failed, collapse = ", "))
  if (length(res$notes)) message(paste0("  ", res$notes, collapse = "\n"))
  message("Bioconductor builds need a C/C++ toolchain and their system libraries; ",
          "run BiocManager::install(c(\"", paste(res$failed, collapse = "\", \""),
          "\"), dependencies = TRUE) to see the full build log.")
  quit(save = "no", status = 1L)
}
message("\nAll packages available. Use sessionInfo() to record the exact versions.")
