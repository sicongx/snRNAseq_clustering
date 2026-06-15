# Resolve repository root (folder containing scripts/) from getwd() or Rmd location.
find_repo_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(6)) {
    if (file.exists(file.path(d, "scripts", "snRNAseq_qc_filter_helpers.R"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) {
      break
    }
    d <- parent
  }
  stop(
    "Could not find repo root. Open RStudio project at snRNAseq_median5sd_clustering/ ",
    "or setwd() to that folder before running."
  )
}
