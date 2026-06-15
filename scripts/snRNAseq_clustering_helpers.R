# Normalization + Harmony + UMAP + clustering helpers (sourced from *\_full_pipeline.Rmd)

snRNAseq_clustering_init <- function(repo_root = "C:/Yale") {
  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(harmony)
    library(future)
  })
  options(future.globals.maxSize = 8 * 1024^3)
  plan("sequential")

  helpers <- file.path(repo_root, "scripts", "snRNAseq_qc_filter_helpers.R")
  if (!exists("threshold_folder_label", mode = "function")) {
    source(helpers)
  }

  assign("use_glmGamPoi", requireNamespace("glmGamPoi", quietly = TRUE), envir = .GlobalEnv)
  if (use_glmGamPoi) {
    message("glmGamPoi detected — SCTransform will use vst.flavor = 'v2'.")
  } else {
    message("glmGamPoi not found — SCTransform will use vst.flavor = 'v1' (install glmGamPoi for speed).")
  }
  invisible(TRUE)
}

sct_params <- function() {
  list(
    vst.flavor = if (exists("use_glmGamPoi") && use_glmGamPoi) "v2" else "v1",
    variable.features.n = 3000,
    verbose = FALSE
  )
}

save_umap_plots <- function(obj, plot_dir, prefix, cluster_col = "seurat_clusters", res_tag = "") {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- if (nzchar(res_tag)) paste0("_", res_tag) else ""

  p_clusters <- DimPlot(obj, reduction = "umap", group.by = cluster_col, label = TRUE) +
    NoLegend() +
    ggtitle(paste(prefix, "— clusters", if (nzchar(res_tag)) paste("res", res_tag) else ""))
  ggsave(
    file.path(plot_dir, paste0("umap_clusters", suffix, ".png")),
    p_clusters, width = 9, height = 7, dpi = 150
  )

  p_sample <- DimPlot(obj, reduction = "umap", group.by = "orig.ident") +
    ggtitle(paste(prefix, "— sample"))
  ggsave(
    file.path(plot_dir, paste0("umap_by_sample", suffix, ".png")),
    p_sample, width = 10, height = 7, dpi = 150
  )

  if ("genotype" %in% colnames(obj@meta.data)) {
    p_genotype <- DimPlot(obj, reduction = "umap", group.by = "genotype") +
      ggtitle(paste(prefix, "— genotype"))
    ggsave(
      file.path(plot_dir, paste0("umap_by_genotype", suffix, ".png")),
      p_genotype, width = 9, height = 7, dpi = 150
    )
  }
}

get_cluster_column <- function(obj) {
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    return("seurat_clusters")
  }
  cols <- grep("_res\\.", colnames(obj@meta.data), value = TRUE)
  if (length(cols) == 0) stop("No cluster column found in metadata.")
  cols[length(cols)]
}

add_genotype_metadata <- function(obj, prefix) {
  genotype <- sub(".*_(wt|het|hom)_.*", "\\1", obj$orig.ident)
  genotype <- toupper(genotype)
  obj$genotype <- factor(genotype, levels = c("WT", "HET", "HOM"))
  obj
}

normalize_summary_df <- function(df) {
  required <- c(
    "object_prefix", "threshold_label", "normalization",
    "resolution", "n_cells", "n_clusters", "filtered_rds", "clustered_rds"
  )
  if (!"normalization" %in% names(df)) {
    df$normalization <- NA_character_
  }
  df[required]
}

run_clustering_pipeline <- function(
  raw_rds,
  th,
  out_dir,
  object_prefix,
  mt_pattern = "^mt-",
  harmony_dims = 1:30,
  resolutions = c(0.5, 2.5)
) {
  label <- threshold_folder_label(th)
  cl_dir <- file.path(out_dir, "clustering", label)
  dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(file.path(cl_dir, "clustering_summary.csv"))) {
    message("\n=== ", object_prefix, " | ", label, " — SKIP (already done) ===")
    return(normalize_summary_df(read.csv(
      file.path(cl_dir, "clustering_summary.csv"),
      stringsAsFactors = FALSE
    )))
  }

  message("\n=== ", object_prefix, " | ", label, " ===")

  seurat_all <- readRDS(raw_rds)
  if (!"percent.mt" %in% colnames(seurat_all@meta.data)) {
    seurat_all[["percent.mt"]] <- PercentageFeatureSet(seurat_all, pattern = mt_pattern)
  }

  seurat_filt <- filter_seurat_by_thresholds(seurat_all, th)
  seurat_filt <- add_genotype_metadata(seurat_filt, object_prefix)
  rm(seurat_all)
  gc()
  message("Cells after filter: ", ncol(seurat_filt))

  filtered_rds <- file.path(cl_dir, paste0(object_prefix, "_filtered.rds"))
  saveRDS(seurat_filt, filtered_rds)

  if (packageVersion("Seurat") >= "5.0.0") {
    seurat_filt <- JoinLayers(seurat_filt)
  }

  norm_method <- if (exists("use_glmGamPoi") && use_glmGamPoi) "SCTransform (glmGamPoi)" else "SCTransform (native)"
  seurat_filt <- tryCatch(
    {
      message("  SCTransform (", if (exists("use_glmGamPoi") && use_glmGamPoi) "glmGamPoi" else "native", ")...")
      do.call(SCTransform, c(list(object = seurat_filt), sct_params()))
    },
    error = function(e) {
      norm_method <<- "LogNormalize"
      message("  SCTransform failed; using LogNormalize. Reason: ", conditionMessage(e))
      seurat_filt <- NormalizeData(seurat_filt, verbose = FALSE)
      seurat_filt <- FindVariableFeatures(seurat_filt, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
      seurat_filt <- ScaleData(seurat_filt, verbose = FALSE)
      seurat_filt
    }
  )
  gc()

  seurat_filt <- RunPCA(seurat_filt, verbose = FALSE)
  seurat_filt <- RunHarmony(seurat_filt, group.by.vars = "orig.ident", verbose = FALSE)
  seurat_filt <- RunUMAP(seurat_filt, reduction = "harmony", dims = harmony_dims)
  seurat_filt <- FindNeighbors(seurat_filt, reduction = "harmony", dims = harmony_dims)

  summary_rows <- list()

  for (res in resolutions) {
    res_tag <- gsub("\\.", "", as.character(res))
    obj <- FindClusters(seurat_filt, resolution = res, verbose = FALSE)
    rds_path <- file.path(cl_dir, paste0(object_prefix, "_harmony_res", res, ".rds"))
    saveRDS(obj, rds_path)

    cluster_col <- get_cluster_column(obj)
    save_umap_plots(obj, cl_dir, label, cluster_col = cluster_col, res_tag = res_tag)

    n_clusters <- length(unique(obj@meta.data[[cluster_col]]))
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      object_prefix = object_prefix,
      threshold_label = label,
      normalization = norm_method,
      resolution = res,
      n_cells = ncol(obj),
      n_clusters = n_clusters,
      filtered_rds = filtered_rds,
      clustered_rds = rds_path,
      stringsAsFactors = FALSE
    )
    message("  res ", res, ": ", n_clusters, " clusters -> ", basename(rds_path))

    cluster_counts <- as.data.frame(table(obj@meta.data[[cluster_col]]))
    colnames(cluster_counts) <- c("cluster", "n_cells")
    write.csv(
      cluster_counts,
      file.path(cl_dir, paste0("cluster_sizes_res", res, ".csv")),
      row.names = FALSE
    )
  }

  summary_df <- do.call(rbind, summary_rows)
  summary_df <- normalize_summary_df(summary_df)
  write.csv(summary_df, file.path(cl_dir, "clustering_summary.csv"), row.names = FALSE)

  invisible(summary_df)
}
