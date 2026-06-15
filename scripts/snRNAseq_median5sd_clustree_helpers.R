# Per-sample nFeature max (median + k * SD), SCT-only clustering, and clustree.
# Used by *_full_pipeline_median5sd_clustree.Rmd (does not alter original helpers).

compute_per_sample_nfeature_max <- function(
  md,
  nFeature_min = 200,
  sd_multiplier = 5,
  sample_col = "orig.ident",
  metric = "nFeature_RNA"
) {
  per_sample <- do.call(rbind, lapply(split(md, md[[sample_col]], drop = TRUE), function(sub) {
    xs <- sub[[metric]]
    med <- stats::median(xs, na.rm = TRUE)
    sdv <- stats::sd(xs, na.rm = TRUE)
    if (is.na(sdv)) {
      sdv <- 0
    }
    nFeature_max <- med + sd_multiplier * sdv
    data.frame(
      sample = unique(as.character(sub[[sample_col]])),
      nFeature_median = med,
      nFeature_sd = sdv,
      sd_multiplier = sd_multiplier,
      nFeature_max = nFeature_max,
      stringsAsFactors = FALSE
    )
  }))
  per_sample <- per_sample[order(per_sample$sample), , drop = FALSE]
  if (any(per_sample$nFeature_max <= nFeature_min)) {
    warning(
      "Some per-sample nFeature_max (median + ", sd_multiplier, " * SD) are <= nFeature_min (",
      nFeature_min, "). Those samples may retain very few cells."
    )
  }
  per_sample
}

make_median5sd_cell_filter <- function(
  mt_max,
  per_sample_thresholds,
  nFeature_min = 200,
  sd_multiplier = 5
) {
  if (is.null(mt_max) || length(mt_max) != 1L || is.na(mt_max) || !is.numeric(mt_max)) {
    stop("Set numeric mt_max before running.")
  }
  if (is.null(per_sample_thresholds) || !is.data.frame(per_sample_thresholds)) {
    stop("per_sample_thresholds must be a data.frame from compute_per_sample_nfeature_max().")
  }
  required <- c("sample", "nFeature_max")
  missing <- setdiff(required, colnames(per_sample_thresholds))
  if (length(missing) > 0) {
    stop("per_sample_thresholds missing columns: ", paste(missing, collapse = ", "))
  }
  list(
    mode = "median5sd",
    nFeature_min = as.numeric(nFeature_min),
    nFeature_max = NA_real_,
    sd_multiplier = as.numeric(sd_multiplier),
    mt_max = as.numeric(mt_max),
    per_sample_max = per_sample_thresholds
  )
}

threshold_folder_label_median5sd <- function(th) {
  paste0(
    "nFmin", th$nFeature_min,
    "_nFmax_median", th$sd_multiplier, "SD",
    "_mtmax", th$mt_max
  )
}

filter_seurat_by_median5sd <- function(seurat_all, th) {
  md <- seurat_all@meta.data
  if (is.null(th$per_sample_max)) {
    stop("cell_filter$per_sample_max is required for median5sd filtering.")
  }
  max_by_cell <- th$per_sample_max$nFeature_max[
    match(as.character(md$orig.ident), th$per_sample_max$sample)
  ]
  if (any(is.na(max_by_cell))) {
    missing_samples <- setdiff(unique(as.character(md$orig.ident)), th$per_sample_max$sample)
    stop("Missing per-sample nFeature_max for: ", paste(missing_samples, collapse = ", "))
  }
  keep <- md$nFeature_RNA > th$nFeature_min &
    md$nFeature_RNA < max_by_cell &
    md$percent.mt < th$mt_max &
    md$scDblFinder.class == "singlet"
  seurat_all[, rownames(md)[keep]]
}

summarize_median5sd_filter_impact <- function(md, th) {
  max_by_cell <- th$per_sample_max$nFeature_max[
    match(as.character(md$orig.ident), th$per_sample_max$sample)
  ]
  build_row <- function(sub, scope_label, sample_label, sample_max) {
    xs <- sub$nFeature_RNA
    in_win <- xs > th$nFeature_min & xs < sample_max
    pass_all_s <- in_win & sub$percent.mt < th$mt_max & sub$scDblFinder.class == "singlet"
    data.frame(
      scope = scope_label,
      sample = sample_label,
      nFeature_min = th$nFeature_min,
      nFeature_max = sample_max,
      mt_max = th$mt_max,
      n_cells_total = nrow(sub),
      n_cells_pass_nFeature_only = sum(in_win),
      n_cells_pass_all_filters = sum(pass_all_s),
      pct_pass_all = round(100 * sum(pass_all_s) / nrow(sub), 2),
      n_excluded_low_nFeature = sum(xs <= th$nFeature_min),
      n_excluded_high_nFeature = sum(xs >= sample_max),
      stringsAsFactors = FALSE
    )
  }

  overall_in_win <- md$nFeature_RNA > th$nFeature_min & md$nFeature_RNA < max_by_cell
  overall_pass <- overall_in_win &
    md$percent.mt < th$mt_max &
    md$scDblFinder.class == "singlet"

  overall_row <- data.frame(
    scope = "overall",
    sample = "all",
    nFeature_min = th$nFeature_min,
    nFeature_max = NA_real_,
    mt_max = th$mt_max,
    n_cells_total = nrow(md),
    n_cells_pass_nFeature_only = sum(overall_in_win),
    n_cells_pass_all_filters = sum(overall_pass),
    pct_pass_all = round(100 * sum(overall_pass) / nrow(md), 2),
    n_excluded_low_nFeature = sum(md$nFeature_RNA <= th$nFeature_min),
    n_excluded_high_nFeature = sum(md$nFeature_RNA >= max_by_cell),
    stringsAsFactors = FALSE
  )

  per_sample <- do.call(rbind, lapply(split(md, md$orig.ident, drop = TRUE), function(sub) {
    sample_label <- unique(as.character(sub$orig.ident))
    sample_max <- th$per_sample_max$nFeature_max[
      match(sample_label, th$per_sample_max$sample)
    ]
    build_row(sub, "per_sample", sample_label, sample_max)
  }))

  rbind(overall_row, per_sample)
}

run_median5sd_threshold_diagnostics <- function(md, out_dir, cell_filter) {
  diag_dir <- file.path(out_dir, "threshold_diagnostics")
  if (dir.exists(diag_dir)) {
    unlink(diag_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

  distribution <- summarize_nfeature_distribution(md)
  sd_mad_bounds <- summarize_nfeature_data_driven_candidates(md, multipliers = 2:5, floor_value = 200)
  round_grid <- summarize_nfeature_round_grid(md)

  write.csv(distribution, file.path(diag_dir, "nfeature_distribution.csv"), row.names = FALSE)
  write.csv(sd_mad_bounds, file.path(diag_dir, "nfeature_sd_mad_bounds.csv"), row.names = FALSE)
  write.csv(round_grid, file.path(diag_dir, "nfeature_round_grid.csv"), row.names = FALSE)

  write.csv(
    cell_filter$per_sample_max,
    file.path(diag_dir, "per_sample_nfeature_max_median5sd.csv"),
    row.names = FALSE
  )

  user_impact <- summarize_median5sd_filter_impact(md, cell_filter)
  write.csv(user_impact, file.path(diag_dir, "user_filter_impact.csv"), row.names = FALSE)

  list(
    diag_dir = diag_dir,
    distribution = distribution,
    sd_mad_bounds = sd_mad_bounds,
    round_grid = round_grid,
    user_impact = user_impact,
    per_sample_max = cell_filter$per_sample_max
  )
}

prepare_per_sample_max_for_facet <- function(per_sample_max, sample_levels) {
  out <- per_sample_max
  out$orig.ident <- factor(out$sample, levels = sample_levels)
  out
}

per_sample_max_box_segments <- function(per_sample_max, sample_levels, half_width = 0.45) {
  ps <- prepare_per_sample_max_for_facet(per_sample_max, sample_levels)
  idx <- as.numeric(ps$orig.ident)
  data.frame(
    orig.ident = ps$orig.ident,
    xmin = idx - half_width,
    xmax = idx + half_width,
    nFeature_max = ps$nFeature_max,
    stringsAsFactors = FALSE
  )
}

plot_qc_density_per_sample_max <- function(
  data,
  metric,
  title,
  xlab,
  per_sample_max,
  vline_min = NULL,
  xlim = NULL
) {
  data$orig.ident <- factor(data$orig.ident, levels = sort(unique(as.character(data$orig.ident))))

  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[metric]])) +
    ggplot2::geom_density(fill = "steelblue", alpha = 0.35, color = "steelblue", linewidth = 0.4) +
    ggplot2::facet_wrap(~orig.ident, scales = "free_y", ncol = 3) +
    ggplot2::labs(title = title, x = xlab, y = "Density") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))

  if (!is.null(vline_min)) {
    p <- p + ggplot2::geom_vline(
      xintercept = vline_min,
      color = "red",
      linetype = "dashed",
      linewidth = 0.7
    )
  }

  if (!is.null(per_sample_max) && nrow(per_sample_max) > 0) {
    ps_max <- prepare_per_sample_max_for_facet(per_sample_max, levels(data$orig.ident))
    p <- p + ggplot2::geom_vline(
      data = ps_max,
      ggplot2::aes(xintercept = nFeature_max),
      inherit.aes = FALSE,
      color = "darkred",
      linetype = "dotted",
      linewidth = 0.7
    )
  }

  if (!is.null(xlim) && length(xlim) == 2L) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }
  p
}

plot_qc_box_per_sample_max <- function(
  data,
  metric,
  title,
  ylab,
  per_sample_max,
  hline_min = NULL
) {
  data$orig.ident <- factor(data$orig.ident, levels = sort(unique(as.character(data$orig.ident))))

  p <- ggplot2::ggplot(data, ggplot2::aes(x = orig.ident, y = .data[[metric]], fill = orig.ident)) +
    ggplot2::geom_boxplot(outlier.size = 0.3) +
    ggplot2::labs(title = title, x = NULL, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip()

  if (!is.null(hline_min)) {
    p <- p + ggplot2::geom_hline(yintercept = hline_min, color = "red", linetype = "dashed", linewidth = 0.7)
  }

  if (!is.null(per_sample_max) && nrow(per_sample_max) > 0) {
    max_segments <- per_sample_max_box_segments(per_sample_max, levels(data$orig.ident))
    p <- p + ggplot2::geom_segment(
      data = max_segments,
      ggplot2::aes(x = xmin, xend = xmax, y = nFeature_max, yend = nFeature_max),
      inherit.aes = FALSE,
      color = "darkred",
      linetype = "dotted",
      linewidth = 0.7
    )
  }
  p
}

plot_qc_scatter_per_sample_max <- function(data, title, per_sample_max, nFeature_min = NULL) {
  data$orig.ident <- factor(data$orig.ident, levels = sort(unique(as.character(data$orig.ident))))

  p <- ggplot2::ggplot(data, ggplot2::aes(x = nFeature_RNA, y = nCount_RNA)) +
    ggplot2::geom_point(color = "steelblue", alpha = 0.08, size = 0.25) +
    ggplot2::facet_wrap(~orig.ident, scales = "free", ncol = 3) +
    ggplot2::labs(title = title, x = "nFeature_RNA", y = "nCount_RNA") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))

  if (!is.null(nFeature_min)) {
    p <- p + ggplot2::geom_vline(xintercept = nFeature_min, color = "red", linetype = "dashed", linewidth = 0.7)
  }

  if (!is.null(per_sample_max) && nrow(per_sample_max) > 0) {
    ps_max <- prepare_per_sample_max_for_facet(per_sample_max, levels(data$orig.ident))
    p <- p + ggplot2::geom_vline(
      data = ps_max,
      ggplot2::aes(xintercept = nFeature_max),
      inherit.aes = FALSE,
      color = "darkred",
      linetype = "dotted",
      linewidth = 0.7
    )
  }
  p
}

save_qc_exploratory_plots_median5sd <- function(
  md,
  th,
  out_dir,
  gene_label = "QC",
  density_xlim_nFeature = c(0, 6000),
  density_xlim_nCount = c(0, 15000)
) {
  plot_dir <- file.path(out_dir, "qc_plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  md$orig.ident <- factor(md$orig.ident, levels = sort(unique(md$orig.ident)))
  per_sample_max <- th$per_sample_max
  per_sample_max$orig.ident <- factor(per_sample_max$sample, levels = levels(md$orig.ident))

  save_gg(
    plot_qc_density_per_sample_max(
      md, "nFeature_RNA", "nFeature — density (zoom)", "nFeature_RNA",
      per_sample_max = per_sample_max,
      vline_min = th$nFeature_min,
      xlim = density_xlim_nFeature
    ),
    file.path(plot_dir, "density_nFeature"), width = 14, height = 10
  )
  save_gg(
    plot_qc_density(md, "nCount_RNA", "nCount — density (zoom)", "nCount_RNA",
                    NULL, NULL, xlim = density_xlim_nCount),
    file.path(plot_dir, "density_nCount"), width = 14, height = 10
  )
  save_gg(
    plot_qc_box_per_sample_max(
      md, "nFeature_RNA", "nFeature — boxplot", "nFeature_RNA",
      per_sample_max = per_sample_max,
      hline_min = th$nFeature_min
    ),
    file.path(plot_dir, "box_nFeature"), width = 10, height = 6
  )
  save_gg(
    plot_qc_box(md, "nCount_RNA", "nCount — boxplot", "nCount_RNA", NULL, NULL),
    file.path(plot_dir, "box_nCount"), width = 10, height = 6
  )
  save_gg(
    plot_qc_scatter_per_sample_max(
      md, "nFeature vs nCount — scatter",
      per_sample_max = per_sample_max,
      nFeature_min = th$nFeature_min
    ),
    file.path(plot_dir, "scatter_nFeature_nCount"), width = 14, height = 10
  )
  invisible(plot_dir)
}

plot_qc_box_zoom_vertical_per_sample_max <- function(
  data,
  metric,
  title,
  ylab,
  per_sample_max,
  y_min = 200,
  y_max = 3000,
  points_per_sample = 1500,
  hline_min = NULL
) {
  data$orig.ident <- factor(data$orig.ident, levels = sort(unique(as.character(data$orig.ident))))
  md_points <- subsample_meta_for_points(data, points_per_sample = points_per_sample)

  use_sina <- requireNamespace("ggforce", quietly = TRUE)
  if (use_sina) {
    point_layer <- ggforce::geom_sina(
      data = md_points,
      ggplot2::aes(x = orig.ident, y = .data[[metric]], color = orig.ident),
      alpha = 0.35,
      size = 0.35,
      maxwidth = 0.45
    )
  } else {
    point_layer <- ggplot2::geom_jitter(
      data = md_points,
      ggplot2::aes(x = orig.ident, y = .data[[metric]], color = orig.ident),
      width = 0.18,
      height = 0,
      alpha = 0.3,
      size = 0.3
    )
  }

  n_samples <- length(levels(data$orig.ident))
  sample_colors <- grDevices::hcl.colors(n_samples, palette = "Set2")
  per_sample_max$orig.ident <- factor(per_sample_max$sample, levels = levels(data$orig.ident))

  p <- ggplot2::ggplot(data, ggplot2::aes(x = orig.ident, y = .data[[metric]], fill = orig.ident)) +
    ggplot2::geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      alpha = 0.25,
      linewidth = 0.45,
      color = "grey30"
    ) +
    point_layer +
    ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +
    ggplot2::labs(title = title, x = NULL, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = sample_colors) +
    ggplot2::scale_color_manual(values = sample_colors)

  if (!is.null(hline_min)) {
    p <- p + ggplot2::geom_hline(yintercept = hline_min, color = "red", linetype = "dashed", linewidth = 0.7)
  }

  if (metric == "nFeature_RNA" && !is.null(per_sample_max) && nrow(per_sample_max) > 0) {
    max_segments <- per_sample_max_box_segments(
      per_sample_max,
      levels(data$orig.ident),
      half_width = 0.275
    )
    p <- p + ggplot2::geom_segment(
      data = max_segments,
      ggplot2::aes(x = xmin, xend = xmax, y = nFeature_max, yend = nFeature_max),
      inherit.aes = FALSE,
      color = "darkred",
      linetype = "dotted",
      linewidth = 0.7
    )
  }
  p
}

save_qc_zoom_box_plots_median5sd <- function(
  md,
  out_dir,
  gene_label,
  th,
  y_min = 200,
  y_max = 3000
) {
  plot_dir <- file.path(out_dir, "qc_plots", "zoom")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  md$orig.ident <- factor(md$orig.ident, levels = sort(unique(as.character(md$orig.ident))))

  p_nFeature <- plot_qc_box_zoom_vertical_per_sample_max(
    md,
    "nFeature_RNA",
    paste0(gene_label, " — nFeature (zoom ", y_min, "–", y_max, ")"),
    "nFeature_RNA",
    per_sample_max = th$per_sample_max,
    y_min = y_min,
    y_max = y_max,
    hline_min = th$nFeature_min
  )
  p_nCount <- plot_qc_box_zoom_vertical(
    md,
    "nCount_RNA",
    paste0(gene_label, " — nCount (zoom ", y_min, "–", y_max, ")"),
    "nCount_RNA",
    y_min = y_min,
    y_max = y_max
  )

  suffix <- paste0("_vertical_", y_min, "_", y_max)
  save_gg(p_nFeature, file.path(plot_dir, paste0("box_nFeature", suffix)), width = 11, height = 7)
  save_gg(p_nCount, file.path(plot_dir, paste0("box_nCount", suffix)), width = 11, height = 7)
  invisible(plot_dir)
}

save_filter_qc_plots_median5sd <- function(
  md,
  th,
  out_dir,
  stage,
  seurat_filt = NULL,
  n_before = NULL,
  md_before = NULL,
  density_xlim_nFeature = c(0, 6000),
  density_xlim_nCount = c(0, 15000),
  qc_zoom_y_min = 200,
  qc_zoom_y_max = 3000
) {
  label <- threshold_folder_label_median5sd(th)
  plot_dir <- file.path(out_dir, "qc_filter_plots", label, stage)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  md$orig.ident <- factor(md$orig.ident, levels = sort(unique(md$orig.ident)))
  per_sample_max <- th$per_sample_max
  per_sample_max$orig.ident <- factor(per_sample_max$sample, levels = levels(md$orig.ident))

  save_gg(
    plot_qc_density_per_sample_max(
      md, "nFeature_RNA", paste(stage, "nFeature density (zoom)"), "nFeature_RNA",
      per_sample_max = per_sample_max,
      vline_min = th$nFeature_min,
      xlim = density_xlim_nFeature
    ),
    file.path(plot_dir, "density_nFeature"), width = 14, height = 10
  )
  save_gg(
    plot_qc_density(
      md, "nCount_RNA", paste(stage, "nCount density (zoom)"), "nCount_RNA",
      NULL, NULL, xlim = density_xlim_nCount
    ),
    file.path(plot_dir, "density_nCount"), width = 14, height = 10
  )

  p_nFeature <- plot_qc_box_zoom_vertical_per_sample_max(
    md,
    "nFeature_RNA",
    paste(stage, "nFeature (zoom ", qc_zoom_y_min, "–", qc_zoom_y_max, ")", sep = ""),
    "nFeature_RNA",
    per_sample_max = per_sample_max,
    y_min = qc_zoom_y_min,
    y_max = qc_zoom_y_max,
    hline_min = th$nFeature_min
  )
  p_nCount <- plot_qc_box_zoom_vertical(
    md,
    "nCount_RNA",
    paste(stage, "nCount (zoom ", qc_zoom_y_min, "–", qc_zoom_y_max, ")", sep = ""),
    "nCount_RNA",
    y_min = qc_zoom_y_min,
    y_max = qc_zoom_y_max
  )
  suffix <- paste0("_vertical_", qc_zoom_y_min, "_", qc_zoom_y_max)
  save_gg(p_nFeature, file.path(plot_dir, paste0("box_nFeature", suffix)), width = 11, height = 7)
  save_gg(p_nCount, file.path(plot_dir, paste0("box_nCount", suffix)), width = 11, height = 7)

  save_gg(
    plot_qc_scatter_per_sample_max(
      md, paste(stage, "nFeature vs nCount"),
      per_sample_max = per_sample_max,
      nFeature_min = th$nFeature_min
    ),
    file.path(plot_dir, "scatter_nFeature_nCount"), width = 14, height = 10
  )

  filter_summary <- NULL
  if (identical(stage, "after_filter") && !is.null(seurat_filt) && !is.null(n_before)) {
    filter_summary <- save_filter_summary_median5sd(
      seurat_filt, th, n_before, out_dir, md_before = md_before
    )
  }

  list(plot_dir = plot_dir, filter_summary = filter_summary)
}

save_filter_summary_median5sd <- function(seurat_filt, th, n_before, out_dir, md_before = NULL) {
  label <- threshold_folder_label_median5sd(th)
  summary_dir <- file.path(out_dir, "filter_summary", label)
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

  cells_per_sample <- as.data.frame(table(seurat_filt$orig.ident))
  colnames(cells_per_sample) <- c("sample", "n_cells")
  write.csv(
    cells_per_sample,
    file.path(summary_dir, "cells_per_sample_after_filter.csv"),
    row.names = FALSE
  )

  write.csv(
    th$per_sample_max,
    file.path(summary_dir, "per_sample_nfeature_max_median5sd.csv"),
    row.names = FALSE
  )

  overall <- data.frame(
    threshold_label = label,
    mode = th$mode,
    nFeature_min = th$nFeature_min,
    nFeature_max_rule = paste0("median + ", th$sd_multiplier, " * SD (per sample)"),
    sd_multiplier = th$sd_multiplier,
    mt_max = th$mt_max,
    n_before = n_before,
    n_after = ncol(seurat_filt),
    removed = n_before - ncol(seurat_filt),
    retention_pct = round(100 * ncol(seurat_filt) / n_before, 1),
    stringsAsFactors = FALSE
  )
  write.csv(
    overall,
    file.path(summary_dir, "filter_summary_overall.csv"),
    row.names = FALSE
  )

  list(
    summary_dir = summary_dir,
    cells_per_sample = cells_per_sample,
    overall = overall
  )
}

run_median5sd_threshold_trial <- function(
  seurat_all,
  th,
  out_dir,
  md_before = NULL,
  density_xlim_nFeature = c(0, 6000),
  density_xlim_nCount = c(0, 15000),
  qc_zoom_y_min = 200,
  qc_zoom_y_max = 3000
) {
  label <- threshold_folder_label_median5sd(th)
  message("  Trial: ", label)

  if (is.null(md_before)) {
    md_before <- seurat_all@meta.data
  }

  plot_args <- list(
    density_xlim_nFeature = density_xlim_nFeature,
    density_xlim_nCount = density_xlim_nCount,
    qc_zoom_y_min = qc_zoom_y_min,
    qc_zoom_y_max = qc_zoom_y_max
  )

  do.call(
    save_filter_qc_plots_median5sd,
    c(list(md = md_before, th = th, out_dir = out_dir, stage = "before_filter"), plot_args)
  )
  n_before <- ncol(seurat_all)
  seurat_filt <- filter_seurat_by_median5sd(seurat_all, th)

  after <- do.call(
    save_filter_qc_plots_median5sd,
    c(
      list(
        md = seurat_filt@meta.data,
        th = th,
        out_dir = out_dir,
        stage = "after_filter",
        seurat_filt = seurat_filt,
        n_before = n_before,
        md_before = md_before
      ),
      plot_args
    )
  )

  list(
    label = label,
    n_before = n_before,
    n_after = ncol(seurat_filt),
    filter_summary = after$filter_summary
  )
}

detect_seurat_clustree_prefix <- function(obj) {
  cols <- grep("_res\\.", colnames(obj@meta.data), value = TRUE)
  if (length(cols) == 0) {
    stop("No resolution cluster columns found for clustree (expected names like RNA_snn_res.0.1).")
  }
  sub("[0-9]+(\\.[0-9]+)?$", "", cols[[1]])
}

save_clustree_plot <- function(obj, plot_dir, prefix = NULL, ...) {
  if (!requireNamespace("clustree", quietly = TRUE)) {
    stop("Install clustree before running Step 6b: install.packages('clustree')")
  }
  if (is.null(prefix)) {
    prefix <- detect_seurat_clustree_prefix(obj)
  }
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  p <- clustree::clustree(obj, prefix = prefix, ...)
  ggplot2::ggsave(
    file.path(plot_dir, "clustree.png"),
    p,
    width = 12,
    height = 10,
    dpi = 150,
    limitsize = FALSE
  )
  p
}

clustering_dir_median5sd <- function(out_dir, th) {
  file.path(out_dir, "clustering", threshold_folder_label_median5sd(th))
}

clustering_rds_paths <- function(cl_dir, object_prefix) {
  list(
    filtered = file.path(cl_dir, paste0(object_prefix, "_filtered.rds")),
    sct_per_sample_dir = file.path(cl_dir, "sct_per_sample"),
    harmony_prepped = file.path(cl_dir, paste0(object_prefix, "_harmony_prepped.rds")),
    clustree = file.path(cl_dir, paste0(object_prefix, "_clustree.rds"))
  )
}

sct_per_sample_rds_path <- function(sct_dir, object_prefix, sample_id) {
  file.path(sct_dir, paste0(object_prefix, "_sct_", sample_id, ".rds"))
}

sct_norm_method_label <- function() {
  if (exists("use_glmGamPoi") && use_glmGamPoi) {
    "SCTransform (glmGamPoi)"
  } else {
    "SCTransform (native)"
  }
}

slim_sct_seurat <- function(obj) {
  DefaultAssay(obj) <- "SCT"
  if ("RNA" %in% Assays(obj)) {
    obj[["RNA"]] <- NULL
  }
  gc()
  obj
}

sct_params_per_sample <- function() {
  c(
    sct_params(),
    list(
      # Store residuals for variable genes only — much lower RAM than all ~24k genes
      return.only.var.genes = TRUE
    )
  )
}

sct_one_sample <- function(obj, verbose = TRUE) {
  if (packageVersion("Seurat") >= "5.0.0") {
    obj <- JoinLayers(obj)
  }
  params <- c(list(object = obj), sct_params_per_sample())
  if (verbose) {
    params$verbose <- TRUE
  }
  message("    SCTransform running (return.only.var.genes = TRUE) ...")
  obj <- do.call(SCTransform, params)
  rm(params)
  gc()
  slim_sct_seurat(obj)
}

list_expected_sct_samples <- function(sample_ids, samples_to_run = NULL) {
  sample_ids <- sort(unique(as.character(sample_ids)))
  if (!is.null(samples_to_run)) {
    missing <- setdiff(samples_to_run, sample_ids)
    if (length(missing) > 0) {
      stop("sct_samples_to_run not found in filtered object: ", paste(missing, collapse = ", "))
    }
    sample_ids <- intersect(sample_ids, samples_to_run)
  }
  sample_ids
}

peek_filtered_seurat_meta <- function(rds_path) {
  obj <- readRDS(rds_path)
  on.exit({
    rm(obj)
    gc()
  }, add = TRUE)
  list(
    n_cells = ncol(obj),
    sample_ids = sort(unique(as.character(obj$orig.ident)))
  )
}

ensure_sct_filter_rds <- function(
  raw_rds,
  th,
  object_prefix,
  paths,
  mt_pattern = "^mt-",
  filtered_rds = NULL,
  force = FALSE
) {
  if (!is.null(filtered_rds) && nzchar(filtered_rds) && file.exists(filtered_rds)) {
    message("Using filtered cells from: ", filtered_rds)
    return(filtered_rds)
  }
  if (file.exists(paths$filtered) && !force) {
    message("Using existing filtered RDS: ", paths$filtered)
    return(paths$filtered)
  }
  message("Building filtered RDS from raw merge (one-time)...")
  seurat_filt <- prep_filtered_seurat(raw_rds, th, object_prefix, mt_pattern)
  saveRDS(seurat_filt, paths$filtered)
  rm(seurat_filt)
  gc()
  paths$filtered
}

ensure_filtered_sample_splits <- function(filter_rds, cl_dir, sample_ids, force = FALSE) {
  split_dir <- file.path(cl_dir, "filtered_per_sample")
  dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)

  split_paths <- file.path(split_dir, paste0(sample_ids, ".rds"))
  if (all(file.exists(split_paths)) && !force) {
    message("Using cached per-sample filtered RDS in ", split_dir)
    return(split_dir)
  }

  message("Splitting filtered object into per-sample RDS (one-time)...")
  obj <- readRDS(filter_rds)
  splits <- SplitObject(obj, split.by = "orig.ident")
  rm(obj)
  gc()

  for (sid in names(splits)) {
    saveRDS(splits[[sid]], file.path(split_dir, paste0(sid, ".rds")))
  }
  rm(splits)
  gc()
  split_dir
}

load_filtered_sample <- function(split_dir, sample_id) {
  sample_path <- file.path(split_dir, paste0(sample_id, ".rds"))
  if (!file.exists(sample_path)) {
    stop("Missing per-sample filtered RDS: ", sample_path, ". Re-run Step 6a with step6_force_rerun = TRUE.")
  }
  obj <- readRDS(sample_path)
  list(object = obj, n_cells = ncol(obj))
}

build_cell_filter_from_raw <- function(
  raw_rds,
  nFeature_min,
  sd_multiplier,
  mt_max,
  mt_pattern = "^mt-"
) {
  seurat_all <- readRDS(raw_rds)
  if (!"percent.mt" %in% colnames(seurat_all@meta.data)) {
    seurat_all[["percent.mt"]] <- PercentageFeatureSet(seurat_all, pattern = mt_pattern)
  }
  md <- seurat_all@meta.data
  per_sample_thresholds <- compute_per_sample_nfeature_max(
    md,
    nFeature_min = nFeature_min,
    sd_multiplier = sd_multiplier
  )
  make_median5sd_cell_filter(
    mt_max = mt_max,
    per_sample_thresholds = per_sample_thresholds,
    nFeature_min = nFeature_min,
    sd_multiplier = sd_multiplier
  )
}

prep_filtered_seurat <- function(raw_rds, th, object_prefix, mt_pattern = "^mt-") {
  seurat_all <- readRDS(raw_rds)
  if (!"percent.mt" %in% colnames(seurat_all@meta.data)) {
    seurat_all[["percent.mt"]] <- PercentageFeatureSet(seurat_all, pattern = mt_pattern)
  }
  seurat_filt <- filter_seurat_by_median5sd(seurat_all, th)
  seurat_filt <- add_genotype_metadata(seurat_filt, object_prefix)
  rm(seurat_all)
  gc()
  seurat_filt
}

run_step6a_sct_per_sample <- function(
  raw_rds,
  th,
  out_dir,
  object_prefix,
  mt_pattern = "^mt-",
  samples_to_run = NULL,
  filtered_rds = NULL,
  force = FALSE
) {
  label <- threshold_folder_label_median5sd(th)
  cl_dir <- clustering_dir_median5sd(out_dir, th)
  dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- clustering_rds_paths(cl_dir, object_prefix)
  dir.create(paths$sct_per_sample_dir, recursive = TRUE, showWarnings = FALSE)

  filter_rds <- ensure_sct_filter_rds(
    raw_rds = raw_rds,
    th = th,
    object_prefix = object_prefix,
    paths = paths,
    mt_pattern = mt_pattern,
    filtered_rds = filtered_rds,
    force = force
  )

  meta <- peek_filtered_seurat_meta(filter_rds)
  message("\n=== Step 6a per-sample SCT | ", object_prefix, " | ", label, " ===")
  message("Cells after filter: ", meta$n_cells)
  message("Memory mode: split filtered RDS per sample; one sample + SCT in RAM at a time.")

  sample_ids <- list_expected_sct_samples(meta$sample_ids, samples_to_run)
  rm(meta)
  gc()

  split_dir <- ensure_filtered_sample_splits(
    filter_rds = filter_rds,
    cl_dir = cl_dir,
    sample_ids = sample_ids,
    force = force
  )

  norm_method <- sct_norm_method_label()
  summary_rows <- list()

  for (sample_id in sample_ids) {
    out_path <- sct_per_sample_rds_path(paths$sct_per_sample_dir, object_prefix, sample_id)
    if (file.exists(out_path) && !force) {
      message("  SKIP ", sample_id, " (", basename(out_path), " exists)")
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        sample = sample_id,
        n_cells = NA_integer_,
        sct_rds = out_path,
        status = "skipped_existing",
        stringsAsFactors = FALSE
      )
      next
    }

    sample_load <- load_filtered_sample(split_dir, sample_id)
    message(
      "  SCTransform: ", sample_id, " (", sample_load$n_cells,
      " cells; loaded from per-sample filtered RDS)"
    )
    obj <- sct_one_sample(sample_load$object)
    rm(sample_load)
    gc()
    saveRDS(obj, out_path)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      sample = sample_id,
      n_cells = ncol(obj),
      sct_rds = out_path,
      status = "completed",
      stringsAsFactors = FALSE
    )
    rm(obj)
    gc()
  }

  per_sample_summary <- do.call(rbind, summary_rows)
  write.csv(
    per_sample_summary,
    file.path(cl_dir, "sct_per_sample_summary.csv"),
    row.names = FALSE
  )

  list(
    cl_dir = cl_dir,
    paths = paths,
    label = label,
    norm_method = norm_method,
    per_sample_summary = per_sample_summary,
    skipped = FALSE
  )
}

run_step6a_merge_harmony <- function(
  th,
  out_dir,
  object_prefix,
  harmony_dims = 1:30,
  graph_name = "RNA_snn",
  integration_features = 3000,
  force = FALSE
) {
  label <- threshold_folder_label_median5sd(th)
  cl_dir <- clustering_dir_median5sd(out_dir, th)
  paths <- clustering_rds_paths(cl_dir, object_prefix)

  if (file.exists(paths$harmony_prepped) && !force) {
    message(
      "\n=== Step 6a-merge | ", object_prefix, " | ", label,
      " — SKIP (", basename(paths$harmony_prepped), " exists; set step6_force_rerun = TRUE to redo) ==="
    )
    return(list(
      cl_dir = cl_dir,
      paths = paths,
      label = label,
      norm_method = sct_norm_method_label(),
      skipped = TRUE
    ))
  }

  sample_paths <- sort(Sys.glob(
    file.path(paths$sct_per_sample_dir, paste0(object_prefix, "_sct_*.rds"))
  ))
  if (length(sample_paths) == 0) {
    stop("Run Step 6a first (no per-sample SCT RDS in ", paths$sct_per_sample_dir, ").")
  }

  message("\n=== Step 6a-merge Harmony/UMAP | ", object_prefix, " | ", label, " ===")
  message("  Loading ", length(sample_paths), " per-sample SCT objects...")
  obj_list <- lapply(sample_paths, readRDS)
  features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = integration_features)

  merged <- obj_list[[1]]
  if (length(obj_list) > 1) {
    for (i in 2:length(obj_list)) {
      merged <- merge(merged, y = obj_list[[i]])
      rm(obj_list[[i]])
      gc()
    }
  }
  rm(obj_list)
  gc()

  VariableFeatures(merged) <- features
  DefaultAssay(merged) <- "SCT"

  merged <- RunPCA(merged, assay = "SCT", verbose = FALSE)
  merged <- RunHarmony(merged, group.by.vars = "orig.ident", verbose = FALSE)
  merged <- RunUMAP(merged, reduction = "harmony", dims = harmony_dims)
  merged <- FindNeighbors(
    merged,
    reduction = "harmony",
    dims = harmony_dims,
    graph.name = graph_name
  )

  norm_method <- sct_norm_method_label()
  saveRDS(merged, paths$harmony_prepped)
  write.csv(
    data.frame(
      step = "6a_merge_harmony",
      threshold_label = label,
      normalization = norm_method,
      n_cells = ncol(merged),
      n_samples = length(sample_paths),
      integration_features = length(features),
      filtered_rds = paths$filtered,
      harmony_prepped_rds = paths$harmony_prepped,
      stringsAsFactors = FALSE
    ),
    file.path(cl_dir, "sctransform_summary.csv"),
    row.names = FALSE
  )

  list(
    cl_dir = cl_dir,
    paths = paths,
    label = label,
    norm_method = norm_method,
    n_cells = ncol(merged),
    n_samples = length(sample_paths),
    skipped = FALSE
  )
}

run_step6a_sctransform <- function(
  raw_rds,
  th,
  out_dir,
  object_prefix,
  mt_pattern = "^mt-",
  harmony_dims = 1:30,
  graph_name = "RNA_snn",
  force = FALSE
) {
  label <- threshold_folder_label_median5sd(th)
  cl_dir <- clustering_dir_median5sd(out_dir, th)
  dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- clustering_rds_paths(cl_dir, object_prefix)

  if (file.exists(paths$harmony_prepped) && !force) {
    message(
      "\n=== Step 6a | ", object_prefix, " | ", label,
      " — SKIP (", basename(paths$harmony_prepped), " exists; set step6_force_rerun = TRUE to redo) ==="
    )
    sct_summary_path <- file.path(cl_dir, "sctransform_summary.csv")
    norm_method <- if (file.exists(sct_summary_path)) {
      read.csv(sct_summary_path, stringsAsFactors = FALSE)$normalization[1]
    } else {
      sct_norm_method_label()
    }
    return(list(
      cl_dir = cl_dir,
      paths = paths,
      label = label,
      norm_method = norm_method,
      skipped = TRUE
    ))
  }

  message("\n=== Step 6a SCTransform + Harmony | ", object_prefix, " | ", label, " ===")
  seurat_filt <- prep_filtered_seurat(raw_rds, th, object_prefix, mt_pattern)
  message("Cells after filter: ", ncol(seurat_filt))
  saveRDS(seurat_filt, paths$filtered)

  if (packageVersion("Seurat") >= "5.0.0") {
    seurat_filt <- JoinLayers(seurat_filt)
  }

  norm_method <- sct_norm_method_label()
  seurat_filt <- tryCatch(
    {
      message(
        "  SCTransform (",
        if (exists("use_glmGamPoi") && use_glmGamPoi) "glmGamPoi" else "native",
        ") on merged filtered object..."
      )
      do.call(SCTransform, c(list(object = seurat_filt), sct_params()))
    },
    error = function(e) {
      norm_method <<- "LogNormalize"
      message("  SCTransform failed; using LogNormalize. Reason: ", conditionMessage(e))
      seurat_filt <- NormalizeData(seurat_filt, verbose = FALSE)
      seurat_filt <- FindVariableFeatures(
        seurat_filt,
        selection.method = "vst",
        nfeatures = 3000,
        verbose = FALSE
      )
      seurat_filt <- ScaleData(seurat_filt, verbose = FALSE)
      seurat_filt
    }
  )
  gc()

  seurat_filt <- RunPCA(seurat_filt, verbose = FALSE)
  seurat_filt <- RunHarmony(seurat_filt, group.by.vars = "orig.ident", verbose = FALSE)
  seurat_filt <- RunUMAP(seurat_filt, reduction = "harmony", dims = harmony_dims)
  seurat_filt <- FindNeighbors(
    seurat_filt,
    reduction = "harmony",
    dims = harmony_dims,
    graph.name = graph_name
  )

  saveRDS(seurat_filt, paths$harmony_prepped)
  write.csv(
    data.frame(
      step = "6a_sctransform",
      threshold_label = label,
      normalization = norm_method,
      n_cells = ncol(seurat_filt),
      filtered_rds = paths$filtered,
      harmony_prepped_rds = paths$harmony_prepped,
      stringsAsFactors = FALSE
    ),
    file.path(cl_dir, "sctransform_summary.csv"),
    row.names = FALSE
  )

  list(
    cl_dir = cl_dir,
    paths = paths,
    label = label,
    norm_method = norm_method,
    n_cells = ncol(seurat_filt),
    skipped = FALSE
  )
}

run_step6b_clustree <- function(
  th,
  out_dir,
  object_prefix,
  clustree_resolutions = seq(0, 1, 0.1),
  graph_name = "RNA_snn",
  force = FALSE
) {
  label <- threshold_folder_label_median5sd(th)
  cl_dir <- clustering_dir_median5sd(out_dir, th)
  paths <- clustering_rds_paths(cl_dir, object_prefix)
  clustree_png <- file.path(cl_dir, "clustree.png")

  if (!file.exists(paths$harmony_prepped)) {
    stop("Run Step 6a first (missing ", paths$harmony_prepped, ").")
  }

  if (file.exists(clustree_png) && file.exists(paths$clustree) && !force) {
    message(
      "\n=== Step 6b | ", object_prefix, " | ", label,
      " — SKIP (clustree outputs exist; set step6_force_rerun = TRUE to redo) ==="
    )
    return(list(
      cl_dir = cl_dir,
      paths = paths,
      label = label,
      clustree_prefix = detect_seurat_clustree_prefix(readRDS(paths$clustree)),
      skipped = TRUE
    ))
  }

  message("\n=== Step 6b clustree | ", object_prefix, " | ", label, " ===")
  seurat_filt <- readRDS(paths$harmony_prepped)

  message("  FindClusters for clustree (", length(clustree_resolutions), " resolutions)...")
  for (res in clustree_resolutions) {
    seurat_filt <- FindClusters(
      seurat_filt,
      resolution = res,
      graph.name = graph_name,
      verbose = FALSE
    )
  }

  clustree_prefix <- detect_seurat_clustree_prefix(seurat_filt)
  message("  clustree prefix: ", clustree_prefix)
  save_clustree_plot(seurat_filt, cl_dir, prefix = clustree_prefix)
  saveRDS(seurat_filt, paths$clustree)

  res_cols <- grep(
    paste0("^", gsub("\\.", "\\\\.", clustree_prefix)),
    colnames(seurat_filt@meta.data),
    value = TRUE
  )
  write.csv(
    data.frame(resolution_column = res_cols, stringsAsFactors = FALSE),
    file.path(cl_dir, "clustree_resolution_columns.csv"),
    row.names = FALSE
  )

  list(
    cl_dir = cl_dir,
    paths = paths,
    label = label,
    clustree_prefix = clustree_prefix,
    skipped = FALSE
  )
}

run_step6c_clustering <- function(
  th,
  out_dir,
  object_prefix,
  resolutions,
  graph_name = "RNA_snn",
  force = FALSE
) {
  if (length(resolutions) == 0) {
    stop("Set clustering_resolutions in Pipeline options before running Step 6c.")
  }

  label <- threshold_folder_label_median5sd(th)
  cl_dir <- clustering_dir_median5sd(out_dir, th)
  paths <- clustering_rds_paths(cl_dir, object_prefix)
  input_rds <- if (file.exists(paths$clustree)) paths$clustree else paths$harmony_prepped

  if (!file.exists(input_rds)) {
    stop("Run Step 6a first (missing ", paths$harmony_prepped, ").")
  }

  summary_path <- file.path(cl_dir, "clustering_summary.csv")
  if (file.exists(summary_path) && !force) {
    existing <- read.csv(summary_path, stringsAsFactors = FALSE)
    if (all(resolutions %in% existing$resolution)) {
      message(
        "\n=== Step 6c | ", object_prefix, " | ", label,
        " — SKIP (clustering_summary.csv exists; set step6_force_rerun = TRUE to redo) ==="
      )
      return(existing)
    }
  }

  message("\n=== Step 6c clustering + UMAP plots | ", object_prefix, " | ", label, " ===")
  seurat_filt <- readRDS(input_rds)
  sct_summary_path <- file.path(cl_dir, "sctransform_summary.csv")
  norm_method <- if (file.exists(sct_summary_path)) {
    read.csv(sct_summary_path, stringsAsFactors = FALSE)$normalization[1]
  } else {
    sct_norm_method_label()
  }

  summary_rows <- list()
  for (res in resolutions) {
    res_tag <- gsub("\\.", "", as.character(res))
    rds_path <- file.path(cl_dir, paste0(object_prefix, "_harmony_res", res, ".rds"))

    if (file.exists(rds_path) && !force) {
      message("  res ", res, " — loading existing ", basename(rds_path))
      obj <- readRDS(rds_path)
    } else {
      obj <- FindClusters(
        seurat_filt,
        resolution = res,
        graph.name = graph_name,
        verbose = FALSE
      )
      saveRDS(obj, rds_path)
      seurat_filt <- obj
    }

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
      filtered_rds = paths$filtered,
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
  write.csv(summary_df, summary_path, row.names = FALSE)
  write.csv(summary_df, file.path(cl_dir, "clustering_run_summary.csv"), row.names = FALSE)
  invisible(summary_df)
}

run_clustering_pipeline_sct_clustree <- function(
  raw_rds,
  th,
  out_dir,
  object_prefix,
  mt_pattern = "^mt-",
  harmony_dims = 1:30,
  resolutions = c(0.5, 2.5),
  clustree_resolutions = seq(0, 1, 0.1),
  graph_name = "RNA_snn",
  force = FALSE
) {
  run_step6a_sctransform(
    raw_rds = raw_rds,
    th = th,
    out_dir = out_dir,
    object_prefix = object_prefix,
    mt_pattern = mt_pattern,
    harmony_dims = harmony_dims,
    graph_name = graph_name,
    force = force
  )
  run_step6b_clustree(
    th = th,
    out_dir = out_dir,
    object_prefix = object_prefix,
    clustree_resolutions = clustree_resolutions,
    graph_name = graph_name,
    force = force
  )
  run_step6c_clustering(
    th = th,
    out_dir = out_dir,
    object_prefix = object_prefix,
    resolutions = resolutions,
    graph_name = graph_name,
    force = force
  )
}

snRNAseq_median5sd_clustree_init <- function(repo_root = "C:/Yale") {
  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(harmony)
    library(future)
  })
  options(future.globals.maxSize = 8 * 1024^3)
  plan("sequential")

  qc_helpers <- file.path(repo_root, "scripts", "snRNAseq_qc_filter_helpers.R")
  median_helpers <- file.path(repo_root, "scripts", "snRNAseq_median5sd_clustree_helpers.R")
  cluster_helpers <- file.path(repo_root, "scripts", "snRNAseq_clustering_helpers.R")

  if (!exists("save_gg", mode = "function")) {
    source(qc_helpers)
  }
  if (!exists("compute_per_sample_nfeature_max", mode = "function")) {
    source(median_helpers)
  }
  if (!exists("sct_params", mode = "function")) {
    source(cluster_helpers)
  }

  assign("use_glmGamPoi", requireNamespace("glmGamPoi", quietly = TRUE), envir = .GlobalEnv)
  if (use_glmGamPoi) {
    message("glmGamPoi detected — SCTransform will use vst.flavor = 'v2'.")
  } else {
    message("glmGamPoi not found — SCTransform will use vst.flavor = 'v1' (install glmGamPoi for speed).")
  }
  invisible(TRUE)
}
