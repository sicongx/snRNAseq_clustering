# Helpers for scn1lab / dyrk1a full pipeline (QC plots + single-user-defined filter)

save_gg <- function(p, path_base, width = 12, height = 7) {
  ggplot2::ggsave(paste0(path_base, ".png"), p, width = width, height = height, dpi = 150, limitsize = FALSE)
}

plot_qc_density <- function(
  data,
  metric,
  title,
  xlab,
  vline_min = NULL,
  vline_max = NULL,
  xlim = NULL
) {
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[metric]])) +
    ggplot2::geom_density(fill = "steelblue", alpha = 0.35, color = "steelblue", linewidth = 0.4) +
    ggplot2::facet_wrap(~orig.ident, scales = "free_y", ncol = 3) +
    ggplot2::labs(title = title, x = xlab, y = "Density") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
  if (!is.null(vline_min)) {
    p <- p + ggplot2::geom_vline(xintercept = vline_min, color = "red", linetype = "dashed", linewidth = 0.7)
  }
  if (!is.null(vline_max)) {
    p <- p + ggplot2::geom_vline(xintercept = vline_max, color = "darkred", linetype = "dotted", linewidth = 0.7)
  }
  if (!is.null(xlim) && length(xlim) == 2L) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }
  p
}

plot_qc_box <- function(data, metric, title, ylab, hline_min = NULL, hline_max = NULL) {
  p <- ggplot2::ggplot(data, ggplot2::aes(x = orig.ident, y = .data[[metric]], fill = orig.ident)) +
    ggplot2::geom_boxplot(outlier.size = 0.3) +
    ggplot2::labs(title = title, x = NULL, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip()
  if (!is.null(hline_min)) {
    p <- p + ggplot2::geom_hline(yintercept = hline_min, color = "red", linetype = "dashed", linewidth = 0.7)
  }
  if (!is.null(hline_max)) {
    p <- p + ggplot2::geom_hline(yintercept = hline_max, color = "darkred", linetype = "dotted", linewidth = 0.7)
  }
  p
}

plot_qc_scatter <- function(data, title, nFeature_min = NULL, nFeature_max = NULL) {
  p <- ggplot2::ggplot(data, ggplot2::aes(x = nFeature_RNA, y = nCount_RNA)) +
    ggplot2::geom_point(color = "steelblue", alpha = 0.08, size = 0.25) +
    ggplot2::facet_wrap(~orig.ident, scales = "free", ncol = 3) +
    ggplot2::labs(title = title, x = "nFeature_RNA", y = "nCount_RNA") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
  if (!is.null(nFeature_min)) {
    p <- p + ggplot2::geom_vline(xintercept = nFeature_min, color = "red", linetype = "dashed", linewidth = 0.7)
  }
  if (!is.null(nFeature_max)) {
    p <- p + ggplot2::geom_vline(xintercept = nFeature_max, color = "darkred", linetype = "dotted", linewidth = 0.7)
  }
  p
}

subsample_meta_for_points <- function(md, points_per_sample = 1500, seed = 1) {
  set.seed(seed)
  do.call(
    rbind,
    lapply(split(md, md$orig.ident, drop = TRUE), function(d) {
      if (nrow(d) > points_per_sample) {
        d[sample.int(nrow(d), points_per_sample), , drop = FALSE]
      } else {
        d
      }
    })
  )
}

plot_qc_box_zoom_vertical <- function(
  data,
  metric,
  title,
  ylab,
  y_min = 200,
  y_max = 3000,
  points_per_sample = 1500,
  hline_min = NULL,
  hline_max = NULL
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
  if (!is.null(hline_max)) {
    p <- p + ggplot2::geom_hline(yintercept = hline_max, color = "darkred", linetype = "dotted", linewidth = 0.7)
  }
  p
}

save_qc_zoom_box_plots <- function(md, out_dir, gene_label, y_min = 200, y_max = 3000, th = NULL) {
  plot_dir <- file.path(out_dir, "qc_plots", "zoom")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  md$orig.ident <- factor(md$orig.ident, levels = sort(unique(as.character(md$orig.ident))))

  nF_min <- if (!is.null(th)) th$nFeature_min else NULL
  nF_max <- if (!is.null(th)) th$nFeature_max else NULL

  p_nFeature <- plot_qc_box_zoom_vertical(
    md,
    "nFeature_RNA",
    paste0(gene_label, " — nFeature (zoom ", y_min, "–", y_max, ")"),
    "nFeature_RNA",
    y_min = y_min,
    y_max = y_max,
    hline_min = nF_min,
    hline_max = nF_max
  )
  p_nCount <- plot_qc_box_zoom_vertical(
    md,
    "nCount_RNA",
    paste0(gene_label, " — nCount (zoom ", y_min, "–", y_max, ")"),
    "nCount_RNA",
    y_min = y_min,
    y_max = y_max
  )

  save_gg(
    p_nFeature,
    file.path(plot_dir, paste0("box_nFeature_vertical_", y_min, "_", y_max)),
    width = 11,
    height = 7
  )
  save_gg(
    p_nCount,
    file.path(plot_dir, paste0("box_nCount_vertical_", y_min, "_", y_max)),
    width = 11,
    height = 7
  )

  list(
    plot_dir = plot_dir,
    nFeature = p_nFeature,
    nCount = p_nCount
  )
}

save_qc_exploratory_plots <- function(
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

  feature_lines <- list(min = th$nFeature_min, max = th$nFeature_max)

  save_gg(
    plot_qc_density(
      md, "nFeature_RNA", "nFeature — density (zoom)", "nFeature_RNA",
      feature_lines$min, feature_lines$max, xlim = density_xlim_nFeature
    ),
    file.path(plot_dir, "density_nFeature"), width = 14, height = 10
  )
  save_gg(
    plot_qc_density(
      md, "nCount_RNA", "nCount — density (zoom)", "nCount_RNA",
      NULL, NULL, xlim = density_xlim_nCount
    ),
    file.path(plot_dir, "density_nCount"), width = 14, height = 10
  )
  save_gg(
    plot_qc_box(md, "nFeature_RNA", "nFeature — boxplot", "nFeature_RNA", feature_lines$min, feature_lines$max),
    file.path(plot_dir, "box_nFeature"), width = 10, height = 6
  )
  save_gg(
    plot_qc_box(md, "nCount_RNA", "nCount — boxplot", "nCount_RNA", NULL, NULL),
    file.path(plot_dir, "box_nCount"), width = 10, height = 6
  )
  save_gg(
    plot_qc_scatter(md, "nFeature vs nCount — scatter", feature_lines$min, feature_lines$max),
    file.path(plot_dir, "scatter_nFeature_nCount"), width = 14, height = 10
  )
  invisible(plot_dir)
}

save_doublet_plots <- function(seurat_objs, out_dir, doublet_summary) {
  plot_dir <- file.path(out_dir, "doublet_plots", "all_samples")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  all_md <- do.call(rbind, lapply(names(seurat_objs), function(id) {
    md <- seurat_objs[[id]]@meta.data
    md$sample <- id
    md
  }))

  p_rate <- ggplot2::ggplot(doublet_summary, ggplot2::aes(x = reorder(sample, -doublet_rate), y = doublet_rate)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::geom_text(ggplot2::aes(label = paste0(doublet_rate, "%")), vjust = -0.3, size = 2.8) +
    ggplot2::labs(title = "Doublet rate per sample", x = NULL, y = "Doublet rate (%)") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p_class <- ggplot2::ggplot(all_md, ggplot2::aes(x = sample, fill = scDblFinder.class)) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
    ggplot2::labs(title = "Singlet vs doublet fraction", x = NULL, y = "Fraction") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  save_gg(p_rate, file.path(plot_dir, "doublet_rate_bar"), width = 10, height = 5)
  save_gg(p_class, file.path(plot_dir, "class_fraction"), width = 10, height = 5)
  invisible(plot_dir)
}

make_cell_filter <- function(nFeature_min, nFeature_max, mt_max) {
  bad <- character(0)
  if (is.null(nFeature_min) || length(nFeature_min) != 1L || is.na(nFeature_min) || !is.numeric(nFeature_min)) {
    bad <- c(bad, "nFeature_min")
  }
  if (is.null(mt_max) || length(mt_max) != 1L || is.na(mt_max) || !is.numeric(mt_max)) {
    bad <- c(bad, "mt_max")
  }
  no_upper <- is.null(nFeature_max) || (length(nFeature_max) == 1L && is.na(nFeature_max))
  if (!no_upper) {
    if (length(nFeature_max) != 1L || is.na(nFeature_max) || !is.numeric(nFeature_max)) {
      bad <- c(bad, "nFeature_max")
    }
  }
  if (length(bad) > 0) {
    stop(
      "Set numeric values in 'User filter settings' before running: ",
      paste(unique(bad), collapse = ", ")
    )
  }
  nFeature_min <- as.numeric(nFeature_min)
  mt_max <- as.numeric(mt_max)
  nFeature_max_out <- if (no_upper) NA_real_ else as.numeric(nFeature_max)
  if (!no_upper && nFeature_min >= nFeature_max_out) {
    stop("nFeature_min must be less than nFeature_max.")
  }
  list(
    mode = "fixed",
    nFeature_min = nFeature_min,
    nFeature_max = nFeature_max_out,
    mt_max = mt_max
  )
}

threshold_folder_label <- function(th) {
  max_label <- if (is.null(th$nFeature_max) || is.na(th$nFeature_max)) {
    "none"
  } else {
    th$nFeature_max
  }
  paste0(
    "nFmin", th$nFeature_min,
    "_nFmax", max_label,
    "_mtmax", th$mt_max
  )
}

filter_seurat_by_thresholds <- function(seurat_all, th) {
  md <- seurat_all@meta.data
  keep <- md$nFeature_RNA > th$nFeature_min &
    md$percent.mt < th$mt_max &
    md$scDblFinder.class == "singlet"
  if (!is.null(th$nFeature_max) && !is.na(th$nFeature_max)) {
    keep <- keep & md$nFeature_RNA < th$nFeature_max
  }
  seurat_all[, rownames(md)[keep]]
}

save_filter_zoom_box_plots <- function(
  md,
  th,
  plot_dir,
  stage,
  gene_label = "",
  y_min = 200,
  y_max = 3000
) {
  nF_min <- th$nFeature_min
  nF_max <- if (is.null(th$nFeature_max) || is.na(th$nFeature_max)) NULL else th$nFeature_max
  tag <- if (nzchar(gene_label)) paste0(gene_label, " — ") else ""

  p_nFeature <- plot_qc_box_zoom_vertical(
    md,
    "nFeature_RNA",
    paste0(tag, stage, " nFeature (zoom ", y_min, "–", y_max, ")"),
    "nFeature_RNA",
    y_min = y_min,
    y_max = y_max,
    hline_min = nF_min,
    hline_max = nF_max
  )
  p_nCount <- plot_qc_box_zoom_vertical(
    md,
    "nCount_RNA",
    paste0(tag, stage, " nCount (zoom ", y_min, "–", y_max, ")"),
    "nCount_RNA",
    y_min = y_min,
    y_max = y_max
  )

  suffix <- paste0("_vertical_", y_min, "_", y_max)
  save_gg(p_nFeature, file.path(plot_dir, paste0("box_nFeature", suffix)), width = 11, height = 7)
  save_gg(p_nCount, file.path(plot_dir, paste0("box_nCount", suffix)), width = 11, height = 7)
  invisible(plot_dir)
}

save_filter_qc_plots <- function(
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
  plot_dir <- file.path(out_dir, "qc_filter_plots", threshold_folder_label(th), stage)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  md$orig.ident <- factor(md$orig.ident, levels = sort(unique(md$orig.ident)))

  nF_min <- th$nFeature_min
  nF_max <- if (is.null(th$nFeature_max) || is.na(th$nFeature_max)) NULL else th$nFeature_max

  save_gg(
    plot_qc_density(
      md, "nFeature_RNA", paste(stage, "nFeature density (zoom)"), "nFeature_RNA",
      nF_min, nF_max, xlim = density_xlim_nFeature
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
  save_filter_zoom_box_plots(md, th, plot_dir, stage, y_min = qc_zoom_y_min, y_max = qc_zoom_y_max)
  save_gg(
    plot_qc_scatter(md, paste(stage, "nFeature vs nCount"), nF_min, nF_max),
    file.path(plot_dir, "scatter_nFeature_nCount"), width = 14, height = 10
  )

  filter_summary <- NULL
  if (identical(stage, "after_filter") && !is.null(seurat_filt) && !is.null(n_before)) {
    filter_summary <- save_filter_summary(seurat_filt, th, n_before, out_dir, md_before = md_before)
  }

  list(plot_dir = plot_dir, filter_summary = filter_summary)
}

save_filter_summary <- function(seurat_filt, th, n_before, out_dir, md_before = NULL) {
  label <- threshold_folder_label(th)
  summary_dir <- file.path(out_dir, "filter_summary", label)
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

  cells_per_sample <- as.data.frame(table(seurat_filt$orig.ident))
  colnames(cells_per_sample) <- c("sample", "n_cells")
  write.csv(
    cells_per_sample,
    file.path(summary_dir, "cells_per_sample_after_filter.csv"),
    row.names = FALSE
  )

  overall <- data.frame(
    threshold_label = label,
    mode = if (!is.null(th$mode)) th$mode else "fixed",
    nFeature_min = if (!is.null(th$nFeature_min)) th$nFeature_min else NA,
    nFeature_max = if (!is.null(th$nFeature_max)) th$nFeature_max else NA,
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

run_threshold_trial <- function(
  seurat_all,
  th,
  out_dir,
  md_before = NULL,
  density_xlim_nFeature = c(0, 6000),
  density_xlim_nCount = c(0, 15000),
  qc_zoom_y_min = 200,
  qc_zoom_y_max = 3000
) {
  label <- threshold_folder_label(th)
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

  do.call(save_filter_qc_plots, c(list(md = md_before, th = th, out_dir = out_dir, stage = "before_filter"), plot_args))
  n_before <- ncol(seurat_all)
  seurat_filt <- filter_seurat_by_thresholds(seurat_all, th)

  after <- do.call(
    save_filter_qc_plots,
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

nfeature_quantiles <- function(x) {
  stats::quantile(x, probs = c(0.25, 0.5, 0.75, 0.9), na.rm = TRUE)
}

round_nfeature_min <- function(x, floor_value = 200) {
  if (is.na(x) || !is.finite(x)) {
    return(floor_value)
  }
  x <- max(floor_value, x)
  step <- if (x < 500) {
    50
  } else if (x < 1500) {
    100
  } else if (x < 3000) {
    250
  } else {
    500
  }
  as.integer(ceiling(x / step) * step)
}

round_nfeature_max <- function(x) {
  if (is.na(x) || !is.finite(x)) {
    return(NA_integer_)
  }
  step <- if (x < 1500) {
    100
  } else if (x < 3500) {
    250
  } else {
    500
  }
  as.integer(max(step, floor(x / step) * step))
}

nfeature_bounds_for_method <- function(median_value, spread_value, multiplier, floor_value = 200) {
  raw_lower <- max(floor_value, median_value - multiplier * spread_value)
  raw_upper <- median_value + multiplier * spread_value
  list(
    raw_lower = raw_lower,
    raw_upper = raw_upper,
    round_min = round_nfeature_min(raw_lower, floor_value = floor_value),
    round_max = round_nfeature_max(raw_upper)
  )
}

nfeature_window_retention <- function(md, nFeature_min, nFeature_max, apply_singlet_mt = FALSE, mt_max = NULL) {
  x <- md$nFeature_RNA
  keep <- x > nFeature_min & x < nFeature_max
  if (isTRUE(apply_singlet_mt)) {
    if ("percent.mt" %in% colnames(md) && !is.null(mt_max)) {
      keep <- keep & md$percent.mt < mt_max
    }
    if ("scDblFinder.class" %in% colnames(md)) {
      keep <- keep & md$scDblFinder.class == "singlet"
    }
  }
  n_kept <- sum(keep)
  n_total <- length(x)
  list(
    n_kept = n_kept,
    n_total = n_total,
    pct_kept = round(100 * n_kept / n_total, 2)
  )
}

summarize_nfeature_distribution <- function(md) {
  build_row <- function(sub, scope_label, sample_label) {
    xs <- sub$nFeature_RNA
    qs <- nfeature_quantiles(xs)
    data.frame(
      scope = scope_label,
      sample = sample_label,
      n_cells_total = nrow(sub),
      nFeature_median = unname(qs["50%"]),
      nFeature_Q25 = unname(qs["25%"]),
      nFeature_Q75 = unname(qs["75%"]),
      nFeature_Q90 = unname(qs["90%"]),
      nFeature_sd = stats::sd(xs),
      nFeature_mad = stats::mad(xs),
      n_below_200 = sum(xs < 200),
      pct_below_200 = round(100 * sum(xs < 200) / nrow(sub), 2),
      n_at_or_above_Q75 = sum(xs >= unname(qs["75%"])),
      pct_at_or_above_Q75 = round(100 * sum(xs >= unname(qs["75%"])) / nrow(sub), 2),
      stringsAsFactors = FALSE
    )
  }

  overall <- build_row(md, "overall", "all")
  per_sample <- do.call(rbind, lapply(split(md, md$orig.ident, drop = TRUE), function(sub) {
    build_row(sub, "per_sample", unique(as.character(sub$orig.ident)))
  }))
  rbind(overall, per_sample)
}

summarize_nfeature_data_driven_candidates <- function(md, multipliers = 2:5, floor_value = 200) {
  build_rows <- function(sub, scope_label, sample_label) {
    xs <- sub$nFeature_RNA
    med <- stats::median(xs)
    sdv <- stats::sd(xs)
    madv <- stats::mad(xs)
    rows <- list()

    for (k in multipliers) {
      for (method_info in list(
        list(method = "median -/+ k * SD", spread = sdv),
        list(method = "median -/+ k * MAD", spread = madv)
      )) {
        bounds <- nfeature_bounds_for_method(med, method_info$spread, k, floor_value = floor_value)
        raw_ret <- nfeature_window_retention(sub, bounds$raw_lower, bounds$raw_upper)
        round_ret <- nfeature_window_retention(sub, bounds$round_min, bounds$round_max)
        rows[[length(rows) + 1L]] <- data.frame(
          scope = scope_label,
          sample = sample_label,
          method = method_info$method,
          multiplier = k,
          nFeature_median = med,
          spread = method_info$spread,
          raw_lower = round(bounds$raw_lower, 1),
          raw_upper = round(bounds$raw_upper, 1),
          round_min = bounds$round_min,
          round_max = bounds$round_max,
          n_cells_raw_window = raw_ret$n_kept,
          pct_raw_window = raw_ret$pct_kept,
          n_cells_round_window = round_ret$n_kept,
          pct_round_window = round_ret$pct_kept,
          stringsAsFactors = FALSE
        )
      }
    }

    rows[[length(rows) + 1L]] <- data.frame(
      scope = scope_label,
      sample = sample_label,
      method = "count nFeature below floor",
      multiplier = NA_real_,
      nFeature_median = med,
      spread = NA_real_,
      raw_lower = NA_real_,
      raw_upper = floor_value,
      round_min = NA_integer_,
      round_max = floor_value,
      n_cells_raw_window = sum(xs < floor_value),
      pct_raw_window = round(100 * sum(xs < floor_value) / length(xs), 2),
      n_cells_round_window = sum(xs < floor_value),
      pct_round_window = round(100 * sum(xs < floor_value) / length(xs), 2),
      stringsAsFactors = FALSE
    )

    do.call(rbind, rows)
  }

  overall <- build_rows(md, "overall", "all")
  per_sample <- do.call(rbind, lapply(split(md, md$orig.ident, drop = TRUE), function(sub) {
    build_rows(sub, "per_sample", unique(as.character(sub$orig.ident)))
  }))
  rbind(overall, per_sample)
}

summarize_nfeature_round_grid <- function(
  md,
  min_candidates = c(200, 300, 400, 500, 600, 800, 1000),
  max_candidates = c(2000, 2500, 3000, 3500, 4000, 5000, 6000)
) {
  rows <- list()
  for (nmin in min_candidates) {
    for (nmax in max_candidates) {
      if (nmin >= nmax) {
        next
      }
      overall <- nfeature_window_retention(md, nmin, nmax)
      per_sample_pct <- vapply(split(md, md$orig.ident, drop = FALSE), function(s) {
        nfeature_window_retention(s, nmin, nmax)$pct_kept
      }, numeric(1))
      rows[[length(rows) + 1L]] <- data.frame(
        nFeature_min = nmin,
        nFeature_max = nmax,
        n_cells_overall = overall$n_kept,
        pct_retained_overall = overall$pct_kept,
        min_pct_retained_any_sample = round(min(per_sample_pct), 2),
        max_pct_retained_any_sample = round(max(per_sample_pct), 2),
        stringsAsFactors = FALSE
      )
    }
  }
  grid <- do.call(rbind, rows)
  grid[order(grid$nFeature_min, grid$nFeature_max), , drop = FALSE]
}

summarize_user_nfeature_filter <- function(md, th) {
  x <- md$nFeature_RNA
  n_total <- length(x)
  in_nfeature_window <- x > th$nFeature_min & x < th$nFeature_max
  pass_mt <- md$percent.mt < th$mt_max
  pass_singlet <- md$scDblFinder.class == "singlet"
  pass_all <- in_nfeature_window & pass_mt & pass_singlet

  build_row <- function(sub, scope_label, sample_label) {
    xs <- sub$nFeature_RNA
    in_win <- xs > th$nFeature_min & xs < th$nFeature_max
    pass_all_s <- in_win & sub$percent.mt < th$mt_max & sub$scDblFinder.class == "singlet"
    data.frame(
      scope = scope_label,
      sample = sample_label,
      nFeature_min = th$nFeature_min,
      nFeature_max = th$nFeature_max,
      mt_max = th$mt_max,
      n_cells_total = nrow(sub),
      n_cells_pass_nFeature_only = sum(in_win),
      n_cells_pass_all_filters = sum(pass_all_s),
      pct_pass_all = round(100 * sum(pass_all_s) / nrow(sub), 2),
      n_excluded_low_nFeature = sum(xs <= th$nFeature_min),
      n_excluded_high_nFeature = sum(xs >= th$nFeature_max),
      stringsAsFactors = FALSE
    )
  }

  rbind(
    build_row(md, "overall", "all"),
    do.call(rbind, lapply(split(md, md$orig.ident, drop = TRUE), function(sub) {
      build_row(sub, "per_sample", unique(as.character(sub$orig.ident)))
    }))
  )
}

run_nfeature_threshold_diagnostics <- function(md, out_dir, user_filter = NULL) {
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

  user_impact <- NULL
  if (!is.null(user_filter)) {
    user_impact <- summarize_user_nfeature_filter(md, user_filter)
    write.csv(user_impact, file.path(diag_dir, "user_filter_impact.csv"), row.names = FALSE)
  }

  list(
    diag_dir = diag_dir,
    distribution = distribution,
    sd_mad_bounds = sd_mad_bounds,
    round_grid = round_grid,
    user_impact = user_impact
  )
}

