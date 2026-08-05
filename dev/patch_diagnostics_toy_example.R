# Real-data example for comparing patch diagnostics across npatches settings.
#
# This reproduces the mouse-colon vignette's "Immune Surveillance" setup on
# the data distributed with SpaceMosaic. Run from the repository root with:
#   Rscript dev/patch_diagnostics_toy_example.R
#
# When sourced in an interactive R session, the objects created below can also
# be passed to runInteractivePlotter(); see the final commented example.


required_packages <- c(
  "cli", "FNN", "ggplot2", "HDF5Array", "igraph", "Matrix",
  "patchwork", "SingleCellExperiment", "SpatialExperiment"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install required packages first: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

if (!file.exists("DESCRIPTION") ||
    read.dcf("DESCRIPTION", fields = "Package")[[1L]] != "SpaceMosaic") {
  stop("Run this script from the SpaceMosaic repository root.", call. = FALSE)
}

# Source development versions rather than an installed package.
source("R/neighbors.R")
source("R/getPatches.R")
source("R/plotting.R")
source("R/interactive_plotting.R")
source("R/patchDiagnostics.R")

# Parameters intended for experimentation. Five iterations keep this real-data
# development example reasonably quick; use 15 for the package default.
npatches_grid <- c(50L, 100L, 150L)
n_iters <- 5L
diagnostic_k <- 10L
strict_diagnostic_k <- 5L
seed <- 123L

# Load the same packaged data used by vignettes/mouse_colon.Rmd.
spe <- HDF5Array::loadHDF5SummarizedExperiment(
  file.path("inst", "extdata", "mouse_colon")
)
xy <- SpatialExperiment::spatialCoords(spe)
rownames(xy) <- spe$cell_id
cell_types <- spe$tier1

# Reproduce the vignette's local immune-cell density predictor.
is_epithelial <- cell_types == "Epithelial"
is_immune <- cell_types == "Immune"
k_neighbors <- 100L
nn_indices <- FNN::get.knn(as.matrix(xy), k = k_neighbors)$nn.index
immune_neighbors <- rowSums(matrix(is_immune[nn_indices], nrow = nrow(xy)))

# Reproduce the vignette's cellular-neighborhood embedding.
Z <- embedCellNeighborhoods(
  mat = SingleCellExperiment::reducedDim(spe, "PCA"),
  xy = xy,
  ks = c(5, 50)
)

xy_use <- xy[is_epithelial, , drop = FALSE]
X_use <- cbind(immune_neighbors = immune_neighbors[is_epithelial])
rownames(X_use) <- rownames(xy_use)
Z_use <- Z[is_epithelial, , drop = FALSE]

# Fit and diagnose the same biological problem at several requested patch
# counts. log_iters = TRUE is needed for membership_stability.
fits <- vector("list", length(npatches_grid))
diagnostics <- vector("list", length(npatches_grid))
names(fits) <- names(diagnostics) <- paste0("npatches_", npatches_grid)

for (i in seq_along(npatches_grid)) {
  requested <- npatches_grid[i]
  message("Running getPatches() with npatches = ", requested, " ...")
  set.seed(seed)
  fits[[i]] <- getPatches(
    xy = xy_use,
    X = X_use,
    Z = Z_use,
    npatches = requested,
    n_iters = n_iters,
    log_iters = TRUE,
    verbose = FALSE
  )
  diagnostics[[i]] <- getPatchDiagnostics(
    xy = xy_use,
    X = X_use,
    patches = fits[[i]],
    k = diagnostic_k,
    strict_k = strict_diagnostic_k
  )
}

# Combine results into tables convenient for parameter comparison.
patch_diagnostics <- do.call(rbind, lapply(seq_along(diagnostics), function(i) {
  out <- diagnostics[[i]]$patch_diagnostics
  out$npatches_requested <- npatches_grid[i]
  out
}))
rownames(patch_diagnostics) <- NULL

assignment_summary <- do.call(rbind, lapply(seq_along(diagnostics), function(i) {
  out <- diagnostics[[i]]$assignment_summary
  out$npatches_requested <- npatches_grid[i]
  out
}))
rownames(assignment_summary) <- NULL

connectivity_curve <- do.call(rbind, lapply(seq_along(diagnostics), function(i) {
  out <- diagnostics[[i]]$connectivity_curve
  out$npatches_requested <- npatches_grid[i]
  out
}))
rownames(connectivity_curve) <- NULL

polygon_data <- do.call(rbind, lapply(seq_along(fits), function(i) {
  patch_data <- diagnostics[[i]]$patch_diagnostics
  out <- getPatchPolys(
    xy = xy_use,
    patch = fits[[i]]$patch,
    patch_data = patch_data
  )
  out$npatches_requested <- npatches_grid[i]
  out
}))
rownames(polygon_data) <- NULL

message("Parameter comparison:")
print(assignment_summary, row.names = FALSE)

# Save reusable numerical results without adding artifacts to the repository.
# Set SPACEMOSAIC_PATCH_DIAGNOSTICS_DIR to choose a permanent output folder.
output_dir <- Sys.getenv(
  "SPACEMOSAIC_PATCH_DIAGNOSTICS_DIR",
  unset = file.path(
    Sys.getenv("TMPDIR", unset = tempdir()),
    "SpaceMosaic_patch_diagnostics"
  )
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results_path <- file.path(output_dir, "mouse_colon_patch_diagnostics.rds")
saveRDS(
  list(
    npatches_grid = npatches_grid,
    n_iters = n_iters,
    diagnostic_k = diagnostic_k,
    strict_diagnostic_k = strict_diagnostic_k,
    fits = fits,
    patch_diagnostics = patch_diagnostics,
    assignment_summary = assignment_summary,
    connectivity_curve = connectivity_curve,
    polygon_data = polygon_data
  ),
  results_path
)

# 1. Spatial view: where the patches lie relative to immune-cell exposure.
cell_plot_data <- data.frame(
  x = xy_use[, 1],
  y = xy_use[, 2],
  immune_neighbors = X_use[, "immune_neighbors"]
)
context_plot_data <- data.frame(
  x = xy[, 1],
  y = xy[, 2],
  is_epithelial = is_epithelial,
  is_immune = is_immune
)
spatial_plot <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = context_plot_data[!context_plot_data$is_epithelial, ],
    mapping = ggplot2::aes(x = x, y = y),
    color = "grey80",
    size = 0.035
  ) +
  ggplot2::geom_point(
    data = context_plot_data[context_plot_data$is_immune, ],
    mapping = ggplot2::aes(x = x, y = y),
    color = "#2166AC",
    size = 0.05
  ) +
  ggplot2::geom_point(
    data = cell_plot_data,
    mapping = ggplot2::aes(x = x, y = y, color = immune_neighbors),
    size = 0.045
  ) +
  ggplot2::geom_polygon(
    data = polygon_data,
    mapping = ggplot2::aes(x = x, y = y, group = patch),
    inherit.aes = FALSE,
    fill = NA,
    color = "#39FF14",
    linewidth = 0.22
  ) +
  ggplot2::facet_wrap(~npatches_requested, nrow = 1) +
  ggplot2::coord_fixed() +
  ggplot2::scale_color_viridis_c(option = "plasma") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::labs(
    title = "Mouse-colon epithelial patches across npatches settings",
    subtitle = paste0(
      "Green lines are patch hulls; immune cells are blue; n_iters = ",
      n_iters
    ),
    color = "Immune cells\namong 100-NN",
    x = NULL,
    y = NULL
  )

# 2. Core cross-patch diagnostic: size versus within-patch X variation. Color
# reveals final-iteration stability and shape flags spatial fragmentation.
diagnostic_plot_data <- patch_diagnostics
diagnostic_plot_data$connectivity <- factor(
  ifelse(
    diagnostic_plot_data$strict_component_fraction < 1,
    "Fragmented",
    "Connected"
  ),
  levels = c("Connected", "Fragmented")
)
fragmentation_labels <- vapply(npatches_grid, function(requested) {
  setting_data <- diagnostic_plot_data[
    diagnostic_plot_data$npatches_requested == requested,
  ]
  n_fragmented <- sum(setting_data$connectivity == "Fragmented")
  sprintf(
    "npatches = %d\nFragmented at k = %d: %d/%d (%.0f%%)",
    requested,
    strict_diagnostic_k,
    n_fragmented,
    nrow(setting_data),
    100 * n_fragmented / nrow(setting_data)
  )
}, character(1L))
names(fragmentation_labels) <- as.character(npatches_grid)

diagnostic_plot <- ggplot2::ggplot(
  diagnostic_plot_data,
  ggplot2::aes(x = n_cells, y = x_sd)
) +
  ggplot2::geom_point(
    ggplot2::aes(
      color = membership_stability,
      shape = connectivity
    ),
    size = 2
  ) +
  ggplot2::facet_wrap(
    ~npatches_requested,
    scales = "free_x",
    labeller = ggplot2::as_labeller(fragmentation_labels)
  ) +
  ggplot2::scale_color_viridis_c(option = "viridis", limits = c(0, 1)) +
  ggplot2::scale_shape_manual(
    values = c(Connected = 16, Fragmented = 4)
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::labs(
    title = "Core patch diagnostics",
    subtitle = "Color shows final-iteration stability; crosses mark fragmented patches",
    x = "Cells per patch",
    y = "Within-patch SD of immune-neighbor count",
    color = "Final membership\nstability",
    shape = "Spatial connectivity"
  )

# 3. Connectivity curves: faint lines are individual patches and the dark line
# is the median within each parameter setting.
connectivity_curve_plot <- ggplot2::ggplot(
  connectivity_curve,
  ggplot2::aes(x = k, y = component_fraction, group = patch)
) +
  ggplot2::geom_line(color = "grey35", alpha = 0.25, linewidth = 0.35) +
  ggplot2::stat_summary(
    mapping = ggplot2::aes(group = 1),
    fun = stats::median,
    geom = "line",
    color = "#0072B2",
    linewidth = 1.2
  ) +
  ggplot2::geom_vline(
    xintercept = strict_diagnostic_k,
    linetype = "dashed",
    color = "#D55E00"
  ) +
  ggplot2::facet_wrap(~npatches_requested, nrow = 1) +
  ggplot2::scale_x_continuous(breaks = seq_len(diagnostic_k)) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::labs(
    title = "Patch connectivity across k",
    subtitle = paste0(
      "Grey lines are patches; blue is the median; dashed line marks strict k = ",
      strict_diagnostic_k
    ),
    x = "Number of spatial neighbors (k)",
    y = "Largest-component fraction"
  )

# 4. Distribution view: useful for identifying tails and comparing parameter
# settings without imposing universal pass/fail thresholds.
distribution_data <- rbind(
  data.frame(
    npatches_requested = patch_diagnostics$npatches_requested,
    metric = "Cells per patch",
    value = patch_diagnostics$n_cells
  ),
  data.frame(
    npatches_requested = patch_diagnostics$npatches_requested,
    metric = "SD of immune-neighbor count",
    value = patch_diagnostics$x_sd
  ),
  data.frame(
    npatches_requested = patch_diagnostics$npatches_requested,
    metric = "Membership stability",
    value = patch_diagnostics$membership_stability
  )
)
distribution_data$npatches_requested <- factor(
  distribution_data$npatches_requested,
  levels = npatches_grid
)
distribution_colors <- grDevices::hcl.colors(
  length(npatches_grid),
  palette = "Dark 3"
)
names(distribution_colors) <- as.character(npatches_grid)
distribution_linetypes <- rep(
  c("solid", "dashed", "dotdash", "dotted", "longdash", "twodash"),
  length.out = length(npatches_grid)
)
names(distribution_linetypes) <- as.character(npatches_grid)

distribution_plot <- ggplot2::ggplot(
  distribution_data,
  ggplot2::aes(
    x = value,
    color = npatches_requested,
    linetype = npatches_requested
  )
) +
  ggplot2::geom_density(
    linewidth = 1,
    na.rm = TRUE
  ) +
  ggplot2::facet_wrap(~metric, scales = "free", ncol = 3) +
  ggplot2::scale_color_manual(
    values = distribution_colors
  ) +
  ggplot2::scale_linetype_manual(
    values = distribution_linetypes
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::labs(
    title = "Diagnostic distributions across parameter settings",
    x = NULL,
    y = "Density",
    color = "Requested patches",
    linetype = "Requested patches"
  )

# 5. Tuning view: requested versus surviving non-empty patches, with assignment
# rate shown explicitly instead of hidden in a pass/fail label.
tuning_plot <- ggplot2::ggplot(
  assignment_summary,
  ggplot2::aes(x = npatches_requested, y = n_nonempty_patches)
) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_line(color = "#2C7FB8") +
  ggplot2::geom_point(
    color = "#2C7FB8",
    size = 4
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f%% assigned", 100 * fraction_assigned)),
    vjust = -0.8,
    size = 3.5
  ) +
  ggplot2::scale_x_continuous(
    breaks = npatches_grid,
    expand = ggplot2::expansion(mult = 0.18)
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.08, 0.18))
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::labs(
    title = "Requested vs non-empty patches",
    x = "npatches requested",
    y = "Non-empty patches"
  )

# Keep the overview plots compact: the two detailed patch plots share one page,
# while the tuning plot occupies a smaller centered panel below the densities.
diagnostic_connectivity_page <- patchwork::wrap_plots(
  diagnostic_plot,
  connectivity_curve_plot,
  ncol = 1,
  heights = c(1, 1.05)
)
compact_tuning_row <- patchwork::wrap_plots(
  patchwork::plot_spacer(),
  tuning_plot,
  patchwork::plot_spacer(),
  nrow = 1,
  widths = c(0.45, 1, 0.45)
)
distribution_tuning_page <- patchwork::wrap_plots(
  distribution_plot,
  compact_tuning_row,
  ncol = 1,
  heights = c(2.1, 1)
)

plot_path <- file.path(output_dir, "mouse_colon_patch_diagnostics.pdf")
grDevices::pdf(plot_path, width = 13, height = 9, onefile = TRUE)
print(spatial_plot)
print(diagnostic_connectivity_page)
print(distribution_tuning_page)
grDevices::dev.off()

message("Results written to: ", normalizePath(results_path))
message("Three-page diagnostic plot written to: ", normalizePath(plot_path))

# Select the middle setting for optional interactive exploration. This is not
# launched automatically because runInteractivePlotter() starts a blocking
# Shiny session. After source("dev/patch_diagnostics_toy_example.R"), run:
#
# selected_name <- paste0("npatches_", npatches_grid[ceiling(length(npatches_grid) / 2)])
# runInteractivePlotter(
#   spe = spe,
#   patches = fits[[selected_name]]$patch,
#   patch_data = diagnostics[[selected_name]]$patch_diagnostics
# )

stopifnot(
  nrow(assignment_summary) == length(npatches_grid),
  all(is.finite(patch_diagnostics$membership_stability)),
  nrow(connectivity_curve) == nrow(patch_diagnostics) * diagnostic_k,
  nrow(polygon_data) > 0L
)

message("Mouse-colon patch diagnostic example completed successfully.")
