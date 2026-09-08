# Minimal end-to-end example for residual Moran diagnostics.
#
# Run from the repository root:
#   Rscript dev/moran_toy_example.R


required_packages <- c("FNN", "Matrix")
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

# Source the development versions rather than an installed package.
source("R/DE.R")
source("R/getPatches.R")
source("R/moran.R")

set.seed(42)

n_per_patch <- 12L
patch <- rep(c("patch_1", "patch_2"), each = n_per_patch)
exposure <- rep(seq(-1, 1, length.out = n_per_patch), 2)

xy <- cbind(
  x = rep(seq_len(n_per_patch), 2),
  y = rep(c(0, 5), each = n_per_patch)
)

# The three neighboring blocks create nonlinear spatial structure that is not
# explained by the linear exposure term.
spatial_signal <- rep(
  c(rep(2, 4), rep(-2, 4), rep(2, 4)),
  2
)

y <- cbind(
  spatial_gene = 10 + 2 * exposure + spatial_signal,
  noise_gene = 10 - exposure + stats::rnorm(length(exposure))
)

fit <- patchDE(
  y = y,
  df = data.frame(exposure),
  patch = patch,
  return_residuals = TRUE,
  verbose = FALSE
)

set.seed(42)
moran_result <- moranTest(
  residuals = fit$residuals,
  xy = xy,
  patch = patch,
  k = 2,
  n_permutations = 999,
  alternative = "greater",
  p_adjust_method = "BH",
  adjustment_scope = "global"
)

print(
  moran_result[, c(
    "patch", "gene", "observed",
    "p_value", "p_adjusted", "status"
  )],
  row.names = FALSE
)

# Create a compact visual explanation of the simulated layout, residuals, and
# resulting Moran statistics. Override the destination with
# SPACEMOSAIC_MORAN_PLOT if desired.
plot_path <- Sys.getenv(
  "SPACEMOSAIC_MORAN_PLOT",
  unset = file.path("dev", "runs", "moran_toy_example.pdf")
)
dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)

residual_palette <- grDevices::colorRampPalette(
  c("#2166AC", "white", "#B2182B")
)(101)
residual_limit <- max(abs(fit$residuals))

residual_colors <- function(values) {
  color_index <- round((values + residual_limit) /
                         (2 * residual_limit) * 100) + 1L
  color_index <- pmax(1L, pmin(101L, color_index))
  residual_palette[color_index]
}

plot_residuals <- function(values, title) {
  graphics::plot(
    xy,
    type = "n",
    asp = 1,
    xlab = "x",
    ylab = "y",
    main = title
  )
  graphics::points(
    xy,
    pch = 21,
    bg = residual_colors(values),
    col = "grey30",
    cex = 1.8
  )
  graphics::text(
    xy,
    labels = seq_len(nrow(xy)),
    pos = 3,
    cex = 0.55
  )
}

grDevices::pdf(plot_path, width = 11, height = 8.5)
graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

patch_colors <- c(patch_1 = "#1B9E77", patch_2 = "#D95F02")
graphics::plot(
  xy,
  type = "n",
  asp = 1,
  xlab = "x",
  ylab = "y",
  main = "Cells, patches, and k-NN graph"
)
for (patch_name in unique(patch)) {
  cell_index <- which(patch == patch_name)
  W <- .buildContiguityGraph(xy[cell_index, , drop = FALSE], k = 2)
  edges <- Matrix::summary(W)
  edges <- edges[edges$i < edges$j, , drop = FALSE]
  graphics::segments(
    x0 = xy[cell_index[edges$i], 1],
    y0 = xy[cell_index[edges$i], 2],
    x1 = xy[cell_index[edges$j], 1],
    y1 = xy[cell_index[edges$j], 2],
    col = "grey75"
  )
}
graphics::points(
  xy,
  pch = 21,
  bg = patch_colors[patch],
  col = "grey20",
  cex = 1.8
)
graphics::text(xy, labels = seq_len(nrow(xy)), pos = 3, cex = 0.55)
graphics::legend(
  "topright",
  legend = names(patch_colors),
  pt.bg = patch_colors,
  pch = 21,
  bty = "n"
)

plot_residuals(
  fit$residuals[, "spatial_gene"],
  "Residuals: spatial gene"
)
plot_residuals(
  fit$residuals[, "noise_gene"],
  "Residuals: noise gene"
)

observed_matrix <- stats::xtabs(observed ~ gene + patch, data = moran_result)
adjusted_matrix <- stats::xtabs(p_adjusted ~ gene + patch, data = moran_result)
expected_i <- unique(moran_result$expected)
bar_limits <- range(c(observed_matrix, expected_i, 0))
bar_positions <- graphics::barplot(
  observed_matrix,
  beside = TRUE,
  col = c("#7570B3", "#E7298A"),
  ylim = bar_limits + c(-0.15, 0.2),
  ylab = "Moran's I",
  main = "Observed Moran's I by gene and patch",
  legend.text = rownames(observed_matrix),
  args.legend = list(x = "bottomright", bty = "n", cex = 0.75)
)
graphics::abline(h = expected_i, lty = 2, col = "grey30")
graphics::text(
  bar_positions,
  observed_matrix,
  labels = sprintf("q=%.3f", adjusted_matrix),
  pos = ifelse(observed_matrix >= 0, 3, 1),
  cex = 0.6
)

grDevices::dev.off()
message("Diagnostic plot written to: ", normalizePath(plot_path))

spatial_rows <- moran_result$gene == "spatial_gene"
noise_rows <- moran_result$gene == "noise_gene"

# These checks make the example useful as a quick manual smoke test while the
# formal regression tests remain in tests/testthat/.
stopifnot(
  all(moran_result$status == "ok"),
  all(moran_result$observed[spatial_rows] > moran_result$observed[noise_rows]),
  all(moran_result$p_adjusted[spatial_rows] < 0.05)
)

message("Moran diagnostic toy example completed successfully.")
