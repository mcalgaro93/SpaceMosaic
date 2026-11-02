#' Plot cell-level Moran diagnostics
#'
#' Create residual maps, local Moran maps, quadrant maps and Moran scatterplots from
#' [getMoranLocal()] results for one complete gene-patch combination. These
#' functions do not rebuild the spatial graph or recalculate local statistics.
#' Each returns a `ggplot` object that can be extended with layers, scales,
#' labels and themes. The optional package `ggplot2` is required only when
#' plotting.
#'
#' @param data Data frame returned by [getMoranLocal()], subset to one gene and
#'   one complete patch. All rows must have status `"ok"`. Do not filter out
#'   individual cells: their lag values refer to the complete original patch.
#' @param standardize Logical; divide centered residuals by
#'   `sqrt(mean(residual_centered^2))` for map colors or scatterplot coordinates
#'   only. In the scatterplot, divide the spatial lag by the same residual
#'   scale. For `plotMoranLocal()`, instead use `lag_scaled` to express local
#'   statistics on a scale whose within-patch mean equals global Moran's I.
#'   Default FALSE.
#' @param size Non-negative numeric point size, passed to `ggplot2`.
#' @param alpha Numeric point opacity between zero and one.
#' @param background Optional data frame with finite numeric columns `x`, `y`.
#'   Drawn in grey behind the selected cells. Additional columns are ignored.
#' @param background_size Non-negative numeric point size for background
#'   cells, independent of `size`. Default 0.5.
#' @param background_alpha Numeric background point opacity between zero and
#'   one, independent of `alpha`. Default 0.3.
#' @param polygons Optional output of [getPatchPolys()], with columns `x`, `y`
#'   and `patch`. Only vertices matching the selected patch are drawn.
#' @param aspect Spatial map aspect ratio: `"fixed"` (default) preserves equal
#'   distances on both axes; `"free"` allows the panel to fill the available
#'   space, potentially stretching the spatial geometry.
#'
#' @details
#' `plotMoranResiduals()` maps centered residuals with a diverging fill scale
#' centered at zero. With `standardize = TRUE`, only the plotting variable is
#' rescaled; local Moran statistics and quadrants are unchanged. Plot data
#' retain the input columns and include the additional column `residual_plot`.
#'
#' `plotMoranLocal()` maps the unmodified `local_i` values by default, with a diverging
#' fill scale centered at zero and symmetric limits. Positive values indicate
#' concordant residual and spatial lag signs (High-High or Low-Low); negative
#' values indicate discordant signs (High-Low or Low-High). Colors describe
#' local association, not statistical significance. Values are not restricted
#' to the interval from -1 to 1. With the binary, non-row-standardized weights used by
#' [getMoranLocal()], magnitudes also depend on the number of graph neighbors.
#' With `standardize = TRUE`, plot
#' `residual_centered * lag_scaled / mean(residual_centered^2)` instead.
#' This multiplies `local_i` by `n / S0`, the inverse mean number of graph
#' neighbors, so the mean displayed value equals the observed global Moran's I.
#' It does not normalize weights separately for each cell or remove differences
#' in graph geometry, and displayed values can still exceed the interval from
#' -1 to 1. Plot data retain the input columns and include `local_i_plot`,
#' the displayed values; input local statistics and quadrants are unchanged.
#' Use common fill limits when comparing maps across genes or patches.
#'
#' `plotMoranQuadrants()` maps descriptive quadrants with the same palette used
#' by `plotMoranScatter()`. Cells with a missing quadrant (on a zero axis) are
#' grey. Quadrant colors do not indicate statistical significance.
#'
#' `plotMoranScatter()` plots `residual_centered` against `lag_scaled`, with
#' reference axes at zero and a line through the origin. Its slope is
#' `sum(residual_centered * lag_scaled) / sum(residual_centered^2)`, the observed
#' global Moran's I for the supplied complete patch. A complete graph can be
#' plotted, but its global Moran's I cannot distinguish spatial arrangements.
#' With `standardize = TRUE`, both axes are divided by
#' `sqrt(mean(residual_centered^2))`. This preserves the slope, quadrants and
#' local statistics while expressing coordinates in residual standard deviation
#' units for comparisons across genes or patches. The lag is not divided by its
#' own standard deviation. Plot data retain the input columns and include
#' `residual_plot` and `lag_plot`, the displayed coordinates. For visual
#' comparisons, also use common axis limits and panel proportions.
#'
#' Spatial maps use a fixed aspect ratio by default; set `aspect = "free"` to
#' allow flexible panel proportions. Both options zoom to the selected cells,
#' even when a larger background is supplied. All four plots use point fill,
#' leaving point color available for borders. Replace the fill scale to customize the
#' palette. Multiple panels can be composed externally; `patchwork` is not
#' required by these functions. The scatterplot leaves the aspect ratio free.
#' No automatic titles or subtitles are added; use `ggplot2::labs()` to supply
#' them when needed. Axis and legend labels are provided by default.
#'
#' @return A `ggplot` object. The input data frame is not modified.
#' @seealso [getMoranLocal()], [getPatchPolys()]
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   xy <- cbind(x = c(0, 1, 3, 7), y = c(0, 1, 0, 1))
#'   local <- getMoranLocal(c(-3, -1, 1, 3), xy, k = 1, verbose = FALSE)
#'   polygons <- getPatchPolys(xy, rep("all", 4))
#'   plotMoranResiduals(local, standardize = TRUE, polygons = polygons)
#'   plotMoranLocal(local, polygons = polygons)
#'   plotMoranLocal(local, standardize = TRUE, polygons = polygons)
#'   plotMoranQuadrants(local)
#'   plotMoranQuadrants(local, aspect = "free")
#'   plotMoranScatter(local) + ggplot2::labs(title = "Residual spatial structure")
#'   plotMoranScatter(local, standardize = TRUE)
#' }
#' @export
plotMoranResiduals <- function(data, standardize = FALSE, size = 1, alpha = 1,
                               background = NULL, polygons = NULL,
                               background_size = 0.5, background_alpha = 0.3,
                               aspect = c("fixed", "free")) {
  aspect <- match.arg(aspect)
  data <- .prepareMoranPlotData(data, c("x", "y"), size, alpha)
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) {
    stop("standardize must be TRUE or FALSE.")
  }
  data$residual_plot <- data$residual_centered
  if (standardize) {
    data$residual_plot <- data$residual_plot / sqrt(mean(data$residual_plot^2))
  }
  limit <- max(abs(data$residual_plot))

  .moranMapBase(data, background, polygons, background_size, background_alpha, aspect) +
    ggplot2::geom_point(
      ggplot2::aes(fill = !!as.name("residual_plot")),
      shape = 21, color = "grey30", stroke = 0.15, size = size, alpha = alpha
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-limit, limit)
    ) +
    ggplot2::labs(
      fill = if (standardize) "Scaled centered residual" else "Centered residual"
    )
}


#' @rdname plotMoranResiduals
#' @export
plotMoranLocal <- function(data, standardize = FALSE, size = 1, alpha = 1,
                           background = NULL, polygons = NULL,
                           background_size = 0.5, background_alpha = 0.3,
                           aspect = c("fixed", "free")) {
  aspect <- match.arg(aspect)
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) {
    stop("standardize must be TRUE or FALSE.")
  }
  columns <- c("x", "y", "local_i", if (standardize) "lag_scaled")
  data <- .prepareMoranPlotData(data, columns, size, alpha)
  data$local_i_plot <- data$local_i
  if (standardize) {
    m2 <- mean(data$residual_centered^2)
    data$local_i_plot <- data$residual_centered * data$lag_scaled / m2
  }
  limit <- max(abs(data$local_i_plot))
  # Keep a non-degenerate scale when every local statistic is zero.
  if (limit == 0) limit <- 1

  .moranMapBase(data, background, polygons, background_size, background_alpha, aspect) +
    ggplot2::geom_point(
      ggplot2::aes(fill = !!as.name("local_i_plot")),
      shape = 21, color = "grey30", stroke = 0.15, size = size, alpha = alpha
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-limit, limit)
    ) +
    ggplot2::labs(
      fill = if (standardize) "Scaled local Moran's I" else "Local Moran's I"
    )
}


#' @rdname plotMoranResiduals
#' @export
plotMoranQuadrants <- function(data, size = 1, alpha = 1,
                               background = NULL, polygons = NULL,
                               background_size = 0.5, background_alpha = 0.3,
                               aspect = c("fixed", "free")) {
  aspect <- match.arg(aspect)
  data <- .prepareMoranPlotData(data, c("x", "y", "quadrant"), size, alpha)
  colors <- .moranQuadrantColors()

  .moranMapBase(data, background, polygons, background_size, background_alpha, aspect) +
    ggplot2::geom_point(
      ggplot2::aes(fill = !!as.name("quadrant")),
      shape = 21, color = "grey30", stroke = 0.15, size = size, alpha = alpha
    ) +
    ggplot2::scale_fill_manual(
      values = colors, limits = names(colors), breaks = names(colors),
      drop = FALSE, na.value = "grey70", 
      guide = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1), 
        theme = ggplot2::theme(legend.margin = ggplot2::margin(0, 0, 0, 0)))
    ) +
    ggplot2::labs(fill = "Residual-lag quadrant")
}


#' @rdname plotMoranResiduals
#' @export
plotMoranScatter <- function(data, standardize = FALSE, size = 1, alpha = 0.7) {
  data <- .prepareMoranPlotData(data, c("lag_scaled", "quadrant"), size, alpha)
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) {
    stop("standardize must be TRUE or FALSE.")
  }
  observed <- sum(data$residual_centered * data$lag_scaled) /
    sum(data$residual_centered^2)
  data$residual_plot <- data$residual_centered
  data$lag_plot <- data$lag_scaled
  if (standardize) {
    residual_scale <- sqrt(mean(data$residual_centered^2))
    data$residual_plot <- data$residual_plot / residual_scale
    data$lag_plot <- data$lag_plot / residual_scale
  }
  colors <- .moranQuadrantColors()

  ggplot2::ggplot(data, ggplot2::aes(
    x = !!as.name("residual_plot"), y = !!as.name("lag_plot")
  )) +
    ggplot2::geom_hline(yintercept = 0, color = "grey80", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
    ggplot2::geom_abline(intercept = 0, slope = observed, color = "grey40") +
    ggplot2::geom_point(
      ggplot2::aes(fill = !!as.name("quadrant")),
      shape = 21, color = "grey30", stroke = 0.15, size = size, alpha = alpha
    ) +
    ggplot2::scale_fill_manual(
      values = colors, limits = names(colors), breaks = names(colors),
      drop = FALSE, na.value = "grey70", 
      guide = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1), 
        theme = ggplot2::theme(legend.margin = ggplot2::margin(0, 0, 0, 0)))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = if (standardize) "Standardized residual" else "Centered residual",
      y = "Rescaled spatial lag",
      fill = "Residual-lag quadrant"
    )
}


.prepareMoranPlotData <- function(data, columns, size, alpha) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install the optional package 'ggplot2' to plot Moran diagnostics.")
  }
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("data must be a non-empty data frame from getMoranLocal().")
  }
  required <- c("cell_index", "patch", "gene", "n_cells", "status",
                "residual_centered", columns)
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Missing columns in data: ", paste(missing, collapse = ", "), ".")
  }
  if (anyNA(data$patch) || anyNA(data$gene) ||
      length(unique(data$patch)) != 1L || length(unique(data$gene)) != 1L) {
    stop("data must contain exactly one gene and one assigned patch.")
  }
  if (anyNA(data$status) || any(data$status != "ok")) {
    stop("All rows must have status 'ok'; this gene-patch combination cannot be plotted.")
  }
  if (anyNA(data$n_cells) || any(data$n_cells != nrow(data)) ||
      nrow(data) < 3L || anyNA(data$cell_index) || anyDuplicated(data$cell_index)) {
    stop("data must contain every cell exactly once for the complete selected patch.")
  }
  if ("k" %in% names(data) && (anyNA(data$k) || length(unique(data$k)) != 1L)) {
    stop("data must come from a single spatial graph configuration.")
  }
  numeric_columns <- intersect(c("x", "y", "residual_centered", "lag_scaled", "local_i"), required)
  for (column in numeric_columns) {
    if (!is.numeric(data[[column]]) || any(!is.finite(data[[column]]))) {
      stop("Column '", column, "' must contain finite numeric values.")
    }
  }
  if (sum(data$residual_centered^2) == 0) {
    stop("Centered residuals must have positive variation.")
  }
  if ("quadrant" %in% columns &&
      any(!is.na(data$quadrant) & !data$quadrant %in% names(.moranQuadrantColors()))) {
    stop("quadrant contains labels not produced by getMoranLocal().")
  }
  if (!is.numeric(size) || length(size) != 1L || !is.finite(size) || size < 0) {
    stop("size must be a non-negative finite number.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must be a finite number between zero and one.")
  }
  as.data.frame(data)
}


.moranMapBase <- function(data, background, polygons, background_size,
                           background_alpha, aspect) {
  if (!is.numeric(background_size) || length(background_size) != 1L ||
      !is.finite(background_size) || background_size < 0) {
    stop("background_size must be a non-negative finite number.")
  }
  if (!is.numeric(background_alpha) || length(background_alpha) != 1L ||
      !is.finite(background_alpha) || background_alpha < 0 || background_alpha > 1) {
    stop("background_alpha must be a finite number between zero and one.")
  }
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = !!as.name("x"), y = !!as.name("y")))
  if (!is.null(background)) {
    background <- .prepareMoranMapCoordinates(background, "background")
    plot <- plot + ggplot2::geom_point(
      data = background, color = "grey85",
      size = background_size, alpha = background_alpha
    )
  }
  if (!is.null(polygons)) {
    polygons <- .prepareMoranMapCoordinates(polygons, "polygons")
    if (!"patch" %in% names(polygons)) stop("polygons must contain a 'patch' column.")
    polygons <- polygons[!is.na(polygons$patch) &
                           as.character(polygons$patch) == as.character(data$patch[1]), , drop = FALSE]
    if (!nrow(polygons)) stop("polygons contains no vertices for the selected patch.")
    plot <- plot + ggplot2::geom_polygon(
      data = polygons, ggplot2::aes(group = !!as.name("patch")),
      fill = NA, color = "grey30", linewidth = 0.2
    )
  }
  coordinates <- if (aspect == "fixed") {
    ggplot2::coord_fixed(xlim = range(data$x), ylim = range(data$y))
  } else {
    ggplot2::coord_cartesian(xlim = range(data$x), ylim = range(data$y))
  }
  plot + coordinates + ggplot2::theme_void()
}


.prepareMoranMapCoordinates <- function(data, label) {
  if (!is.data.frame(data) || !all(c("x", "y") %in% names(data)) ||
      !is.numeric(data$x) || !is.numeric(data$y) ||
      any(!is.finite(data$x)) || any(!is.finite(data$y))) {
    stop(label, " must be a data frame with finite numeric columns 'x' and 'y'.")
  }
  as.data.frame(data)
}


.moranQuadrantColors <- function() {
  c("High-High" = "#D73027", "Low-Low" = "#4575B4",
    "High-Low" = "#FDAE61", "Low-High" = "#66C2A5")
}

