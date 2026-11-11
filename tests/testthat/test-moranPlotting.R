skip_if_not_installed("ggplot2")

make_moran_plot_data <- function() {
  xy <- cbind(x = c(0, 1, 3, 7), y = c(0, 1, 0, 1))
  getMoranLocal(c(-3, -1, 1, 3), xy, k = 1, verbose = FALSE)
}

test_that("residual maps rescale only colors and remain extensible ggplots", {
  d <- make_moran_plot_data()
  original <- d
  plain <- plotMoranResiduals(d)
  scaled <- plotMoranResiduals(d, standardize = TRUE, size = 3, alpha = 0.5)
  expect_s3_class(plain, "ggplot")
  expect_equal(plain$data$residual_plot, d$residual_centered)
  expect_equal(scaled$data$residual_plot, d$residual_centered / sqrt(5))
  expect_equal(scaled$data$local_i, plain$data$local_i)
  expect_identical(d, original)

  built <- ggplot2::ggplot_build(scaled)
  expect_equal(built$data[[1]]$x, d$x)
  expect_equal(built$data[[1]]$y, d$y)
  expect_equal(built$data[[1]]$size, rep(3, 4))
  expect_equal(built$data[[1]]$alpha, rep(0.5, 4))

  customized <- suppressMessages(
    plain + ggplot2::scale_fill_gradient2(limits = c(-6, 6))
  ) + ggplot2::labs(title = "Custom title") +
    ggplot2::geom_text(ggplot2::aes(label = cell_id))
  expect_equal(customized$labels$title, "Custom title")
  expect_silent(ggplot2::ggplot_build(customized))
})

test_that("local maps preserve statistic magnitudes and use a symmetric diverging scale", {
  # Map the original local values, including negative association.
  xy <- cbind(x = 0:3, y = c(0, 1, 0, 1))
  d <- getMoranLocal(c(-3, 1, -1, 3), xy, k = 1, verbose = FALSE)
  original <- d
  p <- plotMoranLocal(d, size = 3, alpha = 0.5)
  expect_s3_class(p, "ggplot")
  expect_identical(p$data[names(d)], d)
  expect_equal(p$data$local_i_plot, d$local_i)
  expect_identical(d, original)
  expect_true(any(d$local_i < 0))
  scale <- p$scales$get_scales("fill")
  limit <- max(abs(d$local_i))
  expect_equal(scale$limits, c(-limit, limit))
  expect_equal(scale$map(c(-limit, 0, limit)), c("#2166AC", "#FFFFFF", "#B2182B"))
  built <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(built$fill, scale$map(d$local_i))
  expect_equal(built$x, d$x)
  expect_equal(built$y, d$y)
  expect_equal(built$size, rep(3, nrow(d)))
  expect_equal(built$alpha, rep(0.5, nrow(d)))
  expect_equal(p$labels$fill, "Local Moran's I")

  clustered <- getMoranLocal(c(-3, -3, -3, 3, 3, 3),
                            cbind(x = 0:5, y = 0), k = 1, verbose = FALSE)
  expect_gt(max(clustered$local_i), 1)
  clustered_plot <- plotMoranLocal(clustered)
  expect_equal(clustered_plot$scales$get_scales("fill")$limits, c(-2, 2))
  expect_false(anyNA(ggplot2::ggplot_build(clustered_plot)$data[[1]]$fill))

  # A nonconstant residual vector can have zero local statistics everywhere.
  zero <- getMoranLocal(c(-1, 0, 1), xy[1:3, ], k = 1, verbose = FALSE)
  expect_equal(zero$local_i, rep(0, 3))
  expect_equal(ggplot2::ggplot_build(plotMoranLocal(zero))$data[[1]]$fill,
               rep("#FFFFFF", 3))
})

test_that("standardized local maps average to global Moran and match scatter products", {
  set.seed(44)
  xy <- matrix(runif(100), ncol = 2)
  r <- rnorm(50)
  d <- getMoranLocal(r, xy, k = 5, verbose = FALSE)
  original <- d
  p <- plotMoranLocal(d, standardize = TRUE)
  global <- getMoranGlobal(r, xy, k = 5, verbose = FALSE)
  expect_equal(mean(p$data$local_i_plot), global$observed)
  scatter <- plotMoranScatter(d, standardize = TRUE)$data
  expect_equal(p$data$local_i_plot, scatter$residual_plot * scatter$lag_plot)
  expect_identical(p$data[names(d)], d)
  expect_identical(d, original)
  expect_equal(p$labels$fill, "Scaled local Moran's I")
  limit <- max(abs(p$data$local_i_plot))
  expect_equal(p$scales$get_scales("fill")$limits, c(-limit, limit))
  expect_equal(ggplot2::ggplot_build(p)$data[[1]]$fill,
               p$scales$get_scales("fill")$map(p$data$local_i_plot))

  transformed <- getMoranLocal(10 + 4 * r, xy, k = 5, verbose = FALSE)
  expect_equal(plotMoranLocal(transformed, standardize = TRUE)$data$local_i_plot,
               p$data$local_i_plot)
  zero <- getMoranLocal(c(-1, 0, 1), cbind(x = 0:2, y = 0), k = 1, verbose = FALSE)
  expect_equal(ggplot2::ggplot_build(plotMoranLocal(zero, standardize = TRUE))$data[[1]]$fill,
               rep("#FFFFFF", 3))
})

test_that("local maps require finite numeric local statistics", {
  d <- make_moran_plot_data()
  expect_error(plotMoranLocal(d[, setdiff(names(d), "local_i")]),
               "Missing columns in data: local_i", fixed = TRUE)
  for (invalid in list(NA_real_, NaN, Inf, -Inf, "1")) {
    d$local_i <- invalid
    expect_error(plotMoranLocal(d), "Column 'local_i' must contain finite numeric values",
                 fixed = TRUE)
  }
  d <- make_moran_plot_data()
  for (invalid in list(NA, NULL, c(TRUE, FALSE), 1, "TRUE")) {
    expect_error(plotMoranLocal(d, standardize = invalid),
                 "standardize must be TRUE or FALSE", fixed = TRUE)
  }
  without_lag <- d[, setdiff(names(d), "lag_scaled")]
  expect_s3_class(plotMoranLocal(without_lag), "ggplot")
  expect_error(plotMoranLocal(without_lag, standardize = TRUE),
               "Missing columns in data: lag_scaled", fixed = TRUE)
  d$lag_scaled[1] <- Inf
  expect_error(plotMoranLocal(d, standardize = TRUE),
               "Column 'lag_scaled' must contain finite numeric values", fixed = TRUE)
})

test_that("maps place context and selected patch outlines behind the cells", {
  d <- make_moran_plot_data()
  background <- data.frame(x = c(-100, 100), y = c(-100, 100))
  polygons <- getPatchPolys(as.matrix(d[c("x", "y")]), d$patch)
  other <- polygons
  other$patch <- "other"
  polygons <- rbind(polygons, other)
  for (fun in list(plotMoranResiduals, plotMoranLocal, plotMoranQuadrants)) {
    p <- fun(d, background = background, polygons = polygons,
             size = 3, alpha = 0.9, background_size = 0.2, background_alpha = 0.1)
    built <- ggplot2::ggplot_build(p)
    expect_equal(built$data[[1]]$x, background$x)
    expect_equal(built$data[[1]]$size, rep(0.2, nrow(background)))
    expect_equal(built$data[[1]]$alpha, rep(0.1, nrow(background)))
    expect_equal(nrow(built$data[[2]]), sum(polygons$patch == "all"))
    expect_equal(built$data[[3]]$x, d$x)
    expect_equal(built$data[[3]]$size, rep(3, nrow(d)))
    expect_equal(built$data[[3]]$alpha, rep(0.9, nrow(d)))
    expect_equal(p$coordinates$limits$x, range(d$x))
    expect_equal(p$coordinates$limits$y, range(d$y))
    expect_equal(p$coordinates$ratio, 1)

    defaults <- ggplot2::ggplot_build(fun(d, background = background))
    expect_equal(defaults$data[[1]]$size, rep(0.5, nrow(background)))
    expect_equal(defaults$data[[1]]$alpha, rep(0.3, nrow(background)))
  }
})

test_that("Moran plots leave titles to the user and scatterplot proportions free", {
  d <- make_moran_plot_data()
  for (fun in list(plotMoranResiduals, plotMoranLocal, plotMoranQuadrants, plotMoranScatter)) {
    p <- fun(d)
    expect_null(p$labels$title)
    expect_null(p$labels$subtitle)
    custom <- p + ggplot2::labs(title = "Fos - patch 27", subtitle = "Selected patch")
    expect_identical(custom$labels$title, "Fos - patch 27")
    expect_identical(custom$labels$subtitle, "Selected patch")
  }
  expect_null(plotMoranScatter(d)$coordinates$ratio)
})

test_that("spatial maps allow fixed or free proportions with identical data and limits", {
  d <- make_moran_plot_data()
  background <- data.frame(x = c(-100, 100), y = c(-100, 100))
  polygons <- getPatchPolys(as.matrix(d[c("x", "y")]), d$patch)
  for (fun in list(plotMoranResiduals, plotMoranLocal, plotMoranQuadrants)) {
    fixed <- fun(d, background = background, polygons = polygons, aspect = "fixed")
    free <- fun(d, background = background, polygons = polygons, aspect = "free")
    expect_equal(fun(d)$coordinates$ratio, 1)
    expect_equal(fixed$coordinates$ratio, 1)
    expect_null(free$coordinates$ratio)
    expect_equal(free$coordinates$limits, fixed$coordinates$limits)
    expect_equal(ggplot2::ggplot_build(free)$data,
                 ggplot2::ggplot_build(fixed)$data)
    expect_error(fun(d, aspect = "invalid"), "arg.*should be one of")
  }
})

test_that("quadrant maps and scatterplots share colors including zero-axis cells", {
  xy <- cbind(x = c(0, 1, 3, 7), y = c(0, 1, 0, 1))
  d <- getMoranLocal(c(-3, 1, -1, 3), xy, k = 1, verbose = FALSE)
  map <- ggplot2::ggplot_build(plotMoranQuadrants(d))
  scatter <- ggplot2::ggplot_build(plotMoranScatter(d))
  expect_equal(map$data[[1]]$fill, scatter$data[[4]]$fill)

  # All lag values are zero for the centered vector on a three-node chain.
  d <- getMoranLocal(c(-1, 0, 1), xy[1:3, ], k = 1, verbose = FALSE)
  map <- ggplot2::ggplot_build(plotMoranQuadrants(d))
  scatter <- ggplot2::ggplot_build(plotMoranScatter(d))
  expect_identical(map$data[[1]]$fill, rep("grey70", 3))
  expect_identical(scatter$data[[4]]$fill, rep("grey70", 3))
})

test_that("scatterplot line reproduces the independently calculated global statistic", {
  set.seed(44)
  xy <- matrix(runif(100), ncol = 2)
  r <- rnorm(50)
  d <- getMoranLocal(r, xy, k = 5, verbose = FALSE)
  global <- getMoranGlobal(r, xy, k = 5, verbose = FALSE)
  p <- plotMoranScatter(d)
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$data[[3]]$slope, global$observed)
  expect_equal(built$data[[3]]$intercept, 0)
  expect_equal(built$data[[4]]$x, d$residual_centered)
  expect_equal(built$data[[4]]$y, d$lag_scaled)
})

test_that("standardized scatterplots remove residual scale and preserve Moran diagnostics", {
  d <- make_moran_plot_data()
  original <- d
  plain <- ggplot2::ggplot_build(plotMoranScatter(d))
  p <- plotMoranScatter(d, standardize = TRUE)
  scaled <- ggplot2::ggplot_build(p)
  expect_equal(scaled$data[[4]]$x, d$residual_centered / sqrt(5))
  expect_equal(scaled$data[[4]]$y, d$lag_scaled / sqrt(5))
  expect_equal(mean(scaled$data[[4]]$x^2), 1)
  expect_equal(scaled$data[[3]]$slope, plain$data[[3]]$slope)
  expect_equal(sum(scaled$data[[4]]$x * scaled$data[[4]]$y) /
                 sum(scaled$data[[4]]$x^2), scaled$data[[3]]$slope)
  expect_equal(scaled$data[[4]]$fill, plain$data[[4]]$fill)
  expect_equal(p$data[names(d)], d)
  expect_identical(d, original)
  expect_match(p$labels$x, "Standardized")
  expect_match(p$labels$y, "residual SD units", fixed = TRUE)

  # Positive rescaling of the original residuals must leave both displayed
  # coordinates unchanged, not just the slope or the quadrant assignments.
  xy <- as.matrix(d[c("x", "y")])
  transformed <- getMoranLocal(10 + 4 * d$residual, xy, k = 1, verbose = FALSE)
  transformed_plot <- ggplot2::ggplot_build(plotMoranScatter(transformed, standardize = TRUE))
  expect_equal(transformed_plot$data[[4]]$x, scaled$data[[4]]$x)
  expect_equal(transformed_plot$data[[4]]$y, scaled$data[[4]]$y)

  for (invalid in list(NA, NULL, c(TRUE, FALSE), 1, "TRUE")) {
    expect_error(plotMoranScatter(d, standardize = invalid),
                 "standardize must be TRUE or FALSE", fixed = TRUE)
  }
})

test_that("Moran plots reject mixed, incomplete and uncomputable results", {
  d <- make_moran_plot_data()
  invalid <- d
  invalid$status <- "constant residuals"
  mixed <- d
  mixed$gene[1] <- "another_gene"
  duplicated <- d
  duplicated$cell_index[1] <- duplicated$cell_index[2]
  for (fun in list(plotMoranResiduals, plotMoranLocal, plotMoranQuadrants, plotMoranScatter)) {
    expect_error(fun(mixed), "exactly one gene", fixed = TRUE)
    expect_error(fun(d[-1, ]), "complete selected patch", fixed = TRUE)
    expect_error(fun(duplicated), "complete selected patch", fixed = TRUE)
    expect_error(fun(invalid), "status 'ok'", fixed = TRUE)
    expect_error(fun(d[0, ]), "non-empty data frame", fixed = TRUE)
    expect_error(fun(d, size = -1), "size must be", fixed = TRUE)
    expect_error(fun(d, alpha = 2), "alpha must be", fixed = TRUE)
  }
  expect_error(plotMoranResiduals(d, standardize = NA), "standardize must be", fixed = TRUE)
  expect_error(plotMoranScatter(d[, setdiff(names(d), "lag_scaled")]),
               "Missing columns", fixed = TRUE)
  expect_error(plotMoranQuadrants(d, background = data.frame(x = 1)),
               "finite numeric columns", fixed = TRUE)
  expect_error(plotMoranQuadrants(d, polygons = data.frame(x = 1, y = 1, patch = "other")),
               "no vertices", fixed = TRUE)
  for (fun in list(plotMoranResiduals, plotMoranLocal, plotMoranQuadrants)) {
    expect_error(fun(d, background_size = -1), "background_size must be", fixed = TRUE)
    expect_error(fun(d, background_alpha = NA), "background_alpha must be", fixed = TRUE)
    expect_error(fun(d, background_alpha = 2), "background_alpha must be", fixed = TRUE)
  }
})

test_that("local maps combine complete patches with independent standardization", {
  xy <- cbind(x = c(0, 1, 3, 7, 20, 21, 23, 27, 28),
              y = c(0, 1, 0, 1, 10, 11, 10, 11, 10))
  r <- c(-3, -1, 1, 3, -20, 5, -5, 8, 12)
  patches <- rep(c("a", "b"), c(4, 5))
  d <- getMoranLocal(r, xy, patch = patches, k = 1, verbose = FALSE)
  # Accept factors with unused levels and preserve interleaved input ordering.
  d$patch <- factor(d$patch, levels = c("b", "unused", "a"))
  d <- d[c(5, 1, 6, 2, 7, 3, 8, 4, 9), ]
  original <- d
  polygons <- getPatchPolys(xy, patches)
  other <- polygons[polygons$patch == "a", ]
  other$patch <- "other"
  polygons <- rbind(polygons, other)
  background <- data.frame(x = c(-100, 100), y = c(-100, 100))
  global <- getMoranGlobal(r, xy, patch = patches, k = 1, verbose = FALSE)

  for (standardize in c(FALSE, TRUE)) {
    p <- plotMoranLocal(d, standardize = standardize, polygons = polygons,
                        background = background)
    expect_identical(p$data[names(d)], d)
    for (patch in c("a", "b")) {
      rows <- which(d$patch == patch)
      separate <- plotMoranLocal(d[rows, ], standardize = standardize)
      expect_equal(p$data$local_i_plot[rows], separate$data$local_i_plot)
      if (standardize) {
        expect_equal(mean(p$data$local_i_plot[rows]),
                     global$observed[global$patch == patch])
      }
    }
    limit <- max(abs(p$data$local_i_plot))
    expect_equal(p$scales$get_scales("fill")$limits, c(-limit, limit))
    expect_equal(p$coordinates$limits$x, range(d$x))
    expect_equal(p$coordinates$limits$y, range(d$y))
    built <- ggplot2::ggplot_build(p)
    expect_equal(length(unique(built$data[[2]]$group)), 2L)
    expect_equal(nrow(built$data[[2]]), sum(polygons$patch %in% c("a", "b")))
    expect_equal(built$data[[3]]$x, d$x)
    expect_equal(built$data[[3]]$fill,
                 p$scales$get_scales("fill")$map(p$data$local_i_plot))
    free <- plotMoranLocal(d, standardize = standardize, aspect = "free")
    expect_null(free$coordinates$ratio)
    expect_equal(free$coordinates$limits, p$coordinates$limits)
  }
  expect_identical(d, original)
})

test_that("multi-patch maps validate every patch and retain single-gene scope", {
  xy <- cbind(x = c(0:3, 10:13), y = rep(c(0, 1), 4))
  d <- getMoranLocal(c(-3, -1, 1, 3, -30, -10, 10, 30), xy,
                     patch = rep(c("a", "b"), each = 4), k = 1, verbose = FALSE)
  for (fun in list(plotMoranResiduals, plotMoranLocal, plotMoranQuadrants)) {
    expect_error(fun(d[-8, ]), "complete selected patch", fixed = TRUE)
    invalid <- d
    invalid$n_cells[8] <- 8
    expect_error(fun(invalid), "complete selected patch", fixed = TRUE)
    invalid <- d
    invalid$cell_index[8] <- invalid$cell_index[7]
    expect_error(fun(invalid), "complete selected patch", fixed = TRUE)
    invalid <- d
    invalid$residual_centered[5:8] <- 0
    expect_error(fun(invalid), "positive variation", fixed = TRUE)
    invalid <- d
    invalid$status[8] <- "constant residuals"
    expect_error(fun(invalid), "status 'ok'", fixed = TRUE)
    invalid <- d
    invalid$gene[5:8] <- "other"
    expect_error(fun(invalid), "exactly one gene", fixed = TRUE)
    invalid <- d
    invalid$patch[8] <- NA_character_
    expect_error(fun(invalid), "assigned patches", fixed = TRUE)
  }
  expect_error(plotMoranScatter(d), "one assigned patch", fixed = TRUE)
})

test_that("residual and quadrant maps combine patches without changing patch diagnostics", {
  xy <- cbind(x = c(0, 1, 3, 7, 20, 21, 23, 27, 28),
              y = c(0, 1, 0, 1, 10, 11, 10, 11, 10))
  patches <- rep(c("a", "b"), c(4, 5))
  d <- getMoranLocal(c(-3, -1, 1, 3, -20, 5, -5, 8, 12), xy,
                     patch = patches, k = 1, verbose = FALSE)
  d$patch <- factor(d$patch, levels = c("b", "unused", "a"))
  d <- d[c(5, 1, 6, 2, 7, 3, 8, 4, 9), ]
  original <- d
  polygons <- getPatchPolys(xy, patches)
  for (standardize in c(FALSE, TRUE)) {
    p <- plotMoranResiduals(d, standardize = standardize, polygons = polygons)
    expect_identical(p$data[names(d)], d)
    for (patch in c("a", "b")) {
      rows <- which(d$patch == patch)
      separate <- plotMoranResiduals(d[rows, ], standardize = standardize)
      expect_equal(p$data$residual_plot[rows], separate$data$residual_plot)
      if (standardize) expect_equal(mean(p$data$residual_plot[rows]^2), 1)
    }
    limit <- max(abs(p$data$residual_plot))
    expect_equal(p$scales$get_scales("fill")$limits, c(-limit, limit))
    built <- ggplot2::ggplot_build(p)
    expect_equal(length(unique(built$data[[1]]$group)), 2L)
    expect_equal(built$data[[2]]$fill,
                 p$scales$get_scales("fill")$map(p$data$residual_plot))
  }
  p <- plotMoranQuadrants(d, polygons = polygons)
  expect_identical(p$data, d)
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$data[[1]]$group)), 2L)
  for (patch in c("a", "b")) {
    rows <- which(d$patch == patch)
    separate <- ggplot2::ggplot_build(plotMoranQuadrants(d[rows, ]))
    expect_equal(built$data[[2]]$fill[rows], separate$data[[1]]$fill)
  }
  for (fun in list(plotMoranResiduals, plotMoranQuadrants)) {
    fixed <- fun(d)
    free <- fun(d, aspect = "free")
    expect_equal(fixed$coordinates$limits$x, range(d$x))
    expect_equal(fixed$coordinates$limits$y, range(d$y))
    expect_equal(free$coordinates$limits, fixed$coordinates$limits)
    expect_null(free$coordinates$ratio)
  }
  expect_identical(d, original)
})
