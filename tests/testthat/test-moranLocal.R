test_that("local Moran matches a hand-calculated chain graph", {
  # Unequal distances avoid ambiguous nearest-neighbor ties.
  xy <- cbind(x = c(0, 1, 3, 7), y = rep(0, 4))
  r <- c(-3, -1, 1, 3)
  result <- getMoranLocal(r + 10, xy, k = 1, verbose = FALSE)

  expect_equal(result$residual, r + 10)
  expect_equal(result$residual_centered, r)
  expect_equal(result$lag_centered, c(-1, -2, 2, 1))
  expect_equal(result$local_i, c(0.6, 0.4, 0.4, 0.6))
  expect_equal(result$local_i_scaled, c(6, 4, 4, 6) / 15)
  expect_equal(result$lag_scaled, (4 / 6) * result$lag_centered)
  expect_equal(with(result, sum(residual_centered * lag_scaled) / sum(residual_centered^2)), 1 / 3)
  expect_equal(sum(result$local_i) / 6, 1 / 3)
  expect_identical(result$quadrant, c("Low-Low", "Low-Low", "High-High", "High-High"))
  expect_true(all(result$status == "ok"))
  expect_false(any(result$complete_graph))
})

test_that("local statistics and scatterplot slopes reproduce global Moran", {
  set.seed(72)
  xy <- matrix(runif(60), ncol = 2)
  residuals <- matrix(rnorm(90), nrow = 30,
                      dimnames = list(paste0("cell_", 1:30), c("a", "b", "c")))
  patch <- rep(c("B", "A", NA_character_), times = c(14, 15, 1))
  for (k in c(1, 3, 8)) {
    local <- getMoranLocal(residuals, xy, patch, k = k, verbose = FALSE)
    global <- getMoranGlobal(residuals, xy, patch, k = k, verbose = FALSE)
    for (i in seq_len(nrow(global))) {
      d <- local[!is.na(local$patch) & local$patch == global$patch[i] &
                   local$gene == global$gene[i], ]
      W <- SpaceMosaic:::.buildContiguityGraph(xy[d$cell_index, , drop = FALSE], k = k)
      expect_equal(sum(d$local_i) / sum(W), global$observed[i])
      expect_equal(mean(d$local_i_scaled), global$observed[i])
      # Fitted slope in the displayed scatterplot has the same value.
      expect_equal(sum(d$residual_centered * d$lag_scaled) / sum(d$residual_centered^2),
                   global$observed[i])
    }
    expect_equal(nrow(local), length(residuals))
    expect_identical(local$cell_index, rep(1:30, 3))
    expect_identical(local$cell_id, rep(rownames(residuals), 3))
    expect_equal(local$residual, as.vector(residuals))
    expect_true(all(local$status[is.na(local$patch)] == "unassigned cell"))
    expect_true(all(is.na(local$local_i[is.na(local$patch)])))
  }
})

test_that("local Moran is scale invariant while centered residuals and lags retain units", {
  set.seed(23)
  xy <- matrix(runif(80), ncol = 2)
  r <- rnorm(40)
  original <- getMoranLocal(r, xy, k = 5, verbose = FALSE)
  transformed <- getMoranLocal(10 + 4 * r, xy, k = 5, verbose = FALSE)

  expect_equal(transformed$residual_centered, 4 * original$residual_centered)
  expect_equal(transformed$lag_centered, 4 * original$lag_centered)
  expect_equal(transformed$lag_scaled, 4 * original$lag_scaled)
  expect_equal(transformed$local_i, original$local_i)
  expect_equal(transformed$local_i_scaled, original$local_i_scaled)
  expect_identical(transformed$quadrant, original$quadrant)

  # Equivalent to the former standardized-residual formulation.
  z <- r - mean(r)
  z_standardized <- z / sqrt(mean(z^2))
  W <- SpaceMosaic:::.buildContiguityGraph(xy, k = 5)
  expect_equal(original$local_i, z_standardized * as.numeric(W %*% z_standardized))
})

test_that("local Moran retains invalid combinations without changing the graph", {
  xy <- cbind(x = 1:8, y = rep(0, 8))
  r <- cbind(valid = c(-2, 0, 2, 1, 2, 3, 4, 5), constant = rep(0, 8),
             invalid = c(1, NA, 3, 4, 5, 6, 7, 8))
  patch <- c("A", "A", "A", "B", "B", "C", NA, NA)
  result <- getMoranLocal(r, xy, patch, k = 1, verbose = FALSE)

  expect_equal(nrow(result), 24)
  a <- result[!is.na(result$patch) & result$patch == "A", ]
  expect_true(all(a$status[a$gene == "valid"] == "ok"))
  expect_true(all(a$status[a$gene == "constant"] == "constant residuals"))
  expect_true(all(a$status[a$gene == "invalid"] == "non-finite residuals"))
  expect_true(all(result$status[result$cell_index %in% 4:6] == "fewer than three cells"))
  invalid <- result$status != "ok"
  for (field in c("residual_centered", "lag_centered", "lag_scaled", "local_i", "local_i_scaled")) {
    expect_true(all(is.na(result[[field]][invalid])))
  }
  expect_true(all(is.na(result$quadrant[invalid])))
})

test_that("quadrants distinguish discordance and leave zero axes unclassified", {
  xy <- cbind(x = c(0, 1, 3, 7), y = rep(0, 4))
  result <- getMoranLocal(c(-3, 1, -1, 3), xy, k = 1, verbose = FALSE)
  expect_identical(result$quadrant, c("Low-High", "High-Low", "Low-High", "High-Low"))
  expect_true(all(result$local_i < 0))

  result <- getMoranLocal(c(-1, 0, 1), xy[1:3, ], k = 1, verbose = FALSE)
  expect_true(all(is.na(result$quadrant)))
  expect_equal(result$local_i, rep(0, 3))
  expect_equal(result$local_i_scaled, rep(0, 3))
})

test_that("complete graphs are flagged without discarding valid calculations", {
  xy <- cbind(x = c(0, 1, 4), y = rep(0, 3))
  expect_warning(
    result <- getMoranLocal(c(-1, 0, 1), xy, k = 10, verbose = FALSE),
    "using k = 2", fixed = TRUE
  )
  expect_true(all(result$complete_graph))
  expect_equal(mean(result$local_i_scaled), -1 / 2)
  expect_true(all(result$status == "ok"))
  expect_equal(with(result, sum(residual_centered * lag_scaled) / sum(residual_centered^2)), -1 / 2)
  expect_identical(result$effective_k, rep(2L, 3))

  result <- getMoranLocal(c(-1, 0, 1), xy, k = 2, verbose = FALSE)
  expect_true(all(result$complete_graph))
  expect_equal(mean(result$local_i_scaled), -1 / 2)
  expect_equal(with(result, sum(residual_centered * lag_scaled) / sum(residual_centered^2)), -1 / 2)
})

test_that("local Moran preserves identifiers and shared input validation", {
  xy <- cbind(x = c(0, 1, 3, 7), y = rep(0, 4))
  rownames(xy) <- letters[1:4]
  r <- c(-3, -1, 1, 3)
  expect_identical(getMoranLocal(r, xy, k = 1, verbose = FALSE)$cell_id, letters[1:4])
  names(r) <- LETTERS[1:4]
  expect_identical(getMoranLocal(r, xy, k = 1, verbose = FALSE)$cell_id, LETTERS[1:4])
  expect_error(getMoranLocal(r[-1], xy), "nrow(residuals) must equal nrow(xy)", fixed = TRUE)
  expect_error(getMoranLocal(r, xy, k = 0), "positive integer", fixed = TRUE)
  expect_error(getMoranLocal(r, xy, verbose = NA), "verbose must be TRUE or FALSE", fixed = TRUE)
  expect_error(getMoranLocal(r, xy, patch = rep(NA, 4)), "non-missing patch ID", fixed = TRUE)
  expect_error(getMoranLocal(matrix(1:8, 4, dimnames = list(NULL, c("g", "g"))), xy),
               "uniquely named column", fixed = TRUE)
  expect_silent(getMoranLocal(r, xy, k = 1, verbose = FALSE))
  expect_message(getMoranLocal(r, xy, k = 1, verbose = TRUE), "getMoranLocal")
})

test_that("local Moran SPE wrapper matches aligned matrix inputs", {
  xy <- cbind(x = c(0, 1, 3, 7), y = rep(0, 4))
  r <- matrix(c(-3, -1, 1, 3), nrow = 1,
               dimnames = list("Fos", paste0("cell_", 1:4)))
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(custom = r), spatialCoords = xy
  )
  spe$region <- "A"
  expected <- getMoranLocal(t(r), xy, rep("A", 4), k = 1, verbose = FALSE)
  result <- getMoranLocal.spe(spe, assay_name = "custom", patch_column = "region",
                             k = 1, verbose = FALSE)
  expect_identical(result, expected)
  expect_error(getMoranLocal.spe(r), "SpatialExperiment object", fixed = TRUE)
  expect_error(getMoranLocal.spe(spe), "not found in assay names", fixed = TRUE)
  expect_error(getMoranLocal.spe(spe, assay_name = "custom"), "not found in colData", fixed = TRUE)
})
