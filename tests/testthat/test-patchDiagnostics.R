test_that("getPatchDiagnostics calculates core diagnostics and assignment summary", {
  xy <- cbind(x = 0:4, y = rep(0, 5))
  rownames(xy) <- paste0("cell", seq_len(nrow(xy)))
  X <- c(1, 3, 9, 2, 6)
  patch <- c("A", "A", NA, "B", "B")
  names(patch) <- rownames(xy)

  membership_log <- cbind(
    iter1 = c("A", "A", NA, "A", "B"),
    iter2 = patch
  )
  rownames(membership_log) <- rownames(xy)
  result <- list(patch = patch, membership_log = membership_log)

  diagnostics <- getPatchDiagnostics(
    xy = xy,
    X = X,
    patches = result,
    k = 1
  )

  expect_named(
    diagnostics,
    c("patch_diagnostics", "assignment_summary", "connectivity_curve")
  )
  expect_equal(diagnostics$patch_diagnostics$n_cells, c(2, 2))
  expect_equal(diagnostics$patch_diagnostics$x_sd, c(sqrt(2), sqrt(8)))
  expect_equal(
    diagnostics$patch_diagnostics$membership_stability,
    c(1, 0.5)
  )

  expect_equal(diagnostics$assignment_summary$n_analyzed_cells, 5)
  expect_equal(diagnostics$assignment_summary$n_assigned_cells, 4)
  expect_equal(diagnostics$assignment_summary$fraction_assigned, 0.8)
  expect_equal(diagnostics$assignment_summary$n_unassigned_cells, 1)
  expect_equal(diagnostics$assignment_summary$fraction_unassigned, 0.2)
  expect_equal(diagnostics$assignment_summary$n_nonempty_patches, 2)
  expect_equal(diagnostics$assignment_summary$connectivity_k, 1)
  expect_equal(diagnostics$assignment_summary$strict_connectivity_k, 1)
})


test_that("result supplies final patches and all logged iterations", {
  xy <- cbind(x = 1:4, y = 0)
  final_patch <- c("A", "A", "B", "B")
  result <- list(
    patch = final_patch,
    membership_log = cbind(
      iter1 = c("A", "B", "B", "B"),
      iter2 = c("A", "A", "B", "A"),
      iter3 = final_patch
    )
  )

  diagnostics <- getPatchDiagnostics(xy, 1:4, patches = result, k = 1)

  expect_equal(
    diagnostics$patch_diagnostics$membership_stability,
    c(1, 0.5)
  )
})


test_that("a final patch vector works without iteration logs", {
  xy <- cbind(x = 1:3, y = 0)
  patch <- c("A", "A", "B")

  diagnostics <- getPatchDiagnostics(
    xy,
    1:3,
    patches = patch,
    k = 1
  )

  expect_true(all(is.na(
    diagnostics$patch_diagnostics$membership_stability
  )))
})


test_that("multivariable X produces one standard deviation per variable", {
  xy <- cbind(x = 1:6, y = 0)
  patch <- rep(c("A", "B"), each = 3)
  X <- cbind(
    exposure = c(1, 2, 3, 10, 12, 14),
    depth = c(4, 4, 4, 2, 3, 4)
  )

  diagnostics <- getPatchDiagnostics(xy, X, patches = patch, k = 1)
  patch_diagnostics <- diagnostics$patch_diagnostics

  expect_true(all(c(
    "x_sd_exposure", "x_sd_depth", "x_n_exposure", "x_n_depth"
  ) %in% names(patch_diagnostics)))
  expect_equal(patch_diagnostics$x_sd_exposure, c(1, 2))
  expect_equal(patch_diagnostics$x_sd_depth, c(0, 1))
  expect_true(all(is.na(
    patch_diagnostics$membership_stability
  )))
})


test_that("connectivity uses the global graph and detects separated fragments", {
  xy <- matrix(
    c(0, 0,
      0, 1,
      10, 0,
      10, 1),
    ncol = 2,
    byrow = TRUE
  )

  diagnostics <- getPatchDiagnostics(
    xy = xy,
    X = 1:4,
    patches = rep("A", 4),
    k = 1
  )

  expect_equal(
    diagnostics$patch_diagnostics$strict_component_fraction,
    0.5
  )
})


test_that("multi-k connectivity distinguishes strict and permissive graphs", {
  xy <- cbind(x = c(0, 1, 10, 11), y = 0)

  result <- getPatchDiagnostics(
    xy = xy,
    X = 1:4,
    patches = rep("A", 4),
    k = 3,
    strict_k = 1
  )
  diagnostics <- result$patch_diagnostics

  expect_equal(diagnostics$strict_component_fraction, 0.5)
  expect_equal(diagnostics$min_connectivity_k, 2L)
  expect_equal(result$connectivity_curve$component_fraction, c(0.5, 1, 1))
})


test_that("single-cell patches are connected without graph edges", {
  result <- getPatchDiagnostics(
    xy = matrix(c(0, 0), nrow = 1),
    X = 1,
    patches = "A"
  )
  diagnostics <- result$patch_diagnostics

  expect_equal(diagnostics$strict_component_fraction, 1)
  expect_equal(diagnostics$min_connectivity_k, 0L)
  expect_equal(nrow(result$connectivity_curve), 0L)

  result_with_other_cells <- getPatchDiagnostics(
    xy = cbind(x = 0:2, y = 0),
    X = 1:3,
    patches = c("A", "B", "B"),
    k = 1
  )
  singleton <- result_with_other_cells$patch_diagnostics$patch == "A"
  expect_equal(
    result_with_other_cells$patch_diagnostics$min_connectivity_k[singleton],
    0L
  )
})


test_that("unavailable X SD retains finite-value counts", {
  xy <- cbind(x = 1:3, y = 0)
  X <- c(1, Inf, 3)
  patch <- c("A", "B", "B")

  diagnostics <- getPatchDiagnostics(xy, X, patches = patch, k = 1)
  patch_diagnostics <- diagnostics$patch_diagnostics

  expect_true(is.na(patch_diagnostics$x_sd[patch_diagnostics$patch == "A"]))
  expect_equal(
    patch_diagnostics$strict_component_fraction[
      patch_diagnostics$patch == "A"
    ],
    1
  )
  expect_equal(patch_diagnostics$x_n, c(1L, 1L))
  expect_false("diagnostic_status" %in% names(patch_diagnostics))
})


test_that("ordinary missing X values are counted and excluded from SD", {
  xy <- cbind(x = 1:4, y = 0)
  X <- c(1, NA, 3, 5)

  diagnostics <- getPatchDiagnostics(
    xy,
    X,
    patches = rep("A", 4),
    k = 1
  )$patch_diagnostics

  expect_equal(diagnostics$x_n, 3)
  expect_equal(diagnostics$x_sd, 2)
})


test_that("all-unassigned input returns an empty patch table", {
  xy <- cbind(x = 1:3, y = 0)
  diagnostics <- getPatchDiagnostics(
    xy = xy,
    X = 1:3,
    patches = rep(NA_character_, 3),
    k = 1
  )

  expect_equal(nrow(diagnostics$patch_diagnostics), 0)
  expect_equal(diagnostics$assignment_summary$n_assigned_cells, 0)
  expect_equal(diagnostics$assignment_summary$fraction_unassigned, 1)
  expect_equal(diagnostics$assignment_summary$n_nonempty_patches, 0)
  expect_equal(nrow(diagnostics$connectivity_curve), 0)
})


test_that("diagnostic table is compatible with getPatchPolys", {
  xy <- matrix(
    c(0, 0, 1, 0, 0, 1,
      3, 3, 4, 3, 3, 4),
    ncol = 2,
    byrow = TRUE
  )
  patch <- rep(c("A", "B"), each = 3)
  diagnostics <- getPatchDiagnostics(xy, 1:6, patches = patch, k = 1)

  polygons <- getPatchPolys(
    xy,
    patch,
    patch_data = diagnostics$patch_diagnostics
  )

  expect_true(all(c(
    "n_cells",
    "x_sd",
    "strict_component_fraction",
    "min_connectivity_k",
    "membership_stability"
  ) %in% names(polygons)))
})


test_that("alignment and input errors are informative", {
  xy <- cbind(x = 1:3, y = 0)
  rownames(xy) <- c("a", "b", "c")

  expect_error(
    getPatchDiagnostics(xy, 1:3, patches = NULL, k = 1),
    "patch vector or the output of getPatches",
    fixed = TRUE
  )
  expect_error(
    getPatchDiagnostics(xy, factor(c("low", "mid", "high")),
                        patches = c("A", "A", "B"), k = 1),
    "numeric vector or matrix",
    fixed = TRUE
  )
  expect_error(
    getPatchDiagnostics(xy, 1:2, patches = c("A", "B", "B"), k = 1),
    "nrow(X) must equal nrow(xy)",
    fixed = TRUE
  )
  expect_error(
    getPatchDiagnostics(
      xy,
      1:3,
      patches = c("A", "A", "B"),
      k = 2,
      strict_k = 3
    ),
    "strict_k must be a positive integer no greater than k",
    fixed = TRUE
  )

  misordered_X <- c(b = 1, a = 2, c = 3)
  expect_error(
    getPatchDiagnostics(
      xy,
      misordered_X,
      patches = c("A", "A", "B"),
      k = 1
    ),
    "rownames(X) must match rownames(xy)",
    fixed = TRUE
  )

  named_patch <- c(a = "A", c = "B", b = "B")
  expect_error(
    getPatchDiagnostics(xy, 1:3, patches = named_patch, k = 1),
    "names(patches) must match rownames(xy)",
    fixed = TRUE
  )
})


test_that("numeric-looking patch identifiers use numeric order", {
  xy <- cbind(x = 1:6, y = 0)
  patch <- c("10", "10", "2", "2", "1", "1")

  diagnostics <- getPatchDiagnostics(xy, 1:6, patches = patch, k = 1)

  expect_equal(diagnostics$patch_diagnostics$patch, c("1", "2", "10"))
})
