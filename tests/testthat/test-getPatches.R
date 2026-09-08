test_that(".buildContiguityGraph creates Moran-compatible weights", {
  xy <- matrix(
    c(0, 0,
      1, 0,
      2, 0,
      3, 0,
      4, 0),
    ncol = 2,
    byrow = TRUE
  )

  W <- SpaceMosaic:::.buildContiguityGraph(xy, k = 2)

  expect_s4_class(W, "sparseMatrix")
  expect_true(isSymmetric(as.matrix(W)))
  expect_true(all(W@x == 1))
  expect_true(all(Matrix::diag(W) == 0))
  expect_true(all(Matrix::rowSums(W) >= 2))
})

test_that(".buildContiguityGraph adapts k to small patches", {
  xy <- matrix(
    c(0, 0,
      1, 0,
      0, 1),
    ncol = 2,
    byrow = TRUE
  )

  # k is reduced from 10 to n - 1, yielding the complete three-node graph.
  expect_warning(
    W <- SpaceMosaic:::.buildContiguityGraph(xy, k = 10),
    "k = 10 exceeds the maximum allowed for 3 observations; using k = 2.",
    fixed = TRUE
  )

  expect_equal(as.matrix(W), matrix(1, 3, 3) - diag(3))
})

test_that(".buildContiguityGraph validates coordinates and k", {
  valid_xy <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)

  expect_error(
    SpaceMosaic:::.buildContiguityGraph(as.data.frame(valid_xy)),
    "numeric matrix",
    fixed = TRUE
  )
  expect_error(
    SpaceMosaic:::.buildContiguityGraph(matrix(1:9, ncol = 3)),
    "exactly two columns",
    fixed = TRUE
  )

  non_finite_xy <- valid_xy
  non_finite_xy[1, 1] <- NA_real_
  expect_error(
    SpaceMosaic:::.buildContiguityGraph(non_finite_xy),
    "finite coordinates",
    fixed = TRUE
  )
  expect_error(
    SpaceMosaic:::.buildContiguityGraph(valid_xy[1, , drop = FALSE]),
    "at least two observations",
    fixed = TRUE
  )
  expect_error(
    SpaceMosaic:::.buildContiguityGraph(valid_xy, k = 0),
    "positive integer",
    fixed = TRUE
  )
  expect_error(
    SpaceMosaic:::.buildContiguityGraph(valid_xy, k = 1.5),
    "positive integer",
    fixed = TRUE
  )
})
