test_that("moranTest builds the graph and reports its settings", {
  xy <- matrix(
    c(0, 0,
      1, 0,
      2, 0,
      3, 0,
      4, 0),
    ncol = 2,
    byrow = TRUE
  )
  residuals <- c(-2, -1, 0, 1, 2)

  set.seed(11)
  result <- moranTest(
    residuals, xy, k = 2, n_permutations = 19, alternative = "greater"
  )
  W <- SpaceMosaic:::.buildContiguityGraph(xy, k = 2)

  expect_equal(result$observed, SpaceMosaic:::.moranI(residuals, W))
  expect_identical(result$k, 2L)
  expect_identical(result$effective_k, 2L)
  expect_identical(result$n_permutations, 19L)
})

test_that("moranTest reports an adjusted k", {
  xy <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)

  set.seed(3)
  expect_warning(
    result <- moranTest(c(-1, 0, 1), xy, k = 10, n_permutations = 9),
    "using k = 2",
    fixed = TRUE
  )

  expect_identical(result$k, 10L)
  expect_identical(result$effective_k, 2L)
})

test_that("moranTest requires residuals and coordinates to be aligned", {
  xy <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)

  expect_error(
    moranTest(c(-1, 1), xy, k = 1, n_permutations = 9),
    "nrow(residuals) must equal nrow(xy)",
    fixed = TRUE
  )
})

test_that("patchDE residuals feed reproducibly into moranTest", {
  xy <- cbind(x = 0:5, y = rep(0, 6))
  df <- data.frame(treatment = rep(c(0, 1), 3))
  spatial_component <- c(-2, -1, -0.5, 0.5, 1, 2)
  y <- cbind(gene_a = 2 + 3 * df$treatment + spatial_component)
  patch <- rep("one", nrow(y))

  fit <- patchDE(
    y, df, patch,
    return_residuals = TRUE,
    verbose = FALSE
  )

  set.seed(99)
  first <- moranTest(
    fit$residuals, xy, patch = patch, k = 2, n_permutations = 19
  )
  set.seed(99)
  second <- moranTest(
    fit$residuals, xy, patch = patch, k = 2, n_permutations = 19
  )

  expect_identical(first, second)
  expect_equal(
    first$observed,
    SpaceMosaic:::.moranI(
      fit$residuals[, "gene_a"],
      SpaceMosaic:::.buildContiguityGraph(xy, k = 2)
    )
  )
})

test_that("moranTest handles multiple genes and patches", {
  xy <- cbind(x = rep(0:3, 2), y = rep(c(0, 10), each = 4))
  patch <- rep(c("one", "two"), each = 4)
  residuals <- cbind(
    gene_a = c(-2, -1, 1, 2, 2, 1, -1, -2),
    gene_b = c(1, 3, 2, 8, 4, 1, 7, 2)
  )

  set.seed(5)
  result <- moranTest(
    residuals, xy, patch = patch,
    k = 2, n_permutations = 19,
    p_adjust_method = "BH",
    adjustment_scope = "global"
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)
  expect_equal(result$patch, rep(c("one", "two"), each = 2))
  expect_equal(result$gene, rep(c("gene_a", "gene_b"), 2))
  expect_true(all(result$status == "ok"))
  expect_equal(result$p_adjusted, stats::p.adjust(result$p_value, "BH"))
})

test_that("moranTest reports non-testable genes without stopping a batch", {
  xy <- cbind(x = 0:3, y = rep(0, 4))
  residuals <- cbind(
    informative = c(-2, -1, 1, 2),
    constant = rep(0, 4),
    non_finite = c(1, 2, NA, 4)
  )

  set.seed(8)
  result <- moranTest(residuals, xy, k = 2, n_permutations = 9)

  expect_equal(
    result$status,
    c("ok", "constant residuals", "non-finite residuals")
  )
  expect_true(is.finite(result$p_value[result$gene == "informative"]))
  expect_true(all(is.na(result$p_value[result$gene != "informative"])))
})

test_that("moranTest can adjust p-values within patches", {
  p <- c(0.01, 0.04, 0.02, 0.5)
  patch <- c("one", "one", "two", "two")
  gene <- rep(c("a", "b"), 2)

  adjusted <- SpaceMosaic:::.adjustMoranPValues(
    p, patch, gene, method = "BH", scope = "patch"
  )

  expect_equal(
    adjusted,
    c(stats::p.adjust(p[1:2], "BH"), stats::p.adjust(p[3:4], "BH"))
  )
})

test_that(".moranI identifies perfectly grouped values", {
  # The only edges join observations with equal centered values, producing
  # perfect positive spatial autocorrelation.
  W <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4),
    j = c(2, 1, 4, 3),
    x = 1,
    dims = c(4, 4)
  )

  expect_equal(SpaceMosaic:::.moranI(c(0, 0, 1, 1), W), 1)
})

test_that(".moranI matches the complete-graph expectation", {
  # On a complete graph without self-edges, Moran's I is -1 / (n - 1) for
  # every non-constant vector.
  W_dense <- matrix(1, nrow = 4, ncol = 4)
  diag(W_dense) <- 0
  W_sparse <- Matrix::Matrix(W_dense, sparse = TRUE)

  expected <- -1 / 3
  x <- c(1, 4, 2, 8)

  expect_equal(SpaceMosaic:::.moranI(x, W_dense), expected)
  expect_equal(SpaceMosaic:::.moranI(x, W_sparse), expected)
})

test_that(".moranI is invariant to shifts and non-zero rescaling", {
  W <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4),
    j = c(2, 1, 3, 2, 4, 3),
    x = 1,
    dims = c(4, 4)
  )
  x <- c(1, 2, 5, 9)
  observed <- SpaceMosaic:::.moranI(x, W)

  expect_equal(SpaceMosaic:::.moranI(x + 100, W), observed)
  expect_equal(SpaceMosaic:::.moranI(-3 * x, W), observed)
})

test_that(".moranI validates its inputs", {
  valid_W <- matrix(
    c(0, 1, 0,
      1, 0, 1,
      0, 1, 0),
    nrow = 3,
    byrow = TRUE
  )

  expect_error(
    SpaceMosaic:::.moranI(c(1, 1, 1), valid_W),
    "constant vector",
    fixed = TRUE
  )
  expect_error(
    SpaceMosaic:::.moranI(c(1, 2, 3), matrix(0, 3, 3)),
    "no edges",
    fixed = TRUE
  )
  expect_error(
    SpaceMosaic:::.moranI(c(1, 2, 3), matrix(0, 2, 2)),
    "matching length",
    fixed = TRUE
  )

  negative_W <- valid_W
  negative_W[1, 2] <- negative_W[2, 1] <- -1
  expect_error(
    SpaceMosaic:::.moranI(c(1, 2, 3), negative_W),
    "non-negative",
    fixed = TRUE
  )

  diagonal_W <- valid_W
  diag(diagonal_W) <- 1
  expect_error(
    SpaceMosaic:::.moranI(c(1, 2, 3), diagonal_W),
    "zero diagonal",
    fixed = TRUE
  )

  asymmetric_W <- valid_W
  asymmetric_W[1, 2] <- 0
  expect_error(
    SpaceMosaic:::.moranI(c(1, 2, 3), asymmetric_W),
    "symmetric",
    fixed = TRUE
  )
})

test_that(".moranPermutationTest matches explicitly calculated permutations", {
  W <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3, 3, 4),
    j = c(2, 1, 3, 2, 4, 3),
    x = 1,
    dims = c(4, 4)
  )
  x <- c(1, 2, 5, 9)

  set.seed(42)
  result <- SpaceMosaic:::.moranPermutationTest(
    x, W, n_permutations = 20, alternative = "greater"
  )

  # Recreate the same permutations and calculate their statistics through the
  # public-facing internal helper to check the optimized calculation.
  set.seed(42)
  simulated <- replicate(
    20,
    SpaceMosaic:::.moranI(sample(x), W)
  )
  expected_p <- (sum(simulated >= SpaceMosaic:::.moranI(x, W)) + 1) / 21

  expect_equal(result$observed, SpaceMosaic:::.moranI(x, W))
  expect_equal(result$expected, -1 / 3)
  expect_equal(result$p_value, expected_p)
  expect_identical(result$n_permutations, 20L)
  expect_identical(result$alternative, "greater")
})

test_that(".moranPermutationTest supports directional alternatives", {
  W <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4),
    j = c(2, 1, 4, 3),
    x = 1,
    dims = c(4, 4)
  )
  x <- c(0, 0, 1, 1)

  set.seed(7)
  greater <- SpaceMosaic:::.moranPermutationTest(
    x, W, n_permutations = 19, alternative = "greater"
  )
  set.seed(7)
  less <- SpaceMosaic:::.moranPermutationTest(
    x, W, n_permutations = 19, alternative = "less"
  )
  set.seed(7)
  two_sided <- SpaceMosaic:::.moranPermutationTest(
    x, W, n_permutations = 19, alternative = "two.sided"
  )

  expect_equal(greater$observed, 1)
  expect_true(greater$p_value <= less$p_value)
  expect_true(two_sided$p_value >= 1 / 20)
  expect_true(two_sided$p_value <= 1)
})

test_that(".moranPermutationTest validates the number of permutations", {
  W <- matrix(
    c(0, 1, 0,
      1, 0, 1,
      0, 1, 0),
    nrow = 3,
    byrow = TRUE
  )

  expect_error(
    SpaceMosaic:::.moranPermutationTest(c(1, 2, 3), W, 0),
    "positive integer",
    fixed = TRUE
  )
  expect_error(
    SpaceMosaic:::.moranPermutationTest(c(1, 2, 3), W, 1.5),
    "positive integer",
    fixed = TRUE
  )
})
