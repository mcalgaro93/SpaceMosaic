test_that("hastyDE returns OLS residuals only when requested", {
  df <- data.frame(
    treatment = c(0, 0, 1, 1, 2, 2),
    batch = factor(c("a", "b", "a", "b", "a", "b"))
  )
  y <- cbind(
    gene_a = c(2, 4, 5, 8, 9, 13),
    gene_b = c(10, 8, 9, 5, 6, 2)
  )
  rownames(y) <- paste0("cell_", seq_len(nrow(y)))

  default_result <- hastyDE(y, df)
  residual_result <- hastyDE(y, df, return_residuals = TRUE)

  expected <- vapply(
    seq_len(ncol(y)),
    function(j) stats::residuals(stats::lm(y[, j] ~ ., data = df)),
    numeric(nrow(y))
  )
  dimnames(expected) <- dimnames(y)

  expect_false("residuals" %in% names(default_result))
  expect_equal(residual_result$residuals, expected)
  expect_identical(dimnames(residual_result$residuals), dimnames(y))
})

test_that("hastyDE returns intercept-only residuals for constant predictors", {
  y <- cbind(gene_a = c(1, 3, 5), gene_b = c(8, 5, 2))
  rownames(y) <- paste0("cell_", seq_len(nrow(y)))
  df <- data.frame(group = rep(1, nrow(y)))

  result <- hastyDE(y, df, return_residuals = TRUE)

  expect_equal(result$residuals, sweep(y, 2, colMeans(y), FUN = "-"))
  expect_true(all(is.na(result$effect)))
})

test_that("hastyDE validates return_residuals", {
  y <- matrix(1:6, ncol = 2)
  df <- data.frame(group = c(0, 1, 2))

  expect_error(
    hastyDE(y, df, return_residuals = NA),
    "TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    hastyDE(y, df, return_residuals = c(TRUE, FALSE)),
    "TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("pearsonResiduals leaves all-zero genes at zero", {
  y <- cbind(
    zero_gene = c(0, 0, 0, 0),
    expressed_gene = c(1, 2, 4, 8)
  )
  totals <- c(10, 20, 30, 40)

  residuals <- pearsonResiduals(y, totals)

  expect_equal(residuals[, "zero_gene"], rep(0, nrow(y)))
  expect_true(all(is.finite(residuals)))
})

test_that("pearsonResiduals is invariant to filtering other genes", {
  y <- cbind(
    zero_gene = c(0, 0, 0, 0),
    gene_a = c(1, 2, 4, 8),
    gene_b = c(7, 3, 2, 1)
  )
  totals <- c(10, 20, 30, 40)
  full <- pearsonResiduals(y, totals)

  for (gene in colnames(y)) {
    filtered <- pearsonResiduals(y[, gene, drop = FALSE], totals)
    expect_equal(filtered[, gene], full[, gene])
  }
})

test_that("pearsonResiduals retains the standard formula for expressed genes", {
  y <- cbind(gene_a = c(1, 2, 4, 8), gene_b = c(7, 3, 2, 1))
  totals <- c(10, 20, 30, 40)
  gene_scale <- colMeans(y)
  expected <- outer(totals, gene_scale) / mean(totals)
  standard_residuals <- (y - expected) / sqrt(expected)

  expect_equal(pearsonResiduals(y, totals), standard_residuals)
})

test_that("patchDE returns residuals aligned with its filtered inputs", {
  y <- cbind(
    gene_a = c(1, 3, 4, 8, 9, 12, 14, 20),
    gene_b = c(10, 8, 7, 3, 9, 7, 4, 2)
  )
  rownames(y) <- paste0("cell_", seq_len(nrow(y)))
  df <- data.frame(treatment = rep(0:3, 2))
  patch <- rep(c("one", "two"), each = 4)
  selected_cells <- patch == "two"

  result <- patchDE(
    y[selected_cells, "gene_b", drop = FALSE],
    df[selected_cells, , drop = FALSE],
    patch[selected_cells],
    return_residuals = TRUE,
    verbose = FALSE
  )
  expected <- hastyDE(
    y[selected_cells, "gene_b", drop = FALSE],
    df[selected_cells, , drop = FALSE],
    return_residuals = TRUE
  )

  expect_named(result, c("de", "residuals"))
  expect_equal(result$residuals, expected$residuals)
  expect_identical(dimnames(result$residuals), dimnames(expected$residuals))
  expect_identical(rownames(result$de$treatment$pvals), "gene_b")
  expect_identical(colnames(result$de$treatment$pvals), "two")
})

test_that("patchDE leaves residuals missing for cells without a patch", {
  y <- cbind(gene_a = c(1, 3, 6, 5, 7, 8, 12))
  rownames(y) <- paste0("cell_", seq_len(nrow(y)))
  df <- data.frame(treatment = c(0, 1, 2, 0, 0, 1, 2))
  patch <- c("one", "one", "one", NA, "two", "two", "two")

  result <- patchDE(
    y, df, patch,
    return_residuals = TRUE,
    verbose = FALSE
  )

  expect_true(is.na(result$residuals[4, "gene_a"]))
  expect_true(all(is.finite(result$residuals[-4, "gene_a"])))
  expect_identical(dimnames(result$residuals), dimnames(y))
})
