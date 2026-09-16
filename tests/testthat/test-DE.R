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

test_that("limmaDE returns coherent moderated inference", {
  df <- data.frame(
    treatment = c(0, 0, 1, 1, 2, 2),
    batch = factor(c("a", "b", "a", "b", "a", "b"))
  )
  y <- cbind(
    gene_a = c(2, 4, 5, 8, 9, 13),
    gene_b = c(10, 8, 9, 5, 6, 2),
    gene_c = c(4, 5, 3, 7, 8, 11)
  )
  rownames(y) <- paste0("cell_", seq_len(nrow(y)))
  design <- stats::model.matrix(~ ., df)
  expected <- limma::eBayes(
    limma::lmFit(t(y), design), trend = FALSE, robust = TRUE
  )

  result <- limmaDE(
    y, df, return_residuals = TRUE, trend = FALSE, robust = TRUE
  )
  expected_se <- sweep(
    expected$stdev.unscaled[, -1, drop = FALSE],
    1, sqrt(expected$s2.post), `*`
  )
  expected_residuals <- t(
    limma::residuals.MArrayLM(expected, y = t(y))
  )

  expect_equal(result$effect, expected$coefficients[, -1, drop = FALSE])
  expect_equal(result$se, expected_se)
  expect_equal(result$p, expected$p.value[, -1, drop = FALSE])
  expect_equal(unname(result$sigma2), unname(expected$sigma^2))
  expect_equal(unname(result$sigma2_post), unname(expected$s2.post))
  expect_equal(unname(result$df_resid), unname(expected$df.residual))
  expect_equal(unname(result$df_total), unname(expected$df.total))
  expect_equal(result$residuals, expected_residuals)
  expect_identical(dimnames(result$residuals), dimnames(y))
})

test_that("limmaDE computes residuals only when requested", {
  y <- cbind(gene_a = c(1, 3, 5, 8), gene_b = c(8, 6, 4, 1))
  df <- data.frame(group = c(0, 0, 1, 1))

  result <- limmaDE(y, df, trend = FALSE)

  expect_false("residuals" %in% names(result))
})

test_that("limmaDE preserves constant predictors as NA", {
  y <- cbind(
    gene_a = c(1, 3, 5, 8, 9, 12),
    gene_b = c(8, 6, 4, 3, 2, 1),
    gene_c = c(2, 3, 2, 5, 6, 8)
  )
  df <- data.frame(constant = 1, varying = c(0, 0, 1, 1, 2, 2))

  expect_warning(
    result <- limmaDE(y, df, trend = FALSE),
    "Dropped zero-variance predictors: constant",
    fixed = TRUE
  )

  expect_true(all(is.na(result$effect[, "constant"])))
  expect_true(all(is.finite(result$effect[, "varying"])))
  expect_identical(colnames(result$effect), c("constant", "varying"))
})

test_that("limmaDE validates logical options", {
  y <- cbind(gene_a = c(1, 2, 4, 8), gene_b = c(8, 4, 2, 1))
  df <- data.frame(group = c(0, 0, 1, 1))

  expect_error(limmaDE(y, df, return_residuals = NA), "TRUE or FALSE")
  expect_error(limmaDE(y, df, trend = 1), "trend must be TRUE or FALSE")
  expect_error(limmaDE(y, df, robust = c(TRUE, FALSE)), "TRUE or FALSE")
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

test_that("patchDE dispatches to limmaDE", {
  y <- cbind(
    gene_a = c(1, 3, 4, 8, 9, 12, 14, 20),
    gene_b = c(10, 8, 7, 3, 9, 7, 4, 2),
    gene_c = c(2, 4, 3, 7, 6, 8, 10, 13)
  )
  rownames(y) <- paste0("cell_", seq_len(nrow(y)))
  df <- data.frame(treatment = rep(0:3, 2))
  patch <- rep(c("one", "two"), each = 4)

  result <- patchDE(
    y, df, patch,
    method = "limma",
    return_residuals = TRUE,
    verbose = FALSE
  )
  expected_one <- limmaDE(
    y[patch == "one", , drop = FALSE],
    df[patch == "one", , drop = FALSE],
    return_residuals = TRUE,
    trend = TRUE,
    robust = TRUE
  )

  expect_equal(result$de$treatment$ests[, "one"], expected_one$effect[, "treatment"])
  expect_equal(result$de$treatment$ses[, "one"], expected_one$se[, "treatment"])
  expect_equal(result$de$treatment$pvals[, "one"], expected_one$p[, "treatment"])
  expect_equal(result$residuals[patch == "one", ], expected_one$residuals)
})

test_that("patchDE disables the limma trend for Pearson residuals", {
  y <- cbind(
    gene_a = c(1, 3, 4, 8, 9, 12),
    gene_b = c(10, 8, 7, 3, 4, 2),
    gene_c = c(2, 4, 3, 7, 6, 9)
  )
  totals <- rowSums(y) + 10
  df <- data.frame(treatment = 0:5)
  patch <- rep("one", nrow(y))
  transformed <- pearsonResiduals(y, totals)

  result <- patchDE(
    y, df, patch,
    method = "limma",
    pearson = TRUE,
    tot = totals,
    verbose = FALSE
  )
  expected <- limmaDE(transformed, df, trend = FALSE, robust = TRUE)

  expect_equal(result$treatment$ests[, "one"], expected$effect[, "treatment"])
  expect_equal(result$treatment$ses[, "one"], expected$se[, "treatment"])
  expect_equal(result$treatment$pvals[, "one"], expected$p[, "treatment"])
})

test_that("patchDE validates its DE method", {
  y <- cbind(gene_a = c(1, 2, 4), gene_b = c(4, 2, 1))
  df <- data.frame(group = 0:2)

  expect_error(
    patchDE(y, df, rep("one", 3), method = "unknown", verbose = FALSE),
    "'arg' should be one of",
    fixed = TRUE
  )
})
