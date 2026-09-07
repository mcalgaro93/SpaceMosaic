#' Test residual spatial autocorrelation with Moran's I
#'
#' Build a symmetric binary k-nearest-neighbor graph within each patch and test
#' global Moran's I by randomly permuting the supplied residuals. A vector
#' performs one test; a cells by genes matrix performs one test per gene and
#' patch. Residuals may come from [patchDE()] with
#' `return_residuals = TRUE`, from [hastyDE()], or from any other fitted model.
#'
#' Moran's I measures whether observations connected in the spatial graph tend
#' to have similar values. Values larger than expected indicate that similar
#' residuals cluster in space; values smaller than expected indicate that
#' neighboring residuals tend to be dissimilar. A value near its expectation
#' under random spatial arrangement provides little evidence of residual
#' spatial autocorrelation.
#'
#' The observed statistic is compared with statistics obtained by randomly
#' permuting the residuals while keeping the spatial graph fixed. The empirical
#' p-value is calculated as `(b + 1) / (n_permutations + 1)`, where `b` is the
#' number of permuted statistics at least as extreme as the observed one. Thus,
#' the smallest attainable p-value is `1 / (n_permutations + 1)`.
#'
#' @param residuals Numeric vector with one residual per observation, or a
#'   numeric matrix with observations in rows and genes in columns.
#' @param xy Numeric matrix with two columns containing spatial coordinates.
#'   Rows must be in the same order as `residuals`.
#' @param patch Optional vector of patch IDs, with one value per observation.
#'   Moran's I is calculated independently within each non-missing patch. NULL
#'   treats all observations as a single patch named `"all"`.
#' @param k Number of nearest neighbors used to construct the graph. If `k` is
#'   not smaller than the number of observations, it is reduced to `n - 1` with
#'   a warning. Default 10.
#' @param n_permutations Number of random permutations used to calculate the
#'   empirical p-value. Default 999.
#' @param alternative Direction of the alternative hypothesis: `"greater"`
#'   tests for positive spatial autocorrelation, `"less"` for negative spatial
#'   autocorrelation, and `"two.sided"` for either direction.
#' @param p_adjust_method Method passed to [stats::p.adjust()] for multiple-test
#'   correction. Use `"none"` to leave p-values unadjusted. Default `"BH"`.
#' @param adjustment_scope Scope of the multiple-test correction: `"global"`
#'   across every gene-patch test, `"patch"` separately within each patch, or
#'   `"gene"` separately within each gene. Default `"global"`.
#'
#' @return A data frame with one row per gene and patch. It contains patch and
#'   gene identifiers, number of cells, observed and expected Moran's I, raw and
#'   adjusted empirical p-values, graph and permutation settings, and `status`.
#'   A status of `"ok"` indicates a completed test; non-testable inputs such as
#'   constant or non-finite residuals return missing statistics with an
#'   explanatory status rather than stopping the remaining tests.
#'
#' @examples
#' xy <- cbind(x = 1:5, y = rep(0, 5))
#' residuals <- c(-2, -1, 0, 1, 2)
#' set.seed(1)
#' moranTest(residuals, xy, k = 2, n_permutations = 19)
#'
#' @export
moranTest <- function(residuals, xy, patch = NULL, k = 10L,
                      n_permutations = 999L,
                      alternative = c("greater", "less", "two.sided"),
                      p_adjust_method = "BH",
                      adjustment_scope = c("global", "patch", "gene")) {
  if (is.numeric(residuals) && is.null(dim(residuals))) {
    residuals <- matrix(residuals, ncol = 1L)
    colnames(residuals) <- "residual"
  } else if (is.matrix(residuals)) {
    if (!is.numeric(residuals)) {
      stop("residuals must be a numeric vector or matrix.")
    }
    residuals <- as.matrix(residuals)
    if (is.null(colnames(residuals))) {
      colnames(residuals) <- paste0("gene_", seq_len(ncol(residuals)))
    }
  } else {
    stop("residuals must be a numeric vector or matrix.")
  }
  if (ncol(residuals) < 1L || anyDuplicated(colnames(residuals))) {
    stop("residuals must have at least one uniquely named column.")
  }
  if (!is.matrix(xy) || !is.numeric(xy) || ncol(xy) != 2L ||
      any(!is.finite(xy))) {
    stop("xy must be a numeric matrix with exactly two finite columns.")
  }
  if (nrow(residuals) != nrow(xy)) {
    stop("nrow(residuals) must equal nrow(xy).")
  }
  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) || k < 1 ||
      k != floor(k)) {
    stop("k must be a positive integer.")
  }
  if (length(n_permutations) != 1L || !is.numeric(n_permutations) ||
      !is.finite(n_permutations) || n_permutations < 1 ||
      n_permutations != floor(n_permutations)) {
    stop("n_permutations must be a positive integer.")
  }
  k <- as.integer(k)
  n_permutations <- as.integer(n_permutations)

  if (is.null(patch)) {
    patch <- rep("all", nrow(residuals))
  } else if (length(patch) != nrow(residuals)) {
    stop("length(patch) must equal nrow(residuals).")
  }
  patch_names <- unique(as.character(patch[!is.na(patch)]))
  if (length(patch_names) == 0L) {
    stop("patch must contain at least one non-missing patch ID.")
  }

  alternative <- match.arg(alternative)
  adjustment_scope <- match.arg(adjustment_scope)
  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L ||
      !p_adjust_method %in% c("none", stats::p.adjust.methods)) {
    stop("p_adjust_method must be 'none' or a method supported by p.adjust().")
  }

  results <- vector("list", length(patch_names) * ncol(residuals))
  result_index <- 0L

  for (patch_name in patch_names) {
    cell_index <- which(!is.na(patch) & as.character(patch) == patch_name)
    n_cells <- length(cell_index)
    patch_residuals <- residuals[cell_index, , drop = FALSE]

    # A graph cannot support global Moran's I with fewer than three nodes.
    if (n_cells >= 3L) {
      W <- .buildContiguityGraph(xy[cell_index, , drop = FALSE], k = k)
      effective_k <- min(k, n_cells - 1L)
      s0 <- sum(W)
    } else {
      W <- NULL
      effective_k <- NA_integer_
      s0 <- NULL
    }

    for (gene in colnames(residuals)) {
      result_index <- result_index + 1L
      values <- patch_residuals[, gene]
      status <- "ok"
      test_result <- NULL

      if (n_cells < 3L) {
        status <- "fewer than three cells"
      } else if (any(!is.finite(values))) {
        status <- "non-finite residuals"
      } else if (all(values == values[1L])) {
        status <- "constant residuals"
      } else {
        test_result <- .moranPermutationTest(
          x = values,
          W = W,
          n_permutations = n_permutations,
          alternative = alternative,
          s0 = s0
        )
      }

      results[[result_index]] <- data.frame(
        patch = patch_name,
        gene = gene,
        n_cells = n_cells,
        observed = if (is.null(test_result)) NA_real_ else test_result$observed,
        expected = if (is.null(test_result)) NA_real_ else test_result$expected,
        p_value = if (is.null(test_result)) NA_real_ else test_result$p_value,
        p_adjusted = NA_real_,
        k = k,
        effective_k = effective_k,
        n_permutations = n_permutations,
        alternative = alternative,
        status = status,
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out$p_adjusted <- .adjustMoranPValues(
    out$p_value,
    patch = out$patch,
    gene = out$gene,
    method = p_adjust_method,
    scope = adjustment_scope
  )
  out
}


#' Adjust p-values from multiple Moran tests
#'
#' Apply a multiple-testing correction globally or independently within each
#' patch or gene while preserving missing values from non-testable inputs.
#'
#' @param p Numeric vector of raw p-values.
#' @param patch,gene Character identifiers aligned with `p`.
#' @param method Adjustment method accepted by [stats::p.adjust()], or `"none"`.
#' @param scope One of `"global"`, `"patch"`, or `"gene"`.
#'
#' @return Numeric vector of adjusted p-values aligned with `p`.
#'
#' @noRd
.adjustMoranPValues <- function(p, patch, gene, method, scope) {
  if (method == "none") {
    return(p)
  }
  adjusted <- rep(NA_real_, length(p))
  groups <- switch(
    scope,
    global = rep("all", length(p)),
    patch = patch,
    gene = gene
  )
  for (group in unique(groups)) {
    use <- groups == group & !is.na(p)
    adjusted[use] <- stats::p.adjust(p[use], method = method)
  }
  adjusted
}


#' Calculate global Moran's I
#'
#' Internal helper that calculates global Moran's I for a numeric vector and a
#' user-supplied spatial weights matrix. The weights matrix is expected to be
#' symmetric, non-negative, and diagonal-free. Weight normalization is handled
#' by the usual `n / S0` factor, where `S0` is the sum of all weights.
#'
#' @param x Numeric vector containing one value per spatial observation.
#' @param W Square dense or sparse matrix of spatial weights. Its dimensions
#'   must match `length(x)`.
#'
#' @return A numeric scalar containing global Moran's I.
#'
#' @noRd
.moranI <- function(x, W) {
  moran_data <- .prepareMoran(x, W)

  .moranIFromCentered(
    z = moran_data$z,
    W = W,
    scale = moran_data$scale
  )
}


#' Prepare inputs for Moran's I calculations
#'
#' Validate the observations and spatial weights, then pre-compute the centered
#' observations and the constant scale factor shared by the observed statistic
#' and all its permutations.
#'
#' @inheritParams .moranI
#'
#' @return A list containing the centered values in `z` and the Moran scale
#'   factor in `scale`.
#'
#' @noRd
.prepareMoran <- function(x, W, s0 = NULL) {
  n <- length(x)

  if (!is.numeric(x) || n < 3L || any(!is.finite(x))) {
    stop("x must contain at least three finite numeric values.")
  }
  if (!is.matrix(W) && !inherits(W, "Matrix")) {
    stop("W must be a matrix or a Matrix object.")
  }
  if (nrow(W) != n || ncol(W) != n) {
    stop("W must be a square matrix matching length(x).")
  }
  if (is.null(s0)) {
    if (any(!is.finite(W))) {
      stop("W must contain only finite weights.")
    }
    if (any(W < 0)) {
      stop("W must contain non-negative weights.")
    }
    if (any(Matrix::diag(W) != 0)) {
      stop("W must have a zero diagonal.")
    }
    if (!isSymmetric(W)) {
      stop("W must be symmetric.")
    }
    s0 <- sum(W)
  }

  # Moran's I is calculated on values centered around their sample mean.
  z <- x - mean(x)
  denominator <- sum(z^2)
  if (denominator == 0) {
    stop("Moran's I is undefined for a constant vector.")
  }

  # S0 is the total weight in the spatial network. A zero value means that
  # none of the observations are connected.
  if (length(s0) != 1L || !is.finite(s0) || s0 <= 0) {
    stop("Moran's I is undefined when W has no edges.")
  }

  list(z = z, scale = n / (s0 * denominator))
}


#' Calculate Moran's I from pre-processed values
#'
#' Fast internal calculation used after validation and centering have already
#' been performed. Keeping the constant scale factor outside this helper avoids
#' recalculating unchanged quantities during permutation tests.
#'
#' @param z Centered numeric observation vector.
#' @inheritParams .moranI
#' @param scale Pre-computed value of `n / (S0 * sum(z^2))`.
#'
#' @return A numeric scalar containing global Moran's I.
#'
#' @noRd
.moranIFromCentered <- function(z, W, scale) {
  # Matrix multiplication preserves sparsity in W and avoids constructing the
  # dense outer product of all pairs of observations.
  numerator <- as.numeric(Matrix::crossprod(z, W %*% z))

  as.numeric(scale * numerator)
}


#' Permutation test for global Moran's I
#'
#' Calculate an empirical p-value by permuting the observations while keeping
#' the spatial weights fixed. Inputs and quantities that remain constant across
#' permutations are prepared only once.
#'
#' @inheritParams .moranI
#' @param n_permutations Number of random permutations. Default 999.
#' @param alternative Direction of the alternative hypothesis: `"greater"`,
#'   `"less"`, or `"two.sided"`.
#' @param s0 Optional pre-computed sum of weights. Used by batch calculations
#'   after the shared graph has already been validated.
#'
#' @return A list with the observed and expected Moran's I, the empirical
#'   p-value, the alternative hypothesis, and the number of permutations.
#'
#' @noRd
.moranPermutationTest <- function(x, W, n_permutations = 999L,
                                  alternative = c("greater", "less",
                                                  "two.sided"),
                                  s0 = NULL) {
  if (length(n_permutations) != 1L || !is.numeric(n_permutations) ||
      !is.finite(n_permutations) || n_permutations < 1 ||
      n_permutations != floor(n_permutations)) {
    stop("n_permutations must be a positive integer.")
  }
  alternative <- match.arg(alternative)
  n_permutations <- as.integer(n_permutations)

  moran_data <- .prepareMoran(x, W, s0 = s0)
  z <- moran_data$z
  scale <- moran_data$scale
  observed <- .moranIFromCentered(z, W, scale)

  # Permutation preserves the mean and sum of squared centered values, so only
  # the spatial cross-product needs to be recalculated on each iteration.
  simulated <- vapply(
    seq_len(n_permutations),
    function(i) .moranIFromCentered(sample(z), W, scale),
    numeric(1)
  )

  expected <- -1 / (length(x) - 1)
  extreme <- switch(
    alternative,
    greater = simulated >= observed,
    less = simulated <= observed,
    two.sided = abs(simulated - expected) >= abs(observed - expected)
  )

  # Adding one to numerator and denominator prevents a zero Monte Carlo
  # p-value and includes the observed arrangement in the reference set.
  p_value <- (sum(extreme) + 1) / (n_permutations + 1)

  list(
    observed = observed,
    expected = expected,
    p_value = p_value,
    alternative = alternative,
    n_permutations = n_permutations
  )
}
