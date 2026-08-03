
#' Very fast DE without best practices
#' Runs simple OLS regression on all genes at once using matrix algebra
#' @param y the normalized counts matrix (cells x genes)
#' @param df a data frame of the variables to be modeled
#' @param return_residuals Logical; if TRUE, include the OLS residual matrix
#'   (cells x genes) in the returned list. Default FALSE.
#' @return A list with:
#'   - effect: matrix of effect sizes (genes x predictors)
#'   - se:     matrix of standard errors
#'   - p:      matrix of p-values
#'   - df_resid: residual degrees of freedom
#'   - residuals: optional matrix of OLS residuals (cells x genes)
#' @export
hastyDE <- function(y, df, return_residuals = FALSE) {
  if (!is.logical(return_residuals) || length(return_residuals) != 1L ||
      is.na(return_residuals)) {
    stop("return_residuals must be TRUE or FALSE.")
  }
  if (!is.matrix(y) && !inherits(y, "Matrix")) {
    y <- as.matrix(y)
  }
  X <- model.matrix(~ ., df)  # intercept + predictors
  if (nrow(X) != nrow(y)) stop("nrow(y) must equal nrow(df).")
  
  n <- nrow(X); p <- ncol(X)
  all_coef_names <- colnames(X)[-1]  # drop intercept
  if (length(all_coef_names) == 0) stop("No predictors found.")
  G <- ncol(y)
  
  # Drop zero-variance columns (keep intercept always)
  col_vars <- apply(X[, -1, drop = FALSE], 2, var)
  keep <- col_vars > 0
  if (!any(keep)) {
    # All predictors are constant — return all NA
    na_mat <- matrix(NA_real_, nrow = G, ncol = length(all_coef_names))
    colnames(na_mat) <- all_coef_names
    rownames(na_mat) <- colnames(y)
    sigma2_out <- rep(NA_real_, G)
    names(sigma2_out) <- colnames(y)
    result <- list(
      effect = na_mat,
      se = na_mat,
      p = na_mat,
      sigma2 = sigma2_out,
      df_resid = n - 1L
    )
    if (return_residuals) {
      # With no varying predictors, the estimable model contains only an
      # intercept, so its residuals are the column-centered observations.
      result$residuals <- as.matrix(
        sweep(y, 2, Matrix::colMeans(y), FUN = "-")
      )
      dimnames(result$residuals) <- dimnames(y)
    }
    return(result)
  }
  dropped <- all_coef_names[!keep]
  if (length(dropped) > 0) {
    warning("Dropped zero-variance predictors: ", paste(dropped, collapse = ", "))
    X <- X[, c(TRUE, keep), drop = FALSE]  # keep intercept + non-constant cols
    p <- ncol(X)
  }
  coef_names <- colnames(X)[-1]
  
  # Crossproducts
  XtX <- Matrix::crossprod(X)
  XtY <- Matrix::crossprod(X, y)
  
  # Inverse (with fallback)
  XtX_inv <- tryCatch(
    chol2inv(chol(XtX)),
    error = function(e) {
      tryCatch(solve(XtX), error = function(e2) NULL)
    }
  )
  if (is.null(XtX_inv)) {
    # Degenerate even after dropping zero-variance columns — return all NA
    na_mat <- matrix(NA_real_, nrow = G, ncol = length(all_coef_names))
    colnames(na_mat) <- all_coef_names
    rownames(na_mat) <- colnames(y)
    sigma2_out <- rep(NA_real_, G)
    names(sigma2_out) <- colnames(y)
    result <- list(
      effect = na_mat,
      se = na_mat,
      p = na_mat,
      sigma2 = sigma2_out,
      df_resid = n - p
    )
    if (return_residuals) {
      # No unique coefficient solution is available for this design, so fitted
      # values and residuals cannot be reported reliably.
      result$residuals <- matrix(
        NA_real_,
        nrow = nrow(y),
        ncol = ncol(y),
        dimnames = dimnames(y)
      )
    }
    return(result)
  }
  
  # Coefficients: p x G
  B <- XtX_inv %*% XtY
  rownames(B) <- colnames(X)
  
  # Residual variance
  yty   <- Matrix::colSums(y * y)
  RSS   <- yty - Matrix::colSums(B * XtY)
  dfres <- n - p
  sigma2 <- RSS / dfres  # length G
  
  # Standard errors: for each predictor j, sqrt(sigma2 * V_jj)
  Vdiag <- diag(XtX_inv)
  SE <- vapply((2:p), function(j) sqrt(sigma2 * Vdiag[j]), numeric(length(sigma2)))
  if (is.matrix(SE)) {
    colnames(SE) <- coef_names
  } else {
    SE <- matrix(SE, ncol = 1)
    colnames(SE) <- coef_names
  }
  
  # Effects and p-values
  Effect <- as.matrix(Matrix::t(B[-1, , drop = FALSE]))   # G x k
  Tstat  <- as.matrix(Effect / SE)
  Pval   <- 2 * pt(abs(Tstat), df = dfres, lower.tail = FALSE)
  
  # Reinsert NA columns for dropped predictors
  if (length(dropped) > 0) {
    Effect <- .reinsertNA(Effect, all_coef_names, coef_names, G, colnames(y))
    SE     <- .reinsertNA(SE, all_coef_names, coef_names, G, colnames(y))
    Pval   <- .reinsertNA(Pval, all_coef_names, coef_names, G, colnames(y))
  }
  
  # Named residual MSE vector (per gene)
  names(sigma2) <- colnames(y)
  
  result <- list(
    effect     = Effect,
    se         = SE,
    p          = Pval,
    sigma2     = sigma2,
    df_resid   = dfres
  )

  if (return_residuals) {
    # Residuals are materialized only on request because the cells-by-genes
    # matrix can be substantially larger than the DE summary statistics.
    result$residuals <- as.matrix(y - X %*% B)
    dimnames(result$residuals) <- dimnames(y)
  }

  result
}

# Helper: reinsert NA columns for dropped predictors
.reinsertNA <- function(mat, all_names, kept_names, G, gene_names) {
  full <- matrix(NA_real_, nrow = G, ncol = length(all_names))
  colnames(full) <- all_names
  rownames(full) <- gene_names
  full[, kept_names] <- mat
  full
}



#' Convert a raw counts matrix to Pearson residuals
#' @param y Raw counts matrix (cells x genes)
#' @param tot Numeric vector of total counts per cell (length = nrow(y))
#'   containing finite, non-negative values.
#' @return Matrix of Pearson residuals (same dimensions as y)
#' @export
pearsonResiduals <- function(y, tot) {
  if (length(tot) != nrow(y)) stop("length(tot) must equal nrow(y).")
  if (any(!is.finite(tot)) || any(tot < 0)) {
    stop("tot must contain finite, non-negative values.")
  }

  y <- as.matrix(y)
  zero_residuals <- matrix(
    0,
    nrow = nrow(y),
    ncol = ncol(y),
    dimnames = dimnames(y)
  )
  mean_tot <- mean(tot)
  genescale <- colMeans(y)

  if (mean_tot == 0) {
    return(zero_residuals)
  }

  expected <- outer(tot, genescale) / mean_tot
  residuals <- (y - expected) / sqrt(expected)
  # A gene with mean zero has zero expected expression in every cell. Its
  # Pearson residual is defined as zero rather than introducing an arbitrary
  # positive gene scale solely to avoid division by zero.
  residuals[expected == 0] <- 0
  residuals
}


#' patchDE: run DE over all patches
#' @param y Expression matrix, cells * genes
#' @param df Data frame to be used as DE predictors
#' @param patch Vector of patch IDs
#' @param pearson Logical; if TRUE, transform y to Pearson residuals before DE
#' @param tot Numeric vector of total counts per cell (required if pearson = TRUE)
#' @param resid_mse Logical; if TRUE, include per-gene residual MSE in output
#' @param return_residuals Logical; if TRUE, return an OLS residual matrix
#'   aligned with the rows and columns of `y`. Default FALSE.
#' @param verbose Show progress. Default TRUE.
#' @return A list keyed by model variable. Each element contains `pvals`,
#'   `ests`, and `ses` matrices with genes in rows and patches in columns.
#'   If `resid_mse = TRUE`, the list also contains a `resid_mse` matrix with
#'   the same orientation. If `return_residuals = TRUE`, a list with components
#'   `de` (the usual result) and `residuals` is returned instead. `residuals` is
#'   a cells by genes matrix aligned with `y`; rows whose patch is missing
#'   contain `NA`.
#' @examples
#' y <- cbind(gene_a = c(2, 4, 5, 8, 9, 13))
#' df <- data.frame(treatment = c(0, 0, 1, 1, 2, 2))
#' patch <- rep("patch_1", nrow(y))
#'
#' fit <- patchDE(
#'   y, df, patch,
#'   return_residuals = TRUE,
#'   verbose = FALSE
#' )
#' fit$de$treatment$ests
#' fit$residuals
#'
#' @export
patchDE <- function(y, df, patch, pearson = FALSE, tot = NULL,
                    resid_mse = FALSE, return_residuals = FALSE,
                    verbose = TRUE) {
  if (length(patch) != nrow(y) || nrow(df) != nrow(y)) {
    stop("nrow(y), nrow(df), and length(patch) must be equal.")
  }
  if (!is.logical(return_residuals) || length(return_residuals) != 1L ||
      is.na(return_residuals)) {
    stop("return_residuals must be TRUE or FALSE.")
  }
  if (pearson) {
    if (is.null(tot)) {
      stop("tot must be provided when pearson = TRUE.")
    }
  }

  # get DE results per patch:
  results <- list()
  if (return_residuals) {
    residual_matrix <- matrix(
      NA_real_,
      nrow = nrow(y),
      ncol = ncol(y),
      dimnames = dimnames(y)
    )
  }
  patches <- unique(patch[!is.na(patch)])
  if (length(patches) == 0L) {
    stop("patch must contain at least one non-missing patch ID.")
  }
  if (verbose) {
    cli::cli_progress_bar("patchDE", total = length(patches))
  }
  for (patchid in patches) {
    cell_index <- which(!is.na(patch) & patch == patchid)
    ysub <- y[cell_index, , drop = FALSE]
    if (pearson) {
      ysub <- pearsonResiduals(ysub, tot = tot[cell_index])
    }
    patch_name <- as.character(patchid)
    results[[patch_name]] <- hastyDE(
      y = ysub,
      df = df[cell_index, , drop = FALSE],
      return_residuals = return_residuals
    )
    if (return_residuals) {
      residual_matrix[cell_index, ] <- results[[patch_name]]$residuals
      # The residual matrix is returned separately and is not needed while the
      # DE summaries are reformatted below.
      results[[patch_name]]$residuals <- NULL
    }
    if (verbose) cli::cli_progress_update()
  }
  if (verbose) cli::cli_progress_done()
  # reformat to a per-variable list:
  variables <- colnames(results[[1]][[1]])
  out <- list()
  for (varname in variables) {
    out[[varname]] <- list()
    out[[varname]]$pvals <- do.call(
      cbind,
      lapply(results, function(tmp) tmp$p[, varname])
    )
    out[[varname]]$ests <- do.call(
      cbind,
      lapply(results, function(tmp) tmp$effect[, varname])
    )
    out[[varname]]$ses <- do.call(
      cbind,
      lapply(results, function(tmp) tmp$se[, varname])
    )
    colnames(out[[varname]]$pvals) <-
      colnames(out[[varname]]$ests) <-
      colnames(out[[varname]]$ses) <- names(results)
    rownames(out[[varname]]$pvals) <-
      rownames(out[[varname]]$ests) <-
      rownames(out[[varname]]$ses) <- rownames(results[[1]][[1]])
  }
  if (resid_mse) {
    out[["resid_mse"]] <- do.call(
      cbind,
      lapply(results, function(tmp) tmp$sigma2)
    )
    colnames(out[["resid_mse"]]) <- names(results)
  }
  if (return_residuals) {
    return(list(de = out, residuals = residual_matrix))
  }
  out
}
