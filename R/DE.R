
#' Very fast DE without best practices
#' Runs simple OLS regression on all genes at once using matrix algebra
#' @param y the normalized counts matrix (cells x genes)
#' @param df a data frame of the variables to be modeled
#' @return A list with:
#'   - effect: matrix of effect sizes (genes x predictors)
#'   - se:     matrix of standard errors
#'   - p:      matrix of p-values
#'   - df_resid: residual degrees of freedom
hastyDE <- function(y, df) {
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
    return(list(effect = na_mat, se = na_mat, p = na_mat, sigma2 = sigma2_out, df_resid = n - 1L))
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
    return(list(effect = na_mat, se = na_mat, p = na_mat, sigma2 = sigma2_out, df_resid = n - p))
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
  
  list(
    effect     = Effect,
    se         = SE,
    p          = Pval,
    sigma2     = sigma2,
    df_resid   = dfres
  )
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
#' @return Matrix of Pearson residuals (same dimensions as y)
#' @export
pearsonResiduals <- function(y, tot) {
  if (length(tot) != nrow(y)) stop("length(tot) must equal nrow(y).")
  genescale <- colMeans(y)
  genescale[genescale == 0] <- min(genescale[genescale > 0], na.rm = TRUE) # avoid division by zero
  expected <- outer(tot, genescale) / mean(tot)
  (as.matrix(y) - expected) / sqrt(expected)
}


#' patchDE: run DE over all patches
#' @param y Expression matrix, cells * genes
#' @param df Data frame to be used as DE predictors
#' @param patch Vector of patch IDs
#' @param pearson Logical; if TRUE, transform y to Pearson residuals before DE
#' @param tot Numeric vector of total counts per cell (required if pearson = TRUE)
#' @param resid_mse Logical; if TRUE, include per-gene residual MSE in output
#' @param verbose Show progress. Default TRUE.
#' @export
patchDE <- function(y, df, patch, pearson = FALSE, tot = NULL, resid_mse = FALSE, verbose = TRUE) {
  if (pearson) {
    if (is.null(tot)) {
      stop("tot must be provided when pearson = TRUE.")
    }
  }
  # get DE results per patch:
  results <- list()
  patches <- setdiff(unique(patch), NA)
  if (verbose) cli::cli_progress_bar("patchDE", total = length(patches))
  for (patchid in patches) {
    patchinds <- (patch == patchid) & !is.na(patch)
    ysub <- y[patchinds, , drop = FALSE]
    if (pearson) {
      ysub <- pearsonResiduals(ysub, tot = tot[patchinds])
    }
    results[[patchid]] <- hastyDE(y = ysub, df = df[patchinds, , drop = FALSE])
    if (verbose) cli::cli_progress_update()
  }
  if (verbose) cli::cli_progress_done()
  # reformat to a per-variable list:
  variables <- colnames(results[[1]][[1]])
  out <- list()
  for (varname in variables) {
    out[[varname]] <- list()
    out[[varname]]$pvals <- sapply(results, function(tmp) {
      tmp$p[, varname]
    })
    out[[varname]]$ests <- sapply(results, function(tmp) {
      tmp$effect[, varname]
    })
    out[[varname]]$ses <- sapply(results, function(tmp) {
      tmp$se[, varname]
    })
   rownames(out[[varname]]$pvals) <- rownames(out[[varname]]$ests) <- rownames(out[[varname]]$ses) <- rownames(results[[1]][[1]])
  }
  if (resid_mse) {
    out[["resid_mse"]] <- sapply(results, function(tmp) tmp$sigma2)
  }
  return(out)
}