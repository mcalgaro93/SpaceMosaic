
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
  coef_names <- colnames(X)[-1]  # drop intercept
  k <- length(coef_names)
  if (k == 0) stop("No predictors found.")
  
  # Crossproducts
  XtX <- Matrix::crossprod(X)
  XtY <- Matrix::crossprod(X, y)
  
  # Inverse (with fallback)
  XtX_inv <- tryCatch(
    chol2inv(chol(XtX)),
    error = function(e) solve(XtX)
  )
  
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
  colnames(SE) <- coef_names
  
  # Effects and p-values
  Effect <- as.matrix(Matrix::t(B[-1, , drop = FALSE]))   # G x k
  Tstat  <- as.matrix(Effect / SE)
  Pval   <- 2 * pt(abs(Tstat), df = dfres, lower.tail = FALSE)
  
  list(
    effect     = Effect,
    se         = SE,
    p          = Pval,
    df_resid   = dfres
  )
}



#' patchDE: run DE over all patches
#' @param y Expression matrix, cells * genes
#' @param df Data frame to be used as DE predictors
#' @param patch Vector of patch IDs
#' @export
patchDE <- function(y, df, patch) {
  # get DE results per patch:
  results <- list()
  patches <- setdiff(unique(patch), NA)
  for (patchid in patches) {
    patchinds <- (patch == patchid) & !is.na(patch)
    results[[patchid]] <- hastyDE(y = y[patchinds, , drop = FALSE],
                                  df = df[patchinds, , drop = FALSE])
  }
  # reformat to a per-variable list:
  variables <- colnames(results[[1]][[1]])
  out <- list()
  for (varname in variables) {
    out[[varname]] <- list()
    out[[varname]]$pvals <- sapply(results, function(tmp){tmp$p})
    out[[varname]]$ests <- sapply(results, function(tmp){tmp$effect})
    out[[varname]]$ses <- sapply(results, function(tmp){tmp$se})
    rownames(out[[varname]]$pvals) <- rownames(out[[varname]]$ests) <- rownames(out[[varname]]$ses) <- rownames(results[[1]][[1]])
  }
  return(out)
}