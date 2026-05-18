#' DE based on limma's linear model fitting
#'
#' This function performs a linear regression of each gene on the predictors in `df`, returning effect sizes, standard errors, p-values, and residuals.
#' @param y the normalized counts matrix (cells x genes)
#' @param df a data frame of the variables to be modeled
#' @return A list with:
#'   - effect: matrix of effect sizes (genes x predictors)
#'   - se:     matrix of standard errors
#'   - p:      matrix of p-values
#'   - residuals: matrix of residuals (cells x genes)
#' 
#' @importFrom limma lmFit eBayes
#'
#' @rdname limmaDE
#'
#' @export
limmaDE <- function(y, df) {
  if (!is.matrix(y) && !inherits(y, "Matrix")) {
    y <- as.matrix(y)
  }
  X <- model.matrix(~ ., df)  # intercept + predictors
  if (nrow(X) != nrow(y)) stop("nrow(y) must equal nrow(df).")
  
  # Fit the linear model using limma
  fit <- limma::lmFit(t(y), X)
  # Apply empirical Bayes moderation of the standard errors
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  
  # Effects and p-values
  effect <- fit$coefficients[, -1, drop = FALSE]  # G x k
  se     <- fit$stdev.unscaled[, -1, drop = FALSE] * sqrt(fit$s2.post) 
  p      <- fit$p.value[, -1, drop = FALSE]

  # Residuals: cells x genes
  residuals <- t(residuals(fit, y = t(y))) 
  
  list(
    effect = effect, 
    se = se, 
    p = p, 
    residuals = residuals
  )
}