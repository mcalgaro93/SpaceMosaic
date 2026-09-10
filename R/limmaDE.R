#' Differential expression with limma
#'
#' Fits the same linear model to every gene and applies empirical-Bayes
#' moderation to the gene-wise residual variances.
#'
#' @param y Numeric expression matrix with cells in rows and genes in columns.
#' @param df Data frame containing the model predictors, with one row per cell.
#' @param return_residuals Logical; if `TRUE`, include the model residuals
#'   (cells by genes). Default `FALSE`.
#' @param trend Logical passed to [limma::eBayes()]. Use `FALSE` for Pearson
#'   residuals and consider `TRUE` for log-expression retaining a
#'   mean--variance trend. Default `TRUE`.
#' @param robust Logical passed to [limma::eBayes()] to use robust empirical
#'   Bayes hyperparameter estimation. Default `TRUE`.
#'
#' @return A list containing gene-by-predictor matrices `effect`, `se`, and
#'   `p`; gene-wise vectors `sigma2` (raw residual variance), `sigma2_post`
#'   (moderated posterior residual variance), `df_resid`, and `df_total`; and,
#'   when requested, a cells-by-genes `residuals` matrix.
#'
#' @importFrom limma eBayes lmFit residuals.MArrayLM
#' @export
limmaDE <- function(y, df, return_residuals = FALSE, trend = TRUE,
                    robust = TRUE) {
  logical_args <- list(
    return_residuals = return_residuals,
    trend = trend,
    robust = robust
  )
  for (arg_name in names(logical_args)) {
    value <- logical_args[[arg_name]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(arg_name, " must be TRUE or FALSE.")
    }
  }

  if (!is.matrix(y) && !inherits(y, "Matrix")) {
    y <- as.matrix(y)
  }
  X <- stats::model.matrix(~ ., df)
  if (nrow(X) != nrow(y)) stop("nrow(y) must equal nrow(df).")

  all_coef_names <- colnames(X)[-1L]
  if (length(all_coef_names) == 0L) stop("No predictors found.")

  # Within small patches, factor levels or numeric predictors can be constant.
  # Fit the estimable part of the model and reinsert those terms as NA.
  varying <- vapply(
    seq.int(2L, ncol(X)),
    function(j) {
      values <- X[, j]
      all(is.finite(values)) && length(unique(values)) > 1L
    },
    logical(1)
  )
  dropped <- all_coef_names[!varying]
  if (length(dropped) > 0L) {
    warning("Dropped zero-variance predictors: ", paste(dropped, collapse = ", "))
  }
  X_fit <- X[, c(TRUE, varying), drop = FALSE]

  fit <- limma::lmFit(t(y), X_fit)
  fit <- limma::eBayes(fit, trend = trend, robust = robust)

  kept_coef_names <- colnames(X_fit)[-1L]
  gene_names <- colnames(y)
  n_genes <- ncol(y)
  effect <- matrix(
    NA_real_, nrow = n_genes, ncol = length(all_coef_names),
    dimnames = list(gene_names, all_coef_names)
  )
  se <- p <- effect
  if (length(kept_coef_names) > 0L) {
    effect[, kept_coef_names] <- fit$coefficients[, -1L, drop = FALSE]
    # The standard errors are the unscaled standard errors multiplied by the
    # posterior residual standard deviation.
    se[, kept_coef_names] <- sweep(
      fit$stdev.unscaled[, -1L, drop = FALSE],
      1L, sqrt(fit$s2.post), `*`
    )
    p[, kept_coef_names] <- fit$p.value[, -1L, drop = FALSE]
  }

  sigma2 <- fit$sigma^2
  names(sigma2) <- gene_names
  sigma2_post <- fit$s2.post
  names(sigma2_post) <- gene_names
  df_resid <- fit$df.residual
  names(df_resid) <- gene_names
  df_total <- fit$df.total
  names(df_total) <- gene_names

  result <- list(
    effect = effect,
    se = se,
    p = p,
    sigma2 = sigma2, # raw residual variance
    sigma2_post = sigma2_post, # moderated posterior residual variance
    df_resid = df_resid, # n - rank(X) degrees of freedom
    df_total = df_total # df_resid + df_prior, df_prior = from empirical Bayes
  )

  if (return_residuals) {
    result$residuals <- t(limma::residuals.MArrayLM(fit, y = t(y)))
    dimnames(result$residuals) <- dimnames(y)
  }

  result
}
