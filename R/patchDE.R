
#' Very fast DE without best practices
#' Runs simple OLS regression on all genes at once using matrix algebra
#' @param y the normalized counts matrix (cells x genes)
#' @param df a data frame of the variables to be modeled
#' @return A list with:
#'   - effect: matrix of effect sizes (genes x predictors)
#'   - se:     matrix of standard errors
#'   - p:      matrix of p-values
#'   - residuals: matrix of residuals (cells x genes)
#' 
#' @importFrom Matrix colSums crossprod t
#' @importFrom stats model.matrix pt
#'
#' @rdname hastyDE
#'
#' @export
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
  
  # Residuals: cells x genes
  Fitted <- X %*% B
  Residuals <- as.matrix(y - Fitted)
  
  list(
    effect     = Effect,
    se         = SE,
    p          = Pval,
    residuals  = Residuals
  )
}

#' patchDE: run DE over all patches
#' @param spe An
#' [`SpatialExperiment`][SpatialExperiment::SpatialExperiment-class] object.
#'
#' @param patchesID The name of the reduced dimension in `spe` that contains the patch assignments. 
#' This should be the same name obtained in `getPatches()`.
#' @param method The DE method to use. Currently only "hasty" and "limma" are implemented. Future versions may include additional methods. See details for more information.
#' @param verbose If TRUE, print additional information about the DE process for each patch.
#' 
#' @details
#' Available DE methods:
#' - "hasty": A very fast DE method that performs ordinary least squares regression for each gene and predictor using matrix algebra.
#' - "limma": A DE method that uses limma's linear model fitting and empirical Bayes moderation.
#' 
#' @examples
#'
#' ## load a subset of 1k cells of the first sample and slice of the MERFISH
#' ## data available in MerfishData::MouseColonIbdCadinu2024()
#' suppressPackageStartupMessages({
#'     library(HDF5Array)
#'     library(SpatialExperiment)
#' })
#'
#' fname <- system.file("extdata", "MerfishData1k", package="SpaceMosaic")
#' spe <- loadHDF5SummarizedExperiment(fname)
#' 
#' spe <- getPatches(spe, cell_type_column="tier1",
#'                   response_cell_type="Epithelial",
#'                   explanatory_cell_type="Immune",
#'                   mmxpixel=0.000109,
#'                   npatches=10, n_iters=5)
#' de_results <- patchDE(spe, "patches_Epithelial_Immune")
#'
#' @return A list of patch-level differential expression results. Each element
#'   corresponds to a model term and contains a nested list with `pvals`,
#'   `ests`, and `ses` matrices indexed by gene and patch.
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom SummarizedExperiment colData "colData<-"
#' @rdname patchDE
#' 
#'
#' @export

patchDE <- function(spe, patchesID, method = "hasty", verbose = FALSE) {

  if (!patchesID %in% reducedDimNames(spe)) {
    msg <- "Patches ID '{patchesID}' not found in reducedDimNames(spe). Run getPatches() first."
    cli_abort(msg)
  }

  if (!"logcounts" %in% assayNames(spe)) {
    msg <- "Assay 'logcounts' not found in spe. Ensure normalized log-transformed counts are available."
    cli_abort(msg)
  }
   
  patch_df <- reducedDim(spe, patchesID)
  use <- which(as.vector(patch_df$response) == TRUE &
               !is.na(patch_df$X) &
               !is.na(as.vector(patch_df$patch)))
  y <- t(assay(spe, "logcounts")[, use, drop = FALSE])
  df <- patch_df[use,"X",drop=FALSE]
  patches <- as.vector(patch_df[use,"patch"])

  # get DE results per patch:
  results <- list()
  patch_names <- unique(patches)

  for (patchid in patch_names) {
    # Print some patch-level statistics for debugging
    if(verbose) {
      cli::cli_inform("Processing patch: {patchid}")
      cli::cli_inform("Number of cells in patch: {sum(patches == patchid, na.rm = TRUE)}")
      cli::cli_inform("Number of predictors in df: {ncol(df)}")
    }
    patchinds <- which(patches == patchid & !is.na(patches))
    switch(method,
           hasty = {
             results[[patchid]] <- hastyDE(y = y[patchinds ,, drop = FALSE],
                                          df = df[patchinds, , drop = FALSE])
           },
           limma = {
             results[[patchid]] <- limmaDE(y = y[patchinds ,, drop = FALSE],
                                          df = df[patchinds, , drop = FALSE])
           },
           {
             msg <- "Unsupported method '{method}'. Choose 'hasty' or 'limma'."
             cli_abort(msg)
           }
    )
  }
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
