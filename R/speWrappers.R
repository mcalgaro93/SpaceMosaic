#' Assign cells to patches using elliptical Gaussian assignments
#'
#' Iterative EM-like algorithm (per iteration):
#'   (a) Estimate per-patch centroid + covariance from xy.
#'   (b) Assign cells using spatial fit + X-diversity boost + Z-penalty + hunger.
#'   (c) Re-estimate per-patch centroid + covariance.
#'   (d) Assign cells using spatial fit only (no X/Z/hunger).
#'   (e) Contiguity check: set orphans to NA.
#'
#' Uses candidate filtering (only evaluates nearby patches per cell) and
#' mutual k-NN contiguity for scalability to 100k+ cells.
#'
#' @param spe SpatialExperiment object. Must have rownames.
#' @param X Design variables for each spatial unit. Either a character vector of
#'   column names in `colData(spe)`, or a numeric matrix/vector aligned to the
#'   rows of `spatialCoords(spe)` (one row per spatial unit). Each column is
#'   scaled to unit SD.
#' @param npatches Number of patches to create.
#' @param Z Optional per-cell context embeddings (cells x features). Either a
#'   single character string naming an entry in `reducedDims(spe)`, or a numeric
#'   matrix aligned to the rows of `spatialCoords(spe)`. If supplied, patches
#'   will prefer Z-coherent regions. NULL disables.
#' @param alpha Weight of the Z penalty. Typical range 0.2--1.0; default 0.5.
#'   Higher values force patches to respect microenvironment boundaries at the
#'   cost of spatial compactness. Set to 0 to ignore Z entirely.
#' @param beta Weight of the X diversity boost. Typical range 0.5/K--3/K where
#'   K = ncol(X); default 1 (appropriate for single-column X). For multi-column
#'   X, scale down proportionally (e.g. beta = 0.2 for K = 5). Raise if patches
#'   are too homogeneous; lower if patches are spatially fragmented.
#' @param hunger_weight Controls how aggressively low-variance patches grab cells.
#'   Typical range 0.3--0.7; default 0.5. 0 = all patches equally hungry
#'   (uniform sizes), 1 = hunger proportional to 1/totvar (maximizes variance
#'   reduction but allows extreme size imbalance). Lower if some patches shrink
#'   to nothing; raise if variance reduction is insufficient.
#' @param max_elongation Maximum ratio of largest to smallest eigenvalue of a
#'   patch covariance matrix. Typical range 2--8; default 4. Lower values force
#'   rounder patches; raise if tissue structures are genuinely elongated.
#' @param max_radius Maximum Euclidean distance from a patch centroid for
#'   assignment. Cells beyond this get zero spatial score. Default NULL (auto:
#'   3x the expected patch radius assuming uniform circular patches). Override
#'   if patches span very different density regions.
#' @param mahal_radius Maximum Mahalanobis radius for assignment. Cells beyond
#'   this (in each patch's own coordinate system) get zero spatial score.
#'   Typical range 2--4; default 3. Lower values make tighter patches with more
#'   unassigned cells; Inf disables the cutoff entirely.
#' @param n_candidates Number of nearest patch centroids to evaluate per cell.
#'   Typical range 10--50; default 20. Higher values are more accurate but
#'   slower. Rarely needs tuning unless npatches is very large (>1000).
#' @param n_iters Number of outer iterations. Typical range 10--30; default 15.
#'   Convergence is usually reached by 10--15; raise if patches are still
#'   shifting at the final iteration (check membership_log).
#' @param init_method Initialization method for patch seeds. "kmeans" uses
#'   isotropic k-means seeding (default). "gradient_ellipse" orients initial
#'   ellipses along the local spatial gradient of X. The gradient_ellipse
#'   option currently requires ncol(X) = 1.
#' @param init_gradient_k Number of nearest neighbors used to estimate local
#'   spatial gradients of X when init_method = "gradient_ellipse". Must be a
#'   finite integer greater than or equal to 3.
#' @param init_gradient_elongation Target initial ellipse elongation ratio
#'   (major/minor eigenvalue ratio) for init_method = "gradient_ellipse". Must
#'   be finite and greater than or equal to 1.
#' @param x_weighted_ellipse_second_pass Logical; if TRUE, the ellipse re-fit in
#'   step (c) up-weights cells whose X is farther from their patch mean. This
#'   can encourage elongated patches along smooth X gradients while keeping step
#'   (d) spatial-only. Applied only when ncol(X) = 1.
#' @param x_ellipse_gamma Strength of X-based up-weighting in second-pass
#'   ellipse fitting. 0 disables weighting; typical range 0.5--1.5.
#' @param x_ellipse_wmax Cap on standardized X-deviation used in the weighting
#'   rule to limit outlier influence. Typical range 2--4.
#' @param log_iters If TRUE, return a list with patch assignments plus
#'   per-iteration diagnostics (SS per patch and membership). Default TRUE.
#' @param patch_column Column name in colData(spe) where to store patch assignments.
#'   Also used to namespace `membership_log`/`ss_log` under `metadata(spe)`.
#'   Default = 'patch'.
#' @param verbose Show progress. Default TRUE.
#' @return The input `spe`, with results attached:
#'   \itemize{
#'     \item `colData(spe)[[patch_column]]`: named vector/factor of final patch assignments.
#'     \item `metadata(spe)[[patch_column]]`: (only if `log_iters = TRUE`) a list with
#'       `membership_log` (n_cells x n_iters matrix of patch assignments per
#'       iteration) and `ss_log` (npatches x n_iters matrix of per-patch
#'       sum-of-squares of X across iterations).
#'   }
#' @export
spe_getPatches <- function(spe, X, npatches,
                            Z = NULL,
                            alpha = 0.5,
                            beta = 1,
                            hunger_weight = 0.5,
                            max_elongation = 4,
                            max_radius = NULL,
                            mahal_radius = 3,
                            n_candidates = 20,
                            n_iters = 15,
                            init_method = c("kmeans", "gradient_ellipse"),
                            init_gradient_k = 30,
                            init_gradient_elongation = 4,
                            x_weighted_ellipse_second_pass = FALSE,
                            x_ellipse_gamma = 1,
                            x_ellipse_wmax = 3,
                            log_iters = TRUE,
                            patch_column = "patch",
                            verbose = TRUE) {

  xy <- SpatialExperiment::spatialCoords(spe)

  if (missing(X) || is.null(X) || length(X) == 0) {
    stop("`X` must be a non-empty character vector of colData names, or a matrix.")
  }

  X_mat <- .resolve_feature_matrix(spe, X, "X", source = "colData")

  Z_mat <- NULL
  if (!is.null(Z)) {
    Z_mat <- .resolve_feature_matrix(spe, Z, "Z", source = "reducedDim")
  }

  init_method <- match.arg(init_method)

  patchResult <- getPatches(xy, X_mat, npatches,
             Z_mat,
             alpha = alpha,
             beta = beta,
             hunger_weight = hunger_weight,
             max_elongation = max_elongation,
             max_radius = max_radius,
             mahal_radius = mahal_radius,
             n_candidates = n_candidates,
             n_iters = n_iters,
             init_method = init_method,
             init_gradient_k = init_gradient_k,
             init_gradient_elongation = init_gradient_elongation,
             x_weighted_ellipse_second_pass = x_weighted_ellipse_second_pass,
             x_ellipse_gamma = x_ellipse_gamma,
             x_ellipse_wmax = x_ellipse_wmax,
             log_iters = log_iters,
             verbose = verbose)

    if (log_iters) {
        patch_vec      <- patchResult$patch
        membership_log <- patchResult$membership_log[colnames(spe), , drop = FALSE]
        ss_log         <- patchResult$ss_log
    } else {
        patch_vec      <- patchResult
        membership_log <- NULL
        ss_log         <- NULL
    }

    patch_vec <- patch_vec[colnames(spe)]
    SummarizedExperiment::colData(spe)[[patch_column]] <- patch_vec

    if (log_iters) {
        S4Vectors::metadata(spe)[[patch_column]] <- list(
        membership_log = membership_log,
        ss_log = ss_log
        )
    }

    spe
}


# Internal helper: resolve either a matrix passed directly, or column/reducedDim
# names to pull from colData(spe) / reducedDim(spe, name).
.resolve_feature_matrix <- function(spe, value, value_name_arg, source = c("colData", "reducedDim")) {
  source <- match.arg(source)

  if (is.matrix(value)) {
    # user passed a matrix directly — use as-is, just sanity check dims
    if (nrow(value) != ncol(spe)) {
      stop("`", value_name_arg, "` matrix must have one row per column of `spe` ",
           "(", ncol(spe), " expected, got ", nrow(value), ").")
    }
    return(value)
  }

  if (is.character(value)) {
    if (source == "colData") {
      cd <- SummarizedExperiment::colData(spe)
      missing_cols <- setdiff(value, colnames(cd))
      if (length(missing_cols) > 0) {
        stop("The following `", value_name_arg, "` are not columns of colData(spe): ",
             paste(missing_cols, collapse = ", "))
      }
      return(as.matrix(cd[, value, drop = FALSE]))
    } else {
      if (length(value) != 1) {
        stop("`", value_name_arg, "` must be a single character string when naming a reducedDim.")
      }
      available_reddims <- SingleCellExperiment::reducedDimNames(spe)
      if (!value %in% available_reddims) {
        stop("`", value_name_arg, "` = '", value, "' not found in reducedDimNames(spe). ",
             "Available: ", paste(available_reddims, collapse = ", "))
      }
      return(SingleCellExperiment::reducedDim(spe, value))
    }
  }

  stop("`", value_name_arg, "` must be either a character vector/string of names, ",
       "or a matrix with ", ncol(spe), " rows (one per spatial unit).")
}


#' Embed cellular neighborhoods from single cell embeddings and positions
#'
#' Creates a neighborhood embedding by averaging a cell embedding matrix over
#' spatial neighbor networks at multiple scales, and stores the result back
#' onto the object as a new reduced dimension.
#'
#' @param spe A SpatialExperiment object.
#' @param embedding Either a matrix of single cell embeddings (cells x
#'   features, with one row per column of \code{spe}), or a single character
#'   string giving the name of a reduced dimension already stored in
#'   \code{reducedDim(spe)}.
#' @param ks Vector giving the number of nearest neighbors for each scale.
#'   Default \code{c(5, 50)}.
#' @param tissue Optional vector giving tissue IDs to prevent cross-tissue
#'   neighbor edges. Default NULL.
#' @param name Character string giving the name under which the result is
#'   stored via \code{reducedDim(spe, name)}. Default
#'   \code{"Z"}.
#' @return \code{spe} with a new reduced dimension, named according to
#'   \code{name}, added via \code{reducedDim(spe, name)}. This matrix has
#'   dimensions n cells x (ncol(embedding_mat) * length(ks)).
#' @export
spe_embedCellNeighborhoods <- function(spe, embedding, ks = c(5, 50), tissue = NULL,
                                       name = "Z") {
    embedding_mat <- .resolve_feature_matrix(spe, embedding, 'embedding', source = "reducedDim")
    reducedDim(spe, name) <- embedCellNeighborhoods(embedding_mat, spatialCoords(spe), ks, tissue)
    spe
}




#' spe_patchDE: run DE over all patches
#'
#' Runs differential expression across spatial patches, using a chosen assay
#' from a SpatialExperiment object and covariates resolved from its colData.
#'
#' @param spe A SpatialExperiment object.
#' @param df Either a character vector naming column(s) of \code{colData(spe)}
#'   to use as DE predictors, or a matrix with one row per column of
#'   \code{spe} (cells), used as-is.
#' @param assay_name Character string giving the name of the assay in
#'   \code{spe} to use as the expression matrix (cells x genes after
#'   transposition). Default "logcounts".
#' @param patch_column Character string giving the column of
#'   \code{colData(spe)} that holds patch IDs for each cell. Default "patch".
#' @param pearson Logical; if TRUE, transform y to Pearson residuals before DE
#' @param tot Numeric vector of total counts per cell (required if pearson = TRUE)
#' @param resid_mse Logical; if TRUE, include per-gene residual MSE in output
#' @param verbose Show progress. Default TRUE.
#' @return A list keyed by model variable. Each element contains `pvals`,
#'   `ests`, and `ses` matrices with genes in rows and patches in columns.
#'   If `resid_mse = TRUE`, the list also contains a `resid_mse` matrix with
#'   the same orientation.
#' @export
spe_patchDE <- function(spe, df, assay_name = "logcounts", patch_column = "patch", pearson = FALSE, tot = NULL, resid_mse = FALSE, verbose = TRUE){
        y <- t(assay(spe,assay_name))

        df <- as.data.frame(.resolve_feature_matrix(spe, df, 'df', source = 'colData'))

        patchDE(y, df, colData(spe)[,patch_column], pearson = pearson, tot = tot, resid_mse = resid_mse, verbose = verbose)
}


#' Compute patch attributes matrix W from a SpatialExperiment
#'
#' Extracts a reducedDim matrix and a colData patch assignment column
#' from a SpatialExperiment, then computes
#' per-patch mean attributes via `getPatchAttributes()`.
#'
#' @param spe A SpatialExperiment object.
#' @param dimred Character or integer scalar specifying which entry of
#'   `reducedDims(spe)` to use as Z. Default `"Z"`.
#' @param patch_col Character scalar naming the column of `colData(spe)`
#'   containing patch assignments. Default `"patch"`.
#' @return Matrix (npatches x features) of per-patch mean attributes.
#' @export
spe_getPatchAttributes <- function(spe, dimred = "Z", patch_col = "patch") {
  if (!dimred %in% SingleCellExperiment::reducedDimNames(spe)) {
    stop(sprintf("'%s' not found in reducedDims(spe).", dimred))
  }
  if (!patch_col %in% colnames(SummarizedExperiment::colData(spe))) {
    stop(sprintf("'%s' not found in colData(spe).", patch_col))
  }

  Z <- SingleCellExperiment::reducedDim(spe, dimred)
  patch <- SummarizedExperiment::colData(spe)[[patch_col]]

  getPatchAttributes(Z, patch)
}