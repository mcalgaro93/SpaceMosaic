#' Embed cellular neighborhoods in a SpatialExperiment
#'
#' Creates a neighborhood embedding for a \code{SpatialExperiment} object by
#' calling \code{\link{embedCellNeighborhoods}} on a stored reduced dimension
#' and the object's spatial coordinates, then storing the result as a new
#' reduced dimension.
#'
#' @param spe A SpatialExperiment object.
#' @param embedding A single character string giving the name of a reduced
#'   dimension stored in \code{reducedDim(spe)} with a 
#'   single cell embeddings matrix (cells x features).
#' @param ks Vector giving the number of nearest neighbors for each scale.
#'   Default \code{c(5, 50)}.
#' @param tissue Optional character string giving the name of a column in
#'   \code{colData(spe)} containing tissue IDs to prevent cross-tissue
#'   neighbor edges. Default NULL.
#' @param name Character string giving the name under which the result is
#'   stored via \code{reducedDim(spe, name)}. Default \code{"Z"}.
#' @return \code{spe} with a new reduced dimension, named according to
#'   \code{name}, added via \code{reducedDim(spe, name)}. This matrix has
#'   dimensions n cells x (ncol(embedding_mat) * length(ks)).
#' @seealso \code{\link{embedCellNeighborhoods}}, which does the multi-scale
#'   neighborhood averaging on a plain matrix.
#'
#' @importFrom SingleCellExperiment reducedDim reducedDim<- reducedDimNames
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
#' @examples 
#' library(SpatialExperiment)
#' spe <- readRDS(system.file("extdata", "cosmx_carcinoma.rds", package = "SpaceMosaic"))
#' spe <- embedCellNeighborhoods.spe(spe, embedding = "PCA", ks = c(5, 50), tissue = 'sample_id')
#' reducedDim(spe, "Z")[1:3, 1:4]
#' @export
embedCellNeighborhoods.spe <- function(spe, embedding, ks = c(5, 50), tissue = NULL,
                                       name = "Z") {

    if (!methods::is(spe, "SpatialExperiment")) {
      stop("`spe` must be a SpatialExperiment object.")
    }

    if (!embedding %in% reducedDimNames(spe)) {
      stop("`embedding` = '", embedding, "' not found in reducedDimNames(spe). ",
           "Available: ", paste(reducedDimNames(spe), collapse = ", "))
    }
    if (!is.null(tissue) && !tissue %in% colnames(colData(spe))) {
      stop("`tissue` = '", tissue, "' not found in colData(spe). ",
           "Available: ", paste(colnames(colData(spe)), collapse = ", "))
    }
    embedding_mat <- SingleCellExperiment::reducedDim(spe, embedding)
    reducedDim(spe, name) <- embedCellNeighborhoods(embedding_mat, spatialCoords(spe), ks, tissue = colData(spe)[[tissue]])
    spe
}


#' Run the SpaceMosaic patching workflow on a `SpatialExperiment`
#'
#' Convenience wrapper for the SpaceMosaic patching workflow on a
#' `SpatialExperiment`. This wrapper contains three steps:
#' `getPatches()` identifies spatial patches, `getPatchDiagnostics()` computes
#' patch-level information from the resulting patch assignments, and
#' `getPatchPolys()` constructs patch polygons using the patch assignments and
#' diagnostic information. Results are stored directly in the
#' `SpatialExperiment` object.
#'
#' Spatial coordinates are taken from `spatialCoords(spe)`, design variables
#' from `colData(spe)`, and (optionally) context embeddings from
#' `reducedDims(spe)`. All results are returned attached to `spe`.
#'
#' @param spe A `SpatialExperiment` object. Must have `colnames` (cell
#'   identifiers), non-empty `spatialCoords()`, and the columns named in `X`
#'   present in `colData(spe)`.
#' @param X Design variables for each spatial unit. A character vector of
#'   column names in `colData(spe)`. Each column is scaled to unit SD by
#'   `getPatches()`.
#' @param npatches Number of patches to create.
#' @param Z Character string or `NULL`. Name of an entry in
#'   `reducedDimNames(spe)` holding per-cell context embeddings
#'   (cells x features). If supplied, patches will prefer Z-coherent regions.
#'   Default `NULL` disables this behaviour.
#' @param alpha Weight of the Z penalty. Typical range 0.2--1.0; default 0.5.
#'   Higher values force patches to respect microenvironment boundaries at the
#'   cost of spatial compactness. Set to 0 to ignore Z entirely.
#' @param beta Weight of the X diversity boost. Typical range 0.5/K--3/K where
#'   K = `length(X)`; default 1 (appropriate for a single design variable). For
#'   multiple design variables, scale down proportionally (e.g. `beta = 0.2`
#'   for K = 5). Raise if patches are too homogeneous; lower if patches are
#'   spatially fragmented.
#' @param hunger_weight Controls how aggressively low-variance patches grab
#'   cells. Typical range 0.3--0.7; default 0.5. 0 = all patches equally hungry
#'   (uniform sizes), 1 = hunger proportional to 1/totvar (maximizes variance
#'   reduction but allows extreme size imbalance).
#' @param max_elongation Maximum ratio of largest to smallest eigenvalue of a
#'   patch covariance matrix. Typical range 2--8; default 4. Lower values force
#'   rounder patches.
#' @param max_radius Maximum Euclidean distance from a patch centroid for
#'   assignment. Cells beyond this get zero spatial score. Default `NULL`
#'   (auto: 3x the expected patch radius assuming uniform circular patches).
#' @param mahal_radius Maximum Mahalanobis radius for assignment. Cells beyond
#'   this (in each patch's own coordinate system) get zero spatial score.
#'   Typical range 2--4; default 3. `Inf` disables the cutoff entirely.
#' @param n_candidates Number of nearest patch centroids to evaluate per cell.
#'   Typical range 10--50; default 20. Higher values are more accurate but
#'   slower.
#' @param n_iters Number of outer iterations. Typical range 10--30; default 15.
#' @param init_method Initialization method for patch seeds. `"kmeans"`
#'   (default) uses isotropic k-means seeding; `"gradient_ellipse"` orients
#'   initial ellipses along the local spatial gradient of X and currently
#'   requires `length(X) == 1`.
#' @param init_gradient_k Number of nearest neighbors used to estimate local
#'   spatial gradients of X when `init_method = "gradient_ellipse"`. A finite
#'   integer >= 3; default 30.
#' @param init_gradient_elongation Target initial ellipse elongation ratio
#'   (major/minor eigenvalue ratio) for `init_method = "gradient_ellipse"`.
#'   Finite and >= 1; default 4.
#' @param x_weighted_ellipse_second_pass Logical; if `TRUE`, the ellipse re-fit
#'   up-weights cells whose X is farther from their patch mean, which can
#'   encourage elongated patches along smooth X gradients. Applied only when
#'   `length(X) == 1`. Default `FALSE`.
#' @param x_ellipse_gamma Strength of X-based up-weighting in second-pass
#'   ellipse fitting. 0 disables weighting; typical range 0.5--1.5; default 1.
#' @param x_ellipse_wmax Cap on standardized X-deviation used in the weighting
#'   rule, to limit outlier influence. Typical range 2--4; default 3.
#' @param log_iters Logical; if `TRUE` (default), per-iteration diagnostics are
#'   retained and stored in
#'   `metadata(spe)$SpaceMosaic$patch_iterations`.
#' @param patch_diagnostics Logical; if `TRUE` (default), compute patch-level
#'   diagnostics with `getPatchDiagnostics()` and merge its elements
#'   (`patch_diagnostics`, `assignment_summary`, `connectivity_curve`) into
#'   `metadata(spe)$SpaceMosaic`.
#' @param k Maximum number of spatial neighbors used by `getPatchDiagnostics()`.
#'   Default 10. Limited internally to `n - 1` for small datasets. Ignored when
#'   `patch_diagnostics = FALSE`.
#' @param strict_k Neighbors used for `strict_component_fraction` in
#'   `getPatchDiagnostics()`. Must be no greater than `k`; `NULL` (default)
#'   uses `min(5, k)`. Ignored when `patch_diagnostics = FALSE`.
#' @param patch_polys Logical; if `TRUE` (default), compute patch polygons with
#'   `getPatchPolys()` and store the result in
#'   `metadata(spe)$SpaceMosaic$patch_polys`. When `patch_diagnostics = FALSE`,
#'   polygons are built without diagnostic information (`patch_data = NULL`).
#' @param patch_column Column name in `colData(spe)` where patch assignments are
#'   stored. Default `"patch"`. An existing column of the same name is
#'   overwritten.
#' @param verbose Logical; show progress. Default `TRUE`.
#'
#' @return The input `spe`, with results attached:
#'   \itemize{
#'     \item `colData(spe)[[patch_column]]`: named vector/factor of final patch
#'       assignments, ordered to match `colnames(spe)`. `NA` denotes an
#'       unassigned cell.
#'     \item `metadata(spe)$SpaceMosaic$patch_iterations`: if
#'       `log_iters = TRUE`, a list with `patch` (final assignments),
#'       `membership_log` (n_cells x n_iters matrix of patch assignments per
#'       iteration) and `ss_log` (npatches x n_iters matrix of per-patch
#'       sum-of-squares of X across iterations).
#'     \item `metadata(spe)$SpaceMosaic$patch_diagnostics`,
#'       `$assignment_summary` and `$connectivity_curve`: the three data frames
#'       returned by `getPatchDiagnostics()`, if `patch_diagnostics = TRUE`.
#'     \item `metadata(spe)$SpaceMosaic$patch_polys`: patch polygons, if
#'       `patch_polys = TRUE`.
#'   }
#'
#' @seealso [getPatches()], [getPatchDiagnostics()], [getPatchPolys()]
#'
#' @importFrom SingleCellExperiment reducedDimNames
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom S4Vectors metadata metadata<-
#'
#' @examples
#' library(SpatialExperiment)
#' spe <- readRDS(system.file("extdata", "cosmx_carcinoma.rds", package = "SpaceMosaic"))
#' spe <- getPatches.spe(
#'   spe = spe,
#'   X = "distance_to_fibroblast",
#'   npatches = 50
#' )
#' head(spe$patch)
#' head(metadata(spe)$SpaceMosaic$patch_diagnostics)
#'
#' @export getPatches.spe
getPatches.spe <- function(spe, X, npatches,
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
                            patch_diagnostics = TRUE,
                            k = 10L, strict_k = NULL,
                            patch_polys = TRUE,
                            patch_column = "patch",
                            verbose = TRUE) {

  if (!methods::is(spe, "SpatialExperiment")) {
    stop("`spe` must be a SpatialExperiment object.")
  }

  xy <- SpatialExperiment::spatialCoords(spe)
  if (is.null(rownames(xy))) {
    rownames(xy) <- colnames(spe)
  }

  if (missing(X) || is.null(X) || length(X) == 0) {
    stop("`X` must be a non-empty character vector of colData names.")
  }

  missing_cols <- setdiff(X, colnames(SummarizedExperiment::colData(spe)))
    if (length(missing_cols) > 0) {
      stop("The following `X` are not columns of colData(spe): ",
           paste(missing_cols, collapse = ", "))
    }

  X_mat <- as.matrix(SummarizedExperiment::colData(spe)[, X, drop = FALSE])

  Z_mat <- NULL
  if (!is.null(Z)) {
    if (!is.character(Z) || length(Z) != 1L || is.na(Z)) {
        stop("`Z` must be NULL or a single character string naming a reduced dimension.")
    }

    if (!Z %in% reducedDimNames(spe)) {
        stop(
        "`Z` = '", Z, "' not found in reducedDimNames(spe). ",
        "Available: ", paste(reducedDimNames(spe), collapse = ", ")
        )
    }

  Z_mat <- SingleCellExperiment::reducedDim(spe, Z)
}

  init_method <- match.arg(init_method)

  patchResult <- getPatches(xy = xy, 
             X = X_mat, 
             npatches = npatches,
             Z = Z_mat,
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
          patch_vec      <- patchResult$patch[colnames(spe)]
          membership_log <- patchResult$membership_log[colnames(spe), , drop = FALSE]
          ss_log         <- patchResult$ss_log
      } else {
          patch_vec <- patchResult[colnames(spe)]
      }

      SummarizedExperiment::colData(spe)[[patch_column]] <- patch_vec

      if (log_iters) {
          S4Vectors::metadata(spe)[["SpaceMosaic"]][["patch_iterations"]] <- list(
              patch = patch_vec,
              membership_log = membership_log,
              ss_log = ss_log
          )
          patch_metadata <- S4Vectors::metadata(spe)[["SpaceMosaic"]][["patch_iterations"]]
      } else {
          patch_metadata <- patch_vec
      }

    if(patch_diagnostics) {
      patch_diagnostics_list <- getPatchDiagnostics(
            xy = spatialCoords(spe),
            X = X_mat,
            patches = patch_metadata,
            k = k,
            strict_k = strict_k
        )
      metadata(spe)$SpaceMosaic[names(patch_diagnostics_list)] <- patch_diagnostics_list
      patch_data <- patch_diagnostics_list$patch_diagnostics
    } else{
      patch_data <- NULL
    }

    if (patch_polys) {
      patch_polys <- getPatchPolys(
        xy = spatialCoords(spe),
        patch = patch_vec,
        patch_data = patch_data
      )
      metadata(spe)$SpaceMosaic$patch_polys <- patch_polys
    }
    spe
}



#' Run patch-level differential expression on a SpatialExperiment
#'
#' Runs \code{patchDE()} on the expression data and patch assignments stored
#' in a \code{SpatialExperiment} object and returns the raw results alongside
#' the object they were computed from.
#'
#' @param spe A \code{SpatialExperiment} object. Must contain the requested
#'   expression assay and patch assignments in \code{colData(spe)}.
#' @param predictor_cols Character vector of column names in
#'   \code{colData(spe)} to use as predictors in \code{patchDE()}. Must be
#'   non-empty.
#' @param assay Character string naming the assay in \code{spe} that contains
#'   the expression matrix. Default \code{"logcounts"}.
#' @param patch_column Character string naming the column in
#'   \code{colData(spe)} that contains patch assignments. Cells with \code{NA}
#'   here are dropped by \code{patchDE()}. Default \code{"patch"}.
#' @param method Differential-expression backend. \code{"hasty"} uses ordinary
#'   least squares and \code{"limma"} uses empirical-Bayes moderated inference.
#'   Default \code{"hasty"} for backward compatibility.
#' @param pearson Logical; if TRUE, expression values are converted to Pearson
#'   residuals within each patch before the model is fitted. Requires
#'   \code{tot}. When \code{method = "limma"} this also disables the
#'   mean-variance trend. Default FALSE.
#' @param tot Optional character string naming a column in
#'   \code{colData(spe)} that holds per-cell total counts. Required when
#'   \code{pearson = TRUE}; ignored otherwise. Default NULL.
#' @param resid_mse Logical; passed to \code{patchDE()} to request per-gene
#'   residual mean squared errors, which are carried through in the \code{de}
#'   element. Default FALSE.
#' @param return_residuals Logical; if TRUE, the cells-by-genes residual matrix
#'   produced by \code{patchDE()} is transposed to genes-by-cells and stored
#'   as assay \code{residual_assay} in the returned \code{spe}. Default FALSE.
#' @param residual_assay Character string naming the assay under which
#'   residuals are stored in \code{spe} when \code{return_residuals = TRUE}.
#'   Matches the default \code{assay_name} expected by
#'   \code{\link{moranTest.spe}}. Default \code{"residuals"}.
#' @param verbose Logical; show a progress bar over patches. Default TRUE.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{de}}{The differential-expression results exactly as returned
#'       by \code{patchDE()}: one element per entry of \code{predictor_cols},
#'       each holding the \code{pvals}, \code{ests} and \code{ses} matrices
#'       (genes by patches), plus anything else the backend produced, such as
#'       residual mean squared errors when \code{resid_mse = TRUE}.}
#'     \item{\code{spe}}{The input \code{SpatialExperiment}. When
#'       \code{return_residuals = TRUE}, this includes a new
#'       \code{residual_assay} assay (genes by cells) holding the residual
#'       matrix from \code{patchDE()}; otherwise it is returned unchanged.}
#'   }
#'
#' @details
#' The expression matrix is extracted from \code{assay(spe, assay)} and
#' transposed from genes-by-cells to cells-by-genes before being passed to
#' \code{patchDE()}. Predictor variables are taken from \code{colData(spe)},
#' and patch assignments from \code{colData(spe)[[patch_column]]}.
#'
#' Differential expression is run on every non-missing patch in \code{spe},
#' and the patch columns of the returned matrices follow the order
#' \code{patchDE()} produced.
#'
#' \code{patchDE()} returns one residual matrix per fit, not one per
#' predictor, so cells with a missing patch assignment (and therefore no fit)
#' get \code{NA} residuals in \code{residual_assay}.
#'
#' @seealso
#' \code{\link{patchDE}},
#' \code{\link{patchMetaAnalysis.spe}},
#' \code{\link{moranTest.spe}}
#'
#' @examples
#' 
#' library(SpatialExperiment)
#' spe <- readRDS(system.file("extdata", "cosmx_carcinoma.rds", package = "SpaceMosaic"))
#' spe_use <- spe[, spe$celltype == "Cancer.cells"]
#' spe_use <- getPatches.spe(
#'   spe = spe_use,
#'   X = "distance_to_fibroblast",
#'   npatches = 50
#' )
#' de_result <- patchDE.spe(spe_use,
#'   predictor_cols = "distance_to_fibroblast"
#' )
#' names(de_result)
#' head(de_result$de$distance_to_fibroblast$pvals)
#' 
#' @export


patchDE.spe <- function(
    spe,
    predictor_cols,
    assay = "logcounts",
    patch_column = "patch",
    method = c("hasty", "limma"),
    pearson = FALSE,
    tot = NULL,
    resid_mse = FALSE,
    return_residuals = FALSE,
    residual_assay = "residuals",
    verbose = TRUE
) {

    method <- match.arg(method)

    if (!methods::is(spe, "SpatialExperiment")) {
        stop("`spe` must be a SpatialExperiment object.")
    }

    if (missing(predictor_cols) ||
        is.null(predictor_cols) ||
        length(predictor_cols) == 0) {
        stop(
            "`predictor_cols` must be a non-empty character vector ",
            "of colData names."
        )
    }

    missing_cols <- setdiff(
        predictor_cols,
        colnames(SummarizedExperiment::colData(spe))
    )

    if (length(missing_cols) > 0) {
        stop(
            "The following `predictor_cols` are not columns of colData(spe): ",
            paste(missing_cols, collapse = ", ")
        )
    }

    if (!patch_column %in%
        colnames(SummarizedExperiment::colData(spe))) {
        stop(sprintf(
            "'%s' not found in colData(spe).",
            patch_column
        ))
    }

    if (!assay %in% SummarizedExperiment::assayNames(spe)) {
        stop(sprintf(
            "'%s' not found in assay names of `spe`. Available assays: %s",
            assay, paste(SummarizedExperiment::assayNames(spe), collapse = ", ")
        ))
    }

    if (pearson && is.null(tot)) {
        stop("`tot` must be supplied when `pearson = TRUE`.")
    }

    # Expression matrix: genes x cells -> cells x genes
    y <- t(SummarizedExperiment::assay(spe, assay))

    df <- SummarizedExperiment::colData(
        spe
    )[, predictor_cols, drop = FALSE]

    if (pearson && !is.null(tot)) {

        if (!tot %in%
            colnames(SummarizedExperiment::colData(spe))) {
            stop(sprintf(
                "'%s' not found in colData(spe).",
                tot
            ))
        }

        tot <- SummarizedExperiment::colData(spe)[[tot]]
    }

    # Run patchDE
    de_res <- patchDE(
        y,
        df,
        SummarizedExperiment::colData(spe)[[patch_column]],
        method = method,
        pearson = pearson,
        tot = tot,
        resid_mse = resid_mse,
        return_residuals = return_residuals,
        verbose = verbose
    )

    # patchDE() wraps its output when residuals are requested
    if (return_residuals) {
        SummarizedExperiment::assay(spe, residual_assay) <- t(de_res$residuals)
        out <- list(de = de_res$de)
    } else {
        out <- list(de = de_res)
    }

    out$spe <- spe

    out

}

#' Meta-analyse patch-level differential expression across an embedding
#'
#' Takes the output of \code{patchDE.spe()}, derives patch-level attributes
#' from a reduced-dimension embedding, and runs \code{patchMetaAnalysis()} to
#' combine the per-patch fits.
#'
#' @param patchDE_result The list returned by \code{patchDE.spe()}. Must
#'   contain a \code{de} element holding the per-predictor differential
#'   expression results and a \code{spe} element holding the
#'   \code{SpatialExperiment} they were computed from.
#' @param embedding_name Character string naming the entry of
#'   \code{reducedDims(spe)} used to derive patch attributes. Default
#'   \code{"Z"}.
#' @param patch_column Character string naming the column in
#'   \code{colData(spe)} that contains patch assignments. Default
#'   \code{"patch"}.
#' @param summarize_subgroups Logical; if TRUE, \code{summarizeSubgroups()} is
#'   run on the meta-analysis output and the result is added to the returned
#'   list. Default FALSE.
#' @param cellmeta_cols Optional character vector of \code{colData(spe)}
#'   columns passed to \code{summarizeSubgroups()} as cell metadata. Used only
#'   when \code{summarize_subgroups = TRUE}. Default NULL.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{meta}}{The meta-analysis results as returned by
#'       \code{patchMetaAnalysis()}: one element per predictor, each holding
#'       \code{pvals}, \code{ests} and \code{ses}, and \code{subgroups} where
#'       the backend produced it.}
#'     \item{\code{W}}{The patches-by-attributes matrix from
#'       \code{getPatchAttributes()}, with rows in the same order as the patch
#'       columns of the differential expression results.}
#'     \item{\code{subgroups_summary}}{The output of
#'       \code{summarizeSubgroups()}. Present only when
#'       \code{summarize_subgroups = TRUE}.}
#'   }
#'
#' @details
#' Patch attributes are computed from the full embedding and then reordered to
#' match the patch columns of the differential expression matrices, so the two
#' line up when passed to \code{patchMetaAnalysis()}. The patches represented
#' in \code{colData(spe)[[patch_column]]} and those in the differential
#' expression results must be the same set; their order need not agree.
#'
#' @seealso
#' \code{\link{patchDE.spe}},
#' \code{\link{patchMetaAnalysis}},
#' \code{\link{getPatchAttributes}},
#' \code{\link{summarizeSubgroups}}
#' 
#' @examples
#' library(SpatialExperiment)
#' spe <- readRDS(system.file("extdata", "cosmx_carcinoma.rds", package = "SpaceMosaic"))
#' spe <- embedCellNeighborhoods.spe(spe, embedding = "PCA", ks = c(5, 50), tissue = 'sample_id')
#' spe_use <- spe[,spe$celltype == "Cancer.cells"]
#' spe_use <- getPatches.spe(
#'   spe = spe_use,
#'   X = "distance_to_fibroblast",
#'   npatches = 50
#' )
#' de_result <- patchDE.spe(spe_use,
#'   predictor_cols = "distance_to_fibroblast"
#' )
#' meta_result <- patchMetaAnalysis.spe(de_result, embedding_name = "Z")
#' head(meta_result$meta$distance_to_fibroblast$pvals)
#'
#' @export

patchMetaAnalysis.spe <- function(
    patchDE_result,
    embedding_name = "Z",
    patch_column = "patch",
    summarize_subgroups = FALSE,
    cellmeta_cols = NULL
) {

    if (!is.list(patchDE_result) ||
        !all(c("de", "spe") %in% names(patchDE_result))) {
        stop(
            "`patchDE_result` must be a list with `de` and `spe` elements, ",
            "as returned by `patchDE.spe()`."
        )
    }

    spe <- patchDE_result$spe
    de_res <- patchDE_result$de

    if (length(de_res) == 0 || is.null(names(de_res))) {
        stop(
            "`patchDE_result$de` must be a named list with one element ",
            "per predictor."
        )
    }

    if (!patch_column %in%
        colnames(SummarizedExperiment::colData(spe))) {
        stop(sprintf(
            "'%s' not found in colData(spe).",
            patch_column
        ))
    }

    if (!embedding_name %in%
        SingleCellExperiment::reducedDimNames(spe)) {
        stop(sprintf(
            "Embedding '%s' not found in reducedDims(spe).",
            embedding_name
        ))
    }

    # Every predictor must carry the statistics the meta-analysis needs
    for (predictor in names(de_res)) {

        if (!all(
            c("pvals", "ests", "ses") %in%
            names(de_res[[predictor]])
        )) {
            stop(
                "Expected `patchDE_result$de[[predictor]]` to contain ",
                "`pvals`, `ests`, and `ses` for predictor '",
                predictor,
                "'."
            )
        }
    }

    patch <- SummarizedExperiment::colData(spe)[[patch_column]]

    # Extract embedding
    Z <- SingleCellExperiment::reducedDim(
        spe,
        embedding_name
    )

    # Calculate patch-level embedding attributes
    W <- getPatchAttributes(
        Z,
        patch
    )

    # Run meta-analysis
    de_res_meta <- patchMetaAnalysis(
        de_res,
        W
    )

    out <- list(
        meta = de_res_meta,
        W = W
    )

    if (summarize_subgroups) {

        if (!is.null(cellmeta_cols)) {

            missing_cols <- setdiff(
                cellmeta_cols,
                colnames(SummarizedExperiment::colData(spe))
            )

            if (length(missing_cols) > 0) {
                stop(
                    "The following `cellmeta_cols` are not columns of ",
                    "colData(spe): ",
                    paste(missing_cols, collapse = ", ")
                )
            }
        }

        out$subgroups_summary <- summarizeSubgroups(
            de_res_meta,
            patch = patch,
            cellmeta = SummarizedExperiment::colData(
                spe
            )[, cellmeta_cols, drop = FALSE]
        )
    }

    out
}

#'   moranTest() method for \code{SpatialExperiment} objects. Extracts
#'   residuals from the specified assay and spatial coordinates from
#'   \code{spatialCoords(spe)}, then dispatches to \code{moranTest}.
#'
#' @param spe A \code{SpatialExperiment} object.
#' @param assay_name Character; name of the assay in \code{spe} containing
#'   the values to test (e.g. Pearson residuals). Default \code{"residuals"}.
#' @param patch_column Name of the column in `colData(spe)` where patch assignments are
#'   stored. Default `"patch"`.
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
#' @return A data frame with one row per gene and patch. It contains patch and
#'   gene identifiers, number of cells, observed and expected Moran's I, raw and
#'   adjusted empirical p-values, graph and permutation settings, and `status`.
#'   A status of `"ok"` indicates a completed test; non-testable inputs such as
#'   constant or non-finite residuals return missing statistics with an
#'   explanatory status rather than stopping the remaining tests.
#'
#' @importFrom methods is
#' @importFrom SummarizedExperiment assayNames assay colData
#' @importFrom SpatialExperiment spatialCoords
#'
#' @examples
#' library(SpatialExperiment)
#' spe <- readRDS(system.file("extdata", "cosmx_carcinoma.rds", package = "SpaceMosaic"))
#' spe_use <- spe[,spe$celltype == "Cancer.cells"]
#' spe_use <- getPatches.spe(
#'   spe = spe_use,
#'   X = "distance_to_fibroblast",
#'   npatches = 50
#' )
#' de_result <- patchDE.spe(spe_use,
#'   predictor_cols = "distance_to_fibroblast",
#'   return_residuals = TRUE
#' )
#' moran_result <- moranTest.spe(de_result$spe[1:10,])
#' @export

moranTest.spe <- function(spe, assay_name = 'residuals' , patch_column = "patch", k = 10L,
                      n_permutations = 999L,
                      alternative = c("greater", "less", "two.sided"),
                      p_adjust_method = "BH",
                      adjustment_scope = c("global", "patch", "gene")){



                  if (!is(spe, "SpatialExperiment")) {
                    stop("'spe' must be a SpatialExperiment object.")
                  }

                  if (!is.character(assay_name) || length(assay_name) != 1) {
                    stop("'assay_name' must be a single character string.")
                  }

                  if (!assay_name %in% assayNames(spe)) {
                    stop(sprintf("'%s' not found in assay names of 'spe'. Available assays: %s",
                                assay_name, paste(assayNames(spe), collapse = ", ")))
                  }

                  if (!patch_column %in% colnames(colData(spe))) {
                    stop(sprintf("'%s' not found in colData(spe). Available: %s",
                                patch_column, paste(colnames(colData(spe)), collapse = ", ")))
                  }

                  res <- moranTest( residuals = t(assay(spe,assay_name)),
                                          xy = spatialCoords(spe),
                                          patch = colData(spe)[[patch_column]],
                                          k = k,
                                          n_permutations = n_permutations,
                                          alternative = alternative,
                                          p_adjust_method = p_adjust_method,
                                          adjustment_scope = adjustment_scope)
                  res
                      }