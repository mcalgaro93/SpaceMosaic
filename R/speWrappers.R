#' @describeIn getPatches Method for \code{SpatialExperiment} objects.
#'   Resolves `X` from `colData(spe)` and `Z` from `reducedDims(spe)`, runs
#'   `getPatches()` on `spatialCoords(spe)`, and attaches the result back onto
#'   `spe` rather than returning it standalone. See `getPatches` for the
#'   underlying algorithm and the meaning/defaults of all tuning parameters
#'   (`alpha`, `beta`, `hunger_weight`, `max_elongation`, `max_radius`,
#'   `mahal_radius`, `n_candidates`, `n_iters`, `init_method` and friends,
#'   `log_iters`, `verbose`).
#'
#' @param spe SpatialExperiment object. Must have rownames.
#' @param X Design variables for each spatial unit. Either a character vector of
#'   column names in `colData(spe)`, or a numeric matrix/vector aligned to the
#'   rows of `spatialCoords(spe)` (one row per spatial unit). Each column is
#'   scaled to unit SD.
#' @param Z Optional per-cell context embeddings (cells x features). Either a
#'   single character string naming an entry in `reducedDims(spe)`, or a numeric
#'   matrix aligned to the rows of `spatialCoords(spe)`. If supplied, patches
#'   will prefer Z-coherent regions. NULL disables.
#' @param patch_column Column name in colData(spe) where to store patch assignments.
#'   Default = 'patch'.
#' @return The input `spe`, with results attached:
#'   \itemize{
#'     \item `colData(spe)[[patch_column]]`: named vector/factor of final patch assignments.
#'     \item `metadata(spe)[['SpaceMosaic']][['patch_iterations']]`: (only if `log_iters = TRUE`) a list with
#'       `membership_log` (n_cells x n_iters matrix of patch assignments per
#'       iteration) and `ss_log` (npatches x n_iters matrix of per-patch
#'       sum-of-squares of X across iterations).
#'   }
#' @export
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
        S4Vectors::metadata(spe)[["SpaceMosaic"]][["patch_iterations"]] <- list(
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


#' @describeIn embedCellNeighborhoods Method for \code{SpatialExperiment}
#'   objects. Resolves `embedding` from `reducedDims(spe)` (or takes it as a
#'   matrix directly), runs `embedCellNeighborhoods()` against
#'   `spatialCoords(spe)`, and stores the result back as a new reduced
#'   dimension on `spe`. See `embedCellNeighborhoods` for how the
#'   multi-scale neighborhood averaging works.
#'
#' @param spe A SpatialExperiment object.
#' @param embedding Either a matrix of single cell embeddings (cells x
#'   features, with one row per column of \code{spe}), or a single character
#'   string giving the name of a reduced dimension already stored in
#'   \code{reducedDim(spe)}.
#' @param name Character string giving the name under which the result is
#'   stored via \code{reducedDim(spe, name)}. Default \code{"Z"}.
#' @return \code{spe} with a new reduced dimension, named according to
#'   \code{name}, added via \code{reducedDim(spe, name)}. This matrix has
#'   dimensions n cells x (ncol(embedding_mat) * length(ks)).
#' @export
embedCellNeighborhoods.spe <- function(spe, embedding, ks = c(5, 50), tissue = NULL,
                                       name = "Z") {
    embedding_mat <- .resolve_feature_matrix(spe, embedding, 'embedding', source = "reducedDim")
    reducedDim(spe, name) <- embedCellNeighborhoods(embedding_mat, spatialCoords(spe), ks, tissue)
    spe
}




#' @describeIn patchDE Method for \code{SpatialExperiment} objects. Extracts
#'   the expression matrix from the chosen assay and resolves `df` from
#'   `colData(spe)`, then dispatches to `patchDE()`. See `patchDE` for
#'   details of the DE model and the meaning of `pearson`, `tot`, and
#'   `resid_mse`.
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
#' @export
patchDE.spe <- function(spe, df, assay_name = "logcounts", patch_column = "patch", pearson = FALSE, tot = NULL, resid_mse = FALSE, verbose = TRUE){
        y <- t(assay(spe,assay_name))

        df <- as.data.frame(.resolve_feature_matrix(spe, df, 'df', source = 'colData'))

        patchDE(y, df, colData(spe)[,patch_column], pearson = pearson, tot = tot, resid_mse = resid_mse, verbose = verbose)

}


#' @describeIn getPatchAttributes Method for \code{SpatialExperiment} objects.
#'   Extracts a reducedDim matrix and a colData patch assignment column from
#'   `spe`, then dispatches to `getPatchAttributes()`. `spe`'s `dimred` and
#'   `patch_col` correspond to the `Z` and `patch` arguments of
#'   `getPatchAttributes()` (renamed here since they now identify columns
#'   rather than being passed as data directly); see `getPatchAttributes`
#'   for how the per-patch means are computed.
#'
#' @param spe A SpatialExperiment object.
#' @param dimred Character or integer scalar specifying which entry of
#'   `reducedDims(spe)` to use as Z. Default `"Z"`.
#' @param patch_col Character scalar naming the column of `colData(spe)`
#'   containing patch assignments. Default `"patch"`.
#' @export
getPatchAttributes.spe <- function(spe, dimred = "Z", patch_col = "patch") {
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

#' @describeIn moranTest Method for \code{SpatialExperiment} objects. Extracts
#'   residuals from the specified assay and spatial coordinates from
#'   \code{spatialCoords(spe)}, then dispatches to \code{moranTest}. See
#'   \code{moranTest} for details of the permutation test and the meaning of
#'   \code{k}, \code{n_permutations}, \code{alternative},
#'   \code{p_adjust_method}, and \code{adjustment_scope}.
#'
#' @param spe A \code{SpatialExperiment} object.
#' @param assay_name Character; name of the assay in \code{spe} containing
#'   the values to test (e.g. Pearson residuals). Default \code{"residuals"}.
#'
#' @export

moranTest.spe <- function(spe, assay_name = 'residuals' , patch = NULL, k = 10L,
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

                  res <- moranTest( residuals = t(assay(spe,assay_name)),
                                          xy = spatialCoords(spe),
                                          patch = patch,
                                          k = k,
                                          n_permutations = n_permutations,
                                          alternative = alternative,
                                          p_adjust_method = p_adjust_method,
                                          adjustment_scope = adjustment_scope)
                  res
                      }

#' @describeIn getPatchDiagnostics Method for \code{SpatialExperiment}
#'   objects. Resolves \code{X} from \code{colData(spe)}, pulls the current
#'   patch assignment from \code{colData(spe)$patch} (merging into
#'   \code{metadata(spe)$patch} if present), and dispatches to
#'   \code{getPatchDiagnostics} using \code{spatialCoords(spe)}. See
#'   \code{getPatchDiagnostics} for the diagnostics returned and the meaning
#'   of \code{k} and \code{strict_k}.
#'
#' @param spe SpatialExperiment object. Must have rownames.
#' @param X Design variables for each spatial unit. Either a character vector of
#'   column names in `colData(spe)`, or a numeric matrix/vector aligned to the
#'   rows of `spatialCoords(spe)` (one row per spatial unit). Each column is
#'   scaled to unit SD.
#'
#' @export

getPatchDiagnostics.spe <- function(spe, X, k = 10L, strict_k = NULL) {



                  if (!is(spe, "SpatialExperiment")) {
                    stop("'spe' must be a SpatialExperiment object.")
                  }

                  X_mat <- .resolve_feature_matrix(spe, X, "X", source = "colData")

                  patch_metadata <- metadata(spe)$SpaceMosaic$patch_iterations
                  if (is.null(patch_metadata)) {
                  patch <- colData(spe)$patch
                  } else {
                  patch_metadata$patch <- colData(spe)$patch
                  }

                  metadata(spe)$SpaceMosaic$patch_diagnostics <- getPatchDiagnostics(
                                          xy = spatialCoords(spe),
                                          X = X_mat,
                                          patch = patch_metadata,
                                          k = k,
                                          strict_k = strict_k)
                  spe
                      }

patchDEWorkflow <- function(spe, predictor_cols, assay = 'logcounts', patch_column = "patch",
                            embedding_name = 'Z', metaanalysis = TRUE, pearson = TRUE,
                            tot = NULL, resid_mse = FALSE, verbose = TRUE) {

    y <- t(assay(spe, assay))
    df <- as.data.frame(.resolve_feature_matrix(spe, predictor_cols, 'predictor_cols', source = 'colData'))

    de_res <- patchDE(y, df, colData(spe)[, patch_column], pearson = pearson, tot = tot,
                      resid_mse = resid_mse, verbose = verbose)

    patch_ids <- colData(spe)$patch
    patch_ids <- as.character(sort(unique(as.numeric(patch_ids[!is.na(patch_ids)]))))

    patchDE_object <- SingleCellExperiment(
        rowData = DataFrame(gene_id = rownames(spe)),
        colData = DataFrame(patch = unique(patch_ids))
    )

    rownames(patchDE_object) <- rownames(spe)
    colnames(patchDE_object) <- patch_ids

    for (predictor in predictor_cols) {
        de_res_predictor <- lapply(de_res[[predictor]], function(x) {
            x[, colnames(patchDE_object), drop = FALSE]
        })
        names(de_res_predictor) <- paste0(predictor, "_", names(de_res_predictor))
        assays(patchDE_object) <- de_res_predictor
    }

    if (metaanalysis) {
        W <- getPatchAttributes.spe(
            spe,
            dimred = embedding_name,
            patch_col = patch_column
        )

        W <- W[colnames(patchDE_object), , drop = FALSE]

        reducedDim(patchDE_object, "W") <- W

        de_res_meta <- patchMetaAnalysis(de_res, reducedDim(patchDE_object, "W"))

        for (predictor in predictor_cols) {
            de_res_meta_predictor <- lapply(de_res_meta[[predictor]], function(x) {
                x[, colnames(patchDE_object), drop = FALSE]
            })
            names(de_res_meta_predictor) <- paste0(predictor, "_meta_", names(de_res_meta_predictor))
            assays(patchDE_object) <- c(
                assays(patchDE_object),
                de_res_meta_predictor
            )
        }
    }

    return(patchDE_object)
}