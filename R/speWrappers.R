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
                            patch_diagnostics = TRUE,
                            k = 10L, strict_k = NULL,
                            patch_polys = TRUE,
                            patch_column = "patch",
                            verbose = TRUE) {

  xy <- SpatialExperiment::spatialCoords(spe)

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
    if (!Z %in% reducedDimNames(spe)) {
      stop("`Z` = '", Z, "' not found in reducedDimNames(spe). ",
           "Available: ", paste(reducedDimNames(spe), collapse = ", "))
    }
    Z_mat <- SingleCellExperiment::reducedDim(spe, Z)

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
       patch_metadata <- list(
           membership_log = membership_log,
           ss_log = ss_log
        )
        S4Vectors::metadata(spe)[["SpaceMosaic"]][["patch_iterations"]] <- patch_metadata
        patch_metadata$patch <- patch_vec


    } else {
      patch_metadata <- patch_vec
    }

    if(patch_diagnostics) {
      patch_data <- getPatchDiagnostics(
            xy = spatialCoords(spe),
            X = X_mat,
            patch = patch_metadata,
            k = k,
            strict_k = strict_k
        )
      metadata(spe)$SpaceMosaic$patch_diagnostics <- patch_data 
      patch_data <- patch_data$patch_diagnostics
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

    if (!embedding %in% reducedDimNames(spe)) {
      stop("`embedding` = '", embedding, "' not found in reducedDimNames(spe). ",
           "Available: ", paste(reducedDimNames(spe), collapse = ", "))
    }
    embedding_mat <- SingleCellExperiment::reducedDim(spe, embedding)
    reducedDim(spe, name) <- embedCellNeighborhoods(embedding_mat, spatialCoords(spe), ks, tissue)
    spe
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


patchDEWorkflow <- function(
    spe,
    predictor_cols,
    assay = "logcounts",
    patch_column = "patch",
    embedding_name = "Z",
    metaanalysis = TRUE,
    pearson = FALSE,
    tot = NULL,
    resid_mse = FALSE,
    verbose = TRUE
) {

    if (missing(predictor_cols) ||
        is.null(predictor_cols) ||
        length(predictor_cols) == 0) {
        stop("`predictor_cols` must be a non-empty character vector of colData names.")
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

    if (!patch_column %in% colnames(SummarizedExperiment::colData(spe))) {
        stop(sprintf(
            "'%s' not found in colData(spe).",
            patch_column
        ))
    }

    # Expression matrix: genes x cells -> cells x genes
    y <- t(SummarizedExperiment::assay(spe, assay))

    df <- SummarizedExperiment::colData(
        spe
    )[, predictor_cols, drop = FALSE]

    # Run patchDE
    de_res <- patchDE(
        y,
        df,
        SummarizedExperiment::colData(spe)[[patch_column]],
        pearson = pearson,
        tot = tot,
        resid_mse = resid_mse,
        verbose = verbose
    )

    # Patch IDs
    patch_ids <- SummarizedExperiment::colData(spe)[[patch_column]]
    patch_ids <- as.character(
        sort(unique(as.numeric(patch_ids[!is.na(patch_ids)])))
    )

    # Output SCE
    patchDE_object <- SingleCellExperiment::SingleCellExperiment(
        rowData = S4Vectors::DataFrame(
            gene_id = rownames(spe)
        ),
        colData = S4Vectors::DataFrame(
            patch = patch_ids
        )
    )

    rownames(patchDE_object) <- rownames(spe)
    colnames(patchDE_object) <- patch_ids

    # ------------------------------------------------------------------
    # Store p-values, estimates and SEs in metadata, and z-scores in
    # one assay per predictor.
    # ------------------------------------------------------------------

    predictor_metadata <- list()
    z_assays <- list()

    for (predictor in predictor_cols) {

        # Extract results for this predictor and retain the output
        # corresponding to the patches in patchDE_object.
        de_res_predictor <- lapply(
            de_res[[predictor]],
            function(x) {
                x[, colnames(patchDE_object), drop = FALSE]
            }
        )


        pvals <- de_res_predictor[["pvals"]]
        ests    <- de_res_predictor[["ests"]]
        ses     <- de_res_predictor[["ses"]]

        # Z score = estimate / standard error
        z <- ests / ses

        # Store z-score as the assay for this predictor
        z_assays[[predictor]] <- z

        # Store the three original quantities in metadata
        predictor_metadata[[predictor]] <- list(
            pvals = pvals,
            ests = ests,
            ses = ses
        )
    }

    # Add all predictor z-score assays at once
    SummarizedExperiment::assays(patchDE_object) <- z_assays

    # Store pvalues / estimates / SEs
    metadata(patchDE_object) <- predictor_metadata

    # ------------------------------------------------------------------
    # Meta-analysis
    # ------------------------------------------------------------------

    if (metaanalysis) {

        if (!embedding_name %in%
            SingleCellExperiment::reducedDimNames(spe)) {
            stop(sprintf(
                "Embedding '%s' not found in reducedDims(spe).",
                embedding_name
            ))
        }

        Z <- SingleCellExperiment::reducedDim(
            spe,
            embedding_name
        )

        patch <- SummarizedExperiment::colData(
            spe
        )[[patch_column]]

        W <- getPatchAttributes(Z, patch)

        W <- W[
            colnames(patchDE_object),
            ,
            drop = FALSE
        ]

        SingleCellExperiment::reducedDim(
            patchDE_object,
            "W"
        ) <- W

        de_res_meta <- patchMetaAnalysis(
            de_res,
            SingleCellExperiment::reducedDim(
                patchDE_object,
                "W"
            )
        )

        # Store meta-analysis z-scores as additional assays.
        meta_z_assays <- list()

        for (predictor in predictor_cols) {

    de_res_meta_predictor <- lapply(
        de_res_meta[[predictor]],
        function(x) {
            x[, colnames(patchDE_object), drop = FALSE]
        }
    )

    if (!all(c("pvals", "ests", "ses") %in%
             names(de_res_meta_predictor))) {
        stop(
            "Expected `de_res_meta[[predictor]]` to contain ",
            "`pvals`, `ests`, and `ses` for predictor '",
            predictor,
            "'."
        )
    }

    meta_ests <- de_res_meta_predictor[["ests"]]
    meta_ses  <- de_res_meta_predictor[["ses"]]

    meta_z_assays[[paste0(predictor, "_meta")]] <-
        meta_ests / meta_ses

    # Store meta-analysis results separately in metadata
    metadata(patchDE_object)[[paste0(predictor, "_meta")]] <- list(
        pvals = de_res_meta_predictor[["pvals"]],
        ests = meta_ests,
        ses = meta_ses
    )
}

        # Add meta-analysis z-score assays
        SummarizedExperiment::assays(patchDE_object) <-
            c(
                SummarizedExperiment::assays(patchDE_object),
                meta_z_assays
            )
    }

    return(patchDE_object)
}