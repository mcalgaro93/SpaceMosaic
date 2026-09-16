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
#' @export
embedCellNeighborhoods.spe <- function(spe, embedding, ks = c(5, 50), tissue = NULL,
                                       name = "Z") {

    if (!embedding %in% reducedDimNames(spe)) {
      stop("`embedding` = '", embedding, "' not found in reducedDimNames(spe). ",
           "Available: ", paste(reducedDimNames(spe), collapse = ", "))
    }
    if(!tissue %in% colnames(colData(spe)) && !is.null(tissue)){
      stop("`tissue` = '", tissue, "' not found in colData(spe). ",
           "Available: ", paste(colnames(colData(spe)), collapse = ", "))
    }
    embedding_mat <- SingleCellExperiment::reducedDim(spe, embedding)
    reducedDim(spe, name) <- embedCellNeighborhoods(embedding_mat, spatialCoords(spe), ks, tissue = colData(spe)[[tissue]])
    spe
}


#' @describeIn getPatches Method for `SpatialExperiment` objects.
#'
#' Convenience wrapper for the SpaceMosaic patching workflow on a
#' `SpatialExperiment`. The workflow consists of three steps:
#' `getPatches()` identifies spatial patches, `patchDiagnostics()` computes
#' patch-level information from the resulting patch assignments, and
#' `getPatchPolys()` constructs patch polygons using the patch assignments
#' and diagnostic information. Results are stored directly in the
#' `SpatialExperiment` object.
#' 
#' @param spe SpatialExperiment object. Must have rownames.
#' @param X Design variables for each spatial unit. A character vector of
#'   column names in `colData(spe)`. Each column is scaled to unit SD.
#' @param Z Optional per-cell context embeddings (cells x features). A single
#'   character string naming an entry in `reducedDims(spe)`.
#'   If supplied, patches will prefer Z-coherent regions. NULL disables.
#' @param patch_column Column name in `colData(spe)` where to store patch
#'   assignments. Default = 'patch'.
#' @param patch_diagnostics Logical; if TRUE, compute patch-level diagnostics
#'   using `patchDiagnostics()` and store the result in
#'   `metadata(spe)$SpaceMosaic$patch_diagnostics`.
#' @param patch_polys Logical; if TRUE, compute patch polygons using
#'   `patchPolys()` and store the result in
#'   `metadata(spe)$SpaceMosaic$patch_polys`.
#' @return The input `spe`, with results attached:
#'   \itemize{
#'     \item `colData(spe)[[patch_column]]`: named vector/factor of final patch
#'       assignments.
#'     \item `metadata(spe)$SpaceMosaic$patch_diagnostics`: patch-level
#'       diagnostics, if `patch_diagnostics = TRUE`.
#'     \item `metadata(spe)$SpaceMosaic$patch_polys`: patch polygons, if
#'       `patch_polys = TRUE`.
#'     \item `metadata(spe)$SpaceMosaic$patch_iterations`: if `log_iters = TRUE`,
#'       a list containing `membership_log` (n_cells x n_iters matrix of patch
#'       assignments per iteration) and `ss_log` (npatches x n_iters matrix of
#'       per-patch sum-of-squares of X across iterations).
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
            patch = patch_metadata,
            k = k,
            strict_k = strict_k
        )
       metadata(spe)$SpaceMosaic <- c(
        metadata(spe)$SpaceMosaic,
        patch_diagnostics_list
    )
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




#' Run patch-level differential expression
#'
#' Runs \code{patchDE()} on the expression data and patch assignments stored
#' in a \code{SpatialExperiment} object and packages the results into a
#' \code{SingleCellExperiment}.
#'
#' @param spe A \code{SpatialExperiment} object. Must contain the requested
#'   expression assay and patch assignments in \code{colData(spe)}.
#' @param predictor_cols Character vector of column names in
#'   \code{colData(spe)} to use as predictors in \code{patchDE()}.
#' @param assay Character string specifying the assay in \code{spe} containing
#'   the expression matrix. Default is \code{"logcounts"}.
#' @param patch_column Character string specifying the column in
#'   \code{colData(spe)} containing patch assignments. Default is
#'   \code{"patch"}.
#' @param method Differential-expression backend. `"hasty"` uses ordinary
#'   least squares and `"limma"` uses empirical-Bayes moderated inference.
#'   Default `"hasty"` for backward compatibility.
#' @param pearson Logical; passed to \code{patchDE()} to control whether
#'   the Pearson-based association method is used. Default is FALSE.
#' @param tot Optional character string specifying the column in
#'   \code{colData(spe)} containing per-cell total counts. If NULL, total
#'   counts are handled by \code{patchDE()}. Default is NULL.
#' @param resid_mse Logical; passed to \code{patchDE()} to control whether
#'   residual mean squared errors are returned. Default is FALSE.
#' @param return_residuals Logical; if TRUE, return an OLS residual matrix
#'   aligned with the rows and columns of `y`. Default FALSE.
#' @param verbose Show progress. Default TRUE.
#' @param verbose Logical; passed to \code{patchDE()} to control verbosity.
#'   Default is TRUE.
#'
#' @return A \code{SingleCellExperiment} object containing one column per
#'   patch and one row per gene. The object contains one assay per predictor,
#'   containing z-scores, and metadata entries containing the p-values,
#'   estimates, and standard errors returned by \code{patchDE()}.
#'
#' @details
#' The expression matrix is extracted from \code{assay(spe, assay)} and
#' transposed from genes-by-cells to cells-by-genes before being passed to
#' \code{patchDE()}. Predictor variables are taken from \code{colData(spe)},
#' and patch assignments are taken from \code{colData(spe)[[patch_column]]}.
#'
#' @seealso
#' \code{\link{patchDE}},
#' \code{\link{patchMetaAnalysis.spe}}
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
    verbose = TRUE
) {

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

    # Expression matrix: genes x cells -> cells x genes
    y <- t(SummarizedExperiment::assay(spe, assay))

    df <- SummarizedExperiment::colData(
        spe
    )[, predictor_cols, drop = FALSE]

    if (!is.null(tot)) {

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

    method <- match.arg(method)

    if(return_residuals){
        residuals <- de_res$residuals
        de_res <- de_res$de
    }



    # Patch IDs
    patch_ids <- SummarizedExperiment::colData(
        spe
    )[[patch_column]]

    
    # Output SCE
    patch_diagnostics <- metadata(spe)$SpaceMosaic$patch_diagnostics

    if (!is.null(patch_diagnostics)) {
        
        patch_ids <- intersect(patch_diagnostics$patch, as.character(patch_ids))
        colData <- patch_diagnostics[patch_diagnostics$patch %in% patch_ids,]

    } else {
        patch_ids <- as.character(
                sort(unique(as.numeric(
                    patch_ids[!is.na(patch_ids)]
                )))
            )
        colData <- S4Vectors::DataFrame(
            patch = patch_ids
        )
    }

    patchDE_object <- SingleCellExperiment::SingleCellExperiment(
        rowData = S4Vectors::DataFrame(
            gene_id = rownames(spe)
        ),
        colData = colData
    )
    
    rownames(patchDE_object) <- rownames(spe)
    colnames(patchDE_object) <- patch_ids

    # Store p-values, estimates, SEs, and z-scores
    predictor_metadata <- list()
    de_assays <- list()
    for (predictor in predictor_cols) {

        de_res_predictor <- lapply(
            de_res[[predictor]],
            function(x) {
                x[, colnames(patchDE_object), drop = FALSE]
            }
        )

        de_assays[[paste0(method,"_",predictor,"_pvals")]] <- de_res_predictor[["pvals"]] 
        de_assays[[paste0(method,"_",predictor,"_ests")]] <- de_res_predictor[["ests"]] 
        de_assays[[paste0(method,"_",predictor,"_ses")]] <- de_res_predictor[["ses"]]


        if(return_residuals){
            predictor_metadata[["residuals"]][[method]][[predictor]] <- residuals
        }
    }

  
    SummarizedExperiment::assays(patchDE_object) <- de_assays

    S4Vectors::metadata(patchDE_object) <- predictor_metadata


    patchDE_object
}

#' Run patch-level meta-analysis
#'
#' Performs patch-level meta-analysis on the results generated by
#' \code{patchDE.spe()}, using embedding-derived patch attributes.
#'
#' @param spe A \code{SpatialExperiment} object containing the embedding
#'   and patch assignments used to calculate patch-level attributes.
#' @param patchDE_object A \code{SingleCellExperiment} returned by
#'   \code{patchDE.spe()}.
#' @param embedding_name Character string specifying the entry in
#'   \code{reducedDims(spe)} to use for meta-analysis. Default is \code{"Z"}.
#' @param patch_column Character string specifying the column in
#'   \code{colData(spe)} containing patch assignments. Default is
#'   \code{"patch"}.
#'
#' @return A \code{SingleCellExperiment} containing the original
#'   patch-level differential expression results together with additional
#'   meta-analysis z-score assays, metadata entries, and the embedding-derived
#'   patch attributes stored as reduced dimension \code{"W"}.
#'
#' @details
#' The embedding specified by \code{embedding_name} is extracted from
#' \code{reducedDims(spe)}. Patch-level attributes are calculated using
#' \code{getPatchAttributes()} and supplied to \code{patchMetaAnalysis()}
#' together with the patch-level differential expression results.
#'
#' @seealso
#' \code{\link{patchDE.spe}},
#' \code{\link{patchMetaAnalysis.spe}},
#' \code{\link{getPatchAttributes}}
#'
#' @export
patchMetaAnalysis.spe <- function(
    spe,
    patchDE_object,
    embedding_name = "Z",
    patch_column = "patch",
    summarize_subgroups = FALSE,
    cellmeta_cols = NULL
) {

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

    # Check that patch IDs match
    spe_patches <- as.character(
        sort(unique(as.numeric(
            SummarizedExperiment::colData(spe)[[patch_column]]
            [!is.na(SummarizedExperiment::colData(spe)[[patch_column]])]
        )))
    )

    sce_patches <- colnames(patchDE_object)

    if (!identical(spe_patches, sce_patches)) {
        stop(
            "Patch IDs in `spe` do not match the columns of ",
            "`patchDE_object`."
        )
    }

    # Extract embedding
    Z <- SingleCellExperiment::reducedDim(
        spe,
        embedding_name
    )

    patch <- SummarizedExperiment::colData(
        spe
    )[[patch_column]]

    # Calculate patch-level embedding attributes
    W <- getPatchAttributes(
        Z,
        patch
    )

    W <- W[
        colnames(patchDE_object),
        ,
        drop = FALSE
    ]

    # Store patch attributes
    SingleCellExperiment::reducedDim(
        patchDE_object,
        "W"
    ) <- W

    # Run meta-analysis
    #
    # IMPORTANT:
    # patchMetaAnalysis() needs the original patchDE result, not the
    # SCE containing the formatted results, unless its API explicitly
    # accepts the latter.
    #
    # Therefore, see note below.

    DEObj <- list()
    assaysPatchDE <- assays(patchDE_object)
    for (nm in names(assaysPatchDE)) {
      parts <- strsplit(nm, "_")[[1]]
      
      var  <- parts[2]
      stat <- parts[3]
      
      if (is.null(DEObj[[var]])) {
        DEObj[[var]] <- list()
      }
      
      DEObj[[var]][[stat]] <- assaysPatchDE[[nm]]
    }
    de_res_meta <- patchMetaAnalysis(
        DEObj,
        W
    )

    # Store meta-analysis results
    meta_de_assays <- list()

    for (predictor in names(de_res_meta)) {

        de_res_meta_predictor <- lapply(
            de_res_meta[[predictor]],
            function(x) {
                x[, colnames(patchDE_object), drop = FALSE]
            }
        )

        if (!all(
            c("pvals", "ests", "ses") %in%
            names(de_res_meta_predictor)
        )) {
            stop(
                "Expected `de_res_meta[[predictor]]` to contain ",
                "`pvals`, `ests`, and `ses` for predictor '",
                predictor,
                "'."
            )
        }

        subgroups_mat <- de_res_meta_predictor[["subgroups"]]
        meta_de_assays[[paste0(parts[1],"_",predictor,"_meta_pvals")]] <- de_res_meta_predictor[["pvals"]] 
        meta_de_assays[[paste0(parts[1],"_",predictor,"_meta_ests")]] <- de_res_meta_predictor[["ests"]] 
        meta_de_assays[[paste0(parts[1],"_",predictor,"_meta_ses")]] <- de_res_meta_predictor[["ses"]]
        meta_de_assays[[paste0(parts[1],"_",predictor,"_meta_subgroups")]] <- subgroups_mat

    }



    assays <- SummarizedExperiment::assays(patchDE_object)
    assays[names(meta_de_assays)] <- meta_de_assays
    SummarizedExperiment::assays(patchDE_object) <- assays
    
    if(summarize_subgroups){
      metadata(patchDE_object)$subgroups_summary <- summarizeSubgroups(
        de_res_meta, patch = colData(spe)$patch, cellmeta = colData(spe)[, cellmeta_cols, drop = FALSE])
    }

    patchDE_object
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