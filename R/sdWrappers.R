# Internal helpers -----------------------------------------------------------

# Check that `sd` is a SpatialData object and that the packages needed to
# work with it are available.
.checkSd <- function(sd) {
    for (pkg in c("spatialdataR", "sf", "duckspatial")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop("Package '", pkg, "' is required for SpatialData objects. ",
                 "Please install it.")
        }
    }
    if (!methods::is(sd, "SpatialData")) {
        stop("`sd` must be a SpatialData object.")
    }
    invisible(TRUE)
}

# Resolve `table` (a name or index into `tables(sd)`) to a table name.
.sdTableName <- function(sd, table) {
    nms <- spatialdataR::tableNames(sd)
    if (length(nms) == 0) {
        stop("`sd` has no tables.")
    }
    if (!(is.character(table) || is.numeric(table)) || length(table) != 1L ||
        is.na(table)) {
        stop("`table` must be a single table name or index.")
    }
    if (is.numeric(table)) {
        if (table < 1 || table > length(nms)) {
            stop("`table` = ", table, " is out of range; `sd` has ",
                 length(nms), " table(s).")
        }
        table <- nms[table]
    }
    if (!table %in% nms) {
        stop("`table` = '", table, "' not found in tableNames(sd). ",
             "Available: ", paste(nms, collapse = ", "))
    }
    table
}

# Extract the SingleCellExperiment stored in `tables(sd)[[table]]`.
.sdTable <- function(sd, table) {
    sce <- spatialdataR::tables(sd)[[table]]
    if (!methods::is(sce, "SingleCellExperiment")) {
        stop("Table '", table, "' must be a SingleCellExperiment object.")
    }
    sce
}

# Extract a cells x 2 matrix of spatial coordinates, with rownames matching
# `colnames(sce)`. Coordinates come from the shapes element annotated by
# `sce` (its `region`), or from `shapes(sd)[[1]]` if `sce` has no region.
# Cells are matched to geometries by instance ID when `sce` has an instance
# key, and by position otherwise. Non-point geometries (e.g. cell polygons)
# are reduced to their centroids.
.sdCoords <- function(sd, sce) {
    if (length(spatialdataR::shapes(sd)) == 0) {
        stop("`shapes(sd)` is empty; cannot extract spatial coordinates.")
    }
    region <- unique(spatialdataR::region(sce))
    if (is.null(region)) {
        region <- spatialdataR::shapeNames(sd)[1]
    } else if (length(region) != 1L) {
        stop("The table annotates several elements (",
             paste(region, collapse = ", "), "); only tables annotating a ",
             "single shapes element are supported.")
    } else if (!region %in% spatialdataR::shapeNames(sd)) {
        stop("The table annotates '", region, "', which is not a shapes ",
             "element of `sd`.")
    }
    shape <- spatialdataR::shapes(sd)[[region]]

    geom <- spatialdataR::data(shape) |>
        sf::st_as_sf()
    if (!all(sf::st_geometry_type(geom) == "POINT")) {
        geom <- sf::st_centroid(sf::st_geometry(geom))
    }
    xy <- sf::st_coordinates(geom)[, c("X", "Y"), drop = FALSE]

    ik <- tryCatch(spatialdataR::instance_key(sce), error = function(e) NULL)
    if (!is.null(ik)) {
        shape_ids <- if (ik %in% names(shape)) {
            spatialdataR::pull(shape, !!ik)
        } else {
            seq_len(nrow(xy))
        }
        idx <- match(spatialdataR::instances(sce), shape_ids)
        if (anyNA(idx)) {
            stop(sum(is.na(idx)), " cell(s) in the table have no matching ",
                 "geometry in shapes element '", region, "'.")
        }
        xy <- xy[idx, , drop = FALSE]
    } else if (nrow(xy) != ncol(sce)) {
        stop("Shapes element '", region, "' has ", nrow(xy), " geometries ",
             "but the table has ", ncol(sce), " cells.")
    }
    rownames(xy) <- colnames(sce)
    xy
}

# Store `sce` back as `tables(sd)[[table]]`.
.sdSetTable <- function(sd, sce, table) {
    spatialdataR::`table<-`(sd, i = table, value = sce)
}

# Wrap an sf object as a SpatialDataShape. spatialdataR::SpatialDataShape()
# writes every data.frame/sf it receives to the same DuckDB table
# ("sdShape"), overwriting the geometries of all previously created shapes.
# Writing `x` to a uniquely named table ourselves and passing the resulting
# duckspatial_df avoids this.
.sdShape <- function(x) {
    conn <- suppressMessages(duckspatial::ddbs_create_conn())
    table_id <- basename(tempfile("sdShape_"))
    suppressWarnings(suppressMessages(
        duckspatial::ddbs_write_table(
            conn = conn, data = x, name = table_id,
            overwrite = FALSE, temp_view = FALSE)
    ))
    x_duck <- duckspatial::as_duckspatial_df(
        x = table_id, conn = conn, crs = NA_character_,
        geom_col = attr(x, "sf_column"))
    spatialdataR::SpatialDataShape(x_duck)
}

# Convert the vertex data frame returned by getPatchPolys() (columns `x`, `y`,
# `patch` and optional per-patch annotations) to an sf object with one
# polygon per patch and the annotations as attributes.
.patchPolysToSf <- function(polys) {
    ids <- unique(polys$patch)
    geoms <- lapply(ids, function(id) {
        v <- as.matrix(polys[polys$patch == id, c("x", "y")])
        sf::st_polygon(list(rbind(v, v[1, ])))
    })
    attrs <- polys[!duplicated(polys$patch),
                   setdiff(names(polys), c("x", "y")), drop = FALSE]
    rownames(attrs) <- NULL
    sf::st_sf(attrs, geometry = sf::st_sfc(geoms))
}


#' Embed cellular neighborhoods in a SpatialData object
#'
#' Creates a neighborhood embedding for a \code{SpatialData} object by
#' calling \code{\link{embedCellNeighborhoods}} on a reduced dimension stored
#' in \code{tables(sd)[[table]]} and the spatial coordinates of the cells, then storing
#' the result as a new reduced dimension of \code{tables(sd)[[table]]}.
#'
#' Spatial coordinates are taken from the shapes element annotated by the
#' table (its \code{region}), falling back to \code{shapes(sd)[[1]]} when
#' the table has no region. Cells are matched to geometries by instance ID
#' when the table has an instance key, and by position otherwise. Non-point
#' geometries are reduced to their centroids.
#'
#' @param sd A SpatialData object whose \code{tables(sd)[[table]]} is a
#'   SingleCellExperiment.
#' @param table Name or index of the table in \code{tables(sd)} to use.
#'   Default \code{1} (the first table).
#' @param embedding A single character string giving the name of a reduced
#'   dimension stored in \code{reducedDim(tables(sd)[[table]])} with a
#'   single cell embeddings matrix (cells x features).
#' @param ks Vector giving the number of nearest neighbors for each scale.
#'   Default \code{c(5, 50)}.
#' @param tissue Optional character string giving the name of a column in
#'   \code{colData(tables(sd)[[table]])} containing tissue IDs to prevent cross-tissue
#'   neighbor edges. Default NULL.
#' @param name Character string giving the name under which the result is
#'   stored via \code{reducedDim(tables(sd)[[table]], name)}. Default \code{"Z"}.
#' @return \code{sd} with a new reduced dimension, named according to
#'   \code{name}, added to \code{tables(sd)[[table]]}. This matrix has
#'   dimensions n cells x (ncol(embedding_mat) * length(ks)).
#' @seealso \code{\link{embedCellNeighborhoods}}, which does the multi-scale
#'   neighborhood averaging on a plain matrix, and
#'   \code{\link{embedCellNeighborhoods.spe}}.
#'
#' @importFrom SingleCellExperiment reducedDim reducedDim<- reducedDimNames
#' @importFrom SummarizedExperiment colData
#' @examplesIf requireNamespace("spatialdataR", quietly = TRUE)
#' library(SpatialExperiment)
#' library(spatialdataR)
#' spe <- readRDS(system.file("extdata", "cosmx_carcinoma.rds", package = "SpaceMosaic"))
#' sd <- spe_to_spatialdata(spe)
#' sd <- embedCellNeighborhoods.sd(sd, embedding = "PCA", ks = c(5, 50), tissue = 'sample_id')
#' reducedDim(table(sd), "Z")[1:3, 1:4]
#' @export
embedCellNeighborhoods.sd <- function(sd, embedding, ks = c(5, 50), tissue = NULL,
                                      name = "Z", table = 1) {

    .checkSd(sd)
    table <- .sdTableName(sd, table)
    sce <- .sdTable(sd, table)

    if (!embedding %in% reducedDimNames(sce)) {
      stop("`embedding` = '", embedding, "' not found in reducedDimNames(tables(sd)[[table]]). ",
           "Available: ", paste(SingleCellExperiment::reducedDimNames(sce), collapse = ", "))
    }
    if (!is.null(tissue) && !tissue %in% colnames(SummarizedExperiment::colData(sce))) {
      stop("`tissue` = '", tissue, "' not found in colData(tables(sd)[[table]]). ",
           "Available: ", paste(colnames(SummarizedExperiment::colData(sce)), collapse = ", "))
    }

    tissue_ids <- if (is.null(tissue)) {
      NULL
    } else {
      SummarizedExperiment::colData(sce)[[tissue]]
    }

    embedding_mat <- SingleCellExperiment::reducedDim(sce, embedding)
    SingleCellExperiment::reducedDim(sce, name) <- embedCellNeighborhoods(embedding_mat, .sdCoords(sd, sce), ks, tissue = tissue_ids)
    return(.sdSetTable(sd, sce, table))
}


#' Run the SpaceMosaic patching workflow on a `SpatialData` object
#'
#' Convenience wrapper for the SpaceMosaic patching workflow on a
#' `SpatialData` object. This wrapper contains three steps:
#' `getPatches()` identifies spatial patches, `getPatchDiagnostics()` computes
#' patch-level information from the resulting patch assignments, and
#' `getPatchPolys()` constructs patch polygons using the patch assignments and
#' diagnostic information. Patch assignments and diagnostics are stored in
#' `tables(sd)[[table]]`, and patch polygons as a new shapes element of `sd`.
#'
#' Spatial coordinates are taken from the shapes element annotated by the
#' table (see [embedCellNeighborhoods.sd()] for details), design variables from `colData(tables(sd)[[table]])`, and (optionally) context
#' embeddings from `reducedDims(tables(sd)[[table]])`.
#'
#' @param sd A `SpatialData` object. `tables(sd)[[table]]` must be a
#'   `SingleCellExperiment` with `colnames` (cell identifiers) and the columns
#'   named in `X` present in its `colData`, and annotate a shapes element
#'   holding the cell geometries.
#' @param X Design variables for each spatial unit. A character vector of
#'   column names in `colData(tables(sd)[[table]])`. Each column is scaled to unit SD by
#'   `getPatches()`.
#' @param npatches Number of patches to create.
#' @param Z Character string or `NULL`. Name of an entry in
#'   `reducedDimNames(tables(sd)[[table]])` holding per-cell context embeddings
#'   (cells x features). If supplied, patches will prefer Z-coherent regions.
#'   Default `NULL` disables this behaviour.
#' @inheritParams embedCellNeighborhoods.sd
#' @inheritParams getPatches.spe
#' @param log_iters Logical; if `TRUE` (default), per-iteration diagnostics are
#'   retained and stored in
#'   `metadata(tables(sd)[[table]])$SpaceMosaic$patch_iterations`.
#' @param patch_diagnostics Logical; if `TRUE` (default), compute patch-level
#'   diagnostics with `getPatchDiagnostics()` and merge its elements
#'   (`patch_diagnostics`, `assignment_summary`, `connectivity_curve`) into
#'   `metadata(tables(sd)[[table]])$SpaceMosaic`.
#' @param patch_polys Logical; if `TRUE` (default), compute patch polygons with
#'   `getPatchPolys()` and store them as the shapes element `patch_shape` of
#'   `sd`. When `patch_diagnostics = FALSE`, polygons are built without
#'   diagnostic information (`patch_data = NULL`).
#' @param patch_shape Name of the shapes element of `sd` where patch polygons
#'   are stored when `patch_polys = TRUE`. Default `"patch_polys"`. An existing
#'   shapes element of the same name is overwritten.
#' @param patch_column Column name in `colData(tables(sd)[[table]])` where patch
#'   assignments are stored. Default `"patch"`. An existing column of the same
#'   name is overwritten.
#'
#' @return The input `sd`, with results attached:
#'   \itemize{
#'     \item `colData(tables(sd)[[table]])[[patch_column]]`: named vector/factor of final
#'       patch assignments, ordered to match `colnames(tables(sd)[[table]])`. `NA`
#'       denotes an unassigned cell.
#'     \item `metadata(tables(sd)[[table]])$SpaceMosaic$patch_iterations`: if
#'       `log_iters = TRUE`, a list with `patch` (final assignments),
#'       `membership_log` (n_cells x n_iters matrix of patch assignments per
#'       iteration) and `ss_log` (npatches x n_iters matrix of per-patch
#'       sum-of-squares of X across iterations).
#'     \item `metadata(tables(sd)[[table]])$SpaceMosaic$patch_diagnostics`,
#'       `$assignment_summary` and `$connectivity_curve`: the three data frames
#'       returned by `getPatchDiagnostics()`, if `patch_diagnostics = TRUE`.
#'     \item `shapes(sd)[[patch_shape]]`: if `patch_polys = TRUE`, one polygon
#'       per patch (the convex hull of its cells), with a `patch` column and,
#'       if `patch_diagnostics = TRUE`, the columns of `patch_diagnostics`.
#'       Patches with fewer than three cells have no polygon.
#'   }
#'
#' @seealso [getPatches()], [getPatchDiagnostics()], [getPatchPolys()],
#'   [getPatches.spe()]
#'
#' @importFrom SingleCellExperiment reducedDimNames
#' @importFrom S4Vectors metadata metadata<-
#'
#' @examplesIf requireNamespace("spatialdataR", quietly = TRUE)
#' library(SpatialExperiment)
#' library(spatialdataR)
#' spe <- readRDS(system.file("extdata", "cosmx_carcinoma.rds", package = "SpaceMosaic"))
#' sd <- spe_to_spatialdata(spe)
#' sd <- getPatches.sd(
#'   sd = sd,
#'   X = "distance_to_fibroblast",
#'   npatches = 50
#' )
#' head(table(sd)$patch)
#' head(metadata(table(sd))$SpaceMosaic$patch_diagnostics)
#' shapes(sd)$patch_polys
#'
#' @export getPatches.sd
getPatches.sd <- function(sd, X, npatches,
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
                           patch_shape = "patch_polys",
                           patch_column = "patch",
                           table = 1,
                           verbose = TRUE) {

  .checkSd(sd)
  table <- .sdTableName(sd, table)
  sce <- .sdTable(sd, table)
  xy <- .sdCoords(sd, sce)

  if (missing(X) || is.null(X) || length(X) == 0) {
    stop("`X` must be a non-empty character vector of colData names.")
  }

  missing_cols <- setdiff(X, colnames(SummarizedExperiment::colData(sce)))
    if (length(missing_cols) > 0) {
      stop("The following `X` are not columns of colData(tables(sd)[[table]]): ",
           paste(missing_cols, collapse = ", "))
    }

  X_mat <- as.matrix(SummarizedExperiment::colData(sce)[, X, drop = FALSE])

  Z_mat <- NULL
  if (!is.null(Z)) {
    if (!is.character(Z) || length(Z) != 1L || is.na(Z)) {
        stop("`Z` must be NULL or a single character string naming a reduced dimension.")
    }

    if (!Z %in% SingleCellExperiment::reducedDimNames(sce)) {
        stop(
        "`Z` = '", Z, "' not found in reducedDimNames(tables(sd)[[table]]). ",
        "Available: ", paste(SingleCellExperiment::reducedDimNames(sce), collapse = ", ")
        )
    }

  Z_mat <- SingleCellExperiment::reducedDim(sce, Z)
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
          patch_vec      <- patchResult$patch[colnames(sce)]
          membership_log <- patchResult$membership_log[colnames(sce), , drop = FALSE]
          ss_log         <- patchResult$ss_log
      } else {
          patch_vec <- patchResult[colnames(sce)]
      }

      SummarizedExperiment::colData(sce)[[patch_column]] <- patch_vec

      if (log_iters) {
          S4Vectors::metadata(sce)[["SpaceMosaic"]][["patch_iterations"]] <- list(
              patch = patch_vec,
              membership_log = membership_log,
              ss_log = ss_log
          )
          patch_metadata <- S4Vectors::metadata(sce)[["SpaceMosaic"]][["patch_iterations"]]
      } else {
          patch_metadata <- patch_vec
      }

    if(patch_diagnostics) {
      patch_diagnostics_list <- getPatchDiagnostics(
            xy = xy,
            X = X_mat,
            patches = patch_metadata,
            k = k,
            strict_k = strict_k
        )
      metadata(sce)$SpaceMosaic[names(patch_diagnostics_list)] <- patch_diagnostics_list
      patch_data <- patch_diagnostics_list$patch_diagnostics
    } else{
      patch_data <- NULL
    }

    sd <- .sdSetTable(sd, sce, table)

    if (patch_polys) {
      polys <- getPatchPolys(
        xy = xy,
        patch = patch_vec,
        patch_data = patch_data
      )
      if (nrow(polys) > 0) {
        sd <- spatialdataR::`shape<-`(
          sd, patch_shape, value = .sdShape(.patchPolysToSf(polys))
        )
      }
    }
    return(sd)
}



#' Run patch-level differential expression on a SpatialData object
#'
#' Runs \code{patchDE()} on the expression data and patch assignments stored
#' in \code{tables(sd)[[table]]} and returns the raw results alongside the object they
#' were computed from.
#'
#' @param sd A \code{SpatialData} object. \code{tables(sd)[[table]]} must contain the
#'   requested expression assay and patch assignments in its \code{colData}.
#' @param predictor_cols Character vector of column names in
#'   \code{colData(tables(sd)[[table]])} to use as predictors in \code{patchDE()}. Must
#'   be non-empty.
#' @param assay Character string naming the assay in \code{tables(sd)[[table]]} that
#'   contains the expression matrix. Default \code{"logcounts"}.
#' @param patch_column Character string naming the column in
#'   \code{colData(tables(sd)[[table]])} that contains patch assignments. Cells with
#'   \code{NA} here are dropped by \code{patchDE()}. Default \code{"patch"}.
#' @param tot Optional character string naming a column in
#'   \code{colData(tables(sd)[[table]])} that holds per-cell total counts. Required when
#'   \code{pearson = TRUE}; ignored otherwise. Default NULL.
#' @param return_residuals Logical; if TRUE, the cells-by-genes residual matrix
#'   produced by \code{patchDE()} is transposed to genes-by-cells and stored
#'   as assay \code{residual_assay} in \code{tables(sd)[[table]]}. Default FALSE.
#' @param residual_assay Character string naming the assay under which
#'   residuals are stored in \code{tables(sd)[[table]]} when
#'   \code{return_residuals = TRUE}. Matches the default \code{assay_name}
#'   expected by \code{\link{moranTest.sd}}. Default \code{"residuals"}.
#' @inheritParams embedCellNeighborhoods.sd
#' @inheritParams patchDE.spe
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{de}}{The differential-expression results exactly as returned
#'       by \code{patchDE()}: one element per entry of \code{predictor_cols},
#'       each holding the \code{pvals}, \code{ests} and \code{ses} matrices
#'       (genes by patches), plus anything else the backend produced, such as
#'       residual mean squared errors when \code{resid_mse = TRUE}.}
#'     \item{\code{sd}}{The input \code{SpatialData} object. When
#'       \code{return_residuals = TRUE}, \code{tables(sd)[[table]]} includes a new
#'       \code{residual_assay} assay (genes by cells) holding the residual
#'       matrix from \code{patchDE()}; otherwise it is returned unchanged.}
#'     \item{\code{table}}{The name of the table that was used, so that
#'       \code{patchMetaAnalysis.sd()} picks the same one by default.}
#'   }
#'
#' @details
#' The expression matrix is extracted from \code{assay(tables(sd)[[table]], assay)} and
#' transposed from genes-by-cells to cells-by-genes before being passed to
#' \code{patchDE()}. Predictor variables are taken from
#' \code{colData(tables(sd)[[table]])}, and patch assignments from
#' \code{colData(tables(sd)[[table]])[[patch_column]]}.
#'
#' Differential expression is run on every non-missing patch, and the patch
#' columns of the returned matrices follow the order \code{patchDE()}
#' produced.
#'
#' \code{patchDE()} returns one residual matrix per fit, not one per
#' predictor, so cells with a missing patch assignment (and therefore no fit)
#' get \code{NA} residuals in \code{residual_assay}.
#'
#' @seealso
#' \code{\link{patchDE}},
#' \code{\link{patchDE.spe}},
#' \code{\link{patchMetaAnalysis.sd}},
#' \code{\link{moranTest.sd}}
#'
#' @examplesIf requireNamespace("spatialdataR", quietly = TRUE)
#' library(SpatialExperiment)
#' library(spatialdataR)
#' spe <- readRDS(system.file("extdata", "cosmx_carcinoma.rds", package = "SpaceMosaic"))
#' spe <- spe[, spe$celltype == "Cancer.cells"]
#' sd <- spe_to_spatialdata(spe)
#' sd <- getPatches.sd(
#'   sd = sd,
#'   X = "distance_to_fibroblast",
#'   npatches = 50
#' )
#' de_result <- patchDE.sd(sd,
#'   predictor_cols = "distance_to_fibroblast"
#' )
#' names(de_result)
#' head(de_result$de$distance_to_fibroblast$pvals)
#'
#' @export
patchDE.sd <- function(
    sd,
    predictor_cols,
    assay = "logcounts",
    patch_column = "patch",
    method = c("hasty", "limma"),
    pearson = FALSE,
    tot = NULL,
    resid_mse = FALSE,
    return_residuals = FALSE,
    residual_assay = "residuals",
    table = 1,
    verbose = TRUE
) {

    method <- match.arg(method)

    .checkSd(sd)
    table <- .sdTableName(sd, table)
    sce <- .sdTable(sd, table)

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
        colnames(SummarizedExperiment::colData(sce))
    )

    if (length(missing_cols) > 0) {
        stop(
            "The following `predictor_cols` are not columns of colData(tables(sd)[[table]]): ",
            paste(missing_cols, collapse = ", ")
        )
    }

    if (!patch_column %in%
        colnames(SummarizedExperiment::colData(sce))) {
        stop(sprintf(
            "'%s' not found in colData(tables(sd)[[table]]).",
            patch_column
        ))
    }

    if (!assay %in% SummarizedExperiment::assayNames(sce)) {
        stop(sprintf(
            "'%s' not found in assay names of `tables(sd)[[table]]`. Available assays: %s",
            assay, paste(SummarizedExperiment::assayNames(sce), collapse = ", ")
        ))
    }

    if (pearson && is.null(tot)) {
        stop("`tot` must be supplied when `pearson = TRUE`.")
    }

    # Expression matrix: genes x cells -> cells x genes
    y <- t(SummarizedExperiment::assay(sce, assay))

    df <- SummarizedExperiment::colData(
        sce
    )[, predictor_cols, drop = FALSE]

    if (pearson && !is.null(tot)) {

        if (!tot %in%
            colnames(SummarizedExperiment::colData(sce))) {
            stop(sprintf(
                "'%s' not found in colData(tables(sd)[[table]]).",
                tot
            ))
        }

        tot <- SummarizedExperiment::colData(sce)[[tot]]
    }

    # Run patchDE
    de_res <- patchDE(
        y,
        df,
        SummarizedExperiment::colData(sce)[[patch_column]],
        method = method,
        pearson = pearson,
        tot = tot,
        resid_mse = resid_mse,
        return_residuals = return_residuals,
        verbose = verbose
    )

    # patchDE() wraps its output when residuals are requested
    if (return_residuals) {
        SummarizedExperiment::assay(sce, residual_assay) <- t(de_res$residuals)
        sd <- .sdSetTable(sd, sce, table)
        out <- list(de = de_res$de)
    } else {
        out <- list(de = de_res)
    }

    out$sd <- sd
    out$table <- table
    return(out)
}

#' Meta-analyse patch-level differential expression across an embedding
#'
#' Takes the output of \code{patchDE.sd()}, derives patch-level attributes
#' from a reduced-dimension embedding stored in \code{tables(sd)[[table]]}, and runs
#' \code{patchMetaAnalysis()} to combine the per-patch fits.
#'
#' @param patchDE_result The list returned by \code{patchDE.sd()}. Must
#'   contain a \code{de} element holding the per-predictor differential
#'   expression results and an \code{sd} element holding the
#'   \code{SpatialData} object they were computed from.
#' @param embedding_name Character string naming the entry of
#'   \code{reducedDims(tables(sd)[[table]])} used to derive patch attributes. Default
#'   \code{"Z"}.
#' @param patch_column Character string naming the column in
#'   \code{colData(tables(sd)[[table]])} that contains patch assignments. Default
#'   \code{"patch"}.
#' @param cellmeta_cols Optional character vector of
#'   \code{colData(tables(sd)[[table]])} columns passed to \code{summarizeSubgroups()}
#'   as cell metadata. Used only when \code{summarize_subgroups = TRUE}.
#'   Default NULL.
#' @param table Name or index of the table in \code{tables(sd)} to use.
#'   Default \code{NULL} uses the table \code{patchDE.sd()} was run on.
#' @inheritParams patchMetaAnalysis.spe
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
#' in \code{colData(tables(sd)[[table]])[[patch_column]]} and those in the differential
#' expression results must be the same set; their order need not agree.
#'
#' @seealso
#' \code{\link{patchDE.sd}},
#' \code{\link{patchMetaAnalysis}},
#' \code{\link{patchMetaAnalysis.spe}},
#' \code{\link{getPatchAttributes}},
#' \code{\link{summarizeSubgroups}}
#'
#' @examplesIf requireNamespace("spatialdataR", quietly = TRUE)
#' library(SpatialExperiment)
#' library(spatialdataR)
#' spe <- readRDS(system.file("extdata", "cosmx_carcinoma.rds", package = "SpaceMosaic"))
#' # Embed neighborhoods on all cells, then keep cancer cells only
#' spe <- embedCellNeighborhoods.spe(spe, embedding = "PCA", ks = c(5, 50), tissue = 'sample_id')
#' spe <- spe[, spe$celltype == "Cancer.cells"]
#' sd <- spe_to_spatialdata(spe)
#' sd <- getPatches.sd(
#'   sd = sd,
#'   X = "distance_to_fibroblast",
#'   npatches = 50
#' )
#' de_result <- patchDE.sd(sd,
#'   predictor_cols = "distance_to_fibroblast"
#' )
#' meta_result <- patchMetaAnalysis.sd(de_result, embedding_name = "Z")
#' head(meta_result$meta$distance_to_fibroblast$pvals)
#'
#' @export
patchMetaAnalysis.sd <- function(
    patchDE_result,
    embedding_name = "Z",
    patch_column = "patch",
    summarize_subgroups = FALSE,
    cellmeta_cols = NULL,
    table = NULL
) {

    if (!is.list(patchDE_result) ||
        !all(c("de", "sd") %in% names(patchDE_result))) {
        stop(
            "`patchDE_result` must be a list with `de` and `sd` elements, ",
            "as returned by `patchDE.sd()`."
        )
    }

    .checkSd(patchDE_result$sd)
    if (is.null(table)) {
        table <- if (is.null(patchDE_result$table)) 1 else patchDE_result$table
    }
    table <- .sdTableName(patchDE_result$sd, table)
    sce <- .sdTable(patchDE_result$sd, table)
    de_res <- patchDE_result$de

    if (length(de_res) == 0 || is.null(names(de_res))) {
        stop(
            "`patchDE_result$de` must be a named list with one element ",
            "per predictor."
        )
    }

    if (!patch_column %in%
        colnames(SummarizedExperiment::colData(sce))) {
        stop(sprintf(
            "'%s' not found in colData(tables(sd)[[table]]).",
            patch_column
        ))
    }

    if (!embedding_name %in%
        SingleCellExperiment::reducedDimNames(sce)) {
        stop(sprintf(
            "Embedding '%s' not found in reducedDims(tables(sd)[[table]]).",
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

    patch <- SummarizedExperiment::colData(sce)[[patch_column]]

    # Extract embedding
    Z <- SingleCellExperiment::reducedDim(
        sce,
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
                colnames(SummarizedExperiment::colData(sce))
            )

            if (length(missing_cols) > 0) {
                stop(
                    "The following `cellmeta_cols` are not columns of ",
                    "colData(tables(sd)[[table]]): ",
                    paste(missing_cols, collapse = ", ")
                )
            }
        }

        out$subgroups_summary <- summarizeSubgroups(
            de_res_meta,
            patch = patch,
            cellmeta = SummarizedExperiment::colData(
                sce
            )[, cellmeta_cols, drop = FALSE]
        )
    }
    return(out)
}

#' Moran's I test on a SpatialData object
#'
#' moranTest() method for \code{SpatialData} objects. Extracts residuals from
#' the specified assay of \code{tables(sd)[[table]]} and spatial coordinates from
#' the shapes element it annotates, then dispatches to \code{moranTest}.
#'
#' @param sd A \code{SpatialData} object.
#' @param assay_name Character; name of the assay in \code{tables(sd)[[table]]}
#'   containing the values to test (e.g. Pearson residuals). Default
#'   \code{"residuals"}.
#' @param patch_column Name of the column in `colData(tables(sd)[[table]])` where patch
#'   assignments are stored. Default `"patch"`.
#' @inheritParams embedCellNeighborhoods.sd
#' @inheritParams moranTest.spe
#' @return A data frame with one row per gene and patch. It contains patch and
#'   gene identifiers, number of cells, observed and expected Moran's I, raw and
#'   adjusted empirical p-values, graph and permutation settings, and `status`.
#'   A status of `"ok"` indicates a completed test; non-testable inputs such as
#'   constant or non-finite residuals return missing statistics with an
#'   explanatory status rather than stopping the remaining tests.
#'
#' @seealso \code{\link{moranTest}}, \code{\link{moranTest.spe}}
#'
#' @importFrom SummarizedExperiment assayNames assay colData
#'
#' @examplesIf requireNamespace("spatialdataR", quietly = TRUE)
#' library(SpatialExperiment)
#' library(spatialdataR)
#' spe <- readRDS(system.file("extdata", "cosmx_carcinoma.rds", package = "SpaceMosaic"))
#' spe <- spe[1:10, spe$celltype == "Cancer.cells"]
#' sd <- spe_to_spatialdata(spe)
#' sd <- getPatches.sd(
#'   sd = sd,
#'   X = "distance_to_fibroblast",
#'   npatches = 50
#' )
#' de_result <- patchDE.sd(sd,
#'   predictor_cols = "distance_to_fibroblast",
#'   return_residuals = TRUE
#' )
#' moran_result <- moranTest.sd(de_result$sd)
#' head(moran_result)
#' @export
moranTest.sd <- function(sd, assay_name = 'residuals', patch_column = "patch", k = 10L,
                         n_permutations = 999L,
                         alternative = c("greater", "less", "two.sided"),
                         p_adjust_method = "BH",
                         adjustment_scope = c("global", "patch", "gene"),
                         table = 1) {

    .checkSd(sd)
    table <- .sdTableName(sd, table)
    sce <- .sdTable(sd, table)

    if (!is.character(assay_name) || length(assay_name) != 1) {
      stop("'assay_name' must be a single character string.")
    }

    if (!assay_name %in% assayNames(sce)) {
      stop(sprintf("'%s' not found in assay names of 'tables(sd)[[table]]'. Available assays: %s",
                   assay_name, paste(assayNames(sce), collapse = ", ")))
    }

    if (!patch_column %in% colnames(colData(sce))) {
      stop(sprintf("'%s' not found in colData(tables(sd)[[table]]). Available: %s",
                   patch_column, paste(colnames(colData(sce)), collapse = ", ")))
    }

    res <- moranTest(residuals = t(assay(sce, assay_name)),
                     xy = .sdCoords(sd, sce),
                     patch = colData(sce)[[patch_column]],
                     k = k,
                     n_permutations = n_permutations,
                     alternative = alternative,
                     p_adjust_method = p_adjust_method,
                     adjustment_scope = adjustment_scope)
    return(res)
}

