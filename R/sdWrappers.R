# Internal helpers -----------------------------------------------------------

# spatialdataR is only suggested: require it when a function needs it.
.require_spatialdataR <- function() {
    if (!requireNamespace("spatialdataR", quietly = TRUE)) {
        msg <- paste("The 'spatialdataR' package is required to use this function.",
                     "Please install it with 'BiocManager::install(\"spatialdataR\")'.")
        cli::cli_abort(c("x" = msg))
    }
}

# Check that `sd` is a SpatialData object and that the other packages needed
# to work with it are available.
.checkSd <- function(sd) {
    for (pkg in c("sf", "duckspatial")) {
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

    # The shapes are a lazy DuckDB table whose row order is not stable across
    # queries, so instance IDs are read from the same collected sf object as
    # the geometries.
    geom <- spatialdataR::data(shape) |>
        sf::st_as_sf()
    ik <- tryCatch(spatialdataR::instance_key(sce), error = function(e) NULL)
    shape_ids <- if (!is.null(ik) && ik %in% names(geom)) geom[[ik]] else NULL
    if (!all(sf::st_geometry_type(geom) == "POINT")) {
        geom <- sf::st_centroid(sf::st_geometry(geom))
    }
    xy <- sf::st_coordinates(geom)[, c("X", "Y"), drop = FALSE]

    if (!is.null(ik)) {
        if (is.null(shape_ids)) {
            shape_ids <- seq_len(nrow(xy))
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
# spatialdataR subsets shapes by position (filter(row_number() %in% i)), which
# only selects the right rows if DuckDB returns them in a stable order; with
# several threads it does not, so the connection is made single-threaded.
# `instance_key` names the column of `x` holding instance IDs; without one,
# spatialdataR uses row positions, which tables cannot be matched against
# reliably.
.sdShape <- function(x, instance_key = NULL) {
    conn <- suppressMessages(duckspatial::ddbs_create_conn())
    DBI::dbExecute(conn, "SET threads = 1")
    table_id <- basename(tempfile("sdShape_"))
    suppressWarnings(suppressMessages(
        duckspatial::ddbs_write_table(
            conn = conn, data = x, name = table_id,
            overwrite = FALSE, temp_view = FALSE)
    ))
    x_duck <- duckspatial::as_duckspatial_df(
        x = table_id, conn = conn, crs = NA_character_,
        geom_col = attr(x, "sf_column"))
    meta <- spatialdataR::SpatialDataAttrs(type = "shape")
    if (!is.null(instance_key)) {
        spatialdataR::instance_key(meta) <- instance_key
    }
    spatialdataR::SpatialDataShape(x_duck, meta = meta)
}

# Make the table `y` annotate the shapes element `region`: column i of `y` is
# the geometry whose `instance_key` column equals `ids[i]`. spatialdataR does
# not export `region_key<-`, and `regions<-` fails on a table without a
# region key, so the key is written to the table's attributes directly.
.sdAnnotate <- function(y, region, instance_key, ids) {
    attrs <- SingleCellExperiment::int_metadata(y)$spatialdata_attrs
    attrs$region_key <- "region"
    SingleCellExperiment::int_metadata(y)$spatialdata_attrs <- attrs
    spatialdataR::instance_key(y) <- instance_key
    spatialdataR::regions(y) <- region
    spatialdataR::instances(y) <- ids
    y
}

# Unique non-missing patch IDs, in numeric order when all are numbers and in
# lexical order otherwise, as in getPatchDiagnostics().
.sortedPatchIds <- function(patch) {
    patch <- as.character(patch)
    patch_ids <- unique(patch[!is.na(patch)])
    numeric_patch_ids <- suppressWarnings(as.numeric(patch_ids))
    if (all(is.finite(numeric_patch_ids))) {
        patch_ids[order(numeric_patch_ids)]
    } else {
        sort(patch_ids)
    }
}

# Build the patch table: a SingleCellExperiment with one column per patch in
# `patch_ids` and no rows. Patch diagnostics go to its colData, the
# connectivity curve (patches x k) and the per-iteration sum of squares
# (patches x iterations) to its reducedDims, and the assignment summary to
# its metadata. Any of `diagnostics` and `ss_log` may be NULL.
.patchTable <- function(patch_ids, diagnostics = NULL, ss_log = NULL) {
    y <- SingleCellExperiment::SingleCellExperiment(
        colData = S4Vectors::DataFrame(row.names = patch_ids)
    )
    if (!is.null(diagnostics)) {
        pd <- diagnostics$patch_diagnostics
        SummarizedExperiment::colData(y) <- S4Vectors::DataFrame(
            pd[match(patch_ids, pd$patch), setdiff(names(pd), "patch"),
               drop = FALSE],
            row.names = patch_ids, check.names = FALSE
        )
        cc <- diagnostics$connectivity_curve
        ks <- sort(unique(cc$k))
        curve <- matrix(NA_real_, nrow = length(patch_ids), ncol = length(ks),
                        dimnames = list(patch_ids, paste0("k", ks)))
        curve[cbind(match(cc$patch, patch_ids), match(cc$k, ks))] <-
            cc$component_fraction
        SingleCellExperiment::reducedDim(y, "connectivity_curve") <- curve
        S4Vectors::metadata(y)$SpaceMosaic$assignment_summary <-
            diagnostics$assignment_summary
    }
    if (!is.null(ss_log)) {
        SingleCellExperiment::reducedDim(y, "ss_log") <-
            ss_log[patch_ids, , drop = FALSE]
    }
    y
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

    .require_spatialdataR()
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
#' `getPatchPolys()` constructs patch polygons. Patch assignments are stored in
#' `tables(sd)[[table]]`, patch polygons as a new shapes element of `sd`, and
#' patch-level results in a new table that annotates those polygons, with one
#' column per patch. Subsetting that table to some patches therefore also
#' drops the polygons of the other patches from `sd`.
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
#' @param log_iters Logical; if `TRUE` (default), per-iteration patch
#'   assignments and sums of squares are retained (see Value).
#' @param patch_diagnostics Logical; if `TRUE` (default), compute patch-level
#'   diagnostics with `getPatchDiagnostics()` and store them in
#'   `tables(sd)[[patch_table]]`.
#' @param patch_shape Name of the shapes element of `sd` where patch polygons
#'   are stored. Default `"patch_polys"`. An existing shapes element of the
#'   same name is overwritten.
#' @param patch_table Name of the table of `sd` where patch-level results are
#'   stored. Default `"patches"`. An existing table of the same name is
#'   overwritten.
#' @param patch_column Column name in `colData(tables(sd)[[table]])` where patch
#'   assignments are stored. Default `"patch"`. An existing column of the same
#'   name is overwritten.
#'
#' @return The input `sd`, with results attached:
#'   \itemize{
#'     \item `colData(tables(sd)[[table]])[[patch_column]]`: named vector/factor of final
#'       patch assignments, ordered to match `colnames(tables(sd)[[table]])`. `NA`
#'       denotes an unassigned cell.
#'     \item `reducedDim(tables(sd)[[table]], paste0(patch_column, "_iterations"))`:
#'       if `log_iters = TRUE`, a cells x iterations character matrix of patch
#'       assignments per iteration.
#'     \item `shapes(sd)[[patch_shape]]`: one polygon per patch (the convex hull of its cells), with a `patch` column
#'       holding the patch ID, which is also its instance ID. Patches with
#'       fewer than three cells have no polygon.
#'     \item `tables(sd)[[patch_table]]`: a `SingleCellExperiment` with no rows and one column per non-empty
#'       patch, named by patch ID, that annotates `shapes(sd)[[patch_shape]]`
#'       with `patch` as instance key. It holds:
#'       \itemize{
#'         \item in its `colData`, if `patch_diagnostics = TRUE`, the columns
#'           of the `patch_diagnostics` data frame returned by
#'           `getPatchDiagnostics()`;
#'         \item `reducedDim(, "connectivity_curve")`, if
#'           `patch_diagnostics = TRUE`, a patches x k matrix of the fraction
#'           of patch cells in the largest connected component at each k;
#'         \item `reducedDim(, "ss_log")`, if `log_iters = TRUE`, a patches x
#'           iterations matrix of per-patch sum-of-squares of X;
#'         \item `metadata()$SpaceMosaic$assignment_summary`, if
#'           `patch_diagnostics = TRUE`, the assignment summary returned by
#'           `getPatchDiagnostics()`.
#'       }
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
#' patches <- table(sd, "patches")
#' head(colData(patches))
#' shapes(sd)$patch_polys
#'
#' # Keeping only patches with at least 120 cells also drops the polygons of
#' # the other patches.
#' table(sd, "patches") <- patches[, patches$n_cells >= 120]
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
                           patch_shape = "patch_polys",
                           patch_table = "patches",
                           patch_column = "patch",
                           table = 1,
                           verbose = TRUE) {

  .require_spatialdataR()
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
        SingleCellExperiment::reducedDim(
            sce, paste0(patch_column, "_iterations")
        ) <- membership_log
        patch_info <- list(patch = patch_vec, membership_log = membership_log)
    } else {
        patch_vec  <- patchResult[colnames(sce)]
        patch_info <- patch_vec
    }

    SummarizedExperiment::colData(sce)[[patch_column]] <- patch_vec
    sd <- .sdSetTable(sd, sce, table)

    polys <- getPatchPolys(xy = xy, patch = patch_vec)
    if (nrow(polys) == 0) {
        warning("No patch has three or more cells, so no patch polygons ",
                "or patch table were stored.")
        return(sd)
    }
    sd <- spatialdataR::`shape<-`(
        sd, patch_shape,
        value = .sdShape(.patchPolysToSf(polys), instance_key = "patch")
    )

    diagnostics <- NULL
    if (patch_diagnostics) {
        diagnostics <- getPatchDiagnostics(
            xy = xy,
            X = X_mat,
            patches = patch_info,
            k = k,
            strict_k = strict_k
        )
    }
    patch_ids <- .sortedPatchIds(patch_vec)
    patch_sce <- .patchTable(
        patch_ids,
        diagnostics = diagnostics,
        ss_log = if (log_iters) patchResult$ss_log
    )
    patch_sce <- .sdAnnotate(patch_sce, patch_shape, "patch", patch_ids)
    sd <- .sdSetTable(sd, patch_sce, patch_table)

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

    .require_spatialdataR()
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

    .require_spatialdataR()
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

    .require_spatialdataR()
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

