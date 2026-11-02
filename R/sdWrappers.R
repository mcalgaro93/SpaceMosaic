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

# Patch table `patch_table` of `sd`, as created by getPatches.sd().
.sdPatchTable <- function(sd, patch_table) {
    if (!patch_table %in% spatialdataR::tableNames(sd)) {
        stop("`sd` has no table '", patch_table, "'; run getPatches.sd() ",
             "first or set `patch_table`.")
    }
    spatialdataR::tables(sd)[[patch_table]]
}

# Rebuild the patch table `y` with one row per gene in `genes`, dropping the
# results of patchDE.sd(), patchMetaAnalysis.sd() and moranTest.sd() (assays,
# the "W" reducedDim, `moran_` colData columns and their metadata) while
# keeping those of getPatches.sd().
.patchTableReset <- function(y, genes) {
    cd <- SummarizedExperiment::colData(y)
    out <- SingleCellExperiment::SingleCellExperiment(
        rowData = S4Vectors::DataFrame(row.names = genes),
        colData = cd[, !startsWith(names(cd), "moran_"), drop = FALSE]
    )
    rd <- SingleCellExperiment::reducedDims(y)
    SingleCellExperiment::reducedDims(out) <- rd[names(rd) != "W"]
    md <- S4Vectors::metadata(y)
    md$SpaceMosaic[c("patchDE", "subgroups_summary", "moranTest")] <- NULL
    S4Vectors::metadata(out) <- md
    .sdAnnotate(out, spatialdataR::region(y), spatialdataR::instance_key(y),
                spatialdataR::instances(y))
}

# Remove the assays of `y` whose names start with `prefix`.
.dropAssays <- function(y, prefix) {
    keep <- !startsWith(as.character(SummarizedExperiment::assayNames(y)),
                        prefix)
    SummarizedExperiment::assays(y) <- SummarizedExperiment::assays(y)[keep]
    y
}

# Store the genes x patches matrices in `mats` as assays of the patch table
# `y`, matched to its rows and columns by name. Entries of `y` without a
# match are NA; genes and patches of `mats` missing from `y` are dropped
# with a warning.
.patchTableAssays <- function(y, mats) {
    dropped_genes <- unique(unlist(lapply(mats, function(m)
        setdiff(rownames(m), rownames(y)))))
    dropped_patches <- unique(unlist(lapply(mats, function(m)
        setdiff(colnames(m), colnames(y)))))
    if (length(dropped_genes) > 0) {
        warning(length(dropped_genes), " gene(s) are not rows of the patch ",
                "table and were dropped: ",
                paste(utils::head(dropped_genes, 5), collapse = ", "),
                if (length(dropped_genes) > 5) ", ...")
    }
    if (length(dropped_patches) > 0) {
        warning(length(dropped_patches), " patch(es) are not columns of the ",
                "patch table and were dropped: ",
                paste(dropped_patches, collapse = ", "))
    }
    for (nm in names(mats)) {
        m <- mats[[nm]]
        full <- matrix(NA, nrow(y), ncol(y), dimnames = dimnames(y))
        storage.mode(full) <- storage.mode(m)
        i <- match(rownames(m), rownames(y))
        j <- match(colnames(m), colnames(y))
        full[i[!is.na(i)], j[!is.na(j)]] <-
            m[!is.na(i), !is.na(j), drop = FALSE]
        SummarizedExperiment::assay(y, nm) <- full
    }
    y
}

# Store the convex hulls of the patches in `patch` (cells with coordinates
# `xy`) as the shapes element `patch_shape` of `sd`, and a patch table
# annotating them (see .patchTable()) as the table `patch_table`.
.sdAddPatchTable <- function(sd, xy, patch, patch_shape, patch_table,
                             diagnostics = NULL, ss_log = NULL) {
    polys <- getPatchPolys(xy = xy, patch = patch)
    if (nrow(polys) == 0) {
        warning("No patch has three or more cells, so no patch polygons ",
                "or patch table were stored.")
        return(sd)
    }
    sd <- spatialdataR::`shape<-`(
        sd, patch_shape,
        value = .sdShape(.patchPolysToSf(polys), instance_key = "patch")
    )
    patch_ids <- .sortedPatchIds(patch)
    patch_sce <- .patchTable(patch_ids, diagnostics = diagnostics,
                             ss_log = ss_log)
    patch_sce <- .sdAnnotate(patch_sce, patch_shape, "patch", patch_ids)
    .sdSetTable(sd, patch_sce, patch_table)
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
    .sdAddPatchTable(
        sd, xy, patch_vec,
        patch_shape = patch_shape,
        patch_table = patch_table,
        diagnostics = diagnostics,
        ss_log = if (log_iters) patchResult$ss_log
    )
}



#' Run patch-level differential expression on a SpatialData object
#'
#' Runs \code{patchDE()} on the expression data and patch assignments stored
#' in \code{tables(sd)[[table]]} and stores the results in the patch table
#' created by \code{getPatches.sd()}.
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
#' @param patch_table Name of the patch table of \code{sd} created by
#'   \code{getPatches.sd()}, where results are stored. Default
#'   \code{"patches"}.
#' @inheritParams embedCellNeighborhoods.sd
#' @inheritParams patchDE.spe
#'
#' @return The input \code{sd}, with results attached:
#'   \itemize{
#'     \item \code{tables(sd)[[patch_table]]} is rebuilt with one row per gene,
#'       and every earlier result of \code{patchDE.sd()},
#'       \code{patchMetaAnalysis.sd()} and \code{moranTest.sd()} is removed
#'       from it. The results of \code{getPatches.sd()} are kept. For each
#'       predictor \code{<p>} of \code{patchDE()} (one per entry of
#'       \code{predictor_cols}, or per model-matrix column for factors) it
#'       gains the genes-by-patches assays \code{est_<p>}, \code{se_<p>} and
#'       \code{pval_<p>}, in that order, plus \code{resid_mse} when
#'       \code{resid_mse = TRUE}. Patches without a fit, such as those
#'       whose cells were all set to \code{NA}, get \code{NA}.
#'       \code{metadata(tables(sd)[[patch_table]])$SpaceMosaic$patchDE}
#'       records the cell table used, the predictors, and the fitted patches
#'       in the order \code{patchDE()} produced, for
#'       \code{patchMetaAnalysis.sd()}.
#'     \item \code{assay(tables(sd)[[table]], residual_assay)}: if
#'       \code{return_residuals = TRUE}, the residual matrix from
#'       \code{patchDE()} (genes by cells).
#'   }
#'
#' @details
#' The expression matrix is extracted from \code{assay(tables(sd)[[table]], assay)} and
#' transposed from genes-by-cells to cells-by-genes before being passed to
#' \code{patchDE()}. Predictor variables are taken from
#' \code{colData(tables(sd)[[table]])}, and patch assignments from
#' \code{colData(tables(sd)[[table]])[[patch_column]]}.
#'
#' Differential expression is run on every non-missing patch. Results are
#' matched to the columns of the patch table by patch ID; fitted patches that
#' are not columns of the patch table are dropped with a warning.
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
#' sd <- patchDE.sd(sd,
#'   predictor_cols = "distance_to_fibroblast"
#' )
#' patches <- table(sd, "patches")
#' assayNames(patches)
#' assay(patches, "pval_distance_to_fibroblast")[1:5, 1:5]
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
    patch_table = "patches",
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

    patch_sce <- .sdPatchTable(sd, patch_table)

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
        de_res <- de_res$de
    }

    # One element per predictor, next to backend extras such as resid_mse.
    predictors <- names(de_res)[vapply(de_res, function(x)
        is.list(x) && "ests" %in% names(x), logical(1))]
    ests <- de_res[[predictors[1]]]$ests

    # One assay per statistic and predictor, estimates first.
    stats <- c(est = "ests", se = "ses", pval = "pvals")
    mats <- list()
    for (stat in names(stats)) {
        for (predictor in predictors) {
            mats[[paste0(stat, "_", predictor)]] <-
                de_res[[predictor]][[stats[[stat]]]]
        }
    }
    if (resid_mse) {
        mats$resid_mse <- de_res$resid_mse
    }

    patch_sce <- .patchTableReset(patch_sce, rownames(ests))
    patch_sce <- .patchTableAssays(patch_sce, mats)
    S4Vectors::metadata(patch_sce)$SpaceMosaic$patchDE <- list(
        table = table,
        predictors = predictors,
        patches = colnames(ests)
    )
    .sdSetTable(sd, patch_sce, patch_table)
}

#' Meta-analyse patch-level differential expression across an embedding
#'
#' Takes the results of \code{patchDE.sd()} from the patch table, derives
#' patch-level attributes from a reduced-dimension embedding stored in
#' \code{tables(sd)[[table]]}, runs \code{patchMetaAnalysis()} to combine the
#' per-patch fits, and stores the results in the patch table.
#'
#' @param sd A \code{SpatialData} object on which \code{patchDE.sd()} has
#'   been run.
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
#' @param patch_table Name of the patch table of \code{sd} holding the
#'   results of \code{patchDE.sd()}, where results are stored. Default
#'   \code{"patches"}.
#' @inheritParams patchMetaAnalysis.spe
#'
#' @return The input \code{sd}, with results attached to
#'   \code{tables(sd)[[patch_table]]}. Earlier results of
#'   \code{patchMetaAnalysis.sd()} are replaced.
#'   \itemize{
#'     \item For each predictor \code{<p>}, the genes-by-patches assays
#'       \code{meta_est_<p>}, \code{meta_se_<p>} and \code{meta_pval_<p>},
#'       and \code{meta_subgroups_<p>} where the backend produced it, as
#'       returned by \code{patchMetaAnalysis()}.
#'     \item \code{reducedDim(, "W")}: the patches-by-attributes matrix from
#'       \code{getPatchAttributes()}. Patches without cells get \code{NA}.
#'     \item \code{metadata()$SpaceMosaic$subgroups_summary}: the output of
#'       \code{summarizeSubgroups()}, if \code{summarize_subgroups = TRUE}.
#'   }
#'
#' @details
#' The per-patch fits are read from the \code{est_<p>}, \code{se_<p>} and
#' \code{pval_<p>} assays of the patch table, restricted to the patches
#' \code{patchDE()} fitted, in the order it produced. Patch attributes are
#' computed from \code{colData(tables(sd)[[table]])[[patch_column]]}; the
#' patches represented there and those fitted by \code{patchDE.sd()} must be
#' the same set; their order need not agree.
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
#' sd <- patchDE.sd(sd,
#'   predictor_cols = "distance_to_fibroblast"
#' )
#' sd <- patchMetaAnalysis.sd(sd, embedding_name = "Z")
#' patches <- table(sd, "patches")
#' assay(patches, "meta_pval_distance_to_fibroblast")[1:5, 1:5]
#'
#' @export
patchMetaAnalysis.sd <- function(
    sd,
    embedding_name = "Z",
    patch_column = "patch",
    summarize_subgroups = FALSE,
    cellmeta_cols = NULL,
    table = NULL,
    patch_table = "patches"
) {

    .require_spatialdataR()
    .checkSd(sd)
    patch_sce <- .sdPatchTable(sd, patch_table)
    de_info <- S4Vectors::metadata(patch_sce)$SpaceMosaic$patchDE
    if (is.null(de_info)) {
        stop("No patchDE.sd() results found in tables(sd)[[patch_table]]; ",
             "run patchDE.sd() first.")
    }
    if (is.null(table)) {
        table <- de_info$table
    }
    table <- .sdTableName(sd, table)
    sce <- .sdTable(sd, table)

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

    # Rebuild the patchDE() output from the assays, with the fitted patches
    # in the order patchDE() produced.
    fitted <- de_info$patches[de_info$patches %in% colnames(patch_sce)]
    de_res <- lapply(stats::setNames(nm = de_info$predictors), function(p) {
        stat <- function(prefix) {
            SummarizedExperiment::assay(
                patch_sce, paste0(prefix, "_", p)
            )[, fitted, drop = FALSE]
        }
        list(pvals = stat("pval"), ests = stat("est"), ses = stat("se"))
    })

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

    # One assay per statistic and predictor, replacing earlier ones.
    stats <- c(ests = "est", ses = "se", pvals = "pval")
    mats <- list()
    for (predictor in names(de_res_meta)) {
        for (stat in names(de_res_meta[[predictor]])) {
            prefix <- if (stat %in% names(stats)) stats[[stat]] else stat
            mats[[paste0("meta_", prefix, "_", predictor)]] <-
                de_res_meta[[predictor]][[stat]]
        }
    }
    patch_sce <- .patchTableAssays(.dropAssays(patch_sce, "meta_"), mats)

    W_full <- matrix(NA_real_, ncol(patch_sce), ncol(W),
                     dimnames = list(colnames(patch_sce), colnames(W)))
    i <- match(rownames(W), colnames(patch_sce))
    W_full[i[!is.na(i)], ] <- W[!is.na(i), , drop = FALSE]
    SingleCellExperiment::reducedDim(patch_sce, "W") <- W_full

    subgroups_summary <- NULL
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

        subgroups_summary <- summarizeSubgroups(
            de_res_meta,
            patch = patch,
            cellmeta = SummarizedExperiment::colData(
                sce
            )[, cellmeta_cols, drop = FALSE]
        )
    }
    S4Vectors::metadata(patch_sce)$SpaceMosaic$subgroups_summary <-
        subgroups_summary

    .sdSetTable(sd, patch_sce, patch_table)
}

#' Moran's I test on a SpatialData object
#'
#' moranTest() method for \code{SpatialData} objects. Extracts residuals from
#' the specified assay of \code{tables(sd)[[table]]} and spatial coordinates from
#' the shapes element it annotates, dispatches to \code{moranTest}, and stores
#' the results in the patch table created by \code{getPatches.sd()}.
#'
#' @param sd A \code{SpatialData} object.
#' @param assay_name Character; name of the assay in \code{tables(sd)[[table]]}
#'   containing the values to test (e.g. Pearson residuals). Default
#'   \code{"residuals"}.
#' @param patch_column Name of the column in `colData(tables(sd)[[table]])` where patch
#'   assignments are stored. Default `"patch"`.
#' @param patch_table Name of the patch table of `sd` created by
#'   `getPatches.sd()`, where results are stored. Default `"patches"`.
#' @inheritParams embedCellNeighborhoods.sd
#' @inheritParams moranTest.spe
#' @return The input `sd`, with the results of `moranTest()` attached to
#'   `tables(sd)[[patch_table]]`, replacing earlier ones:
#'   \itemize{
#'     \item the genes-by-patches assays `moran_observed`, `moran_expected`,
#'       `moran_p_value`, `moran_p_adjusted` (observed and expected Moran's I,
#'       raw and adjusted empirical p-values) and `moran_status`. A status of
#'       `"ok"` indicates a completed test; non-testable inputs such as
#'       constant or non-finite residuals have missing statistics and an
#'       explanatory status. Genes or patches that were not tested are `NA`.
#'     \item the `colData` columns `moran_n_cells` and `moran_effective_k`:
#'       the number of cells and the number of neighbors used per patch.
#'     \item `metadata()$SpaceMosaic$moranTest`: the cell table and assay
#'       used and the graph, permutation and adjustment settings.
#'   }
#'   If the patch table has no rows yet (`patchDE.sd()` has not been run),
#'   it is given one row per tested gene.
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
#' sd <- patchDE.sd(sd,
#'   predictor_cols = "distance_to_fibroblast",
#'   return_residuals = TRUE
#' )
#' sd <- moranTest.sd(sd)
#' patches <- table(sd, "patches")
#' assay(patches, "moran_p_adjusted")[1:5, 1:5]
#' @export
moranTest.sd <- function(sd, assay_name = 'residuals', patch_column = "patch", k = 10L,
                         n_permutations = 999L,
                         alternative = c("greater", "less", "two.sided"),
                         p_adjust_method = "BH",
                         adjustment_scope = c("global", "patch", "gene"),
                         table = 1,
                         patch_table = "patches") {

    .require_spatialdataR()
    .checkSd(sd)
    table <- .sdTableName(sd, table)
    sce <- .sdTable(sd, table)
    alternative <- match.arg(alternative)
    adjustment_scope <- match.arg(adjustment_scope)

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

    patch_sce <- .sdPatchTable(sd, patch_table)

    res <- moranTest(residuals = t(assay(sce, assay_name)),
                     xy = .sdCoords(sd, sce),
                     patch = colData(sce)[[patch_column]],
                     k = k,
                     n_permutations = n_permutations,
                     alternative = alternative,
                     p_adjust_method = p_adjust_method,
                     adjustment_scope = adjustment_scope)

    if (nrow(patch_sce) == 0) {
        patch_sce <- .patchTableReset(patch_sce, unique(res$gene))
    }

    # Long results (one row per gene and patch) -> genes x patches matrices
    genes <- unique(res$gene)
    patches <- unique(res$patch)
    idx <- cbind(match(res$gene, genes), match(res$patch, patches))
    mats <- list()
    for (stat in c("observed", "expected", "p_value", "p_adjusted", "status")) {
        m <- matrix(NA, length(genes), length(patches),
                    dimnames = list(genes, patches))
        storage.mode(m) <- storage.mode(res[[stat]])
        m[idx] <- res[[stat]]
        mats[[paste0("moran_", stat)]] <- m
    }
    patch_sce <- .patchTableAssays(.dropAssays(patch_sce, "moran_"), mats)

    # Per-patch values
    first <- res[!duplicated(res$patch), , drop = FALSE]
    i <- match(colnames(patch_sce), first$patch)
    patch_sce$moran_n_cells <- first$n_cells[i]
    patch_sce$moran_effective_k <- first$effective_k[i]

    S4Vectors::metadata(patch_sce)$SpaceMosaic$moranTest <- list(
        table = table,
        assay_name = assay_name,
        k = k,
        n_permutations = n_permutations,
        alternative = alternative,
        p_adjust_method = p_adjust_method,
        adjustment_scope = adjustment_scope
    )
    .sdSetTable(sd, patch_sce, patch_table)
}

