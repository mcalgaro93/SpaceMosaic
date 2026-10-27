#' Convert a SpatialExperiment object to a SpatialData object
#'
#' Stores the spatial coordinates of a \code{SpatialExperiment} as a points
#' shapes element and everything else (assays, \code{colData},
#' \code{reducedDims}, ...) as a \code{SingleCellExperiment} table annotating
#' that shapes element. If \code{spe} has an image in \code{imgData(spe)}, the
#' first one is added as an image element.
#'
#' @param spe A SpatialExperiment object.
#' @param sd Optional, a SpatialData object where to add the SpatialExperiment
#'   object content. Default \code{NULL} creates a new one.
#' @param shape_name Name of the shapes element that holds the cell
#'   coordinates. Default \code{"shape1"}.
#' @param table_name Name of the table element that holds the
#'   SingleCellExperiment. Default \code{"table"}.
#' @param add_cell_id If TRUE, stores the cell IDs (colnames) in the shapes
#'   so each table row can be matched to its geometry. Default TRUE.
#' @return A SpatialData object.
#'
#' @seealso \code{\link{getPatches.sd}} and the other \code{.sd} wrappers.
#'
#' @examplesIf requireNamespace("spatialdataR", quietly = TRUE)
#' library(SpatialExperiment)
#' spe <- readRDS(system.file("extdata", "cosmx_carcinoma.rds", package = "SpaceMosaic"))
#' sd <- spe_to_spatialdata(spe)
#' sd
#' @export
spe_to_spatialdata <- function(spe, sd = NULL,
                               shape_name  = "shape1",
                               table_name  = "table",
                               add_cell_id = TRUE) {

  # ---- Input checks ----
  for (pkg in c("spatialdataR", "sf", "duckspatial")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required. Please install it.")
    }
  }
  if (!methods::is(spe, "SpatialExperiment")) {
    stop("`spe` must be a SpatialExperiment object.")
  }
  if (add_cell_id && is.null(colnames(spe))) {
    stop("`spe` has no colnames; set them or use `add_cell_id = FALSE`.")
  }
  coords <- SpatialExperiment::spatialCoords(spe)

  # ---- Table: keep only the SingleCellExperiment part ----
  sce <- methods::as(spe, "SingleCellExperiment")

  # ---- Coordinates -> duckspatial points ----
  coords_df <- as.data.frame(coords)
  colnames(coords_df) <- c("x", "y")
  if (add_cell_id) {
    coords_df$cell_id <- colnames(spe)
  }

  # .sdShape() avoids spatialdataR's shared DuckDB table (see its comment).
  coords_sf <- sf::st_as_sf(coords_df, coords = c("x", "y"), crs = NA)
  sd_shapes <- .sdShape(coords_sf)

  # ---- Build SpatialData and attach the table ----
  if (is.null(sd)) {
    shapes_list <- stats::setNames(list(sd_shapes), shape_name)
    sd <- spatialdataR::SpatialData(shapes = shapes_list)
  } else {
    sd <- spatialdataR::`shape<-`(sd, shape_name, value = sd_shapes)
  }
  sd <- spatialdataR::setTable(sd, shape_name, sce, name = table_name)

  # ---- Image (if any) ----
  if (nrow(SpatialExperiment::imgData(spe)) > 0) {
    if (!requireNamespace("png", quietly = TRUE)) {
      stop("Package 'png' is required to transfer the image of `spe`.")
    }
    img <- png::readPNG(SpatialExperiment::imgSource(SpatialExperiment::getImg(spe)))
    spatialdataR::images(sd)$image1 <- spatialdataR::SpatialDataImage(img)
  }
  sd
}
