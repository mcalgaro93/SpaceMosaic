#' Convert a SpatialExperiment object to a SpatialData object
#'
#' Stores the spatial coordinates of a \code{SpatialExperiment} as a points
#' shapes element and everything else (assays, \code{colData},
#' \code{reducedDims}, ...) as a \code{SingleCellExperiment} table annotating
#' that shapes element. Every image in \code{imgData(spe)} is added as an
#' image element.
#'
#' Images are stored at their own resolution, in pixels, without a coordinate
#' transformation. The cell coordinates are in full-resolution pixels, so an
#' image only lines up with the cells if its \code{scaleFactor} is 1;
#' otherwise multiply the image coordinates by \code{1 / scaleFactor} to
#' match them to the cells.
#'
#' @param spe A SpatialExperiment object.
#' @param shape_name Name of the shapes element that holds the cell
#'   coordinates. Default \code{"shape1"}.
#' @param table_name Name of the table element that holds the
#'   SingleCellExperiment. Default \code{"table"}.
#' @param add_cell_id If TRUE, stores the cell IDs (colnames) in the shapes
#'   so each table row can be matched to its geometry. Default TRUE.
#' @param image_names Names of the image elements, one per row of
#'   \code{imgData(spe)}. Default \code{NULL} uses
#'   \code{"<sample_id>_<image_id>"}, e.g. \code{"section1_lowres"}.
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
spe_to_spatialdata <- function(spe,
                               shape_name  = "shape1",
                               table_name  = "table",
                               add_cell_id = TRUE,
                               image_names = NULL) {

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
  img_data <- SpatialExperiment::imgData(spe)
  if (is.null(image_names)) {
    image_names <- paste(img_data$sample_id, img_data$image_id, sep = "_")
  }
  if (length(image_names) != nrow(img_data) || anyDuplicated(image_names)) {
    stop("`image_names` must have one unique name per row of imgData(spe) (",
         nrow(img_data), ").")
  }
  coords <- SpatialExperiment::spatialCoords(spe)

  # ---- Table: keep only the SingleCellExperiment part ----
  sce <- methods::as(spe, "SingleCellExperiment")

  # ---- Coordinates -> duckspatial points ----
  coords_df <- as.data.frame(coords)
  colnames(coords_df) <- c("x", "y")
  # Explicit instance IDs: the shapes are stored in DuckDB, which does not
  # preserve row order, so cells cannot be matched to geometries by position.
  coords_df$instance_id <- seq_len(nrow(coords_df))
  if (add_cell_id) {
    coords_df$cell_id <- colnames(spe)
  }

  # .sdShape() avoids spatialdataR's shared DuckDB table (see its comment).
  coords_sf <- sf::st_as_sf(coords_df, coords = c("x", "y"), crs = NA)
  sd_shapes <- .sdShape(coords_sf)

  # ---- Build SpatialData and attach the table ----
  shapes_list <- stats::setNames(list(sd_shapes), shape_name)
  sd <- spatialdataR::SpatialData(shapes = shapes_list)
  sd <- spatialdataR::setTable(sd, shape_name, sce, name = table_name)
  # setTable() copies the shapes' instance IDs onto the table by position;
  # reset them so table column i is matched to the geometry of cell i.
  sce <- spatialdataR::tables(sd)[[table_name]]
  spatialdataR::instances(sce) <- seq_len(ncol(sce))
  sd <- .sdSetTable(sd, sce, table_name)

  # ---- Images (if any) ----
  for (i in seq_len(nrow(img_data))) {
    # imgRaster() reads stored, loaded and remote images alike.
    img <- SpatialExperiment::imgRaster(
      spe, sample_id = img_data$sample_id[i], image_id = img_data$image_id[i]
    )
    # Raster (y, x) of colours -> numeric array (c, y, x) in [0, 1].
    rgb <- grDevices::col2rgb(as.matrix(img)) / 255
    img <- array(rgb, dim = c(3, nrow(img), ncol(img)))
    spatialdataR::images(sd)[[image_names[i]]] <-
      spatialdataR::SpatialDataImage(img)
  }
  sd
}
