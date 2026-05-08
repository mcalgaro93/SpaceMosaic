#' Find rough polygon boundaries of patches for visualizations
#' @param xy Cells' xy positions
#' @param patch Vector of patch assignments, aligned to the rows of xy
#' @return A named list of polygons, one per patch
#' @export
getPatchPolys <- function(xy, patch) {
  polys <- list()
  cluster_levels <- unique(patch)
  for (i in seq_along(unique(patch))) {
    k <- cluster_levels[i]
    idx <- which(patch == k)
    
    # Only attempt hull if >= 3 points
    if (length(idx) >= 3) {
      pts_k <- xy[idx, , drop = FALSE]       # M_k × 2 matrix of points in cluster k
      hull_indices <- chull(pts_k)           # indices (1..M_k) along the convex hull
      polys[[i]]  <- pts_k[hull_indices, ]       # hull vertices, in order
      names(polys)[i] <- unique(patch)[i]
    }
  }
  return(polys)
}