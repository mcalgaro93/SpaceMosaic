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


#' Plot patch assignments at each iteration
#'
#' Given the output of \code{getPatches(log_iters = TRUE)} and the xy coordinates,
#' draws a series of plots showing patch evolution across iterations.
#'
#' @param xy Matrix of cells' xy positions (same as passed to getPatches).
#' @param result Output of \code{getPatches(..., log_iters = TRUE)}.
#' @param iters Which iterations to plot. Default NULL = all.
#' @param cols Optional vector of colors (one per unique patch ID). If NULL,
#'   uses a default palette.
#' @param cex Point size. Default 0.5.
#' @param ask If TRUE, prompt between plots (default: interactive sessions only).
#' @export
plotPatchIterations <- function(xy, result, iters = NULL,
                                cols = NULL, cex = 0.5, ask = interactive()) {
  stopifnot(is.list(result), "membership_log" %in% names(result))
  membership_log <- result$membership_log
  n_iters <- ncol(membership_log)

  if (is.null(iters)) iters <- seq_len(n_iters)
  iters <- iters[iters >= 1 & iters <= n_iters]

  ## build color palette if not supplied
  all_ids <- sort(unique(as.vector(membership_log[, iters])))
  all_ids <- all_ids[!is.na(all_ids)]
  if (is.null(cols)) {
    np <- length(all_ids)
    cols <- grDevices::hcl.colors(np, palette = "Set 3")
  }
  names(cols) <- all_ids[seq_along(cols)]

  oldpar <- graphics::par(ask = ask)
  on.exit(graphics::par(oldpar))

  for (it in iters) {
    patches <- membership_log[, it]
    patch_factor <- as.numeric(as.factor(patches))
    cellcols <- cols[as.character(patches)]
    cellcols[is.na(patches)] <- "grey80"

    graphics::plot(xy, pch = 16, cex = cex, col = cellcols,
                   main = paste0("Iteration ", it),
                   xlab = "", ylab = "", asp = 1)
    polys <- getPatchPolys(xy, patches)
    for (p in polys) {
      if (!is.null(p)) graphics::polygon(p[, 1], p[, 2], border = "black", lwd = 1.5)
    }
  }
}