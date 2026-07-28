#' Get patch polygons for visualization
#'
#' Computes the convex hull of every patch containing at least three cells and
#' returns the vertices in a tidy data frame ready for
#' `ggplot2::geom_polygon()`. Optional patch-level annotations (for example a
#' cluster, group, size, or a user-computed diagnostic) are joined to every
#' vertex by patch identifier.
#'
#' @param xy Cells' xy positions. The first two columns are used as coordinates.
#' @param patch Vector of patch assignments, aligned to the rows of `xy`.
#' @param patch_data Optional data frame containing one row per patch. It must
#'   contain a `patch` column, or have patch identifiers as row names. Any
#'   additional columns are included in the result.
#' @return A data frame with columns `x`, `y`, and `patch`, followed by any
#'   columns supplied in `patch_data`. Patches with fewer than three cells are
#'   omitted.
#' @examples
#' xy <- matrix(c(0, 0, 1, 0, 0, 1, 3, 3, 4, 3, 3, 4), ncol = 2, byrow = TRUE)
#' patch <- rep(c("A", "B"), each = 3)
#' info <- data.frame(patch = c("A", "B"), cluster = c("Cluster 1", "Cluster 2"))
#' getPatchPolys(xy, patch, patch_data = info)
#' @export
getPatchPolys <- function(xy, patch, patch_data = NULL) {
  xy <- as.matrix(xy)
  if (ncol(xy) < 2L) stop("xy must contain at least two coordinate columns.")
  if (!is.numeric(xy[, 1]) || !is.numeric(xy[, 2])) {
    stop("The first two columns of xy must be numeric.")
  }
  if (nrow(xy) != length(patch)) {
    stop("length(patch) must equal nrow(xy).")
  }

  patch <- as.character(patch)
  patch_ids <- unique(patch[!is.na(patch)])
  polygons <- lapply(patch_ids, function(id) {
    idx <- which(!is.na(patch) & patch == id)
    if (length(idx) < 3L) return(NULL)

    points <- xy[idx, 1:2, drop = FALSE]
    hull <- grDevices::chull(points[, 1], points[, 2])
    data.frame(
      x = points[hull, 1],
      y = points[hull, 2],
      patch = id,
      stringsAsFactors = FALSE
    )
  })
  polygons <- Filter(Negate(is.null), polygons)
  if (length(polygons) == 0L) {
    out <- data.frame(x = numeric(), y = numeric(), patch = character())
  } else {
    out <- do.call(rbind, polygons)
    rownames(out) <- NULL
  }

  if (is.null(patch_data)) return(out)
  patch_data <- .normalizePatchData(patch_data)
  annotation_names <- setdiff(names(patch_data), "patch")
  if (length(annotation_names) == 0L) return(out)
  if (any(annotation_names %in% c("x", "y"))) {
    stop("patch_data cannot contain columns named 'x' or 'y'.")
  }

  annotation_rows <- match(out$patch, patch_data$patch)
  out <- cbind(out, patch_data[annotation_rows, annotation_names, drop = FALSE])
  rownames(out) <- NULL
  out
}

.normalizePatchData <- function(patch_data) {
  patch_data <- as.data.frame(patch_data, stringsAsFactors = FALSE)
  if (!"patch" %in% names(patch_data)) {
    ids <- rownames(patch_data)
    default_ids <- identical(ids, as.character(seq_len(nrow(patch_data))))
    if (is.null(ids) || default_ids) {
      stop("patch_data must contain a 'patch' column or patch identifiers as row names.")
    }
    patch_data <- data.frame(
      patch = ids,
      patch_data,
      row.names = NULL,
      check.names = FALSE
    )
  }

  patch_data$patch <- as.character(patch_data$patch)
  if (anyNA(patch_data$patch) || any(!nzchar(patch_data$patch))) {
    stop("Patch identifiers in patch_data cannot be missing or empty.")
  }
  if (anyDuplicated(patch_data$patch)) {
    stop("patch_data must contain at most one row per patch.")
  }
  patch_data
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
    polygons <- getPatchPolys(xy, patches)
    for (id in unique(polygons$patch)) {
      p <- polygons[polygons$patch == id, , drop = FALSE]
      graphics::polygon(p$x, p$y, border = "black", lwd = 1.5)
    }
  }
}
