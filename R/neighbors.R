#' Create spatial network from N nearest neighbors
#'
#' For each cell identify \code{N} nearest neighbors in Euclidean space and
#' create an edge between them in graph structure.
#'
#' Edges will only be created for cells that have the same \code{subset} value,
#' usually the slide column id but could also be a slide plus FOV id to only
#' create edges within an FOV.
#'
#' @param x Spatial coordinate.
#' @param y Spatial coordinate.
#' @param N Positive integer giving the number of nearest neighbors. Every
#'   subset must contain at least `N + 1` cells.
#' @param subset Either a scalar or a vector with the same length as `x` and
#'   `y`. Only cells sharing a subset value become neighbors.
#' @return Sparse adjacency matrix with distances.
#' @importFrom data.table data.table rbindlist
#' @importFrom spatstat.geom nnwhich nndist
#' @importFrom Matrix sparseMatrix
#' @export
nearestNeighborGraph <- function(x, y, N, subset = 1) {
  if (length(x) != length(y)) {
    stop("x and y must have the same length.")
  }
  if (length(N) != 1L || !is.finite(N) || N < 1 || N != floor(N)) {
    stop("N must be a positive integer.")
  }
  N <- as.integer(N)
  if (length(subset) == 1L) {
    subset <- rep(subset, length(x))
  }
  if (length(subset) != length(x)) {
    stop("subset must be a scalar or have the same length as x and y.")
  }
  if (anyNA(x) || anyNA(y) || anyNA(subset)) {
    stop("x, y, and subset must not contain missing values.")
  }

  DT <- data.table::data.table(x = x, y = y, subset = subset)
  nearestNeighbor <- function(i) {
    subset_dt <- DT[DT[["subset"]] == i]
    idx <- which(DT[["subset"]] == i)
    if (nrow(subset_dt) <= N) {
      stop(
        sprintf(
          "Each subset must contain at least N + 1 cells; subset '%s' contains %d.",
          as.character(i),
          nrow(subset_dt)
        )
      )
    }
    ndist <- spatstat.geom::nndist(subset_dt[, list(x, y)], k = 1:N)
    nwhich <- spatstat.geom::nnwhich(subset_dt[, list(x, y)], k = 1:N)
    ij <- data.table::data.table(
      i = rep(idx, times = N),
      j = idx[as.vector(nwhich)],
      x = as.vector(ndist)
    )
    return(ij)
  }
  ij <- data.table::rbindlist(lapply(unique(subset), nearestNeighbor))
  adj.m <- Matrix::sparseMatrix(i = ij$i, j = ij$j, x = ij$x,
                                dims = c(nrow(DT), nrow(DT)))
  return(adj.m)
}


#' Column-wise neighbor means
#'
#' For each cell, compute the mean of each column of \code{x} across its
#' neighbors.
#'
#' @param x A numeric matrix.
#' @param neighbors A sparse adjacency matrix.
#' @return Matrix of the same dimensions as \code{x}. Rows without neighbors
#'   are returned as zero rows.
#' @importFrom Matrix rowSums Diagonal
#' @export
neighbor_colMeans <- function(x, neighbors) {
  if (nrow(neighbors) != ncol(neighbors)) {
    stop("neighbors must be a square adjacency matrix.")
  }
  if (nrow(x) != nrow(neighbors)) {
    stop("nrow(x) must equal nrow(neighbors).")
  }

  binary_neighbors <- neighbors
  binary_neighbors@x <- rep(1, length(binary_neighbors@x))
  neighbor_counts <- Matrix::rowSums(binary_neighbors)
  inverse_counts <- numeric(length(neighbor_counts))
  has_neighbors <- neighbor_counts > 0
  inverse_counts[has_neighbors] <- 1 / neighbor_counts[has_neighbors]
  normalized_neighbors <-
    Matrix::Diagonal(x = inverse_counts) %*% binary_neighbors
  out <- normalized_neighbors %*% x
  return(out)
}


#' Embed cellular neighborhoods from single cell embeddings and positions
#'
#' Creates a neighborhood embedding by averaging a cell embedding matrix over
#' spatial neighbor networks at multiple scales.
#'
#' @param mat Single cell embeddings matrix (cells x features).
#' @param xy 2-column matrix of cells' positions.
#' @param ks Vector giving the number of nearest neighbors for each scale.
#'   Default \code{c(5, 50)}.
#' @param tissue Optional vector giving tissue IDs to prevent cross-tissue
#'   neighbor edges. Default NULL.
#' @return A matrix of cellular neighborhood embeddings,
#'   dimensions n cells x (ncol(mat) * length(ks)).
#' @export
embedCellNeighborhoods <- function(mat, xy, ks = c(5, 50), tissue = NULL) {
  if (is.null(tissue)) tissue <- 1
  neighborslist <- list()
  for (k in ks) {
    neighborslist[[paste0("k", k)]] <-
      nearestNeighborGraph(xy[, 1], xy[, 2], N = k, subset = tissue)
  }
  env <- c()
  for (name in names(neighborslist)) {
    env <- cbind(env, as.matrix(neighbor_colMeans(mat, neighbors = neighborslist[[name]])))
    newinds <- (ncol(env) - ncol(mat) + 1):ncol(env)
    colnames(env)[newinds] <- paste0(name, colnames(env)[newinds])
  }
  return(as.matrix(env))
}
