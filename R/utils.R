#' Adapted from InSituCor:::nearestNeighborGraph
#' For each cell identify \code{N} nearest neighbors in Euclidean space and
#' create an edge between them in graph structure
#'
#' Edges will only be created for cells that have the same \code{subset} value,
#' usually the slide column id but could also be a slide plus FOV id to only
#' create edges within an FOV.
#'
#' @param x spatial coordinate
#' @param y spatial coordinate
#' @param N number of nearest neighbors
#'
#' @return sparse adjacency matrix with distances
#' @importFrom spatstat.geom nnwhich
#' @importFrom spatstat.geom nndist
#' @importFrom Matrix sparseMatrix
.nearestNeighborGraph <- function(xy, N) {
    ndist <- spatstat.geom::nndist(xy, k=1:N)
    nwhich <- spatstat.geom::nnwhich(xy, k=1:N)
    A <- sparseMatrix(i = rep(1:nrow(xy), times=N), j = as.vector(nwhich),
                      x = as.vector(ndist), dims = c(nrow(xy), nrow(xy)))
  A
}


#' Adapted from InSituCor:::nearestNeighborGraph
#' for each cell, tabulate the distinct values of x over its neighbors:
#' @param x A vector of categorical values
#' @param neighbors A (probably sparse) adjacency matrix
.neighbor_tabulate <- function(x, neighbors) {
  uniquevals <- unique(x)
  sapply(uniquevals, function(val){
    vec <- (x == val) * 1
    .neighbor_sum(vec, neighbors)
  })
}

#' Adapted from InSituCor:::nearestNeighborGraph
#' for each cell, get the sum of x's values over its neighbors:
#' @param x A numeric vector
#' @param neighbors A (probably sparse) adjacency matrix
#' @importFrom Matrix rowSums t
.neighbor_sum <- function(x, neighbors) {
  rowSums(t(t(1*(neighbors != 0)) * x))
}

#' Adapted from InSituCor:::nearestNeighborGraph
#' for each cell, get the mean of x's values over its neighbors:
#' @param x A numeric vector
#' @param neighbors A (probably sparse) adjacency matrix
#' @importFrom Matrix rowSums
.neighbor_mean <- function(x, neighbors) {
  .neighbor_sum(x, neighbors) / pmax(rowSums(neighbors != 0), 1)
}

