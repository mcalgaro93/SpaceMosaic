
#' Assign cells to patches using elliptical Gaussian assignments
#'
#' Iterative EM-like algorithm:
#'   (a) Estimate per-patch centroid + covariance from xy.
#'   (b) Assign cells using spatial fit + X-diversity boost + Z-penalty + hunger.
#'   (c) Re-estimate per-patch centroid + covariance.
#'   (d) Assign cells using spatial fit only (no X/Z/hunger).
#'   (e) Contiguity check: set orphans to NA.
#'
#' @param xy Matrix of cells' xy positions. Must have rownames.
#' @param X Numeric vector or matrix of design variables, aligned to rows of xy.
#'   If a vector, treated as a single column. Each column is scaled to unit SD.
#' @param npatches Number of patches to create.
#' @param Z Optional matrix of per-cell context embeddings (cells x features).
#' @param alpha Weight of the Z penalty. Default 0.5.
#' @param beta Weight of the X diversity boost. Default 1.
#' @param hunger_weight Controls how aggressively low-variance patches grab cells.
#'   0 = uniform hunger, 1 = hunger proportional to 1/totvar. Default 0.5.
#' @param max_elongation Maximum ratio of largest to smallest eigenvalue of a
#'   patch covariance matrix. Constrains how narrow patches can be. Default 4.
#' @param max_radius Maximum Mahalanobis radius for assignment. Cells beyond
#'   this threshold get zero spatial score for the patch. Default 3.
#' @param n_iters Number of outer iterations. Default 15.
#' @param log_iters If TRUE, return a list with patch assignments plus
#'   per-iteration diagnostics (SS per patch and membership). Default TRUE.
#' @param verbose Show progress. Default TRUE.
#' @return If log_iters = FALSE, a named vector of patch assignments.
#'   If log_iters = TRUE, a list with:
#'   \item{patch}{Named vector of final patch assignments.}
#'   \item{ss_log}{Matrix (npatches x n_iters) of per-patch sum-of-squares of X.}
#'   \item{membership_log}{Matrix (n_cells x n_iters) of patch assignments per iteration.}
#' @export
getPatches <- function(xy, X, npatches,
                             Z = NULL,
                             alpha = 0.5,
                             beta = 1,
                             hunger_weight = 0.5,
                             max_elongation = 4,
                             max_radius = 3,
                             n_iters = 15,
                             log_iters = TRUE,
                             verbose = TRUE) {

  ## coerce X to matrix
  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  X <- as.matrix(X)
  stopifnot(nrow(xy) == nrow(X))
  stopifnot(!is.null(rownames(xy)))
  n <- nrow(xy)

  ## scale X column-wise
  X <- scale(X, center = FALSE, scale = apply(X, 2, sd, na.rm = TRUE))

  ## build contiguity graph: intersection of Delaunay + 10-NN
  spatial_nn <- .buildContiguityGraph(xy)

  ## prepare Z
  z_info <- .prepareZ(Z, xy)

  ## initialize with kmeans on xy
  patch <- .initKmeans(xy, npatches)

  ## set up iteration logs
  if (log_iters) {
    membership_log <- matrix(NA_character_, nrow = n, ncol = n_iters)
    ss_log <- matrix(NA_real_, nrow = npatches, ncol = n_iters)
    rownames(ss_log) <- as.character(seq_len(npatches))
  }

  if (verbose) cli::cli_progress_bar("Ellipse iterations", total = n_iters)

  for (iter in seq_len(n_iters)) {

    ## (a) estimate ellipse parameters
    params <- .estimateEllipses(xy, patch, max_elongation)

    ## (b) assign using spatial + X + Z + hunger
    patch <- .assignCells(xy, X, patch, params,
                          z_info = z_info, alpha = alpha, beta = beta,
                          hunger_weight = hunger_weight,
                          max_radius = max_radius,
                          use_xz = TRUE)

    ## (c) re-estimate ellipse parameters
    params <- .estimateEllipses(xy, patch, max_elongation)

    ## (d) assign using spatial only
    patch <- .assignCells(xy, X, patch, params,
                          z_info = z_info, alpha = alpha, beta = beta,
                          hunger_weight = hunger_weight,
                          max_radius = max_radius,
                          use_xz = FALSE)

    ## (e) contiguity check: set orphans to NA
    patch <- .dropOrphans(patch, spatial_nn)

    ## log iteration state
    if (log_iters) {
      membership_log[, iter] <- patch
      pnames <- unique(patch[!is.na(patch)])
      for (p in pnames) {
        cells <- which(patch == p)
        if (length(cells) > 1) {
          ss_log[p, iter] <- sum(scale(X[cells, , drop = FALSE], scale = FALSE)^2)
        } else {
          ss_log[p, iter] <- 0
        }
      }
    }

    if (verbose) cli::cli_progress_update()
  }

  if (verbose) cli::cli_progress_done()

  names(patch) <- rownames(xy)

  if (log_iters) {
    colnames(membership_log) <- paste0("iter", seq_len(n_iters))
    rownames(membership_log) <- rownames(xy)
    colnames(ss_log) <- paste0("iter", seq_len(n_iters))
    return(list(patch = patch, ss_log = ss_log, membership_log = membership_log))
  }
  return(patch)
}


#' Build contiguity graph: intersection of Delaunay triangulation and 10-NN
#' @param xy n x 2 coordinate matrix.
#' @return Sparse adjacency matrix (n x n).
.buildContiguityGraph <- function(xy) {
  dd <- deldir::deldir(xy[, 1], xy[, 2])
  edges <- dd$delsgs
  delaunay_nn <- Matrix::sparseMatrix(
    i = c(edges$ind1, edges$ind2),
    j = c(edges$ind2, edges$ind1),
    x = rep(1, 2 * nrow(edges)),
    dims = c(nrow(xy), nrow(xy))
  )
  knn_graph <- nearestNeighborGraph(xy[, 1], xy[, 2], N = 10)
  delaunay_nn * knn_graph
}


#' Prepare Z: spatially smooth and compute z_scale
#' @param Z Raw Z matrix (cells x features), or NULL.
#' @param xy n x 2 coordinate matrix.
#' @return List with Z_smooth and z_scale, or list(NULL, NULL).
.prepareZ <- function(Z, xy) {
  if (is.null(Z)) return(list(Z_smooth = NULL, z_scale = NULL))
  Z <- as.matrix(Z)
  stopifnot(nrow(Z) == nrow(xy))
  z_neighbors <- nearestNeighborGraph(xy[, 1], xy[, 2], N = 10)
  Z_smooth <- as.matrix(neighbor_colMeans(Z, z_neighbors))
  ## z_scale: mean Z-distance between spatial neighbors
  zn_summary <- Matrix::summary(z_neighbors)
  diffs <- Z_smooth[zn_summary$i, , drop = FALSE] - Z_smooth[zn_summary$j, , drop = FALSE]
  z_scale <- mean(sqrt(rowSums(diffs^2)))
  if (z_scale == 0) z_scale <- 1
  list(Z_smooth = Z_smooth, z_scale = z_scale)
}


#' Initialize patches with kmeans on xy
#' @param xy n x 2 coordinate matrix.
#' @param npatches Number of patches.
#' @return Character vector of patch IDs.
.initKmeans <- function(xy, npatches) {
  km <- stats::kmeans(xy, centers = npatches, nstart = 5, iter.max = 50)
  as.character(km$cluster)
}


#' Estimate per-patch centroid and covariance from xy
#' @param xy n x 2 coordinate matrix.
#' @param patch Character vector of patch assignments (may contain NA).
#' @param max_elongation Max eigenvalue ratio for covariance regularization.
#' @return List with components:
#'   \item{centroids}{Named list of 2-vectors.}
#'   \item{covmats}{Named list of 2x2 covariance matrices.}
#'   \item{inv_covmats}{Named list of 2x2 inverse covariance matrices.}
#'   \item{log_det}{Named numeric vector of log-determinants.}
.estimateEllipses <- function(xy, patch, max_elongation) {
  pnames <- sort(unique(patch[!is.na(patch)]))
  centroids <- list()
  covmats <- list()
  inv_covmats <- list()
  log_det <- numeric(length(pnames))
  names(log_det) <- pnames

  for (p in pnames) {
    cells <- which(patch == p)
    centroids[[p]] <- colMeans(xy[cells, , drop = FALSE])
    if (length(cells) < 3) {
      ## too few cells: use identity scaled to global variance
      covmats[[p]] <- diag(2) * mean(apply(xy, 2, var))
    } else {
      covmats[[p]] <- stats::cov(xy[cells, , drop = FALSE])
    }
    ## regularize elongation
    covmats[[p]] <- .regularizeCov(covmats[[p]], max_elongation)
    inv_covmats[[p]] <- solve(covmats[[p]])
    log_det[p] <- log(det(covmats[[p]]))
  }
  list(centroids = centroids, covmats = covmats,
       inv_covmats = inv_covmats, log_det = log_det,
       pnames = pnames)
}


#' Regularize a 2x2 covariance matrix to limit elongation
#' @param S 2x2 covariance matrix.
#' @param max_elongation Max ratio of eigenvalues.
#' @return Regularized 2x2 covariance matrix.
.regularizeCov <- function(S, max_elongation) {
  e <- eigen(S, symmetric = TRUE)
  vals <- pmax(e$values, 1e-10)
  ratio <- vals[1] / vals[2]
  if (ratio > max_elongation) {
    vals[2] <- vals[1] / max_elongation
  }
  e$vectors %*% diag(vals) %*% t(e$vectors)
}


#' Assign cells to patches based on scores
#'
#' Score = spatial_score + beta * X_boost - alpha * Z_penalty + log(hunger)
#' When use_xz = FALSE, only spatial_score is used.
#'
#' @param xy n x 2 coordinate matrix.
#' @param X Numeric matrix of design variables (n x p).
#' @param patch Current patch assignments.
#' @param params Output of .estimateEllipses.
#' @param z_info Output of .prepareZ.
#' @param alpha Z penalty weight.
#' @param beta X diversity boost weight.
#' @param hunger_weight Hunger interpolation weight.
#' @param max_radius Max Mahalanobis distance for assignment.
#' @param use_xz Whether to include X, Z, and hunger terms.
#' @return Updated patch vector.
.assignCells <- function(xy, X, patch, params,
                         z_info, alpha, beta, hunger_weight,
                         max_radius, use_xz) {
  n <- nrow(xy)
  pnames <- params$pnames
  np <- length(pnames)

  ## compute Mahalanobis distances & spatial scores (vectorized)
  scores <- matrix(-Inf, nrow = n, ncol = np)
  colnames(scores) <- pnames

  for (j in seq_along(pnames)) {
    p <- pnames[j]
    d <- xy - matrix(params$centroids[[p]], nrow = n, ncol = 2, byrow = TRUE)
    mah2 <- rowSums((d %*% params$inv_covmats[[p]]) * d)
    spatial <- -0.5 * mah2 - 0.5 * params$log_det[p]
    ## mask cells beyond max_radius
    spatial[mah2 > max_radius^2] <- -Inf
    scores[, j] <- spatial
  }

  if (use_xz) {
    ## X diversity boost: sum of squared deviations from patch mean across columns
    patch_mean_X <- list()
    for (p in pnames) {
      cells <- which(patch == p)
      patch_mean_X[[p]] <- colMeans(X[cells, , drop = FALSE], na.rm = TRUE)
    }
    for (j in seq_along(pnames)) {
      p <- pnames[j]
      dX <- X - matrix(patch_mean_X[[p]], nrow = n, ncol = ncol(X), byrow = TRUE)
      scores[, j] <- scores[, j] + beta * rowSums(dX^2)
    }

    ## Z penalty
    if (!is.null(z_info$Z_smooth)) {
      patch_mean_Z <- list()
      for (p in pnames) {
        cells <- which(patch == p)
        patch_mean_Z[[p]] <- colMeans(z_info$Z_smooth[cells, , drop = FALSE])
      }
      for (j in seq_along(pnames)) {
        p <- pnames[j]
        dz <- z_info$Z_smooth - matrix(patch_mean_Z[[p]], nrow = n,
                                       ncol = ncol(z_info$Z_smooth), byrow = TRUE)
        z_dist2 <- rowSums(dz^2) / (z_info$z_scale^2)
        scores[, j] <- scores[, j] - alpha * z_dist2
      }
    }

    ## hunger: variance-based (trace of within-patch covariance * n)
    patch_totvar <- numeric(np)
    names(patch_totvar) <- pnames
    for (p in pnames) {
      cells <- which(patch == p)
      if (length(cells) > 1) {
        patch_totvar[p] <- sum(apply(X[cells, , drop = FALSE], 2, var, na.rm = TRUE)) *
                           length(cells)
      }
    }
    mean_totvar <- mean(patch_totvar, na.rm = TRUE)
    if (!is.finite(mean_totvar) || mean_totvar == 0) mean_totvar <- 1
    hunger_denom <- (1 - hunger_weight) * mean_totvar + hunger_weight * patch_totvar
    hunger <- ifelse(hunger_denom > 0, 1 / hunger_denom, 0)
    total_hunger <- sum(hunger)
    if (total_hunger > 0) hunger <- hunger / total_hunger
    ## add log-hunger (avoid log(0))
    log_hunger <- ifelse(hunger > 0, log(hunger), -Inf)
    scores <- scores + matrix(log_hunger, nrow = n, ncol = np, byrow = TRUE)
  }

  ## assign to best patch
  new_patch <- pnames[max.col(scores, ties.method = "first")]
  ## cells with all -Inf scores stay NA
  all_neg_inf <- apply(scores, 1, function(r) all(is.infinite(r) & r < 0))
  new_patch[all_neg_inf] <- NA
  new_patch
}


#' Drop orphan cells that aren't contiguous with the main body of their patch
#' @param patch Character vector of patch assignments.
#' @param spatial_nn Sparse adjacency matrix (n x n).
#' @return Updated patch vector with orphans set to NA.
.dropOrphans <- function(patch, spatial_nn) {
  n <- length(patch)
  trip <- Matrix::summary(spatial_nn)
  same_patch <- (!is.na(patch[trip$i])) & (!is.na(patch[trip$j])) &
                (patch[trip$i] == patch[trip$j])
  masked <- Matrix::sparseMatrix(
    i = trip$i[same_patch], j = trip$j[same_patch],
    x = rep(1, sum(same_patch)), dims = c(n, n)
  )
  g <- igraph::graph_from_adjacency_matrix(masked, mode = "max", weighted = NULL)
  comp <- igraph::components(g)

  for (p in unique(patch[!is.na(patch)])) {
    pcells <- which(patch == p)
    if (length(pcells) <= 1) next
    comp_ids <- comp$membership[pcells]
    tab <- table(comp_ids)
    largest_comp <- as.integer(names(tab)[which.max(tab)])
    orphans <- pcells[comp_ids != largest_comp]
    if (length(orphans) > 0) {
      patch[orphans] <- NA
    }
  }
  patch
}
