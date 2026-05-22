#' Assign cells to patches using elliptical Gaussian assignments
#'
#' Iterative EM-like algorithm (per iteration):
#'   (a) Estimate per-patch centroid + covariance from xy.
#'   (b) Assign cells using spatial fit + X-diversity boost + Z-penalty + hunger.
#'   (c) Re-estimate per-patch centroid + covariance.
#'   (d) Assign cells using spatial fit only (no X/Z/hunger).
#'   (e) Contiguity check: set orphans to NA.
#'
#' Uses candidate filtering (only evaluates nearby patches per cell) and
#' mutual k-NN contiguity for scalability to 100k+ cells.
#'
#' @param xy Matrix of cells' xy positions. Must have rownames.
#' @param X Numeric vector or matrix of design variables, aligned to rows of xy.
#'   If a vector, treated as a single column. Each column is scaled to unit SD.
#' @param npatches Number of patches to create.
#' @param Z Optional matrix of per-cell context embeddings (cells x features).
#'   If supplied, patches will prefer Z-coherent regions. NULL disables.
#' @param alpha Weight of the Z penalty. Typical range 0.2--1.0; default 0.5.
#'   Higher values force patches to respect microenvironment boundaries at the
#'   cost of spatial compactness. Set to 0 to ignore Z entirely.
#' @param beta Weight of the X diversity boost. Typical range 0.5/K--3/K where
#'   K = ncol(X); default 1 (appropriate for single-column X). For multi-column
#'   X, scale down proportionally (e.g. beta = 0.2 for K = 5). Raise if patches
#'   are too homogeneous; lower if patches are spatially fragmented.
#' @param hunger_weight Controls how aggressively low-variance patches grab cells.
#'   Typical range 0.3--0.7; default 0.5. 0 = all patches equally hungry
#'   (uniform sizes), 1 = hunger proportional to 1/totvar (maximizes variance
#'   reduction but allows extreme size imbalance). Lower if some patches shrink
#'   to nothing; raise if variance reduction is insufficient.
#' @param max_elongation Maximum ratio of largest to smallest eigenvalue of a
#'   patch covariance matrix. Typical range 2--8; default 4. Lower values force
#'   rounder patches; raise if tissue structures are genuinely elongated.
#' @param max_radius Maximum Euclidean distance from a patch centroid for
#'   assignment. Cells beyond this get zero spatial score. Default NULL (auto:
#'   3x the expected patch radius assuming uniform circular patches). Override
#'   if patches span very different density regions.
#' @param mahal_radius Maximum Mahalanobis radius for assignment. Cells beyond
#'   this (in each patch's own coordinate system) get zero spatial score.
#'   Typical range 2--4; default 3. Lower values make tighter patches with more
#'   unassigned cells; Inf disables the cutoff entirely.
#' @param n_candidates Number of nearest patch centroids to evaluate per cell.
#'   Typical range 10--50; default 20. Higher values are more accurate but
#'   slower. Rarely needs tuning unless npatches is very large (>1000).
#' @param n_iters Number of outer iterations. Typical range 10--30; default 15.
#'   Convergence is usually reached by 10--15; raise if patches are still
#'   shifting at the final iteration (check membership_log).
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
                       max_radius = NULL,
                       mahal_radius = 3,
                       n_candidates = 20,
                       n_iters = 15,
                       log_iters = TRUE,
                       verbose = TRUE) {

  ## coerce X to matrix
  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  X <- as.matrix(X)
  stopifnot(nrow(xy) == nrow(X))
  stopifnot(!is.null(rownames(xy)))
  n <- nrow(xy)
  n_candidates <- min(n_candidates, npatches)

  ## scale X column-wise
  X <- scale(X, center = FALSE, scale = apply(X, 2, sd, na.rm = TRUE))

  ## auto-compute Euclidean radius cap if not supplied
  if (is.null(max_radius)) {
    hull <- grDevices::chull(xy)
    hx <- xy[hull, 1]; hy <- xy[hull, 2]
    hull_area <- 0.5 * abs(sum(hx * c(hy[-1], hy[1]) - c(hx[-1], hx[1]) * hy))
    max_radius <- 3 * sqrt(hull_area / npatches / pi)
  }

  ## build contiguity graph: mutual k-NN
  if (verbose) cli::cli_alert_info("Building contiguity graph...")
  spatial_nn <- .buildContiguityGraph(xy)

  ## prepare Z
  if (verbose) cli::cli_alert_info("Preparing Z...")
  z_info <- .prepareZ(Z, xy)

  ## initialize with kmeans on xy
  if (verbose) cli::cli_alert_info("Initializing patches...")
  patch <- .initKmeans(xy, npatches)

  ## set up iteration logs
  if (log_iters) {
    membership_log <- matrix(NA_character_, nrow = n, ncol = n_iters)
    ss_log <- matrix(NA_real_, nrow = npatches, ncol = n_iters)
    rownames(ss_log) <- as.character(seq_len(npatches))
  }

  if (verbose) cli::cli_progress_bar("Patch iterations", total = n_iters)

  for (iter in seq_len(n_iters)) {

    ## (a) estimate ellipse parameters
    params <- .estimateEllipses(xy, patch, max_elongation)

    ## (b) assign using spatial + X + Z + hunger
    patch <- .assignCells(xy, X, patch, params,
                          z_info = z_info, alpha = alpha, beta = beta,
                          hunger_weight = hunger_weight,
                          max_radius = max_radius,
                          mahal_radius = mahal_radius,
                          n_candidates = n_candidates,
                          use_xz = TRUE)

    ## (c) re-estimate ellipse parameters
    params <- .estimateEllipses(xy, patch, max_elongation)

    ## (d) assign using spatial only
    patch <- .assignCells(xy, X, patch, params,
                          z_info = z_info, alpha = alpha, beta = beta,
                          hunger_weight = hunger_weight,
                          max_radius = max_radius,
                          mahal_radius = mahal_radius,
                          n_candidates = n_candidates,
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


#' Prepare Z: spatially smooth and compute z_scale
#' @param Z Raw Z matrix (cells x features), or NULL.
#' @param xy n x 2 coordinate matrix.
#' @return List with Z_smooth and z_scale, or list(NULL, NULL).
.prepareZ <- function(Z, xy) {
  if (is.null(Z)) return(list(Z_smooth = NULL, z_scale = NULL))
  Z <- as.matrix(Z)
  stopifnot(nrow(Z) == nrow(xy))
  n <- nrow(xy)
  ## use FNN for neighbor smoothing
  knn <- FNN::get.knn(xy, k = 10)
  ## smooth Z: for each cell, average Z over its 10 nearest neighbors
  ## vectorized via index matrix
  nn_idx <- knn$nn.index  # n x 10
  Z_smooth <- matrix(0, nrow = n, ncol = ncol(Z))
  for (k in seq_len(10)) {
    Z_smooth <- Z_smooth + Z[nn_idx[, k], , drop = FALSE]
  }
  Z_smooth <- Z_smooth / 10
  colnames(Z_smooth) <- colnames(Z)

  ## z_scale: mean Z-distance between spatial neighbors (sampled for speed)
  sample_n <- min(n, 10000)
  sample_idx <- sample.int(n, sample_n)
  dists <- numeric(sample_n * 10)
  for (s in seq_len(sample_n)) {
    idx <- sample_idx[s]
    nbrs <- nn_idx[idx, ]
    diffs <- Z_smooth[nbrs, , drop = FALSE] -
             matrix(Z_smooth[idx, ], nrow = 10, ncol = ncol(Z_smooth), byrow = TRUE)
    dists[((s - 1) * 10 + 1):(s * 10)] <- sqrt(rowSums(diffs^2))
  }
  z_scale <- mean(dists)
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
#'   \item{centroids}{Matrix (npatches x 2) of patch centroids.}
#'   \item{inv_covmats}{Named list of 2x2 inverse covariance matrices.}
#'   \item{log_det}{Named numeric vector of log-determinants.}
#'   \item{pnames}{Sorted patch names.}
.estimateEllipses <- function(xy, patch, max_elongation) {
  pnames <- sort(unique(patch[!is.na(patch)]))
  np <- length(pnames)
  ## pre-split indices
  cell_lists <- split(seq_len(length(patch)), patch)
  cell_lists <- cell_lists[pnames]

  centroids_mat <- matrix(NA_real_, nrow = np, ncol = 2)
  rownames(centroids_mat) <- pnames
  inv_covmats <- vector("list", np)
  names(inv_covmats) <- pnames
  log_det <- numeric(np)
  names(log_det) <- pnames
  global_var <- mean(apply(xy, 2, var))

  for (j in seq_len(np)) {
    cells <- cell_lists[[j]]
    centroids_mat[j, ] <- colMeans(xy[cells, , drop = FALSE])
    if (length(cells) < 3) {
      S <- diag(2) * global_var
    } else {
      S <- stats::cov(xy[cells, , drop = FALSE])
    }
    S <- .regularizeCov(S, max_elongation)
    inv_covmats[[pnames[j]]] <- solve(S)
    log_det[pnames[j]] <- log(det(S))
  }
  list(centroids = centroids_mat, inv_covmats = inv_covmats,
       log_det = log_det, pnames = pnames)
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


#' Assign cells to patches using candidate-filtered scoring
#'
#' For each cell, only evaluates the n_candidates nearest patch centroids.
#' Score = spatial_score + beta * X_boost - alpha * Z_penalty + log(hunger)
#'
#' @param xy n x 2 coordinate matrix.
#' @param X Numeric matrix of design variables (n x p).
#' @param patch Current patch assignments.
#' @param params Output of .estimateEllipses.
#' @param z_info Output of .prepareZ.
#' @param alpha Z penalty weight.
#' @param beta X diversity boost weight.
#' @param hunger_weight Hunger interpolation weight.
#' @param max_radius Max Euclidean distance for assignment.
#' @param mahal_radius Max Mahalanobis distance for assignment (Inf to disable).
#' @param n_candidates Number of candidate patches per cell.
#' @param use_xz Whether to include X, Z, and hunger terms.
#' @return Updated patch vector.
.assignCells <- function(xy, X, patch, params,
                         z_info, alpha, beta, hunger_weight,
                         max_radius, mahal_radius, n_candidates, use_xz) {
  n <- nrow(xy)
  pnames <- params$pnames
  np <- length(pnames)
  n_candidates <- min(n_candidates, np)

  ## find candidate patches for each cell (nearest centroids)
  candidates <- FNN::get.knnx(params$centroids, xy, k = n_candidates)
  cand_idx <- candidates$nn.index  # n x n_candidates

  ## precompute patch-level statistics (only when use_xz = TRUE)
  if (use_xz) {
    cell_lists <- split(seq_len(n), patch)

    ## patch mean X (np x ncol(X))
    patch_mean_X <- matrix(0, nrow = np, ncol = ncol(X))
    rownames(patch_mean_X) <- pnames
    for (j in seq_len(np)) {
      cells <- cell_lists[[pnames[j]]]
      if (!is.null(cells) && length(cells) > 0) {
        patch_mean_X[j, ] <- colMeans(X[cells, , drop = FALSE], na.rm = TRUE)
      }
    }

    ## patch mean Z (if applicable)
    patch_mean_Z <- NULL
    if (!is.null(z_info$Z_smooth)) {
      nz <- ncol(z_info$Z_smooth)
      patch_mean_Z <- matrix(0, nrow = np, ncol = nz)
      rownames(patch_mean_Z) <- pnames
      for (j in seq_len(np)) {
        cells <- cell_lists[[pnames[j]]]
        if (!is.null(cells) && length(cells) > 0) {
          patch_mean_Z[j, ] <- colMeans(z_info$Z_smooth[cells, , drop = FALSE])
        }
      }
    }

    ## hunger: variance-based (trace of within-patch covariance * n)
    patch_totvar <- numeric(np)
    for (j in seq_len(np)) {
      cells <- cell_lists[[pnames[j]]]
      if (!is.null(cells) && length(cells) > 1) {
        patch_totvar[j] <- sum(apply(X[cells, , drop = FALSE], 2, var, na.rm = TRUE)) *
                           length(cells)
      }
    }
    mean_totvar <- mean(patch_totvar, na.rm = TRUE)
    if (!is.finite(mean_totvar) || mean_totvar == 0) mean_totvar <- 1
    hunger_denom <- (1 - hunger_weight) * mean_totvar + hunger_weight * patch_totvar
    hunger <- ifelse(hunger_denom > 0, 1 / hunger_denom, 0)
    total_hunger <- sum(hunger)
    if (total_hunger > 0) hunger <- hunger / total_hunger
    log_hunger <- ifelse(hunger > 0, log(hunger), -Inf)
  }

  ## score each cell against its candidates
  best_patch_idx <- integer(n)
  best_score <- rep(-Inf, n)

  for (c_rank in seq_len(n_candidates)) {
    ## which patch index does each cell consider at this rank?
    pidx <- cand_idx[, c_rank]  # length n, values 1..np

    ## spatial score: Mahalanobis distance
    ## group cells by their candidate patch for batch computation
    centroid_xy <- params$centroids[pidx, , drop = FALSE]  # n x 2
    d <- xy - centroid_xy
    spatial_score <- rep(-Inf, n)
    patch_groups <- split(seq_len(n), pidx)
    for (pg in names(patch_groups)) {
      idx <- patch_groups[[pg]]
      j <- as.integer(pg)
      p <- pnames[j]
      d_sub <- d[idx, , drop = FALSE]
      mah2 <- rowSums((d_sub %*% params$inv_covmats[[p]]) * d_sub)
      euc2 <- rowSums(d_sub^2)
      s <- -0.5 * mah2 - 0.5 * params$log_det[p]
      s[euc2 > max_radius^2 | mah2 > mahal_radius^2] <- -Inf
      spatial_score[idx] <- s
    }

    ## X diversity boost, Z penalty, hunger (only when use_xz = TRUE)
    if (use_xz) {
      dX <- X - patch_mean_X[pidx, , drop = FALSE]
      x_boost <- beta * rowSums(dX^2)

      z_penalty <- 0
      if (!is.null(patch_mean_Z)) {
        dZ <- z_info$Z_smooth - patch_mean_Z[pidx, , drop = FALSE]
        z_penalty <- alpha * rowSums(dZ^2) / (z_info$z_scale^2)
      }

      score <- spatial_score + x_boost - z_penalty + log_hunger[pidx]
    } else {
      score <- spatial_score
    }

    ## update best
    improved <- score > best_score
    best_score[improved] <- score[improved]
    best_patch_idx[improved] <- pidx[improved]
  }

  ## assign
  best_patch_idx[best_patch_idx == 0L] <- NA_integer_
  new_patch <- pnames[best_patch_idx]
  ## cells with all -Inf scores stay NA
  new_patch[best_score == -Inf] <- NA
  new_patch
}


#' Build contiguity graph using symmetric k-NN
#'
#' Two cells are connected if either is among the other's k nearest neighbors.
#' Much faster than Delaunay triangulation for large n.
#'
#' @param xy n x 2 coordinate matrix.
#' @param k Number of neighbors for k-NN. Default 10.
#' @return Sparse binary adjacency matrix (n x n).
.buildContiguityGraph <- function(xy, k = 10) {
  n <- nrow(xy)
  knn <- FNN::get.knn(xy, k = k)
  i_vec <- rep(seq_len(n), each = k)
  j_vec <- as.vector(t(knn$nn.index))
  directed <- Matrix::sparseMatrix(
    i = i_vec, j = j_vec,
    x = rep(1, length(i_vec)),
    dims = c(n, n)
  )
  ## symmetric: connect if either direction exists
  sym <- directed + Matrix::t(directed)
  sym@x[] <- 1
  sym
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
  ## build graph from edge list (avoids sparse->dense coercion)
  keep <- same_patch & (trip$i < trip$j)
  g <- igraph::make_empty_graph(n = n, directed = FALSE)
  if (any(keep)) {
    g <- igraph::add_edges(g, as.vector(t(cbind(trip$i[keep], trip$j[keep]))))
  }
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
