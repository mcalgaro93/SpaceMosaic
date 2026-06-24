
#' Compute patch attributes matrix W from Z and patch assignments
#'
#' Takes colMeans of Z for each patch.
#' @param Z Matrix or data.frame of per-cell features (cells x features).
#'   If a data.frame is supplied, it is converted to a numeric design matrix
#'   using \\code{model.matrix(~ . - 1, data = Z)}.
#' @param patch Character vector of patch assignments (may contain NA).
#' @return Matrix (npatches x features) of per-patch mean attributes.
#' @export
getPatchAttributes <- function(Z, patch) {
  if (is.data.frame(Z)) {
    Z <- stats::model.matrix(~ . - 1, data = Z)
  } else {
    Z <- as.matrix(Z)
    storage.mode(Z) <- "numeric"
  }

  if (length(patch) != nrow(Z)) stop("length(patch) must equal nrow(Z).")

  pnames <- sort(unique(patch[!is.na(patch)]))
  W <- t(sapply(pnames, function(p) colMeans(Z[which(patch == p), , drop = FALSE])))
  rownames(W) <- pnames
  W
}


#' Use local meta-analysis to strengthen per-patch DE results
#'
#' For each patch, finds k nearest neighbors in W-space, computes an
#' inverse-variance-weighted prior from those neighbors' estimates,
#' then performs a Bayesian normal-normal update with the patch's own estimate.
#'
#' @param DEobj Results of patchDE. A list keyed by variable name, each with
#'   $ests, $ses, $pvals matrices (genes x patches).
#' @param W Matrix of patch attributes (patches x features). Rownames must
#'   match the column names of the DE matrices.
#' @param k Number of nearest neighbors to use. Default 15.
#' @param min_effect Minimum absolute posterior estimate for a patch to be
#'   considered significant in subgroup detection. Default 0.
#' @param max_pval Maximum posterior p-value for significance. Default 0.05.
#' @param min_patches Minimum number of patches in a subgroup. Smaller
#'   connected components are set to 0. Default 3.
#' @return A list in the same structure as DEobj, with updated estimates,
#'   SEs, and p-values derived from the local meta-analysis, plus a
#'   \code{$subgroups} matrix (genes x patches) with positive cluster IDs
#'   (1, 2, ...) for significant positive regions, negative IDs (-1, -2, ...)
#'   for significant negative regions, and 0 elsewhere.
#' @export
patchMetaAnalysis <- function(DEobj, W, k = 15,
                              min_effect = 0, max_pval = 0.05,
                              min_patches = 3) {
  ## restrict to patches present in DE results
  de_patches <- colnames(DEobj[[1]]$ests)
  shared_patches <- intersect(rownames(W), de_patches)
  W <- W[shared_patches, , drop = FALSE]

  ## cap k

  k <- min(k, nrow(W) - 1)

  ## find k nearest neighbors for each patch in W-space
  patch_names <- rownames(W)
  nn <- .patchNeighborNetwork(W, k)

  ## apply meta-analysis per variable
  out <- list()
  for (varname in names(DEobj)) {
    ests <- DEobj[[varname]]$ests   # genes x patches
    ses <- DEobj[[varname]]$ses

    ## align patch order
    shared <- intersect(colnames(ests), patch_names)
    ests <- ests[, shared, drop = FALSE]
    ses <- ses[, shared, drop = FALSE]
    nn_aligned <- nn[shared, ]

    post <- .bayesianUpdate(ests, ses, nn_aligned, shared)

    subgroups <- .findSubgroups(post$ests, post$pvals, nn_aligned,
                                min_effect = min_effect,
                                max_pval = max_pval,
                                min_patches = min_patches)

    out[[varname]] <- list(
      ests = post$ests,
      ses = post$ses,
      pvals = post$pvals,
      subgroups = subgroups
    )
  }
  out
}


#' Build k-nearest-neighbor network among patches based on W
#' @param W Patches x features matrix.
#' @param k Number of neighbors.
#' @return Integer matrix (npatches x k) of neighbor indices (column indices into W).
.patchNeighborNetwork <- function(W, k) {
  ## Euclidean distance, find k nearest (excluding self)
  n <- nrow(W)
  k <- min(k, n - 1)
  nn <- FNN::get.knn(W, k = k)$nn.index
  rownames(nn) <- rownames(W)
  nn
}


#' Perform Bayesian normal-normal update using neighbor prior
#' @param ests Genes x patches matrix of estimates.
#' @param ses Genes x patches matrix of standard errors.
#' @param nn Integer matrix (patches x k) of neighbor indices.
#' @param patch_names Character vector of patch names (column order).
#' @return List with $ests, $ses, $pvals matrices (genes x patches).
.bayesianUpdate <- function(ests, ses, nn, patch_names) {
  ngenes <- nrow(ests)
  np <- ncol(ests)
  k <- ncol(nn)

  post_ests <- matrix(NA_real_, nrow = ngenes, ncol = np)
  post_ses <- matrix(NA_real_, nrow = ngenes, ncol = np)
  rownames(post_ests) <- rownames(post_ses) <- rownames(ests)
  colnames(post_ests) <- colnames(post_ses) <- patch_names

  for (j in seq_len(np)) {
    ## neighbor indices for this patch
    nbr_idx <- nn[j, ]

    ## inverse-variance-weighted prior from neighbors
    nbr_ests <- ests[, nbr_idx, drop = FALSE]  # genes x k
    nbr_ses <- ses[, nbr_idx, drop = FALSE]
    nbr_prec <- 1 / (nbr_ses^2)
    ## handle infinite/NA precision
    nbr_prec[!is.finite(nbr_prec)] <- 0

    prior_prec <- rowSums(nbr_prec)
    prior_mean <- ifelse(prior_prec > 0,
                         rowSums(nbr_ests * nbr_prec) / prior_prec,
                         0)
    prior_se <- ifelse(prior_prec > 0, 1 / sqrt(prior_prec), Inf)

    ## patch's own data
    patch_est <- ests[, j]
    patch_se <- ses[, j]
    patch_prec <- 1 / (patch_se^2)
    patch_prec[!is.finite(patch_prec)] <- 0

    ## Bayesian update: normal-normal conjugate
    post_prec <- prior_prec + patch_prec
    post_ests[, j] <- ifelse(post_prec > 0,
                             (prior_mean * prior_prec + patch_est * patch_prec) / post_prec,
                             patch_est)
    post_ses[, j] <- ifelse(post_prec > 0, 1 / sqrt(post_prec), patch_se)
  }

  ## two-sided p-values from posterior
  z <- post_ests / post_ses
  post_pvals <- 2 * stats::pnorm(-abs(z))
  post_pvals[!is.finite(z)] <- 1

  list(ests = post_ests, ses = post_ses, pvals = post_pvals)
}


#' Find subgroups of patches with coherent significant effects
#'
#' Thresholds posterior estimates by effect size and p-value, then finds
#' connected components on the W-neighbor graph among significant patches.
#'
#' @param ests Genes x patches matrix of posterior estimates.
#' @param pvals Genes x patches matrix of posterior p-values.
#' @param nn Integer matrix (patches x k) of W-neighbor indices.
#' @param min_effect Minimum absolute estimate for significance.
#' @param max_pval Maximum p-value for significance.
#' @param min_patches Minimum cluster size; smaller clusters set to 0.
#' @return Matrix (genes x patches) of subgroup labels. Positive clusters
#'   labeled 1, 2, ...; negative clusters labeled -1, -2, ...; non-significant = 0.
.findSubgroups <- function(ests, pvals, nn, min_effect, max_pval, min_patches) {
  ngenes <- nrow(ests)
  np <- ncol(ests)
  pnames <- colnames(ests)

  ## build W-neighbor adjacency graph (undirected) from nn matrix
  from <- rep(seq_len(np), each = ncol(nn))
  to <- as.vector(t(nn))
  g <- igraph::make_empty_graph(n = np, directed = FALSE)
  g <- igraph::add_edges(g, as.vector(rbind(from, to)))
  g <- igraph::simplify(g)

  result <- matrix(0L, nrow = ngenes, ncol = np)
  rownames(result) <- rownames(ests)
  colnames(result) <- pnames

  for (i in seq_len(ngenes)) {
    est_i <- ests[i, ]
    pval_i <- pvals[i, ]

    ## positive significant patches
    pos_idx <- which(est_i >= min_effect & pval_i <= max_pval)
    pos_label <- 0L
    if (length(pos_idx) > 0) {
      sub_g <- igraph::induced_subgraph(g, pos_idx)
      comp <- igraph::components(sub_g)
      for (cl in seq_len(comp$no)) {
        members <- pos_idx[comp$membership == cl]
        if (length(members) >= min_patches) {
          pos_label <- pos_label + 1L
          result[i, members] <- pos_label
        }
      }
    }

    ## negative significant patches
    neg_idx <- which(est_i <= -min_effect & pval_i <= max_pval)
    neg_label <- 0L
    if (length(neg_idx) > 0) {
      sub_g <- igraph::induced_subgraph(g, neg_idx)
      comp <- igraph::components(sub_g)
      for (cl in seq_len(comp$no)) {
        members <- neg_idx[comp$membership == cl]
        if (length(members) >= min_patches) {
          neg_label <- neg_label - 1L
          result[i, members] <- neg_label
        }
      }
    }
  }
  result
}


#' Summarize subgroups with enrichment of cell metadata
#'
#' For each gene x subgroup combination, computes the proportion of patches
#' in the subgroup, the mean posterior effect size, and enrichment scores
#' for each metadata variable.
#'
#' Enrichment for numeric variables: (mean_in_subgroup - mean_overall) / sd_overall.
#' Enrichment for categorical levels: log2(proportion_in_subgroup / proportion_overall).
#'
#' @param meta Output of \code{patchMetaAnalysis}. A list keyed by variable name,
#'   each with \code{$ests}, \code{$subgroups} matrices (genes x patches).
#' @param patch Character vector of patch assignments (may contain NA), aligned to
#'   rows of \code{cellmeta}.
#' @param cellmeta Data frame of cell-level metadata, rows aligned to \code{patch}.
#' @return Data frame with one row per gene x subgroup. Columns include
#'   \code{variable}, \code{gene}, \code{subgroup_id}, \code{n_patches},
#'   \code{prop_patches}, \code{mean_effect}, and one enrichment column per
#'   metadata variable/level.
#' @export
summarizeSubgroups <- function(meta, patch, cellmeta) {
  cellmeta <- as.data.frame(cellmeta)
  all_patches <- unique(patch[!is.na(patch)])
  n_total_patches <- length(all_patches)

  ## precompute per-patch metadata summaries
  patch_meta <- .patchMetaSummary(patch, cellmeta, all_patches)

  rows <- list()
  for (varname in names(meta)) {
    sg_mat <- meta[[varname]]$subgroups
    ests_mat <- meta[[varname]]$ests
    pnames <- colnames(sg_mat)
    gene_names <- rownames(sg_mat)

    for (i in seq_len(nrow(sg_mat))) {
      labels <- sg_mat[i, ]
      unique_labels <- setdiff(unique(labels), 0)
      if (length(unique_labels) == 0) next

      for (lab in unique_labels) {
        sg_patches <- pnames[labels == lab]
        n_sg <- length(sg_patches)

        row <- data.frame(
          variable = varname,
          gene = gene_names[i],
          subgroup_id = lab,
          n_patches = n_sg,
          prop_patches = n_sg / n_total_patches,
          mean_effect = mean(ests_mat[i, sg_patches]),
          stringsAsFactors = FALSE
        )

        enrich <- .computeEnrichment(sg_patches, patch_meta)
        row <- cbind(row, enrich)
        rows <- c(rows, list(row))
      }
    }
  }

  if (length(rows) == 0) {
    return(data.frame(variable = character(), gene = character(),
                      subgroup_id = integer(), n_patches = integer(),
                      prop_patches = numeric(), mean_effect = numeric(),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}


#' Compute per-patch summaries of cell metadata
#' @param patch Patch assignment vector.
#' @param cellmeta Cell metadata data frame.
#' @param all_patches Character vector of all patch names.
#' @return List of precomputed summary statistics.
.patchMetaSummary <- function(patch, cellmeta, all_patches) {
  numeric_vars <- names(cellmeta)[sapply(cellmeta, is.numeric)]
  cat_vars <- names(cellmeta)[sapply(cellmeta, function(x) is.character(x) || is.factor(x))]
  assigned <- !is.na(patch)

  ## numeric: per-patch means, overall mean/sd
  numeric_means <- NULL
  numeric_overall_mean <- NULL
  numeric_overall_sd <- NULL
  if (length(numeric_vars) > 0) {
    numeric_means <- matrix(NA_real_, nrow = length(all_patches),
                            ncol = length(numeric_vars))
    rownames(numeric_means) <- all_patches
    colnames(numeric_means) <- numeric_vars
    for (p in all_patches) {
      idx <- which(patch == p)
      if (length(idx) > 0)
        numeric_means[p, ] <- colMeans(cellmeta[idx, numeric_vars, drop = FALSE],
                                       na.rm = TRUE)
    }
    numeric_overall_mean <- colMeans(cellmeta[assigned, numeric_vars, drop = FALSE],
                                     na.rm = TRUE)
    numeric_overall_sd <- sapply(cellmeta[assigned, numeric_vars, drop = FALSE],
                                 sd, na.rm = TRUE)
    numeric_overall_sd[numeric_overall_sd == 0] <- 1
  }

  ## categorical: per-patch proportions, overall proportions
  cat_props <- list()
  cat_overall_props <- list()
  for (v in cat_vars) {
    vals <- as.character(cellmeta[[v]])
    levs <- sort(unique(vals[!is.na(vals)]))
    overall_tab <- table(factor(vals[assigned], levels = levs))
    cat_overall_props[[v]] <- as.numeric(overall_tab) / sum(overall_tab)
    names(cat_overall_props[[v]]) <- levs
    prop_mat <- matrix(0, nrow = length(all_patches), ncol = length(levs))
    rownames(prop_mat) <- all_patches
    colnames(prop_mat) <- levs
    for (p in all_patches) {
      idx <- which(patch == p)
      if (length(idx) > 0) {
        tab <- table(factor(vals[idx], levels = levs))
        prop_mat[p, ] <- as.numeric(tab) / sum(tab)
      }
    }
    cat_props[[v]] <- prop_mat
  }

  list(numeric_means = numeric_means,
       numeric_overall_mean = numeric_overall_mean,
       numeric_overall_sd = numeric_overall_sd,
       numeric_vars = numeric_vars,
       cat_props = cat_props,
       cat_overall_props = cat_overall_props,
       cat_vars = cat_vars)
}


#' Compute enrichment scores for a set of patches relative to overall
#' @param sg_patches Character vector of patch names in the subgroup.
#' @param patch_meta Output of .patchMetaSummary.
#' @return Single-row data frame of enrichment columns.
.computeEnrichment <- function(sg_patches, patch_meta) {
  result <- list()

  ## numeric: (mean_in_subgroup - overall_mean) / overall_sd
  if (length(patch_meta$numeric_vars) > 0) {
    sg_means <- colMeans(patch_meta$numeric_means[sg_patches, , drop = FALSE],
                         na.rm = TRUE)
    for (v in patch_meta$numeric_vars) {
      result[[v]] <- (sg_means[v] - patch_meta$numeric_overall_mean[v]) /
                     patch_meta$numeric_overall_sd[v]
    }
  }

  ## categorical: log2(prop_in_subgroup / prop_overall) per level
  for (v in patch_meta$cat_vars) {
    sg_prop <- colMeans(patch_meta$cat_props[[v]][sg_patches, , drop = FALSE],
                        na.rm = TRUE)
    overall_prop <- patch_meta$cat_overall_props[[v]]
    for (lev in names(overall_prop)) {
      col_name <- paste0(v, "_", lev)
      p_sg <- sg_prop[lev]
      p_all <- overall_prop[lev]
      if (p_all > 0 && p_sg > 0) {
        result[[col_name]] <- log2(p_sg / p_all)
      } else if (p_sg == 0) {
        result[[col_name]] <- -Inf
      } else {
        result[[col_name]] <- Inf
      }
    }
  }

  as.data.frame(result, stringsAsFactors = FALSE)
}
