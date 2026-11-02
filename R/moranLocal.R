#' Calculate local Moran diagnostics for regression residuals
#'
#' Build the same symmetric binary within-patch k-nearest-neighbor graphs as
#' [getMoranGlobal()] and return cell-level quantities for Moran scatterplots
#' and spatial maps. A vector describes one gene; a cells by genes matrix
#' describes multiple genes. This function is descriptive: it does not perform
#' local significance tests or classify cells as significant spatial clusters.
#'
#' @inheritParams getMoranGlobal
#' @param verbose Logical; report the workload and show progress over patches.
#'
#' @details
#' Within each gene and patch, write the residual vector as `r` and its
#' centered values as `z = r - mean(r)`, exactly as in [getMoranGlobal()].
#' Let `W` be the symmetric binary spatial weights matrix, `S0 = sum(W)` and
#' `m2 = sum(z^2) / n`. The local statistic is
#' \deqn{\mathbf{I}_{\mathrm{local}} =
#'   \frac{n}{\mathbf{z}^{\top}\mathbf{z}}
#'   [\mathbf{z} \odot (W\mathbf{z})],}
#' where \eqn{\odot} denotes element-wise multiplication. Equivalently,
#' \eqn{I_i^{\mathrm{local}} = z_i (W\mathbf{z})_i / m_2}.
#' The returned quantities are:
#' \describe{
#'   \item{`residual_centered`}{`z`, the residual minus its within-patch mean.}
#'   \item{`lag_centered`}{`W %*% z`, the sum of centered residuals of linked
#'     cells. Note that binary weights are not normalized by row, so this is not a
#'     neighbor average.}
#'   \item{`local_i`}{`z * lag_centered / m2`, local Moran's I 
#'     (this corresponds to the `mlvar = TRUE` convention in `spdep`).}
#'   \item{`lag_scaled`}{`(n / S0) * lag_centered`. A scatterplot with
#'     `residual_centered` on the x-axis and `lag_scaled` on the y-axis has a
#'     regression slope through the origin equal to the observed global
#'     Moran's I.}
#'   \item{`quadrant`}{`"High-High"`, `"Low-Low"`, `"High-Low"`, or
#'     `"Low-High"`, using the signs of `residual_centered` and `lag_centered`.
#'     Cells exactly on either zero axis have a missing quadrant. These labels
#'     describe signs, not statistical significance.}
#' }
#'
#' The local statistics satisfy `sum(local_i) / S0 = I_global`. Global Moran's
#' I can also be recovered as `sum(residual_centered * lag_scaled) /
#' sum(residual_centered^2)`. 
#'
#' A graph is built once per patch and reused across genes. Subset genes and
#' cells before calling this function to limit the size of its cell-level
#' output. Subsetting cells changes the graph: retain every cell of a selected
#' patch to reproduce the corresponding [getMoranGlobal()] result.
#'
#' Complete graphs are flagged by `complete_graph = TRUE` even when they arise
#' through symmetrization. Their calculations have status `"ok"`, but global
#' Moran's I cannot distinguish spatial arrangements on a complete graph.
#'
#' @return A data frame with one row per input cell and gene, ordered by gene
#'   and then original cell position. Columns include `cell_index` (input row
#'   number), `cell_id` (residual row names, coordinate row names, or input row
#'   numbers as character strings), `patch`, `gene`, `x`, `y`, `n_cells`,
#'   `residual`, `residual_centered`, `lag_centered`, `lag_scaled`, `local_i`,
#'   `quadrant`, `k`, `effective_k`, `complete_graph`,
#'   and `status`.
#'
#'   Unassigned cells are retained with status `"unassigned cell"`. Patch-gene
#'   combinations with fewer than three cells, non-finite residuals or constant
#'   residuals use the same statuses as [getMoranGlobal()]. Their calculated
#'   quantities and quadrants are missing; original residuals and coordinates
#'   are preserved. All cells of a non-finite patch-gene combination are flagged
#'   rather than dropping individual cells and changing its graph.
#'
#' @seealso [getMoranGlobal()], [plotMoranResiduals()], [plotMoranLocal()], [plotMoranQuadrants()],
#'   [plotMoranScatter()]
#' @examples
#' xy <- cbind(x = c(0, 1, 3, 7), y = rep(0, 4))
#' residuals <- c(-3, -1, 1, 3)
#' local <- getMoranLocal(residuals, xy, k = 1, verbose = FALSE)
#' global <- getMoranGlobal(residuals, xy, k = 1, verbose = FALSE)
#' with(local, sum(residual_centered * lag_scaled) / sum(residual_centered^2))
#' global$observed
#' plot(local$residual_centered, local$lag_scaled,
#'      xlab = "Centered residual", ylab = "Globally scaled spatial lag")
#' abline(a = 0, b = global$observed)
#' @export
getMoranLocal <- function(residuals, xy, patch = NULL, k = 10L,
                          verbose = TRUE) {
  # Preserve vector names before the shared preparation converts to a matrix.
  cell_ids <- if (is.null(dim(residuals))) names(residuals) else rownames(residuals)
  inputs <- .prepareMoranInputs(residuals, xy, patch, k, verbose)
  residuals <- inputs$residuals
  patch <- as.character(inputs$patch)
  k <- inputs$k
  n <- nrow(residuals)
  if (is.null(cell_ids)) cell_ids <- rownames(xy)
  if (is.null(cell_ids)) cell_ids <- as.character(seq_len(n))

  patch_names <- inputs$patch_names
  if (anyNA(patch)) patch_names <- c(patch_names, NA_character_)
  results <- vector("list", length(patch_names) * ncol(residuals))
  result_index <- 0L

  if (verbose) {
    cli::cli_alert_info(sprintf(
      "getMoranLocal: %d cells x %d genes = %s cell-gene rows (no local tests).",
      n, ncol(residuals),
      format(as.double(n) * ncol(residuals), big.mark = ",", trim = TRUE)
    ))
    cli::cli_progress_bar("Local Moran diagnostics", total = length(patch_names))
  }

  for (patch_name in patch_names) {
    unassigned <- is.na(patch_name)
    cell_index <- if (unassigned) which(is.na(patch)) else which(patch == patch_name)
    n_cells <- length(cell_index)
    W <- NULL
    effective_k <- NA_integer_
    complete_graph <- NA
    if (!unassigned && n_cells >= 3L) {
      W <- .buildContiguityGraph(xy[cell_index, , drop = FALSE], k = k)
      effective_k <- min(k, n_cells - 1L)
      s0 <- sum(W)
      complete_graph <- s0 == as.double(n_cells) * (n_cells - 1L)
    }

    for (gene_index in seq_len(ncol(residuals))) {
      values <- residuals[cell_index, gene_index]
      status <- "ok"
      if (unassigned) {
        status <- "unassigned cell"
      } else if (n_cells < 3L) {
        status <- "fewer than three cells"
      } else if (any(!is.finite(values))) {
        status <- "non-finite residuals"
      } else if (all(values == values[1L])) {
        status <- "constant residuals"
      }

      z <- lag_centered <- lag_scaled <- local_i <- rep(NA_real_, n_cells)
      quadrant <- rep(NA_character_, n_cells)
      if (status == "ok") {
        prepared <- .prepareMoran(values, W, s0 = s0)
        z <- prepared$z
        m2 <- sum(z^2) / n_cells
        lag_centered <- as.numeric(W %*% z)
        lag_scaled <- (n_cells / s0) * lag_centered
        local_i <- z * lag_centered / m2
        quadrant[z > 0 & lag_centered > 0] <- "High-High"
        quadrant[z < 0 & lag_centered < 0] <- "Low-Low"
        quadrant[z > 0 & lag_centered < 0] <- "High-Low"
        quadrant[z < 0 & lag_centered > 0] <- "Low-High"
      }

      result_index <- result_index + 1L
      results[[result_index]] <- data.frame(
        cell_index = cell_index, cell_id = cell_ids[cell_index],
        patch = patch_name, gene = colnames(residuals)[gene_index],
        x = xy[cell_index, 1], y = xy[cell_index, 2],
        n_cells = if (unassigned) NA_integer_ else n_cells,
        residual = values, residual_centered = z, lag_centered = lag_centered,
        lag_scaled = lag_scaled, local_i = local_i,
        quadrant = quadrant,
        k = k, effective_k = effective_k, complete_graph = complete_graph,
        status = status, stringsAsFactors = FALSE
      )
    }
    if (verbose) cli::cli_progress_update()
  }

  if (verbose) cli::cli_progress_done()
  out <- do.call(rbind, results)
  out <- out[order(match(out$gene, colnames(residuals)), out$cell_index), ]
  rownames(out) <- NULL
  out
}
