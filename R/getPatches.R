#' @title Get patches
#'
#' @description Assign cells to "patches" for stratified differential expression analysis.
#' Suggested tuning parameters to use: 
#'  1. Control granularity: npatches.
#'  2. Control learning rate: bizesize, then possibly alpha. 
#'  3. Patch size and shape: maxradius, roundness
#'  4. Patch position: initwithhotspots
#'
#' @param spe An
#' [`SpatialExperiment`][SpatialExperiment::SpatialExperiment-class] object.
#'
#' @param cell_type_column A character string specifying the column name in
#' `colData(spe)` that contains the cell type annotations for all cell types.
#'
#' @param response_cell_type A character string specifying the response cell
#' type.
#'
#' @param explanatory_cell_type A character string specifying the explanatory
#' cell type.
#'
#' @param mmxpixel Conversion factor from pixels to millimeters to apply to the
#' spatial coordinates. If NULL (default), no conversion is applied.
#' 
#' @param npatches A positive integer specifying the number of patches to create.
#'
#' @param bitesize How far a patch can expand its borders in any iteration.
#' Controls the "learning rate". Default value is 0.1.
#' @param maxradius Cells beyond this distance from a patch centroid can't be
#' assigned to it. Used to prevent excessively expansive patches. Default value
#' is 0.5
#' @param roundness Tuning parameter from 0-1 guiding patches' tendency to have
#' smooth borders. Lower values get better performance (more consistent var(X)
#' per patch) but produce more uneven borders. Default value is 0.5.
#' @param initwithhotspots TRUE or FALSE for whether to initialize patches
#' around var(X) hotspots (TRUE, default value) or just by clustering cells'
#' positions (FALSE). Anecdotally, TRUE performs better.
#' @param n_iters How many iterations to run. Default value is 25.
#' @param alpha Number from 0-1. How aggressively patches with low var(X) will
#' grab cells from bigger neighboring patches. Default value is 0.5.
#' @param effectivezerodist Lower threshold for cells' distance to patch
#' boundaries. Lower values make it less likely that cells near patch boundaries
#' will be assigned to the smaller patch, which can help prevent patches from
#' getting too big. Default value is 0.001.
#' @param plotprogress Logical for whether to show progress graphically over
#' iterations. Default value is FALSE.
#'
#' @return A character vector of patch assignments for each cell of the type
#' specified by `response_cell_type`. The vector is named with the cell IDs.
#'
#' @examples
#'
#' library(MerfishData)
#'
#' spe <- MouseColonIbdCadinu2024()
#'
#' ## subset to sample 1 slice 1
#' mask <- spe$sample_id == 1 & spe$slice_id == 1
#' spe_sub <- spe[, mask]
#' 
#' patches <- getPatches(spe, cell_type_column="tier1",
#'                       response_cell_type="Epithelial",
#'                       explanatory_cell_type="Immune",
#'                       mmxpixel=0.000109,
#'                       npatches=40, n_iters=10)
#'
#' @importFrom cli cli_abort
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom Matrix rowSums
#'
#' @rdname getPatches
#'
#' @export
getPatches <- function(
    spe,
    cell_type_column,
    response_cell_type,
    explanatory_cell_type,
    mmxpixel=NULL,
    npatches,
    bitesize = 0.1,
    maxradius = 0.5,
    roundness = 0.5,
    initwithhotspots = TRUE,
    n_iters = 25,
    alpha = 0.5,
    effectivezerodist = 0.001,
    plotprogress = FALSE
  ) {
  
    if (!cell_type_column %in% colnames(colData(spe)))
        cli_abort("The specified cell type column '{cell_type_column}' is not present in colData(spe).")
    if (sum(colData(spe)[, cell_type_column] == response_cell_type) == 0)
      cli_abort("No cells of the specified response cell type '{response_cell_type}' were found in colData(spe).")
    if (sum(colData(spe)[, cell_type_column] == explanatory_cell_type) == 0)
      cli_abort("No cells of the specified explanatory cell type '{explanatory_cell_type}' were found in colData(spe).")

    xy <- spatialCoords(spe)
    if (is.null(colnames(spe)))
      rownames(xy) <- paste0("cell", seq_len(nrow(xy)))
    else
      rownames(xy) <- colnames(spe)

    neighbors <- .nearestNeighborGraph(xy, N=50)
    env <- .neighbor_tabulate(as.vector(colData(spe)[, cell_type_column]),
                              neighbors)
    rownames(env) <- rownames(xy)
    if (!is.null(mmpixel))
      xy <- xy * mmxpixel

    mask <- colData(spe)[, cell_type_column] == response_cell_type
    xy <- xy[mask, ]
    X <- env[mask, explanatory_cell_type]

    ## scale X
    X <- X / sd(X, na.rm=TRUE)

    ## get initial patch assignments
    seeds <- .initializePatches(xy, X, neighbors, npatches, initwithhotspots)

    ## set up data frame to track cells
    celldf <- data.frame(X=X)
    rownames(celldf) <- rownames(xy)
    celldf$patch <- rep(NA, nrow(celldf))
    celldf$patch <- seeds

    ## initialize patch data frame
    uniqueseeds <- setdiff(seeds, NA)
    ## for patches with fewer than 2 cells 'totvar' and 'hunger' will be zero
    patchdf <- data.frame(totvar = rep(0, length(uniqueseeds)),
                          hunger = rep(0, length(uniqueseeds)))
    rownames(patchdf) <- uniqueseeds

    ## get patch centroids
    centroids <- c()
    for (name in unique(seeds))
        centroids <- rbind(centroids, colMeans(xy[(seeds == name) & !is.na(seeds), , drop=FALSE]))
    rownames(centroids) <- unique(seeds)

    ## get each cell's distance to each patch:
    cell2patchproximity <- .getproximity2patch(xy, centroids, celldf$patch,
                                               bitesize = bitesize,
                                               effectivezerodist = effectivezerodist,
                                               maxradius = maxradius, roundness = roundness)

    ## reduce patches to those with proximity info 
    patchdf <- patchdf[rownames(patchdf) %in% colnames(cell2patchproximity), , drop=FALSE]

    stopifnot(nrow(patchdf) == ncol(cell2patchproximity)) ## QC

    ## which cells are just too far to assign
    celldf$toofar <- Matrix::rowSums(cell2patchproximity != 0) == 0

    cli_progress_bar("Computing patches", total = n_iters)
    for (iter in seq_len(n_iters)) {

        ## get patch stats
        for (name in uniqueseeds) {
            ## need at least 2 cells to get a variance estimate
            mask <- celldf$patch == name
            if (sum(mask, na.rm = TRUE) > 1)
                patchdf[name, "totvar"] <- var(celldf$X[mask], na.rm = TRUE) *
                                           sum(mask, na.rm = T)
        }

        ## discard NA values from mean calculation and guard against undefined
        ## mean and zero division
        mean_totvar <- mean(patchdf[, "totvar"], na.rm=TRUE)
        if (!is.finite(mean_totvar))
            mean_totvar <- 0

        hunger_denom <- (1 - alpha) * mean_totvar + alpha * patchdf[, "totvar"]
        patchdf[, "hunger"] <- ifelse(is.finite(hunger_denom) & hunger_denom > 0,
                                      1 / hunger_denom, 0)

        ## discard NA values from sum calculation and guard against zero division
        total_hunger <- sum(patchdf$hunger, na.rm=TRUE)
        if (is.finite(total_hunger) && total_hunger > 0)
            patchdf$hunger <- patchdf$hunger / total_hunger

        ## reduce patches to those with proximity info
        mask <- rownames(patchdf) %in% colnames(cell2patchproximity)
        stopifnot(all(mask)) ## QC
        patchdf <- patchdf[mask, , drop=FALSE]

        idx <- apply(cell2patchproximity[, rownames(patchdf)] %*% diag(patchdf$hunger), 1,
                     which.max)

        ## add cells to nearest (hunger-weighted) patch:
        celldf$patch <- rownames(patchdf)[idx]
        celldf$patch[celldf$toofar] = NA

        ## update patch centroids
        centroids <- c()
        patchnames <- setdiff(unique(celldf$patch), NA)
        for (name in patchnames) {
          mask <- (celldf$patch == name) & !is.na(celldf$patch)
          stopifnot(sum(mask) > 0) ## QC
          centroids <- rbind(centroids, colMeans(xy[mask, , drop=FALSE]))
        }
        rownames(centroids) <- patchnames

        ## update each cell's distance to each patch
        cell2patchproximity <- .getproximity2patch(xy, centroids,
                                                   patch = celldf$patch,
                                                   bitesize = bitesize,
                                                   effectivezerodist = effectivezerodist,
                                                   maxradius = maxradius,
                                                   roundness = roundness)
        celldf$toofar <- rowSums(cell2patchproximity != 0) == 0

        stopifnot(nrow(patchdf) == ncol(cell2patchproximity)) ## QC

        cli_progress_update()
    }
    cli_progress_done()

    return(list(celldf=celldf, patchdf=patchdf))
}

#' @importFrom spdep knearneigh knn2nb nb2listw localG
#' @importFrom stats kmeans
.initializePatches <- function(xy, X, neighbors, npatches, initwithhotspots) {
    seeds <- rep(NA, nrow(xy))
    if (initwithhotspots) {
        ## get local X, X^2, var
        cnX <- .neighbor_mean(neighbors=neighbors, x=X)
        cnX2 <- .neighbor_mean(neighbors=neighbors, x=X^2)
        cnvar <- cnX2 - cnX^2

        ## get variance hotspot stats
        knn <- knearneigh(xy, k=8)
        nb <- knn2nb(knn)
        ## convert to a spatial-weights object (row-standardized)
        lw <- nb2listw(nb, style="W", zero.policy=TRUE)
        ## compute the Getis–Ord G* statistic (Gi*) for each cell
        gi <- localG(cnvar, lw, zero.policy=TRUE)

        ## call hotspots, and link them with dbscan
        ishot - gi >= 1
        ## cluster hotspot points to get patch seeds
        hotspotcluster <- kmeans(xy[ishot, ], centers=npatches, iter.max=30)$cluster
        seeds[ishot] <- paste0("patch", hotspotcluster)
    } else {
        ## simple initial clustering: Mclust on xy alone
        simplecluster <- kmeans(xy, centers = npatches, iter.max = 30)$cluster
        seeds <- paste0("patch", simplecluster)
    }

    return(seeds)
}

## Score every cell's proximity to every patch.
## "Proximity" is roughly an inverse distance. Cells sufficiently far from a
## patch get proximity = 0. Distance is calculated as a weighted geometric mean
## of a cell's distance to a patch's center and to its nearest cell within the
## patch (i.e. to its border).
##
## xy Matrix of cells' xy positions
## centroids Matrix of patch centroids (patches * xy position)
## patch Vector of cells' patch assignments. Aligned to the rows to xy.
## bitesize How far a patch can expand its borders in any iteration. Controls the "learning rate"
## maxradius Cells beyond this distance from a patch centroid can't be assigned
## to it. Used to prevent excessively expansive patches.
## roundness Tuning parameter from 0-1 guiding patches' tendency to have smooth
## borders. Lower values get better performance (more consistent var(X) per
## patch) but produce more uneven borders.
## effectivezerodist Lower threshold for cells' distance to patch boundaries.
##
## It returns a matrix of cell * patch proximity scores, usually zero for cells
## that are too far from a given patch.
.getproximity2patch <- function(xy, centroids, patch,
                                bitesize, maxradius, roundness, effectivezerodist) {
    patchnames <- setdiff(unique(patch), NA)
    roundness <- max(min(roundness, 1), 0)
  
    # dist to centroid:  
    centroiddistmat <- matrix(NA, nrow(xy), length(patchnames),
                              dimnames = list(rownames(xy), patchnames))
    for (name in patchnames) {
        ## clamp distances to prevent zeros, which would result below in infinite proximities
        edist <- sqrt((xy[, 1] - centroids[name, 1])^2 + (xy[, 2] - centroids[name, 2])^2)
        centroiddistmat[, name] <- pmax(edist, effectivezerodist)
      }

    # dist to any cells in patch:
    anydistmat <- matrix(NA, nrow(xy), length(patchnames),
                         dimnames = list(rownames(xy), patchnames))
    for (name in patchnames) {
      anydistmat[, name] <- FNN::get.knnx(
        data  = as.matrix(xy[(patch == name) & !is.na(patch), , drop=FALSE]),
        query = as.matrix(xy),
        k = 1
      )$nn.dist[, 1]
    }
    anydistmat[anydistmat < effectivezerodist] <- effectivezerodist

    # weighted average (geomean) of distances to border and centroid:
    bigdist <- exp((1 - roundness) * log2(anydistmat) + roundness *
                   log2(centroiddistmat[rownames(anydistmat), colnames(anydistmat)]))

    # convert to proximity: 
    proxmat <- 1 / bigdist

    # censor cells that are too far:
    proxmat[anydistmat > bitesize] <- 0
    proxmat[centroiddistmat[rownames(anydistmat), colnames(anydistmat)] > maxradius] <- 0

    return(proxmat)
}
