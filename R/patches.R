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
#' @param verbose Show progress on the computations. Default value is TRUE.
#'
#' @return The input
#' [`SpatialExperiment`][SpatialExperiment::SpatialExperiment-class] object with
#' a new reduced dimension data frame containing the patch assignments for the
#' cells annotated to the given cell type response. The reduced dimension data
#' frame is stored as a reduced dimension under the name
#' `patches_<response_cell_type>_<explanatory_cell_type>`. This data frame is
#' in fact a [`DataFrame`][S4Vectors::DataFrame-class] object
#' with a metadata slot where the patch-level information is stored.
#'
#' @examples
#'
#' ## load a subset of 1k cells of the first sample and slice of the MERFISH
#' ## data available in MerfishData::MouseColonIbdCadinu2024()
#' suppressPackageStartupMessages({
#'     library(HDF5Array)
#'     library(SpatialExperiment)
#' })
#'
#' fname <- system.file("extdata", "MerfishData1k", package="SpaceMosaic")
#' spe <- loadHDF5SummarizedExperiment(fname)
#' 
#' spe <- getPatches(spe, cell_type_column="tier1",
#'                   response_cell_type="Epithelial",
#'                   explanatory_cell_type="Immune",
#'                   mmxpixel=0.000109,
#'                   npatches=10, n_iters=5)
#' spe
#'
#' ## a new reduced dimension data frame called "patches_Epithelial_Immune"
#' ## has been added to the input SpatialExperiment object
#' reducedDimNames(spe)
#'
#' ## using the 'reducedDim()' accessor, we can fetch the patch assignments
#' patches <- reducedDim(spe, "patches_Epithelial_Immune")
#' patches
#'
#' ## the patch-level information is in the metadata slot of that data frame
#' metadata(patches)
#'
#' @importFrom cli cli_abort cli_progress_bar cli_progress_update
#' @importFrom cli cli_progress_done
#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment reducedDimNames reducedDim "reducedDim<-"
#' @importFrom Matrix rowSums
#' @importFrom S4Vectors DataFrame Rle "metadata<-"
#' @importFrom BiocGenerics sd var
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
    verbose = TRUE
  ) {
  
    if (!cell_type_column %in% colnames(colData(spe))) {
        msg <- "Cell type column '{cell_type_column}' missing from colData(spe)."
        cli_abort(msg)
    }
    if (sum(colData(spe)[, cell_type_column] == response_cell_type) == 0) {
      msg <- paste("No cells of the specified response cell type",
                  "'{response_cell_type}' were found in colData(spe).")
      cli_abort(msg)
    }
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
    if (!is.null(mmxpixel))
      xy <- xy * mmxpixel

    mask_response_cells <- colData(spe)[, cell_type_column] == response_cell_type
    xy <- xy[mask_response_cells, ]
    X <- env[mask_response_cells, explanatory_cell_type]

    ## scale X
    X <- X / sd(X, na.rm=TRUE)

    ## get initial patch assignments
    seeds <- .initializePatches(xy, X, npatches, initwithhotspots)

    ## set up data frame to track cells
    celldf <- data.frame(xy=xy, X=X)
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

    if (verbose)
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

        ppxhunger <- cell2patchproximity[, rownames(patchdf)] %*% diag(patchdf$hunger)
        idx <- apply(ppxhunger, 1, which.max)

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

        if (verbose)
            cli_progress_update()
    }
    if (verbose)
        cli_progress_done()

    ## put the results into the original SPE object as a reduced dimension data frame
    celldfall <- DataFrame(response=Rle(mask_response_cells),
                           X=NA_real_,
                           patch=Rle(NA_character_),
                           toofar=Rle(NA))
    rownames(celldfall) <- colnames(spe)
    celldfall[mask_response_cells, "X"] <- celldf$X
    celldfall[mask_response_cells, "patch"] <- Rle(factor(celldf$patch))
    celldfall[mask_response_cells, "toofar"] <- Rle(factor(celldf$toofar))
    md <- list(mmxpixel=mmxpixel, patchdf=DataFrame(patchdf))
    metadata(celldfall) <- md
    patchesID <- sprintf("patches_%s_%s", response_cell_type,
                         explanatory_cell_type)
    reducedDim(spe, patchesID) <- celldfall

    return(spe)
}

#' @importFrom spdep knearneigh knn2nb nb2listw localG
#' @importFrom stats kmeans
.initializePatches <- function(xy, X, npatches, initwithhotspots) {
    seeds <- rep(NA, nrow(xy))
    if (initwithhotspots) {
        ## get local X, X^2, var
        neighbors <- .nearestNeighborGraph(xy, N=50)
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
        ishot <- gi >= 1
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

#' @importFrom FNN get.knnx
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
      anydistmat[, name] <- get.knnx(
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

#' @title Plot patches
#'
#' @description Visualize the spatial distribution of the patches computed by
#' `getPatches()` and their statistical power.
#'
#' @param spe An
#' [`SpatialExperiment`][SpatialExperiment::SpatialExperiment-class] object.
#'
#' @param patchesID A character string specifying the name of the reduced
#' dimension data frame containing the patch assignments returned by
#' `getPatches()`. This identifier should be one of the character strings
#' returned by `reducedDimNames(spe)`, and should be of the form
#' `patches_<response_cell_type>_<explanatory_cell_type>`.
#'
#' @param what A character string specifying what to plot: `assignments` plots
#' the spatial distribution of the patch assignments to cells, `power` to plot
#' the predictor sum of squares for each patch, and `both` to plot both of
#' these things. Default value is `both`.
#'
#' @param with A character string specifying which plotting system to use,
#' `ggplot` (default) or `base`. If `ggplot` is chosen but packages
#' [`ggplot2`][ggplot2::ggplot2-package] and
#' [`cowplot`][cowplot::cowplot-package] are not installed, the function will
#' fall back to `base` plotting and issue a message to the user.
#'
#' @return If `with="ggplot"`, the default value, a `ggplot` object is returned.
#' If `with="base"`, the function does not return anything, and it produces the
#' plot.
#'
#' @examples
#'
#'
#' ## load a subset of 1k cells of the first sample and slice of the MERFISH
#' ## data available in MerfishData::MouseColonIbdCadinu2024()
#' library(HDF5Array)
#' library(SpatialExperiment)
#'
#' fname <- system.file("extdata", "MerfishData1k", package="SpaceMosaic")
#' spe <- loadHDF5SummarizedExperiment(fname)
#' 
#' spe <- getPatches(spe, cell_type_column="tier1",
#'                   response_cell_type="Epithelial",
#'                   explanatory_cell_type="Immune",
#'                   mmxpixel=0.000109,
#'                   npatches=10, n_iters=5)
#'
#' \dontrun{
#' plotPatches(spe, patchesID="patches_Epithelial_Immune")
#' }
#'
#' @importFrom utils installed.packages
#' @importFrom cli cli_abort cli_alert_info
#' @importFrom S4Vectors decode metadata
#' @importFrom SingleCellExperiment reducedDimNames reducedDim "reducedDim<-"
#'
#' @rdname plotPatches
#'
#' @export
plotPatches <- function(
    spe,
    patchesID,
    what=c("both", "assignments", "power"),
    with=c("ggplot", "base")) {

    if (!patchesID %in% reducedDimNames(spe)) {
        msg <- "Patches ID '{patchesID}' not found in reducedDimNames(spe)."
        cli_abort(c("x"=msg))
    }
    what <- match.arg(what)
    with <- match.arg(with)

    if (with == "ggplot") {
        instpkgs <- installed.packages(noCache=TRUE)[, "Package"]
        if (!"ggplot2" %in% instpkgs) {
            cli_alert_info("ggplot2 is not installed.")
            with <- "base"
        }
        if (!"cowplot" %in% instpkgs) {
            cli_alert_info("cowplot is not installed.")
            with <- "base"
        }
        if (with == "base")
            cli_alert_info("Falling back to base plotting.")
    }

    celldf <- reducedDim(spe, patchesID)
    mmxpixel <- metadata(celldf)$mmxpixel
    xy <- spatialCoords(spe)[decode(celldf$response), , drop=FALSE]
    if (!is.null(mmxpixel))
        xy <- xy * mmxpixel

    celldf <- celldf[celldf$response, , drop=FALSE]
    if (nrow(celldf) == 0)
        cli_abort(c("x"="No rows in reducedDim(spe, '{patchesID}')."))

    patchdf <- metadata(celldf)$patchdf

    if (with == "ggplot")
        .plot_patches_ggplot(celldf, xy, patchdf, what)
      else
        .plot_patches_base(celldf, xy, patchdf, what)
}

#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
.plot_patches_ggplot <- function(celldf, xy, patchdf, what) {
    if (!.isPackageLoaded("ggplot2")) {
        loaded <- suppressPackageStartupMessages(requireNamespace("ggplot2",
                                                                  quietly=TRUE))
        if (!loaded)
            cli_abort(c("x"="ggplot2 could not be loaded"))
    }
    if (!.isPackageLoaded("cowplot")) {
        loaded <- suppressPackageStartupMessages(requireNamespace("cowplot",
                                                                  quietly=TRUE))
        if (!loaded)
            cli_abort(c("x"="cowplot could not be loaded"))
    }

    ## order and color patches by power for plotting
    ord <- order(patchdf$totvar)
    patchesbytotvar <- factor(rownames(patchdf), levels=rownames(patchdf)[ord])
    patchcols <- colorRampPalette(brewer.pal(8, "Set1"))(nrow(patchdf))
    names(patchcols) <- levels(patchesbytotvar)

    .data <- ggplot2::.data
    p1 <- ggplot2::ggplot(data=celldf, ggplot2::aes(x=xy[, "x"], y=xy[, "y"],
                                               color=.data$patch)) +
              ggplot2::geom_point(size=0.5) +
              ggplot2::scale_color_manual(values=patchcols) + 
              ggplot2::theme_bw() +
              ggplot2::theme(## aspect.ratio=1,
                             panel.grid.major=ggplot2::element_blank(),
                             panel.grid.minor=ggplot2::element_blank(),
                             axis.title=ggplot2::element_text(size=14),
                             axis.text=ggplot2::element_text(size=12),
                             legend.position="none") +
              ggplot2::labs(x="X", y="Y")

    p2 <- ggplot2::ggplot(data=patchdf, ggplot2::aes(x=patchesbytotvar,
                                                     y=.data$totvar,
                                                fill=patchesbytotvar)) +
              ggplot2::geom_bar(stat="identity") +
              ggplot2::scale_fill_manual(values=patchcols) +
              ggplot2::theme_bw() +
              ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank(),
                             panel.border = ggplot2::element_blank(),
                             axis.line.y = ggplot2::element_line(color="black"),
                             axis.text.x = ggplot2::element_blank(),
                             axis.ticks.x = ggplot2::element_blank(),
                             axis.title.x = ggplot2::element_text(size=14,
                                                  margin=ggplot2::margin(t=20)),
                             axis.title.y = ggplot2::element_text(size=14),
                             axis.text = ggplot2::element_text(size=12),
                             legend.position = "none",
                             plot.margin = ggplot2::margin(20, 10, 20, 10)) +
                             ggplot2::scale_y_continuous(expand = c(0, 0)) +
              ggplot2::labs(x="Patches", y="Predictor sum of squares")

    if (what == "assignments")
        return(p1)
    else if (what == "power")
        return(p2)

    ## both
    p <- cowplot::plot_grid(p1, p2, align="v", ncol=1)

    return(p)
}

#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics par points barplot mtext
.plot_patches_base <- function(celldf, xy, patchdf, what) {
    ## order and color patches by power for plotting
    patchdf <- patchdf[order(patchdf$totvar), ]
    patchcols <- colorRampPalette(brewer.pal(8, "Set1"))(nrow(patchdf))
    names(patchcols) <- rownames(patchdf)

    if (what == "both")
        par(mfrow=c(2, 1))

    if (what %in% c("both", "assignments")) {
        if (what == "assignments")
            par(mfrow=c(1, 1))
        par(mar=c(4, 4, 2, 1))
        plot(xy, asp=1, pch=16, cex=0.1, col="grey80", main="", las=1,
             cex.axis=1.2, cex.lab=1.5, xlab="X", ylab="Y")
        points(xy, pch=16, cex=0.4, col=patchcols[decode(celldf$patch)])
    }

    if (what %in% c("both", "power")) {
        if (what == "power")
            par(mfrow=c(1, 1))
        par(mar=c(3, 5, 2, 1))
        barplot(patchdf$totvar, col=patchcols[rownames(patchdf)], las=1,
                ylab="Predictor sum of squares", xlab="", cex.axis=1.2,
                cex.lab=1.5)
        mtext("Patches", side=1, line=1.5, cex=1.5)
    }
}
