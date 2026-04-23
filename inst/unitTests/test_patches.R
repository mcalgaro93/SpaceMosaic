test_getPatches <- function() {
    message("Running unit tests for getPatches()")

    suppressPackageStartupMessages({
        library(HDF5Array)
        library(SpatialExperiment)
    })

    fname <- system.file("extdata", "MerfishData1k", package="SpaceMosaic")
    spe <- loadHDF5SummarizedExperiment(fname)
 
    spe <- getPatches(spe, cell_type_column="tier1",
                      response_cell_type="Epithelial",
                      explanatory_cell_type="Immune",
                      mmxpixel=0.000109,
                      npatches=10, n_iters=5)
    checkTrue("patches_Epithelial_Immune" %in% reducedDimNames(spe))
}

test_plotPatches <- function() {
    message("Running unit tests for plotPatches()")

    suppressPackageStartupMessages({
        library(HDF5Array)
        library(SpatialExperiment)
        library(ggplot2)
    })

    fname <- system.file("extdata", "MerfishData1k", package="SpaceMosaic")
    spe <- loadHDF5SummarizedExperiment(fname)
 
    spe <- getPatches(spe, cell_type_column="tier1",
                      response_cell_type="Epithelial",
                      explanatory_cell_type="Immune",
                      mmxpixel=0.000109,
                      npatches=10, n_iters=5)

    p <- plotPatches(spe, "patches_Epithelial_Immune")
    checkTrue(is(p, "ggplot"))
}
