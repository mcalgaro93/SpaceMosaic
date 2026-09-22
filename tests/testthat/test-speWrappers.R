# These tests check that every `.spe`-suffixed wrapper in R/speWrappers.R
# produces results identical to calling the corresponding matrix-based
# function directly on the pieces extracted from a `SpatialExperiment`. The
# bundled carcinoma dataset is subset to a small, arbitrary set of cells and
# genes purely to keep the tests fast; the subset is not meant to be a
# biologically meaningful patch layout.

spe_path <- system.file("extdata", "cosmx_carcinoma.rds", package = "SpaceMosaic")

spe_full <- readRDS(spe_path)
set.seed(42)
cells_use <- sample(colnames(spe_full), 400)
genes_use <- sample(rownames(spe_full), 30)
spe <- spe_full[genes_use, cells_use]

test_that("embedCellNeighborhoods.spe matches embedCellNeighborhoods on the matrix workflow", {
  expected <- embedCellNeighborhoods(
    mat = SingleCellExperiment::reducedDim(spe, "PCA"),
    xy = SpatialExperiment::spatialCoords(spe),
    ks = c(5, 20),
    tissue = SummarizedExperiment::colData(spe)[["sample_id"]]
  )

  result <- embedCellNeighborhoods.spe(
    spe, embedding = "PCA", ks = c(5, 20), tissue = "sample_id", name = "Z_test"
  )
  result_Z <- SingleCellExperiment::reducedDim(result, "Z_test")

  expect_equal(unname(result_Z), unname(expected))
  expect_identical(colnames(result_Z), colnames(expected))
  expect_identical(rownames(result_Z), colnames(spe))
})

test_that("getPatches.spe matches getPatches/getPatchDiagnostics/getPatchPolys on the matrix workflow", {
  X_col <- "distance_to_fibroblast"
  npatches <- 4

  xy <- SpatialExperiment::spatialCoords(spe)
  X_mat <- as.matrix(SummarizedExperiment::colData(spe)[, X_col, drop = FALSE])

  set.seed(123)
  expected_patches <- getPatches(
    xy = xy, X = X_mat, npatches = npatches, n_iters = 3, verbose = FALSE
  )

  set.seed(123)
  result <- getPatches.spe(
    spe, X = X_col, npatches = npatches, n_iters = 3, verbose = FALSE
  )

  expect_equal(
    unname(SummarizedExperiment::colData(result)$patch),
    unname(expected_patches$patch[colnames(spe)])
  )

  expected_diagnostics <- getPatchDiagnostics(
    xy = xy, X = X_mat, patches = expected_patches, k = 10
  )
  result_metadata <- S4Vectors::metadata(result)$SpaceMosaic

  expect_equal(
    result_metadata$patch_diagnostics, expected_diagnostics$patch_diagnostics
  )
  expect_equal(
    result_metadata$assignment_summary, expected_diagnostics$assignment_summary
  )
  expect_equal(
    result_metadata$connectivity_curve, expected_diagnostics$connectivity_curve
  )

  expected_polys <- getPatchPolys(
    xy = xy,
    patch = expected_patches$patch[colnames(spe)],
    patch_data = expected_diagnostics$patch_diagnostics
  )
  expect_equal(result_metadata$patch_polys, expected_polys)
})

test_that("patchDE.spe matches patchDE on the matrix workflow", {
  patch <- rep(c("1", "2"), length.out = ncol(spe))
  spe_patch <- spe
  SummarizedExperiment::colData(spe_patch)$patch <- patch

  y <- t(SummarizedExperiment::assay(spe_patch, "logcounts"))
  df <- SummarizedExperiment::colData(spe_patch)[, "distance_to_fibroblast", drop = FALSE]

  expected <- patchDE(
    y, df, patch, method = "limma", return_residuals = TRUE, verbose = FALSE
  )

  result <- patchDE.spe(
    spe_patch,
    predictor_cols = "distance_to_fibroblast",
    method = "limma",
    return_residuals = TRUE,
    verbose = FALSE
  )

  expect_equal(result$de, expected$de)
  expect_equal(
    unname(SummarizedExperiment::assay(result$spe, "residuals")),
    unname(t(expected$residuals))
  )
})

test_that("patchMetaAnalysis.spe matches getPatchAttributes/patchMetaAnalysis on the matrix workflow", {
  patch <- as.character(rep(1:3, length.out = ncol(spe)))
  spe_patch <- spe
  SummarizedExperiment::colData(spe_patch)$patch <- patch

  de_result <- patchDE.spe(
    spe_patch,
    predictor_cols = "distance_to_fibroblast",
    method = "limma",
    verbose = FALSE
  )

  result <- patchMetaAnalysis.spe(de_result, embedding_name = "Z")

  Z <- SingleCellExperiment::reducedDim(spe_patch, "Z")
  W <- getPatchAttributes(Z, patch)
  de_patches <- colnames(de_result$de[[1]]$pvals)
  W <- W[de_patches, , drop = FALSE]
  expected_meta <- patchMetaAnalysis(de_result$de, W)

  expect_equal(result$meta, expected_meta)
  expect_equal(result$W, W)
})

test_that("moranTest.spe matches moranTest on the matrix workflow", {
  patch <- rep(c("1", "2"), length.out = ncol(spe))
  spe_patch <- spe
  SummarizedExperiment::colData(spe_patch)$patch <- patch

  de_result <- patchDE.spe(
    spe_patch,
    predictor_cols = "distance_to_fibroblast",
    method = "limma",
    return_residuals = TRUE,
    verbose = FALSE
  )

  set.seed(99)
  result <- moranTest.spe(de_result$spe, k = 5, n_permutations = 19)

  expected_residuals <- t(SummarizedExperiment::assay(de_result$spe, "residuals"))
  set.seed(99)
  expected <- moranTest(
    residuals = expected_residuals,
    xy = SpatialExperiment::spatialCoords(de_result$spe),
    patch = SummarizedExperiment::colData(de_result$spe)$patch,
    k = 5,
    n_permutations = 19
  )
  expect_equal(result, expected)
})

test_that("getPatches.spe on a further cell-subsetted SpatialExperiment matches the matrix workflow", {
  # Subsets spe by a boolean mask (as in the getPatches.spe/patchDE.spe
  # examples), rather than relying on the fixed sample already used to build
  # `spe`, to check that patching still lines up correctly on an arbitrary
  # subset of cells.
  spe_sub <- spe[, spe$celltype == "Cancer.cells"]
  npatches <- 4

  set.seed(7)
  result <- getPatches.spe(
    spe_sub, X = "distance_to_fibroblast", npatches = npatches, n_iters = 3,
    verbose = FALSE
  )

  xy <- SpatialExperiment::spatialCoords(spe_sub)
  X_mat <- as.matrix(
    SummarizedExperiment::colData(spe_sub)[, "distance_to_fibroblast", drop = FALSE]
  )
  set.seed(7)
  expected <- getPatches(
    xy = xy, X = X_mat, npatches = npatches, n_iters = 3, verbose = FALSE
  )

  expect_equal(
    unname(SummarizedExperiment::colData(result)$patch),
    unname(expected$patch[colnames(spe_sub)])
  )
})

test_that("moranTest.spe on a gene-subsetted SpatialExperiment matches the matrix workflow", {
  # Mirrors the row-subsetting shown in moranTest.spe's own @examples
  # (`de_result$spe[1:10, ]`), checking that subsetting genes after patchDE.spe
  # correctly subsets the residuals assay passed on to moranTest.spe.
  patch <- rep(c("1", "2"), length.out = ncol(spe))
  spe_patch <- spe
  SummarizedExperiment::colData(spe_patch)$patch <- patch

  de_result <- patchDE.spe(
    spe_patch,
    predictor_cols = "distance_to_fibroblast",
    method = "limma",
    return_residuals = TRUE,
    verbose = FALSE
  )
  spe_gene_sub <- de_result$spe[1:10, ]

  set.seed(99)
  result <- moranTest.spe(spe_gene_sub, k = 5, n_permutations = 19)

  expected_residuals <- t(SummarizedExperiment::assay(spe_gene_sub, "residuals"))
  set.seed(99)
  expected <- moranTest(
    residuals = expected_residuals,
    xy = SpatialExperiment::spatialCoords(spe_gene_sub),
    patch = SummarizedExperiment::colData(spe_gene_sub)$patch,
    k = 5,
    n_permutations = 19
  )
  expect_equal(result, expected)
})

test_that("getPatches.spe skips metadata storage when log_iters/patch_diagnostics/patch_polys are FALSE", {
  npatches <- 4

  set.seed(123)
  result <- getPatches.spe(
    spe, X = "distance_to_fibroblast", npatches = npatches, n_iters = 3,
    log_iters = FALSE, patch_diagnostics = FALSE, patch_polys = FALSE,
    verbose = FALSE
  )

  xy <- SpatialExperiment::spatialCoords(spe)
  X_mat <- as.matrix(SummarizedExperiment::colData(spe)[, "distance_to_fibroblast", drop = FALSE])
  set.seed(123)
  expected_patch <- getPatches(
    xy = xy, X = X_mat, npatches = npatches, n_iters = 3, log_iters = FALSE,
    verbose = FALSE
  )

  expect_equal(
    unname(SummarizedExperiment::colData(result)$patch),
    unname(expected_patch[colnames(spe)])
  )
  expect_null(S4Vectors::metadata(result)$SpaceMosaic)
})

test_that("patchDE.spe with return_residuals = FALSE returns unwrapped `de` and leaves assays untouched", {
  patch <- rep(c("1", "2"), length.out = ncol(spe))
  spe_patch <- spe
  SummarizedExperiment::colData(spe_patch)$patch <- patch

  # method defaults to "hasty", also exercising that untested backend.
  result <- patchDE.spe(
    spe_patch, predictor_cols = "distance_to_fibroblast", verbose = FALSE
  )

  y <- t(SummarizedExperiment::assay(spe_patch, "logcounts"))
  df <- SummarizedExperiment::colData(spe_patch)[, "distance_to_fibroblast", drop = FALSE]
  expected <- patchDE(y, df, patch, method = "hasty", verbose = FALSE)

  expect_equal(result$de, expected)
  expect_identical(
    SummarizedExperiment::assayNames(result$spe),
    SummarizedExperiment::assayNames(spe_patch)
  )
})

test_that("patchDE.spe with pearson = TRUE matches patchDE on the matrix workflow", {
  patch <- rep(c("1", "2"), length.out = ncol(spe))
  spe_patch <- spe
  SummarizedExperiment::colData(spe_patch)$patch <- patch

  result <- patchDE.spe(
    spe_patch,
    predictor_cols = "distance_to_fibroblast",
    assay = "counts",
    pearson = TRUE,
    tot = "total",
    verbose = FALSE
  )

  y <- t(SummarizedExperiment::assay(spe_patch, "counts"))
  df <- SummarizedExperiment::colData(spe_patch)[, "distance_to_fibroblast", drop = FALSE]
  expected <- patchDE(
    y, df, patch,
    pearson = TRUE,
    tot = SummarizedExperiment::colData(spe_patch)$total,
    verbose = FALSE
  )

  expect_equal(result$de, expected)
})

test_that("patchMetaAnalysis.spe with summarize_subgroups = TRUE matches summarizeSubgroups on the matrix workflow", {
  patch <- as.character(rep(1:3, length.out = ncol(spe)))
  spe_patch <- spe
  SummarizedExperiment::colData(spe_patch)$patch <- patch

  de_result <- patchDE.spe(
    spe_patch,
    predictor_cols = "distance_to_fibroblast",
    method = "limma",
    verbose = FALSE
  )

  result <- patchMetaAnalysis.spe(
    de_result,
    embedding_name = "Z",
    summarize_subgroups = TRUE,
    cellmeta_cols = "distance_to_fibroblast"
  )

  Z <- SingleCellExperiment::reducedDim(spe_patch, "Z")
  W <- getPatchAttributes(Z, patch)
  de_patches <- colnames(de_result$de[[1]]$pvals)
  W <- W[de_patches, , drop = FALSE]
  expected_meta <- patchMetaAnalysis(de_result$de, W)
  expected_subgroups <- summarizeSubgroups(
    expected_meta,
    patch = patch,
    cellmeta = SummarizedExperiment::colData(spe_patch)[, "distance_to_fibroblast", drop = FALSE]
  )

  expect_equal(result$subgroups_summary, expected_subgroups)
})
