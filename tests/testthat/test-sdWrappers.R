# These tests check that every `.sd`-suffixed wrapper in R/sdWrappers.R
# produces results identical to the corresponding `.spe`-suffixed wrapper in
# R/speWrappers.R when both are run on the same data: a `SpatialExperiment`
# and the `SpatialData` object built from it with `spe_to_spatialdata()`. As
# in test-speWrappers.R, the bundled carcinoma dataset is subset to a small,
# arbitrary set of cells and genes purely to keep the tests fast.

skip_if_not_installed("spatialdataR")
skip_if_not_installed("sf")

spe_path <- system.file("extdata", "cosmx_carcinoma.rds", package = "SpaceMosaic")

spe_full <- readRDS(spe_path)
set.seed(42)
cells_use <- sample(colnames(spe_full), 400)
genes_use <- sample(rownames(spe_full), 30)
spe <- spe_full[genes_use, cells_use]
sd <- spe_to_spatialdata(spe)

# Table of a SpatialData object, by name or index.
sd_table <- function(sd, table = 1) spatialdataR::tables(sd)[[table]]

# Add a patch column to both the SpatialExperiment and its SpatialData copy.
with_patch <- function(patch) {
  spe_patch <- spe
  SummarizedExperiment::colData(spe_patch)$patch <- patch
  list(spe = spe_patch, sd = spe_to_spatialdata(spe_patch))
}

test_that("spe_to_spatialdata keeps the table and coordinates of the SpatialExperiment", {
  sce <- sd_table(sd)

  expect_s4_class(sd, "SpatialData")
  expect_identical(spatialdataR::tableNames(sd), "table")
  expect_identical(spatialdataR::shapeNames(sd), "shape1")
  expect_identical(colnames(sce), colnames(spe))
  expect_identical(rownames(sce), rownames(spe))
  expect_identical(
    SummarizedExperiment::assayNames(sce), SummarizedExperiment::assayNames(spe)
  )
  expect_equal(
    SummarizedExperiment::assay(sce, "logcounts"),
    SummarizedExperiment::assay(spe, "logcounts")
  )
  expect_equal(
    SummarizedExperiment::colData(sce), SummarizedExperiment::colData(spe)
  )
  expect_equal(
    SingleCellExperiment::reducedDims(sce), SingleCellExperiment::reducedDims(spe)
  )
  expect_equal(
    unname(SpaceMosaic:::.sdCoords(sd, sce)),
    unname(SpatialExperiment::spatialCoords(spe))
  )
  expect_identical(rownames(SpaceMosaic:::.sdCoords(sd, sce)), colnames(spe))
})

test_that("coordinates are matched to cells by instance ID, not by position", {
  set.seed(1)
  perm <- sample(ncol(spe))
  sd_perm <- sd
  table_perm <- sd_table(sd_perm)[, perm]
  spatialdataR::table(sd_perm) <- table_perm

  expect_equal(
    unname(SpaceMosaic:::.sdCoords(sd_perm, sd_table(sd_perm))),
    unname(SpatialExperiment::spatialCoords(spe)[perm, ])
  )
})

test_that("coordinates match the SpatialExperiment on the full dataset", {
  # DuckDB does not preserve row order on larger tables, so cells must be
  # matched to geometries by instance ID rather than by position.
  sd_full <- spe_to_spatialdata(spe_full)
  for (i in 1:2) {
    expect_equal(
      unname(SpaceMosaic:::.sdCoords(sd_full, sd_table(sd_full))),
      unname(SpatialExperiment::spatialCoords(spe_full))
    )
  }
})

test_that("creating a new SpatialData object leaves earlier ones untouched", {
  # spatialdataR::SpatialDataShape() writes every data.frame to one shared
  # DuckDB table; spe_to_spatialdata() must not be affected by this.
  spe_sub <- spe[, spe$celltype == "Cancer.cells"]
  sd_a <- spe_to_spatialdata(spe)
  sd_b <- spe_to_spatialdata(spe_sub)
  invisible(spatialdataR::SpatialDataShape(data.frame(x = 1:3, y = 1:3)))

  expect_equal(
    unname(SpaceMosaic:::.sdCoords(sd_a, sd_table(sd_a))),
    unname(SpatialExperiment::spatialCoords(spe))
  )
  expect_equal(
    unname(SpaceMosaic:::.sdCoords(sd_b, sd_table(sd_b))),
    unname(SpatialExperiment::spatialCoords(spe_sub))
  )

  # Also when both are added to the same SpatialData object.
  sd_two <- spe_to_spatialdata(
    spe_sub, sd = sd_a, shape_name = "shape2", table_name = "table2"
  )
  expect_equal(
    unname(SpaceMosaic:::.sdCoords(sd_two, sd_table(sd_two, "table"))),
    unname(SpatialExperiment::spatialCoords(spe))
  )
  expect_equal(
    unname(SpaceMosaic:::.sdCoords(sd_two, sd_table(sd_two, "table2"))),
    unname(SpatialExperiment::spatialCoords(spe_sub))
  )
})

test_that("embedCellNeighborhoods.sd matches embedCellNeighborhoods.spe", {
  expected <- embedCellNeighborhoods.spe(
    spe, embedding = "PCA", ks = c(5, 20), tissue = "sample_id", name = "Z_test"
  )
  result <- embedCellNeighborhoods.sd(
    sd, embedding = "PCA", ks = c(5, 20), tissue = "sample_id", name = "Z_test"
  )

  expect_s4_class(result, "SpatialData")
  expect_equal(
    SingleCellExperiment::reducedDim(sd_table(result), "Z_test"),
    SingleCellExperiment::reducedDim(expected, "Z_test")
  )
})

test_that("getPatches.sd matches getPatches.spe", {
  set.seed(123)
  expected <- getPatches.spe(
    spe, X = "distance_to_fibroblast", npatches = 4, n_iters = 3, verbose = FALSE
  )
  set.seed(123)
  result <- getPatches.sd(
    sd, X = "distance_to_fibroblast", npatches = 4, n_iters = 3, verbose = FALSE
  )
  result_table <- sd_table(result)

  expect_identical(result_table$patch, expected$patch)

  # Everything but the polygons is stored in the table's metadata ...
  expected_metadata <- S4Vectors::metadata(expected)$SpaceMosaic
  expect_equal(
    S4Vectors::metadata(result_table)$SpaceMosaic,
    expected_metadata[setdiff(names(expected_metadata), "patch_polys")]
  )

  # ... while the polygons are stored as a shapes element, one per patch.
  expected_polys <- expected_metadata$patch_polys
  result_polys <- spatialdataR::data(spatialdataR::shapes(result)$patch_polys) |>
    sf::st_as_sf()
  expect_identical(result_polys$patch, unique(expected_polys$patch))
  for (id in result_polys$patch) {
    vertices <- sf::st_coordinates(result_polys[result_polys$patch == id, ])
    hull <- as.matrix(expected_polys[expected_polys$patch == id, c("x", "y")])
    expect_equal(unname(vertices[, c("X", "Y")]), unname(rbind(hull, hull[1, ])))
  }
  annotations <- setdiff(names(expected_polys), c("x", "y"))
  expect_equal(
    as.data.frame(sf::st_drop_geometry(result_polys))[, annotations],
    `rownames<-`(expected_polys[!duplicated(expected_polys$patch), annotations], NULL)
  )
})

test_that("getPatches.sd stores patch polygons under `patch_shape`", {
  set.seed(123)
  result <- getPatches.sd(
    sd, X = "distance_to_fibroblast", npatches = 4, n_iters = 3,
    patch_shape = "my_patches", verbose = FALSE
  )

  expect_true("my_patches" %in% spatialdataR::shapeNames(result))
  expect_false("patch_polys" %in% spatialdataR::shapeNames(result))
  # Cell coordinates are untouched by adding the polygons.
  expect_equal(
    unname(SpaceMosaic:::.sdCoords(result, sd_table(result))),
    unname(SpatialExperiment::spatialCoords(spe))
  )
})

test_that("getPatches.sd with Z matches getPatches.spe", {
  set.seed(123)
  expected <- getPatches.spe(
    spe, X = "distance_to_fibroblast", npatches = 4, Z = "Z", n_iters = 3,
    verbose = FALSE
  )
  set.seed(123)
  result <- getPatches.sd(
    sd, X = "distance_to_fibroblast", npatches = 4, Z = "Z", n_iters = 3,
    verbose = FALSE
  )

  expect_identical(sd_table(result)$patch, expected$patch)
})

test_that("getPatches.sd on a cell-subsetted object matches getPatches.spe", {
  spe_sub <- spe[, spe$celltype == "Cancer.cells"]
  sd_sub <- spe_to_spatialdata(spe_sub)

  set.seed(7)
  expected <- getPatches.spe(
    spe_sub, X = "distance_to_fibroblast", npatches = 4, n_iters = 3,
    verbose = FALSE
  )
  set.seed(7)
  result <- getPatches.sd(
    sd_sub, X = "distance_to_fibroblast", npatches = 4, n_iters = 3,
    verbose = FALSE
  )

  expect_identical(sd_table(result)$patch, expected$patch)
})

test_that("getPatches.sd skips metadata storage when log_iters/patch_diagnostics/patch_polys are FALSE", {
  set.seed(123)
  expected <- getPatches.spe(
    spe, X = "distance_to_fibroblast", npatches = 4, n_iters = 3,
    log_iters = FALSE, patch_diagnostics = FALSE, patch_polys = FALSE,
    verbose = FALSE
  )
  set.seed(123)
  result <- getPatches.sd(
    sd, X = "distance_to_fibroblast", npatches = 4, n_iters = 3,
    log_iters = FALSE, patch_diagnostics = FALSE, patch_polys = FALSE,
    verbose = FALSE
  )
  result_table <- sd_table(result)

  expect_identical(result_table$patch, expected$patch)
  expect_null(S4Vectors::metadata(result_table)$SpaceMosaic)
  expect_identical(spatialdataR::shapeNames(result), "shape1")
})

test_that("patchDE.sd matches patchDE.spe", {
  obj <- with_patch(rep(c("1", "2"), length.out = ncol(spe)))

  expected <- patchDE.spe(
    obj$spe, predictor_cols = "distance_to_fibroblast", method = "limma",
    return_residuals = TRUE, verbose = FALSE
  )
  result <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", method = "limma",
    return_residuals = TRUE, verbose = FALSE
  )

  expect_named(result, c("de", "sd", "table"))
  expect_identical(result$table, "table")
  expect_equal(result$de, expected$de)
  expect_equal(
    SummarizedExperiment::assay(sd_table(result$sd), "residuals"),
    SummarizedExperiment::assay(expected$spe, "residuals")
  )
})

test_that("patchDE.sd with return_residuals = FALSE returns unwrapped `de` and leaves assays untouched", {
  obj <- with_patch(rep(c("1", "2"), length.out = ncol(spe)))

  # method defaults to "hasty".
  expected <- patchDE.spe(
    obj$spe, predictor_cols = "distance_to_fibroblast", verbose = FALSE
  )
  result <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", verbose = FALSE
  )

  expect_equal(result$de, expected$de)
  expect_identical(
    SummarizedExperiment::assayNames(sd_table(result$sd)),
    SummarizedExperiment::assayNames(obj$spe)
  )
})

test_that("patchDE.sd with pearson = TRUE matches patchDE.spe", {
  obj <- with_patch(rep(c("1", "2"), length.out = ncol(spe)))

  expected <- patchDE.spe(
    obj$spe, predictor_cols = "distance_to_fibroblast", assay = "counts",
    pearson = TRUE, tot = "total", verbose = FALSE
  )
  result <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", assay = "counts",
    pearson = TRUE, tot = "total", verbose = FALSE
  )

  expect_equal(result$de, expected$de)
})

test_that("patchMetaAnalysis.sd matches patchMetaAnalysis.spe", {
  obj <- with_patch(as.character(rep(1:3, length.out = ncol(spe))))

  de_spe <- patchDE.spe(
    obj$spe, predictor_cols = "distance_to_fibroblast", method = "limma",
    verbose = FALSE
  )
  de_sd <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", method = "limma",
    verbose = FALSE
  )

  expected <- patchMetaAnalysis.spe(de_spe, embedding_name = "Z")
  result <- patchMetaAnalysis.sd(de_sd, embedding_name = "Z")

  expect_equal(result, expected)
})

test_that("patchMetaAnalysis.sd with summarize_subgroups = TRUE matches patchMetaAnalysis.spe", {
  obj <- with_patch(as.character(rep(1:3, length.out = ncol(spe))))

  de_spe <- patchDE.spe(
    obj$spe, predictor_cols = "distance_to_fibroblast", method = "limma",
    verbose = FALSE
  )
  de_sd <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", method = "limma",
    verbose = FALSE
  )

  expected <- patchMetaAnalysis.spe(
    de_spe, embedding_name = "Z", summarize_subgroups = TRUE,
    cellmeta_cols = "distance_to_fibroblast"
  )
  result <- patchMetaAnalysis.sd(
    de_sd, embedding_name = "Z", summarize_subgroups = TRUE,
    cellmeta_cols = "distance_to_fibroblast"
  )

  expect_equal(result$subgroups_summary, expected$subgroups_summary)
})

test_that("moranTest.sd matches moranTest.spe", {
  obj <- with_patch(rep(c("1", "2"), length.out = ncol(spe)))

  de_spe <- patchDE.spe(
    obj$spe, predictor_cols = "distance_to_fibroblast", method = "limma",
    return_residuals = TRUE, verbose = FALSE
  )
  de_sd <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", method = "limma",
    return_residuals = TRUE, verbose = FALSE
  )

  set.seed(99)
  expected <- moranTest.spe(de_spe$spe, k = 5, n_permutations = 19)
  set.seed(99)
  result <- moranTest.sd(de_sd$sd, k = 5, n_permutations = 19)

  expect_equal(result, expected)
})

test_that("moranTest.sd on a gene-subsetted table matches moranTest.spe", {
  obj <- with_patch(rep(c("1", "2"), length.out = ncol(spe)))

  de_spe <- patchDE.spe(
    obj$spe, predictor_cols = "distance_to_fibroblast", method = "limma",
    return_residuals = TRUE, verbose = FALSE
  )
  de_sd <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", method = "limma",
    return_residuals = TRUE, verbose = FALSE
  )
  sd_gene_sub <- de_sd$sd
  spatialdataR::table(sd_gene_sub) <- sd_table(sd_gene_sub)[1:10, ]

  set.seed(99)
  expected <- moranTest.spe(de_spe$spe[1:10, ], k = 5, n_permutations = 19)
  set.seed(99)
  result <- moranTest.sd(sd_gene_sub, k = 5, n_permutations = 19)

  expect_equal(result, expected)
})

test_that("the `table` argument selects which table is used and modified", {
  spe_sub <- spe[, spe$celltype == "Cancer.cells"]
  sd_two <- spe_to_spatialdata(
    spe_sub, sd = sd, shape_name = "shape2", table_name = "table2"
  )
  expect_identical(spatialdataR::tableNames(sd_two), c("table", "table2"))

  set.seed(7)
  expected <- getPatches.spe(
    spe_sub, X = "distance_to_fibroblast", npatches = 4, n_iters = 3,
    verbose = FALSE
  )
  set.seed(7)
  by_name <- getPatches.sd(
    sd_two, X = "distance_to_fibroblast", npatches = 4, n_iters = 3,
    table = "table2", verbose = FALSE
  )
  set.seed(7)
  by_index <- getPatches.sd(
    sd_two, X = "distance_to_fibroblast", npatches = 4, n_iters = 3,
    table = 2, verbose = FALSE
  )

  expect_identical(sd_table(by_name, "table2")$patch, expected$patch)
  expect_identical(sd_table(by_index, "table2")$patch, expected$patch)
  expect_null(sd_table(by_name, "table")$patch)

  # patchDE.sd records the table, and patchMetaAnalysis.sd reuses it.
  de_sd <- patchDE.sd(
    by_name, predictor_cols = "distance_to_fibroblast", table = "table2",
    verbose = FALSE
  )
  de_spe <- patchDE.spe(
    expected, predictor_cols = "distance_to_fibroblast", verbose = FALSE
  )
  expect_identical(de_sd$table, "table2")
  expect_equal(de_sd$de, de_spe$de)
  expect_equal(
    patchMetaAnalysis.sd(de_sd, embedding_name = "Z"),
    patchMetaAnalysis.spe(de_spe, embedding_name = "Z")
  )
})

test_that(".sd wrappers give informative errors on invalid input", {
  expect_error(getPatches.sd(spe, X = "distance_to_fibroblast", npatches = 4),
               "must be a SpatialData object")
  expect_error(getPatches.sd(sd, X = "distance_to_fibroblast", npatches = 4,
                             table = "nope"),
               "not found in tableNames")
  expect_error(getPatches.sd(sd, X = "distance_to_fibroblast", npatches = 4,
                             table = 5),
               "out of range")
  expect_error(getPatches.sd(sd, X = "nope", npatches = 4),
               "not columns of colData")
  expect_error(embedCellNeighborhoods.sd(sd, embedding = "nope"),
               "not found in reducedDimNames")
  expect_error(patchDE.sd(sd, predictor_cols = "distance_to_fibroblast",
                          patch_column = "nope"),
               "not found in colData")
  expect_error(patchMetaAnalysis.sd(list(de = list())),
               "must be a list with `de` and `sd` elements")
  expect_error(moranTest.sd(sd), "not found in assay names")
})
