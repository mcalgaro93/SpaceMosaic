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

# Add a patch column to both the SpatialExperiment and its SpatialData copy,
# and give the latter the patch polygons and patch table that
# getPatches.sd() would create.
with_patch <- function(patch) {
  spe_patch <- spe
  SummarizedExperiment::colData(spe_patch)$patch <- patch
  sd_patch <- spe_to_spatialdata(spe_patch)
  xy <- SpaceMosaic:::.sdCoords(sd_patch, sd_table(sd_patch))
  sd_patch <- SpaceMosaic:::.sdAddPatchTable(
    sd_patch, xy, patch, patch_shape = "patch_polys", patch_table = "patches"
  )
  list(spe = spe_patch, sd = sd_patch)
}

# Assay `name` of the patch table, restricted to the genes and patches of
# `like` in the same order.
patch_assay <- function(sd, name, like) {
  SummarizedExperiment::assay(sd_table(sd, "patches"), name)[
    rownames(like), colnames(like), drop = FALSE]
}

# Check the patch-table assays written by patchDE.sd() against the `de`
# element returned by patchDE.spe().
expect_de_assays <- function(sd, de) {
  for (predictor in names(de)) {
    if (predictor == "resid_mse") next
    stats <- c(est = "ests", se = "ses", pval = "pvals")
    for (stat in names(stats)) {
      like <- de[[predictor]][[stats[[stat]]]]
      expect_equal(patch_assay(sd, paste0(stat, "_", predictor), like), like)
    }
  }
}

# Check the patch-table assays and W written by patchMetaAnalysis.sd()
# against the list returned by patchMetaAnalysis.spe().
expect_meta_results <- function(sd, expected) {
  prefixes <- c(ests = "est", ses = "se", pvals = "pval", subgroups = "subgroups")
  for (predictor in names(expected$meta)) {
    for (stat in names(expected$meta[[predictor]])) {
      like <- expected$meta[[predictor]][[stat]]
      expect_equal(
        patch_assay(sd, paste0("meta_", prefixes[[stat]], "_", predictor), like),
        like
      )
    }
  }
  W <- SingleCellExperiment::reducedDim(sd_table(sd, "patches"), "W")
  expect_equal(W[rownames(expected$W), , drop = FALSE], expected$W)
}

# Check the patch-table assays written by moranTest.sd() against the data
# frame returned by moranTest.spe().
expect_moran_results <- function(sd, expected) {
  patches <- sd_table(sd, "patches")
  idx <- cbind(expected$gene, expected$patch)
  for (stat in c("observed", "expected", "p_value", "p_adjusted", "status")) {
    expect_equal(
      SummarizedExperiment::assay(patches, paste0("moran_", stat))[idx],
      expected[[stat]]
    )
  }
  first <- expected[!duplicated(expected$patch), ]
  expect_equal(patches[, first$patch]$moran_n_cells, first$n_cells)
  expect_equal(patches[, first$patch]$moran_effective_k, first$effective_k)
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
  sd_two <- spatialdataR::combine(
    sd_a,
    spe_to_spatialdata(spe_sub, shape_name = "shape2", table_name = "table2")
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

  # Per-iteration assignments are stored as a reducedDim of the cell table,
  # which keeps no SpaceMosaic metadata ...
  expected_metadata <- S4Vectors::metadata(expected)$SpaceMosaic
  expect_identical(
    SingleCellExperiment::reducedDim(result_table, "patch_iterations"),
    expected_metadata$patch_iterations$membership_log
  )
  expect_null(S4Vectors::metadata(result_table)$SpaceMosaic)

  # ... patch-level results in a table annotating the polygons ...
  patches <- sd_table(result, "patches")
  diagnostics <- expected_metadata$patch_diagnostics
  expect_identical(colnames(patches), diagnostics$patch)
  expect_identical(spatialdataR::region(patches), "patch_polys")
  expect_identical(spatialdataR::instance_key(patches), "patch")
  expect_identical(spatialdataR::instances(patches), diagnostics$patch)
  expect_equal(
    as.data.frame(SummarizedExperiment::colData(patches)),
    `rownames<-`(diagnostics[, names(diagnostics) != "patch"], diagnostics$patch)
  )
  curve <- SingleCellExperiment::reducedDim(patches, "connectivity_curve")
  expected_curve <- expected_metadata$connectivity_curve
  expect_equal(
    curve[cbind(expected_curve$patch, paste0("k", expected_curve$k))],
    expected_curve$component_fraction
  )
  expect_equal(
    SingleCellExperiment::reducedDim(patches, "ss_log"),
    expected_metadata$patch_iterations$ss_log[diagnostics$patch, ]
  )
  expect_equal(
    S4Vectors::metadata(patches)$SpaceMosaic$assignment_summary,
    expected_metadata$assignment_summary
  )

  # ... and the polygons as a shapes element, one per patch, keyed by patch.
  expected_polys <- expected_metadata$patch_polys
  polys_shape <- spatialdataR::shapes(result)$patch_polys
  result_polys <- spatialdataR::data(polys_shape) |> sf::st_as_sf()
  expect_identical(names(sf::st_drop_geometry(result_polys)), "patch")
  expect_setequal(result_polys$patch, unique(expected_polys$patch))
  expect_identical(spatialdataR::instances(polys_shape), result_polys$patch)
  for (id in result_polys$patch) {
    vertices <- sf::st_coordinates(result_polys[result_polys$patch == id, ])
    hull <- as.matrix(expected_polys[expected_polys$patch == id, c("x", "y")])
    expect_equal(unname(vertices[, c("X", "Y")]), unname(rbind(hull, hull[1, ])))
  }
})

test_that("subsetting the patch table drops the polygons of excluded patches", {
  set.seed(123)
  result <- getPatches.sd(
    sd, X = "distance_to_fibroblast", npatches = 4, n_iters = 3, verbose = FALSE
  )
  patches <- sd_table(result, "patches")
  keep <- colnames(patches)[c(3, 1)]
  result <- spatialdataR::`table<-`(result, i = "patches", value = patches[, keep])

  polys <- spatialdataR::data(spatialdataR::shapes(result)$patch_polys) |>
    sf::st_as_sf()
  expect_setequal(polys$patch, keep)
})

test_that("getPatches.sd stores patch results under `patch_shape` and `patch_table`", {
  set.seed(123)
  result <- getPatches.sd(
    sd, X = "distance_to_fibroblast", npatches = 4, n_iters = 3,
    patch_shape = "my_patches", patch_table = "my_patch_table", verbose = FALSE
  )

  expect_true("my_patches" %in% spatialdataR::shapeNames(result))
  expect_false("patch_polys" %in% spatialdataR::shapeNames(result))
  expect_identical(spatialdataR::tableNames(result), c("table", "my_patch_table"))
  expect_identical(
    spatialdataR::region(sd_table(result, "my_patch_table")), "my_patches"
  )
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

test_that("getPatches.sd stores only polygons and patch IDs when log_iters/patch_diagnostics are FALSE", {
  set.seed(123)
  expected <- getPatches.spe(
    spe, X = "distance_to_fibroblast", npatches = 4, n_iters = 3,
    log_iters = FALSE, patch_diagnostics = FALSE, verbose = FALSE
  )
  set.seed(123)
  result <- getPatches.sd(
    sd, X = "distance_to_fibroblast", npatches = 4, n_iters = 3,
    log_iters = FALSE, patch_diagnostics = FALSE, verbose = FALSE
  )
  result_table <- sd_table(result)

  expect_identical(result_table$patch, expected$patch)
  expect_null(S4Vectors::metadata(result_table)$SpaceMosaic)
  expect_identical(SingleCellExperiment::reducedDimNames(result_table),
                   SingleCellExperiment::reducedDimNames(sd_table(sd)))
  expect_identical(spatialdataR::shapeNames(result), c("shape1", "patch_polys"))

  patches <- sd_table(result, "patches")
  expect_setequal(colnames(patches), unique(stats::na.omit(expected$patch)))
  expect_identical(ncol(SummarizedExperiment::colData(patches)), 0L)
  expect_identical(SingleCellExperiment::reducedDimNames(patches), character())
  expect_null(S4Vectors::metadata(patches)$SpaceMosaic)
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

  expect_s4_class(result, "SpatialData")
  patches <- sd_table(result, "patches")
  expect_identical(rownames(patches), rownames(spe))
  expect_identical(
    SummarizedExperiment::assayNames(patches),
    paste0(c("est_", "se_", "pval_"), "distance_to_fibroblast")
  )
  expect_de_assays(result, expected$de)
  de_info <- S4Vectors::metadata(patches)$SpaceMosaic$patchDE
  expect_identical(de_info$table, "table")
  expect_identical(de_info$predictors, "distance_to_fibroblast")
  expect_identical(
    de_info$patches, colnames(expected$de$distance_to_fibroblast$ests)
  )
  expect_equal(
    SummarizedExperiment::assay(sd_table(result), "residuals"),
    SummarizedExperiment::assay(expected$spe, "residuals")
  )
})

test_that("patchDE.sd with return_residuals = FALSE leaves the cell table's assays untouched", {
  obj <- with_patch(rep(c("1", "2"), length.out = ncol(spe)))

  # method defaults to "hasty".
  expected <- patchDE.spe(
    obj$spe, predictor_cols = "distance_to_fibroblast", verbose = FALSE
  )
  result <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", verbose = FALSE
  )

  expect_de_assays(result, expected$de)
  expect_identical(
    SummarizedExperiment::assayNames(sd_table(result)),
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

  expect_de_assays(result, expected$de)
})

test_that("patchDE.sd fills patches without a fit with NA", {
  obj <- with_patch(as.character(rep(1:3, length.out = ncol(spe))))
  sce <- sd_table(obj$sd)
  sce$patch[sce$patch == "3"] <- NA
  spatialdataR::table(obj$sd, "table") <- sce

  result <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", verbose = FALSE
  )
  ests <- SummarizedExperiment::assay(
    sd_table(result, "patches"), "est_distance_to_fibroblast"
  )

  expect_identical(colnames(ests), c("1", "2", "3"))
  expect_true(all(is.na(ests[, "3"])))
  expect_false(anyNA(ests[, c("1", "2")]))
})

test_that("patchMetaAnalysis.sd matches patchMetaAnalysis.spe", {
  obj <- with_patch(as.character(rep(1:3, length.out = ncol(spe))))

  de_spe <- patchDE.spe(
    obj$spe, predictor_cols = "distance_to_fibroblast", method = "limma",
    verbose = FALSE
  )
  sd_de <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", method = "limma",
    verbose = FALSE
  )

  expected <- patchMetaAnalysis.spe(de_spe, embedding_name = "Z")
  result <- patchMetaAnalysis.sd(sd_de, embedding_name = "Z")

  expect_meta_results(result, expected)
  expect_null(
    S4Vectors::metadata(sd_table(result, "patches"))$SpaceMosaic$subgroups_summary
  )
})

test_that("patchMetaAnalysis.sd with summarize_subgroups = TRUE matches patchMetaAnalysis.spe", {
  obj <- with_patch(as.character(rep(1:3, length.out = ncol(spe))))

  de_spe <- patchDE.spe(
    obj$spe, predictor_cols = "distance_to_fibroblast", method = "limma",
    verbose = FALSE
  )
  sd_de <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", method = "limma",
    verbose = FALSE
  )

  expected <- patchMetaAnalysis.spe(
    de_spe, embedding_name = "Z", summarize_subgroups = TRUE,
    cellmeta_cols = "distance_to_fibroblast"
  )
  result <- patchMetaAnalysis.sd(
    sd_de, embedding_name = "Z", summarize_subgroups = TRUE,
    cellmeta_cols = "distance_to_fibroblast"
  )

  expect_meta_results(result, expected)
  expect_equal(
    S4Vectors::metadata(sd_table(result, "patches"))$SpaceMosaic$subgroups_summary,
    expected$subgroups_summary
  )
})

test_that("patchMetaAnalysis.sd works after patchDE.sd with resid_mse = TRUE", {
  obj <- with_patch(as.character(rep(1:3, length.out = ncol(spe))))

  de_spe <- patchDE.spe(
    obj$spe, predictor_cols = "distance_to_fibroblast", verbose = FALSE
  )
  sd_de <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", resid_mse = TRUE,
    verbose = FALSE
  )
  expect_true(
    "resid_mse" %in% SummarizedExperiment::assayNames(sd_table(sd_de, "patches"))
  )

  expected <- patchMetaAnalysis.spe(de_spe, embedding_name = "Z")
  result <- patchMetaAnalysis.sd(sd_de, embedding_name = "Z")

  expect_meta_results(result, expected)
})

test_that("moranTest.sd matches moranTest.spe", {
  obj <- with_patch(rep(c("1", "2"), length.out = ncol(spe)))

  de_spe <- patchDE.spe(
    obj$spe, predictor_cols = "distance_to_fibroblast", method = "limma",
    return_residuals = TRUE, verbose = FALSE
  )
  sd_de <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", method = "limma",
    return_residuals = TRUE, verbose = FALSE
  )

  set.seed(99)
  expected <- moranTest.spe(de_spe$spe, k = 5, n_permutations = 19)
  set.seed(99)
  result <- moranTest.sd(sd_de, k = 5, n_permutations = 19)

  expect_moran_results(result, expected)
  moran_info <- S4Vectors::metadata(sd_table(result, "patches"))$SpaceMosaic$moranTest
  expect_identical(moran_info$assay_name, "residuals")
  expect_identical(moran_info$k, 5)
  expect_identical(moran_info$alternative, "greater")
})

test_that("moranTest.sd on a gene-subsetted table matches moranTest.spe", {
  obj <- with_patch(rep(c("1", "2"), length.out = ncol(spe)))

  de_spe <- patchDE.spe(
    obj$spe, predictor_cols = "distance_to_fibroblast", method = "limma",
    return_residuals = TRUE, verbose = FALSE
  )
  sd_gene_sub <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", method = "limma",
    return_residuals = TRUE, verbose = FALSE
  )
  spatialdataR::table(sd_gene_sub, "table") <- sd_table(sd_gene_sub)[1:10, ]

  set.seed(99)
  expected <- moranTest.spe(de_spe$spe[1:10, ], k = 5, n_permutations = 19)
  set.seed(99)
  result <- moranTest.sd(sd_gene_sub, k = 5, n_permutations = 19)

  # The patch table keeps every gene of patchDE.sd(); untested ones are NA.
  expect_moran_results(result, expected)
  observed <- SummarizedExperiment::assay(sd_table(result, "patches"), "moran_observed")
  expect_identical(nrow(observed), nrow(spe))
  expect_true(all(is.na(observed[-(1:10), ])))
})

test_that("moranTest.sd gives the patch table rows when patchDE.sd has not been run", {
  obj <- with_patch(rep(c("1", "2"), length.out = ncol(spe)))
  sce <- sd_table(obj$sd)
  SummarizedExperiment::assay(sce, "residuals") <-
    as.matrix(SummarizedExperiment::assay(sce, "logcounts"))
  spatialdataR::table(obj$sd, "table") <- sce
  spe_res <- obj$spe
  SummarizedExperiment::assay(spe_res, "residuals") <-
    as.matrix(SummarizedExperiment::assay(spe_res, "logcounts"))

  set.seed(99)
  expected <- moranTest.spe(spe_res, k = 5, n_permutations = 19)
  set.seed(99)
  result <- moranTest.sd(obj$sd, k = 5, n_permutations = 19)

  expect_identical(rownames(sd_table(result, "patches")), rownames(spe))
  expect_moran_results(result, expected)
})

test_that("patchDE.sd replaces earlier DE, meta-analysis and Moran results", {
  obj <- with_patch(rep(c("1", "2"), length.out = ncol(spe)))

  set.seed(1)
  sd_all <- patchDE.sd(
    obj$sd, predictor_cols = "distance_to_fibroblast", resid_mse = TRUE,
    return_residuals = TRUE, verbose = FALSE
  ) |>
    patchMetaAnalysis.sd(embedding_name = "Z", summarize_subgroups = TRUE,
                         cellmeta_cols = "distance_to_fibroblast") |>
    moranTest.sd(k = 5, n_permutations = 19)
  before <- sd_table(sd_all, "patches")
  expect_true(all(c("resid_mse", "moran_status") %in%
                    SummarizedExperiment::assayNames(before)))

  result <- patchDE.sd(
    sd_all, predictor_cols = "total", verbose = FALSE
  )
  patches <- sd_table(result, "patches")

  expect_identical(
    SummarizedExperiment::assayNames(patches),
    paste0(c("est_", "se_", "pval_"), "total")
  )
  expect_false("W" %in% SingleCellExperiment::reducedDimNames(patches))
  expect_false(any(startsWith(names(SummarizedExperiment::colData(patches)), "moran_")))
  expect_identical(names(S4Vectors::metadata(patches)$SpaceMosaic), "patchDE")
})

test_that("the `table` argument selects which table is used and modified", {
  spe_sub <- spe[, spe$celltype == "Cancer.cells"]
  sd_two <- spatialdataR::combine(
    sd,
    spe_to_spatialdata(spe_sub, shape_name = "shape2", table_name = "table2")
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
  sd_de <- patchDE.sd(
    by_name, predictor_cols = "distance_to_fibroblast", table = "table2",
    verbose = FALSE
  )
  de_spe <- patchDE.spe(
    expected, predictor_cols = "distance_to_fibroblast", verbose = FALSE
  )
  expect_identical(
    S4Vectors::metadata(sd_table(sd_de, "patches"))$SpaceMosaic$patchDE$table,
    "table2"
  )
  expect_de_assays(sd_de, de_spe$de)
  expect_meta_results(
    patchMetaAnalysis.sd(sd_de, embedding_name = "Z"),
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
  sd_no_patch_table <- spe_to_spatialdata(with_patch(rep("1", ncol(spe)))$spe)
  expect_error(patchDE.sd(sd_no_patch_table,
                          predictor_cols = "distance_to_fibroblast"),
               "has no table 'patches'")
  expect_error(patchMetaAnalysis.sd(list(de = list())),
               "must be a SpatialData object")
  expect_error(patchMetaAnalysis.sd(with_patch(rep("1", ncol(spe)))$sd),
               "run patchDE.sd\\(\\) first")
  expect_error(moranTest.sd(sd), "not found in assay names")
})
