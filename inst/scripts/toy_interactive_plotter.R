# Minimal mouse-colon example for the SpaceMosaic interactive plotter

devtools::load_all()
library(HDF5Array)
library(SpatialExperiment)

# Load the dataset used in vignettes/mouse_colon.Rmd.
h5_path <- system.file("extdata", "mouse_colon", package = "SpaceMosaic")
spe <- HDF5Array::loadHDF5SummarizedExperiment(h5_path)
xy <- spatialCoords(spe)
cell_ids <- as.character(spe$cell_id)
rownames(xy) <- cell_ids
rownames(spatialCoords(spe)) <- cell_ids

# Immune-surveillance toy question: how many immune cells surround each cell?
is_epithelial <- spe$tier1 == "Epithelial"
is_immune <- spe$tier1 == "Immune"
# Use a 50-nearest-neighbor graph to count the number of immune cells in each cell's neighborhood.
nn <- FNN::get.knn(as.matrix(xy), k = 50)$nn.index
# Count the number of immune cells in each cell's neighborhood and store it as a new column in `spe`.
spe$local_immune_cells <- rowSums(matrix(is_immune[nn], nrow = nrow(nn)))

# Encode the cellular neighborhoods and identify epithelial patches.
Z <- embedCellNeighborhoods(
  mat = reducedDim(spe, "PCA"),
  xy = xy,
  ks = c(5, 50)
)

set.seed(42)
patch_fit <- getPatches(
  xy = xy[is_epithelial, ],
  X = spe$local_immune_cells[is_epithelial],
  Z = Z[is_epithelial, ],
  npatches = 40,
  n_iters = 6,
  verbose = FALSE
)
epithelial_patches <- patch_fit$patch

# runInteractivePlotter expects assignments aligned to all cells in `spe`.
patches <- rep(NA_character_, ncol(spe))
names(patches) <- cell_ids
patches[is_epithelial] <- epithelial_patches

# Add a few real, interpretable patch-level annotations.
cell_groups <- split(which(is_epithelial), epithelial_patches)
most_common <- function(x) names(which.max(table(x)))
patch_data <- data.frame(
  patch = names(cell_groups),
  n_cells = lengths(cell_groups),
  mean_immune_cells = vapply(
    cell_groups, function(i) mean(spe$local_immune_cells[i]), numeric(1)
  ),
  dominant_subtype = vapply(
    cell_groups, function(i) most_common(spe$tier2[i]), character(1)
  ),
  neighborhood = vapply(
    cell_groups, function(i) most_common(spe$leiden_neigh[i]), character(1)
  )
)

# Test the association between expression and local immune abundance within
# each patch, then borrow information from patches with similar neighborhoods.
counts <- as.matrix(assay(spe, "counts")[, is_epithelial])
de_res <- patchDE(
  y = t(counts),
  df = data.frame(X = spe$local_immune_cells[is_epithelial]),
  patch = epithelial_patches,
  pearson = TRUE,
  tot = colSums(counts),
  verbose = FALSE
)
W <- getPatchAttributes(
  Z = Z[is_epithelial, , drop = FALSE],
  patch = epithelial_patches
)
meta_res <- patchMetaAnalysis(DEobj = de_res, W = W, k = 10)

# Gene-by-patch posterior Z-scores consumed by the interactive plotter.
metats <- meta_res$X$ests / meta_res$X$ses
# Get the top 10 genes with the largest number of patches with |Z| > 2.
top_genes <- rownames(metats)[order(rowMeans(abs(metats) > 2), decreasing = TRUE)[1:10]]

# Explore cell annotations, switch the patch fill, reorder layers and export
# the resulting plot (or its reproducible R code) directly from the app.
runInteractivePlotter(
  spe = spe,
  patches = patches,
  metats = metats[top_genes, , drop = FALSE],
  patch_data = patch_data,
  meaningful_vars = c(
    "tier1", "tier2", "leiden_neigh", "local_immune_cells"
  )
)
