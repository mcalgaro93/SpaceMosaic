# Generate the mouse colon example dataset distributed in
# inst/extdata/mouse_colon.
#
# Data source
# -----------
# The complete dataset is obtained from the Bioconductor ExperimentHub package
# MerfishData with MerfishData::MouseColonIbdCadinu2024(). The source dataset
# contains MERFISH measurements from the onset and recovery of mouse colitis.
#
# Primary reference
# -----------------
# Cadinu P, Sivanathan KN, Misra A, et al. (2024). Charting the cellular
# biogeography in colitis reveals fibroblast trajectories and coordinated
# spatial remodeling. Cell 187(8):2010-2028.e30.
# DOI: 10.1016/j.cell.2024.03.013
# PMID: 38569542
#
# Processing performed here
# -------------------------
# * retain sample_id == "1" and slice_id == "1";
# * retain only the "counts" assay;
# * materialize counts as a sparse dgCMatrix to reduce the example size;
# * add the epithelial-subtype and cellular-neighborhood color palettes used by
#   SpaceMosaic examples;
# * save the resulting SpatialExperiment as an HDF5-backed
#   SummarizedExperiment.
#
# Run this script from the root of the SpaceMosaic source repository. The
# source object is large and may need to be downloaded into the local
# ExperimentHub cache. Existing output is replaced only after all validation
# checks below pass.
#
# Last validated on 2026-07-14 with R 4.6.0, Bioconductor 3.23,
# MerfishData 1.14.1, Matrix 1.7.5, HDF5Array 1.40.0,
# SummarizedExperiment 1.42.0, and S4Vectors 0.50.1. The exact environment used
# for each future build is written to sessionInfo.txt in the output directory.

required_packages <- c(
    "HDF5Array",
    "Matrix",
    "MerfishData",
    "S4Vectors",
    "SummarizedExperiment"
)

missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
    stop(
        "Install the packages required to generate the dataset: ",
        paste(missing_packages, collapse = ", "),
        call. = FALSE
    )
}

description <- "DESCRIPTION"
if (!file.exists(description) ||
        read.dcf(description, fields = "Package")[[1L]] != "SpaceMosaic") {
    stop("Run this script from the SpaceMosaic repository root.", call. = FALSE)
}

source_spe <- MerfishData::MouseColonIbdCadinu2024()

# Isolate the first sample and slice for a reasonably sized example dataset.
keep <- as.character(source_spe$sample_id) == "1" &
    as.character(source_spe$slice_id) == "1"
if (anyNA(keep) || sum(keep) != 25119L) {
    stop(
        "The expected sample 1 / slice 1 subset has changed; review the source ",
        "dataset before regenerating the example.",
        call. = FALSE
    )
}
mouse_colon <- source_spe[, keep]
rm(source_spe)

# Retain only raw counts and materialize them as a column-compressed sparse
# matrix. This intentionally drops logcounts and other assay representations.
SummarizedExperiment::assays(mouse_colon) <-
    SummarizedExperiment::assays(mouse_colon)["counts"]
dense_counts <- SummarizedExperiment::assay(mouse_colon, "counts")
sparse_counts <- methods::as(dense_counts, "dgCMatrix")
SummarizedExperiment::assay(
    mouse_colon,
    "counts",
    withDimnames = FALSE
) <- sparse_counts
rm(dense_counts, sparse_counts)

epithelial_subtype_colors <- c(
    "TA" = "#9b8fd6",
    "Stem cells" = "#889a50",
    "Goblet 1" = "#d8c66b",
    "Colonocytes" = "#61429e",
    "Repair associated  (Arg1+)" = "#e47b98",
    "EEC" = "#e4be4f",
    "Goblet 2" = "#c2d096",
    "IAE 3" = "#92478e",
    "IAE 2" = "#a7493a",
    "IAE 1" = "#d8d2f4"
)

neighbor_colors <- c(
    "ME1" = "#4daf4a",
    "SM1" = "#984ea3",
    "MU2" = "#f781bf",
    "MU1" = "#a65628",
    "ME2" = "#377eb8",
    "MU3" = "#a4b862",
    "MES1" = "#b84120",
    "MES3" = "#381981",
    "MU4" = "#6fd3e0",
    "MES2" = "#a28b40"
)

epithelial_subtypes <- unique(as.character(
    mouse_colon$tier2[mouse_colon$tier1 == "Epithelial"]
))
neighbors <- unique(as.character(mouse_colon$leiden_neigh))

if (!setequal(epithelial_subtypes, names(epithelial_subtype_colors))) {
    stop(
        "Epithelial subtypes have changed; update the color palette explicitly.",
        call. = FALSE
    )
}
if (!setequal(neighbors, names(neighbor_colors))) {
    stop(
        "Cellular neighborhoods have changed; update the color palette explicitly.",
        call. = FALSE
    )
}

object_metadata <- S4Vectors::metadata(mouse_colon)
object_metadata$epithelial_subtype_colors <- epithelial_subtype_colors
object_metadata$neighbor_colors <- neighbor_colors
object_metadata$space_mosaic_provenance <- list(
    generated_on = as.character(Sys.Date()),
    source = "MerfishData::MouseColonIbdCadinu2024()",
    source_package_version = as.character(utils::packageVersion("MerfishData")),
    reference_doi = "10.1016/j.cell.2024.03.013",
    reference_pmid = "38569542",
    filters = c(sample_id = "1", slice_id = "1"),
    retained_assays = "counts",
    counts_class = "dgCMatrix"
)
S4Vectors::metadata(mouse_colon) <- object_metadata

stopifnot(
    identical(dim(mouse_colon), c(943L, 25119L)),
    identical(SummarizedExperiment::assayNames(mouse_colon), "counts"),
    methods::is(SummarizedExperiment::assay(mouse_colon), "dgCMatrix")
)

output_dir <- file.path("inst", "extdata", "mouse_colon")
HDF5Array::saveHDF5SummarizedExperiment(
    mouse_colon,
    dir = output_dir,
    replace = TRUE
)

writeLines(
    capture.output(utils::sessionInfo()),
    file.path(output_dir, "sessionInfo.txt")
)
