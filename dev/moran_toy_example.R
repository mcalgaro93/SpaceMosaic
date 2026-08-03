# Minimal end-to-end example for residual Moran diagnostics.
#
# Run from the repository root:
#   Rscript dev/moran_toy_example.R


required_packages <- c("FNN", "Matrix")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install required packages first: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

if (!file.exists("DESCRIPTION") ||
    read.dcf("DESCRIPTION", fields = "Package")[[1L]] != "SpaceMosaic") {
  stop("Run this script from the SpaceMosaic repository root.", call. = FALSE)
}

# Source the development versions rather than an installed package.
source("R/DE.R")
source("R/getPatches.R")
source("R/moran.R")

set.seed(42)

n_per_patch <- 12L
patch <- rep(c("patch_1", "patch_2"), each = n_per_patch)
exposure <- rep(seq(-1, 1, length.out = n_per_patch), 2)

xy <- cbind(
  x = rep(seq_len(n_per_patch), 2),
  y = rep(c(0, 5), each = n_per_patch)
)

# The three neighboring blocks create nonlinear spatial structure that is not
# explained by the linear exposure term.
spatial_signal <- rep(
  c(rep(2, 4), rep(-2, 4), rep(2, 4)),
  2
)

y <- cbind(
  spatial_gene = 10 + 2 * exposure + spatial_signal,
  noise_gene = 10 - exposure + stats::rnorm(length(exposure))
)

fit <- patchDE(
  y = y,
  df = data.frame(exposure),
  patch = patch,
  return_residuals = TRUE,
  verbose = FALSE
)

set.seed(42)
moran_result <- moranTest(
  residuals = fit$residuals,
  xy = xy,
  patch = patch,
  k = 2,
  n_permutations = 999,
  alternative = "greater",
  p_adjust_method = "BH",
  adjustment_scope = "global"
)

print(
  moran_result[, c(
    "patch", "gene", "observed",
    "p_value", "p_adjusted", "status"
  )],
  row.names = FALSE
)

spatial_rows <- moran_result$gene == "spatial_gene"
noise_rows <- moran_result$gene == "noise_gene"

# These checks make the example useful as a quick manual smoke test while the
# formal regression tests remain in tests/testthat/.
stopifnot(
  all(moran_result$status == "ok"),
  all(moran_result$observed[spatial_rows] > moran_result$observed[noise_rows]),
  all(moran_result$p_adjusted[spatial_rows] < 0.05)
)

message("Moran diagnostic toy example completed successfully.")
