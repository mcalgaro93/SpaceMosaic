# Ad hoc tuning sandbox for getPatches.
#
# Run from repository root:
#   Rscript dev/getPatches_tuning_sandbox.R


runname <- "oldpass2-alpha0beta1hunger2mahal12"


npatches <- 80
n_iters <- 12
n_candidates <- 20
max_elongation <- 20 #3
alpha <- 0           #.3
beta <- 1            #.5
hunger_weight <- 2   #.3
mahal_radius <- 12    #3
x_weighted_ellipse_second_pass <- FALSE
x_ellipse_gamma <- 1
x_ellipse_wmax <- 3


required_packages <- c(
  "FNN",
  "Matrix",
  "igraph",
  "cli"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages) > 0) {
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

# Use source so local edits in R/getPatches.R are picked up immediately.
source("R/getPatches.R")
source("R/plotting.R")

mini_path <- file.path("inst", "extdata", "miniCRC.RDS")
if (!file.exists(mini_path)) {
  stop("Could not find example dataset: ", mini_path, call. = FALSE)
}
miniCRC <- readRDS(mini_path)

if (!is.list(miniCRC) || !all(c("xy", "X") %in% names(miniCRC))) {
  stop("miniCRC must be a list containing at least components 'xy' and 'X'.", call. = FALSE)
}

xy <- as.matrix(miniCRC$xy)
if (ncol(xy) < 2) {
  stop("miniCRC$xy must have at least 2 columns for spatial coordinates.", call. = FALSE)
}
xy <- xy[, 1:2, drop = FALSE]
mode(xy) <- "numeric"

X <- miniCRC$X
if (is.null(dim(X))) {
  X <- matrix(X, ncol = 1)
}
X <- as.matrix(X)
mode(X) <- "numeric"

if (nrow(xy) != nrow(X)) {
  stop("miniCRC$xy and miniCRC$X must have the same number of rows.", call. = FALSE)
}

cell_names <- rownames(xy)
if (is.null(cell_names)) {
  if (!is.null(names(miniCRC$X)) && length(miniCRC$X) == nrow(xy)) {
    cell_names <- names(miniCRC$X)
  } else {
    cell_names <- as.character(seq_len(nrow(xy)))
  }
}
rownames(xy) <- cell_names
rownames(X) <- cell_names

complete <- stats::complete.cases(xy) & stats::complete.cases(X)
xy <- xy[complete, , drop = FALSE]
X <- X[complete, , drop = FALSE]

n_cells <- nrow(xy)
message("Cells used: ", format(n_cells, big.mark = ","))

message(
  sprintf(
    paste0(
      "Running getPatches | alpha=%.2f beta=%.2f hunger=%.2f mahal=%.2f ",
      "xw2=%s gamma=%.2f wmax=%.2f"
    ),
    alpha,
    beta,
    hunger_weight,
    mahal_radius,
    as.character(x_weighted_ellipse_second_pass),
    x_ellipse_gamma,
    x_ellipse_wmax
  )
)

set.seed(0)
t0 <- proc.time()[[3L]]
result <- getPatches(
  xy = xy,
  X = X,
  npatches = npatches,
  Z = NULL,
  alpha = alpha,
  beta = beta,
  hunger_weight = hunger_weight,
  max_elongation = max_elongation,
  mahal_radius = mahal_radius,
  n_candidates = n_candidates,
  n_iters = n_iters,
  x_weighted_ellipse_second_pass = x_weighted_ellipse_second_pass,
  x_ellipse_gamma = x_ellipse_gamma,
  x_ellipse_wmax = x_ellipse_wmax,
  log_iters = TRUE,
  verbose = TRUE
)
elapsed <- proc.time()[[3L]] - t0
message("Elapsed seconds: ", round(elapsed, 2))

get_x_signal <- function(X_mat) {
  if (ncol(X_mat) == 1L) return(X_mat[, 1])
  message("X has multiple columns; using the first column for coloring.")
  X_mat[, 1]
}

plot_patch_iterations_with_x <- function(xy_mat, membership_log, x_signal,
                                         point_cex = 0.6,
                                         boundary_col = "black",
                                         boundary_lwd = 1.2,
                                         out_file = NULL) {
  stopifnot(nrow(xy_mat) == nrow(membership_log), length(x_signal) == nrow(xy_mat))

  n_it <- ncol(membership_log)
  n_col <- 1#max(1L, min(4L, n_it))
  n_row <- 1#ceiling(n_it / n_col)

  finite_vals <- x_signal[is.finite(x_signal)]
  if (length(finite_vals) == 0L) {
    point_cols <- rep("grey70", length(x_signal))
  } else {
    rng <- range(finite_vals)
    if (rng[1] == rng[2]) {
      point_cols <- rep("#2c7fb8", length(x_signal))
    } else {
      pal <- grDevices::colorRampPalette(c("#0d0887", "#f0f921"))(200)
      scaled <- (x_signal - rng[1]) / (rng[2] - rng[1])
      idx <- pmax(1L, pmin(200L, floor(scaled * 199) + 1L))
      point_cols <- pal[idx]
      point_cols[!is.finite(scaled)] <- "grey70"
    }
  }

  if (!is.null(out_file)) {
    grDevices::pdf(out_file, width = 6, height = 6)#width = 4.5 * n_col, height = 4.5 * n_row, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mfrow = c(n_row, n_col), mar = c(0,0,1,0), xaxs = "i", yaxs = "i")

  for (it in seq_len(n_it)) {
    patches <- membership_log[, it]

    graphics::plot(
      xy_mat[, 1],
      xy_mat[, 2],
      col = point_cols,
      pch = 16,
      cex = point_cex,
      xlab = "",
      ylab = "",
      main = paste0("Iteration ", it),
      asp = 1
    )

    polygons <- getPatchPolys(xy_mat, patches)
    if (nrow(polygons) > 0L) {
      for (pid in unique(polygons$patch)) {
        p <- polygons[polygons$patch == pid, , drop = FALSE]
        graphics::polygon(p$x, p$y, border = boundary_col, lwd = boundary_lwd)
      }
    }
  }
}

x_signal <- get_x_signal(X)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
plot_file <- file.path("dev", "runs", paste0("getPatches_iterations_", timestamp, "_", runname, ".pdf"))
plot_patch_iterations_with_x(
  xy_mat = xy,
  membership_log = result$membership_log,
  x_signal = x_signal,
  out_file = plot_file
)
message("Saved iteration plots to: ", plot_file)
