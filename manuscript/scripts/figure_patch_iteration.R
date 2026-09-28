#!/usr/bin/env Rscript

# Conceptual figure for patch initialization and one refinement iteration in
# SpaceMosaic. Synthetic data are used to make the sequence visually explicit.

required_packages <- c("ggplot2", "patchwork", "svglite", "ragg")
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
    stop(
        "Install the following packages before running this script: ",
        paste(missing_packages, collapse = ", "),
        call. = FALSE
    )
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
    stop("Run this file with Rscript.", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", script_arg))
manuscript_dir <- normalizePath(file.path(dirname(script_file), ".."))
figure_dir <- file.path(manuscript_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(28092026)

palette <- list(
    patch1 = "#0072B2",
    patch2 = "#E69F00",
    patch3 = "#009E73",
    patch4 = "#CC79A7",
    focal = "#D55E00",
    unassigned = "#B8B8B8",
    ink = "#222222"
)
patch_colours <- c(
    `1` = palette$patch1,
    `2` = palette$patch2,
    `3` = palette$patch3,
    `4` = palette$patch4,
    Unassigned = palette$unassigned
)

theme_methods <- function(base_size = 10.5) {
    ggplot2::theme_minimal(base_size = base_size) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold", colour = palette$ink, size = base_size + 1,
                margin = ggplot2::margin(b = 5)
            ),
            plot.subtitle = ggplot2::element_text(
                colour = "#4D4D4D", size = base_size - 0.5,
                margin = ggplot2::margin(b = 7)
            ),
            plot.tag = ggplot2::element_text(
                face = "bold", colour = palette$ink, size = base_size + 2
            ),
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            legend.title = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(8, 8, 8, 8)
        )
}

# Irregular tissue-like point cloud with a continuous design variable X.
n_cells <- 260L
theta <- stats::runif(n_cells, 0, 2 * pi)
radius <- sqrt(stats::runif(n_cells, 0.02, 1))
cells <- data.frame(
    id = seq_len(n_cells),
    x = 6.3 * radius * cos(theta) + 0.65 * sin(2 * theta) +
        stats::rnorm(n_cells, sd = 0.12),
    y = 4.1 * radius * sin(theta) + 0.35 * cos(3 * theta) +
        stats::rnorm(n_cells, sd = 0.12)
)
cells <- subset(cells, !(x < -3.7 & y > 1.4))
cells$X <- as.numeric(scale(
    0.42 * cells$x - 0.20 * cells$y +
        1.15 * sin((cells$x + cells$y) / 2.3) +
        stats::rnorm(nrow(cells), sd = 0.35)
))

ellipse_path <- function(center, covariance, level = 2, n = 180L) {
    eigen_decomp <- eigen(covariance, symmetric = TRUE)
    eigenvalues <- pmax(eigen_decomp$values, 0.06)
    angle <- seq(0, 2 * pi, length.out = n)
    circle <- cbind(cos(angle), sin(angle))
    xy <- level * circle %*% diag(sqrt(eigenvalues)) %*%
        t(eigen_decomp$vectors)
    data.frame(x = xy[, 1] + center[1], y = xy[, 2] + center[2])
}

summarize_patches <- function(labels, level = 2) {
    active <- sort(unique(labels[!is.na(labels)]))
    ellipses <- lapply(active, function(p) {
        idx <- which(labels == p)
        xy <- as.matrix(cells[idx, c("x", "y")])
        covariance <- if (nrow(xy) >= 3L) stats::cov(xy) else diag(0.2, 2)
        out <- ellipse_path(colMeans(xy), covariance, level = level)
        out$patch <- as.character(p)
        out
    })
    do.call(rbind, ellipses)
}

spatial_score <- function(labels, beta = 0.24) {
    active <- sort(unique(labels))
    score <- matrix(-Inf, nrow = nrow(cells), ncol = length(active))
    for (j in seq_along(active)) {
        p <- active[j]
        idx <- which(labels == p)
        xy_p <- as.matrix(cells[idx, c("x", "y")])
        center <- colMeans(xy_p)
        covariance <- stats::cov(xy_p) + diag(0.08, 2)
        inverse <- solve(covariance)
        delta <- sweep(as.matrix(cells[, c("x", "y")]), 2, center)
        mahalanobis2 <- rowSums((delta %*% inverse) * delta)
        spatial <- -0.5 * mahalanobis2 - 0.5 * log(det(covariance))
        x_mean <- mean(cells$X[idx])
        score[, j] <- spatial + beta * (cells$X - x_mean)^2
    }
    active[max.col(score, ties.method = "first")]
}

spatial_only_assignment <- function(labels) {
    active <- sort(unique(labels))
    score <- matrix(-Inf, nrow = nrow(cells), ncol = length(active))
    for (j in seq_along(active)) {
        p <- active[j]
        idx <- which(labels == p)
        xy_p <- as.matrix(cells[idx, c("x", "y")])
        center <- colMeans(xy_p)
        covariance <- stats::cov(xy_p) + diag(0.08, 2)
        inverse <- solve(covariance)
        delta <- sweep(as.matrix(cells[, c("x", "y")]), 2, center)
        mahalanobis2 <- rowSums((delta %*% inverse) * delta)
        score[, j] <- -0.5 * mahalanobis2 - 0.5 * log(det(covariance))
    }
    active[max.col(score, ties.method = "first")]
}

# A: inputs; B: spatial k-means initialization.
set.seed(28092026)
initial <- stats::kmeans(
    as.matrix(cells[, c("x", "y")]), centers = 4,
    nstart = 5, iter.max = 50
)$cluster

# Relabel clusters from left to right only to keep the visual stable.
cluster_centers <- aggregate(cells[, c("x", "y")], list(initial), mean)
old_levels <- cluster_centers$Group.1[order(cluster_centers$x)]
initial <- match(initial, old_levels)
initial_ellipses <- summarize_patches(initial, level = 1.75)

# C: one joint assignment and ellipse refit. Black outlines mark cells whose
# provisional label differs from the k-means label.
joint <- spatial_score(initial)
joint_ellipses <- summarize_patches(joint, level = 1.75)

# D: spatial-only labels. A small, deliberately separated island is shown to
# explain the action of the subsequent largest-component filter.
spatial_only <- spatial_only_assignment(joint)
island_patch <- spatial_only[which.max(cells$x)]
island_candidates <- order(cells$x + 0.4 * cells$y, decreasing = FALSE)[1:4]
spatial_with_island <- spatial_only
spatial_with_island[island_candidates] <- island_patch
final <- spatial_with_island
final[island_candidates] <- NA_integer_

plot_cells <- function(labels) {
    data.frame(
        cells,
        patch = factor(
            ifelse(is.na(labels), "Unassigned", as.character(labels)),
            levels = c("1", "2", "3", "4", "Unassigned")
        )
    )
}

p_a <- ggplot2::ggplot(cells, ggplot2::aes(x, y)) +
    ggplot2::geom_point(
        ggplot2::aes(fill = X), shape = 21, size = 2.4,
        colour = "white", stroke = 0.25
    ) +
    ggplot2::scale_fill_gradient2(
        low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426",
        midpoint = 0, name = "Scaled X"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
        title = "Aligned target-cell inputs",
        subtitle = "All target cells enter together; no hotspot seeds are selected"
    ) +
    theme_methods() +
    ggplot2::theme(legend.position = "bottom")

initial_df <- plot_cells(initial)
p_b <- ggplot2::ggplot(initial_df, ggplot2::aes(x, y)) +
    ggplot2::geom_point(
        ggplot2::aes(fill = patch), shape = 21, size = 2.35,
        colour = "white", stroke = 0.25
    ) +
    ggplot2::geom_path(
        data = initial_ellipses,
        ggplot2::aes(x, y, group = patch, colour = patch),
        linewidth = 0.8, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = patch_colours, drop = FALSE) +
    ggplot2::scale_colour_manual(values = patch_colours, guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(
        title = "Spatial k-means initialization",
        subtitle = "A complete starting partition defines initial centroids and ellipses"
    ) +
    theme_methods() +
    ggplot2::theme(legend.position = "none")

joint_df <- plot_cells(joint)
joint_df$changed <- joint != initial
p_c <- ggplot2::ggplot(joint_df, ggplot2::aes(x, y)) +
    ggplot2::geom_point(
        ggplot2::aes(fill = patch), shape = 21, size = 2.35,
        colour = "white", stroke = 0.25
    ) +
    ggplot2::geom_point(
        data = subset(joint_df, changed),
        shape = 21, size = 3.15, fill = NA,
        colour = palette$ink, stroke = 0.55
    ) +
    ggplot2::geom_path(
        data = joint_ellipses,
        ggplot2::aes(x, y, group = patch, colour = patch),
        linewidth = 0.8, inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = patch_colours, drop = FALSE) +
    ggplot2::scale_colour_manual(values = patch_colours, guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(
        title = "Joint assignment and ellipse refit",
        subtitle = "Spatial, X, Z, and hunger terms can move boundary cells"
    ) +
    theme_methods() +
    ggplot2::theme(legend.position = "none")

final_df <- plot_cells(final)
island_df <- plot_cells(spatial_with_island)[island_candidates, ]
p_d <- ggplot2::ggplot(final_df, ggplot2::aes(x, y)) +
    ggplot2::geom_point(
        ggplot2::aes(fill = patch), shape = 21, size = 2.35,
        colour = "white", stroke = 0.25
    ) +
    ggplot2::geom_point(
        data = island_df, shape = 4, size = 3.5,
        colour = palette$ink, stroke = 0.9
    ) +
    ggplot2::annotate(
        "label",
        x = mean(island_df$x) + 0.7,
        y = mean(island_df$y) + 0.75,
        label = "disconnected island\nset to NA",
        size = 3.0, linewidth = 0.2, fill = "white", colour = palette$ink
    ) +
    ggplot2::scale_fill_manual(values = patch_colours, drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(
        title = "Spatial-only assignment and contiguity",
        subtitle = "Only the largest connected component of each label is retained"
    ) +
    theme_methods() +
    ggplot2::theme(legend.position = "bottom")

figure <- (p_a | p_b) / (p_c | p_d) +
    patchwork::plot_annotation(
        title = "Patch initialization and one iterative refinement cycle",
        subtitle = paste(
            "Patches begin as a complete spatial partition and are repeatedly reshaped,",
            "reassigned, and filtered for contiguity"
        ),
        tag_levels = "A",
        theme = ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold", size = 16, colour = palette$ink
            ),
            plot.subtitle = ggplot2::element_text(
                size = 11, colour = "#4D4D4D",
                margin = ggplot2::margin(b = 8)
            )
        )
    ) &
    ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )

base_name <- file.path(figure_dir, "patch_iteration")
ggplot2::ggsave(
    paste0(base_name, ".svg"), figure,
    width = 13, height = 9.2, units = "in", bg = "white",
    device = svglite::svglite
)
ggplot2::ggsave(
    paste0(base_name, ".png"), figure,
    width = 13, height = 9.2, units = "in", dpi = 320, bg = "white",
    device = ragg::agg_png
)
ggplot2::ggsave(
    paste0(base_name, ".pdf"), figure,
    width = 13, height = 9.2, units = "in", bg = "white"
)

message("Wrote figure files to: ", figure_dir)
