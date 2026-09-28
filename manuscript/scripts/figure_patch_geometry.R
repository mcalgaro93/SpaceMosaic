#!/usr/bin/env Rscript

# Conceptual figure for the elliptical patch model and assignment score used by
# SpaceMosaic. The values are synthetic and chosen to expose the geometry.

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

set.seed(25092026)

palette <- list(
    focal = "#D55E00",
    patch = "#0072B2",
    local = "#E69F00",
    context = "#009E73",
    rejected = "#CC79A7",
    ink = "#222222",
    grey = "#777777",
    pale = "#F3F5F7"
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
            axis.title = ggplot2::element_text(colour = palette$ink),
            axis.text = ggplot2::element_text(colour = palette$ink),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(colour = "#ECECEC", linewidth = 0.3),
            legend.title = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(8, 8, 8, 8)
        )
}

# Define an elongated patch covariance.
mu <- c(0, 0)
angle <- 30 * pi / 180
rotation <- matrix(
    c(cos(angle), sin(angle), -sin(angle), cos(angle)),
    nrow = 2
)
eigenvalues <- c(4.0, 0.65)
sigma <- rotation %*% diag(eigenvalues) %*% t(rotation)
sigma_inv <- solve(sigma)
major_direction <- rotation[, 1]
minor_direction <- rotation[, 2]

# Draw current patch cells from the same covariance.
n_patch_cells <- 58L
standard_points <- matrix(stats::rnorm(n_patch_cells * 2L), ncol = 2)
patch_xy <- standard_points %*% diag(sqrt(eigenvalues)) %*% t(rotation)
patch_cells <- data.frame(x = patch_xy[, 1], y = patch_xy[, 2])

ellipse_path <- function(radius, n = 240L) {
    theta <- seq(0, 2 * pi, length.out = n)
    unit_circle <- cbind(cos(theta), sin(theta))
    xy <- radius * unit_circle %*% diag(sqrt(eigenvalues)) %*% t(rotation)
    data.frame(x = xy[, 1], y = xy[, 2], radius = factor(radius))
}
ellipses <- do.call(rbind, lapply(c(1, 2, 3), ellipse_path))

euclidean_cap <- 4.4
circle_theta <- seq(0, 2 * pi, length.out = 300L)
radius_circle <- data.frame(
    x = euclidean_cap * cos(circle_theta),
    y = euclidean_cap * sin(circle_theta)
)

# A and B have the same Euclidean distance but lie on different ellipse axes.
# C lies along the major axis but exceeds the Euclidean radius cap.
candidate_distance <- 3.1
candidates <- rbind(
    A = candidate_distance * major_direction,
    B = candidate_distance * minor_direction,
    C = 5.15 * major_direction
)
candidates <- data.frame(
    candidate = rownames(candidates),
    x = candidates[, 1],
    y = candidates[, 2],
    stringsAsFactors = FALSE
)
candidates$euclidean <- sqrt(candidates$x^2 + candidates$y^2)
candidates$mahalanobis <- sqrt(vapply(
    seq_len(nrow(candidates)),
    function(i) {
        d <- as.numeric(candidates[i, c("x", "y")])
        t(d) %*% sigma_inv %*% d
    },
    numeric(1)
))
candidates$status <- c(
    "Admissible",
    "Beyond Mahalanobis radius",
    "Beyond Euclidean radius"
)

candidate_segments <- data.frame(
    x = 0, y = 0,
    xend = candidates$x,
    yend = candidates$y,
    candidate = candidates$candidate
)

# Panel A: full spatial geometry.
p_a <- ggplot2::ggplot() +
    ggplot2::geom_path(
        data = radius_circle, ggplot2::aes(x, y),
        colour = palette$rejected, linewidth = 0.75, linetype = "22"
    ) +
    ggplot2::geom_path(
        data = ellipses,
        ggplot2::aes(x, y, group = radius, colour = radius),
        linewidth = 0.75
    ) +
    ggplot2::geom_point(
        data = patch_cells, ggplot2::aes(x, y),
        shape = 21, size = 2.1, fill = "#B9DBEC",
        colour = "white", stroke = 0.25
    ) +
    ggplot2::geom_segment(
        data = candidate_segments,
        ggplot2::aes(x, y, xend = xend, yend = yend, linetype = candidate),
        colour = palette$grey, linewidth = 0.55,
        arrow = grid::arrow(length = grid::unit(0.11, "inches"))
    ) +
    ggplot2::geom_point(
        ggplot2::aes(x = mu[1], y = mu[2]),
        shape = 4, size = 5.2, stroke = 1.3, colour = palette$ink
    ) +
    ggplot2::geom_point(
        data = candidates,
        ggplot2::aes(x, y, fill = status),
        shape = 21, size = 4.2, colour = palette$ink, stroke = 0.8
    ) +
    ggplot2::geom_label(
        data = candidates,
        ggplot2::aes(x, y, label = candidate),
        nudge_x = c(0.35, 0.35, 0.35),
        nudge_y = c(0.35, 0.35, 0.35),
        size = 3.2, fontface = "bold", linewidth = 0.2,
        fill = "white", colour = palette$ink
    ) +
    ggplot2::annotate(
        "label", x = -0.2, y = 0.45, label = "centroid ~~ mu[p]",
        hjust = 1, size = 3.0, linewidth = 0.2, fill = "white",
        parse = TRUE
    ) +
    ggplot2::annotate(
        "text", x = -4.7, y = 4.15,
        label = "dashed circle: Euclidean cap",
        hjust = 0, colour = palette$rejected, size = 3.1
    ) +
    ggplot2::scale_colour_manual(
        values = c("1" = "#9ECAE1", "2" = "#4292C6", "3" = "#08519C"),
        labels = c("Mahalanobis radius 1", "Mahalanobis radius 2", "Mahalanobis radius 3")
    ) +
    ggplot2::scale_fill_manual(values = c(
        "Admissible" = palette$context,
        "Beyond Mahalanobis radius" = palette$local,
        "Beyond Euclidean radius" = palette$rejected
    )) +
    ggplot2::scale_linetype_manual(values = c(A = "solid", B = "dashed", C = "dotdash"), guide = "none") +
    ggplot2::coord_equal(xlim = c(-5.3, 5.7), ylim = c(-4.7, 4.7)) +
    ggplot2::labs(
        title = "Elliptical spatial model",
        subtitle = "Patch covariance determines orientation, extent, and admissible assignments"
    ) +
    theme_methods() +
    ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.position = "bottom"
    )

# Panel B: demonstrate why Euclidean and Mahalanobis distances are not
# interchangeable for an elongated patch.
distance_df <- rbind(
    data.frame(
        candidate = candidates$candidate[1:2],
        distance = candidates$euclidean[1:2],
        measure = "Euclidean distance"
    ),
    data.frame(
        candidate = candidates$candidate[1:2],
        distance = candidates$mahalanobis[1:2],
        measure = "Mahalanobis distance"
    )
)
distance_df$candidate <- factor(distance_df$candidate, levels = c("A", "B"))

p_b <- ggplot2::ggplot(
    distance_df,
    ggplot2::aes(candidate, distance, fill = candidate)
) +
    ggplot2::geom_hline(
        data = data.frame(measure = "Mahalanobis distance", y = 3),
        ggplot2::aes(yintercept = y),
        linetype = "22", linewidth = 0.65, colour = palette$rejected
    ) +
    ggplot2::geom_col(width = 0.62, colour = "white") +
    ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", distance)),
        vjust = -0.45, size = 3.2, colour = palette$ink
    ) +
    ggplot2::facet_wrap(~measure, nrow = 1) +
    ggplot2::scale_fill_manual(values = c(A = palette$context, B = palette$local)) +
    ggplot2::scale_y_continuous(
        limits = c(0, 4.5), expand = ggplot2::expansion(mult = c(0, 0.04))
    ) +
    ggplot2::labs(
        title = "Same Euclidean distance, different spatial fit",
        subtitle = "Candidate B crosses the Mahalanobis cutoff because it lies across the narrow axis",
        x = "Candidate cell", y = "Distance"
    ) +
    theme_methods() +
    ggplot2::theme(
        strip.text = ggplot2::element_text(face = "bold", colour = palette$ink),
        legend.position = "none"
    )

# Panel C: signed contributions to the joint score. Values are illustrative;
# the signs and interpretation match the implemented score.
score_df <- data.frame(
    component = factor(
        c("Spatial fit", "X diversity", "Z mismatch", "Patch hunger"),
        levels = rev(c("Spatial fit", "X diversity", "Z mismatch", "Patch hunger"))
    ),
    contribution = c(-1.45, 1.10, -0.72, -0.24),
    explanation = c(
        "G[ap]",
        "+ beta %.% group('||', x[a]^'*' - bar(x)[p]^'*', '||')^2",
        "- alpha %.% group('||', z[a]^'*' - bar(z)[p]^'*', '||')^2 / sigma[Z]^2",
        "+ log(h[p])"
    ),
    stringsAsFactors = FALSE
)
score_df$direction <- ifelse(score_df$contribution >= 0, "Increases score", "Decreases score")

p_c <- ggplot2::ggplot(
    score_df,
    ggplot2::aes(contribution, component, fill = direction)
) +
    ggplot2::geom_vline(xintercept = 0, colour = palette$ink, linewidth = 0.45) +
    ggplot2::geom_col(width = 0.58, colour = "white") +
    ggplot2::geom_text(
        ggplot2::aes(
            label = explanation,
            x = ifelse(contribution >= 0, contribution + 0.08, contribution - 0.08),
            hjust = ifelse(contribution >= 0, 0, 1)
        ),
        size = 3.0, colour = palette$ink, parse = TRUE
    ) +
    ggplot2::scale_fill_manual(values = c(
        "Increases score" = palette$context,
        "Decreases score" = palette$rejected
    )) +
    ggplot2::scale_x_continuous(
        limits = c(-2.5, 2.5),
        breaks = c(-2, -1, 0, 1, 2)
    ) +
    ggplot2::labs(
        title = "Signed components of the joint assignment score",
        subtitle = "The candidate is assigned to the patch with the largest total score",
        x = expression("Contribution to " * Q[ap]), y = NULL
    ) +
    theme_methods() +
    ggplot2::theme(legend.position = "bottom")

figure <- p_a / (p_b | p_c) +
    patchwork::plot_layout(heights = c(1.28, 1)) +
    patchwork::plot_annotation(
        title = "Patch geometry and joint cell-to-patch assignment",
        subtitle = paste(
            "Elliptical spatial fit is combined with design contrast, context coherence,",
            "and patch hunger"
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
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))

base_name <- file.path(figure_dir, "patch_geometry")
ggplot2::ggsave(
    paste0(base_name, ".svg"), figure,
    width = 13, height = 10, units = "in", bg = "white",
    device = svglite::svglite
)
ggplot2::ggsave(
    paste0(base_name, ".png"), figure,
    width = 13, height = 10, units = "in", dpi = 320, bg = "white",
    device = ragg::agg_png
)
ggplot2::ggsave(
    paste0(base_name, ".pdf"), figure,
    width = 13, height = 10, units = "in", bg = "white"
)

message("Wrote figure files to: ", figure_dir)
