#!/usr/bin/env Rscript

# Conceptual figure for the multiscale cellular-context embedding used by
# SpaceMosaic. The synthetic data are deterministic and are used only to
# illustrate the sequence of operations.

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

set.seed(24092026)

palette <- list(
    focal = "#D55E00",
    target = "#0072B2",
    context = "#9A9A9A",
    local = "#E69F00",
    broad = "#56B4E9",
    ink = "#222222",
    grid = "#D9D9D9",
    pale = "#F4F4F4"
)

theme_methods <- function(base_size = 10) {
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
            legend.position = "bottom",
            plot.margin = ggplot2::margin(8, 8, 8, 8)
        )
}

# Generate an irregular tissue-like point pattern with one focal target cell.
n_background <- 84L
theta <- stats::runif(n_background, 0, 2 * pi)
radius <- sqrt(stats::runif(n_background, 0.05, 1))
cells <- data.frame(
    cell = paste0("c", seq_len(n_background)),
    x = 5.8 * radius * cos(theta) + stats::rnorm(n_background, sd = 0.22),
    y = 4.0 * radius * sin(theta) + stats::rnorm(n_background, sd = 0.18),
    stringsAsFactors = FALSE
)
cells <- rbind(
    data.frame(cell = "focal", x = 0, y = 0, stringsAsFactors = FALSE),
    cells
)

target_probability <- stats::plogis(0.45 - 0.12 * cells$x + 0.08 * cells$y)
cells$is_target <- stats::runif(nrow(cells)) < target_probability
cells$is_target[1] <- TRUE
cells$population <- ifelse(cells$is_target, "Target cell", "Other cell")

# Synthetic PCA scores. They vary smoothly in space plus cell-level noise.
cells$PC1 <- as.numeric(scale(0.65 * cells$x + sin(cells$y / 1.8) +
    stats::rnorm(nrow(cells), sd = 0.75)))
cells$PC2 <- as.numeric(scale(-0.55 * cells$y + cos(cells$x / 1.6) +
    stats::rnorm(nrow(cells), sd = 0.70)))
cells$PC3 <- as.numeric(scale(sin(cells$x / 2.1) - cos(cells$y / 1.5) +
    stats::rnorm(nrow(cells), sd = 0.65)))
embedding <- as.matrix(cells[, c("PC1", "PC2", "PC3")])

nearest_indices <- function(index, candidates, k) {
    candidates <- setdiff(candidates, index)
    distance2 <- (cells$x[candidates] - cells$x[index])^2 +
        (cells$y[candidates] - cells$y[index])^2
    candidates[order(distance2)[seq_len(min(k, length(candidates)))]]
}

neighbor_mean <- function(index, k) {
    idx <- nearest_indices(index, seq_len(nrow(cells)), k)
    colMeans(embedding[idx, , drop = FALSE])
}

# Multiscale Z is computed for every cell before restricting to target cells.
z5 <- t(vapply(seq_len(nrow(cells)), neighbor_mean, numeric(3), k = 5L))
z50 <- t(vapply(seq_len(nrow(cells)), neighbor_mean, numeric(3), k = 50L))
colnames(z5) <- paste0("PC", 1:3, "_k5")
colnames(z50) <- paste0("PC", 1:3, "_k50")
z_all <- cbind(z5, z50)

focal_index <- 1L
neighbors_5 <- nearest_indices(focal_index, seq_len(nrow(cells)), 5L)
neighbors_50 <- nearest_indices(focal_index, seq_len(nrow(cells)), 50L)
target_indices <- which(cells$is_target)
target_neighbors_10 <- nearest_indices(focal_index, target_indices, 10L)
z_focal_smoothed <- colMeans(z_all[target_neighbors_10, , drop = FALSE])

# Panel A: the starting cell-level embedding.
p_a <- ggplot2::ggplot(cells, ggplot2::aes(x, y)) +
    ggplot2::geom_point(
        ggplot2::aes(fill = PC1, shape = population),
        size = 2.4, colour = "white", stroke = 0.35
    ) +
    ggplot2::geom_point(
        data = cells[focal_index, ], shape = 21, size = 5.0,
        fill = palette$focal, colour = palette$ink, stroke = 1.0
    ) +
    ggplot2::annotate(
        "label", x = 0.35, y = 0.55, label = "focal cell  i",
        hjust = 0, size = 3.1, linewidth = 0.2,
        colour = palette$ink, fill = "white"
    ) +
    ggplot2::scale_fill_gradient2(
        low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426",
        midpoint = 0, name = "PC1 score"
    ) +
    ggplot2::scale_shape_manual(values = c("Target cell" = 21, "Other cell" = 24)) +
    ggplot2::coord_equal() +
    ggplot2::labs(
        title = "Cell-level embedding",
        subtitle = "Each cell carries a vector of PCA scores"
    ) +
    theme_methods()

# Panel B: show the two neighborhood scales as facets.
make_scale_frame <- function(k, indices, label) {
    out <- cells
    out$scale <- label
    out$membership <- "Outside neighborhood"
    out$membership[indices] <- paste0(k, " nearest cells")
    out$membership[focal_index] <- "Focal cell"
    out
}
scale_cells <- rbind(
    make_scale_frame(5L, neighbors_5, "Local scale  k = 5"),
    make_scale_frame(50L, neighbors_50, "Broader scale  k = 50")
)
scale_cells$scale <- factor(
    scale_cells$scale,
    levels = c("Local scale  k = 5", "Broader scale  k = 50")
)

p_b <- ggplot2::ggplot(scale_cells, ggplot2::aes(x, y)) +
    ggplot2::geom_segment(
        data = subset(scale_cells, membership != "Outside neighborhood" &
            membership != "Focal cell"),
        ggplot2::aes(x = 0, y = 0, xend = x, yend = y),
        colour = "#BEBEBE", linewidth = 0.25, alpha = 0.6,
        inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
        ggplot2::aes(fill = membership), shape = 21, size = 1.9,
        colour = "white", stroke = 0.25
    ) +
    ggplot2::geom_point(
        data = subset(scale_cells, membership == "Focal cell"),
        shape = 21, size = 4.2, fill = palette$focal,
        colour = palette$ink, stroke = 0.9
    ) +
    ggplot2::facet_wrap(~scale, nrow = 1) +
    ggplot2::scale_fill_manual(values = c(
        "5 nearest cells" = palette$local,
        "50 nearest cells" = palette$broad,
        "Focal cell" = palette$focal,
        "Outside neighborhood" = "#E6E6E6"
    )) +
    ggplot2::coord_equal() +
    ggplot2::labs(
        title = "Multiscale spatial averaging",
        subtitle = "PCA scores are averaged separately at each scale"
    ) +
    theme_methods() +
    ggplot2::theme(
        strip.text = ggplot2::element_text(face = "bold", colour = palette$ink),
        legend.position = "none"
    )

# Panel C: make the two concatenated blocks, and the later smoothed vector,
# visible as a compact heat map.
heat_values <- rbind(
    "Concatenated focal Z" = c(z5[focal_index, ], z50[focal_index, ]),
    "Smoothed target Z*" = z_focal_smoothed
)
heat_df <- do.call(
    rbind,
    lapply(seq_len(nrow(heat_values)), function(i) {
        data.frame(
            operation = rownames(heat_values)[i],
            column = seq_len(ncol(heat_values)),
            value = as.numeric(heat_values[i, ]),
            stringsAsFactors = FALSE
        )
    })
)
heat_df$operation <- factor(
    heat_df$operation,
    levels = rev(rownames(heat_values))
)
heat_df$block <- rep(c("k = 5", "k = 50"), each = 3L, times = 2L)

p_c <- ggplot2::ggplot(
    heat_df,
    ggplot2::aes(column, operation, fill = value)
) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", value)),
        size = 3.0, colour = palette$ink
    ) +
    ggplot2::annotate(
        "segment", x = 3.5, xend = 3.5, y = 0.5, yend = 2.5,
        linewidth = 0.9, colour = palette$ink
    ) +
    ggplot2::annotate("text", x = 2, y = 2.62, label = "z(i, k = 5)",
        fontface = "bold", size = 3.2) +
    ggplot2::annotate("text", x = 5, y = 2.62, label = "z(i, k = 50)",
        fontface = "bold", size = 3.2) +
    ggplot2::scale_x_continuous(
        breaks = seq_len(6L), labels = rep(c("PC1", "PC2", "PC3"), 2),
        expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_fill_gradient2(
        low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426",
        midpoint = 0, guide = "none"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
        title = "Concatenation and smoothing",
        subtitle = "The two scales remain separate blocks in Z"
    ) +
    theme_methods() +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(colour = palette$ink, size = 9),
        axis.text.y = ggplot2::element_text(colour = palette$ink, size = 8.5),
        plot.margin = ggplot2::margin(16, 8, 8, 8)
    )

# Panel D: the second averaging step uses the 10 nearest target cells only.
smooth_cells <- cells
smooth_cells$status <- ifelse(
    seq_len(nrow(cells)) %in% target_neighbors_10,
    "10 nearest target cells",
    ifelse(cells$is_target, "Other target cells", "Non-target cells")
)
smooth_cells$status[focal_index] <- "Focal target cell"

p_d <- ggplot2::ggplot(smooth_cells, ggplot2::aes(x, y)) +
    ggplot2::geom_segment(
        data = smooth_cells[target_neighbors_10, ],
        ggplot2::aes(x = 0, y = 0, xend = x, yend = y),
        colour = palette$target, linewidth = 0.45, alpha = 0.55,
        inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
        ggplot2::aes(fill = status), shape = 21, size = 2.1,
        colour = "white", stroke = 0.3
    ) +
    ggplot2::geom_point(
        data = smooth_cells[focal_index, ], shape = 21, size = 4.8,
        fill = palette$focal, colour = palette$ink, stroke = 1
    ) +
    ggplot2::annotate(
        "label", x = 2.25, y = 2.85,
        label = "z*(a) = mean Z of\n10 nearest target cells",
        size = 3.0, linewidth = 0.25, fill = "white",
        colour = palette$ink
    ) +
    ggplot2::scale_fill_manual(values = c(
        "Focal target cell" = palette$focal,
        "10 nearest target cells" = palette$target,
        "Other target cells" = "#A7CBE2",
        "Non-target cells" = "#E3E3E3"
    )) +
    ggplot2::coord_equal() +
    ggplot2::labs(
        title = "Target-cell smoothing",
        subtitle = "Multiscale Z vectors are averaged over nearby targets"
    ) +
    theme_methods() +
    ggplot2::theme(legend.position = "none")

figure <- (p_a | p_b) / (p_c | p_d) +
    patchwork::plot_annotation(
        title = "Construction of the multiscale cellular-context representation",
        subtitle = paste(
            "Single-cell embeddings are averaged at two spatial scales, concatenated,",
            "and regularized across nearby target cells"
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

base_name <- file.path(figure_dir, "multiscale_embedding")
ggplot2::ggsave(
    paste0(base_name, ".svg"), figure,
    width = 13, height = 8.5, units = "in", bg = "white",
    device = svglite::svglite
)
ggplot2::ggsave(
    paste0(base_name, ".png"), figure,
    width = 13, height = 8.5, units = "in", dpi = 320, bg = "white",
    device = ragg::agg_png
)
ggplot2::ggsave(
    paste0(base_name, ".pdf"), figure,
    width = 13, height = 8.5, units = "in", bg = "white"
)

message("Wrote figure files to: ", figure_dir)
