# Minimal parameter sweep for Scenario 2 (Stem Cell Niche).
#
# Run from the repository root:
#   Rscript dev/tune_stem_cell_niche_patches.R
#
# For a quick smoke test without changing the grid:
#   SPACEMOSAIC_MAX_RUNS=20 Rscript dev/tune_stem_cell_niche_patches.R
# For a reproducible random sample of the full canonical grid:
#   SPACEMOSAIC_SAMPLE_RUNS=200 Rscript dev/tune_stem_cell_niche_patches.R
# Control parallel workers (default: up to 4):
#   SPACEMOSAIC_WORKERS=6 Rscript dev/tune_stem_cell_niche_patches.R

CONFIG <- list(
    seed = 123L,
    n_iters = 10L,
    top_n_plots = 4L,
    near_optimal_fraction = 0.10,
    min_cells = c("Stem cells" = 10L, "Colonocytes" = 10L, "TA" = 10L),
    grid = expand.grid(
        init_method = c("kmeans", "gradient_ellipse"),
        n_candidates = 20,
        npatches = 80,
        max_elongation = c(2, 4, 8),
        mahal_radius = c(2, 3, 4),
        beta = 1,
        hunger_weight = c(0.3, 0.5, 0.7),
        init_gradient_elongation = 4,
        init_gradient_k = 30,
        x_weighted_ellipse_second_pass = c(FALSE, TRUE),
        x_ellipse_gamma = 1,
        x_ellipse_wmax = 3,
        stringsAsFactors = FALSE
    )
)

required <- c("FNN", "HDF5Array", "SpatialExperiment", "rpart")
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) {
    stop("Install required packages first: ", paste(missing, collapse = ", "))
}
if (!file.exists("DESCRIPTION") || !file.exists("R/getPatches.R")) {
    stop("Run this script from the SpaceMosaic repository root.")
}
source("R/getPatches.R")

# Read defaults directly from getPatches() so report annotations cannot drift
# from the implementation. init_method uses the first match.arg choice.
method_formals <- formals(getPatches)
default_names <- c(
    "alpha", "beta", "hunger_weight", "max_elongation", "max_radius",
    "mahal_radius", "n_candidates", "n_iters", "init_method",
    "init_gradient_k", "init_gradient_elongation",
    "x_weighted_ellipse_second_pass", "x_ellipse_gamma", "x_ellipse_wmax"
)
METHOD_DEFAULTS <- setNames(lapply(default_names, function(name) {
    value <- eval(method_formals[[name]], envir = baseenv())
    if (name == "init_method") value[1L] else value
}), default_names)

data_path <- file.path("inst", "extdata", "mouse_colon")
spe <- HDF5Array::loadHDF5SummarizedExperiment(data_path)
xy_all <- as.matrix(SpatialExperiment::spatialCoords(spe))[, 1:2, drop = FALSE]
tier1 <- as.character(spe$tier1)
tier2_all <- as.character(spe$tier2)
use <- tier1 == "Epithelial" & stats::complete.cases(xy_all) & !is.na(tier2_all)
xy <- xy_all[use, , drop = FALSE]
tier2 <- tier2_all[use]
if (is.null(rownames(xy))) rownames(xy) <- as.character(which(use))

is_stem_all <- tier2_all == "Stem cells" & stats::complete.cases(xy_all)
if (!any(is_stem_all)) stop("No tier2 == 'Stem cells' cells found.")
X <- log1p(FNN::get.knnx(
    data = xy_all[is_stem_all, , drop = FALSE], query = xy, k = 1L
)$nn.dist[, 1L])

targets <- CONFIG$min_cells
if (!all(names(targets) %in% unique(tier2))) {
    stop("Missing tier2 labels: ",
         paste(setdiff(names(targets), unique(tier2)), collapse = ", "))
}

# The HDF5-backed object is no longer needed. Dropping it before forking avoids
# retaining its handles and reduces copy-on-write memory pressure.
rm(spe, xy_all, tier1, tier2_all, use, is_stem_all)
invisible(gc())

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path("dev", "runs", paste0("stem_cell_niche_", timestamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# gradient-specific parameters do not affect kmeans initialization. Canonicalize
# them before unique() so equivalent kmeans configurations run only once.
grid <- CONFIG$grid
is_kmeans <- grid$init_method == "kmeans"
grid$init_gradient_k[is_kmeans] <- NA
grid$init_gradient_elongation[is_kmeans] <- NA
grid <- unique(grid)
sample_runs <- suppressWarnings(as.integer(
    Sys.getenv("SPACEMOSAIC_SAMPLE_RUNS", "")
))
if (is.finite(sample_runs) && sample_runs > 0L && sample_runs < nrow(grid)) {
    set.seed(CONFIG$seed)
    grid <- grid[sample.int(nrow(grid), sample_runs), , drop = FALSE]
}
max_runs <- suppressWarnings(as.integer(Sys.getenv("SPACEMOSAIC_MAX_RUNS", "")))
if (is.finite(max_runs) && max_runs > 0L) grid <- utils::head(grid, max_runs)
grid$run_id <- sprintf("run_%03d", seq_len(nrow(grid)))

detected_cores <- parallel::detectCores(logical = FALSE)
if (!is.finite(detected_cores)) detected_cores <- 2L
default_workers <- max(1L, min(4L, detected_cores - 1L))
workers <- suppressWarnings(as.integer(Sys.getenv(
    "SPACEMOSAIC_WORKERS", as.character(default_workers)
)))
if (!is.finite(workers) || workers < 1L) workers <- 1L
workers <- min(workers, nrow(grid))
if (.Platform$OS.type == "windows" && workers > 1L) {
    warning("Fork parallelism is unavailable on Windows; using one worker.")
    workers <- 1L
}

get_param <- function(row, name, default) {
    if (!name %in% names(row) || is.na(row[[name]][1L])) default
    else row[[name]][1L]
}

summarize_patches <- function(patch, run_id) {
    ids <- sort(unique(patch[!is.na(patch)]))
    rows <- lapply(ids, function(id) {
        idx <- which(patch == id)
        counts <- table(factor(tier2[idx], levels = names(targets)))
        fractions <- pmin(as.numeric(counts) / as.numeric(targets), 1)
        data.frame(
            run_id = run_id, patch = as.character(id), size = length(idx),
            stem_cells = unname(counts["Stem cells"]),
            colonocytes = unname(counts["Colonocytes"]),
            ta_cells = unname(counts["TA"]),
            balance = min(fractions),
            good = all(as.numeric(counts) >= as.numeric(targets)),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

run_configuration <- function(i) {
    p <- grid[i, , drop = FALSE]
    Sys.setenv(
        OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1"
    )
    set.seed(CONFIG$seed)
    started <- proc.time()[[3L]]
    fit <- try(getPatches(
        xy = xy, X = X, Z = NULL, npatches = p$npatches,
        alpha = 0, beta = p$beta, hunger_weight = p$hunger_weight,
        max_elongation = p$max_elongation, mahal_radius = p$mahal_radius,
        n_candidates = p$n_candidates, n_iters = CONFIG$n_iters,
        init_method = p$init_method,
        init_gradient_k = get_param(p, "init_gradient_k", 30L),
        init_gradient_elongation = get_param(p, "init_gradient_elongation", 4),
        x_weighted_ellipse_second_pass = p$x_weighted_ellipse_second_pass,
        x_ellipse_gamma = p$x_ellipse_gamma,
        x_ellipse_wmax = p$x_ellipse_wmax,
        log_iters = FALSE, verbose = FALSE
    ), silent = TRUE)
    elapsed <- proc.time()[[3L]] - started
    if (inherits(fit, "try-error")) {
        return(list(
            result = data.frame(
                run_id = p$run_id, good_fraction = NA_real_,
                mean_balance = NA_real_,
                n_observed_patches = NA_integer_,
                elapsed_seconds = elapsed, error = as.character(fit)
            ),
            detail = NULL
        ))
    }
    metrics <- summarize_patches(fit, p$run_id)
    good_fraction <- mean(metrics$good)
    mean_balance <- mean(metrics$balance)
    list(
        result = data.frame(
            run_id = p$run_id, good_fraction = good_fraction,
            mean_balance = mean_balance, n_observed_patches = nrow(metrics),
            elapsed_seconds = elapsed, error = NA_character_
        ),
        detail = list(patch = fit, metrics = metrics)
    )
}

message(sprintf(
    "Running %d configurations with %d worker%s...",
    nrow(grid), workers, if (workers == 1L) "" else "s"
))
sweep_started <- proc.time()[[3L]]
if (workers > 1L) {
    run_outputs <- parallel::mclapply(
        seq_len(nrow(grid)), run_configuration,
        mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE
    )
} else {
    run_outputs <- lapply(seq_len(nrow(grid)), run_configuration)
}
sweep_elapsed <- proc.time()[[3L]] - sweep_started
results <- lapply(run_outputs, `[[`, "result")
details <- lapply(run_outputs, `[[`, "detail")
message(sprintf("Parameter sweep completed in %.1f seconds.", sweep_elapsed))

summary_df <- merge(grid, do.call(rbind, results), by = "run_id", all.x = TRUE,
                    sort = FALSE)
summary_df <- summary_df[
    order(-summary_df$good_fraction, -summary_df$mean_balance, na.last = TRUE),
]
valid <- summary_df[
    is.finite(summary_df$good_fraction) & is.finite(summary_df$mean_balance),
    , drop = FALSE
]
if (!nrow(valid)) stop("All parameter combinations failed; inspect run_summary.csv.")

good_cutoff <- unname(stats::quantile(
    valid$good_fraction, probs = 1 - CONFIG$near_optimal_fraction,
    type = 1, na.rm = TRUE
))
balance_cutoff <- unname(stats::quantile(
    valid$mean_balance, probs = 1 - CONFIG$near_optimal_fraction,
    type = 1, na.rm = TRUE
))
valid$near_optimal <- valid$good_fraction >= good_cutoff &
    valid$mean_balance >= balance_cutoff
optimal <- valid[1L, , drop = FALSE]
summary_df$near_optimal <- valid$near_optimal[
    match(summary_df$run_id, valid$run_id)
]

parameter_columns <- c(
    "init_method", "n_candidates", "npatches", "max_elongation",
    "mahal_radius", "beta", "hunger_weight", "init_gradient_elongation",
    "init_gradient_k", "x_weighted_ellipse_second_pass",
    "x_ellipse_gamma", "x_ellipse_wmax"
)
varying_parameters <- parameter_columns[
    vapply(valid[parameter_columns], function(x) length(unique(x[!is.na(x)])) > 1L,
           logical(1L))
]

# One interpretable classification tree: which parameter rules identify runs
# in the configured top fraction for both biological metrics?
tree_data <- valid[, c(varying_parameters, "near_optimal"), drop = FALSE]
tree_data$near_optimal <- factor(
    ifelse(tree_data$near_optimal, "near-optimal", "other"),
    levels = c("other", "near-optimal")
)
tree_control <- rpart::rpart.control(
    cp = 0.005, minbucket = max(3L, floor(nrow(tree_data) / 50L)), xval = 10L
)
optimal_tree <- if (nlevels(droplevels(tree_data$near_optimal)) > 1L) {
    rpart::rpart(
        near_optimal ~ ., data = tree_data, method = "class",
        control = tree_control
    )
} else {
    NULL
}

importance <- if (is.null(optimal_tree$variable.importance)) {
    data.frame(parameter = character(), importance = numeric())
} else {
    x <- optimal_tree$variable.importance
    data.frame(parameter = names(x), importance = unname(x) / sum(x))
}
importance <- importance[order(importance$importance, decreasing = TRUE), ]

near_ranges <- do.call(rbind, lapply(parameter_columns, function(parameter) {
    x <- valid[[parameter]][valid$near_optimal]
    x <- x[!is.na(x)]
    data.frame(
        parameter = parameter,
        values = paste(sort(unique(x)), collapse = " | "),
        stringsAsFactors = FALSE
    )
}))

all_patch_metrics <- do.call(rbind, lapply(details, function(x) {
    if (!is.null(x)) x$metrics
}))
utils::write.csv(summary_df, file.path(out_dir, "run_summary.csv"), row.names = FALSE)
utils::write.csv(all_patch_metrics, file.path(out_dir, "patch_metrics.csv"),
                 row.names = FALSE)
utils::write.csv(optimal, file.path(out_dir, "optimal_parameter_set.csv"),
                 row.names = FALSE)
utils::write.csv(near_ranges, file.path(out_dir, "near_optimal_ranges.csv"),
                 row.names = FALSE)
utils::write.csv(importance, file.path(out_dir, "parameter_importance.csv"),
                 row.names = FALSE)
saveRDS(list(
    config = CONFIG, summary = summary_df, details = details, xy = xy, X = X,
    tier2 = tier2, optimal = optimal, near_optimal_ranges = near_ranges,
    optimal_tree = optimal_tree
), file.path(out_dir, "results.rds"))

format_value <- function(x) {
    if (is.null(x)) return("NULL")
    if (is.logical(x)) return(ifelse(x, "TRUE", "FALSE"))
    if (is.numeric(x)) return(format(x, trim = TRUE, scientific = FALSE))
    as.character(x)
}

parameter_annotation <- function(row) {
    actual <- list(
        npatches = row$npatches, alpha = 0, beta = row$beta,
        hunger_weight = row$hunger_weight,
        max_elongation = row$max_elongation, max_radius = NULL,
        mahal_radius = row$mahal_radius, n_candidates = row$n_candidates,
        n_iters = CONFIG$n_iters, init_method = row$init_method,
        init_gradient_k = if (row$init_method == "kmeans") NA else row$init_gradient_k,
        init_gradient_elongation = if (row$init_method == "kmeans") NA
            else row$init_gradient_elongation,
        x_weighted_ellipse_second_pass = row$x_weighted_ellipse_second_pass,
        x_ellipse_gamma = if (!row$x_weighted_ellipse_second_pass) NA
            else row$x_ellipse_gamma,
        x_ellipse_wmax = if (!row$x_weighted_ellipse_second_pass) NA
            else row$x_ellipse_wmax
    )
    default_text <- character()
    changed_text <- character()
    ignored_text <- character()
    for (name in names(actual)) {
        value <- actual[[name]]
        if (length(value) == 1L && is.na(value)) {
            default_suffix <- if (name %in% names(METHOD_DEFAULTS)) {
                paste0(" [default=", format_value(METHOD_DEFAULTS[[name]]), "]")
            } else {
                ""
            }
            ignored_text <- c(
                ignored_text, paste0(name, "=n/a", default_suffix)
            )
        } else if (name %in% names(METHOD_DEFAULTS) &&
                   identical(value, METHOD_DEFAULTS[[name]])) {
            default_text <- c(
                default_text,
                paste0(name, "=", format_value(value), " [default]")
            )
        } else {
            suffix <- if (name %in% names(METHOD_DEFAULTS)) {
                paste0(" [default=", format_value(METHOD_DEFAULTS[[name]]), "]")
            } else {
                " [required]"
            }
            changed_text <- c(
                changed_text,
                paste0(name, "=", format_value(value), suffix)
            )
        }
    }
    list(
        defaults = paste(default_text, collapse = ", "),
        changed = paste(changed_text, collapse = ", "),
        ignored = paste(ignored_text, collapse = ", ")
    )
}

add_parameter_caption <- function(row) {
    annotation <- parameter_annotation(row)
    default_label <- paste(
        strwrap(paste("DEFAULT:", annotation$defaults), width = 145),
        collapse = "\n"
    )
    changed_label <- paste(
        strwrap(paste("CHANGED:", annotation$changed), width = 145),
        collapse = "\n"
    )
    mtext(default_label, side = 1, line = 4.3, outer = TRUE,
          cex = 0.55, col = "grey45")
    mtext(changed_label, side = 1, line = 2.5, outer = TRUE,
          cex = 0.55, col = "#D55E00")
    if (nzchar(annotation$ignored)) {
        ignored_label <- paste(
            strwrap(paste("IGNORED:", annotation$ignored), width = 145),
            collapse = "\n"
        )
        mtext(ignored_label, side = 1, line = 0.8, outer = TRUE,
              cex = 0.55, col = "grey60")
    }
}

draw_interpretable_tree <- function(tree) {
    frame <- tree$frame
    node_ids <- as.integer(row.names(frame))
    is_leaf <- frame$var == "<leaf>"
    split_questions <- setNames(rep("", nrow(frame)), node_ids)
    split_cursor <- 1L
    xlevels <- attr(tree, "xlevels")

    for (i in seq_len(nrow(frame))) {
        if (is_leaf[i]) next
        variable <- frame$var[i]
        split <- tree$splits[split_cursor, ]
        if (split[["ncat"]] < 0) {
            if (variable == "x_weighted_ellipse_second_pass") {
                question <- "X-weighted second pass = TRUE?"
            } else {
                readable_name <- unname(c(
                    max_elongation = "Maximum elongation",
                    mahal_radius = "Mahalanobis radius",
                    hunger_weight = "Hunger weight"
                )[variable])
                if (is.na(readable_name)) readable_name <- variable
                question <- sprintf("%s >= %s?", readable_name,
                                    format(split[["index"]], trim = TRUE))
            }
        } else {
            codes <- tree$csplit[split[["index"]], ]
            right_levels <- xlevels[[variable]][codes == 3L]
            readable_variable <- if (variable == "init_method") {
                "Initialization"
            } else {
                variable
            }
            readable_levels <- gsub("_", " ", right_levels, fixed = TRUE)
            if (length(right_levels) == 1L) {
                question <- paste0(
                    readable_variable, " = ", readable_levels, "?"
                )
            } else {
                question <- paste0(
                    readable_variable, " in {",
                    paste(readable_levels, collapse = ", "), "}?"
                )
            }
        }
        split_questions[as.character(node_ids[i])] <- question
        split_cursor <- split_cursor + 1L + frame$ncompete[i] +
            frame$nsurrogate[i]
    }

    depths <- floor(log(node_ids, base = 2))
    leaf_ids <- node_ids[is_leaf]
    leaf_spacing <- 1.3
    x_position <- setNames(seq_along(leaf_ids) * leaf_spacing, leaf_ids)
    position_node <- function(node) {
        key <- as.character(node)
        if (!is.na(x_position[key])) return(unname(x_position[key]))
        left <- position_node(node * 2L)
        right <- position_node(node * 2L + 1L)
        x_position[key] <<- mean(c(left, right))
        unname(x_position[key])
    }
    position_node(1L)
    x <- unname(x_position[as.character(node_ids)])
    y <- max(depths) - depths

    graphics::plot.new()
    graphics::plot.window(
        xlim = c(0.5, length(leaf_ids) * leaf_spacing + 0.8),
        ylim = c(-0.65, max(y) + 0.65)
    )

    for (i in which(!is_leaf)) {
        parent <- node_ids[i]
        for (is_true in c(FALSE, TRUE)) {
            child <- parent * 2L + as.integer(is_true)
            child_i <- match(child, node_ids)
            if (is.na(child_i)) next
            graphics::arrows(
                x[i], y[i] - 0.22, x[child_i], y[child_i] + 0.22,
                length = 0.07, lwd = 1.1, col = "grey40"
            )
            graphics::text(
                0.42 * x[i] + 0.58 * x[child_i],
                0.42 * y[i] + 0.58 * y[child_i],
                if (is_true) "TRUE" else "FALSE",
                cex = 0.68,
                col = if (is_true) "#0072B2" else "grey35"
            )
        }
    }

    for (i in seq_len(nrow(frame))) {
        near_count <- frame$yval2[i, 3L]
        total <- frame$n[i]
        pct <- 100 * near_count / total
        predicted_near <- frame$yval[i] == 2L
        heading <- if (is_leaf[i]) {
            if (predicted_near) "NEAR-OPTIMAL" else "OTHER"
        } else {
            split_questions[as.character(node_ids[i])]
        }
        label <- paste0(
            heading, "\nNear-optimal: ", near_count, "/", total,
            " (", sprintf("%.0f", pct), "%)"
        )
        fill <- if (predicted_near) "#F6C89F" else "#E6E6E6"
        graphics::rect(
            x[i] - 0.58, y[i] - 0.22, x[i] + 0.58, y[i] + 0.22,
            col = fill, border = if (predicted_near) "#D55E00" else "grey45",
            lwd = 1.2
        )
        graphics::text(x[i], y[i], label, cex = 0.63)
    }
}

# Plot the global optimum and the best configuration for each initialization
# method. Fill remaining slots with the best distinct near-optimal runs.
best_by_init <- do.call(rbind, lapply(split(valid, valid$init_method), function(x) {
    x[1L, , drop = FALSE]
}))
plot_ids <- unique(c(
    optimal$run_id, best_by_init$run_id,
    head(valid$run_id[valid$near_optimal], CONFIG$top_n_plots)
))
plot_ids <- head(plot_ids, CONFIG$top_n_plots)

pdf(file.path(out_dir, "sweet_spot.pdf"), width = 11, height = 8.5, onefile = TRUE)
palette <- hcl.colors(101, "Inferno")

# Two-metric overview.
par(mar = c(6, 6, 6, 2), oma = c(6, 0, 0, 0))
balance_pct <- 100 * valid$mean_balance
good_pct <- 100 * valid$good_fraction
balance_cutoff_pct <- 100 * balance_cutoff
good_cutoff_pct <- 100 * good_cutoff
plot(
     balance_pct, good_pct, type = "n",
     xlab = paste0(
         "Average target balance across patches (%)\n",
         "100% means that all three cell-count targets are met"
     ),
     ylab = "Patches meeting all three cell-count targets (%)",
     main = "Biological patch quality across getPatches configurations"
)
usr <- par("usr")
rect(
    balance_cutoff_pct, good_cutoff_pct, usr[2L], usr[4L],
    col = grDevices::adjustcolor("#D55E00", alpha.f = 0.10),
    border = NA
)
abline(v = balance_cutoff_pct, h = good_cutoff_pct,
       lty = 2, col = "#D55E00")
points(
    balance_pct, good_pct,
    col = ifelse(valid$near_optimal, "#D55E00", "#80808070"),
    pch = ifelse(valid$init_method == "kmeans", 17, 19)
)
points(100 * optimal$mean_balance, 100 * optimal$good_fraction,
       pch = 21, cex = 2.2, lwd = 2)
text(100 * optimal$mean_balance, 100 * optimal$good_fraction,
     paste("optimal", optimal$run_id), pos = 2)
legend("bottomright",
       legend = c("kmeans", "gradient ellipse", "near-optimal region"),
       pch = c(17, 19, 15), col = c("grey30", "grey30", "#D55E00"),
       pt.cex = c(1, 1, 1.5), bty = "n")
mtext(sprintf(
    "Near-optimal configurations rank in the top %.0f%% for both complementary quality criteria",
    100 * CONFIG$near_optimal_fraction
), side = 3, line = 1.2, cex = 0.85)
mtext(sprintf(
    "Per-patch targets: Stem cells >= %d, Colonocytes >= %d, TA >= %d",
    targets[["Stem cells"]], targets[["Colonocytes"]], targets[["TA"]]
), side = 3, line = 0.1, cex = 0.75, col = "grey35")
add_parameter_caption(optimal)

# Decision rules for near-optimal parameter regions.
if (!is.null(optimal_tree)) {
    if (nrow(optimal_tree$frame) > 1L) {
        par(mar = c(2, 2, 5, 2), oma = c(6, 0, 0, 0), xpd = TRUE)
        draw_interpretable_tree(optimal_tree)
        title("Decision tree: parameter rules identifying near-optimal runs")
    } else {
        par(mar = c(2, 2, 5, 2), oma = c(6, 0, 0, 0))
        plot.new()
        title("Decision tree: no reliable split found")
        text(
            0.5, 0.55,
            paste(
                "The tested subset does not contain enough evidence for a",
                "stable parameter rule.\nRun a larger portion of the grid or",
                "increase near_optimal_fraction."
            ),
            cex = 1.05
        )
    }
    add_parameter_caption(optimal)
}

# Compact explicit statement of the best configuration and acceptable ranges.
par(mar = c(2, 2, 4, 2), oma = c(0, 0, 0, 0))
plot.new()
title("Recommended parameter set")
annotation <- parameter_annotation(optimal)
text(0.02, 0.88, paste("Best observed run:", optimal$run_id),
     adj = 0, cex = 1.25)
text(0.02, 0.79, sprintf(
    "good patches %.1f%% | mean balance %.3f",
    100 * optimal$good_fraction, optimal$mean_balance
), adj = 0)
text(0.02, 0.66, "Values equal to getPatches defaults", adj = 0,
     col = "grey45")
text(0.04, 0.59, paste(strwrap(annotation$defaults, width = 125),
                       collapse = "\n"),
     adj = c(0, 1), cex = 0.75, col = "grey45")
text(0.02, 0.46, "Values changed from getPatches defaults", adj = 0,
     col = "#D55E00")
text(0.04, 0.39, paste(strwrap(annotation$changed, width = 125),
                       collapse = "\n"),
     adj = c(0, 1), cex = 0.75, col = "#D55E00")
if (nzchar(annotation$ignored)) {
    text(0.02, 0.26, "Parameters ignored by the selected configuration", adj = 0,
         col = "grey60")
    text(0.04, 0.19, annotation$ignored, adj = 0, cex = 0.75, col = "grey60")
}

for (id in plot_ids) {
    row <- valid[match(id, valid$run_id), , drop = FALSE]
    original_i <- match(id, grid$run_id)
    patch <- details[[original_i]]$patch
    metrics <- details[[original_i]]$metrics
    quality <- setNames(metrics$balance, metrics$patch)
    cols <- rep("grey90", length(patch))
    ok <- !is.na(patch)
    cols[ok] <- palette[1L + round(100 * quality[as.character(patch[ok])])]
    par(mar = c(1, 1, 4, 1), oma = c(7, 0, 0, 0))
    plot(xy, asp = 1, pch = 16, cex = 0.35, col = cols, axes = FALSE,
         xlab = "", ylab = "",
         main = sprintf(
             "%s | %s | good %.1f%% | balance %.3f",
             id, row$init_method, 100 * row$good_fraction,
             row$mean_balance
         ))
    for (pid in unique(patch[ok])) {
        idx <- which(patch == pid)
        if (length(idx) >= 3L) {
            hull <- grDevices::chull(xy[idx, , drop = FALSE])
            polygon(xy[idx[hull], , drop = FALSE],
                    border = "#20202080", lwd = 0.5)
        }
    }
    legend(
        "topright",
        title = "Patch target balance",
        legend = c(
            "1.00 (all targets met)", "0.75", "0.50", "0.25", "0.00",
            "Unassigned cell"
        ),
        fill = c(
            palette[c(101L, 76L, 51L, 26L, 1L)],
            "grey90"
        ),
        border = NA,
        bg = grDevices::adjustcolor("white", alpha.f = 0.88),
        box.col = "grey70",
        cex = 0.72
    )
    add_parameter_caption(row)
}
dev.off()

message("Optimal parameter set:")
print(optimal[, c(parameter_columns, "good_fraction", "mean_balance")],
      row.names = FALSE)
message("Equivalent kmeans combinations removed: gradient initialization ",
        "parameters are ignored when init_method = 'kmeans'.")
message("Results written to: ", out_dir)
