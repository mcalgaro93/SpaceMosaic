# Benchmark patch-level parallelization in patchDE() and moranTest().
#
# Run from the repository root:
#   Rscript dev/benchmark_patch_parallelization.R
#
# patchDE covers many small patches and fewer larger patches, using both DE
# methods. moranTest uses the many-small-patches scenario.
#
# Optional environment variables:
#   SPACEMOSAIC_BENCH_WORKERS=1,2,4
#   SPACEMOSAIC_BENCH_REPEATS=3
#   SPACEMOSAIC_BENCH_SMALL_CELLS=10000
#   SPACEMOSAIC_BENCH_SMALL_GENES=500
#   SPACEMOSAIC_BENCH_SMALL_PATCHES=100
#   SPACEMOSAIC_BENCH_LARGE_CELLS=20000
#   SPACEMOSAIC_BENCH_LARGE_GENES=1000
#   SPACEMOSAIC_BENCH_LARGE_PATCHES=20
#   SPACEMOSAIC_BENCH_MORAN_GENES=40
#   SPACEMOSAIC_BENCH_PERMUTATIONS=99
#   SPACEMOSAIC_BENCH_OUTPUT=dev/runs/parallel_patch_benchmark.csv
#   SPACEMOSAIC_BENCH_PLOT=dev/runs/parallel_patch_benchmark.pdf

required_packages <- c("BiocParallel", "FNN", "ggplot2", "limma", "Matrix")
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

source("R/DE.R")
source("R/limmaDE.R")
source("R/getPatches.R")
source("R/moran.R")

read_integer <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = default)))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
  }
  value
}

worker_text <- Sys.getenv("SPACEMOSAIC_BENCH_WORKERS", unset = "1,2,4")
workers <- suppressWarnings(
  as.integer(strsplit(worker_text, ",", fixed = TRUE)[[1L]])
)
if (length(workers) == 0L || anyNA(workers) || any(workers < 1L)) {
  stop("SPACEMOSAIC_BENCH_WORKERS must contain positive integers.", call. = FALSE)
}
workers <- unique(workers)
n_repeats <- read_integer("SPACEMOSAIC_BENCH_REPEATS", 3L)
n_moran_genes <- read_integer("SPACEMOSAIC_BENCH_MORAN_GENES", 40L)
n_permutations <- read_integer("SPACEMOSAIC_BENCH_PERMUTATIONS", 99L)
output_path <- Sys.getenv(
  "SPACEMOSAIC_BENCH_OUTPUT",
  unset = file.path("dev", "runs", "parallel_patch_benchmark.csv")
)
plot_path <- Sys.getenv(
  "SPACEMOSAIC_BENCH_PLOT",
  unset = file.path("dev", "runs", "parallel_patch_benchmark.pdf")
)

scenarios <- data.frame(
  scenario = c("many_small", "fewer_large"),
  n_cells = c(
    read_integer("SPACEMOSAIC_BENCH_SMALL_CELLS", 10000L, 6L),
    read_integer("SPACEMOSAIC_BENCH_LARGE_CELLS", 20000L, 6L)
  ),
  n_genes = c(
    read_integer("SPACEMOSAIC_BENCH_SMALL_GENES", 500L),
    read_integer("SPACEMOSAIC_BENCH_LARGE_GENES", 1000L)
  ),
  n_patches = c(
    read_integer("SPACEMOSAIC_BENCH_SMALL_PATCHES", 100L, 2L),
    read_integer("SPACEMOSAIC_BENCH_LARGE_PATCHES", 20L, 2L)
  ),
  stringsAsFactors = FALSE
)
if (any(scenarios$n_patches > scenarios$n_cells %/% 3L)) {
  stop("Each patch must contain at least three cells.", call. = FALSE)
}

simulate_data <- function(n_cells, n_genes, n_patches, seed) {
  set.seed(seed)
  patch <- rep(sprintf("patch_%03d", seq_len(n_patches)), length.out = n_cells)
  patch <- sample(patch)
  exposure <- stats::rnorm(n_cells)
  batch <- factor(sample(c("A", "B"), n_cells, replace = TRUE))
  df <- data.frame(exposure = exposure, batch = batch)

  gene_effect <- stats::rnorm(n_genes, sd = 0.35)
  y <- matrix(stats::rnorm(n_cells * n_genes), nrow = n_cells)
  y <- y + tcrossprod(exposure, gene_effect)
  colnames(y) <- sprintf("gene_%04d", seq_len(n_genes))
  rownames(y) <- sprintf("cell_%06d", seq_len(n_cells))

  patch_number <- match(patch, unique(patch))
  patch_centers <- cbind(
    x = 10 * cos(2 * pi * seq_len(n_patches) / n_patches),
    y = 10 * sin(2 * pi * seq_len(n_patches) / n_patches)
  )
  xy <- patch_centers[patch_number, , drop = FALSE] +
    matrix(stats::rnorm(n_cells * 2L), ncol = 2L)

  list(y = y, df = df, patch = patch, xy = xy)
}

make_backend <- function(n_workers, seed = NULL) {
  if (n_workers == 1L) {
    BiocParallel::SerialParam(RNGseed = seed)
  } else {
    BiocParallel::SnowParam(workers = n_workers, RNGseed = seed)
  }
}

time_patch_de <- function(data, method, n_workers) {
  elapsed <- system.time({
    value <- patchDE(
      y = data$y, df = data$df, patch = data$patch,
      method = method, verbose = FALSE,
      BPPARAM = make_backend(n_workers)
    )
  })[["elapsed"]]
  list(elapsed = unname(elapsed), value = value)
}

time_moran <- function(data, residuals, n_workers) {
  elapsed <- system.time({
    value <- moranTest(
      residuals = residuals, xy = data$xy, patch = data$patch,
      n_permutations = n_permutations,
      BPPARAM = make_backend(n_workers, seed = 20260916L)
    )
  })[["elapsed"]]
  list(elapsed = unname(elapsed), value = value)
}

rows <- list()
row_index <- 0L

for (scenario_index in seq_len(nrow(scenarios))) {
  scenario <- scenarios[scenario_index, ]
  message(sprintf(
    "Preparing %s: %d cells, %d genes, %d patches",
    scenario$scenario, scenario$n_cells, scenario$n_genes, scenario$n_patches
  ))
  data <- simulate_data(
    scenario$n_cells, scenario$n_genes, scenario$n_patches,
    seed = 20260916L + scenario_index
  )

  for (method in c("hasty", "limma")) {
    reference <- NULL
    for (repeat_index in seq_len(n_repeats)) {
      for (n_workers in workers) {
        message(sprintf(
          "patchDE/%s/%s: repeat %d/%d, workers = %d",
          method, scenario$scenario, repeat_index, n_repeats, n_workers
        ))
        result <- time_patch_de(data, method, n_workers)
        if (is.null(reference)) {
          reference <- result$value
        } else if (!identical(result$value, reference)) {
          stop(
            "patchDE/", method, " returned different results with ",
            n_workers, " worker(s).", call. = FALSE
          )
        }
        row_index <- row_index + 1L
        rows[[row_index]] <- data.frame(
          scenario = scenario$scenario, operation = "patchDE", method = method,
          workers = n_workers, repetition = repeat_index,
          elapsed_seconds = result$elapsed, n_cells = scenario$n_cells,
          n_genes = scenario$n_genes, n_patches = scenario$n_patches,
          n_permutations = NA_integer_, stringsAsFactors = FALSE
        )
      }
    }
  }

  if (scenario$scenario == "many_small") {
    moran_gene_count <- min(n_moran_genes, scenario$n_genes)
    residual_fit <- patchDE(
      y = data$y[, seq_len(moran_gene_count), drop = FALSE],
      df = data$df, patch = data$patch,
      return_residuals = TRUE, verbose = FALSE
    )
    reference <- NULL
    for (repeat_index in seq_len(n_repeats)) {
      for (n_workers in workers) {
        message(sprintf(
          "moranTest/%s: repeat %d/%d, workers = %d",
          scenario$scenario, repeat_index, n_repeats, n_workers
        ))
        result <- time_moran(data, residual_fit$residuals, n_workers)
        if (is.null(reference)) {
          reference <- result$value
        } else if (!identical(result$value, reference)) {
          stop(
            "moranTest returned different results with ", n_workers,
            " worker(s).", call. = FALSE
          )
        }
        row_index <- row_index + 1L
        rows[[row_index]] <- data.frame(
          scenario = scenario$scenario, operation = "moranTest",
          method = "permutation", workers = n_workers,
          repetition = repeat_index, elapsed_seconds = result$elapsed,
          n_cells = scenario$n_cells, n_genes = moran_gene_count,
          n_patches = scenario$n_patches, n_permutations = n_permutations,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  rm(data)
  invisible(gc())
}

results <- do.call(rbind, rows)
benchmark_key <- interaction(
  results$scenario, results$operation, results$method,
  drop = TRUE, lex.order = TRUE
)
serial_median <- tapply(
  results$elapsed_seconds[results$workers == 1L],
  benchmark_key[results$workers == 1L], stats::median
)
results$speedup_vs_serial <- serial_median[as.character(benchmark_key)] /
  results$elapsed_seconds

summary_table <- aggregate(
  cbind(elapsed_seconds, speedup_vs_serial) ~
    scenario + operation + method + workers,
  data = results, FUN = stats::median
)
names(summary_table)[5:6] <- c("median_seconds", "median_speedup")

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, output_path, row.names = FALSE)
print(summary_table, row.names = FALSE)
message("Raw benchmark results written to ", output_path)

summary_table$benchmark <- ifelse(
  summary_table$operation == "patchDE",
  paste0("patchDE (", summary_table$method, ")"),
  "moranTest"
)
scenario_labels <- c(
  fewer_large = "Fewer, larger patches",
  many_small = "Many small patches"
)
plot_data <- rbind(
  transform(
    summary_table,
    metric = "Elapsed time (seconds)", value = median_seconds
  ),
  transform(
    summary_table,
    metric = "Speedup vs serial", value = median_speedup
  )
)
plot_data$metric <- factor(
  plot_data$metric,
  levels = c("Elapsed time (seconds)", "Speedup vs serial")
)
baseline <- data.frame(
  metric = factor("Speedup vs serial", levels = levels(plot_data$metric)),
  value = 1
)

performance_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = factor(workers), y = value,
    color = benchmark, group = benchmark
  )
) +
  ggplot2::geom_hline(
    data = baseline,
    ggplot2::aes(yintercept = value),
    inherit.aes = FALSE,
    linetype = "dashed", color = "grey50"
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::facet_grid(
    metric ~ scenario,
    scales = "free_y",
    labeller = ggplot2::labeller(
      scenario = ggplot2::as_labeller(scenario_labels)
    )
  ) +
  ggplot2::labs(
    title = "Patch-level parallel performance",
    subtitle = paste(
      "Median of three runs; speedup values above 1 favor parallel execution"
    ),
    x = "Workers", y = NULL, color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
grDevices::pdf(plot_path, width = 9, height = 7)
print(performance_plot)
grDevices::dev.off()
message("Benchmark plots written to ", plot_path)
