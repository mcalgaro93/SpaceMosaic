# Compare the current descriptive Moran implementation with a block prototype.
# Run from the repository root: Rscript dev/benchmark_moran_blocks.R
# This experiment does not change the package implementation.
# Five repeated timings use one simulated matrix per scenario. The kernel
# comparison excludes graph/output construction; end-to-end also vectorizes
# result-table assembly. Neither comparison includes fitting the OLS model.
# Rprofmem reports cumulative R allocations, NOT peak RAM or all native memory.

required <- c("Matrix", "FNN", "cli")
stopifnot(all(vapply(required, requireNamespace, logical(1), quietly = TRUE)))
source("R/getPatches.R")
source("R/moran.R")

output_dir <- file.path("dev", "runs", "moran_blocks_benchmark")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
repetitions <- 5L
block_sizes <- c(64L, 256L, 1024L)
k <- 10L

# Same statistic and invalid-input statuses as the current descriptive path.
# Graph construction is outside this kernel, as it is for the current gene loop.
block_kernel <- function(residuals, W, block_size) {
  n <- nrow(residuals)
  s0 <- sum(W)
  observed <- rep(NA_real_, ncol(residuals))
  status <- rep("ok", ncol(residuals))
  for (start in seq.int(1L, ncol(residuals), by = block_size)) {
    idx <- seq.int(start, min(start + block_size - 1L, ncol(residuals)))
    values <- residuals[, idx, drop = FALSE]
    finite <- colSums(!is.finite(values)) == 0L
    status[idx[!finite]] <- "non-finite residuals"
    if (!any(finite)) next
    valid_idx <- idx[finite]
    z <- values[, finite, drop = FALSE]
    z <- sweep(z, 2L, colMeans(z), "-")
    denominator <- colSums(z * z)
    variable <- denominator > 0
    status[valid_idx[!variable]] <- "constant residuals"
    if (!any(variable)) next
    z <- z[, variable, drop = FALSE]
    lag <- W %*% z
    observed[valid_idx[variable]] <- as.numeric(
      (n / s0) * Matrix::colSums(z * lag) / denominator[variable])
  }
  list(observed = observed, status = status)
}

current_kernel <- function(residuals, W) {
  s0 <- sum(W)
  vapply(seq_len(ncol(residuals)), function(g) {
    .moranPermutationTest(residuals[, g], W, n_permutations = 0,
                          alternative = "greater", s0 = s0)$observed
  }, numeric(1))
}

# Single-patch equivalent of getMoranGlobal(), including graph building and output.
# Multiple patches repeat this work independently; timings are per patch.
block_test <- function(residuals, xy, block_size) {
  n <- nrow(residuals)
  W <- .buildContiguityGraph(xy, k = k)
  result <- block_kernel(residuals, W, block_size)
  ok <- result$status == "ok"
  data.frame(patch = "all", gene = colnames(residuals), n_cells = n,
             observed = result$observed,
             expected = ifelse(ok, -1 / (n - 1), NA_real_),
             p_value = NA_real_, p_adjusted = NA_real_, k = k,
             effective_k = min(k, n - 1L), n_permutations = 0L,
             alternative = "greater", status = result$status,
             stringsAsFactors = FALSE)
}

# Accumulate enough work for the timer's millisecond resolution. Calibration
# and garbage collection are outside each measured batch. Timings include any
# garbage collection triggered by the computations themselves.
batch_iterations <- function(fun) {
  elapsed <- system.time(fun())[["elapsed"]]
  as.integer(min(50L, max(1L, ceiling(0.1 / max(elapsed, 0.001)))))
}
measure <- function(fun, iterations) {
  gc()
  unname(system.time(for (i in seq_len(iterations)) fun())[["elapsed"]]) /
    iterations
}
allocated_bytes <- function(fun) {
  if (!isTRUE(capabilities("profmem"))) return(NA_real_)
  path <- tempfile()
  gc()
  Rprofmem(path)
  tryCatch(fun(), finally = Rprofmem(NULL))
  lines <- readLines(path, warn = FALSE)
  sizes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  unlink(path)
  sum(sizes, na.rm = TRUE)
}

# Invalid genes must not affect neighboring columns in a block.
set.seed(20261013)
check_xy <- matrix(runif(200), ncol = 2)
check_y <- matrix(rnorm(100 * 8), nrow = 100,
                  dimnames = list(NULL, paste0("gene_", 1:8)))
check_y[, 2] <- 0
check_y[1, 3] <- NA_real_
check_y[2, 4] <- Inf
reference <- getMoranGlobal(check_y, check_xy, verbose = FALSE)
for (size in c(1L, 3L, 64L)) {
  candidate <- block_test(check_y, check_xy, size)
  stopifnot(identical(reference$status, candidate$status),
            isTRUE(all.equal(reference$observed, candidate$observed,
                             tolerance = 1e-12)),
            isTRUE(all.equal(reference$expected, candidate$expected)))
}

grid <- expand.grid(n_cells = c(100L, 300L, 1000L),
                    n_genes = c(100L, 1000L, 5000L))
raw <- allocations <- accuracy <- list()
for (case in seq_len(nrow(grid))) {
  n <- grid$n_cells[case]
  g <- grid$n_genes[case]
  set.seed(1000 + case)
  xy <- matrix(runif(n * 2), ncol = 2)
  W <- .buildContiguityGraph(xy, k = k)
  # Synthetic regression residuals: independent and spatially structured genes.
  design <- cbind(1, xy[, 1])
  y <- matrix(rnorm(n * g), nrow = n)
  spatial <- sin(8 * xy[, 1]) * cos(8 * xy[, 2])
  y[, seq_len(ceiling(g / 2))] <- y[, seq_len(ceiling(g / 2)), drop = FALSE] +
    spatial
  residuals <- qr.resid(qr(design), y)
  colnames(residuals) <- paste0("gene_", seq_len(g))
  reference <- getMoranGlobal(residuals, xy, k = k, verbose = FALSE)
  for (size in block_sizes) {
    candidate <- block_test(residuals, xy, size)
    error <- max(abs(reference$observed - candidate$observed))
    stopifnot(error < 1e-12, identical(reference$status, candidate$status))
    accuracy[[length(accuracy) + 1L]] <- data.frame(
      n_cells = n, n_genes = g, block_size = size, max_absolute_error = error)
  }
  for (scope in c("kernel", "end_to_end")) {
    funs <- list(current = if (scope == "kernel") {
      function() current_kernel(residuals, W)
    } else {
      function() getMoranGlobal(residuals, xy, k = k, verbose = FALSE)
    })
    for (size in block_sizes) {
      funs[[paste0("block_", size)]] <- local({
        block_size <- size
        if (scope == "kernel") {
          function() block_kernel(residuals, W, block_size)
        } else {
          function() block_test(residuals, xy, block_size)
        }
      })
    }
    iterations <- vapply(funs, batch_iterations, integer(1))
    for (rep in seq_len(repetitions)) {
      # Randomized order prevents always assigning a cold/warm machine to one method.
      for (method in sample(names(funs))) {
        elapsed <- measure(funs[[method]], iterations[[method]])
        raw[[length(raw) + 1L]] <- data.frame(n_cells = n, n_genes = g,
          scope = scope, method = method, repetition = rep,
          batch_iterations = iterations[[method]], elapsed_seconds = elapsed)
      }
    }
    for (method in names(funs)) {
      allocations[[length(allocations) + 1L]] <- data.frame(
        n_cells = n, n_genes = g, scope = scope, method = method,
        allocated_mib = allocated_bytes(funs[[method]]) / 1024^2)
    }
  }
  message(sprintf("Completed %d/%d: %d cells x %d genes", case, nrow(grid), n, g))
}
raw <- do.call(rbind, raw)
allocations <- do.call(rbind, allocations)
accuracy <- do.call(rbind, accuracy)
keys <- c("n_cells", "n_genes", "scope", "method")
summary <- aggregate(raw["elapsed_seconds"], raw[keys], median)
names(summary)[names(summary) == "elapsed_seconds"] <- "median_seconds"
spread <- aggregate(raw["elapsed_seconds"], raw[keys],
                    function(x) diff(range(x)))
names(spread)[names(spread) == "elapsed_seconds"] <- "range_seconds"
summary <- merge(merge(summary, spread, by = keys), allocations, by = keys)
baseline <- summary[summary$method == "current",
                    c("n_cells", "n_genes", "scope", "median_seconds")]
names(baseline)[4] <- "current_seconds"
summary <- merge(summary, baseline, by = c("n_cells", "n_genes", "scope"))
summary$speedup <- summary$current_seconds / summary$median_seconds
write.csv(raw, file.path(output_dir, "timings_raw.csv"), row.names = FALSE)
write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
write.csv(accuracy, file.path(output_dir, "accuracy.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

png(file.path(output_dir, "benchmark.png"), width = 1500, height = 900, res = 140)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
colors <- c(current = "grey35", block_64 = "#0072B2",
            block_256 = "#D55E00", block_1024 = "#009E73")
for (scope in c("kernel", "end_to_end")) {
  for (n in c(100L, 300L, 1000L)) {
    d <- summary[summary$scope == scope & summary$n_cells == n, ]
    plot(NA, xlim = range(grid$n_genes), ylim = range(d$median_seconds),
         log = "xy", xlab = "Genes per patch", ylab = "Seconds per patch (median)",
         main = paste(if (scope == "kernel") "Statistic only" else "Graph + output", n, "cells"))
    for (method in names(colors)) {
      s <- d[d$method == method, ]
      s <- s[order(s$n_genes), ]
      lines(s$n_genes, s$median_seconds, type = "b", pch = 16, col = colors[method])
    }
    if (n == 100L) legend("topleft", legend = names(colors), col = colors,
                          lty = 1, pch = 16, bty = "n", cex = 0.7)
  }
}
invisible(dev.off())
print(summary[summary$n_genes == 5000,
              c("n_cells", "scope", "method", "median_seconds", "speedup",
                "allocated_mib")], row.names = FALSE)
cat("Maximum numerical difference:", max(accuracy$max_absolute_error), "\n")
cat("Results:", normalizePath(output_dir), "\n")
