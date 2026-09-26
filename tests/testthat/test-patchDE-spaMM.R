.patchSpatialExample <- function() {
    set.seed(72)
    n <- 80L
    df <- data.frame(X = rnorm(n))
    xy <- matrix(runif(n * 2), ncol = 2)
    tot <- sample(500:1500, n, replace = TRUE)
    y <- cbind(gene_a = rnbinom(n, mu = tot / 1000 * exp(2 + df$X / 2), size = 3),
               zero = 0)
    list(y = y, df = df, patch = rep(c("second", "first"), n / 2),
         xy = xy, tot = tot, method = "spaMM", verbose = FALSE,
         return_diagnostics = TRUE)
}

test_that("spatial patch fitting aligns interleaved cells with direct fits", {
    skip_if_not_installed("spaMM")
    ex <- .patchSpatialExample()
    ex$patch[1] <- NA
    ex$spatial_control <- list(nu = 0.7, control = list())
    fit <- do.call(patchDE, ex)
    for (id in c("first", "second")) {
        idx <- which(ex$patch == id)
        direct <- spaMMDE(ex$y[idx, ], ex$df[idx, , drop = FALSE],
                          ex$xy[idx, ], ex$tot[idx], nu = 0.7)
        expect_equal(fit$de$X$ests[, id], direct$effect[, "X"])
        expect_equal(fit$de$X$ses[, id], direct$se[, "X"])
        expect_equal(fit$de$X$pvals[, id], direct$p[, "X"])
        expect_identical(fit$diagnostics$status[fit$diagnostics$patch == id],
                         direct$diagnostics$status)
    }
    expect_identical(colnames(fit$de$X$ests), c("first", "second"))
    expect_identical(attr(fit$de, "effect_scale"), "log")
    expect_identical(attr(fit$de, "method"), "spaMM")
    expect_equal(fit$diagnostics$nu, rep(0.7, 4))
    ex$return_diagnostics <- FALSE
    expect_equal(do.call(patchDE, ex), fit$de)
})

test_that("spatial fits agree between serial and socket workers", {
    skip_if_not_installed("spaMM")
    ex <- .patchSpatialExample()
    serial <- do.call(patchDE, ex)
    parallel <- do.call(patchDE, c(ex, list(BPPARAM = BiocParallel::SnowParam(2))))
    expect_equal(parallel$de, serial$de, tolerance = 1e-6)
    serial$diagnostics$elapsed <- parallel$diagnostics$elapsed <- NULL
    expect_equal(parallel$diagnostics, serial$diagnostics, tolerance = 1e-6)
})

test_that("character contrasts remain consistent across disjoint patch levels", {
    skip_if_not_installed("spaMM")
    ex <- .patchSpatialExample()
    ex$df$group <- ifelse(ex$patch == "first", "a", "b")
    # Each patch observes a single group level, which every backend warns about.
    fit <- suppressWarnings(do.call(patchDE, ex))
    expect_true("groupb" %in% names(fit$de))
    expect_true(all(is.na(fit$de$groupb$ests)))
    expect_true(all(vapply(fit$diagnostics$dropped_predictors,
                           function(x) "groupb" %in% x, logical(1))))
})

test_that("spatial patch options fail explicitly when incompatible", {
    ex <- .patchSpatialExample()
    for (opt in c("pearson", "resid_mse", "return_residuals")) {
        args <- ex; args[[opt]] <- TRUE
        expect_error(do.call(patchDE, args), "must be FALSE")
    }
    args <- ex; args$tot <- NULL
    expect_error(do.call(patchDE, args), "positive finite")
    args <- ex; args$xy <- ex$xy[-1, ]
    expect_error(do.call(patchDE, args), "one row per cell")
    args <- ex; args$spatial_control <- list(typo = 1)
    expect_error(do.call(patchDE, args), "only nu")
    args <- ex; args$method <- "hasty"
    expect_error(do.call(patchDE, args), "only for method")
})

test_that("spatial results feed meta-analysis while preserving missing fits", {
    skip_if_not_installed("spaMM")
    ex <- .patchSpatialExample()
    fit <- do.call(patchDE, ex)
    expect_true(all(is.finite(fit$de$X$ests["gene_a", ])))
    W <- matrix(c(0, 1), ncol = 1, dimnames = list(c("first", "second"), "context"))
    meta <- patchMetaAnalysis(fit$de, W, k = 1, min_patches = 1)
    expect_true(all(is.finite(meta$X$ests["gene_a", ])))
    expect_true(all(is.na(meta$X$ests["zero", ])))
    expect_true(all(is.na(meta$X$pvals["zero", ])))
    expect_identical(attr(meta, "effect_scale"), "log")
})
