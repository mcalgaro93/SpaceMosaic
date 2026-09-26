.spatialDEExample <- function() {
    set.seed(42)
    df <- data.frame(X = rnorm(60))
    xy <- matrix(runif(120), ncol = 2)
    y <- cbind(gene_a = rnbinom(60, mu = exp(2 + 0.5 * df$X), size = 3))
    list(y = y, df = df, xy = xy, tot = rep(1000, 60))
}

test_that("spaMMDE agrees with a direct NB spatial fit", {
    skip_if_not_installed("spaMM")
    ex <- .spatialDEExample()
    result <- do.call(spaMMDE, ex)
    dat <- data.frame(response = ex$y[, 1], X = ex$df$X,
                      coord_x = ex$xy[, 1], coord_y = ex$xy[, 2],
                      log_total = log(ex$tot))
    fit <- spaMM::fitme(
        response ~ X + offset(log_total) + Matern(1 | coord_x + coord_y),
        data = dat, family = spaMM::negbin(), fixed = list(nu = 0.5), method = "ML"
    )
    expected <- unname(spaMM::fixef(fit)["X"])
    expected_se <- sqrt(stats::vcov(fit)["X", "X"])
    expect_equal(result$effect[1, 1], expected, tolerance = 1e-6)
    expect_equal(result$se[1, 1], unname(expected_se), tolerance = 1e-6)
    expect_equal(result$p[1, 1], 2 * pnorm(-abs(expected / expected_se)),
                 ignore_attr = TRUE)
    expect_identical(dimnames(result$effect), list("gene_a", "X"))
    expect_identical(result$diagnostics$status, "ok")
    expect_equal(result$diagnostics$nb_shape, unname(spaMM::get_fittedPars(fit)$NB_shape))
    expect_identical(result$effect_scale, "log")
})

test_that("sparse counts and constant columns preserve the estimable result", {
    skip_if_not_installed("spaMM")
    ex <- .spatialDEExample()
    base <- do.call(spaMMDE, ex)
    ex$y <- Matrix::Matrix(cbind(ex$y, zero = 0), sparse = TRUE)
    ex$df$constant <- 1
    ex$df$single_level <- factor(rep("a", 60))
    expect_warning(res <- do.call(spaMMDE, ex),
                   "Dropped zero-variance predictors: constant, single_level")
    expect_equal(res$effect[1, "X"], base$effect[1, "X"])
    expect_true(all(is.na(res$effect[, c("constant", "single_level")])))
    expect_identical(res$diagnostics$status, c("ok", "all_zero"))
    expect_true(all(is.na(res$p[2, ])))
    expect_equal(res$diagnostics$elapsed[2], 0)
})

test_that("non-estimable designs return explicit diagnostics", {
    ex <- .spatialDEExample()
    ex$df$duplicate <- ex$df$X
    expect_identical(do.call(spaMMDE, ex)$diagnostics$status, "rank_deficient")
    ex$df <- data.frame(X = rep(1, 60))
    expect_warning(res <- do.call(spaMMDE, ex),
                   "Dropped zero-variance predictors: X")
    expect_identical(res$diagnostics$status, "no_varying_predictors")
    expect_true(all(is.na(res$effect)))
    ex <- list(y = matrix(c(1, 2), ncol = 1), df = data.frame(X = 1:2),
               xy = cbind(1:2, 1:2), tot = c(10, 10))
    expect_identical(do.call(spaMMDE, ex)$diagnostics$status, "insufficient_cells")
})

test_that("input validation rejects invalid counts and alignment", {
    ex <- .spatialDEExample()
    check <- function(field, value, pattern) {
        args <- ex
        args[[field]] <- value
        expect_error(do.call(spaMMDE, args), pattern)
    }
    check("y", ex$y + 0.5, "integer")
    check("y", -ex$y, "non-negative")
    bad <- ex$y; bad[1, 1] <- NA_real_
    check("y", bad, "finite")
    check("xy", ex$xy[-1, ], "one row per cell")
    check("xy", matrix(1, 60, 2), "distinct locations")
    check("tot", rep(0, 60), "positive finite")
    check("tot", rep(1000, 59), "one positive")
    check("df", data.frame(X = rep(NA_real_, 60)), "without NA")
    check("df", ex$df[-1, , drop = FALSE], "one row per cell")
    expect_error(do.call(spaMMDE, c(ex, list(nu = -1))), "nu must")
    expect_error(do.call(spaMMDE, c(ex, list(control = 1))), "control must")
})

test_that("failed genes do not interrupt the remaining genes", {
    skip_if_not_installed("spaMM")
    ex <- .spatialDEExample()
    ex$y <- cbind(failed = ex$y[, 1], valid = ex$y[, 1])
    original <- .fitSpaMMGene
    calls <- 0L
    local_mocked_bindings(.fitSpaMMGene = function(...) {
        calls <<- calls + 1L
        if (calls == 1L) {
            return(list(status = "fit_error", message = "test failure", elapsed = 0))
        }
        original(...)
    })
    result <- do.call(spaMMDE, ex)
    expect_identical(result$diagnostics$status, c("fit_error", "ok"))
    expect_true(is.na(result$effect[1, 1]))
    expect_true(is.finite(result$effect[2, 1]))
    expect_identical(result$diagnostics$message[1], "test failure")
})

test_that("spaMM fitting errors and warnings are captured", {
    skip_if_not_installed("spaMM")
    ex <- .spatialDEExample()
    local_mocked_bindings(fitme = function(...) stop("optimizer failed"),
                          .package = "spaMM")
    result <- do.call(spaMMDE, ex)
    expect_identical(result$diagnostics$status, "fit_error")
    expect_match(result$diagnostics$message, "optimizer failed")
    expect_true(all(is.na(result$effect)))
})

test_that("factor contrasts and non-syntactic names are preserved", {
    skip_if_not_installed("spaMM")
    ex <- .spatialDEExample()
    ex$df <- data.frame(`cell group` = factor(rep(c("a", "b"), 30)),
                        `coord_x` = ex$df$X, check.names = FALSE)
    result <- do.call(spaMMDE, ex)
    X <- model.matrix(~ ., ex$df)
    expected_names <- colnames(X)[-1]
    explicit <- ex
    explicit$df <- as.data.frame(X[, -1, drop = FALSE])
    names(explicit$df) <- c("group_b", "exposure")
    other <- do.call(spaMMDE, explicit)
    expect_identical(colnames(result$effect), expected_names)
    expect_equal(unname(result$effect), unname(other$effect))
    expect_equal(unname(result$se), unname(other$se))
    expect_identical(result$diagnostics$status, "ok")
})

test_that("signaled and stored warnings are flagged without discarding the fit", {
    skip_if_not_installed("spaMM")
    ex <- .spatialDEExample()
    clean <- do.call(spaMMDE, ex)
    original <- spaMM::fitme
    local_mocked_bindings(fitme = function(...) {
        fit <- original(...)
        warning("test signaled warning")
        fit$warnings <- list(test = "test stored warning")
        fit
    }, .package = "spaMM")
    result <- do.call(spaMMDE, ex)
    expect_identical(result$diagnostics$status, "fit_warning")
    expect_match(result$diagnostics$message, "test signaled warning")
    expect_match(result$diagnostics$message, "test stored warning")
    # spaMM warns for reasons ranging from fatal to cosmetic, so usable
    # estimates are reported and the caller filters on status.
    expect_equal(result$effect, clean$effect)
    expect_equal(result$se, clean$se)
    expect_equal(result$p, clean$p)
    expect_equal(result$diagnostics$nb_shape, clean$diagnostics$nb_shape)
})

test_that("unusable inference outranks the warning label", {
    skip_if_not_installed("spaMM")
    ex <- .spatialDEExample()
    original_fitme <- spaMM::fitme
    original_fixef <- spaMM::fixef
    local_mocked_bindings(
        fitme = function(...) {
            fit <- original_fitme(...)
            warning("test signaled warning")
            fit
        },
        fixef = function(...) {
            coefs <- original_fixef(...)
            coefs[] <- NaN
            coefs
        },
        .package = "spaMM"
    )
    result <- do.call(spaMMDE, ex)
    expect_identical(result$diagnostics$status, "invalid_inference")
    expect_match(result$diagnostics$message, "Non-finite coefficients")
    expect_match(result$diagnostics$message, "test signaled warning")
    expect_true(all(is.na(result$effect)))
    expect_true(all(is.na(result$se)))
    expect_true(all(is.na(result$p)))
})

test_that("counts are validated through the values a matrix stores", {
    ex <- .spatialDEExample()
    sparse <- Matrix::Matrix(ex$y, sparse = TRUE)
    sparse[1, 1] <- 0.5
    ex$y <- sparse
    expect_error(do.call(spaMMDE, ex), "integer")
    # Reading the x slot keeps a sparse matrix sparse during validation.
    expect_identical(.storedValues(Matrix::Matrix(c(0, 3, 0, 5), 2, sparse = TRUE)),
                     c(3, 5))
    expect_identical(.storedValues(matrix(1:4, 2)), 1:4)
})
