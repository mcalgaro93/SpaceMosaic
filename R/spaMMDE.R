#' Spatial differential expression within a single patch
#'
#' Fit gene-wise negative binomial models with a log link, library-size offset
#' and an isotropic Matern spatial random effect using \pkg{spaMM}.
#'
#' @param y Non-negative integer counts, cells by genes (matrix or Matrix).
#' @param df Data frame of additive fixed-effect predictors, aligned with y.
#'   Factors use the usual R contrasts. An intercept is included automatically.
#' @param xy Numeric matrix of two spatial coordinates, aligned with y. Use
#'   consistent physical units; coordinates are not rescaled internally.
#' @param tot Positive, finite total counts per cell, calculated before gene
#'   filtering. These enter the model as offset(log(tot)).
#' @param nu Positive fixed Matern smoothness. Default 0.5 gives exponential
#'   correlation. Correlation scale and spatial variance are estimated per gene.
#' @param control List passed to [spaMM::fitme()] as its control argument.
#'
#' @return A list containing gene-by-predictor matrices `effect`, `se`, and `p`,
#'   and a `diagnostics` data frame with one row per gene. Coefficients and SEs
#'   are on the natural-log scale; p-values are unadjusted two-sided normal Wald
#'   approximations. Diagnostics include `status`, `message`, `n_cells`,
#'   `n_nonzero`, `elapsed`, `nb_shape`, `spatial_variance`, and `rho` (spaMM's
#'   inverse spatial scale, not a distance), plus `nu`. The result also records
#'   `effect_scale`, `method`, and `dropped_predictors`.
#'
#' @details This function analyzes one patch and runs sequentially. Inputs must
#'   have identical row order; row names are not used to reorder cells. It does
#'   not adjust for expression in neighboring cells or return residuals.
#'
#'   Constant design columns are omitted, reported as NA and signaled with a
#'   warning, as in [hastyDE()] and [limmaDE()]. A remaining rank-deficient
#'   design is not fitted, to avoid reporting arbitrary aliased coefficients.
#'
#'   Statuses record how far a gene got. `ok` means no warning was detected and
#'   finite coefficients with positive SEs were obtained. `fit_warning` means
#'   spaMM warned but the estimates are still finite with positive SEs; those
#'   estimates are returned, because spaMM warns for reasons that range from
#'   fatal to cosmetic and only the caller can decide which matter. Filter on
#'   `status` to discard them. Both signaled warnings and warnings stored by
#'   spaMM are recorded in `message`. `invalid_inference` (non-finite estimates
#'   or non-positive SEs), `fit_error`, `all_zero` and the design statuses
#'   return NA inference. Neither `ok` nor `fit_warning` certifies convergence
#'   to a global optimum or statistical calibration. Wald inference needs
#'   validation for small patches. No non-spatial fallback is used.
#'
#'   One negative binomial Matern model is fitted per gene, and every fit
#'   re-estimates the spatial scale, the spatial variance and the
#'   overdispersion from scratch. Cost grows linearly in the number of genes and
#'   steeply in the number of cells per patch, so restricting `y` to genes of
#'   interest beforehand is a requirement in practice, not an optimization.
#'
#' @importFrom methods slotNames
#' @examples
#' if (requireNamespace("spaMM", quietly = TRUE)) {
#'     set.seed(1)
#'     xy <- matrix(runif(80), ncol = 2)
#'     df <- data.frame(exposure = rnorm(40))
#'     y <- cbind(gene_a = rnbinom(40, mu = exp(2 + df$exposure), size = 3))
#'     fit <- spaMMDE(y, df, xy, tot = rep(1000, 40))
#'     fit$effect
#'     fit$diagnostics
#' }
#' @export

spaMMDE <- function(y, df, xy, tot, nu = 0.5, control = list()) {
    # Check inputs
    if (!(is.matrix(y) || inherits(y, "Matrix")) ||
        nrow(y) < 2L || ncol(y) < 1L) {
        stop("y must be a numeric cells-by-genes count matrix with at least two cells.")
    }
    # Validate the stored values instead of y itself: testing a sparse matrix
    # element-wise would expand it to a dense cells-by-genes copy only to check
    # that it holds counts.
    stored_counts <- .storedValues(y)
    if (!is.numeric(stored_counts) && !is.logical(stored_counts)) {
        stop("y must be a numeric cells-by-genes count matrix with at least two cells.")
    }
    if (any(!is.finite(stored_counts)) || any(stored_counts < 0) ||
        any(stored_counts != floor(stored_counts))) {
        stop("y must contain finite, non-negative integer counts.")
    }
    n <- nrow(y)
    if (!is.data.frame(df) || nrow(df) != n || ncol(df) < 1L ||
        is.null(names(df)) || anyNA(names(df)) || any(!nzchar(names(df))) ||
        anyDuplicated(names(df))) {
        stop("df must have one row per cell and uniquely named predictors.")
    }
    if (anyNA(df) || !all(vapply(df, function(z) {
        (is.numeric(z) && all(is.finite(z))) || is.factor(z) ||
            is.character(z) || is.logical(z)
    }, logical(1)))) {
        stop("df predictors must be finite numeric, logical, character or factor values without NA.")
    }
    if (!is.matrix(xy) || !is.numeric(xy) ||
        !identical(dim(xy), c(n, 2L)) || any(!is.finite(xy))) {
        stop("xy must be a finite numeric matrix with two columns and one row per cell.")
    }
    if (nrow(unique(xy)) < 2L) stop("xy must contain at least two distinct locations.")
    if (!is.numeric(tot) || !is.null(dim(tot)) || length(tot) != n ||
        any(!is.finite(tot)) || any(tot <= 0)) {
        stop("tot must contain one positive finite total per cell.")
    }
    if (!is.numeric(nu) || length(nu) != 1L || !is.finite(nu) || nu <= 0) {
        stop("nu must be a positive finite number.")
    }
    if (!is.list(control)) stop("control must be a list.")
    # Check gene names
    genes <- colnames(y)
    if (is.null(genes)) genes <- paste0("gene_", seq_len(ncol(y)))
    if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
        stop("y must have unique, non-empty gene names when supplied.")
    }
    # model.matrix cannot contrast a factor with only one defined level.
    for (j in seq_along(df)) {
        if (is.character(df[[j]])) df[[j]] <- factor(df[[j]])
        if (is.factor(df[[j]]) && nlevels(df[[j]]) < 2L) df[[j]] <- rep(0, n)
    }
    # Construct model matrix and identify varying predictors
    design <- stats::model.matrix(~ ., data = df)
    if (any(!is.finite(design))) stop("The model matrix must be finite.")
    terms <- colnames(design)[-1L]
    varying <- vapply(seq_along(terms), function(j) {
        length(unique(design[, j + 1L])) > 1L
    }, logical(1))
    X <- design[, c(TRUE, varying), drop = FALSE]
    if (any(!varying)) {
        warning("Dropped zero-variance predictors: ",
                paste(terms[!varying], collapse = ", "))
    }
    # Initialize outputs
    effect <- matrix(NA_real_, ncol(y), length(terms),
                     dimnames = list(genes, terms))
    se <- p <- effect
    # Initialize diagnostics data frame
    diagnostics <- data.frame(
        gene = genes, n_cells = n, n_nonzero = as.integer(Matrix::colSums(y > 0)),
        status = "not_fitted", message = "", elapsed = 0,
        nb_shape = NA_real_, spatial_variance = NA_real_, rho = NA_real_, nu = nu,
        stringsAsFactors = FALSE
    )
    # Check for design problems before fitting
    design_problem <- if (!any(varying)) {
        "no_varying_predictors"
    } else if (qr(X)$rank < ncol(X)) {
        "rank_deficient"
    } else if (n <= ncol(X)) {
        "insufficient_cells"
    } else NULL
    # Fit each gene
    if (!is.null(design_problem)) {
        diagnostics$status <- design_problem
        diagnostics$message <- "The requested fixed effects cannot be estimated in this patch."
    } else {
        if (!requireNamespace("spaMM", quietly = TRUE)) {
            stop("Install 'spaMM' to use spaMMDE().")
        }
        # Safe internal names avoid collisions with user predictors and preserve
        # the original model.matrix contrasts, including non-syntactic names.
        predictor_names <- paste0("predictor", seq_len(ncol(X) - 1L))
        dat <- as.data.frame(X[, -1L, drop = FALSE])
        names(dat) <- predictor_names
        dat$coord_x <- xy[, 1L]
        dat$coord_y <- xy[, 2L]
        dat$log_total <- log(tot)
        formula <- stats::reformulate(
            c(predictor_names, "offset(log_total)", "Matern(1 | coord_x + coord_y)"),
            response = "response"
        )
        # Fit each gene in a tryCatch to isolate errors and warnings
        for (g in seq_len(ncol(y))) {
            counts <- as.numeric(y[, g])
            if (all(counts == 0)) {
                diagnostics$status[g] <- "all_zero"
                diagnostics$message[g] <- "Gene has zero counts in every cell; fitting was skipped."
                next
            }
            dat$response <- counts
            result <- .fitSpaMMGene(formula, dat, nu, control, predictor_names)
            diagnostics$status[g] <- result$status
            diagnostics$message[g] <- result$message
            diagnostics$elapsed[g] <- result$elapsed
            if (result$status %in% c("ok", "fit_warning")) {
                effect[g, varying] <- result$effect
                se[g, varying] <- result$se
                p[g, varying] <- 2 * stats::pnorm(-abs(result$effect / result$se))
                diagnostics$nb_shape[g] <- result$nb_shape
                diagnostics$spatial_variance[g] <- result$spatial_variance
                diagnostics$rho[g] <- result$rho
            }
        }
    }
    list(effect = effect, se = se, p = p, diagnostics = diagnostics,
         effect_scale = "log", method = "NB-Matern-ML",
         dropped_predictors = terms[!varying])
}

# Return the values a matrix actually stores, without densifying a sparse one.
.storedValues <- function(y) {
    if (isS4(y) && "x" %in% methods::slotNames(y)) return(y@x)
    as.vector(y)
}

# Isolate errors and warnings for each gene; never abort the remaining genes.
.fitSpaMMGene <- function(formula, dat, nu, control, predictor_names) {
    warnings <- character()
    started <- proc.time()[["elapsed"]]
    result <- tryCatch(withCallingHandlers({
        fit <- spaMM::fitme(formula, data = dat, family = spaMM::negbin(),
                            fixed = list(nu = nu), method = "ML", control = control)
        stored <- unlist(fit$warnings, use.names = FALSE)
        warnings <- c(warnings, as.character(stored))
        beta <- spaMM::fixef(fit)[predictor_names]
        V <- stats::vcov(fit)[predictor_names, predictor_names, drop = FALSE]
        ses <- sqrt(diag(V))
        pars <- spaMM::get_fittedPars(fit)
        scalar <- function(x) if (length(x) == 1L) as.numeric(x) else NA_real_
        list(status = "ok", message = "", effect = beta, se = ses,
             nb_shape = scalar(pars$NB_shape),
             spatial_variance = scalar(pars$lambda),
             rho = scalar(pars$corrPars[[1L]]$rho))
    }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) list(status = "fit_error", message = conditionMessage(e)))
    if (result$status == "ok") {
        # Usability is checked before the warning label, never after: a
        # fit_warning keeps its estimates, so the label may only be applied once
        # the estimates are known to be finite with positive standard errors.
        if (any(!is.finite(result$effect)) || any(!is.finite(result$se)) ||
            any(result$se <= 0)) {
            result$status <- "invalid_inference"
            result$message <- "Non-finite coefficients or non-positive/non-finite standard errors."
        } else if (length(warnings)) {
            result$status <- "fit_warning"
        }
    }
    result$message <- paste(unique(c(result$message[nzchar(result$message)], warnings)), collapse = "; ")
    result$elapsed <- proc.time()[["elapsed"]] - started
    result
}
