#' Post-hoc enrichment of DE-significant patch sets
#'
#' For each DE variable in a `patchMetaAnalysis()` result, this function builds
#' direction-specific binary outcomes per gene (up: z > z_cutoff,
#' down: z < -z_cutoff), then runs logistic regression for each attribute
#' column to quantify enrichment in significant patches.
#'
#' The tested attribute is fit one-at-a-time, optionally adjusting for
#' user-specified confounders from `attrib`.
#'
#' @param meta Output of `patchMetaAnalysis()`: named list where each entry has
#'   `ests` and `ses` matrices (genes x patches).
#' @param attrib Matrix/data.frame of patch attributes (patches x variables),
#'   typically from `getPatchAttributes()`. Rownames must be patch IDs.
#' @param z_cutoff Z-stat threshold for significance. Defaults to 2.
#' @param meta_variables Optional character vector selecting which entries in
#'   `names(meta)` to analyze. If NULL (default), analyze all.
#' @param gene_subset Optional character vector of genes to analyze. If NULL
#'   (default), analyze all genes available for each DE variable.
#' @param adjust_for Optional character vector of `attrib` column names to include
#'   as adjustment covariates in every model.
#' @param min_pos Minimum number of positive patches required to attempt a model
#'   for a given gene-direction combination. Defaults to 5.
#' @param p_adjust_method Method passed to `stats::p.adjust()` for q-values.
#'   Defaults to "BH".
#' @param verbose Show progress. Default TRUE.
#' @return A list with:
#'   \item{results}{Data frame with one row per
#'     DE-variable x gene x direction x tested-attribute regression. Columns:
#'     `de_variable`, `gene`, `direction`, `attribute`, `n_patches_used`,
#'     `n_sig`, `n_nonsig`, `log_odds`, `se`, `p_value`, `q_value`, `converged`,
#'     `fit_note`.}
#'   \item{skipped}{Data frame with skipped model attempts and reason.}
#'   \item{params}{Named list with function arguments used.}
#' @export
posthocEnrichment <- function(meta, attrib,
                              z_cutoff = 2,
                              meta_variables = NULL,
                              gene_subset = NULL,
                              adjust_for = NULL,
                              min_pos = 5,
                              p_adjust_method = "BH",
                              verbose = TRUE) {
  if (!is.list(meta) || length(meta) == 0) {
    stop("meta must be a non-empty list from patchMetaAnalysis().", call. = FALSE)
  }
  if (is.null(names(meta)) || any(!nzchar(names(meta)))) {
    stop("meta must be a named list.", call. = FALSE)
  }

  attrib_df <- as.data.frame(attrib, stringsAsFactors = FALSE)
  if (is.null(rownames(attrib_df))) {
    stop("attrib must have rownames corresponding to patch IDs.", call. = FALSE)
  }

  if (is.null(meta_variables)) {
    meta_variables <- names(meta)
  } else {
    missing_meta <- setdiff(meta_variables, names(meta))
    if (length(missing_meta) > 0) {
      stop(
        "meta_variables not found in meta: ",
        paste(missing_meta, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (!is.null(gene_subset)) {
    gene_subset <- unique(as.character(gene_subset))
    gene_subset <- gene_subset[nzchar(gene_subset)]
    if (length(gene_subset) == 0) {
      stop("gene_subset was provided but contains no valid gene names.", call. = FALSE)
    }
  }

  if (is.null(adjust_for)) adjust_for <- character(0)
  missing_adjust <- setdiff(adjust_for, colnames(attrib_df))
  if (length(missing_adjust) > 0) {
    stop(
      "adjust_for columns not found in attrib: ",
      paste(missing_adjust, collapse = ", "),
      call. = FALSE
    )
  }

  tested_attributes <- setdiff(colnames(attrib_df), adjust_for)
  if (length(tested_attributes) == 0) {
    stop("No testable attrib columns remain after excluding adjust_for.", call. = FALSE)
  }

  rows <- vector("list", 0)
  skipped <- vector("list", 0)

  total_steps <- 0L
  for (v in meta_variables) {
    .checkMetaEntry(meta[[v]], v)
    ests <- meta[[v]]$ests
    if (is.null(gene_subset)) {
      n_gene_use <- nrow(ests)
    } else {
      if (is.null(rownames(ests))) {
        stop(
          "gene_subset requires rownames on meta[['", v, "']]$ests.",
          call. = FALSE
        )
      }
      n_gene_use <- sum(rownames(ests) %in% gene_subset)
    }
    total_steps <- total_steps + (2L * n_gene_use * length(tested_attributes))
  }
  use_progress <- verbose && total_steps > 0L
  if (use_progress) {
    cli::cli_progress_bar("posthocEnrichment", total = total_steps)
  } else if (verbose) {
    cli::cli_alert_info("No analyzable tests after applying current filters.")
  }

  for (varname in meta_variables) {
    ests <- as.matrix(meta[[varname]]$ests)
    ses <- as.matrix(meta[[varname]]$ses)
    zmat <- ests / ses
    zmat[!is.finite(zmat)] <- NA_real_

    if (!is.null(gene_subset)) {
      gene_hits <- intersect(gene_subset, rownames(zmat))
      if (length(gene_hits) == 0) next
      zmat <- zmat[gene_hits, , drop = FALSE]
    }

    shared_patches <- intersect(colnames(ests), rownames(attrib_df))
    if (length(shared_patches) == 0) {
      stop(
        "No shared patch IDs between meta[['", varname, "']] and attrib.",
        call. = FALSE
      )
    }

    zmat <- zmat[, shared_patches, drop = FALSE]
    W <- attrib_df[shared_patches, , drop = FALSE]

    for (g in seq_len(nrow(zmat))) {
      gene <- rownames(zmat)[g]
      if (is.null(gene) || !nzchar(gene)) gene <- as.character(g)
      z <- zmat[g, ]

      outcomes <- list(
        up = z > z_cutoff,
        down = z < -z_cutoff
      )

      for (dirn in names(outcomes)) {
        y <- outcomes[[dirn]]

        for (attr in tested_attributes) {
          dat <- data.frame(
            y = y,
            attr = W[[attr]],
            stringsAsFactors = FALSE
          )
          if (length(adjust_for) > 0) {
            dat <- cbind(dat, W[, adjust_for, drop = FALSE])
          }

          # Complete-case model frame for this specific test.
          keep <- stats::complete.cases(dat)
          dat <- dat[keep, , drop = FALSE]

          if (nrow(dat) == 0) {
            skipped[[length(skipped) + 1L]] <- data.frame(
              de_variable = varname,
              gene = gene,
              direction = dirn,
              attribute = attr,
              reason = "no_complete_cases",
              stringsAsFactors = FALSE
            )
            if (use_progress) cli::cli_progress_update()
            next
          }

          n_sig <- sum(dat$y)
          n_nonsig <- nrow(dat) - n_sig

          if (n_sig < min_pos || n_nonsig < min_pos) {
            skipped[[length(skipped) + 1L]] <- data.frame(
              de_variable = varname,
              gene = gene,
              direction = dirn,
              attribute = attr,
              reason = "min_pos_not_met",
              stringsAsFactors = FALSE
            )
            if (use_progress) cli::cli_progress_update()
            next
          }

          if (!is.numeric(dat$attr)) {
            skipped[[length(skipped) + 1L]] <- data.frame(
              de_variable = varname,
              gene = gene,
              direction = dirn,
              attribute = attr,
              reason = "non_numeric_attribute",
              stringsAsFactors = FALSE
            )
            if (use_progress) cli::cli_progress_update()
            next
          }

          if (stats::sd(dat$attr) == 0) {
            skipped[[length(skipped) + 1L]] <- data.frame(
              de_variable = varname,
              gene = gene,
              direction = dirn,
              attribute = attr,
              reason = "zero_variance_attribute",
              stringsAsFactors = FALSE
            )
            if (use_progress) cli::cli_progress_update()
            next
          }

          fit <- .fitPatchEnrichmentModel(dat)
          rows[[length(rows) + 1L]] <- data.frame(
            de_variable = varname,
            gene = gene,
            direction = dirn,
            attribute = attr,
            n_patches_used = nrow(dat),
            n_sig = n_sig,
            n_nonsig = n_nonsig,
            log_odds = fit$estimate,
            se = fit$se,
            p_value = fit$p_value,
            converged = fit$converged,
            fit_note = fit$note,
            stringsAsFactors = FALSE
          )

          if (use_progress) cli::cli_progress_update()
        }
      }
    }
  }

  if (use_progress) cli::cli_progress_done()

  if (length(rows) == 0) {
    results <- data.frame(
      de_variable = character(),
      gene = character(),
      direction = character(),
      attribute = character(),
      n_patches_used = integer(),
      n_sig = integer(),
      n_nonsig = integer(),
      log_odds = numeric(),
      se = numeric(),
      p_value = numeric(),
      q_value = numeric(),
      converged = logical(),
      fit_note = character(),
      stringsAsFactors = FALSE
    )
  } else {
    results <- do.call(rbind, rows)
    results$q_value <- NA_real_
    ok <- is.finite(results$p_value)
    results$q_value[ok] <- stats::p.adjust(results$p_value[ok], method = p_adjust_method)
  }

  if (length(skipped) == 0) {
    skipped_df <- data.frame(
      de_variable = character(),
      gene = character(),
      direction = character(),
      attribute = character(),
      reason = character(),
      stringsAsFactors = FALSE
    )
  } else {
    skipped_df <- do.call(rbind, skipped)
  }

  list(
    results = results,
    skipped = skipped_df,
    params = list(
      z_cutoff = z_cutoff,
      meta_variables = meta_variables,
      gene_subset = gene_subset,
      adjust_for = adjust_for,
      min_pos = min_pos,
      p_adjust_method = p_adjust_method
    )
  )
}


.checkMetaEntry <- function(x, name) {
  if (!is.list(x) || !all(c("ests", "ses") %in% names(x))) {
    stop(
      "meta[['", name, "']] must contain ests and ses matrices.",
      call. = FALSE
    )
  }
  ests <- x$ests
  ses <- x$ses
  if (is.null(dim(ests)) || is.null(dim(ses)) || !all(dim(ests) == dim(ses))) {
    stop(
      "meta[['", name, "']] ests/ses must have identical matrix dimensions.",
      call. = FALSE
    )
  }
  if (is.null(colnames(ests))) {
    stop("meta[['", name, "']] ests must have patch column names.", call. = FALSE)
  }
}


.fitPatchEnrichmentModel <- function(dat) {
  out <- list(
    estimate = NA_real_,
    se = NA_real_,
    p_value = NA_real_,
    converged = FALSE,
    note = "fit_failed"
  )

  fit <- tryCatch(
    stats::glm(y ~ ., data = dat, family = stats::binomial()),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    out$note <- "glm_error"
    return(out)
  }

  out$converged <- isTRUE(fit$converged)
  coefs <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(coefs) || !("attr" %in% rownames(coefs))) {
    out$note <- "coef_missing"
    return(out)
  }

  out$estimate <- coefs["attr", "Estimate"]
  out$se <- coefs["attr", "Std. Error"]
  out$p_value <- coefs["attr", "Pr(>|z|)"]
  out$note <- if (out$converged) "ok" else "not_converged"
  out
}
