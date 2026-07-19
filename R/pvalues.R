#' Compute robust *p*-values for one or two `lavaan` objects
#'
#' Compute one or several *p*-values for a fitted `lavaan` model. Available
#' methods include penalized eigenvalue block averaging, penalized regression,
#' and familiar robust corrections. The defaults select the recommended method
#' for the fitted model.
#'
#' The `tests` argument is a character vector. Each element has one of the forms
#' `TEST`, `TEST_UG`, `TEST_ML`, `TEST_RLS`, `TEST_UG_ML`, or `TEST_UG_RLS`.
#' For example, `SB_UG_RLS`.
#'
#' * The first part names the test. The available names are listed below.
#' * Add `UG` to use the unbiased estimator of the fourth-order moment matrix
#'   from Du and Bentler (2022). Leave it out to use the standard biased matrix.
#'   Performance depends on the data-generating process, so neither gamma choice
#'   dominates in every setting.
#' * The final part chooses the chi-square statistic. `ML` uses the
#'   normal-discrepancy statistic (Bollen, 1989). `RLS` uses Browne's (1974)
#'   reweighted least squares statistic.
#'
#' The `peba` method is the recommended default. It partitions the eigenvalues
#' into `j` equally sized sets (if not possible, the smallest set is incomplete),
#' shrinks them towards their common mean, and averages within each set. Provide
#' a positive integer `j` no greater than the test df. Values from `2` through
#' `6` are usually good candidates. The method was introduced by Foldnes, Moss,
#' and Grønneberg
#' (2025).
#'
#' `pols` is a penalized regression method with a finite positive penalization
#' term. Foldnes, Moss, and Grønneberg (2025) studied `pols=2`, which has
#' good performance in a variety of contexts.
#'
#' `pall` penalizes all eigenvalues in ugamma, while `all` uses all eigenvalues
#' without penalization. `pall` is the recommended option for nested models, for
#' which the penalized methods were extended and evaluated by Foldnes,
#' Grønneberg, and Moss (2026).
#'
#' The `eba` method is the unpenalized predecessor of `peba` (Foldnes and
#' Grønneberg, 2018). It averages within the eigenvalue blocks without
#' shrinkage. `peba` usually performs better. `eba` remains available for
#' comparisons, and values from `j=2` through `j=4` tend to work best.
#'
#' Familiar corrections are available too:
#'
#' * `std` uses the ordinary chi-square reference distribution.
#' * `sb` gives the Satorra--Bentler *p*-value (Satorra and Bentler, 1994).
#' * `ss` gives the scaled and shifted *p*-value (Asparouhov and Muthén, 2010).
#' * `sf` gives the scaled *F* *p*-value (Wu and Lin, 2016).
#'
#' ## Estimators and data types
#'
#' [semTests-support] (`?semTests-support`) gives the complete matrix of
#' supported estimators, data types, and configurations. Here is the short
#' version.
#'
#' The limiting null law of the test statistic is a weighted sum of
#' chi-squares for any minimum-discrepancy estimator. `pvalues()` supports
#' ML/MLM/MLR, GLS, ULS,
#' FIML (missing data), and categorical DWLS/ULS families, with single- and
#' multi-group continuous, ordered, and mixed-indicator fits.
#' `pvalues_nested()` supports continuous estimators and categorical
#' Satorra-2000 comparisons with a delta restriction map. The fitted model must
#' expose the asymptotic moment covariance. Robust choices such as
#' `test = "satorra.bentler"` or `estimator = "MLM"`, `"MLR"`, or `"DWLS"`
#' provide it. The RLS statistic (`browne.residual.nt.model`) and the unbiased
#' (`UG`) Du-Bentler gamma are defined for the classical continuous,
#' complete-data, random-x ML case. Other families use the standard statistic
#' and biased gamma. Fixed or conditional observed exogenous predictors are
#' currently refused. Full WLS/ADF is refused outright: its weight is already
#' the inverse moment covariance, so the correction is exactly the identity and
#' every test would equal the ordinary chi-square.
#'
#' GLS, ULS, categorical DWLS/ULS, FIML missing data, and nested FIML or
#' categorical comparison are validated to numerical tolerance against an
#' independent implementation and are marked stable; see the Stability section
#' in [semTests-support].
#'
#' Both entry points require a converged fit and warn if lavaan reports an
#' inadmissible solution. Nested comparisons require the same requested
#' estimator, fitting conventions, variables, groups, sample sizes, raw data,
#' and missingness mask. `semTests` checks comparability. The substantive
#' nesting argument remains part of the analysis.
#'
#' For FIML, `fiml.convention = "observed"` (the default) uses observed
#' saturated and model information throughout. An independent implementation
#' validates this construction. The `"lavaan"` option reproduces lavaan
#' 0.7-2's robust-test construction. For a single model it uses lavaan's
#' inspected `UGamma`. For nested models it uses lavaan's H1 weight convention
#' and selected model information. The returned object records the convention,
#' which keeps this inferential choice visible in saved results.
#'
#' For a single FIML model, `"lavaan"` means the eigenvalue spectrum returned by
#' `lavInspect(fit, "UGamma")`. The scalar Yuan--Bentler--Mplus test stored by a
#' default MLR fit is a different correction and may give a different result.
#' For nested FIML models, `"lavaan"` reproduces
#' `lavTestLRT(..., method = "satorra.2000")`.
#'
#' @param object,m0,m1 One or two `lavaan` objects. `pvalues` does
#'   goodness-of-fit testing on one object, while `pvalues_nested` does
#'   hypothesis testing on two nested models.
#' @param tests A non-empty character vector of tests. Each element uses one of
#'   `TEST`, `TEST_UG`, `TEST_ML`, `TEST_RLS`, `TEST_UG_ML`, or
#'   `TEST_UG_RLS`. The defaults are the recommended options. `EBA` and `pEBA`
#'   take a positive integer number of blocks (for example, `"pEBA4"`), no
#'   greater than the test df. `pOLS` takes a finite positive penalty (for
#'   example, `"pOLS2"`).
#' @param method Nested reduction method. Only Satorra's `"2000"` construction
#'   (the paper-recommended default) is available. The Satorra--Bentler `"2001"`
#'   construction has been withdrawn because it performs poorly, and requesting
#'   it points you to `"2000"`.
#' @param A.method For nested FIML or categorical models, choose `"exact"` for
#'   the literal parameter restriction map or `"delta"` for the local
#'   moment-tangent restriction map. `"delta"` is the default and applies to a
#'   wider range of equivalent model parameterizations. Categorical nested
#'   models support `"delta"` only.
#' @param fiml.convention For FIML fits, use the fully observed-information
#'   convention (`"observed"`, the default) or reproduce lavaan 0.7-2's
#'   inspected robust-test spectrum (`"lavaan"`).
#' @name pvalues
#' @export
#' @examples
#' library("semTests")
#' library("lavaan")
#' model <- "visual  =~ x1 + x2 + x3
#'           textual =~ x4 + x5 + x6
#'           speed   =~ x7 + x8 + x9"
#' object <- cfa(model, HolzingerSwineford1939, estimator = "MLM")
#' pvalues(object)
#'
#' # For the pEBA6 method with biased gamma and ML chisq statistic:
#' pvalues(object, "pEBA6_ML")
#'
#' # Nested model comparison (constrain the textual loadings to be equal):
#' constrained <- "visual  =~ x1 + x2 + x3
#'                 textual =~ a*x4 + a*x5 + a*x6
#'                 speed   =~ x7 + x8 + x9"
#' m1 <- cfa(model, HolzingerSwineford1939, estimator = "MLM")
#' m0 <- cfa(constrained, HolzingerSwineford1939, estimator = "MLM")
#' pvalues_nested(m0, m1)
#'
#' @return A named numeric vector of p-values, of class `semTests_pvalues`,
#'   carrying an `"semtests"` attribute that records the requested tests,
#'   requested and base estimator, base-statistic and gamma choices,
#'   information type, data type, parameterization where applicable, nested
#'   method and restriction map where applicable, and degrees of freedom.
#'
#' @seealso [semTests-support] for the full list of supported configurations.
#'
#' @references
#'
#' Foldnes, N., Moss, J., & Grønneberg, S. (2025). Improved goodness of fit
#' procedures for structural equation models. *Structural Equation Modeling: A
#' Multidisciplinary Journal*, 32(1), 1--13.
#' \doi{10.1080/10705511.2024.2372028}
#'
#' Foldnes, N., Grønneberg, S., & Moss, J. (2026). Penalized eigenvalue block
#' averaging: Extension to nested model comparison and Monte Carlo evaluations.
#' *Behavior Research Methods*, 58, article 107.
#' \doi{10.3758/s13428-026-02968-4}
#'
#' Satorra, A. (2000). Scaled and adjusted restricted tests in multi-sample
#' analysis of moment structures. In R. D. H. Heijmans, D. S. G. Pollock, &
#' A. Satorra (Eds.), *Innovations in Multivariate Statistical Analysis*
#' (pp. 233--247). Kluwer Academic.
#' \doi{10.1007/978-1-4615-4603-0_17}
#'
#' Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test
#' statistic for moment structure analysis. *Psychometrika*, 66(4), 507--514.
#' \doi{10.1007/BF02296192}
#'
#' Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and
#' standard errors in covariance structure analysis. In A. von Eye &
#' C. C. Clogg (Eds.), *Latent Variables Analysis: Applications for
#' Developmental Research* (pp. 399--419). Sage.
#'
#' Asparouhov, T., & Muthén, B. O. (2010). *Simple second order chi-square
#' correction*. Mplus Technical Appendix.
#' https://www.statmodel.com/download/WLSMV_new_chi21.pdf
#'
#' Wu, H., & Lin, J. (2016). A Scaled F Distribution as an Approximation to the
#' Distribution of Test Statistics in Covariance Structure Analysis.
#' *Structural Equation Modeling*.
#' \doi{10.1080/10705511.2015.1057733}
#'
#' Foldnes, N., & Grønneberg, S. (2018). Approximating Test Statistics Using
#' Eigenvalue Block Averaging. *Structural Equation Modeling*, 25(1), 101--114.
#' \doi{10.1080/10705511.2017.1373021}
#'
#' Du, H., & Bentler, P. M. (2022). 40-Year Old Unbiased Distribution Free
#' Estimator Reliably Improves SEM Statistics for Nonnormal Data.
#' *Structural Equation Modeling: A Multidisciplinary Journal*, 29(6),
#' 872--887. \doi{10.1080/10705511.2022.2063870}
#'
#' Kenward, M. G., & Molenberghs, G. (1998). Likelihood based frequentist
#' inference when data are missing at random. *Statistical Science*, 13(3),
#' 236--247. \doi{10.1214/ss/1028905886}
#'
#' Bollen, K. A. (1989). *Structural Equations with Latent Variables*.
#' John Wiley & Sons. \doi{10.1002/9781118619179}
#'
#' Browne, M. W. (1974). Generalized least squares estimators in the analysis
#' of covariance structures. *South African Statistical Journal*, 8, 1--24.
pvalues <- function(object,
                    tests = if (is_classic_nt(object)) "pEBA4_RLS" else "pEBA4",
                    fiml.convention = c("observed", "lavaan")) {
  # The default is fit-appropriate: the historical RLS default for the classical
  # normal-theory ML case, the suffix-free form (standard statistic, biased
  # gamma) otherwise.
  check_lavaan(object, "object")
  options <- parse_tests(tests)
  check_supported(object)
  fiml.convention <- match.arg(fiml.convention)
  result <- unlist(lapply(options, function(option) {
    do.call(
      compute_pvalues,
      c(
        list(m0 = object),
        option,
        list(fiml.convention = fiml.convention)
      )
    )
  }))
  provenance_fiml <- if (is_fiml(object)) fiml.convention else NA
  as_semtests(
    result,
    fit_provenance(
      object,
      nested = FALSE,
      fiml.convention = provenance_fiml,
      df = fit_df(object),
      tests = tests,
      parsed_options = options
    )
  )
}

#' @rdname pvalues
#' @export
pvalues_nested <- function(m0, m1, method = c("2000", "2001"),
                           tests = if (is_classic_nt(m0)) "PALL_UG_ML" else "PALL",
                           A.method = c("delta", "exact"),
                           fiml.convention = c("observed", "lavaan")) {
  method <- match.arg(method)
  A.method <- match.arg(A.method)
  fiml.convention <- match.arg(fiml.convention)

  # Class gate first: the df computation below reads m0@test / m1@test, which
  # would otherwise die on a cryptic S4 slot error for a non-lavaan argument.
  check_lavaan(m0, "m0")
  check_lavaan(m1, "m1")
  options <- parse_tests(tests)

  m <- fit_df(m0) - fit_df(m1)

  # Check for identical df setting
  if (m == 0L) {
    semtests_abort(
      "Nested models must have different degrees of freedom.",
      "semTests_error_incompatible_models"
    )
  } else if (m < 0) {
    semtests_warn(
      paste0(
        "`m0` has fewer degrees of freedom than `m1`; treating `m1` as the ",
        "constrained model and swapping the inputs."
      ),
      "semTests_warning_model_order"
    )
    tmp <- m1
    m1 <- m0
    m0 <- tmp
  }
  df <- fit_df(m0) - fit_df(m1)

  check_supported_nested(m0, m1, method, A.method)

  result <- unlist(lapply(options, function(option) {
    do.call(
      compute_pvalues,
      c(
        list(m0 = m0, m1 = m1),
        option,
        list(
          A.method = A.method,
          fiml.convention = fiml.convention
        )
      )
    )
  }))
  provenance_A_method <- if (is_fiml(m0) || is_fiml(m1) ||
    isTRUE(m0@Model@categorical) ||
    isTRUE(m1@Model@categorical)) {
    A.method
  } else {
    NA
  }
  provenance_fiml <- if (is_fiml(m0) || is_fiml(m1)) fiml.convention else NA
  as_semtests(
    result,
    fit_provenance(
      m0,
      nested = TRUE,
      method = method,
      A.method = provenance_A_method,
      fiml.convention = provenance_fiml,
      df = df,
      tests = tests,
      parsed_options = options
    )
  )
}

#' Provenance of a `semTests_pvalues` object.
#'
#' Records the fit-level options actually used to compute the p-values, so the
#' returned object is self-describing across estimators and data types.
#' @keywords internal
fit_provenance <- function(fit, nested, method = NA, A.method = NA,
                           fiml.convention = NA, df = NULL, tests = NULL,
                           parsed_options = NULL) {
  estimator_requested <- fit@Options$estimator.orig
  if (is.null(estimator_requested) || !length(estimator_requested) ||
    is.na(estimator_requested[1])) {
    estimator_requested <- fit@Options$estimator
  }
  out <- list(
    estimator           = fit@Options$estimator,
    estimator_requested = estimator_requested,
    lavaan_test         = names(fit@test),
    information         = fit@Options$information,
    data_type           = if (isTRUE(fit@Model@categorical)) "categorical" else "continuous",
    missing             = fit@Options$missing,
    df                  = if (is.null(df)) unname(lavaan::fitmeasures(fit, "df")) else df,
    nested              = nested,
    method              = method
  )
  if (isTRUE(fit@Model@categorical)) {
    out$parameterization <- fit@Options$parameterization
    out$spectrum_source <- "lavaan UGamma"
  }
  if (!is.na(A.method[1])) out$A.method <- A.method
  if (!is.na(fiml.convention[1])) out$fiml.convention <- fiml.convention
  if (!is.null(tests) && !is.null(parsed_options)) {
    resolved_statistic <- vapply(parsed_options, function(option) {
      if (option$chisq == "auto") {
        if (is_classic_nt(fit)) "rls" else "ml"
      } else {
        option$chisq
      }
    }, character(1))
    gamma_type <- vapply(parsed_options, function(option) {
      if (identical(option$trad, "std")) {
        "not used"
      } else if (option$unbiased == 2L) {
        "unbiased"
      } else {
        "biased"
      }
    }, character(1))
    out$requested_tests <- as.character(tests)
    out$statistic <- stats::setNames(resolved_statistic, tests)
    out$gamma <- stats::setNames(gamma_type, tests)
  }
  out
}

#' @keywords internal
as_semtests <- function(x, info) {
  attr(x, "semtests") <- info
  class(x) <- c("semTests_pvalues", class(x))
  x
}

#' Print method for p-values from [pvalues()] / [pvalues_nested()].
#'
#' Prints the p-values, then a one-line provenance footer (estimator, data type,
#' information, df) recording the options used.
#' @param x A `semTests_pvalues` object.
#' @param ... Passed to the default print method.
#' @return `x`, invisibly.
#' @exportS3Method
print.semTests_pvalues <- function(x, ...) {
  info <- attr(x, "semtests")
  y <- unclass(x)
  attr(y, "semtests") <- NULL
  print.default(y, ...)
  if (!is.null(info)) {
    fiml <- !is.null(info$missing) && info$missing[1] %in% c("ml", "fiml")
    requested <- if (!is.null(info$estimator_requested) &&
      info$estimator_requested[1] != info$estimator[1]) {
      sprintf(" (%s)", info$estimator_requested[1])
    } else {
      ""
    }
    # lavaan stores `information` as a length-2 vector; show the (single) choice.
    cat(sprintf(
      "estimator: %s%s%s | data: %s | information: %s | df: %s%s%s\n",
      info$estimator[1], requested, if (fiml) " (FIML)" else "",
      info$data_type[1], info$information[1], info$df[1],
      if (!is.null(info$fiml.convention)) {
        sprintf(" | FIML convention: %s", info$fiml.convention[1])
      } else {
        ""
      },
      if (isTRUE(info$nested)) {
        a_method <- if (!is.null(info$A.method) &&
          !is.na(info$A.method[1])) {
          sprintf(", A.method %s", info$A.method[1])
        } else {
          ""
        }
        sprintf(" | nested (method %s%s)", info$method[1], a_method)
      } else {
        ""
      }
    ))
  }
  invisible(x)
}

#' Calculate traditional pvalues.
#' @param df,chisq,lambdas,type Parameters needed to calculate the p-values.
#' @returns Traditional p-values.
#' @keywords internal
trad_pvalue <- function(df, chisq, lambdas, type = c("std", "sf", "ss", "sb", "pall", "all")) {
  type <- match.arg(type)
  if (type == "std") {
    return(1 - stats::pchisq(chisq, df))
  }
  if (type == "sf") {
    return(scaled_f(chisq, lambdas))
  }
  if (type == "ss") {
    return(scaled_and_shifted(chisq, lambdas))
  }
  if (type == "sb") {
    m <- length(lambdas)
    return(as.numeric(1 - stats::pchisq(chisq * m / sum(lambdas), df = m)))
  }

  if (type == "all") {
    return(as.numeric(pvalue_all(chisq, lambdas)))
  }

  as.numeric(pall(chisq, lambdas))
}

#' Compute p-values for one parsed test specification.
#' @keywords internal
compute_pvalues <- function(m0, unbiased, trad, eba, peba, pols,
                            chisq = c("ml", "rls"), m1 = NULL,
                            A.method = "delta",
                            fiml.convention = "observed") {
  use_trad <- setdiff(trad, "std")
  if (is.null(m1)) {
    df <- fit_df(m0)
    chisqs <- make_chisqs(chisq, m0)
    if (is_fiml(m0)) {
      if (unbiased != 1) {
        stop("The unbiased (Du-Bentler) gamma is not defined for FIML; ",
          "drop `UG` from the test name to use the biased FIML spectrum. ",
          "See `?semTests-support`.",
          call. = FALSE
        )
      }
      lambdas_list <- fiml_lambdas(m0, df, fiml.convention)
    } else if (unbiased == 1) {
      # lavaan already exposes the biased single-model UGamma, including all
      # threshold, group-weight, constraint, and pairwise-missingness details.
      # Use it directly instead of rebuilding U %*% Gamma.
      lambdas_list <- lavaan_lambdas(m0, df)
    } else {
      ug_list <- ugamma(m0, unbiased)
      lambdas_list <- lapply(ug_list, function(ug) {
        ugamma_eigenvalues(
          ug, df,
          context = "single-model unbiased UGamma"
        )
      })
    }
  } else {
    # The nested reduction (method 2000) is estimator-agnostic for continuous
    # data -- the estimator enters only through WLS.V and the information matrix.
    # Fit-shape support (categorical / missing / FIML) is enforced upstream by
    # check_supported_nested(); here we branch on the FIML vs complete-data
    # computation and reject the FIML-incompatible UG gamma (a per-test option).
    df <- fit_df(m0) - fit_df(m1)
    chisqs <- make_chisqs(chisq, m0, m1)
    if (is_fiml(m0) || is_fiml(m1)) {
      if (unbiased != 1) {
        stop("The unbiased (Du-Bentler) gamma is not defined for FIML; ",
          "drop `UG` from the test name to use the biased FIML spectrum. ",
          "See `?semTests-support`.",
          call. = FALSE
        )
      }
      lambdas_list <- fiml_lambdas_nested(
        m0, m1, df,
        A.method = A.method,
        fiml.convention = fiml.convention
      )
    } else {
      lambdas_list <- lambdas_nested(m0, m1, unbiased, df)
    }

    if (min(unlist(lambdas_list)) < 0) {
      semtests_abort(
        paste0(
          "Nested comparison produced materially negative leading ",
          "eigenvalues. The fits may be unstable or not genuinely nested; ",
          "inspect them before trusting the p-values."
        ),
        "semTests_error_unstable_spectrum"
      )
    }
  }

  block_options <- list(EBA = eba, pEBA = peba)
  for (label in names(block_options)) {
    blocks <- block_options[[label]]
    if (!is.null(blocks) && any(blocks > df)) {
      semtests_abort(
        paste0(
          label, " cannot use more blocks than the test degrees ",
          "of freedom (", df, ")."
        ),
        "semTests_error_invalid_tests"
      )
    }
  }

  return_value <- c()
  for (i in seq_along(chisqs)) {
    chisq <- chisqs[i]
    result <- unlist(lapply(seq_along(lambdas_list), function(j) {
      lambdas <- lambdas_list[[j]]

      if (!is.null(peba)) {
        ppeba <- sapply(peba, function(k) peba_pvalue(chisq, lambdas, k))
        names(ppeba) <- paste0("peba", peba)
      } else {
        ppeba <- NULL
      }

      if (!is.null(eba)) {
        peba <- sapply(eba, function(k) eba_pvalue(chisq, lambdas, k))
        names(peba) <- paste0("eba", eba)
      } else {
        peba <- NULL
      }

      if (!is.null(pols)) {
        ppols <- sapply(pols, function(gamma) pols_pvalue(chisq, lambdas, gamma))
        names(ppols) <- paste0("pols", pols)
      } else {
        ppols <- NULL
      }

      ptrad <- sapply(use_trad, function(x) trad_pvalue(df, chisq, lambdas, x))
      names(ptrad) <- use_trad

      out <- pmin(pmax(c(ptrad, peba, ppeba, ppols), 0), 1)
      name <- if (names(lambdas_list)[[j]] == "ug_biased") "" else "_ug"
      name <- paste0(name, "_", names(chisqs)[i])

      if (length(out) != 0) {
        names(out) <- paste0(names(out), name)
      }
      out
    }))

    if ("std" %in% trad) {
      std <- c(trad_pvalue(df, chisq, NULL, "std"))
      names(std) <- paste0("std_", names(chisqs)[i])
      result <- c(std, result)
    }

    return_value <- c(return_value, result)
  }

  return_value
}
