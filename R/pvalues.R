#' Calculate *p*-values for one or two `lavaan` objects.
#'
#' Calculate *p*-values for a `lavaan` object using several methods,
#' including penalized eigenvalue block-averaging and penalized regression
#' estimators. The recommended choices of *p*-values are included as default
#' values. Multiple *p*-values can be returned simultaneously.
#'
#' The test argument is a list of character strings on the form
#' (test)_(ug?)_(ml?), for instance, `SB_UG_RLS`.
#'
#' * The first part of the string specifies the desired test. The supported tests are listed below.
#' * If `UG` is included in the string the unbiased estimator of the
#' fourth order moment matrix (Du, Bentler, 2022) is used. If not, the
#' standard biased matrix is used. There is no simple relationship between
#' *p*-value performance and the choice of `unbiased`.
#' * The final part of specifies the chi square statistic. The `ML`
#' choice uses the chi square based on the normal discrepancy function (Bollen, 2014).
#' The `RLS` choice (default) uses the reweighted least squares statistic of Browne (1974).
#'
#' The `eba` method partitions the eigenvalues into `j` equally sized sets
#' (if not possible, the smallest set is incomplete), and takes the mean
#' eigenvalue of these sets. Provide a list of integers `j` to partition
#' with respect to. The method was proposed by Foldnes & Grønneberg (2018).
#' `eba` with `j=2` -- `j=4` appear to work best.
#'
#' The `peba` method is a penalized variant of `eba`, described in
#' (Foldnes, Moss, Grønneberg, 2024). It typically outperforms `eba`, and
#' the best choice of `j` are typically about `2`--`6`.
#'
#' `pols` is a penalized regression method with a penalization term from ranging
#' from 0 to infinity. Foldnes, Moss, Grønneberg (2024) studied `pols=2`, which
#' has good performance in a variety of contexts.
#'
#' In addition, you may specify a
#' * `std` the standard *p*-value where the choice of `chisq` is approximated by a chi square distribution.
#' * `sb` Satorra-Bentler *p*-value. The *p*-value proposed by Satorra and Bentler (1994).
#' * `ss` The scaled and shifted *p*-value proposed by Asparouhov & Muthén (2010).
#' * `sf` The scaled *F* *p*-value proposed by Wu and Lin (2016).
#'
#' The `unbiased` argument is `TRUE` if the the unbiased estimator of the
#' fourth order moment matrix (Du, Bentler, 2022) is used. If `FALSE`, the
#' standard biased matrix is used. There is no simple relationship between
#' p-value performance and the choice of `unbiased`.
#'
#' The `chisq` argument controls which basic test statistic is used. The `ml`
#' choice uses the chi square based on the normal discrepancy function (Bollen, 2014).
#' The `rls` choice uses the reweighted least squares statistic of Browne (1974).
#'
#' @param object,m0,m1 One or two `lavaan` objects. `pvalues` does goodness-of-fit testing on one object,
#'    `pvalues_nested` does hypothesis testing on two nested models.
#' @param tests A list of tests to evaluate on the
#'    form `"(test)_(ug?)_(rls?)"`; see the default arguments and details below..
#' @param method For nested models, choose between `2000` and `2001`. Note: `2001` and Satorra-Bentler will not correspond with the variant in the paper.
#' @name pvalues
#' @export
#' @examples
#' library("semTests")
#' model <- "A =~ A1+A2+A3+A4+A5;
#'           C =~ C1+C2+C3+C4+C5"
#' n <- 200
#' object <- lavaan::sem(model, psych::bfi[1:n, 1:10], estimator = "MLM")
#' pvalues(object)
#'
#' # For the pEBA6 method with biased gamma and ML chisq statistic:
#' pvalues(object, "pEBA6_ML")
#'
#' @return A named vector of p-values.
#'
#' @references
#'
#' Foldnes, N., Moss, J., & Grønneberg, S. (2024). Improved goodness of fit procedures for structural equation models. Structural Equation Modeling: A Multidisciplinary Journal, 1-13. https://doi.org/10.1080/10705511.2024.2372028
#'
#' Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and standard errors in covariance structure analysis. https://psycnet.apa.org/record/1996-97111-016
#'
#' Asparouhov, & Muthén. (2010). Simple second order chi-square correction. Mplus Technical Appendix. https://www.statmodel.com/download/WLSMV_new_chi21.pdf
#'
#' Wu, H., & Lin, J. (2016). A Scaled F Distribution as an Approximation to the Distribution of Test Statistics in Covariance Structure Analysis. Structural Equation Modeling. https://doi.org/10.1080/10705511.2015.1057733
#'
#' Foldnes, N., & Grønneberg, S. (2018). Approximating Test Statistics Using Eigenvalue Block Averaging. Structural Equation Modeling, 25(1), 101-114. https://doi.org/10.1080/10705511.2017.1373021
#'
#' Du, H., & Bentler, P. M. (2022). 40-Year Old Unbiased Distribution Free Estimator Reliably Improves SEM Statistics for Nonnormal Data. Structural Equation Modeling: A Multidisciplinary Journal, 29(6), 872-887. https://doi.org/10.1080/10705511.2022.2063870
#'
#' Bollen, K. A. (2014). Structural Equations with Latent Variables (Vol. 210). John Wiley & Sons. https://doi.org/10.1002/9781118619179
#'
#' Browne. (1974). Generalized least squares estimators in the analysis of covariance structures. South African Statistical Journal. https://doi.org/10.10520/aja0038271x_175
pvalues <- \(object, tests = c("SB_UG_RLS", "pEBA2_UG_RLS", "pEBA4_RLS", "pEBA6_RLS", "pOLS_RLS")) {
  pvalues_internal(object, tests)
}

#' @keywords internal
pvalues_internal <- \(object, tests = c("SB_UG_RLS", "pEBA2_UG_RLS", "pEBA4_RLS", "pEBA6_RLS", "pOLS_RLS"), trad = NULL, eba = NULL, peba = NULL, pols = NULL, unbiased = 1, chisq = c("ml"), extras = FALSE) {
  if (is.null(tests) && is.null(trad) && is.null(eba) && is.null(peba) && is.null(pols)) {
    stop("Please provide some p-values to calculate.")
  }
  if (is.null(tests)) {
    pvalues_(object, unbiased = unbiased, trad = trad, eba = eba, peba = peba, pols = pols, chisq = chisq, extras = extras)
  } else {
    options <- lapply(tests, \(test) split_input(test))
    sapply(options, \(option) do.call(pvalues_, c(object, option, extras = extras)))
  }
}

#' @rdname pvalues
#' @export
pvalues_nested <- \(m0, m1, method = c("2000", "2001"), tests = c("SB_UG_RLS", "pEBA2_UG_RLS", "pEBA4_RLS", "pEBA6_RLS", "pOLS_RLS")) {
  pvalues_nested_internal(m0, m1, method = method, tests = tests)
}

#' @keywords internal
pvalues_nested_internal <- \(m0, m1, method = c("2000", "2001"), tests = c("SB_UG_RLS", "pEBA2_UG_RLS", "pEBA4_RLS", "pEBA6_RLS", "pOLS_RLS"), trad = NULL, eba = NULL, peba = NULL, pols = NULL, unbiased = 1, chisq = "ml", extras = FALSE) {
  method <- match.arg(method)

  if (is.null(tests) && is.null(trad) && is.null(eba) && is.null(peba) && is.null(pols)) {
    stop("Please provide some p-values to calculate.")
  }

  m <- m0@test[[1]]$df - m1@test[[1]]$df

  # Check for identical df setting
  if (m == 0L) {
    stop("Cannot test models with the same degree of freedom.")
  } else if (m < 0) {
    m = m1
    m1 = m0
    m0 = m
  }

  if (is.null(tests)) {
    pvalues_(m0, m1, unbiased = unbiased, trad = trad, eba = eba, peba = peba, pols = pols, chisq = chisq, extras = extras, method = method)
  } else {
    options <- lapply(tests, \(test) split_input(test))
    sapply(options, \(option) do.call(pvalues_, c(m0, m1, option, extras = extras, method = method)))
  }
}

#' P value function for one and two arguments.
#'
#' @keywords internal
#' @name pvalue_internal
#' @return pvalues.
NULL

#' Calculate traditional pvalues.
#' @param df,chisq,lambdas,type Parameters needed to calculate the p-values.
#' @returns Traditional p-values.
#' @keywords internal
trad_pvalue <- \(df, chisq, lambdas, type = c("std", "sf", "ss", "sb", "pall", "all")) {
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

  if (type == "pall") {
    return(as.numeric(pall(chisq, lambdas)))
  }

}

#' @rdname pvalue_internal
pvalues_ <- \(m0, m1, unbiased, trad, eba, peba, pols, chisq = c("ml", "rls"), extras = FALSE, method) {
  use_trad <- setdiff(trad, "std")
  bad_2001 <- FALSE
  if (missing(m1)) {
    df <- lavaan::fitmeasures(m0, "df")
    chisqs <- make_chisqs(chisq, m0)
    ug_list <- ugamma(m0, unbiased)
    lambdas_list <- lapply(ug_list, \(ug) Re(eigen(ug, only.values = TRUE)$values)[seq(df)])
  } else {
    if (m0@Options$estimator != "ML" || m1@Options$estimator != "ML") {
      stop("Only the 'ML' estimator supported.")
    }
    df <- lavaan::fitmeasures(m0, "df") - lavaan::fitmeasures(m1, "df")
    chisqs <- make_chisqs(chisq, m0, m1)
    ug_list <- ugamma_nested(m0, m1, method, unbiased)
    ug_list <- lapply(ug_list, sparsify)
    lambdas_list <- lapply(ug_list, \(ug) Re(RSpectra::eigs(ug, k = df, which = "LR", opts = list(retvec = FALSE))$values))

    if(min(unlist(lambdas_list)) < 0) {
      warning("Negative eigenvalues encountered in the first df eigenvalues of UGamma, defaulting to method = '2000'.")
      ug_list <- ugamma_nested(m0, m1, "2000", unbiased)
      lambdas_list <- lapply(ug_list, \(ug) Re(RSpectra::eigs(ug, k = df, which = "LR", opts = list(retvec = FALSE))$values))
      bad_2001 <- TRUE
    }
  }


  return_value <- c()
  for (i in seq_along(chisqs)) {
    chisq <- chisqs[i]
    result <- unlist(lapply(seq_along(ug_list), \(j) {
      ug <- ug_list[[j]]
      lambdas <- lambdas_list[[j]]

      if (!is.null(peba)) {
        ppeba <- sapply(peba, \(k) peba_pvalue(chisq, lambdas, k))
        names(ppeba) <- paste0("peba", peba)
      } else {
        ppeba <- NULL
      }

      if (!is.null(eba)) {
        peba <- sapply(eba, \(k) eba_pvalue(chisq, lambdas, k))
        names(peba) <- paste0("eba", eba)
      } else {
        peba <- NULL
      }

      if (!is.null(pols)) {
        ppols <- sapply(pols, \(gamma) pols_pvalue(chisq, lambdas, gamma))
        names(ppols) <- paste0("pols", pols)
      } else {
        ppols <- NULL
      }

      ptrad <- sapply(use_trad, \(x) trad_pvalue(df, chisq, lambdas, x))
      names(ptrad) <- use_trad

      out <- pmax(c(ptrad, peba, ppeba, ppols), 0)
      name <- if (names(ug_list)[[j]] == "ug_biased") "" else "_ug"
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

  if(!missing(m1)) {
    attr(return_value, "bad_2001") = bad_2001
    attr(return_value, method) = method
  }

  if (extras) {
    n <- length(lambdas_list)
    names(lambdas_list) <- rep("lambda", n)
    if (unbiased == 2) {
      names(lambdas_list) <- c(rep("lambda_biased", n / 2), rep("lambda_unbiased", n / 2))
    }
    c(return_value, chisqs, lambdas_list)
  } else {
    return_value
  }
}
