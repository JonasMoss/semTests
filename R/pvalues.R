#' Calculate p-values for one or two lavaan objects.
#'
#' Calculate p-values for a `lavaan` object using several methods,
#' including penalized eigenvalue block-averaging and penalized regression
#' estimators. The choice `peba=4` together with `chisq = "rls"` and `ub`
#' is recommended. Multiple p-values can be returned simultaneously.
#'
#' The traditional methods include:
#' * `pstd` the standard *p*-value where the choice of `chisq` is approximated by a chi square distribution.
#' * `psb` Satorra-Bentler *p*-value. The *p*-value proposed by Satorra and Bentler (1994).
#' * `pss` The scaled and shifted *p*-value proposed by Asparouhov & Muthén (2010).
#' * `pcf` The Scaled F *p*-value proposed by Wu and Lin (2016).
#' * `pfull` *p*-value based on all eigenvalues of the asymptotic covariance matrix matrix.
#'
#' The `eba` method partitions the eigenvalues into `j` equally sized sets
#' (if not possible, the smallest set is incomplete), and takes the mean
#' eigenvalue of these sets. Provide a list of integers `j` to partition
#' with respect to. The method was proposed by Foldnes & Grønneberg (2018).
#' `eba` with `j=2` or `j=4` appear to work best.
#'
#' The `peba` method is a penalized variant of `eba`, described in
#' (Foldnes, Moss, Grønneberg, WIP). It typically outperforms `eba`, and
#' the best choice of `j` is typically `6`.
#'
#' `pols` is a penalized regression method with a penalization term from ranging
#' from 0 to infitity. Foldnes, Moss, Grønneberg (WIP) studied `pols=2`, which
#' has good performance in a variety of contexts.
#'
#' The `unbiased` argument is `TRUE` if the the unbiased estimator of the
#' fourth order moment matrix (Du, Bentler, 2022) is used. If `FALSE`, the
#' standard biased matrix is used. There is no simple relationship between
#' p-value performance and the choice of `unbiased`.
#'
#' The `chisq` argument controls which basic test statistic is used. The `trad`
#' choice uses the chi square based on the normal discrepancy function (Bollen, 2014).
#' The `rls` choice uses the reweighted least squares statistic of Browne (1974).
#'
#' @param object A `lavaan` object.
#' @param trad List of traditional p-values to calculate.
#'    Not calculated if `NULL.`
#' @param eba List of which `eba` p-values to calculate.
#'    Not calculated if `NULL.`
#' @param peba List of which `peba` p-values to calculate.
#'    Not calculated if `NULL.`
#' @param pols List of penalization parameters to use in the penalized
#'    OLS p-value. Not calculated if `NULL.`
#' @param unbiased A number between 1 and 3. 1: Calculate using the biased
#'    gamma matrix (default). 2: Calculate using the unbiased gamma matrix.
#'    3: Calculate using both gammas.
#' @param chisq Which chi-square statistic to base the calculations on.
#' @param extras Returns the estimated eigenvalues and basic test statistics
#'    if checked.
#' @export
#' @return A named vector of p-values.
#'
#' @references
#' Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and standard errors in covariance structure analysis. https://psycnet.apa.org/record/1996-97111-016
#'
#' Asparouhov, & Muthén. (2010). Simple second order chi-square correction. Mplus Technical Appendix. https://www.statmodel.com/download/WLSMV_new_chi21.pdf
#'
#' Wu, H., & Lin, J. (2016). A Scaled F Distribution as an Approximation to the Distribution of Test Statistics in Covariance Structure Analysis. Structural Equation Modeling. https://doi.org/10.1080/10705511.2015.1057733
#'
#' Foldnes, N., & Grønneberg, S. (2018). Approximating Test Statistics Using Eigenvalue Block Averaging. Structural Equation Modeling, 25(1), 101–114. https://doi.org/10.1080/10705511.2017.1373021
#'
#' Du, H., & Bentler, P. M. (2022). 40-Year Old Unbiased Distribution Free Estimator Reliably Improves SEM Statistics for Nonnormal Data. Structural Equation Modeling: A Multidisciplinary Journal, 29(6), 872–887. https://doi.org/10.1080/10705511.2022.2063870
#'
#' Bollen, K. A. (2014). Structural Equations with Latent Variables (Vol. 210). John Wiley & Sons. https://doi.org/10.1002/9781118619179
#'
#' Browne. (1974). Generalized least squares estimators in the analysis of covariance structures. South African Statistical Journal. https://doi.org/10.10520/aja0038271x_175
pvalues <- \(object, trad = NULL, eba = NULL, peba = c(2, 4), pols = 2, unbiased = 1, chisq = c("rls", "trad"), extras = FALSE) {
  if (is.null(trad) && is.null(eba) && is.null(peba) && is.null(pols)) {
    stop("Please provide some p-values to calculate.")
  }
  pvalues_one(object, unbiased = unbiased, trad = trad, eba = eba, peba = peba, pols = pols, chisq = chisq, extras = extras)
}

#' P value function for one and two arguments.
#'
#' @keywords internal
#' @param object,m0,m1 `lavaan` objects.
#' @param trad List of traditional p-values to calculate.
#' @param eba List of which eba p-values to calculate.
#' @name pvalue_internal
#' @return pvalues.
NULL

#' Calculate traditional pvalues.
#' @param df,chisq,lambda,type Parameters needed to calculate the p-values.
#' @returns Traditional p-values.
#' @keywords internal
trad_pvalue <- \(df, chisq, lambdas, type = c("pstd", "psf", "pss", "psb", "pfull")) {
  type <- match.arg(type)
  if (type == "pstd") {
    return(1 - stats::pchisq(chisq, df))
  }
  if (type == "psf") {
    return(scaled_f(chisq, lambdas))
  }
  if (type == "pss") {
    return(scaled_and_shifted(chisq, lambdas))
  }
  if (type == "pfull") {
    return(CompQuadForm::imhof(chisq, lambdas)$Qq)
  }
  if (type == "psb") {
    m <- length(lambdas)
    return(as.numeric(1 - stats::pchisq(chisq * m / sum(lambdas), df = m)))
  }
}

#' Calculate rls
#' @param object lavaan object.
#' @returns RLS.
#' @keywords internal
rls <- \(object) {
  s_inv <- chol2inv(chol(lavaan::lavInspect(object, "sigma.hat")))
  residuals <- lavaan::lavInspect(object, "residuals")
  mm <- residuals[[1]] %*% s_inv
  n <- lavaan::lavInspect(object, "nobs")
  .5 * n * sum(mm * t(mm))
}

#' @keywords internal
make_chisqs <- \(chisq, object) {
  chisqs <- c()
  if ("trad" %in% chisq) chisqs["trad"] <- lavaan::fitmeasures(object, "chisq")
  if ("rls" %in% chisq) chisqs["rls"] <- rls(object)
  chisqs
}

#' @rdname pvalue_internal
pvalues_one <- \(object, unbiased, trad, eba, peba, pols, chisq = c("trad", "rls"), extras = FALSE) {
  df <- lavaan::fitmeasures(object, "df")
  chisqs <- make_chisqs(chisq, object)
  use_trad <- setdiff(trad, "pstd")
  ug_list <- ugamma_no_groups(object, unbiased)
  lambdas_list <- lapply(ug_list, \(ug) Re(eigen(ug)$values)[seq(df)])
  return_value <- c()
  for (i in seq_along(chisqs)) {
    chisq <- chisqs[i]
    result <- unlist(lapply(seq_along(ug_list), \(j) {
      ug <- ug_list[[j]]
      lambdas <- lambdas_list[[j]]

      if (!is.null(peba)) {
        ppeba <- sapply(peba, \(k) peba_pvalue(chisq, lambdas, k))
        names(ppeba) <- paste0("ppeba", peba)
      } else {
        ppeba <- NULL
      }

      if (!is.null(eba)) {
        peba <- sapply(eba, \(k) eba_pvalue(chisq, lambdas, k))
        names(peba) <- paste0("peba", eba)
      } else {
        peba <- NULL
      }

      if (!is.null(pols)) {
        ppols <- sapply(pols, \(gamma) pols_pvalue(chisq, lambdas, gamma))
        names(ppols) <- paste0("pols_", pols)
      } else {
        ppols <- NULL
      }

      ptrad <- sapply(use_trad, \(x) trad_pvalue(df, chisq, lambdas, x))
      names(ptrad) <- use_trad

      out <- pmax(c(ptrad, peba, ppeba, ppols), 0)
      name <- if (names(ug_list)[[j]] == "ug_biased") "" else "_ub"
      name <- paste0(name, "_", names(chisqs)[i])

      if (length(out) != 0) {
        names(out) <- paste0(names(out), name)
      }
      out
    }))

    if ("pstd" %in% trad) {
      pstd <- c(trad_pvalue(df, chisq, NULL, "pstd"))
      names(pstd) <- paste0("pstd_", names(chisqs)[i])
      result <- c(pstd, result)
    }

    return_value <- c(return_value, result)
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

#' Calculate the scaled and shifted / the mean-variance adjusted p-value
#'
#' @param chisq Chi-square fit value from a lavaan object.
#' @param lambdas Eigenvalues of UG matrix.
#' @name laavan_tests
#' @keywords internal
#' @return The scaled and shifted p-value or the mean-variance adjusted p-value.
NULL

#' @rdname laavan_tests
scaled_and_shifted <- \(chisq, lambdas) {
  df <- length(lambdas)
  tr_ug <- sum(lambdas)
  tr_ug2 <- sum(lambdas^2)
  a <- sqrt(df / tr_ug2)
  b <- df - sqrt(df * tr_ug^2 / tr_ug2)
  t3 <- unname(chisq * a + b)
  1 - stats::pchisq(t3, df = df)
}

#' Calculate the scaled_f p-value.
#' @param chisq Chi-square fit value from a lavaan object.
#' @param eig eig of UG matrix.
#' @return scaled f p-value.
#' @keywords internal
scaled_f <- \(chisq, eig) {
  s1 <- sum(eig)
  s2 <- sum(eig^2)
  s3 <- sum(eig^3)

  denom <- 2 * s1 * s2^2 - s1^2 * s3 + 2 * s2 * s3
  if (denom > 0) {
    d1f3 <- s1 * (s1^2 * s2 - 2 * s2^2 + 4 * s1 * s3) / denom
    d2f3 <- (s1^2 * s2 + 2 * s2^2) / (s3 * s1 - s2^2) + 6
    if (d2f3 < 6) d2f3 <- Inf
    cf3 <- s1 * (s1^2 * s2 - 2 * s2^2 + 4 * s1 * s3) /
      (s1^2 * s2 - 4 * s2^2 + 6 * s1 * s3)
  } else {
    d1f3 <- Inf
    d2f3 <- s1^2 / s2 + 4
    cf3 <- s1 * (s1^2 + 2 * s2) / (s1^2 + 4 * s2)
  }

  unname(1 - stats::pf(chisq / cf3, d1f3, d2f3))
}

#' Calculate the jth eba pvalue.
#' @keywords internal
eba_pvalue <- \(chisq, lambdas, j) {
  m <- length(lambdas)
  k <- ceiling(m / j)
  eig <- lambdas
  eig <- c(eig, rep(NA, k * j - length(eig)))
  dim(eig) <- c(k, j)
  eig_means <- colMeans(eig, na.rm = TRUE)
  repeated <- rep(eig_means, each = k)[seq(m)]
  CompQuadForm::imhof(chisq, repeated)$Qq
}

#' Calculate the jth eba pvalue.
#' @keywords internal
peba_pvalue <- \(chisq, lambdas, j) {
  m <- length(lambdas)
  k <- ceiling(m / j)
  eig <- lambdas
  eig <- c(eig, rep(NA, k * j - length(eig)))
  dim(eig) <- c(k, j)
  eig_means <- colMeans(eig, na.rm = TRUE)
  eig_mean <- mean(lambdas)
  repeated <- rep(eig_means, each = k)[seq(m)]
  CompQuadForm::imhof(chisq, (repeated + eig_mean) / 2)$Qq
}

#' Calculate penalized OLS pvalue.
#' @keywords internal
pols_pvalue <- \(chisq, lambdas, gamma) {
  x <- seq_along(lambdas)
  beta1_hat <- 1 / gamma * stats::cov(x, lambdas) / stats::var(x)
  beta0_hat <- mean(lambdas) - beta1_hat * mean(x)
  lambda_hat <- pmax(beta0_hat + beta1_hat * x, 0)
  CompQuadForm::imhof(chisq, lambda_hat)$Qq
}

#' Calculate non-nested gamma without mean structure
#' @keywords internal
ugamma_no_groups <- \(object, unbiased = 1) {
  u <- lavaan::lavInspect(object, "U")
  object@Options$gamma.unbiased <- FALSE
  gamma <- lavaan::lavInspect(object, "gamma")

  out <- list()

  if (unbiased == 1 || unbiased == 3) {
    out <- list(ug_biased = u %*% gamma)
  }

  if (unbiased == 2 || unbiased == 3) {
    gamma_unb <- gamma_unbiased(object, gamma)
    out <- c(out, list(ug_unbiased = u %*% gamma_unb))
  }

  out
}

