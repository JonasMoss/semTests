#' Get gamma from a model.
#' @param m1 Model to extract gamma from.
#' @param unbiased Biased (1), unbiased (2), or both (3).
#' @param m0 Optional second model, used if `m1` does not work.
#' @return List of (un)biased gammas.
#' @keywords internal

gamma <- function(m1, unbiased = 1, m0 = NULL) {

  stopifnot(unbiased %in% c(1,2,3))
  if (unbiased %in% c(2, 3) && !is_classic_nt(m1)) {
    stop("The unbiased (Du-Bentler) gamma is only defined for continuous, ",
         "complete-data ML; drop `UG` from the test name to use the biased gamma. ",
         "See `?semTests-support`.",
         call. = FALSE)
  }
  gamma_biased <- gamma_from_lavaan(m1, m0)

  gamma_list = list()

  if (unbiased %in% c(2, 3)) {
    gamma_unbiased = gamma_to_gamma_unbiased(gamma_biased, m1)
    gamma_unbiased_rescaled = gamma_rescale(gamma_unbiased, m1)
    gamma_list <- c(gamma_list, list(ug_unbiased = gamma_unbiased_rescaled))
  }

  if (unbiased %in% c(1, 3)) {
    gamma_biased_rescaled = gamma_rescale(gamma_biased, m1)
    gamma_list <- c(gamma_list, list(ug_biased = gamma_biased_rescaled))
  }

  gamma_list
}

#' @keywords internal
gamma_rescale <- function(gamma_list, m1) {
  if (length(gamma_list) > 0) {
    fg <- unlist(m1@SampleStats@nobs) / m1@SampleStats@ntotal
    for (i in (seq_along(gamma_list))) {
      gamma_list[[i]] <- 1 / fg[i] * gamma_list[[i]]
    }
    lavaan::lav_matrix_bdiag(gamma_list)
  }
}

#' @keywords internal
gamma_from_lavaan <- function(m1, m0 = NULL) {
  inspect_gamma <- function(fit) {
    tryCatch(
      lavaan::lavInspect(fit, "gamma"),
      error = function(e) NULL
    )
  }
  gamma_biased <- inspect_gamma(m1)

  if(is.null(gamma_biased)) {
    if (!is.null(m0)) {
      gamma_biased <- inspect_gamma(m0)
      if (is.null(gamma_biased)) {
        semtests_abort(
          paste0("lavaan did not expose a Gamma/NACOV matrix for `m0` or `m1`. ",
                 "Refit both models with a compatible robust test or estimator ",
                 "(for example, `estimator = \"MLM\"` for ML, or ",
                 "`test = \"satorra.bentler\"` where supported)."),
          "semTests_error_missing_gamma"
        )
      }
    } else {
      semtests_abort(
        paste0("lavaan did not expose a Gamma/NACOV matrix for `object`. Refit ",
               "the model with a compatible robust test or estimator (for ",
               "example, `estimator = \"MLM\"` for ML, or ",
               "`test = \"satorra.bentler\"` where supported)."),
        "semTests_error_missing_gamma"
      )
    }
  }

  if(!is.list(gamma_biased)) {
    list(gamma_biased)
  } else {
    gamma_biased
  }
}

#' Calculate unbiased gamma from gamma and object.
#' @param gammas List of gammas for each group.
#' @param object `lavaan` object that corresponds to gamma.
#' @return List of unbiased gammas.
#' @keywords internal
gamma_to_gamma_unbiased <- function(gammas, object) {

  p <- length(object@Data@ov$name)
  meanstructure <- object@Options$meanstructure


  groups <- (object@Data@ngroups > 1)
  access <- function(x, g) {
    if(groups) {
      x[[g]]
    } else {
      x
    }
  }

  gamma_list = list()
  for(g in seq(object@Data@ngroups)) {
    N <- lavaan::lavInspect(object, what = "nobs")[[g]]
    Gamma <- gammas[[g]]
    COV <- access(lavaan::lavInspect(object, what = "samplestats"), g)$cov
    cov.vech <- lavaan::lav_matrix_vech(COV)
    GammaNT.cov <- 2 * lavaan::lav_matrix_duplication_ginv_pre_post(COV %x% COV)
    Gamma.cov <- Gamma
    if (meanstructure) {
      Gamma.cov <- Gamma[-(1:p), -(1:p), drop = FALSE]
      Gamma.mean.cov <- Gamma[1:p, -(1:p), drop = FALSE]
    } else {
      Gamma.cov <- Gamma
    }

    Gamma.u <- (N * (N - 1) / (N - 2) / (N - 3) * Gamma.cov -
                  N / (N - 2) / (N - 3) * (GammaNT.cov -
                                             2 / (N - 1) * tcrossprod(cov.vech)))
    if (meanstructure) {
      Gamma <- lavaan::lav_matrix_bdiag(COV, Gamma.u)
      Gamma[1:p, (p + 1):ncol(Gamma)] <- Gamma.mean.cov * N / (N - 2)
      Gamma[(p + 1):ncol(Gamma), 1:p] <- t(Gamma.mean.cov * N / (N - 2))
    } else {
      Gamma <- Gamma.u
    }
    gamma_list = c(gamma_list, list(Gamma))
  }

  gamma_list
}

#' Calculate nested ugamma.
#'
#' This can also be used with restrictions.
#'
#' @param m0,m1 Two nested `lavaan` objects.
#' @param gamma Gamma weighted by groups.
#' @param a The `A` restriction matrix. If `NULL`, it is calculated from the
#'   two models' public inspected delta matrices.
#' @keywords internal
#' @return Ugamma for nested object.
lav_ugamma_nested_2000 <- function(m0, m1, gamma, a = NULL) {
  wls_v <- lavaan::lavTech(m1, "WLS.V")
  pi <- lavaan::lavInspect(m1, "delta")

  p_inv <- lavaan::lavInspect(m1, what = "inverted.information")

  if (is.null(a)) {
    a <- do.call(get_a_matrix, list(m1, m0))
    if (m1@Model@eq.constraints) {
      a <- a %*% t(m1@Model@eq.constraints.K)
    }
  }

  paapaap <- p_inv %*% t(a) %*% generalized_inverse(a %*% p_inv %*% t(a)) %*% a %*% p_inv

  # Compute scaling factor
  fg <- unlist(m1@SampleStats@nobs) / m1@SampleStats@ntotal

  # We need the global gamma, cf. eq.~(10) in Satorra (2000).
  # And the global V.
  # Note that the calculation is not the same as in lavaan(2000), as we
  # need more than the trace to be correct. (The simplications used in lavaan
  # do not hold for the whole ugamma.)

  for (i in (seq_along(wls_v))) {
    wls_v[[i]] <- fg[i] * wls_v[[i]]
  }

  v_global <- lavaan::lav_matrix_bdiag(wls_v)
  pi_global <- if (is.list(pi)) {
    do.call(rbind, pi)
  } else {
    pi
  }

  # U global version, eq.~(22) in Satorra (2000).
  u_global <- v_global %*% pi_global %*% paapaap %*% t(pi_global) %*% v_global
  return(u_global %*% gamma)
}


#' Calculate non-nested gamma
#' @keywords internal
ugamma <- function(object, unbiased = 1) {
  u0 <- lavaan::lavInspect(object, "U")
  gamma_list <- gamma(object, unbiased)
  lapply(gamma_list, function(gamma) u0 %*% gamma)
}

#' Extract a stable leading UGamma spectrum.
#'
#' The theoretical spectrum is real and non-negative. Small negative or
#' imaginary components arise from the nonsymmetric numerical representation
#' of a matrix that is similar to a symmetric positive-semidefinite product.
#' Numerical negatives are truncated; a materially negative leading value is
#' either retained for the method-2001 fallback or rejected.
#' @keywords internal
ugamma_eigenvalues <- function(ugamma, df, context = "UGamma",
                               allow_negative = FALSE) {
  ugamma <- as.matrix(ugamma)
  if (!is.numeric(ugamma) || length(dim(ugamma)) != 2L ||
      nrow(ugamma) != ncol(ugamma)) {
    stop(context, " must be a square numeric matrix.", call. = FALSE)
  }
  if (df <= 0L) return(numeric(0L))
  if (df > nrow(ugamma)) {
    stop(context, " dimension is smaller than the model degrees of freedom.",
         call. = FALSE)
  }
  if (any(!is.finite(ugamma))) {
    stop(context, " contains non-finite entries.", call. = FALSE)
  }

  values <- eigen(ugamma, only.values = TRUE)$values
  scale <- max(1, max(Mod(values)))
  tolerance <- 1e-8 * scale
  if (any(!is.finite(values))) {
    stop(context, " contains non-finite eigenvalues.", call. = FALSE)
  }
  if (max(abs(Im(values))) > tolerance) {
    stop(context, " has materially complex eigenvalues.", call. = FALSE)
  }

  lambdas <- sort(Re(values), decreasing = TRUE)[seq_len(df)]
  numerical_negative <- lambdas < 0 & lambdas >= -tolerance
  lambdas[numerical_negative] <- 0
  if (!allow_negative && any(lambdas < -tolerance)) {
    stop(context, " has materially negative leading eigenvalues; the fit may ",
         "be numerically unstable or the requested models may not be nested.",
         call. = FALSE)
  }
  if (!allow_negative && sum(lambdas > tolerance) < df) {
    stop(context, " rank is smaller than the model degrees of freedom.",
         call. = FALSE)
  }
  lambdas
}

#' Biased single-model spectrum exposed by lavaan.
#' @keywords internal
lavaan_lambdas <- function(object, df) {
  ugamma <- tryCatch(
    lavaan::lavInspect(object, "UGamma"),
    error = function(e) {
      stop("lavaan could not construct UGamma for this fit: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  list(
    ug_biased = ugamma_eigenvalues(
      ugamma, df,
      context = "lavaan's single-model UGamma"
    )
  )
}

#' Calculate nested gamma
#' @keywords internal
ugamma_nested <- function(m0, m1, method = c("2000", "2001"), unbiased = 1) {
  method <- match.arg(method)
  gamma_list <- gamma(m1, unbiased, m0)
  f <- function(gamma) {
    if (method == "2001") {
      u0 <- lavaan::lavInspect(m0, "U")
      u1 <- lavaan::lavInspect(m1, "U")
      (u0 - u1) %*% gamma
    } else {
      lav_ugamma_nested_2000(m0, m1, gamma)
    }
  }

  lapply(gamma_list, f)
}

#' Gamma-free factors of the reduced nested spectrum (Satorra 2000).
#'
#' Returns the restriction-space factor `D` and companion `C` such that the `m`
#' nonzero eigenvalues of the full `q x q` UGamma equal those of
#' `C^{-1} D' Gamma D`, an `m x m` problem with `m` the number of restrictions.
#' This is the materialised reduction of Moss (2026): the full `q x q` U matrix
#' and its eigendecomposition are never formed. With `U = D C^+ D'` one has
#' `D = V Delta P^+ A'` and `C = A P^+ A'`, where `V` is the (group-weighted)
#' weight, `Delta` the Jacobian, `P^+` the inverted information, and `A` the
#' restriction matrix.
#' @keywords internal
nested_factor_2000 <- function(m0, m1, a = NULL) {
  wls_v <- lavaan::lavTech(m1, "WLS.V")
  pi <- lavaan::lavInspect(m1, "delta")
  p_inv <- lavaan::lavInspect(m1, what = "inverted.information")

  if (is.null(a)) {
    a <- do.call(get_a_matrix, list(m1, m0))
    if (m1@Model@eq.constraints) {
      a <- a %*% t(m1@Model@eq.constraints.K)
    }
  }

  fg <- unlist(m1@SampleStats@nobs) / m1@SampleStats@ntotal
  for (i in seq_along(wls_v)) wls_v[[i]] <- fg[i] * wls_v[[i]]
  v_global <- lavaan::lav_matrix_bdiag(wls_v)
  pi_global <- if (is.list(pi)) do.call(rbind, pi) else pi

  list(
    D = v_global %*% pi_global %*% p_inv %*% t(a),
    C = a %*% p_inv %*% t(a)
  )
}

#' Reference spectrum of the nested test, without forming the full UGamma.
#'
#' For Satorra's (2000) method the nonzero eigenvalues are those of the `m x m`
#' matrix `C^{-1} D' Gamma D` (see [nested_factor_2000]); the `q x q` UGamma and
#' its eigendecomposition are avoided. The (2001) method has no such reduction,
#' so the top-`df` eigenvalues of the full `(U0 - U1) Gamma` are returned.
#' @param m0,m1 Two nested `lavaan` objects.
#' @param method Either `"2000"` or `"2001"`.
#' @param unbiased Biased (1), unbiased (2), or both (3) gamma.
#' @param df Number of restrictions (degrees-of-freedom difference).
#' @return A list of eigenvalue vectors, one per gamma estimate.
#' @keywords internal
lambdas_nested <- function(m0, m1, method = c("2000", "2001"), unbiased = 1, df) {
  method <- match.arg(method)
  gamma_list <- gamma(m1, unbiased, m0)

  if (method == "2001") {
    du <- lavaan::lavInspect(m0, "U") - lavaan::lavInspect(m1, "U")
    lapply(gamma_list, function(g) {
      ugamma_eigenvalues(
        du %*% g, df,
        context = "nested method-2001 UGamma",
        allow_negative = TRUE
      )
    })
  } else {
    factors <- nested_factor_2000(m0, m1)
    c_inv <- generalized_inverse(factors$C)
    lapply(gamma_list, function(g) {
      m_phi <- t(factors$D) %*% g %*% factors$D
      ugamma_eigenvalues(
        c_inv %*% m_phi, df,
        context = "nested method-2000 restriction UGamma"
      )
    })
  }
}
