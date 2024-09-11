#' Get gamma from a model.
#' @param m1 Model to extract gamma from.
#' @param unbiased Biased (1), unbiased (2), or both (3).
#' @param m0 Optional second model, used if `m1` does not work.
#' @return List of (un)biased gammas.
#' @keywords internal

gamma <- \(m1, unbiased = 1, m0 = NULL) {

  stopifnot(unbiased %in% c(1,2,3))
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
gamma_rescale <- \(gamma_list, m1) {
  if (length(gamma_list) > 0) {
    fg <- unlist(m1@SampleStats@nobs) / m1@SampleStats@ntotal
    for (i in (seq_along(gamma_list))) {
      gamma_list[[i]] <- 1 / fg[i] * gamma_list[[i]]
    }
    lavaan::lav_matrix_bdiag(gamma_list)
  }
}

#' @keywords internal
gamma_from_lavaan <- \(m1, m0 = NULL) {
  gamma_biased <- lavaan::lavInspect(m1, "gamma")

  if(is.null(gamma_biased)) {
    if (!is.null(m0)) {
      gamma_biased <- lavaan::lavInspect(m0, "gamma")
      if (is.null(gamma_biased)) {
        stop("Could not calculate the gamma matrix from `m0` or `m1`. Use either `estimator = \"MLM\"' or `test=\"satorra.bentler\"' when fitting your lavaan model.")
      }
    } else {
      stop("Could not calculate the gamma matrix from the lavaan object. Use either `estimator = \"MLM\"' or `test=\"satorra.bentler\"' when fitting your lavaan model.")
    }
  }

  if(!is.list(gamma_biased)) {
    list(gamma_biased)
  } else {
    gamma_biased
  }
}

#' Calculate unbiase gamma from gamma and object.
#' @param gammas List of gammas for each group.
#' @param object `lavaan` object that corresponds to gamma.
#' @return List of unbiased gammas.
#' @keywords internal
gamma_to_gamma_unbiased <- \(gammas, object) {

  p <- length(object@Data@ov$name)
  meanstructure <- object@Options$meanstructure


  groups <- (object@Data@ngroups > 1)
  access <- \(x, g) {
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
#' @param a The `A` matrix. If if `NULL`, gets calculated by
#'    `lavaan:::lav_test_diff_A` with `method = method`.
#' @param method Method passed to `lavaan:::lav_test_diff_A`.
#' @keywords internal
#' @return Ugamma for nested object.
lav_ugamma_nested_2000 <- \(m0, m1, gamma, a = NULL, method = "delta") {
  wls_v <- lavaan::lavTech(m1, "WLS.V")
  pi <- lavaan::lavInspect(m1, "delta")

  p_inv <- lavaan::lavInspect(m1, what = "inverted.information")

  if (is.null(a)) {
    a <- do.call(lavaan:::lav_test_diff_A, list(m1, m0, method = "delta", reference = "H1"))
    if (m1@Model@eq.constraints) {
      a <- a %*% t(m1@Model@eq.constraints.K)
    }
  }

  paapaap <- p_inv %*% t(a) %*% MASS::ginv(a %*% p_inv %*% t(a)) %*% a %*% p_inv

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
ugamma <- \(object, unbiased = 1) {
  u0 <- lavaan::lavInspect(object, "U")
  gamma_list <- gamma(object, unbiased)
  lapply(gamma_list, \(gamma) u0 %*% gamma)
}

#' Calculate nested gamma
#' @keywords internal
ugamma_nested <- \(m0, m1, method = c("2000", "2001"), unbiased = 1) {
  method <- match.arg(method)
  gamma_list <- gamma(m1, unbiased, m0)
  f <- \(gamma) {
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

