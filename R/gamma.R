#' Get gamma from a model.
#' @param m1 Model to extract gamma from.
#' @param unbiased Biased (1), unbiased (2), or both (3).
#' @param m0 Optional second model, used if `m1` does not work.
#' @return List of (un)biased gammas.
#' @keywords internal

gamma_matrices <- function(m1, unbiased = 1, m0 = NULL) {
  stopifnot(unbiased %in% c(1, 2, 3))
  if (unbiased %in% c(2, 3) && !is_classic_nt(m1)) {
    stop("The unbiased (Du-Bentler) gamma is only defined for continuous, ",
      "complete-data ML with random observed exogenous covariates ",
      "(`fixed.x = FALSE`, `conditional.x = FALSE`). Drop `UG` from the ",
      "test name to use the biased gamma. ",
      "See `?semTests-support`.",
      call. = FALSE
    )
  }
  gamma_biased <- gamma_from_lavaan(m1, m0)

  gamma_list <- list()

  if (unbiased %in% c(2, 3)) {
    gamma_unbiased <- gamma_to_gamma_unbiased(gamma_biased, m1)
    gamma_unbiased_rescaled <- gamma_rescale(gamma_unbiased, m1)
    gamma_list <- c(gamma_list, list(ug_unbiased = gamma_unbiased_rescaled))
  }

  if (unbiased %in% c(1, 3)) {
    gamma_biased_rescaled <- gamma_rescale(gamma_biased, m1)
    gamma_list <- c(gamma_list, list(ug_biased = gamma_biased_rescaled))
  }

  gamma_list
}

#' @keywords internal
gamma_rescale <- function(gamma_list, m1) {
  if (!length(gamma_list)) {
    return(NULL)
  }
  group_weights <- unlist(m1@SampleStats@nobs) / m1@SampleStats@ntotal
  gamma_list <- Map(
    function(gamma_matrix, weight) gamma_matrix / weight,
    gamma_list,
    group_weights
  )
  lavaan::lav_matrix_bdiag(gamma_list)
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

  if (is.null(gamma_biased)) {
    if (!is.null(m0)) {
      gamma_biased <- inspect_gamma(m0)
      if (is.null(gamma_biased)) {
        semtests_abort(
          paste0(
            "lavaan did not expose a Gamma/NACOV matrix for `m0` or `m1`. ",
            "Refit both models with a compatible robust test or estimator ",
            "(for example, `estimator = \"MLM\"` for ML, or ",
            "`test = \"satorra.bentler\"` where supported)."
          ),
          "semTests_error_missing_gamma"
        )
      }
    } else {
      semtests_abort(
        paste0(
          "lavaan did not expose a Gamma/NACOV matrix for `object`. Refit ",
          "the model with a compatible robust test or estimator (for ",
          "example, `estimator = \"MLM\"` for ML, or ",
          "`test = \"satorra.bentler\"` where supported)."
        ),
        "semTests_error_missing_gamma"
      )
    }
  }

  if (!is.list(gamma_biased)) {
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
    if (groups) {
      x[[g]]
    } else {
      x
    }
  }

  sample_sizes <- lavaan::lavInspect(object, what = "nobs")
  sample_stats <- lavaan::lavInspect(object, what = "samplestats")
  gamma_list <- vector("list", object@Data@ngroups)
  for (g in seq_len(object@Data@ngroups)) {
    sample_size <- sample_sizes[[g]]
    gamma_matrix <- gammas[[g]]
    covariance <- access(sample_stats, g)$cov
    covariance_vech <- lavaan::lav_matrix_vech(covariance)
    gamma_nt_covariance <- 2 *
      lavaan::lav_matrix_duplication_ginv_pre_post(
        covariance %x% covariance
      )
    if (meanstructure) {
      gamma_covariance <- gamma_matrix[-seq_len(p), -seq_len(p), drop = FALSE]
      gamma_mean_covariance <- gamma_matrix[
        seq_len(p),
        -seq_len(p),
        drop = FALSE
      ]
    } else {
      gamma_covariance <- gamma_matrix
    }

    gamma_unbiased <- (
      sample_size * (sample_size - 1) /
        (sample_size - 2) / (sample_size - 3) * gamma_covariance -
        sample_size / (sample_size - 2) / (sample_size - 3) *
          (
            gamma_nt_covariance -
              2 / (sample_size - 1) * tcrossprod(covariance_vech)
          )
    )
    if (meanstructure) {
      gamma_matrix <- lavaan::lav_matrix_bdiag(covariance, gamma_unbiased)
      covariance_columns <- (p + 1L):ncol(gamma_matrix)
      correction <- gamma_mean_covariance *
        sample_size / (sample_size - 2)
      gamma_matrix[seq_len(p), covariance_columns] <- correction
      gamma_matrix[covariance_columns, seq_len(p)] <- t(correction)
    } else {
      gamma_matrix <- gamma_unbiased
    }
    gamma_list[[g]] <- gamma_matrix
  }

  gamma_list
}

#' Calculate nested ugamma.
#'
#' This can also be used with restrictions.
#'
#' @param m0,m1 Two nested `lavaan` objects.
#' @param gamma_matrix Gamma weighted by groups.
#' @param a The `A` restriction matrix. If `NULL`, it is calculated from the
#'   two models' public inspected delta matrices.
#' @keywords internal
#' @return Ugamma for nested object.
lav_ugamma_nested_2000 <- function(m0, m1, gamma_matrix, a = NULL) {
  factors <- nested_factor_2000(m0, m1, a)
  factors$D %*%
    generalized_inverse(factors$C) %*%
    t(factors$D) %*%
    gamma_matrix
}


#' Calculate non-nested gamma
#' @keywords internal
ugamma <- function(object, unbiased = 1) {
  u0 <- lavaan::lavInspect(object, "U")
  gamma_list <- gamma_matrices(object, unbiased)
  lapply(gamma_list, function(gamma_matrix) u0 %*% gamma_matrix)
}

#' Extract a stable leading UGamma spectrum.
#'
#' The theoretical spectrum is real and non-negative. Small negative or
#' imaginary components arise from the nonsymmetric numerical representation
#' of a matrix that is similar to a symmetric positive-semidefinite product.
#' Numerical negatives are truncated; materially negative values are rejected.
#' @keywords internal
ugamma_eigenvalues <- function(ugamma, df, context = "UGamma") {
  ugamma <- as.matrix(ugamma)
  if (!is.numeric(ugamma) || length(dim(ugamma)) != 2L ||
    nrow(ugamma) != ncol(ugamma)) {
    stop(context, " must be a square numeric matrix.", call. = FALSE)
  }
  if (df <= 0L) {
    return(numeric(0L))
  }
  if (df > nrow(ugamma)) {
    stop(context, " dimension is smaller than the model degrees of freedom.",
      call. = FALSE
    )
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
  if (any(lambdas < -tolerance)) {
    stop(context, " has materially negative leading eigenvalues; the fit may ",
      "be numerically unstable or the requested models may not be nested.",
      call. = FALSE
    )
  }
  if (sum(lambdas > tolerance) < df) {
    stop(context, " rank is smaller than the model degrees of freedom.",
      call. = FALSE
    )
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
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  list(
    ug_biased = ugamma_eigenvalues(
      ugamma, df,
      context = "lavaan's single-model UGamma"
    )
  )
}

#' Full nested UGamma reference used to validate the reduced implementation.
#' @keywords internal
ugamma_nested_reference <- function(m0, m1, unbiased = 1) {
  gamma_list <- gamma_matrices(m1, unbiased, m0)
  lapply(
    gamma_list,
    function(gamma_matrix) {
      lav_ugamma_nested_2000(m0, m1, gamma_matrix)
    }
  )
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
#' matrix `C^{-1} D' Gamma D` (see [nested_factor_2000]). This avoids the
#' `q x q` UGamma and its eigendecomposition.
#' @param m0,m1 Two nested `lavaan` objects.
#' @param unbiased Biased (1), unbiased (2), or both (3) gamma.
#' @param df Number of restrictions (degrees-of-freedom difference).
#' @return A list of eigenvalue vectors, one per gamma estimate.
#' @keywords internal
lambdas_nested <- function(m0, m1, unbiased = 1, df) {
  gamma_list <- gamma_matrices(m1, unbiased, m0)
  factors <- nested_factor_2000(m0, m1)
  c_inv <- generalized_inverse(factors$C)
  lapply(gamma_list, function(gamma_matrix) {
    m_phi <- t(factors$D) %*% gamma_matrix %*% factors$D
    ugamma_eigenvalues(
      c_inv %*% m_phi, df,
      context = "nested method-2000 restriction UGamma"
    )
  })
}
