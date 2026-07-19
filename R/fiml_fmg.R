#' Is this a continuous FIML/missing-data lavaan fit?
#' @keywords internal
is_fiml <- function(fit) {
  !isTRUE(fit@Model@categorical) &&
    fit@Options$missing %in% c("ml", "fiml", "ml.x", "fiml.x", "direct")
}

#' @keywords internal
fiml_check_supported <- function(fit, context = "FIML FMG") {
  if (isTRUE(fit@Model@categorical)) {
    stop(context, " is only implemented for continuous lavaan fits. ",
      "See `?semTests-support`.",
      call. = FALSE
    )
  }
  if (uses_fixed_or_conditional_x(fit)) {
    stop(context, " is not implemented for models with fixed or conditional ",
      "observed exogenous covariates. Refit with `fixed.x = FALSE` and ",
      "`conditional.x = FALSE` for joint random-x inference, or see ",
      "`?semTests-support`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @keywords internal
fiml_sym <- function(x) {
  (x + t(x)) / 2
}

# Convert a lavaan inspection to one matrix per group.
fiml_group_matrices <- function(x, fit, what) {
  blocks <- if (is.list(x)) x else list(x)
  n_groups <- fit@Data@ngroups
  if (length(blocks) != n_groups) {
    stop("Expected ", n_groups, " ", what, " ",
      if (n_groups == 1L) "matrix" else "matrices",
      " for a ", n_groups, "-group FIML fit.",
      call. = FALSE
    )
  }
  lapply(blocks, as.matrix)
}

# Group proportions used by lavaan's global moment-space convention.
fiml_group_weights <- function(fit) {
  weights <- as.numeric(unlist(fit@SampleStats@nobs)) /
    fit@SampleStats@ntotal
  if (length(weights) != fit@Data@ngroups ||
    any(!is.finite(weights)) || any(weights <= 0) ||
    abs(sum(weights) - 1) > sqrt(.Machine$double.eps)) {
    stop("Could not construct valid FIML group weights.", call. = FALSE)
  }
  weights
}

# Assemble independent group blocks in the global moment space. Projector
# weights use n_g/N; Gamma uses its reciprocal so U Gamma keeps the intended
# group scaling.
fiml_group_bdiag <- function(x, fit, what,
                             scale = c("weight", "inverse_weight")) {
  scale <- match.arg(scale)
  blocks <- fiml_group_matrices(x, fit, what)
  weights <- fiml_group_weights(fit)
  multipliers <- switch(scale,
    weight = weights,
    inverse_weight = 1 / weights
  )
  blocks <- Map(function(block, multiplier) {
    block * multiplier
  }, blocks, multipliers)
  lavaan::lav_matrix_bdiag(blocks)
}

#' @keywords internal
fiml_data_blocks <- function(fit) {
  fiml_check_supported(fit)
  blocks <- fiml_group_matrices(
    lavaan::lavInspect(fit, "data"), fit, "raw-data"
  )
  ov <- lavaan::lavNames(fit, "ov")
  lapply(blocks, function(x) {
    missing_names <- setdiff(ov, colnames(x))
    if (length(missing_names)) {
      stop("FIML raw data do not contain all observed variables.",
        call. = FALSE
      )
    }
    x[, ov, drop = FALSE]
  })
}

#' Verify that two FIML fits use identical observations.
#' @keywords internal
fiml_check_same_data <- function(m0, m1) {
  X0 <- fiml_data_blocks(m0)
  X1 <- fiml_data_blocks(m1)
  same_values <- identical(m0@Data@group.label, m1@Data@group.label) &&
    length(X0) == length(X1)
  if (same_values) {
    same_values <- all(vapply(seq_along(X0), function(group) {
      x0 <- X0[[group]]
      x1 <- X1[[group]]
      same_block <- identical(dim(x0), dim(x1)) &&
        identical(colnames(x0), colnames(x1)) &&
        identical(is.na(x0), is.na(x1))
      if (!same_block) {
        return(FALSE)
      }
      observed <- !is.na(x0)
      isTRUE(all.equal(
        unname(x0[observed]), unname(x1[observed]),
        tolerance = 0, check.attributes = FALSE
      ))
    }, logical(1)))
  }
  if (!same_values) {
    stop("Nested FIML fits must use the same raw data and missingness mask.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stable identities for lavaan's full free parameters.
#' @keywords internal
fiml_parameter_keys <- function(fit) {
  pt <- lavaan::parTable(fit)
  pt <- pt[pt$free > 0L, , drop = FALSE]
  pt <- pt[order(pt$free), , drop = FALSE]
  if (nrow(pt) != fit@Model@nx.free ||
    !identical(as.integer(pt$free), seq_len(fit@Model@nx.free))) {
    stop("Could not identify lavaan's full free-parameter ordering.",
      call. = FALSE
    )
  }
  fields <- intersect(
    c("lhs", "op", "rhs", "group", "level", "block"),
    names(pt)
  )
  do.call(paste, c(unname(pt[fields]), sep = "\r"))
}

#' Align an H0 parameter basis to H1's full-parameter ordering.
#' @keywords internal
fiml_align_parameter_basis <- function(m0, m1, basis0) {
  key0 <- fiml_parameter_keys(m0)
  key1 <- fiml_parameter_keys(m1)
  if (length(key0) != length(key1) ||
    anyDuplicated(key0) || anyDuplicated(key1) ||
    !setequal(key0, key1)) {
    stop("FIML nested exact restriction requires the same underlying full ",
      "free-parameter set in both models. Use `A.method = \"delta\"` for ",
      "models represented by different parameter sets.",
      call. = FALSE
    )
  }
  basis0[match(key1, key0), , drop = FALSE]
}

#' @keywords internal
fiml_validate_restriction <- function(restriction, df, npar) {
  restriction <- as.matrix(restriction)
  if (nrow(restriction) != df) {
    stop("FIML nested restriction rank (", nrow(restriction),
      ") does not match df difference (", df, ").",
      call. = FALSE
    )
  }
  if (ncol(restriction) != npar) {
    stop("FIML restriction matrix does not match H1's effective parameter ",
      "space.",
      call. = FALSE
    )
  }
  restriction
}

#' @keywords internal
fiml_restriction_exact <- function(m0, m1, df) {
  basis1 <- model_parameter_basis(m1)
  basis0 <- fiml_align_parameter_basis(
    m0,
    m1,
    model_parameter_basis(m0)
  )
  map <- generalized_inverse(basis1) %*% basis0
  restriction <- t(get_orthogonal_complement(map))
  fiml_validate_restriction(restriction, df, ncol(basis1))
}

#' @keywords internal
fiml_restriction_delta <- function(m0, m1, df) {
  delta1 <- model_effective_delta(m1)
  delta0 <- model_effective_delta(m0)
  map <- generalized_inverse(delta1) %*% delta0
  restriction <- t(get_orthogonal_complement(map))
  fiml_validate_restriction(
    restriction,
    df,
    ncol(model_parameter_basis(m1))
  )
}

#' Saturated observed-data FIML information in moment space.
#'
#' lavaan 0.7-2 exposes the corrected FIML H1 information. Setting
#' `h1.information` on a local copy requests the unstructured, observed
#' saturated information without refitting or changing the fitted parameters.
#' @keywords internal
fiml_h1_information_observed <- function(fit) {
  object <- fit
  object@Options$h1.information[] <- "unstructured"
  fiml_group_bdiag(
    lavaan::lavTech(object, "h1.information.observed"),
    fit, "observed H1 information",
    scale = "weight"
  )
}

#' @keywords internal
fiml_gamma_matrix <- function(fit) {
  fiml_group_bdiag(
    lavaan::lavTech(fit, "Gamma"),
    fit, "Gamma",
    scale = "inverse_weight"
  )
}

#' @keywords internal
fiml_delta_effective <- function(fit) {
  model_effective_delta(fit)
}

# Expected H1 weight in lavaan's grouped moment-space convention.
fiml_h1_information_lavaan <- function(fit) {
  fiml_group_bdiag(
    lavaan::lavTech(fit, "WLS.V"),
    fit, "WLS.V",
    scale = "weight"
  )
}

#' Return the leading symmetric-product eigenvalues without forming AB.
#' @keywords internal
fiml_sandwich_eigenvalues <- function(U, Gamma, df) {
  Gamma <- fiml_sym(Gamma)
  eg <- eigen(Gamma, symmetric = TRUE)
  root <- eg$vectors %*%
    diag(sqrt(pmax(eg$values, 0)), nrow(Gamma)) %*%
    t(eg$vectors)
  reduced <- fiml_sym(root %*% fiml_sym(U) %*% root)
  ev <- eigen(reduced, symmetric = TRUE, only.values = TRUE)$values
  sort(Re(ev), decreasing = TRUE)[seq_len(df)]
}

#' FIML goodness-of-fit eigenvalues.
#'
#' The observed convention uses observed saturated information, matching the
#' convention independently implemented and validated in magmaan. The lavaan
#' convention returns lavaan's own inspected `UGamma` spectrum.
#' @keywords internal
fiml_lambdas <- function(fit, df,
                         fiml.convention = c("observed", "lavaan")) {
  fiml.convention <- match.arg(fiml.convention)
  if (df <= 0L) {
    return(list(ug_biased = numeric(0L)))
  }
  fiml_check_supported(fit)

  lambdas <- if (fiml.convention == "lavaan") {
    ugamma_eigenvalues(
      lavaan::lavInspect(fit, "UGamma"), df,
      context = "lavaan's FIML UGamma"
    )
  } else {
    H <- fiml_h1_information_observed(fit)
    Gamma <- fiml_gamma_matrix(fit)
    Delta <- fiml_delta_effective(fit)
    if (nrow(Delta) != nrow(H) || !identical(dim(H), dim(Gamma))) {
      stop("FIML H1 information, Gamma, and delta dimensions do not agree.",
        call. = FALSE
      )
    }
    HD <- H %*% Delta
    P <- fiml_sym(crossprod(Delta, HD))
    U <- H - HD %*% generalized_inverse(P) %*% t(HD)
    fiml_sandwich_eigenvalues(U, Gamma, df)
  }

  list(ug_biased = lambdas)
}

#' Gamma-free nested FIML ingredients in H1's effective parameter space.
#' @keywords internal
fiml_nested_ingredients <- function(m1,
                                    fiml.convention = c("observed", "lavaan")) {
  fiml.convention <- match.arg(fiml.convention)
  basis1 <- model_parameter_basis(m1)
  delta <- model_effective_delta(m1)

  if (fiml.convention == "observed") {
    V <- fiml_h1_information_observed(m1)
    information <- as.matrix(lavaan::lavTech(m1, "information.observed"))
  } else {
    V <- fiml_h1_information_lavaan(m1)
    information <- as.matrix(lavaan::lavTech(m1, "information"))
  }

  list(
    A1 = fiml_sym(crossprod(basis1, information %*% basis1)),
    B1 = fiml_sym(crossprod(
      V %*% delta,
      fiml_gamma_matrix(m1) %*% (V %*% delta)
    ))
  )
}

#' FIML Satorra-2000 nested restriction eigenvalues.
#' @keywords internal
fiml_lambdas_nested <- function(m0, m1, df,
                                A.method = c("delta", "exact"),
                                fiml.convention = c("observed", "lavaan")) {
  A.method <- match.arg(A.method)
  fiml.convention <- match.arg(fiml.convention)
  fiml_check_supported(m0, "Nested FIML FMG")
  fiml_check_supported(m1, "Nested FIML FMG")
  fiml_check_same_data(m0, m1)

  A <- switch(A.method,
    delta = fiml_restriction_delta(m0, m1, df),
    exact = fiml_restriction_exact(m0, m1, df)
  )
  ingredients <- fiml_nested_ingredients(m1, fiml.convention)
  if (ncol(A) != nrow(ingredients$A1)) {
    stop("FIML restriction matrix and H1 information dimensions do not agree.",
      call. = FALSE
    )
  }

  Y <- generalized_inverse(ingredients$A1) %*% t(A)
  C <- fiml_sym(A %*% Y)
  S <- fiml_sym(t(Y) %*% ingredients$B1 %*% Y)
  ev <- Re(eigen(generalized_inverse(C) %*% S,
    only.values = TRUE
  )$values)
  list(ug_biased = sort(ev, decreasing = TRUE))
}
