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
         "See `?semTests-support`.", call. = FALSE)
  }
  if (fit@Data@ngroups != 1L) {
    stop(context, " is currently implemented for single-group FIML fits. ",
         "See `?semTests-support`.", call. = FALSE)
  }
  if (isTRUE(fit@Model@conditional.x) || any(fit@Model@nexo > 0L)) {
    stop(context, " is not implemented for models with fixed exogenous ",
         "covariates. See `?semTests-support`.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
fiml_sym <- function(x) {
  (x + t(x)) / 2
}

#' Extract the single matrix returned by a one-group lavaan inspection.
#' @keywords internal
fiml_one_group_matrix <- function(x, what) {
  if (is.list(x)) {
    if (length(x) != 1L) {
      stop("Expected one ", what, " matrix for a single-group FIML fit.",
           call. = FALSE)
    }
    x <- x[[1L]]
  }
  as.matrix(x)
}

#' @keywords internal
fiml_data_matrix <- function(fit) {
  fiml_check_supported(fit)
  x <- as.data.frame(lavaan::lavInspect(fit, "data"))
  ov <- lavaan::lavNames(fit, "ov")
  as.matrix(x[, ov, drop = FALSE])
}

#' Verify that two FIML fits use identical observations.
#' @keywords internal
fiml_check_same_data <- function(m0, m1) {
  X0 <- fiml_data_matrix(m0)
  X1 <- fiml_data_matrix(m1)
  same_values <- identical(dim(X0), dim(X1)) &&
    identical(colnames(X0), colnames(X1)) &&
    identical(is.na(X0), is.na(X1))
  if (same_values) {
    observed <- !is.na(X0)
    same_values <- isTRUE(all.equal(
      unname(X0[observed]), unname(X1[observed]),
      tolerance = 0, check.attributes = FALSE
    ))
  }
  if (!same_values) {
    stop("Nested FIML fits must use the same raw data and missingness mask.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' A basis mapping effective parameters to lavaan's full free parameters.
#'
#' This handles lavaan's explicit equality-constraint basis, simple equality
#' constraints, and equality constraints combined with inequalities or bounds.
#' @keywords internal
fiml_K_matrix <- function(fit) {
  model <- fit@Model
  if (isTRUE(model@eq.constraints)) {
    return(as.matrix(model@eq.constraints.K))
  }
  if (methods::.hasSlot(model, "ceq.simple.only") &&
      isTRUE(model@ceq.simple.only)) {
    K <- as.matrix(model@ceq.simple.K)
    if (!nrow(K)) return(diag(model@nx.free))
    return(qr.Q(qr(K)))
  }
  if (methods::.hasSlot(model, "ceq.JAC") && nrow(model@ceq.JAC) > 0L) {
    return(get_orthogonal_complement(t(as.matrix(model@ceq.JAC))))
  }
  diag(model@nx.free)
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
         call. = FALSE)
  }
  fields <- intersect(
    c("lhs", "op", "rhs", "group", "level", "block"),
    names(pt)
  )
  do.call(paste, c(unname(pt[fields]), sep = "\r"))
}

#' Align an H0 full-parameter basis to H1's parameter ordering.
#' @keywords internal
fiml_align_K0 <- function(m0, m1, K0) {
  key0 <- fiml_parameter_keys(m0)
  key1 <- fiml_parameter_keys(m1)
  if (length(key0) != length(key1) ||
      anyDuplicated(key0) || anyDuplicated(key1) ||
      !setequal(key0, key1)) {
    stop("FIML nested exact restriction requires the same underlying full ",
         "free-parameter set in both models. Use `A.method = \"delta\"` for ",
         "models represented by different parameter sets.", call. = FALSE)
  }
  K0[match(key1, key0), , drop = FALSE]
}

#' @keywords internal
fiml_validate_A <- function(A, df, npar) {
  A <- as.matrix(A)
  if (nrow(A) != df) {
    stop("FIML nested restriction rank (", nrow(A),
         ") does not match df difference (", df, ").", call. = FALSE)
  }
  if (ncol(A) != npar) {
    stop("FIML restriction matrix does not match H1's effective parameter ",
         "space.", call. = FALSE)
  }
  A
}

#' @keywords internal
fiml_A_exact <- function(m0, m1, df) {
  K1 <- fiml_K_matrix(m1)
  K0 <- fiml_align_K0(m0, m1, fiml_K_matrix(m0))
  H <- generalized_inverse(K1) %*% K0
  A <- t(get_orthogonal_complement(H))
  fiml_validate_A(A, df, ncol(K1))
}

#' @keywords internal
fiml_A_delta <- function(m0, m1, df) {
  K1 <- fiml_K_matrix(m1)
  K0 <- fiml_K_matrix(m0)
  Delta1 <- fiml_one_group_matrix(lavaan::lavInspect(m1, "delta"), "delta")
  Delta0 <- fiml_one_group_matrix(lavaan::lavInspect(m0, "delta"), "delta")
  D1 <- Delta1 %*% K1
  D0 <- Delta0 %*% K0
  H <- generalized_inverse(D1) %*% D0
  A <- t(get_orthogonal_complement(H))
  fiml_validate_A(A, df, ncol(K1))
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
  fiml_one_group_matrix(
    lavaan::lavTech(object, "h1.information.observed"),
    "observed H1 information"
  )
}

#' @keywords internal
fiml_gamma_matrix <- function(fit) {
  fiml_one_group_matrix(lavaan::lavTech(fit, "Gamma"), "Gamma")
}

#' @keywords internal
fiml_delta_effective <- function(fit) {
  Delta <- fiml_one_group_matrix(lavaan::lavInspect(fit, "delta"), "delta")
  Delta %*% fiml_K_matrix(fit)
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
  if (df <= 0L) return(list(ug_biased = numeric(0L)))
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
           call. = FALSE)
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
  K1 <- fiml_K_matrix(m1)
  Delta <- fiml_one_group_matrix(lavaan::lavInspect(m1, "delta"), "delta")
  Delta <- Delta %*% K1

  if (fiml.convention == "observed") {
    V <- fiml_h1_information_observed(m1)
    information <- as.matrix(lavaan::lavTech(m1, "information.observed"))
  } else {
    V <- fiml_one_group_matrix(lavaan::lavTech(m1, "WLS.V"), "WLS.V")
    information <- as.matrix(lavaan::lavTech(m1, "information"))
  }

  list(
    A1 = fiml_sym(crossprod(K1, information %*% K1)),
    B1 = fiml_sym(crossprod(V %*% Delta,
                            fiml_gamma_matrix(m1) %*% (V %*% Delta)))
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

  A <- switch(
    A.method,
    delta = fiml_A_delta(m0, m1, df),
    exact = fiml_A_exact(m0, m1, df)
  )
  ingredients <- fiml_nested_ingredients(m1, fiml.convention)
  if (ncol(A) != nrow(ingredients$A1)) {
    stop("FIML restriction matrix and H1 information dimensions do not agree.",
         call. = FALSE)
  }

  Y <- generalized_inverse(ingredients$A1) %*% t(A)
  C <- fiml_sym(A %*% Y)
  S <- fiml_sym(t(Y) %*% ingredients$B1 %*% Y)
  ev <- Re(eigen(generalized_inverse(C) %*% S,
                 only.values = TRUE)$values)
  list(ug_biased = sort(ev, decreasing = TRUE))
}
