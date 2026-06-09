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
fiml_vech_index <- function(p, r, c) {
  if (r < c) {
    tmp <- r
    r <- c
    c <- tmp
  }
  as.integer((c - 1L) * p - (c - 1L) * (c - 2L) / 2L + (r - c + 1L))
}

#' @keywords internal
fiml_sym <- function(x) {
  (x + t(x)) / 2
}

#' @keywords internal
fiml_sym_solve <- function(A, B = diag(nrow(A))) {
  solve(fiml_sym(A), B)
}

#' @keywords internal
fiml_observed_indices <- function(x) {
  which(!is.na(x))
}

#' @keywords internal
fiml_data_matrix <- function(fit) {
  fiml_check_supported(fit)
  x <- as.data.frame(lavaan::lavInspect(fit, "data"))
  ov <- lavaan::lavNames(fit, "ov")
  as.matrix(x[, ov, drop = FALSE])
}

#' @keywords internal
fiml_saturated_scores <- function(X, mu, Sigma) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  pstar <- p * (p + 1L) / 2L
  scores <- matrix(0, n, p + pstar)

  for (row in seq_len(n)) {
    obs <- fiml_observed_indices(X[row, ])
    if (!length(obs)) stop("FIML row has no observed variables.", call. = FALSE)

    x_o <- X[row, obs]
    Sigma_o <- Sigma[obs, obs, drop = FALSE]
    d <- x_o - mu[obs]
    A <- solve(Sigma_o)
    z <- drop(A %*% d)
    G <- fiml_sym(A - tcrossprod(z))

    scores[row, obs] <- -2 * z

    q <- length(obs)
    for (cj in seq_len(q)) {
      full_c <- obs[cj]
      for (ri in cj:q) {
        full_r0 <- obs[ri]
        full_r <- max(full_r0, full_c)
        full_l <- min(full_r0, full_c)
        idx <- p + fiml_vech_index(p, full_r, full_l)
        scores[row, idx] <- scores[row, idx] +
          if (ri == cj) G[ri, cj] else 2 * G[ri, cj]
      }
    }
  }

  scores
}

#' @keywords internal
fiml_basis_trace_adad <- function(A, left, right) {
  lu <- c(left$a, left$b)
  lv <- c(left$b, left$a)
  ru <- c(right$a, right$b)
  rv <- c(right$b, right$a)
  ln <- if (left$a == left$b) 1L else 2L
  rn <- if (right$a == right$b) 1L else 2L

  out <- 0
  for (i in seq_len(ln)) {
    for (j in seq_len(rn)) {
      out <- out + A[rv[j], lu[i]] * A[lv[i], ru[j]]
    }
  }
  out
}

#' @keywords internal
fiml_saturated_hessian_analytic <- function(X, mu, Sigma) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  pstar <- p * (p + 1L) / 2L
  H <- matrix(0, p + pstar, p + pstar)

  for (row in seq_len(n)) {
    obs <- fiml_observed_indices(X[row, ])
    if (!length(obs)) stop("FIML row has no observed variables.", call. = FALSE)

    x_o <- X[row, obs]
    Sigma_o <- Sigma[obs, obs, drop = FALSE]
    d <- x_o - mu[obs]
    A <- solve(Sigma_o)
    z <- drop(A %*% d)
    q <- length(obs)

    H[obs, obs] <- H[obs, obs] + 2 * A

    basis <- vector("list", q * (q + 1L) / 2L)
    kk <- 0L
    for (cj in seq_len(q)) {
      full_c <- obs[cj]
      for (ri in cj:q) {
        full_r0 <- obs[ri]
        full_r <- max(full_r0, full_c)
        full_l <- min(full_r0, full_c)
        Dz <- numeric(q)
        if (ri == cj) {
          Dz[ri] <- z[ri]
        } else {
          Dz[ri] <- z[cj]
          Dz[cj] <- z[ri]
        }
        ADz <- drop(A %*% Dz)
        full_index <- p + fiml_vech_index(p, full_r, full_l)

        h <- 2 * ADz
        H[obs, full_index] <- H[obs, full_index] + h
        H[full_index, obs] <- H[full_index, obs] + h

        kk <- kk + 1L
        basis[[kk]] <- list(
          full_index = full_index,
          a = ri,
          b = cj,
          Dz = Dz,
          ADz = ADz
        )
      }
    }

    for (k in seq_along(basis)) {
      bk <- basis[[k]]
      for (l in seq_len(k)) {
        bl <- basis[[l]]
        h <- -fiml_basis_trace_adad(A, bl, bk) +
          sum(bl$Dz * bk$ADz) + sum(bk$Dz * bl$ADz)
        H[bk$full_index, bl$full_index] <-
          H[bk$full_index, bl$full_index] + h
        if (k != l) {
          H[bl$full_index, bk$full_index] <-
            H[bl$full_index, bk$full_index] + h
        }
      }
    }
  }

  fiml_sym(H / n)
}

#' Saturated observed-data FIML eta-space pieces.
#' @keywords internal
fiml_saturated_moments <- function(fit) {
  X <- fiml_data_matrix(fit)
  ov <- lavaan::lavNames(fit, "ov")
  ss <- lavaan::lavInspect(fit, "sampstat")
  mu <- as.numeric(ss$mean[ov])
  Sigma <- as.matrix(ss$cov[ov, ov])

  scores_dev <- fiml_saturated_scores(X, mu, Sigma)
  H_dev_mean <- fiml_saturated_hessian_analytic(X, mu, Sigma)
  n <- nrow(X)
  H <- (n / 2) * H_dev_mean
  J <- 0.25 * crossprod(scores_dev)
  Hinv <- fiml_sym_solve(H)
  Gamma <- Hinv %*% J %*% Hinv

  list(
    mean = mu,
    cov = Sigma,
    H = fiml_sym(H),
    J = fiml_sym(J),
    Gamma = fiml_sym(Gamma),
    data = X
  )
}

#' FIML goodness-of-fit eigenvalues in saturated eta-space.
#' @keywords internal
fiml_lambdas <- function(fit, df) {
  if (df <= 0L) return(list(ug_biased = numeric(0L)))
  sm <- fiml_saturated_moments(fit)
  H <- sm$H
  Gamma <- sm$Gamma
  Delta <- as.matrix(lavaan::lavInspect(fit, "delta"))
  if (nrow(Delta) != nrow(H)) {
    stop("FIML delta rows do not match saturated moment dimension.",
         call. = FALSE)
  }

  HD <- H %*% Delta
  P <- fiml_sym(crossprod(Delta, HD))
  U <- H - HD %*% generalized_inverse(P) %*% t(HD)
  U <- fiml_sym(U)

  Gs <- fiml_sym(Gamma)
  reduced <- tryCatch({
    R <- chol(Gs)
    R %*% U %*% t(R)
  }, error = function(e) {
    eg <- eigen(Gs, symmetric = TRUE)
    sq <- eg$vectors %*% diag(sqrt(pmax(eg$values, 0)), nrow(Gs)) %*%
      t(eg$vectors)
    sq %*% U %*% sq
  })
  reduced <- fiml_sym(reduced)
  ev <- eigen(reduced, symmetric = TRUE, only.values = TRUE)$values
  lambdas <- sort(Re(ev), decreasing = TRUE)[seq_len(df)]
  list(ug_biased = lambdas)
}

#' @keywords internal
fiml_K_matrix <- function(fit) {
  model <- fit@Model
  if (isTRUE(model@eq.constraints)) {
    return(model@eq.constraints.K)
  }
  if (methods::.hasSlot(model, "ceq.simple.only") &&
      isTRUE(model@ceq.simple.only)) {
    return(t(model@ceq.simple.K))
  }
  diag(model@nx.free)
}

#' @keywords internal
fiml_A_exact <- function(m0, m1, df) {
  K1 <- fiml_K_matrix(m1)
  K0 <- fiml_K_matrix(m0)
  if (nrow(K1) != nrow(K0)) {
    stop("FIML nested exact restriction requires matching unconstrained ",
         "parameter spaces.", call. = FALSE)
  }
  H <- generalized_inverse(K1) %*% K0
  A <- t(get_orthogonal_complement(H))
  if (nrow(A) != df) {
    stop("FIML nested exact restriction rank (", nrow(A),
         ") does not match df difference (", df, ").", call. = FALSE)
  }
  A
}

#' @keywords internal
fiml_A_delta <- function(m0, m1, df) {
  A <- get_a_matrix(m1, m0)
  if (m1@Model@eq.constraints) {
    A <- A %*% t(m1@Model@eq.constraints.K)
  }
  if (nrow(A) != df) {
    stop("FIML nested delta restriction rank (", nrow(A),
         ") does not match df difference (", df, ").", call. = FALSE)
  }
  A
}

#' FIML Satorra-2000 nested restriction eigenvalues.
#' @keywords internal
fiml_lambdas_nested <- function(m0, m1, df,
                                A.method = c("exact", "delta")) {
  A.method <- match.arg(A.method)
  fiml_check_supported(m0, "Nested FIML FMG")
  fiml_check_supported(m1, "Nested FIML FMG")

  X0 <- fiml_data_matrix(m0)
  X1 <- fiml_data_matrix(m1)
  if (!identical(dim(X0), dim(X1)) ||
      !identical(colnames(X0), colnames(X1)) ||
      !identical(is.na(X0), is.na(X1)) ||
      any(abs(X0[!is.na(X0)] - X1[!is.na(X1)]) > 0)) {
    stop("Nested FIML fits must use the same raw data and missingness mask.",
         call. = FALSE)
  }

  sm <- fiml_saturated_moments(m1)
  H <- sm$H
  Gamma <- sm$Gamma
  Delta1 <- as.matrix(lavaan::lavInspect(m1, "delta"))
  if (nrow(Delta1) != nrow(H)) {
    stop("FIML H1 delta rows do not match saturated moment dimension.",
         call. = FALSE)
  }

  A <- switch(A.method,
              exact = fiml_A_exact(m0, m1, df),
              delta = fiml_A_delta(m0, m1, df))
  if (ncol(A) != ncol(Delta1)) {
    stop("FIML restriction matrix does not match H1 free parameter space.",
         call. = FALSE)
  }

  HD <- H %*% Delta1
  P <- fiml_sym(crossprod(Delta1, HD))
  P_inv <- generalized_inverse(P)
  Y <- P_inv %*% t(A)
  C <- fiml_sym(A %*% Y)
  mid <- t(HD) %*% Gamma %*% HD
  S <- fiml_sym(t(Y) %*% mid %*% Y)
  ev <- Re(eigen(generalized_inverse(C) %*% S, only.values = TRUE)$values)
  list(ug_biased = sort(ev, decreasing = TRUE))
}
