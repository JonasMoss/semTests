# These functions are copied from lavaan: https://github.com/yrosseel/lavaan


get_orthogonal_complement <- function(mat) {
  qr_decomp <- qr(mat)
  q_mat <- qr.Q(qr_decomp, complete = TRUE)
  q_mat[, -seq_len(qr_decomp$rank), drop = FALSE]
}

get_i_minus_b_inv <- function(MLIST) {
  beta <- MLIST$beta
  nr <- nrow(MLIST$psi)

  if (!is.null(beta)) {
    tmp <- -beta
    tmp[lav_matrix_diag_idx(nr)] <- 1
    return(solve(tmp))
  }

  diag(nr)

}

lav_matrix_vech_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  if (n < 100L) {
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n), n, n, byrow = TRUE)
    if (diagonal) which(ROW >= COL) else which(ROW > COL)
  } else {
    if (diagonal) {
      unlist(lapply(
        seq_len(n),
        function(j) (j - 1L) * n + seq.int(j, n)
      ))
    } else {
      unlist(lapply(
        seq_len(n - 1L),
        function(j) (j - 1L) * n + seq.int(j + 1L, n)
      ))
    }
  }
}

lav_matrix_vechru_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  if (n < 100L) {
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n), n, n, byrow = TRUE)
    tmp <- matrix(seq_len(n * n), n, n, byrow = TRUE)
    if (diagonal) tmp[ROW >= COL] else tmp[ROW > COL]
  } else {
    if (diagonal) {
      unlist(lapply(
        seq_len(n),
        function(j) seq.int(j - 1, n - 1) * n + j
      ))
    } else {
      unlist(lapply(
        seq_len(n - 1L),
        function(j) seq.int(j, n - 1) * n + j
      ))
    }
  }
}

lav_matrix_diag_idx <- function(n = 1L) {
  n <- as.integer(n)
  if (n < 1L) return(integer(0L))
  1L + (seq_len(n) - 1L) * (n + 1L)
}

lav_matrix_diagh_idx <- function(n = 1L) {
  n <- as.integer(n)
  if (n < 1L) return(integer(0L))
  if (n == 1L) return(1L)
  c(1L, cumsum(n:2L) + 1L)
}

lav_matrix_vec <- function(x) as.vector(x)

lav_matrix_vech <- function(S, diagonal = TRUE) {
  ROW <- row(S)
  COL <- col(S)
  if (diagonal) S[ROW >= COL] else S[ROW > COL]
}

lav_matrix_vechu_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  # if (lav_use_lavaanC() && n > 1L) {
  #   return(lavaanC::m_vechu_idx(n, diagonal))
  # }
  if (n < 100L) {
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n), n, n, byrow = TRUE)
    if (diagonal) which(ROW <= COL) else which(ROW < COL)
  } else {
    if (diagonal) {
      unlist(lapply(seq_len(n), function(j) seq.int(j) + (j - 1) * n))
    } else {
      unlist(lapply(seq_len(n - 1L), function(j) seq.int(j) + j * n))
    }
  }
}

lav_func_jacobian_complex <- function(func, x,
                                      h = .Machine$double.eps, ...,
                                      fallback.simple = TRUE) {
  f0 <- try(func(x * (0 + 1i), ...), silent = TRUE)
  if (!is.complex(f0)) {
    if (fallback.simple) {
      dx <- lav_func_jacobian_simple(func = func, x = x, h = sqrt(h), ...)
      return(dx)
    } else {
      stop("function does not return a complex value")
    }
  }
  if (inherits(f0, "try-error")) {
    if (fallback.simple) {
      dx <- lav_func_jacobian_simple(func = func, x = x, h = sqrt(h), ...)
      return(dx)
    } else {
      stop("Function does not support complex arguments.")
    }
  }
  nres <- length(f0)
  nvar <- length(x)

  h <- pmax(h, abs(h * x))

  tmp <- x + h
  h <- (tmp - x)

  dx <- matrix(as.numeric(NA), nres, nvar)
  for (p in seq_len(nvar)) {
    dx[, p] <- Im(func(x + h * 1i * (seq.int(nvar) == p), ...)) / h[p]
  }

  dx
}

lav_matrix_vech_reverse <- lav_matrix_vechru_reverse <- lav_matrix_upper2full <-
  function(x, diagonal = TRUE) {
    if (diagonal) {
      p <- (sqrt(1 + 8 * length(x)) - 1) / 2
    } else {
      p <- (sqrt(1 + 8 * length(x)) + 1) / 2
    }

    S <- numeric(p * p)
    S[lav_matrix_vech_idx(p, diagonal = diagonal)] <- x
    S[lav_matrix_vechru_idx(p, diagonal = diagonal)] <- x

    attr(S, "dim") <- c(p, p)
    S
  }


computeSigmaHat.LISREL <- function(MLIST = NULL,
                                                        delta = TRUE) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  PSI <- MLIST$psi
  THETA <- MLIST$theta
  BETA <- MLIST$beta
  WMAT <- MLIST$wmat

  # standard: no composites
  if (is.null(WMAT)) {
    # beta?
    if (is.null(BETA)) {
      LAMBDA..IB.inv <- LAMBDA
    } else {
      IB.inv <- get_i_minus_b_inv(MLIST = MLIST)
      LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    }

    # compute V(Y*|x_i)
    VYx <- tcrossprod(LAMBDA..IB.inv %*% PSI, LAMBDA..IB.inv) + THETA

    # composites, or mix of composites and latent variables
  } else {
    # - first join LAMBDA and WMAT
    # - create 'T' matrix: - identity for regular lv's,
    #                      - THETA block-diagonal for composites
    # - create C_0: VETA, but zero diagonal elements for composites
    cov.idx <- which(apply(LAMBDA, 1L,
                           function(x) sum(x == 0) == ncol(LAMBDA)))
    clv.idx <- which(apply(LAMBDA, 2L,
                           function(x) sum(x == 0) == nrow(LAMBDA)))
    # regular latent variables
    rlv.idx <- seq_len(ncol(LAMBDA))[-clv.idx]

    # combine LAMBDA and WMAT
    LW <- LAMBDA + WMAT

    Tmat <- diag(nrow(LAMBDA))
    Tmat[cov.idx, cov.idx] <- THETA[cov.idx, cov.idx]
    wtw <- t(LW[,clv.idx, drop = FALSE]) %*% Tmat %*% LW[,clv.idx, drop = FALSE]
    wtw.inv <- solve(wtw)
    WTW.inv <- diag(ncol(LAMBDA))
    WTW.inv[clv.idx, clv.idx] <- wtw.inv

    if (is.null(BETA)) {
      IB.inv <- diag(nrow(PSI))
    } else {
      IB.inv <- get_i_minus_b_inv(MLIST = MLIST)
    }
    VETA <- IB.inv %*% PSI %*% t(IB.inv)
    C0 <- VETA; diag(C0)[clv.idx] <- 0

    VYx <- Tmat %*% LW %*% WTW.inv %*% C0 %*% t(WTW.inv) %*% t(LW) %*% Tmat + THETA
  }

  # if delta, scale
  if (delta && !is.null(MLIST$delta)) {
    DELTA <- diag(MLIST$delta[, 1L], nrow = nvar, ncol = nvar)
    VYx <- DELTA %*% VYx %*% DELTA
  }

  VYx
}

lav_func_jacobian_simple <- function(func, x,
                                     h = sqrt(.Machine$double.eps), ...) {
  f0 <- func(x, ...)
  nres <- length(f0)
  nvar <- length(x)

  # determine 'h' per element of x
  h <- pmax(h, abs(h * x))

  # get exact h, per x
  tmp <- x + h
  h <- (tmp - x)

  # simple 'forward' method
  dx <- matrix(as.numeric(NA), nres, nvar)
  for (p in seq_len(nvar)) {
    dx[, p] <- (func(x + h * (seq.int(nvar) == p), ...) - func(x, ...)) / h[p]
  }

  dx
}
