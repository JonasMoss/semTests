# These functions are copied from lavaan: https://github.com/yrosseel/lavaan


derivative.sigma.LISREL <- function(m = "lambda",
                                    idx = seq_len(length(MLIST[[m]])),
                                    MLIST = NULL,
                                    vech = TRUE,
                                    delta = TRUE) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  PSI <- MLIST$psi
  WMAT <- MLIST$wmat

  # for composites (vec version)
  compute.sigma <- function(x, mm = "wmat", MLIST = NULL) {
    mlist <- MLIST
    if (mm %in% c("psi", "theta")) {
      mlist[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist[[mm]][, ] <- x
    }
    lav_matrix_vec(computeSigmaHat.LISREL(mlist))
  }

  composites <- FALSE
  if (!is.null(WMAT)) {
    composites <- TRUE
  }

  # only lower.tri part of sigma (not same order as elimination matrix?)
  v.idx <- lav_matrix_vech_idx(nvar)
  pstar <- nvar * (nvar + 1) / 2

  # shortcut for gamma, nu, alpha, tau,.... : empty matrix
  if (m == "nu" || m == "alpha" || m == "tau" || m == "gamma" ||
      m == "gw" || m == "cov.x" || m == "mean.x") {
    return(matrix(0.0, nrow = pstar, ncol = length(idx)))
  }

  # Delta?
  delta.flag <- FALSE
  if (delta && !is.null(MLIST$delta)) {
    DELTA <- MLIST$delta
    delta.flag <- TRUE
  } else if (m == "delta") { # modindices?
    return(matrix(0.0, nrow = pstar, ncol = length(idx)))
  }

  # beta?
  if (!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  } else {
    IB.inv <- get_i_minus_b_inv(MLIST = MLIST)
  }

  # pre
  # if(m == "lambda" || m == "beta")
  #    IK <- diag(nvar*nvar) + lav_matrix_commutation(nvar, nvar)
  if (m == "lambda" || m == "beta") {
    L1 <- LAMBDA %*% IB.inv %*% PSI %*% t(IB.inv)
  }
  if (m == "beta" || m == "psi") {
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }

  # here we go:
  if (m == "lambda") {
    KOL.idx <- matrix(1:(nvar * nfac), nvar, nfac, byrow = TRUE)[idx]
    DX <- (L1 %x% diag(nvar))[, idx, drop = FALSE] +
      (diag(nvar) %x% L1)[, KOL.idx, drop = FALSE]
  } else if (m == "beta") {
    if (composites) {
      DX <- lav_func_jacobian_complex(func = compute.sigma,
                                      x = lav_matrix_vec(MLIST$beta),
                                      mm = "beta", MLIST = MLIST)
      DX <- DX[, idx, drop = FALSE]
    } else {
      KOL.idx <- matrix(1:(nfac * nfac), nfac, nfac, byrow = TRUE)[idx]
      DX <- (L1 %x% LAMBDA..IB.inv)[, idx, drop = FALSE] +
        (LAMBDA..IB.inv %x% L1)[, KOL.idx, drop = FALSE]
      # this is not really needed (because we select idx=m.el.idx)
      # but just in case we need all elements of beta...
      DX[, which(idx %in% lav_matrix_diag_idx(nfac))] <- 0.0
    }
  } else if (m == "psi") {
    if (composites) {
      tmp <- lav_func_jacobian_complex(func = compute.sigma,
                                       x = lav_matrix_vech(MLIST$psi),
                                       mm = "psi", MLIST = MLIST)
      DX <- matrix(0, nrow = nrow(tmp), ncol = length(PSI))
      DX[, lav_matrix_vech_idx(nrow(PSI))] <- tmp
      DX[, lav_matrix_vechu_idx(nrow(PSI), diagonal = FALSE)] <-
        DX[, lav_matrix_vech_idx(nrow(PSI), diagonal = FALSE), drop = FALSE]
      DX <- DX[, idx, drop = FALSE]
    } else {
      DX <- (LAMBDA..IB.inv %x% LAMBDA..IB.inv)
      # symmetry correction, but keeping all duplicated elements
      # since we depend on idx=m.el.idx
      lower.idx <- lav_matrix_vech_idx(nfac, diagonal = FALSE)
      upper.idx <- lav_matrix_vechru_idx(nfac, diagonal = FALSE)
      offdiagSum <- DX[, lower.idx] + DX[, upper.idx]
      DX[, c(lower.idx, upper.idx)] <- cbind(offdiagSum, offdiagSum)
      DX <- DX[, idx, drop = FALSE]
    }
  } else if (m == "theta") {
    # DX <- diag(nvar*nvar) # very sparse...
    DX <- matrix(0, nvar * nvar, length(idx))
    DX[cbind(idx, seq_along(idx))] <- 1
    # symmetry correction not needed, since all off-diagonal elements
    # are zero?
  } else if (m == "delta") {
    Omega <- computeSigmaHat.LISREL(MLIST, delta = FALSE)
    DD <- diag(DELTA[, 1], nvar, nvar)
    DD.Omega <- (DD %*% Omega)
    A <- DD.Omega %x% diag(nvar)
    B <- diag(nvar) %x% DD.Omega
    DX <- A[, lav_matrix_diag_idx(nvar), drop = FALSE] +
      B[, lav_matrix_diag_idx(nvar), drop = FALSE]
    DX <- DX[, idx, drop = FALSE]
  } else if (m == "wmat") {
    # just a dummy to get us going
    DX <- lav_func_jacobian_complex(func = compute.sigma,
                                    x = lav_matrix_vec(WMAT),
                                    mm = "wmat", MLIST = MLIST)
    DX <- DX[, idx, drop = FALSE]

  } else {
    stop()
  }

  if (delta.flag && !m == "delta") {
    DX <- DX * as.vector(DELTA %x% DELTA)
  }

  # vech?
  if (vech) {
    DX <- DX[v.idx, , drop = FALSE]
  }

  DX
}

derivative.mu.LISREL <- function(m = "alpha",
                                 # all model matrix elements, or only a few?
                                 idx = seq_len(length(MLIST[[m]])),
                                 MLIST = NULL) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  WMAT <- MLIST$wmat

  # shortcut for empty matrices
  if (m == "gamma" || m == "psi" || m == "theta" || m == "tau" ||
      m == "delta" || m == "gw" || m == "cov.x" || m == "mean.x") {
    return(matrix(0.0, nrow = nvar, ncol = length(idx)))
  }

  # missing alpha
  if (is.null(MLIST$alpha)) {
    ALPHA <- matrix(0, nfac, 1L)
  } else {
    ALPHA <- MLIST$alpha
  }


  # beta?
  if (!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  } else {
    IB.inv <- get_i_minus_b_inv(MLIST = MLIST)
  }

  if (m == "nu") {
    DX <- diag(nvar)
  } else if (m == "lambda") {
    DX <- t(IB.inv %*% ALPHA) %x% diag(nvar)
  } else if (m == "wmat") {
    # dummy, just to get us going
    DX <- t(IB.inv %*% ALPHA) %x% diag(nvar)
  } else if (m == "beta") {
    DX <- t(IB.inv %*% ALPHA) %x% (LAMBDA %*% IB.inv)
    # this is not really needed (because we select idx=m.el.idx)
    DX[, lav_matrix_diag_idx(nfac)] <- 0.0
  } else if (m == "alpha") {
    DX <- LAMBDA %*% IB.inv
  } else {
    stop()
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}

derivative.th.LISREL <- function(m = "tau",
                                 # all model matrix elements, or only a few?
                                 idx = seq_len(length(MLIST[[m]])),
                                 th.idx = NULL,
                                 MLIST = NULL,
                                 delta = TRUE) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  TAU <- MLIST$tau
  nth <- nrow(TAU)

  # missing alpha
  if (is.null(MLIST$alpha)) {
    ALPHA <- matrix(0, nfac, 1L)
  } else {
    ALPHA <- MLIST$alpha
  }

  # missing nu
  if (is.null(MLIST$nu)) {
    NU <- matrix(0, nvar, 1L)
  } else {
    NU <- MLIST$nu
  }

  # Delta?
  delta.flag <- FALSE
  if (delta && !is.null(MLIST$delta)) {
    DELTA <- MLIST$delta
    delta.flag <- TRUE
  }

  if (is.null(th.idx)) {
    th.idx <- seq_len(nth)
    nlev <- rep(1L, nvar)
    K_nu <- diag(nvar)
  } else {
    nlev <- tabulate(th.idx, nbins = nvar)
    nlev[nlev == 0L] <- 1L
    K_nu <- matrix(0, sum(nlev), nvar)
    K_nu[cbind(seq_len(sum(nlev)), rep(seq_len(nvar), times = nlev))] <- 1.0
  }

  # shortcut for empty matrices
  if (m == "gamma" || m == "psi" || m == "theta" || m == "gw" ||
      m == "cov.x" || m == "mean.x") {
    return(matrix(0.0, nrow = length(th.idx), ncol = length(idx)))
  }

  # beta?
  if (!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  } else {
    IB.inv <- get_i_minus_b_inv(MLIST = MLIST)
  }

  if (m == "tau") {
    DX <- matrix(0, nrow = length(th.idx), ncol = nth)
    DX[th.idx > 0L, ] <- diag(nth)
    if (delta.flag) {
      DX <- DX * as.vector(K_nu %*% DELTA)
    }
  } else if (m == "nu") {
    DX <- (-1) * K_nu
    if (delta.flag) {
      DX <- DX * as.vector(K_nu %*% DELTA)
    }
  } else if (m == "lambda") {
    DX <- (-1) * t(IB.inv %*% ALPHA) %x% diag(nvar)
    DX <- K_nu %*% DX
    if (delta.flag) {
      DX <- DX * as.vector(K_nu %*% DELTA)
    }
  } else if (m == "beta") {
    DX <- (-1) * t(IB.inv %*% ALPHA) %x% (LAMBDA %*% IB.inv)
    # this is not really needed (because we select idx=m.el.idx)
    DX[, lav_matrix_diag_idx(nfac)] <- 0.0
    DX <- K_nu %*% DX
    if (delta.flag) {
      DX <- DX * as.vector(K_nu %*% DELTA)
    }
  } else if (m == "alpha") {
    DX <- (-1) * LAMBDA %*% IB.inv
    DX <- K_nu %*% DX
    if (delta.flag) {
      DX <- DX * as.vector(K_nu %*% DELTA)
    }
  } else if (m == "delta") {
    DX1 <- matrix(0, nrow = length(th.idx), ncol = 1)
    DX1[th.idx > 0L, ] <- TAU
    DX2 <- NU + LAMBDA %*% IB.inv %*% ALPHA
    DX2 <- K_nu %*% DX2
    DX <- K_nu * as.vector(DX1 - DX2)
  } else {
    stop()
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}

derivative.pi.LISREL <- function(m = "lambda",
                                 # all model matrix elements, or only a few?
                                 idx = seq_len(length(MLIST[[m]])),
                                 MLIST = NULL) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  GAMMA <- MLIST$gamma
  nexo <- ncol(GAMMA)

  # Delta?
  delta.flag <- FALSE
  if (!is.null(MLIST$delta)) {
    DELTA.diag <- MLIST$delta[, 1L]
    delta.flag <- TRUE
  }

  # shortcut for empty matrices
  if (m == "tau" || m == "nu" || m == "alpha" || m == "psi" ||
      m == "theta" || m == "gw" || m == "cov.x" || m == "mean.x") {
    return(matrix(0.0, nrow = nvar * nexo, ncol = length(idx)))
  }

  # beta?
  if (!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  } else {
    IB.inv <- get_i_minus_b_inv(MLIST = MLIST)
  }

  if (m == "lambda") {
    DX <- t(IB.inv %*% GAMMA) %x% diag(nvar)
    if (delta.flag) {
      DX <- DX * DELTA.diag
    }
  } else if (m == "beta") {
    DX <- t(IB.inv %*% GAMMA) %x% (LAMBDA %*% IB.inv)
    # this is not really needed (because we select idx=m.el.idx)
    DX[, lav_matrix_diag_idx(nfac)] <- 0.0
    if (delta.flag) {
      DX <- DX * DELTA.diag
    }
  } else if (m == "gamma") {
    DX <- diag(nexo) %x% (LAMBDA %*% IB.inv)
    if (delta.flag) {
      DX <- DX * DELTA.diag
    }
  } else if (m == "delta") {
    PRE <- rep(1, nexo) %x% diag(nvar)
    DX <- PRE * as.vector(LAMBDA %*% IB.inv %*% GAMMA)
  } else {
    stop()
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}

derivative.gw.LISREL <- function(m = "gw",
                                 idx = seq_len(length(MLIST[[m]])),
                                 MLIST = NULL) {
  # shortcut for empty matrices
  if (m != "gw") {
    return(matrix(0.0, nrow = 1L, ncol = length(idx)))
  } else {
    # m == "gw"
    DX <- matrix(1.0, 1, 1)
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}

derivative.theta.LISREL <- function(m = "theta",
                                    # all model matrix elements, or only a few?
                                    idx = seq_len(length(MLIST[[m]])),
                                    MLIST = NULL) {
  THETA <- MLIST$theta
  nvar <- nrow(THETA)
  v.idx <- lav_matrix_vech_idx(nvar)

  # shortcut for empty matrices
  if (m != "theta") {
    DX <- matrix(0.0, nrow = length(THETA), ncol = length(idx))
    return(DX[v.idx, , drop = FALSE])
  } else {
    # m == "theta"
    DX <- diag(1, nrow = length(THETA), ncol = length(THETA))
  }

  DX <- DX[v.idx, idx, drop = FALSE]
  DX
}
