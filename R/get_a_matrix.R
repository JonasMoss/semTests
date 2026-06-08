# These functions are copied from lavaan: https://github.com/yrosseel/lavaan

get_a_matrix <- function(m1, m0) {

  delta <- function(m) {
    delta_list <- computeDelta(m@Model)
    delta <- do.call(rbind, delta_list)
    if (m@Model@eq.constraints) {
      return(delta %*% m@Model@eq.constraints.K)
    }

    if (methods::.hasSlot(m@Model, "ceq.simple.only") && m@Model@ceq.simple.only) {
      return(delta %*% t(m@Model@ceq.simple.K))
    }

    delta
  }

  t(get_orthogonal_complement(generalized_inverse(delta(m1)) %*% delta(m0)))
}

computeDelta <- function(lavmodel = NULL, GLIST. = NULL,
                        m.el.idx. = NULL, x.el.idx. = NULL,
                        ceq.simple = FALSE,
                        force.conditional.x.false = FALSE) {

  representation <- lavmodel@representation
  categorical <- lavmodel@categorical
  if (methods::.hasSlot(lavmodel, "correlation")) {
    correlation <- lavmodel@correlation
  } else {
    correlation <- FALSE
  }
  conditional.x <- lavmodel@conditional.x
  group.w.free <- lavmodel@group.w.free
  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks
  nvar <- lavmodel@nvar
  num.idx <- lavmodel@num.idx
  th.idx <- lavmodel@th.idx
  nexo <- lavmodel@nexo
  parameterization <- lavmodel@parameterization

  if (parameterization == "theta") {
    stop("Parameterization = 'theta' is not supported.")
  }

  if (categorical) {
    stop("Categorical is not supported.")
  }

  # number of thresholds per group (if any)
  nth <- sapply(th.idx, function(x) sum(x > 0L))

  # state or final?
  if (is.null(GLIST.)) {
    GLIST <- lavmodel@GLIST
  } else {
    GLIST <- GLIST.
  }

  # type = "free" or something else?
  type <- "nonfree"
  m.el.idx <- m.el.idx.
  x.el.idx <- x.el.idx.
  if (is.null(m.el.idx) && is.null(x.el.idx)) {
    type <- "free"
  }

  # number of rows in DELTA.group
  pstar <- integer(nblocks)
  for (g in 1:nblocks) {
    pstar[g] <- as.integer(nvar[g] * (nvar[g] + 1) / 2)
    if (lavmodel@meanstructure) {
      pstar[g] <- nvar[g] + pstar[g] # first the means, then sigma
    }
    if (correlation) {
      pstar[g] <- pstar[g] - nvar[g] # remove variances
    }
    if (conditional.x && nexo[g] > 0L) {
      pstar[g] <- pstar[g] + (nvar[g] * nexo[g]) # add slopes
    }
    if (group.w.free) {
      pstar[g] <- pstar[g] + 1L # add group weight
    }
  }


  # number of columns in DELTA + m.el.idx/x.el.idx
  if (type == "free") {
    if (methods::.hasSlot(lavmodel, "ceq.simple.only") && lavmodel@ceq.simple.only) {
      NCOL <- lavmodel@nx.unco
    } else {
      NCOL <- lavmodel@nx.free
    }
    m.el.idx <- x.el.idx <- vector("list", length = length(GLIST))
    for (mm in 1:length(GLIST)) {
      m.el.idx[[mm]] <- lavmodel@m.free.idx[[mm]]
      if (methods::.hasSlot(lavmodel, "ceq.simple.only") &&
          lavmodel@ceq.simple.only) {
        x.el.idx[[mm]] <- lavmodel@x.unco.idx[[mm]]
      } else {
        x.el.idx[[mm]] <- lavmodel@x.free.idx[[mm]]
      }
      # handle symmetric matrices
      if (lavmodel@isSymmetric[mm]) {
        # since we use 'x.free.idx', only symmetric elements
        # are duplicated (not the equal ones, only in x.free.free)
        dix <- duplicated(x.el.idx[[mm]])
        if (any(dix)) {
          m.el.idx[[mm]] <- m.el.idx[[mm]][!dix]
          x.el.idx[[mm]] <- x.el.idx[[mm]][!dix]
        }
      }
    }
  } else {
    NCOL <- sum(unlist(lapply(x.el.idx, function(x) length(unique(x)))))
  }


  # compute Delta
  Delta <- vector("list", length = nblocks)
  for (g in 1:nblocks) {
    Delta.group <- matrix(0, nrow = pstar[g], ncol = NCOL)
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]

    for (mm in mm.in.group) {
      mname <- names(lavmodel@GLIST)[mm]

      if (!length(m.el.idx[[mm]])) next

      if (representation == "LISREL") {
        # Sigma
        DELTA <- dxSigma <-
          derivative.sigma.LISREL(
            m = mname,
            idx = m.el.idx[[mm]],
            MLIST = GLIST[mm.in.group],
            delta = parameterization == "delta"
          )

        if (!categorical && correlation) {
          rm.idx <- lav_matrix_diagh_idx(nvar[g])
          DELTA <- DELTA[-rm.idx, , drop = FALSE]
        }

        if (!categorical) {
          if (conditional.x) {
            DELTA.mu <- derivative.mu.LISREL(
              m = mname,
              idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
            )

            if (lavmodel@nexo[g] > 0L) {
              DELTA.pi <- derivative.pi.LISREL(
                m = mname,
                idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
              )

              if (lavmodel@multilevel) {
                DELTA <- rbind(DELTA.mu, DELTA.pi, DELTA)
              } else {
                nEls <- NROW(DELTA.mu) + NROW(DELTA.pi)
                tmp <- rbind(DELTA.mu, DELTA.pi)
                row.idx <- as.vector(matrix(seq.int(nEls),
                                                 nrow = lavmodel@nexo[g] + 1L,
                                                 ncol = lavmodel@nvar[g], byrow = TRUE
                ))
                DELTA.beta <- tmp[row.idx, , drop = FALSE]
                DELTA <- rbind(DELTA.beta, DELTA)
              }
            } else {
              DELTA <- rbind(DELTA.mu, DELTA)
            }
          } else if (!conditional.x && lavmodel@meanstructure) {
            DELTA.mu <- derivative.mu.LISREL(
              m = mname,
              idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
            )
            DELTA <- rbind(DELTA.mu, DELTA)
          }
        }
        if (group.w.free) {
          DELTA.gw <- derivative.gw.LISREL(
            m = mname,
            idx = m.el.idx[[mm]],
            MLIST = GLIST[mm.in.group]
          )
          DELTA <- rbind(DELTA.gw, DELTA)
        }
      } else if (representation == "RAM") {
        stop("RAM representation not supported")
      } else {
        stop()
      }
      Delta.group[, x.el.idx[[mm]]] <- DELTA
    } # mm

    # if type == "free" take care of equality constraints
    if (type == "free" && ceq.simple &&
        methods::.hasSlot(lavmodel, "ceq.simple.only") && lavmodel@ceq.simple.only) {
      Delta.group <- Delta.group %*% lavmodel@ceq.simple.K
    }

    Delta[[g]] <- Delta.group
  } # g

  # if multilevel, rbind levels within group
  if (methods::.hasSlot(lavmodel, "multilevel") && lavmodel@multilevel) {
    DELTA <- vector("list", length = lavmodel@ngroups)
    for (g in 1:lavmodel@ngroups) {
      DELTA[[g]] <- rbind(
        Delta[[(g - 1) * 2 + 1]],
        Delta[[(g - 1) * 2 + 2]]
      )
    }
    Delta <- DELTA
  }

  Delta
}
