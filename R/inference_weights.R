# Estimated-weight first stages. Every block carries jointly computed moment
# influences Z and weight influences contracted with the current residual.
# Centering is groupwise because group sample sizes are conditioned upon.

continuous_first_stage <- function(object, g, weight_effect) {
  X <- object@Data@X[[g]]
  n <- nrow(X)
  p <- ncol(X)
  means <- isTRUE(object@Options$meanstructure)
  X <- sweep(X,2,colMeans(X),"-")
  S <- crossprod(X)/n
  ij <- which(lower.tri(S,diag=TRUE),arr.ind=TRUE)
  Zc <- sweep(X[,ij[,1],drop=FALSE]*X[,ij[,2],drop=FALSE],2,S[ij],"-")
  Z <- if (means) cbind(X,Zc) else Zc
  scale_cov <- if (isTRUE(object@Options$sample.cov.rescale)) 1 else n/(n-1)
  moment_IF <- if (means) cbind(X,Zc*scale_cov) else Zc*scale_cov
  Gamma <- crossprod(Z)/n
  V <- as.matrix(fiml_group_matrices(lavaan::lavInspect(object,"wls.v"),object,"weights")[[g]])
  est <- object@Options$estimator
  expected_weight <- switch(est,
    WLS=inference_solve(Gamma),
    DWLS=diag(1/diag(inference_solve(gaussian_moment_geometry(S*scale_cov,means)$information))),
    GLS=gaussian_moment_geometry(S*scale_cov,means)$information,
    ULS=diag(ncol(Z)))
  if (!isTRUE(all.equal(as.numeric(V),as.numeric(expected_weight),tolerance=1e-7,check.attributes=FALSE)))
    stop("The continuous fitted weight does not match the supported empirical weight map.",call.=FALSE)
  # Derivative of Gamma with respect to the estimated centering mean. The
  # derivative with respect to S is constant in cases and E[Z]=0, so vanishes.
  C <- lapply(seq_len(p),function(j) {
    Dz <- -sweep(X[,ij[,1],drop=FALSE],2,as.numeric(ij[,2]==j),"*") -
           sweep(X[,ij[,2],drop=FALSE],2,as.numeric(ij[,1]==j),"*")
    if (means) {
      Dm <- matrix(0,n,p)
      Dm[,j] <- -1
      Dz <- cbind(Dm,Dz)
    }
    A <- crossprod(Z,Dz)/n
    A+t(A)
  })
  correction <- function(residual) {
    if (weight_effect=="fixed" || est=="ULS") return(matrix(0,n,ncol(Z)))
    if (est=="WLS") {
      v <- drop(V%*%residual)
      dGv <- sweep(Z,1,drop(Z%*%v),"*") -
        matrix(drop(Gamma%*%v),n,ncol(Z),byrow=TRUE)
      for (j in seq_len(p)) dGv <- dGv + tcrossprod(X[,j],drop(C[[j]]%*%v))
      return(-dGv%*%V)
    }
    # GLS and continuous DWLS use the normal-theory weight functional.
    deriv <- matrix(0,ncol(Z),nrow(ij))
    for (j in seq_len(nrow(ij))) {
      h <- 1e-5*max(1,abs((S*scale_cov)[ij[j,,drop=FALSE]]))
      E <- matrix(0,p,p)
      E[ij[j,1],ij[j,2]] <- E[ij[j,2],ij[j,1]] <- h
      plus <- gaussian_moment_geometry(S*scale_cov+E,means)$information
      minus <- gaussian_moment_geometry(S*scale_cov-E,means)$information
      if (est=="DWLS") {
        plus <- diag(1/diag(inference_solve(plus)))
        minus <- diag(1/diag(inference_solve(minus)))
      }
      deriv[,j] <- (plus-minus)%*%residual/(2*h)
    }
    (Zc*scale_cov)%*%t(deriv)
  }
  list(observed=as.numeric(object@SampleStats@WLS.obs[[g]]),weight=V,
       influence=moment_IF,correction=correction,nobs=n)
}

# The ordinal first stage: the moment influence of the thresholds and
# polychorics and, for estimated weights, the influence of the weight matrix,
# all in closed form (polychoric_influence.R).
ordinal_first_stage <- function(object,g,weight_effect,step) {
  X <- object@Data@X[[g]]
  n <- nrow(X)
  p <- ncol(X)
  if (any(!is.finite(X))) stop("Ordinal weight inference requires complete observations.",call.=FALSE)
  levels <- vapply(seq_len(p),function(j) max(X[,j]),numeric(1))
  if (any(vapply(seq_len(p),function(j) !identical(sort(unique(X[,j])),as.numeric(seq_len(levels[j]))),logical(1))))
    stop("Ordinal categories must be consecutive and nonempty.",call.=FALSE)
  # lavaan smooths empty cells of 2x2 tables (zero.add), which the closed-form
  # first stage does not follow, so require every pairwise table to be full.
  for (u in seq_len(p - 1L)) for (v in seq(u + 1L, p)) {
    cells <- tabulate(X[, u] + levels[u] * (X[, v] - 1), levels[u] * levels[v])
    if (any(cells == 0))
      stop("Ordinal inference requires nonempty pairwise cells; ",
           object@Data@ov.names[[g]][u], " and ", object@Data@ov.names[[g]][v],
           " have an empty one.",
           call. = FALSE)
  }
  opt <- object@Options
  observed <- as.numeric(object@SampleStats@WLS.obs[[g]])
  V <- as.matrix(fiml_group_matrices(lavaan::lavInspect(object,"wls.v"),object,"weights")[[g]])

  # Moment influence, in closed form. Thresholds and polychorics are stacked in
  # that order in WLS.obs, so split them back out to index by variable.
  sizes <- levels - 1
  offsets <- c(0, cumsum(sizes))
  thresholds <- lapply(seq_len(p), function(j) observed[offsets[j] + seq_len(sizes[j])])
  Z <- ordinal_moment_influence(X, thresholds, observed[-seq_len(sum(sizes))])

  # Fixed weights need no derivative of the weight functional, and the moment
  # influence above is already exact, so nothing below is reached.
  if (weight_effect != "estimated") {
    return(list(observed=observed, weight=V, influence=Z,
                correction=function(residual) matrix(0,n,length(residual)), nobs=n))
  }
  # ULS holds the weight fixed at the identity, so it has no weight influence
  # either, whatever weight_effect says.
  if (opt$estimator == "ULS") {
    return(list(observed=observed, weight=V, influence=Z,
                correction=function(residual) matrix(0,n,length(residual)), nobs=n))
  }
  # DWLS weights by the reciprocal of Gamma's diagonal, so its weight influence
  # needs only that diagonal. Full WLS inverts all of Gamma, so it needs every
  # entry: dW_i = -W dGamma_i W. Both are closed form.
  wi <- ordinal_weight_influence(X, thresholds, observed[-seq_len(sum(sizes))],
                                 step = sqrt(step),
                                 full = opt$estimator == "WLS")
  correction <- if (opt$estimator == "DWLS") {
    scale <- diag(wi$gamma)^2
    function(residual) {
      sweep(sweep(-wi$diagonal, 2L, scale, "/"), 2L, residual, "*")
    }
  } else {
    function(residual) -wi$derivative(drop(V %*% residual)) %*% V
  }
  list(observed = observed, weight = V, influence = Z,
       correction = function(residual) {
         out <- correction(residual)
         sweep(out, 2L, colMeans(out), "-")
       },
       nobs = n)
}

ls_first_stage <- function(object,weight_effect,step) {
  if (isTRUE(object@Options$gamma.n.minus.one) || isTRUE(object@Options$gamma.unbiased))
    stop("Estimated-weight inference does not support alternative Gamma normalization.",call.=FALSE)
  if (length(step)!=1L || !is.finite(step) || step<=0 || step>.01)
    stop("derivative_step must be in (0, .01].",call.=FALSE)
  categorical <- isTRUE(object@Model@categorical)
  blocks <- lapply(seq_len(object@Data@ngroups),function(g) {
    if (categorical) ordinal_first_stage(object,g,weight_effect,step) else
      continuous_first_stage(object,g,weight_effect)
  })
  list(blocks=blocks,method=if (categorical) "closed-form ordinal first-stage influence" else
    "empirical continuous moments with analytic Gamma influence")
}

ls_parameter_gradient <- function(object,basis,first_stage) {
  delta <- inference_delta(object)
  fitted <- inference_implied(object)
  out <- numeric(ncol(basis))
  for (g in seq_along(first_stage$blocks)) {
    b <- first_stage$blocks[[g]]
    D <- delta[[g]]%*%basis
    residual <- b$observed-fitted[[g]]
    out <- out+(b$nobs-1)*drop(crossprod(D,b$weight%*%residual))
  }
  out
}

ls_parameter_evaluation <- function(object,basis,first_stage) {
  delta <- inference_delta(object)
  fitted <- inference_implied(object)
  k <- ncol(basis)
  J <- E <- matrix(0,k,k)
  score <- numeric(k)
  n <- numeric(length(first_stage$blocks))
  for (g in seq_along(first_stage$blocks)) {
    b <- first_stage$blocks[[g]]
    n[g] <- b$nobs
    D <- delta[[g]]%*%basis
    residual <- b$observed-fitted[[g]]
    score <- score+(n[g]-1)*drop(crossprod(D,b$weight%*%residual))
    Y <- ((n[g]-1)/n[g])*(b$influence%*%b$weight+b$correction(residual))%*%D
    J <- J+crossprod(Y)
    E <- E+(n[g]-1)*crossprod(D,b$weight%*%D)
  }
  list(score=score,score_covariance=symmetrize(J),expected=symmetrize(E),nobs=sum(n),group_nobs=n)
}
