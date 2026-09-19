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

# Differentiate the native ordinal frequency functional, including first-stage
# derivative matrices and their effect on Gamma. This deliberately slow oracle
# path collapses duplicate response patterns and does not refit the SEM.
ordinal_first_stage <- function(object,g,weight_effect,step) {
  X <- object@Data@X[[g]]
  n <- nrow(X)
  p <- ncol(X)
  if (any(!is.finite(X))) stop("Ordinal weight inference requires complete observations.",call.=FALSE)
  levels <- vapply(seq_len(p),function(j) max(X[,j]),numeric(1))
  if (any(vapply(seq_len(p),function(j) !identical(sort(unique(X[,j])),as.numeric(seq_len(levels[j]))),logical(1))))
    stop("Ordinal categories must be consecutive and nonempty.",call.=FALSE)
  opt <- object@Options
  # Tighten only this local numerical evaluator, leaving the lavaan namespace
  # and the fitted object untouched. Default pair optimization is too coarse
  # for differentiating the first-stage functional.
  ns <- asNamespace("lavaan")
  # These five are not part of lavaan's public interface. There is no exported
  # route to "recompute the first-stage polychoric map under perturbed sampling
  # weights with a tighter root", which is what differentiating the first-stage
  # functional needs. Check them up front so a lavaan change fails here, with a
  # message naming the cause, rather than deep inside a numerical derivative.
  #
  # weight_effect does not avoid this. Only the weight derivative dW below is
  # conditional on it; the moment influence Z is differentiated either way, and
  # that is what needs the map. The fallback is pvalues_nested(), which reaches
  # categorical nested comparisons analytically through Satorra 2000.
  required <- c("lav_bvord_cor_twostep_fit", "lav_bvord_cache_from_args",
                "lav_bvord_min_fns", "lav_samp_step2", "muthen1984")
  missing <- required[!vapply(required, exists, logical(1),
                              envir = ns, inherits = FALSE)]
  if (length(missing)) {
    stop("Ordinal first-stage inference relies on lavaan internals that this ",
         "lavaan version does not provide: ", paste(missing, collapse = ", "),
         ". Use pvalues_nested() for categorical nested comparisons.",
         call. = FALSE)
  }
  local <- new.env(parent=ns)
  pair_fit <- get("lav_bvord_cor_twostep_fit",ns)
  # The refined pair fit reaches muthen1984 only through lexical scoping: it
  # calls lav_samp_step2, which calls lav_bvord_cor_twostep_fit, and both are
  # re-parented into `local` below. If lavaan ever resolves either by a route
  # that does not consult `local`, the shadow is bypassed and the derivative is
  # silently computed at lavaan's coarse tolerance instead. Count the calls so
  # that failure is loud. The existing agreement check against WLS.obs cannot
  # catch it: both tolerances agree at the unperturbed weights.
  local$refined_calls <- 0L
  local$lav_bvord_cor_twostep_fit <- function(fit_y1,fit_y2,wt=NULL,
                                                zero_cell_flag=FALSE,...) {
    local$refined_calls <- local$refined_calls + 1L
    initial <- pair_fit(fit_y1=fit_y1,fit_y2=fit_y2,wt=wt,
                        zero_cell_flag=zero_cell_flag,...)
    cache <- get("lav_bvord_cache_from_args",ns)(fit_y1=fit_y1,fit_y2=fit_y2,wt=wt)
    if (any(cache$freq==0)) stop("Ordinal weight derivatives require nonempty pairwise cells.",call.=FALSE)
    fns <- get("lav_bvord_min_fns",ns)()
    gradient <- function(rho) {
      fns$objective(rho,cache)
      as.numeric(fns$gradient(rho,cache))
    }
    root <- stats::uniroot(gradient,c(max(-.999,initial-.05),min(.999,initial+.05)),tol=1e-12)$root
    if (zero_cell_flag) attr(root,"zero.cell.flag") <- FALSE
    root
  }
  step2 <- get("lav_samp_step2",ns)
  environment(step2) <- local
  local$lav_samp_step2 <- step2
  first_stage_map <- get("muthen1984",ns)
  environment(first_stage_map) <- local
  map <- function(w) {
    a <- first_stage_map(data_1=X,ov_names=colnames(X),
      ov_types=rep("ordered",p),ov_levels=levels,wt=w,sampling_weights_type="frequency",
      zero_add=opt$zero.add,zero_keep_margins=opt$zero.keep.margins,
      zero_cell_warn=FALSE,zero_cell_tables=TRUE)
    if (nrow(a$zero.cell.tables)>0L) stop("Ordinal weight derivatives require nonempty pairwise cells; smoothing is not differentiable here.",call.=FALSE)
    Gamma <- a$WLS.W*n
    V <- switch(opt$estimator,DWLS=diag(1/diag(Gamma)),WLS=inference_solve(Gamma),ULS=diag(ncol(Gamma)))
    list(moment=c(unlist(a$TH),a$COR[lower.tri(a$COR)]),weight=V)
  }
  base <- map(rep(1,n))
  if (local$refined_calls == 0L) {
    stop("The refined polychoric root was never reached, so first-stage ",
         "derivatives would be taken at lavaan's default pair tolerance. ",
         "This lavaan version resolves its two-step fit differently. ",
         "Use pvalues_nested() for categorical nested comparisons.",
         call. = FALSE)
  }
  observed <- as.numeric(object@SampleStats@WLS.obs[[g]])
  V <- as.matrix(fiml_group_matrices(lavaan::lavInspect(object,"wls.v"),object,"weights")[[g]])
  if (!isTRUE(all.equal(base$moment,observed,tolerance=2e-6,check.attributes=FALSE)) ||
      !isTRUE(all.equal(as.numeric(base$weight),as.numeric(V),tolerance=2e-6,check.attributes=FALSE)))
    stop("The fitted ordinal first stage does not match the native frequency weight map.",call.=FALSE)
  patterns <- split(seq_len(n),apply(X,1,paste,collapse=":"))
  Z <- matrix(0,n,length(observed))
  dW <- vector("list",length(patterns))
  for (j in seq_along(patterns)) {
    rows <- patterns[[j]]
    # A larger probability step plus Richardson extrapolation suppresses
    # the polychoric optimizer noise amplified by very small perturbations.
    h <- sqrt(step)*length(rows)/n
    direction <- rep(-1,n)
    direction[rows] <- n/length(rows)-1
    plus <- map(1+h*direction)
    minus <- map(1-h*direction)
    plus_half <- map(1+h/2*direction)
    minus_half <- map(1-h/2*direction)
    coarse <- (plus$moment-minus$moment)/(2*h)
    fine <- (plus_half$moment-minus_half$moment)/h
    Z[rows,] <- matrix((4*fine-coarse)/3,length(rows),length(observed),byrow=TRUE)
    if (weight_effect=="estimated") {
      coarse_W <- (plus$weight-minus$weight)/(2*h)
      fine_W <- (plus_half$weight-minus_half$weight)/h
      if (max(abs(fine_W-coarse_W)) > .005*max(1,abs(fine_W)))
        stop("Ordinal weight derivatives are unstable at this step; inspect the first-stage fit.",call.=FALSE)
      derivative <- (4*fine_W-coarse_W)/3
      dW[[j]] <- if (opt$estimator=="DWLS") diag(derivative) else derivative
    }
  }
  Z <- sweep(Z,2,colMeans(Z),"-")
  correction <- function(residual) {
    out <- matrix(0,n,length(residual))
    if (weight_effect=="estimated") for (j in seq_along(patterns)) {
      rows <- patterns[[j]]
      value <- if (opt$estimator=="DWLS") dW[[j]]*residual else drop(dW[[j]]%*%residual)
      out[rows,] <- matrix(value,length(rows),length(residual),byrow=TRUE)
    }
    sweep(out,2,colMeans(out),"-")
  }
  list(observed=observed,weight=V,influence=Z,correction=correction,nobs=n)
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
  list(blocks=blocks,method=if (categorical) "lavaan ordinal frequency-functional numerical IJ" else
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
