# Shared profile geometry in summed criterion units.
inference_solve <- function(A, B = diag(nrow(A))) {
  root <- inference_chol(fiml_sym(A), "Sensitivity/metric")
  backsolve(root, forwardsolve(t(root), B))
}

inference_complement <- function(K, dimension = nrow(K), tolerance = 1e-10) {
  if (!ncol(K)) return(diag(dimension))
  norms <- sqrt(colSums(K^2))
  if (any(norms == 0)) stop("Nuisance basis contains zero columns.", call. = FALSE)
  q <- qr(sweep(K, 2, norms, "/"), tol = tolerance)
  if (q$rank != ncol(K)) stop("Nuisance basis is rank deficient.", call. = FALSE)
  if (ncol(K) >= dimension) stop("There are no testable directions.", call. = FALSE)
  tail <- diag(dimension)[, seq.int(ncol(K) + 1L, dimension), drop = FALSE]
  qr.qy(q, tail)
}

#' Profile a score over nuisance directions
#'
#' All inputs use summed criterion units. Sensitivity is the negative derivative
#' of the estimating score, including residual/model curvature when applicable.
#' Covariance must include all first-stage contributions for an estimated
#' criterion. This function removes nuisance estimation before applying the
#' score quadratic; it does not estimate any of these inputs.
#'
#' @param score Score vector at the restricted estimate.
#' @param curvature Symmetric sensitivity matrix in the alternative space.
#' @param score_covariance Covariance of the summed score in that space.
#' @param nuisance Full-column-rank basis for the null tangent, with zero columns
#'   when there are no nuisance parameters.
#' @param metric Positive metric in the alternative space, defaulting to
#'   sensitivity; `NULL` fully studentizes the profiled score.
#' @param nobs Total number of independent observations, for provenance.
#' @param hypothesis Description of the null and its population target.
#' @param keep Retain inputs and nuisance projection.
#' @param tolerance Matrix validation tolerance.
#' @return A `semTests_test`, usable with [pvalues()].
#' @keywords internal
profile_score_test_matrices <- function(score, curvature, score_covariance, nuisance,
                                      metric = curvature, nobs = NULL,
                                      hypothesis = "Supplied restricted pseudo-target",
                                      keep = FALSE, tolerance = 1e-10) {
  inference_options(hypothesis,keep,tolerance)
  score <- inference_vector(score, "score")
  p <- length(score)
  H <- inference_symmetric(curvature, p, "curvature", tolerance)
  J <- inference_symmetric(score_covariance, p, "score_covariance", tolerance)
  ev <- eigen(J,symmetric=TRUE,only.values=TRUE)$values
  if (min(ev) < -tolerance*max(abs(ev)))
    stop("score_covariance must be positive semidefinite.",call.=FALSE)
  if (!is.matrix(nuisance) || !is.numeric(nuisance) || nrow(nuisance) != p ||
      any(!is.finite(nuisance))) stop("nuisance must be a finite conformable matrix.", call. = FALSE)
  D <- inference_complement(nuisance, p, tolerance)
  G <- D
  if (ncol(nuisance)) {
    G <- D - nuisance %*% inference_solve(crossprod(nuisance, H %*% nuisance),
                                         crossprod(nuisance, H %*% D))
  }
  M <- if (is.null(metric)) NULL else {
    metric <- inference_symmetric(metric, p, "metric", tolerance)
    crossprod(G, metric %*% G)
  }
  out <- score_test_matrices(drop(crossprod(G, score)), M,
    crossprod(G, J %*% G), nobs = nobs, hypothesis = hypothesis,
    keep = keep, tolerance = tolerance)
  out$diagnostics$nuisance_score <- drop(crossprod(nuisance, score))
  if (keep) out$ingredients$projection <- G
  out
}

# Negative Hessian in saturated mean/vech coordinates. X is centered at the
# fitted mean; T includes mean residuals, not just the sample covariance.
gaussian_observed_information <- function(X, geometry) {
  K <- geometry$precision
  q <- X %*% K
  T <- crossprod(q)
  n <- nrow(X)
  i <- geometry$index[, 1L]
  j <- geometry$index[, 2L]
  H <- (K[i,i,drop=FALSE] * T[j,j,drop=FALSE] +
        K[i,j,drop=FALSE] * T[j,i,drop=FALSE] +
        K[j,i,drop=FALSE] * T[i,j,drop=FALSE] +
        K[j,j,drop=FALSE] * T[i,i,drop=FALSE] -
        n * (K[i,i,drop=FALSE] * K[j,j,drop=FALSE] +
             K[i,j,drop=FALSE] * K[j,i,drop=FALSE])) * tcrossprod(geometry$weight)
  if (geometry$meanstructure) {
    v <- colSums(q)
    C <- sweep(K[,i,drop=FALSE], 2, v[j], "*") +
         sweep(K[,j,drop=FALSE], 2, v[i], "*")
    C <- sweep(C, 2, geometry$weight, "*")
    H <- rbind(cbind(n*K, C), cbind(t(C), H))
  }
  fiml_sym(H)
}

# Pattern records retain only observed coordinates. No n by full-moment score
# array is constructed; projection takes place before score cross-products.
gaussian_pattern_parts <- function(object) {
  X <- lavaan::lavInspect(object, "data")
  if (!is.list(X)) X <- list(X)
  moments <- lavaan::lavInspect(object, "implied")
  if (!is.null(moments$cov)) moments <- list(moments)
  means <- isTRUE(object@Options$meanstructure)
  records <- list()
  H <- O <- list()
  ns <- numeric(length(X))
  offset <- 0L
  for (g in seq_along(X)) {
    x <- X[[g]]
    p <- ncol(x)
    if (!identical(colnames(x), colnames(moments[[g]]$cov)) || any(is.infinite(x)))
      stop("Raw observations and implied-moment ordering do not agree.", call. = FALSE)
    # lavaan excludes entirely empty cases from nobs and its likelihood.
    x <- x[rowSums(!is.na(x)) > 0L, , drop = FALSE]
    ns[g] <- nrow(x)
    if (ns[g] != object@SampleStats@nobs[[g]]) stop("Observed case counts do not agree.", call. = FALSE)
    center <- if (means) as.numeric(moments[[g]]$mean) else colMeans(x)
    x <- sweep(x, 2, center, "-")
    fullij <- which(lower.tri(matrix(0,p,p),diag=TRUE), arr.ind=TRUE)
    lookup <- matrix(0L,p,p)
    lookup[fullij] <- seq_len(nrow(fullij)) + if (means) p else 0L
    m <- nrow(fullij) + if (means) p else 0L
    Hg <- Og <- matrix(0,m,m)
    patterns <- split(seq_len(nrow(x)), apply(!is.na(x),1,paste0,collapse=""))
    for (rows in patterns) {
      observed <- which(!is.na(x[rows[1L],]))
      z <- x[rows,observed,drop=FALSE]
      geo <- gaussian_moment_geometry(moments[[g]]$cov[observed,observed,drop=FALSE], means)
      localij <- cbind(observed[geo$index[,1]],observed[geo$index[,2]])
      ix <- c(if (means) observed else integer(), lookup[localij])
      Hg[ix,ix] <- Hg[ix,ix] + nrow(z)*geo$information
      Og[ix,ix] <- Og[ix,ix] + gaussian_observed_information(z,geo)
      records[[length(records)+1L]] <- list(X=z, geometry=geo, index=offset+ix, group=g,
        observed=observed, center=center[observed], nobs=nrow(z),
        total=colSums(z), crossproduct=crossprod(z))
    }
    H[[g]] <- Hg
    O[[g]] <- Og
    offset <- offset+m
  }
  list(records=records, expected=H, observed=O, group_nobs=ns, nobs=sum(ns))
}

pattern_score_products <- function(parts, G, score_covariance, tile_rows=256L) {
  score <- numeric(ncol(G))
  score_cov <- matrix(0,ncol(G),ncol(G))
  group_score <- matrix(0,length(parts$group_nobs),ncol(G))
  for (rec in parts$records) {
    for (first in seq.int(1L,nrow(rec$X),by=tile_rows)) {
      rows <- first:min(nrow(rec$X),first+tile_rows-1L)
      Y <- gaussian_score_rows(rec$X[rows,,drop=FALSE],rec$geometry) %*% G[rec$index,,drop=FALSE]
      group_score[rec$group,] <- group_score[rec$group,]+colSums(Y)
      score_cov <- score_cov+crossprod(Y)
    }
  }
  score <- colSums(group_score)
  if (score_covariance == "centered") score_cov <- score_cov-crossprod(group_score/sqrt(parts$group_nobs))
  list(score=score,score_covariance=fiml_sym(score_cov))
}

gaussian_global_score_parts <- function(object, score_covariance, curvature, metric="information") {
  parts <- gaussian_pattern_parts(object)
  E <- lavaan::lav_matrix_bdiag(parts$expected)
  H <- if (curvature == "expected") E else lavaan::lav_matrix_bdiag(parts$observed)
  K <- model_effective_delta(object)
  D <- inference_complement(K)
  G <- D - K %*% inference_solve(crossprod(K,H%*%K),crossprod(K,H%*%D))
  M <- if (metric=="observed") lavaan::lav_matrix_bdiag(parts$observed) else E
  products <- pattern_score_products(parts,G,score_covariance)
  # Model-based information in the covariance slot. Distinct from the empirical
  # products only where the information equality fails, i.e. under missing data.
  if (score_covariance == "expected") {
    products$score_covariance <- fiml_sym(crossprod(G,E%*%G))
  }
  nuisance <- pattern_score_products(parts,K,"uncentered")$score
  c(products,list(metric=crossprod(G,M%*%G),nobs=parts$nobs,projection=G,
    tangent=K,information=parts$expected,nuisance_gradient=nuisance/parts$nobs,
    group_nobs=parts$group_nobs))
}

# Only these local functions depend on lavaan's internal model evaluator. No
# refitting, calls to the original data expression, or changes to user fits.
inference_model_at <- function(object, theta) {
  object@Model <- get("lav_model_set_parameters",asNamespace("lavaan"))(object@Model,theta)
  object
}

inference_implied <- function(object) {
  get("lav_model_wls_est",asNamespace("lavaan"))(object@Model)
}

inference_delta <- function(object) {
  get("lav_model_delta",asNamespace("lavaan"))(object@Model)
}

# Evaluate pattern scores at a local Model, not the cached fitted @implied slot.
gaussian_parameter_evaluation <- function(object, basis, score_covariance="centered") {
  # lavInspect("implied") reads @implied; refresh only the local copy.
  object@implied <- get("lav_model_implied",asNamespace("lavaan"))(object@Model)
  parts <- gaussian_pattern_parts(object)
  delta <- do.call(rbind,inference_delta(object)) %*% basis
  out <- pattern_score_products(parts,delta,score_covariance)
  E <- matrix(0,ncol(basis),ncol(basis))
  offset <- 0L
  for (H in parts$expected) {
    ix <- offset+seq_len(nrow(H))
    E <- E+crossprod(delta[ix,,drop=FALSE],H%*%delta[ix,,drop=FALSE])
    offset <- offset+nrow(H)
  }
  c(out,list(expected=E,nobs=parts$nobs,group_nobs=parts$group_nobs))
}

criterion_hessian <- function(gradient, theta, basis, step) {
  if (length(step)!=1L || !is.finite(step) || step <= 0 || step > .01)
    stop("derivative_step must be in (0, .01].",call.=FALSE)
  k <- ncol(basis)
  out <- matrix(0,k,k)
  for (j in seq_len(k)) {
    h <- step*max(1,sqrt(sum((basis[,j]*theta)^2)))
    shift <- h*basis[,j]
    out[,j] <- -(gradient(theta+shift)-gradient(theta-shift))/(2*h)
  }
  asymmetry <- max(abs(out-t(out)))/max(1,max(abs(out)))
  if (asymmetry > 1e-3) stop("Numerical criterion curvature is asymmetric; reduce derivative_step.",call.=FALSE)
  list(matrix=fiml_sym(out),asymmetry=asymmetry)
}

check_criterion_adapter <- function(object) {
  check_lavaan(object)
  if (object@Options$estimator == "ML") {
    check_ml_inference_adapter(object,"Criterion inference",allow_fiml=TRUE)
    return(invisible(TRUE))
  }
  check_fit_quality(object)
  if (!object@Options$estimator %in% c("DWLS","WLS","GLS","ULS") ||
      uses_fixed_or_conditional_x(object) || object@Data@nlevels!=1L ||
      !identical(object@Data@data.type,"full") || length(object@Data@cluster) ||
      isTRUE(object@Data@std.ov) || any(lengths(object@Data@weights)>0L) ||
      length(object@Data@sampling.weights) || !identical(object@Options$missing,"listwise") ||
      !is.null(object@call$NACOV) || !is.null(object@call$WLS.V) ||
      !is.null(object@call$sample.cov) || isTRUE(object@Options$correlation) ||
      isTRUE(object@Options$group.w.free))
    stop("Criterion inference requires unweighted, single-level raw-data ML or complete-data LS with native weights and joint random-x inference.",call.=FALSE)
  if (length(object@Model@ceq.nonlinear.idx) || nrow(object@Model@cin.JAC)>0L ||
      any(object@ParTable$op %in% c("<",">")) || any(is.finite(object@ParTable$lower)) ||
      any(is.finite(object@ParTable$upper)))
    stop("Criterion inference supports affine equality constraints only.",call.=FALSE)
  if (isTRUE(object@Model@categorical) &&
      (length(lavaan::lavNames(object,"ov.num")) || !object@Options$estimator %in% c("DWLS","WLS","ULS")))
    stop("The categorical first-stage adapter currently requires all variables ordered.",call.=FALSE)
  invisible(TRUE)
}

parameter_inference_parts <- function(object, curvature="observed",score_covariance="centered",
                                      weight_effect="estimated",derivative_step=1e-4,
                                      basis=model_parameter_basis(object),
                                      theta=as.numeric(lavaan::coef(object)),first_stage=NULL) {
  if (object@Options$estimator == "ML") {
    evaluate <- function(x) gaussian_parameter_evaluation(inference_model_at(object,x),basis,score_covariance)
    base <- evaluate(theta)
    at <- inference_model_at(object,theta)
    at@implied <- get("lav_model_implied",asNamespace("lavaan"))(at@Model)
    records <- gaussian_pattern_parts(at)$records
    gradient <- function(x) gaussian_parameter_gradient(inference_model_at(object,x),basis,records)
  } else {
    if (is.null(first_stage)) first_stage <- ls_first_stage(object,weight_effect,derivative_step)
    evaluate <- function(x) ls_parameter_evaluation(inference_model_at(object,x),basis,first_stage)
    base <- evaluate(theta)
    gradient <- function(x) ls_parameter_gradient(inference_model_at(object,x),basis,first_stage)
  }
  if (score_covariance == "expected") base$score_covariance <- base$expected
  curvature <- if (curvature == "observed") criterion_hessian(gradient,theta,basis,derivative_step) else
    list(matrix=base$expected,asymmetry=0)
  c(base,list(curvature=curvature$matrix,curvature_asymmetry=curvature$asymmetry,
    provenance=list(curvature=curvature,score_covariance=score_covariance,
      weight_effect=if (object@Options$estimator=="ML") "none" else weight_effect,
      derivative_step=derivative_step,first_stage=if (is.null(first_stage)) "Gaussian observed-data likelihood" else first_stage$method)))
}

# The Hessian differentiates aggregate scores. Reuse each pattern's sufficient
# statistics, avoiding repeated casewise covariance construction at every
# parameter perturbation. The empirical score_covariance is evaluated only at the estimate.
gaussian_parameter_gradient <- function(object,basis,records) {
  implied <- get("lav_model_implied",asNamespace("lavaan"))(object@Model)
  delta <- do.call(rbind,inference_delta(object))%*%basis
  total_score <- numeric(nrow(delta))
  means <- isTRUE(object@Options$meanstructure)
  for (rec in records) {
    obs <- rec$observed
    K <- inference_solve(implied$cov[[rec$group]][obs,obs,drop=FALSE])
    shift <- if (means) rec$center-as.numeric(implied$mean[[rec$group]])[obs] else rep(0,length(obs))
    total <- rec$total+rec$nobs*shift
    M <- rec$crossproduct+tcrossprod(rec$total,shift)+tcrossprod(shift,rec$total)+rec$nobs*tcrossprod(shift)
    cov_score <- (K%*%M%*%K-rec$nobs*K)[rec$geometry$index]*rec$geometry$weight
    score <- c(if (means) drop(K%*%total) else numeric(),cov_score)
    total_score[rec$index] <- total_score[rec$index]+score
  }
  drop(crossprod(delta,total_score))
}
