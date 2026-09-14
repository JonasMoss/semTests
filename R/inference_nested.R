# Exact affine embedding in the alternative's full free-parameter coordinates.
# Different latent-variable parameterizations require their own map.
nested_criterion_embedding <- function(object, alternative, tolerance) {
  check_lavaan(alternative,"alternative")
  # Reuse the comprehensive data/criterion check without requiring an optimized
  # alternative. Information/SE choices do not define the criterion here.
  if (!identical(object@Options$estimator,alternative@Options$estimator))
    stop("Nested models must use the same fitting criterion.",call.=FALSE)
  comparable <- alternative
  comparable@Options$estimator.orig <- requested_estimator(object)
  for (nm in c("information","h1.information","observed.information"))
    comparable@Options[[nm]] <- object@Options[[nm]]
  check_nested_pair(object,comparable)
  if (isTRUE(object@Model@categorical)) {
    if (!identical(object@Options$parameterization,alternative@Options$parameterization))
      stop("Nested models must share the categorical parameterization.",call.=FALSE)
  }
  if (length(alternative@Model@ceq.nonlinear.idx) || nrow(alternative@Model@cin.JAC)>0L ||
      any(alternative@ParTable$op %in% c("<",">")) ||
      any(is.finite(alternative@ParTable$lower)) || any(is.finite(alternative@ParTable$upper)))
    stop("The alternative must have only affine equality constraints.",call.=FALSE)
  keys <- function(pt) do.call(paste,c(unname(pt[intersect(c("lhs","op","rhs","group","level","block"),names(pt))]),sep="\r"))
  pt0 <- lavaan::parTable(object)
  pt1 <- lavaan::parTable(alternative)
  full1 <- pt1[pt1$free>0,,drop=FALSE]
  full1 <- full1[order(full1$free),,drop=FALSE]
  if (anyDuplicated(full1$free) || !identical(as.integer(full1$free),seq_len(alternative@Model@nx.free)))
    stop("Alternative free-parameter ordering is ambiguous.",call.=FALSE)
  match0 <- match(keys(full1),keys(pt0))
  if (anyNA(match0)) stop("The restricted parameter table must contain every alternative parameter, including fixed zeros.",call.=FALSE)
  theta <- pt0$est[match0]
  free0 <- pt0$free[match0]
  map <- matrix(0,nrow(full1),object@Model@nx.free)
  ix <- which(free0>0)
  map[cbind(ix,free0[ix])] <- 1
  K0 <- map%*%model_parameter_basis(object)
  K1 <- model_parameter_basis(alternative)
  # K bases may be nonorthonormal (lavaan explicit equalities).
  A <- inference_solve(crossprod(K1),crossprod(K1,K0))
  if (max(abs(K1%*%A-K0)) > tolerance*100*max(1,abs(K0)))
    stop("The restricted tangent is not nested in the alternative.",call.=FALSE)
  # Check the affine offset as well as its tangent.
  # Unoptimized starting values need not satisfy retained affine constraints.
  # Use the constraint offset itself, never the template's initial estimates.
  origin <- if (isTRUE(alternative@Model@eq.constraints))
    alternative@Model@eq.constraints.k0 else numeric(length(theta))
  displacement <- theta-origin
  projected <- K1%*%inference_solve(crossprod(K1),crossprod(K1,displacement))
  if (max(abs(displacement-projected)) > 1e-6*max(1,abs(theta)))
    stop("The restricted estimate does not satisfy alternative affine constraints.",call.=FALSE)
  at <- inference_model_at(alternative,theta)
  expected <- inference_implied(object)
  actual <- inference_implied(at)
  if (!isTRUE(all.equal(expected,actual,tolerance=1e-7,check.attributes=FALSE)))
    stop("The alternative at the restricted estimate has different implied moments; parameterizations do not align.",call.=FALSE)
  D0 <- do.call(rbind,inference_delta(object))%*%model_parameter_basis(object)
  D1 <- do.call(rbind,inference_delta(at))%*%K1
  if (max(abs(D1%*%A-D0)) > 1e-6*max(1,abs(D0)))
    stop("Nested moment tangents do not align.",call.=FALSE)
  list(theta=theta,basis=K1,nuisance=A,object=at)
}

#' Joint nested profile-score test at a restricted fit
#'
#' Tests the common pseudo-target restriction under one shared criterion.
#' The observed sensitivity retains full model curvature; the score covariance
#' includes fitted-mean effects for ML and joint moment/weight effects for LS.
#' The alternative supplies the parameterization and may be created with
#' `do.fit = FALSE`: it is evaluated at the restricted estimates without refitting.
#' Use the same latent-variable scaling and include released paths explicitly
#' as fixed parameters in the restricted model. Only affine constraints and
#' exactly aligned parameter tables are supported.
#'
#' @param object Converged restricted lavaan fit.
#' @param alternative Nested alternative lavaan fit or unoptimized template,
#'   using the same observations, estimator and moment conventions.
#' @param metric `"information"` uses the selected sensitivity for the quadratic,
#'   `"expected"` uses Fisher/Gauss--Newton curvature while retaining the selected
#'   nuisance sensitivity, and `"sandwich"` fully studentizes the profile score.
#' @param curvature `"observed"` (default) or `"expected"` comparator.
#' @param score_covariance Group-centered or uncentered Gaussian score products.
#'   LS influence functions are centered by definition and require `"centered"`.
#' @param weight_effect `"estimated"` includes weight estimation; `"fixed"`
#'   supplies a comparator for LS. Ignored for ML.
#' @param derivative_step Relative central-difference step for full curvature.
#'   Ordinal frequency perturbations use its square root times pattern frequency
#'   with Richardson extrapolation; see [wald_test()].
#' @param hypothesis Description of the common pseudo-target null.
#' @param keep Retain profile matrices and estimating-equation ingredients.
#' @param tolerance Matrix validation tolerance.
#' @return A `semTests_test`. Expected curvature or fixed weights alone do not
#'   establish validity under fixed structural misspecification.
#' @keywords internal
score_test_nested <- function(object,alternative,
                              metric=c("information","expected","sandwich"),
                              curvature=c("observed","expected"),
                              score_covariance=c("centered","uncentered","expected"),
                              weight_effect=c("estimated","fixed"),
                              derivative_step=1e-4,
                              hypothesis="Common criterion pseudo-target satisfies the restricted model",
                              keep=FALSE,tolerance=1e-10) {
  check_criterion_adapter(object)
  metric <- match.arg(metric)
  curvature <- match.arg(curvature)
  score_covariance <- match.arg(score_covariance)
  weight_effect <- match.arg(weight_effect)
  inference_options(hypothesis,keep,tolerance)
  if (object@Options$estimator!="ML" && score_covariance!="centered")
    stop("LS uses centered first-stage influence functions; choose score_covariance = 'centered'.",call.=FALSE)
  embed <- nested_criterion_embedding(object,alternative,tolerance)
  first <- if (object@Options$estimator=="ML") NULL else ls_first_stage(object,weight_effect,derivative_step)
  parts <- parameter_inference_parts(embed$object,curvature,score_covariance,weight_effect,
    derivative_step,basis=embed$basis,theta=embed$theta,first_stage=first)
  M <- switch(metric,information=parts$curvature,expected=parts$expected,sandwich=NULL)
  out <- profile_score_test_matrices(parts$score,parts$curvature,parts$score_covariance,
    embed$nuisance,M,parts$nobs,hypothesis,keep,tolerance)
  out$provenance <- c(out$provenance,parts$provenance,list(
    lavaan_version=as.character(utils::packageVersion("lavaan")),
    estimator=requested_estimator(object),metric=metric,group_nobs=parts$group_nobs))
  out$provenance$backend <- "lavaan restricted criterion"
  out$diagnostics$nuisance_stationarity_norm <- max(abs(out$diagnostics$nuisance_score),0)/parts$nobs
  out$diagnostics$curvature_asymmetry <- parts$curvature_asymmetry
  if (keep) out$ingredients$estimating_equations <- parts
  out
}

#' Nested criterion-difference test with a complete profile reference
#'
#' Uses twice the difference of the summed fitting criteria, with the spectrum
#' computed at the restricted fit using full profile curvature. The null is a
#' common population minimizer under the shared criterion. Both fits must
#' converge and use identical first-stage moments and weights. This is distinct
#' from the existing exact-fit nested corrections in [pvalues_nested()].
#'
#' @inheritParams score_test_nested
#' @return A `semTests_test` with a weighted chi-square reference. The reference
#'   is evaluated at the restricted fit; negative criterion differences fail.
#' @keywords internal
criterion_test_nested <- function(object,alternative,
                                  curvature=c("observed","expected"),
                                  score_covariance=c("centered","uncentered","expected"),
                                  weight_effect=c("estimated","fixed"),
                                  derivative_step=1e-4,
                                  hypothesis="Equal population minimum under the shared criterion",
                                  keep=FALSE,tolerance=1e-10) {
  check_lavaan(object)
  check_lavaan(alternative,"alternative")
  check_fit_quality(alternative,"alternative")
  if (object@Options$estimator!="ML") {
    for (slot in c("WLS.obs","WLS.V","WLS.VD")) {
      if (!isTRUE(all.equal(methods::slot(object@SampleStats,slot),
                           methods::slot(alternative@SampleStats,slot),tolerance=0)))
        stop("Nested criterion differences require identical first-stage moments and weights.",call.=FALSE)
    }
  }
  out <- score_test_nested(object,alternative,metric="information",
    curvature=match.arg(curvature),score_covariance=match.arg(score_covariance),
    weight_effect=match.arg(weight_effect),derivative_step=derivative_step,
    hypothesis=hypothesis,keep=keep,tolerance=tolerance)
  statistic <- 2*out$nobs*(as.numeric(lavaan::fitMeasures(object,"fmin"))-
                          as.numeric(lavaan::fitMeasures(alternative,"fmin")))
  if (!is.finite(statistic) || statistic < -tolerance)
    stop("The restricted fitting criterion is smaller than the alternative; check optimization.",call.=FALSE)
  out$type <- "criterion difference"
  out$statistic <- max(0,statistic)
  out$p_value <- pvalue_all(out$statistic,out$eigenvalues)
  out$provenance$reference_fit <- "restricted"
  out
}
