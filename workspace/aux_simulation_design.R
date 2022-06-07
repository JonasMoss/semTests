#generate model syntax based on MEade 2008, but with random loadings
# unit variance in all indicators


#'
#' @param num_fac
#' @param num_indices
#' @param corr
#' @param min
#' @param max
#' @param no_start
#'
#'
#' Meade, A. W., Johnson, E. C., & Braddy, P. W. (2008).
#'   Power and sensitivity of alternative fit indices in tests of measurement
#'   invariance. Journal of Applied Psychology, 93(3), 568–592.
#'   https://doi.org/10.1037/0021-9010.93.3.568

get_model <- function(numFac, numIndic, corr=.3, min=.3, max=.8, nostart=FALSE){
  loadings <- runif(numFac*numIndic, min, max)
  resvars <- paste("xy~~start(",1-loadings^2, ")*xy")
  resvars <- stringr::str_replace_all(resvars, "y", as.character(1:(numFac*numIndic)))
  resvars <- paste(resvars, collapse="\n")
  if(nostart==F){
    model <- sapply(1:numFac, function(k){
      tmp <- paste0("start(", loadings[(k-1)*numIndic+1:numIndic], ")*x", (k-1)*numIndic+1:numIndic)
      paste0("F", k, "=~", paste0(tmp, collapse="+"))
  })
    model <- paste(paste(model, collapse="\n"), "\n", resvars)
    #add correlations
    for(i in 1:(numFac-1))
      model <- paste(model, "\n", paste0("F",i,"~~", paste0("start(",corr,")*F",(i+1):numFac, collapse="+")))

  } else {
    model <- sapply(1:numFac, function(k){
      tmp <- paste0("x", (k-1)*numIndic+1:numIndic)
      paste0("F", k, "=~", paste0(tmp, collapse="+"), ";")
    })
  }

  return(paste(model, collapse="\n"))

}

#test all ones on diagonal
#mymodel <- get_model(4,5)
#lavInspect(sem(mymodel, data=NULL), "sigma.hat")
