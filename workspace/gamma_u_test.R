### =========================================================
### Verify that our unbiased gamma agree with lavaan and
###   Han Du.
### =========================================================

library("lavaan")

# https://github.com/hduquant/Distribution-Free-Estimator/blob/main/supplemental%20material/simulation%20R%20code/simulationhoff%20(skew2).R
handu <- \(x) {
  n <- nrow(x)
  p <- ncol(x)
  ps <- p * (p + 1) / 2

  mean_xt <- apply(x, 2, mean)
  x_c <- x - matrix(1, n, 1) %*% mean_xt
  k <- 1
  sigmaele <- matrix(NA, nrow = ps, ncol = n)

  for (i in 1:p) {
    for (j in i:p) {
      sigmaele[k, ] <- x_c[, i] * x_c[, j]
      k <- k + 1
    }
  }

  kk <- 1
  index <- matrix(NA, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in i:p) {
      index[i, j] <- index[j, i] <- kk
      kk <- kk + 1
    }
  }


  gammaadf.u <- matrix(NA, nrow = ps, ncol = ps)
  for (i in 1:p) {
    for (j in i:p) {
      for (k in 1:p) {
        for (l in 1:p) {
          if (index[k, l] >= index[i, j] & l >= k) {
            #    print(c(i,j,k,l))
            gammaadf.u[index[k, l], index[i, j]] <- n * (n - 1) / (n - 2) / (n - 3) * (sum(sigmaele[index[k, l], ] * sigmaele[index[i, j], ]) / n -
                                                                                         sum(sigmaele[index[k, l], ]) * sum(sigmaele[index[i, j], ]) / n^2) -
              n / (n - 2) / (n - 3) * (sum(sigmaele[index[k, i], ]) * sum(sigmaele[index[l, j], ]) / n^2 +
                                         sum(sigmaele[index[k, j], ]) * sum(sigmaele[index[i, l], ]) / n^2 -
                                         2 / (n - 1) * sum(sigmaele[index[k, l], ]) * sum(sigmaele[index[i, j], ]) / n^2)
            gammaadf.u[index[i, j], index[k, l]] <- gammaadf.u[index[k, l], index[i, j]]
          }
        }
      }
    }
  }
  gammaadf.u
}

# obtained from Njål Foldnes
gammau <- \(Data, meanstructure = FALSE, cov.biased = TRUE) {

  Y <- unname(as.matrix(Data)); N <- nrow(Y); p <- ncol(Y)

  # 'biased' covariance matrix
  COV <- COV.unbiased <- cov(Y)
  if(cov.biased) {
    COV <- COV * (N-1)/N
  }
  cov.vech <- lav_matrix_vech(COV)

  # compute 'biased' Gamma
  Gamma <- lavaan:::lav_samplestats_Gamma(Data, meanstructure = meanstructure)

  # if meanstructure = TRUE, partition Gamma in mean/cov parts
  if(meanstructure) {
    Gamma.cov <- Gamma[-(1:p), -(1:p), drop = FALSE]
    Gamma.mean.cov <- Gamma[1:p, -(1:p), drop = FALSE]
  } else {
    Gamma.cov <- Gamma
  }

  # normal-theory Gamma (cov only)
  GammaNT.cov <- 2 * lav_matrix_duplication_ginv_pre_post(COV %x% COV)

  # Browne's unbiased DF estimator (COV part)
  Gamma1u <- ( N*(N-1)/(N-2)/(N-3) * Gamma.cov -
                 N/(N-2)/(N-3) * ( GammaNT.cov -
                                     2/(N-1) * tcrossprod(cov.vech) )  )

  if(!meanstructure) {
    return(Gamma1u)
  }

  # gammaadf2u is Du Han's unbiased DF estimator (MEAN/COV part)
  gammaadf2u <- Gamma.mean.cov * N/(N-2)

  # assemble gamma
  gammaadfu <- lav_matrix_bdiag(COV.unbiased,  Gamma1u)
  gammaadfu[1:p,(p+1):ncol(gammaadfu)] <- gammaadf2u
  gammaadfu[(p+1):ncol(gammaadfu),1:p] <- t(gammaadf2u)

  gammaadfu
}

# Now we can test them.
n <- nrow(x)
gamma_lavaan <- gammau(x, cov.biased = TRUE)
gamma_our <- gamma_est_unbiased(x, cov(x))
gamma_du <- handu(x)

# Comparisons
sum(abs(gamma_lavaan - gamma_our))
sum(abs(gamma_lavaan - gamma_du))
sum(abs(gamma_our - gamma_du))
