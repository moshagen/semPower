##########################  determine ncp, chi-square from Fmin   #####################


#' getNCP
#'
#' Computes the non-centrality parameter from the population minimum of the fit-function 
#' (dividing by N - 1 following the Wishart likelihood): `ncp = (N - 1) * F0`.
#'
#' @param Fmin population minimum of the fit-function (can be a list for multiple group models).
#' @param n number of observations (can be a list for multiple group models).
#' @return Returns the implied NCP.
getNCP <- function(Fmin, n){
  NCP <- unlist(Fmin) * (unlist(n) - 1)
  sum(NCP)
}


#' getChiSquare.NCP
#'
#' Computes chi-square from the non-centrality parameter: `chi-square = ncp + df`.
#'
#' @param NCP non-centrality parameter
#' @param df model degrees of freedom
#' @return Returns chi-square
getChiSquare.NCP <- function(NCP, df){
  chiSquare <- NCP + df
  chiSquare
}

#' getChiSquare.F
#'
#' Computes the (Wishart-) chi-square from the population minimum of the fit-function: 
#' `chi-square = (N - 1) * F0 + df = ncp + df`. Note that F0 is the population minimum. 
#' Using F_hat would give `chi-square = (N - 1) * F_hat`.
#'
#' @param Fmin population minimum of the fit-function (can be a list for multiple group models).
#' @param n number of observations  (can be a list for multiple group models).
#' @param df model degrees of freedom
#' @return Returns chi-square
getChiSquare.F <- function(Fmin, n, df){
  chiSquare <- getNCP(Fmin, n) + df
  chiSquare
}


##########################  determine Fmin from RMSEA , Mc , GFI , AGFI   #####################


#' getF
#' 
#' Computes the minimum of the ML-fit-function from known fit indices.
#' 
#' @param effect magnitude of effect
#' @param effect.measure measure of effect, one of `'fmin'`, `'rmsea'`, `'agfi'`, `'gfi'`, `'mc'`
#' @param df model degrees of freedom
#' @param p number of observed variables
#' @param SigmaHat model implied covariance matrix
#' @param Sigma observed (or population) covariance matrix
#' @param muHat model implied mean vector
#' @param mu observed (or population) mean vector
#' @return Returns Fmin
getF <- function(effect, effect.measure, df = NULL, p = NULL, SigmaHat = NULL, Sigma = NULL, muHat = NULL, mu = NULL){
  fmin <- effect
  if(is.null(SigmaHat)){ # sufficient to check for on NULL matrix; primary validity check is in validateInput
    switch(effect.measure,
           "RMSEA" = fmin <- getF.RMSEA(effect, df),
           "MC" = fmin <- getF.Mc(effect),
           "GFI" = fmin <- getF.GFI(effect, p),
           "AGFI" = fmin <- getF.AGFI(effect, df, p)
    )
  }else{
    fmin <- effect <- getF.Sigma(SigmaHat, Sigma, muHat, mu)
  }
  fmin
}


#' getF.RMSEA
#'
#' Computes the minimum of the ML-fit-function from RMSEA:
#' `F_min = rmsea^2 * df`.
#'
#' @param RMSEA RMSEA
#' @param df model degrees of freedom
#' @return Returns Fmin
getF.RMSEA <- function(RMSEA, df){
  fmin <- RMSEA^2 * df
  fmin
}


#' getF.Mc
#'
#' Computes the minimum of the ML-fit-function from Mc.
#'
#' @param Mc Mc
#' @return Returns Fmin
getF.Mc <- function(Mc){
  fmin <- -2 * log(Mc)
  fmin
}

#' getF.GFI
#'
#' Computes the minimum of the ML-fit-function from GFI.
#'
#' @param GFI GFI
#' @param p number of observed variables
#' @return Returns Fmin
getF.GFI <- function(GFI, p){
  fmin <- -(GFI - 1) * p / (2 * GFI)
  fmin
}


#' getF.AGFI
#'
#' Computes the minimum of the ML-fit-function from AGFI.
#'
#' @param AGFI AGFI
#' @param df model degrees of freedom
#' @param p number of observed variables
#' @return Returns Fmin
getF.AGFI <- function(AGFI, df, p){
  fmin <- -(AGFI - 1) * df * p / (p * p + p + (2 * AGFI - 2) * df)
  fmin
}


##########################  determine RMSEA Mc GFI AGFI from Fmin   #####################

#' getIndices.F
#'
#' Computes known indices from the minimum of the ML-fit-function.
#'
#' @param fmin minimum of the ML-fit-function
#' @param df model degrees of freedom
#' @param p number of observed variables
#' @param SigmaHat model implied covariance matrix
#' @param Sigma population covariance matrix
#' @param muHat model implied means
#' @param mu population means
#' @param N list of sample weights
#' @return list of indices
getIndices.F <- function(fmin, df, p = NULL, SigmaHat = NULL, Sigma = NULL, muHat = NULL, mu = NULL, N = NULL){
  # may happen in simulated power
  if(fmin <= 0){
    fit <- list(
      rmsea = 0,
      mc = 1,
      gfi = 1, 
      agfi = 1
    )
  }else{
    fit <- list(
      rmsea = getRMSEA.F(fmin, df, nGroups = ifelse(length(N) > 1, length(N), 1)),
      mc = getMc.F(fmin),
      gfi = NULL,
      agfi = NULL
    )
    if(!is.null(p)){
      fit[['gfi']] <- getGFI.F(fmin, p)
      fit[['agfi']] <- getAGFI.F(fmin, df, p)
    }
  }
  if(!is.null(SigmaHat)){
    if(length(N) > 1){
      fit[['srmr']] <- getSRMR.Sigma.mgroups(SigmaHat, Sigma, muHat, mu, N)
      fit[['cfi']] <- getCFI.Sigma.mgroups(SigmaHat, Sigma, muHat, mu, N)
    }else{
      fit[['srmr']] <- getSRMR.Sigma(SigmaHat, Sigma, muHat, mu)
      fit[['cfi']] <- getCFI.Sigma(SigmaHat, Sigma, muHat, mu)
    }
  }
  fit
}


#' getRMSEA.F
#'
#' Computes RMSEA from the minimum of the ML-fit-function
#' `F_min = rmsea^2 * df`.
#'
#' @param Fmin minimum of the ML-fit-function
#' @param df model degrees of freedom
#' @param nGroups the number of groups
#' @return Returns RMSEA
getRMSEA.F <- function(Fmin, df, nGroups = 1){
  RMSEA <- sqrt(Fmin / df) * sqrt(nGroups)
  RMSEA
}

#' getMc.F
#'
#' Computes Mc from the minimum of the ML-fit-function.
#'
#' @param Fmin minimum of the ML-fit-function
#' @return Returns Mc
getMc.F <- function(Fmin){
  Mc <- exp(-.5 * Fmin)
  Mc
}


#' getGFI.F
#'
#' Computes GFI from the minimum of the ML-fit-function.
#'
#' @param Fmin minimum of the ML-fit-function
#' @param p number of observed variables
#' @return Returns GFI
getGFI.F <- function(Fmin, p){
  GFI <- p / (p + 2 * Fmin)
  GFI
}


#' getAGFI.F
#'
#' Computes AGFI from the minimum of the ML-fit-function.
#'
#' @param Fmin minimum of the ML-fit-function
#' @param df model degrees of freedom
#' @param p number of observed variables
#' @return Returns AGFI
getAGFI.F <- function(Fmin, df, p){
  AGFI <- -(Fmin * p * p + (Fmin - df) * p - 2 * df * Fmin) / (df * p + 2 * df * Fmin)
  AGFI
}



##########################  calculate Fmin RMSEA SRMR CFI from covariance matrix #####################


#' getF.Sigma
#'
#' Computes the minimum of the ML-fit-function given the model-implied and the observed (or population) covariance matrix:
#' `F_min = tr(S %*% SigmaHat^-1) - p + ln(det(SigmaHat)) - ln(det(S))`. When a meanstructure is included, 
#' `(mu - muHat)' SigmaHat^-1 (mu - muHat)` is added.
#'
#' @param SigmaHat model implied covariance matrix
#' @param S observed (or population) covariance matrix
#' @param muHat model implied mean vector
#' @param mu observed (or population) mean vector
#' @return Returns Fmin
getF.Sigma <- function(SigmaHat, S, muHat = NULL, mu = NULL){
  checkPositiveDefinite(SigmaHat)
  checkPositiveDefinite(S)
  fmin <- sum(diag(S %*% solve(SigmaHat))) + log(det(SigmaHat)) - log(det(S)) - ncol(S)
  if(!is.null(mu)){
    fmean <- t(c(mu - muHat)) %*% solve(SigmaHat) %*% c(mu - muHat)
    fmin <- fmin + fmean  
  }
  fmin
}

#' getSRMR.Sigma
#'
#' Computes SRMR given the model-implied and the observed (or population) covariance matrix, 
#' using the Hu & Bentler approach to standardization.
#'
#' @param SigmaHat model implied covariance matrix
#' @param S observed (or population) covariance matrix
#' @param muHat model implied mean vector
#' @param mu observed (or population) mean vector
#' @return Returns SRMR
getSRMR.Sigma <- function(SigmaHat, S, muHat = NULL, mu = NULL){
  checkPositiveDefinite(SigmaHat)
  checkPositiveDefinite(S)
  p <- ncol(S)
  
  # bollen approach to standardization
  # m <- cov2cor(S) - cov2cor(SigmaHat)
  
  # hu+bentler approach to std
  sqrt.d <- 1 / sqrt(diag(S))
  D <- diag(sqrt.d, ncol = length(sqrt.d))
  m <- D %*% (S - SigmaHat) %*% D
  
  fols <- sum(m[lower.tri(m, diag = TRUE)]^2)
  
  # mplus variant
  #fols <- (sum(m[lower.tri(m, diag=F)]^2)  +  sum(((diag(S) - diag(SigmaHat))/diag(S))^2)) 
  
  enum <- fols
  denum <- (p * (p + 1) / 2)
  
  if(!is.null(mu)){
    # mplus / bollen approach
    #stdMeanResid <- sum( mu/sqrt(diag(S)) - muHat/sqrt(diag(SigmaHat)) )^2
    
    # hu+bentler approach
    stdMeanResid <- sum( (D %*% (mu - muHat))^2 )
    
    enum <- enum + stdMeanResid
    denum <- denum + p
  }
  
  srmr <- sqrt(enum / denum)    
  srmr
}

#' getSRMR.Sigma.mgroups 
#'
#' Computes SRMR given the model-implied and the observed (or population) covariance matrix for multiple group models
#' using the Hu & Bentler approach to standardization and the MPlus approach to multiple group sampling weights 
#' (weight squared sums of residuals).
#'
#' @param SigmaHat a list of model implied covariance matrices
#' @param S a list of observed (or population) covariance matrices
#' @param muHat model implied mean vector
#' @param mu observed (or population) mean vector
#' @param N a list of group weights
#' @return Returns SRMR
getSRMR.Sigma.mgroups <- function(SigmaHat, S, muHat = NULL, mu = NULL, N){
  if(is.null(mu)){
    srmrs <- sapply(seq_along(SigmaHat), function(x) getSRMR.Sigma(SigmaHat[[x]], S[[x]]))
  }else{
    srmrs <- sapply(seq_along(SigmaHat), function(x) getSRMR.Sigma(SigmaHat[[x]], S[[x]], muHat[[x]], mu[[x]]))
  }
  
  # lavaan approach: apply sample weights to srmr
  # srmr <- (sum(srmrs*N)/sum(N))
  # mplus approach: apply sample weights to squared sums of res
  srmr <- sqrt( sum(unlist(srmrs)^2 * unlist(N)) / sum(unlist(N)) )
  srmr
}


#' getCFI.Sigma
#'
#' Computes CFI given the model-implied and the observed (or population) covariance matrix: 
#' `CFI = (F_null - F_hyp) / F_null`.
#'
#' @param SigmaHat model implied covariance matrix
#' @param S observed (or population) covariance matrix
#' @param muHat model implied mean vector
#' @param mu observed (or population) mean vector
#' @return Returns CFI
getCFI.Sigma <- function(SigmaHat, S, muHat = NULL, mu = NULL){
  checkPositiveDefinite(SigmaHat)
  checkPositiveDefinite(S)
  fm <- getF.Sigma(SigmaHat, S, muHat, mu)
  SigmaHatNull <- diag(ncol(S))
  diag(SigmaHatNull) <- diag(S)
  muHatNull <- mu # as in mplus: baseline model has unrestricted means
  f0 <- getF.Sigma(SigmaHatNull, S, muHatNull, mu)
  cfi <- (f0-fm)/f0
  cfi
}

#' getCFI.Sigma.mgroups
#'
#' Computes CFI given the model-implied and the observed (or population) covariance matrix for multiple group models.
#' `CFI = (F_null - F_hyp) / F_null` applying multiple group sampling weights to `F_hyp` and `F_null`. 
#'
#' @param SigmaHat a list of model implied covariance matrix
#' @param S a list of observed (or population) covariance matrix
#' @param muHat model implied mean vector
#' @param mu observed (or population) mean vector
#' @param N a list of group weights
#' @return Returns CFI
getCFI.Sigma.mgroups <- function(SigmaHat, S, muHat = NULL, mu = NULL, N){
  N <- unlist(N)
  
  if(is.null(mu)){
    fmin.g <- sapply(seq_along(S), function(x){getF.Sigma(SigmaHat[[x]], S[[x]])})
    fnull.g <- sapply(seq_along(S), function(x){
      SigmaHatNull <- diag(ncol(S[[x]]))
      diag(SigmaHatNull) <- diag(S[[x]])
      getF.Sigma(SigmaHatNull, S[[x]])
    })
  }else{
    fmin.g <- sapply(seq_along(S), function(x){getF.Sigma(SigmaHat[[x]], S[[x]], muHat[[x]], mu[[x]])})
    fnull.g <- sapply(seq_along(S), function(x){
      SigmaHatNull <- diag(ncol(S[[x]]))
      diag(SigmaHatNull) <- diag(S[[x]])
      muHatNull <- mu[[x]]   # as in mplus: baseline model has unrestricted means
      getF.Sigma(SigmaHatNull, S[[x]], muHatNull[[x]], mu[[x]])
    })
  }
  
  
  # approach A: apply sampling weights to CFI
  #cfi.g <- (fnull.g - fmin.g) / fnull.g
  #cfi <-  sum(cfi.g * N) / sum(N)
  
  # approach B: apply sampling weights to fmin and fnull
  fmin <- sum(fmin.g * N) / sum(N)
  fnull <- sum(fnull.g * N) / sum(N)
  cfi <- (fnull - fmin) / fnull
  
  return(cfi)
}
