##########################  determine ncp, chi-square from Fmin   #####################


#' getNCP
#'
#' calculates non-centrality parameter from the population minimum of the fit-function
#'
#' ncp = (n-1) * F
#'
#' @param Fmin population minimum of the fit-function
#' @param n number of observations
#' @return NCP
getNCP <- function(Fmin, n){
  NCP <- unlist(Fmin) * (unlist(n)-1)
  sum(NCP)
}


#' getChiSquare.NCP
#'
#' calculates chi-square from NCP
#'
#' chi = ncp + df
#'
#' @param NCP non-centrality parameter
#' @param df model degrees of freedom
#' @return chiSquare
getChiSquare.NCP <- function(NCP, df){
  chiSquare <- NCP + df
  chiSquare
}

#' getChiSquare.F
#'
#' calculates chis-square from the population minimum of the fit-function
#'
#' chi = (n-1)*F +  df = ncp + df
#'
#' note that F is the population minimum; using F_hat would give chi = (n-1)*F_hat
#'
#' @param Fmin population minimum of the fit-function
#' @param n number of observations
#' @param df model degrees of freedom
#' @return NCP
getChiSquare.F <- function(Fmin, n, df){
  chiSquare <- getNCP(Fmin, n) + df
  chiSquare
}


##########################  determine Fmin from RMSEA , Mc , GFI , AGFI   #####################


#' getF
#' calculates minimum of the ML-fit-function from known fit indices
#' @param effect magnitude of effect
#' @param effect.measure measure of effect, one of 'fmin','rmsea','agfi','gfi','mc'
#' @param df model degrees of freedom
#' @param p number of observed varaibles
#' @param SigmaHat model implied covariance matrix
#' @param Sigma population covariance matrix
#' @return Fmin
getF <- function(effect, effect.measure, df = NULL, p = NULL, SigmaHat = NULL, Sigma = NULL){
  fmin <- effect
  if(is.null(SigmaHat)){ # sufficient to check for on NULL matrix; primary validity check is in validateInput
    switch(effect.measure,
           "RMSEA" = fmin <- getF.RMSEA(effect, df),
           "MC" = fmin <- getF.Mc(effect),
           "GFI" = fmin <- getF.GFI(effect, p),
           "AGFI" = fmin <- getF.AGFI(effect, df, p)
    )
  }else{
    fmin <- effect <- getF.Sigma(SigmaHat, Sigma)
  }
  fmin
}


#' getF.RMSEA
#'
#' calculates minimum of the ML-fit-function from RMSEA
#'
#' F_min = rmsea^2 * df
#'
#' @param RMSEA RMSEA
#' @param df model degrees of freedom
#' @return Fmin
getF.RMSEA <- function(RMSEA, df){
  fmin <- RMSEA^2 * df
  fmin
}


#' getF.Mc
#'
#' calculates minimum of the ML-fit-function from Mc
#'
#'
#' @param Mc Mc
#' @return Fmin
getF.Mc <- function(Mc){
  fmin <- -2 * log(Mc)
  fmin
}

#' getF.GFI
#'
#' calculates minimum of the ML-fit-function from AGFI
#'
#'
#' @param GFI GFI
#' @param p number of observed variables
#' @return Fmin
getF.GFI <- function(GFI, p){
  fmin <- -(GFI-1)*p/(2*GFI)
  fmin
}


#' getF.AGFI
#'
#' calculates minimum of the ML-fit-function from AGFI
#'
#' F_min = rmsea^2 * df
#'
#' @param AGFI AGFI
#' @param df model degrees of freedom
#' @param p number of observed variables
#' @return Fmin
getF.AGFI <- function(AGFI, df, p){
  fmin <- -(AGFI-1)*df*p/(p*p+p+(2*AGFI-2)*df)
  fmin
}


##########################  determine RMSEA Mc GFI AGFI from Fmin   #####################

#' getIndices.F
#'
#' calculates known indices from minimum of the ML-fit-function
#'
#'
#' @param fmin minimum of the ML-fit-function
#' @param df model degrees of freedom
#' @param p number of observed variables
#' @param SigmaHat model implied covariance matrix
#' @param Sigma population covariance matrix
#' @param N list of sample weights
#' @return list of indices
getIndices.F <- function(fmin, df, p = NULL, SigmaHat = NULL, Sigma = NULL, N = NULL){
  fit <- list(
    rmsea = getRMSEA.F(fmin, df, nGroups=ifelse(length(N)>1, length(N), 1)),
    mc = getMc.F(fmin),
    gfi = NULL,
    agfi = NULL
  )
  if(!is.null(p)){
    fit$gfi <- getGFI.F(fmin, p)
    fit$agfi <- getAGFI.F(fmin, df, p)
  }
  if(!is.null(SigmaHat)){
    if(length(N) > 1){
      fit$srmr <- getSRMR.Sigma.mgroups(SigmaHat, Sigma, N)
      fit$cfi <- getCFI.Sigma.mgroups(SigmaHat, Sigma, N)
    }else{
      fit$srmr <- getSRMR.Sigma(SigmaHat, Sigma)
      fit$cfi <- getCFI.Sigma(SigmaHat, Sigma)
    }
  }
  fit
}



#' getRMSEA.F
#'
#' calculates RMSEA from minimum of the ML-fit-function
#'
#' F_min = rmsea^2 * df
#'
#' @param Fmin minimum of the ML-fit-function
#' @param df model degrees of freedom
#' @param nGroups the number of groups
#' @return RMSEA
getRMSEA.F <- function(Fmin, df, nGroups = 1){
  RMSEA <- sqrt(Fmin/df) * sqrt(nGroups)
  RMSEA
}

#' getMc.F
#'
#' calculates Mc from minimum of the ML-fit-function
#'
#'
#' @param Fmin minimum of the ML-fit-function
#' @return Mc
getMc.F <- function(Fmin){
  Mc <- exp(-.5*Fmin)
  Mc
}


#' getGFI.F
#'
#' calculates GFI from minimum of the ML-fit-function
#'
#'
#' @param Fmin minimum of the ML-fit-function
#' @param p number of observed variables
#' @return GFI
getGFI.F <- function(Fmin, p){
  GFI <- p/(p + 2*Fmin)
  GFI
}


#' getAGFI.F
#'
#' calculates AGFI from minimum of the ML-fit-function
#'
#'
#' @param Fmin minimum of the ML-fit-function
#' @param df model degrees of freedom
#' @param p number of observed variables
#' @return AGFI
getAGFI.F <- function(Fmin, df, p){
  AGFI <- -(Fmin*p*p+(Fmin-df)*p-2*df*Fmin)/(df*p+2*df*Fmin)
  AGFI
}



##########################  calculate Fmin RMSEA SRMR CFI from covariance matrix #####################


#' getF.Sigma
#'
#' calculates minimum of the ML-fit-function given model-implied and observed covariance matrix.
#'
#' F_min = tr(S %*% SigmaHat^-1) - p + ln(det(SigmaHat)) - ln(det(S))
#'
#' @param SigmaHat model implied covariance matrix
#' @param S observed (or population) covariance matrix
#' @return Fmin
getF.Sigma <- function(SigmaHat, S){
  checkPositiveDefinite(SigmaHat)
  checkPositiveDefinite(S)
  fmin <- sum(diag(S %*% solve(SigmaHat))) + log(det(SigmaHat)) - log(det(S)) - ncol(S)
  fmin
}

#' getSRMR.Sigma
#'
#' calculates SRMR given model-implied and observed covariance matrix.
#'
#'
#' @param SigmaHat model implied covariance matrix
#' @param S observed (or population) covariance matrix
#' @return SRMR
getSRMR.Sigma <- function(SigmaHat, S){
  checkPositiveDefinite(SigmaHat)
  checkPositiveDefinite(S)
  p <- ncol(S)

  # bollen approach to standardization
  # m <- cov2cor(S) - cov2cor(SigmaHat)
  
  # hu+bentler approach to std
  sqrt.d <- 1/sqrt(diag(S))
  D <- diag(sqrt.d, ncol=length(sqrt.d))
  m <- D %*% (S - SigmaHat) %*% D
  
  fols <- sum(m[lower.tri(m, diag=T)]^2)
  # mplus variant
  #fols <- (sum(m[lower.tri(m, diag=F)]^2)  +  sum(((diag(S) - diag(SigmaHat))/diag(S))^2)) 
  
  srmr <- sqrt(fols / (p*(p+1)/2))
  srmr
}

#' getSRMR.Sigma.mgroups 
#'
#' calculates SRMR given model-implied and observed covariance matrix for multiple group models
#'
#' @param SigmaHat a list of model implied covariance matrices
#' @param S a list of observed (or population) covariance matrices
#' @param N a list of group weights
#' @return SRMR
getSRMR.Sigma.mgroups <- function(SigmaHat, S, N){
  srmrs <- sapply(seq_along(SigmaHat), function(x) getSRMR.Sigma(SigmaHat[[x]], S[[x]]))
  # lavaan approach: apply sample weights to srmr
  # srmr <- (sum(srmrs*N)/sum(N))
  # mplus approach: apply sample weights to squared sums of res
  srmr <- sqrt(sum(unlist(srmrs)^2*unlist(N))/sum(unlist(N)))
  srmr
}




#' getCFI.Sigma
#'
#' calculates CFI given model-implied and observed covariance matrix.
#'
#' cfi= (f_null - f_hyp) / f_null
#'
#' @param SigmaHat model implied covariance matrix
#' @param S observed (or population) covariance matrix
#' @return CFI
getCFI.Sigma <- function(SigmaHat, S){
  checkPositiveDefinite(SigmaHat)
  checkPositiveDefinite(S)
  fm <- getF.Sigma(SigmaHat, S)
  SigmaHatNull <- diag(ncol(S))
  diag(SigmaHatNull) <- diag(S)
  f0 <- getF.Sigma(SigmaHatNull, S)
  cfi <- (f0-fm)/f0
  cfi
}

#' getCFI.Sigma.mgroups
#'
#' calculates CFI given model-implied and observed covariance matrix for multiple group models.
#'
#' cfi= (f_null - f_hyp) / f_null
#'
#' @param SigmaHat a list of model implied covariance matrix
#' @param S a list of observed (or population) covariance matrix
#' @param N a list of group weights
#' @return CFI
getCFI.Sigma.mgroups <- function(SigmaHat, S, N){

  fmin.g <- sapply(seq_along(S), function(x){getF.Sigma(SigmaHat[[x]], S[[x]])})
  fmin <- sum(unlist(fmin.g) * (unlist(N)-1))
  
  fnull.g <- sapply(seq_along(S), function(x){
    SigmaHatNull <- diag(ncol(S[[x]]))
    diag(SigmaHatNull) <- diag(S[[x]])
    getF.Sigma(SigmaHatNull, S[[x]])
    })
  fnull <- sum(unlist(fnull.g) * (unlist(N)-1))
  
  cfi <- (fnull - fmin)/fnull
  return(cfi)
}



##########################  output and formatting #####################


#' getFormattedResults
#'
#' returned dataframe containing formatted results
#' 
#' @param type type of power analysis
#' @param result result object (list)
#' @param digits number of significant digits
#' @return data.frame
getFormattedResults <- function(type, result, digits = 6){

  ########### common across power types

  if(!is.null(result$srmr) && !is.null(result$gfi)){
    rows <- c('F0','RMSEA','SRMR','Mc','GFI','AGFI','CFI', '')
    head <- data.frame(rows)
    values <- c(result$fmin, result$rmsea, result$srmr, result$mc, result$gfi, result$agfi, result$cfi)
  }else if(!is.null(result$gfi)){
    rows <- c('F0','RMSEA','Mc','GFI','AGFI', '')
    head <- data.frame(rows)
    values <- c(result$fmin, result$rmsea, result$mc, result$gfi, result$agfi)
  }else{
    rows <- c('F0','RMSEA','Mc', '')
    head <- data.frame(rows)
    values <- c(result$fmin, result$rmsea, result$mc)
  }

  head$values <- c(formatC(values, format='f', digits = digits),'')


  ########### a-priori

  if(type == 'a-priori'){

    rows <- c('df','Required Num Observations','')
    body <- data.frame(rows)
    body$values <- c(result$df, result$requiredN, '')

    rows <- c('Critical Chi-Square', 'NCP', 'Alpha', 'Beta', 'Power (1-beta)', 'Implied Alpha/Beta Ratio')
    foot <- data.frame(rows)

    v1 <- formatC(c(result$chiCrit, result$impliedNCP), format = 'f', digits = digits)
    v1 <- sapply(v1, substr, 1, digits+2)

    v2 <- c(result$alpha, result$impliedBeta, result$impliedPower, result$impliedAbratio)
    # determine whether to use float or scientific number format
    v2.f <- rep('f', length(v2))
    v2.f[v2 < 1e-5 | v2 > 1e5] <- 'e'
    v2 <- sapply(seq_along(v2),  function(y, z, i) { formatC(x=v2[i], format=v2.f[i], digits = digits)}, y=v2, z= v2.f)
    foot$values <- c(v1, v2)

    # manually correct some quirks
    if(result$impliedBeta < 1e-5){
      foot$values[foot$rows=='Power (1-beta)'] <- '> 0.9999'
    }


    out <- rbind(head, body, foot)
    rownames(out) <- colnames(out) <- NULL

  }


  ########### post-hoc

  if(type == 'post-hoc'){

    ifelse(!is.list(result$N), rows <- c('df','Num Observations'), rows <- c('df','Num Observations',' '))
    body <- data.frame(rows)
    ifelse(!is.list(result$N), 
           body$values <- c(result$df, result$N), 
           body$values <- c(result$df, sum(unlist(result$N)), paste0('(',paste(result$N, collapse = ', '),')'))
          )

    rows <- c('NCP', '', 'Critical Chi-Square')
    body2 <- data.frame(rows)
    v1 <- c(formatC(result$ncp, format = 'f', digits = digits), '',
                      formatC(result$chiCrit, format = 'f', digits = digits)
                      )
    body2$values <- sapply(v1, substr, 1, (digits+2))

    rows <- c('Alpha', 'Beta', 'Power (1-beta)', 'Implied Alpha/Beta Ratio')
    foot <- data.frame(rows)
    v <- c(result$alpha, result$beta, result$power, result$impliedAbratio)
    # determine whether to use float or scientific number format
    v.f <- rep('f', length(v))
    v.f[v < 1e-5 | v > 1e5] <- 'e'
    foot$values <- sapply(seq_along(v),  function(y, z, i) { formatC(x=v[i], format=v.f[i], digits = digits)}, y=v, z= v.f)

    # manually correct some quirks
    if(result$beta < 1e-5){
      foot$values[foot$rows=='Power (1-beta)'] <- '> 0.9999'
    }
    if(result$beta == 0){
      foot$values[foot$rows=='Beta'] <- '< 1.00e-320'
      foot$values[foot$rows=='Implied Alpha/Beta Ratio'] <- '> 1.00e-320'
    }

    out <- rbind(head, body, body2, foot)
    rownames(out) <- colnames(out) <- NULL

  }

  ########### compromise

  if(type == 'compromise'){

    rows <- c('df','Num Observations', 'Desired Alpha/Beta Ratio', '','Critical Chi-Square')
    body <- data.frame(rows)
    if(!result$bPrecisionWarning){
      body$values <- c(result$df, result$N, formatC(result$desiredAbratio, format = 'f', digits = digits), '',
                       substr(formatC(result$chiCrit, format = 'f', digits = digits), 1, digits+2))
    }else{
      smax <- substr(formatC(result$max, format = 'f', digits = digits), 1, digits+2)
      smin <- substr(formatC(result$min, format = 'f', digits = digits), 1, digits+2)
      sChiCrit <- paste(smax,'< Chi-Square < ', smin)
      body$values <- c(result$df, result$N, formatC(result$desiredAbratio, format = 'f', digits = digits), '',
                       sChiCrit)
    }

    rows <- c('Implied Alpha', 'Implied Beta', 'Implied Power (1-beta)', 'Actual Alpha/Beta Ratio')
    foot <- data.frame(rows)
    v <- c(result$impliedAlpha, result$impliedBeta, result$impliedPower, result$impliedAbratio)
    # determine whether to use float or scientific number format
    v.f <- rep('f', length(v))
    v.f[v < 1e-5 | v > 1e5] <- 'e'
    foot$values <- sapply(seq_along(v),  function(y, z, i) { formatC(x=v[i], format=v.f[i], digits = digits)}, y=v, z= v.f)

    # manually correct some quirks
    if(result$impliedBeta < 1e-5){
      foot$values[foot$rows=='Implied Power (1-beta)'] <- '> 0.9999'
    }
    if(result$bPrecisionWarning){
      foot$values[foot$rows=='Implied Beta'] <- '< 1.00e-240'
      foot$values[foot$rows=='Implied Alpha'] <- '< 1.00e-320'
      foot$values[foot$rows=='Actual Alpha/Beta Ratio'] <- ' '
    }

    out <- rbind(head, body, foot)
    rownames(out) <- colnames(out) <- NULL


  }

  out
}

