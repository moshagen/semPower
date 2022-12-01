##########################  helper functions  #####################

#' semPower.genSigma
#'
#' Generate a covariance matrix and associated lavaan model strings based on defined model features.
#' This requires the lavaan package.
#' 
#' @param Phi factor correlation (or covariance) matrix or single number giving correlation between all factors or NULL for a uncorrelated factors. 
#' @param Lambda factor loading matrix. 
#' @param tau intercepts. If NULL and alpha is set, these are assumed to be zero. 
#' @param Alpha factor means. If NUll and tau is set, these are assumed to be zero. 
#' @param loadings a list providing the standardized factor loadings by factor. Must not contain secondary loadings.   
#' @param nIndicator vector indicating the number of indicators for each factor, e.g. c(4, 6) to define two factors with 4 and 6 indicators, respectively 
#' @param loadM vector giving mean loadings for each factor or single number to use for every loading
#' @param loadSD vector giving the standard deviation of loadings for each factor for use in conjunction with loadM. When NULL, SDs are set to zero.
#' @param loadMinMax list giving the minimum and maximum loading for each factor or vector to apply to all factors 
#' @param useReferenceIndicator whether to identify factors in accompanying true model string by a reference indicator (TRUE) or by setting their variance to 1 (FALSE). ()giving the minimum and maximum loading for each factor or vector to apply to all factors 
#' @param metricInvariance a list containing the factor indices for which the analysis model should apply metric invariance labels, e.g. list(c(1,2), c(3,4)) to assume invariance for f1 and f2 as well as f3 and f4. 
#' @return a list containing the implied covariance matrix (Sigma),  the implied loading (Lambda) and factor-covariance matrix (Phi), the implied indicator means (mu), intercepts (tau), and latent means (alpha), as well as the associated lavaan model string defining the population (modelPop) and a lavaan model string defining a corresponding true cfa analysis model (modelTrue) 
#' @examples
#' \dontrun{
#' # Provide factor correlation for a two-factor model, the number of indicators by factor, 
#' # and a single loading which is equal for all indicators
#' genSigma <- semPower.genSigma(phi = .2, nIndicator = c(5, 6), loadM = .5)
#' 
#' # Provide factor correlation matrix and loading matrix
#' Phi <- matrix(c(
#'                 c(1.0, 0.2, 0.5),
#'                 c(0.2, 1.0, 0.3),
#'                 c(0.5, 0.3, 1.0)
#'                ), byrow = TRUE, ncol = 3)
#' Lambda <- matrix(c(
#'                 c(0.5, 0.0, 0.0),
#'                 c(0.4, 0.0, 0.0),
#'                 c(0.3, 0.0, 0.0),
#'                 c(0.0, 0.7, 0.0),
#'                 c(0.0, 0.8, 0.0),
#'                 c(0.0, 0.5, 0.0),
#'                 c(0.0, 0.0, 0.5),
#'                 c(0.0, 0.0, 0.4),
#'                 c(0.0, 0.0, 0.6),
#'                ), byrow = TRUE, ncol = 3)
#'                
#' genSigma <- semPower.genSigma(Phi = Phi, Lambda = Lambda)
#' 
#' # same as above, but providing reduced loading matrix, i.e..
#' # only defining primary loadings; all secondary loadings are zero.
#' loadings <- list(
#'                c(0.4, 0.5, 0.8),
#'                c(0.7, 0.6, 0.5, 0.4, 0.5),
#'                c(0.5, 0.5, 0.6, 0.8)
#'                )
#'                
#' genSigma <- semPower.genSigma(Phi = Phi, loadings = loadings)
#'   
#' # same as above, but providing
#' # the number of indicators by factor 
#' # and min-max loading for all factors (sampling from a uniform)
#' genSigma <- semPower.genSigma(Phi = Phi, nIndicator = c(3, 5, 4), 
#'                               loadMinMax = c(.3, .8))
#'
#' # same as above, but providing mean and sd loading for all factors
#' genSigma <- semPower.genSigma(Phi = Phi, nIndicator = c(3, 5, 4), 
#'                               loadM = .5, loadSD = .1)
#'
#'
#' # same as above, but providing mean and sd of loadings for each factor
#' genSigma <- semPower.genSigma(Phi = Phi, nIndicator = c(3, 5, 4), 
#'                               loadM = c(.5, .6, .7), loadSD = c(0, .05, .01))
#'
#' # same as above, but using min-max loadings for each factor
#' loadMinMax <- list(
#'                    c(.4, .6),
#'                    c(.5, .8),
#'                    c(.3, .7)
#'                    )
#'
#' genSigma <- semPower.genSigma(Phi = Phi, nIndicator = c(3, 5, 4), 
#'                               loadMinMax = loadMinMax)
#'                               
#' }
#' @importFrom stats rnorm runif 
#' @importFrom utils installed.packages
#' @export
semPower.genSigma <- function(Phi = NULL, 
                              Lambda = NULL,
                              tau = NULL,
                              Alpha = NULL,  # capital Alpha, so to distinguish from alpha error
                              loadings = NULL, 
                              nIndicator = NULL, 
                              loadM = NULL, 
                              loadSD = NULL, 
                              loadMinMax = NULL,
                              useReferenceIndicator = FALSE,
                              metricInvariance = NULL,
                              ...){
  
  # check whether lavaan is available
  if(!'lavaan' %in% rownames(installed.packages())) stop('This function depends on the lavaan package, so install lavaan first.')
  
  # validate input
  if(is.null(nIndicator) & is.null(loadings) & is.null(Lambda)) stop('Either provide Labmda, loadings, or number of indicators')
  if(!is.null(nIndicator) & !is.null(loadings)) stop('Either provide loadings or number of indicators, but not both.')
  if(!is.null(nIndicator) & !is.null(Lambda)) stop('Either provide Lambda or number of indicators, but not both.')
  if(!is.null(Lambda) & !is.null(loadings)) stop('Either provide Lambda or loadings, but not both.')
  if(is.null(Lambda) & is.null(nIndicator) & !is.list(loadings)) stop('loadings must be a list')
  
  if(!is.null(Lambda)){
    nfac <- ncol(Lambda)
  }else{
    nfac <- ifelse(is.null(loadings), length(nIndicator), length(loadings))
  }
  
  if(is.null(Phi)) Phi <- diag(nfac)
  if(length(Phi) == 1){
    Phi <- matrix(Phi, ncol = nfac, nrow = nfac)
    diag(Phi) <- 1
  } 
  checkPositiveDefinite(Phi)
  if(ncol(Phi) != nfac) stop('Phi must have the same number of rows/columns as the number of factors.') 

  if(is.null(loadings) && is.null(Lambda)){
    if(any(!sapply(nIndicator, function(x) x %% 1 == 0))) stop('Number of indicators must be a integer')
    invisible(sapply(nIndicator, function(x) checkBounded(x, 'Number of indicators ', bound = c(1, 10000), inclusive = TRUE)))
    if(is.null(loadM) & is.null(loadMinMax)) stop('Either mean loading or min-max loading need to be defined')
    if(is.null(loadMinMax) & length(loadM) == 1) loadM <- rep(loadM, nfac)
    if(is.null(loadMinMax) & length(loadM) != nfac) stop('Nindicator and mean loading must of same size')
    if(is.null(loadMinMax)) invisible(sapply(loadM, function(x) checkBounded(x, 'All loadings', bound = c(-1, 1), inclusive = TRUE)))
    
    if(!is.null(loadMinMax) && (!is.null(loadM) || !is.null(loadSD))) stop('Either specify mean and SD of loadings or specify min-max loading, both not both.')
    if(!is.null(loadMinMax)) invisible(sapply(unlist(loadMinMax), function(x) checkBounded(x, 'All loadings', bound = c(-1, 1), inclusive = TRUE)))
    if(length(loadMinMax) == 2) loadMinMax <- lapply(1:nfac, function(x) loadMinMax)
    if(is.null(loadMinMax) && is.null(loadSD)) loadSD <- rep(0, nfac)
    if(is.null(loadMinMax) && length(loadSD) == 1) loadSD <- rep(loadSD,nfac)
    if(is.null(loadMinMax) && !is.null(loadSD)) invisible(sapply(loadSD, function(x) checkBounded(x, 'Standard deviations', bound = c(0, .5), inclusive = TRUE)))
  }else{
    if(is.null(Lambda)){
      invisible(lapply(loadings, function(x) lapply(x, function(x) checkBounded(x, 'All loadings', bound = c(-1, 1), inclusive = TRUE))))
      nIndicator <- unlist(lapply(loadings, length))  # crossloadings are disallowed
    }else{
      invisible(apply(Lambda, c(1, 2), function(x) checkBounded(x, 'All loadings', bound = c(-1, 1), inclusive = TRUE)))
      if(any(apply(Lambda, 1, function(x) sum(x^2)) > 1)) stop('Loadings imply negative residual variance(s). Note that loadings must be standardized.')
      nIndicator <- apply(Lambda, 2, function(x) sum(x != 0))
    }
  }
  
  if(!is.null(tau)){
    if(length(tau) != sum(nIndicator)) stop('Intercepts (tau) must be of same length as the number of indicators')
  }
  if(!is.null(Alpha)){
    if(length(Alpha) != nfac) stop('Latent means (Alpha) must be of same length as the number of factors')
  }
  if(!is.null(tau) & is.null(Alpha)) Alpha <- rep(0, nfac)
  if(!is.null(Alpha) & is.null(tau)) tau <- rep(0, sum(nIndicator))
  
  if(!is.null(metricInvariance)){
    if(!is.list(metricInvariance)) stop('metricInvariance must be a list')
    if(any(unlist(lapply(metricInvariance, function(x) length(x))) < 2)) stop('each list entry in metricInvariance must involve at least two factors')
    if(max(unlist(metricInvariance)) > nfac | min(unlist(metricInvariance)) <= 0) stop('factor index < 1 or > nfactors in metricInvariance')
    if(any(unlist(lapply(metricInvariance, function(x) var(nIndicator[x]))) != 0)) stop('factors in metriInvariance must have the same number of indicators')
    metricInvarianceLabels <- lapply(1:length(metricInvariance), function(x) paste0('l', x, 1:nIndicator[metricInvariance[[x]][1]], '*')) 
  }


  # define lambda 
  if(is.null(Lambda)){
    Lambda <- matrix(0, ncol = nfac, nrow = sum(nIndicator)) # store loading matrix
    sidx <- 1
    for(f in 1:nfac){
      eidx <- sidx + (nIndicator[f] - 1)
      if(!is.null(loadM)){
        cload <- round(rnorm(nIndicator[f], loadM[f], loadSD[f]), 2)
        if(any(cload < -1) | any(cload > 1)) warning('Sampled loadings outside [-1, 1] were set to -1 or 1.')
        cload[cload < -1] <- -1
        cload[cload > 1] <- 1
      }else if(!is.null(loadMinMax)){
        cload <- round(runif(nIndicator[f], loadMinMax[[f]][1], loadMinMax[[f]][2]), 2)
      }else if(!is.null(loadings)){
        cload <- loadings[[f]]  # loadings only contains primary loadings
      }else{
        stop('loading not found')  
      }
      Lambda[sidx:eidx, f] <- cload
      sidx <- eidx + 1
    }
  }
  
  # define theta
  dTheta <- diag(1 - Lambda %*% Phi %*% t(Lambda))
  
  ### we could now do Sigma  = Lambda %*% Phi %*% t(Lambda) + diag(Theta), 
  ### but we let lav do this anyway, also because it's nice to have a 
  ### population model string
  tok <- list()
  for(f in 1:ncol(Lambda)){
    iIdx <- which(Lambda[, f] != 0)
    cload <- Lambda[iIdx, f]
    tok[f] <- 
      paste(
        paste0('f', f, ' =~ ', paste0(cload, '*', paste0('x', iIdx), collapse = ' + ')),
        paste0('f', f, ' ~~ ', Phi[f, f],'*', 'f', f),
        sep='\n')
  }
  # residuals
  for(i in 1:nrow(Lambda)){
    tok <- append(tok, paste0(paste0('x', i), ' ~~ ', dTheta[i], '*', paste0('x', i), collapse = '\n'))
  }
  # define factor cor
  if(nfac > 1){
    for(f in 1:(nfac - 1)){
      for(ff in (f + 1):nfac){
        tok <- append(tok, paste0('f', f, ' ~~ ', Phi[f, ff], '*f', ff))
      }
    }
  }
  
  # add means
  if(!is.null(tau)){
    tok <- append(tok, paste0(paste0('x', 1:nrow(Lambda)), ' ~ ', tau, '*1', collapse = '\n'))
  }
  if(!is.null(Alpha)){
    tok <- append(tok, paste0(paste0('f', 1:ncol(Lambda)), ' ~ ', Alpha, '*1', collapse = '\n'))
  }
  
  modelPop <- paste(unlist(tok), collapse = '\n')

  # get Sigma (and mu). 
  Sigma <- lavaan::fitted(lavaan::sem(modelPop))$cov
  # checkPositiveDefinite(Sigma) # redundant, lav checks this
  mu <- NULL
  if(!is.null(tau)) mu <- lavaan::fitted(lavaan::sem(modelPop))$mean
  
  
  ### also not needed, but since we are at in anyway, we 
  ### can also provide a (correct) cfa analysis model string 
  tok <- list()
  sidx <- 1
  for(f in 1:ncol(Lambda)){
    iIdx <- which(Lambda[, f] != 0)
    # add invariance constrains
    if(any(unlist(lapply(metricInvariance, function(x) f %in% x)))){
      labelIdx <- which(unlist(lapply(metricInvariance, function(x) f %in% x)))
      clabel <- metricInvarianceLabels[[labelIdx]]
      if(useReferenceIndicator){
        # scale by first loading instead of 1
        tok[f] <- paste0('f', f, ' =~ ', Lambda[iIdx[1], f], '*x', iIdx[1], ' + ', paste0(clabel, 'x', iIdx, collapse = ' + '))
      }else{
        tok[f] <- paste0('f', f, ' =~ NA*x', iIdx[1],' + ', paste0(clabel, 'x', iIdx, collapse = ' + '))
      }
    }else{
      if(useReferenceIndicator){
        # scale by first loading instead of 1
        tok[f] <- paste0('f', f, ' =~ ', Lambda[iIdx[1], f], '*', paste0('x', iIdx, collapse = ' + '))
      }else{
        tok[f] <- paste0('f', f, ' =~ NA*', paste0('x', iIdx, collapse = ' + '))
      }
    }
  }
  modelTrue <- paste(c(unlist(tok)), collapse = '\n')
  if(!useReferenceIndicator){
    modelTrue <- paste(c(modelTrue,
                       # factor variances are always 1, regardless of phi
                       sapply(1:nfac, function(f) paste0('f', f, ' ~~ 1*f', f))),
                       collapse = '\n')
  }

  
  list(Sigma = Sigma, 
       mu = mu,
       Lambda = Lambda, 
       Phi = Phi,
       tau = tau,
       Alpha = Alpha,
       modelPop = modelPop, 
       modelTrue = modelTrue)
}


#' getPhi.B
#'
#' Computes implied correlations from Beta matrix (using all-y notation), disallowing recursive paths.
#' 
#' @param B matrix of regression coefficients (all-y notation). Must only contain non-zero lower-triangular elements, so the first row only includes zeros. 
#' @param lPsi matrix of residual correlations. This is not the Psi matrix, but a lesser version ignoring all variances and containing correlations (when standardized = TRUE) off the diagonal. Can be omitted for no correlations beyond those implied by B. 
#' @param standardized whether B and lPsi shall be interpreted as standardized coefficients (TRUE) or as unstandardized coefficients (FALSE). 
#' @return the implied correlation matrix
#' @examples
#' \dontrun{
#' # mediation model
#' B <- matrix(c(
#' c(.00, .00, .00),
#' c(.10, .00, .00),
#' c(.20, .30, .00)
#' ), byrow = TRUE, ncol = 3)
#' Phi <- getPhi.B(B)
#' 
#' # clpm type model with residual correlations at wave 2 + 3
#' B <- matrix(c(
#'   c(.0, .0, .0, .0, 0, 0),  # X1
#'   c(.0, .0, .0, .0, 0, 0),  # Y1
#'   c(.7, .1, .0, .0, 0, 0),  # X2
#'   c(.2, .8, .0, .0, 0, 0),  # Y2
#'   c(.0, .0, .7, .1, 0, 0),  # X3
#'   c(.0, .0, .2, .8, 0, 0)   # Y3
#' ), byrow = TRUE, ncol = 6)
# 
#' lPsi <- matrix(0, ncol = ncol(B), nrow = nrow(B))
#' lPsi[3,4] <- lPsi[4,3] <- .2
#' lPsi[5,6] <- lPsi[6,5] <- .3
#' 
#' Phi <- getPhi.B(B, lPsi, standardized = FALSE)
#' 
#' }
getPhi.B <- function(B, lPsi = NULL, standardized = TRUE){
  
  checkSquare(B)
  if(any(B[upper.tri(B, diag = TRUE)] != 0)) stop('B may not contain any non-zero values on or upper the diagonal.')
  if(any(rowSums(B^2) > 1)) stop('B implies negative residual variances.')
  invisible(lapply(B, function(x) lapply(x, function(x) checkBounded(x, 'All elements in B', bound = c(-1, 1)))))
  if(!is.null(lPsi)){
    checkSymmetricSquare(lPsi)
    if(ncol(lPsi) != ncol(B)) stop('lPsi must be of same dimension as B')
    invisible(lapply(lPsi, function(x) lapply(x, function(x) checkBounded(x, 'All elements in lPsi', bound = c(-1, 1)))))
    diag(lPsi) <- 0
  }else{
    lPsi <- diag(ncol(B))
  }
  
  if(!standardized){
    
    ## this treats B and lPsi as unstandardized parameters
    ## and thus yields Phi as variance/covariance matrix
    invIB <- solve(diag(ncol(B)) - B)
    Psi <- diag(ncol(B)) + lPsi 
    Phi <- invIB %*% Psi %*% t(invIB) 

  }else{
    
    ## for std, exploit the structure of B to build Phi recursively
    ## there must be a simpler way to do this...
    exog <- apply(B, 1, function(x) !any(x != 0))
    Be <- B[!exog, ]
    if(!is.matrix(Be)) Be <- t(matrix(Be))
    
    Phi <- diag(ncol(B)) 
    for(i in 1:nrow(Be)){
      idx <- i + sum(exog) - 1
      cb <- matrix(Be[i, 1:idx])
      cr <- Phi[1:idx, 1:idx]
      predR <- (cr %*% cb)
      # add residual covariances
      if(!is.null(lPsi) & any(lPsi[idx, 1:(idx + 1)] != 0)){
        cR <- rbind(cr, t(predR))
        cR <- cbind(cR, t(t(c((predR),1))))
        cB <- B[1:(idx + 1), 1:(idx + 1)]
        rootResidualVar <-  sqrt(diag(1 - diag(cB %*% cR %*% t(cB)))) 
        cPsi <- lPsi[1:(idx + 1), 1:(idx + 1)]
        corR <- cR + rootResidualVar %*% cPsi %*% t(rootResidualVar) 
        predR <- corR[(idx + 1), 1:idx]
      }
      Phi[(idx + 1), 1:idx] <- t(predR)
      Phi[1:idx, (idx + 1)] <- predR
    }
    
  }

  checkPositiveDefinite(Phi)
  
  Phi
}


#' semPower.getDf
#'
#' Convenience function to determine the degrees of freedom of a given model provided as lavaan model string. 
#' This requires the lavaan package.
#' 
#' @param lavModel the lavaan model string 
#' @param nGroups for multigroup models: the number of groups 
#' @param group.equal for multigroup models: type of group equality constraints (loadings, intercepts, means, residuals, residual.covariances, lv.variances, lv.covariances, regressions)
#' @return df
#' @examples
#' \dontrun{
#' lavModel <- '
#' f1 =~ x1 + x2 + x3 + x4
#' f2 =~ x5 + x6 + x7 + x8
#' f3 =~ y1 + y2 + y3
#' f3 ~ f1 + f2
#' '
#' semPower.getDf(lavModel)
#' 
#' # multigroup version
#' semPower.getDf(lavModel, nGroups = 3)  
#' semPower.getDf(lavModel, nGroups = 3, group.equal = c('loadings'))
#' semPower.getDf(lavModel, nGroups = 3, group.equal = c('loadings', 'intercepts'))
#' }
#' @importFrom utils installed.packages
#' @export
semPower.getDf <- function(lavModel, nGroups = NULL, group.equal = NULL){
  # check whether lavaan is available
  if(!'lavaan' %in% rownames(installed.packages())) stop('This function depends on the lavaan package, so install lavaan first.')
  
  # Fit model to dummy covariance matrix instead of counting parameters. 
  # This should also account for parameter restrictions and other intricacies
  # Model fitting will give warnings we just can ignore
  tryCatch({
    params <- lavaan::sem(lavModel)
    dummyS <- diag(params@Model@nvar)
    rownames(dummyS) <- params@Model@dimNames[[1]][[1]]
    if(is.null(nGroups) || nGroups == 1){
      dummyFit <- suppressWarnings(lavaan::sem(lavModel, sample.cov = dummyS, sample.nobs = 1000, warn = FALSE))
    }else{
      if(is.null(group.equal)){
        dummyFit <- suppressWarnings(lavaan::sem(lavModel, sample.cov = lapply(1:nGroups, function(x) dummyS), sample.nobs = rep(1000, nGroups), warn = FALSE))
      }else{
        dummyFit <- suppressWarnings(lavaan::sem(lavModel, sample.cov = lapply(1:nGroups, function(x) dummyS), sample.nobs = rep(1000, nGroups), group.equal = group.equal, warn = FALSE))
      }
    }
    dummyFit@test$standard$df
  }, 
  warning = function(w){
    warning(w)
  }, 
  error = function(e){
    stop(e)
  }
  )
}