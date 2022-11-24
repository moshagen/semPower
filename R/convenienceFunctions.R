
#' semPower.genSigma
#'
#' Convenience function to generate a covariance matrix and associated lavaan model strings based on defined model features.
#' This requires the lavaan package.
#' 
#' @param Phi factor correlation matrix or single number giving correlation between all factors or NULL for a uncorrelated factors. 
#' @param Lambda factor loading matrix (standardized). 
#' @param tau intercepts. If NULL and alpha is set, these are assumed to be zero. 
#' @param Alpha factor means. If NUll and tau is set, these are assumed to be zero. 
#' @param loadings a list providing the standardized factor loadings by factor. Must not contain secondary loadings.   
#' @param nIndicator vector indicating the number of indicators for each factor, e.g. c(4, 6) to define two factors with 4 and 6 indicators, respectively 
#' @param loadM vector giving mean loadings for each factor or single number to use for every loading
#' @param loadSD vector giving the standard deviation of loadings for each factor for use in conjunction with loadM. When NULL, SDs are set to zero.
#' @param loadMinMax list giving the minimum and maximum loading for each factor or vector to apply to all factors 
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
                              ...){
  
  # check whether lavaan is available
  if(!'lavaan' %in% rownames(installed.packages())) stop('This function depends on the lavaan package, so install lavaan first.')
  
  # validate input
  if(is.null(nIndicator) & is.null(loadings) & is.null(Lambda)) stop('Either provide Labmda, loadings, or number of indicators')
  if(!is.null(nIndicator) & !is.null(loadings)) stop('Either provide loadings or number of indicators, but not both.')
  if(!is.null(nIndicator) & !is.null(Lambda)) stop('Either provide Lambda or number of indicators, but not both.')
  if(!is.null(Lambda) & !is.null(loadings)) stop('Either provide Lambda or loadings, but not both.')
  if(is.null(nIndicator) & !is.list(loadings)) stop('loadings must be a list')
  
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
  invisible(apply(Phi, c(1, 2), function(x) checkBounded(x, 'All elements in Phi', bound = c(-1, 1), inclusive = TRUE)))
  
  if(is.null(loadings)){
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
      invisible(lapply(loadings, function(x) lapply(x, function(x) checkBounded(x, 'All loadings', bound = c(-1, 1)))))
      nIndicator <- unlist(lapply(loadings, length))  # crossloadings are disallowed
    }else{
      invisible(apply(Lambda, c(1, 2), function(x) checkBounded(x, 'All loadings', bound = c(-1, 1))))
      if(any(apply(Lambda, 1, function(x) sum(x^2)) > 1)) stop('Loadings imply negative residual variance(s). Note that loadings must be standardized.')
      nIndicator <- apply(Lambda, 2, function(x) length(x))
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
  
  
  ### pop model
  # define factors
  tok <- list()
  iLambda <- matrix(0, ncol = nfac, nrow = sum(nIndicator)) # store loading matrix
  sidx <- 1
  for(f in 1:nfac){
    eidx <- sidx + (nIndicator[f] - 1)
    if(!is.null(loadM)){
      cload <- round(rnorm(nIndicator[f], loadM[f], loadSD[f]), 2)
      if(any(cload <= -1) | any(cload >= 1)) warning('Sampled loadings outside [-1, 1] were set to -.99 or .99.')
      cload[cload <= -1] <- -.99
      cload[cload >= 1] <- .99
    }else if(!is.null(loadMinMax)){
      cload <- round(runif(nIndicator[f], loadMinMax[[f]][1], loadMinMax[[f]][2]), 2)
    }else if(!is.null(loadings)){
      cload <- loadings[[f]]  # loadings only contains primary loadings
    }else{
      cload <- Lambda[sidx:eidx, f]  # loadings only contains primary loadings
    }
    iLambda[sidx:eidx, f] <- cload
    tok[f] <- 
      paste(
        paste0('f', f, ' =~ ', paste0(cload, '*', paste0('x', sidx:eidx), collapse = ' + ')),
        paste0('f', f, ' ~~ 1*', 'f', f),
        paste0(paste0('x', sidx:eidx), ' ~~ ', 1 - cload^2, '*', paste0('x', sidx:eidx), collapse = '\n'),
        sep='\n')
    sidx <- eidx + 1
  }
  
  # define factor cor
  for(f in 1:(nfac - 1)){
    for(ff in (f + 1):nfac){
      tok <- append(tok, paste0('f', f, ' ~~ ', Phi[f, ff], '*f', ff))
    }
  }
  
  # add means
  if(!is.null(tau)){
    tok <- append(tok, paste0(paste0('x', 1:sum(nIndicator)), ' ~ ', tau, '*1', collapse = '\n'))
  }
  if(!is.null(Alpha)){
    tok <- append(tok, paste0(paste0('f', 1:nfac), ' ~ ', Alpha, '*1', collapse = '\n'))
  }
  
  modelPop <- paste(unlist(tok), collapse = '\n')
  
  
  ### create true cfa analysis model string 
  tok <- list()
  sidx <- 1
  for(f in 1:nfac){
    eidx <- sidx + (nIndicator[f] - 1)
    tok[f] <- paste0('f', f, ' =~ NA*', paste0('x', sidx:eidx, collapse = ' + '))
    sidx <- eidx + 1
  }
  modelTrue <- paste(c(
    unlist(tok),
    sapply(1:nfac, function(x) paste0('f', x, ' ~~ 1*f', x))), 
    collapse = '\n')
  
  
  # get Sigma (and mu)
  Sigma <- lavaan::fitted(lavaan::sem(modelPop))$cov
  mu <- NULL
  if(!is.null(tau)) mu <- lavaan::fitted(lavaan::sem(modelPop))$mean
  
  list(Sigma = Sigma, 
       mu = mu,
       Lambda = iLambda, 
       Phi = Phi,
       tau = tau,
       Alpha = Alpha,
       modelPop = modelPop, 
       modelTrue = modelTrue)
  
}


#' semPower.powerLav
#'
#' Perform power analysis on population and model-implied Sigmas as defined through lavaan model strings
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of 'a-priori', 'post-hoc', 'compromise'
#' @param modelPop lavaan model string defining the true model.
#' @param modelH0 lavaan model string defining the (incorrect) analysis model
#' @param modelH1 lavaan model string defining the comparison model. If omitted, the saturated model is the comparison model.
#' @param Sigma population covariance matrix (if modelPop is not set).
#' @param mu population means (if modelPop is not set).
#' @param ... other parameters related to the specific type of power analysis requested
#' @return a list containing the results of the power analysis, the population covariance matrix Sigma, the H0 implied matrix SigmaHat, and the H1 matrix SigmaH1, as well as various lavaan model strings
#' @examples
#' \dontrun{
#' ## a priori power analysis for the null hypothesis that the correlation between 
#' ## two cfa factors with a true correlation of .2 differs from zero
#' # define population model 
#' mPop = '
#'   f1 =~ .5*x1 + .6*x2 + .4*x3
#'   f2 =~ .7*x4 + .8*x5 + .3*x6
#'   x1 ~~ .75*x1
#'   x2 ~~ .64*x2
#'   x3 ~~ .84*x3
#'   x4 ~~ .51*x4
#'   x5 ~~ .36*x5
#'   x6 ~~ .91*x6
#'   f1 ~~ 1*f1
#'   f2 ~~ 1*f2
#'   f1 ~~ .2*f2
#' '
#' # define analysis model (restricting the factor correlation to zero) 
#' mAna = '
#'   f1 =~ x1 + x2 + x3
#'   f2 =~ x4 + x5 + x6
#'   f1 ~~ 0*f2
#' '
#' # do a priori power analsis
#' lavpower.ap <- semPower.powerLav(type = 'a-priori', 
#'                                  modelPop = mPop, modelH0 = mAna,
#'                                  alpha = .05, beta = .05)
#' summary(lavpower.ap$power)
#'
#' }
#' @importFrom utils installed.packages
#' @export
semPower.powerLav <- function(type, 
                              modelPop = NULL, modelH0 = NULL, modelH1 = NULL, 
                              Sigma = NULL, mu = NULL, ...){

  # check whether lavaan is available
  if(!'lavaan' %in% rownames(installed.packages())) stop('This function depends on the lavaan package, so install lavaan first.')
  # validate power type
  type <- checkPowerTypes(type)
  # validate input
  if(is.null(modelH0)) stop('Provide a lavaan model string defining the analysis (H0) model.')
  if(is.null(modelPop) & is.null(Sigma)) stop('Either provide a lavaan model string defining the population model or provide the population covariance matrix Sigma.')
  
  
  # determine population Sigma and mu
  SigmaH1 <- Sigma
  muH1 <- mu
  if(is.null(Sigma)){
    SigmaH1 <- Sigma <- lavaan::fitted(lavaan::sem(modelPop))$cov
    muH1 <- mu <- lavaan::fitted(lavaan::sem(modelPop))$mean
  }
  # override H1 sigma when comparison model differs from the saturated model 
  if(!is.null(modelH1)){
    h1.model <- lavaan::sem(modelH1, sample.cov = Sigma, sample.mean = mu, sample.nobs = 1000, likelihood = 'wishart', sample.cov.rescale = FALSE)
    if(!h1.model@optim$converged) stop('The H1 model did not converge.')
    SigmaH1 <- lavaan::fitted(h1.model)$cov
    if(!is.null(mu)) muH1 <- lavaan::fitted(h1.model)$mean
  }
  
  # get H0 sigmaHat
  h0.model <- lavaan::sem(modelH0, sample.cov = Sigma,  sample.mean = mu, sample.nobs = 1000, likelihood = 'wishart', sample.cov.rescale = FALSE)
  if(!h0.model@optim$converged) stop('The H0 model did not converge.')
  SigmaHat <- lavaan::fitted(h0.model)$cov
  muHat <- lavaan::fitted(h0.model)$mean

  # determine df
  df.h0 <- df <- h0.model@test$standard$df
  if(!is.null(modelH1)){
    df.h1 <- h1.model@test$standard$df
    df <- df.h0 - df.h1
  }
  
  # do power analysis
  power <- semPower(type = type, 
                    SigmaHat = SigmaHat, Sigma = SigmaH1, 
                    muHat = muHat, mu = muH1, 
                    df = df, 
                    ...)
  
  list(power = power, 
       SigmaHat = SigmaHat, SigmaHatH1 = SigmaH1, Sigma = Sigma,
       muHat = muHat, muHatH1 = muH1, mu = mu,
       modelPop = modelPop, modelH0 = modelH0, modelH1 = modelH1)  
  
}



#' semPower.powerCFA
#'
#' Convenience function for performing power analysis for simple CFA models involving one hypothesized zero correlation between factors.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of 'a-priori', 'post-hoc', 'compromise'
#' @param comparison comparison model, one of 'saturated' or 'restricted'. This determines the df for power analyses. 'Saturated' provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. 'Restricted' provides power to reject the model when compared to a model that just restricts the parameter defined by nullCor to zero, so the df are always 1.
#' @param nullCor vector of size 2 indicating which factor correlation in phi is hypothesized to equal zero, e.g. c(1, 2) to refer to the correlation between first and second factor. Can be omitted for two-factor models.
#' @param ... other parameters specifying the factor model (see [semPower.genSigma()]) and the type of power analyses 
#' @return a list containing the results of the power analysis, Sigma and SigmaHat, the implied loading matrix (lambda), as well as several lavaan model strings (modelPop, modelTrue, and modelAna) 
#' @examples
#' \dontrun{
#' # a priori power analysis only providing the number of indicators to define 
#' # two factors with correlation of phi and same loading for all indicators
#' cfapower.ap <- semPower.powerCFA(type = 'a-priori', 
#'                                  Phi = .2, nIndicator = c(5, 6), loadM = .5,
#'                                  alpha = .05, beta = .05)
#'                                  summary(cfapower.ap$power)
#'
#' # sanity check: fit true model to population Sigma to evaluate everything was set up as intended
#' summary(lavaan::sem(cfapower.ap$modelTrue, sample.cov = cfapower.ap$Sigma,
#'                     sample.nobs = 1000, likelihood = 'wishart', sample.cov.rescale = FALSE),
#'         stand = TRUE)
#'
#' # peek into lavaan model strings:
#' # population model
#' cfapower.ap$modelPop
#' # (incorrect) analysis model
#' cfapower.ap$modelAna
#' 
#' # or plug the population Sigma and model-implied SigmaHat 
#' # into a regular power analysis command 
#' ph <- semPower.aPriori(SigmaHat = cfapower.ap$SigmaHat, Sigma = cfapower.ap$Sigma, 
#'                        df = 1, alpha = .05, beta = .05)
#'                        summary(ph)
#'
#' # same as above, but compare to the saturated model 
#' # (rather than to the less restricted model)
#' #' cfapower.ap <- semPower.powerCFA(type = 'a-priori', comparison = 'saturated', 
#'                                  Phi = .2, nIndicator = c(5, 6), loadM = .5,
#'                                  alpha = .05, beta = .05)
#'
#' # same as above, but request a compromise power analysis
#' cfapower.cp <- semPower.powerCFA(type = 'compromise',
#'                                  Phi = .2, nIndicator = c(5, 6), loadM = .5,
#'                                  abratio = 1, N = 200)
#'
#' # same as above, but request a post-hoc power analysis
#' cfapower.ph <- semPower.powerCFA(type = 'post-hoc', 
#'                                  Phi = .2, nIndicator = c(5, 6), loadM = .5,
#'                                  alpha = .05, N = 200)
#'
#' # post-hoc power analysis providing factor correlation matrix 
#' # and reduced loading matrix 
#' Phi <- matrix(c(
#'                 c(1.0, 0.1),
#'                 c(0.1, 1.0)
#'               ), byrow = T, ncol = 2)
#'
#' # loadings: only define primary loadings 
#' # must not contain secondary loadings
#' loadings <- list(
#'                c(0.4, 0.5, 0.8),
#'                c(0.7, 0.6, 0.5, 0.4)
#'                )
#'
#' cfapower <- semPower.powerCFA(type = 'post-hoc',
#'                               Phi = phi, loadings = loadings,
#'                               alpha = .05, N = 250)
#' 
#' }
#' @seealso [semPower.genSigma()]
#' @export
semPower.powerCFA <- function(type, comparison = 'restricted', nullCor = NULL, ...){

  # validate input
  comparison <- tolower(comparison)
  if(!comparison %in% c('saturated', 'restricted')) stop('Comparison model must be one of "saturated" or "restricted"')
  
  # generate sigma 
  generated <- semPower.genSigma(...)

  ### now validate nullCor
  if(is.null(nullCor) & ncol(generated$Phi) == 2) nullCor <- c(1, 2)
  if(is.null(nullCor)) stop('nullCor must be defined')
  if(length(nullCor) != 2) stop('nullCor must be a vector of size 2')
  if(nullCor[1] == nullCor[2]) stop('nullCor may not refer to variances.')
  
  
  ### ana model 
  modelAna <- paste(c(
    generated$modelTrue,
    paste0('f', nullCor, collapse = ' ~~ 0*')),
    collapse = '\n')
  
  
  ### plug model strings into lavpower
  modelH1 <- NULL
  if(comparison == 'restricted') modelH1 <- generated$modelTrue
  lavpower <- semPower.powerLav(type = type,
                                modelPop = generated$modelPop,
                                modelH0 = modelAna,
                                modelH1 = modelH1,
                                ...)
  
  append(lavpower, generated)
  
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
      dummyFit <- suppressWarnings(lavaan::sem(lavModel, sample.cov = dummyS, sample.nobs = 1000, warn = F))
    }else{
      if(is.null(group.equal)){
        dummyFit <- suppressWarnings(lavaan::sem(lavModel, sample.cov = lapply(1:nGroups, function(x) dummyS), sample.nobs = rep(1000, nGroups), warn = F))
      }else{
        dummyFit <- suppressWarnings(lavaan::sem(lavModel, sample.cov = lapply(1:nGroups, function(x) dummyS), sample.nobs = rep(1000, nGroups), group.equal = group.equal, warn = F))
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
