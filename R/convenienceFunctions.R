#' semPower.powerCFA
#'
#' Convenience function for performing power analysis for simple CFA models involving one hypothesized zero correlation between factors.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of 'a-priori', 'post-hoc', 'compromise'
#' @param comparison comparison model, one of 'saturated' or 'restricted'. This determines the df for power analyses. 'Saturated' provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. 'Restricted' provides power to reject the model when compared to a model that just restricts the parameter defined by nullCor to zero, so the df are always 1.
#' @param phi factor correlation matrix or single number giving correlation between all factors 
#' @param nullCor vector of size 2 indicating which factor correlation in phi is hypothesized to equal zero, e.g. c(1, 2) to refer to the correlation between first and second factor
#' @param loadings a list providing factor loadings by factor. Must not contain secondary loadings.   
#' @param nIndicator vector indicating the number of indicators for each factor, e.g. c(4, 6) to define two factors with 4 and 6 indicators, respectively 
#' @param loadM vector giving mean loadings for each factor or single number to use for every loading
#' @param loadSD vector giving the standard deviation of loadings for each factor for use in conjunction with loadM. When NULL, SDs are set to zero.
#' @param loadMinMax list giving the minimum and maximum loading for each factor or vector to apply to all factors 
#' @return a list containing the results of the power analysis, Sigma and SigmaHat, the implied loading matrix (lambda), as well as several lavaan model strings (modelPop, modelTrue, and modelAna) 
#' @param ... other parameters related to the specific type of power analysis requested
#' @examples
#' \dontrun{
#' # a priori power analysis only providing the number of indicators to define 
#' # two factors with correlation of phi and same loading for all indicators
#' cfapower.ap <- semPower.powerCFA(type = 'a-priori', 
#'                                  phi = .2, nIndicator = c(5, 6), loadM = .5,
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
#'                                  phi = .2, nIndicator = c(5, 6), loadM = .5,
#'                                  alpha = .05, beta = .05)
#'
#' # same as above, but request a compromise power analysis
#' cfapower.cp <- semPower.powerCFA(type = 'compromise',
#'                                  phi = .2, nIndicator = c(5, 6), loadM = .5,
#'                                  abratio = 1, N = 200)
#'
#' # same as above, but request a post-hoc power analysis
#' cfapower.ph <- semPower.powerCFA(type = 'post-hoc', 
#'                                  phi = .2, nIndicator = c(5, 6), loadM = .5,
#'                                  alpha = .05, N = 200)
#'
#' # post-hoc power analysis providing factor correlation matrix 
#' # and reduced loading matrix 
#' phi <- matrix(c(
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
#'                               phi = phi, nullCor = c(1, 2), loadings = loadings,
#'                               alpha = .05, N = 250)
#'
#' # post-hoc power analysis providing factor correlation matrix, 
#' # number of indicators by factor, and min-max loading for all factors 
#' phi <- matrix(c(
#'                 c(1.0, 0.2, 0.5),
#'                 c(0.2, 1.0, 0.3),
#'                 c(0.5, 0.3, 1.0)
#'                ), byrow = TRUE, ncol = 3)
#'
#' cfapower <- semPower.powerCFA(type = 'post-hoc',
#'                               phi = phi, nullCor = c(1, 2), 
#'                               nIndicator = c(6, 5, 4), loadMinMax = c(.3, .8),
#'                               alpha = .05, N = 250)
#'
#' # same as above, but providing mean and sd loading for all factors
#' cfapower <- semPower.powerCFA(type = 'post-hoc',
#'                               phi = phi, nullCor = c(1, 2), 
#'                               nIndicator = c(6, 5, 4), loadM = .5, loadSD = .1, 
#'                               alpha = .05, N = 250)
#'
#' # same as above, but hypothesizing zero correlation between factors 2 and 3
#' cfapower <- semPower.powerCFA(type = 'post-hoc',
#'                               phi = phi, nullCor = c(2, 3), 
#'                               nIndicator = c(6, 5, 4), loadM = .5, loadSD = .1,
#'                               alpha = .05, N = 250)
#'
#' # same as above, but providing mean and sd of loadings for each factor
#' cfapower <- semPower.powerCFA(type = 'post-hoc',
#'                               phi = phi, nullCor = c(1, 2), 
#'                               nIndicator = c(3, 6, 5), 
#'                               loadM = c(.5, .6, .7), loadSD = c(.1, .05, 0),
#'                               alpha = .05, N = 250)
#'
#' # same as above, but using min-max loadings for each factor
#' loadMinMax <- list(
#'                    c(.4, .6),
#'                    c(.5, .8),
#'                    c(.3, .7)
#'                    )
#'
#' cfapower <- semPower.powerCFA(type = 'post-hoc',
#'                               phi = phi, nullCor = c(1, 2), nIndicator = c(3, 6, 5), 
#'                               loadMinMax = loadMinMax,
#'                               alpha = .05, N = 250)
#' 
#' }
#' @importFrom stats rnorm runif 
#' @importFrom utils installed.packages
#' @export
semPower.powerCFA <- function(type, comparison = 'restricted',
                              phi, nullCor = NULL, loadings = NULL, 
                              nIndicator = NULL, 
                              loadM = NULL, loadSD = NULL, 
                              loadMinMax = NULL, 
                              ...){
  
  # check whether lavaan is availabe
  if(!'lavaan' %in% rownames(installed.packages())) stop('This function depends on the lavaan package, so install lavaan first.')
  
  # validate power type
  type <- checkPowerTypes(type)
  
  # validate input
  comparison <- tolower(comparison)
  if(!comparison %in% c('saturated', 'restricted')) stop('Comparison model must be one of "saturated" or "restricted"')
  
  if(is.null(nIndicator) & is.null(loadings)) stop('Either provide loadings or number of indicators')
  if(!is.null(nIndicator) & !is.null(loadings)) stop('Either provide loadings or number of indicators, but not both.')
  if(is.null(nIndicator) & !is.list(loadings)) stop('loadings must be a list')
  
  nfac <- ifelse(is.null(loadings), length(nIndicator), length(loadings))
  if(nfac < 2) stop('At least two factors are required')
  
  if(is.null(nullCor) & length(phi) == 1) nullCor <- c(1, 2)
  if(is.null(nullCor)) stop('nullCor must be defined')
  if(length(nullCor) != 2) stop('nullCor must be a vector of size 2')
  if(nullCor[1] == nullCor[2]) stop('nullCor may not refer to variances.')
  
  if(length(phi) == 1){
    phi <- matrix(phi, ncol = nfac, nrow = nfac)
    diag(phi) <- 1
  } 
  checkPositiveDefinite(phi)
  if(ncol(phi) != nfac) stop('phi must have the same number of rows/columns as the number of factors') 
  invisible(apply(phi, c(1, 2), function(x) checkBounded(x, 'All elements in phi', bound = c(-1, 1), inclusive = TRUE)))
  if(nullCor[1] > nrow(phi) | nullCor[2] > ncol(phi)) stop('nullCor does not refer to a valid correlation in phi')
  
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
    invisible(lapply(loadings, function(x) lapply(x, function(x) checkBounded(x, 'All loadings', bound = c(-1, 1)))))
    nIndicator <- unlist(lapply(loadings, length))  # crossloadings are disallowed
  }

  ### pop model
  # define factors
  tok <- list()
  lambda <- matrix(0, ncol = nfac, nrow = sum(nIndicator)) # store loading matrix
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
    }else{
      cload <- loadings[[f]]  # loadings only contains primary loadings
    }
    lambda[sidx:eidx, f] <- cload
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
      tok <- append(tok, paste0('f', f, ' ~~ ', phi[f, ff], '*f', ff))
    }
  }
  
  modelPop <- paste(unlist(tok), collapse = '\n')
  
  
  ### ana model 
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
  
  modelAna <- paste(c(
    modelTrue,
    paste0('f', nullCor, collapse = ' ~~ 0*')),
    collapse = '\n')
  
  
  # get sigmas
  Sigma <- lavaan::fitted(lavaan::sem(modelPop))$cov
  hyp.model <- lavaan::sem(modelAna, sample.cov = Sigma, sample.nobs = 1000, likelihood = 'wishart', sample.cov.rescale = FALSE)
  SigmaHat <- lavaan::fitted(hyp.model)$cov
  
  # do power analysis
  df <- ifelse(comparison == 'saturated', hyp.model@test$standard$df, 1)
  power <- semPower(type = type, Sigma = Sigma, SigmaHat = SigmaHat, df = df, ...)
  
  list(power = power, Sigma = Sigma, SigmaHat = SigmaHat, modelPop = modelPop, modelTrue = modelTrue, modelAna = modelAna, lambda = lambda, phi = phi)  
  
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
