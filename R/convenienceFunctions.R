#' semPower.powerLav
#'
#' Perform power analysis on population and model-implied Sigmas as defined through lavaan model strings
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of 'a-priori', 'post-hoc', 'compromise'.
#' @param modelPop lavaan model string defining the true model.
#' @param modelH0 lavaan model string defining the (incorrect) analysis model.
#' @param modelH1 lavaan model string defining the comparison model. If omitted, the saturated model is the comparison model.
#' @param fitH1model whether to fit the H1 model. If FALSE, the H1 model is assumed to show the same fit as the saturated model, and only the delta df are computed.
#' @param Sigma population covariance matrix (if modelPop is not set).
#' @param mu population means (if modelPop is not set).
#' @param simulatedPower whether to perform a simulated (TRUE) (rather than analytical, FALSE) power analysis.
#' @param lavOptions a list of additional options passed to lavaan, e.g., list(estimator = 'mlm') to request robust ML estimation. Probably mostly useful in conjunction with simulatedPower. 
#' @param ... other parameters related to the specific type of power analysis requested
#' @return a list containing the results of the power analysis, the population covariance matrix Sigma and mean vector mu, the H0 implied matrix SigmaHat and mean vector muHat, as well as various lavaan model strings (modelPop, modelH0, and modelH1)
#' @examples
#' \dontrun{
#' ## a priori power analysis for the null hypothesis that the correlation between 
#' ## two cfa factors with a true correlation of .2 differs from zero
#' # define population model 
#' mPop <- '
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
#' mAna <- '
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
                              modelPop = NULL, modelH0 = NULL, modelH1 = NULL, fitH1model = TRUE, 
                              Sigma = NULL, mu = NULL, 
                              simulatedPower = FALSE, 
                              lavOptions = NULL, ...){
  
  # check whether lavaan is available
  if(!'lavaan' %in% rownames(installed.packages())) stop('This function depends on the lavaan package, so install lavaan first.')
  # validate input
  type <- checkPowerTypes(type)
  if(is.null(modelH0)) stop('Provide a lavaan model string defining the analysis (H0) model.')
  if(is.null(modelPop) && is.null(Sigma)) stop('Either provide a lavaan model string defining the population model or provide the population covariance matrix Sigma.')
  if(!is.null(modelPop) && !is.null(Sigma)) stop('Either provide a lavaan model string defining the population model or provide the population covariance matrix Sigma, but not both.')
  if(simulatedPower && type == 'compromise') stop('Simulated power is not available for compromise power analysis, because this would require a vast (infeasible) number of simulation runs to yield reliable results.')

  # determine population Sigma / mu
  if(is.null(Sigma)){
    Sigma <- lavaan::fitted(lavaan::sem(modelPop))[['cov']]
    mu <- lavaan::fitted(lavaan::sem(modelPop))[['mean']]
  }
  
  # analytical power
  if(!simulatedPower){
    
    # we need to call lavaan() directly with defaults as defined in sem()
    lavOptions <- getLavOptions(lavOptions)
    if(!is.null(lavOptions[['estimator']]) && toupper(lavOptions[['estimator']]) != "ML") stop('Analytical power is only available with ML estimation. Note that power based on ML derivatives (mlm etc) is asymptotically identical.')

    # get H0 sigmaHat / muHat
    modelH0Fit <- do.call(lavaan::lavaan,
                          append(list(model = modelH0,
                                      sample.cov = Sigma,
                                      sample.mean = mu),
                                 lavOptions))
    if(!modelH0Fit@optim[['converged']]) stop('The H0 model did not converge.')
    SigmaHat <- lavaan::fitted(modelH0Fit)[['cov']]
    muHat <- lavaan::fitted(modelH0Fit)[['mean']]
    df <- dfH0 <- modelH0Fit@test[['standard']][['df']]
    
    # get H1 comparison model and deltaF
    if(!is.null(modelH1) && fitH1model){
      modelH1Fit <- do.call(lavaan::lavaan, 
                            append(list(model = modelH1, 
                                        sample.cov = Sigma,
                                        sample.mean = mu),
                                   lavOptions))
      if(!modelH1Fit@optim[['converged']]) stop('The H1 model did not converge.')
      dfH1 <- modelH1Fit@test[['standard']][['df']]
      if(dfH1 >= dfH0) stop('The df of the H1 model are not larger than the df of the H0 model, as they should be.')
      # get delta F
      fminH1 <- getF.Sigma(lavaan::fitted(modelH1Fit)[['cov']], Sigma, lavaan::fitted(modelH1Fit)[['mean']], mu)
      fminH0 <- getF.Sigma(lavaan::fitted(modelH0Fit)[['cov']], Sigma, lavaan::fitted(modelH0Fit)[['mean']], mu)
      deltaF <- fminH0 - fminH1
      df <- (dfH0 - dfH1)
    }else if (!is.null(modelH1) && !fitH1model){
      df <- df - semPower.getDf(modelH1)
    }

    # we use sigma for the comparison with the saturated model (so we also get additional fitindices) 
    # but delta f for the comparison with an explicit h1 model.
    if(is.null(modelH1) || !fitH1model){
      power <- semPower(type = type, 
                        SigmaHat = SigmaHat, Sigma = Sigma, 
                        muHat = muHat, mu = mu, 
                        df = df, 
                        ...)    
    }else{
      power <- semPower(type = type, 
                        effect = deltaF, effect.measure = "F0", 
                        df = df, 
                        ...)    
    }    
        
  # simulated power
  }else{
    power <- semPower(type = type, 
                      Sigma = Sigma, mu = mu, 
                      modelH0 = modelH0, modelH1 = modelH1, fitH1model = fitH1model,
                      simulatedPower = simulatedPower, 
                      lavOptions = lavOptions,   # simulate() takes care of proper lavOptions
                      ...)
    SigmaHat <- muHat <- NULL
  }

  list(power = power, 
       SigmaHat = SigmaHat, Sigma = Sigma,
       muHat = muHat, mu = mu,
       modelPop = modelPop, modelH0 = modelH0, modelH1 = modelH1)
}



#' semPower.powerCFA
#'
#' Convenience function for performing power analysis for simple CFA models involving one hypothesized zero correlation between factors.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of 'a-priori', 'post-hoc', 'compromise'
#' @param comparison comparison model, one of 'saturated' or 'restricted'. This determines the df for power analyses. 'Saturated' provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. 'Restricted' provides power to reject the model when compared to a model that just restricts the parameter defined by nullCor to zero, so the df are always 1.
#' @param Phi either a single number defining the correlation between exactly two factors or the factor correlation matrix.
#' @param nullEffect defines the hypothesis of interest. Valid are 'cor = 0' (the default) and 'corX = corZ' to test for the equality of correlations. Define the correlations to be set to equality in nullWhich 
#' @param nullWhich vector of size 2 indicating which factor correlation in phi is hypothesized to equal zero when nullEffect = 'cor = 0' or list of vectors defining which correlations to restrict to equality when nullEffect = 'corX = corZ'. Can also contain more than two correlations, e.g., list(c(1,2), c(1,3), c(2,3)) to set phi[1,2] = phi[1,3] = phi[2,3]
#' @param ... other parameters specifying the factor model (see [semPower.genSigma()]) and the type of power analysis 
#' @return a list containing the results of the power analysis, the population covariance matrix Sigma and mean vector mu, the H0 implied matrix SigmaHat and mean vector muHat, as well as various lavaan model strings (modelH0, and modelH1)
#' @examples
#' \dontrun{
#' # a priori power analysis only providing the number of indicators to define 
#' # two factors with correlation of phi and same loading for all indicators
#' cfapower.ap <- semPower.powerCFA(type = 'a-priori',
#'                                  nullWhich = c(1, 2), 
#'                                  Phi = .2, nIndicator = c(5, 6), loadM = .5,
#'                                  alpha = .05, beta = .05)
#'                                  summary(cfapower.ap$power)
#'
#' # same as above, but compare to the saturated model 
#' # (rather than to the less restricted model)
#' cfapower.ap <- semPower.powerCFA(type = 'a-priori', comparison = 'saturated', 
#'                                  nullWhich = c(1, 2), 
#'                                  Phi = .2, nIndicator = c(5, 6), loadM = .5,
#'                                  alpha = .05, beta = .05)
#'
#' # same as above, but request a compromise power analysis
#' cfapower.cp <- semPower.powerCFA(type = 'compromise',
#'                                  nullWhich = c(1, 2), 
#'                                  Phi = .2, nIndicator = c(5, 6), loadM = .5,
#'                                  abratio = 1, N = 200)
#'
#' # same as above, but request a post-hoc power analysis
#' cfapower.ph <- semPower.powerCFA(type = 'post-hoc', 
#'                                  nullWhich = c(1, 2), 
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
#'                               nullWhich = c(1, 2), 
#'                               Phi = Phi, loadings = loadings,
#'                               alpha = .05, N = 250)
#' 
#' }
#' @seealso [semPower.genSigma()]
#' @export
semPower.powerCFA <- function(type, comparison = 'restricted', 
                              Phi = NULL,
                              nullEffect = 'cor = 0',
                              nullWhich = NULL, ...){
  
  # validate input
  comparison <- checkComparisonModel(comparison)
  if(is.null(Phi)) stop('Phi must be defined')
  if(is.null(nullEffect)) stop('nullEffect must be defined.')
  if(length(nullEffect) > 1) stop('nullEffect must contain a single hypothesis')
  nullEffect <- unlist(lapply(nullEffect, function(x) tolower(trimws(x))))
  nullEffect <- gsub(" ", "", nullEffect, fixed = TRUE)
  if(any(unlist(lapply(nullEffect, function(x) !x %in% c('cor=0', 'corx=corz'))))) stop('nullEffect must be either cor=0 or corx=corz')
  
  # generate sigma 
  generated <- semPower.genSigma(Phi = Phi, ...)
  
  ### now do validation of nullWhich, since we now know Phi
  if(is.null(nullWhich) && ncol(generated[['Phi']]) == 2) nullWhich <- c(1, 2)
  if(is.null(nullWhich)) stop('nullWhich must be defined.')
  if(!is.list(nullWhich)) nullWhich <- list(nullWhich)
  if(any(unlist(lapply(nullWhich, function(x) length(x) != 2)))) stop('nullWhich may only contain vectors of size two.')
  if(any(unlist(lapply(nullWhich, function(x) x[1] == x[2])))) stop('elements in nullWhich may not refer to variances.')
  if(any(unlist(lapply(nullWhich, function(x) (x[1] < 1 | x[2] < 1 | x[1] > ncol(generated[['Phi']]) | x[2] > ncol(generated[['Phi']])))))) stop('At least one element in nullWhich is an out of bounds index concerning Phi.')
  if(length(nullWhich) > 1){
    for(i in 1:(length(nullWhich) - 1)){
      for(j in (i + 1):length(nullWhich)){
        if(nullWhich[[i]][1] == nullWhich[[j]][1] && nullWhich[[i]][2] == nullWhich[[j]][2]) stop('elements in nullWhich may not refer to the same correlation')
      }
    }
  }
  
  ### H0 model 
  modelH0 <- generated[['modelTrueCFA']]
  if(nullEffect == 'cor=0'){
    modelH0 <- paste(c(modelH0,
      paste0('f', nullWhich[[1]], collapse = ' ~~ 0*')),
      collapse = '\n')
  }else{
    labs <- list()
    tok <- ''
    for(i in 1:length(nullWhich)){
      cl <- paste0('pf',paste0(formatC(nullWhich[[i]], width = 2, flag = 0), collapse = ''))
      tok <- paste(tok, paste0('f', nullWhich[[i]][1], ' ~~ ', cl, '*f', nullWhich[[i]][2]), sep = '\n')
      labs <- append(labs, cl)
    }
    labs <- unlist(labs)
    for(i in 1:(length(labs) - 1)){
      for(j in (i + 1):length(labs)){
        tok <- paste(tok, paste(labs[i], ' == ', labs[j]), sep = '\n')
      }
    }
    modelH0 <- paste(c(modelH0, tok), collapse = '\n')
  }
  
  modelH1 <- NULL
  if(comparison == 'restricted'){
    # h1 model always fits perfectly, only needed for delta df
    modelH1 <- generated[['modelTrueCFA']]
    fitH1model <- FALSE 
  } 
  
  semPower.powerLav(type = type,
                    Sigma = generated[['Sigma']],
                    modelH0 = modelH0,
                    modelH1 = modelH1,
                    fitH1model = fitH1model,
                    ...)
  
}

#' semPower.powerRegression
#'
#' Convenience function for performing power analysis on slope(s) in a latent regression.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of 'a-priori', 'post-hoc', 'compromise'
#' @param comparison comparison model, one of 'saturated' or 'restricted'. This determines the df for power analyses. 'Saturated' provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. 'Restricted' provides power to reject the model when compared to a model that just restricts the indirect effect to zero, so the df are always 1.
#' @param slopes vector of standardized slopes (or a single number for a single slope) of the k predictors for Y. 
#' @param nullEffect defines the hypothesis of interest. Valid are 'slope = 0' (the default) and 'slopeX = slopeZ' to test for the equality of slopes. Define the slopes to set to equality in nullWhich 
#' @param nullWhich single number indicating which slope is hypothesized to equal zero when nullEffect = 'slope = 0' or vector defines which slopes to restrict to equality when nullEffect = 'slopeX = slopeZ'. Can also contain more than two slopes.
#' @param corXX correlation(s) between the k predictors (X). Either NULL, a single number (for k = 2 predictors), or a matrix. If NULL, the predictors are uncorrelated. 
#' @param ... other parameters specifying the factor model (see [semPower.genSigma()]), where the first factor corresponds to Y (so you need nfactors = slopes + 1), and the type of power analysis 
#' @return a list containing the results of the power analysis, the population covariance matrix Sigma and mean vector mu, the H0 implied matrix SigmaHat and mean vector muHat, as well as various lavaan model strings (modelH0, and modelH1)
#' @examples
#' \dontrun{
#' # latent regression of the form Y = .2*X1 + .3*X2, where X1 and X2 correlate by .4
#' # request power for the hypothesis that the slope of X1 is zero. 
#' # providing the number of indicators by factor (Y, X1, X2) each loading by the same magnitude on its designed factor.
#' regPower <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), corXX = .4, nullWhich = 1, 
#'                                      nIndicator = c(3, 5, 4),
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#' summary(regPower$power)
#' 
#' # same as above, but ask for power to detect the  slope of X2
#' regPower <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), corXX = .4, nullWhich = 2, 
#'                                      nIndicator = c(3, 5, 4),
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#' 
#' # latent regression with three predictors, providing the predictor intercorrelation matrix
#' corXX <- matrix(c(
#'   c(1.00, 0.20, 0.30),
#'   c(0.20, 1.00, 0.10),
#'   c(0.30, 0.10, 1.00)
#' ), ncol = 3,byrow = TRUE)
#' regPower <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3, .4), corXX = corXX, 
#'                                      nullWhich = 1,
#'                                      nIndicator = c(4, 3, 5, 4),
#'                                      loadM = c(.5, .5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'
#' # same as above, but testing the equality of the first and second slope
#' regPower <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3, .4), corXX = corXX, 
#'                                      nullEffect = 'slopeX = slopeZ', 
#'                                      nullWhich = c(1, 2),
#'                                      nIndicator = c(4, 3, 5, 4),
#'                                      loadM = c(.5, .5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'                                      
#' # same as above, but testing the equality of all three slopes
#' regPower <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3, .4), corXX = corXX, 
#'                                      nullEffect = 'slopeX = slopeZ', 
#'                                      nullWhich = c(1, 2, 3),
#'                                      nIndicator = c(4, 3, 5, 4),
#'                                      loadM = c(.5, .5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'                                      
#' }
#' @seealso [semPower.genSigma()]
#' @export
semPower.powerRegression <- function(type, comparison = 'restricted',
                                     slopes = NULL, 
                                     corXX = NULL, 
                                     nullEffect = 'slope = 0',
                                     nullWhich = NULL,
                                     ...){
  
  comparison <- checkComparisonModel(comparison)
  
  # we override Phi and Sigma later, so let's make sure it is not set in ellipsis argument
  if('Phi' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Phi, because the factor correlations depend on corXX and the slopes.')
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma, because Sigma is determined as function of corXX and the slopes.')
  
  # validate input
  if(is.null(slopes)) stop('slopes cannot be NULL.')
  if(!is.vector(slopes) && !is.matrix(slopes)) stop('slopes must be a single number of a vector')
  invisible(lapply(slopes, function(x) checkBounded(x, 'All slopes ', bound = c(-1, 1), inclusive = TRUE)))
  if(sum(slopes^2) > 1) stop('slopes imply a negative residual variance for Y, make sure that the sum of the squared slopes is < 1')
  if(!is.matrix(slopes)) slopes <- matrix(slopes, nrow = length(slopes))
  
  if(is.null(corXX)) corXX <- diag(nrow(slopes)) 
  if(is.vector(corXX) && length(corXX) > 1) stop('corXX must be a single number or a matrix') 
  if(!is.matrix(corXX)){
    corXX <- matrix(corXX, nrow = 2, ncol = 2) 
    diag(corXX) <- 1
  } 
  checkPositiveDefinite(corXX)
  if(ncol(corXX) != nrow(slopes)) stop('Dimension of corXX does not match number of predictors.')
  
  if(is.null(nullEffect)) stop('nullEffect must be defined.')
  if(length(nullEffect) > 1) stop('nullEffect must contain a single hypothesis')
  nullEffect <- unlist(lapply(nullEffect, function(x) tolower(trimws(x))))
  nullEffect <- gsub(" ", "", nullEffect, fixed = TRUE)
  if(any(unlist(lapply(nullEffect, function(x) !x %in% c('slope=0', 'slopex=slopez'))))) stop('nullEffect must be either slope=0 or slopex=slopez')
  
  if(is.null(nullWhich)) stop('nullWhich must be defined.')
  if(any(nullWhich < 1) || any(nullWhich > nrow(slopes))) stop('nullWhich is invalid.')
  if(nullEffect == 'slopex=slopez'){
    if(length(nullWhich) < 2 || length(nullWhich) > nrow(slopes)) stop('nullWhich must contain at least two slopes when nullEffect is slopex=slopez, but not more slopes than available')
  }else{
    if(length(nullWhich) > 1) stop('nullWhich must be a single number when nullEffect is slope=0')
  }
  nullWhich <- nullWhich + 1 # because first factor is criterion
  
  # calc implied sigma. we do this here, because this is a special case and simpler than defining B and calling getPhi.B  
  corXY <- (corXX %*% slopes)
  Phi <- t(c(1, corXY))
  Phi <- rbind(Phi, cbind(corXY, corXX))
  generated <- semPower.genSigma(Phi = Phi, useReferenceIndicator = TRUE, ...)

  ### create ana model string
  # add regressions 
  model <- paste(generated[['modelTrueCFA']], 
                 paste0('f1 ~ ', paste0(paste0('pf',(1 + 1:ncol(corXX))), '*f',(1 + 1:ncol(corXX)), collapse = ' + ')), 
                 sep = '\n')
  modelTrue <- model
  
  if(nullEffect == 'slopex=slopez'){
    tok <- ''
    for(i in 1:(length(nullWhich) - 1)){
      for(j in (i + 1):length(nullWhich)){
        tok <- paste(tok, paste0('pf', nullWhich[i], ' == ', 'pf', nullWhich[j]), sep = '\n')
      }
    }
    modelH0 <- paste(model, '\n', tok)  
  }else{
    modelH0 <- paste(model, '\n', paste0('pf', nullWhich,' == 0'))  
  }
  
  modelH1 <- NULL
  if(comparison == 'restricted'){
    # h1 model always fits perfectly, only needed for delta df
    modelH1 <- modelTrue
    fitH1model <- FALSE 
  } 
  
  semPower.powerLav(type, 
                    Sigma = generated[['Sigma']], 
                    modelH0 = modelH0, 
                    modelH1 = modelH1, 
                    fitH1model = fitH1model,
                    ...)
}


#' semPower.powerMediation
#'
#' Convenience function for performing power analysis concerning indirect effect(s) in a mediation model.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of 'a-priori', 'post-hoc', 'compromise'
#' @param comparison comparison model, one of 'saturated' or 'restricted'. This determines the df for power analyses. 'Saturated' provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. 'Restricted' provides power to reject the model when compared to a model that just restricts the indirect effect to zero, so the df are always 1.
#' @param bYX the standardized slope (direct effect) for X -> Y 
#' @param bMX the standardized slope for X -> M
#' @param bYM the standardized slope for M -> Y
#' @param Beta matrix of regression weights connecting the latent factors, akin to all-Y notation.
#' @param indirect a list of indices indicating the elements of B that define the indirect effect of interest, e.g. list(c(2,1),c(3,2)).
#' @param ... other parameters specifying the factor model (see [semPower.genSigma()]) and the type of power analysis 
#' @return a list containing the results of the power analysis, the population covariance matrix Sigma and mean vector mu, the H0 implied matrix SigmaHat and mean vector muHat, as well as various lavaan model strings (modelH0, and modelH1)
#' @examples
#' \dontrun{
#' # simple case of X -> M -> Y mediation
#' # X -- .30 --> M -- .40 --> Y 
#' # X --------- .25 --------> Y 
#' # providing the number of indicators by factor (X, M, Y) each loading by the same magnitude on its designed factor.
#' medPower <- semPower.powerMediation(type = 'a-priori', 
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(3, 5, 4), 
#'                                     loadM = c(.5, .6, .7),                                     
#'                                     alpha = .05, beta = .05)
#' summary(medPower$power)
#'                                     
#' # same mediation model as above, but assuming single loading of 1 for all factors (=> mediation model with manifest variables)
#' medPower <- semPower.powerMediation(type = 'a-priori', 
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(1, 1, 1), loadM = 1,
#'                                     alpha = .05, beta = .05)
#' 
#' # same latent mediation model as above, but specifying loadings through Lambda 
#' Lambda <- matrix(c(
#'                 c(0.5, 0.0, 0.0),    # X, M, Y
#'                 c(0.4, 0.0, 0.0),
#'                 c(0.3, 0.0, 0.0),
#'                 c(0.0, 0.7, 0.0),
#'                 c(0.0, 0.8, 0.0),
#'                 c(0.0, 0.5, 0.0),
#'                 c(0.0, 0.0, 0.5),
#'                 c(0.0, 0.0, 0.4),
#'                 c(0.0, 0.0, 0.6)
#'                ), byrow = TRUE, ncol = 3)
#' medPower <- semPower.powerMediation(type = 'a-priori', 
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     Lambda = Lambda,                                     
#'                                     alpha = .05, beta = .05)
#'
#' # same mediation model as above, but specifying Beta 
#' B <- matrix(c(
#'               c(.00, .00, .00),
#'               c(.30, .00, .00),
#'               c(.25, .40, .00)
#'               ), byrow = TRUE, ncol = 3)
#' medPower <- semPower.powerMediation(type = 'a-priori', 
#'                                     Beta = B, indirect = list(c(2,1), c(3,2)),
#'                                     Lambda = Lambda,
#'                                     alpha = .05, beta = .05)
#' 
#' # Beta for a more complex mediation hypothesis
#' # of the form X -> M1 -> M2 -> Y 
#' # and using a reduced loading matrix
#' B <- matrix(c(
#'               c(.00, .00, .00, .00),
#'               c(.20, .00, .00, .00),
#'               c(.00, .30, .00, .00),
#'               c(.00, .00, .40, .00)
#'               ), byrow = TRUE, ncol = 4)
#' # only define primary loadings by factor
#' loadings <- list(
#'                c(0.4, 0.5, 0.8),
#'                c(0.7, 0.6, 0.5, 0.8),
#'                c(0.5, 0.6, 0.3, 0.4, 0.6),
#'                c(0.6, 0.7, 0.8)
#'                )
#'
#' medPower <- semPower.powerMediation(type = 'a-priori', 
#'                                     Beta = B, indirect = list(c(2,1), c(3,2), c(4,3)),
#'                                     loadings = loadings,
#'                                     alpha = .05, beta = .05)
#'  
#' }
#' @seealso [semPower.genSigma()]
#' @export
semPower.powerMediation <- function(type, comparison = 'restricted',
                                    bYX = NULL, bMX = NULL, bYM = NULL,
                                    Beta = NULL, indirect = NULL, ...){
  
  comparison <- checkComparisonModel(comparison)
  
  # we override Phi and Sigma later, so let's make sure it is not set in ellipsis argument
  if('Phi' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Phi, because the factor correlations depend on Beta (or the slopes).')
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma, because Sigma is determined as function of Beta (or the slopes).')
  
  # validate input
  if(!is.null(Beta) && (!is.null(bYX) || !is.null(bMX) || !is.null(bYM))) stop('Either provide bYX, bMX, and bYM or provide Beta, but not both.')
  if(is.null(Beta)){
    if(is.null(bYX) || is.null(bMX) || is.null(bYM)) stop('Provide bYX, bYM, and bYM or provide Beta')
    if(length(bYX) != 1) stop('bYX must be a single slope (X -> Y)')
    if(length(bMX) != 1) stop('bMX must be a single slope (X -> M)')
    if(length(bYM) != 1) stop('bYM must be a single slope (M -> Y)')
    invisible(lapply(c(bYX, bMX, bYM), function(x) checkBounded(x, 'All slopes ', bound = c(-1, 1), inclusive = TRUE)))
    if((bYX^2 + bYM^2) > 1) stop('bYX and bYM imply a negative residual variance for Y, make sure that the sum of the squared slopes on Y is < 1')
    if(bMX == 0 || bYM == 0) stop('One of bMX and bYM is zero, implying the indirect effect is zero. The indirect effect must differ from zero.')
    indirect <- list(c(2, 1), c(3, 2))
  }
  
  if(!is.null(Beta)){
    if(is.null(indirect)) stop('indirect must not be NULL when Beta is defined.')
    if(isSymmetric(Beta)) stop('Beta must be symmetric.')
    invisible(apply(Beta, c(1, 2), function(x) checkBounded(x, 'All elements in Beta', bound = c(-1, 1), inclusive = TRUE)))
    if(any(diag(Beta) != 0)) stop('All diagonal elements of Beta must be zero.')
    # negative implied residual variances are checked in getPhi.B
    if(any(lapply(indirect, function(x) length(x)) != 2)) stop('Indirect must be a list containing vectors of size two each')
    if(any(unlist(lapply(indirect, function(x) any(x > ncol(Beta)))))) stop('At least one element in indirect is an out of bounds index concerning B')
    if(any(unlist(lapply(indirect, function(x) Beta[x[1], x[2]])) == 0)) stop('Beta and indirect imply an indirect effect of zero. The indirect effect must differ from zero.')
  }
  
  B <- Beta
  if(is.null(B)){
    B <- matrix(c(
      c(0, 0, 0),      # X
      c(bMX, 0, 0),    # M
      c(bYX, bYM, 0)   # Y
    ), byrow = TRUE, ncol = 3)
  }
  
  ### get Sigma
  # we want the completely standardized slopes, so transform to standard cfa model by converting B to implied phi
  Phi <- getPhi.B(B) 
  generated <- semPower.genSigma(Phi = Phi, useReferenceIndicator = TRUE, ...)

  ### create model strings
  model <- generated[['modelTrueCFA']]
  # add mediation structure
  for(f in 1:ncol(B)){
    fidx <- which(B[f, ] != 0)
    if(length(fidx) != 0){
      tok <- paste0('f', f, ' ~ ', paste(paste0('pf', paste0(formatC(f, width = 2, flag = 0), formatC(fidx, width = 2, flag = 0)), '*'), paste0('f', fidx), sep = '', collapse = ' + '))
      model <- paste(model, tok, sep='\n')
    }
  }
  # add indirect effects
  ind <- unlist(lapply(indirect, function(x) paste0('pf', paste0(formatC(x, width = 2, flag = 0), collapse = ''))))
  modelTrue <- paste(model, '\n', 'ind := ', paste(ind, collapse = '*'))
  
  # lav doesn't like constraining the indirect effect to zero. 
  # we instead constrain the smallest of the contained direct effects, so this
  # actually gives power for a single slope. however, this seems to closely reflect
  # power for the indirect effect (and indeed works much better than using ind=0 as 
  # comparison model, which grossly overestimates the true effect)
  # modelH0 <- paste(modelTrue, '\n', 'ind == 0')  
  cs <- indirect[[which.min(unlist(lapply(indirect, function(x) B[x[1], x[2]])))]]
  mb <- paste0('pf', paste(formatC(cs, width = 2, flag = 0), collapse = ''))
  modelH0 <- paste(modelTrue, '\n', paste0(mb,' == 0'))  
  
  modelH1 <- NULL
  if(comparison == 'restricted'){
    # h1 model always fits perfectly, only needed for delta df
    modelH1 <- modelTrue
    fitH1model <- FALSE 
  } 
  
  semPower.powerLav(type, 
                    Sigma = generated[['Sigma']], 
                    modelH0 = modelH0, 
                    modelH1 = modelH1, 
                    fitH1model = fitH1model,
                    ...)
}


#' semPower.powerCLPM
#'
#' Convenience function for performing power analysis on effects in a cross-lagged panel model (CLPM).
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of 'a-priori', 'post-hoc', 'compromise'
#' @param comparison comparison model, one of 'saturated' or 'restricted'. This determines the df for power analyses. 'Saturated' provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. 'Restricted' provides power to reject the model when compared to a model that just restricts the effect of interest to zero, so the df are always 1.
#' @param nWaves number of waves, must be >= 2.
#' @param autoregEffects vector of the autoregressive effects of X and Y (constant across waves), or a list of vectors of autoregressive effects for X and Y from wave to wave, e.g. list(c(.7, .6), c(.5, .5)) for a autoregressive effect of .7 for x1->x2 and .6 for x2->x3
#' @param crossedEffects vector of crossed effects of x on y (X -> Y) and vice versa (both constant across waves), or a list of vectors of crossed effects giving the crossed effect of x on y (and vice versa) for each wave, e.g. list(c(.2, .3), c(.1, .1)) for x1->y2 = .2 and x2->y3 = .3.
#' @param rXY vector of (residual-)correlations between X and Y for each wave. If NULL, all (residual-)correlations are zero. 
#' @param waveEqual parameters that are assumed to be equal across waves in both the H0 and the H1 model. Valid are 'autoregX' and 'autoregY' for autoregressive effects, 'crossedX' and 'crossedY' for crossed effects, 'corXY' for residual correlations, or NULL for none (so that all parameters are freely estimated, subject to the constraints defined in nullEffect). 
#' @param nullEffect defines the hypothesis of interest. Valid are the same arguments as in waveEqual and 'autoregX = 0', 'autoregY = 0', 'crossedX = 0', 'crossedY = 0' to constrain the X or Y autoregressive effects or the crossed effects to zero, 'autoregX = autoregY' and 'crossedX = crossedY' to constrain them to be equal for X and Y.
#' @param nullWhich used in conjunction with nullEffect to identify which parameter to constrain when there are > 2 waves and parameters are not constant across waves. For example, nullEffect = 'autoregX = 0' with nullWhich = 2 would constrain the second autoregressive effect for X to zero.    
#' @param metricInvariance whether metric invariance over waves is assumed (TRUE) or not (FALSE). Whereas this does not change the difference in the df for comparisons to the restricted model, power might be affected.
#' @param ... other parameters specifying the factor model (see [semPower.genSigma()]) and the type of power analysis 
#' @return a list containing the results of the power analysis, the population covariance matrix Sigma and mean vector mu, the H0 implied matrix SigmaHat and mean vector muHat, as well as various lavaan model strings (modelH0, and modelH1)
#' @examples
#' \dontrun{
#'  
#' }
#' @seealso [semPower.genSigma()]
#' @export
semPower.powerCLPM <- function(type, comparison = 'restricted',
                               nWaves = NULL, 
                               autoregEffects = NULL, crossedEffects = NULL, 
                               rXY = NULL,
                               waveEqual = NULL, 
                               nullEffect = NULL, nullWhich = NULL,
                               metricInvariance = TRUE,
                               ...){
  
  # TODO: lagged effects would be nice
  # TODO: do we need autocorrelated residuals?
  
  comparison <- checkComparisonModel(comparison)
  
  # we override Beta and Sigma later, so let's make sure it is not set in ellipsis argument
  if('Beta' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Beta.')
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma.')
  
  #validate input
  if(is.null(autoregEffects) ||  is.null(crossedEffects)) stop('autoregEffects and crossedEffects may not be NULL.')
  if(is.null(nWaves) || is.na(nWaves) || nWaves < 2) stop('nWaves must be >= 2.')
  if(is.null(rXY)) rXY <- rep(0, nWaves)
  if(length(rXY) != nWaves) stop('rXY must be of length nWaves')
  invisible(lapply(rXY, function(x) checkBounded(x, 'All rXY ', bound = c(-1, 1), inclusive = FALSE)))
  if(!is.list(autoregEffects)) autoregEffects <- list(rep(autoregEffects[1], (nWaves - 1)), rep(autoregEffects[2], (nWaves - 1)))
  if(!is.list(crossedEffects)) crossedEffects <- list(rep(crossedEffects[1], (nWaves - 1)), rep(crossedEffects[2], (nWaves - 1)))
  invisible(lapply(autoregEffects, function(x) lapply(x, function(x) checkBounded(x, 'All autoregEffects ', bound = c(-1, 1), inclusive = FALSE))))
  invisible(lapply(crossedEffects, function(x) lapply(x, function(x) checkBounded(x, 'All autoregEffects ', bound = c(-1, 1), inclusive = FALSE))))
  if(length(autoregEffects) != length(crossedEffects) || (length(crossedEffects) != 2 && length(crossedEffects) != (nWaves - 1))) stop('autoregEffects and crossedEffects must be of length nWaves - 1 or be of length 2.')
  
  if(!is.null(waveEqual)){
    waveEqual <- unlist(lapply(waveEqual, function(x) tolower(trimws(x))))
    if(any(unlist(lapply(waveEqual, function(x) !x %in% c('autoregx', 'autoregy', 'crossedx', 'crossedy', 'corxy'))))) stop('waveEqual may only contain autoregX, autoregY, crossedX, crossedY, corXY')
  }
  
  if(is.null(nullEffect)) stop('nullEffect must be defined.')
  # we do not allow stacking of hypotheses. there might be a use case for this,
  # but this would complicate defining the relevant parameter when these vary across waves. 
  if(length(nullEffect) > 1) stop('nullEffect must contain a single hypothesis')
  nullEffect <- unlist(lapply(nullEffect, function(x) tolower(trimws(x))))
  nullEffect <- gsub(" ", "", nullEffect, fixed = TRUE)
  nullValid <- c('autoregx', 'autoregy', 'crossedx', 'crossedy', 'corxy',
                 'autoregx=0', 'autoregy=0', 'crossedx=0', 'crossedy=0',
                 'autoregx=autoregy', 'crossedx=crossedy', 'corxy=0')
  if(any(unlist(lapply(nullEffect, function(x) !x %in% nullValid)))) stop('Unknown value for nullEffect')
  if(any(nullEffect %in% waveEqual)) stop('You cannot set the same parameters in nullEffect and waveEqual')
  if(nWaves == 2 && nullEffect %in% c('autoregx', 'autoregy','crossedx', 'crossedy', 'corxy')) stop('for two waves, there is only one crossedX and crossedY effect, only one autoregressive effect each, and only one X-Y residual correlation. Did you mean crossedX = 0 or autoregX = 0?')
  
  if(is.null(nullWhich) && nWaves == 2) nullWhich <- 1
  if(is.null(nullWhich) && nWaves > 2){
    msg <- 'nullWhich must be defined when there are more than 2 waves and relevant parameters are not constant across waves'
    if(is.null(waveEqual) && !nullEffect %in% c('autoregx', 'autoregy', 'crossedx', 'crossedy')) stop(msg) 
    if(!'autoregx' %in% waveEqual && nullEffect %in% c('autoregx=0', 'autoregx=autoregy')) stop(msg) 
    if(!'autoregy' %in% waveEqual && nullEffect %in% c('autoregy=0', 'autoregx=autoregy')) stop(msg) 
    if(!'crossedx' %in% waveEqual && nullEffect %in% c('crossedx=0', 'crossedx=crossedy')) stop(msg) 
    if(!'crossedy' %in% waveEqual && nullEffect %in% c('crossedy=0', 'crossedx=crossedy')) stop(msg) 
    if(!'corxy' %in% waveEqual && nullEffect %in% c('corxy=0')) stop(msg) 
    nullWhich <- 1 # this should be the proper default for all remaining cases
  }
  if(!is.null(nullWhich)){
    if(!is.numeric(nullWhich) || length(nullWhich) > 1) stop('nullWhich must be a single number.')
    if(nullWhich < 1 || (nullEffect != 'corxy=0' && nullWhich > (nWaves - 1))) stop('nullWhich must lie between 1 and nWaves - 1.')
  }
  
  ### create B
  B <- matrix(0, ncol = 2*nWaves, nrow = 2*nWaves)
  # add autoregEffects and crossed-effects
  for(i in 1:(nWaves - 1)){
    xidx <- 2 + 2*(i - 1) + 1
    yidx <- xidx + 1
    # autoregEffects
    B[xidx, (xidx - 2)] <- autoregEffects[[1]][i]
    B[yidx, (yidx - 2)] <- autoregEffects[[2]][i]
    # crossed effects
    B[yidx, (xidx - 2)] <- crossedEffects[[1]][i]
    B[xidx, (yidx - 2)] <- crossedEffects[[2]][i]
  }
  
  ### create Psi
  Psi <- diag(ncol(B))
  if(any(rXY != 0)){
    for(i in 1:nWaves){
      Psi[2*i, (2*i - 1)] <- Psi[(2*i - 1), 2*i] <- rXY[i]
    }
  }
  
  # add metric invariance constrains to analysis model
  metricInvarianceList <- NULL
  if(metricInvariance){
    metricInvarianceList <- list(
      seq(1, 2*nWaves, 2),
      seq(2, 2*nWaves, 2)  
    )
  }
  
  ### get Sigma
  generated <- semPower.genSigma(Beta = B, Psi = Psi, 
                                 useReferenceIndicator = TRUE, 
                                 metricInvariance = metricInvarianceList, 
                                 ...)

  ### create model strings
  model <- generated[['modelTrueCFA']]
  
  # add CLPM structure 
  for(f in 3:ncol(B)){     # omit rows 1:2
    fidx <- which(B[f, ] != 0)
    if(length(fidx) != 0){
      tok <- paste0('f', f, ' ~ ', paste(paste0('pf', paste0(formatC(f, width = 2, flag = 0), formatC(fidx, width = 2, flag = 0)), '*'), paste0('f', fidx), sep = '', collapse = ' + '))
      model <- paste(model, tok, sep='\n')
    }
  }
  # add (residual) correlations 
  model <- paste(model, 'f1 ~~ pf0201*f2', sep='\n')
  for(i in 2:nWaves){
    tok <- paste0('f',(2*i - 1),' ~~ ', paste0('pf', paste0(formatC(2*i, width = 2, flag = 0), formatC(2*i - 1, width = 2, flag = 0)), '*'), 'f', 2*i)
    model <- paste(model, tok, sep='\n')
  }
  
  ### define H1 and H0 model
  # first get constraints that may be part of either model
  tokStabX <- tokStabY <- tokCrossedX <- tokCrossedY <- tokCorXY <- ''
  # we also do this for autoregx=0 and autoregx=autoregy, because we need pStabX later; tokStabX is only used for autoregx 
  if('autoregx' %in% waveEqual || nullEffect %in% c('autoregx', 'autoregx=0', 'autoregx=autoregy')){
    xw <- seq(2*nWaves - 1, 2, -2)
    pStabX <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(xw - 2, width = 2, flag = 0))
    for(i in 1:(length(pStabX) - 1)){
      for(j in (i + 1):length(pStabX)){
        tokStabX <- paste(tokStabX, paste0(pStabX[i], '==', pStabX[j]), sep = '\n')
      }  
    }
    pStabX <- pStabX[order(pStabX)]
  }
  if('autoregy' %in% waveEqual || nullEffect %in% c('autoregy', 'autoregy=0', 'autoregx=autoregy')){
    yw <- seq(2*nWaves, 3, -2)
    pStabY <- paste0('pf', formatC(yw, width = 2, flag = 0), formatC(yw - 2, width = 2, flag = 0))
    for(i in 1:(length(pStabY) - 1)){
      for(j in (i + 1):length(pStabY)){
        tokStabY <- paste(tokStabY, paste0(pStabY[i], '==', pStabY[j]), sep = '\n')
      }  
    }
    pStabY <- pStabY[order(pStabY)]
  }
  # we also do this for crossedx=0 and crossedx=crossedy, because we need pCrossedX later; tokCrossedX is only used for crossedx 
  if('crossedx' %in% waveEqual || nullEffect %in% c('crossedx', 'crossedx=0', 'crossedx=crossedy')){  
    xw <- seq(2*nWaves - 3, 0, -2)
    yw <- seq(2*nWaves, 3, -2)
    pCrossedX <- paste0('pf', formatC(yw, width = 2, flag = 0), formatC(xw, width = 2, flag = 0))
    for(i in 1:(length(pCrossedX) - 1)){
      for(j in (i + 1):length(pCrossedX)){
        tokCrossedX <- paste(tokCrossedX, paste0(pCrossedX[i], '==', pCrossedX[j]), sep = '\n')
      }  
    }
    pCrossedX <- pCrossedX[order(pCrossedX)]
  }
  if('crossedy' %in% waveEqual || nullEffect %in% c('crossedy', 'crossedy=0', 'crossedx=crossedy')){
    xw <- seq(2*nWaves - 1, 2, -2)
    yw <- seq(2*nWaves - 2, 1, -2)
    pCrossedY <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(yw, width = 2, flag = 0))
    for(i in 1:(length(pCrossedY) - 1)){
      for(j in (i + 1):length(pCrossedY)){
        tokCrossedY <- paste(tokCrossedY, paste0(pCrossedY[i], '==', pCrossedY[j]), sep = '\n')
      }  
    }
    pCrossedY <- pCrossedY[order(pCrossedY)]
  }
  if('corxy' %in% waveEqual || nullEffect %in% c('corxy', 'corxy=0')){
    xw <- seq(2*nWaves - 1, 2, -2)
    yw <- seq(2*nWaves, 3, -2)
    pCorXY <- paste0('pf', formatC(yw, width = 2, flag = 0), formatC(xw, width = 2, flag = 0))
    for(i in 1:(length(pCorXY) - 1)){
      for(j in (i + 1):length(pCorXY)){
        tokCorXY <- paste(tokCorXY, paste0(pCorXY[i], '==', pCorXY[j]), sep = '\n')
      }  
    }
    pCorXY <- pCorXY[order(pCorXY)]
  }
  
  # add constraints to H1 model
  modelH1 <- model
  if(!is.null(waveEqual)){
    if('autoregx' %in% waveEqual) modelH1 <- paste(modelH1, tokStabX, sep = '\n')
    if('autoregy' %in% waveEqual) modelH1 <- paste(modelH1, tokStabY, sep = '\n')
    if('crossedx' %in% waveEqual) modelH1 <- paste(modelH1, tokCrossedX, sep = '\n')
    if('crossedy' %in% waveEqual) modelH1 <- paste(modelH1, tokCrossedY, sep = '\n')
    if('corxy' %in% waveEqual) modelH1 <- paste(modelH1, tokCorXY, sep = '\n')
  }
  
  ## add constraints to H0 model
  modelH0 <- modelH1  
  # modelH1 constraints are not in nullEffect, so ask again for each type: 
  if('autoregx' %in% nullEffect) modelH0 <- paste(modelH0, tokStabX, sep = '\n')
  if('autoregy' %in% nullEffect) modelH0 <- paste(modelH0, tokStabY, sep = '\n')
  if('crossedx' %in% nullEffect) modelH0 <- paste(modelH0, tokCrossedX, sep = '\n')
  if('crossedy' %in% nullEffect) modelH0 <- paste(modelH0, tokCrossedY, sep = '\n')
  if('corxy' %in% nullEffect) modelH0 <- paste(modelH0, tokCorXY, sep = '\n')
  if('autoregx=0' %in% nullEffect){
    tok <- paste0(pStabX[nullWhich], ' == 0')
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('autoregy=0' %in% nullEffect){
    tok <- paste0(pStabY[nullWhich], ' == 0')
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('crossedx=0' %in% nullEffect){
    tok <- paste0(pCrossedX[nullWhich], ' == 0')
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('crossedy=0' %in% nullEffect){
    tok <- paste0(pCrossedY[nullWhich], ' == 0')
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('autoregx=autoregy' %in% nullEffect){
    tok <- paste0(pStabX[nullWhich], ' == ', pStabY[nullWhich])
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('crossedx=crossedy' %in% nullEffect){
    tok <- paste0(pCrossedX[nullWhich], ' == ', pCrossedY[nullWhich])
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('corxy=0' %in% nullEffect){
    pCorXY <- c('pf0201', pCorXY)   # add exog cor
    tok <- paste0(pCorXY[nullWhich], ' == 0')
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  
  # here we actually fit modelH1 in case of a restricted comparison
  # because we cannot be sure that user input yields perfectly fitting h1 models, 
  # when there are additional constraints (waveequal or invariance)
  # maybe it makes sense to throw a warning if the h1 model yields f > 0 
  if(comparison == 'saturated') modelH1 <- NULL
  
  semPower.powerLav(type, 
                    Sigma = generated[['Sigma']],
                    modelH0 = modelH0, 
                    modelH1 = modelH1, 
                    ...)
}


