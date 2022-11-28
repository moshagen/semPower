##########################  convenience functions  #####################

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
#' @param ... other parameters specifying the factor model (see [semPower.genSigma()]) and the type of power analysis 
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
  comparison <- checkComparisonModel(comparison)
  
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

#' semPower.powerRegression
#'
#' Convenience function for performing power analysis on slope(s) in a latent regression.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of 'a-priori', 'post-hoc', 'compromise'
#' @param comparison comparison model, one of 'saturated' or 'restricted'. This determines the df for power analyses. 'Saturated' provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. 'Restricted' provides power to reject the model when compared to a model that just restricts the indirect effect to zero, so the df are always 1.
#' @param slope vector of standardized slopes (or a single number for a single predictor) predicting Y. 
#' @param nullSlope single number indicating which of the slope(s) is hypothesized to equal zero, defaults to 1. 
#' @param corXX correlation(s) between the k predictors (X). Either NULL, a single number (for k = 2 predictors), or a matrix. If NULL, the predictors are uncorrelated. 
#' @param ... other parameters specifying the factor model (see [semPower.genSigma()]), where the first factor corresponds to Y, and the type of power analysis 
#' @return a list containing the results of the power analysis, Sigma and SigmaHat, the implied loading matrix (lambda), as well as several lavaan model strings (modelPop, modelTrue, and modelAna) 
#' @examples
#' \dontrun{
#' # latent regression of the form Y = .2*X1 + .3*X2, where X1 and X2 correlate by .4
#' # request power for the hypothesis that the slope of X1 ist zero. 
#' # providing the number of indicators by factor (Y, X1, X2) each loading by the same magnitude on its designed factor.
#' regPower <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), corXX = .4,
#'                                      nIndicator = c(3, 5, 4),
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#' summary(regPower$power)
#' 
#' # same as above, but ask for power to detect the  slope of X2
#' regPower <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), nullSlope = 2, corXX = .4,
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
#'                                      nIndicator = c(3, 5, 4),
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'                                      
#' }
#' @seealso [semPower.genSigma()]
#' @export
semPower.powerRegression <- function(type, comparison = 'restricted',
                                     slopes = NULL, nullSlope = 1, corXX = NULL, 
                                     ...){
  
  comparison <- checkComparisonModel(comparison)
  
  # we override Phi and Sigma later, so let's make sure it is not set in ellipsis argument
  if('Phi' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Phi, because the factor correlations depend on corXX and the slopes.')
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma, because Sigma is determined as function of corXX and the slopes.')
  
  # validate input
  if(is.null(slopes)) stop('slopes cannot be NULL.')
  if(!is.vector(slopes) & !is.matrix(slopes)) stop('slopes must be a single number of a vector')
  invisible(lapply(slopes, function(x) checkBounded(x, 'All slopes ', bound = c(-1, 1), inclusive = TRUE)))
  if(sum(slopes^2) > 1) stop('slopes imply a negative residual variance for Y, make sure that the sum of the squared slopes is < 1')
  if(!is.matrix(slopes)) slopes <- matrix(slopes, nrow = length(slopes))
  if(nullSlope < 1 | nullSlope > nrow(slopes)) stop('nullSlope is invalid.')
  
  if(is.null(corXX)) corXX <- diag(nrow(slopes)) 
  if(is.vector(corXX) & length(corXX) > 1) stop('corXX must be a single number or a matrix') 
  if(!is.matrix(corXX)){
    corXX <- matrix(corXX, nrow = 2, ncol = 2) 
    diag(corXX) <- 1
  } 
  checkPositiveDefinite(corXX)
  if(ncol(corXX) != nrow(slopes)) stop('Dimension of corXX does not match number of predictors.')
  
  # calc implied sigma  
  corXY <- (corXX %*% slopes)
  Phi <- t(c(1, corXY))
  Phi <- rbind(Phi, cbind(corXY, corXX))
  generated <- semPower.genSigma(Phi = Phi, useReferenceIndicator = TRUE, ...)
  Sigma <- generated$Sigma
  
  ### create ana model string
  # add regressions 
  model <- paste(generated$modelTrue, 
                 paste0('f1 ~ ', paste0(paste0('pf',(1 + 1:ncol(corXX))), '*f',(1 + 1:ncol(corXX)), collapse = '+')), 
                 sep = '\n')
  modelTrue <- model
  modelAna <- paste(model, '\n', paste0('pf', (nullSlope + 1),' == 0'))  
  
  modelH1 <- NULL
  if(comparison == 'restricted') modelH1 <- modelTrue
  
  semPower.powerLav(type, modelH0 = modelAna, modelH1 = modelH1, Sigma = Sigma, ...)
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
#' @return a list containing the results of the power analysis, Sigma and SigmaHat, the implied loading matrix (lambda), as well as several lavaan model strings (modelPop, modelTrue, and modelAna) 
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
#'                 c(0.0, 0.0, 0.6),
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
  if(!is.null(Beta) & (!is.null(bYX) | !is.null(bMX) | !is.null(bYM))) stop('Either provide bYX, bMX, and bYM or provide Beta, but not both.')
  if(is.null(Beta)){
    if(is.null(bYX) | is.null(bMX) | is.null(bYM)) stop('Provide bYX, bYM, and bYM or provide Beta')
    if(length(bYX) != 1) stop('bYX must be a single slope (X -> Y)')
    if(length(bMX) != 1) stop('bMX must be a single slope (X -> M)')
    if(length(bYM) != 1) stop('bYM must be a single slope (M -> Y)')
    invisible(lapply(c(bYX, bMX, bYM), function(x) checkBounded(x, 'All slopes ', bound = c(-1, 1), inclusive = TRUE)))
    if((bYX^2 + bYM^2) > 1) stop('bYX and bYM imply a negative residual variance for Y, make sure that the sum of the squared slopes on Y is < 1')
    if(bMX == 0 | bYM == 0) stop('One of bMX and bYM is zero, implying the indirect effect is zero. The indirect effect must differ from zero.')
    indirect <- list(c(2,1), c(3,2))
  }
  
  if(!is.null(Beta)){
    if(is.null(indirect)) stop('indirect must not be NULL when Beta is defined.')
    if(isSymmetric(Beta)) stop('Beta must be symmetric.')
    invisible(apply(Beta, c(1, 2), function(x) checkBounded(x, 'All elements in Beta', bound = c(-1, 1), inclusive = TRUE)))
    if(any(diag(Beta) != 0)) stop('All diagonal elements of Beta must be zero.')
    # negative implied residual variances are checked later
    if(any(lapply(indirect, function(x) length(x)) != 2)) stop('Indirect must be a list containing vectors of size two each')
    if(any(unlist(lapply(indirect, function(x) any(x > ncol(B)))))) stop('At least one element in indirect is an out of bounds index concerning B')
    if(any(unlist(lapply(indirect, function(x) B[x[1], x[2]])) == 0)) stop('Beta and indirect imply an indirect effect of zero. The indirect effect must differ from zero.')
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
  # transform to standard cfa model by converting B to implied phi
  Phi <- getPhi.B(B) 
  generated <- semPower.genSigma(Phi = Phi, useReferenceIndicator = TRUE, ...)
  Sigma <- generated$Sigma
  
  ### create ana model string
  model <- generated$modelTrue
  # add mediation structure
  for(f in 1:ncol(B)){
    fidx <- which(B[f, ] != 0)
    if(length(fidx) != 0){
      tok <- paste0('f', f, ' ~ ', paste(paste0('pf', paste0(f, fidx), '*'), paste0('f', fidx), sep = '', collapse = ' + '))
      model <- paste(model, tok, sep='\n')
    }
  }
  # add indirect effects
  ind <- unlist(lapply(indirect, function(x) paste0('pf', paste0(x, collapse = ''))))
  modelTrue <- paste(model, '\n', 'ind := ', paste(ind, collapse = '*'))
  
  # lav doesn't like constraining the indirect effect to zero. 
  # we instead constrain the smallest of the contained direct effects, so this
  # actually gives power for a single slope. however, this seems to closely reflect
  # power for the indirect effect (and indeed works much better than using ind=0 as 
  # comparison model, which grossly overestimates the true effect)
  #modelAna <- paste(modelTrue, '\n', 'ind == 0')  
  mb <- paste0('pf', paste(which(B == min(B[B != 0]), arr.ind = TRUE)[1, ], collapse = ''))
  modelAna <- paste(modelTrue, '\n', paste0(mb,' == 0'))  

  modelH1 <- NULL
  if(comparison == 'restricted') modelH1 <- modelTrue
  
  semPower.powerLav(type, modelH0 = modelAna, modelH1 = modelH1, Sigma = Sigma, ...)
}



