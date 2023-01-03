#' semPower.powerLav
#'
#' Perform a power analysis given `lavaan` model strings defining the H0 and the H1 model based on either 
#' a `lavaan` model string defining the population model or the population covariance matrix Sigma and the population means mu.
#' This requires the `lavaan` package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param modelPop `lavaan` model string defining the true model. Can be omitted when `Sigma` is set.
#' @param modelH0 `lavaan` model string defining the (incorrect) analysis model.
#' @param modelH1 `lavaan` model string defining the comparison model. If omitted, the saturated model is the comparison model.
#' @param fitH1model whether to fit the H1 model. If `FALSE`, the H1 model is assumed to show the same fit as the saturated model, and only the delta df are computed.
#' @param Sigma can be used instead of `modelPop`: population covariance matrix. A list for multiple group models.
#' @param mu can be used instead of `modelPop`: vector of population means. Can be omitted for no meanstructure. A list for multiple group models.
#' @param simulatedPower whether to perform a simulated (`TRUE`, rather than analytical, `FALSE`) power analysis. See [simulate()] for additional options.
#' @param lavOptions a list of additional options passed to `lavaan`, e. g., `list(estimator = 'mlm')` to request robust ML estimation. Mostly useful in conjunction with `simulatedPower`. 
#' @param lavOptionsH1 alternative options passed to `lavaan` that are only used for the H1 model. If `NULL`, identical to `lavOptions`. Probably only useful for multigroup models.
#' @param ... other parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()].
#' @return A list containing the following components is returned:
#' \item{`power`}{the results of the power analysis. Use the `summary` method to obtain formatted results.}
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
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
#' mH0 <- '
#'   f1 =~ x1 + x2 + x3
#'   f2 =~ x4 + x5 + x6
#'   f1 ~~ 0*f2
#' '
#' # do a priori power analsis
#' lavpower.ap <- semPower.powerLav(type = 'a-priori', 
#'                                  modelPop = mPop, modelH0 = mH0,
#'                                  alpha = .05, beta = .05)
#' summary(lavpower.ap$power)
#'
#' }
#' @seealso [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @importFrom utils installed.packages
#' @export
semPower.powerLav <- function(type, 
                              modelPop = NULL, 
                              modelH0 = NULL, modelH1 = NULL, fitH1model = TRUE, 
                              Sigma = NULL, mu = NULL, 
                              simulatedPower = FALSE, 
                              lavOptions = NULL, lavOptionsH1 = lavOptions, 
                              ...){
  
  # check whether lavaan is available
  if(!'lavaan' %in% rownames(installed.packages())) stop('This function depends on the lavaan package, so install lavaan first.')
  # validate input
  type <- checkPowerTypes(type)
  if(is.null(modelH0)) stop('Provide a lavaan model string defining the analysis (H0) model.')
  if(is.null(modelPop) && is.null(Sigma)) stop('Either provide a lavaan model string defining the population model or provide the population covariance matrix Sigma.')
  if(!is.null(modelPop) && !is.null(Sigma)) stop('Either provide a lavaan model string defining the population model or provide the population covariance matrix Sigma, but not both.')
  if(simulatedPower && type == 'compromise') stop('Simulated power is not available for compromise power analysis, because this would require a vast (infeasible) number of simulation runs to yield reliable results.')
  if(!is.null(modelPop) && !is.list(modelPop)) modelPop <- list(modelPop)
  if(!is.null(Sigma) && !is.list(Sigma)) Sigma <- list(Sigma)

  # determine population Sigma / mu
  if(is.null(Sigma)){
    Sigma <- lapply(modelPop, function(x) lavaan::fitted(lavaan::sem(x))[['cov']])
    mu <- lapply(modelPop, function(x) lavaan::fitted(lavaan::sem(x))[['mean']])
  }
  
  # analytical power
  if(!simulatedPower){
    
    # we need to call lavaan() directly with defaults as defined in sem()
    lavOptions <- getLavOptions(lavOptions, nGroups = length(Sigma))
    if(!is.null(lavOptions[['estimator']]) && toupper(lavOptions[['estimator']]) != "ML") stop('Analytical power is only available with ML estimation. Note that power based on ML derivatives (mlm etc) is asymptotically identical.')

    # get H0 sigmaHat / muHat
    modelH0Fit <- do.call(lavaan::lavaan,
                          append(list(model = modelH0,
                                      sample.cov = if(length(Sigma) > 1) Sigma else Sigma[[1]],
                                      sample.mean = if(length(Sigma) > 1) mu else mu[[1]]),
                                 lavOptions))
    if(!modelH0Fit@optim[['converged']]) stop('The H0 model did not converge.')
    if(length(Sigma) > 1){
      # multigroup case
      SigmaHat <- lapply(1:length(Sigma), function(x) lavaan::fitted(modelH0Fit)[[x]][['cov']])
      muHat <- lapply(1:length(Sigma), function(x) lavaan::fitted(modelH0Fit)[[x]][['mean']])
    }else{
      # single group case
      SigmaHat <- list(lavaan::fitted(modelH0Fit)[['cov']])
      muHat <- list(lavaan::fitted(modelH0Fit)[['mean']])
    }
    df <- dfH0 <- modelH0Fit@test[['standard']][['df']]  # this is probably invalid for estm with adjusted df
    
    # get H1 comparison model and deltaF
    if(!is.null(modelH1) && fitH1model){
      lavOptionsH1 <- getLavOptions(lavOptionsH1, nGroups = length(Sigma))
      modelH1Fit <- do.call(lavaan::lavaan,
                            append(list(model = modelH1,
                                        sample.cov = if(length(Sigma) > 1) Sigma else Sigma[[1]],
                                        sample.mean = if(length(Sigma) > 1) mu else mu[[1]]),
                                   lavOptionsH1))
      if(!modelH1Fit@optim[['converged']]) stop('The H1 model did not converge.')
      dfH1 <- modelH1Fit@test[['standard']][['df']]
      if(dfH1 >= dfH0) stop('The df of the H0 model are not larger than the df of the H1 model, as they should be.')
      # get delta F
      if(length(Sigma) > 1){
        # multigroup case
        fminH0 <- lapply(1:length(Sigma), 
                         function(x) getF.Sigma(lavaan::fitted(modelH0Fit)[[x]][['cov']], Sigma[[x]], 
                                                lavaan::fitted(modelH0Fit)[[x]][['mean']], mu[[x]]))
        fminH1 <- lapply(1:length(Sigma), 
                         function(x) getF.Sigma(lavaan::fitted(modelH1Fit)[[x]][['cov']], Sigma[[x]], 
                                                lavaan::fitted(modelH1Fit)[[x]][['mean']], mu[[x]]))
        deltaF <- lapply(1:length(Sigma), function(x) fminH0[[x]] - fminH1[[x]]) # result must be a list
      }else{
        # single group case
        fminH0 <- getF.Sigma(lavaan::fitted(modelH0Fit)[['cov']], Sigma[[1]], 
                             lavaan::fitted(modelH0Fit)[['mean']], mu[[1]])
        fminH1 <- getF.Sigma(lavaan::fitted(modelH1Fit)[['cov']], Sigma[[1]], 
                             lavaan::fitted(modelH1Fit)[['mean']], mu[[1]])
        deltaF <- fminH0 - fminH1
      }
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
                      # simulate() takes care of proper lavOptions
                      lavOptions = lavOptions,  
                      lavOptionsH1 = lavOptionsH1, 
                      ...)
    SigmaHat <- muHat <- NULL
  }

  # remove list structure for single group models
  if(length(Sigma) == 1){
    Sigma <- Sigma[[1]]
    SigmaHat <- SigmaHat[[1]]
    mu <- mu[[1]]
    muHat <- muHat[[1]]
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
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'Saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'Restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param Phi either a single number defining the correlation between exactly two factors or the factor correlation matrix. A list for multiple group models.
#' @param nullEffect defines the hypothesis of interest, must be one of `'cor = 0'` (the default) to test whether a correlation is zero, `'corX = corZ'` to test for the equality of correlations, and `'corA = corB'` to test for the equality of a correlation across groups. Define the correlations to be set to equality in `nullWhich` and the groups in `nullWhichGroups`. 
#' @param nullWhich vector of size 2 indicating which factor correlation in `Phi` is hypothesized to equal zero when `nullEffect = 'cor = 0'`, or to restrict to equality across groups when `nullEffect = 'corA = corB'`, or list of vectors defining which correlations to restrict to equality when `nullEffect = 'corX = corZ'`. Can also contain more than two correlations, e.g., `list(c(1, 2), c(1, 3), c(2, 3))` to set `Phi[1, 2] = Phi[1, 3] = Phi[2, 3]`. If omitted, the correlation between the first and the second factor is targeted, i. e. `nullWhich = c(1, 2)`.
#' @param nullWhichGroups for `nullEffect = 'corA = corB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param ... other parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and specifying the factor model (see [semPower.genSigma()]). See details.
#' @return A list containing the following components is returned:
#' \item{`power`}{the results of the power analysis. Use the `summary` method to obtain formatted results.}
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' @details 
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the *definition of the factor model*:
#' * `Lambda`: The factor loading matrix (with the number of columns equaling the number of factors).
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading in the second factor by .8, .6, .6, and .4..
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`, defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Defines the mean loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the mean loadings of the first factor to equal .5 and those of the second factor do equal .6.
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` need to be defined. 
#' If the model contains observed variables only, use `Lambda = diag(x)` where `x` is the number of variables.
#' 
#' Additional arguments related to the requested type of *power analysis*:
#' * `alpha`: The alpha error probability. Required for `type = 'a-priori'` and `type = 'post-hoc'`.
#' * Either `beta` or `power`: The beta error probability and the statistical power (1 - beta), respectively. Only for `type = 'a-priori'`.
#' * `N`: The sample size. Always required for `type = 'post-hoc'` and `type = 'compromise'`. For `type = 'a-priori'` and multiple group analysis, `N` is a list of group weights.
#' * `abratio`: The ratio of alpha to beta. Only for `type = 'compromise'`. 
#' 
#' Optional arguments if a *simulated power analysis* (`simulatedPower = TRUE`) is requested:
#' * `nReplications`: The number of simulation runs. Defaults to 250, but larger numbers greatly improve accuracy at the expense of increased computation time.
#' * `minConvergenceRate`: The required minimum convergence rate. Defaults to .50.
#' 
#' @examples
#' \dontrun{
#' # a priori power analysis only providing the number of indicators to define 
#' # two factors with correlation of phi and same loading for all indicators
#' cfapower.ap <- semPower.powerCFA(type = 'a-priori',
#'                                  nullWhich = c(1, 2), 
#'                                  Phi = .2, nIndicator = c(5, 6), loadM = .5,
#'                                  alpha = .05, beta = .05)
#' summary(cfapower.ap$power)
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
#' cfapower.ph <- semPower.powerCFA(type = 'post-hoc',
#'                                  nullWhich = c(1, 2), 
#'                                  Phi = Phi, loadings = loadings,
#'                                  alpha = .05, N = 250)
#' 
#' # multigroup case to test that there are no group differences 
#' # concerning the correlation between two factors   
#' cfapower.ph <- semPower.powerCFA(type = 'post-hoc', comparison = 'restricted', 
#'                                  nullEffect = 'corA=corB',
#'                                  Phi = list(.2, .1), loadM = .5, 
#'                                  nIndicator = c(3, 3), 
#'                                  alpha = .05, N = c(250, 250))
#' 
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerCFA <- function(type, comparison = 'restricted', 
                              Phi = NULL,
                              nullEffect = 'cor = 0',
                              nullWhich = NULL, 
                              nullWhichGroups = NULL, 
                              ...){
  
  # validate input
  comparison <- checkComparisonModel(comparison)
  if(is.null(Phi)) stop('Phi must be defined')
  nullEffect <- checkNullEffect(nullEffect, c('cor=0', 'corx=corz', 'cora=corb'))
  if(!is.null(nullWhichGroups) && !is.list(Phi)) stop('Phi must be provided for each group.')
  if(nullEffect == 'cora=corb' && !is.list(Phi)) stop('corA=corB refers to muligroup analysis, so Phi must be a list.')
  if(is.list(Phi) && !is.null(nullWhichGroups)) lapply(as.list(nullWhichGroups), function(x) checkBounded(x, bound(1, length(Phi)))) 
    
  # generate sigma 
  generated <- semPower.genSigma(Phi = Phi, ...)
  
  ### now do validation of nullWhich, since we now know Phi
  isMultigroup <- is.list(Phi)
  if(isMultigroup) nfac <- ncol(generated[[1]][['Phi']]) else nfac <- ncol(generated[['Phi']]) 
  if(is.null(nullWhich) && nfac == 2) nullWhich <- c(1, 2)
  if(is.null(nullWhich)) stop('nullWhich must be defined.')
  if(!is.list(nullWhich)) nullWhich <- list(nullWhich)
  if(any(unlist(lapply(nullWhich, function(x) length(x) != 2)))) stop('nullWhich may only contain vectors of size two.')
  if(any(unlist(lapply(nullWhich, function(x) x[1] == x[2])))) stop('elements in nullWhich may not refer to variances.')
  if(any(unlist(lapply(nullWhich, function(x) (x[1] < 1 | x[2] < 1 | x[1] > nfac | x[2] > nfac))))) stop('At least one element in nullWhich is an out of bounds index concerning Phi.')
  if(length(nullWhich) > 1){
    for(i in 1:(length(nullWhich) - 1)){
      for(j in (i + 1):length(nullWhich)){
        if(nullWhich[[i]][1] == nullWhich[[j]][1] && nullWhich[[i]][2] == nullWhich[[j]][2]) stop('elements in nullWhich may not refer to the same correlation')
      }
    }
  }
  
  ### H0 model
  if(isMultigroup) modelH0 <- generated[[1]][['modelTrueCFA']] else modelH0 <- generated[['modelTrueCFA']] 
  if(nullEffect == 'cor=0'){
    modelH0 <- paste(c(modelH0,
      paste0('f', nullWhich[[1]], collapse = ' ~~ 0*')),
      collapse = '\n')
  }else if(nullEffect == 'corx=corz'){
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
  }else if(nullEffect == 'cora=corb'){
    if(is.null(nullWhichGroups)) nullWhichGroups <- 1:length(Phi)
    lab <- rep('NA', length(Phi))
    lab[nullWhichGroups] <- 'pf1'
    lab <- paste0('c(', paste(lab, collapse = ','), ')*')
    modelH0 <- paste(c(modelH0,
                       paste0('f', nullWhich[[1]], collapse = paste0(' ~~ ', lab))),
                       collapse = '\n')
  }else{
    stop('nullEffect not defined')
  }

  # we always enforce invariance constraints in the multigroup case
  lavOptions <- NULL
  if(isMultigroup) lavOptions <- list(group.equal = c('loadings', 'lv.variances'))

  modelH1 <- NULL
  if(comparison == 'restricted'){
    if(isMultigroup) modelH1 <- generated[[1]][['modelTrueCFA']] else modelH1 <- generated[['modelTrueCFA']] 
    # single group case: the h1 model always fits perfectly
    # multigroup case: we cannot be sure that user input yields a perfectly fitting model
    fitH1model <- isMultigroup 
  } 
  
  if(isMultigroup) Sigma <- lapply(generated, '[[', 'Sigma') else Sigma <- generated[['Sigma']] 
  
  semPower.powerLav(type = type,
                    Sigma = Sigma,
                    modelH0 = modelH0,
                    modelH1 = modelH1,
                    fitH1model = fitH1model,
                    lavOptions = lavOptions,
                    ...)
  
}

#' semPower.powerRegression
#'
#' Convenience function for performing power analysis on slope(s) in a latent regression of the form Y = XB.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'Saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'Restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param slopes vector of standardized slopes (or a single number for a single slope) of the k predictors for Y. A list of slopes for multigroup models.
#' @param corXX correlation(s) between the k predictors (X). Either `NULL` for uncorrelated predictors, a single number (for k = 2 predictors), or a matrix. Can also be a list for multigroup models providing the correlations by group of matrices (otherwise, the same correlations are used in all groups). 
#' @param nullEffect defines the hypothesis of interest, must be one of `'slope = 0'` (the default) to test whether a slope is zero, `'slopeX = slopeZ'` to test for the equality of slopes, or `'slopeA = slopeB'` to test for the equality of slopes across groups. Define the slopes to set to equality in `nullWhich`.
#' @param nullWhich single number indicating which slope is hypothesized to equal zero when `nullEffect = 'slope = 0'`, or indicating which slope to restrict to equality across groups when `nullEffect = 'slopeA = slopeB'`, or vector defining the slopes to restrict to equality when `nullEffect = 'slopeX = slopeZ'`. Can also contain more than two slopes, e.g. `c(1, 2, 3)` to constrain the first three slopes to equality.
#' @param nullWhichGroups for `nullEffect = 'slopeA = slopeB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param ... other parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and specifying the factor model (see [semPower.genSigma()]). Note the first factor is treated as Y and the subsequent factors as the predictors X_k.
#' @return A list containing the following components is returned:
#' \item{`power`}{the results of the power analysis. Use the `summary` method to obtain formatted results.}
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
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
#' # multigroup example                                      
#' regPower <- semPower.powerRegression(type = 'post-hoc',
#'                                      slopes = list(c(.2, .3, .4), c(.2, .0, .4)), 
#'                                      corXX = corXX, 
#'                                      nullEffect = 'slopeA = slopeB', 
#'                                      nullWhich = 2,
#'                                      nIndicator = c(4, 3, 5, 4),
#'                                      loadM = c(.5, .5, .6, .7), 
#'                                      N = list(250, 250),
#'                                      alpha = .05, beta = .05)
#'                                      
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerRegression <- function(type, comparison = 'restricted',
                                     slopes = NULL, 
                                     corXX = NULL, 
                                     nullEffect = 'slope = 0',
                                     nullWhich = NULL,
                                     nullWhichGroups = NULL,
                                     ...){
  
  comparison <- checkComparisonModel(comparison)
  
  # we override Phi and Sigma later, so let's make sure it is not set in ellipsis argument
  if('Phi' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Phi, because the factor correlations depend on corXX and the slopes.')
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma, because Sigma is determined as function of corXX and the slopes.')
  
  # validate input
  nullEffect <- checkNullEffect(nullEffect, c('slope=0', 'slopex=slopez', 'slopea=slopeb'))
  if(is.null(slopes)) stop('slopes cannot be NULL.')
  if(nullEffect == 'slopea=slobeb' && !is.list(slopes)) stop('slopes must be a list when a multiple group analysis is requested.')
  if(is.list(slopes)) if(length(unique(unlist(lapply(slopes, length)))) > 1) stop('the same number of slopes must be provided for each group.')
  if(!is.list(slopes)) slopes <- list(slopes)
  if(!is.null(corXX) && !is.list(corXX)) corXX <- list(corXX)

  if(any(!unlist(lapply(slopes, is.vector))) && any(!unlist(lapply(slopes, is.matrix)))) stop('slopes must be a single number or a vector')
  lapply(slopes, function(y) invisible(lapply(y, function(x) checkBounded(x, 'All slopes ', bound = c(-1, 1), inclusive = TRUE))))
  if(any(unlist(lapply(slopes, function(x) sum(x^2) > 1)))) stop('slopes imply a negative residual variance for Y, make sure that the sum of the squared slopes is < 1')
  slopes <- lapply(slopes, function(x) if(!is.matrix(x)) matrix(x, nrow = length(x)) else x)
  isMultigroup <- length(slopes) > 1

  if(is.null(corXX)) corXX <- lapply(slopes, function(x) diag(nrow(x))) 
  if(any(unlist(lapply(corXX, is.vector)) && unlist(lapply(corXX, length)) > 1)) stop('corXX must be a single number or a matrix') 
  if(isMultigroup && length(corXX) == 1) corXX <- rep(corXX, length(slopes)) # assume same corXX for all groups
  corXX <- lapply(corXX, function(x) {
    if(!is.matrix(x)){
      xx <- matrix(x, nrow = 2, ncol = 2) 
      diag(xx) <- 1
      xx
    }else{
      x
    }
  })
  invisible(lapply(corXX, checkPositiveDefinite))
  if(any(unlist(lapply(seq_along(corXX), function(x) ncol(corXX[[x]]) != nrow(slopes[[x]]))))) stop('Dimension of corXX does not match number of predictors.')
  
  if(is.null(nullWhich)) stop('nullWhich must be defined.')
  if(any(nullWhich < 1) || any(nullWhich > nrow(slopes[[1]]))) stop('nullWhich is invalid.')
  if(nullEffect == 'slopex=slopez'){
    if(length(nullWhich) < 2 || length(nullWhich) > nrow(slopes[[1]])) stop('nullWhich must contain at least two slopes when nullEffect is slopex=slopez, but not more slopes than available')
  }else{
    if(length(nullWhich) > 1) stop('nullWhich must be a single number when nullEffect is slope=0')
  }
  nullWhich <- nullWhich + 1 # because first factor is criterion
  
  # calc implied sigma. 
  # for single groups, phi is defined directly, because this special case is simpler than defining B and calling getPhi.B  
  if(!isMultigroup){
    corXY <- (corXX[[1]] %*% slopes[[1]])
    Phi <- t(c(1, corXY))
    Phi <- rbind(Phi, cbind(corXY, corXX[[1]]))
    generated <- semPower.genSigma(Phi = Phi, useReferenceIndicator = TRUE, ...)
    # for multigroup, the approach does not work, because we need to go unstandardized 
  }else{
    B <- lapply(slopes, function(x){
      cB <- matrix(0, ncol = (length(x) + 1), nrow = (length(x) + 1))
      cB[1, ] <- c(0, x)  # we want the first factor endogenous
      cB
    })
    Psi <- lapply(seq_along(slopes), function(x){
      cPsi <- diag(ncol(B[[x]]))
      cPsi[2:nrow(cPsi), 2:ncol(cPsi)] <- corXX[[x]]
      cPsi
    })
    generated <- semPower.genSigma(B = B, Psi = Psi, 
                                   useReferenceIndicator = TRUE, ...)
  }
  

  ### create ana model string
  # add regressions 
  if(isMultigroup) model <- generated[[1]][['modelTrueCFA']] else model <- generated[['modelTrueCFA']]
  np <- (1 + 1:ncol(corXX[[1]]))
  
  if(nullEffect == 'slopex=slopez'){
    tok <- ''
    for(i in 1:(length(nullWhich) - 1)){
      for(j in (i + 1):length(nullWhich)){
        tok <- paste(tok, paste0('pf', nullWhich[i], ' == ', 'pf', nullWhich[j]), sep = '\n')
      }
    }
    modelH0 <- paste(model, 
                     paste0('f1 ~ ', paste0(paste0('pf', np), '*f', np, collapse = ' + ')),
                     tok,
                     sep = '\n')
  }else if(nullEffect == 'slope=0'){
    modelH0 <- paste(model, 
                     paste0('f1 ~ ', paste0(paste0('pf', np), '*f', np, collapse = ' + ')),
                     paste0('pf', nullWhich,' == 0'),
                     sep = '\n')
  }else if(nullEffect == 'slopea=slopeb'){
    if(is.null(nullWhichGroups)) nullWhichGroups <- seq_along(slopes)
    lab <- rep('NA', length(slopes))
    lab[nullWhichGroups] <- 'pf1'
    lab <- paste0('c(', paste(lab, collapse = ','), ')*')
    modelH0 <- paste(model, 
                     paste0('f1 ~ ', paste0(lab, 'f', np[(nullWhich - 1)], ' + '), paste0('f',np[-(nullWhich - 1)], collapse = ' + ')),
                     sep = '\n')
  }else{
    stop('nullEffect not defined.')
  }

  # we always enforce invariance constraints in the multigroup case
  lavOptions <- NULL
  if(isMultigroup) lavOptions <- list(group.equal = c('loadings', 'lv.variances'))
  
  modelH1 <- NULL
  if(comparison == 'restricted'){
    if(!isMultigroup){
      modelH1 <- paste(model, 
                       paste0('f1 ~ ', paste0(paste0('pf',(1 + 1:ncol(corXX[[1]]))), '*f',(1 + 1:ncol(corXX[[1]])), collapse = ' + ')),
                       sep = '\n')
    }else{
      # no slope labels in multigroup case
      modelH1 <- paste(model, 
                       paste0('f1 ~ ', paste0('f',(1 + 1:ncol(corXX[[1]])), collapse = ' + ')),
                       sep = '\n')
      
    }
    # single group case: the h1 model always fits perfectly
    # multigroup case: we cannot be sure that user input yields a perfectly fitting model
    fitH1model <- isMultigroup 
  } 

  if(isMultigroup) Sigma <- lapply(generated, '[[', 'Sigma') else Sigma <- generated[['Sigma']] 
  
  semPower.powerLav(type, 
                    Sigma = Sigma, 
                    modelH0 = modelH0, 
                    modelH1 = modelH1, 
                    fitH1model = fitH1model,
                    lavOptions = lavOptions,
                    ...)
}


#' semPower.powerMediation
#'
#' Convenience function for performing power analysis concerning indirect effect(s) in a mediation model.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'Saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'Restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param bYX the standardized slope (direct effect) for X -> Y. A list for multiple group models. Can be `NULL` if `Beta` is set.
#' @param bMX the standardized slope for X -> M. A list for multiple group models. Can be `NULL` if `Beta` is set.
#' @param bYM the standardized slope for M -> Y. A list for multiple group models. Can be `NULL` if `Beta` is set.
#' @param Beta can be used instead of `bYX`, `bMX`, and `bYM`: matrix of regression weights connecting the latent factors (all-Y notation). Exogenous variables must be in the first row(s), so the upper triangular of Beta must be zero. A list for multiple group models.
#' @param indirect `NULL` unless `Beta` is set. Otherwise a list of vectors of size 2 indicating the elements of `Beta` that define the indirect effect of interest, e.g. `list(c(2, 1), c(3, 2))`. See details.
#' @param nullEffect defines the hypothesis of interest, must be one of `'ind = 0'` (the default) to test whether the indirect effect is zero or `'indA = indB'` to test for the equality of indirect effects across groups. See details.
#' @param nullWhichGroups for `nullEffect = 'indA = indB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param ... other parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and specifying the factor model (see [semPower.genSigma()]).  Note that in case of a simple mediation, the order of factors is X, M, Y. 
#' @return A list containing the following components is returned:
#' \item{`power`}{the results of the power analysis. Use the `summary` method to obtain formatted results.}
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' @details
#' Notes on implementation:
#' * For models without latent variables, `nullEffect = 'ind = 0'` and `nullEffect = 'indA = indB'` constrain the indirect effect to zero and to equality, respectively.
#' * For models with latent variables and `nullEffect = 'ind = 0'`, power is approximated by constraining the smallest slope contained in the indirect effect to zero. 
#' * For models with latent variables multiple groups (i.e., `nullEffect = 'indA = indB'`), there is currently no way to determine power, 
#' because implementing equality constrains on the indirect effects leads to non-convergence and 
#' the approach to constrain a single or all slopes to equality across groups misrepresents the actual hypothesis of interest. 
#' 
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
#' # same mediation model as above, but assuming single loading of 1 for all factors 
#' # (=> mediation model with manifest variables)
#' medPower <- semPower.powerMediation(type = 'a-priori', 
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(1, 1, 1), loadM = 1,
#'                                     alpha = .05, beta = .05)
#' 
#' # same latent mediation model as above, but specifying loadings through Lambda 
#' Lambda <- matrix(c(
#'                  #  X,   M,   Y
#'                 c(0.5, 0.0, 0.0),    
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
#'               c(.00, .00, .00, .00),       # X
#'               c(.20, .00, .00, .00),       # M1
#'               c(.00, .30, .00, .00),       # M2
#'               c(.00, .00, .40, .00)        # Y
#'               ), byrow = TRUE, ncol = 4)
#' # only define primary loadings by factor
#' loadings <- list(
#'                c(0.4, 0.5, 0.8),           # X
#'                c(0.7, 0.6, 0.5, 0.8),      # M1
#'                c(0.5, 0.6, 0.3, 0.4, 0.6), # M2
#'                c(0.6, 0.7, 0.8)            # Y
#'                )
#'
#' medPower <- semPower.powerMediation(type = 'a-priori', 
#'                                     Beta = B, indirect = list(c(2,1), c(3,2), c(4,3)),
#'                                     loadings = loadings,
#'                                     alpha = .05, beta = .05)
#'  
#' # multigroup example
#' medPower <- semPower.powerMediation(type = 'a-priori', 
#'                                     nullEffect = 'indA = indB',
#'                                     bYX = list(.25, .25), bMX = list(.3, .3), bYM = list(.4, .5),
#'                                     Lambda = diag(3), N = list(1, 1),
#'                                     alpha = .05, beta = .05)
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerMediation <- function(type, comparison = 'restricted',
                                    bYX = NULL, bMX = NULL, bYM = NULL,
                                    Beta = NULL, indirect = NULL, 
                                    nullEffect = 'ind = 0',
                                    nullWhichGroups = NULL, 
                                    ...){
  
  comparison <- checkComparisonModel(comparison)
  
  # we override Phi and Sigma later, so let's make sure it is not set in ellipsis argument
  if('Phi' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Phi, because the factor correlations depend on Beta (or the slopes).')
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma, because Sigma is determined as function of Beta (or the slopes).')
  
  # validate input
  nullEffect <- checkNullEffect(nullEffect, c('ind=0', 'inda=indb'))
  if(!is.null(Beta) && (!is.null(bYX) || !is.null(bMX) || !is.null(bYM))) stop('Either provide bYX, bMX, and bYM or provide Beta, but not both.')
  if(is.null(Beta)){
    if(is.null(bYX) || is.null(bMX) || is.null(bYM)) stop('Provide bYX, bYM, and bYM or provide Beta')
    isMultigroup <- is.list(bYX) && length(bYX) > 1
    if(!is.list(bYX)) bYX <- list(bYX)
    if(!is.list(bMX)) bMX <- list(bMX)
    if(!is.list(bYM)) bYM <- list(bYM)
    if(length(unique(unlist(lapply(list(bYX, bMX, bYM), length)))) != 1) stop('bYX, bYM, and bYM must be of same lenght in multiple group analysis.')
    if(any(unlist(lapply(bYX, length)) != 1)) stop('Each bYX must be a single slope (X -> Y)')
    if(any(unlist(lapply(bMX, length)) != 1)) stop('Each bMX must be a single slope (X -> M)')
    if(any(unlist(lapply(bYM, length)) != 1)) stop('Each bYX must be a single slope (M -> Y)')
    if(!isMultigroup) invisible(lapply(c(bYX, bMX, bYM), function(x) checkBounded(x, 'All slopes ', bound = c(-1, 1), inclusive = TRUE)))
    if(any(unlist(bYX)^2 + unlist(bYM)^2 > 1)) stop('bYX and bYM imply a negative residual variance for Y, make sure that the sum of the squared slopes on Y is < 1')
    if(!isMultigroup && (bMX == 0 || bYM == 0)) stop('One of bMX and bYM is zero, implying the indirect effect is zero. The indirect effect must differ from zero.')
    indirect <- list(c(2, 1), c(3, 2))
  }
  
  if(!is.null(Beta)){
    isMultigroup <- is.list(Beta) && length(Beta) > 1
    if(!is.list(Beta)) Beta <- list(Beta)
    if(any(unlist(lapply(Beta, function(x) any(diag(x) != 0) )))) stop('All diagonal elements of Beta must be zero.')
    lapply(Beta, function(y) invisible(apply(y, c(1, 2), function(x) checkBounded(x, 'All elements in Beta', bound = c(-1, 1), inclusive = TRUE))))
    if(isMultigroup && (length(unique(unlist(lapply(Beta, ncol)))) > 1 || length(unique(unlist(lapply(Beta, nrow)))) > 1)) stop('Beta must be of same dimension for all groups') 
    # negative implied residual variances are checked in getPhi.B
    if(is.null(indirect)) stop('indirect must not be NULL when Beta is defined.')
    if(any(lapply(indirect, function(x) length(x)) != 2)) stop('Indirect must be a list containing vectors of size two each')
    if(any(unlist(lapply(indirect, function(x) any(x > ncol(Beta)))))) stop('At least one element in indirect is an out of bounds index concerning B')
    if(any(unlist(lapply(indirect, function(x) Beta[[1]][x[1], x[2]])) == 0)) stop('Beta and indirect imply an indirect effect of zero. The indirect effect must differ from zero.')
  }
  if(nullEffect == 'inda=indb' && !isMultigroup) stop('bYX, bYM, and bYM or Beta must be a list for multiple group comparisons.')
  if(nullEffect == 'ind=0' && isMultigroup) stop('Multiple group models are only valid for nullEffect = inda=indb.')
  
  B <- Beta
  if(is.null(B)){
    B <- lapply(seq_along(bMX), function(x){
      matrix(c(
        c(0, 0, 0),                # X
        c(bMX[[x]], 0, 0),         # M
        c(bYX[[x]], bYM[[x]], 0)   # Y
      ), byrow = TRUE, ncol = 3)
    })
  }
  
  if(!is.null(nullWhichGroups)) lapply(nullWhichGroups, function(x) checkBounded(x, 'All elements in nullWhichGroups'), bound = c(1, length(B)), inclusive = TRUE)
  
    
  ### get Sigma
  # we want completely standardized slopes, so transform to standard cfa model by converting B to implied phi
  Phi <- lapply(B, getPhi.B)
  if(!isMultigroup) Phi <- Phi[[1]]
  
  generated <- semPower.genSigma(Phi = Phi, useReferenceIndicator = TRUE, ...)
  if(!isMultigroup) isObserved <- ncol(generated[['Sigma']]) == ncol(Phi) else isObserved <- ncol(generated[[1]][['Sigma']]) == ncol(Phi[[1]])
  
  ### create model strings
  if(!isMultigroup) model <- generated[['modelTrueCFA']] else model <- generated[[1]][['modelTrueCFA']]
  # add mediation structure
  for(f in 1:ncol(B[[1]])){
    fidx <- unique(unlist(lapply(B, function(x) which(x[f, ] != 0))))
    if(length(fidx) != 0){
      clab <- paste0('pf', paste0(formatC(f, width = 2, flag = 0), formatC(fidx, width = 2, flag = 0)))
      if(isMultigroup) clab <- unlist(lapply(clab, function(x) paste0('c(', paste0(x, 'g', seq_along(B), collapse = ', '), ')')))
      tok <- paste0('f', f, ' ~ ', paste(clab, paste0('*f', fidx), sep = '', collapse = ' + '))
      model <- paste(model, tok, sep='\n')
    }
  }
  # add indirect effects
  ind <- unlist(lapply(indirect, function(x) paste0('pf', paste0(formatC(x, width = 2, flag = 0), collapse = ''))))
  if(!isMultigroup){
    model <- paste(model, '\n', 'ind := ', paste(ind, collapse = '*'))
  }else{
    ind <- lapply(seq_along(B), function(x) paste0(ind, 'g', x))
    tok <- paste0('ind', seq_along(B), ' := ', unlist(lapply(ind, function(x) paste(x, collapse = '*'))), collapse = '\n')
    model <- paste(model, '\n', tok)
  }

  # lav doesn't like constraining the indirect effect to zero or to equality for latent variable models, 
  # so we apply different approaches depending on whether there are latent variables
  if(nullEffect == 'ind=0'){
    if(isObserved){
      # constrain indirect effect
      modelH0 <- paste(model, '\n', 'ind == 0')
    }else{
      # constrain the smallest of the contained direct effects, so this
      # actually gives power for a single slope. however, this seems to closely reflect
      # power for the indirect effect (and indeed works much better than using ind=0 as 
      # comparison model, which grossly overestimates the true effect)
      cs <- indirect[[which.min(unlist(lapply(indirect, function(x) B[[1]][x[1], x[2]])))]]
      mb <- paste0('pf', paste(formatC(cs, width = 2, flag = 0), collapse = ''))
      modelH0 <- paste(model, '\n', paste0(mb,' == 0'))  
    }
  }else if(nullEffect == 'inda=indb'){
    if(isObserved){
      # constrain indirect effect
      if(is.null(nullWhichGroups)) nullWhichGroups <- seq_along(B)
      indeffects <- paste0('ind', nullWhichGroups)
      tok <- list()
      for(i in 1:(length(indeffects) - 1)){
        for(j in (i + 1):length(indeffects)){
          tok <- append(tok, paste0(indeffects[i], ' == ', indeffects[j]))
        }
      }
      modelH0 <- paste(c(model, unlist(tok)), collapse = '\n')
    }else{
      # setting indirect effects to equality does not work 
      stop('Multigroup comparisons of indirect effects are not supported for latent variable models.')
    }
  }else{
    stop('nullEffect not defined.')
  }

  # enforce invariance constraints in the multigroup case
  lavOptions <- NULL
  if(isMultigroup) lavOptions <- list(group.equal = c('loadings', 'lv.variances'))
  
  modelH1 <- NULL
  if(comparison == 'restricted'){
    # h1 model always fits perfectly, only needed for delta df
    modelH1 <- model
    # single group case: the h1 model always fits perfectly
    # multigroup case: we cannot be sure that user input yields a perfectly fitting model
    fitH1model <- isMultigroup 
  } 
  
  if(isMultigroup) Sigma <- lapply(generated, '[[', 'Sigma') else Sigma <- generated[['Sigma']] 
  
  semPower.powerLav(type, 
                    Sigma = Sigma, 
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
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'Saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'Restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param nWaves number of waves, must be >= 2.
#' @param autoregEffects vector of the autoregressive effects of X and Y (constant across waves), or a list of vectors of autoregressive effects for X and Y from wave to wave, e.g. `list(c(.7, .6), c(.5, .5))` for a autoregressive effect of .7 for x1->x2 and .6 for x2->x3
#' @param crossedEffects vector of crossed effects of X on Y (X -> Y) and vice versa (both constant across waves), or a list of vectors of crossed effects giving the crossed effect of X on Y (and vice versa) for each wave, e. g. `list(c(.2, .3), c(.1, .1))` for x1->y2 = .2 and x2->y3 = .3.
#' @param rXY vector of (residual-)correlations between X and Y for each wave. If NULL, all (residual-)correlations are zero. 
#' @param waveEqual parameters that are assumed to be equal across waves in both the H0 and the H1 model. Valid are `'autoregX'` and `'autoregY'` for autoregressive effects, `'crossedX'` and `'crossedY'` for crossed effects, `'corXY'` for residual correlations, or `NULL` for none (so that all parameters are freely estimated, subject to the constraints defined in `nullEffect`). 
#' @param nullEffect defines the hypothesis of interest. Valid are the same arguments as in `waveEqual` and additionally `'autoregX = 0'`, `'autoregY = 0'`, `'crossedX = 0'`, `'crossedY = 0'` to constrain the X or Y autoregressive effects or the crossed effects to zero, `'autoregX = autoregY'` and `'crossedX = crossedY'` to constrain them to be equal for X and Y.
#' @param nullWhich used in conjunction with `nullEffect` to identify which parameter to constrain when there are > 2 waves and parameters are not constant across waves. For example, `nullEffect = 'autoregX = 0'` with `nullWhich = 2` would constrain the second autoregressive effect for X to zero.    
#' @param metricInvariance whether metric invariance over waves is assumed (`TRUE`, the default) or not (`FALSE`). This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model, where the df are not affected  by invariance constraints).
#' @param ... other parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and specifying the factor model (see [semPower.genSigma()]). Note that the order of factors is (x1, y1, x2, y2, ..., x_nWaves, y_nWaves). 
#' @return A list containing the following components is returned:
#' \item{`power`}{the results of the power analysis. Use the `summary` method to obtain formatted results.}
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' @examples
#' \dontrun{
#'  
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
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
  if(is.null(autoregEffects) || is.null(crossedEffects)) stop('autoregEffects and crossedEffects may not be NULL.')
  if(is.null(nWaves) || is.na(nWaves) || nWaves < 2) stop('nWaves must be >= 2.')
  if(is.null(rXY)) rXY <- rep(0, nWaves)
  if(length(rXY) != nWaves) stop('rXY must be of length nWaves')
  invisible(lapply(rXY, function(x) checkBounded(x, 'All rXY ', bound = c(-1, 1), inclusive = FALSE)))
  if(!is.list(autoregEffects)) autoregEffects <- list(rep(autoregEffects[1], (nWaves - 1)), rep(autoregEffects[2], (nWaves - 1)))
  if(!is.list(crossedEffects)) crossedEffects <- list(rep(crossedEffects[1], (nWaves - 1)), rep(crossedEffects[2], (nWaves - 1)))
  invisible(lapply(autoregEffects, function(x) lapply(x, function(x) checkBounded(x, 'All autoregEffects ', bound = c(-1, 1), inclusive = FALSE))))
  invisible(lapply(crossedEffects, function(x) lapply(x, function(x) checkBounded(x, 'All autoregEffects ', bound = c(-1, 1), inclusive = FALSE))))
  if(length(autoregEffects) != length(crossedEffects) || (length(crossedEffects) != 2 && length(crossedEffects) != (nWaves - 1))) stop('autoregEffects and crossedEffects must be of length nWaves - 1 or be of length 2.')
  if(is.list(autoregEffects)) if(length(autoregEffects[[1]]) != length(autoregEffects[[2]])) stop('autoregEffects for X and Y must be of equal length.')
  if(is.list(autoregEffects)) if(length(crossedEffects[[1]]) != length(crossedEffects[[2]])) stop('CrossedEffects for X and Y must be of equal length.')
  if(is.list(autoregEffects)) if(length(autoregEffects[[1]]) != length(crossedEffects[[2]])) stop('autoregEffects and crossedEffects must be of equal length.')  
  if(is.list(autoregEffects)) if(length(autoregEffects[[1]]) != (nWaves - 1)) stop('autoregEffects must be of length nWaves - 1.')  
  if(is.list(autoregEffects)) if(length(crossedEffects[[1]]) != (nWaves - 1)) stop('crossedEffects must be of length nWaves - 1.')   
  
  if(!is.null(waveEqual)){
    waveEqual <- unlist(lapply(waveEqual, function(x) tolower(trimws(x))))
    if(any(unlist(lapply(waveEqual, function(x) !x %in% c('autoregx', 'autoregy', 'crossedx', 'crossedy', 'corxy'))))) stop('waveEqual may only contain autoregX, autoregY, crossedX, crossedY, corXY')
  }
  
  # we do not allow stacking of hypotheses. there might be a use case for this,
  # but this would complicate defining the relevant parameter when these vary across waves. 
  nullValid <- c('autoregx', 'autoregy', 'crossedx', 'crossedy', 'corxy',
                 'autoregx=0', 'autoregy=0', 'crossedx=0', 'crossedy=0',
                 'autoregx=autoregy', 'crossedx=crossedy', 'corxy=0')
  nullEffect <- checkNullEffect(nullEffect, nullValid)

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
  tokAutoregX <- tokAutoregY <- tokCrossedX <- tokCrossedY <- tokCorXY <- ''
  # we also do this for autoregx=0 and autoregx=autoregy, because we need pAutoregX later; tokAutoregX is only used for autoregx 
  if('autoregx' %in% waveEqual || nullEffect %in% c('autoregx', 'autoregx=0', 'autoregx=autoregy')){
    xw <- seq(2*nWaves - 1, 2, -2)
    pAutoregX <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(xw - 2, width = 2, flag = 0))
    for(i in 1:(length(pAutoregX) - 1)){
      for(j in (i + 1):length(pAutoregX)){
        tokAutoregX <- paste(tokAutoregX, paste0(pAutoregX[i], '==', pAutoregX[j]), sep = '\n')
      }  
    }
    pAutoregX <- pAutoregX[order(pAutoregX)]
  }
  if('autoregy' %in% waveEqual || nullEffect %in% c('autoregy', 'autoregy=0', 'autoregx=autoregy')){
    yw <- seq(2*nWaves, 3, -2)
    pAutoregY <- paste0('pf', formatC(yw, width = 2, flag = 0), formatC(yw - 2, width = 2, flag = 0))
    for(i in 1:(length(pAutoregY) - 1)){
      for(j in (i + 1):length(pAutoregY)){
        tokAutoregY <- paste(tokAutoregY, paste0(pAutoregY[i], '==', pAutoregY[j]), sep = '\n')
      }  
    }
    pAutoregY <- pAutoregY[order(pAutoregY)]
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
    if('autoregx' %in% waveEqual) modelH1 <- paste(modelH1, tokAutoregX, sep = '\n')
    if('autoregy' %in% waveEqual) modelH1 <- paste(modelH1, tokAutoregY, sep = '\n')
    if('crossedx' %in% waveEqual) modelH1 <- paste(modelH1, tokCrossedX, sep = '\n')
    if('crossedy' %in% waveEqual) modelH1 <- paste(modelH1, tokCrossedY, sep = '\n')
    if('corxy' %in% waveEqual) modelH1 <- paste(modelH1, tokCorXY, sep = '\n')
  }
  
  ## add constraints to H0 model
  modelH0 <- modelH1  
  # modelH1 constraints are not in nullEffect, so ask again for each type: 
  if('autoregx' %in% nullEffect) modelH0 <- paste(modelH0, tokAutoregX, sep = '\n')
  if('autoregy' %in% nullEffect) modelH0 <- paste(modelH0, tokAutoregY, sep = '\n')
  if('crossedx' %in% nullEffect) modelH0 <- paste(modelH0, tokCrossedX, sep = '\n')
  if('crossedy' %in% nullEffect) modelH0 <- paste(modelH0, tokCrossedY, sep = '\n')
  if('corxy' %in% nullEffect) modelH0 <- paste(modelH0, tokCorXY, sep = '\n')
  if('autoregx=0' %in% nullEffect){
    tok <- paste0(pAutoregX[nullWhich], ' == 0')
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('autoregy=0' %in% nullEffect){
    tok <- paste0(pAutoregY[nullWhich], ' == 0')
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
    tok <- paste0(pAutoregX[nullWhich], ' == ', pAutoregY[nullWhich])
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



#' semPower.powerRICLPM
#'
#' Convenience function for performing power analysis on effects in a random intercept cross-lagged panel model (RI-CLPM).
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'Saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'Restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param nWaves number of waves, must be >= 3.
#' @param autoregEffects vector of the autoregressive effects of X and Y (constant across waves), or a list of vectors of autoregressive effects for X and Y from wave to wave, e.g. `list(c(.7, .6), c(.5, .5))` for an autoregressive effect of .7 for X1->X2 and .6 for X2->X3 and autoregressive effects of .5 for Y1->Y2 and Y2 -> Y3
#' @param crossedEffects vector of crossed effects of X on Y (X -> Y) and vice versa (both constant across waves), or a list of vectors of crossed effects giving the crossed effect of X on Y (and vice versa) for each wave, e.g. `list(c(.2, .3), c(.1, .1))` for X1->Y2 = .2, X2->Y3 = .3, Y1->Y2 = .1, and Y2->Y3 = .1.
#' @param rXY vector of (residual-)correlations between X and Y for each wave. If `NULL`, all (residual-)correlations are zero. 
#' @param rBXBY correlation between random intercept factors.
#' @param waveEqual parameters that are assumed to be equal across waves in both the H0 and the H1 model. Valid are `'autoregX'` and `'autoregY'` for autoregressive effects, `'crossedX'` and `'crossedY'` for crossed effects, `'corXY'` for residual correlations, or `NULL` for none (so that all parameters are freely estimated). 
#' @param nullEffect defines the hypothesis of interest. Valid are the same arguments as in `waveEqual` and additionally `'autoregX = 0'`, `'autoregY = 0'`, `'crossedX = 0'`, `'crossedY = 0'` to constrain the X or Y autoregressive effects or the crossed effects to zero, `'corBXBY = 0'` to constrain the correlation between the random intercepts to zero, and `'autoregX = autoregY'` and `'crossedX = crossedY'` to constrain them to be equal for X and Y.
#' @param nullWhich used in conjunction with `nullEffect` to identify which parameter to constrain when there are > 2 waves and parameters are not constant across waves. For example, `nullEffect = 'autoregX = 0'` with `nullWhich = 2` would constrain the second autoregressive effect for X to zero.    
#' @param metricInvariance whether metric invariance over waves is assumed (`TRUE`, the default) or not (`FALSE`). This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model, where the df are not affected by invariance constraints).
#' @param ... other parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and specifying the factor model (see [semPower.genSigma()]). Note that the order of factors is (x1, y1, x2, y2, ..., x_nWaves, y_nWaves), see details. 
#' @details 
#' Specification of the factor model assumes the following ordering:
#' * `Lambda`: Columns should be in order X1, Y1, X2, Y2, ..., X_nWaves, Y_nWaves. 
#' * `loadings`: List of vectors providing the factor loadings for each factor ordered by wave, e.g., list(c(.2, .2, .2), c(.4, .4, .4, .4), c(.2, .2, .2), c(.4, .4, .4, .4), c(.2, .2, .2), c(.4, .4, .4, .4)) to define loadings of .2 for the three indicators of X at waves 1-3 and loadings of .4 for the four indicators of Y at waves 1-3. Must not contain secondary loadings.   
#' * `nIndicator` Vector indicating the number of indicators for each factor ordered by wave, e.g. c(3, 4, 3, 4, 3, 4) to define three indicators for factor X at waves 1-3 and four indicators for factor Y at waves 1-3.
#' * `loadM` Vector giving mean loadings for each factor ordered by wave, e.g., c(.5, .6, .5, .6, .5, .6) to define loadings of .5 for X at waves 1-3 and loadings of .6 for Y at waves 1-3; or single number to use for every loading.
#' 
#' @return A list containing the following components is returned:
#' \item{`power`}{the results of the power analysis. Use the `summary` method to obtain formatted results.}
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' @examples
#' \dontrun{
#'  
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerRICLPM <- function(type, comparison = 'restricted',
                                 nWaves = NULL, 
                                 autoregEffects = NULL, crossedEffects = NULL, 
                                 rXY = NULL,
                                 rBXBY = NULL,
                                 waveEqual = NULL, 
                                 nullEffect = NULL, nullWhich = NULL,
                                 metricInvariance = TRUE,
                                 ...){
  
  # TODO: do we need autocorrelated residuals?
  
  comparison <- checkComparisonModel(comparison)
  
  # we override Beta and Sigma later, so let's make sure it is not set in ellipsis argument
  if('Beta' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Beta.')
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma.')
  
  # validate input
  if(is.null(autoregEffects) ||  is.null(crossedEffects)) stop('autoregEffects and crossedEffects may not be NULL.')
  if(is.null(nWaves) | is.na(nWaves) | nWaves < 3) stop('nWaves must be >= 3.')
  if(is.null(rXY)) rXY <- rep(0, nWaves)
  if(length(rXY) != nWaves) stop('rXY must be of length nWaves')
  invisible(lapply(rXY, function(x) checkBounded(x, 'All rXY ', bound = c(-1, 1), inclusive = FALSE)))
  if(is.null(rBXBY)) rBXBY <- 0
  if(length(rBXBY) != 1) stop('rBXBY must contain a single number.')
  checkBounded(rBXBY, 'rBXBY ', bound = c(-1, 1), inclusive = FALSE)
  if(!is.list(autoregEffects)) autoregEffects <- list(rep(autoregEffects[1], (nWaves - 1)), rep(autoregEffects[2], (nWaves - 1)))
  if(!is.list(crossedEffects)) crossedEffects <- list(rep(crossedEffects[1], (nWaves - 1)), rep(crossedEffects[2], (nWaves - 1)))
  invisible(lapply(autoregEffects, function(x) lapply(x, function(x) checkBounded(x, 'All autoregEffects ', bound = c(-1, 1), inclusive = FALSE))))
  invisible(lapply(crossedEffects, function(x) lapply(x, function(x) checkBounded(x, 'All autoregEffects ', bound = c(-1, 1), inclusive = FALSE))))
  if(length(autoregEffects) != length(crossedEffects) || (length(crossedEffects) != 2 && length(crossedEffects) != (nWaves - 1))) stop('autoregEffects and crossedEffects must be of length nWaves - 1 or be of length 2.')
  if(is.list(autoregEffects)) if(length(autoregEffects[[1]]) != length(autoregEffects[[2]])) stop('autoregEffects for X and Y must be of equal length.')
  if(is.list(autoregEffects)) if(length(crossedEffects[[1]]) != length(crossedEffects[[2]])) stop('CrossedEffects for X and Y must be of equal length.')
  if(is.list(autoregEffects)) if(length(autoregEffects[[1]]) != length(crossedEffects[[2]])) stop('autoregEffects and crossedEffects must be of equal length.')  
  if(is.list(autoregEffects)) if(length(autoregEffects[[1]]) != (nWaves - 1)) stop('autoregEffects must be of length nWaves - 1.')  
  if(is.list(autoregEffects)) if(length(crossedEffects[[1]]) != (nWaves - 1)) stop('crossedEffects must be of length nWaves - 1.')   
  
  if(!is.null(waveEqual)){
    waveEqual <- unlist(lapply(waveEqual, function(x) tolower(trimws(x))))
    if(any(unlist(lapply(waveEqual, function(x) !x %in% c('autoregx', 'autoregy', 'crossedx', 'crossedy', 'corxy'))))) stop('waveEqual may only contain autoregX, autoregY, crossedX, crossedY, corXY')
  }
  
  # we do not allow stacking of hypotheses. there might be a use case for this,
  # but this would complicate defining the relevant parameter when these vary across waves. 
  nullValid <- c('autoregx', 'autoregy', 'crossedx', 'crossedy', 'corxy',
                 'autoregx=0', 'autoregy=0', 'crossedx=0', 'crossedy=0',
                 'autoregx=autoregy', 'crossedx=crossedy', 'corxy=0', 'corbxby=0')
  nullEffect <- checkNullEffect(nullEffect, nullValid)
  if(any(nullEffect %in% waveEqual)) stop('You cannot set the same parameters in nullEffect and waveEqual.')
  
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
    if(nullEffect == 'corbxby=0' && nullWhich != 1) stop('If nullEffect is "corBXBY = 0", nullWhich must be 1.')
    if(nullWhich < 1 || (nullEffect != 'corxy=0' && nullWhich > (nWaves - 1))) stop('nullWhich must lie between 1 and nWaves - 1.')
  }

  
  ### create Lambda 
  args <- list(...)
  Lambda  <- args$Lambda
  if(is.null(Lambda)){
    Lambda <- genLambda(args[['loadings']], args[['nIndicator']],
                        args[['loadM']], args[['loadSD']], args[['loadMinMax']])
  }
  if(ncol(Lambda) != 2*nWaves) stop('Number of factors must be 2*nWaves.')
  
  # modify Lambda according to RI-CLPM structure
  Lambda <- cbind(matrix(0, nrow = nrow(Lambda), ncol = (2*nWaves + 2)), Lambda) # add between + within factors
  # cols: Bx, By, Wx_1, Wy_1,..., Wx_nWaves, Wy_nWaves, Fx_1, Fy_1, ..., Fx_nWaves, Fy_nWaves
  
  
  ### create B
  B <- matrix(0, ncol = (4*nWaves + 2), nrow = (4*nWaves + 2)) 
  # Bx, By, Wx_1, Wy_1,..., Wx_nWaves, Wy_nWaves, Fx_1, Fy_1, ..., Fx_nWaves, Fy_nWaves
  
  # define between factors (random intercepts)
  B[seq((3 + 2*nWaves), (4*nWaves + 2), 2), 1] <- 1 # BX
  B[seq((4 + 2*nWaves), (4*nWaves + 2), 2), 2] <- 1 # BY
  
  # define within factors
  diag(B[((3 + 2*nWaves):(4*nWaves + 2)), (3:(2 + 2*nWaves))]) <- 1 
  
  # add autoregressive effects and crossed-effects
  for(i in 1:(nWaves - 1)){
    xidx <- 2 + 2*(i - 1) + 3
    yidx <- xidx + 1
    # autoregressive effects
    B[xidx, (xidx - 2)] <- autoregEffects[[1]][i]
    B[yidx, (yidx - 2)] <- autoregEffects[[2]][i]
    # crossed effects
    B[yidx, (xidx - 2)] <- crossedEffects[[1]][i]
    B[xidx, (yidx - 2)] <- crossedEffects[[2]][i]
  }
  
  
  ### create Psi
  Psi <- diag(ncol(B))
  # Bx, By, Wx_1, Wy_1,..., Wx_nWaves, Wy_nWaves, Fx_1, Fy_1, ..., Fx_nWaves, Fy_nWaves
  
  # add cor between random intercepts
  Psi[2,1] <- Psi[1,2] <- rBXBY
  
  # set residual variance of Fx_1, ..., Fy_nWaves to 0
  diag(Psi[(3 + 2*nWaves):(4*nWaves + 2), (3 + 2*nWaves):(4*nWaves + 2)]) <- 0 
  
  # add (residual) correlations between within-factors
  if(any(rXY != 0)){
    for(i in 1:nWaves){
      Psi[(2*i + 2), (2*i + 1)] <- Psi[(2*i + 1), (2*i + 2)] <- rXY[i]
    }
  }
  
  
  # add metric invariance constrains
  metricInvarianceList <- NULL
  if(metricInvariance){
    metricInvarianceList <- list(
      seq(3 + 2*nWaves, 2 + 4*nWaves, 2),
      seq(4 + 2*nWaves, 2 + 4*nWaves, 2)  
    )
  }
  
  
  ### get model-implied sigma
  if(!is.null(args[['Lambda']])) args[['Lambda']] <- NULL # delete user-provided Lambda 
  
  generated <- do.call(what = semPower.genSigma, 
                       args = append(list(Lambda = Lambda, Beta = B, Psi = Psi,
                                          useReferenceIndicator = TRUE,
                                          metricInvariance = metricInvarianceList), args))
  
  Sigma <- generated$Sigma
  
  
  ### create ana model string
  
  # define random intercept (between) factors
  tok1 <- paste0('f1 =~ ', paste0('1*', paste0('f', seq((3 + 2*nWaves), ncol(B), 2)), collapse = ' + '))
  tok2 <- paste0('f2 =~ ', paste0('1*', paste0('f', seq((4 + 2*nWaves), ncol(B), 2)), collapse = ' + '))
  model <- paste(tok1, tok2, sep='\n')
  
  # define residualized (within) factors
  for(i in 1:(2*nWaves)){
    widx <- seq(3, 2*nWaves + 2)[i]
    fidx <- seq(2*nWaves + 3, 4*nWaves + 2)[i]
    model <- paste(model, paste0('f', widx, ' =~ 1*f', fidx), sep = '\n')
  }
  
  # add unresidualized factors
  model <- paste(model, generated$modelTrueCFA, sep='\n')
  
  # set residual variance of unresidualized factors to 0
  for(f in (2*nWaves + 3):(4*nWaves + 2)){
    tok <- paste0('f', f, ' ~~ 0*', 'f', f)
    model <- paste(model, tok, sep='\n')
  }
  
  # estimate correlation between random intercepts
  model <- paste(model, 'f1 ~~ pf0201*f2', sep='\n')
  
  # set correlation between random intercepts and residualized factors at wave 1 to 0
  model <- paste(model, 'f1 + f2 ~~ 0*f3 + 0*f4', sep='\n')
  
  # add autoregressive and cross-lagged effects
  for(f in 3:(2*nWaves + 2)){ # omit rows for random intercepts and unresidualized factors
    fidx <- which(B[f, ] != 0)
    if(length(fidx) != 0){
      tok <- paste0('f', f, ' ~ ', paste(paste0('pf', paste0(formatC(f, width = 2, flag = 0), formatC(fidx, width = 2, flag = 0)), '*'), paste0('f', fidx), sep = '', collapse = ' + '))
      model <- paste(model, tok, sep='\n')
    }
  }
  
  # add (residual) correlations 
  for(i in 1:nWaves){
    tok <- paste0('f',(1 + 2*i),' ~~ ', paste0('pf', paste0(formatC(2 + 2*i, width = 2, flag = 0), formatC(1 + 2*i, width = 2, flag = 0)), '*'), 'f', (2 + 2*i))
    model <- paste(model, tok, sep='\n')
  }


  ### define H1 and ana model
  
  # first get constraints that may be part of either model
  tok.autoregx <- tok.autoregy <- tok.crossedx <- tok.crossedy <- tok.corxy <- ''
  
  # we also do this for autoregx=0 and autoregx=autoregy, because we need p.autoregx later; tok.autoregx is only used for autoregx 
  if('autoregx' %in% waveEqual || nullEffect %in% c('autoregx', 'autoregx=0', 'autoregx=autoregy')){
    xw <- seq(2 + 2*nWaves - 1, 5, -2)
    p.autoregx <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(xw - 2, width = 2, flag = 0))
    for(i in 1:(length(p.autoregx) - 1)){
      for(j in (i + 1):length(p.autoregx)){
        tok.autoregx <- paste(tok.autoregx, paste0(p.autoregx[i], '==', p.autoregx[j]), sep = '\n')
      }  
    }
    p.autoregx <- p.autoregx[order(p.autoregx)]
  }
  if('autoregy' %in% waveEqual || nullEffect %in% c('autoregy', 'autoregy=0', 'autoregx=autoregy')){
    yw <- seq(2 + 2*nWaves, 6, -2)
    p.autoregy <- paste0('pf', formatC(yw, width = 2, flag = 0), formatC(yw - 2, width = 2, flag = 0))
    for(i in 1:(length(p.autoregy) - 1)){
      for(j in (i + 1):length(p.autoregy)){
        tok.autoregy <- paste(tok.autoregy, paste0(p.autoregy[i], '==', p.autoregy[j]), sep = '\n')
      }  
    }
    p.autoregy <- p.autoregy[order(p.autoregy)]
  }
  # we also do this for crossedx=0 and crossedx=crossedy, because we need p.crossedx later; tok.crossedX is only used for crossedx 
  if('crossedx' %in% waveEqual || nullEffect %in% c('crossedx', 'crossedx=0', 'crossedx=crossedy')){  
    xw <- seq(2 + 2*nWaves - 3, 3, -2)
    yw <- seq(2 + 2*nWaves, 6, -2)
    p.crossedx <- paste0('pf', formatC(yw, width = 2, flag = 0), formatC(xw, width = 2, flag = 0))
    for(i in 1:(length(p.crossedx) - 1)){
      for(j in (i + 1):length(p.crossedx)){
        tok.crossedx <- paste(tok.crossedx, paste0(p.crossedx[i], '==', p.crossedx[j]), sep = '\n')
      }  
    }
    p.crossedx <- p.crossedx[order(p.crossedx)]
  }
  if('crossedy' %in% waveEqual || nullEffect %in% c('crossedy', 'crossedy=0', 'crossedx=crossedy')){
    xw <- seq(2 + 2*nWaves - 1, 5, -2)
    yw <- seq(2 + 2*nWaves - 2, 4, -2)
    p.crossedy <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(yw, width = 2, flag = 0))
    for(i in 1:(length(p.crossedy) - 1)){
      for(j in (i + 1):length(p.crossedy)){
        tok.crossedy <- paste(tok.crossedy, paste0(p.crossedy[i], '==', p.crossedy[j]), sep = '\n')
      }  
    }
    p.crossedy <- p.crossedy[order(p.crossedy)]
  }
  if('corxy' %in% waveEqual || nullEffect %in% c('corxy', 'corxy=0')){
    xw <- seq(2+ 2*nWaves - 1, 5, -2)
    yw <- seq(2 + 2*nWaves, 6, -2)
    p.corxy <- paste0('pf', formatC(yw, width = 2, flag = 0), formatC(xw, width = 2, flag = 0))
    for(i in 1:(length(p.corxy) - 1)){
      for(j in (i + 1):length(p.corxy)){
        tok.corxy <- paste(tok.corxy, paste0(p.corxy[i], '==', p.corxy[j]), sep = '\n')
      }  
    }
    p.corxy <- p.corxy[order(p.corxy)]
  }
  
  ## add constraints to H1 model
  modelH1 <- model
  if(!is.null(waveEqual)){
    if('autoregx' %in% waveEqual) modelH1 <- paste(modelH1, tok.autoregx, sep = '\n')
    if('autoregy' %in% waveEqual) modelH1 <- paste(modelH1, tok.autoregy, sep = '\n')
    if('crossedx' %in% waveEqual) modelH1 <- paste(modelH1, tok.crossedx, sep = '\n')
    if('crossedy' %in% waveEqual) modelH1 <- paste(modelH1, tok.crossedy, sep = '\n')
    if('corxy' %in% waveEqual) modelH1 <- paste(modelH1, tok.corxy, sep = '\n')
  }
  
  ## add constraints to ana model
  modelH0 <- modelH1  
  # modelH1 constraints are not in nullEffect, so ask again for each type: 
  if('autoregx' %in% nullEffect) modelH0 <- paste(modelH0, tok.autoregx, sep = '\n')
  if('autoregy' %in% nullEffect) modelH0 <- paste(modelH0, tok.autoregy, sep = '\n')
  if('crossedx' %in% nullEffect) modelH0 <- paste(modelH0, tok.crossedx, sep = '\n')
  if('crossedy' %in% nullEffect) modelH0 <- paste(modelH0, tok.crossedy, sep = '\n')
  if('corxy' %in% nullEffect) modelH0 <- paste(modelH0, tok.corxy, sep = '\n')
  if('autoregx=0' %in% nullEffect){
    tok <- paste0(p.autoregx[nullWhich], ' == 0')
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('autoregy=0' %in% nullEffect){
    tok <- paste0(p.autoregy[nullWhich], ' == 0')
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('crossedx=0' %in% nullEffect){
    tok <- paste0(p.crossedx[nullWhich], ' == 0')
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('crossedy=0' %in% nullEffect){
    tok <- paste0(p.crossedy[nullWhich], ' == 0')
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('autoregx=autoregy' %in% nullEffect){
    tok <- paste0(p.autoregx[nullWhich], ' == ', p.autoregy[nullWhich])
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('crossedx=crossedy' %in% nullEffect){
    tok <- paste0(p.crossedx[nullWhich], ' == ', p.crossedy[nullWhich])
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('corxy=0' %in% nullEffect){
    p.corxy <- c('pf0403', p.corxy)   # add exog cor
    tok <- paste0(p.corxy[nullWhich], ' == 0')
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  if('corbxby=0' %in% nullEffect){
    tok <- paste0('pf0201', ' == 0')
    modelH0 <- paste(modelH0, tok, sep = '\n')
  } 
  
  # here we actually fit modelH1 in case of a restricted comparison
  # because we cannot be sure that user input yields perfectly fitting h1 models 
  # when there are additional constraints (waveequal or invariance)
  # maybe it makes sense to throw a warning if the h1 model yields f > 0 
  if(comparison == 'saturated') modelH1 <- NULL
  
  semPower.powerLav(type, 
                    modelH0 = modelH0, 
                    modelH1 = modelH1, 
                    Sigma = Sigma,
                    ...)
}

#' semPower.powerMI
#'
#' Convenience function for performing power analysis for multigroup measurement invariance models.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, either `'saturated'` or one of `'configural'`, `'metric'`, `'scalar'`, or a vector of restrictions in `lavaan` format (with `'none'` for no restrictions). See details.
#' @param nullEffect defines the hypothesis of interest. One of `'metric'`, `'scalar'`, `'residual'`, or a vector of restrictions in `lavaan` format. See details.   
#' @param ... other parameters specifying the factor model (see [semPower.genSigma()]) and the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()].
#' @return A list containing the following components is returned:
#' \item{`power`}{the results of the power analysis. Use the `summary` method to obtain formatted results.}
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' @details
#' Multigroup invariance models fit the specified model simultaneously to various groups and place increasingly
#' restrictive cross-group equality constrains on the model parameters. The typical - but not in all parts necessary -
#' sequence is (a) configural, (b) metric, (c) scalar, and (d) residual invariance, where each level of invariance is
#' compared against the previous level (e.g., scalar vs. metric). In the context, power analysis asks for the 
#' power to reject a particular level of invariance. 
#'  
#' The models defined in the `comparison` and the `nullEffect` arguments can be specified as follows:
#' \itemize{
#' \item `configural`: no invariance constraints. Shows the same fit as the saturated model, so only the delta df differ. 
#' \item `metric`: all loadings are restricted to equality. 
#' \item `scalar`: all loadings and (item-)intercepts are restricted to equality. 
#' \item `residual`: all loadings, (item-)intercepts, and (item-)residuals are restricted to equality.
#' }
#' For greater flexibility, the models can also be defined using lavaan style group.equal restrictions as a vector: 
#' \itemize{
#' \item `'none'`: no invariance constraints. Shows the same fit as the saturated model, so only the delta df differ. 
#' \item `c('loadings')`: all loadings are restricted to equality. 
#' \item `c('loadings', 'intercepts')`: all loadings and (item-)intercepts are restricted to equality. 
#' \item `c('loadings', 'intercepts', 'residuals')`: all loadings, (item-)intercepts, and (item-)residuals are restricted to equality.
#' \item `c('loadings', 'residuals')`: all loadings and (item-)residuals are restricted to equality.
#' \item `c('loadings', 'intercepts', 'means')`: all loadings, (item-)intercepts, and latent means are restricted to equality.
#' }
#' Note that variance scaling is used, so invariance of variances (`'lv.variances'`) is always met. 
#' @examples
#' \dontrun{
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerMI <- function(type, 
                             comparison = NULL,
                             nullEffect = NULL,
                             ...){
  
  # validate input
  lavGroupStrings <- c('loadings', 'intercepts', 'residuals', 'residual.covariances', 'lv.variances', 'lv.covariances','regressions', 'means')
  useLavOptions <- any(grepl(paste(lavGroupStrings, collapse = '|'), comparison)) || any(grepl(paste(lavGroupStrings, collapse = '|'), nullEffect))
  # we only check typos etc when not using lavstrings
  if(!useLavOptions){
    comparison <- checkNullEffect(comparison, c('saturated', 'configural', 'metric', 'scalar'))
    nullEffect <- checkNullEffect(nullEffect, c('metric', 'scalar', 'residual'))
    if(which(c('saturated', 'configural', 'metric', 'scalar') %in% comparison) >= 
       (2 + which(c('metric', 'scalar', 'residuals') %in% nullEffect))) stop('Model defined in nullEffect is not nested in comparison model.')
  }else{
    if(!any(c('saturated', 'none') %in% comparison) && !all(comparison %in% nullEffect)) stop('Comparison model is not nested in hypothesized model; all restrictions in comparison must also be present in nullEffect.')
  }
  
  ### generate sigmas
  # we use variance scaling. If using reference indicators instead, 
  # the model string requires adaptation (in the likely case of unequal loadings), 
  # because currently only the generated modelstring from the first group (generated[[1]]) 
  # is used to define modelH0 and modelH1
  # the downside is that lv.variances is always true.
  generated <- semPower.genSigma(...)   
  
  # more input validations
  if(!is.list(generated[[1]])) stop('Loadings, Phi, Beta, etc. must be provided as a list for each group.')
  if(is.null(generated[[1]][['mu']])){
    inv <- FALSE
    if(useLavOptions){
      inv <- any(grepl('intercepts|means', comparison)) || any(grepl('intercepts|means', nullEffect))
    }else{
      inv <- (nullEffect == 'scalar' || nullEffect == 'residuals')
    }
    if(inv) stop('The models imply a meanstructure, so tau and/or Alpha need to be defined.')
  }

  # models are the same, the only difference pertains to lavOptions
  modelH0 <- modelH1 <- generated[[1]][['modelTrueCFA']]
  
  # set proper lavOptions
  lavOptionsH1 <- NULL
  if(!useLavOptions){
    lavOptionsH0 <- list(group.equal = switch(nullEffect,
                                              'metric' = c('loadings'),
                                              'scalar' = c('loadings', 'intercepts'),
                                              'residual' = c('loadings', 'intercepts', 'residuals')
    ))
    if(comparison %in% c('metric', 'scalar', 'residuals')){
      lavOptionsH1 <- list(group.equal = switch(comparison,
                                                'metric' = c('loadings'),
                                                'scalar' = c('loadings', 'intercepts'),
                                                'residual' = c('loadings', 'intercepts', 'residuals')
      ))
    }
  }else{
    lavOptionsH0 <- list(group.equal = nullEffect)
    if(!any(c('saturated', 'none') %in% comparison)){
      lavOptionsH1 <- list(group.equal = comparison)
    }
  }
  
  if('saturated' %in% comparison) modelH1 <- NULL
  
  Sigma <- lapply(generated, '[[', 'Sigma')
  mu <- NULL
  if(!useLavOptions){
    if(nullEffect == 'scalar' || nullEffect == 'residuals')
      mu <- lapply(generated, '[[', 'mu')
  }else{
    if(any(grepl('intercepts|means', comparison)) || any(grepl('intercepts|means', nullEffect)))
      mu <- lapply(generated, '[[', 'mu')
  }
  
  semPower.powerLav(type = type,
                    Sigma = Sigma,
                    mu = mu,
                    modelH0 = modelH0,
                    modelH1 = modelH1,
                    fitH1model = TRUE,
                    lavOptions = lavOptionsH0,
                    lavOptionsH1 = lavOptionsH1,
                    ...)
  
}

