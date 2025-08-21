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
#' @param fittingFunction one of `'ML'` (default), `'WLS'`, `'DWLS'`, `'ULS'`. Defines the fitting function used to obtain SigmaHat in analytical power analyses. This also implies a certain discrepancy function used to obtain Fmin.
#' @param simulatedPower whether to perform a simulated (`TRUE`, rather than analytical, `FALSE`) power analysis. See [simulate()] for additional options.
#' @param lavOptions a list of additional options passed to `lavaan`, e. g., `list(estimator = 'mlm')` to request robust ML estimation. Mostly useful in conjunction with `simulatedPower`. 
#' @param lavOptionsH1 alternative options passed to `lavaan` that are only used for the H1 model. If `NULL`, identical to `lavOptions`. Probably only useful for multigroup models.
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()]. See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details
#' Generic function to perform a power analysis based  on a true population covariance matrix Sigma 
#' and a model implied covariance matrix SigmaHat (and optionally the associated mean vectors), 
#' where SigmaHat (and muHat) is determined by fitting a respective H0 model using `lavaan`, 
#' and Sigma (and mu) can also be provided through a corresponding `lavaan` model string.
#' 
#' All `semPower` convenience functions internally call this function.
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the requested type of **power analysis**:
#' * `alpha`: The alpha error probability. Required for `type = 'a-priori'` and `type = 'post-hoc'`.
#' * Either `beta` or `power`: The beta error probability and the statistical power (1 - beta), respectively. Only for `type = 'a-priori'`.
#' * `N`: The sample size. Always required for `type = 'post-hoc'` and `type = 'compromise'`. For `type = 'a-priori'` and multiple group analysis, `N` is a list of group weights.
#' * `abratio`: The ratio of alpha to beta. Only for `type = 'compromise'`. 
#' 
#' If a **simulated power analysis** (`simulatedPower = TRUE`) is requested, optional arguments can be provided as a list to `simOptions`:
#' * `nReplications`: The targeted number of simulation runs. Defaults to 250, but larger numbers greatly improve accuracy at the expense of increased computation time.
#' * `minConvergenceRate`:  The minimum convergence rate required, defaults to .5. The maximum actual simulation runs are increased by a factor of 1/minConvergenceRate.
#' * `type`: specifies whether the data should be generated from a population assuming multivariate normality (`'normal'`; the default), or based on an approach generating non-normal data (`'IG'`, `'mnonr'`, `'RC'`, or `'VM'`). 
#' The approaches generating non-normal data require additional arguments detailed below.
#' * `missingVars`: vector specifying the variables containing missing data (defaults to NULL).
#' * `missingVarProp`: can be used instead of `missingVars`: The proportion of variables containing missing data (defaults to zero).
#' * `missingProp`: The proportion of missingness for variables containing missing data (defaults to zero), either a single value or a vector giving the probabilities for each variable.
#' * `missingMechanism`: The missing data mechanism, one of `MCAR` (the default), `MAR`, or `NMAR`.
#' * `nCores`: The number of cores to use for parallel processing. Defaults to 1 (= no parallel processing). This requires the `doSNOW` package.
#' 
#' `type = 'IG'` implements the independent generator approach (IG, Foldnes & Olsson, 2016) approach 
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors. This requires the `covsim` package.
#' 
#' `type = 'mnonr'` implements the approach suggested by Qu, Liu, & Zhang (2020) and requires provision of  Mardia's multivariate skewness (`skewness`)  and kurtosis  (`kurtosis`), where 
#' skewness must be non-negative and kurtosis must be at least 1.641 skewness + p (p + 0.774), where p is the number of variables. This requires the `mnonr` package.
#' 
#' `type = 'RK'` implements the approach suggested by Ruscio & Kaczetow (2008) and requires provision of the population distributions
#'  of each variable (`distributions`). `distributions` must be a list (if all variables shall be based on the same population distribution) or a list of lists. 
#'  Each component must specify the population distribution (e.g. `rchisq`) and additional arguments (`list(df = 2)`).
#' 
#' `type = 'VM'` implements the third-order polynomial method (Vale & Maurelli, 1983) 
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.
#' 
#' @examples
#' \dontrun{
#' # set up two CFA factors with a true correlation of .2
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
#' # define the H0 analysis model (restricting the factor correlation to zero) 
#' mH0 <- '
#'   f1 =~ x1 + x2 + x3
#'   f2 =~ x4 + x5 + x6
#'   f1 ~~ 0*f2
#' '
#' # determine N to reject the H0 that the correlation is zero 
#' # with a power of 95% on alpha = .05
#' powerLav <- semPower.powerLav(type = 'a-priori', 
#'                               modelPop = mPop, modelH0 = mH0,
#'                               alpha = .05, beta = .05)
#' summary(powerLav)
#' 
#' # same as above, but also define an H1 comparison model 
#' mH1 <- '
#'   f1 =~ x1 + x2 + x3
#'   f2 =~ x4 + x5 + x6
#'   f1 ~~ f2
#' '
#' powerLav <- semPower.powerLav(type = 'a-priori', 
#'                               modelPop = mPop, modelH0 = mH0, modelH1 = mH1,
#'                               alpha = .05, beta = .05)
#' 
#' # same as above, but use covariance matrix input instead of modelPop
#' gen <- semPower.genSigma(Phi = .2, 
#'                          loadings = list(c(.5, .6, .4), c(.7, .8, .3)))
#' Sigma <- gen$Sigma
#' powerLav <- semPower.powerLav(type = 'a-priori', 
#'                               Sigma = Sigma, modelH0 = mH0,
#'                               alpha = .05, beta = .05)
#'
#' # note all of the above is identical to the output provided by the semPower.powerCFA function
#' powerCFA <- semPower.powerCFA(type = 'a-priori',
#'                               comparison = 'saturated',
#'                               Phi = .2, 
#'                               loadings = list(c(.5, .6, .4), c(.7, .8, .3)), 
#'                               alpha = .05, beta = .05)
#' 
#' # same as above, but perform simulated power analysis
#' # with 250 replications using a robust ML test-statistic
#' set.seed(300121)
#' powerLav <- semPower.powerLav(type = 'a-priori', 
#'                               Sigma = Sigma, modelH0 = mH0,
#'                               alpha = .05, beta = .05, 
#'                               simulatedPower = TRUE,
#'                               simOptions = list(nReplications = 250)
#'                               lavOptions = list(estimator = 'MLM'))
#' }
#' @seealso [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @importFrom utils installed.packages
#' @export
semPower.powerLav <- function(type, 
                              modelPop = NULL, 
                              modelH0 = NULL, modelH1 = NULL, fitH1model = TRUE, 
                              Sigma = NULL, mu = NULL,
                              fittingFunction = 'ML', 
                              simulatedPower = FALSE, 
                              lavOptions = NULL, lavOptionsH1 = lavOptions, 
                              ...){
  
  # check whether lavaan is available
  if(!'lavaan' %in% rownames(installed.packages())) stop('This function depends on the lavaan package, so install lavaan first.')
  # validate input
  type <- checkPowerTypes(type)
  knownFittingFunctions <- c('ML', 'WLS', 'DWLS', 'ULS')
  fittingFunction <- toupper(fittingFunction)
  if(!fittingFunction %in% knownFittingFunctions){
    stop(paste("Fitting function must be one of", paste(knownFittingFunctions, collapse = ', '), '.'))
  }
  
  ### actually, estimators such as DWLS and ULS use a reduced weight matrix for parameter estimation, 
  ### but the full weight matrix for model tests. This does not match what is reported by lavaan, however, 
  ### so for now we use DWLS/ULS for both, estimating SigmaHat and computing fmin. 
  #discrepancyFunction <- getDiscrepancyFunctionFromFittingFunction(fittingFunction)
  discrepancyFunction <- fittingFunction
  
  if(is.null(modelH0)) stop('Provide a lavaan model string defining the analysis (H0) model.')
  if(is.null(modelPop) && is.null(Sigma)) stop('Either provide a lavaan model string defining the population model or provide the population covariance matrix Sigma.')
  if(!is.null(modelPop) && !is.null(Sigma)) stop('Either provide a lavaan model string defining the population model or provide the population covariance matrix Sigma, but not both.')
  if(simulatedPower && type == 'compromise') stop('Simulated power is not available for compromise power analysis, because this would require a vast (infeasible) number of simulation runs to yield reliable results.')
  if(!is.null(modelPop) && !is.list(modelPop)) modelPop <- list(modelPop)
  if(!is.null(Sigma) && !is.list(Sigma)) Sigma <- list(Sigma)
  if(!is.null(mu) && !is.list(mu)) mu <- list(mu)
  

  # determine population Sigma / mu
  if(is.null(Sigma)){
    Sigma <- lapply(modelPop, function(x) orderLavCov(lavaan::fitted(lavaan::sem(x))[['cov']]))
    mu <- lapply(modelPop, function(x) orderLavMu(lavaan::fitted(lavaan::sem(x))[['mean']]))
  }
  # make sure sigma/mu are properly ordered
  Sigma <- lapply(Sigma, function(x) orderLavCov(x))
  if(!is.null(mu)) mu <- lapply(mu, function(x) orderLavMu(x))
  
  nGroups <- length(Sigma)
  

  # analytical power  (also called once when simulated power is requested to obtain analytic comparison)
  if(!simulatedPower){
    
    # we need to call lavaan() directly with defaults as defined in sem()
    lavOptions <- getLavOptions(lavOptions, nGroups = nGroups)
    # make sure proper estimator is set
    if(!is.null(lavOptions[['estimator']]) && 
       !toupper(lavOptions[['estimator']]) %in% knownFittingFunctions) stop('Analytical power is only available using ML, WLS, DWLS, or ULS estimation. Note that most robust derivatives (mlm etc) are asymptotically identical.')
    lavOptions[['estimator']] <- fittingFunction
    if(fittingFunction %in% c('WLS','DWLS','ULS')) lavOptions[['likelihood']] <- NULL
    if(fittingFunction == 'WLS'){
      lavOptions[['WLS.V']] <- lapply(seq(nGroups), function(x) getWLSv(Sigma[[x]], mu[[x]]))
      if(nGroups == 1) lavOptions[['WLS.V']] <- lavOptions[['WLS.V']][[1]]
    } 
    if(fittingFunction == 'DWLS'){
      lavOptions[['WLS.V']] <- lapply(seq(nGroups), function(x) getWLSv(Sigma[[x]], mu[[x]], diag = TRUE))
      if(nGroups == 1) lavOptions[['WLS.V']] <- lavOptions[['WLS.V']][[1]]
    } 

    # get H0 sigmaHat / muHat
    modelH0Fit <- suppressWarnings(  # suppress here, throw warnings in second try
      do.call(lavaan::lavaan,
                          append(list(model = modelH0,
                                      sample.cov = if(nGroups > 1) Sigma else Sigma[[1]],  # lav complains when providing a list of length 1
                                      sample.mean = if(nGroups > 1) mu else mu[[1]]),
                                 lavOptions))
      )
    if(!modelH0Fit@optim[['converged']]){
      # try again using starting values from previous run
      modelH0Fit <- do.call(lavaan::lavaan,
                            append(list(model = modelH0,
                                        sample.cov = if(nGroups > 1) Sigma else Sigma[[1]],
                                        sample.mean = if(nGroups > 1) mu else mu[[1]], 
                                        start = modelH0Fit),
                                   lavOptions))
      if(!modelH0Fit@optim[['converged']]) stop('The H0 model did not converge.')
    } 
    if(nGroups > 1){
      # multigroup case
      SigmaHat <- lapply(seq(nGroups), function(x) orderLavCov(lavaan::fitted(modelH0Fit)[[x]][['cov']]))
      muHat <- lapply(seq(nGroups), function(x) orderLavMu(lavaan::fitted(modelH0Fit)[[x]][['mean']]))
    }else{
      # single group case (lav doesnt return a list...)
      SigmaHat <- list(orderLavCov(lavaan::fitted(modelH0Fit)[['cov']]))
      muHat <- list(orderLavMu(lavaan::fitted(modelH0Fit)[['mean']]))
    }
    df <- dfH0 <- modelH0Fit@test[['standard']][['df']]  # this is probably invalid for estm with adjusted df
    
    # get H1 comparison model and deltaF
    if(!is.null(modelH1) && fitH1model){
      
      # we need to call lavaan() directly with defaults as defined in sem()
      lavOptionsH1 <- getLavOptions(lavOptionsH1, nGroups = nGroups)
      # make sure proper estimator is set
      lavOptionsH1[['estimator']] <- fittingFunction
      if(fittingFunction %in% c('WLS','DWLS','ULS')) lavOptionsH1[['likelihood']] <- NULL
      if(fittingFunction == 'WLS'){
        lavOptionsH1[['WLS.V']] <- lapply(seq(nGroups), function(x) getWLSv(Sigma[[x]], mu[[x]]))
        if(nGroups == 1) lavOptionsH1[['WLS.V']] <- lavOptionsH1[['WLS.V']][[1]]
      } 
      if(fittingFunction == 'DWLS'){
        lavOptionsH1[['WLS.V']] <- lapply(seq(nGroups), function(x) getWLSv(Sigma[[x]], mu[[x]], diag = TRUE))
        if(nGroups == 1) lavOptionsH1[['WLS.V']] <- lavOptionsH1[['WLS.V']][[1]]
      } 

      modelH1Fit <- suppressWarnings(  # suppress here, throw warnings in second try
        do.call(lavaan::lavaan,
                append(list(model = modelH1,
                            sample.cov = if(length(Sigma) > 1) Sigma else Sigma[[1]],
                            sample.mean = if(length(Sigma) > 1) mu else mu[[1]]),
                       lavOptionsH1))        
      )
      if(!modelH1Fit@optim[['converged']]){
        # try again using starting values from previous run
        modelH1Fit <- do.call(lavaan::lavaan,
                              append(list(model = modelH1,
                                          sample.cov = if(nGroups > 1) Sigma else Sigma[[1]],
                                          sample.mean = if(nGroups > 1) mu else mu[[1]], 
                                          start = modelH1Fit),
                                     lavOptionsH1))
        
        if(!modelH1Fit@optim[['converged']]) stop('The H1 model did not converge.')
      }
      dfH1 <- modelH1Fit@test[['standard']][['df']]
      if(dfH1 >= dfH0) stop('The df of the H0 model are not larger than the df of the H1 model, as they should be.')
      # get delta F
      if(nGroups > 1){
        # multigroup case
        fminH0 <- lapply(seq(nGroups), 
                         function(x) getF.Sigma(
                           orderLavCov(lavaan::fitted(modelH0Fit)[[x]][['cov']]), Sigma[[x]],
                           orderLavMu(lavaan::fitted(modelH0Fit)[[x]][['mean']]), mu[[x]],
                           fittingFunction = discrepancyFunction))
        fminH1 <- lapply(seq(nGroups), 
                         function(x) getF.Sigma(
                           orderLavCov(lavaan::fitted(modelH1Fit)[[x]][['cov']]), Sigma[[x]], 
                           orderLavMu(lavaan::fitted(modelH1Fit)[[x]][['mean']]), mu[[x]],
                           fittingFunction = discrepancyFunction))
        deltaF <- lapply(seq(nGroups), function(x) fminH0[[x]] - fminH1[[x]]) # result must be a list
      }else{
        # single group case
        fminH0 <- getF.Sigma(orderLavCov(lavaan::fitted(modelH0Fit)[['cov']]), Sigma[[1]], 
                             orderLavMu(lavaan::fitted(modelH0Fit)[['mean']]), mu[[1]],
                             fittingFunction = discrepancyFunction)
        fminH1 <- getF.Sigma(orderLavCov(lavaan::fitted(modelH1Fit)[['cov']]), Sigma[[1]], 
                             orderLavMu(lavaan::fitted(modelH1Fit)[['mean']]), mu[[1]],
                             fittingFunction = discrepancyFunction)
        deltaF <- fminH0 - fminH1
      }
      if(sum(unlist(deltaF)) < 1e-10) warning(paste0('The H0 model shows the same discrepancy as the H1 model. This usually happens when the H0 model contains restrictions that are valid in the population. Check the definition of the population values and the H0 constraints.'))
      if(any(fminH1 > 1e-6)) warning(paste0('H1 model yields imperfect fit (F0 = ', round(unlist(fminH1)[which(unlist(fminH1) > 1e-6)[1]], 6), '). This may happen if the H1 model contains restrictions on parameters (such as invariance constraints) that actually differ in the population. Verify that this is intended.'))
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
                        fittingFunction = discrepancyFunction,
                        ...)    
    }else{
      power <- semPower(type = type, 
                        effect = deltaF, effect.measure = "F0", 
                        df = df, 
                        fittingFunction = discrepancyFunction, # actually irrelevant
                        ...)    
    }    
        
  # simulated power
  }else{
    power <- semPower(type = type, 
                      Sigma = Sigma, mu = mu, 
                      modelH0 = modelH0, modelH1 = modelH1, fitH1model = fitH1model,
                      fittingFunction = fittingFunction,
                      simulatedPower = simulatedPower, 
                      # simulate() takes care of proper lavOptions
                      lavOptions = lavOptions,  
                      lavOptionsH1 = lavOptionsH1, 
                      ...)
    SigmaHat <- muHat <- NULL
  }

  # remove list structure for single group models
  if(nGroups == 1){
    Sigma <- Sigma[[1]]
    SigmaHat <- SigmaHat[[1]]
    mu <- mu[[1]]
    muHat <- muHat[[1]]
  }

  result <- append(power, 
                list(SigmaHat = SigmaHat, Sigma = Sigma,
                     muHat = muHat, mu = mu,
                     modelPop = modelPop, modelH0 = modelH0, modelH1 = modelH1))
  
  switch(type, 
         'a-priori'= {
           class(result) <- "semPower.aPriori"
         },
         'post-hoc' = {
           class(result) <- "semPower.postHoc"
         },
         'compromise' = {
           class(result) <- "semPower.compromise"
         }
  )

  result
}
