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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.  This requires the `semTools` package.
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
  if(!is.null(mu) && !is.list(mu)) mu <- list(mu)
  
  # the following is probably no longer needed, as we now always supply lav-friendly restrictions
  
  # # lav doesn't like both equality constrains and value constrains on the same parameters, so
  # # transform this by dropping equality constrains and assign value constrains to the affected parameters
  # modelH0 <- makeRestrictionsLavFriendly(modelH0)
  # if(!is.null(modelH1)) modelH1 <- makeRestrictionsLavFriendly(modelH1)

  # determine population Sigma / mu
  if(is.null(Sigma)){
    Sigma <- lapply(modelPop, function(x) orderLavCov(lavaan::fitted(lavaan::sem(x))[['cov']]))
    mu <- lapply(modelPop, function(x) orderLavMu(lavaan::fitted(lavaan::sem(x))[['mean']]))
  }
  
  # analytical power
  if(!simulatedPower){
    
    # we need to call lavaan() directly with defaults as defined in sem()
    lavOptions <- getLavOptions(lavOptions, nGroups = length(Sigma))
    if(!is.null(lavOptions[['estimator']]) && toupper(lavOptions[['estimator']]) != "ML") stop('Analytical power is only available with ML estimation. Note that power based on ML derivatives (mlm etc) is asymptotically identical.')

    # get H0 sigmaHat / muHat
    modelH0Fit <- suppressWarnings(  # suppress here, throw warnings in second try
      do.call(lavaan::lavaan,
                          append(list(model = modelH0,
                                      sample.cov = if(length(Sigma) > 1) Sigma else Sigma[[1]],
                                      sample.mean = if(length(Sigma) > 1) mu else mu[[1]]),
                                 lavOptions))
      )
    if(!modelH0Fit@optim[['converged']]){
      # try again using starting values from previous run
      modelH0Fit <- do.call(lavaan::lavaan,
                            append(list(model = modelH0,
                                        sample.cov = if(length(Sigma) > 1) Sigma else Sigma[[1]],
                                        sample.mean = if(length(Sigma) > 1) mu else mu[[1]], 
                                        start = modelH0Fit),
                                   lavOptions))
      if(!modelH0Fit@optim[['converged']]) stop('The H0 model did not converge.')
    } 
    if(length(Sigma) > 1){
      # multigroup case
      SigmaHat <- lapply(1:length(Sigma), function(x) orderLavCov(lavaan::fitted(modelH0Fit)[[x]][['cov']]))
      muHat <- lapply(1:length(Sigma), function(x) orderLavMu(lavaan::fitted(modelH0Fit)[[x]][['mean']]))
    }else{
      # single group case
      SigmaHat <- list(orderLavCov(lavaan::fitted(modelH0Fit)[['cov']]))
      muHat <- list(orderLavMu(lavaan::fitted(modelH0Fit)[['mean']]))
    }
    df <- dfH0 <- modelH0Fit@test[['standard']][['df']]  # this is probably invalid for estm with adjusted df
    
    # get H1 comparison model and deltaF
    if(!is.null(modelH1) && fitH1model){
      lavOptionsH1 <- getLavOptions(lavOptionsH1, nGroups = length(Sigma))
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
                                          sample.cov = if(length(Sigma) > 1) Sigma else Sigma[[1]],
                                          sample.mean = if(length(Sigma) > 1) mu else mu[[1]], 
                                          start = modelH1Fit),
                                     lavOptionsH1))
        
        if(!modelH1Fit@optim[['converged']]) stop('The H1 model did not converge.')
      }
      dfH1 <- modelH1Fit@test[['standard']][['df']]
      if(dfH1 >= dfH0) stop('The df of the H0 model are not larger than the df of the H1 model, as they should be.')
      # get delta F
      if(length(Sigma) > 1){
        # multigroup case
        fminH0 <- lapply(1:length(Sigma), 
                         function(x) getF.Sigma(orderLavCov(lavaan::fitted(modelH0Fit)[[x]][['cov']]), Sigma[[x]], 
                                                orderLavMu(lavaan::fitted(modelH0Fit)[[x]][['mean']]), mu[[x]]))
        fminH1 <- lapply(1:length(Sigma), 
                         function(x) getF.Sigma(orderLavCov(lavaan::fitted(modelH1Fit)[[x]][['cov']]), Sigma[[x]], 
                                                orderLavMu(lavaan::fitted(modelH1Fit)[[x]][['mean']]), mu[[x]]))
        deltaF <- lapply(1:length(Sigma), function(x) fminH0[[x]] - fminH1[[x]]) # result must be a list
      }else{
        # single group case
        fminH0 <- getF.Sigma(orderLavCov(lavaan::fitted(modelH0Fit)[['cov']]), Sigma[[1]], 
                             orderLavMu(lavaan::fitted(modelH0Fit)[['mean']]), mu[[1]])
        fminH1 <- getF.Sigma(orderLavCov(lavaan::fitted(modelH1Fit)[['cov']]), Sigma[[1]], 
                             orderLavMu(lavaan::fitted(modelH1Fit)[['mean']]), mu[[1]])
        deltaF <- fminH0 - fminH1
      }
      if(sum(unlist(deltaF)) < 1e-10) warning(paste0('The H0 model fits shows the same discrepancy as the H1 model. This usually happens when the H0 model contains restrictions that are valid in the population. Check the definition of the population values and the H0 constraints.'))
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



#' semPower.powerCFA
#'
#' Convenience function for performing power analyses for CFA models to reject one of the following hypotheses: 
#' (a) a zero correlation between two factors, (b) the equality of two correlations between factors,
#' or (c) the equality of a correlation between two factors across two or more groups. 
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param Phi either a single number defining the correlation between exactly two factors or the factor correlation matrix. A list for multiple group models.
#' @param nullEffect defines the hypothesis of interest, must be one of `'cor = 0'` (the default) to test whether a correlation is zero, `'corX = corZ'` to test for the equality of correlations, `'corA = corB'` to test for the equality of a correlation across groups, and `loading = 0` to test whether a loading is zero. Define the correlations to be set to equality in `nullWhich` and the groups in `nullWhichGroups`. 
#' @param nullWhich vector of size 2 indicating which element in Lambda should equal zero when `nullEffect = 'loading = 0'`, or which factor correlation in `Phi` is hypothesized to equal zero when `nullEffect = 'cor = 0'`, or to restrict to equality across groups when `nullEffect = 'corA = corB'`, or list of vectors defining which correlations to restrict to equality when `nullEffect = 'corX = corZ'`. Can also contain more than two correlations, e.g., `list(c(1, 2), c(1, 3), c(2, 3))` to set `Phi[1, 2] = Phi[1, 3] = Phi[2, 3]`. If omitted, the correlation between the first and the second factor is targeted, i. e., `nullWhich = c(1, 2)`.
#' @param nullWhichGroups for `nullEffect = 'corA = corB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details 
#' 
#' This function performs a power analysis to reject various hypotheses arising
#' in standard CFA models:
#' * `nullEffect = 'loading = 0'`: Tests the hypothesis that a loading is zero.
#' * `nullEffect = 'cor = 0'`: Tests the hypothesis that the correlation between two factors is zero. 
#' * `nullEffect = 'corX = corZ'`: Tests the hypothesis that two or more correlations between three or more factors are equal to each other.
#' * `nullEffect = 'corA = corB'`: Tests the hypothesis that the correlation between two factors is equal in two or more groups (always assuming metric invariance).
#' 
#' For hypotheses regarding regression relationships between factors, see [semPower.powerRegression()].
#' For hypotheses regarding mediation effects, see [semPower.powerMediation()].
#' For hypotheses regarding measurement invariance, see [semPower.powerMI()].
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model**:
#' * `Lambda`: The factor loading matrix (with the number of columns equaling the number of factors).
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading on the second factor by .8, .6, .6, and .4.
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the loadings of the first factor to equal .5 and those of the second factor do equal .6.
#' 
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` need to be defined.
#' If the model contains observed variables only, use `Lambda = diag(x)` where `x` is the number of variables.
#' 
#' Additional arguments related to the requested type of **power analysis**:
#' * `alpha`: The alpha error probability. Required for `type = 'a-priori'` and `type = 'post-hoc'`.
#' * Either `beta` or `power`: The beta error probability and the statistical power (1 - beta), respectively. Only for `type = 'a-priori'`.
#' * `N`: The sample size. Always required for `type = 'post-hoc'` and `type = 'compromise'`. For `type = 'a-priori'` and multiple group analysis, `N` is a list of group weights.
#' * `abratio`: The ratio of alpha to beta. Only for `type = 'compromise'`. 
#' 
#' If a **simulated power analysis** (`simulatedPower = TRUE`) is requested, optional arguments can be provided as a list to `simOptions`:
#' * `nReplications`: The targeted number of simulation runs. Defaults to 250, but larger numbers greatly improve accuracy at the expense of increased computation time.
#' * `minConvergenceRate`:  The minimum convergence rate required, defaults to .5. The maximum actual simulation runs are increased by a factor of 1/minConvergenceRate.
#' * `type`: specifies whether the data should be generated from a population assuming multivariate normality (`'normal'`; the default), or based on an approach generating non-normal data (`'IG'`, `'mnonr'`, `'RC'`, or `'VM'`). 
#' * `missingVars`: vector specifying the variables containing missing data (defaults to NULL).
#' * `missingVarProp`: can be used instead of `missingVars`: The proportion of variables containing missing data (defaults to zero).
#' * `missingProp`: The proportion of missingness for variables containing missing data (defaults to zero), either a single value or a vector giving the probabilities for each variable.
#' * `missingMechanism`: The missing data mechanism, one of `MCAR` (the default), `MAR`, or `NMAR`.
#' The approaches generating non-normal data require additional arguments detailed below.
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.  This requires the `semTools` package.
#' 
#' @examples
#' \dontrun{
#' # get required N to detect a correlation of >= .2 between two factors
#' # with a power of 95% on alpha = 5%, where the factors are  
#' # measured by 5 and 6 indicators, respectively, and all loadings are equal to .5
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               Phi = .2, 
#'                               nIndicator = c(5, 6), loadM = .5,
#'                               alpha = .05, beta = .05)
#' # show summary
#' summary(powercfa)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powercfa$modelH1, sample.cov = powercfa$Sigma, 
#' sample.nobs = powercfa$requiredN, sample.cov.rescale = FALSE)
#' lavaan::sem(powercfa$modelH0, sample.cov = powercfa$Sigma, 
#' sample.nobs = powercfa$requiredN, sample.cov.rescale = FALSE)
#' 
#' # same as above, but determine power with N = 500 on alpha = .05
#' powercfa <- semPower.powerCFA(type = 'post-hoc',
#'                               Phi = .2, 
#'                               nIndicator = c(5, 6), loadM = .5,
#'                               alpha = .05, N = 500)
#' 
#' # same as above, but determine the critical chi-square with N = 500 so that alpha = beta
#' powercfa <- semPower.powerCFA(type = 'compromise',
#'                               Phi = .2, 
#'                               nIndicator = c(5, 6), loadM = .5,
#'                               abratio = 1, N = 500)
#'                               
#' # same as above, but compare to the saturated model 
#' # (rather than to the less restricted model)
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               comparison = 'saturated',
#'                               Phi = .2, 
#'                               nIndicator = c(5, 6), loadM = .5,
#'                               alpha = .05, beta = .05)
#'                               
#' # same as above, but provide a reduced loading matrix defining
#' # three indicators with loadings of .7, .6, and .5 on the first factor and
#' # four indicators with loadings of .5, .6, .4, .8 on the second factor 
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               Phi = .2, 
#'                               loadings = list(c(.7, .6, .5), 
#'                                               c(.5, .6, .4, .8)),
#'                               alpha = .05, beta = .05)
#'
#' # detect that the loading of indicator 4 on the first factor differs from zero
#' Lambda <- matrix(c(
#'   c(.8, 0),
#'   c(.4, 0),
#'   c(.6, 0),
#'   c(.1, .5),
#'   c(0, .6),
#'   c(0, .7)
#' ), ncol = 2, byrow = TRUE)
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               Phi = .2,
#'                               nullEffect = 'loading = 0', 
#'                               nullWhich = c(4, 1), 
#'                               Lambda = Lambda,
#'                               alpha = .05, beta = .05)
#'
#'
#' # get required N to detect a correlation of >= .3 between factors 1 and 3  
#' # in a three factor model. Factors are measured by 3 indicators each, and all loadings 
#' # on the first, second, and third factor are .5, .6, and .7, respectively.
#' Phi <- matrix(c(
#'   c(1.00, 0.20, 0.30),
#'   c(0.20, 1.00, 0.10),
#'   c(0.30, 0.10, 1.00)
#' ), ncol = 3,byrow = TRUE)
#' 
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               Phi = Phi,
#'                               nullWhich = c(1, 3), 
#'                               nIndicator = c(3, 3, 3), loadM = c(.5, .6, .7),
#'                               alpha = .05, beta = .05)
#' 
#' # same as above, but ask for N to detect that 
#' # the correlation between factors 1 and 2 (of r = .2) differs from
#' # the correlation between factors 2 and 3 (of r = .3).
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               Phi = Phi,
#'                               nullEffect = 'corX = corZ',
#'                               nullWhich = list(c(1, 2), c(1, 3)), 
#'                               nIndicator = c(3, 3, 3), loadM = c(.5, .6, .7),
#'                               alpha = .05, beta = .05)
#'                               
#' # same as above, but ask for N to detect that all three correlations are unequal
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               Phi = Phi,
#'                               nullEffect = 'corX = corZ',
#'                               nullWhich = list(c(1, 2), c(1, 3), c(2, 3)), 
#'                               nIndicator = c(3, 3, 3), loadM = c(.5, .6, .7),
#'                               alpha = .05, beta = .05)
#'                               
#' # get required N to detect that the correlation between two factors
#' # in group 1 (of r = .2) differs from the one in group 2 (of r = .4). 
#' # The measurement model is identical for both groups:
#' # The first factor is measured by 3 indicators loading by .7 each, 
#' # the second factor is measured by 6 indicators loading by .5 each.
#' # Both groups are sized equally (N = list(1, 1)).
#' powercfa <- semPower.powerCFA(type = 'a-priori', 
#'                               nullEffect = 'corA = corB',
#'                               Phi = list(.2, .4), 
#'                               loadM = c(.7, .5), 
#'                               nIndicator = c(3, 6), 
#'                               alpha = .05, beta = .05, N = list(1, 1))
#'
#' # request a simulated post-hoc power analysis with 500 replications.
#' set.seed(300121)
#' powercfa <- semPower.powerCFA(type = 'post-hoc',
#'                               Phi = .2, 
#'                               nIndicator = c(5, 6), loadM = .5,
#'                               alpha = .05, N = 500, 
#'                               simulatedPower = TRUE, 
#'                               simOptions = list(nReplications = 500))
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
  args <- list(...)
  
  # validate input
  checkEllipsis(...)
  comparison <- checkComparisonModel(comparison)
  nullEffect <- checkNullEffect(nullEffect, c('cor=0', 'corx=corz', 'cora=corb', 'loading=0'))
  if(nullEffect == 'loading=0' && is.null(Phi)) Phi <- diag(1)
  if(is.null(Phi)) stop('Phi must be defined')
  if(!is.null(nullWhichGroups) && !is.list(Phi)) stop('Phi must be provided for each group.')
  if(nullEffect == 'cora=corb' && !is.list(Phi)) stop('corA=corB refers to multigroup analysis, so Phi must be a list.')
  if(nullEffect == 'corx=corz' && (!is.matrix(Phi) || ncol(Phi) <= 2)) stop('corx=corz compares two correlations and thus requires at least three factors and a correlation matrix for Phi.')
  if(is.list(Phi) && !is.null(nullWhichGroups)) lapply(as.list(nullWhichGroups), function(x) checkBounded(x, 'Each element in nullWhichGroups', bound = c(1, length(Phi)), inclusive = TRUE)) 
    
  # generate sigma 
  generated <- semPower.genSigma(Phi = Phi, ...)
  
  ### now do validation of nullWhich, since we now know Phi
  isMultigroup <- is.list(Phi)
  if(isMultigroup) nfac <- ncol(generated[[1]][['Phi']]) else nfac <- ncol(generated[['Phi']])
  if(nullEffect != 'loading=0'){
    if(is.null(nullWhich) && nfac == 2) nullWhich <- c(1, 2)
    if(is.null(nullWhich)) stop('nullWhich must be defined.')
    if(!is.list(nullWhich)) nullWhich <- list(nullWhich)
    if(any(unlist(lapply(nullWhich, function(x) length(x) != 2)))) stop('nullWhich may only contain vectors of size two.')
    if(any(unlist(lapply(nullWhich, function(x) x[1] == x[2])))) stop('elements in nullWhich may not refer to variances.')
    if(any(unlist(lapply(nullWhich, function(x) (x[1] < 1 || x[2] < 1 || x[1] > nfac || x[2] > nfac))))) stop('At least one element in nullWhich is an out of bounds index concerning Phi.')
    if(length(nullWhich) > 1){
      for(i in 1:(length(nullWhich) - 1)){
        for(j in (i + 1):length(nullWhich)){
          if(nullWhich[[i]][1] == nullWhich[[j]][1] && nullWhich[[i]][2] == nullWhich[[j]][2]) stop('elements in nullWhich may not refer to the same correlation')
        }
      }
    }
  }else{
    if(is.null(nullWhich)) stop('nullWhich must be defined.')
    if(length(nullWhich) != 2) stop('nullWhich must be a vector with two elements.')
    if(nullWhich[1] > nrow(generated[['Lambda']]) || nullWhich[2] > ncol(generated[['Lambda']])  ) stop('nullWhich refers to an invalid element in Lambda. The first entry must be <= nIndicator, the second <= nFactors.')
    if(generated[['Lambda']][nullWhich[1], nullWhich[2]] == 0) stop('nullWhich refers to a loading with a population value of zero. The loading referred by nullWhich must differ from zero.')
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
    lab <- paste0('ff', 1:length(Phi))
    lab[nullWhichGroups] <- 'pf1'
    lab <- paste0('c(', paste(lab, collapse = ','), ')*')
    modelH0 <- paste(c(modelH0,
                       paste0('f', nullWhich[[1]], collapse = paste0(' ~~ ', lab))),
                       collapse = '\n')
  }else if(nullEffect == 'loading=0'){
    tInd <- paste0('x', nullWhich[1])
    tFac <- paste0('f', nullWhich[2])
    tok <- strsplit(modelH0, '\n', fixed = TRUE)[[1]]
    tok <- lapply(tok, function(x){
      if(startsWith(x, tFac) && grepl('=~', x)){
        t <- lapply(strsplit(strsplit(x, '=~', fixed = TRUE)[[1]][2], '+', fixed = TRUE)[[1]], trimws)
        idx <- grep(tInd, unlist(t))
        if(grepl('NA', t[[idx]])){
          t[[idx]] <- sub('NA', '0', t[[idx]])
        }else{
          t[[idx]] <- paste0('0*', t[[idx]])
        }
        paste0(tFac, ' =~ ', paste0(unlist(t), collapse = ' + '))
      }else{
        x
      }
    })
    modelH0 <- paste(unlist(tok), collapse = '\n')
  }else{
    stop('nullEffect not defined')
  }

  # we always enforce invariance constraints in the multigroup case
  if(isMultigroup){
    args[['lavOptions']] <- append(args[['lavOptions']], list(group.equal = c('loadings', 'lv.variances')))
  } 

  modelH1 <- NULL
  fitH1model <- FALSE
  if(comparison == 'restricted'){
    if(isMultigroup) modelH1 <- generated[[1]][['modelTrueCFA']] else modelH1 <- generated[['modelTrueCFA']] 
    # single group case: the h1 model always fits perfectly
    # multigroup case: we cannot be sure that user input yields a perfectly fitting model
    fitH1model <- isMultigroup 
  } 
  
  if(isMultigroup) Sigma <- lapply(generated, '[[', 'Sigma') else Sigma <- generated[['Sigma']] 
  
  do.call(semPower.powerLav, append(list(
                    type = type,
                    Sigma = Sigma,
                    modelH0 = modelH0,
                    modelH1 = modelH1,
                    fitH1model = fitH1model),
                    args)
          )

}

#' semPower.powerRegression
#'
#' Convenience function for performing power analysis on slope(s) in a latent regression of the form Y = XB.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param slopes vector of slopes (or a single number for a single slope) of the k predictors for Y. A list of slopes for multigroup models.
#' @param corXX correlation(s) between the k predictors (X). Either `NULL` for uncorrelated predictors, a single number (for k = 2 predictors), or a matrix. Can also be a list for multigroup models providing the correlations by group of matrices (otherwise, the same correlations are used in all groups). 
#' @param nullEffect defines the hypothesis of interest, must be one of `'slope = 0'` (the default) to test whether a slope is zero, `'slopeX = slopeZ'` to test for the equality of slopes, or `'slopeA = slopeB'` to test for the equality of slopes across groups. Define the slopes to set to equality in `nullWhich`.
#' @param nullWhich single number indicating which slope is hypothesized to equal zero when `nullEffect = 'slope = 0'`, or indicating which slope to restrict to equality across groups when `nullEffect = 'slopeA = slopeB'`, or vector defining the slopes to restrict to equality when `nullEffect = 'slopeX = slopeZ'`. Can also contain more than two slopes, e.g. `c(1, 2, 3)` to constrain the first three slopes to equality.
#' @param nullWhichGroups for `nullEffect = 'slopeA = slopeB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param standardized whether all parameters should be standardized (`TRUE`, the default). If `FALSE`, all regression relations are unstandardized.
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. The first factor is treated as Y and the subsequent factors as the predictors X_k. See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details 
#' 
#' This function performs a power analysis to reject various hypotheses arising
#' in SEM models involving a simple regression relation of the form `Y = b_1*X_1 + ... + b_k*X_k` between the factors:
#' * `nullEffect = 'slope = 0'`: Tests the hypothesis that the slope for a predictor is zero. 
#' * `nullEffect = 'slopeX = slopeZ'`: Tests the hypothesis that two or more slopes are equal to each other.
#' * `nullEffect = 'slopeA = slopeB'`: Tests the hypothesis that the slope for a predictor is equal in two or more groups (always assuming metric invariance).
#' 
#' For hypotheses regarding mediation effects, see [semPower.powerMediation()]. For hypothesis in autoregressive models, see  [semPower.powerAutoregressive()].
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model**:
#' * `Lambda`: The factor loading matrix (with the number of columns equaling the number of factors).
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading on the second factor by .8, .6, .6, and .4.
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the loadings of the first factor to equal .5 and those of the second factor do equal .6.
#' 
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` need to be defined.
#' If the model contains observed variables only, use `Lambda = diag(x)` where `x` is the number of variables.
#' 
#' Note that the first factor acts as the criterion Y, the subsequent factors as predictors X_1 to X_k.
#' 
#' Additional arguments related to the requested type of **power analysis**:
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.  This requires the `semTools` package.
#' 
#' @examples
#' \dontrun{
#' # latent regression of the form `Y = .2*X1 + .3*X2`, where X1 and X2 correlate by .4
#' # obtain required N to reject the hypothesis that the slope of X1 is zero 
#' # with a power of 95% on alpha = 5%,   
#' # where Y is measured by 3 indicators loading by .5 each,
#' # X1 by 5 indicators loading by .6 each, and
#' # X2 by 4 indicators loading by .7 each. 
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), corXX = .4, 
#'                                      nullWhich = 1, 
#'                                      nIndicator = c(3, 5, 4), 
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#' # show summary
#' summary(powerReg)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerReg$modelH1, sample.cov = powerReg$Sigma, 
#' sample.nobs = powerReg$requiredN, sample.cov.rescale = FALSE)
#' lavaan::sem(powerReg$modelH0, sample.cov = powerReg$Sigma, 
#' sample.nobs = powerReg$requiredN, sample.cov.rescale = FALSE)
#' 
#' # same as above, but determine power with N = 500 on alpha = .05 
#' powerReg <- semPower.powerRegression(type = 'post-hoc',
#'                                      slopes = c(.2, .3), corXX = .4, 
#'                                      nullWhich = 1, 
#'                                      nIndicator = c(3, 5, 4), 
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, N = 500)
#'                                      
#' # same as above, but determine the critical chi-square with N = 500 so that alpha = beta 
#' powerReg <- semPower.powerRegression(type = 'compromise',
#'                                      slopes = c(.2, .3), corXX = .4, 
#'                                      nullWhich = 1, 
#'                                      nIndicator = c(3, 5, 4), 
#'                                      loadM = c(.5, .6, .7),
#'                                      abratio = .05, N = 500)
#'                                      
#' # same as above, but ask for the required N to detect that the slope of X2 is zero
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), corXX = .4, 
#'                                      nullWhich = 2, 
#'                                      nIndicator = c(3, 5, 4), 
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'
#' # same as above, but define unstandardized slopes
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), corXX = .4,
#'                                      nullWhich = 2, 
#'                                      standardized = FALSE,
#'                                      nIndicator = c(3, 5, 4), 
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'                                      
#' # same as above, but compare to the saturated model
#' # (rather than to the less restricted model)
#' powerReg <- semPower.powerRegression(type = 'a-priori', 
#'                                      comparison = 'saturated',
#'                                      slopes = c(.2, .3), corXX = .4, 
#'                                      nullWhich = 2, 
#'                                      nIndicator = c(3, 5, 4), 
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'                                      
#' # same as above, but provide a reduced loading matrix defining
#' # three indicators with loadings of .7, .6, .5 on the first factor (Y),
#' # four indicators with loadings of .5, .6, .4, .8 on the second factor (X1), and
#' # three indicators with loadings of .8, .7, .8 on the third factor (X2).
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), corXX = .4, 
#'                                      nullWhich = 2, 
#'                                      loadings = list(
#'                                        c(.7, .6, .5), 
#'                                        c(.5, .6, .4, .8),
#'                                        c(.8, .7, .8)),
#'                                      alpha = .05, beta = .05)
#'                               
#' # latent regression of the form `Y = .2*X1 + .3*X2 + .4*X3`, 
#' # providing the predictor intercorrelation matrix,
#' # and ask for the required N to detect that the first slope differs from zero.
#' corXX <- matrix(c(
#'   #   X1    X2    X3
#'   c(1.00, 0.20, 0.30),  # X1
#'   c(0.20, 1.00, 0.10),  # X2
#'   c(0.30, 0.10, 1.00)   # X3
#' ), ncol = 3,byrow = TRUE)
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3, .4), corXX = corXX, 
#'                                      nullWhich = 1,
#'                                      nIndicator = c(4, 3, 5, 4),
#'                                      loadM = c(.5, .5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'
#' # same as above, but ask for the required N to detect that 
#' # the slope for X1 (b = .2) and the slope for X2 (b = .3) differ from each other
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3, .4), corXX = corXX, 
#'                                      nullEffect = 'slopeX = slopeZ', 
#'                                      nullWhich = c(1, 2),
#'                                      nIndicator = c(4, 3, 5, 4),
#'                                      loadM = c(.5, .5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'                                      
#' # same as above, but ask for the required N to reject the hypothesis that 
#' # all three slopes are equal to each other
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3, .4), corXX = corXX, 
#'                                      nullEffect = 'slopeX = slopeZ', 
#'                                      nullWhich = c(1, 2, 3),
#'                                      nIndicator = c(4, 3, 5, 4),
#'                                      loadM = c(.5, .5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'     
#' # get required N to detect that 
#' # the slope for X2 group 1 (of b2 = .3) differs from the slope for X2 in group 2 (of b = .0). 
#' # The remaining slopes are equal in both groups (b1 = .2, b3 = .4).
#' # The measurement model is identical in both groups:
#' # The criterion (Y) is measured by 4 indicators loading by .5 each, 
#' # Predictors X1 and X3 are both measured by 5 indicators loading by .6 each,
#' # Predictor X2 is measured by 3 indicators loading by .7 each.
#' # Both groups are sized equally (N = list(1, 1)).
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = list(c(.2, .3, .4), 
#'                                      c(.2, .0, .4)), 
#'                                      corXX = corXX, 
#'                                      nullEffect = 'slopeA = slopeB', 
#'                                      nullWhich = 2,
#'                                      nIndicator = list(
#'                                         c(4, 5, 3, 5), 
#'                                         c(4, 5, 3, 5)),
#'                                      loadM = list(
#'                                         c(.5, .6, .7, .6),
#'                                         c(.5, .6, .7, .6)), 
#'                                      alpha = .05, beta = .05, 
#'                                      N = list(1, 1))
#'
#'# request a simulated post-hoc power analysis with 500 replications 
#'# to detect that the slope of X1 differs from zero.
#' set.seed(300121)
#' powerReg <- semPower.powerRegression(type = 'post-hoc',
#'                                      slopes = c(.2, .1), 
#'                                      nullWhich = 1,
#'                                      nIndicator = c(4, 3, 3), loadM = .5,
#'                                      alpha = .05, N = 500, 
#'                                      simulatedPower = TRUE, 
#'                                      simOptions = list(nReplications = 500))
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerRegression <- function(type, comparison = 'restricted',
                                     slopes = NULL, 
                                     corXX = NULL, 
                                     nullEffect = 'slope = 0',
                                     nullWhich = NULL,
                                     nullWhichGroups = NULL,
                                     standardized = TRUE,
                                     ...){
  
  args <- list(...)
  comparison <- checkComparisonModel(comparison)
  checkEllipsis(...)
  
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
  if(any(unlist(lapply(corXX, is.vector))) && length(unique(unlist(lapply(corXX, length)))) > 1) stop('corXX must be a single number or a matrix') 
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
    if(length(nullWhich) > 1) stop('nullWhich must be a single number when nullEffect is slope=0 or slopeA=slopeB')
  }
  nullWhich <- nullWhich + 1 # because first factor is criterion
  
  # get temporary Lambda so that we can check whether number of factors matches slopes + 1
  tLambda  <- args[['Lambda']]
  if(is.null(tLambda)){
    if(!isMultigroup){
      tLambda <- genLambda(args[['loadings']], args[['nIndicator']],
                           args[['loadM']], args[['loadSD']], args[['loadMinMax']])
    }else{
      tLambda <- genLambda(args[['loadings']][[1]], args[['nIndicator']][[1]],
                           args[['loadM']][[1]], args[['loadSD']][[1]], args[['loadMinMax']][[1]])
      
    }
  }
  if(ncol(tLambda) != (1 + nrow(slopes[[1]]))) stop('The number of factors does not match the number of slopes + 1. Remember to define a measurement model including both the criterion (Y) and all predictors (X).')
  
  ### calc implied sigma. 
  # standardized case: transform B and Psi to Phi
  if(standardized){
    B <- lapply(slopes, function(x){
      cB <- matrix(0, ncol = (length(x) + 1), nrow = (length(x) + 1))
      cB[nrow(cB), ] <- c(x, 0)
      cB
    })
    Psi <- lapply(seq_along(slopes), function(x){
      cPsi <- diag(ncol(B[[x]]))
      cPsi[1:(nrow(cPsi) - 1), 1:(ncol(cPsi) - 1)] <- corXX[[x]]
      cPsi
    })
    cPhi <- lapply(seq_along(B), function(x) getPhi.B(B[[x]], Psi[[x]]))
    # change order so that Y is first factor
    Phi <- lapply(cPhi, function(x){
      nf <- ncol(x)
      nPhi <- diag(nf)
      nPhi[2:nf, 2:nf] <- x[1:(nf - 1), 1:(nf - 1)]
      nPhi[1, 1:nf] <- c(1, x[nf, 1:(nf-1)])
      nPhi[1:nf, 1] <- c(1, x[1:(nf-1), nf])
      nPhi
    })
    
    generated <- semPower.genSigma(Phi = if(length(Phi) > 1) Phi else Phi[[1]], 
                                   useReferenceIndicator = TRUE, ...)  
    
  # unstandardized case
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
    generated <- semPower.genSigma(Beta = if(length(B) > 1) B else B[[1]], 
                                   Psi = if(length(Psi) > 1) Psi else Psi[[1]],
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
    lab <- paste0('ff', seq_along(slopes))
    lab[nullWhichGroups] <- 'pf1'
    lab <- paste0('c(', paste(lab, collapse = ','), ')*')
    modelH0 <- paste(model, 
                     paste0('f1 ~ ', paste0(lab, 'f', np[(nullWhich - 1)], ' + '), paste0('f',np[-(nullWhich - 1)], collapse = ' + ')),
                     sep = '\n')
  }else{
    stop('nullEffect not defined.')
  }

  # we always enforce metric invariance in the multigroup case
  if(isMultigroup){
    args[['lavOptions']] <- append(args[['lavOptions']], list(group.equal = c('loadings')))
  } 
  
  modelH1 <- NULL
  fitH1model <- FALSE
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
  
  do.call(semPower.powerLav, append(list(
    type = type,
    Sigma = Sigma,
    modelH0 = modelH0,
    modelH1 = modelH1,
    fitH1model = fitH1model),
    args)
  )
}


#' semPower.powerMediation
#'
#' Convenience function for performing power analysis concerning indirect effect(s) in a mediation model.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param bYX the slope (direct effect) for X -> Y. A list for multiple group models. Can be `NULL` if `Beta` is set.
#' @param bMX the slope for X -> M. A list for multiple group models. Can be `NULL` if `Beta` is set.
#' @param bYM the slope for M -> Y. A list for multiple group models. Can be `NULL` if `Beta` is set.
#' @param Beta can be used instead of `bYX`, `bMX`, and `bYM`: matrix of regression weights connecting the latent factors (all-Y notation). Exogenous variables must be in the first row(s), so the upper triangular of Beta must be zero. A list for multiple group models.
#' @param indirect `NULL` unless `Beta` is set. Otherwise a list of vectors of size 2 indicating the elements of `Beta` that define the indirect effect of interest, e.g. `list(c(2, 1), c(3, 2))`. See details.
#' @param nullEffect defines the hypothesis of interest, must be one of `'ind = 0'` (the default) to test whether the indirect effect is zero or `'indA = indB'` to test for the equality of indirect effects across groups. See details.
#' @param nullWhichGroups for `nullEffect = 'indA = indB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param standardized whether all parameters should be standardized (`TRUE`, the default). If `FALSE`, all regression relations are unstandardized.
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. In case of a simple mediation, the order of factors is X, M, Y. See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details
#' This function performs a power analysis to reject various hypotheses arising
#' in the context of mediation:
#' * `nullEffect = 'ind = 0'`: Tests the hypothesis that an indirect effect is zero. 
#' * `nullEffect = 'indA = indB'`: Tests the hypothesis that an indirect effect is equal in two or more groups. This is currently only possible for models without latent variables.
#' 
#' The indirect effect of interest can be specified in two ways:
#' * If a simple mediation involving three variables of the form `X -> M -> Y` is assumed, the arguments
#' `bYX`, `bMX`, and `bYM` are used to define the respective slopes, e. g.  `bYX = .4`, `bMX = .5`, and `bYM = .3` translates to
#' `X -- .5 --> M -- .3 --> Y` and  `X -- .4 --> Y`.
#' * More complex mediation structures can be defined by providing the `Beta` matrix along with `indirect` specifying which paths define the indirect effect. See examples below.
#' 
#' Notes on implementation:
#' * For models without latent variables, `nullEffect = 'ind = 0'` and `nullEffect = 'indA = indB'` constrain the indirect effect to zero and to equality, respectively, yielding the test described in Tofighi & Kelley (2020).
#' * For models with latent variables and `nullEffect = 'ind = 0'`, power is (sometimes roughly) approximated by constraining the smallest slope contained in the indirect effect to zero.
#' * For models with latent variables multiple groups (i. e., `nullEffect = 'indA = indB'`), there is currently no way to determine power. 
#' 
#' Tofighi, D., & Kelley, K. (2020). Improved inference in mediation analysis: Introducing the model-based constrained optimization procedure. *Psychological Methods, 25(4)*, 496–515. https://doi.org/10.1037/met0000259
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model**:
#' * `Lambda`: The factor loading matrix (with the number of columns equaling the number of factors).
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading on the second factor by .8, .6, .6, and .4.
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the loadings of the first factor to equal .5 and those of the second factor do equal .6.
#' 
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` need to be defined.
#' If the model contains observed variables only, use `Lambda = diag(x)` where `x` is the number of variables.
#'
#' Note that in case of a simple mediation model involving three variables, the order of the factors is X, M, Y, i. e., the first factor is treated as X, the second as M, and the thrird as Y. In case of a more complex mediation defined via the `Beta` matrix, the order of factors matches the order of `Beta`. 
#' 
#' Additional arguments related to the requested type of **power analysis**:
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.  This requires the `semTools` package.
#' 
#' @examples
#' \dontrun{
#' # simple case of X -> M -> Y mediation in the form of
#' # X -- .30 --> M -- .40 --> Y
#' # X --------- .25 --------> Y
#' # determine the required N to detect the indirect effect of >= .12 (= .3 * .4) 
#' # with a power of 95% on alpha = 5%, where   
#' # X is measured by 3 indicators loading by .5 each, 
#' # M is measured by 5 indicators loading by .6 each, 
#' # Y is measured by 4 indicators loading by .7 each.
#' powerMed <- semPower.powerMediation(type = 'a-priori',
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(3, 5, 4),
#'                                     loadM = c(.5, .6, .7),
#'                                     alpha = .05, beta = .05)
#' # show summary
#' summary(powerMed)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerMed$modelH1, sample.cov = powerMed$Sigma,
#' sample.nobs = powerMed$requiredN, sample.cov.rescale = FALSE)
#' lavaan::sem(powerMed$modelH0, sample.cov = powerMed$Sigma,
#' sample.nobs = powerMed$requiredN, sample.cov.rescale = FALSE)
#' 
#' # same as above, but determine power with N = 500 on alpha = .05
#' powerMed <- semPower.powerMediation(type = 'post-hoc',
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(3, 5, 4),
#'                                     loadM = c(.5, .6, .7),
#'                                     alpha = .05, N = 500)
#' 
#' # same as above, but determine the critical chi-square with N = 500 so that alpha = beta
#' powerMed <- semPower.powerMediation(type = 'compromise',
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(3, 5, 4),
#'                                     loadM = c(.5, .6, .7),
#'                                     abratio = 1, N = 500)
#' 
#' # same as above, but compare to the saturated model
#' # (rather than to the less restricted model)
#' powerMed <- semPower.powerMediation(type = 'a-priori',
#'                                     comparison = 'saturated',
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(3, 5, 4),
#'                                     loadM = c(.5, .6, .7),
#'                                     alpha = .05, beta = .05)
#' 
#' # same as above, but assuming observed variables only (Lambda = diag(3))
#' powerMed <- semPower.powerMediation(type = 'a-priori',
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     Lambda = diag(3),
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # same mediation model as above, but specifying Beta and indirect
#' Beta <- matrix(c(
#'   #   X    M    Y
#'   c(.00, .00, .00),    # X
#'   c(.30, .00, .00),    # M
#'   c(.25, .40, .00)     # Y
#' ), byrow = TRUE, ncol = 3)
#' powerMed <- semPower.powerMediation(type = 'a-priori',
#'                                     Beta = Beta, 
#'                                     indirect = list(c(2, 1), c(3, 2)),
#'                                     nIndicator = c(3, 5, 4),
#'                                     loadM = c(.5, .6, .7),
#'                                     alpha = .05, beta = .05)
#' 
#' # Beta for a more complex mediation hypothesis
#' # of the form X -- .2 --> M1 -- .3 --> M2 -- .40 -> Y 
#' # (and all other effects being zero)
#' # using a reduced loading matrix to define that
#' # X is measured by 3 indicators loading by .4, .5, .8 
#' # M1 is measured by 4 indicators loading by .7, .6, .5, .8
#' # M2 is measured by 5 indicators loading by .5, .6, .3, .4, .6 
#' # Y is measured by 4 indicators loading by .6, .7, .8
#' Beta <- matrix(c(
#'   c(.00, .00, .00, .00),       # X
#'   c(.20, .00, .00, .00),       # M1
#'   c(.00, .30, .00, .00),       # M2
#'   c(.00, .00, .40, .00)        # Y
#' ), byrow = TRUE, ncol = 4)
#' loadings <- list(
#'   c(0.4, 0.5, 0.8),           # X
#'   c(0.7, 0.6, 0.5, 0.8),      # M1
#'   c(0.5, 0.6, 0.3, 0.4, 0.6), # M2
#'   c(0.6, 0.7, 0.8)            # Y
#' )
#' powerMed <- semPower.powerMediation(type = 'a-priori',
#'                                     Beta = B, 
#'                                     indirect = list(c(2, 1), 
#'                                                     c(3, 2), 
#'                                                     c(4, 3)),
#'                                     loadings = loadings,
#'                                     alpha = .05, beta = .05)
#' 
#' # Determine required N to detect that the indirect effect 
#' # in group 1 (of .2 * .3 = .09) differs from the indirect effect 
#' # in group 2 (of .3 * .5 = .15).
#' # The direct effect of X on Y is .25 in both groups.  
#' # The model is based on observed variables only (Lambda = diag(3))
#' # Both groups are sized equally (N = list(1, 1)).
#' powerMed <- semPower.powerMediation(type = 'a-priori',
#'                                     nullEffect = 'indA = indB',
#'                                     bYX = list(.25, .25), 
#'                                     bMX = list(.2, .3), 
#'                                     bYM = list(.3, .5),
#'                                     Lambda = diag(3),
#'                                     alpha = .05, beta = .05, N = list(1, 1))
#' 
#' # same as above, but specifying Beta 
#' Beta1 <- matrix(c(
#'   c(.00, .00, .00),    # X
#'   c(.20, .00, .00),    # M
#'   c(.25, .30, .00)     # Y
#' ), byrow = TRUE, ncol = 3)
#' Beta2 <- matrix(c(
#'   c(.00, .00, .00),    # X
#'   c(.30, .00, .00),    # M
#'   c(.25, .50, .00)     # Y
#' ), byrow = TRUE, ncol = 3)
#' powerMed <- semPower.powerMediation(type = 'a-priori',
#'                                     nullEffect = 'indA = indB',
#'                                     Beta = list(Beta1, Beta2), 
#'                                     indirect = list(c(2, 1), c(3, 2)),
#'                                     Lambda = diag(3),
#'                                     alpha = .05, beta = .05, N = list(1, 1))
#' 
#' # request a simulated post-hoc power analysis with 500 replications.
#' set.seed(300121)
#' powerMed <- semPower.powerMediation(type = 'post-hoc',
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(3, 5, 4),
#'                                     loadM = c(.5, .6, .7),
#'                                     alpha = .05, N = 500,
#'                                     simulatedPower = TRUE, 
#'                                     simOptions = list(nReplications = 500))
#'}
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerMediation <- function(type, comparison = 'restricted',
                                    bYX = NULL, bMX = NULL, bYM = NULL,
                                    Beta = NULL, indirect = NULL, 
                                    nullEffect = 'ind = 0',
                                    nullWhichGroups = NULL,
                                    standardized = TRUE,
                                    ...){
  
  args <- list(...)
  comparison <- checkComparisonModel(comparison)
  checkEllipsis(...)
  
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
  if(standardized){
    Phi <- lapply(B, getPhi.B)
    generated <- semPower.genSigma(Phi = if(!isMultigroup) Phi[[1]] else Phi, 
                                   useReferenceIndicator = TRUE, ...)
  }else{
    generated <- semPower.genSigma(Beta = if(!isMultigroup) B[[1]] else B, 
                                   useReferenceIndicator = TRUE, ...)
  }
  if(!isMultigroup) isObserved <- ncol(generated[['Sigma']]) == ncol(B[[1]]) else isObserved <- ncol(generated[[1]][['Sigma']]) == ncol(B[[1]])  
  
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
  # so as a temporary workaround we apply different approaches depending on whether there are latent variables
  if(nullEffect == 'ind=0'){
    if(isObserved){
      # constrain indirect effect
      modelH0 <- paste(model, '\n', 'ind == 0')
    }else{
      # constrain the smallest of the contained direct effects, so this
      # actually gives power for a single slope. whereas this works much better than using 
      # ind=0 as comparison model, power is only approximated.
      # this approach should be replaced by the proper way of constraining 
      # the indirect effect to zero, once lav supports a suited optimizer
      # such as NPSOL or SLSQP.
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
  if(isMultigroup){
    args[['lavOptions']] <- append(args[['lavOptions']], list(group.equal = c('loadings')))
  } 
  
  modelH1 <- NULL
  fitH1model <- FALSE
  if(comparison == 'restricted'){
    # h1 model always fits perfectly, only needed for delta df
    modelH1 <- model
    # single group case: the h1 model always fits perfectly
    # multigroup case: we cannot be sure that user input yields a perfectly fitting model
    fitH1model <- isMultigroup 
  } 
  
  if(isMultigroup) Sigma <- lapply(generated, '[[', 'Sigma') else Sigma <- generated[['Sigma']] 
  
  do.call(semPower.powerLav, append(list(
    type = type,
    Sigma = Sigma,
    modelH0 = modelH0,
    modelH1 = modelH1,
    fitH1model = fitH1model),
    args)
  )

}


#' semPower.powerCLPM
#'
#' Convenience function for performing power analysis on effects in a cross-lagged panel model (CLPM).
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param nWaves number of waves, must be >= 2.
#' @param autoregEffects vector of the autoregressive effects of X and Y (constant across waves), or a list of vectors of autoregressive effects for X and Y from wave to wave, e.g. `list(c(.7, .6), c(.5, .5))` for a autoregressive effect of .7 for `X1 -> X2` and .6 for `X2 -> X3` and autoregressive effects of .5 for `Y1 -> Y2` and `Y2 -> Y3`. Must be a list of lists for multiple groups models. If the list structure is omitted, no group differences are assumed. 
#' @param crossedEffects vector of crossed effects of X on Y `(X -> Y)` and vice versa (both constant across waves), or a list of vectors of crossed effects giving the crossed effect of X on Y (and vice versa) for each wave, e.g. `list(c(.2, .3), c(.1, .1))` for `X1 - > Y2` = .2, `X2 -> Y3` = .3, `Y1 -> Y2` = .1, and `Y2 -> Y3` = .1. Must be a list of lists for multiple groups models. If the list structure is omitted, no group differences are assumed.
#' @param rXY vector of (residual-)correlations between X and Y for each wave. If `NULL`, all (residual-)correlations are zero. Can be a list for multiple groups models, otherwise no no group differences are assumed.
#' @param waveEqual parameters that are assumed to be equal across waves in both the H0 and the H1 model. Valid are `'autoregX'` and `'autoregY'` for autoregressive effects, `'crossedX'` and `'crossedY'` for crossed effects, `'corXY'` for residual correlations, or `NULL` for none (so that all parameters are freely estimated, subject to the constraints defined in `nullEffect`). 
#' @param nullEffect defines the hypothesis of interest. Valid are the same arguments as in `waveEqual` and additionally `'autoregX = 0'`, `'autoregY = 0'`, `'crossedX = 0'`, `'crossedY = 0'` to constrain the X or Y autoregressive effects or the crossed effects to zero, `'autoregX = autoregY'` and `'crossedX = crossedY'` to constrain them to be equal for X and Y, and `'autoregXA = autoregXB'`, `'autoregYA = autoregYB'`, `'crossedXA = crossedXB'`, `'crossedYA = crossedYB'` to constrain them to be equal across groups. 
#' @param nullWhich used in conjunction with `nullEffect` to identify which parameter to constrain when there are > 2 waves and parameters are not constant across waves. For example, `nullEffect = 'autoregX = 0'` with `nullWhich = 2` would constrain the second autoregressive effect for X to zero.    
#' @param nullWhichGroups for hypothesis involving cross-groups comparisons, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param standardized whether all parameters should be standardized (`TRUE`, the default). If `FALSE`, all regression relations are unstandardized.
#' @param metricInvariance whether metric invariance over waves is assumed (`TRUE`, the default) or not (`FALSE`). This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model).
#' @param autocorResiduals whether the residuals of the indicators of latent variables are autocorrelated over waves (`TRUE`, the default) or not (`FALSE`). This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model).
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. The order of factors is (X1, Y1, X2, Y2, ..., X_nWaves, Y_nWaves). See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details
#' This function performs a power analysis to reject various hypotheses arising
#' in crossed-lagged panel models (CLPM). In a standard CLPM implemented here, 
#' two variables X and Y are repeatedly assessed at two or more different time points (`nWaves`), yielding 
#' autoregressive effects (stabilities; X1 -> X2 and Y1 -> Y2),
#' synchronous effects (X1 <-> Y1, X2 <-> Y2), and cross-lagged effects (X1 -> Y2 and Y1 -> X2). 
#' CLPM including more than two waves are typically implemented assuming that the parameters are constant across waves (`waveEqual`), and usually omit lag-2 effects (e.g., X1 -> Y3). 
#' CLPM based on latent factors usually assume at least metric invariance of the factors over waves (`metricInvariance`).
#'  
#' Relevant hypotheses in arising in a CLPM are:
#' * `autoregX = 0` and `autoregY = 0`: Tests the hypothesis that the autoregressive effect of X and Y, respectively, is zero. 
#' * `crossedX = 0` and `crossedY = 0`: Tests the hypothesis that the crossed effect of X on Y (`crossedX`) and of Y on X (`crossedY`), respectively, is zero. 
#' * `autoregX = autoregY`: Tests the hypothesis that the autoregressive effect of X and Y are equal.
#' * `crossedX = crossedY`: Tests the hypothesis that the crossed effect of X on Y (`crossedX`) and of Y on X (`crossedY`) are equal.
#' * `autoregX` and `autoregY`: Tests the hypothesis that the autoregressive effect of X and Y, respectively, is equal across waves. 
#' * `crossedX` and `crossedY`: Tests the hypothesis that the crossed effect of X on Y (`crossedX`) and of Y on X (`crossedY`), respectively, is equal across waves. 
#' * `corXY`: Tests the hypothesis that the (residual-)correlations between X and Y are equal across waves. 
#' * `autoregXA = autoregXB` and `autoregYA = autoregYB`: Tests the hypothesis that the autoregressive effect of either X or Y are equal across groups.
#' * `crossedXA = crossedXB` and `crossedYA = crossedYB`: Tests the hypothesis that the crossed effect of X on Y (`crossedX`) or of Y on X (`crossedY`), respectively, is equal across groups.
#' 
#' For hypotheses regarding the random-intercept CLPM, see [semPower.powerRICLPM()]. For hypothesis in autoregressive models, see [semPower.powerAutoregressive()].
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model**:
#' * `Lambda`: The factor loading matrix (with the number of columns equaling the number of factors).
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading on the second factor by .8, .6, .6, and .4.
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the loadings of the first factor to equal .5 and those of the second factor do equal .6.
#' 
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` need to be defined.
#' If the model contains observed variables only, use `Lambda = diag(x)` where `x` is the number of variables.
#'
#' Note that the order of the factors is (X1, Y1, X2, Y2, ..., X_nWaves, Y_nWaves), i. e., the first factor is treated as the first measurement of X, the second as the first measurement of Y, the third as the second measurement of X, etc.. 
#' 
#' Additional arguments related to the requested type of **power analysis**:
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.  This requires the `semTools` package.
#' 
#' @examples
#' \dontrun{
#' # Determine required N in a 2-wave CLPM
#' # to detect a crossed-effect of X (X1 -> Y2) of >= .2 
#' # with a power of 95% on alpha = 5%, where
#' # X1 and X2 are measured by 5 indicators loading by .5 each, and
#' # Y1 and Y2 are measured by 3 indicators loading by .6 each, and
#' # there is no synchronous correlation between X and Y (rXY = NULL), 
#' # the stability of X is .8,
#' # the stability of Y is .7, and
#' # the crossed-effect of Y (Y1 -> X2) is .1.
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'crossedX = 0',
#'                                 nIndicator = c(5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#' 
#' # show summary
#' summary(powerCLPM)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerCLPM$modelH1, sample.cov = powerCLPM$Sigma,
#'             sample.nobs = powerCLPM$requiredN, 
#'             sample.cov.rescale = FALSE)
#' lavaan::sem(powerCLPM$modelH0, sample.cov = powerCLPM$Sigma,
#'             sample.nobs = powerCLPM$requiredN, 
#'             sample.cov.rescale = FALSE)
#' 
#' # same as above, but determine power with N = 500 on alpha = .05
#' powerCLPM <- semPower.powerCLPM(type = 'post-hoc',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'crossedX = 0',
#'                                 nIndicator = c(5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6),
#'                                 alpha = .05, N = 500)
#' 
#' # same as above, but determine the critical chi-square with N = 500 so that alpha = beta
#' powerCLPM <- semPower.powerCLPM(type = 'compromise',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'crossedX = 0',
#'                                 nIndicator = c(5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6),
#'                                 abratio = 1, N = 500)
#' 
#' # same as above, but compare to the saturated model
#' # (rather than to the less restricted model)
#' powerCLPM <- semPower.powerCLPM(type = 'compromise',
#'                                 comparison = 'saturated',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'crossedX = 0',
#'                                 nIndicator = c(5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6),
#'                                 abratio = 1, N = 500)
#' 
#' # same as above, but assume only observed variables
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'crossedX = 0',
#'                                 Lambda = diag(4),
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but provide reduced loadings matrix to define that
#' # X1 and X2 are measured by 5 indicators each loading by .4, .5, .6, .5, .4 
#' # Y1 and Y2 are measured by 3 indicators each loading by .8, .6, .7
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'crossedX = 0',
#'                                 loadings = list(
#'                                   c(.4, .5, .6, .5, .4),    # X1
#'                                   c(.8, .6, .7),            # Y1
#'                                   c(.4, .5, .6, .5, .4),    # X2
#'                                   c(.8, .6, .7)             # Y2
#'                                 ),
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but do not assume metric invariance across waves
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'crossedX = 0',
#'                                 nIndicator = c(5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6),
#'                                 metricInvariance = FALSE,
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but determine N to detect that the crossed-effect of Y (Y1 -> X2) is >= .1.
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'crossedY = 0',
#'                                 nIndicator = c(5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but determine N to detect that the stability of X (X1 -> X2) is >= .8.
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'autoregX = 0',
#'                                 nIndicator = c(5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but determine N to detect that the stability of Y (Y1 -> Y2) is >= .7.
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'autoregY = 0',
#'                                 nIndicator = c(5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but determine N to detect that 
#' # the crossed effect of X (X1 -> Y2) of .2 differs from 
#' # the crossed effect of Y (Y1 -> X2) of .1  
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'crossedX = crossedY',
#'                                 nIndicator = c(5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but determine N to detect that 
#' # the autoregressive effect of X (X1 -> X2) of .8 differs from 
#' # the autoregressive effect of Y (Y1 -> Y2) of .7  
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'autoregX = autoregY',
#'                                 nIndicator = c(5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but assume that the synchronous correlation between X and Y 
#' # is .3 at the first wave, and the respective residual correlation is .2 at the second wave, 
#' # and determine N to detect that synchronous residual correlation (at wave 2) is => .2.
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = c(.3, .2),
#'                                 nullEffect = 'corXY = 0',
#'                                 nIndicator = c(5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#' 
#' # Determine required N in a 3-wave CLPM
#' # to detect a crossed-effect of X (X1 -> Y2 and X2 -> Y3) of >= .2 
#' # with a power of 95% on alpha = 5%, where
#' # the crossed, autoregressive, and synchronous effects of X and Y are equal over waves,
#' # X1, X2, and X3 are measured by 5 indicators loading by .5 each, and
#' # Y1, Y2, and Y3 are measured by 3 indicators loading by .6 each, and
#' # the synchronous correlation between X and Y are .2 across all three waves, and
#' # the stability of X is .8 across all three waves,
#' # the stability of Y is .7 across all three waves, and
#' # the crossed-effect of Y (Y1 -> X2, and Y2 -> Y3) is .1.
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 3,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = c(.2, .2, .2),
#'                                 waveEqual = c('autoregX', 'autoregY', 
#'                                               'crossedX', 'crossedY'),
#'                                 nullEffect = 'crossedX = 0',
#'                                 nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#' 
#' # Determine required N in a 3-wave CLPM to detect that 
#' # the crossed-effect of X in wave 1 (X1 -> Y2) of .20 is equal to the 
#' # the crossed-effect of X in wave 2 (X2 -> Y3) of .10 
#' # with a power of 95% on alpha = 5%, where
#' # the autoregressive effects of X and Y are equal over waves,
#' # X1, X2, and X3 are measured by 5 indicators loading by .5 each, and
#' # Y1, Y2, and Y3 are measured by 3 indicators loading by .6 each, and
#' # the synchronous correlation between X and Y are .2, .3, and .4 at the first, 
#' # second, and third wave, and
#' # the stability of X is .8 across all three waves,
#' # the stability of Y is .7 across all three waves, and
#' # the crossed-effects of Y (Y1 -> X2, and Y2 -> X3) are both .1 
#' # (but freely estimated for each wave).
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 3,
#'                                 autoregEffects = c(.8, .7), 
#'                                 crossedEffects = list(
#'                                   c(.20, .10),   # X1 -> Y2, X2 -> Y3
#'                                   c(.05, .10)),  # Y1 -> X2, Y2 -> X3
#'                                 rXY = c(.2, .3, .4),
#'                                 nullEffect = 'crossedX',
#'                                 waveEqual = c('autoregX', 'autoregY'),
#'                                 nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but determine N to detect that 
#' # the crossed-effect of X at wave 2 is >= .10.
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 3,
#'                                 autoregEffects = c(.8, .7), 
#'                                 crossedEffects = list(
#'                                   c(.20, .10),   # X1 -> Y2, X2 -> Y3
#'                                   c(.05, .10)),  # Y1 -> X2, Y2 -> X3
#'                                 rXY = c(.2, .3, .4),
#'                                 nullEffect = 'crossedX',
#'                                 nullWhich = 2,
#'                                 waveEqual = c('autoregX', 'autoregY'),
#'                                 nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but determine N to detect that 
#' # the residual correlation between X and Y at wave 2 (of .3) differs from 
#' # the residual correlation between X and Y at wave 3 (of .4)
#' # and define unstandardized parameters
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 3,
#'                                 autoregEffects = c(.8, .7), 
#'                                 crossedEffects = list(
#'                                   c(.20, .10),   # X1 -> Y2, X2 -> Y3
#'                                   c(.05, .10)),  # Y1 -> X2, Y2 -> X3
#'                                 rXY = c(.2, .3, .4),
#'                                 nullEffect = 'corXY',
#'                                 waveEqual = c('autoregX', 'autoregY'),
#'                                 standardized = FALSE,
#'                                 nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#'
#'
#' # multiple group example
#' # determine power in a 3-wave CLPM to detect that 
#' # the autoregressive effect of X in group 1 (of .8) differs from the 
#' # autoregressive effect of X in group 2 (of .6)
#' # with a 500 observations in both groups on alpha = 5%, where
#' # the autoregressive effects of X and Y are equal over waves (but not across groups),
#' # the cross-lagged effects of X and Y are equal over waves (and also across groups),
#' # X1, X2, and X3 are measured by 5 indicators loading by .5 each, and
#' # Y1, Y2, and Y3 are measured by 3 indicators loading by .4 each, and
#' # there are no synchronous correlation between X and Y.
#' powerCLPM <- semPower.powerCLPM(type = 'post-hoc', alpha = .05, N = list(500, 500),
#'                                 nWaves = 3,
#'                                 autoregEffects = list(
#'                                 # group 1
#'                                 list(c(.8, .8),    # X1 -> X2, X2 -> X3 
#'                                      c(.7, .7)),   # Y1 -> Y2, Y2 -> Y3
#'                                 # group 2
#'                                 list(c(.6, .6),    # X1 -> X2, X2 -> X3 
#'                                      c(.7, .7))    # Y1 -> Y2, Y2 -> Y3
#'                                 ),
#'                                 crossedEffects = c(.2, .1),
#'                                 waveEqual = c('autoregX', 'autoregY', 'crossedX', 'crossedY'),
#'                                 rXY = NULL,
#'                                 nullEffect = 'autoregxa=autoregxb',
#'                                 nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                 loadM = c(.5, .4, .5, .4, .5, .4))
#'                                 
#' 
#' # request a simulated post-hoc power analysis with 500 replications.
#' set.seed(300121)
#' powerCLPM <- semPower.powerCLPM(type = 'post-hoc',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'crossedX = 0',
#'                                 Lambda = diag(4),
#'                                 alpha = .05, N = 500, 
#'                                 simulatedPower = TRUE, 
#'                                 simOptions = list(nReplications = 500))
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerCLPM <- function(type, comparison = 'restricted',
                               nWaves = NULL, 
                               autoregEffects = NULL, 
                               crossedEffects = NULL, 
                               rXY = NULL,
                               waveEqual = NULL, 
                               nullEffect = NULL, 
                               nullWhich = NULL,
                               nullWhichGroups = NULL,
                               standardized = TRUE,
                               metricInvariance = TRUE,
                               autocorResiduals = TRUE,
                               ...){
  
  # TODO: lagged effects would be nice

  comparison <- checkComparisonModel(comparison)
  checkEllipsis(...)
  
  # we override Beta and Sigma later, so let's make sure it is not set in ellipsis argument
  if('Beta' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Beta.')
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma.')
  
  # validate input
  if(is.null(autoregEffects) || is.null(crossedEffects)) stop('autoregEffects and crossedEffects may not be NULL.')
  if(is.null(nWaves) || is.na(nWaves) || nWaves < 2) stop('nWaves must be >= 2.')

  # we do not allow stacking of hypotheses. there might be a use case for this,
  # but this would complicate defining the relevant parameter when these vary across waves. 
  nullValid <- c('autoregx', 'autoregy', 'crossedx', 'crossedy', 'corxy',
                 'autoregx=0', 'autoregy=0', 'crossedx=0', 'crossedy=0',
                 'autoregx=autoregy', 'crossedx=crossedy', 'corxy=0',
                 'autoregxa=autoregxb', 'autoregya=autoregyb', 'crossedxa=crossedxb', 'crossedya=crossedyb'
  )
  nullEffect <- checkNullEffect(nullEffect, nullValid)
  
  # create list structure for autoregEffects, crossedEffects, and corXY
  # assume multigroup when list structure is present for either autoreg or crossed effects
  ngA <- ifelse(is.list(autoregEffects[[1]]), length(autoregEffects), 1)
  ngX <- ifelse(is.list(crossedEffects[[1]]), length(crossedEffects), 1)
  if(sum(c(ngA, ngX) > 1) > 1  && ngA != ngX) stop('Specify the same number of groups for both autoregEffects and crossedEffects.')
  nGroups <- max(c(ngA, ngX))
  isMultigroup <- nGroups > 1
  
  if(isMultigroup && !nullEffect %in% c('autoregxa=autoregxb', 'autoregya=autoregyb', 'crossedxa=crossedxb', 'crossedya=crossedyb')) stop('Multigroup analysis are only supported for nullEffect = autoregxa=autoregxb, autoregya=autoregyb, crossedxa=crossedxb, crossedya=crossedyb')
  if(!isMultigroup && nullEffect %in% c('autoregxa=autoregxb', 'autoregya=autoregyb', 'crossedxa=crossedxb', 'crossedya=crossedyb')) stop('nullEffect = autoregxa=autoregxb, autoregya=autoregyb, crossedxa=crossedxb, crossedya=crossedyb imply multigroup analyses, but no list structure for any relevant parameter provided.')
  if(isMultigroup && is.null(nullWhichGroups)) nullWhichGroups <- seq(nGroups)

  # [[groups]][[X, Y]][[waves]]
  if(!is.list(autoregEffects)) autoregEffects <- list(rep(autoregEffects[[1]], (nWaves - 1)), rep(autoregEffects[[2]], (nWaves - 1)))
  if(!is.list(crossedEffects)) crossedEffects <- list(rep(crossedEffects[[1]], (nWaves - 1)), rep(crossedEffects[[2]], (nWaves - 1)))
  if(!is.list(autoregEffects[[1]])) autoregEffects <- rep(list(autoregEffects), nGroups)
  if(!is.list(crossedEffects[[1]])) crossedEffects <- rep(list(crossedEffects), nGroups)
  if(length(autoregEffects[[1]][[1]]) == 1) autoregEffects <- lapply(autoregEffects, function(x) lapply(x, function(y) rep(y, nWaves - 1)))
  if(length(crossedEffects[[1]][[1]]) == 1) crossedEffects <- lapply(crossedEffects, function(x) lapply(x, function(y) rep(y, nWaves - 1)))

  if(is.null(rXY)) rXY <- rep(0, nWaves)
  if(is.list(rXY) && length(rXY) != nGroups) stop('corXY implies a different number of groups as autoregEffects or crossedEffects.')
  if(!is.list(rXY)) rXY <- rep(list(rXY), nGroups)

  if(any(unlist(lapply(rXY, function(x) length(x) != nWaves)))) stop('rXY must be of length nWaves')
  invisible(lapply(rXY, function(y) lapply(y, function(x) checkBounded(x, 'All rXY ', bound = c(-1, 1), inclusive = FALSE))))

  if(any(unlist(lapply(autoregEffects, function(x) length(x) != 2)))) stop('Provide autoregEffects for X and Y.')
  if(any(unlist(lapply(crossedEffects, function(x) length(x) != 2)))) stop('Provide crossedEffects for X and Y..')
  if(any(unlist(lapply(autoregEffects, function(x) length(x[[1]]) != length(x[[2]]))))) stop('autoregEffects for X and Y must be of equal length.')
  if(any(unlist(lapply(crossedEffects, function(x) length(x[[1]]) != length(x[[2]]))))) stop('crossedEffects for X and Y must be of equal length.')
  if(any(unlist(lapply(autoregEffects, function(x) length(x[[1]]) != (nWaves - 1))))) stop('autoregEffects must be of length nWaves - 1.')
  if(any(unlist(lapply(crossedEffects, function(x) length(x[[1]]) != (nWaves - 1))))) stop('crossedEffects must be of length nWaves - 1.')
  invisible(lapply(autoregEffects, function(y) lapply(y, function(x) lapply(x, function(x) checkBounded(x, 'All autoregressive effects ', bound = c(-1, 1), inclusive = FALSE)))))
  invisible(lapply(crossedEffects, function(y) lapply(y, function(x) lapply(x, function(x) checkBounded(x, 'All crossed effects ', bound = c(-1, 1), inclusive = FALSE)))))
  
  if(!is.null(waveEqual)){
    waveEqual <- unlist(lapply(waveEqual, function(x) tolower(trimws(x))))
    if(any(unlist(lapply(waveEqual, function(x) !x %in% c('autoregx', 'autoregy', 'crossedx', 'crossedy', 'corxy'))))) stop('waveEqual may only contain autoregX, autoregY, crossedX, crossedY, corXY')
  }
  
  if(any(unlist(lapply(nullEffect, function(x) !x %in% nullValid)))) stop('Unknown value for nullEffect.')
  if(any(nullEffect %in% waveEqual)) stop('You cannot set the same parameters in nullEffect and waveEqual.')
  if(nWaves == 2 && nullEffect %in% c('autoregx', 'autoregy','crossedx', 'crossedy', 'corxy')) stop('For two waves, there is only one crossedX and crossedY effect, only one autoregressive effect each, and only one X-Y residual correlation. Did you mean crossedX = 0 or autoregX = 0?')
  
  if(is.null(nullWhich) && nWaves == 2) nullWhich <- 1
  if(is.null(nullWhich) && nWaves > 2){
    msg <- 'nullWhich must be defined when there are more than 2 waves and relevant parameters are not constant across waves.'
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
  if('corxy' %in% nullEffect && standardized) stop('Power analysis for nullEffect == "corxy" can only be performed on unstandardized parameters. Repeat with standardized = FALSE.') 
  
  
  ### create B
  Beta <- lapply(seq(nGroups), function(x){
    B <- matrix(0, ncol = 2*nWaves, nrow = 2*nWaves)
    # add autoregEffects and crossed-effects
    for(i in 1:(nWaves - 1)){
      xidx <- 2 + 2*(i - 1) + 1
      yidx <- xidx + 1
      # autoregEffects
      B[xidx, (xidx - 2)] <- autoregEffects[[x]][[1]][i]
      B[yidx, (yidx - 2)] <- autoregEffects[[x]][[2]][i]
      # crossed effects
      B[yidx, (xidx - 2)] <- crossedEffects[[x]][[1]][i]
      B[xidx, (yidx - 2)] <- crossedEffects[[x]][[2]][i]
    }
    B
  })
  
  ### create Psi
  Psi <- lapply(seq(nGroups), function(x){
    P <- diag(ncol(Beta[[1]]))
    if(any(rXY[[x]] != 0)){
      for(i in 1:nWaves){
        P[2*i, (2*i - 1)] <- P[(2*i - 1), 2*i] <- rXY[[x]][i]
      }
    }
    P
  })
  
  # add metric invariance constraints to analysis model
  metricInvarianceFactors <- NULL
  if(metricInvariance){
    metricInvarianceFactors <- list(
      seq(1, 2*nWaves, 2),
      seq(2, 2*nWaves, 2)  
    )
  }
  
  ### get Sigma
  if(standardized){
    Phi <- lapply(seq_along(Beta), function(x) getPhi.B(Beta[[x]], Psi[[x]]))
    generated <- semPower.genSigma(Phi = if(!isMultigroup) Phi[[1]] else Phi, 
                                   useReferenceIndicator = TRUE,
                                   metricInvariance = metricInvarianceFactors, 
                                   nGroups = nGroups,
                                   ...)
  }else{  
    generated <- semPower.genSigma(Beta = if(!isMultigroup) Beta[[1]] else Beta, 
                                   Psi = if(!isMultigroup) Psi[[1]] else Psi, 
                                   useReferenceIndicator = TRUE,
                                   metricInvariance = metricInvarianceFactors, 
                                   nGroups = nGroups,
                                   ...)
  }
  
  ### create model strings
  if(!isMultigroup) model <- generated[['modelTrueCFA']] else model <- generated[[1]][['modelTrueCFA']]
  
  # add CLPM structure 
  for(f in 3:ncol(Beta[[1]])){     # omit rows 1:2
    fidx <- which(Beta[[1]][f, ] != 0)
    if(length(fidx) != 0){
      tok <- paste0('f', f, ' ~ ', paste(paste0('pf', paste0(formatC(f, width = 2, flag = 0), formatC(fidx, width = 2, flag = 0)), '*'), paste0('f', fidx), sep = '', collapse = ' + '))
      model <- paste(model, tok, sep='\n')
    }
  }
  # add (residual) correlations 
  # (note that lav estimates these by default)
  model <- paste(model, 'f1 ~~ pf0201*f2', sep='\n')
  for(i in 2:nWaves){
    tok <- paste0('f',(2*i - 1),' ~~ ', paste0('pf', paste0(formatC(2*i, width = 2, flag = 0), formatC(2*i - 1, width = 2, flag = 0)), '*'), 'f', 2*i)
    model <- paste(model, tok, sep='\n')
  }
  
  # add autocorrelated residuals
  if(autocorResiduals){
    if(!isMultigroup) Lambda <- generated[['Lambda']] else Lambda <- generated[[1]][['Lambda']]
    # do this only when there is at least one latent variable
    if(nrow(Lambda) > 2*nWaves){
      autocorResidualsFactors <- metricInvarianceFactors  # same structure 
      tok <- list()
      for(x in seq_along(autocorResidualsFactors)){
        ci <- lapply(autocorResidualsFactors[[x]], function(f) paste0('x', which(Lambda[, f] != 0)))
        if(length(ci[[1]]) > 1){
          for(i in 1:(length(ci) - 1)){
            for(j in (i + 1) : length(ci)){
              tok <- append(tok, paste(ci[[i]], '~~', ci[[j]]))
            }
          }
        }
      }
      model <- paste(model, paste(tok, collapse = '\n'), sep ='\n')
    }
  }
  
  ### define H1 and H0 model
  
  # first get relevant parameter labels that may be subject to waveEqual or nullEffect constraints 
  xw <- seq(2*nWaves - 1, 2, -2)
  pAutoregX <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(xw - 2, width = 2, flag = 0))
  pAutoregX <- pAutoregX[order(pAutoregX)]
  
  yw <- seq(2*nWaves, 3, -2)
  pAutoregY <- paste0('pf', formatC(yw, width = 2, flag = 0), formatC(yw - 2, width = 2, flag = 0))
  pAutoregY <- pAutoregY[order(pAutoregY)]
  
  xw <- seq(2*nWaves - 3, 0, -2)
  yw <- seq(2*nWaves, 3, -2)
  pCrossedX <- paste0('pf', formatC(yw, width = 2, flag = 0), formatC(xw, width = 2, flag = 0))
  pCrossedX <- pCrossedX[order(pCrossedX)]
  
  xw <- seq(2*nWaves - 1, 2, -2)
  yw <- seq(2*nWaves - 2, 1, -2)
  pCrossedY <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(yw, width = 2, flag = 0))
  pCrossedY <- pCrossedY[order(pCrossedY)]
  
  xw <- seq(2*nWaves - 1, 2, -2)
  yw <- seq(2*nWaves, 3, -2)
  pCorXY <- paste0('pf', formatC(yw, width = 2, flag = 0), formatC(xw, width = 2, flag = 0))
  pCorXY <- pCorXY[order(pCorXY)]
  
  # multigroup case
  if(isMultigroup){
    # remove group specific labels from measurement part to enforce metric invariance
    model <- gsub(paste0('_g', seq(nGroups), collapse = '|'), '_gc', model)
    # assign group labels to all structural parameters
    patt <- 'pf0201'
    repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
    model <- gsub(patt, repl, model)
    for(pp in seq(nWaves - 1)){
      patt <- pAutoregX[pp]
      repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      model <- gsub(patt, repl, model)
      patt <- pAutoregY[pp]
      repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      model <- gsub(patt, repl, model)
      patt <- pCrossedX[pp]
      repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      model <- gsub(patt, repl, model)
      patt <- pCrossedY[pp]
      repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      model <- gsub(patt, repl, model)
      patt <- pCorXY[pp]
      repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      model <- gsub(patt, repl, model)
    }
  }
  
  # add constraints to H1 model
  modelH1 <- model
  if(!is.null(waveEqual)){
    if('autoregx' %in% waveEqual) modelH1 <- gsub(paste(pAutoregX, collapse = '|'), paste(pAutoregX, collapse = ''), modelH1)
    if('autoregy' %in% waveEqual) modelH1 <- gsub(paste(pAutoregY, collapse = '|'), paste(pAutoregY, collapse = ''), modelH1)
    if('crossedx' %in% waveEqual) modelH1 <- gsub(paste(pCrossedX, collapse = '|'), paste(pCrossedX, collapse = ''), modelH1)
    if('crossedy' %in% waveEqual) modelH1 <- gsub(paste(pCrossedY, collapse = '|'), paste(pCrossedY, collapse = ''), modelH1)
    if('corxy' %in% waveEqual) modelH1 <- gsub(paste(pCorXY, collapse = '|'), paste(pCorXY, collapse = ''), modelH1)
  }
  
  # add additional constraints to H0 model
  modelH0 <- modelH1
  
  # wave equal constraints not included in modelH1:
  if('autoregx' %in% nullEffect) modelH0 <- gsub(paste(pAutoregX, collapse = '|'), paste(pAutoregX, collapse = ''), modelH0)
  if('autoregy' %in% nullEffect) modelH0 <- gsub(paste(pAutoregY, collapse = '|'), paste(pAutoregY, collapse = ''), modelH0)
  if('crossedx' %in% nullEffect) modelH0 <- gsub(paste(pCrossedX, collapse = '|'), paste(pCrossedX, collapse = ''), modelH0)
  if('crossedy' %in% nullEffect) modelH0 <- gsub(paste(pCrossedY, collapse = '|'), paste(pCrossedY, collapse = ''), modelH0)
  if('corxy' %in% nullEffect) modelH0 <- gsub(paste(pCorXY, collapse = '|'), paste(pCorXY, collapse = ''), modelH0)
  
  # zero and equality constraints:
  if('autoregx=0' %in% nullEffect){
    if('autoregx' %in% waveEqual){
      modelH0 <- gsub(paste(pAutoregX, collapse = ''), '0', modelH0)
    }else{
      modelH0 <- gsub(pAutoregX[nullWhich], '0', modelH0)
    }
  }
  if('autoregy=0' %in% nullEffect){
    if('autoregy' %in% waveEqual){
      modelH0 <- gsub(paste(pAutoregY, collapse = ''), '0', modelH0)
    }else{
      modelH0 <- gsub(pAutoregY[nullWhich], '0', modelH0)
    }
  }
  if('crossedx=0' %in% nullEffect){
    if('crossedx' %in% waveEqual){
      modelH0 <- gsub(paste(pCrossedX, collapse = ''), '0', modelH0)
    }else{
      modelH0 <- gsub(pCrossedX[nullWhich], '0', modelH0)
    }
  }
  if('crossedy=0' %in% nullEffect){
    if('crossedy' %in% waveEqual){
      modelH0 <- gsub(paste(pCrossedY, collapse = ''), '0', modelH0)
    }else{
      modelH0 <- gsub(pCrossedY[nullWhich], '0', modelH0)
    }
  }
  if('autoregx=autoregy' %in% nullEffect){
    if('autoregx' %in% waveEqual && 'autoregy' %in% waveEqual){
      patt <- paste(c(paste(pAutoregX, collapse = ''), paste(pAutoregY, collapse = '')), collapse = '|')
    }else if('autoregx' %in% waveEqual){
      patt <- paste(c(paste(pAutoregX, collapse = ''), pAutoregY[nullWhich]), collapse = '|')
    }else if('autoregy' %in% waveEqual){
      patt <- paste(c(pAutoregX[nullWhich], paste(pAutoregY, collapse = '')), collapse = '|')
    }else{
      patt <- paste(c(pAutoregX[nullWhich], pAutoregY[nullWhich]), collapse = '|')
    }
    repl <-  gsub('\\|', '', patt)
    modelH0 <- gsub(patt, repl, modelH0)
  }
  if('crossedx=crossedy' %in% nullEffect){
    if('crossedx' %in% waveEqual && 'crossedy' %in% waveEqual){
      patt <- paste(c(paste(pCrossedX, collapse = ''), paste(pCrossedY, collapse = '')), collapse = '|')
    }else if('crossedx' %in% waveEqual){
      patt <- paste(c(paste(pCrossedX, collapse = ''), pCrossedY[nullWhich]), collapse = '|')
    }else if('crossedy' %in% waveEqual){
      patt <- paste(c(pCrossedX[nullWhich], paste(pCrossedY, collapse = '')), collapse = '|')
    }else{
      patt <- paste(c(pCrossedX[nullWhich], pCrossedY[nullWhich]), collapse = '|')
    }
    repl <-  gsub('\\|', '', patt)
    modelH0 <- gsub(patt, repl, modelH0)
  }
  if('corxy=0' %in% nullEffect){
    if('corxy' %in% waveEqual){
      modelH0 <- gsub(paste(pCorXY, collapse = ''), '0', modelH0)
    }else{
      pCorXY <- c('pf0201', pCorXY)   # add exog cor
      modelH0 <- gsub(pCorXY[nullWhich], '0', modelH0)
    }
  }
  # multigroup cases
  if('autoregxa=autoregxb' %in% nullEffect){
    if('autoregx' %in% waveEqual){
      patt <- paste0(paste(pAutoregX, collapse = ''), '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(paste(pAutoregX, collapse = ''), '_gc')
    }else{
      patt <- paste0(pAutoregX[nullWhich], '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(pAutoregX[nullWhich], '_gc')
    }
    modelH0 <- gsub(patt, repl, modelH0)
  }
  if('autoregya=autoregyb' %in% nullEffect){
    if('autoregy' %in% waveEqual){
      patt <- paste0(paste(pAutoregY, collapse = ''), '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(paste(pAutoregY, collapse = ''), '_gc')
    }else{
      patt <- paste0(pAutoregY[nullWhich], '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(pAutoregY[nullWhich], '_gc')
    }
    modelH0 <- gsub(patt, repl, modelH0)
  }
  if('crossedxa=crossedxb' %in% nullEffect){
    if('crossedx' %in% waveEqual){
      patt <- paste0(paste(pCrossedX, collapse = ''), '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(paste(pCrossedX, collapse = ''), '_gc')
    }else{
      patt <- paste0(pCrossedX[nullWhich], '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(pCrossedX[nullWhich], '_gc')
    }
    modelH0 <- gsub(patt, repl, modelH0)
  }
  if('crossedya=crossedyb' %in% nullEffect){
    if('crossedy' %in% waveEqual){
      patt <- paste0(paste(pCrossedY, collapse = ''), '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(paste(pCrossedY, collapse = ''), '_gc')
    }else{
      patt <- paste0(pCrossedY[nullWhich], '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(pCrossedY[nullWhich], '_gc')
    }
    modelH0 <- gsub(patt, repl, modelH0)
  }
  
  # here we actually fit modelH1 in case of a restricted comparison
  # because we cannot be sure that user input yields perfectly fitting h1 models, 
  # when there are additional constraints (waveequal or invariance)
  if(comparison == 'saturated') modelH1 <- NULL
  
  if(isMultigroup) Sigma <- lapply(generated, '[[', 'Sigma') else Sigma <- generated[['Sigma']] 
  
  semPower.powerLav(type, 
                    Sigma = Sigma,
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
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param nWaves number of waves, must be >= 3.
#' @param autoregEffects vector of the autoregressive effects of X and Y (constant across waves), or a list of vectors of autoregressive effects for X and Y from wave to wave, e.g. `list(c(.7, .6), c(.5, .5))` for an autoregressive effect of .7 for X1->X2 and .6 for X2->X3 and autoregressive effects of .5 for Y1->Y2 and Y2->Y3.
#' @param crossedEffects vector of crossed effects of X on Y (X -> Y) and vice versa (both constant across waves), or a list of vectors of crossed effects giving the crossed effect of X on Y (and vice versa) for each wave, e.g. `list(c(.2, .3), c(.1, .1))` for X1->Y2 = .2, X2->Y3 = .3, Y1->Y2 = .1, and Y2->Y3 = .1.
#' @param rXY vector of (residual-)correlations between X and Y for each wave. If `NULL`, all (residual-)correlations are zero. 
#' @param rBXBY correlation between random intercept factors. If `NULL`, the correlation is zero. 
#' @param waveEqual parameters that are assumed to be equal across waves in both the H0 and the H1 model. Valid are `'autoregX'` and `'autoregY'` for autoregressive effects, `'crossedX'` and `'crossedY'` for crossed effects, `'corXY'` for residual correlations, or `NULL` for none (so that all parameters are freely estimated, subject to the constraints defined in `nullEffect`). 
#' @param nullEffect defines the hypothesis of interest. Valid are the same arguments as in `waveEqual` and additionally `'autoregX = 0'`, `'autoregY = 0'`, `'crossedX = 0'`, `'crossedY = 0'` to constrain the X or Y autoregressive effects or the crossed effects to zero, `'corBXBY = 0'` to constrain the correlation between the random intercepts to zero, and `'autoregX = autoregY'` and `'crossedX = crossedY'` to constrain them to be equal for X and Y, and `'autoregXA = autoregXB'`, `'autoregYA = autoregYB'`, `'crossedXA = crossedXB'`, `'crossedYA = crossedYB'`, and `corBXBYA = corBXBYB` to constrain them to be equal across groups.
#' @param nullWhich used in conjunction with `nullEffect` to identify which parameter to constrain when there are > 2 waves and parameters are not constant across waves. For example, `nullEffect = 'autoregX = 0'` with `nullWhich = 2` would constrain the second autoregressive effect for X to zero.    
#' @param nullWhichGroups for hypothesis involving cross-groups comparisons, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param metricInvariance whether metric invariance over waves is assumed (`TRUE`, the default) or not (`FALSE`). This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model, where the df are not affected by invariance constraints).
#' @param autocorResiduals whether the residuals of the indicators of latent variables are autocorrelated over waves (`TRUE`, the default) or not (`FALSE`). This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model).
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. The order of factors is (X1, Y1, X2, Y2, ..., X_nWaves, Y_nWaves). See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details 
#' This function performs a power analysis to reject various hypotheses arising in a random intercept crossed-lagged panel model (RI-CLPM). 
#' In a standard RI-CLPM implemented here, two variables X and Y are repeatedly assessed at three or more different time points (`nWaves`), 
#' yielding autoregressive effects (`X1 -> X2`, `X2 -> X3`, `Y1 -> Y2`, `Y2 -> Y3`), synchronous effects (`X1 <-> Y1`, `X2 <-> Y2`, `X3 <-> Y3`), and cross-lagged effects (`X1 -> Y2`, `X2 -> Y3`, `Y1 -> X2`, `Y2 -> X3`). 
#' RI-CLPMs are typically implemented assuming that the parameters are constant across waves (`waveEqual`), and usually omit lag-2 effects (e.g., `X1 -> Y3`). 
#' RI-CLPMs based on latent factors usually assume at least metric invariance of the factors over waves (`metricInvariance`).
#'  
#' Relevant hypotheses in arising in a RI-CLPM are:
#' * `autoregX = 0` and `autoregY = 0`: Tests the hypothesis that the autoregressive effect of X and Y, respectively, is zero. 
#' * `crossedX = 0` and `crossedY = 0`: Tests the hypothesis that the crossed effect of X on Y (`crossedX`) and Y on X (`crossedY`), respectively, is zero. 
#' * `autoregX = autoregY`: Tests the hypothesis that the autoregressive effect of X and Y are equal.
#' * `crossedX = crossedY`: Tests the hypothesis that the crossed effect of X on Y (`crossedX`) and Y on X (`crossedY`) are equal.
#' * `autoregX` and `autoregY`: Tests the hypothesis that the autoregressive effect of X and Y, respectively, is equal across waves. 
#' * `crossedX` and `crossedY`: Tests the hypothesis that the crossed effect of X on Y (`crossedX`) and Y on X (`crossedY`), respectively, is equal across waves. 
#' * `corXY`: Tests the hypothesis that the (residual-)correlations between X and Y are equal across waves. 
#' * `corBXBY = 0`: Tests the hypothesis that the correlation between the random intercept factors of X and Y is zero.
#' * `autoregXA = autoregXB` and `autoregYA = autoregYB`: Tests the hypothesis that the autoregressive effect of either X or Y are equal across groups.
#' * `crossedXA = crossedXB` and `crossedYA = crossedYB`: Tests the hypothesis that the crossed effect of X on Y (`crossedX`) or of Y on X (`crossedY`), respectively, is equal across groups.
#' * `corBXBYA = corBXBYB`: Tests the hypothesis that the correlation between the random intercept factors is equal across groups.
#' 
#' For hypotheses regarding the traditional CLPM, see [semPower.powerCLPM()].
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model**:
#' * `Lambda`: The factor loading matrix (with the number of columns equaling 2 times the number of waves). Columns should be in order X1, Y1, X2, Y2, ..., X_nWaves, Y_nWaves.
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure ordered by wave, e.g., list(c(.2, .2, .2), c(.4, .4, .4, .4), c(.2, .2, .2), c(.4, .4, .4, .4), c(.2, .2, .2), c(.4, .4, .4, .4)) defines loadings of .2 for the three indicators of X at waves 1-3 and loadings of .4 for the four indicators of Y at waves 1-3. Must not contain secondary loadings.   
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators for each factor ordered by wave, e.g. c(3, 4, 3, 4, 3, 4) defines three indicators for X at waves 1-3 and four indicators for Y at waves 1-3. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor at each wave (if a vector is provided), e. g. `loadM = c(.5, .6, .5, .6, .5, .6)` defines mean loadings of .5 for X at waves 1-3 and mean loadings of .6 for Y at waves 1-3.
#' 
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` need to be defined.
#' If the model contains observed variables only, use `Lambda = diag(x)` where `x` is the number of variables.
#'
#' Note that the order of the factors is (X1, Y1, X2, Y2, ..., X_nWaves, Y_nWaves), i. e., the first factor is treated as the first measurement of X, the second as the first measurement of Y, the third as the second measurement of X, etc.. 
#' 
#' Additional arguments related to the requested type of **power analysis**:
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.  This requires the `semTools` package.
#' 
#' @examples
#' \dontrun{
#' # Determine required N in a 3-wave RI-CLPM
#' # to detect crossed effects of X (X1 -> Y2 and X2 -> Y3) of >= .2
#' # with a power of 95% on alpha = 5%, where
#' # X1, X2, and X3 are measured by 5 indicators loading by .5 each, and
#' # Y1, Y2, and Y3 are measured by 3 indicators loading by .4 each, and
#' # there is no synchronous correlation between X and Y (rXY = NULL),
#' # the correlation between the random intercept factors of X and Y (rBXBY) is .1,
#' # the autoregressive effects of X are .8 (equal across waves),
#' # the autoregressive effects of Y are .7 (equal across waves), and
#' # the crossed effects of Y (Y1 -> X2 and Y2 -> X3) are .1 (equal across waves).
#' 
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = .1,
#'                                     nullEffect = 'crossedX = 0',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     alpha = .05, beta = .05)
#' 
#' # show summary
#' summary(powerRICLPM)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerRICLPM$modelH1, sample.cov = powerRICLPM$Sigma,
#'             sample.nobs = powerRICLPM$requiredN, sample.cov.rescale = FALSE)
#' lavaan::sem(powerRICLPM$modelH0, sample.cov = powerRICLPM$Sigma,
#'             sample.nobs = powerRICLPM$requiredN, sample.cov.rescale = FALSE)
#' 
#' 
#' # same as above, but determine power with N = 500 on alpha = .05
#' powerRICLPM <- semPower.powerRICLPM(type = 'post-hoc',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = .1,
#'                                     nullEffect = 'crossedX = 0',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     alpha = .05, N = 500)
#' 
#' 
#' # same as above, but determine the critical chi-square with N = 500 so that alpha = beta
#' powerRICLPM <- semPower.powerRICLPM(type = 'compromise',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = .1,
#'                                     nullEffect = 'crossedX = 0',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     abratio = 1, N = 500)
#' 
#' 
#' # same as above, but compare to the saturated model
#' # (rather than to the less restricted model)
#' powerRICLPM <- semPower.powerRICLPM(type = 'compromise',
#'                                     comparison = 'saturated',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = .1,
#'                                     nullEffect = 'crossedX = 0',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     abratio = 1, N = 500)
#' 
#' 
#' # same as above, but assume only observed variables
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = .1,
#'                                     nullEffect = 'crossedX = 0',
#'                                     Lambda = diag(6),
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # same as above, but provide reduced loadings matrix to define that
#' # X1, X2, and X3 are measured by 5 indicators each loading by .5, .4, .5, .4, .3
#' # Y1, Y2, and Y3 are measured by 3 indicators each loading by .4, .3, .2
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = .1,
#'                                     nullEffect = 'crossedX = 0',
#'                                     loadings = list(
#'                                       c(.5, .4, .5, .4, .3),    # X1
#'                                       c(.4, .3, .2),            # Y1
#'                                       c(.5, .4, .5, .4, .3),    # X2
#'                                       c(.4, .3, .2),            # Y2
#'                                       c(.5, .4, .5, .4, .3),    # X3
#'                                       c(.4, .3, .2)             # Y3
#'                                     ),
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # same as above, but do not assume metric invariance across waves
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = .1,
#'                                     nullEffect = 'crossedX = 0',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     metricInvariance = FALSE,
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # same as above, but determine N to detect that the crossed effect of Y 
#' # (Y1 -> X2 and Y2 -> X3) is >= .1.
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = .1,
#'                                     nullEffect = 'crossedY = 0',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # same as above, but determine N to detect that the autoregressive effect 
#' # of X (X1 -> X2 and X2 -> X3) is >= .8.
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = .1,
#'                                     nullEffect = 'autoregX = 0',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # same as above, but determine N to detect that the autoregressive effect 
#' # of Y (Y1 -> Y2) is >= .7.
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = .1,
#'                                     nullEffect = 'autoregY = 0',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # same as above, but determine N to detect that
#' # the crossed effect of X (X1 -> Y2) of .2 differs from
#' # the crossed effect of Y (Y1 -> X2) of .1
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = .1,
#'                                     nullEffect = 'crossedX = crossedY',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # same as above, but determine N to detect that
#' # the autoregressive effect of X (X1 -> X2) of .8 differs from
#' # the autoregressive effect of Y (Y1 -> Y2) of .7
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = .1,
#'                                     nullEffect = 'autoregX = autoregY',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # same as above, but determine N to detect that the correlation between the 
#' # random intercept factors is >= .1
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = .1,
#'                                     nullEffect = 'corBXBY = 0',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # same as above, but assume that the synchronous (residual-)correlations between
#' #  X and Y are equal across waves, 
#' # namely a synchronous correlation of .05 at the first wave and residual correlations 
#' # of .05 at the second and third wave,
#' # and determine N to detect a crossed effect of X (X1 -> Y2 and X2 -> Y3) of >= .2
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY', 
#'                                                   'corXY'),
#'                                     rXY = c(.05, .05, .05),
#'                                     rBXBY = .1,
#'                                     nullEffect = 'crossedX = 0',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # same as above, but assume that the synchronous correlation between X and Y
#' # is .3 at the first wave, and the respective residual correlations are .2 at 
#' # the second wave and .3 at the third wave,
#' # and determine N to detect that the synchronous residual correlation at wave 2 is => .2.
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = c(.3, .2, .3),
#'                                     rBXBY = .1,
#'                                     nullEffect = 'corXY = 0',
#'                                     nullWhich = 2,
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # Determine required N in a 3-wave RI-CLPM to detect that
#' # the crossed effect of X at wave 1 (X1 -> Y2) of .20 is equal to the
#' # the crossed effect of X at wave 2 (X2 -> Y3) of .05
#' # with a power of 95% on alpha = 5%, where
#' # the autoregressive effects of X and Y are equal over waves,
#' # X1, X2, and X3 are measured by 5 indicators loading by .5 each, and
#' # Y1, Y2, and Y3 are measured by 3 indicators loading by .4 each, and
#' # the synchronous correlation between X and Y are .2, .3, and .4 at the first, 
#' # second, and third wave, 
#' # the correlation between the random intercept factors of X and Y is .1, and
#' # the autoregressive effect of X is .8 across all three waves,
#' # the autoregressive effect of Y is .7 across all three waves, and
#' # the crossed effects of Y (Y1 -> X2, and Y2 -> Y3) are both .1 
#' # (but freely estimated for each wave).
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = list(
#'                                       # X   Y
#'                                       c(.20, .10),  # wave 1 -> wave 2
#'                                       c(.05, .10)), # wave 2 -> wave 3
#'                                     waveEqual = c('autoregX', 'autoregY'),
#'                                     rXY = c(.2, .3, .4),
#'                                     rBXBY = .1,
#'                                     nullEffect = 'crossedX',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # same as above, but determine N to detect that
#' # the crossed effect of X at wave 2 is >= .05.
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = list(
#'                                       # X   Y
#'                                       c(.20, .10),  # wave 1 -> wave 2
#'                                       c(.05, .10)), # wave 2 -> wave 3
#'                                     waveEqual = c('autoregX', 'autoregY'),
#'                                     rXY = c(.2, .3, .4),
#'                                     rBXBY = .1,
#'                                     nullEffect = 'crossedX = 0',
#'                                     nullWhich = 2,
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     alpha = .05, beta = .05)
#' 
#' 
#' # same as above, but determine N to detect that
#' # the residual correlation between X and Y at wave 2 (of .3) differs from
#' # the residual correlation between X and Y at wave 3 (of .4).
#' powerRICLPM <- semPower.powerRICLPM(type = 'a-priori',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = list(
#'                                       # X   Y
#'                                       c(.20, .10),  # wave 1 -> wave 2
#'                                       c(.05, .10)), # wave 2 -> wave 3
#'                                     waveEqual = c('autoregX', 'autoregY'),
#'                                     rXY = c(.2, .3, .4),
#'                                     rBXBY = .1,
#'                                     nullEffect = 'corXY',
#'                                     nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                     loadM = c(.5, .4, .5, .4, .5, .4),
#'                                     alpha = .05, beta = .05)
#'
#'
#' # multigroup example
#' # Determine the achieved power N in a 3-wave RI-CLPM to detect that
#' # the crossed effect of X at wave 1 (X1 -> Y2) in group 1 of .25 differs
#' # from the crossed effect of X at wave 1 (X1 -> Y2) in group 2 of .15,
#' # where both groups comprise 500 observations and alpha = 5%, and
#' # the measurement model is equal for both groups, and
#' # the crossed effects of X (X1 -> Y2, and X2 -> Y3) are .25 and .10 in the first group, 
#' # the crossed effects of X (X1 -> Y2, and X2 -> Y3) are .15 and .05 in the second group, 
#' # the crossed effects of Y (Y1 -> X2, and Y2 -> X3) are .05 and .15 in the first group, 
#' # the crossed effects of Y (Y1 -> X2, and Y2 -> X3) are .01 and .10 in the second group, and
#' # the autoregressive effects of X (of .5) and Y (of .4) are equal over waves and over groups 
#' # (but freely estimated in each group).
#' powerRICLPM <- semPower.powerRICLPM(type = 'post-hoc', alpha = .05, N = list(500, 500),
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.5, .4), # group and wave constant 
#'                                     crossedEffects = list(
#'                                       # group 1
#'                                       list(
#'                                         c(.25, .10),   # X
#'                                         c(.05, .15)    # Y 
#'                                       ),
#'                                       # group 2
#'                                       list(
#'                                         c(.15, .05),   # X
#'                                         c(.01, .10)    # Y 
#'                                       )
#'                                     ),
#'                                     rXY = NULL,        # identity
#'                                     rBXBY = NULL,      # identity 
#'                                     nullEffect = 'crossedXA = crossedXB',
#'                                     nullWhich = 1,
#'                                     nIndicator = rep(3, 6), 
#'                                     loadM = c(.5, .6, .5, .6, .5, .6),
#'                                     metricInvariance = TRUE,
#'                                     waveEqual = c('autoregX', 'autoregY')
#'                                     )
#' 
#' 
#' # Request a simulated post-hoc power analysis with 500 replications
#' # to detect crossed effects of X (X1 -> Y2 and X2 -> Y3) of >= .2
#' # with a power of 95% on alpha = 5% in a RI-CLPM with 3 waves, 
#' # where there are only observed variables and 
#' # there is no synchronous correlation between X and Y (rXY = NULL),
#' # and no correlation between the random intercept factors of X and Y (rBXBY = NULL),
#' # the autoregressive effects of X are .8 (equal across waves),
#' # the autoregressive effects of Y are .7 (equal across waves), and
#' # the crossed effects of Y (Y1 -> X2 and Y2 -> X3) are .1 (equal across waves).
#' set.seed(300121)
#' powerRICLPM <- semPower.powerRICLPM(type = 'post-hoc',
#'                                     nWaves = 3,
#'                                     autoregEffects = c(.8, .7),
#'                                     crossedEffects = c(.2, .1),
#'                                     waveEqual = c('autoregX', 'autoregY', 
#'                                                   'crossedX', 'crossedY'),
#'                                     rXY = NULL,
#'                                     rBXBY = NULL,
#'                                     nullEffect = 'crossedX = 0',
#'                                     Lambda = diag(6),
#'                                     alpha = .05, N = 500,
#'                                     simulatedPower = TRUE, 
#'                                     simOptions = list(nReplications = 500))
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerRICLPM <- function(type, comparison = 'restricted',
                                 nWaves = NULL, 
                                 autoregEffects = NULL, 
                                 crossedEffects = NULL, 
                                 rXY = NULL,
                                 rBXBY = NULL,
                                 waveEqual = NULL, 
                                 nullEffect = NULL,
                                 nullWhichGroups = NULL,
                                 nullWhich = NULL,
                                 metricInvariance = TRUE,
                                 autocorResiduals = TRUE,
                                 ...){
  
  comparison <- checkComparisonModel(comparison)
  checkEllipsis(...)
  
  # validate input
  if('standardized' %in% names(list(...)) && list(...)[['standardized']]) stop('Standardized is not available for RICLPM.')
  if(is.null(autoregEffects) ||  is.null(crossedEffects)) stop('autoregEffects and crossedEffects may not be NULL.')
  if(is.null(nWaves) || is.na(nWaves) || nWaves < 3) stop('nWaves must be >= 3.')
  
  # we do not allow stacking of hypotheses. there might be a use case for this,
  # but this would complicate defining the relevant parameter when these vary across waves. 
  nullValid <- c('autoregx', 'autoregy', 'crossedx', 'crossedy', 'corxy',
                 'autoregx=0', 'autoregy=0', 'crossedx=0', 'crossedy=0',
                 'autoregx=autoregy', 'crossedx=crossedy', 'corxy=0', 'corbxby=0',
                 'autoregxa=autoregxb', 'autoregya=autoregyb', 
                 'crossedxa=crossedxb', 'crossedya=crossedyb',
                 'corbxbya=corbxbyb')
  nullEffect <- checkNullEffect(nullEffect, nullValid)
  
  # create list structure for autoregEffects, crossedEffects, and corXY
  # assume multigroup when list structure is present for either autoreg or crossed effects
  ngA <- ifelse(is.list(autoregEffects[[1]]), length(autoregEffects), 1)
  ngX <- ifelse(is.list(crossedEffects[[1]]), length(crossedEffects), 1)
  if(sum(c(ngA, ngX) > 1) > 1  && ngA != ngX) stop('Specify the same number of groups for both autoregEffects and crossedEffects.')
  nGroups <- max(c(ngA, ngX))
  isMultigroup <- nGroups > 1
  
  if(isMultigroup && !nullEffect %in% c('autoregxa=autoregxb', 'autoregya=autoregyb', 'crossedxa=crossedxb', 'crossedya=crossedyb','corbxbya=corbxbyb')) stop('Multigroup analysis are only supported for nullEffect = autoregxa=autoregxb, autoregya=autoregyb, crossedxa=crossedxb, crossedya=crossedyb, corbxbya=corbxbyb')
  if(!isMultigroup && nullEffect %in% c('autoregxa=autoregxb', 'autoregya=autoregyb', 'crossedxa=crossedxb', 'crossedya=crossedyb','corbxbya=corbxbyb')) stop('nullEffect = autoregxa=autoregxb, autoregya=autoregyb, crossedxa=crossedxb, crossedya=crossedyb, corbxbya=corbxbyb imply multigroup analyses, but no list structure for any relevant parameter provided.')
  if(isMultigroup && is.null(nullWhichGroups)) nullWhichGroups <- seq(nGroups)
  
  # [[groups]][[X, Y]][[waves]]
  if(!is.list(autoregEffects)) autoregEffects <- list(rep(autoregEffects[[1]], (nWaves - 1)), rep(autoregEffects[[2]], (nWaves - 1)))
  if(!is.list(crossedEffects)) crossedEffects <- list(rep(crossedEffects[[1]], (nWaves - 1)), rep(crossedEffects[[2]], (nWaves - 1)))
  if(!is.list(autoregEffects[[1]])) autoregEffects <- rep(list(autoregEffects), nGroups)
  if(!is.list(crossedEffects[[1]])) crossedEffects <- rep(list(crossedEffects), nGroups)
  if(length(autoregEffects[[1]][[1]]) == 1) autoregEffects <- lapply(autoregEffects, function(x) lapply(x, function(y) rep(y, nWaves - 1)))
  if(length(crossedEffects[[1]][[1]]) == 1) crossedEffects <- lapply(crossedEffects, function(x) lapply(x, function(y) rep(y, nWaves - 1)))
  
  if(any(unlist(lapply(autoregEffects, function(x) length(x) != 2)))) stop('Provide autoregEffects for X and Y.')
  if(any(unlist(lapply(crossedEffects, function(x) length(x) != 2)))) stop('Provide crossedEffects for X and Y..')
  if(any(unlist(lapply(autoregEffects, function(x) length(x[[1]]) != length(x[[2]]))))) stop('autoregEffects for X and Y must be of equal length.')
  if(any(unlist(lapply(crossedEffects, function(x) length(x[[1]]) != length(x[[2]]))))) stop('crossedEffects for X and Y must be of equal length.')
  if(any(unlist(lapply(autoregEffects, function(x) length(x[[1]]) != (nWaves - 1))))) stop('autoregEffects must be of length nWaves - 1.')
  if(any(unlist(lapply(crossedEffects, function(x) length(x[[1]]) != (nWaves - 1))))) stop('crossedEffects must be of length nWaves - 1.')
  invisible(lapply(autoregEffects, function(y) lapply(y, function(x) lapply(x, function(x) checkBounded(x, 'All autoregressive effects ', bound = c(-1, 1), inclusive = FALSE)))))
  invisible(lapply(crossedEffects, function(y) lapply(y, function(x) lapply(x, function(x) checkBounded(x, 'All crossed effects ', bound = c(-1, 1), inclusive = FALSE)))))
  
  if(is.null(rXY)) rXY <- rep(0, nWaves)
  if(is.list(rXY) && length(rXY) != nGroups) stop('corXY implies a different number of groups as autoregEffects or crossedEffects.')
  if(!is.list(rXY)) rXY <- rep(list(rXY), nGroups)
  if(any(unlist(lapply(rXY, function(x) length(x) != nWaves)))) stop('rXY must be of length nWaves')
  invisible(lapply(rXY, function(y) lapply(y, function(x) checkBounded(x, 'All rXY ', bound = c(-1, 1), inclusive = FALSE))))
  
  if(nullEffect %in% c('corbxbya=corbxbyb') && (is.null(rBXBY) || length(rBXBY) == 1)) stop('nullEffect corbxbya=corbxbyb requires that rBXBY is specified for each group.')
  if(is.null(rBXBY)) rBXBY <- rep(list(0), nGroups)
  if(isMultigroup && length(rBXBY) == 1) rBXBY <- rep(list(rBXBY), nGroups)
  if(any(unlist(lapply(rBXBY, function(x) length(x) != 1)))) stop('rBXBY must contain a single number or be a list of single numbers')
  invisible(lapply(rBXBY, function(x) checkBounded(x, 'All rBXBY ', bound = c(-1, 1), inclusive = FALSE)))
  
  if(!is.null(waveEqual)){
    waveEqual <- unlist(lapply(waveEqual, function(x) tolower(trimws(x))))
    if(any(unlist(lapply(waveEqual, function(x) !x %in% c('autoregx', 'autoregy', 'crossedx', 'crossedy', 'corxy'))))) stop('waveEqual may only contain autoregX, autoregY, crossedX, crossedY, corXY')
  }
  
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
  Lambda  <- args[['Lambda']]
  if(is.null(Lambda)){
    Lambda <- genLambda(args[['loadings']], args[['nIndicator']],
                        args[['loadM']], args[['loadSD']], args[['loadMinMax']])
  }
  if(ncol(Lambda) != 2*nWaves) stop('Number of factors must be 2*nWaves.')
  
  # modify Lambda according to RI-CLPM structure
  Lambda <- cbind(matrix(0, nrow = nrow(Lambda), ncol = (2*nWaves + 2)), Lambda) # add between + within factors
  # cols: Bx, By, Wx_1, Wy_1,..., Wx_nWaves, Wy_nWaves, Fx_1, Fy_1, ..., Fx_nWaves, Fy_nWaves
  
  Lambda <- rep(list(Lambda), nGroups)  # require same measurement model across groups
  
  ### create Beta
  Beta <- lapply(seq(nGroups), function(g){
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
      B[xidx, (xidx - 2)] <- autoregEffects[[g]][[1]][i] 
      B[yidx, (yidx - 2)] <- autoregEffects[[g]][[2]][i]
      # crossed effects
      B[yidx, (xidx - 2)] <- crossedEffects[[g]][[1]][i]
      B[xidx, (yidx - 2)] <- crossedEffects[[g]][[2]][i]
    }
    
    B
  })
  
  ### create Psi
  Psi <- lapply(seq(nGroups), function(g){
    P <- diag(ncol(Beta[[1]]))
    # Bx, By, Wx_1, Wy_1,..., Wx_nWaves, Wy_nWaves, Fx_1, Fy_1, ..., Fx_nWaves, Fy_nWaves
    
    # add cor between random intercepts
    P[2,1] <- P[1,2] <- rBXBY[[g]]
    
    # set residual variance of Fx_1, ..., Fy_nWaves to 0
    diag(P[(3 + 2*nWaves):(4*nWaves + 2), (3 + 2*nWaves):(4*nWaves + 2)]) <- 0 
    
    # add (residual) correlations between within-factors
    if(any(rXY[[g]] != 0)){
      for(i in 1:nWaves){
        P[(2*i + 2), (2*i + 1)] <- P[(2*i + 1), (2*i + 2)] <- rXY[[g]][i]
      }
    }
    
    P
  })
  
  # add metric invariance constrains
  metricInvarianceFactors <- NULL
  if(metricInvariance){
    metricInvarianceFactors <- list(
      seq(3 + 2*nWaves, 2 + 4*nWaves, 2),
      seq(4 + 2*nWaves, 2 + 4*nWaves, 2)  
    )
  }
  
  ### get model-implied sigma
  generated <- semPower.genSigma(Beta = if(!isMultigroup) Beta[[1]] else Beta, 
                                 Psi = if(!isMultigroup) Psi[[1]] else Psi, 
                                 Lambda = if(!isMultigroup) Lambda[[1]] else Lambda, 
                                 useReferenceIndicator = TRUE,
                                 metricInvariance = metricInvarianceFactors,
                                 nGroups = nGroups)
  
  
  ### create ana model string
  if(!isMultigroup) modelCFA <- generated[['modelTrueCFA']] else modelCFA <- generated[[1]][['modelTrueCFA']]
  
  # define random intercept (between) factors
  tok1 <- paste0('f1 =~ ', paste0('1*', paste0('f', seq((3 + 2*nWaves), ncol(Beta[[1]]), 2)), collapse = ' + '))
  tok2 <- paste0('f2 =~ ', paste0('1*', paste0('f', seq((4 + 2*nWaves), ncol(Beta[[1]]), 2)), collapse = ' + '))
  model <- paste(tok1, tok2, sep='\n')
  
  # define residualized (within) factors
  for(i in 1:(2*nWaves)){
    widx <- seq(3, 2*nWaves + 2)[i]
    fidx <- seq(2*nWaves + 3, 4*nWaves + 2)[i]
    model <- paste(model, paste0('f', widx, ' =~ 1*f', fidx), sep = '\n')
  }
  
  # add unresidualized factors
  model <- paste(model, modelCFA, sep='\n')
  
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
    fidx <- which(Beta[[1]][f, ] != 0)
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
  
  
  # add autocorrelated residuals
  if(autocorResiduals){
    if(!isMultigroup) Lambda <- generated[['Lambda']] else Lambda <- generated[[1]][['Lambda']]
    # do this only when there is at least one latent variable
    if(nrow(Lambda) > 2*nWaves){
      autocorResidualsFactors <- metricInvarianceFactors  # same structure 
      tok <- list()
      for(x in seq_along(autocorResidualsFactors)){
        ci <- lapply(autocorResidualsFactors[[x]], function(f) paste0('x', which(Lambda[, f] != 0)))
        if(length(ci[[1]]) > 1){
          for(i in 1:(length(ci) - 1)){
            for(j in (i + 1) : length(ci)){
              tok <- append(tok, paste(ci[[i]], '~~', ci[[j]]))
            }
          }
        }
      }
      model <- paste(model, paste(tok, collapse = '\n'), sep ='\n')
    }
  }
  
  
  ### define H1 and H0 model
  
  # first get relevant parameter labels that may be subject to waveEqual or nullEffect constraints 
  xw <- seq(2 + 2*nWaves - 1, 5, -2)
  pAutoregX <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(xw - 2, width = 2, flag = 0))
  pAutoregX <- pAutoregX[order(pAutoregX)]
  
  yw <- seq(2 + 2*nWaves, 6, -2)
  pAutoregY <- paste0('pf', formatC(yw, width = 2, flag = 0), formatC(yw - 2, width = 2, flag = 0))
  pAutoregY <- pAutoregY[order(pAutoregY)]
  
  xw <- seq(2 + 2*nWaves - 3, 3, -2)
  yw <- seq(2 + 2*nWaves, 6, -2)
  pCrossedX <- paste0('pf', formatC(yw, width = 2, flag = 0), formatC(xw, width = 2, flag = 0))
  pCrossedX <- pCrossedX[order(pCrossedX)]
  
  xw <- seq(2 + 2*nWaves - 1, 5, -2)
  yw <- seq(2 + 2*nWaves - 2, 4, -2)
  pCrossedY <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(yw, width = 2, flag = 0))
  pCrossedY <- pCrossedY[order(pCrossedY)]
  
  xw <- seq(2+ 2*nWaves - 1, 5, -2)
  yw <- seq(2 + 2*nWaves, 6, -2)
  pCorXY <- paste0('pf', formatC(yw, width = 2, flag = 0), formatC(xw, width = 2, flag = 0))
  pCorXY <- pCorXY[order(pCorXY)]
  
  # multigroup case
  if(isMultigroup){
    # remove group specific labels from measurement part to enforce metric invariance
    model <- gsub(paste0('_g', seq(nGroups), collapse = '|'), '_gc', model)
    # assign group labels to all structural parameters (measurement part is held equal across groups)
    patt <- 'pf0201' #  RI cor
    repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
    model <- gsub(patt, repl, model)
    patt <- 'pf0403' # exog x,y cor
    repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
    model <- gsub(patt, repl, model)
    for(pp in seq(nWaves - 1)){
      patt <- pAutoregX[pp]
      repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      model <- gsub(patt, repl, model)
      patt <- pAutoregY[pp]
      repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      model <- gsub(patt, repl, model)
      patt <- pCrossedX[pp]
      repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      model <- gsub(patt, repl, model)
      patt <- pCrossedY[pp]
      repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      model <- gsub(patt, repl, model)
      patt <- pCorXY[pp]
      repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      model <- gsub(patt, repl, model)
    }
  }
  
  # add constraints to H1 model
  modelH1 <- model
  if(!is.null(waveEqual)){
    if('autoregx' %in% waveEqual) modelH1 <- gsub(paste(pAutoregX, collapse = '|'), paste(pAutoregX, collapse = ''), modelH1)
    if('autoregy' %in% waveEqual) modelH1 <- gsub(paste(pAutoregY, collapse = '|'), paste(pAutoregY, collapse = ''), modelH1)
    if('crossedx' %in% waveEqual) modelH1 <- gsub(paste(pCrossedX, collapse = '|'), paste(pCrossedX, collapse = ''), modelH1)
    if('crossedy' %in% waveEqual) modelH1 <- gsub(paste(pCrossedY, collapse = '|'), paste(pCrossedY, collapse = ''), modelH1)
    if('corxy' %in% waveEqual) modelH1 <- gsub(paste(pCorXY, collapse = '|'), paste(pCorXY, collapse = ''), modelH1)
  }
  
  # add additional constraints to H0 model
  modelH0 <- modelH1
  
  # wave equal constraints not included in modelH1:
  if('autoregx' %in% nullEffect) modelH0 <- gsub(paste(pAutoregX, collapse = '|'), paste(pAutoregX, collapse = ''), modelH0)
  if('autoregy' %in% nullEffect) modelH0 <- gsub(paste(pAutoregY, collapse = '|'), paste(pAutoregY, collapse = ''), modelH0)
  if('crossedx' %in% nullEffect) modelH0 <- gsub(paste(pCrossedX, collapse = '|'), paste(pCrossedX, collapse = ''), modelH0)
  if('crossedy' %in% nullEffect) modelH0 <- gsub(paste(pCrossedY, collapse = '|'), paste(pCrossedY, collapse = ''), modelH0)
  if('corxy' %in% nullEffect) modelH0 <- gsub(paste(pCorXY, collapse = '|'), paste(pCorXY, collapse = ''), modelH0)
  
  # zero and equality constraints:
  if('autoregx=0' %in% nullEffect){
    if('autoregx' %in% waveEqual){
      modelH0 <- gsub(paste(pAutoregX, collapse = ''), '0', modelH0)
    }else{
      modelH0 <- gsub(pAutoregX[nullWhich], '0', modelH0)
    }
  }
  if('autoregy=0' %in% nullEffect){
    if('autoregy' %in% waveEqual){
      modelH0 <- gsub(paste(pAutoregY, collapse = ''), '0', modelH0)
    }else{
      modelH0 <- gsub(pAutoregY[nullWhich], '0', modelH0)
    }
  }
  if('crossedx=0' %in% nullEffect){
    if('crossedx' %in% waveEqual){
      modelH0 <- gsub(paste(pCrossedX, collapse = ''), '0', modelH0)
    }else{
      modelH0 <- gsub(pCrossedX[nullWhich], '0', modelH0)
    }
  }
  if('crossedy=0' %in% nullEffect){
    if('crossedy' %in% waveEqual){
      modelH0 <- gsub(paste(pCrossedY, collapse = ''), '0', modelH0)
    }else{
      modelH0 <- gsub(pCrossedY[nullWhich], '0', modelH0)
    }
  }
  if('autoregx=autoregy' %in% nullEffect){
    if('autoregx' %in% waveEqual && 'autoregy' %in% waveEqual){
      patt <- paste(c(paste(pAutoregX, collapse = ''), paste(pAutoregY, collapse = '')), collapse = '|')
    }else if('autoregx' %in% waveEqual){
      patt <- paste(c(paste(pAutoregX, collapse = ''), pAutoregY[nullWhich]), collapse = '|')
    }else if('autoregy' %in% waveEqual){
      patt <- paste(c(pAutoregX[nullWhich], paste(pAutoregY, collapse = '')), collapse = '|')
    }else{
      patt <- paste(c(pAutoregX[nullWhich], pAutoregY[nullWhich]), collapse = '|')
    }
    repl <-  gsub('\\|', '', patt)
    modelH0 <- gsub(patt, repl, modelH0)
  }
  if('crossedx=crossedy' %in% nullEffect){
    if('crossedx' %in% waveEqual && 'crossedy' %in% waveEqual){
      patt <- paste(c(paste(pCrossedX, collapse = ''), paste(pCrossedY, collapse = '')), collapse = '|')
    }else if('crossedx' %in% waveEqual){
      patt <- paste(c(paste(pCrossedX, collapse = ''), pCrossedY[nullWhich]), collapse = '|')
    }else if('crossedy' %in% waveEqual){
      patt <- paste(c(pCrossedX[nullWhich], paste(pCrossedY, collapse = '')), collapse = '|')
    }else{
      patt <- paste(c(pCrossedX[nullWhich], pCrossedY[nullWhich]), collapse = '|')
    }
    repl <-  gsub('\\|', '', patt)
    modelH0 <- gsub(patt, repl, modelH0)
  }
  if('corxy=0' %in% nullEffect){
    if('corxy' %in% waveEqual){
      modelH0 <- gsub(paste(pCorXY, collapse = ''), '0', modelH0)
    }else{
      pCorXY <- c('pf0403', pCorXY)   # add exog cor
      modelH0 <- gsub(pCorXY[nullWhich], '0', modelH0)
    }
  }
  if('corbxby=0' %in% nullEffect){
    modelH0 <- gsub('pf0201', '0', modelH0)
  }
  
  # multigroup cases
  if('autoregxa=autoregxb' %in% nullEffect){
    if('autoregx' %in% waveEqual){
      patt <- paste0(paste(pAutoregX, collapse = ''), '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(paste(pAutoregX, collapse = ''), '_gc')
    }else{
      patt <- paste0(pAutoregX[nullWhich], '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(pAutoregX[nullWhich], '_gc')
    }
    modelH0 <- gsub(patt, repl, modelH0)
  }
  if('autoregya=autoregyb' %in% nullEffect){
    if('autoregy' %in% waveEqual){
      patt <- paste0(paste(pAutoregY, collapse = ''), '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(paste(pAutoregY, collapse = ''), '_gc')
    }else{
      patt <- paste0(pAutoregY[nullWhich], '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(pAutoregY[nullWhich], '_gc')
    }
    modelH0 <- gsub(patt, repl, modelH0)
  }
  if('crossedxa=crossedxb' %in% nullEffect){
    if('crossedx' %in% waveEqual){
      patt <- paste0(paste(pCrossedX, collapse = ''), '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(paste(pCrossedX, collapse = ''), '_gc')
    }else{
      patt <- paste0(pCrossedX[nullWhich], '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(pCrossedX[nullWhich], '_gc')
    }
    modelH0 <- gsub(patt, repl, modelH0)
  }
  if('crossedya=crossedyb' %in% nullEffect){
    if('crossedy' %in% waveEqual){
      patt <- paste0(paste(pCrossedY, collapse = ''), '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(paste(pCrossedY, collapse = ''), '_gc')
    }else{
      patt <- paste0(pCrossedY[nullWhich], '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(pCrossedY[nullWhich], '_gc')
    }
    modelH0 <- gsub(patt, repl, modelH0)
  }  
  if('corbxbya=corbxbyb' %in% nullEffect){
    patt <- paste(paste0('pf0201_g', nullWhichGroups), collapse = '|')
    modelH0 <- gsub(patt, 'pf0201_gc', modelH0)
  }
  
  
  # here we actually fit modelH1 in case of a restricted comparison
  # because we cannot be sure that user input yields perfectly fitting h1 models 
  # when there are additional constraints (waveequal or invariance)
  if(comparison == 'saturated') modelH1 <- NULL
  
  if(isMultigroup) Sigma <- lapply(generated, '[[', 'Sigma') else Sigma <- generated[['Sigma']] 
  
  semPower.powerLav(type, 
                    modelH0 = modelH0, 
                    modelH1 = modelH1, 
                    Sigma = Sigma,
                    ...)
}

#' semPower.powerMI
#'
#' Convenience function for performing power analyses for hypothesis arising 
#' in multigroup measurement invariance models concerning a specific level of invariance.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, either `'saturated'` or one of `'configural'`, `'metric'`, `'scalar'`, `'covariances'`, or a vector of restrictions in `lavaan` format (with `'none'` for no restrictions). See details.
#' @param nullEffect defines the hypothesis (i.e., level of invariance) of interest. One of `'metric'`, `'scalar'`, `'residual'`, `'covariances'`, `'means'` or a vector of restrictions in `lavaan` format. See details.   
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details
#' This function performs a power analysis to reject various hypotheses arising
#' in the context of multigroup measurement invariance. Multigroup invariance models 
#' fit the specified model simultaneously to various groups and place increasingly
#' restrictive cross-group equality constrains on the model parameters. The typical - but not in all parts necessary -
#' sequence is (a) configural, (b) metric, (c) scalar, and (d) residual invariance, where each level of invariance is
#' compared against the previous level (e.g., scalar vs. metric). Power analysis provides  
#' the power (or the required N) to reject a particular level of invariance.
#' 
#' For hypotheses regarding longitudinal invariance, see [semPower.powerLI()].
#'  
#' The models defined in the `comparison` and the `nullEffect` arguments can be specified in two ways. Either specify
#' a specific level of invariance that includes all previous levels:
#' \itemize{
#' \item `'configural'`: no invariance constraints. Shows the same fit as the saturated model, so only the delta df differ. 
#' \item `'metric'`: all loadings are restricted to equality. 
#' \item `'scalar'`: all loadings and (indicator-)intercepts are restricted to equality. 
#' \item `'residual'`: all loadings, (indicator-)intercepts, and (indicator-)residuals are restricted to equality.
#' \item `'covariances'`: all loadings, (indicator-)intercepts, and (indicator-)residuals, and latent covariances are restricted to equality.
#' \item `'means'`: all loadings, (indicator-)intercepts, (indicator-)residuals, latent covariances, and latent means are restricted to equality.
#' }
#' 
#' For example, setting `comparison = 'metric'` and `nullEffect = 'scalar'` determines power 
#' to reject the hypothesis that the constraints placed in the scalar invariance model 
#' (restricting loadings and intercepts) over the 
#' metric invariance model (restricting only the loadings) are defensible.
#'  
#' For greater flexibility, the models can also be defined using `lavaan` style `group.equal` restrictions as a vector: 
#' \itemize{
#' \item `'none'`: no invariance constraints and thus representing a configural invariance model. Shows the same fit as the saturated model, so only the delta df differ. 
#' \item `c('loadings')`: all loadings are restricted to equality. 
#' \item `c('loadings', 'intercepts')`: all loadings and (indicator-)intercepts are restricted to equality. 
#' \item `c('loadings', 'intercepts', 'residuals')`: all loadings, (indicator-)intercepts, and (indicator-)residuals are restricted to equality.
#' \item `c('loadings', 'residuals')`: all loadings and (indicator-)residuals are restricted to equality.
#' \item `c('loadings', 'intercepts', 'means')`: all loadings, (indicator-)intercepts, and latent factor means are restricted to equality.
#' }
#' 
#' For example, setting `comparison = c('loadings')` and `nullEffect = 'c('loadings', 'intercepts')'` 
#' determines power to reject the hypothesis that the constraints placed in the scalar invariance model 
#' (restricting loadings and intercepts) over the  metric invariance model (restricting only the loadings) are defensible.
#' Note that variance scaling is used, so invariance of variances (`'lv.variances'`) is always met. 
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model**:
#' * `Lambda`: The factor loading matrix (with the number of columns equaling the number of factors).
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading on the second factor by .8, .6, .6, and .4.
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the loadings of the first factor to equal .5 and those of the second factor do equal .6.
#' * `Theta`: Variance-covariance matrix of the indicator residuals, which should be a diagonal matrix. Required when residual non-invariance is to be detected. When `NULL`, Theta is a diagonal matrix with elements such that all variances are 1. 
#' * `tau`: Defines the item intercepts, required whenever a model involves hypotheses about means (e.g., scalar invariance). If `NULL` and `Alpha` is set, all intercepts are assumed to equal zero.
#' * `Alpha`: Defines the latent means, required whenever a model involves hypotheses about latent means (e.g., latent mean invariance). If `NULL` and `tau` is set, all latent means are assumed to equal zero. Because variance scaling is used so that all factor variances are 1, latent mean differences can be interpreted akin to Cohen's d as standardized mean differences.
#' 
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` always need to be defined, 
#' and `Theta`, `tau` and `Alpha` need to be defined for particular levels of invariance. 
#' As this function operates on multiple groups, either argument is a list whenever there are 
#' group differences in the respective parameters. When no list is provided, the same 
#' parameter values are assumed for all groups.
#' 
#' Additional arguments related to the requested type of **power analysis**:
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.  This requires the `semTools` package.
#' 
#' @examples
#' \dontrun{
#' # obtain the required N to reject the hypothesis of metric invariance
#' # in comparison to the configural invariance model 
#' # with a power of 95% on alpha = 5% 
#' # assuming equally sized groups (N = list(1, 1)) 
#' # for a factor model involving a single factor which 
#' # is measured by 5 indicators (in both groups)
#' # loading by .5 each in the first group and 
#' # loading by .6 each in the second group.
#' powerMI <- semPower.powerMI(type = 'a-priori',
#'                             comparison = 'configural', 
#'                             nullEffect = 'metric',
#'                             nIndicator = list(5, 5),
#'                             loadM = list(.5, .6),
#'                             alpha = .05, beta = .05, N = list(1, 1))
#' 
#' # show summary
#' summary(powerMI)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerMI$modelH1, sample.cov = list(powerMI$Sigma[[1]], powerMI$Sigma[[2]]),
#'             sample.nobs = as.list(powerMI$requiredN.g), sample.cov.rescale = FALSE)
#' lavaan::sem(powerMI$modelH0, sample.cov = list(powerMI$Sigma[[1]], powerMI$Sigma[[2]]),
#'             sample.nobs = as.list(powerMI$requiredN.g), sample.cov.rescale = FALSE)
#' 
#' # same as above, but determine power with N = 500 in each group on alpha = .05
#' powerMI <- semPower.powerMI(type = 'post-hoc',
#'                             comparison = 'configural', 
#'                             nullEffect = 'metric',
#'                             nIndicator = 5,
#'                             loadM = list(.5, .6),
#'                             alpha = .05, N = list(500, 500))
#' 
#' # same as above, but determine the critical chi-square with N = 500 in each 
#' # group so that alpha = beta
#' powerMI <- semPower.powerMI(type = 'compromise',
#'                             comparison = 'configural', 
#'                             nullEffect = 'metric',
#'                             nIndicator = 5,
#'                             loadM = list(.5, .6),
#'                             abratio = 1, N = list(500, 500))
#' 
#' # same as above, but compare to the saturated model
#' # (rather than to the configural invariance model)
#' powerMI <- semPower.powerMI(type = 'a-priori',
#'                             comparison = 'saturated', 
#'                             nullEffect = 'metric',
#'                             nIndicator = 5,
#'                             loadM = list(.5, .6),
#'                             alpha = .05, beta = .05, N = list(1, 1))
#' 
#' # same as above, but provide individual factor loadings by group using a 
#' # reduced loading matrix to define a  single factor model with three indicators
#' # loading by .4, .6, .5 in the first group and 
#' # loading by .5, .6, .7 in the second group
#' powerMI <- semPower.powerMI(type = 'a-priori',
#'                             comparison = 'saturated', 
#'                             nullEffect = 'metric',
#'                             loadings = list(
#'                               list(c(.4, .6, .5)), 
#'                               list(c(.5, .6, .7))),
#'                             alpha = .05, beta = .05, N = list(1, 1))
#' 
#' # same as above, but make first group twice as large as the second group 
#' powerMI <- semPower.powerMI(type = 'a-priori',
#'                             comparison = 'saturated', 
#'                             nullEffect = 'metric',
#'                             loadings = list(
#'                               list(c(.4, .6, .5)), 
#'                               list(c(.5, .6, .7))),
#'                             alpha = .05, beta = .05, N = list(2, 1))
#' 
#' # obtain the required N to reject the hypothesis of scalar invariance
#' # in comparison to the metric invariance model 
#' # with a power of 95% on alpha = 5% 
#' # assuming equally sized groups (N = list(1, 1)) 
#' # for a two factor model, where both factors are  
#' # measured by 3 indicators each and all loadings equal .5 (in both groups),
#' # the factor correlation is .3 in both groups, and the
#' # all intercepts are 0.0 in the first group, but
#' # all intercepts are 0.1 in the second group
#' powerMI <- semPower.powerMI(type = 'a-priori',
#'                             comparison = 'metric', 
#'                             nullEffect = 'scalar',
#'                             Phi = list(.3, .3),
#'                             nIndicator = list(
#'                               c(3, 3), 
#'                               c(3, 3)),
#'                             loadM = .5,
#'                             tau = list(
#'                               rep(0.0, 6), 
#'                               rep(0.1, 6) 
#'                             ),
#'                             alpha = .05, beta = .05, N = list(1, 1))
#' 
#' # same as above, but use lavaan group.equal strings 
#' powerMI <- semPower.powerMI(type = 'a-priori',
#'                             comparison = c('loadings'), 
#'                             nullEffect = c('loadings', 'intercepts'),
#'                             Phi = list(.3, .3),
#'                             nIndicator = list(
#'                               c(3, 3), 
#'                               c(3, 3)),
#'                             loadM = .5,
#'                             tau = list(
#'                               rep(0.0, 6), 
#'                               rep(0.1, 6) 
#'                             ),
#'                             alpha = .05, beta = .05, N = list(1, 1))
#' 
#' # same as above, but
#' # obtain the required N to reject the hypothesis of equal latent means
#' # in comparison to the scalar invariance model;
#' # all intercepts are zero in both groups, 
#' # in the first group, the latent means equal 0.0, 
#' # in the second group, the latent mean of the factors are 0.0 and 0.5
#' powerMI <- semPower.powerMI(type = 'a-priori',
#'                             comparison = c('loadings', 'intercepts'), 
#'                             nullEffect = c('loadings', 'intercepts', 'means'),
#'                             Phi = list(.3, .3),
#'                             nIndicator = list(
#'                               c(3, 3), 
#'                               c(3, 3)),
#'                             loadM = .5,
#'                             tau = list(
#'                               rep(0.0, 6), 
#'                               rep(0.0, 6) 
#'                             ),
#'                             Alpha = list(
#'                               c(0.0, 0.0),
#'                               c(0.0, 0.5)
#'                             ),
#'                             alpha = .05, beta = .05, N = list(1, 1))
#' 
#' # request a simulated post-hoc power analysis with 500 replications
#' # to reject the hypothesis of metric invariance.
#' set.seed(300121)
#' powerMI <- semPower.powerMI(type = 'post-hoc',
#'                             comparison = 'configural', 
#'                             nullEffect = 'metric',
#'                             nIndicator = list(5, 5),
#'                             loadM = list(.5, .6),
#'                             alpha = .05, N = list(500, 500), 
#'                             simulatedPower = TRUE, 
#'                             simOptions = list(nReplications = 500))
#'                              
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerMI <- function(type, 
                             comparison = NULL,
                             nullEffect = NULL,
                             ...){
  
  args <- list(...)
  
  # validate input
  checkEllipsis(...)
  lavGroupStrings <- c('loadings', 'intercepts', 'residuals', 'residual.covariances', 'lv.covariances','regressions')
  useLavOptions <- any(grepl(paste(lavGroupStrings, collapse = '|'), comparison)) || any(grepl(paste(lavGroupStrings, collapse = '|'), nullEffect))
  # we only check typos etc when not using lavstrings
  if(!useLavOptions){
    comparison <- checkNullEffect(comparison, c('saturated', 'configural', 'metric', 'scalar', 'covariances'))
    nullEffect <- checkNullEffect(nullEffect, c('metric', 'scalar', 'residual', 'covariances', 'means'))
    if(which(c('saturated', 'configural', 'metric', 'scalar', 'covariances') %in% comparison) >= 
       (2 + which(c('metric', 'scalar', 'residuals', 'covariances', 'means') %in% nullEffect))) stop('Model defined in nullEffect is not nested in comparison model.')
  }else{
    if('lv.variances' %in% comparison || 'lv.variances' %in% nullEffect) stop('Variance scaling is used, so invariance of latent variances is always met.')
    if(!any(c('saturated', 'none') %in% comparison) && !all(comparison %in% nullEffect)) stop('Comparison model is not nested in hypothesized model; all restrictions in comparison must also be present in nullEffect.')
  }
  
  ### generate sigmas
  # we use variance scaling, so the first loading may also differ across groups.
  # If using reference indicators instead, the first loading must be equal across groups.
  # Not sure whether to expose this to users.
  # the downside is that lv.variances is always true.
  generated <- semPower.genSigma(..., useReferenceIndicator = FALSE)   
  
  # more input validations
  if(!is.list(generated[[1]])) stop('Loadings, Phi, Beta, etc. must be provided as a list for each group.')
  if(is.null(generated[[1]][['mu']])){
    inv <- FALSE
    if(useLavOptions){
      inv <- any(grepl('intercepts|means', comparison)) || any(grepl('intercepts|means', nullEffect))
    }else{
      inv <- any(c('scalar', 'residuals', 'covariances', 'means') %in% nullEffect)
    }
    if(inv) stop('The models imply a meanstructure, so tau and/or Alpha need to be defined.')
  }
  
  # models are the same, the only difference pertains to lavOptions
  modelH0 <- modelH1 <- generated[[1]][['modelTrueCFA']]

  # set proper lavOptions
  lavOptionsH1 <- list()
  if(!useLavOptions){
    lavOptionsH0 <- list(group.equal = switch(nullEffect,
                                              'metric' = c('loadings'),
                                              'scalar' = c('loadings', 'intercepts'),
                                              'residual' = c('loadings', 'intercepts', 'residuals'),
                                              'covariances' = c('loadings', 'intercepts', 'residuals', 'lv.covariances'),
                                              'means' = c('loadings', 'intercepts', 'residuals', 'lv.covariances', 'means')
                                              
    ))
    if(comparison %in% c('metric', 'scalar', 'residuals', 'covariances')){
      lavOptionsH1 <- list(group.equal = switch(comparison,
                                                'metric' = c('loadings'),
                                                'scalar' = c('loadings', 'intercepts'),
                                                'residual' = c('loadings', 'intercepts', 'residuals'),
                                                'covariances' = c('loadings', 'intercepts', 'residuals', 'lv.covariances')
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
    if(any(c('scalar', 'residuals', 'covariances', 'means') %in% nullEffect))
      mu <- lapply(generated, '[[', 'mu')
  }else{
    if(any(grepl('intercepts|means', comparison)) || any(grepl('intercepts|means', nullEffect)))
      mu <- lapply(generated, '[[', 'mu')
  }
  
  args[['lavOptions']] <- append(args[['lavOptions']], lavOptionsH0)
  args[['lavOptionsH1']] <- append(args[['lavOptionsH1']], lavOptionsH1)
  
  do.call(semPower.powerLav, append(list(
    type = type,
    Sigma = Sigma,
    mu = mu,
    modelH0 = modelH0,
    modelH1 = modelH1,
    fitH1model = TRUE),
    args)
  )

}

#' semPower.powerLI
#'
#' Convenience function for performing power analyses for hypothesis arising 
#' in longitudinal measurement invariance models concerning a specific level of invariance.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, either `'saturated'` or one of `'configural'`, `'metric'`, `'scalar'`, `'residual'`, `'covariances'`, `'means'`, or a vector of restrictions in `lavaan` format (with `'none'` for no restrictions). See details.
#' @param nullEffect defines the hypothesis (i.e., level of invariance) of interest. Accepts the same arguments as `comparison`. See details.   
#' @param autocorResiduals whether the residuals of the indicators of latent variables are autocorrelated over waves (`TRUE`, the default) or not (`FALSE`). This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model).
#' @param Phi the factor correlation matrix. Can be `NULL` for uncorrelated factors.
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details
#' This function performs a power analysis to reject various hypotheses arising
#' in the context of longitudinal measurement invariance, where a single attribute is measured repeatedly. The typical - but not in all parts necessary -
#' sequence concerning the measurement part is (a) configural, (b) metric, (c) scalar, (d) residual invariance, 
#' and concerning the structural part  (e) latent covariances, (f) latent means, where each level of invariance is
#' compared against the previous level (e.g., scalar vs. metric). Power analysis provides  
#' the power (or the required N) to reject a particular level of invariance.
#' 
#' For hypotheses regarding multiple group invariance, see [semPower.powerMI()]. For hypotheses regarding autoregressive models, see [semPower.powerAutoregressive()]. For hypotheses in an ARMA model, see [semPower.powerARMA()].
#'  
#' There are two ways to specify the models defined in the `comparison` and the `nullEffect` arguments. Either, one may
#' specify a specific level of invariance that includes all previous levels:
#' \itemize{
#' \item `'configural'`: no invariance constraints. Shows the same fit as the saturated model, so only the delta df differ. 
#' \item `'metric'`: all loadings are restricted to equality over measurement occasions. Note that reference scaling is used, so the first indicator should be invariant.
#' \item `'scalar'`: all loadings and (indicator-)intercepts are restricted to equality. 
#' \item `'residual'`: all loadings, (indicator-)intercepts, and (indicator-)residuals are restricted to equality.
#' \item `'covariances'`: all loadings, (indicator-)intercepts, (indicator-)residuals, and latent covariances are restricted to equality.
#' \item `'means'`: all loadings, (indicator-)intercepts, (indicator-)residuals, latent covariances, and latent means are restricted to equality.
#' }
#' 
#' For example, setting `comparison = 'metric'` and `nullEffect = 'scalar'` determines power 
#' to reject the hypothesis that the constraints placed in the scalar invariance model 
#' (restricting loadings and intercepts) over the 
#' metric invariance model (restricting only the loadings) are defensible.
#'  
#' For greater flexibility, the models can also be defined using `lavaan` style restrictions as a vector, namely
#' `'none'` (no restrictions), `'loadings'` (loadings), `'intercepts'` (intercepts), `'residuals'` (residuals), `'lv.covariances'` (latent covariances), `'means'` (latent means).
#'  For instance: 
#' \itemize{
#' \item `'none'`: no invariance constraints and thus representing a configural invariance model. Shows the same fit as the saturated model, so only the delta df differ. 
#' \item `c('loadings')`: all loadings are restricted to equality. Note that reference scaling is used, so the first indicator should be invariant. 
#' \item `c('loadings', 'intercepts')`: all loadings and (indicator-)intercepts are restricted to equality. 
#' \item `c('loadings', 'intercepts', 'residuals')`: all loadings, (indicator-)intercepts, and (indicator-)residuals are restricted to equality.
#' \item `c('loadings', 'residuals')`: all loadings and (indicator-)residuals are restricted to equality.
#' \item `c('loadings', 'intercepts', 'means')`: all loadings, (indicator-)intercepts, and latent factor means are restricted to equality.
#' \item `c('loadings', 'residuals', 'lv.covariances')`: all loadings, (indicator-)residuals, and latent factor covariances are restricted to equality.
#' }
#' 
#' For example, setting `comparison = c('loadings')` and `nullEffect = 'c('loadings', 'intercepts')'` 
#' determines power to reject the hypothesis that the constraints placed in the scalar invariance model 
#' (restricting loadings and intercepts) over the  metric invariance model (restricting only the loadings) are defensible.
#' Note that variance scaling is used, so invariance of variances (`'lv.variances'`) is always met. Latent means are identified using single occasion identification.
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model**:
#' * `Lambda`: The factor loading matrix (with the number of columns equaling the number of factors).
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading on the second factor by .8, .6, .6, and .4.
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the loadings of the first factor to equal .5 and those of the second factor do equal .6.
#' * `Theta`: Variance-covariance matrix of the indicator residuals, which should be a diagonal matrix. Required when residual non-invariance is to be detected. When `NULL`, Theta is a diagonal matrix with elements such that all variances are 1. 
#' * `tau`: Defines the indicator intercepts, required whenever a model involves hypotheses about means (e.g., scalar invariance). If `NULL` and `Alpha` is set, all intercepts are assumed to equal zero.
#' * `Alpha`: Defines the latent means, required whenever a model involves hypotheses about latent means (e.g., latent mean invariance). If `NULL` and `tau` is set, all latent means are assumed to equal zero. Because variance scaling is used so that all factor variances are 1, latent mean differences can be interpreted akin to Cohen's d as standardized mean differences.
#' 
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` always need to be defined, 
#' and `Theta`, `tau` and `Alpha` need to be defined for particular levels of invariance.
#' 
#' 
#' Additional arguments related to the requested type of **power analysis**:
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.  This requires the `semTools` package.
#' 
#' @examples
#' \dontrun{
#' 
#' 
#' # obtain the required N to reject the hypothesis of metric invariance
#' # in comparison to the configural invariance model
#' # with a power of 80% on alpha = 5%
#' # for amodel involving a two factors (= two measurements) which
#' # is measured by 5 indicators
#' # loading by .5 each at the first measurement occasion
#' # loading by .6 each in the second measurement occasion,
#' # and assuming autocorrelated residuals
#' powerLI <- semPower.powerLI(
#'   type = 'a-priori', alpha = .05, power = .80, 
#'   comparison = 'configural',
#'   nullEffect = 'metric',
#'   nIndicator = c(5, 5),
#'   loadM = c(.5, .6),
#'   autocorResiduals = TRUE
#' )
#' 
#' # show summary
#' summary(powerLI)
#' 
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerLI$modelH1, sample.cov = powerLI$Sigma,
#'             sample.nobs = 1000, sample.cov.rescale = FALSE)
#' lavaan::sem(powerLI$modelH0, sample.cov = powerLI$Sigma,
#'             sample.nobs = 1000, sample.cov.rescale = FALSE)
#' 
#' 
#' 
#' # same as above, but determine power with N = 500 on alpha = .05
#' powerLI <- semPower.powerLI(
#'   type = 'post-hoc', alpha = .05, N = 500, 
#'   comparison = 'configural',
#'   nullEffect = 'metric',
#'   nIndicator = c(5, 5),
#'   loadM = c(.5, .6),
#'   autocorResiduals = TRUE
#' )
#' 
#' 
#' # same as above, but determine the critical chi-square with N = 500 in each
#' # group so that alpha = beta
#' powerLI <- semPower.powerLI(
#'   type = 'compromise', abratio = 1, N = 500, 
#'   comparison = 'configural',
#'   nullEffect = 'metric',
#'   nIndicator = c(5, 5),
#'   loadM = c(.5, .6),
#'   autocorResiduals = TRUE
#' )
#' 
#' 
#' # same as above, but compare to the saturated model
#' # (rather than to the configural invariance model)
#' powerLI <- semPower.powerLI(
#'   type = 'a-priori', alpha = .05, power = .80, 
#'   comparison = 'saturated',
#'   nullEffect = 'metric',
#'   nIndicator = c(5, 5),
#'   loadM = c(.5, .6),
#'   autocorResiduals = TRUE
#' )
#' 
#' 
#' # same as above, but provide individual factor loadings by group using a
#' # reduced loading matrix to define a  single factor model with three indicators
#' # loading by .4, .6, .5 at the first measurement occasion and
#' # loading by .5, .6, .7 at the second measurement occasion 
#' powerLI <- semPower.powerLI(
#'   type = 'a-priori', alpha = .05, power = .80, 
#'   comparison = 'configural',
#'   nullEffect = 'metric',
#'   loadings = list(
#'     c(.4, .6, .5),
#'     c(.5, .6, .7)
#'   ),
#'   autocorResiduals = TRUE
#' )
#' 
#' # obtain the required N to reject the hypothesis of scalar invariance
#' # in comparison to the metric invariance model
#' # with a power of 80% on alpha = 5%
#' # for a two factor model, where both factors are
#' # measured by 3 indicators each and all loadings equal .5 (at both measurements),
#' # all intercepts are 0.0 at the first measurement occasion, but
#' # all intercepts are 0.2 at the second measurement occasion and
#' powerLI <- semPower.powerLI(
#'   type = 'a-priori', alpha = .05, power = .80, 
#'   comparison = 'metric',
#'   nullEffect = 'scalar',
#'   nIndicator = c(5, 5),
#'   loadM = c(.5, .5),
#'   tau = c(0, 0, 0, 0, 0, 
#'           .2, .2, .2, .2, .2),
#'   autocorResiduals = TRUE
#' )
#' 
#' # same as above, but use lavaan strings 
#' powerLI <- semPower.powerLI(
#'   type = 'a-priori', alpha = .05, power = .80, 
#'   comparison = c('loadings'),
#'   nullEffect = c('loadings', 'intercepts'),
#'   nIndicator = c(5, 5),
#'   loadM = c(.5, .5),
#'   tau = c(0, 0, 0, 0, 0, 
#'           .2, .2, .2, .2, .2),
#'   autocorResiduals = TRUE
#' )
#' 
#' 
#' # obtain the required N to reject the hypothesis of equal latent means
#' # in comparison to the scalar invariance model;
#' # all intercepts are zero in both groups,
#' # at the first measurement occasion, the latent mean is 0.0,
#' # at the first measurement occasion, the latent mean is 0.5
#' powerLI <- semPower.powerLI(
#'   type = 'a-priori', alpha = .05, power = .80, 
#'   comparison = c('loadings', 'intercepts'),
#'   nullEffect = c('loadings', 'intercepts', 'means'),
#'   nIndicator = c(5, 5),
#'   loadM = c(.5, .5),
#'   tau = rep(0, 10),
#'   Alpha = c(0, .5),
#'   autocorResiduals = TRUE
#' )
#' 
#' # obtain the required N to reject the hypothesis of equal covariances
#' # in comparison to the residual invariance model;
#'  Phi <- matrix(c(
#'   c(1, .3, .1),
#'    c(.3, 1, .2),
#'    c(.1, .2, 1)
#'  ), nrow=3, byrow = TRUE)
#'  powerLI <- semPower.powerLI(
#'    type = 'a-priori', alpha = .05, power = .80,
#'    comparison = 'residual',
#'    nullEffect = 'covariances',
#'    nIndicator = c(3, 3, 3),
#'    loadM = c(.5, .5, .5),
#'    Phi = Phi,
#'    tau = rep(0, 9)
#' )   
#'  
#' # request a simulated post-hoc power analysis with 250 replications
#' # to reject the hypothesis of equal latent means.
#' set.seed(300121)
#' powerLI <- semPower.powerLI(
#'   type = 'post-hoc', alpha = .05, N = 500, 
#'   comparison = c('loadings', 'intercepts'),
#'   nullEffect = c('loadings', 'intercepts', 'means'),
#'   nIndicator = c(5, 5),
#'   loadM = c(.5, .5),
#'   tau = rep(0, 10),
#'   Alpha = c(0, .5),
#'   autocorResiduals = TRUE,
#'   simulatedPower = TRUE,
#'   simOptions = list(nReplications = 250)  
#' )
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerLI <- function(type, 
                             comparison = NULL,
                             nullEffect = NULL,
                             autocorResiduals = TRUE,
                             Phi = NULL,
                             ...){
  
  args <- list(...)
  
  # validate input
  checkEllipsis(...)
  
  # allow both lavstring and single string
  nullValid <- c('metric', 'scalar', 'residual', 'covariances', 'means')
  nullValidLav <- c('loadings', 'intercepts', 'residuals', 'lv.covariances', 'means')
  
  compValid <- c('saturated', 'configural', 'none', nullValid, nullValidLav)
  comparison <- lapply(comparison, function(x) checkNullEffect(x, compValid, 'comparison'))
  nullEffect <- lapply(nullEffect, function(x) checkNullEffect(x, c(nullValid, nullValidLav), 'nullEffect'))
  
  if(unlist(comparison)[1] == 'configural' || unlist(comparison)[1] == 'none') comparison <- 'configural'
  if(unlist(comparison)[1] == 'saturated') comparison <- 'saturated'
  
  # check and translate to lavstring
  if(length(nullEffect) > 1 || nullEffect == 'loadings'){
    if(any(unlist(lapply(nullEffect, function(x) x %in% nullValid[-5])))) stop('Either use lavaan-type strings or use predefined strings, but do not mix.')
    if(any(unlist(lapply(comparison, function(x) x %in% c(nullValid))))) stop('Either use lavaan-type strings or use predefined strings, but do not mix.')
    if((length(nullEffect) == 1 && nullEffect != 'loadings') && length(nullEffect) <= length(comparison)) stop('The H0 model must contain all restrictions of the comparison model plus one additional restriction.')
    if(!any(c('saturated', 'configural') %in% comparison) && any(unlist(lapply(comparison, function(x) !x %in% nullEffect)))) stop('The H0 model must contain all restrictions of the comparison model plus one additional restriction.')
  }else{
    if(nullEffect %in% nullValidLav[-5]) stop('Either use lavaan-type strings or use predefined strings, but do not mix.')
    if(comparison %in% nullValidLav[-5]) stop('Either use lavaan-type strings or use predefined strings, but do not mix.')
    if(which(c('saturated', 'configural', nullValid) %in% unlist(comparison)) >= (2 + which(nullValid %in% unlist(nullEffect)))) stop('Model defined in nullEffect is not nested in comparison model.')
    
    # translate to lavstring
    lc <- ln <- list()
    idxC <- which(c('saturated', 'configural', nullValid) == comparison)
    idxN <- which(c('saturated', 'configural',  nullValid) == nullEffect)
    if(idxC > 2) lc <- append(lc, 'loadings')
    if(idxN > 2) ln <- append(ln, 'loadings')
    if(idxC > 3) lc <- append(lc, 'intercepts')
    if(idxN > 3) ln <- append(ln, 'intercepts')
    if(idxC > 4) lc <- append(lc, 'residuals')
    if(idxN > 4) ln <- append(ln, 'residuals')
    if(idxC > 5) lc <- append(lc, 'lv.covariances')
    if(idxN > 5) ln <- append(ln, 'lv.covariances')
    if(idxC > 6) lc <- append(lc, 'means')
    if(idxN > 6) ln <- append(ln, 'means')
    nullEffect <- ln
    if(!any(c('saturated', 'configural') %in% comparison)) comparison <- lc
  }
  comparison <- unlist(comparison)
  nullEffect <- unlist(nullEffect)
  
  if('means' %in% nullEffect && !any(c('saturated', 'intercepts') %in% comparison)) stop('Latent means cannot be estimated without constraints on intercepts (scalar invariance).')
  
  
  ### generate sigma
  
  # use variance scaling, so constraints on loadings are always valid.
  generated <- semPower.genSigma(..., Phi = Phi, useReferenceIndicator = FALSE)
  
  if(is.null(generated[['Lambda']])) stop('powerLI does not support multiple group models (do not define lists for measurement parameters).')
  if(ncol(generated[['Lambda']]) == nrow(generated[['Lambda']])) stop('Longitudinal invariance requires latent variables with multiple indicators.')
  if(ncol(generated[['Lambda']]) < 3 && 'lv.covariances' %in% c(comparison, nullEffect)) stop('Equality of covariances is only meaningful for at least 3 factors (measurements).')
  if(!is.null(Phi)){
    if(!is.matrix(Phi) && length(Phi) > 1) stop('If there are more than two factors, Phi must be a matrix.')
    if(is.matrix(Phi) && ncol(Phi) != ncol(generated[['Lambda']])) stop('The dimensions of Phi must equal the number of factors.')
  }
  
  metricInvarianceFactors <- list(seq(ncol(generated[['Lambda']])))  
  
  #### gen model strings
  
  # configural
  modelH1 <- modelH0 <- generated[['modelTrueCFA']]
  
  # loadings
  if('loadings' %in% c(comparison, nullEffect)){
    generatedMetric <- semPower.genSigma(..., Phi = Phi, useReferenceIndicator = FALSE, metricInvariance = metricInvarianceFactors)
    if('loadings' %in% comparison) modelH1 <- generatedMetric[['modelTrueCFA']]
    if('loadings' %in% nullEffect) modelH0 <- generatedMetric[['modelTrueCFA']]
  }
  # intercepts
  if('intercepts' %in% c(comparison, nullEffect)){
    tok <- list()
    for(x in seq_along(metricInvarianceFactors)){
      ci <- lapply(metricInvarianceFactors[[x]], function(f) paste0('x', which(generated[['Lambda']][, f] != 0)))
      lab <- paste0('i', 1:length(ci[[1]]))
      tok <- append(tok, unlist(lapply(ci, function(x) paste0(x, ' ~ ', lab, '*1') )))
    }
    if('intercepts' %in% comparison) modelH1 <- append(modelH1, tok)
    if('intercepts' %in% nullEffect) modelH0 <- append(modelH0, tok)
    ## TODO or shall we estimate latent means when there are constraints on intercepts?
    
    # means
    if('means' %in% nullEffect){
      tok <- 'f1 ~ 0*1'   # single occasion identification
      ci <- paste0('f', 2:ncol(generated[['Lambda']]))
      # estm means in H1 model
      tokH1 <- append(tok, unlist(lapply(ci, function(x) paste0(x, ' ~ 1') )))
      modelH1 <- append(modelH1, tokH1)
      # restrict means in H0 model
      tokH0 <- append(tok, unlist(lapply(ci, function(x) paste0(x, ' ~ 0*1') )))
      modelH0 <- append(modelH0, tokH0)
    }
    
  }
  # residuals
  if('residuals' %in% c(comparison, nullEffect)){
    tok <- list()
    for(x in seq_along(metricInvarianceFactors)){
      ci <- lapply(metricInvarianceFactors[[x]], function(f) paste0('x', which(generated[['Lambda']][, f] != 0)))
      lab <- paste0('r', 1:length(ci[[1]]))
      tok <- append(tok, unlist(lapply(ci, function(x) paste0(x, ' ~~ ', lab, '*', x) )))
    }
    if('residuals' %in% comparison) modelH1 <- append(modelH1, tok)
    if('residuals' %in% nullEffect) modelH0 <- append(modelH0, tok)
  }
  # lv covar
  if('lv.covariances' %in% c(comparison, nullEffect)){
    tok <- list()
    for(i in 1:(ncol(generated[['Lambda']]) - 1)){
      for(j in (i + 1):ncol(generated[['Lambda']])){
        tok <- append(tok, paste0('f', i, '~~ c*f', j))
      }
    }
    if('lv.covariances' %in% comparison) modelH1 <- append(modelH1, tok)
    if('lv.covariances' %in% nullEffect) modelH0 <- append(modelH0, tok)
  }

  
  # add autocorrelated residuals
  if(autocorResiduals){
    # do this only when there is at least one latent variable
    autocorResidualsFactors <- metricInvarianceFactors  # same structure 
    tok <- list()
    for(x in seq_along(autocorResidualsFactors)){
      ci <- lapply(autocorResidualsFactors[[x]], function(f) paste0('x', which(generated[['Lambda']][, f] != 0)))
      if(length(ci[[1]]) > 1){
        for(i in 1:(length(ci) - 1)){
          for(j in (i + 1) : length(ci)){
            tok <- append(tok, paste(ci[[i]], '~~', ci[[j]]))
          }
        }
      }
    }
    modelH1 <- append(modelH1, tok)
    modelH0 <- append(modelH0, tok)
  }

    
  
  modelH0 <- paste(unlist(modelH0), collapse = '\n')
  modelH1 <- paste(unlist(modelH1), collapse = '\n')
  
  if(comparison[1] == 'saturated'){
    modelH1 <- NULL
  } 
  
  mu <- NULL
  if(any(c('intercepts', 'means') %in% c(comparison, nullEffect)))
    mu <- generated[['mu']]
  
  
  do.call(semPower.powerLav, append(list(
    type = type,
    Sigma = generated[['Sigma']],
    mu = mu,
    modelH0 = modelH0,
    modelH1 = modelH1,
    fitH1model = TRUE),
    args)
  )

}


#' semPower.powerPath
#'
#' Convenience function for performing power analyses for hypothesis arising 
#' in a generic path model.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param Beta matrix of regression slopes between latent variables (all-Y notation). A list for multiple group models. Exogenous variables must occupy the first rows in `Beta` when `standardized = TRUE`. See details. 
#' @param Psi variance-covariance matrix of latent (residual) factors. If `standardized = TRUE`, the diagonal is ignored and all off-diagonal elements are treated as correlations. If `NULL`, an identity matrix is assumed. A list for multiple group models. See details.
#' @param nullEffect defines the hypothesis of interest, must be one of `'beta = 0'` (the default) to test whether a regression slope is zero, `'betaX = betaZ'` to test for the equality of slopes, and `'betaX = betaZ'` to test for the equality of a slope across groups. Define the slopes to be set to equality in `nullWhich` and the groups in `nullWhichGroups`. 
#' @param nullWhich vector of size 2 indicating which slope in `Beta` is hypothesized to equal zero when `nullEffect = 'beta = 0'`, or to restrict to equality across groups when `nullEffect = 'betaA = betaB'`, or list of vectors defining which correlations to restrict to equality when `nullEffect = 'betaX = betaZ'`. Can also contain more than two slopes, e.g., `list(c(2, 1), c(3, 1), c(3, 2))` to set `Beta[2, 1] = Beta[3, 1] = Beta[3, 2]`.
#' @param nullWhichGroups for `nullEffect = 'betaA = betaB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param standardized whether all parameters should be standardized (`TRUE`, the default). If `FALSE`, all regression relations are unstandardized.
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details
#' This function performs a power analysis to reject a hypothesis arising
#' in a generic structural equation model specifying regression relations between the factors via the Beta matrix:  
#' * `nullEffect = 'beta = 0'`: Tests the hypothesis that a slope is zero. 
#' * `nullEffect = 'betaX = betaZ'`: Tests the hypothesis that two or more slopes are equal to each other.
#' * `nullEffect = 'betaA = betaB'`: Tests the hypothesis that a slope is equal in two or more groups (always assuming metric invariance).
#' 
#' This function provides a generic way to perform power analyses (as compared to other functions covering special cases in a more accessible manner).
#' 
#' A specific hypothesis is defined by setting `nullEffect` to define the hypothesis type, 
#' `nullWhich` to define the slope(s) that are targeted, and by providing 
#' the `Beta` (and optionally the `Psi`) matrix to define the population structure.
#'  
#' To understand the structure of `Beta` and `Psi`, consider the general structural equation model, 
#' \deqn{\Sigma = \Lambda (I - B)^{-1} \Psi [(I - B)^{-1}]'  \Lambda' + \Theta } 
#' where \eqn{B} is the \eqn{m \cdot m} matrix containing the regression slopes and \eqn{\Psi} is the (residual) variance-covariance matrix of the \eqn{m} factors. 
#' 
#' As an example, suppose there are four factors (X1, X2, X3, X4), and Beta is defined as follows:
#' \eqn{
#' \begin{array}{lrrr} 
#'     & X_1 & X_2 & X_3 & X_4\\ 
#' X_1 & 0.0 & 0.0 & 0.0 & 0.0 \\ 
#' X_2 & 0.0 & 0.0 & 0.0 & 0.0  \\ 
#' X_3 & 0.2 & 0.3 & 0.0 & 0.0  \\ 
#' X_4 & 0.3 & 0.5 & 0.0 & 0.0  \\ 
#' \end{array}
#' }
#' 
#' Each row specifies how a particular factor is predicted by the available factors, 
#' so the above implies the following regression relations:
#' 
#' \eqn{
#' X_1 = 0.0 \cdot X_1 +  0.0 \cdot X_2 + 0.0 \cdot X_3 + 0.0 \cdot X_4 \\
#' X_2 = 0.0 \cdot X_1 +  0.0 \cdot X_2 + 0.0 \cdot X_3 + 0.0 \cdot X_4 \\
#' X_3 = 0.2 \cdot X_1 +  0.3 \cdot X_2 + 0.0 \cdot X_3 + 0.0 \cdot X_4 \\
#' X_4 = 0.3 \cdot X_1 +  0.5 \cdot X_2 + 0.0 \cdot X_3 + 0.0 \cdot X_4 
#' }
#' 
#' which simplifies to
#' 
#' \eqn{
#' X_3 = 0.2 \cdot X_1 + 0.3 \cdot X_2 \\
#' X_4 = 0.3 \cdot X_1 + 0.5 \cdot X_2 
#' }
#' 
#' Further suppose that Psi is
#' \eqn{
#' \begin{array}{lrrr} 
#'     & X_1 & X_2 & X_3 & X_4\\ 
#' X_1 & 1.0 & 0.3 & 0.0 & 0.0 \\ 
#' X_2 & 0.3 & 1.0 & 0.0 & 0.0 \\ 
#' X_3 & 0.0 & 0.0 & 1.0 & 0.2 \\ 
#' X_4 & 0.0 & 0.0 & 0.2 & 1.0 \\ 
#' \end{array}
#' }
#' 
#' which implies a correlation between X1 and X2 of .3 and a residual correlation
#' between X3 and X4 of .2. 
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model**:
#' * `Lambda`: The factor loading matrix (with the number of columns equaling the number of factors).
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading on the second factor by .8, .6, .6, and .4.
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the loadings of the first factor to equal .5 and those of the second factor do equal .6.
#' 
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` always need to be defined. 
#' 
#' Additional arguments related to the requested type of **power analysis**:
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.  This requires the `semTools` package.
#' 
#' @examples
#' \dontrun{
#' # set up pathmodel in the form of
#' # f2 = .2*f1
#' # f3 = .3*f2
#' # f4 = .1*f1 + .4*f3
#' # obtain the required N to detect that the 
#' # slope f1 -> f4 is >= .10 
#' # with a power of 95% on alpha = 5%
#' # where f1 is measured by 3, f2 by 4, f3 by 5, and f4 by 6 indicators, 
#' # and all loadings are .5
#' Beta <- matrix(c(
#'   c(.00, .00, .00, .00),       # f1
#'   c(.20, .00, .00, .00),       # f2
#'   c(.00, .30, .00, .00),       # f3
#'   c(.10, .00, .40, .00)        # f4
#' ), byrow = TRUE, ncol = 4)
#' powerPath <- semPower.powerPath(type = 'a-priori',
#'                                 Beta = Beta,
#'                                 nullWhich = c(4, 1),
#'                                 nIndicator = c(3, 4, 5, 6), 
#'                                 loadM = .5,
#'                                 alpha = .05, beta = .05)
#' # show summary
#' summary(powerPath)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerPath$modelH1, sample.cov = powerPath$Sigma,
#' sample.nobs = powerPath$requiredN, sample.cov.rescale = FALSE)
#' lavaan::sem(powerPath$modelH0, sample.cov = powerPath$Sigma,
#' sample.nobs = powerPath$requiredN, sample.cov.rescale = FALSE)
#' 
#' # same as above, but detect that the slope f3 -> f4 is >= .30 
#' powerPath <- semPower.powerPath(type = 'a-priori',
#'                                 Beta = Beta,
#'                                 nullWhich = c(4, 3),
#'                                 nIndicator = c(3, 4, 5, 6), 
#'                                 loadM = .5,
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but detect that 
#' # the slope f1 -> f2 (of .20) differs from the slope f2 -> f3 (of .30) 
#' powerPath <- semPower.powerPath(type = 'a-priori',
#'                                 Beta = Beta,
#'                                 nullEffect = 'betaX = betaZ',
#'                                 nullWhich = list(c(2, 1), c(3, 2)),
#'                                 nIndicator = c(3, 4, 5, 6), 
#'                                 loadM = .5,
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but consider a multiple group model with equally sized groups, 
#' # and obtain the required N to detect that the slope 
#' # in group 1 (of .20) differs from the one in group 2 (of .40)
#' Beta1 <- Beta2 <- matrix(c(
#'   c(.00, .00, .00, .00),       # f1
#'   c(.20, .00, .00, .00),       # f2
#'   c(.00, .30, .00, .00),       # f3
#'   c(.10, .00, .40, .00)        # f4
#' ), byrow = TRUE, ncol = 4)
#' Beta2[2, 1] <- .40
#' Beta <- list(Beta1, Beta2)
#' powerPath <- semPower.powerPath(type = 'a-priori',
#'                                 Beta = Beta,
#'                                 nullEffect = 'betaA = betaB',
#'                                 nullWhich = c(2, 1),
#'                                 nIndicator = list(
#'                                   c(3, 4, 5, 6), 
#'                                   c(3, 4, 5, 6)), 
#'                                 loadM = .5,
#'                                 alpha = .05, beta = .05, N = list(1, 1))
#' }
#' @seealso [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerPath <- function(type, comparison = 'restricted',
                               Beta,
                               Psi = NULL,
                               nullEffect = 'beta = 0',
                               nullWhich = NULL, 
                               nullWhichGroups = NULL,
                               standardized = TRUE,
                               ...){
  
  args <- list(...)
  comparison <- checkComparisonModel(comparison)
  checkEllipsis(...)
  
  # we override Sigma later, so let's make sure it is not set in ellipsis argument
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma, because Sigma is determined as function of Beta (or the slopes).')
  
  # validate input
  nullEffect <- checkNullEffect(nullEffect, c('beta=0', 'betax=betaz', 'betaa=betab'))
  if(is.null(Beta)) stop('Beta may not be NULL.')
  isMultigroup <- is.list(Beta) && length(Beta) > 1
  if(!is.list(Beta)) Beta <- list(Beta)
  if(any(unlist(lapply(Beta, function(x) any(diag(x) != 0) )))) stop('All diagonal elements of Beta must be zero.')
  if(standardized && any(unlist(lapply(Beta, function(x) any(x[upper.tri(x, diag = TRUE)] != 0))))) stop('All upper triangular elements in Beta must be zero when requesting standardized parameters. Remember exogenous variables must occupy the first rows in Beta.')
  if(isMultigroup && (length(unique(unlist(lapply(Beta, ncol)))) > 1 || length(unique(unlist(lapply(Beta, nrow)))) > 1)) stop('Beta must be of same dimension for all groups') 
  lapply(Beta, function(x) checkSquare(x, 'Beta'))
  if(!is.null(Psi)){
    if(!is.list(Psi)) Psi <- list(Psi)
    if(isMultigroup && (length(unique(unlist(lapply(Psi, ncol)))) > 1 || length(unique(unlist(lapply(Psi, nrow)))) > 1)) stop('Psi must be of same dimension for all groups') 
    lapply(Psi, function(x) checkSymmetricSquare(x, 'Psi'))
    if(any(unlist(lapply(Psi, function(x) any(eigen(x)$values <= 0))))) warning('Phi is not positive definite.')
    if(ncol(Psi[[1]]) != ncol(Beta[[1]])) stop('Beta and Psi must be of same dimension.')
  }
  if(is.null(nullWhich)) stop('nullWhich must not be NULL.')
  if(any(unlist(lapply(nullWhich, function(x) any(x > ncol(Beta[[1]])))))) stop('At least one element in nullWhich is an out of bounds index concerning Beta.')
  if(nullEffect == 'betax=betaz' && any(lapply(nullWhich, function(x) length(x)) != 2)) stop('nullWhich must be a list containing vectors of size two each.')
  if(nullEffect == 'betaa=betab' && !isMultigroup) stop('Beta must be a list for multiple group comparisons.')
  if(nullEffect != 'betaa=betab' && isMultigroup) stop('Multiple group models are only valid for nullEffect = "betaA=betaB".')
  if(nullEffect == 'beta=0' && any(unlist(lapply(Beta, function(x) x[nullWhich[1], nullWhich[2]] == 0)))) stop('nullWhich must not refer to a slope with a population value of zero.')
  if(!is.null(nullWhichGroups)) lapply(nullWhichGroups, function(x) checkBounded(x, 'All elements in nullWhichGroups', bound = c(1, length(Beta)), inclusive = TRUE))
  if(!is.list(nullWhich)) nullWhich <- list(nullWhich)
  if(isMultigroup && is.null(nullWhichGroups)) nullWhichGroups <- seq_along(Beta)
  
  
  ### get Sigma
  if(standardized){
    Phi <- lapply(seq_along(Beta), function(x) getPhi.B(Beta[[x]], Psi[[x]]))
    generated <- semPower.genSigma(Phi = if(!isMultigroup) Phi[[1]] else Phi, 
                                   useReferenceIndicator = TRUE, ...)  
  }else{  
    generated <- semPower.genSigma(Beta = if(!isMultigroup) Beta[[1]] else Beta, 
                                   Psi = if(!isMultigroup) Psi[[1]] else Psi, 
                                   useReferenceIndicator = TRUE, ...)
  }
  
  ### create model strings
  # we need to use modelTrueCFA and define regression relations here,
  # because we need labels for the H0 model and because standardized based on phi yields no regression structure 
  if(!isMultigroup) model <- generated[['modelTrueCFA']] else model <- generated[[1]][['modelTrueCFA']]
  tok <- list()
  tokH1 <- list()
  for(f in seq(ncol(Beta[[1]]))){
    ifelse(isMultigroup, idx <- unique(unlist(lapply(Beta, function(x) which(x[f, ] != 0)))), idx <- which(Beta[[1]][f, ] != 0))
    if(length(idx) > 0){
      # cIdx <- lapply(idx, function(x) c(f, x)) %in% nullWhich   ## not sure why this does not work when looping over f?
      cIdx <- unlist(lapply(lapply(idx, function(x) c(f, x)), function(y) any(unlist(lapply(nullWhich, function(z) all(y == z))))))
      if(nullEffect == 'beta=0'){
        prefix <- rep('', length(idx))
        prefix[cIdx] <- '0*'
      }else if(nullEffect == 'betax=betaz'){
        prefix <- rep('', length(idx))
        prefix[cIdx] <- 'pc*'
      }else if(nullEffect == 'betaa=betab'){
        clab <- paste0('pf', paste0(formatC(f, width = 2, flag = 0), formatC(idx, width = 2, flag = 0)))
        prefix <- lapply(clab, function(x) paste0(x, '_g', seq_along(Beta)))
        prefix[cIdx] <- lapply(prefix[cIdx], function(x) unlist(lapply(x, function(y) gsub(paste0('_g', nullWhichGroups, collapse = '|'), '_gc', y))))
        prefix <- unlist(lapply(prefix, function(x) paste0('c(', paste(x, collapse = ', ') ,')*')))
      }
      tok <- append(tok, paste0('f', f, ' ~ ', paste(prefix, paste0('f', idx), sep = '', collapse = ' + ')))
      tokH1 <- append(tokH1, paste0('f', f, ' ~ ', paste(paste0('f', idx), sep = '', collapse = ' + ')))
    }
  }
  # (residual) correlations, assuming the same Psi across groups
  if(is.null(Psi)) Psi <- list(diag(ncol(Beta[[1]])))
  for(f in 1:(ncol(Beta[[1]]) - 1)){
    for(ff in (f + 1):ncol(Beta[[1]])){
      if(Psi[[1]][f, ff] != 0){
        tok <- append(tok, paste0('f', f, ' ~~ f', ff))
        tokH1 <- append(tokH1, paste0('f', f, ' ~~ f', ff))
      }
    }
  }
  modelH0 <- paste(c(model, unlist(tok)), collapse = '\n')
  modelH1 <- paste(c(model, unlist(tokH1)), collapse = '\n')
  
  # always enforce invariance constraints in the multigroup case
  if(isMultigroup){
    args[['lavOptions']] <- append(args[['lavOptions']], list(group.equal = c('loadings')))
  } 
  
  # always fit H1 model
  fitH1model <- TRUE 
  if(comparison == 'saturated'){
    modelH1 <- NULL
  } 
  
  if(isMultigroup) Sigma <- lapply(generated, '[[', 'Sigma') else Sigma <- generated[['Sigma']] 
  
  do.call(semPower.powerLav, append(list(
    type = type,
    Sigma = Sigma,
    modelH0 = modelH0,
    modelH1 = modelH1,
    fitH1model = fitH1model),
    args)
  )
}

#' semPower.powerBifactor
#'
#' Perform a power analysis for models including one or more bifactors to reject one of the following hypotheses: 
#' (a) a zero correlation between two factors, (b) the equality of two correlations between factors,
#' or (c) the equality of a correlation between two factors across two or more groups. 
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param bfLoadings a single vector or a list containing one or more vectors giving the loadings on each bifactor. For example, list(rep(.6, 10), rep(.6, 10)) defines two bifactors with 10 indicators each, loading by .6 each. Can be a list of lists for multiple group models.
#' @param bfWhichFactors a list containing one or more vectors defining which (specific) factors defined in the respective arguments in ... are part of the bifactor structure. See details.
#' @param Phi either a single number defining the correlation between exactly two factors or the factor correlation matrix. Must only contain the bifactor(s) and the covariate(s). Must be a list for multiple group models. Phi assumes the following order (bifactor_1, bifactor_2, ..., bifactor_j, covariate_1,  covariate_2, ...,  covariate_k). See details.
#' @param nullEffect defines the hypothesis of interest, must be one of `'cor = 0'` (the default) to test whether a correlation is zero, `'corX = corZ'` to test for the equality of correlations, and `'corA = corB'` to test for the equality of a correlation across groups. Define the correlations to be set to equality in `nullWhich` and the groups in `nullWhichGroups`. 
#' @param nullWhich vector of size 2 indicating which factor correlation in `Phi` is hypothesized to equal zero when `nullEffect = 'cor = 0'`, or to restrict to equality across groups when `nullEffect = 'corA = corB'`, or list of vectors defining which correlations to restrict to equality when `nullEffect = 'corX = corZ'`. Can also contain more than two correlations, e.g., `list(c(1, 2), c(1, 3), c(2, 3))` to set `Phi[1, 2] = Phi[1, 3] = Phi[2, 3]`. 
#' @param nullWhichGroups for `nullEffect = 'corA = corB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model concerning the specific factors and the covariate(s). See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details 
#' 
#' This function performs a power analysis to reject various hypotheses arising
#' in a model including a bifactor structure:
#' * `nullEffect = 'cor = 0'`: Tests the hypothesis that the correlation between a bifactor and another factor (which can also be a bifactor) is zero.
#' * `nullEffect = 'corX = corZ'`: Tests the hypothesis that two or more correlations involving one or more bifactors are equal to each other.
#' * `nullEffect = 'corA = corB'`: Tests the hypothesis that the correlation between the bifactor and another factor (which can also be a  bifactor) is equal in two or more groups (always assuming metric invariance).
#' 
#' A bifactor structure is defined by specifying the loadings on the general factor in `bfLoadings`, the comprised specific 
#' factors in `bfWhichFactors`, and the loadings on the specific factors in either `Lambda`, or `loadings`, 
#' or `nIndicator` and `loadM`. The latter arguments also include the loadings defining the 
#' covariate(s).
#' 
#' The correlations betwen the bifactor(s) and the covariate(s) are defined in `Phi`, which 
#' must omit the specific factors and only includes the bifactor(s) and the covariate(s) assuming 
#' the following order: (bifactor_1, bifactor_2, ..., bifactor_j, covariate_1,  covariate_2, ...,  covariate_k).
#' 
#' For example, the following defines a single bifactor with 10 indicators loading by .5 each. 
#' The bifactor structure involves 3 specific factors measured by 3 indicators each, each loading by
#' .3, .2, and .1 on the respective specific factor (in addition to the loadings on the bifactor). 
#' Furthermore, two covariate with 5 indicators each, all loading by .7, are defined. The correlation
#' between the covariates is .5, the one between the bifactor and the first and second covariate are
#' .3 and .2, respectively.
#' 
#' ```
#' bfLoadings <- rep(.5, 10)
#' bfWhichFactors <- c(1, 2, 3)
#' loadings <- list(
#'   rep(.3, 3),   # specific factor 1
#'   rep(.2, 3),   # specific factor 2
#'   rep(.1, 3),   # specific factor 3
#'   rep(.7, 5),   # covariate 1
#'   rep(.7, 5)    # covariate 2
#' )
#' Phi <- matrix(c(
#'   c(1, .3, .2),   # bifactor
#'   c(.3, 1, .5),   # covariate 1
#'   c(.2, .5, 1)   # covariate 2
#' ), ncol = 3, byrow = TRUE) 
#' ```
#'  
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model** concerning the specific factors and the covariate(s):
#' * `Lambda`: The factor loading matrix (with the number of columns equaling the number of factors).
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading on the second factor by .8, .6, .6, and .4.
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the loadings of the first factor to equal .5 and those of the second factor do equal .6.
#' 
#' Additional arguments related to the requested type of **power analysis**:
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.  This requires the `semTools` package.
#' 
#' @examples
#' \dontrun{
#' # get required N to detect a correlation of >= .3 between
#' # a single bifactor with 11 indicators all loadings by .6
#' # spanning the indicators of 3 specific factors
#' # with three indicators each, loading by .2, .15, and .25, respectively
#' # and a covariate measured by 4 indicators loading by .7 each,
#' # with a power of 95% on alpha = 5%
#' bfLoadings <- rep(.6, 11)
#' bfWhichFactors <- c(1, 2, 3)
#' loadings <- list(
#'   # specific factors
#'   rep(.2, 3),
#'   rep(.15, 3),
#'   rep(.25, 3),
#'   # covariate
#'   rep(.7, 4)
#' )
#' Phi <- .3    # bifactor - covariate
#' powerbifactor <- semPower.powerBifactor(type = 'a-priori',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         alpha = .05, beta = .05)
#' # show summary
#' summary(powerbifactor)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerbifactor$modelH1, sample.cov = powerbifactor$Sigma,
#'             sample.nobs = powerbifactor$requiredN, 
#'             sample.cov.rescale = FALSE)
#' lavaan::sem(powerbifactor$modelH0, sample.cov = powerbifactor$Sigma,
#'             sample.nobs = powerbifactor$requiredN, 
#'             sample.cov.rescale = FALSE)
#' 
#' # same as above, but determine power with N = 500 on alpha = .05
#' powerbifactor <- semPower.powerBifactor(type = 'post-hoc',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         alpha = .05, N = 500)
#' 
#' # same as above, but determine the critical chi-square with N = 500 so that alpha = beta
#' powerbifactor <- semPower.powerBifactor(type = 'compromise',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         abratio = 1, N = 500)
#' 
#' # same as above, but compare to the saturated model
#' # (rather than to the less restricted model)
#' powerbifactor <- semPower.powerBifactor(type = 'a-priori',
#'                                         comparison = 'saturated',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         alpha = .05, beta = .05)
#' 
#' # define two bifactors with 10 indicators each, where
#' # all loadings are .6 on the first and .5 on the second bifactor.
#' # the first bifactor spans the indicators of specific factors 1-3,
#' # the second bifactor spans the indicators of specific factors 4-6,
#' # all specific factors are measured by three indicators each,
#' # loadings are .2, .15, .25, .1, .15., and.2, respectively.
#' # define an additional  covariate measured by 4 indicators loading by .6 each.
#' # get required N to detect a correlation of >= .3 between the bifactors
#' # with a power of 95% on alpha = 5%
#' bfLoadings <- list(rep(.6, 10),
#'                    rep(.6, 10))
#' bfWhichFactors <- list(c(1, 2, 3),
#'                        c(4, 5, 6))
#' loadings <- list(
#'   # specific factors for bf1
#'   rep(.2, 3),
#'   rep(.15, 3),
#'   rep(.25, 3),
#'   # specific factors bf2
#'   rep(.1, 3),
#'   rep(.15, 3),
#'   rep(.2, 3),
#'   # covariate
#'   rep(.6, 4)
#' )
#' Phi <- diag(3)
#' Phi[1, 2] <- Phi[2, 1] <- .3    # bifactor1 - bifactor2
#' Phi[1, 3] <- Phi[3, 1] <- .5    # bifactor1 - covariate
#' Phi[2, 3] <- Phi[3, 2] <- .1    # bifactor2 - covariate
#' 
#' powerbifactor <- semPower.powerBifactor(type = 'a-priori',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         alpha = .05, beta = .05)
#' 
#' # same as above, but get required N to detect that
#' # the correlation between the first bifactor and the covariate (of r=.5) differs from
#' # the correlation between the second bifactor and the covariate (of r=.1)
#' powerbifactor <- semPower.powerBifactor(type = 'a-priori',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullEffect = 'corx = corz',
#'                                         nullWhich = list(c(1, 3), c(2, 3)),
#'                                         loadings = loadings,
#'                                         alpha = .05, beta = .05)
#' 
#' # multiple group example: get required N to detect that
#' # the correlation of a bifactor with 10 indicators
#' # spanning three specific factors with 3 indicators each
#' # to a covariate in group 1  (of r = .3)
#' # differs from the one in group 2 (of r = .1)
#' bfLoadings <- rep(.6, 10)
#' bfWhichFactors <- c(1, 2, 3)
#' 
#' loadings <- list(
#'   # specific factors
#'   rep(.2, 3),
#'   rep(.15, 3),
#'   rep(.25, 3),
#'   # covariate
#'   rep(.7, 4)
#' )
#' Phi1 <- Phi2 <- diag(2)
#' Phi1[1, 2] <- Phi1[2, 1] <- .3    # bifactor - covariate
#' Phi2[1, 2] <- Phi2[2, 1] <- .1    # bifactor - covariate
#' Phi <- list(Phi1, Phi2)
#' powerbifactor <- semPower.powerBifactor(type = 'a-priori',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullEffect = 'corA = corB',
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         alpha = .05, beta = .05, 
#'                                         N = list(1, 1))
#'                                         
#' # request a simulated post-hoc power analysis with 500 replications.
#' bfLoadings <- rep(.6, 11)
#' bfWhichFactors <- c(1, 2, 3)
#' loadings <- list(
#'   # specific factors
#'   rep(.2, 3),
#'   rep(.15, 3),
#'   rep(.1, 3),
#'   # covariate
#'   rep(.7, 5)
#' )
#' Phi <- .2  
#' set.seed(300121)
#' powerbifactor <- semPower.powerBifactor(type = 'post-hoc',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         alpha = .05, N = 500, 
#'                                         simulatedPower = TRUE,
#'                                         simOptions = list(nReplications = 500)
#'                                         )
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerBifactor <- function(type, comparison = 'restricted', 
                                   bfLoadings = NULL,
                                   bfWhichFactors = NULL,
                                   Phi = NULL,
                                   nullEffect = 'cor = 0',
                                   nullWhich = NULL, 
                                   nullWhichGroups = NULL, 
                                   ...){
  
  # validate input
  checkEllipsis(...)
  args <- list(...)
  comparison <- checkComparisonModel(comparison)
  if(is.null(Phi)) stop('Phi must be defined')
  nullEffect <- checkNullEffect(nullEffect, c('cor=0', 'corx=corz', 'cora=corb'))
  
  # assume multigroup structure when phi is a list
  isMultigroup <- is.list(Phi)
  if(!is.list(Phi)) Phi <- list(Phi)
  nGroups <- length(Phi)
  Phi <- lapply(Phi, function(x){
    if(length(x) == 1){
      matrix(c(1, x, x, 1), ncol = 2)
    }else{
      x
    } 
  })
  if(length(unique(unlist(lapply(Phi, function(x) ncol(x))))) != 1 || length(unique(unlist(lapply(Phi, function(x) nrow(x))))) != 1) stop('Phi must be of same dimension in all groups.')
  if(nullEffect == 'cora=corb' && !isMultigroup) stop('Phi must be a list when nullEffect = "corA = corB".')
  if(nullEffect != 'cora=corb' && isMultigroup) stop('Multigroups are only supported for nullEffect = "corA = corB"')
  
  if(!is.list(bfLoadings)) bfLoadings <- list(list(bfLoadings)) # (groups)(bifactors)
  if(!is.list(bfLoadings[[1]])) bfLoadings <- list(bfLoadings) # (groups)(bifactors)
  if(!is.list(bfWhichFactors)) bfWhichFactors <- list(list(bfWhichFactors)) # (groups)(bifactors)
  if(!is.list(bfWhichFactors[[1]])) bfWhichFactors <- list(bfWhichFactors) # (groups)(bifactors)
  # assume same measurement model for all groups if not specified otherwise
  if(isMultigroup && length(bfLoadings) == 1) bfLoadings <- lapply(seq(nGroups), function(x) bfLoadings[[1]])
  if(isMultigroup && length(bfWhichFactors) == 1) bfWhichFactors <- lapply(seq(nGroups), function(x) bfWhichFactors[[1]])
  numBifactors <- length(bfLoadings[[1]])
  
  if(nullEffect == 'corx=corz' && !is.list(nullWhich)) stop('nullWhich must be a list when nullEffect = "corx=corz".')
  if(!is.list(nullWhich)) nullWhich <- list(nullWhich)
  if(any(unlist(lapply(nullWhich, function(x) length(x) != 2)))) stop('nullWhich may only contain vectors of size two.')
  if(any(unlist(lapply(nullWhich, function(x) x[1] == x[2])))) stop('elements in nullWhich may not refer to variances.')
  if(max(unlist(nullWhich)) > ncol(Phi[[1]])) stop('At least one element in nullWhich is an out of bounds index concerning Phi. Recall that Phi must only comprise bifactor(s) and covariate(s), but not specific factors.')
  
  ### create Lambda
  # sLambda only contains specific factors and covariate(s)
  sLambda <- args[['Lambda']]
  if(is.null(sLambda)){
    sLambda <- lapply(seq(nGroups), function(x){
      if(is.list(args[['loadings']][[1]]) || 
         any(unlist(lapply(args[names(args) %in% c('loadM', 'nIndicator', 'loadSD', 'loadMinMax')], is.list)))){
        genLambda(args[['loadings']][[x]], args[['nIndicator']][[x]],
                  args[['loadM']][[x]], args[['loadSD']][[x]], args[['loadMinMax']][[x]])
      }else{
        genLambda(args[['loadings']], args[['nIndicator']],
                  args[['loadM']], args[['loadSD']], args[['loadMinMax']])
      }
    })
  }
  if(isMultigroup && length(sLambda) == 1) sLambda <- lapply(seq(nGroups), function(x) rep(sLambda[[1]], nGroups))
  if(any(unlist(lapply(seq(nGroups), function(g) lapply(seq(numBifactors), function(f){
    length(bfLoadings[[g]][[f]]) < sum(sLambda[[g]][ , bfWhichFactors[[g]][[f]]] != 0)
  }))))) stop('Bifactors must have at least the same number of indicators as the comprised by the respective specific factors. If you want an indicator only loading on a specific factor, assign the respective bifactor loading a value of zero.')
  
  # add bifactor structure to Lambda
  Lambda <- lapply(seq(nGroups), function(g){
    csLambda <- sLambda[[g]]
    cBfLoadings <- bfLoadings[[g]]
    cBfWhichFactors <- bfWhichFactors[[g]]
    extraIndicators <- unlist(lapply(seq(numBifactors), function(x){
      length(cBfLoadings[[x]]) - sum(csLambda[, cBfWhichFactors[[x]]] != 0)
    }))
    cLambda <- matrix(0, nrow = (sum(extraIndicators) + nrow(csLambda)), ncol = (numBifactors + ncol(csLambda)))  
    # specific factors + covariate(s)
    cLambda[(sum(extraIndicators) + 1) : nrow(cLambda), (numBifactors + 1) : ncol(cLambda)] <- csLambda
    extraLoadings <- lapply(seq_along(cBfLoadings), function(x) cBfLoadings[[x]][1:extraIndicators[x]])
    for(i in 1:numBifactors){
      # bifactor extra indicators
      if(i == 1) idx <- 1 : extraIndicators[i] else idx <- (1 + extraIndicators[i - 1]) : (extraIndicators[i - 1] + extraIndicators[i])
      cLambda[idx, i] <- extraLoadings[[i]]
      # bifactor specific indicators 
      idx <- length(unlist(extraLoadings)) + unlist(lapply(cBfWhichFactors[[i]], function(x) which(csLambda[, x] != 0)))
      cLambda[idx, i] <- cBfLoadings[[i]][(1 + extraIndicators[i]) : length(cBfLoadings[[i]])]
    }
    cLambda
  })
  
  ### create Phi
  sPhi <- Phi
  nf <- 1:ncol(Lambda[[1]])
  posCov <- nf[!nf %in% c(1 : numBifactors, (numBifactors + unlist(bfWhichFactors[[1]])))]
  numCov <- length(posCov)
  Phi <- lapply(seq(nGroups), function(g){
    csPhi <- sPhi[[g]]
    if(ncol(csPhi) != (numCov + numBifactors)) stop('Incorrect dimensions for Phi. Remember that Phi may only refer to the bifactor(s) and the covariate(s), but must omit all specific factors.')
    cPhi <- diag(ncol(Lambda[[g]]))
    # bf-bf
    cPhi[1:numBifactors, 1:numBifactors] <-  csPhi[1:numBifactors, 1:numBifactors]
    # bf-cov
    if(numCov > 0){
      cPhi[posCov, 1:numBifactors] <- csPhi[(numBifactors + 1) : ncol(csPhi), 1 : numBifactors]
      cPhi[1:numBifactors, posCov] <- csPhi[1 : numBifactors, (numBifactors + 1) : ncol(csPhi)]
      # cov-cov
      cPhi[posCov, posCov] <- csPhi[(numBifactors + 1) : ncol(csPhi), (numBifactors + 1) : ncol(csPhi)]
    }
    cPhi
  })
  
  ### generate sigma 
  if(isMultigroup){
    generated <- semPower.genSigma(Lambda = Lambda, Phi = Phi)
    model <- generated[[1]][['modelTrueCFA']]
  }else{
    generated <- semPower.genSigma(Lambda = Lambda[[1]], Phi = Phi[[1]])
    model <- generated[['modelTrueCFA']] 
  }
  
  ### create model strings 
  # add bifactor orthogonalization constrains
  tok <- lapply(seq(numBifactors), function(x){
    ctok <- list()
    csf <- paste0('f', c(x, numBifactors + bfWhichFactors[[1]][[x]]))
    for(i in 1 : (length(csf) - 1)){
      for(j in (i + 1): length(csf)){
        ctok <- append(ctok, paste0(csf[i], ' ~~ 0*', csf[j]))
      }
    }
    ctok
  })
  model <- paste(c(model, unlist(tok)), collapse = '\n')
  
  ### H0 model
  # adapt nullWhich to account for specific factors
  nullWhich <- lapply(nullWhich, function(x){
    if(x[[1]] > numBifactors) x[[1]] <- posCov[(x[[1]] - numBifactors)]
    if(x[[2]] > numBifactors) x[[2]] <- posCov[(x[[2]] - numBifactors)]
    x
  })  
  if(nullEffect == 'cor=0'){
    modelH0 <- paste(c(model,
                       paste0('f', nullWhich[[1]], collapse = ' ~~ 0*')),
                     collapse = '\n')
  }else if(nullEffect == 'corx=corz'){
    labs <- list()
    tok <- ''
    for(i in seq_along(nullWhich)){
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
    modelH0 <- paste(c(model, tok), collapse = '\n')
  }else if(nullEffect == 'cora=corb'){
    if(is.null(nullWhichGroups)) nullWhichGroups <- seq(nGroups)
    lab <- paste0('ff', seq(nGroups))
    lab[nullWhichGroups] <- 'pf1'
    lab <- paste0('c(', paste(lab, collapse = ','), ')*')
    modelH0 <- paste(c(model,
                       paste0('f', nullWhich[[1]], collapse = paste0(' ~~ ', lab))),
                     collapse = '\n')
  }else{
    stop('nullEffect not defined')
  }
  
  # we always enforce invariance constraints in the multigroup case
  if(isMultigroup){
    args[['lavOptions']] <- append(args[['lavOptions']], list(group.equal = c('loadings', 'lv.variances')))
  } 
  
  modelH1 <- NULL
  fitH1model <- FALSE
  if(comparison == 'restricted'){
    modelH1 <- model 
    # single group case: the h1 model always fits perfectly
    # multigroup case: we cannot be sure that user input yields a perfectly fitting model
    fitH1model <- isMultigroup 
  } 
  
  if(isMultigroup) Sigma <- lapply(generated, '[[', 'Sigma') else Sigma <- generated[['Sigma']] 
  
  do.call(semPower.powerLav, append(list(
    type = type,
    Sigma = Sigma,
    modelH0 = modelH0,
    modelH1 = modelH1,
    fitH1model = fitH1model),
    args)
  )
  
}



#' semPower.powerAutoreg
#'
#' Convenience function for performing power analysis on effects in an autoregressive model.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param nWaves number of waves, must be >= 2.
#' @param autoregEffects vector of the autoregressive effects, e.g. `c(.7, .6)` for  autoregressive effects of .7 for `X1 -> X2` and .6 for `X2 -> X3`. Must be a list for multiple groups models.
#' @param lag1Effects alternative name for `autoregEffects`.
#' @param lag2Effects vector of lag-2 effects, e.g. `c(.2, .1)` for lag-2 effects of .2 for `X1 -> X3` and .1 for `X2 -> X4`.
#' @param lag3Effects vector of lag-3 effects, e.g. `c(.2)` for a lag-3 effect of .2 for `X1 -> X4`.
#' @param means vector of means for `X`. Can be omitted for no meanstructure.
#' @param variances vector of (residual-)variances for `X`. When omitted and `standardized = FALSE`, all (residual-)variances are equal to 1. When omitted and `standardized = TRUE`, the (residual-)variances are determined so that all variances are 1, and will thus typically differ from each other. When provided, `standardized` must be `FALSE`.
#' @param waveEqual parameters that are assumed to be equal across waves in both the H0 and the H1 model. Valid are `'lag1'` (or equivalently `'autoreg'`), `'lag2'`, and  `'lag3'`, or `NULL` for none (so that all parameters are freely estimated, subject to the constraints defined in `nullEffect`). 
#' @param nullEffect defines the hypothesis of interest. Valid are the same arguments as in `waveEqual` and additionally `'lag1 = 0'` (or equivalently `'autoregX = 0'`) `'lag2 = 0'`, `'lag3 = 0'` to constrain the autoregressive, lag-2, or lag-3 effects to zero, and `'autoregA = autoregB'` to the autoregressive effects be equal across groups. 
#' @param nullWhich used in conjunction with `nullEffect` to identify which parameter to constrain when there are > 2 waves and parameters are not constant across waves. For example, `nullEffect = 'lag1 = 0'` with `nullWhich = 2` would constrain the second autoregressive effect for X to zero.    
#' @param nullWhichGroups for hypothesis involving cross-groups comparisons, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param standardized whether all parameters should be standardized (`TRUE`, the default). If `FALSE`, all regression relations are unstandardized.
#' @param invariance whether metric invariance over waves is assumed (`TRUE`, the default) or not (`FALSE`). When means are part of the model, invariant intercepts are also assumed. This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model).
#' @param autocorResiduals whether the residuals of the indicators of latent variables are autocorrelated over waves (`TRUE`, the default) or not (`FALSE`). This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model).
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. The order of factors is (X1, X2, ..., X_nWaves). See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details
#' This function performs a power analysis to reject various hypotheses arising
#' in simple autoregressive (simplex) models, where one variable is repeatedly 
#' assessed at two or more different time points (`nWaves`), yielding 
#' autoregressive effects (aka lag-1 effects or stabilities, ; `X1 -> X2 -> X3`), and optionally 
#' lagged effects (`X1 ->  X3`), variances, and means.
#'  
#' Relevant hypotheses in arising in an autogressive model are:
#' * `autoreg` or `lag1`: Tests the hypothesis that the autoregressive (lag-1) effects are equal across waves (stationarity of autoregressive parameters).
#' * `lag2`: Tests the hypothesis that the lag-2 effects are equal across waves (stationarity of lag-2 effects).
#' * `lag3`: Tests the hypothesis that the lag-3 effects are equal across waves (stationarity of lag-3 effects).
#' * `var`: Tests the hypothesis that the residual-variances of X (i.e., X_2, ..., X_nWaves) are equal across waves (stationarity of variance).
#' * `mean`: Tests the hypothesis that the conditional means of X (i.e., X_2, ..., X_nWaves) are equal across waves (stationarity of means).
#' * `autoreg = 0` or `lag1 = 0`: Tests the hypothesis that the autoregressive (lag-1) effect of X is zero. 
#' * `lag2 = 0` and `lag3 = 0`: Tests the hypothesis that a lag-2 or a lag-3 effect is zero.
#' * `autoregA = autoregB` or `lag1A = lag1B`: : Tests the hypothesis that the autoregressive effect of X is equal across groups.
#' 
#' For hypotheses in an ARMA model, see [semPower.powerARMA()]. For hypotheses regarding a CLPM structure, see [semPower.powerCLPM()]. For hypotheses regarding longitudinal measurement invariance, see [semPower.powerLI()].
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model**:
#' * `Lambda`: The factor loading matrix (with the number of columns equaling the number of factors).
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading on the second factor by .8, .6, .6, and .4.
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the loadings of the first factor to equal .5 and those of the second factor do equal .6.
#' 
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` need to be defined.
#' If the model contains observed variables only, use `Lambda = diag(x)` where `x` is the number of variables.
#'
#' Note that the order of the factors is (X1, X2, ..., X_nWaves).
#' 
#' Additional arguments related to the requested type of **power analysis**:
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.  This requires the `semTools` package.
#' 
#' @examples
#' \dontrun{
#' # Determine required N in a 4-wave autoregressive model
#' # to detect an autoregressive effect between X1 -> X2 of >= .5
#' # with a power of 80% on alpha = 5%, where
#' # X is measured by 3 indicators loading by .5 each (at each wave), and 
#' # the autoregressive effecst are .5 (X1 -> X2), .7 (X2 -> X3), and .6 (X3 -> X4), and
#' # there are no lagged effects, and
#' # metric invariance and autocorrelated residuals are assumed
#' powerAutoreg <- semPower.powerAutoreg(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 4, 
#'   autoregEffects = c(.5, .7, .6),
#'   nullEffect = 'autoreg=0',
#'   nullWhich = 1,
#'   nIndicator = rep(3, 4), loadM = .5,
#'   invariance = TRUE, autocorResiduals = TRUE)
#' 
#' # show summary
#' summary(powerAutoreg)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerAutoreg$modelH1, sample.cov = powerAutoreg$Sigma,
#'             sample.nobs = powerAutoreg$requiredN,
#'             sample.cov.rescale = FALSE)
#' lavaan::sem(powerAutoreg$modelH0, sample.cov = powerAutoreg$Sigma,
#'             sample.nobs = powerAutoreg$requiredN,
#'             sample.cov.rescale = FALSE)
#' 
#' 
#' # same as above, but determine power with N = 250 on alpha = .05
#' powerAutoreg <- semPower.powerAutoreg(
#'   'post-hoc', alpha = .05, N = 250,
#'   nWaves = 4, 
#'   autoregEffects = c(.5, .7, .6),
#'   nullEffect = 'autoreg=0',
#'   nullWhich = 1,
#'   nIndicator = rep(3, 4), loadM = .5,
#'   invariance = TRUE, autocorResiduals = TRUE)
#' 
#' # same as above, but determine the critical chi-square with N = 250 so that alpha = beta
#' powerAutoreg <- semPower.powerAutoreg(
#'   'compromise', abratio = 1, N = 250,
#'   nWaves = 4, 
#'   autoregEffects = c(.5, .7, .6),
#'   nullEffect = 'autoreg=0',
#'   nullWhich = 1,
#'   nIndicator = rep(3, 4), loadM = .5,
#'   invariance = TRUE, autocorResiduals = TRUE)
#' 
#' # same as above, but compare to the saturated model
#' # (rather than to the less restricted model)
#' powerAutoreg <- semPower.powerAutoreg(
#'   'post-hoc', alpha = .05, N = 250,
#'   comparison = 'saturated',
#'   nWaves = 4, 
#'   autoregEffects = c(.5, .7, .6),
#'   nullEffect = 'autoreg=0',
#'   nullWhich = 1,
#'   nIndicator = rep(3, 4), loadM = .5,
#'   invariance = TRUE, autocorResiduals = TRUE)
#' 
#' # same as above, but assume only observed variables
#' powerAutoreg <- semPower.powerAutoreg(
#'   'post-hoc', alpha = .05, N = 250,
#'   nWaves = 4, 
#'   autoregEffects = c(.5, .7, .6),
#'   nullEffect = 'autoreg=0',
#'   nullWhich = 1,
#'   Lambda = diag(4))
#' 
#' # same as above, but provide reduced loadings matrix to define that
#' # X is measured by 3 indicators each loading by .8, .6, .7 (at each wave)
#' powerAutoreg <- semPower.powerAutoreg(
#'   'post-hoc', alpha = .05, N = 250,
#'   nWaves = 4, 
#'   autoregEffects = c(.5, .7, .6),
#'   nullEffect = 'autoreg=0',
#'   nullWhich = 1,
#'   loadings = list(
#'     c(.8, .6, .7),   # X1
#'     c(.8, .6, .7),   # X2
#'     c(.8, .6, .7),   # X3
#'     c(.8, .6, .7)    # X4
#'   ), 
#'   invariance = TRUE, autocorResiduals = TRUE)
#' 
#' # same as above, but assume wave-constant autoregressive effects
#' powerAutoreg <- semPower.powerAutoreg(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 4, 
#'   autoregEffects = c(.6, .6, .6),
#'   waveEqual = c('autoreg'),
#'   nullEffect = 'autoreg=0',
#'   nIndicator = rep(3, 4), loadM = .5,
#'   invariance = TRUE, autocorResiduals = TRUE)
#' 
#' 
#' # same as above, but detect that autoregressive effects are not wave-constant
#' powerAutoreg <- semPower.powerAutoreg(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 4, 
#'   autoregEffects = c(.6, .7, .8),
#'   nullEffect = 'autoreg',
#'   nIndicator = rep(3, 4), loadM = .5,
#'   invariance = TRUE, autocorResiduals = TRUE)
#' 
#' # same as above, but include lag-2 and lag-3 effects
#' powerAutoreg <- semPower.powerAutoreg(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 4, 
#'   autoregEffects = c(.6, .6, .6),
#'   lag2Effects = c(.25, .20),
#'   lag3Effects = c(.15),
#'   waveEqual = c('autoreg'),
#'   nullEffect = 'autoreg=0',
#'   nIndicator = rep(3, 4), loadM = .5,
#'   invariance = TRUE, autocorResiduals = TRUE)
#' 
#' 
#' # same as above, but detect that first lag-2 effect differs from zero
#' powerAutoreg <- semPower.powerAutoreg(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 4, 
#'   autoregEffects = c(.6, .6, .6),
#'   lag2Effects = c(.25, .20),
#'   lag3Effects = c(.15),
#'   waveEqual = c('autoreg'),
#'   nullEffect = 'lag2=0',
#'   nullWhich = 1,
#'   nIndicator = rep(3, 4), loadM = .5,
#'   invariance = TRUE, autocorResiduals = TRUE)
#' 
#' 
#' # same as above, but assume wave-constant lag2 effects
#' powerAutoreg <- semPower.powerAutoreg(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 4, 
#'   autoregEffects = c(.6, .6, .6),
#'   lag2Effects = c(.25, .25),
#'   lag3Effects = c(.15),
#'   waveEqual = c('autoreg', 'lag2'),
#'   nullEffect = 'lag2=0',
#'   nIndicator = rep(3, 4), loadM = .5,
#'   invariance = TRUE, autocorResiduals = TRUE)
#' 
#' 
#' # same as above, but detect that lag3 effect differs from zero
#' powerAutoreg <- semPower.powerAutoreg(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 4, 
#'   autoregEffects = c(.6, .6, .6),
#'   lag2Effects = c(.25, .25),
#'   lag3Effects = c(.15),
#'   waveEqual = c('autoreg', 'lag2'),
#'   nullEffect = 'lag3=0',
#'   nIndicator = rep(3, 4), loadM = .5,
#'   invariance = TRUE, autocorResiduals = TRUE)
#' 
#' 
#' # Determine required N in a 3-wave autoregressive model
#' # assuming wave-constant autoregressive effects 
#' # that the autoregressive effects in group 1
#' # differ from those in group 2
#' # with a power of 80% on alpha = 5%, where
#' # X is measured by 3 indicators loading by .5 each (at each wave and in each group), and 
#' # the autoregressive effect is .7 in group 1 and
#' # the autoregressive effect is .5 in group 2 and
#' # there are no lagged effects, and
#' # metric invariance over both time and groups and autocorrelated residuals are assumed and
#' # the groups are equal-sized
#' powerAutoreg <- semPower.powerAutoreg(
#'   'a-priori', alpha = .05, power = .80, N = list(1, 1),
#'   nWaves = 3, 
#'   autoregEffects = list(
#'     c(.7, .7),
#'     c(.5, .5)
#'   ),
#'   waveEqual = c('autoreg'),
#'   nullEffect = 'autoregA = autoregB',
#'   nullWhich = 1,
#'   nIndicator = rep(3, 3), loadM = .5,
#'   invariance = TRUE, autocorResiduals = TRUE)
#'   
#' # Determine required N in a 4-wave autoregressive model
#' # to detect that the factor residual-variances (X2, X3, X4) differ
#' # with a power of 80% on alpha = 5%, where
#' # the (residual-)variances are 1, .5, 1.5, and 1, respectively,  
#' # X is measured by 3 indicators loading by .5 each (at each wave), and
#' # the autoregressive effects are .6, and
#' # both the H0 and the H1 assume wave-constant autoregressive effects, and
#' # there are no lagged effects, and
#' # metric invariance and autocorrelated residuals are assumed
#' powerAutoreg <- semPower.powerAutoreg(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 4,
#'   autoregEffects = c(.6, .6, .6),
#'   variances = c(1, .5, 1.5, 1),
#'   waveEqual = c('autoreg'),
#'   nullEffect = 'var',
#'   nullWhich = 1,
#'   nIndicator = rep(3, 4), loadM = .5,
#'   standardized = FALSE,
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE)
#' 
#' # same as above, but 
#' # include latent means and 
#' # detect that latent means differ and
#' # assume wave-constant variances and autoregressive parameters for both H0 and H1
#' powerAutoreg <- semPower.powerAutoreg(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 4,
#'   autoregEffects = c(.6, .6, .6),
#'   variances = c(1, 1, 1, 1),
#'   means = c(0, .5, 1, .7),
#'   waveEqual = c('autoreg', 'var'),
#'   nullEffect = 'mean',
#'   nullWhich = 1,
#'   nIndicator = rep(3, 4), loadM = .5,
#'   standardized = FALSE,
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE)
#'   
#' # request a simulated post-hoc power analysis with 500 replications
#' set.seed(300121)
#' powerAutoreg <- semPower.powerAutoreg(
#'   'post-hoc', alpha = .05, N = 500,
#'   nWaves = 3, 
#'   autoregEffects = c(.7, .7),
#'   waveEqual = c('autoreg'),
#'   nullEffect = 'autoreg = 0',
#'   nullWhich = 1,
#'   nIndicator = rep(3, 3), loadM = .5,
#'   invariance = TRUE, autocorResiduals = TRUE, 
#'   simulatedPower = TRUE,
#'   simOptions = list(nReplications = 500)
#'   )
#'   
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerAutoreg <- function(type, comparison = 'restricted',
                                  nWaves = NULL, 
                                  autoregEffects = NULL, 
                                  lag1Effects = autoregEffects, 
                                  lag2Effects = NULL, 
                                  lag3Effects = NULL,
                                  means = NULL,
                                  variances = NULL,
                                  waveEqual = NULL, 
                                  nullEffect = NULL, 
                                  nullWhich = NULL,
                                  nullWhichGroups = NULL,
                                  standardized = TRUE,
                                  invariance = TRUE,
                                  autocorResiduals = TRUE,
                                  ...){
  
  comparison <- checkComparisonModel(comparison)
  checkEllipsis(...)
  
  # we override Beta and Sigma later, so let's make sure it is not set in ellipsis argument
  if('Beta' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Beta.')
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma.')
  
  # validate input
  if(is.null(autoregEffects)) autoregEffects <- lag1Effects
  if(is.null(autoregEffects)) stop('autoregEffects / lag1Effects may not be NULL.')
  if(is.null(nWaves) || is.na(nWaves) || nWaves < 2) stop('nWaves must be >= 2.')
  
  nullValid <- c('autoreg', 'autoreg=0', 'autorega=autoregb',
                 'lag1', 'lag1=0', 'lag1a=lag1b',
                 'lag2', 'lag3', 'lag2=0', 'lag3=0', 'lagged=0',
                 'var', 'mean')
  nullEffect <- checkNullEffect(nullEffect, nullValid)
  
  # [[groups]][[waves]]
  isMultigroup <- is.list(autoregEffects)
  if(!is.list(autoregEffects)) autoregEffects <- list(autoregEffects)
  nGroups <- length(autoregEffects)
  
  if(isMultigroup && !nullEffect %in% c('autorega=autoregb')) stop('Multigroup analysis are only supported for nullEffect = autoregA=autoregB')
  if(!isMultigroup && nullEffect %in% c('autorega=autoregb')) stop('nullEffect = autorega=autoregb imply multigroup analyses, but no list structure for any relevant parameter provided.')
  if(isMultigroup && is.null(nullWhichGroups)) nullWhichGroups <- seq(nGroups)
  
  if(any(unlist(lapply(autoregEffects, function(x) length(x) != (nWaves - 1))))) stop('autoregEffects must be of length nWaves - 1.')
  invisible(lapply(autoregEffects, function(x) lapply(x, function(x) checkBounded(x, 'All autoregressive effects ', bound = c(-1, 1), inclusive = FALSE))))
  
  if(!is.null(lag2Effects) && nWaves < 3) stop('There are no lag-2 effects when nWaves < 3')
  if(!is.null(lag3Effects) && nWaves < 4) stop('There are no lag-2 effects when nWaves < 4')
  if(isMultigroup && !is.null(lag2Effects) && !is.list(lag2Effects)) stop('For multigroup models, lag2Effects must be a list')
  if(isMultigroup && !is.null(lag3Effects) && !is.list(lag3Effects)) stop('For multigroup models, lag3Effects must be a list')
  if(isMultigroup && !is.null(lag2Effects) && length(lag2Effects) != nGroups) stop('lag2Effects must be provided for each group.')
  if(isMultigroup && !is.null(lag3Effects) && length(lag3Effects) != nGroups) stop('lag3Effects must be provided for each group.')

  estimateLag2Effects <- !is.null(lag2Effects)
  if(nWaves > 2){
    if(is.null(lag2Effects)) lag2Effects <- rep(0, nWaves - 2)
    if(!is.list(lag2Effects)) lag2Effects <- rep(list(lag2Effects), nGroups)
    if(any(unlist(lapply(lag2Effects, function(x) length(x) != (nWaves - 2))))) stop('lag2Effects must be of length nWaves - 2.')
    invisible(lapply(lag2Effects, function(x) lapply(x, function(x) checkBounded(x, 'All lag2 effects ', bound = c(-1, 1), inclusive = FALSE))))
  }
  estimateLag3Effects <- !is.null(lag3Effects)
  if(nWaves > 3){
    if(is.null(lag3Effects)) lag3Effects <- rep(0, nWaves - 3)
    if(!is.list(lag3Effects)) lag3Effects <- rep(list(lag3Effects), nGroups)
    if(any(unlist(lapply(lag3Effects, function(x) length(x) != (nWaves - 3))))) stop('lag3Effects must be of length nWaves - 3.')
    invisible(lapply(lag3Effects, function(x) lapply(x, function(x) checkBounded(x, 'All lag3 effects ', bound = c(-1, 1), inclusive = FALSE))))
  }
  if(!is.null(variances) && standardized) stop('When variances are provided, standardized must be FALSE.') 
  if(isMultigroup && !is.null(means) && !is.list(means)) stop('For multigroup models, means must be a list')
  if(isMultigroup && !is.null(variances) && !is.list(variances)) stop('For multigroup models, means must be a list')
  if(!is.null(means) && !is.list(means)) means <- list(means)
  if(!is.null(variances) && !is.list(variances)) variances <- list(variances)
  if(!is.null(means) && any(unlist(lapply(means, function(x) length(x) != nWaves)))) stop('means must be provided for each wave.') 
  if(!is.null(variances) && any(unlist(lapply(variances, function(x) length(x) != nWaves)))) stop('variances must be provided for each wave.') 
  
  if(!is.null(waveEqual)){
    waveEqual <- unlist(lapply(waveEqual, function(x) tolower(trimws(x))))
    if(any(unlist(lapply(waveEqual, function(x) !x %in% c('autoreg', 'lag1', 'lag2', 'lag3', 'var', 'mean'))))) stop('waveEqual may only contain autoreg, lag1, lag2, lag3, mean, or var.')
  }else{
    waveEqual <- ''
  }
  
  if(nullEffect == 'lag1') nullEffect <- 'autoreg'
  if(nullEffect == 'lag1=0') nullEffect <- 'autoreg=0'
  if(nullEffect == 'lag1a=lag1b') nullEffect <- 'autorega=autoregb'
  if('lag1' %in% waveEqual) waveEqual[which('lag1' %in% waveEqual)] <- 'autoreg'
  if('lag1=0' %in% waveEqual) waveEqual[which('lag1' %in% waveEqual)] <- 'autoreg=0'
  if('lag1a=lag1b' %in% waveEqual) waveEqual[which('lag1' %in% waveEqual)] <- 'autorega=autoregb'

  if(nullEffect %in% waveEqual) stop('You cannot set the same parameters in nullEffect and waveEqual')
  if(nWaves == 2 && nullEffect %in% c('autoreg')) stop('For two waves, there is only one one autoregressive effect. Did you mean autoreg = 0?')
  if(nWaves == 3 && nullEffect %in% c('lag2')) stop('For three waves, there is only one lag2 effect. Did you mean lag2 = 0?')
  if(nWaves == 4 && nullEffect %in% c('lag3')) stop('For three waves, there is only one lag3 effect. Did you mean lag3 = 0?')
  if((nullEffect == 'var' || 'var' %in% waveEqual) && is.null(variances)) stop('Either nullEffect or waveEqual refer to variances, but no variances are provided.')
  if((nullEffect == 'mean' || 'mean' %in% waveEqual) && is.null(variances)) stop('Either nullEffect or waveEqual refer to means, but no means are provided.')
  if((nullEffect == 'mean' || 'mean' %in% waveEqual) && !invariance) stop('When hypothesis include latent means, invariance must be TRUE.')
  
  if(is.null(nullWhich) && nWaves == 2) nullWhich <- 1
  if(is.null(nullWhich)){
    if(nWaves > 2 && !c('autoreg') %in% waveEqual && nullEffect %in% c('autoreg=0', 'autorega=autoregb')) stop('nullWhich must be defined when there are more than 2 waves and the autoregressive effects are not constant across waves.') 
    if(nWaves > 3 && !c('lag2') %in% waveEqual && nullEffect %in% c('lag2=0', 'lagged=0')) stop('nullWhich must be defined when there are more than 3 waves and the lag2 effects are not constant across waves.') 
    if(nWaves > 4 && !c('lag3') %in% waveEqual && nullEffect %in% c('lag3=0', 'lagged=0')) stop('nullWhich must be defined when there are more than 4 waves and the lag3 effects are not constant across waves.') 
    nullWhich <- 1 # this should be the proper default for all remaining cases
  }
  if(!is.null(nullWhich)){
    if(!is.numeric(nullWhich) || length(nullWhich) > 1) stop('nullWhich must be a single number.')
    if(nullWhich < 1 || nullWhich > (nWaves - 1)) stop('nullWhich must lie between 1 and nWaves - 1.')
    if('lag2' %in% nullEffect && nullWhich > (nWaves - 2)) stop('For lag2 effects, nullWhich must lie between 1 and nWaves - 2.')
    if('lag3' %in% nullEffect && nullWhich > (nWaves - 3)) stop('For lag3 effects, nullWhich must lie between 1 and nWaves - 3.')
  }
  
  ### create B
  Beta <- lapply(seq(nGroups), function(x){
    B <- matrix(0, ncol = nWaves, nrow = nWaves)
    # add autoregEffects
    for(i in seq(nWaves - 1)){
      B[(i + 1), i] <- autoregEffects[[x]][i]
    }
    # lag-2 effects
    if(nWaves > 2){
      for(i in seq(nWaves - 2)){
        B[(i + 2), i] <- lag2Effects[[x]][i]
      }
    }
    # lag-3 effects
    if(nWaves > 3){
      for(i in seq(nWaves - 3)){
        B[(i + 3), i] <- lag3Effects[[x]][i]
      }
    }
    B
  })
  
  # add metric invariance constraints to analysis model
  metricInvarianceFactors <- NULL
  if(invariance) metricInvarianceFactors <- list(seq(nWaves))

  # define latent means
  Alpha <- NULL
  if(!is.null(means)){
    Alpha <- means
  }
  
  ### get Sigma
  if(standardized){
    Phi <- lapply(seq_along(Beta), function(x) getPhi.B(Beta[[x]]))
    generated <- semPower.genSigma(Phi = if(!isMultigroup) Phi[[1]] else Phi, 
                                   Alpha = if(!isMultigroup) Alpha[[1]] else Alpha,
                                   useReferenceIndicator = TRUE,
                                   metricInvariance = metricInvarianceFactors, 
                                   nGroups = nGroups,
                                   ...)
  }else{
    # Psi
    Psi <- NULL
    if(!is.null(variances)){
      Psi <- lapply(seq(variances), function(x){
        diag(variances[[x]])
      })
    }
    generated <- semPower.genSigma(Beta = if(!isMultigroup) Beta[[1]] else Beta,
                                   Psi = if(!isMultigroup) Psi[[1]] else Psi,
                                   Alpha = if(!isMultigroup) Alpha[[1]] else Alpha,
                                   useReferenceIndicator = TRUE,
                                   metricInvariance = metricInvarianceFactors, 
                                   nGroups = nGroups,
                                   ...)
  }
  
  ### create model strings
  if(!isMultigroup) model <- generated[['modelTrueCFA']] else model <- generated[[1]][['modelTrueCFA']]
  if(!isMultigroup) Lambda <- generated[['Lambda']] else Lambda <- generated[[1]][['Lambda']]
  
  # add autoregressive structure 
  for(f in 2:ncol(Beta[[1]])){     # omit first row
    fidx <- (f - 1)
    if(estimateLag2Effects && f > 2) fidx <- c(fidx, (f - 2)) # estm lag2 effects regardless of these are zero
    if(estimateLag3Effects && f > 3) fidx <- c(fidx, (f - 3)) # estm lag3 effects regardless of these are zero
    if(length(fidx) != 0){
      tok <- paste0('f', f, ' ~ ', paste(paste0('pf', paste0(formatC(f, width = 2, flag = 0), formatC(fidx, width = 2, flag = 0)), '*'), paste0('f', fidx), sep = '', collapse = ' + '))
      model <- paste(model, tok, sep='\n')
    }
  }

  # add variances
  if(!is.null(variances)){
    for(f in 1:ncol(Beta[[1]])){
      tok <- paste0('f', f, ' ~~ ', 'pvf', formatC(f, width = 2, flag = 0), '*f', f)
      model <- paste(model, tok, sep='\n')
    }
  }
  
  # estm means + intercepts
  if(!is.null(means) && invariance){
    # fix first factor mean to zero
    model <- paste(model, 'f1 ~ 0*1', sep='\n')
    for(f in 2:ncol(Beta[[1]])){  # start at 2
      tok <- paste0('f', f, ' ~ ', 'pmf', formatC(f, width = 2, flag = 0), '*1')
      model <- paste(model, tok, sep='\n')
    }
    # add invariance constraints on intercepts
    tok <- list()
    for(x in seq_along(metricInvarianceFactors)){
      ci <- lapply(metricInvarianceFactors[[x]], function(f) paste0('x', which(Lambda[, f] != 0)))
      lab <- paste0('pmi', 1:length(ci[[1]]))
      tok <- append(tok, unlist(lapply(ci, function(x) paste0(x, ' ~ ', lab, '*1') )))
    }
    model <- paste(model, paste(unlist(tok), collapse = '\n'), sep='\n')
  }
  
  # add autocorrelated residuals
  if(autocorResiduals){
    # do this only when there is at least one latent variable
    if(nrow(Lambda) > nWaves){
      autocorResidualsFactors <- metricInvarianceFactors  # same structure 
      tok <- list()
      for(x in seq_along(autocorResidualsFactors)){
        ci <- lapply(autocorResidualsFactors[[x]], function(f) paste0('x', which(Lambda[, f] != 0)))
        if(length(ci[[1]]) > 1){
          for(i in 1:(length(ci) - 1)){
            for(j in (i + 1) : length(ci)){
              tok <- append(tok, paste(ci[[i]], '~~', ci[[j]]))
            }
          }
        }
      }
      model <- paste(model, paste(tok, collapse = '\n'), sep ='\n')
    }
  }
  
  
  ### define H1 and H0 model
  
  # first get relevant parameter labels that may be subject to waveEqual or nullEffect constraints 
  xw <- seq(nWaves, 2)
  pAutoregX <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(xw - 1, width = 2, flag = 0))
  pAutoregX <- pAutoregX[order(pAutoregX)]
  if(nWaves > 2){
    xw <- seq(nWaves, 3)
    pLag2 <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(xw - 2, width = 2, flag = 0))
    pLag2 <- pLag2[order(pLag2)]
  }
  if(nWaves > 3){
    xw <- seq(nWaves, 4)
    pLag3 <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(xw - 3, width = 2, flag = 0))
    pLag3 <- pLag3[order(pLag3)]
  }
  pVariances <- paste0('pvf', formatC(2:ncol(Beta[[1]]), width = 2, flag = 0)) # first factor var is no residual var
  pMeans <- paste0('pmf', formatC(2:ncol(Beta[[1]]), width = 2, flag = 0)) # first factor mean is zero
  
  
  # multigroup case
  if(isMultigroup){
    # remove group specific labels from measurement part to enforce metric invariance
    model <- gsub(paste0('_g', seq(nGroups), collapse = '|'), '_gc', model)
    # assign group labels to all structural parameters
    for(pp in seq(nWaves - 1)){
      patt <- pAutoregX[pp]
      repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      model <- gsub(patt, repl, model)
    }
    if(nWaves > 2){
      for(pp in seq(pLag2)){
        patt <- pLag2[pp]
        repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
        model <- gsub(patt, repl, model)
      }
    }
    if(nWaves > 3){
      for(pp in seq(pLag3)){
        patt <- pLag3[pp]
        repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
        model <- gsub(patt, repl, model)
      }
    }
  }
  
  # add constraints to H1 model
  modelH1 <- model
  if(c('autoreg') %in% waveEqual) modelH1 <- gsub(paste(pAutoregX, collapse = '|'), paste(pAutoregX, collapse = ''), modelH1)
  if(c('lag2') %in% waveEqual) modelH1 <- gsub(paste(pLag2, collapse = '|'), paste(pLag2, collapse = ''), modelH1)
  if(c('lag3') %in% waveEqual) modelH1 <- gsub(paste(pLag3, collapse = '|'), paste(pLag3, collapse = ''), modelH1)
  if(c('mean') %in% waveEqual) modelH1 <- gsub(paste(pMeans, collapse = '|'), paste(pMeans, collapse = ''), modelH1)
  if(c('var') %in% waveEqual) modelH1 <- gsub(paste(pVariances, collapse = '|'), paste(pVariances, collapse = ''), modelH1)
  
  # add additional constraints to H0 model
  modelH0 <- modelH1
  
  # wave equal constraints not included in modelH1:
  if(nullEffect %in% 'autoreg') modelH0 <- gsub(paste(pAutoregX, collapse = '|'), paste(pAutoregX, collapse = ''), modelH0)
  if(nullEffect %in% 'lag2') modelH0 <- gsub(paste(pLag2, collapse = '|'), paste(pLag2, collapse = ''), modelH0)
  if(nullEffect %in% 'lag3') modelH0 <- gsub(paste(pLag3, collapse = '|'), paste(pLag3, collapse = ''), modelH0)
  if(nullEffect %in% 'mean') modelH0 <- gsub(paste(pMeans, collapse = '|'), paste(pMeans, collapse = ''), modelH0)
  if(nullEffect %in% 'var') modelH0 <- gsub(paste(pVariances, collapse = '|'), paste(pVariances, collapse = ''), modelH0)
  
  # zero and equality constraints:
  if(nullEffect %in% 'autoreg=0'){
    if('autoreg' %in% waveEqual){
      modelH0 <- gsub(paste(pAutoregX, collapse = ''), '0', modelH0)
    }else{
      modelH0 <- gsub(pAutoregX[nullWhich], '0', modelH0)
    }
  }
  if(nullEffect %in% c('lag2=0', 'lagged')){
    if('lag2' %in% waveEqual){
      modelH0 <- gsub(paste(pLag2, collapse = ''), '0', modelH0)
    }else{
      modelH0 <- gsub(pLag2[nullWhich], '0', modelH0)
    }
  }
  if(nullEffect %in% c('lag3=0', 'lagged')){
    if('lag3' %in% waveEqual){
      modelH0 <- gsub(paste(pLag3, collapse = ''), '0', modelH0)
    }else{
      modelH0 <- gsub(pLag3[nullWhich], '0', modelH0)
    }
  }
  
  # multigroup cases
  if(nullEffect %in% 'autorega=autoregb'){
    if(c('autoreg') %in% waveEqual){
      patt <- paste0(paste(pAutoregX, collapse = ''), '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(paste(pAutoregX, collapse = ''), '_gc')
    }else{
      patt <- paste0(pAutoregX[nullWhich], '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(pAutoregX[nullWhich], '_gc')
    }
    modelH0 <- gsub(patt, repl, modelH0)
  }
  
  # here we actually fit modelH1 in case of a restricted comparison
  # because we cannot be sure that user input yields perfectly fitting h1 models, 
  # when there are additional constraints (waveequal or invariance)
  if(comparison == 'saturated') modelH1 <- NULL
  
  if(isMultigroup) Sigma <- lapply(generated, '[[', 'Sigma') else Sigma <- generated[['Sigma']] 
  mu <- NULL
  if(!is.null(means)){
    if(isMultigroup) mu <- lapply(generated, '[[', 'mu') else mu <- generated[['mu']] 
  }
  
  semPower.powerLav(type, 
                    Sigma = Sigma,
                    mu = mu,
                    modelH0 = modelH0, 
                    modelH1 = modelH1, 
                    ...)
}




#' semPower.powerARMA
#'
#' Convenience function for performing power analysis on effects in an ARMA model.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param nWaves number of waves, must be >= 2.
#' @param autoregEffects vector of the lag-1 autoregressive effects, e.g. `c(.7, .6)` for  autoregressive effects of .7 for `X1 -> X2` and .6 for `X2 -> X3`. Must be a list for multiple groups models.
#' @param autoregLag1 alternative name for autoregEffects.
#' @param autoregLag2 vector of lag-2 effects, e.g. `c(.2, .1)` for lag-2 effects of .2 for `X1 -> X3` and .1 for `X2 -> X4`.
#' @param autoregLag3 vector of lag-3 effects, e.g. `c(.2)` for a lag-3 effect of .2 for `X1 -> X4`.
#' @param mvAvgLag1 vector of the lag-1 moving average parameters, e.g. `c(.4, .3)` for moving average parameters of .4 for `N1 -> X2` and .3 for `N2 -> X3`. Must be a list for multiple groups models.
#' @param mvAvgLag2 vector of the lag-2 moving average parameters, e.g. `c(.3, .2)` for moving average parameters effects of .2 for `N1 -> X3` and .2 for `N2 -> X4`. Must be a list for multiple groups models.
#' @param mvAvgLag3 vector of the lag-3  moving average parameters, e.g. `c(.2)` for a moving average parameter of .2 for `N1 -> X4`. Must be a list for multiple groups models.
#' @param means vector of means of `X`. May be `NULL` for no meanstructure.
#' @param variances vector of variances of the noise factors `N` (= residual variances of `X`).
#' @param waveEqual parameters that are assumed to be equal across waves in both the H0 and the H1 model. Because ARMA models are likely not identified when no such constraints are imposed, this may not be empty. Valid are `'autoreg'`, `'autoregLag2'`, and  `'autoregLag3'` for autoregressive effects, `'mvAvg'`, `'mvAvgLag2'`, and  `'mvAvgLag3'` for moving average effects, `var` for the variance of the noise factors (starting at wave 2), `mean` for the conditional means of X  (starting at wave 2).
#' @param groupEqual parameters that are restricted across groups in both the H0 and the H1 model, when `nullEffect` implies a multiple group model. Valid are `autoreg`, `mvAvg`, `var`, `mean`.
#' @param nullEffect defines the hypothesis of interest. Valid are the same arguments as in `waveEqual` and additionally `'autoreg = 0'`, `'autoregLag2 = 0'`, `'autoregLag3 = 0'`, `'mvAvg = 0'`, `'mvAvgLag2 = 0'`, `'mvAvgLag3 = 0'`,  to constrain the autoregressive or moving average effects to zero, and `'autoregA = autoregB'`, `'mvAvgA = mvAvgB'`, `'varA = varB'`, `'meanA = meanB'` to constrain the autoregressive (lag-1) effects, moving average (lag-1) effects, variances of the noise factors, or means of the X to be equal across groups. 
#' @param nullWhich used in conjunction with `nullEffect` to identify which parameter to constrain when there are multiple waves and parameters are not constant across waves. For example, `nullEffect = 'autoreg = 0'` with `nullWhich = 2` would constrain the second autoregressive effect for X to zero.    
#' @param nullWhichGroups for hypothesis involving cross-groups comparisons, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param invariance whether metric invariance over waves is assumed (`TRUE`, the default) or not (`FALSE`). When means are part of the model, invariant intercepts are also assumed. This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model).
#' @param autocorResiduals whether the residuals of the indicators of latent variables are autocorrelated over waves (`TRUE`, the default) or not (`FALSE`). This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model).
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. The order of factors is (X1, X2, ..., X_nWaves). See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details
#' This function performs a power analysis to reject various hypotheses arising
#' in models with autoregressive and moving average parameters (ARMA models), where one variable `X` is repeatedly 
#' assessed at different time points (`nWaves`), and autoregressive (lag-1 effects; `X1 -> X2 -> X3`, and optionally lag-2 and lag-3) effects, 
#' and moving average parameters (`N1 ->  X2`, or equivalently for lag-2 and lag-3 effects) are assumed.
#'  
#' Relevant hypotheses in arising in an ARMA model are:
#' * `autoreg`: Tests the hypothesis that the autoregressive lag-1 effects are equal across waves (stationarity of autoregressive lag-1 effects).
#' * `autoregLag2`: Tests the hypothesis that the autoregressive lag-2 effects are equal across waves (stationarity of autoregressive lag-2 effects).
#' * `autoregLag3`: Tests the hypothesis that the autoregressive lag-3 effects are equal across waves (stationarity of autoregressive lag-3 effects).
#' * `mvAvg`: Tests the hypothesis that the moving average lag-1 parameters are equal across waves (stationarity of moving average lag-1 effects).
#' * `mvAvgLag2`: Tests the hypothesis that the moving average lag-2 parameters are equal across waves (stationarity of moving average lag-2 effects).
#' * `mvAvgLag3`: Tests the hypothesis that the moving average lag-3 parameters are equal across waves (stationarity of moving average lag-3 effects).
#' * `var`: Tests the hypothesis that the variances of the noise factors N (= the residual variances of X) are equal across waves 2 to nWaves (stationarity of variance).
#' * `mean`: Tests the hypothesis that the conditional means of X are equal across waves 2 to nWaves (stationarity of means).
#' * `autoreg = 0`, `autoregLag2 = 0`, `autoregLag3 = 0`: Tests the hypothesis that the autoregressive effects of the specified lag is zero. 
#' * `mvAvg = 0`, `mvAvgLag2 = 0`, `mvAvgLag3 = 0`: Tests the hypothesis that the moving average parameter of the specified lag is zero. 
#' * `autoregA = autoregB`: Tests the hypothesis that the autoregressive lag-1 effect is equal across groups.
#' * `mvAvgA = mvAvgB`: Tests the hypothesis that the moving average lag-1 parameter is equal across groups.
#' * `varA = varB`: Tests the hypothesis that the variance of the noise factors are equal across groups.
#' * `meanA = meanB`: Tests the hypothesis that latent means are equal across groups.
#' 
#' For hypotheses regarding a simple autoregression, see [semPower.powerAutoregressive()]. For hypotheses regarding a CLPM structure, see [semPower.powerCLPM()].  For hypotheses regarding longitudinal measurement invariance, see [semPower.powerLI()].
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model**:
#' * `Lambda`: The factor loading matrix (with the number of columns equaling the number of factors).
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading on the second factor by .8, .6, .6, and .4.
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the loadings of the first factor to equal .5 and those of the second factor do equal .6.
#' 
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` need to be defined. Note that neither may contain the noise factors.
#' If the model contains observed variables only, use `Lambda = diag(x)` where `x` is the number of variables.
#'
#' The order of the factors is (X1, X2, ..., X_nWaves).
#' 
#' Additional arguments related to the requested type of **power analysis**:
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.  This requires the `semTools` package.
#' 
#' @examples
#' \dontrun{
#' # Determine required N in a 10-wave ARMA model
#' # to detect that the autoregressive effects differ across waves
#' # with a power of 80% on alpha = 5%, where
#' # X is measured by 3 indicators loading by .5 each (at each wave), and
#' # the autoregressive effects vary between .5 and .7, and
#' # the moving average parameters are .3 at each wave and
#' # are assumed to be constant across waves (in both the H0 and the H1 model) and
#' # there are no lagged effects, and
#' # metric invariance and autocorrelated residuals are assumed
#' powerARMA <- semPower.powerARMA(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 10,
#'   autoregLag1 = c(.5, .7, .6, .5, .7, .6, .6, .5, .6),
#'   mvAvgLag1 = rep(.3, 9),
#'   variances = rep(1, 10),
#'   waveEqual = c('mvAvg'),
#'   nullEffect = 'autoreg',
#'   nIndicator = rep(3, 10), loadM = .5,
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE
#' )
#' 
#' # show summary
#' summary(powerARMA)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerARMA$modelH1, sample.cov = powerARMA$Sigma,
#'             sample.nobs = powerARMA$requiredN,
#'             sample.cov.rescale = FALSE)
#' lavaan::sem(powerARMA$modelH0, sample.cov = powerARMA$Sigma,
#'             sample.nobs = powerARMA$requiredN,
#'             sample.cov.rescale = FALSE)
#' 
#' 
#' # same as above, but determine power with N = 250 on alpha = .05
#' powerARMA <- semPower.powerARMA(
#'   'post-hoc', alpha = .05, N = 250,
#'   nWaves = 10,
#'   autoregLag1 = c(.5, .7, .6, .5, .7, .6, .6, .5, .6),
#'   mvAvgLag1 = rep(.3, 9),
#'   variances = rep(1, 10),
#'   waveEqual = c('mvAvg'),
#'   nullEffect = 'autoreg',
#'   nIndicator = rep(3, 10), loadM = .5,
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE
#' )
#' 
#' # same as above, but determine the critical chi-square with N = 250 so that alpha = beta
#' powerARMA <- semPower.powerARMA(
#'   'compromise', abratio = 1, N = 250,
#'   nWaves = 10,
#'   autoregLag1 = c(.5, .7, .6, .5, .7, .6, .6, .5, .6),
#'   mvAvgLag1 = rep(.3, 9),
#'   variances = rep(1, 10),
#'   waveEqual = c('mvAvg'),
#'   nullEffect = 'autoreg',
#'   nIndicator = rep(3, 10), loadM = .5,
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE
#' )
#'   
#' # same as above, but compare to the saturated model
#' # (rather than to the less restricted model)
#' powerARMA <- semPower.powerARMA(
#'   'a-priori', alpha = .05, power = .80, comparison = 'saturated',
#'   nWaves = 10,
#'   autoregLag1 = c(.5, .7, .6, .5, .7, .6, .6, .5, .6),
#'   mvAvgLag1 = rep(.3, 9),
#'   variances = rep(1, 10),
#'   waveEqual = c('mvAvg'),
#'   nullEffect = 'autoreg',
#'   nIndicator = rep(3, 10), loadM = .5,
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE
#' )
#' 
#' 
#' # same as above, but assume only observed variables
#' powerARMA <- semPower.powerARMA(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 10,
#'   autoregLag1 = c(.5, .7, .6, .5, .7, .6, .6, .5, .6),
#'   mvAvgLag1 = rep(.3, 9),
#'   variances = rep(1, 10),
#'   waveEqual = c('mvAvg'),
#'   nullEffect = 'autoreg',
#'   Lambda = diag(1, 10),
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE
#' )
#' 
#' # same as above, but provide reduced loadings matrix to define that
#' # X is measured by 3 indicators each loading by .5, .6, .4 (at each wave)
#' powerARMA <- semPower.powerARMA(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 10,
#'   autoregLag1 = c(.5, .7, .6, .5, .7, .6, .6, .5, .6),
#'   mvAvgLag1 = rep(.3, 9),
#'   variances = rep(1, 10),
#'   waveEqual = c('mvAvg'),
#'   nullEffect = 'autoreg',
#'   loadings = list(
#'     c(.5, .6, .4),  # X1
#'     c(.5, .6, .4),  # X2
#'     c(.5, .6, .4),  # X3
#'     c(.5, .6, .4),  # X4
#'     c(.5, .6, .4),  # X5
#'     c(.5, .6, .4),  # X6
#'     c(.5, .6, .4),  # X7
#'     c(.5, .6, .4),  # X8
#'     c(.5, .6, .4),  # X9
#'     c(.5, .6, .4)   # X10
#'   ),
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE
#' )
#' 
#' # same as above, but detect that the moving average parameters differ across waves
#' # with a power of 80% on alpha = 5%, where
#' # the moving average parameters vary between .05 and .4, and
#' # the autoregressive effects are .5 at each wave and
#' # are assumed to be constant across waves (in both the H0 and the H1 model)
#' powerARMA <- semPower.powerARMA(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 10,
#'   autoregLag1 = rep(.5, 9),
#'   mvAvgLag1 = c(.1, .05, .2, .1, .1, .3, .4, .4, .4),
#'   variances = rep(1, 10),
#'   waveEqual = c('autoreg'),
#'   nullEffect = 'mvAvg',
#'   nIndicator = rep(3, 10), loadM = .5,
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE
#' )
#' 
#' 
#' # same as above, but detect that the (noise) variances differ across waves
#' # with a power of 80% on alpha = 5%, where
#' # the variances vary between 0.5 and 2, and
#' # the autoregressive effects are .5 at each wave and
#' # the moving average parameters are .3 at each wave and
#' # bothj are assumed to be constant across waves (in both the H0 and the H1 model)
#' powerARMA <- semPower.powerARMA(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 10,
#'   autoregLag1 = rep(.5, 9),
#'   mvAvgLag1 = rep(.3, 9),
#'   variances = c(1, .5, .7, .6, .7, .9, 1.2, 1.7, 2.0, 1.5),
#'   waveEqual = c('autoreg', 'mvAvg'),
#'   nullEffect = 'var',
#'   nIndicator = rep(3, 10), loadM = .5,
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE
#' )
#' 
#' 
#' # same as above, but include a meanstructure and
#' # detect that the means differ across waves
#' # with a power of 80% on alpha = 5%, where
#' # the means vary between 0 and .5, and
#' # the autoregressive effects are .5 at each wave and
#' # the moving average parameters are .3 at each wave and
#' # the variances are 1 at each wave and
#' # all are assumed to be constant across waves (in both the H0 and the H1 model) and
#' # metric and scalar invariance is assumed
#' powerARMA <- semPower.powerARMA(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 10,
#'   autoregLag1 = rep(.5, 9),
#'   mvAvgLag1 = rep(.3, 9),
#'   variances = rep(1, 10),
#'   means = c(0, .1, .2, .3, .4, .5, .3, .4, .5, .5),
#'   waveEqual = c('autoreg', 'mvAvg', 'var'),
#'   nullEffect = 'mean',
#'   nIndicator = rep(3, 10), loadM = .5,
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE
#' )
#' 
#' # Determine required N in a 10-wave ARMA model
#' # to detect that the autoregressive lag-2 effects differ from zero
#' # with a power of 80% on alpha = 5%, where
#' # the lag-2 autoregressive effects are .2 at each wave and 
#' # the lag-2 autoregressive effects are .1 at each wave and
#' # the autoregressive effects are .5 at each wave and
#' # the moving average parameters are .3 at each wave and
#' # the noise variances are equal to 1 in each wave,
#' # and all are assumed to be constant across waves (in both the H0 and the H1 model) and
#' # metric invariance and autocorrelated residuals are assumed, and
#' # the autoregressive lag2- and lag3-effects are estimated 
#' powerARMA <- semPower.powerARMA(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 10,
#'   autoregLag1 = rep(.5, 9),
#'   autoregLag2 = rep(.2, 8),
#'   autoregLag3 = rep(.1, 7),
#'   mvAvgLag1 = rep(.3, 9),
#'   variances = rep(1, 10),
#'   waveEqual = c('mvAvg', 'autoreg', 'var', 'autoreglag2', 'autoreglag3'),
#'   nullEffect = 'autoreglag2 = 0',
#'   nIndicator = rep(3, 10), loadM = .5,
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE
#' )
#' 
#' # similar as above, but get required N to detect that 
#' # lag-2 moving average parameters are constant across waves 
#' powerARMA <- semPower.powerARMA(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 10,
#'   autoregLag1 = rep(.5, 9),
#'   autoregLag2 = rep(.2, 8),
#'   mvAvgLag1 = rep(.3, 9),
#'   mvAvgLag2 = c(.1, .2, .3, .1, .2, .3, .1, .1),
#'   variances = rep(1, 10),
#'   waveEqual = c('mvAvg', 'autoreg', 'var', 'autoreglag2'),
#'   nullEffect = 'mvAvgLag2',
#'   nIndicator = rep(3, 10), loadM = .5,
#'   invariance = TRUE
#' )
#' 
#' 
#' # Determine required N in a 5-wave ARMA model
#' # to detect that the autoregressive effects in group 1
#' # differ from the ones in group 2, where
#' # both groups are equal-sized
#' # with a power of 80% on alpha = 5%, where
#' # X is measured by 3 indicators loading by .5 each (at each wave and in each group), and
#' # the autoregressive effects in group 1 are .5 (constant across waves) and
#' # the autoregressive effects in group 2 are .6 (constant across waves) and
#' # the moving average parameters are .25 at each wave and in both groups and
#' # the variances are 1 at each wave and in both groups and 
#' # all are assumed to be constant across waves (in both the H0 and the H1 model)
#' # metric invariance (across both waves and groups) and
#' # autocorrelated residuals are assumed
#' powerARMA <- semPower.powerARMA(
#'   'a-priori', alpha = .05, power = .80, N = list(1, 1),
#'   nWaves = 5,
#'   autoregLag1 = list(
#'     c(.5, .5, .5, .5),   # group 1
#'     c(.6, .6, .6, .6)),  # group 2
#'   mvAvgLag1 = rep(.25, 4),
#'   variances = rep(1, 5),
#'   waveEqual = c('autoreg', 'var', 'mvavg'),
#'   nullEffect = 'autoregA = autoregB',
#'   nIndicator = rep(3, 5), loadM = .5,
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE
#' )
#' 
#' 
#' # Determine required N in a 5-wave ARMA model
#' # to detect that the means in group 1
#' # differ from the means in group 2, where
#' # both groups are equal-sized
#' # with a power of 80% on alpha = 5%, where
#' # X is measured by 3 indicators loading by .5 each (at each wave and in each group), and
#' # the autoregressive effects are .5 at each wave and in both groups and
#' # the moving average parameters are .25 at each wave and in both groups and
#' # the variances are 1 at each wave and in both groups and 
#' # all are assumed to be constant across waves (in both the H0 and the H1 model) and
#' # invariance of variances, autoregressive effects, and moving average parameters 
#' # across groups as well as
#' # metric and scalar invariance (across both waves and groups) and
#' # autocorrelated residuals are assumed
#' powerARMA <- semPower.powerARMA(
#'   'a-priori', alpha = .05, power = .80, N = list(1, 1),
#'   nWaves = 5,
#'   autoregLag1 = list(
#'     c(.5, .5, .5, .5),   # group 1
#'     c(.5, .5, .5, .5)),  # group 2
#'   mvAvgLag1 = rep(.25, 4),
#'   variances = rep(1, 5),
#'   means = list(
#'     c(0, .1, .1, .1, .1),  # group 1
#'     c(0, .4, .4, .4, .4)   # group 2
#'   ),
#'   waveEqual = c('autoreg', 'var', 'mvavg', 'mean'),
#'   groupEqual = c('var', 'autoreg', 'mvavg'),
#'   nullEffect = 'meanA = meanB',
#'   nIndicator = rep(3, 5), loadM = .5,
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE
#' )
#' 
#' # perform a simulated post-hoc power analysis
#' # with 250 replications
#' set.seed(300121)
#' powerARMA <- semPower.powerARMA(
#'   'post-hoc', alpha = .05, N = 500,
#'   nWaves = 5,
#'   autoregLag1 = c(.3, .7, .6, .3),
#'   mvAvgLag1 = rep(.3, 4),
#'   variances = rep(1, 5),
#'   waveEqual = c('mvAvg'),
#'   nullEffect = 'autoreg',
#'   nIndicator = rep(3, 5), loadM = .5,
#'   invariance = TRUE, 
#'   autocorResiduals = TRUE, 
#'   simulatedPower = TRUE,
#'   simOptions = list(nReplications = 250)
#' )
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerARMA <- function(type, comparison = 'restricted',
                               nWaves = NULL, 
                               autoregEffects = NULL, 
                               autoregLag1 = autoregEffects, 
                               autoregLag2 = NULL, 
                               autoregLag3 = NULL,
                               mvAvgLag1 = NULL,
                               mvAvgLag2 = NULL,
                               mvAvgLag3 = NULL,
                               means = NULL,
                               variances = NULL,
                               waveEqual = NULL, 
                               groupEqual = NULL,
                               nullEffect = NULL,
                               nullWhich = NULL,
                               nullWhichGroups = NULL,
                               invariance = TRUE,
                               autocorResiduals = TRUE,
                               ...){
  
  comparison <- checkComparisonModel(comparison)
  checkEllipsis(...)

  if('standardized' %in% names(list(...)) && list(...)[['standardized']]) stop('Standardized is not available for ARMA.')
  
  # we override Beta and Sigma later, so let's make sure it is not set in ellipsis argument
  if('Beta' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Beta.')
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma.')

  nullValid <- c('autoreg', 'autoreglag2', 'autoreglag3', 
                 'autoreg=0','autoreglag2=0', 'autoreglag3=0',
                 'mvavg', 'mvavglag2', 'mvavglag3',
                 'mvavg=0', 'mvavglag2=0', 'mvavglag3=0',
                 'var', 'mean',
                 'autorega=autoregb', 'mvavga=mvavgb', 
                 'vara=varb','meana=meanb')
  nullEffect <- checkNullEffect(nullEffect, nullValid)
  waveEqualValid <- c('autoreg', 'autoreglag2', 'autoreglag3',
                      'mvavg', 'mvavglag2', 'mvavglag3',
                      'var', 'mean')
  waveEqual <- unlist(lapply(waveEqual, function(x) tolower(trimws(x))))
  if(any(unlist(lapply(waveEqual, function(x) !x %in% waveEqualValid)))) stop(paste('waveEqual may only contain', paste(waveEqualValid, collapse = ' ')))
  # identification requires(?) at least one wave-equal constraint, so make this mandatory
  if(is.null(waveEqual)) stop('ARMA models without wave-equality constraints on any parameter are likely not identified, so specify at least one wave-constant parameter in waveEqual (e.g. autoreg, mvAvg) ')

  if(is.null(nWaves) || is.na(nWaves) || nWaves < 2) stop('nWaves must be >= 2.') 
  if(is.null(autoregEffects)) autoregEffects <- autoregLag1

  # determine whether we have a multigroup model
  nGroups <- 1
  cargs <- list(autoregEffects, mvAvgLag1, means, variances)
  cargs <- cargs[!sapply(cargs, is.null)]
  isMultigroup <- any(unlist(lapply(cargs, is.list)))
  if(isMultigroup){
    cargs <- cargs[sapply(cargs, is.list)]
    ig <- unlist(lapply(cargs, length))
    if(length(unique(ig[ig != 1])) > 1) stop('Non-null list arguments supplied to autoregEffects, mvAvgLag1, means, or variances imply a different number of groups. Make sure all lists have the same length or provide no list for no group differences.')
    nGroups <- max(ig)
    if(is.null(nullWhichGroups)) nullWhichGroups <- seq(nGroups)
  }

  if(isMultigroup && !nullEffect %in% c('autorega=autoregb','mvavga=mvavgb','vara=varb','meana=meanb')) stop('Multigroup analysis are only supported for nullEffect = autoregA=autoregB, mvavgA=mvavgaB, varA=varB, meanA=meanB')
  if(!isMultigroup && nullEffect %in% c('autorega=autoregb','mvavga=mvavgb','vara=varb','meana=meanb')) stop('nullEffect = autoregA=autoregB, mvavgA=mvavgaB, varA=varB, meanA=meanB imply multigroup analyses, but no list structure for any relevant parameter provided.')
  if(isMultigroup && !is.null(groupEqual)){
    groupEqualValid <- c('autoreg', 'mvavg', 'var', 'mean')  
    groupEqual <- unlist(lapply(groupEqual, function(x) tolower(trimws(x))))
    if(any(unlist(lapply(groupEqual, function(x) !x %in% groupEqualValid)))) stop(paste('groupEqual may only contain', paste(groupEqualValid, collapse = ' ')))
  }
  
  if(nullEffect %in% waveEqual) stop('You cannot set the same parameters in nullEffect and waveEqual')
  if(nullEffect %in% groupEqual) stop('You cannot set the same parameters in nullEffect and groupEqual')
  if(is.null(means) && (nullEffect %in% c('meana=meanb', 'mean') || 'mean' %in% waveEqual)) stop('Either nullEffect or waveEqual refer to means, but no means provided.')
  
  estimateAutoregLag2 <- !is.null(autoregLag2)
  estimateAutoregLag3 <- !is.null(autoregLag3)
  estimateMvAvgLag2 <- !is.null(mvAvgLag2)
  estimateMvAvgLag3 <- !is.null(mvAvgLag3)
  
  if(is.null(autoregLag2)) autoregLag2 <- rep(0, nWaves - 2)
  if(is.null(mvAvgLag2)) mvAvgLag2 <- rep(0, nWaves - 2)
  if(is.null(autoregLag3)) if(nWaves > 2) autoregLag3 <- rep(0, nWaves - 3) else autoregLag3 <- 0 # just to init properly, never actually used
  if(is.null(mvAvgLag3)) if(nWaves > 2) mvAvgLag3 <- rep(0, nWaves - 3) else mvAvgLag3 <- 0

  if(!is.list(autoregEffects)) autoregEffects <- rep(list(autoregEffects), nGroups)
  if(!is.list(autoregLag2)) autoregLag2 <- rep(list(autoregLag2), nGroups)
  if(!is.list(autoregLag3)) autoregLag3 <- rep(list(autoregLag3), nGroups)
  if(!is.list(mvAvgLag1)) mvAvgLag1 <- rep(list(mvAvgLag1), nGroups)
  if(!is.list(mvAvgLag2)) mvAvgLag2 <- rep(list(mvAvgLag2), nGroups)
  if(!is.list(mvAvgLag3)) mvAvgLag3 <- rep(list(mvAvgLag3), nGroups)
  if(!is.list(variances)) variances <- rep(list(variances), nGroups)
  if(!is.null(means) && !is.list(means)) means <- rep(list(means), nGroups)
  
  if(any(unlist(lapply(autoregEffects, function(x) length(x) != (nWaves - 1))))) stop('autoregEffects must be of length nWaves - 1.')
  if(any(unlist(lapply(autoregLag2, function(x) length(x) != (nWaves - 2))))) stop('autoregLag2 must be of length nWaves - 2.')
  if(any(unlist(lapply(autoregLag3, function(x) length(x) != (nWaves - 3))))) stop('autoregLag3 must be of length nWaves - 3.')
  if(any(unlist(lapply(mvAvgLag1, function(x) length(x) != (nWaves - 1))))) stop('mvAvgLag1 must be of length nWaves - 1.')
  if(any(unlist(lapply(mvAvgLag2, function(x) length(x) != (nWaves - 2))))) stop('mvAvgLag2 must be of length nWaves - 2.')
  if(any(unlist(lapply(mvAvgLag3, function(x) length(x) != (nWaves - 3))))) stop('mvAvgLag3 must be of length nWaves - 3.')
  if(!is.null(means) && any(unlist(lapply(means, function(x) length(x) != nWaves)))) stop('means must be of length nWaves.')
  if(any(unlist(lapply(variances, function(x) length(x) != nWaves)))) stop('variances must be of length nWaves.')

  if((nullEffect == 'var' || 'var' %in% waveEqual) && !invariance) stop('When nullEffect or waveEqual  refer to variances, invariance must be TRUE.')
  if((nullEffect == 'mean' || 'mean' %in% waveEqual) && !invariance) stop('When nullEffect or waveEqual refer to latent means, invariance must be TRUE.')

  if(is.null(nullWhich)){
    msg <- 'nullWhich must be defined when nullEffect refers to parameters that are not equal across waves.'
    if(nullEffect %in% c('autoreg=0', 'autorega=autoregb') && !'autoreg' %in% waveEqual) stop(msg)
    if(nullEffect %in% 'autoreglag2=0' && !'autoreglag2' %in% waveEqual && nWaves > 3) stop(msg)
    if(nullEffect %in% 'autoreglag3=0' && !'autoreglag3' %in% waveEqual && nWaves > 4) stop(msg)
    if(nullEffect %in% c('mvavg=0', 'mvavga=mvavgb') && !'mvavg' %in% waveEqual) stop(msg)
    if(nullEffect %in% 'mvavglag2=0' && !'mvavglag2' %in% waveEqual && nWaves > 3) stop(msg)
    if(nullEffect %in% 'mvavglag3=0' && !'mvavglag3' %in% waveEqual && nWaves > 4) stop(msg)
    if(nullEffect %in% 'vara=varb' && !'var' %in% waveEqual) stop(msg)
    if(nullEffect %in% 'meana=meanb' && !'var' %in% waveEqual) stop(msg)
    nullWhich <- 1 # proper default for all remaining cases
  }
  if(!is.null(nullWhich)){
    if(!is.numeric(nullWhich) || length(nullWhich) > 1) stop('nullWhich must be a single number.')
    if(nullWhich < 1 || nullWhich > (nWaves - 1)) stop('nullWhich must lie between 1 and nWaves - 1.')
    if(nullEffect %in% c('autoreglag2=0','autoreglag2', 'mvavglag2=0', 'mvavglag2') && nullWhich > (nWaves - 2)) stop('For lag2 effects, nullWhich must lie between 1 and nWaves - 2.')
    if(nullEffect %in% c('autoreglag3=0','autoreglag3', 'mvavglag3=0', 'mvavglag3') && nullWhich > (nWaves - 3)) stop('For lag3 effects, nullWhich must lie between 1 and nWaves - 3.')
  }

    
  ### get and modify lambda
  args <- list(...)
  Lambda  <- args[['Lambda']]
  if(is.null(Lambda)){
    Lambda <- genLambda(args[['loadings']], args[['nIndicator']],
                        args[['loadM']], args[['loadSD']], args[['loadMinMax']])
  }
  if(ncol(Lambda) != nWaves) stop('Number of factors must be equal to nWaves.')
  
  Lambda <- cbind(Lambda, matrix(0, nrow = nrow(Lambda), ncol = nWaves))  # add noise factors
  Lambda <- rep(list(Lambda), nGroups)  # require same measurement model across groups
  
  ### create Beta
  Beta <- lapply(seq(nGroups), function(x){
    B <- matrix(0, ncol = 2*nWaves, nrow = 2*nWaves)
    # add autoregEffects
    for(i in seq(nWaves - 1)){
      B[(i + 1), i] <- autoregEffects[[x]][i]
    }
    # lag-2 effects
    if(nWaves > 2){
      for(i in seq(nWaves - 2)){
        B[(i + 2), i] <- autoregLag2[[x]][i]
      }
    }
    # lag-3 effects
    if(nWaves > 3){
      for(i in seq(nWaves - 3)){
        B[(i + 3), i] <- autoregLag3[[x]][i]
      }
    }
    # lag-1 mov avgs
    for(i in seq(nWaves - 1)){
      idx <- nWaves + i
      B[(i + 1), idx] <- mvAvgLag1[[x]][i]
    }
    # lag-2 mov avgs
    if(nWaves > 2){
      for(i in seq(nWaves - 2)){
        idx <- nWaves + i
        B[(i + 2), idx] <- mvAvgLag2[[x]][i]
      }
    }
    # lag-3 mov avgs
    if(nWaves > 3){
      for(i in seq(nWaves - 3)){
        idx <- nWaves + i
        B[(i + 3), idx] <- mvAvgLag3[[x]][i]
      }
    }
    # noise factors
    for(i in seq(nWaves)){
      idx <- nWaves + i
      B[i, idx] <- 1
    }
    B
  })
  
  
  ### create Psi
  Psi <- lapply(seq(nGroups), function(x){
    P <- matrix(0, ncol = 2*nWaves, nrow = 2*nWaves)
    # var noise factors
    diag(P[seq((nWaves + 1), 2*nWaves), seq((nWaves + 1), 2*nWaves)]) <- variances[[x]]
    P
  })
  
  ### define latent means
  Alpha <- NULL
  if(!is.null(means)){
    Alpha <- lapply(means, function(x) c(x, rep(0, nWaves))) # F  + noise factors
  }
  
  ### add metric invariance constraints to analysis model
  metricInvarianceFactors <- NULL
  if(invariance) metricInvarianceFactors <- list(seq(nWaves))
  
  
  ### get sigma
  generated <- semPower.genSigma(Beta = if(!isMultigroup) Beta[[1]] else Beta,
                                 Psi = if(!isMultigroup) Psi[[1]] else Psi,
                                 Alpha = if(!isMultigroup) Alpha[[1]] else Alpha,
                                 Lambda = if(!isMultigroup) Lambda[[1]] else Lambda,
                                 useReferenceIndicator = TRUE,
                                 metricInvariance = metricInvarianceFactors, 
                                 nGroups = nGroups)
  
  
  ### create model string
  if(!isMultigroup) model <- generated[['modelTrueCFA']] else model <- generated[[1]][['modelTrueCFA']]
  Lambda <- Lambda[[1]]
  
  # resid var X
  tok <- paste0('f', 1:nWaves, ' ~~ 0*', 'f', 1:nWaves)
  model <- paste(model, paste(tok, collapse = '\n'), sep = '\n')
  
  # noise fac
  tok <- paste0('n', 1:nWaves, '=~ 1*', 'f', 1:nWaves)
  model <- paste(model, paste(tok, collapse = '\n'), sep = '\n')
  
  # var noise
  tok <- paste0('n', 1:nWaves, '~~ pvn', formatC(1:nWaves, width = 2, flag = 0), '*n', 1:nWaves)
  model <- paste(model, paste(tok, collapse = '\n'), sep = '\n')
  
  # covar noise
  tok <- list()
  for(i in 1:(nWaves - 1)){
    for(j in (i + 1):nWaves){
      if(i < j){
        tok <- append(tok,list(paste0('n',i, ' ~~ 0*n',j)))
      }
    }
  }
  model <- paste(model, paste(unlist(tok), collapse = '\n'), sep = '\n')
  
  # autoreg
  for(f in 2:nWaves){     # omit first row
    fidx <- (f - 1)
    if(estimateAutoregLag2 && f > 2) fidx <- c(fidx, (f - 2)) # estm lag2 effects regardless of these are zero
    if(estimateAutoregLag3 && f > 3) fidx <- c(fidx, (f - 3)) # estm lag3 effects regardless of these are zero
    tok <- paste0('f', f, ' ~ ', paste(paste0('pf', paste0(formatC(f, width = 2, flag = 0), formatC(fidx, width = 2, flag = 0)), '*'), paste0('f', fidx), sep = '', collapse = ' + '))
    model <- paste(model, tok, sep='\n')
  }
  
  # mv avg
  for(f in 2:nWaves){     # omit first row
    fidx <- (f - 1)
    if(estimateMvAvgLag2 && f > 2) fidx <- c(fidx, (f - 2)) # estm lag2 effects regardless of these are zero
    if(estimateMvAvgLag3 && f > 3) fidx <- c(fidx, (f - 3)) # estm lag3 effects regardless of these are zero
    tok <- paste0('f', f, ' ~ ', paste(paste0('pn', paste0(formatC(f, width = 2, flag = 0), formatC(fidx, width = 2, flag = 0)), '*'), paste0('n', fidx), sep = '', collapse = ' + '))
    model <- paste(model, tok, sep='\n')
  }
  
  # estm means + intercepts
  if(!is.null(means) && invariance){
    # fix first factor mean to zero
    model <- paste(model, 'f1 ~ 0*1', sep='\n')
    for(f in 2:nWaves){  # start at 2
      tok <- paste0('f', f, ' ~ ', 'pmf', formatC(f, width = 2, flag = 0), '*1')
      model <- paste(model, tok, sep='\n')
    }
    # fix means of noise factors to zero
    for(f in 1:nWaves){
      tok <- paste0('n', f, ' ~ 0*1')
      model <- paste(model, tok, sep='\n')
    }
    
    # add invariance constraints on intercepts
    tok <- list()
    for(x in seq_along(metricInvarianceFactors)){
      ci <- lapply(metricInvarianceFactors[[x]], function(f) paste0('x', which(Lambda[, f] != 0)))
      lab <- paste0('pmi', 1:length(ci[[1]]))
      tok <- append(tok, unlist(lapply(ci, function(x) paste0(x, ' ~ ', lab, '*1') )))
    }
    model <- paste(model, paste(unlist(tok), collapse = '\n'), sep='\n')
  }
  
  # add autocorrelated residuals
  if(autocorResiduals){
    # do this only when there is at least one latent variable
    if(nrow(Lambda) > nWaves){
      autocorResidualsFactors <- metricInvarianceFactors  # same structure 
      tok <- list()
      for(x in seq_along(autocorResidualsFactors)){
        ci <- lapply(autocorResidualsFactors[[x]], function(f) paste0('x', which(Lambda[, f] != 0)))
        if(length(ci[[1]]) > 1){
          for(i in 1:(length(ci) - 1)){
            for(j in (i + 1) : length(ci)){
              tok <- append(tok, paste(ci[[i]], '~~', ci[[j]]))
            }
          }
        }
      }
      model <- paste(model, paste(tok, collapse = '\n'), sep ='\n')
    }
  }
  
  
  ### define H1 and H0 model
  
  # first get relevant parameter labels that may be subject to waveEqual or nullEffect constraints 
  xw <- seq(nWaves, 2)
  pAutoregX <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(xw - 1, width = 2, flag = 0))
  pAutoregX <- pAutoregX[order(pAutoregX)]
  pMvAvg <- paste0('pn', formatC(xw, width = 2, flag = 0), formatC(xw - 1, width = 2, flag = 0))
  pMvAvg <- pMvAvg[order(pMvAvg)]
  pNoiseVar <- paste0('pvn', formatC(1:nWaves, width = 2, flag = 0))
  pMeans <- paste0('pmf', formatC(1:nWaves, width = 2, flag = 0))
  if(nWaves > 2){
    xw <- seq(nWaves, 3)
    pAutoregLag2 <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(xw - 2, width = 2, flag = 0))
    pAutoregLag2 <- pAutoregLag2[order(pAutoregLag2)]
    pMvAvgLag2 <- paste0('pn', formatC(xw, width = 2, flag = 0), formatC(xw - 2, width = 2, flag = 0))
    pMvAvgLag2 <- pMvAvgLag2[order(pMvAvgLag2)]
  }
  if(nWaves > 3){
    xw <- seq(nWaves, 4)
    pAutoregLag3 <- paste0('pf', formatC(xw, width = 2, flag = 0), formatC(xw - 3, width = 2, flag = 0))
    pAutoregLag3 <- pAutoregLag3[order(pAutoregLag3)]
    pMvAvgLag3 <- paste0('pn', formatC(xw, width = 2, flag = 0), formatC(xw - 3, width = 2, flag = 0))
    pMvAvgLag3 <- pMvAvgLag3[order(pMvAvgLag3)]
  }
  
  
  # multigroup case
  if(isMultigroup){
    # remove group specific labels from measurement part to enforce metric invariance
    model <- gsub(paste0('_g', seq(nGroups), collapse = '|'), '_gc', model)
    # add equal group labels to intercepts to enforce scalar invariance
    for(x in seq_along(metricInvarianceFactors)){
      ci <- lapply(metricInvarianceFactors[[x]], function(f) paste0('x', which(Lambda[, f] != 0)))
      lab <- paste0('pmi', 1:length(ci[[1]]))
      for(pp in seq(lab)){
        patt <- lab[pp]
        repl <- paste0('c(', paste(paste0(patt, rep('_gc', nGroups)), collapse = ', '), ')')
        model <- gsub(patt, repl, model)
      }
    }  
    # assign group labels to all structural parameters
    for(pp in seq(pAutoregX)){
      patt <- pAutoregX[pp]
      if('autoreg' %in% groupEqual){
        repl <- paste0('c(', paste(paste0(patt, rep('_gc', nGroups)), collapse = ', '), ')')
      }else{
        repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      }
      model <- gsub(patt, repl, model)
    }
    for(pp in seq(pMvAvg)){
      patt <- pMvAvg[pp]
      if('mvavg' %in% groupEqual){
        repl <- paste0('c(', paste(paste0(patt, rep('_gc', nGroups)), collapse = ', '), ')')
      }else{
        repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      }
      model <- gsub(patt, repl, model)
    }
    for(pp in seq(pMeans)){
      patt <- pMeans[pp]
      if('mean' %in% groupEqual){
        repl <- paste0('c(', paste(paste0(patt, rep('_gc', nGroups)), collapse = ', '), ')')
      }else{
        repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      }
      model <- gsub(patt, repl, model)
    }
    for(pp in seq(pNoiseVar)){
      patt <- pNoiseVar[pp]
      if('var' %in% groupEqual){
        repl <- paste0('c(', paste(paste0(patt, rep('_gc', nGroups)), collapse = ', '), ')')
      }else{
        repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
      }
      model <- gsub(patt, repl, model)
    }
    if(nWaves > 2){
      for(pp in seq(pAutoregLag2)){
        patt <- pAutoregLag2[pp]
        repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
        model <- gsub(patt, repl, model)
      }
      for(pp in seq(pMvAvgLag2)){
        patt <- pMvAvgLag2[pp]
        repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
        model <- gsub(patt, repl, model)
      }  
    }
    if(nWaves > 3){
      for(pp in seq(pAutoregLag3)){
        patt <- pAutoregLag3[pp]
        repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
        model <- gsub(patt, repl, model)
      }
      for(pp in seq(pMvAvgLag3)){
        patt <- pMvAvgLag3[pp]
        repl <- paste0('c(', paste(paste0(patt, '_g', seq(nGroups)), collapse = ', '), ')')
        model <- gsub(patt, repl, model)
      }
    }
  }

  # create H1 model
  modelH1 <- model
  if(c('autoreg') %in% waveEqual) modelH1 <- gsub(paste(pAutoregX, collapse = '|'), 'pf', modelH1)
  if(c('autoreglag2') %in% waveEqual) modelH1 <- gsub(paste(pAutoregLag2, collapse = '|'), 'pf2', modelH1)
  if(c('autoreglag3') %in% waveEqual) modelH1 <- gsub(paste(pAutoregLag3, collapse = '|'), 'pf3', modelH1)
  if(c('mvavg') %in% waveEqual) modelH1 <- gsub(paste(pMvAvg, collapse = '|'), 'pn', modelH1)
  if(c('mvavglag2') %in% waveEqual) modelH1 <- gsub(paste(pMvAvgLag2, collapse = '|'), 'pn2', modelH1)
  if(c('mvavglag3') %in% waveEqual) modelH1 <- gsub(paste(pMvAvgLag3, collapse = '|'), 'pn3', modelH1)
  # wave-equal contraints for means and vars do not include first measurement
  if(c('var') %in% waveEqual) modelH1 <- gsub(paste(pNoiseVar[-1], collapse = '|'), 'pvn', modelH1)
  if(c('mean') %in% waveEqual) modelH1 <- gsub(paste(pMeans[-1], collapse = '|'), 'pmf', modelH1)
  
  
  # add additional constraints to H0 model
  modelH0 <- modelH1
  
  # wave equal constraints not included in modelH1:
  if(c('autoreg') %in% nullEffect) modelH0 <- gsub(paste(pAutoregX, collapse = '|'), 'pf', modelH0)
  if(c('autoreglag2') %in% nullEffect) modelH0 <- gsub(paste(pAutoregLag2, collapse = '|'), 'pf2', modelH0)
  if(c('autoreglag3') %in% nullEffect) modelH0 <- gsub(paste(pAutoregLag3, collapse = '|'), 'pf3', modelH0)
  if(c('mvavg') %in% nullEffect) modelH0 <- gsub(paste(pMvAvg, collapse = '|'), 'pn', modelH0)
  if(c('mvavglag2') %in% nullEffect) modelH0 <- gsub(paste(pMvAvgLag2, collapse = '|'), 'pn2', modelH0)
  if(c('mvavglag3') %in% nullEffect) modelH0 <- gsub(paste(pMvAvgLag3, collapse = '|'), 'pn3', modelH0)
  # wave-equal contraints for means and vars do not include first measurement
  if(c('var') %in% nullEffect) modelH0 <- gsub(paste(pNoiseVar[-1], collapse = '|'), 'pvn', modelH0)
  if(c('mean') %in% nullEffect) modelH0 <- gsub(paste(pMeans[-1], collapse = '|'), 'pmf', modelH0)
  
  # zero constraints:
  if('autoreg=0' %in% nullEffect){
    if('autoreg' %in% waveEqual){
      modelH0 <- gsub('pf', '0', modelH0)
    }else{
      modelH0 <- gsub(pAutoregX[nullWhich], '0', modelH0)
    }
  }
  if('autoreglag2=0' %in% nullEffect){
    if('autoreglag2' %in% waveEqual){
      modelH0 <- gsub('pf2', '0', modelH0)
    }else{
      modelH0 <- gsub(pAutoregLag2[nullWhich], '0', modelH0)
    }
  }
  if('autoreglag3=0' %in% nullEffect){
    if('autoreglag3' %in% waveEqual){
      modelH0 <- gsub('pf3', '0', modelH0)
    }else{
      modelH0 <- gsub(pAutoregLag3[nullWhich], '0', modelH0)
    }
  }
  if('mvavg=0' %in% nullEffect){
    if('mvavg' %in% waveEqual){
      modelH0 <- gsub('pn', '0', modelH0)
    }else{
      modelH0 <- gsub(pMvAvg[nullWhich], '0', modelH0)
    }
  }
  if('mvavglag2=0' %in% nullEffect){
    if('mvavglag2' %in% waveEqual){
      modelH0 <- gsub('pn2', '0', modelH0)
    }else{
      modelH0 <- gsub(pMvAvgLag2[nullWhich], '0', modelH0)
    }
  }
  if('mvavglag3=0' %in% nullEffect){
    if('mvavglag3' %in% waveEqual){
      modelH0 <- gsub('pn3', '0', modelH0)
    }else{
      modelH0 <- gsub(pMvAvgLag3[nullWhich], '0', modelH0)
    }
  }
  
  # multigroup constraints
  if('autorega=autoregb' %in% nullEffect){
    if(c('autoreg') %in% waveEqual){
      patt <- paste0('pf_g', nullWhichGroups, collapse = '|')
      repl <-'pf_gc'
    }else{
      patt <- paste0(pAutoregX[nullWhich], '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(pAutoregX[nullWhich], '_gc')
    }
    modelH0 <- gsub(patt, repl, modelH0)
  }
  if('mvavga=mvavgb' %in% nullEffect){
    if(c('mvavg') %in% waveEqual){
      patt <- paste0('pn_g', nullWhichGroups, collapse = '|')
      repl <-'pn_gc'
    }else{
      patt <- paste0(pMvAvg[nullWhich], '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(pMvAvg[nullWhich], '_gc')
    }
    modelH0 <- gsub(patt, repl, modelH0)
  }
  if('vara=varb' %in% nullEffect){
    if(c('var') %in% waveEqual){
      patt <- paste0('pvn_g', nullWhichGroups, collapse = '|')
      repl <-'pvn_gc'
    }else{
      patt <- paste0(pNoiseVar[nullWhich], '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(pNoiseVar[nullWhich], '_gc')
    }
    modelH0 <- gsub(patt, repl, modelH0)
  }
  if('meana=meanb' %in% nullEffect){
    if(c('mean') %in% waveEqual){
      patt <- paste0('pmf_g', nullWhichGroups, collapse = '|')
      repl <-'pmf_gc'
    }else{
      patt <- paste0(pMeans[nullWhich], '_g', nullWhichGroups, collapse = '|')
      repl <- paste0(pMeans[nullWhich], '_gc')
    }
    modelH0 <- gsub(patt, repl, modelH0)
  }

  
  # here we actually fit modelH1 in case of a restricted comparison
  # because we cannot be sure that user input yields perfectly fitting h1 models, 
  # when there are additional constraints (waveequal or invariance)
  if(comparison == 'saturated') modelH1 <- NULL
  
  if(isMultigroup) Sigma <- lapply(generated, '[[', 'Sigma') else Sigma <- generated[['Sigma']] 
  mu <- NULL
  if(!is.null(means)){
    if(isMultigroup) mu <- lapply(generated, '[[', 'mu') else mu <- generated[['mu']] 
  }
  
  semPower.powerLav(type, 
                    Sigma = Sigma,
                    mu = mu,
                    modelH0 = modelH0, 
                    modelH1 = modelH1, 
                    ...)
}

#' semPower.powerLGCM
#'
#' Convenience function for performing power analysis on effects in a latent growth curve model (LGCM).
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param nWaves number of waves, must be >= 3 for linear and > 3 for quadratic trends.
#' @param means vector providing the means of the intercept and the linear slope factor (and the quadratic slope factor, if `quadratic = TRUE`). A list for multiple group models.
#' @param variances vector providing the variances of the intercept and the linear slope factor (and the quadratic slope factor, if `quadratic = TRUE`). Can be omitted, if a matrix is provided to `covariances`. Takes precedence over the diagonal in `covariances` when both are defined. A list for multiple group models.
#' @param covariances either the variance-covariance matrix between the intercept and the slope (and the quadratic slope factor, if `quadratic = TRUE`), or a single number giving the covariance between intercept and slope factor, or `NULL` for orthogonal factors. A list for multiple group models.
#' @param quadratic whether to include a quadratic slope factor in addition to a linear slope factor. Defaults to `FALSE` for no quadratic slope factor.
#' @param timeCodes vector of length `nWaves` defining the loadings on the slope factor. If omitted, the time codes default to (0, 1, ..., (nWaves - 1)).
#' @param ticExogSlopes vector defining the slopes for an exogenous time-invariant covariate in the prediction of the intercept and slope factors (and the quadratic slope factor, if `quadratic = TRUE`). Can be omitted for no covariate.
#' @param ticEndogSlopes vector defining the slopes for the intercept and slope factors (and the quadratic slope factor, if `quadratic = TRUE`) in the prediction of an endogenous time-invariant covariate. Can be omitted for no covariate. 
#' @param groupEqual parameters that are restricted across groups in both the H0 and the H1 model, when `nullEffect` implies a multiple group model. Valid are `'imean'`, `'smean'`, `'s2mean'` to restrict the means of the intercept, linear slope, and quadratic slope factors, and `'ivar'`, `'svar'`, `'s2var'` for their variances, and `'iscov'`, `'is2cov'`, `'ss2cov'` for the covariances between intercept and slope, intercept and quadratic slope, and linear and quadratic slope.
#' @param nullEffect defines the hypothesis of interest. See details for valid arguments. 
#' @param nullWhichGroups for hypothesis involving cross-groups comparisons, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param autocorResiduals whether the residuals of the indicators of latent variables are autocorrelated over waves (`TRUE`, the default) or not (`FALSE`). This is only applied to second order LGCMs. This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model).
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. The order of factors is (X1, X2, ..., X_nWaves, ticExogSlopes, ticEndogSlopes). See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details
#' This function performs a power analysis to reject various hypotheses arising
#' in latent growth curve models (LGCM), where one variable `X` is repeatedly 
#' assessed at different time points (`nWaves`), and a latent intercept and a 
#' linear (and optionally a quadratic) latent slope factor is assumed.
#'  
#' Relevant hypotheses in arising in a LCGM are:
#' * `iMean = 0`, `sMean = 0`, `s2Mean = 0`: Tests the hypothesis that the mean of the intercept, linear slope, and quadratic slope factors, respectively, is zero.
#' * `iVar = 0`, `sVar = 0`, `s2Var = 0`: Tests the hypothesis that the variance of the intercept, linear slope, and quadratic slope factors, respectively, is zero.
#' * `isCov = 0`: Tests the hypothesis that covariance between the intercept and linear slope factor is zero.
#' * `is2Cov = 0`: Tests the hypothesis that covariance between the intercept and quadratic slope factor is zero.
#' * `ss2Cov = 0`: Tests the hypothesis that covariance between the linear and the quadratic slope factor is zero.
#' * `betaIT = 0`, `betaST = 0`, `betaS2T = 0`: Tests the hypothesis that the slope for an exogenous time-invariant covariate in the prediction of the intercept, the linear slope, and the quadratic slope factor, respectively, is zero (`TIC -> I, S, S2`). 
#' * `betaTI = 0`, `betaTS = 0`, `betaTS2 = 0`: Tests the hypothesis that the slope the intercept, the linear slope, and the quadratic slope factor, respectively, in the prediction of an endogenous time-invariant covariate is zero (`I, S, S2 -> TIC`).
#' * `iMeanA = iMeanB`, `sMeanA = sMeanB`, `s2MeanA = s2MeanB`: Tests the hypothesis that the means of the intercept, linear slope, and quadratic slope factors, respectively, are equal across groups.
#' * `iVarA = iVarB`, `sVarA = sVarB`, `s2VarA = s2VarB`: Tests the hypothesis that the variances of the intercept, linear slope, and quadratic slope factors, respectively, are equal across groups.
#' * `isCovA = isCovA`: Tests the hypothesis that covariance between the intercept and linear slope factor is equal across groups.
#' * `is2CovA = is2CovA`: Tests the hypothesis that the covariance between the intercept and quadratic slope factor is equal across groups.
#' * `ss2CovA = ss2CovA`: Tests the hypothesis that the covariance between the linear and quadratic slope factor is equal across groups.
#' * `betaITA = betaITB`, `betaSTA = betaSTB`, `betaS2TA = betaS2TB`: Tests the hypothesis that the slopes for the time-invariant covariate in the prediction of the intercept, the linear slope, and the quadratic slope factor, respectively, are equal across groups (`TIC -> I, S, S2`). 
#' * `betaTIA = betaTIB`, `betaTSA = betaTSB`, `betaTS2A = betaTS2B`: Tests the hypothesis that the slope the intercept, the linear slope, and the quadratic slope factor, respectively, in the prediction of the time-invariant covariate are equal across groups (`I, S, S2 -> TIC`).
#' 
#' For hypotheses regarding longitudinal invariance, see [semPower.powerLI()]. For hypotheses regarding a simple autoregression, see [semPower.powerAutoregressive()]. For hypotheses in an ARMA model, see [semPower.powerARMA()].
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model**:
#' * `Lambda`: The factor loading matrix (with the number of columns equaling the number of factors).
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading on the second factor by .8, .6, .6, and .4.
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the loadings of the first factor to equal .5 and those of the second factor do equal .6.
#' 
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` need to be defined. Neither may contain entries referring to the intercept and slope factors.
#' If the model contains observed variables only, use `Lambda = diag(x)` where `x` is the number of variables.
#'
#' The order of the factors is (X1, X2, ..., X_nWaves, ticExogenous, ticEndogenous). If ticExogenous is undefined, ticEndogenous takes its place.
#' 
#' Additional arguments related to the requested type of **power analysis**:
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.  This requires the `semTools` package.
#' 
#' @examples
#' \dontrun{
#' # Determine required N in a 3-wave LGCM model
#' # to detect that the mean of the slope factor differs from zero
#' # with a power of 80% on alpha = 5%, where
#' # X is measured by 3 indicators loading by .5 each (at each wave), and
#' # the mean of the intercept factor is .5 and
#' # the mean of the slope factor is .2 and
#' # the variance of the intercept factor is 1 and
#' # the variance of the slope factor is .5 and
#' # the intercept-slope covariance is .25 and
#' # autocorrelated residuals are assumed
#' powerLGCM <- semPower.powerLGCM(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 3,
#'   means = c(.5, .2),     # i, s
#'   variances = c(1, .5),  # i, s
#'   covariances = .25,
#'   nullEffect = 'sMean = 0',
#'   nIndicator = rep(3, 3), loadM = .5,
#'   autocorResiduals = TRUE
#' )
#' 
#' # show summary
#' summary(powerLGCM)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerLGCM$modelH1, sample.cov = powerLGCM$Sigma,
#'             sample.mean = powerLGCM$mu,
#'             sample.nobs = powerLGCM$requiredN,
#'             sample.cov.rescale = FALSE)
#' lavaan::sem(powerLGCM$modelH0, sample.cov = powerLGCM$Sigma,
#'             sample.mean = powerLGCM$mu,
#'             sample.nobs = powerLGCM$requiredN,
#'             sample.cov.rescale = FALSE)
#' 
#' 
#' # same as above, but determine power with N = 250 on alpha = .05
#' powerLGCM <- semPower.powerLGCM(
#'   'post-hoc', alpha = .05, N = 250,
#'   nWaves = 3,
#'   means = c(.5, .2),
#'   variances = c(1, .5),
#'   covariances = .25,
#'   nullEffect = 'sMean = 0',
#'   nIndicator = rep(3, 3), loadM = .5,
#'   autocorResiduals = TRUE
#' )
#' 
#' # same as above, but determine the critical chi-square with N = 250 so that alpha = beta
#' powerLGCM <- semPower.powerLGCM(
#'   'compromise', abratio = 1, N = 250,
#'   nWaves = 3,
#'   means = c(.5, .2),
#'   variances = c(1, .5),
#'   covariances = .25,
#'   nullEffect = 'sMean = 0',
#'   nIndicator = rep(3, 3), loadM = .5,
#'   autocorResiduals = TRUE
#' )
#' 
#' # same as above, but compare to the saturated model
#' # (rather than to the less restricted model)
#' powerLGCM <- semPower.powerLGCM(
#'   'a-priori', alpha = .05, power = .80, comparison = 'saturated',
#'   nWaves = 3,
#'   means = c(.5, .2),
#'   variances = c(1, .5),
#'   covariances = .25,
#'   nullEffect = 'sMean = 0',
#'   nIndicator = rep(3, 3), loadM = .5,
#'   autocorResiduals = TRUE
#' )
#' 
#' # same as above, but assume only observed variables
#' powerLGCM <- semPower.powerLGCM(
#'   'a-priori', alpha = .05, power = .80, comparison = 'saturated',
#'   nWaves = 3,
#'   means = c(.5, .2),
#'   variances = c(1, .5),
#'   covariances = .25,
#'   nullEffect = 'sMean = 0',
#'   Lambda = diag(3),
#'   autocorResiduals = TRUE
#' )
#' 
#' # same as above, but provide reduced loadings matrix to define that
#' # X is measured by 3 indicators each loading by .5, .6, .4 (at each wave)
#' powerLGCM <- semPower.powerLGCM(
#'   'a-priori', alpha = .05, power = .80, comparison = 'saturated',
#'   nWaves = 3,
#'   means = c(.5, .2),
#'   variances = c(1, .5),
#'   covariances = .25,
#'   nullEffect = 'sMean = 0',
#'   loadings = list(
#'     c(.5, .6, .4),  # X1
#'     c(.5, .6, .4),  # X2 
#'     c(.5, .6, .4)   # X3
#'   ),
#'   autocorResiduals = TRUE
#' )
#' 
#' # same as above, but get required N to detect that
#' # the variance of the intercept factor differs from zero
#' powerLGCM <- semPower.powerLGCM(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 3,
#'   means = c(.5, .2),
#'   variances = c(1, .5),
#'   covariances = .25,
#'   nullEffect = 'iVar = 0',
#'   nIndicator = rep(3, 3), loadM = .5,
#'   autocorResiduals = TRUE
#' )
#' 
#' # same as above, but get required N to detect that
#' # the intercept-slope covariance differs from zero
#' powerLGCM <- semPower.powerLGCM(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 3,
#'   means = c(.5, .2),
#'   variances = c(1, .5),
#'   covariances = .25,
#'   nullEffect = 'iscov = 0',
#'   nIndicator = rep(3, 3), loadM = .5,
#'   autocorResiduals = TRUE
#' )
#' 
#' 
#' # include a quadratic slope factor
#' # and get required N to detect that
#' # its variance differs from zero.
#' # provide the variance-covariance matrix
#' # between intercept, slope, and quadratic slope factors
#' powerLGCM <- semPower.powerLGCM(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 4,
#'   quadratic = TRUE,
#'   means = c(.5, .2, .1),
#'   covariances = matrix(c(
#'     # i  s   s2
#'     c(1, .3, .2),
#'     c(.3, .5, .01),
#'     c(.2, .01, .1)
#'   ), ncol = 3, byrow = TRUE),
#'   nullEffect = 's2var = 0',
#'   nIndicator = rep(3, 4), loadM = .5,
#'   autocorResiduals = TRUE
#' )
#' 
#' # Determine required N in a 3-wave LGCM model
#' # to detect that the slope of an time-invariant covariate (TIC)
#' # on the slope factor differs from zero. 
#' # The TIC is measured by 4 indicators loading
#' # by .7, .7, .5, and .8. The slope of the TIC in the prediction of
#' # the intercept factor is .5, and in the prediction of the slope factor is .4.
#' powerLGCM <- semPower.powerLGCM(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 3,
#'   means = c(.5, .2),
#'   variances = c(1, .5),
#'   covariances = .25,
#'   ticExogSlopes = c(.5, .4),  # i s
#'   nullEffect = 'betaST = 0',
#'   loadings = list(
#'     c(.5, .6, .4),      # X1
#'     c(.5, .6, .4),      # X2 
#'     c(.5, .6, .4),      # X3
#'     c(.7, .7, .5, .8)   # TIC
#'   ),
#'   autocorResiduals = TRUE
#' )
#' 
#' # Determine required N in a 3-wave LGCM model
#' # to detect that the slope of the slope factor in 
#' # the prediction of a time-invariant covariate (TIC) differs from zero. 
#' # The TIC is measured by 4 indicators loading
#' # by .7, .7, .5, and .8. The slopes of the intercept and the slope factors in
#' # the prediction of the TIC are .1 and .3, respectively.
#' powerLGCM <- semPower.powerLGCM(
#'   'a-priori', alpha = .05, power = .80,
#'   nWaves = 3,
#'   means = c(.5, .2),
#'   variances = c(1, .5),
#'   covariances = .25,
#'   ticEndogSlopes = c(.1, .3),  # i s
#'   nullEffect = 'betaTS = 0',
#'   loadings = list(
#'     c(.5, .6, .4),      # X1
#'     c(.5, .6, .4),      # X2 
#'     c(.5, .6, .4),      # X3
#'     c(.7, .7, .5, .8)   # TIC
#'   ),
#'   autocorResiduals = TRUE
#' )
#' 
#' # Determine required N in a 3-wave LGCM model
#' # to detect that the mean of the slope factor in group 1
#' # differs from the mean of the slope factor in group 2
#' # with a power of 80% on alpha = 5%, where
#' # X is measured by 3 indicators loading by .5 each (at each wave and in each group), and
#' # the means of the intercept factor in group 1 and 2 are .5 and .25
#' # the means of the slope factor in group 1 and 2 are .25 and .4
#' # the variance of the intercept factor is 1 in both groups and
#' # the variance of the slope factor is .5in both groups and
#' # the intercept-slope covariance is .25 in both groups and
#' # autocorrelated residuals are assumed and
#' # the groups are equal-sized
#' powerLGCM <- semPower.powerLGCM(
#'   'a-priori', alpha = .05, power = .80, N = list(1, 1),
#'   nWaves = 3,
#'   means = list(
#'     # i, s
#'     c(.5, .2),     # group 1 
#'     c(.25, .4)),   # group 2
#'   variances = c(1, .5),
#'   covariances = .25,
#'   nullEffect = 'sMeanA = sMeanB',
#'   nIndicator = rep(3, 3), loadM = .5,
#'   autocorResiduals = TRUE
#' )
#' 
#' # similar as above, but get required N to detect that
#' # the intercept-slope covariance differs across groups, 
#' # assuming that intercept and slope variances are equal across groups.
#' powerLGCM <- semPower.powerLGCM(
#'   'a-priori', alpha = .05, power = .80, N = list(1, 1),
#'   nWaves = 3,
#'   means = c(.5, .2),  
#'   variances = c(1, .5),
#'   covariances = list(
#'     c(.25), # group 1 
#'     c(.1)), # group 2
#'   nullEffect = 'isCovA = isCovB',
#'   groupEqual = c('ivar', 'svar'),
#'   nIndicator = rep(3, 3), loadM = .5,
#'   autocorResiduals = TRUE
#' )
#' 
#' # perform a simulated post-hoc power analysis
#' # with 250 replications
#' set.seed(300121)
#' powerLGCM <- semPower.powerLGCM(
#'   'post-hoc', alpha = .05, N = 500,
#'   nWaves = 3,
#'   means = c(.5, .2),     # i, s
#'   variances = c(1, .5),  # i, s
#'   covariances = .25,
#'   nullEffect = 'sMean = 0',
#'   nIndicator = rep(3, 3), loadM = .5,
#'   autocorResiduals = TRUE,
#'   simulatedPower = TRUE,
#'   simOptions = list(nReplications = 250)
#' )
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerLGCM <- function(type, comparison = 'restricted',
                               nWaves = NULL, 
                               means = NULL, 
                               variances = NULL, 
                               covariances = NULL,
                               quadratic = FALSE,
                               timeCodes = NULL,
                               ticExogSlopes = NULL,
                               ticEndogSlopes = NULL,
                               groupEqual = NULL,
                               nullEffect = NULL,
                               nullWhichGroups = NULL,
                               autocorResiduals = TRUE,
                               ...){
  
  comparison <- checkComparisonModel(comparison)
  checkEllipsis(...)
  
  if('standardized' %in% names(list(...)) && list(...)[['standardized']]) stop('Standardized is not available for LGCM')
  
  # we override Beta and Sigma later, so let's make sure it is not set in ellipsis argument
  if('Beta' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Beta.')
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma.')
  
  nullValid <- c('imean=0', 'smean=0', 's2mean=0', 
                 'ivar=0', 'svar=0', 's2var=0', 
                 'iscov=0', 'is2cov=0', 'ss2cov=0',
                 'betait=0', 'betast=0', 'betas2t=0',  # tic -> i, s, s2
                 'betati=0', 'betats=0', 'betats2=0',  # i, s, s2 -> tic
                 'imeana=imeanb', 'smeana=smeanb', 's2meana=s2meanb',
                 'ivara=ivarb', 'svara=svarb', 's2vara=s2varb',
                 'iscova=iscovb', 'is2cova=is2covb', 'ss2cova=ss2covb',
                 'betatia=betatib', 'betatsa=betatsb', 'betats2a=betats2b',
                 'betaita=betaitb', 'betasta=betastb','betas2ta=betas2tb')
  nullEffect <- checkNullEffect(nullEffect, nullValid)
  
  if(is.null(nWaves) || is.na(nWaves) || nWaves < 3) stop('nWaves must be >= 3.') 
  if(quadratic && (is.null(nWaves) || is.na(nWaves) || nWaves < 4)) stop('For quadratic slopes, nWaves must be >= 4.') 
  if(is.null(means)) stop('Means must be provided.') 
  if(is.null(variances) && is.null(covariances)) stop('Variances (or a covariance matrix) must be provided.') 
  if(is.null(timeCodes)) timeCodes <- seq(nWaves) - 1
  
  # determine whether we have a multigroup model
  nGroups <- 1
  cargs <- list(means, variances, covariances, ticExogSlopes, ticEndogSlopes)
  cargs <- cargs[!sapply(cargs, is.null)]
  isMultigroup <- any(unlist(lapply(cargs, is.list)))
  if(isMultigroup){
    cargs <- cargs[sapply(cargs, is.list)]
    ig <- unlist(lapply(cargs, length))
    if(length(unique(ig[ig != 1])) > 1) stop('Non-null list arguments supplied to means, variances, covariances, ticPredictor, ticCriterion imply a different number of groups. Make sure all lists have the same length or provide no list for no group differences.')
    nGroups <- max(ig)
    if(is.null(nullWhichGroups)) nullWhichGroups <- seq(nGroups)
  }
  
  if(isMultigroup && 
     nullEffect %in% c('imean=0', 'smean=0', 's2mean=0', 'iscov=0', 'ivar=0', 'svar=0', 's2var=0', 'is2cov=0', 'ss2cov=0',
                       'betati=0', 'betats=0', 'betats2=0', 'betait=0', 'betast=0', 'betas2t=0')) 
    stop('Multiple groups are not supported for specified nullEffect. Make sure not to define lists for relevant parameters.')
  
  if(!isMultigroup && 
     nullEffect %in% c('imeana=imeanb', 'smeana=smeanb', 's2meana=s2meanb',
                       'ivara=ivarb', 'svara=svarb', 's2vara=s2varb',
                       'iscova=iscovb', 'is2cova=is2covb', 'ss2cova=ss2covb',
                       'betatia=betatib', 'betatsa=betatsb', 'betats2a=betats2b',
                       'betaita=betaitb', 'betasta=betastb','betas2ta=betas2tb')) 
    stop('Specified nullEffect implies multiple groups, but no list structure for a relevant parameter provided.')
  
  if(is.null(ticExogSlopes) && nullEffect %in% c('betait=0', 'betast=0', 'betas2t=0', 'betaita=betaitb', 'betasta=betastb','betas2ta=betas2tb'))
    stop('Specified nullEffect refers to an exogenous time-invariant covariate as predictor, but slopes in ticExogSlopes are not provided.')
  if(is.null(ticEndogSlopes) && nullEffect %in% c('betati=0', 'betats=0', 'betats2=0', 'betatia=betatb', 'betatsa=betatsb', 'betats2a=betats2b'))
    stop('Specified nullEffect refers to an endogenous time-invariant covariate as criterion, but slopes in ticEndogSlopes are not provided.')
  if(grepl('s2', nullEffect) && !quadratic) stop('nulleffect refers to quadratic slope, but no quadratic slope requested.')
  
  if(!is.list(means)) means <- rep(list(means), nGroups)
  if(!is.null(variances) && !is.list(variances)) variances <- rep(list(variances), nGroups)
  if(!is.null(covariances) && !is.list(covariances)) covariances <- rep(list(covariances), nGroups)
  if(!is.null(ticExogSlopes) && !is.list(ticExogSlopes)) ticExogSlopes <- rep(list(ticExogSlopes), nGroups)
  if(!is.null(ticEndogSlopes) && !is.list(ticEndogSlopes)) ticEndogSlopes <- rep(list(ticEndogSlopes), nGroups)
  
  vLength <- 2 + sum(quadratic)
  if(is.null(covariances)) covariances <- rep(list(diag(vLength)), nGroups)
  if(is.null(variances) && !all(unlist(lapply(covariances, is.matrix)))) stop('If variances are not provided, covariances must be a matrix.') 
  if(!any(unlist(lapply(covariances, is.matrix)))){
    if(quadratic) stop('If quadratic slopes are included, covariances must be a matrix.')
    covariances <- lapply(covariances, function(x){
      if(length(x) != 1) stop('If no matrix is provided, covariance argument must contain a single number.')
      m <- diag(2)
      m[1, 2] <- m[2, 1] <- x
      m
    })
  }else{
    lapply(covariances, function(x){
      checkSymmetricSquare(x, 'Matrix provided to covariances argument')
      if(ncol(x) != vLength) stop(paste('Matrix provided to covariances argument must be of dimension ', vLength, ' * ', vLength))
    })
  } 
  if(is.null(variances)) variances <- lapply(covariances, diag)
  if(!is.null(variances)) covariances <- lapply(seq(covariances), function(x){
    diag(covariances[[x]]) <- variances[[x]]
    covariances[[x]]
  })
  
  if(any(unlist(lapply(means, function(x) length(x) != vLength)))) stop(paste('Exactly', vLength, 'means must be provided.'))
  if(any(unlist(lapply(variances, function(x) length(x) != vLength)))) stop(paste('Exactly', vLength, 'variances must be provided.'))
  if(!is.null(ticExogSlopes) && any(unlist(lapply(ticExogSlopes, function(x) length(x) != vLength)))) stop(paste('Exactly', vLength, 'slopes for ticExogSlopes must be provided.'))
  if(!is.null(ticEndogSlopes) && any(unlist(lapply(ticEndogSlopes, function(x) length(x) != vLength)))) stop(paste('Exactly', vLength, 'slopes for ticEndogSlopes must be provided.'))
  
  if(isMultigroup && !is.null(groupEqual)){
    groupEqualValid <- c('imean', 'smean', 's2mean', 'ivar', 'svar', 's2var',
                         'iscov', 'is2cov', 'ss2cov')
    groupEqual <- unlist(lapply(groupEqual, function(x) tolower(trimws(x))))
    if(any(unlist(lapply(groupEqual, function(x) !x %in% groupEqualValid)))) stop(paste('groupEqual may only contain', paste(groupEqualValid, collapse = ' ')))
    if(nullEffect %in% groupEqual) stop('You cannot set the same parameters in nullEffect and groupEqual.')
  }
  
  
  ### fetch and update Lambda
  args <- list(...)
  Lambda  <- args[['Lambda']]
  if(is.null(Lambda)){
    Lambda <- genLambda(args[['loadings']], args[['nIndicator']],
                        args[['loadM']], args[['loadSD']], args[['loadMinMax']])
  }
  if(ncol(Lambda) != (nWaves + sum(!is.null(ticExogSlopes)) + sum(!is.null(ticEndogSlopes)))) stop('A single factor is needed for each wave, optionally plus additional single factors for time-invariant covariates.')
  
  # order of factors is: f1 - f_nwaves, (int, slope, slope2) , ticpred, ticcrit
  m1 <- cbind(Lambda[, seq(nWaves)],
              rep(0, nrow(Lambda)),
              rep(0, nrow(Lambda)))
  if(quadratic) m1 <- cbind(m1, rep(0, nrow(Lambda)))
  if(ncol(Lambda) > nWaves){
    Lambda <- cbind(m1, Lambda[, (nWaves + 1) : ncol(Lambda)])
  }else{
    Lambda <- m1
  }
  
  ### create Beta
  Beta <- lapply(seq(nGroups), function(x){
    dim <- nWaves + 2 + sum(quadratic) + sum(!is.null(ticExogSlopes)) + sum(!is.null(ticEndogSlopes))
    B <- matrix(0, ncol = dim, nrow = dim)
    
    for(i in seq(nWaves)){
      # intercept factor
      B[i, (nWaves + 1)] <- 1
      # slope factor
      B[i, (nWaves + 2)] <- timeCodes[i]
      # slope2 factor
      if(quadratic) B[i, (nWaves + 3)] <- timeCodes[i]^2
    }
    # ticExogSlopes
    if(!is.null(ticExogSlopes)){
      idx <- seq((1 + nWaves), nWaves + 2 + sum(quadratic))
      idx2 <- 1 + nWaves + 2 + sum(quadratic)
      B[idx, idx2] <- ticExogSlopes[[x]] 
    } 
    # ticEndogSlopes
    if(!is.null(ticEndogSlopes)){
      idx <- seq((1 + nWaves), nWaves + 2 + sum(quadratic))
      idx2 <- 1 + nWaves + 2 + sum(quadratic) + sum(!is.null(ticExogSlopes))
      B[idx2, idx] <- ticEndogSlopes[[x]] 
    }
    
    B
  })
  
  ### create Psi
  Psi <- lapply(seq(nGroups), function(x){
    # vars for tics = 1
    P <- diag(ncol(Beta[[x]]))
    # resid var
    cB <- Beta[[x]][seq(nWaves), seq(nWaves + 1, nWaves + 2 + sum(quadratic))]
    lPsi <- 1 / (diag(cB %*% covariances[[x]] %*% t(cB)))
    diag(P[seq(nWaves), seq(nWaves)]) <- lPsi
    # i
    P[nWaves + 1, nWaves + 1] <- variances[[x]][1]
    # s
    P[nWaves + 2, nWaves + 2] <- variances[[x]][2]
    # i-s 
    P[nWaves + 1, nWaves + 2] <- P[nWaves + 2, nWaves + 1] <- covariances[[x]][1, 2]
    if(quadratic){
      # s2
      P[(nWaves + 3), (nWaves + 3)] <- variances[[x]][3]    
      # i-s2
      P[(nWaves + 3), (nWaves + 1)] <- P[(nWaves + 1), (nWaves + 3)] <- covariances[[x]][1, 3]
      # s-s2
      P[(nWaves + 3), (nWaves + 2)] <- P[(nWaves + 2), (nWaves + 3)] <- covariances[[x]][2, 3]
    } 
    P
  })
  
  ### create Alpha
  Alpha <- lapply(seq(nGroups), function(x){
    A <- rep(0, ncol(Beta[[1]]))
    idx <- seq(nWaves + 1, nWaves + 2 + sum(quadratic))
    A[idx] <- means[[x]]
    A
  })
  
  
  ### gen Sigma
  generated <- semPower.genSigma(Beta = if(!isMultigroup) Beta[[1]] else Beta,
                                 Psi = if(!isMultigroup) Psi[[1]] else Psi,
                                 Alpha = if(!isMultigroup) Alpha[[1]] else Alpha,
                                 Lambda = if(!isMultigroup) Lambda else rep(list(Lambda), nGroups),
                                 useReferenceIndicator = TRUE,
                                 metricInvariance = list(seq(nWaves)), 
                                 nGroups = nGroups)
  
  ### define model strings
  if(!isMultigroup) model <- generated[['modelTrueCFA']] else model <- generated[[1]][['modelTrueCFA']]
  
  # i, s, s2
  model <- append(model, paste0('i =~ ', paste0('1*f', seq(nWaves), collapse = ' + ')))
  model <- append(model, paste0('s =~ ', paste0(timeCodes, '*f', seq(nWaves), collapse = ' + ')))
  if(quadratic) model <- append(model, paste0('s2 =~ ', paste0(timeCodes^2, '*f', seq(nWaves), collapse = ' + ')))
  
  # intercepts
  for(i in seq(nWaves)){
    ci <- which(Lambda[ , i] != 0)  
    model <- append(model, paste0('x', ci[1], ' ~ 0*1')) # referent indicator
    if(length(ci) > 1){
      if(isMultigroup){
        # do this here: use same group labels to enforce scalar invariance across groups
        lab <- list()
        for(x in 2:length(ci)){
          lab <- append(lab, paste0('c(', paste(rep(paste0('i', x, '_gc'), nGroups), collapse = ', '), ')'))
        }
        lab <- unlist(lab)
      }else{
        lab <- paste0('i', 2:length(ci))
      }
      model <- append(model, paste0('x', ci[-1], ' ~ ', lab, '*1'))
    }
  }
  
  # var+cov
  model <- append(model, c('i ~~ vi*i', 's ~~ vs*s', 'i ~~ cis*s'))
  if(quadratic) model <- append(model, c('s2 ~~ vs2*s2', 'i ~~ cis2*s2', 's ~~ css2*s2'))
  # means
  model <- append(model, paste0('f', seq(nWaves), ' ~ 0*1'))
  model <- append(model, c('i ~ mi*1', 's ~ ms*1'))
  if(quadratic) model <- append(model, c('s2 ~ ms2*1'))
  
  # ticExogSlopes
  if(!is.null(ticExogSlopes)){
    idx <- 1 + nWaves + 2 + sum(quadratic)
    model <- append(model, c(paste0('i ~ it*f', idx), paste0('s ~ st*f', idx)))
    if(quadratic) model <- append(model, paste0('s2 ~ s2t*f', idx))
  }
  # tic crit
  if(!is.null(ticEndogSlopes)){
    idx <- 1 + nWaves + 2 + sum(quadratic) + sum(!is.null(ticExogSlopes))
    model <- append(model, c(paste0('f', idx,' ~ ti*i'), paste0('f', idx,' ~ ts*s'))) 
    if(quadratic) model <- append(model, paste0('f', idx,' ~ ts2*s2')) 
  }
  # cov tic pred/crit 
  if(!is.null(ticExogSlopes) && !is.null(ticEndogSlopes)){
    idx <- 1 + nWaves + 2 + sum(!is.null(quadratic))
    model <- append(model, paste0('f', idx, ' ~~ f', (idx + 1))) 
  }
  
  # add autocorrelated residuals
  if(autocorResiduals){
    # do this only for multiple indicator latents 
    iLambda <- Lambda[ , seq(nWaves)]
    if(any(apply(iLambda, 2, function(x) sum(x != 0)) > 1)){
      autocorResidualsFactors <- list(seq(nWaves)) 
      tok <- list()
      for(x in seq_along(autocorResidualsFactors)){
        ci <- lapply(autocorResidualsFactors[[x]], function(f) paste0('x', which(Lambda[, f] != 0)))
        if(length(ci[[1]]) > 1){
          for(i in 1:(length(ci) - 1)){
            for(j in (i + 1) : length(ci)){
              tok <- append(tok, paste(ci[[i]], '~~', ci[[j]]))
            }
          }
        }
      }
      model <- append(model, tok)
    }
  }
  
  if(isMultigroup){
    # remove group specific labels from measurement part to enforce metric invariance
    model <- gsub(paste0('_g', seq(nGroups), collapse = '|'), '_gc', model)
    # equal group labels to intercepts were already added above
    
    ### assign group labels to all structural parameters
    # imean
    lab <- paste0('mi_g', seq(nGroups))
    if('imean' %in% groupEqual) lab <- rep('mi_gc', nGroups)
    model <- gsub('mi', paste0('c(', paste(lab, collapse = ', '), ')'), model)
    # smean
    lab <- paste0('ms_g', seq(nGroups))
    if('smean' %in% groupEqual) lab <- rep('ms_gc', nGroups)
                  # so that this does not match ms2
    model <- gsub('ms\\*', paste0('c(', paste(lab, collapse = ', '), ')*'), model)
    # s2mean
    lab <- paste0('ms2_g', seq(nGroups))
    if('s2mean' %in% groupEqual) lab <- rep('ms2_gc', nGroups)
    model <- gsub('ms2', paste0('c(', paste(lab, collapse = ', '), ')'), model)
    # ivar
    lab <- paste0('vi_g', seq(nGroups))
    if('ivar' %in% groupEqual) lab <- rep('vi_gc', nGroups)
    model <- gsub('vi', paste0('c(', paste(lab, collapse = ', '), ')'), model)
    # svar
    lab <- paste0('vs_g', seq(nGroups))
    if('svar' %in% groupEqual) lab <- rep('vs_gc', nGroups)
    model <- gsub('vs\\*', paste0('c(', paste(lab, collapse = ', '), ')*'), model)
    # s2var
    lab <- paste0('vs2_g', seq(nGroups))
    if('s2var' %in% groupEqual) lab <- rep('vs2_gc', nGroups)
    model <- gsub('vs2', paste0('c(', paste(lab, collapse = ', '), ')'), model)
    # iscov
    lab <- paste0('cis_g', seq(nGroups))
    if('iscov' %in% groupEqual) lab <- rep('cis_gc', nGroups)
    model <- gsub('cis\\*', paste0('c(', paste(lab, collapse = ', '), ')*'), model)
    # is2cov
    lab <- paste0('cis2_g', seq(nGroups))
    if('is2cov' %in% groupEqual) lab <- rep('cis2_gc', nGroups)
    model <- gsub('cis2', paste0('c(', paste(lab, collapse = ', '), ')'), model)
    # ss2cov
    lab <- paste0('css2_g', seq(nGroups))
    if('ss2cov' %in% groupEqual) lab <- rep('css2_gc', nGroups)
    model <- gsub('css2', paste0('c(', paste(lab, collapse = ', '), ')'), model)
    
    # tics (no groupequal option, cannot think of a use case)
    model <- gsub('it', paste0('c(', paste(paste0('it_g', seq(nGroups)), collapse = ', '), ')'), model)
    model <- gsub('st', paste0('c(', paste(paste0('st_g', seq(nGroups)), collapse = ', '), ')'), model)
    model <- gsub('s2t', paste0('c(', paste(paste0('s2t_g', seq(nGroups)), collapse = ', '), ')'), model)
    model <- gsub('ti', paste0('c(', paste(paste0('ti_g', seq(nGroups)), collapse = ', '), ')'), model)
    model <- gsub('ts\\*', paste0('c(', paste(paste0('ts_g', seq(nGroups)), collapse = ', '), ')*'), model)
    model <- gsub('ts2', paste0('c(', paste(paste0('ts2_g', seq(nGroups)), collapse = ', '), ')'), model)
  }

  modelH1 <- paste(unlist(model), collapse = '\n')
  
  # add additional constraints to H0 model
  modelH0 <- modelH1
  
  ### zero constraints
  # mean, cov
  if('imean=0' %in% nullEffect) modelH0 <- gsub('mi', '0', modelH0)
  if('smean=0' %in% nullEffect) modelH0 <- gsub('ms\\*', '0*', modelH0)
  if('s2mean=0' %in% nullEffect) modelH0 <- gsub('ms2', '0', modelH0)
  if('ivar=0' %in% nullEffect) modelH0 <- gsub('vi', '0', modelH0)
  if('svar=0' %in% nullEffect) modelH0 <- gsub('vs\\*', '0*', modelH0)
  if('s2var=0' %in% nullEffect) modelH0 <- gsub('vs2', '0', modelH0)
  if('iscov=0' %in% nullEffect) modelH0 <- gsub('cis\\*', '0*', modelH0)
  if('is2cov=0' %in% nullEffect) modelH0 <- gsub('cis2', '0', modelH0)
  if('ss2cov=0' %in% nullEffect) modelH0 <- gsub('css2', '0', modelH0)
  # tic
  if('betait=0' %in% nullEffect) modelH0 <- gsub('it', '0', modelH0)
  if('betast=0' %in% nullEffect) modelH0 <- gsub('st', '0', modelH0)
  if('betas2t=0' %in% nullEffect) modelH0 <- gsub('s2t', '0', modelH0)
  if('betati=0' %in% nullEffect) modelH0 <- gsub('ti', '0', modelH0)
  if('betats=0' %in% nullEffect) modelH0 <- gsub('ts\\*', '0*', modelH0)
  if('betats2=0' %in% nullEffect) modelH0 <- gsub('ts2', '0', modelH0)
  
  ### multigroup constraints
  if(isMultigroup){
    # mean
    if('imeana=imeanb' %in% nullEffect) modelH0 <- gsub(paste0('mi_g', nullWhichGroups, collapse = '|'), 'mi_gc', modelH0)
    if('smeana=smeanb' %in% nullEffect) modelH0 <- gsub(paste0('ms_g', nullWhichGroups, collapse = '|'), 'ms_gc', modelH0)
    if('s2meana=s2meanb' %in% nullEffect) modelH0 <- gsub(paste0('ms2_g', nullWhichGroups, collapse = '|'), 'ms2_gc', modelH0)
    # var
    if('ivara=ivarb' %in% nullEffect) modelH0 <- gsub(paste0('vi_g', nullWhichGroups, collapse = '|'), 'vi_gc', modelH0)
    if('svara=svarb' %in% nullEffect) modelH0 <- gsub(paste0('vs_g', nullWhichGroups, collapse = '|'), 'vs_gc', modelH0)
    if('s2vara=s2varb' %in% nullEffect) modelH0 <- gsub(paste0('vs2_g', nullWhichGroups, collapse = '|'), 'vs2_gc', modelH0)
    # cov
    if('iscova=iscovb' %in% nullEffect) modelH0 <- gsub(paste0('cis_g', nullWhichGroups, collapse = '|'), 'cis_gc', modelH0)
    if('is2cova=is2covb' %in% nullEffect) modelH0 <- gsub(paste0('cis2_g', nullWhichGroups, collapse = '|'), 'cis2_gc', modelH0)
    if('ss2cova=ss2covb' %in% nullEffect) modelH0 <- gsub(paste0('css2_g', nullWhichGroups, collapse = '|'), 'css2_gc', modelH0)
    # tic
    if('betaita=betaitb' %in% nullEffect) modelH0 <- gsub(paste0('it_g', nullWhichGroups, collapse = '|'), 'it_gc', modelH0)
    if('betasta=betastb' %in% nullEffect) modelH0 <- gsub(paste0('st_g', nullWhichGroups, collapse = '|'), 'st_gc', modelH0)
    if('betas2ta=betas2tb' %in% nullEffect) modelH0 <- gsub(paste0('s2t_g', nullWhichGroups, collapse = '|'), 's2t_gc', modelH0)
    if('betatia=betatib' %in% nullEffect) modelH0 <- gsub(paste0('ti_g', nullWhichGroups, collapse = '|'), 'ti_gc', modelH0)
    if('betatsa=betatsb' %in% nullEffect) modelH0 <- gsub(paste0('ts_g', nullWhichGroups, collapse = '|'), 'ts_gc', modelH0)
    if('betats2a=betats2b' %in% nullEffect) modelH0 <- gsub(paste0('ts2_g', nullWhichGroups, collapse = '|'), 'ts2_gc', modelH0)
  }
  
  # here we actually fit modelH1 in case of a restricted comparison
  # because we cannot be sure that user input yields perfectly fitting h1 models, 
  # when there are additional constraints (waveequal or invariance)
  if(comparison == 'saturated') modelH1 <- NULL
  
  if(isMultigroup) Sigma <- lapply(generated, '[[', 'Sigma') else Sigma <- generated[['Sigma']] 
  if(isMultigroup) mu <- lapply(generated, '[[', 'mu') else mu <- generated[['mu']] 
  
  semPower.powerLav(type, 
                    Sigma = Sigma,
                    mu = mu,
                    modelH0 = modelH0, 
                    modelH1 = modelH1, 
                    ...)
  
}
