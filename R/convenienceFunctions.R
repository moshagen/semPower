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
#' @return A list containing the following components is returned:
#' \item{`power`}{the results of the power analysis. Use the `summary` method to obtain formatted results.}
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
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
#' Optional arguments if a **simulated power analysis** (`simulatedPower = TRUE`) is requested:
#' * `nReplications`: The number of simulation runs. Defaults to 250, but larger numbers greatly improve accuracy at the expense of increased computation time.
#' * `minConvergenceRate`: The required minimum convergence rate. Defaults to .50.
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
#' summary(powerLav$power)
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
#' gen <- semPower.genSigma(Phi = .2, loadings = list(c(.5, .6, .4), c(.7, .8, .3)))
#' Sigma <- gen$Sigma
#' powerLav <- semPower.powerLav(type = 'a-priori', 
#'                               Sigma = Sigma, modelH0 = mH0,
#'                               alpha = .05, beta = .05)
#'
#' # note all of the above is identical to the output provided by the powerCFA convenience function
#' powerCFA <- semPower.powerCFA(type = 'a-priori',
#'                               comparison = 'saturated',
#'                               Phi = .2, loadings = list(c(.5, .6, .4), c(.7, .8, .3)), 
#'                               alpha = .05, beta = .05)
#' 
#' # same as above, but use simulated power analysis based on a robust ML estimator
#' powerLav <- semPower.powerLav(type = 'a-priori', 
#'                               Sigma = Sigma, modelH0 = mH0,
#'                               alpha = .05, beta = .05, 
#'                               simulatedPower = TRUE, nReplications = 500,
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
#' Convenience function for performing power analyses for CFA models to reject one of the following hypotheses: 
#' (a) a zero correlation between two factors, (b) the equality of two correlations between factors,
#' or (c) the equality of a correlation between two factors across two or more groups. 
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param Phi either a single number defining the correlation between exactly two factors or the factor correlation matrix. A list for multiple group models.
#' @param nullEffect defines the hypothesis of interest, must be one of `'cor = 0'` (the default) to test whether a correlation is zero, `'corX = corZ'` to test for the equality of correlations, and `'corA = corB'` to test for the equality of a correlation across groups. Define the correlations to be set to equality in `nullWhich` and the groups in `nullWhichGroups`. 
#' @param nullWhich vector of size 2 indicating which factor correlation in `Phi` is hypothesized to equal zero when `nullEffect = 'cor = 0'`, or to restrict to equality across groups when `nullEffect = 'corA = corB'`, or list of vectors defining which correlations to restrict to equality when `nullEffect = 'corX = corZ'`. Can also contain more than two correlations, e.g., `list(c(1, 2), c(1, 3), c(2, 3))` to set `Phi[1, 2] = Phi[1, 3] = Phi[2, 3]`. If omitted, the correlation between the first and the second factor is targeted, i. e., `nullWhich = c(1, 2)`.
#' @param nullWhichGroups for `nullEffect = 'corA = corB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. See details.
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
#' This function performs a power analysis to reject various hypotheses arising
#' in standard CFA models:
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
#' Optional arguments if a **simulated power analysis** (`simulatedPower = TRUE`) is requested:
#' * `nReplications`: The number of simulation runs. Defaults to 250, but larger numbers greatly improve accuracy at the expense of increased computation time.
#' * `minConvergenceRate`: The required minimum convergence rate. Defaults to .50.
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
#' summary(powercfa$power)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powercfa$modelH1, sample.cov = powercfa$Sigma, 
#' sample.nobs = powercfa$power$requiredN, sample.cov.rescale = FALSE)
#' lavaan::sem(powercfa$modelH0, sample.cov = powercfa$Sigma, 
#' sample.nobs = powercfa$power$requiredN, sample.cov.rescale = FALSE)
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
#'                               loadings = list(c(.7, .6, .5), c(.5, .6, .4, .8)),
#'                               alpha = .05, beta = .05)
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
#' powercfa <- semPower.powerCFA(type = 'post-hoc',
#'                               Phi = .2, 
#'                               nIndicator = c(5, 6), loadM = .5,
#'                               alpha = .05, N = 500, 
#'                               simulatedPower = TRUE, nReplications = 500)
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
  checkEllipsis(...)
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
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param slopes vector of standardized slopes (or a single number for a single slope) of the k predictors for Y. A list of slopes for multigroup models.
#' @param corXX correlation(s) between the k predictors (X). Either `NULL` for uncorrelated predictors, a single number (for k = 2 predictors), or a matrix. Can also be a list for multigroup models providing the correlations by group of matrices (otherwise, the same correlations are used in all groups). 
#' @param nullEffect defines the hypothesis of interest, must be one of `'slope = 0'` (the default) to test whether a slope is zero, `'slopeX = slopeZ'` to test for the equality of slopes, or `'slopeA = slopeB'` to test for the equality of slopes across groups. Define the slopes to set to equality in `nullWhich`.
#' @param nullWhich single number indicating which slope is hypothesized to equal zero when `nullEffect = 'slope = 0'`, or indicating which slope to restrict to equality across groups when `nullEffect = 'slopeA = slopeB'`, or vector defining the slopes to restrict to equality when `nullEffect = 'slopeX = slopeZ'`. Can also contain more than two slopes, e.g. `c(1, 2, 3)` to constrain the first three slopes to equality.
#' @param nullWhichGroups for `nullEffect = 'slopeA = slopeB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. The first factor is treated as Y and the subsequent factors as the predictors X_k. See details.
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
#' This function performs a power analysis to reject various hypotheses arising
#' in SEM models involving a simple regression relation of the form `Y = b_1*X_1 + ... + b_k*X_k` between the factors:
#' * `nullEffect = 'slope = 0'`: Tests the hypothesis that the slope for a predictor is zero. 
#' * `nullEffect = 'slopeX = slopeZ'`: Tests the hypothesis that two or more slopes are equal to each other.
#' * `nullEffect = 'slopeA = slopeB'`: Tests the hypothesis that the slope for a predictor is equal in two or more groups (always assuming metric invariance).
#' 
#' For hypotheses regarding mediation effects, see [semPower.powerMediation()].
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
#' Optional arguments if a **simulated power analysis** (`simulatedPower = TRUE`) is requested:
#' * `nReplications`: The number of simulation runs. Defaults to 250, but larger numbers greatly improve accuracy at the expense of increased computation time.
#' * `minConvergenceRate`: The required minimum convergence rate. Defaults to .50.
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
#'                                      slopes = c(.2, .3), corXX = .4, nullWhich = 1, 
#'                                      nIndicator = c(3, 5, 4), loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#' # show summary
#' summary(powerReg$power)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerReg$modelH1, sample.cov = powerReg$Sigma, 
#' sample.nobs = powerReg$power$requiredN, sample.cov.rescale = FALSE)
#' lavaan::sem(powerReg$modelH0, sample.cov = powerReg$Sigma, 
#' sample.nobs = powerReg$power$requiredN, sample.cov.rescale = FALSE)
#' 
#' # same as above, but determine power with N = 500 on alpha = .05 
#' powerReg <- semPower.powerRegression(type = 'post-hoc',
#'                                      slopes = c(.2, .3), corXX = .4, nullWhich = 1, 
#'                                      nIndicator = c(3, 5, 4), loadM = c(.5, .6, .7),
#'                                      alpha = .05, N = 500)
#'                                      
#' # same as above, but determine the critical chi-square with N = 500 so that alpha = beta 
#' powerReg <- semPower.powerRegression(type = 'compromise',
#'                                      slopes = c(.2, .3), corXX = .4, nullWhich = 1, 
#'                                      nIndicator = c(3, 5, 4), loadM = c(.5, .6, .7),
#'                                      abratio = .05, N = 500)
#'                                      
#' # same as above, but ask for the required N to detect that the slope of X2 is zero
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), corXX = .4, nullWhich = 2, 
#'                                      nIndicator = c(3, 5, 4), loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'                                      
#' # same as above, but compare to the saturated model
#' # (rather than to the less restricted model)
#' powerReg <- semPower.powerRegression(type = 'a-priori', comparison = 'saturated',
#'                                      slopes = c(.2, .3), corXX = .4, nullWhich = 2, 
#'                                      nIndicator = c(3, 5, 4), loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'                                      
#' # same as above, but provide a reduced loading matrix defining
#' # three indicators with loadings of .7, .6, .5 on the first factor (Y),
#' # four indicators with loadings of .5, .6, .4, .8 on the second factor (X1), and
#' # three indicators with loadings of .8, .7, .8 on the third factor (X2).
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), corXX = .4, nullWhich = 2, 
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
#'                                      slopes = list(c(.2, .3, .4), c(.2, .0, .4)), 
#'                                      corXX = corXX, 
#'                                      nullEffect = 'slopeA = slopeB', 
#'                                      nullWhich = 2,
#'                                      nIndicator = c(4, 5, 3, 5),
#'                                      loadM = c(.5, .6, .7, .6), 
#'                                      alpha = .05, beta = .05, N = list(1, 1))
#'
#'# request a simulated post-hoc power analysis with 500 replications 
#'# to detect that the slope of X1 differs from zero.
#' powerReg <- semPower.powerRegression(type = 'post-hoc',
#'                                      slopes = c(.2, .1), 
#'                                      nullWhich = 1,
#'                                      nIndicator = c(4, 3, 3), loadM = .5,
#'                                      alpha = .05, N = 500, 
#'                                      simulatedPower = TRUE, nReplications = 500)
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
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param bYX the standardized slope (direct effect) for X -> Y. A list for multiple group models. Can be `NULL` if `Beta` is set.
#' @param bMX the standardized slope for X -> M. A list for multiple group models. Can be `NULL` if `Beta` is set.
#' @param bYM the standardized slope for M -> Y. A list for multiple group models. Can be `NULL` if `Beta` is set.
#' @param Beta can be used instead of `bYX`, `bMX`, and `bYM`: matrix of regression weights connecting the latent factors (all-Y notation). Exogenous variables must be in the first row(s), so the upper triangular of Beta must be zero. A list for multiple group models.
#' @param indirect `NULL` unless `Beta` is set. Otherwise a list of vectors of size 2 indicating the elements of `Beta` that define the indirect effect of interest, e.g. `list(c(2, 1), c(3, 2))`. See details.
#' @param nullEffect defines the hypothesis of interest, must be one of `'ind = 0'` (the default) to test whether the indirect effect is zero or `'indA = indB'` to test for the equality of indirect effects across groups. See details.
#' @param nullWhichGroups for `nullEffect = 'indA = indB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. In case of a simple mediation, the order of factors is X, M, Y. See details.
#' @return A list containing the following components is returned:
#' \item{`power`}{the results of the power analysis. Use the `summary` method to obtain formatted results.}
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' @details
#' This function performs a power analysis to reject various hypotheses arising
#' in the context of mediation:
#' * `nullEffect = 'ind = 0'`: Tests the hypothesis that an indirect effect is zero. 
#' * `nullEffect = 'indA = indB'`: Tests the hypothesis that an indirect effect is equal in two or more groups (always assuming metric invariance). This is currently only possible for models without latent variables.
#' 
#' The indirect effect of interest can be specified in two ways:
#' * If a simple mediation involving three variables of the form `X -> M -> Y` is assumed, the arguments
#' `bYX`, `bMX`, and `bYM` are used to define the respective slopes, e. g.  `bYX = .4`, `bMX = .5`, and `bYM = .3` translates to
#' `X -- .5 --> M -- .3 --> Y` and  `X -- .4 --> Y`.
#' * More complex mediation structures can be defined by providing the `Beta` matrix along with `indirect` specifying which paths define the indirect effect. See examples below.
#' 
#' Notes on implementation:
#' * For models without latent variables, `nullEffect = 'ind = 0'` and `nullEffect = 'indA = indB'` constrain the indirect effect to zero and to equality, respectively.
#' * For models with latent variables and `nullEffect = 'ind = 0'`, power is approximated by constraining the smallest slope contained in the indirect effect to zero. Although this hypothesis is not identical to the hypothesis of the absence of an indirect effect, it provides a close approximation of the resulting power.
#' * For models with latent variables multiple groups (i. e., `nullEffect = 'indA = indB'`), there is currently no way to determine power, 
#' because implementing equality constrains on the indirect effects leads to non-convergence and 
#' the approach to constrain a single or all slopes to equality across groups can severely misrepresent the actual hypothesis of interest. 
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
#' Optional arguments if a **simulated power analysis** (`simulatedPower = TRUE`) is requested:
#' * `nReplications`: The number of simulation runs. Defaults to 250, but larger numbers greatly improve accuracy at the expense of increased computation time.
#' * `minConvergenceRate`: The required minimum convergence rate. Defaults to .50.
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
#' summary(powerMed$power)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerMed$modelH1, sample.cov = powerMed$Sigma,
#' sample.nobs = powerMed$power$requiredN, sample.cov.rescale = FALSE)
#' lavaan::sem(powerMed$modelH0, sample.cov = powerMed$Sigma,
#' sample.nobs = powerMed$power$requiredN, sample.cov.rescale = FALSE)
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
#'                                     indirect = list(c(2, 1), c(3, 2), c(4, 3)),
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
#' powerMed <- semPower.powerMediation(type = 'post-hoc',
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(3, 5, 4),
#'                                     loadM = c(.5, .6, .7),
#'                                     alpha = .05, N = 500,
#'                                     simulatedPower = TRUE, nReplications = 500)
#'}
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerMediation <- function(type, comparison = 'restricted',
                                    bYX = NULL, bMX = NULL, bYM = NULL,
                                    Beta = NULL, indirect = NULL, 
                                    nullEffect = 'ind = 0',
                                    nullWhichGroups = NULL, 
                                    ...){
  
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
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param nWaves number of waves, must be >= 2.
#' @param autoregEffects vector of the autoregressive effects of X and Y (constant across waves), or a list of vectors of autoregressive effects for X and Y from wave to wave, e.g. `list(c(.7, .6), c(.5, .5))` for a autoregressive effect of .7 for x1->x2 and .6 for x2->x3
#' @param crossedEffects vector of crossed effects of X on Y (X -> Y) and vice versa (both constant across waves), or a list of vectors of crossed effects giving the crossed effect of X on Y (and vice versa) for each wave, e. g. `list(c(.2, .3), c(.1, .1))` for x1->y2 = .2 and x2->y3 = .3.
#' @param rXY vector of (residual-)correlations between X and Y for each wave. If NULL, all (residual-)correlations are zero. 
#' @param waveEqual parameters that are assumed to be equal across waves in both the H0 and the H1 model. Valid are `'autoregX'` and `'autoregY'` for autoregressive effects, `'crossedX'` and `'crossedY'` for crossed effects, `'corXY'` for residual correlations, or `NULL` for none (so that all parameters are freely estimated, subject to the constraints defined in `nullEffect`). 
#' @param nullEffect defines the hypothesis of interest. Valid are the same arguments as in `waveEqual` and additionally `'autoregX = 0'`, `'autoregY = 0'`, `'crossedX = 0'`, `'crossedY = 0'` to constrain the X or Y autoregressive effects or the crossed effects to zero, `'autoregX = autoregY'` and `'crossedX = crossedY'` to constrain them to be equal for X and Y.
#' @param nullWhich used in conjunction with `nullEffect` to identify which parameter to constrain when there are > 2 waves and parameters are not constant across waves. For example, `nullEffect = 'autoregX = 0'` with `nullWhich = 2` would constrain the second autoregressive effect for X to zero.    
#' @param metricInvariance whether metric invariance over waves is assumed (`TRUE`, the default) or not (`FALSE`). This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model, where the df are not affected  by invariance constraints).
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. The order of factors is (X1, Y1, X2, Y2, ..., X_nWaves, Y_nWaves). See details.
#' @return A list containing the following components is returned:
#' \item{`power`}{the results of the power analysis. Use the `summary` method to obtain formatted results.}
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
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
#' * `crossedX = 0` and `crossedY = 0`: Tests the hypothesis that the crossed effect of X on Y (`crossedX`) and Y on X (`crossedY`), respectively, is zero. 
#' * `autoregX = autoregY`: Tests the hypothesis that the autoregressive effect of X and Y are equal.
#' * `crossedX = crossedY`: Tests the hypothesis that the crossed effect of X on Y (`crossedX`) and Y on X (`crossedY`) are equal.
#' * `autoregX` and `autoregY`: Tests the hypothesis that the autoregressive effect of X and Y, respectively, is equal across waves. 
#' * `crossedX` and `crossedY`: Tests the hypothesis that the crossed effect of X on Y (`crossedX`) and Y on X (`crossedY`), respectively, is equal across waves. 
#' * `corXY`: Tests the hypothesis that the (residual-)correlations between X and Y are equal across waves. 
#' 
#' For hypotheses regarding the random-intercept CLPM, see [semPower.powerRICLPM()].
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
#' Optional arguments if a **simulated power analysis** (`simulatedPower = TRUE`) is requested:
#' * `nReplications`: The number of simulation runs. Defaults to 250, but larger numbers greatly improve accuracy at the expense of increased computation time.
#' * `minConvergenceRate`: The required minimum convergence rate. Defaults to .50.
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
#' summary(powerCLPM$power)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerCLPM$modelH1, sample.cov = powerCLPM$Sigma,
#'             sample.nobs = powerCLPM$power$requiredN, sample.cov.rescale = FALSE)
#' lavaan::sem(powerCLPM$modelH0, sample.cov = powerCLPM$Sigma,
#'             sample.nobs = powerCLPM$power$requiredN, sample.cov.rescale = FALSE)
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
#' # X1 and X2 are measured by 5 indicators each loading by .4, .5, .6, .5, .4, .6 
#' # Y1 and Y2 are measured by 3 indicators each loading by .8, .6, .7
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'crossedX = 0',
#'                                 loadings = list(
#'                                   c(.4, .5, .6, .5, .4, .6),    # X1
#'                                   c(.8, .6, .7),                # Y1
#'                                   c(.4, .5, .6, .5, .4, .6),    # X2
#'                                   c(.8, .6, .7)                 # Y2
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
#' # Y1, Y2, and X3 are measured by 3 indicators loading by .6 each, and
#' # the synchronous correlation between X and Y are .2 across all three waves, and
#' # the stability of X is .8 across all three waves,
#' # the stability of Y is .7 across all three waves, and
#' # the crossed-effect of Y (Y1 -> X2, and Y2 -> Y3) is .1.
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 3,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = c(.2, .2, .2),
#'                                 waveEqual = c('autoregX', 'autoregY', 'crossedX', 'crossedY', 'corXY'),
#'                                 nullEffect = 'crossedX = 0',
#'                                 nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#' 
#' # Determine required N in a 3-wave CLPM to detect that 
#' # the crossed-effect of X in waves 1 (X1 -> Y2) of .20 is equal to the 
#' # the crossed-effect of X in waves 2 (X2 -> Y3) of .05 
#' # with a power of 95% on alpha = 5%, where
#' # the autoregressive effects of X and Y are equal over waves,
#' # X1, X2, and X3 are measured by 5 indicators loading by .5 each, and
#' # Y1, Y2, and X3 are measured by 3 indicators loading by .6 each, and
#' # the synchronous correlation between X and Y are .2, .3, and .4 at the first, second, and third wave, and
#' # the stability of X is .8 across all three waves,
#' # the stability of Y is .7 across all three waves, and
#' # the crossed-effects of Y (Y1 -> X2, and Y2 -> Y3) are both .1 (but freely estimated for each wave).
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 3,
#'                                 autoregEffects = c(.8, .7), 
#'                                 crossedEffects = list(
#'                                   #  X    Y
#'                                   c(.20, .10),   # wave1 -> wave2
#'                                   c(.05, .10)),  # wave2 -> wave3
#'                                 rXY = c(.2, .3, .4),
#'                                 nullEffect = 'crossedX',
#'                                 waveEqual = c('autoregX', 'autoregY'),
#'                                 nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but determine N to detect that 
#' # the crossed-effect of X at wave 2 is >= .1.
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 3,
#'                                 autoregEffects = c(.8, .7), 
#'                                 crossedEffects = list(
#'                                   #  X    Y
#'                                   c(.20, .10),   # wave1 -> wave2
#'                                   c(.05, .10)),  # wave2 -> wave3
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
#' # the residual correlation between X and Y at wave 3 (of .4) differs from 
#' powerCLPM <- semPower.powerCLPM(type = 'a-priori',
#'                                 nWaves = 3,
#'                                 autoregEffects = c(.8, .7), 
#'                                 crossedEffects = list(
#'                                    #  X    Y
#'                                   c(.20, .10),   # wave1 -> wave2
#'                                   c(.05, .10)),  # wave2 -> wave3
#'                                 rXY = c(.2, .3, .4),
#'                                 nullEffect = 'corXY',
#'                                 waveEqual = c('autoregX', 'autoregY'),
#'                                 nIndicator = c(5, 3, 5, 3, 5, 3),
#'                                 loadM = c(.5, .6, .5, .6, .5, .6),
#'                                 alpha = .05, beta = .05)
#' 
#' # request a simulated post-hoc power analysis with 500 replications.
#' powerCLPM <- semPower.powerCLPM(type = 'post-hoc',
#'                                 nWaves = 2,
#'                                 autoregEffects = c(.8, .7),
#'                                 crossedEffects = c(.2, .1),
#'                                 rXY = NULL,
#'                                 nullEffect = 'crossedX = 0',
#'                                 Lambda = diag(4),
#'                                 alpha = .05, N = 500, 
#'                                 simulatedPower = TRUE, nReplications = 500)
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
  checkEllipsis(...)
  
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
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param nWaves number of waves, must be >= 3.
#' @param autoregEffects vector of the autoregressive effects of X and Y (constant across waves), or a list of vectors of autoregressive effects for X and Y from wave to wave, e.g. `list(c(.7, .6), c(.5, .5))` for an autoregressive effect of .7 for X1->X2 and .6 for X2->X3 and autoregressive effects of .5 for Y1->Y2 and Y2 -> Y3
#' @param crossedEffects vector of crossed effects of X on Y (X -> Y) and vice versa (both constant across waves), or a list of vectors of crossed effects giving the crossed effect of X on Y (and vice versa) for each wave, e.g. `list(c(.2, .3), c(.1, .1))` for X1->Y2 = .2, X2->Y3 = .3, Y1->Y2 = .1, and Y2->Y3 = .1.
#' @param rXY vector of (residual-)correlations between X and Y for each wave. If `NULL`, all (residual-)correlations are zero. 
#' @param rBXBY correlation between random intercept factors.
#' @param waveEqual parameters that are assumed to be equal across waves in both the H0 and the H1 model. Valid are `'autoregX'` and `'autoregY'` for autoregressive effects, `'crossedX'` and `'crossedY'` for crossed effects, `'corXY'` for residual correlations, or `NULL` for none (so that all parameters are freely estimated). 
#' @param nullEffect defines the hypothesis of interest. Valid are the same arguments as in `waveEqual` and additionally `'autoregX = 0'`, `'autoregY = 0'`, `'crossedX = 0'`, `'crossedY = 0'` to constrain the X or Y autoregressive effects or the crossed effects to zero, `'corBXBY = 0'` to constrain the correlation between the random intercepts to zero, and `'autoregX = autoregY'` and `'crossedX = crossedY'` to constrain them to be equal for X and Y.
#' @param nullWhich used in conjunction with `nullEffect` to identify which parameter to constrain when there are > 2 waves and parameters are not constant across waves. For example, `nullEffect = 'autoregX = 0'` with `nullWhich = 2` would constrain the second autoregressive effect for X to zero.    
#' @param metricInvariance whether metric invariance over waves is assumed (`TRUE`, the default) or not (`FALSE`). This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model, where the df are not affected by invariance constraints).
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. The order of factors is (X1, Y1, X2, Y2, ..., X_nWaves, Y_nWaves). See details.
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
  checkEllipsis(...)
  
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


  ### define H1 and H0 model
  
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
  
  ## add constraints to H0 model
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
#' Convenience function for performing power analyses for hypothesis arising 
#' in multigroup measurement invariance models concerning a specific level of invariance.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, either `'saturated'` or one of `'configural'`, `'metric'`, `'scalar'`, or a vector of restrictions in `lavaan` format (with `'none'` for no restrictions). See details.
#' @param nullEffect defines the hypothesis (i.e., level of invariance) of interest. One of `'metric'`, `'scalar'`, `'residual'`, or a vector of restrictions in `lavaan` format. See details.   
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. See details.
#' @return A list containing the following components is returned:
#' \item{`power`}{the results of the power analysis. Use the `summary` method to obtain formatted results.}
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' @details
#' This function performs a power analysis to reject various hypotheses arising
#' in the context of multigroup measurement invariance. Multigroup invariance models 
#' fit the specified model simultaneously to various groups and place increasingly
#' restrictive cross-group equality constrains on the model parameters. The typical - but not in all parts necessary -
#' sequence is (a) configural, (b) metric, (c) scalar, and (d) residual invariance, where each level of invariance is
#' compared against the previous level (e.g., scalar vs. metric). Power analysis provides  
#' the power (or the required N) to reject a particular level of invariance.
#'  
#' The models defined in the `comparison` and the `nullEffect` arguments can be specified as follows:
#' \itemize{
#' \item `configural`: no invariance constraints. Shows the same fit as the saturated model, so only the delta df differ. 
#' \item `metric`: all loadings are restricted to equality. 
#' \item `scalar`: all loadings and (item-)intercepts are restricted to equality. 
#' \item `residual`: all loadings, (item-)intercepts, and (item-)residuals are restricted to equality.
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
#' \item `c('loadings', 'intercepts')`: all loadings and (item-)intercepts are restricted to equality. 
#' \item `c('loadings', 'intercepts', 'residuals')`: all loadings, (item-)intercepts, and (item-)residuals are restricted to equality.
#' \item `c('loadings', 'residuals')`: all loadings and (item-)residuals are restricted to equality.
#' \item `c('loadings', 'intercepts', 'means')`: all loadings, (item-)intercepts, and latent factor means are restricted to equality.
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
#' * `tau`: Defines the item intercepts, required whenever a model involves hypotheses about means (e.g., scalar invariance). If `NULL` and `Alpha` is set, all intercepts are assumed to equal zero.
#' * `Alpha`: Defines the latent means, required whenever a model involves hypotheses about latent means (e.g., latent mean invariance). If `NULL` and `tau` is set, all latent means are assumed to equal zero. Because variance scaling is used so that all factor variances are 1, latent mean differences can be interpreted akin to Cohen's d as standardized mean differences.
#' 
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` always need to be defined, 
#' and `tau` and `Alpha` need to be defined for particular levels of invariance. 
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
#' Optional arguments if a **simulated power analysis** (`simulatedPower = TRUE`) is requested:
#' * `nReplications`: The number of simulation runs. Defaults to 250, but larger numbers greatly improve accuracy at the expense of increased computation time.
#' * `minConvergenceRate`: The required minimum convergence rate. Defaults to .50.
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
#' summary(powerMI$power)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerMI$modelH1, sample.cov = list(powerMI$Sigma[[1]], powerMI$Sigma[[2]]),
#'             sample.nobs = as.list(powerMI$power$requiredN.g), sample.cov.rescale = FALSE)
#' lavaan::sem(powerMI$modelH0, sample.cov = list(powerMI$Sigma[[1]], powerMI$Sigma[[2]]),
#'             sample.nobs = as.list(powerMI$power$requiredN.g), sample.cov.rescale = FALSE)
#' 
#' # same as above, but determine power with N = 500 in each group on alpha = .05
#' powerMI <- semPower.powerMI(type = 'post-hoc',
#'                             comparison = 'configural', 
#'                             nullEffect = 'metric',
#'                             nIndicator = 5,
#'                             loadM = list(.5, .6),
#'                             alpha = .05, N = list(500, 500))
#' 
#' # same as above, but determine the critical chi-square with N = 500 in each group so that alpha = beta
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
#' powerMI <- semPower.powerMI(type = 'post-hoc',
#'                             comparison = 'configural', 
#'                             nullEffect = 'metric',
#'                             nIndicator = list(5, 5),
#'                             loadM = list(.5, .6),
#'                             alpha = .05, N = list(500, 500), 
#'                             simulatedPower = TRUE, nReplications = 500)
#'                              
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerMI <- function(type, 
                             comparison = NULL,
                             nullEffect = NULL,
                             ...){
  
  # validate input
  checkEllipsis(...)
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

#' semPower.powerPath
#'
#' Convenience function for performing power analyses for hypothesis arising 
#' in a generic path model.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param Beta matrix of regression slopes between latent variables (all-y notation). A list for multiple group models. See details. 
#' @param Psi variance-covariance matrix of latent (residual) factors. If `NULL`, a diagonal matrix is assumed. A list for multiple group models. See details.
#' @param nullEffect defines the hypothesis of interest, must be one of `'beta = 0'` (the default) to test whether a regression slope is zero, `'betaX = betaZ'` to test for the equality of slopes, and `'betaX = betaZ'` to test for the equality of a slope across groups. Define the slopes to be set to equality in `nullWhich` and the groups in `nullWhichGroups`. 
#' @param nullWhich vector of size 2 indicating which slope in `Beta` is hypothesized to equal zero when `nullEffect = 'beta = 0'`, or to restrict to equality across groups when `nullEffect = 'betaA = betaB'`, or list of vectors defining which correlations to restrict to equality when `nullEffect = 'betaX = betaZ'`. Can also contain more than two slopes, e.g., `list(c(2, 1), c(3, 1), c(3, 2))` to set `Beta[2, 1] = Beta[3, 1] = Beta[3, 2]`.
#' @param nullWhichGroups for `nullEffect = 'betaA = betaB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. See details.
#' @return A list containing the following components is returned:
#' \item{`power`}{the results of the power analysis. Use the `summary` method to obtain formatted results.}
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
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
#' \deqn{\Sigma = \Lambda (I - \Beta)^{-1} \Psi [(I - \Beta)^{-1}]'  \Lambda' + \Theta } 
#' where \eqn{\Beta} is the \eqn{m \cdot m} matrix containing the regression slopes and \eqn{\Psi} is the (residual) variance-covariance matrix of the \eqn{m} factors. 
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
#' Optional arguments if a **simulated power analysis** (`simulatedPower = TRUE`) is requested:
#' * `nReplications`: The number of simulation runs. Defaults to 250, but larger numbers greatly improve accuracy at the expense of increased computation time.
#' * `minConvergenceRate`: The required minimum convergence rate. Defaults to .50.
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
#' summary(powerPath$power)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerPath$modelH1, sample.cov = powerPath$Sigma,
#' sample.nobs = powerPath$power$requiredN, sample.cov.rescale = FALSE)
#' lavaan::sem(powerPath$modelH0, sample.cov = powerPath$Sigma,
#' sample.nobs = powerPath$power$requiredN, sample.cov.rescale = FALSE)
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
#' Beta1 <- matrix(c(
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
                               ...){
  
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
  if(!is.null(nullWhichGroups)) lapply(nullWhichGroups, function(x) checkBounded(x, 'All elements in nullWhichGroups'), bound = c(1, length(Beta)), inclusive = TRUE)
  
  ### get Sigma
  generated <- semPower.genSigma(Beta = if(!isMultigroup) Beta[[1]] else Beta, 
                                 Psi = if(!isMultigroup) Psi[[1]] else Psi, 
                                 useReferenceIndicator = TRUE, ...)
  
  ### create model strings
  # we need labels, so we use modelTrueCFA and define regression relations here
  if(!isMultigroup) model <- generated[['modelTrueCFA']] else model <- generated[[1]][['modelTrueCFA']]
  tok <- list()
  for(f in 1:ncol(Beta[[1]])){
    if(!isMultigroup){
      idx <- which(Beta[[1]][f, ] != 0)
      if(length(idx) > 0){
        clab <- paste0('pf', paste0(formatC(f, width = 2, flag = 0), formatC(idx, width = 2, flag = 0)))
        tok <- append(tok, paste0('f', f, ' ~ ', paste(clab, paste0('*f', idx), sep = '', collapse = ' + ')))
      }
    }else{
      idx <- unique(unlist(lapply(Beta, function(x) which(x[f, ] != 0))))
      if(length(idx) > 0){
        clab <- paste0('pf', paste0(formatC(f, width = 2, flag = 0), formatC(idx, width = 2, flag = 0)))
        clab <- unlist(lapply(clab, function(x) paste0('c(', paste0(x, 'g', seq_along(Beta), collapse = ', '),')')))
        tok <- append(tok, paste0('f', f, ' ~ ', paste(clab, paste0('*f', idx), sep = '', collapse = ' + ')))
      }
    }
  }
  # (residual) correlations, assuming the same Psi across groups
  if(is.null(Psi)) if(isMultigroup) Psi <- list(generated[[1]][['Psi']]) else Psi <- list(generated[['Psi']])
  for(f in 1:(ncol(Beta[[1]]) - 1)){
    for(ff in (f + 1):ncol(Beta[[1]])){
      if(Psi[[1]][f, ff] != 0)
        tok <- append(tok, paste0('f', f, ' ~~ f', ff))
    }
  }
  modelH1 <- paste(c(model, unlist(tok)), collapse = '\n')

  # define H0 model
  if(nullEffect == 'beta=0'){
    modelH0 <- paste(modelH1, 
                     paste0('pf', paste0(formatC(nullWhich[1], width = 2, flag = 0), formatC(nullWhich[2], width = 2, flag = 0)), ' == 0'),
                     sep = '\n')
  }else if(nullEffect == 'betax=betaz'){
    tok <- paste0('pf', paste0(formatC(nullWhich[[1]][1], width = 2, flag = 0), formatC(nullWhich[[1]][2], width = 2, flag = 0)))
    tok <- paste0(tok, ' == ','pf', paste0(formatC(nullWhich[[2]][1], width = 2, flag = 0), formatC(nullWhich[[2]][2], width = 2, flag = 0)))
    modelH0 <- paste(modelH1, tok, sep = '\n')
  }else if(nullEffect == 'betaa=betab'){
    if(is.null(nullWhichGroups)) nullWhichGroups <- seq_along(Beta)
    clab <- paste0('pf', paste0(formatC(nullWhich[1], width = 2, flag = 0), formatC(nullWhich[2], width = 2, flag = 0)))
    clab <- lapply(seq_along(nullWhichGroups), function(x) paste0(clab, 'g', nullWhichGroups[x]))
    tok <- list()
    for(i in 1:(length(clab) - 1)){
      for(j in (i + 1):length(clab)){
        tok <- append(tok, paste0(clab[i], ' == ', clab[j]))
      }
    }
    modelH0 <- paste(modelH1, tok, sep = '\n')
  }else{
    stop('nullEffect not defined.')
  }

  # always enforce invariance constraints in the multigroup case
  lavOptions <- NULL
  if(isMultigroup) lavOptions <- list(group.equal = c('loadings'))

  # always fit H1 model
  fitH1model <- TRUE 
  if(comparison == 'saturated'){
    modelH1 <- NULL
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