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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.
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
