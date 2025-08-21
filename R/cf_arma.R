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
#' For hypotheses regarding a simple autoregression, see [semPower.powerAutoreg()]. For hypotheses regarding a CLPM structure, see [semPower.powerCLPM()].  For hypotheses regarding longitudinal measurement invariance, see [semPower.powerLI()].
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.
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
  
  Lambda <- cbind(matrix(0, nrow = nrow(Lambda), ncol = nWaves), Lambda)  # add noise factors
  Lambda <- rep(list(Lambda), nGroups)  # require same measurement model across groups
  
  ### create Beta
  # n1 - n_nwaves, f1 - f_nwaves,  
  Beta <- lapply(seq(nGroups), function(x){
    B <- matrix(0, ncol = 2*nWaves, nrow = 2*nWaves)
    # add autoregEffects
    for(i in seq(nWaves - 1)){
      idx <- nWaves + i
      B[(idx + 1), idx] <- autoregEffects[[x]][i]
    }
    # lag-2 effects
    if(nWaves > 2){
      for(i in seq(nWaves - 2)){
        idx <- nWaves + i
        B[(idx + 2), idx] <- autoregLag2[[x]][i]
      }
    }
    # lag-3 effects
    if(nWaves > 3){
      for(i in seq(nWaves - 3)){
        idx <- nWaves + i
        B[(idx + 3), idx] <- autoregLag3[[x]][i]
      }
    }
    # lag-1 mov avgs
    for(i in seq(nWaves - 1)){
      idx <- nWaves + i
      B[(idx + 1), i] <- mvAvgLag1[[x]][i]
    }
    # lag-2 mov avgs
    if(nWaves > 2){
      for(i in seq(nWaves - 2)){
        idx <- nWaves + i
        B[(idx + 2), i] <- mvAvgLag2[[x]][i]
      }
    }
    # lag-3 mov avgs
    if(nWaves > 3){
      for(i in seq(nWaves - 3)){
        idx <- nWaves + i
        B[(idx + 3), i] <- mvAvgLag3[[x]][i]
      }
    }
    # noise factors
    for(i in seq(nWaves)){
      idx <- nWaves + i
      B[idx, i] <- 1
    }
    B
  })
  
  ### create Psi
  Psi <- lapply(seq(nGroups), function(x){
    P <- matrix(0, ncol = 2*nWaves, nrow = 2*nWaves)
    # var noise factors
    diag(P[seq(nWaves), seq(nWaves)]) <- variances[[x]]
    P
  })
  
  ### define latent means
  Alpha <- NULL
  if(!is.null(means)){
    Alpha <- lapply(means, function(x) c(rep(0, nWaves), x)) # noise factors + F
  }
  
  ### add metric invariance constraints to analysis model
  metricInvarianceFactors <- NULL
  if(invariance) metricInvarianceFactors <- list(seq((nWaves + 1), 2*nWaves))
  
  
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
  
  # correct factor labels
  sf <- paste0('f', seq((nWaves + 1), 2*nWaves))
  repl <- paste0('f', seq(nWaves))
  for(i in seq(sf)){
    model <- gsub(sf[i], repl[i], model)
  }
  
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