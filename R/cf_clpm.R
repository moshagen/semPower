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
#' @param rXY vector of (residual-)correlations between X and Y for each wave. If `NULL`, all (residual-)correlations are zero. Can be a list for multiple groups models, otherwise no group differences are assumed.
#' @param waveEqual parameters that are assumed to be equal across waves in both the H0 and the H1 model. Valid are `'autoregX'` and `'autoregY'` for autoregressive effects, `'crossedX'` and `'crossedY'` for crossed effects, `'corXY'` for residual correlations, or `NULL` for none (so that all parameters are freely estimated, subject to the constraints defined in `nullEffect`). 
#' @param nullEffect defines the hypothesis of interest. Valid are the same arguments as in `waveEqual` and additionally `'autoregX = 0'`, `'autoregY = 0'`, `'crossedX = 0'`, `'crossedY = 0'` to constrain the X or Y autoregressive effects or the crossed effects to zero, `'autoregX = autoregY'` and `'crossedX = crossedY'` to constrain them to be equal for X and Y, and `'autoregXA = autoregXB'`, `'autoregYA = autoregYB'`, `'crossedXA = crossedXB'`, `'crossedYA = crossedYB'` to constrain them to be equal across groups. 
#' @param nullWhich used in conjunction with `nullEffect` to identify which parameter to constrain when there are > 2 waves and parameters are not constant across waves. For example, `nullEffect = 'autoregX = 0'` with `nullWhich = 2` would constrain the second autoregressive effect for X to zero.    
#' @param nullWhichGroups for hypothesis involving cross-groups comparisons, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param standardized whether all parameters should be treated as standardized (`TRUE`, the default), implying that unstandardized and standardized regression relations have the same value. If `FALSE`, all regression relations are unstandardized.
#' @param standardizedResidualCovariances whether the residual covariances provided in `rXY` should be interpreted as correlations. When `TRUE` (the default) the unstandardized residual covariances differ from the those provided in `rXY`. When `FALSE`, the values provided in `rXY` are the unstandardized residual covariances, and the standardized residual correlations differ.
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
#' For hypotheses regarding the random-intercept CLPM, see [semPower.powerRICLPM()]. For hypothesis in autoregressive models, see [semPower.powerAutoreg()].
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
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.
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
                               standardizedResidualCovariances = TRUE,
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
  ngA <- ifelse(is.list(autoregEffects[[1]]), length(autoregEffects), 1)
  ngX <- ifelse(is.list(crossedEffects[[1]]), length(crossedEffects), 1)
  ngR <- ifelse(is.list(rXY), length(rXY), 1)
  ig <- unique(c(ngA, ngX, ngR))
  if(length(ig[ig > 1]) > 1) stop('Non-null list arguments supplied to autoregEffects, crossedEffects, or rXY imply a different number of groups. Make sure all lists have the same length or provide no list for no group differences.')
  nGroups <- max(ig)
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
  if('corxy' %in% nullEffect && standardizedResidualCovariances) stop('Power analysis for nullEffect == "corxy" can only be performed when the residual covariances are not standardized. Repeat with standardizedResidualCovariances = FALSE.') 
  
  
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
  
  if(standardized){ 
    Psi <- lapply(seq(Beta), function(x) getPsi.B(Beta[[x]], Psi[[x]], standResCov = standardizedResidualCovariances))
  }
  
  # add metric invariance constraints to analysis model
  metricInvarianceFactors <- NULL
  if(metricInvariance){
    metricInvarianceFactors <- list(
      seq(1, 2*nWaves, 2),
      seq(2, 2*nWaves, 2)  
    )
  }
  
  ### get Sigma
  generated <- semPower.genSigma(Beta = if(!isMultigroup) Beta[[1]] else Beta,
                                 Psi = if(!isMultigroup) Psi[[1]] else Psi,
                                 useReferenceIndicator = TRUE,
                                 metricInvariance = metricInvarianceFactors,
                                 nGroups = nGroups,
                                 ...)
  
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
