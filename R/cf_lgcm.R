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
#' For hypotheses regarding longitudinal invariance, see [semPower.powerLI()]. For hypotheses regarding a simple autoregression, see [semPower.powerAutoreg()]. For hypotheses in an ARMA model, see [semPower.powerARMA()].
#' 
#' Note that power analyses concerning the hypotheses `iVar = 0`, `sVar = 0`, and `s2Var = 0` are only approximate, 
#' because the H0 model involves a parameter constraint on the boundary of the parameter space (a variance of zero), 
#' so that the correct limiting distribution is a mixture of non-central \eqn{\chi^2} distributions 
#' (see Stoel et al., 2006). In effect, power is (slightly) underestimated.
#' 
#' Stoel, R. D., Garre, F. G., Dolan, C., & Van Den Wittenboer, G. (2006). On the likelihood ratio test in structural equation modeling when parameters are subject to boundary constraints. *Psychological Methods, 11*, 439-455.
#'
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
#' * `nCores`: The number of cores to use for parallel processing. Defaults to 1 (= no parallel processing). This requires the `doFuture` package.
#' * `futureStrategy`: A string specifying the strategy used in parallel processing (when `nCores` >  1). Defaults to `'multisession'`. This is passed to the `plan` method of the `doFuture` package. See the `doFuture` package for valid options.
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
  if(any(unlist(lapply(covariances, function(x) any(abs(cov2cor(x)) > 1))))) warning('Provided variances and covariances imply a correlation larger 1 or smaller -1.')
  
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
