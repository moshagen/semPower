#' semPower.powerLI
#'
#' Convenience function for performing power analyses for hypothesis arising 
#' in longitudinal measurement invariance models concerning a specific level of invariance.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, either `'saturated'` or one of `'configural'`, `'metric'`, `'scalar'`, `'residual'`, `'covariances'`, `'means'`, or a vector of restrictions in `lavaan` format (with `'none'` for no restrictions). See details.
#' @param nullEffect defines the hypothesis (i.e., level of invariance) of interest. Accepts the same arguments as `comparison`. See details.   
#' @param useReferenceIndicator whether to identify the factor variances and factor means using a referent indicator (`TRUE`) or using a standardization approach setting the means to 0 and the variances to 1 (`FALSE`, the default). See details.
#' @param singleOccasionIdent when `useReferenceIndicator = FALSE`, whether to maintain the identification constraints only for a single measurement (`TRUE`) or whether to maintain all scaling constraints (`FALSE`). See details.
#' @param autocorResiduals whether the residuals of the indicators of latent variables are autocorrelated over waves (`TRUE`, the default) or not (`FALSE`). This affects the df when the comparison model is the saturated model and generally affects power (also for comparisons to the restricted model).
#' @param Phi the factor covariance matrix. Can be `NULL` for uncorrelated factors and unit variances.
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
#' and concerning the structural part  (e) latent variances, (f) latent covariances, and (g) latent means, where each level of invariance is
#' compared against the previous level (e.g., scalar vs. metric). Power analysis provides  
#' the power (or the required N) to reject a particular level of invariance.
#' 
#' For hypotheses regarding multiple group invariance, see [semPower.powerMI()]. For hypotheses regarding autoregressive models, see [semPower.powerAutoreg()]. For hypotheses in an ARMA model, see [semPower.powerARMA()].
#'  
#' The general model structure posits at least one factor measured at two ore more occasions. The invariance constraints over occasions are always 
#' placed separately for each factor. For example, if two factors are measured twice, equality constraints on loadings are placed over the measurements of each factor separately.
#' There are two ways to specify the models defined in the `comparison` and the `nullEffect` arguments. Either, one may
#' specify a specific level of invariance that includes all previous levels:
#' \itemize{
#' \item `'configural'`: no invariance constraints. Shows the same fit as the saturated model, so only the delta df differ. 
#' \item `'metric'`: all loadings (of a particular attribute) are restricted to equality over measurement occasions. Note that if reference scaling is used, the first indicator should be invariant.
#' \item `'scalar'`: all loadings and (indicator-)intercepts are restricted to equality. Note that if reference scaling is used, the first indicator should be invariant.
#' \item `'residual'`: all loadings, (indicator-)intercepts, and (indicator-)residuals are restricted to equality.
#' \item `'variances'`: all loadings, (indicator-)intercepts, (indicator-)residuals, and latent variances are restricted to equality.
#' \item `'covariances'`: all loadings, (indicator-)intercepts, (indicator-)residuals, latent variances, and latent covariances are restricted to equality. This only refers to the factor covariance(s) at each measurement occasion, but does not involve cross-temporal relations such as stabilities.
#' \item `'means'`: all loadings, (indicator-)intercepts, (indicator-)residuals, latent variances, latent covariances, and latent means are restricted to equality.
#' }
#' 
#' For example, setting `comparison = 'metric'` and `nullEffect = 'scalar'` determines power 
#' to reject the hypothesis that the constraints placed in the scalar invariance model 
#' (restricting loadings and intercepts over measurement occasions) over the 
#' metric invariance model (restricting only the loadings over measurement occasions) are defensible.
#'  
#' For greater flexibility, the models can also be defined using `lavaan` style restrictions as a vector, namely
#' `'none'` (no restrictions), `'loadings'` (loadings), `'intercepts'` (intercepts), `'residuals'` (residuals), `'lv.variances'` (latent variances), `'lv.covariances'` (latent covariances), `'means'` (latent means).
#'  For instance: 
#' \itemize{
#' \item `'none'`: no invariance constraints and thus representing a configural invariance model. Shows the same fit as the saturated model, so only the delta df differ. 
#' \item `c('loadings')`: all loadings (of a particular attribute) are restricted to equality. Note that if reference scaling is used, the first indicator should be invariant. 
#' \item `c('loadings', 'intercepts')`: all loadings and (indicator-)intercepts are restricted to equality. Note that if reference scaling is used, the first indicator should be invariant.
#' \item `c('loadings', 'intercepts', 'residuals')`: all loadings, (indicator-)intercepts, and (indicator-)residuals are restricted to equality.
#' \item `c('loadings', 'residuals')`: all loadings and (indicator-)residuals are restricted to equality.
#' \item `c('loadings', 'intercepts', 'means')`: all loadings, (indicator-)intercepts, and latent factor means are restricted to equality.
#' \item `c('loadings', 'residuals', 'lv.variances')`: all loadings, (indicator-)residuals, and latent factor variances are restricted to equality.
#' }
#' 
#' For example, setting `comparison = c('loadings')` and `nullEffect = 'c('loadings', 'intercepts')'` 
#' determines power to reject the hypothesis that the constraints placed in the scalar invariance model 
#' (restricting loadings and intercepts) over the  metric invariance model (restricting only the loadings) are defensible.
#' 
#' By default, the factors are identified by setting their means to 0 and their variances to 1 (`useReferenceIndicator = FALSE`), so that all loadings and indicator intercepts are freely estimated. 
#' In models involving constraints on loadings or intercepts, a single occasion identification approach is used by default (`singleOccasionIdent = TRUE`), so that the factor variances and means are constrained at only a single occasion, but are freely estimated at all remaining occasions. This mimics the `lavaan` approach to invariance testing.
#' If setting `singleOccasionIdent = FALSE`, the identification constraints are maintained in models assuming metric or scalar invariance, so that the factor means and variances are constrained to 0 and 1, respectively, across all occasions (so invariance of factor variances and factor means is always met). Concerning metric invariance, this was the default behavior in `semPower` versions prior 2.1.4. 
#' A referent scaling approach can be be defined using `useReferenceIndicator = TRUE`, so that the loading and the intercept of the first indicator of each factor are constrained at all occasions, and the factor variances and means are freely estimated. 
#' Reference scaling is always maintained in all invariance models, so `singleOccasionIdent` has no effect (i.e. is always `FALSE`)
#' 
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model**:
#' * `Lambda`: The factor loading matrix, where the order is F1_T1, F2_T1, ..., F1_T2, F2_T2, ..., F1_TZ, ..., FX_TZ
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.5, .7, .7))` defines two factors with three indicators loading on the first factor by .5, , 4., and .6, and three indicators loading on the second factor by .5, .7, .7.
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators by factor, e. g., `nIndicator = c(3, 3)` defines a two factor model with three indicators each, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the loadings of the first factor to equal .5 and those of the second factor do equal .6.
#' * `Theta`: Variance-covariance matrix of the indicator residuals, which should be a diagonal matrix. Required when residual non-invariance is to be detected. When `NULL`, Theta is a diagonal matrix with elements such that all variances are 1. 
#' * `tau`: Defines the indicator intercepts, required whenever a model involves hypotheses about means (e.g., scalar invariance). If `NULL` and `Alpha` is set, all intercepts are assumed to equal zero.
#' * `Alpha`: Defines the latent means, required whenever a model involves hypotheses about latent means (e.g., latent mean invariance). If `NULL` and `tau` is set, all latent means are assumed to equal zero. If variance scaling is used so that all factor variances are 1, latent mean differences can be interpreted akin to Cohen's d as standardized mean differences.
#' 
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` always need to be defined, 
#' and `Theta`, `tau` and `Alpha` need to be defined for particular levels of invariance.
#' 
#' If longitudinal invariance for a single factor measured at T occasions is of interest, `loadings` can be a single list with T elements, and `nIndicator`, `loadM`, `tau` and `Alpha` can be single vectors. 
#' If longitudinal invariance for a multiple factor measured at T occasions is of interest, `loadings`, `nIndicator`, `loadM`, `tau` and `Alpha` must be wrapped in a list with T elements, where each element refers to a single measurement occasion. 
#' `Phi` and `Theta` always must be a single matrix. The order of `Phi` must be F1_T1, F2_T1, ..., F1_T2, F2_T2, ..., F1_TZ, F2_TZ, ..., FX_TZ. The order of `Theta` must match
#' the order of `Phi` and thus first contains the respective indicators of F1 at T1, F2 at T1, ..., F1 at T2, F2 at T2, and so on.
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
#'   loadings = list(
#'     c(.7, .6, .5, .7),   # factor 1 at T1
#'     c(.4, .4, .4, .4)    # factor 1 at T2
#'   ),
#'   autocorResiduals = TRUE,
#'   useReferenceIndicator = FALSE,
#'   singleOccasionIdent = TRUE
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
#'   loadings = list(
#'     c(.7, .6, .5, .7),   # factor 1 at T1
#'     c(.4, .4, .4, .4)    # factor 1 at T2
#'   ),
#'   autocorResiduals = TRUE,
#'   useReferenceIndicator = FALSE,
#'   singleOccasionIdent = TRUE
#' )
#' 
#' 
#' # same as above, but determine the critical chi-square with N = 500 in each
#' # group so that alpha = beta
#' powerLI <- semPower.powerLI(
#'   type = 'compromise', abratio = 1, N = 500,
#'   comparison = 'configural',
#'   nullEffect = 'metric',
#'   loadings = list(
#'     c(.7, .6, .5, .7),   # factor 1 at T1
#'     c(.4, .4, .4, .4)    # factor 1 at T2
#'   ),
#'   autocorResiduals = TRUE,
#'   useReferenceIndicator = FALSE,
#'   singleOccasionIdent = TRUE
#' )
#' 
#' 
#' # same as above, but compare to the saturated model
#' # (rather than to the configural invariance model)
#' powerLI <- semPower.powerLI(
#'   type = 'a-priori', alpha = .05, power = .80,
#'   comparison = 'saturated',
#'   nullEffect = 'metric',
#'   loadings = list(
#'     c(.7, .6, .5, .7),   # factor 1 at T1
#'     c(.4, .4, .4, .4)    # factor 1 at T2
#'   ),
#'   autocorResiduals = TRUE,
#'   useReferenceIndicator = FALSE,
#'   singleOccasionIdent = TRUE
#' )
#' 
#' 
#' # define two factors measured at three occasions using a
#' # reduced loading matrix, and provide factor covariance matrix
#' Phi <- matrix(c(
#'   # F1T1 F2T1 F1T2 F2T2 F1T3 F2T3 
#'   c(1,   .1,  .7,   0,  .5,   0),     # F1_T1 
#'   c(.1,   1,  .0,  .6,   0,  .3),     # F2_T1 
#'   c(.7,  .0,   1,  .1,  .7,   0),     # F1_T2 
#'   c(.0,  .6,  .1,   1,   0,  .6),     # F2_T2 
#'   c(.5,  .0,  .7,   0,   1,  .1),     # F1_T3 
#'   c(.0,  .3,   0,  .6,  .1,   1)      # F2_T3
#' ), nrow=6, ncol=6, byrow = TRUE)
#' powerLI <- semPower.powerLI(
#'   type = 'a-priori', alpha = .05, power = .80,
#'   comparison = 'configural',
#'   nullEffect = 'metric',
#'   loadings = list(
#'     list(
#'       c(.8, .8, .8),   # factor 1 at T1
#'       c(.5, .5, .5)    # factor 2 at T1
#'     ),
#'     list(
#'       c(.7, .6, .7),   # factor 1 at T2
#'       c(.4, .6, .4)    # factor 2 at T2
#'     ),
#'     list(
#'       c(.6, .8, .4),   # factor 1 at T3
#'       c(.5, .5, .7)    # factor 2 at T3
#'     )
#'   ),
#'   Phi = Phi,
#'   autocorResiduals = TRUE,
#'   useReferenceIndicator = FALSE,
#'   singleOccasionIdent = TRUE
#' )
#' 
#' # obtain the required N to reject the hypothesis of scalar invariance
#' # in comparison to the metric invariance model
#' # with a power of 80% on alpha = 5%
#' # all intercepts are 0.0 at the first measurement occasion, but
#' # all intercepts are 0.2 at the second measurement occasion and
#' powerLI <- semPower.powerLI(
#'   type = 'a-priori', alpha = .05, power = .80,
#'   comparison = 'metric',
#'   nullEffect = 'scalar',
#'   loadings = list(
#'     c(.7, .6, .5, .7),   # factor 1 at T1
#'     c(.7, .6, .5, .7)    # factor 1 at T2
#'   ),
#'   tau = c(0, 0, 0, 0,
#'           .1, .2, .3, .2),  
#'   autocorResiduals = TRUE,
#'   useReferenceIndicator = FALSE,
#'   singleOccasionIdent = TRUE
#' )
#' 
#' # same as above, but use lavaan strings to define invariance models
#' powerLI <- semPower.powerLI(
#'   type = 'a-priori', alpha = .05, power = .80,
#'   comparison = c('loadings'),
#'   nullEffect = c('loadings', 'intercepts'),
#'   loadings = list(
#'     c(.7, .6, .5, .7),   # factor 1 at T1
#'     c(.7, .6, .5, .7)    # factor 1 at T2
#'   ),
#'   tau = c(0, 0, 0, 0,
#'           .1, .2, .3, .2),  
#'   autocorResiduals = TRUE,
#'   useReferenceIndicator = FALSE,
#'   singleOccasionIdent = TRUE
#' )
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
#'   loadings = list(
#'     c(.7, .6, .5, .7),   # factor 1 at T1
#'     c(.7, .6, .5, .7)    # factor 1 at T2
#'   ),
#'   tau = c(0, 0, 0, 0,
#'           0, 0, 0, 0),  
#'   Alpha = c(0, .5),
#'   autocorResiduals = TRUE,
#'   useReferenceIndicator = FALSE,
#'   singleOccasionIdent = TRUE
#' )
#'    
#' # request a simulated post-hoc power analysis with 250 replications
#' # to reject the hypothesis of equal latent means.
#' set.seed(300121)
#' powerLI <- semPower.powerLI(
#'   type = 'post-hoc', alpha = .05, N = 500,
#'   comparison = 'configural',
#'   nullEffect = 'metric',
#'   loadings = list(
#'     c(.7, .6, .5, .7),   # factor 1 at T1
#'     c(.4, .4, .4, .4)    # factor 1 at T2
#'   ),
#'   autocorResiduals = TRUE,
#'   useReferenceIndicator = FALSE,
#'   singleOccasionIdent = TRUE,
#'   simulatedPower = TRUE,
#'   simOptions = list(nReplications = 250)
#' )
#' 
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerLI <- function(type, 
                             comparison = NULL,
                             nullEffect = NULL,
                             autocorResiduals = TRUE,
                             Phi = NULL,
                             useReferenceIndicator = FALSE,
                             singleOccasionIdent = !useReferenceIndicator,
                             ...){
  
  args <- list(...)
  
  # validate input
  checkEllipsis(...)
  
  # allow both lavstring and single string
  nullValid <- c('metric', 'scalar', 'residual', 'variances', 'covariances', 'means')
  nullValidLav <- c('loadings', 'intercepts', 'residuals', 'lv.variances', 'lv.covariances', 'means')
  
  compValid <- c('saturated', 'configural', 'none', nullValid, nullValidLav)
  comparison <- lapply(comparison, function(x) checkNullEffect(x, compValid, 'comparison'))
  nullEffect <- lapply(nullEffect, function(x) checkNullEffect(x, c(nullValid, nullValidLav), 'nullEffect'))
  
  if(unlist(comparison)[1] == 'configural' || unlist(comparison)[1] == 'none') comparison <- 'configural'
  if(unlist(comparison)[1] == 'saturated') comparison <- 'saturated'
  
  # check and translate to lavstring
  if(length(nullEffect) > 1 || nullEffect == 'loadings'){
    if(any(unlist(lapply(nullEffect, function(x) x %in% nullValid[-6])))) stop('Either use lavaan-type strings or use predefined strings, but do not mix.')
    if(any(unlist(lapply(comparison, function(x) x %in% c(nullValid))))) stop('Either use lavaan-type strings or use predefined strings, but do not mix.')
    if((length(nullEffect) == 1 && nullEffect != 'loadings') && length(nullEffect) <= length(comparison)) stop('The H0 model must contain all restrictions of the comparison model plus one additional restriction.')
    if(!any(c('saturated', 'configural') %in% comparison) && any(unlist(lapply(comparison, function(x) !x %in% nullEffect)))) stop('The H0 model must contain all restrictions of the comparison model plus one additional restriction.')
  }else{
    if(nullEffect %in% nullValidLav[-6]) stop('Either use lavaan-type strings or use predefined strings, but do not mix.')
    if(comparison %in% nullValidLav[-6]) stop('Either use lavaan-type strings or use predefined strings, but do not mix.')
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
    if(idxC > 5) lc <- append(lc, 'lv.variances')
    if(idxN > 5) ln <- append(ln, 'lv.variances')
    if(idxC > 6) lc <- append(lc, 'lv.covariances')
    if(idxN > 6) ln <- append(ln, 'lv.covariances')
    if(idxC > 7) lc <- append(lc, 'means')
    if(idxN > 7) ln <- append(ln, 'means')
    nullEffect <- ln
    if(!any(c('saturated', 'configural') %in% comparison)) comparison <- lc
  }
  comparison <- unlist(comparison)
  nullEffect <- unlist(nullEffect)
  
  if('means' %in% nullEffect && !any(c('saturated', 'intercepts') %in% comparison)) stop('Latent means cannot be estimated without constraints on intercepts (scalar invariance).')
  if(!useReferenceIndicator && !singleOccasionIdent && ('lv.variances' %in% comparison || 'lv.variances' %in% nullEffect)) stop('Variance scaling without single occasion identification is used, so invariance of latent variances is always met. Either use single occasion identification or referent scaling.')
  if(!useReferenceIndicator && !singleOccasionIdent && ('means' %in% nullEffect)) stop('Mean scaling without single occasion identification is used, so invariance of latent means is always met. Either use single occasion identification or referent scaling.')
  if(useReferenceIndicator && singleOccasionIdent) stop('Identification constraints on referent indicators (useReferenceIndicator = TRUE) must be maintained, so singleOccasionIdent must be FALSE.')
  
  
  ### create lambda: here provided list elements refer to measurement occasions, so lambda[oc][fac]. 
  # we need lambda in order  f1_t1, f2_t1, ..., f1_t2, f2_t2, ... 
  tGenerated <- semPower.genSigma(Lambda = args[['Lambda']], loadings = args[['loadings']], loadM = args[['loadM']], loadSD = args[['loadSD']], nIndicator = args[['nIndicator']])
  if(!is.null(tGenerated[['Lambda']])){
    # if no list structure is provided, assume single factor model at each measurement occasion
    Lambda <- tGenerated[['Lambda']]
    nWaves <- ncol(Lambda)
    nFactors <- 1
  }else{
    # create lambda from provided lists
    tLambda <- lapply(tGenerated, '[[', 'Lambda')
    if(length(unique((unlist(lapply(tLambda, ncol))))) != 1) stop('The number of factors must be identical across measurements.')
    if(length(unique((unlist(lapply(tLambda, nrow))))) != 1) stop('The number of indicators must be identical across measurements.')
    nWaves <- length(tLambda)
    nFactors <- ncol(tLambda[[1]])
    nVar <- nrow(tLambda[[1]])
    Lambda <- matrix(0, ncol = nWaves*nFactors, nrow = nWaves*nVar)
    for(i in seq(nWaves)){
      cIdx <- 1 + nFactors * (i - 1)
      cIdx2 <- cIdx + nFactors - 1
      rIdx <- 1 + nVar * (i - 1)
      rIdyx2 <- rIdx + nVar - 1 
      Lambda[rIdx:rIdyx2, cIdx:cIdx2] <- tLambda[[i]]
    }
  }  
  metricInvarianceFactors <- lapply(seq(nFactors), function(x) seq(x, nFactors*nWaves, nFactors))

  ### check phi
  if(is.null(Phi)){
    Phi <- diag(ncol(Lambda))
  }else{
    if(nFactors == 1){
      if(!'means' %in% c(comparison, nullEffect) && 'lv.covariances' %in% c(comparison, nullEffect)) stop('Equality of covariances is only meaningful when there are at least 2 factors at each measurement.')
      if(!is.matrix(Phi) && nWaves != 2) stop('Phi must be matrix when there are more than 2 measurements.')
      if(!is.matrix(Phi)){
        m <- diag(2)
        m[1,2] <- m[2,1] <- Phi
        Phi <- m
      }
    }
    checkSymmetricSquare(Phi, 'must be a symmetric square matrix.')
    if(ncol(Phi) != ncol(Lambda)) stop('The dimensions of Phi must equal the number of factors times the number of measurements.')
  }
  
  ### check tau and Alpha
  tau  <- args[['tau']]
  Alpha  <- args[['Alpha']]
  if(is.null(tau) && any(c('intercepts', 'means') %in% c(comparison, nullEffect))) stop('Tests involving scalar invariance require specification of the intercepts (tau).')
  if(is.null(Alpha) && 'means' %in% c(comparison, nullEffect)) stop('Tests involving mean invariance require specification of the factor means (Alpha).')
  if(is.null(tau)){
    tau <- rep(0, nrow(Lambda))
  }else{
    if(!is.list(tau) && nFactors > 1) stop('tau must be a list of vectors giving the intercepts for each indicator at each measurement.') 
    if(is.list(tau) && nFactors > 1 && length(unique(lapply(tau, length))) != 1) stop('tau must be a list of vectors giving the intercepts for each indicator at each measurement.') 
    tau <- unlist(tau)
    if(length(tau) != nrow(Lambda)) stop('tau must be of length indicators times measurements.')
  }
  if(is.null(Alpha)){
    Alpha <- rep(0, ncol(Lambda))
  }else{
    if(!is.list(Alpha) && nFactors > 1) stop('Alpha must be a list of vectors giving the factor means at each measurement.') 
    if(is.list(Alpha) && nFactors > 1 && length(unique(lapply(Alpha, length))) != 1) stop('Alpha must be a list of vectors giving the factor means at each measurement.') 
    Alpha <- unlist(Alpha)
    if(length(Alpha) != ncol(Lambda)) stop('Alpha must be of length factors times measurements.')
  }
  
  ### throw warning if identification constraints imply misfit (only works in single factor cases)
  if(useReferenceIndicator){
    lapply(metricInvarianceFactors, function(mif){
      mIdx <- sapply(mif, function(x) (which(Lambda[ , x] != 0)))[1, ] # referent loadings indices
      rIdx <- unlist(lapply(seq(mif), function(r) mIdx[r] + (mif[r]-1) * nrow(Lambda))) # adapt indices to vech
      if(any(length(unique(c(Lambda)[rIdx])) != 1)) warning('The loading of least one reference indicator differs across measurements, violating metric invariance. This is probably not intended. Either define the same loading for the first indicator of each factor or use variance scaling (useReferenceIndicator = FALSE).')
      if(any(length(unique(tau[mIdx])) != 1)) warning('The intercepts of least one reference indicator differs across measurements. This is probably not intended. Either define the same intercept for the first indicator of each factor or use mean scaling (useReferenceIndicator = FALSE).')
    })
  }else{
    if(!singleOccasionIdent){
      if(any(unlist(lapply(metricInvarianceFactors, function(mif) length(unique(diag(Phi)[mif])) != 1)))) warning('Factor variances differ across measurements. This is probably not intended. Either define the same variance for each factor across measurements or use single occasion identification (singleOccasionIdent = TRUE) or use referent scaling (useReferenceIndicator = TRUE).')  
      if(any(unlist(lapply(metricInvarianceFactors, function(mif) length(unique(Alpha[mif])) != 1)))) warning('Factor means differ across measurements. This is probably not intended. Either define the same mean for each factor across measurements or use single occasion identification (singleOccasionIdent = TRUE) or use referent scaling (useReferenceIndicator = TRUE).')  
    }
  }
  
  ### check Theta
  Theta  <- args[['Theta']]
  if(is.null(Theta) && any(c('residuals') %in% c(comparison, nullEffect))) stop('Tests involving residual invariance require specification of the residual variances (Theta).')
  if(!is.null(Theta)){
    checkSymmetricSquare(Theta, 'must be a symmetric square matrix.')
    if(ncol(Theta) != nrow(Lambda)) stop('The dimensions of Theta must equal the number of indicators.')
  }

    
  ### generate sigma
  
  generated <- semPower.genSigma(Lambda = Lambda, Phi = Phi, tau = tau, Alpha = Alpha, Theta = Theta, useReferenceIndicator = useReferenceIndicator)

  
  ### generate model strings

  # configural
  modelH1 <- modelH0 <- generated[['modelTrueCFA']]

  # loadings
  if('loadings' %in% c(comparison, nullEffect)){
    generatedMetric <- semPower.genSigma(Lambda = Lambda, Phi = Phi, tau = tau, Alpha = Alpha, useReferenceIndicator = useReferenceIndicator, metricInvariance = metricInvarianceFactors)
    if('loadings' %in% comparison) modelH1 <- generatedMetric[['modelTrueCFA']]
    if('loadings' %in% nullEffect) modelH0 <- generatedMetric[['modelTrueCFA']]
    # only fix first variance in single occasion ident
    if(singleOccasionIdent){
      for(i in (nFactors + 1):(nFactors*nWaves)){
        if('loadings' %in% comparison){
          modelH1 <- sub(paste0('f', i, ' ~~ ', '1*f', i), paste0('f', i, ' ~~ f', i), modelH1, fixed = TRUE)
        }
        if('loadings' %in% nullEffect){
          modelH0 <- sub(paste0('f', i, ' ~~ ', '1*f', i), paste0('f', i, ' ~~ f', i), modelH0, fixed = TRUE)
        }
      }
    }
    if(useReferenceIndicator){
      modelH0 <- append(modelH0, lapply(seq(ncol(Lambda)), function(f) paste0('f', f, ' ~~ f', f)))
      modelH1 <- append(modelH1, lapply(seq(ncol(Lambda)), function(f) paste0('f', f, ' ~~ f', f)))
    }
  }
  
  # intercepts
  if('intercepts' %in% c(comparison, nullEffect)){
    tok <- list()
    sIdx <- 1
    for(x in seq_along(metricInvarianceFactors)){
      ci <- lapply(metricInvarianceFactors[[x]], function(f) paste0('x', which(Lambda[, f] != 0)))
      eIdx <- (sIdx + length(ci[[1]]) - 1)
      lab <- paste0('i', sIdx:eIdx)
      if(useReferenceIndicator){
        tok <- append(tok, unlist(lapply(ci, function(x) paste0(x[1], ' ~ 0*1') )))
        tok <- append(tok, unlist(lapply(ci, function(x) paste0(x[-1], ' ~ ', lab[-1], '*1') )))
      }else{
        tok <- append(tok, unlist(lapply(ci, function(x) paste0(x, ' ~ ', lab, '*1') )))
      }
      sIdx <- sIdx + length(ci[[1]])
    }
    
    if('intercepts' %in% comparison) modelH1 <- append(modelH1, tok)
    if('intercepts' %in% nullEffect) modelH0 <- append(modelH0, tok)
    
    # estimate means when scalar invariance is set depending on scaling approach 
    if(useReferenceIndicator){
      ci <- paste0('f', seq(ncol(Lambda)), ' ~ 1')   # estimate all means
    }else if(singleOccasionIdent){
      ci <- paste0('f', seq(nFactors), ' ~ 0*1') 
      ci <- paste0('f', (nFactors + 1):ncol(Lambda), ' ~ 1')      # estimate all but first mean
    }else{
      ci <- paste0('f', seq(ncol(Lambda)), ' ~ 0*1')  # constrain all means
    }
    
    if('intercepts' %in% comparison) modelH1 <- append(modelH1, ci)
    if('intercepts' %in% nullEffect && !'means' %in% nullEffect) modelH0 <- append(modelH0, ci)
    # constrain means in h0 model when these are targeted by nulleffect 
    if('means' %in% nullEffect){
      # mean scaling only
      if('intercepts' %in% nullEffect && !useReferenceIndicator && !singleOccasionIdent){
        modelH0 <- append(modelH0, ci) 
      }else{
        # single occ and referent scaling
        tok <- list()
        for(x in seq_along(metricInvarianceFactors)){
          if(singleOccasionIdent){
            # first mean is zero, so all means are zero.
            tok <- append(tok, paste0('f', metricInvarianceFactors[[x]], ' ~ 0*1'))
          }else if(useReferenceIndicator){
            # reference indicator: restrict to equality
            tok <- append(tok, paste0('f', metricInvarianceFactors[[x]], ' ~ ', paste0('m', x), '*1'))
          }
        }
        modelH0 <- append(modelH0, tok)
      }
    }
  }
  
  # residuals
  if('residuals' %in% c(comparison, nullEffect)){
    tok <- list()
    sIdx <- 1
    for(x in seq_along(metricInvarianceFactors)){
      ci <- lapply(metricInvarianceFactors[[x]], function(f) paste0('x', which(Lambda[, f] != 0)))
      eIdx <- (sIdx + length(ci[[1]]) - 1)
      lab <- paste0('r', sIdx:eIdx)
      tok <- append(tok, unlist(lapply(ci, function(x) paste0(x, ' ~~ ', lab, '*', x) )))
      sIdx <- sIdx + length(ci[[1]])
    }
    if('residuals' %in% comparison) modelH1 <- append(modelH1, tok)
    if('residuals' %in% nullEffect) modelH0 <- append(modelH0, tok)
  }
  
  # lv var
  if('lv.variances' %in% c(comparison, nullEffect)){
    for(x in seq(metricInvarianceFactors)){
      for(ff in seq(length(metricInvarianceFactors[[x]]))){
        f <- metricInvarianceFactors[[x]][ff]
        pat <- paste0('f', f, ' ~~ f', f)
        if(useReferenceIndicator){
          rep <- paste0('f', f, ' ~~ v', x, '*f', f)
        }else if(singleOccasionIdent){
          rep <- paste0('f', f, ' ~~ 1*f', f)
        }  
        if('lv.variances' %in% comparison){
          for(i in seq(modelH1)){
            modelH1[[i]] <- sub(pat, rep, modelH1[[i]] , fixed = TRUE)
          }
        }
        if('lv.variances' %in% nullEffect){
          for(i in seq(modelH0)){
            modelH0[[i]] <- sub(pat, rep, modelH0[[i]] , fixed = TRUE)
          }
        }
      }
    }
  }
  
  # lv covar
  if(nFactors > 1 && 'lv.covariances' %in% c(comparison, nullEffect)){
    # constrain factor covariances to be equal across measurements, but do not constrain stabilities to equality 
    tok <- list()
    labs <- paste0('c', 1:(nFactors*(nFactors-1)/2)) 
    covFactors <- lapply(seq(length(metricInvarianceFactors[[1]])), function(x) unlist(lapply(metricInvarianceFactors, '[[', x)))
    for(f in seq(covFactors)){
      idx <- 1 
      x <- covFactors[[f]]
      for(i in min(x):max(x)){
        if((i + 1) <= max(x)){
          for(j in (i + 1):max(x)){
            tok <- append(tok, paste0('f', i, ' ~~ ',labs[idx],'*f', j))
            idx <- idx + 1
          }
        }
      }
    }
    if('lv.covariances' %in% comparison) modelH1 <- append(modelH1, tok)
    if('lv.covariances' %in% nullEffect) modelH0 <- append(modelH0, tok)
  }

  
  # add autocorrelated residuals
  if(autocorResiduals){
    autocorResidualsFactors <- metricInvarianceFactors  # same structure 
    tok <- list()
    for(x in seq_along(autocorResidualsFactors)){
      ci <- lapply(autocorResidualsFactors[[x]], function(f) paste0('x', which(Lambda[, f] != 0)))
      for(i in 1:(length(ci) - 1)){
        for(j in (i + 1) : length(ci)){
          tok <- append(tok, paste(ci[[i]], '~~', ci[[j]]))
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
  if(any(c('intercepts', 'means') %in% c(comparison, nullEffect))){
    mu <- generated[['mu']]
  }
  
  
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