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
#' # the intercepts are .1, .2, .3, .4, .5, .6 in the second group
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
#'                               seq(.1, .6, .1) 
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
#'                               seq(.1, .6, .1) 
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