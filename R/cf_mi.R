#' semPower.powerMI
#'
#' Convenience function for performing power analyses for hypothesis arising 
#' in multigroup measurement invariance models concerning a specific level of invariance.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, either `'saturated'` or one of `'configural'`, `'metric'`, `'scalar'`, `'variances'`, `'covariances'`, or a vector of restrictions in `lavaan` format (with `'none'` for no restrictions). See details.
#' @param nullEffect defines the hypothesis (i.e., level of invariance) of interest. One of `'metric'`, `'scalar'`, `'residual'`, `'variances'`, `'covariances'`, `'means'` or a vector of restrictions in `lavaan` format. See details.   
#' @param useReferenceIndicator whether to identify the factor variances and factor means using a referent indicator (`TRUE`) or using a standardization approach setting the means to 0 and the variances to 1 (`FALSE`, the default). See details.
#' @param singleGroupIdent when `useReferenceIndicator = FALSE`, whether to maintain the identification constraints only for a single group (`TRUE`) or whether to maintain all scaling constraints (`FALSE`). See details.
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
#' \item `'variances'`: all loadings, (indicator-)intercepts, (indicator-)residuals, and latent variances are restricted to equality.
#' \item `'covariances'`: all loadings, (indicator-)intercepts, (indicator-)residuals, latent variances, and latent covariances are restricted to equality.
#' \item `'means'`: all loadings, (indicator-)intercepts, (indicator-)residuals, latent variances, latent covariances, and latent means are restricted to equality.
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
#' 
#' By default, the factors are identified by setting their means to 0 and their variances to 1 (`useReferenceIndicator = FALSE`), so that all loadings and indicator intercepts are freely estimated. 
#' In models involving constraints on loadings or intercepts, a single group identification approach is used by default (`singleGroupIdent = TRUE`), so that the factor variances and means are constrained in only a single group, but are freely estimated in all remaining groups. This mimics the `lavaan` approach to invariance testing.
#' If setting `singleGroupIdent = FALSE`, the identification constraints are maintained in models assuming metric or scalar invariance, so that the factor means and variances are constrained to 0 and 1, respectively, in all groups. Concerning metric invariance, this was the default behavior in `semPower` versions prior 2.1.4. 
#' A referent scaling approach can be be defined using `useReferenceIndicator = TRUE`, so that the loading and the intercept of the first indicator of each factor are constrained in each group, and the factor variances and means are freely estimated. 
#' Reference scaling is always maintained in all invariance models, so `singleGroupIdent` has no effect (i.e. is always `FALSE`)
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
#' # with a power of 80% on alpha = 5% 
#' # assuming equally sized groups (N = list(1, 1)) 
#' # for a factor model involving a single factor which 
#' # is measured by 4 indicators (in both groups)
#' # loading by .9, .8, .7, .6 in the first group and 
#' # loading by .7, .7, .8, .7 in the second group.
#' powerMI <- semPower.powerMI(type = 'a-priori',
#'                             alpha = .05, power = .80, N = list(1, 1),
#'                             comparison = 'configural', 
#'                             nullEffect = 'metric',
#'                             loadings = list(
#'                               list(c(.9, .8, .7, .6)),
#'                               list(c(.7, .7, .8, .7)))
#' )
#' 
#' # show summary
#' summary(powerMI)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerMI$modelH1, sample.cov = list(powerMI$Sigma[[1]], powerMI$Sigma[[2]]),
#'             sample.nobs = as.list(powerMI$requiredN.g), sample.cov.rescale = FALSE)
#' 
#' # same as above, but determine power with N = 500 in each group on alpha = .05
#' powerMI <- semPower.powerMI(type = 'post-hoc',
#'                             alpha = .05, N = list(500, 500),
#'                             comparison = 'configural', 
#'                             nullEffect = 'metric',
#'                             loadings = list(
#'                               list(c(.9, .8, .7, .6)),
#'                               list(c(.7, .7, .8, .7)))
#'                             )
#' 
#' # same as above, but determine the critical chi-square with N = 500 in each 
#' # group so that alpha = beta
#' powerMI <- semPower.powerMI(type = 'compromise',
#'                             abratio = 1, N = list(500, 500),
#'                             comparison = 'configural', 
#'                             nullEffect = 'metric',
#'                             loadings = list(
#'                               list(c(.9, .8, .7, .6)),
#'                               list(c(.7, .7, .8, .7)))
#'                             )
#' 
#' # same as above, but compare to the saturated model
#' # (rather than to the configural invariance model)
#' powerMI <- semPower.powerMI(type = 'a-priori',
#'                             alpha = .05, beta = .05, N = list(1, 1),
#'                             comparison = 'saturated', 
#'                             nullEffect = 'metric',
#'                             loadings = list(
#'                               list(c(.9, .8, .7, .6)),
#'                               list(c(.7, .7, .8, .7)))
#'                             )
#' 
#' # same as above, but make first group twice as large as the second group 
#' powerMI <- semPower.powerMI(type = 'a-priori',
#'                             alpha = .05, beta = .05, N = list(2, 1),
#'                             comparison = 'saturated', 
#'                             nullEffect = 'metric',
#'                             loadings = list(
#'                               list(c(.9, .8, .7, .6)),
#'                               list(c(.7, .7, .8, .7)))
#'                             )
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
#'                             alpha = .05, beta = .05, N = list(1, 1)
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
#'                             ))
#' 
#' # same as above, but use lavaan group.equal strings 
#' powerMI <- semPower.powerMI(type = 'a-priori',
#'                             alpha = .05, beta = .05, N = list(1, 1)
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
#'                             ))
#' 
#' # same as above, but
#' # obtain the required N to reject the hypothesis of equal latent means
#' # in comparison to the scalar invariance model;
#' # all intercepts are zero in both groups, 
#' # in the first group, the latent means equal 0.0, 
#' # in the second group, the latent mean of the factors are 0.0 and 0.5
#' powerMI <- semPower.powerMI(type = 'a-priori',
#'                             alpha = .05, beta = .05, N = list(1, 1),
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
#'                             ))                             
#' 
#' # request a simulated post-hoc power analysis with 500 replications
#' # to reject the hypothesis of metric invariance.
#' set.seed(300121)
#' powerMI <- semPower.powerMI(type = 'post-hoc',
#'                             alpha = .05, N = list(500, 500), 
#'                             comparison = 'configural', 
#'                             nullEffect = 'metric',
#'                             loadings = list(
#'                               list(c(.9, .8, .7, .6)),
#'                               list(c(.7, .7, .8, .7))),
#'                             simulatedPower = TRUE, 
#'                             simOptions = list(nReplications = 500))
#'                              
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerMI <- function(type, 
                             comparison = NULL,
                             nullEffect = NULL,
                             useReferenceIndicator = FALSE,
                             singleGroupIdent = !useReferenceIndicator,
                             ...){
  
  args <- list(...)
  
  # validate input
  checkEllipsis(...)
  lavGroupStrings <- c('loadings', 'intercepts', 'residuals', 'residual.covariances', 'lv.variances', 'lv.covariances','regressions')
  useLavOptions <- any(grepl(paste(lavGroupStrings, collapse = '|'), comparison)) || any(grepl(paste(lavGroupStrings, collapse = '|'), nullEffect))
  # we only check typos etc when not using lavstrings
  if(!useLavOptions){
    comparison <- checkNullEffect(comparison, c('saturated', 'configural', 'metric', 'scalar', 'variances', 'covariances'))
    nullEffect <- checkNullEffect(nullEffect, c('metric', 'scalar', 'residual', 'variances', 'covariances', 'means'))
    if(which(c('saturated', 'configural', 'metric', 'scalar', 'variances', 'covariances') %in% comparison) >= 
       (2 + which(c('metric', 'scalar', 'residuals', 'variances', 'covariances', 'means') %in% nullEffect))) stop('Model defined in nullEffect is not nested in comparison model.')
  }else{
    if(!useReferenceIndicator && !singleGroupIdent && ('lv.variances' %in% comparison || 'lv.variances' %in% nullEffect)) stop('Variance scaling without single group identification is used, so invariance of latent variances is always met. Either use single group identification or referent scaling.')
    if(!useReferenceIndicator && !singleGroupIdent && ('means' %in% nullEffect)) stop('Mean scaling without single group identification is used, so invariance of latent means is always met. Either use single group identification or referent scaling.')
    if(!any(c('saturated', 'none') %in% comparison) && !all(comparison %in% nullEffect)) stop('Comparison model is not nested in hypothesized model; all restrictions in comparison must also be present in nullEffect.')
  }
  if(useReferenceIndicator && singleGroupIdent) stop('Identification constraints on referent indicators (useReferenceIndicator = TRUE) must be maintained, so singleGroupIdent must be FALSE.')
  
  ### generate sigmas
  generated <- semPower.genSigma(..., useReferenceIndicator = useReferenceIndicator)   
  
  # more input validations
  if(!is.list(generated[[1]])) stop('Loadings, Phi, Beta, etc. must be provided as a list for each group.')
  if(is.null(generated[[1]][['mu']])){
    inv <- FALSE
    if(useLavOptions){
      inv <- any(grepl('intercepts|means', comparison)) || any(grepl('intercepts|means', nullEffect))
    }else{
      inv <- any(c('scalar', 'residuals', 'variances', 'covariances', 'means') %in% nullEffect)
    }
    if(inv) stop('The models imply a meanstructure, so tau and/or Alpha need to be defined.')
  }
  # throw warning if identification constraints imply misfit 
  if(useReferenceIndicator){
    lambda <- lapply(generated, '[[', 'Lambda')
    iIdx <- sapply(seq(ncol(lambda[[1]])), function(x) min(which(lambda[[1]][ , x] != 0)))
    rIdx <- unlist(lapply(seq(iIdx), function(x) iIdx[x] + (x-1) * nrow(lambda[[1]]))) # adapt index
    if(any(unlist(lapply(seq(rIdx), function(x) length(unique(sapply(lambda, '[[', rIdx[x]))) != 1)))) warning('Loadings of the referent indicator differ across groups. This is probably not intended. Either define the same loading for the first indicator of each factor or use variance scaling (useReferenceIndicator = FALSE).')
    if(!is.null(generated[[1]][['tau']])){
      tau <- do.call(cbind, lapply(generated, '[[', 'tau'))
      if(any(unlist(lapply(seq(iIdx), function(x) length(unique(tau[iIdx[x], ])) != 1))))  warning('Intercepts of the referent indicator differ across groups. This is probably not intended. Either define the same intercept for the first indicator of each factor or use mean scaling (useReferenceIndicator = FALSE).')
    }
  }else{
    if(!singleGroupIdent){
      phi <- do.call(cbind, lapply(lapply(generated, '[[', 'Phi'), diag))
      if(any(apply(phi, 1, function(x) length(unique(x)) != 1))) warning('Factor variances differ across groups. This is probably not intended. Either define the same variance for each factor across groups or use referent scaling (useReferenceIndicator = TRUE).')  
      if(!is.null(generated[[1]][['Alpha']])){
        alpha <- do.call(cbind, lapply(generated, '[[', 'Alpha'))
        if(any(apply(alpha, 1, function(x) length(unique(x)) != 1))) warning('Factor means differ across groups. This is probably not intended. Either define the same mean for each factor across groups or use referent scaling (useReferenceIndicator = TRUE).')
      }
    }
  }
  
  # set proper lavOptions
  lavOptionsH1 <- list()
  if(!useLavOptions){
    lavOptionsH0 <- list(group.equal = switch(nullEffect,
                                              'metric' = c('loadings'),
                                              'scalar' = c('loadings', 'intercepts'),
                                              'residual' = c('loadings', 'intercepts', 'residuals'),
                                              'variances' = c('loadings', 'intercepts', 'residuals', 'lv.variances'),
                                              'covariances' = c('loadings', 'intercepts', 'residuals', 'lv.variances', 'lv.covariances'),
                                              'means' = c('loadings', 'intercepts', 'residuals', 'lv.variances', 'lv.covariances', 'means')
                                              
    ))
    if(comparison %in% c('metric', 'scalar', 'residuals', 'variances', 'covariances')){
      lavOptionsH1 <- list(group.equal = switch(comparison,
                                                'metric' = c('loadings'),
                                                'scalar' = c('loadings', 'intercepts'),
                                                'residual' = c('loadings', 'intercepts', 'residuals'),
                                                'variances' = c('loadings', 'intercepts', 'residuals', 'lv.variances'),
                                                'covariances' = c('loadings', 'intercepts', 'residuals', 'lv.variances', 'lv.covariances')
      ))
    }
  }else{
    lavOptionsH0 <- list(group.equal = nullEffect)
    if(!any(c('saturated', 'none') %in% comparison)){
      lavOptionsH1 <- list(group.equal = comparison)
    }
  }

  
  # define model strings
  
  # for referent scaling, models are the same, the only difference pertains to lavOptions
  modelH0 <- modelH1 <- generated[[1]][['modelTrueCFA']]

  # add constraints on referent intercept
  if(useReferenceIndicator){
    if("intercepts" %in% lavOptionsH0[['group.equal']]){
      modelH0 <- paste(modelH0, paste0('x', iIdx, ' ~ 0*1', collapse = '\n'), sep = '\n')
      modelH0 <- paste(modelH0, paste0('f', seq(ncol(generated[[1]][['Lambda']])), ' ~ 1', collapse = '\n'), sep = '\n')
    }    
  }
    
  # variance scaling
  if(!useReferenceIndicator){
    # for single group identification, freely estimate latent variances/means in all but one group when loadings/intercepts are equal
    nGroups <- length(generated)
    nFactors <- ncol(generated[[1]][['Lambda']])
    if(singleGroupIdent){
      if("loadings" %in% lavOptionsH0[['group.equal']]){
        replacement <- paste0('c(1,', paste(rep('NA', (nGroups-1)), collapse = ','), ')*')
        for(i in 1:nFactors){
          modelH0 <- sub(paste0('f', i, ' ~~ ', '1*f', i), paste0('f', i, ' ~~ ', replacement, 'f', i), modelH0, fixed = TRUE)
        }
        if("loadings" %in% lavOptionsH1[['group.equal']]){
          modelH1 <- sub(paste0('f', i, ' ~~ ', '1*f', i), paste0('f', i, ' ~~ ', replacement, 'f', i), modelH1, fixed = TRUE)
        }
      }
      if("intercepts" %in% lavOptionsH0[['group.equal']]){
        tok <- paste0('c(0,', paste(rep('NA', (nGroups-1)), collapse = ','), ')*')
        modelH0 <- paste(modelH0, paste0(unlist(lapply(1:nFactors, function(x) paste0('f', x, ' ~ ', tok, '1'))), collapse = '\n'), sep = '\n')
        if("intercepts" %in% lavOptionsH1[['group.equal']]){
          modelH1 <- paste(modelH1, paste0(unlist(lapply(1:nFactors, function(x) paste0('f', x, ' ~ ', tok, '1'))), collapse = '\n'), sep = '\n')
        }
      }
      
    # otherwise retain scaling constraints; also for factor means (!)
    }else{
      if("intercepts" %in% lavOptionsH0[['group.equal']]){
        tok <- paste0('c(', paste(rep('0', nGroups), collapse = ','), ')*')
        modelH0 <- paste(modelH0, paste0(unlist(lapply(1:nFactors, function(x) paste0('f', x, ' ~ ', tok, '1'))), collapse = '\n'), sep = '\n')
        if("intercepts" %in% lavOptionsH1[['group.equal']]){
          modelH1 <- paste(modelH1, paste0(unlist(lapply(1:nFactors, function(x) paste0('f', x, ' ~ ', tok, '1'))), collapse = '\n'), sep = '\n')
        }
      }
    }
  }

  
  
  if('saturated' %in% comparison) modelH1 <- NULL
  
  Sigma <- lapply(generated, '[[', 'Sigma')
  mu <- NULL
  if(!useLavOptions){
    if(any(c('scalar', 'residuals', 'variances', 'covariances', 'means') %in% nullEffect))
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