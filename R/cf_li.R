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
#' For hypotheses regarding multiple group invariance, see [semPower.powerMI()]. For hypotheses regarding autoregressive models, see [semPower.powerAutoreg()]. For hypotheses in an ARMA model, see [semPower.powerARMA()].
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