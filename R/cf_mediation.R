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
#' @param estimateDirectEffects Whether to estimate all direct effects (`TRUE`, the default). If `FALSE`, only direct effects unequal zero in the population are estimated.
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
#' * More complex mediation structures can be defined by providing the `Beta` matrix along with `indirect` specifying which paths define the indirect effect. See examples below. To specify residual correlations, use [semPower.powerPath()] in conjunction with [semPower.powerLav()].
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
#' Note that in case of a simple mediation model involving three variables, the order of the factors is X, M, Y, i. e., the first factor is treated as X, the second as M, and the third as Y. In case of a more complex mediation defined via the `Beta` matrix, the order of factors matches the order of `Beta`. 
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
#' # simple case of X -> M -> Y mediation in the form of
#' # X -- .30 --> M -- .40 --> Y
#' # X --------- .25 --------> Y
#' # determine the required N to detect the indirect effect of >= .12 (= .3 * .4) 
#' # with a power of 95% on alpha = 5%, where   
#' # X is measured by 3 indicators loading by .5 each, 
#' # M is measured by 5 indicators loading by .6 each, 
#' # Y is measured by 4 indicators loading by .7 each.
#' powerMed <- semPower.powerMediation(type = 'a-priori',
#'                                     alpha = .05, beta = .05,
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(3, 5, 4),
#'                                     loadM = c(.5, .6, .7)
#'                                     )
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
#'                                     alpha = .05, N = 500,
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(3, 5, 4),
#'                                     loadM = c(.5, .6, .7))
#' 
#' 
#' # same as above, but determine the critical chi-square with N = 500 so that alpha = beta
#' powerMed <- semPower.powerMediation(type = 'compromise',
#'                                     abratio = 1, N = 500
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(3, 5, 4),
#'                                     loadM = c(.5, .6, .7)
#'                                     )
#' 
#' 
#' # same as above, but compare to the saturated model
#' # (rather than to the less restricted model)
#' powerMed <- semPower.powerMediation(type = 'a-priori',
#'                                     alpha = .05, beta = .05,
#'                                     comparison = 'saturated',
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(3, 5, 4),
#'                                     loadM = c(.5, .6, .7)
#'                                     )
#'                                     
#' 
#' # same as above, but assuming observed variables only (Lambda = diag(3))
#' powerMed <- semPower.powerMediation(type = 'a-priori',
#'                                     alpha = .05, beta = .05,
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     Lambda = diag(3)
#'                                     )
#'                                     
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
#'                                     alpha = .05, beta = .05
#'                                     Beta = Beta, 
#'                                     indirect = list(c(2, 1), c(3, 2)),
#'                                     nIndicator = c(3, 5, 4),
#'                                     loadM = c(.5, .6, .7)
#'                                     )
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
#'                                     alpha = .05, beta = .05,
#'                                     Beta = Beta, 
#'                                     indirect = list(c(2, 1), 
#'                                                     c(3, 2), 
#'                                                     c(4, 3)),
#'                                     loadings = loadings)
#'                                     
#' 
#' # Determine required N to detect that the indirect effect 
#' # in group 1 (of .2 * .3 = .09) differs from the indirect effect 
#' # in group 2 (of .3 * .5 = .15).
#' # The direct effect of X on Y is .25 in both groups.  
#' # The model is based on observed variables only (Lambda = diag(3))
#' # Both groups are sized equally (N = list(1, 1)).
#' powerMed <- semPower.powerMediation(type = 'a-priori',
#'                                     alpha = .05, beta = .05, N = list(1, 1),
#'                                     nullEffect = 'indA = indB',
#'                                     bYX = list(.25, .25), 
#'                                     bMX = list(.2, .3), 
#'                                     bYM = list(.3, .5),
#'                                     Lambda = diag(3)
#'                                     )
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
#'                                     alpha = .05, beta = .05, N = list(1, 1),
#'                                     nullEffect = 'indA = indB',
#'                                     Beta = list(Beta1, Beta2), 
#'                                     indirect = list(c(2, 1), c(3, 2)),
#'                                     Lambda = diag(3)
#'                                     )
#' 
#' # request a simulated post-hoc power analysis with 500 replications.
#' set.seed(300121)
#' powerMed <- semPower.powerMediation(type = 'post-hoc',
#'                                     alpha = .05, N = 500,
#'                                     bYX = .25, bMX = .3, bYM = .4,
#'                                     nIndicator = c(3, 5, 4),
#'                                     loadM = c(.5, .6, .7),
#'                                     simulatedPower = TRUE, 
#'                                     simOptions = list(nReplications = 500))
#'}
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerMediation <- function(type, comparison = 'restricted',
                                    bYX = NULL, bMX = NULL, bYM = NULL,
                                    Beta = NULL, indirect = NULL, 
                                    estimateDirectEffects = TRUE,
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
  if(isObserved) model <- '' # dummy latents don't work with non-linear constraints
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
  # ensure we did not miss direct effects
  if(estimateDirectEffects){
    # possible direct effects
    idx <- cbind(unlist(lapply(2:nrow(B[[1]]), function(x) x:nrow(B[[1]]))), rep(1:(nrow(B[[1]])-1), (nrow(B[[1]])-1):1))
    # skip indirect and beta!=0
    de <- unique(rbind(do.call(rbind, indirect), idx))
    xx <- lapply(B, function(x) x[de] == 0)
    de <- matrix(de[apply(do.call(cbind, xx), 1, function(x) any(x)), ], ncol=2)
    # add to modelstring
    if(nrow(de) > 0){
      for(i in 1:nrow(de)){
        f <- de[i, 1]
        fidx <- de[i, 2]
        clab <- paste0('pf', paste0(formatC(f, width = 2, flag = 0), formatC(fidx, width = 2, flag = 0)))
        if(isMultigroup) clab <- unlist(lapply(clab, function(x) paste0('c(', paste0(x, 'g', seq_along(B), collapse = ', '), ')')))
        tok <- paste0('f', f, ' ~ ', paste(clab, paste0('*f', fidx), sep = '', collapse = ' + '))
        model <- paste(model, tok, sep='\n')
      }
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
  
  # for observed only models, replace factor labels (f) by observed variables as found in sigma
  if(isObserved){
    if(isMultigroup) tSigma <- generated[[1]][['Sigma']] else tSigma <- generated[['Sigma']] 
    ff <- paste0('f', 1:ncol(tSigma))
    for(i in 1:ncol(tSigma)){
      model <- gsub(ff[i], colnames(tSigma)[i], model)  # later reused for modelH1
      modelH0 <- gsub(ff[i], colnames(tSigma)[i], modelH0)
    }
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