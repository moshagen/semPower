#' semPower.powerBifactor
#'
#' Perform a power analysis for models including one or more bifactors to reject one of the following hypotheses: 
#' (a) a zero correlation between two factors, (b) the equality of two correlations between factors,
#' or (c) the equality of a correlation between two factors across two or more groups. 
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param bfLoadings a single vector or a list containing one or more vectors giving the loadings on each bifactor. For example, list(rep(.6, 10), rep(.6, 10)) defines two bifactors with 10 indicators each, loading by .6 each. Can be a list of lists for multiple group models.
#' @param bfWhichFactors a list containing one or more vectors defining which (specific) factors defined in the respective arguments in ... are part of the bifactor structure. See details.
#' @param Phi either a single number defining the correlation between exactly two factors or the factor correlation matrix. Must only contain the bifactor(s) and the covariate(s). Must be a list for multiple group models. Phi assumes the following order (bifactor_1, bifactor_2, ..., bifactor_j, covariate_1,  covariate_2, ...,  covariate_k). See details.
#' @param nullEffect defines the hypothesis of interest, must be one of `'cor = 0'` (the default) to test whether a correlation is zero, `'corX = corZ'` to test for the equality of correlations, and `'corA = corB'` to test for the equality of a correlation across groups. Define the correlations to be set to equality in `nullWhich` and the groups in `nullWhichGroups`. 
#' @param nullWhich vector of size 2 indicating which factor correlation in `Phi` is hypothesized to equal zero when `nullEffect = 'cor = 0'`, or to restrict to equality across groups when `nullEffect = 'corA = corB'`, or list of vectors defining which correlations to restrict to equality when `nullEffect = 'corX = corZ'`. Can also contain more than two correlations, e.g., `list(c(1, 2), c(1, 3), c(2, 3))` to set `Phi[1, 2] = Phi[1, 3] = Phi[2, 3]`. 
#' @param nullWhichGroups for `nullEffect = 'corA = corB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model concerning the specific factors and the covariate(s). See details.
#' @return a list. Use the `summary` method to obtain formatted results. Beyond the results of the power analysis and a number of effect size measures, the list contains the following components:
#' \item{`Sigma`}{the population covariance matrix. A list for multiple group models.}
#' \item{`mu`}{the population mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`SigmaHat`}{the H0 model implied covariance matrix. A list for multiple group models.}
#' \item{`muHat`}{the H0 model implied mean vector or `NULL` when no meanstructure is involved. A list for multiple group models.}
#' \item{`modelH0`}{`lavaan` H0 model string.}
#' \item{`modelH1`}{`lavaan` H1 model string or `NULL` when the comparison refers to the saturated model.}
#' \item{`simRes`}{detailed simulation results when a simulated power analysis (`simulatedPower = TRUE`) was performed.}
#' @details 
#' 
#' This function performs a power analysis to reject various hypotheses arising
#' in a model including a bifactor structure:
#' * `nullEffect = 'cor = 0'`: Tests the hypothesis that the correlation between a bifactor and another factor (which can also be a bifactor) is zero.
#' * `nullEffect = 'corX = corZ'`: Tests the hypothesis that two or more correlations involving one or more bifactors are equal to each other.
#' * `nullEffect = 'corA = corB'`: Tests the hypothesis that the correlation between the bifactor and another factor (which can also be a  bifactor) is equal in two or more groups (always assuming metric invariance).
#' 
#' A bifactor structure is defined by specifying the loadings on the general factor in `bfLoadings`, the comprised specific 
#' factors in `bfWhichFactors`, and the loadings on the specific factors in either `Lambda`, or `loadings`, 
#' or `nIndicator` and `loadM`. The latter arguments also include the loadings defining the 
#' covariate(s).
#' 
#' The correlations betwen the bifactor(s) and the covariate(s) are defined in `Phi`, which 
#' must omit the specific factors and only includes the bifactor(s) and the covariate(s) assuming 
#' the following order: (bifactor_1, bifactor_2, ..., bifactor_j, covariate_1,  covariate_2, ...,  covariate_k).
#' 
#' For example, the following defines a single bifactor with 10 indicators loading by .5 each. 
#' The bifactor structure involves 3 specific factors measured by 3 indicators each, each loading by
#' .3, .2, and .1 on the respective specific factor (in addition to the loadings on the bifactor). 
#' Furthermore, two covariate with 5 indicators each, all loading by .7, are defined. The correlation
#' between the covariates is .5, the one between the bifactor and the first and second covariate are
#' .3 and .2, respectively.
#' 
#' ```
#' bfLoadings <- rep(.5, 10)
#' bfWhichFactors <- c(1, 2, 3)
#' loadings <- list(
#'   rep(.3, 3),   # specific factor 1
#'   rep(.2, 3),   # specific factor 2
#'   rep(.1, 3),   # specific factor 3
#'   rep(.7, 5),   # covariate 1
#'   rep(.7, 5)    # covariate 2
#' )
#' Phi <- matrix(c(
#'   c(1, .3, .2),   # bifactor
#'   c(.3, 1, .5),   # covariate 1
#'   c(.2, .5, 1)   # covariate 2
#' ), ncol = 3, byrow = TRUE) 
#' ```
#'  
#' Beyond the arguments explicitly contained in the function call, additional arguments 
#' are required specifying the factor model and the requested type of power analysis.  
#' 
#' Additional arguments related to the **definition of the factor model** concerning the specific factors and the covariate(s). The loadings on the bifactor must be provided via `bfLoadings`.
#' * `Lambda`: The factor loading matrix (with the number of columns equaling the number of specific factors and covariates).
#' * `loadings`: Can be used instead of `Lambda`: Defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading on the second factor by .8, .6, .6, and .4.
#' * `nIndicator`: Can be used instead of `Lambda`: Used in conjunction with `loadM`. Defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' * `loadM`: Can be used instead of `Lambda`: Used in conjunction with `nIndicator`. Defines the loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the loadings of the first factor to equal .5 and those of the second factor do equal .6.
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
#' # get required N to detect a correlation of >= .3 between
#' # a single bifactor with 11 indicators all loadings by .6
#' # spanning the indicators of 3 specific factors
#' # with three indicators each, loading by .2, .15, and .25, respectively
#' # and a covariate measured by 4 indicators loading by .7 each,
#' # with a power of 95% on alpha = 5%
#' bfLoadings <- rep(.6, 11)
#' bfWhichFactors <- c(1, 2, 3)
#' loadings <- list(
#'   # specific factors
#'   rep(.2, 3),
#'   rep(.15, 3),
#'   rep(.25, 3),
#'   # covariate
#'   rep(.7, 4)
#' )
#' Phi <- .3    # bifactor - covariate
#' powerbifactor <- semPower.powerBifactor(type = 'a-priori',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         alpha = .05, beta = .05)
#' # show summary
#' summary(powerbifactor)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerbifactor$modelH1, sample.cov = powerbifactor$Sigma,
#'             sample.nobs = powerbifactor$requiredN, 
#'             sample.cov.rescale = FALSE)
#' lavaan::sem(powerbifactor$modelH0, sample.cov = powerbifactor$Sigma,
#'             sample.nobs = powerbifactor$requiredN, 
#'             sample.cov.rescale = FALSE)
#' 
#' # same as above, but determine power with N = 500 on alpha = .05
#' powerbifactor <- semPower.powerBifactor(type = 'post-hoc',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         alpha = .05, N = 500)
#' 
#' # same as above, but determine the critical chi-square with N = 500 so that alpha = beta
#' powerbifactor <- semPower.powerBifactor(type = 'compromise',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         abratio = 1, N = 500)
#' 
#' # same as above, but compare to the saturated model
#' # (rather than to the less restricted model)
#' powerbifactor <- semPower.powerBifactor(type = 'a-priori',
#'                                         comparison = 'saturated',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         alpha = .05, beta = .05)
#' 
#' # define two bifactors with 10 indicators each, where
#' # all loadings are .6 on the first and .5 on the second bifactor.
#' # the first bifactor spans the indicators of specific factors 1-3,
#' # the second bifactor spans the indicators of specific factors 4-6,
#' # all specific factors are measured by three indicators each,
#' # loadings are .2, .15, .25, .1, .15., and.2, respectively.
#' # define an additional  covariate measured by 4 indicators loading by .6 each.
#' # get required N to detect a correlation of >= .3 between the bifactors
#' # with a power of 95% on alpha = 5%
#' bfLoadings <- list(rep(.6, 10),
#'                    rep(.6, 10))
#' bfWhichFactors <- list(c(1, 2, 3),
#'                        c(4, 5, 6))
#' loadings <- list(
#'   # specific factors for bf1
#'   rep(.2, 3),
#'   rep(.15, 3),
#'   rep(.25, 3),
#'   # specific factors bf2
#'   rep(.1, 3),
#'   rep(.15, 3),
#'   rep(.2, 3),
#'   # covariate
#'   rep(.6, 4)
#' )
#' Phi <- diag(3)
#' Phi[1, 2] <- Phi[2, 1] <- .3    # bifactor1 - bifactor2
#' Phi[1, 3] <- Phi[3, 1] <- .5    # bifactor1 - covariate
#' Phi[2, 3] <- Phi[3, 2] <- .1    # bifactor2 - covariate
#' 
#' powerbifactor <- semPower.powerBifactor(type = 'a-priori',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         alpha = .05, beta = .05)
#' 
#' # same as above, but get required N to detect that
#' # the correlation between the first bifactor and the covariate (of r=.5) differs from
#' # the correlation between the second bifactor and the covariate (of r=.1)
#' powerbifactor <- semPower.powerBifactor(type = 'a-priori',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullEffect = 'corx = corz',
#'                                         nullWhich = list(c(1, 3), c(2, 3)),
#'                                         loadings = loadings,
#'                                         alpha = .05, beta = .05)
#' 
#' # multiple group example: get required N to detect that
#' # the correlation of a bifactor with 10 indicators
#' # spanning three specific factors with 3 indicators each
#' # to a covariate in group 1  (of r = .3)
#' # differs from the one in group 2 (of r = .1)
#' bfLoadings <- rep(.6, 10)
#' bfWhichFactors <- c(1, 2, 3)
#' loadings <- list(
#'   # specific factors
#'   rep(.2, 3),
#'   rep(.15, 3),
#'   rep(.25, 3),
#'   # covariate
#'   rep(.7, 4)
#' )
#' Phi1 <- Phi2 <- diag(2)
#' Phi1[1, 2] <- Phi1[2, 1] <- .3    # bifactor - covariate
#' Phi2[1, 2] <- Phi2[2, 1] <- .1    # bifactor - covariate
#' Phi <- list(Phi1, Phi2)
#' powerbifactor <- semPower.powerBifactor(type = 'a-priori',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullEffect = 'corA = corB',
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         alpha = .05, beta = .05, 
#'                                         N = list(1, 1))
#'                                         
#' # request a simulated post-hoc power analysis with 500 replications.
#' bfLoadings <- rep(.6, 11)
#' bfWhichFactors <- c(1, 2, 3)
#' loadings <- list(
#'   # specific factors
#'   rep(.2, 3),
#'   rep(.15, 3),
#'   rep(.1, 3),
#'   # covariate
#'   rep(.7, 5)
#' )
#' Phi <- .2  
#' set.seed(300121)
#' powerbifactor <- semPower.powerBifactor(type = 'post-hoc',
#'                                         bfLoadings = bfLoadings,
#'                                         bfWhichFactors = bfWhichFactors,
#'                                         Phi = Phi,
#'                                         nullWhich = c(1, 2),
#'                                         loadings = loadings,
#'                                         alpha = .05, N = 500, 
#'                                         simulatedPower = TRUE,
#'                                         simOptions = list(nReplications = 500)
#'                                         )
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerBifactor <- function(type, comparison = 'restricted', 
                                   bfLoadings = NULL,
                                   bfWhichFactors = NULL,
                                   Phi = NULL,
                                   nullEffect = 'cor = 0',
                                   nullWhich = NULL, 
                                   nullWhichGroups = NULL, 
                                   ...){
  
  # validate input
  checkEllipsis(...)
  args <- list(...)
  comparison <- checkComparisonModel(comparison)
  if(is.null(Phi)) stop('Phi must be defined')
  nullEffect <- checkNullEffect(nullEffect, c('cor=0', 'corx=corz', 'cora=corb'))
  
  # assume multigroup structure when phi is a list
  isMultigroup <- is.list(Phi)
  if(!is.list(Phi)) Phi <- list(Phi)
  nGroups <- length(Phi)
  Phi <- lapply(Phi, function(x){
    if(length(x) == 1){
      matrix(c(1, x, x, 1), ncol = 2)
    }else{
      x
    } 
  })
  if(length(unique(unlist(lapply(Phi, function(x) ncol(x))))) != 1 || length(unique(unlist(lapply(Phi, function(x) nrow(x))))) != 1) stop('Phi must be of same dimension in all groups.')
  if(nullEffect == 'cora=corb' && !isMultigroup) stop('Phi must be a list when nullEffect = "corA = corB".')
  if(nullEffect != 'cora=corb' && isMultigroup) stop('Multigroups are only supported for nullEffect = "corA = corB"')
  
  if(!is.list(bfLoadings)) bfLoadings <- list(list(bfLoadings)) # (groups)(bifactors)
  if(!is.list(bfLoadings[[1]])) bfLoadings <- list(bfLoadings) # (groups)(bifactors)
  if(!is.list(bfWhichFactors)) bfWhichFactors <- list(list(bfWhichFactors)) # (groups)(bifactors)
  if(!is.list(bfWhichFactors[[1]])) bfWhichFactors <- list(bfWhichFactors) # (groups)(bifactors)
  # assume same measurement model for all groups if not specified otherwise
  if(isMultigroup && length(bfLoadings) == 1) bfLoadings <- lapply(seq(nGroups), function(x) bfLoadings[[1]])
  if(isMultigroup && length(bfWhichFactors) == 1) bfWhichFactors <- lapply(seq(nGroups), function(x) bfWhichFactors[[1]])
  numBifactors <- length(bfLoadings[[1]])
  
  if(nullEffect == 'corx=corz' && !is.list(nullWhich)) stop('nullWhich must be a list when nullEffect = "corx=corz".')
  if(!is.list(nullWhich)) nullWhich <- list(nullWhich)
  if(any(unlist(lapply(nullWhich, function(x) length(x) != 2)))) stop('nullWhich may only contain vectors of size two.')
  if(any(unlist(lapply(nullWhich, function(x) x[1] == x[2])))) stop('elements in nullWhich may not refer to variances.')
  if(max(unlist(nullWhich)) > ncol(Phi[[1]])) stop('At least one element in nullWhich is an out of bounds index concerning Phi. Recall that Phi must only comprise bifactor(s) and covariate(s), but not specific factors.')
  
  ### create Lambda
  # sLambda only contains specific factors and covariate(s)
  sLambda <- args[['Lambda']]
  if(!is.null(sLambda)){
    if(!is.list(sLambda)) sLambda <- list(sLambda)
    # assume same measurement model in the multigroup case unless specified otherwise
    if(isMultigroup && length(sLambda) == 1) sLambda <- lapply(seq(nGroups), function(x) sLambda[[1]])
  }
  if(is.null(sLambda)){
    sLambda <- lapply(seq(nGroups), function(x){
      if(is.list(args[['loadings']][[1]]) || 
         any(unlist(lapply(args[names(args) %in% c('loadM', 'nIndicator', 'loadSD', 'loadMinMax')], is.list)))){
        genLambda(args[['loadings']][[x]], args[['nIndicator']][[x]],
                  args[['loadM']][[x]], args[['loadSD']][[x]], args[['loadMinMax']][[x]])
      }else{
        genLambda(args[['loadings']], args[['nIndicator']],
                  args[['loadM']], args[['loadSD']], args[['loadMinMax']])
      }
    })
  }
  if(any(unlist(lapply(seq(nGroups), function(g) lapply(seq(numBifactors), function(f){
    length(bfLoadings[[g]][[f]]) < sum(sLambda[[g]][ , bfWhichFactors[[g]][[f]]] != 0)
  }))))) stop('Bifactors must have at least the same number of indicators as the comprised by the respective specific factors. If you want an indicator only loading on a specific factor, assign the respective bifactor loading a value of zero.')
  
  # add bifactor structure to Lambda
  Lambda <- lapply(seq(nGroups), function(g){
    csLambda <- sLambda[[g]]
    cBfLoadings <- bfLoadings[[g]]
    cBfWhichFactors <- bfWhichFactors[[g]]
    extraIndicators <- unlist(lapply(seq(numBifactors), function(x){
      length(cBfLoadings[[x]]) - sum(csLambda[, cBfWhichFactors[[x]]] != 0)
    }))
    cLambda <- matrix(0, nrow = (sum(extraIndicators) + nrow(csLambda)), ncol = (numBifactors + ncol(csLambda)))  
    # specific factors + covariate(s)
    cLambda[(sum(extraIndicators) + 1) : nrow(cLambda), (numBifactors + 1) : ncol(cLambda)] <- csLambda
    extraLoadings <- lapply(seq_along(cBfLoadings), function(x) cBfLoadings[[x]][1:extraIndicators[x]])
    for(i in 1:numBifactors){
      # bifactor extra indicators
      if(i == 1) idx <- 1 : extraIndicators[i] else idx <- (1 + extraIndicators[i - 1]) : (extraIndicators[i - 1] + extraIndicators[i])
      cLambda[idx, i] <- extraLoadings[[i]]
      # bifactor specific indicators 
      idx <- length(unlist(extraLoadings)) + unlist(lapply(cBfWhichFactors[[i]], function(x) which(csLambda[, x] != 0)))
      cLambda[idx, i] <- cBfLoadings[[i]][(1 + extraIndicators[i]) : length(cBfLoadings[[i]])]
    }
    cLambda
  })
  
  # warn if loadings dont satisfy metric invariance
  if(isMultigroup){
    lambdas <- do.call(cbind, lapply(Lambda, c))
    if(any(apply(lambdas, 1, function(x) length(unique(x)) != 1))) warning('At least one loading differs across groups, violating metric invariance. Verify that this is intended.')
  }
  
  ### create Phi
  sPhi <- Phi
  nf <- 1:ncol(Lambda[[1]])
  posCov <- nf[!nf %in% c(1 : numBifactors, (numBifactors + unlist(bfWhichFactors[[1]])))]
  numCov <- length(posCov)
  Phi <- lapply(seq(nGroups), function(g){
    csPhi <- sPhi[[g]]
    if(ncol(csPhi) != (numCov + numBifactors)) stop('Incorrect dimensions for Phi. Remember that Phi may only refer to the bifactor(s) and the covariate(s), but must omit all specific factors.')
    cPhi <- diag(ncol(Lambda[[g]]))
    # bf-bf
    cPhi[1:numBifactors, 1:numBifactors] <-  csPhi[1:numBifactors, 1:numBifactors]
    # bf-cov
    if(numCov > 0){
      cPhi[posCov, 1:numBifactors] <- csPhi[(numBifactors + 1) : ncol(csPhi), 1 : numBifactors]
      cPhi[1:numBifactors, posCov] <- csPhi[1 : numBifactors, (numBifactors + 1) : ncol(csPhi)]
      # cov-cov
      cPhi[posCov, posCov] <- csPhi[(numBifactors + 1) : ncol(csPhi), (numBifactors + 1) : ncol(csPhi)]
    }
    cPhi
  })
  
  ### generate sigma 
  if(isMultigroup){
    generated <- semPower.genSigma(Lambda = Lambda, Phi = Phi)
    model <- generated[[1]][['modelTrueCFA']]
  }else{
    generated <- semPower.genSigma(Lambda = Lambda[[1]], Phi = Phi[[1]])
    model <- generated[['modelTrueCFA']] 
  }
  
  ### create model strings 
  # add bifactor orthogonalization constrains
  tok <- lapply(seq(numBifactors), function(x){
    ctok <- list()
    csf <- paste0('f', c(x, numBifactors + bfWhichFactors[[1]][[x]]))
    for(i in 1 : (length(csf) - 1)){
      for(j in (i + 1): length(csf)){
        ctok <- append(ctok, paste0(csf[i], ' ~~ 0*', csf[j]))
      }
    }
    ctok
  })
  model <- paste(c(model, unlist(tok)), collapse = '\n')
  
  ### H0 model
  # adapt nullWhich to account for specific factors
  nullWhich <- lapply(nullWhich, function(x){
    if(x[[1]] > numBifactors) x[[1]] <- posCov[(x[[1]] - numBifactors)]
    if(x[[2]] > numBifactors) x[[2]] <- posCov[(x[[2]] - numBifactors)]
    x
  })  
  if(nullEffect == 'cor=0'){
    modelH0 <- paste(c(model,
                       paste0('f', nullWhich[[1]], collapse = ' ~~ 0*')),
                     collapse = '\n')
  }else if(nullEffect == 'corx=corz'){
    labs <- list()
    tok <- ''
    for(i in seq_along(nullWhich)){
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
    modelH0 <- paste(c(model, tok), collapse = '\n')
  }else if(nullEffect == 'cora=corb'){
    if(is.null(nullWhichGroups)) nullWhichGroups <- seq(nGroups)
    lab <- paste0('ff', seq(nGroups))
    lab[nullWhichGroups] <- 'pf1'
    lab <- paste0('c(', paste(lab, collapse = ','), ')*')
    modelH0 <- paste(c(model,
                       paste0('f', nullWhich[[1]], collapse = paste0(' ~~ ', lab))),
                     collapse = '\n')
  }else{
    stop('nullEffect not defined')
  }
  
  # we always enforce invariance constraints in the multigroup case
  if(isMultigroup){
    args[['lavOptions']] <- append(args[['lavOptions']], list(group.equal = c('loadings', 'lv.variances')))
  } 
  
  modelH1 <- NULL
  fitH1model <- FALSE
  if(comparison == 'restricted'){
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