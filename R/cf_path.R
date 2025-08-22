#' semPower.powerPath
#'
#' Convenience function for performing power analyses for hypothesis arising 
#' in a generic path model.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param Beta matrix of regression slopes between latent variables (all-Y notation). A list for multiple group models. Exogenous variables must occupy the first rows in `Beta` when `standardized = TRUE`. See details. 
#' @param Psi variance-covariance matrix of latent (residual) factors. If `standardized = TRUE`, the diagonal is ignored and all off-diagonal elements are treated as correlations. If `NULL`, an identity matrix is assumed. A list for multiple group models. See details.
#' @param nullEffect defines the hypothesis of interest, must be one of `'beta = 0'` (the default) to test whether a regression slope is zero, `'betaX = betaZ'` to test for the equality of slopes, and `'betaX = betaZ'` to test for the equality of a slope across groups. Define the slopes to be set to equality in `nullWhich` and the groups in `nullWhichGroups`. 
#' @param nullWhich vector of size 2 indicating which slope in `Beta` is hypothesized to equal zero when `nullEffect = 'beta = 0'`, or to restrict to equality across groups when `nullEffect = 'betaA = betaB'`, or list of vectors defining which correlations to restrict to equality when `nullEffect = 'betaX = betaZ'`. Can also contain more than two slopes, e.g., `list(c(2, 1), c(3, 1), c(3, 2))` to set `Beta[2, 1] = Beta[3, 1] = Beta[3, 2]`.
#' @param nullWhichGroups for `nullEffect = 'betaA = betaB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param standardized whether all parameters should be standardized (`TRUE`, the default). If `FALSE`, all regression relations are unstandardized.
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
#' This function performs a power analysis to reject a hypothesis arising
#' in a generic structural equation model specifying regression relations between the factors via the Beta matrix:  
#' * `nullEffect = 'beta = 0'`: Tests the hypothesis that a slope is zero. 
#' * `nullEffect = 'betaX = betaZ'`: Tests the hypothesis that two or more slopes are equal to each other.
#' * `nullEffect = 'betaA = betaB'`: Tests the hypothesis that a slope is equal in two or more groups (always assuming metric invariance).
#' 
#' This function provides a generic way to perform power analyses (as compared to other functions covering special cases in a more accessible manner).
#' 
#' A specific hypothesis is defined by setting `nullEffect` to define the hypothesis type, 
#' `nullWhich` to define the slope(s) that are targeted, and by providing 
#' the `Beta` (and optionally the `Psi`) matrix to define the population structure.
#'  
#' To understand the structure of `Beta` and `Psi`, consider the general structural equation model, 
#' \deqn{\Sigma = \Lambda (I - B)^{-1} \Psi [(I - B)^{-1}]'  \Lambda' + \Theta } 
#' where \eqn{B} is the \eqn{m \cdot m} matrix containing the regression slopes and \eqn{\Psi} is the (residual) variance-covariance matrix of the \eqn{m} factors. 
#' 
#' As an example, suppose there are four factors (X1, X2, X3, X4), and Beta is defined as follows:
#' 
#' \eqn{
#' \begin{array}{lrrrr} 
#'     & X_1 & X_2 & X_3 & X_4\\ 
#' X_1 & 0.0 & 0.0 & 0.0 & 0.0 \\ 
#' X_2 & 0.0 & 0.0 & 0.0 & 0.0  \\ 
#' X_3 & 0.2 & 0.3 & 0.0 & 0.0  \\ 
#' X_4 & 0.3 & 0.5 & 0.0 & 0.0  \\ 
#' \end{array}
#' }
#' 
#' Each row specifies how a particular factor is predicted by the available factors, 
#' so the above implies the following regression relations:
#' 
#' \eqn{
#' X_1 = 0.0 \cdot X_1 +  0.0 \cdot X_2 + 0.0 \cdot X_3 + 0.0 \cdot X_4 \\
#' X_2 = 0.0 \cdot X_1 +  0.0 \cdot X_2 + 0.0 \cdot X_3 + 0.0 \cdot X_4 \\
#' X_3 = 0.2 \cdot X_1 +  0.3 \cdot X_2 + 0.0 \cdot X_3 + 0.0 \cdot X_4 \\
#' X_4 = 0.3 \cdot X_1 +  0.5 \cdot X_2 + 0.0 \cdot X_3 + 0.0 \cdot X_4 
#' }
#' 
#' which simplifies to
#' 
#' \eqn{
#' X_3 = 0.2 \cdot X_1 + 0.3 \cdot X_2 \\
#' X_4 = 0.3 \cdot X_1 + 0.5 \cdot X_2 
#' }
#' 
#' Further suppose that Psi is
#' 
#' \eqn{
#' \begin{array}{lrrrr} 
#'     & X_1 & X_2 & X_3 & X_4\\ 
#' X_1 & 1.0 & 0.3 & 0.0 & 0.0 \\ 
#' X_2 & 0.3 & 1.0 & 0.0 & 0.0 \\ 
#' X_3 & 0.0 & 0.0 & 1.0 & 0.2 \\ 
#' X_4 & 0.0 & 0.0 & 0.2 & 1.0 \\ 
#' \end{array}
#' }
#' 
#' which implies a correlation between X1 and X2 of .3 and a residual correlation
#' between X3 and X4 of .2. 
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
#' So either `Lambda`, or `loadings`, or `nIndicator` and `loadM` always need to be defined. 
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
#' # set up pathmodel in the form of
#' # f2 = .2*f1
#' # f3 = .3*f2
#' # f4 = .1*f1 + .4*f3
#' # obtain the required N to detect that the 
#' # slope f1 -> f4 is >= .10 
#' # with a power of 95% on alpha = 5%
#' # where f1 is measured by 3, f2 by 4, f3 by 5, and f4 by 6 indicators, 
#' # and all loadings are .5
#' Beta <- matrix(c(
#'   c(.00, .00, .00, .00),       # f1
#'   c(.20, .00, .00, .00),       # f2
#'   c(.00, .30, .00, .00),       # f3
#'   c(.10, .00, .40, .00)        # f4
#' ), byrow = TRUE, ncol = 4)
#' powerPath <- semPower.powerPath(type = 'a-priori',
#'                                 Beta = Beta,
#'                                 nullWhich = c(4, 1),
#'                                 nIndicator = c(3, 4, 5, 6), 
#'                                 loadM = .5,
#'                                 alpha = .05, beta = .05)
#' # show summary
#' summary(powerPath)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerPath$modelH1, sample.cov = powerPath$Sigma,
#' sample.nobs = powerPath$requiredN, sample.cov.rescale = FALSE)
#' lavaan::sem(powerPath$modelH0, sample.cov = powerPath$Sigma,
#' sample.nobs = powerPath$requiredN, sample.cov.rescale = FALSE)
#' 
#' # same as above, but detect that the slope f3 -> f4 is >= .30 
#' powerPath <- semPower.powerPath(type = 'a-priori',
#'                                 Beta = Beta,
#'                                 nullWhich = c(4, 3),
#'                                 nIndicator = c(3, 4, 5, 6), 
#'                                 loadM = .5,
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but detect that 
#' # the slope f1 -> f2 (of .20) differs from the slope f2 -> f3 (of .30) 
#' powerPath <- semPower.powerPath(type = 'a-priori',
#'                                 Beta = Beta,
#'                                 nullEffect = 'betaX = betaZ',
#'                                 nullWhich = list(c(2, 1), c(3, 2)),
#'                                 nIndicator = c(3, 4, 5, 6), 
#'                                 loadM = .5,
#'                                 alpha = .05, beta = .05)
#' 
#' # same as above, but consider a multiple group model with equally sized groups, 
#' # and obtain the required N to detect that the slope 
#' # in group 1 (of .20) differs from the one in group 2 (of .40)
#' Beta1 <- Beta2 <- matrix(c(
#'   c(.00, .00, .00, .00),       # f1
#'   c(.20, .00, .00, .00),       # f2
#'   c(.00, .30, .00, .00),       # f3
#'   c(.10, .00, .40, .00)        # f4
#' ), byrow = TRUE, ncol = 4)
#' Beta2[2, 1] <- .40
#' Beta <- list(Beta1, Beta2)
#' powerPath <- semPower.powerPath(type = 'a-priori',
#'                                 Beta = Beta,
#'                                 nullEffect = 'betaA = betaB',
#'                                 nullWhich = c(2, 1),
#'                                 nIndicator = c(3, 4, 5, 6), 
#'                                 loadM = .5,
#'                                 alpha = .05, beta = .05, N = list(1, 1))
#' }
#' @seealso [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerPath <- function(type, comparison = 'restricted',
                               Beta,
                               Psi = NULL,
                               nullEffect = 'beta = 0',
                               nullWhich = NULL, 
                               nullWhichGroups = NULL,
                               standardized = TRUE,
                               ...){
  
  args <- list(...)
  comparison <- checkComparisonModel(comparison)
  checkEllipsis(...)
  
  # we override Sigma later, so let's make sure it is not set in ellipsis argument
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma, because Sigma is determined as function of Beta (or the slopes).')
  
  # validate input
  nullEffect <- checkNullEffect(nullEffect, c('beta=0', 'betax=betaz', 'betaa=betab'))
  if(is.null(Beta)) stop('Beta may not be NULL.')
  isMultigroup <- is.list(Beta) && length(Beta) > 1
  if(!is.list(Beta)) Beta <- list(Beta)
  if(any(unlist(lapply(Beta, function(x) any(diag(x) != 0) )))) stop('All diagonal elements of Beta must be zero.')
  if(standardized && any(unlist(lapply(Beta, function(x) any(x[upper.tri(x, diag = TRUE)] != 0))))) stop('All upper triangular elements in Beta must be zero when requesting standardized parameters. Remember exogenous variables must occupy the first rows in Beta.')
  if(isMultigroup && (length(unique(unlist(lapply(Beta, ncol)))) > 1 || length(unique(unlist(lapply(Beta, nrow)))) > 1)) stop('Beta must be of same dimension for all groups') 
  lapply(Beta, function(x) checkSquare(x, 'Beta'))
  if(!is.null(Psi)){
    if(!is.list(Psi)) Psi <- list(Psi)
    if(isMultigroup && (length(unique(unlist(lapply(Psi, ncol)))) > 1 || length(unique(unlist(lapply(Psi, nrow)))) > 1)) stop('Psi must be of same dimension for all groups') 
    lapply(Psi, function(x) checkSymmetricSquare(x, 'Psi'))
    if(any(unlist(lapply(Psi, function(x) any(eigen(x)$values <= 0))))) warning('Phi is not positive definite.')
    if(ncol(Psi[[1]]) != ncol(Beta[[1]])) stop('Beta and Psi must be of same dimension.')
  }
  if(is.null(nullWhich)) stop('nullWhich must not be NULL.')
  if(any(unlist(lapply(nullWhich, function(x) any(x > ncol(Beta[[1]])))))) stop('At least one element in nullWhich is an out of bounds index concerning Beta.')
  if(nullEffect == 'betax=betaz' && any(lapply(nullWhich, function(x) length(x)) != 2)) stop('nullWhich must be a list containing vectors of size two each.')
  if(nullEffect == 'betaa=betab' && !isMultigroup) stop('Beta must be a list for multiple group comparisons.')
  if(nullEffect != 'betaa=betab' && isMultigroup) stop('Multiple group models are only valid for nullEffect = "betaA=betaB".')
  if(nullEffect != 'betax=betaz' && is.list(nullWhich)) stop('nullWhich may only contain a single vector of size two unless nullEffect = "betaX = betaZ".')
  if(nullEffect == 'beta=0' && any(unlist(lapply(Beta, function(x) x[nullWhich[1], nullWhich[2]] == 0)))) stop('nullWhich must not refer to a slope with a population value of zero.')
  if(!is.null(nullWhichGroups)) lapply(nullWhichGroups, function(x) checkBounded(x, 'All elements in nullWhichGroups', bound = c(1, length(Beta)), inclusive = TRUE))
  if(!is.list(nullWhich)) nullWhich <- list(nullWhich)
  if(isMultigroup && is.null(nullWhichGroups)) nullWhichGroups <- seq_along(Beta)
  
  
  ### get Sigma
  if(standardized){
    Phi <- lapply(seq_along(Beta), function(x) getPhi.B(Beta[[x]], Psi[[x]]))
    generated <- semPower.genSigma(Phi = if(!isMultigroup) Phi[[1]] else Phi, 
                                   useReferenceIndicator = TRUE, ...)  
  }else{  
    generated <- semPower.genSigma(Beta = if(!isMultigroup) Beta[[1]] else Beta, 
                                   Psi = if(!isMultigroup) Psi[[1]] else Psi, 
                                   useReferenceIndicator = TRUE, ...)
  }
  
  ### create model strings
  # we need to use modelTrueCFA and define regression relations here,
  # because we need labels for the H0 model and because standardized based on phi yields no regression structure 
  if(!isMultigroup) model <- generated[['modelTrueCFA']] else model <- generated[[1]][['modelTrueCFA']]
  tok <- list()
  tokH1 <- list()
  for(f in seq(ncol(Beta[[1]]))){
    ifelse(isMultigroup, idx <- unique(unlist(lapply(Beta, function(x) which(x[f, ] != 0)))), idx <- which(Beta[[1]][f, ] != 0))
    if(length(idx) > 0){
      # cIdx <- lapply(idx, function(x) c(f, x)) %in% nullWhich   ## not sure why this does not work when looping over f?
      cIdx <- unlist(lapply(lapply(idx, function(x) c(f, x)), function(y) any(unlist(lapply(nullWhich, function(z) all(y == z))))))
      if(nullEffect == 'beta=0'){
        prefix <- rep('', length(idx))
        prefix[cIdx] <- '0*'
      }else if(nullEffect == 'betax=betaz'){
        prefix <- rep('', length(idx))
        prefix[cIdx] <- 'pc*'
      }else if(nullEffect == 'betaa=betab'){
        clab <- paste0('pf', paste0(formatC(f, width = 2, flag = 0), formatC(idx, width = 2, flag = 0)))
        prefix <- lapply(clab, function(x) paste0(x, '_g', seq_along(Beta)))
        prefix[cIdx] <- lapply(prefix[cIdx], function(x) unlist(lapply(x, function(y) gsub(paste0('_g', nullWhichGroups, collapse = '|'), '_gc', y))))
        prefix <- unlist(lapply(prefix, function(x) paste0('c(', paste(x, collapse = ', ') ,')*')))
      }
      tok <- append(tok, paste0('f', f, ' ~ ', paste(prefix, paste0('f', idx), sep = '', collapse = ' + ')))
      tokH1 <- append(tokH1, paste0('f', f, ' ~ ', paste(paste0('f', idx), sep = '', collapse = ' + ')))
    }
  }
  # (residual) correlations, assuming the same Psi across groups
  if(is.null(Psi)) Psi <- list(diag(ncol(Beta[[1]])))
  for(f in 1:(ncol(Beta[[1]]) - 1)){
    for(ff in (f + 1):ncol(Beta[[1]])){
      if(Psi[[1]][f, ff] != 0){
        tok <- append(tok, paste0('f', f, ' ~~ f', ff))
        tokH1 <- append(tokH1, paste0('f', f, ' ~~ f', ff))
      }
    }
  }
  modelH0 <- paste(c(model, unlist(tok)), collapse = '\n')
  modelH1 <- paste(c(model, unlist(tokH1)), collapse = '\n')
  
  # always enforce invariance constraints in the multigroup case
  if(isMultigroup){
    args[['lavOptions']] <- append(args[['lavOptions']], list(group.equal = c('loadings')))
  } 
  
  # always fit H1 model
  fitH1model <- TRUE 
  if(comparison == 'saturated'){
    modelH1 <- NULL
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
