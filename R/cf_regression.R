#' semPower.powerRegression
#'
#' Convenience function for performing power analysis on slope(s) in a latent regression of the form Y = XB.
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param slopes vector of slopes (or a single number for a single slope) of the k predictors for Y. A list of slopes for multigroup models.
#' @param corXX correlation(s) between the k predictors (X). Either `NULL` for uncorrelated predictors, a single number (for k = 2 predictors), or a matrix. Can also be a list for multigroup models providing the correlations by group of matrices (otherwise, the same correlations are used in all groups). 
#' @param nullEffect defines the hypothesis of interest, must be one of `'slope = 0'` (the default) to test whether a slope is zero, `'slopeX = slopeZ'` to test for the equality of slopes, or `'slopeA = slopeB'` to test for the equality of slopes across groups. Define the slopes to set to equality in `nullWhich`.
#' @param nullWhich single number indicating which slope is hypothesized to equal zero when `nullEffect = 'slope = 0'`, or indicating which slope to restrict to equality across groups when `nullEffect = 'slopeA = slopeB'`, or vector defining the slopes to restrict to equality when `nullEffect = 'slopeX = slopeZ'`. Can also contain more than two slopes, e.g. `c(1, 2, 3)` to constrain the first three slopes to equality.
#' @param nullWhichGroups for `nullEffect = 'slopeA = slopeB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
#' @param standardized whether all parameters should be standardized (`TRUE`, the default). If `FALSE`, all regression relations are unstandardized.
#' @param ... mandatory further parameters related to the specific type of power analysis requested, see [semPower.aPriori()], [semPower.postHoc()], and [semPower.compromise()], and parameters specifying the factor model. The first factor is treated as Y and the subsequent factors as the predictors X_k. See details.
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
#' in SEM models involving a simple regression relation of the form `Y = b_1*X_1 + ... + b_k*X_k` between the factors:
#' * `nullEffect = 'slope = 0'`: Tests the hypothesis that the slope for a predictor is zero. 
#' * `nullEffect = 'slopeX = slopeZ'`: Tests the hypothesis that two or more slopes are equal to each other.
#' * `nullEffect = 'slopeA = slopeB'`: Tests the hypothesis that the slope for a predictor is equal in two or more groups (always assuming metric invariance).
#' 
#' For hypotheses regarding mediation effects, see [semPower.powerMediation()]. For hypothesis in autoregressive models, see  [semPower.powerAutoreg()].
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
#' Note that the first factor acts as the criterion Y, the subsequent factors as predictors X_1 to X_k.
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
#' # latent regression of the form `Y = .2*X1 + .3*X2`, where X1 and X2 correlate by .4
#' # obtain required N to reject the hypothesis that the slope of X1 is zero 
#' # with a power of 95% on alpha = 5%,   
#' # where Y is measured by 3 indicators loading by .5 each,
#' # X1 by 5 indicators loading by .6 each, and
#' # X2 by 4 indicators loading by .7 each. 
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), corXX = .4, 
#'                                      nullWhich = 1, 
#'                                      nIndicator = c(3, 5, 4), 
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#' # show summary
#' summary(powerReg)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powerReg$modelH1, sample.cov = powerReg$Sigma, 
#' sample.nobs = powerReg$requiredN, sample.cov.rescale = FALSE)
#' lavaan::sem(powerReg$modelH0, sample.cov = powerReg$Sigma, 
#' sample.nobs = powerReg$requiredN, sample.cov.rescale = FALSE)
#' 
#' # same as above, but determine power with N = 500 on alpha = .05 
#' powerReg <- semPower.powerRegression(type = 'post-hoc',
#'                                      slopes = c(.2, .3), corXX = .4, 
#'                                      nullWhich = 1, 
#'                                      nIndicator = c(3, 5, 4), 
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, N = 500)
#'                                      
#' # same as above, but determine the critical chi-square with N = 500 so that alpha = beta 
#' powerReg <- semPower.powerRegression(type = 'compromise',
#'                                      slopes = c(.2, .3), corXX = .4, 
#'                                      nullWhich = 1, 
#'                                      nIndicator = c(3, 5, 4), 
#'                                      loadM = c(.5, .6, .7),
#'                                      abratio = .05, N = 500)
#'                                      
#' # same as above, but ask for the required N to detect that the slope of X2 is zero
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), corXX = .4, 
#'                                      nullWhich = 2, 
#'                                      nIndicator = c(3, 5, 4), 
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'
#' # same as above, but define unstandardized slopes
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), corXX = .4,
#'                                      nullWhich = 2, 
#'                                      standardized = FALSE,
#'                                      nIndicator = c(3, 5, 4), 
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'                                      
#' # same as above, but compare to the saturated model
#' # (rather than to the less restricted model)
#' powerReg <- semPower.powerRegression(type = 'a-priori', 
#'                                      comparison = 'saturated',
#'                                      slopes = c(.2, .3), corXX = .4, 
#'                                      nullWhich = 2, 
#'                                      nIndicator = c(3, 5, 4), 
#'                                      loadM = c(.5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'                                      
#' # same as above, but provide a reduced loading matrix defining
#' # three indicators with loadings of .7, .6, .5 on the first factor (Y),
#' # four indicators with loadings of .5, .6, .4, .8 on the second factor (X1), and
#' # three indicators with loadings of .8, .7, .8 on the third factor (X2).
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3), corXX = .4, 
#'                                      nullWhich = 2, 
#'                                      loadings = list(
#'                                        c(.7, .6, .5), 
#'                                        c(.5, .6, .4, .8),
#'                                        c(.8, .7, .8)),
#'                                      alpha = .05, beta = .05)
#'                               
#' # latent regression of the form `Y = .2*X1 + .3*X2 + .4*X3`, 
#' # providing the predictor intercorrelation matrix,
#' # and ask for the required N to detect that the first slope differs from zero.
#' corXX <- matrix(c(
#'   #   X1    X2    X3
#'   c(1.00, 0.20, 0.30),  # X1
#'   c(0.20, 1.00, 0.10),  # X2
#'   c(0.30, 0.10, 1.00)   # X3
#' ), ncol = 3,byrow = TRUE)
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3, .4), corXX = corXX, 
#'                                      nullWhich = 1,
#'                                      nIndicator = c(4, 3, 5, 4),
#'                                      loadM = c(.5, .5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'
#' # same as above, but ask for the required N to detect that 
#' # the slope for X1 (b = .2) and the slope for X2 (b = .3) differ from each other
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3, .4), corXX = corXX, 
#'                                      nullEffect = 'slopeX = slopeZ', 
#'                                      nullWhich = c(1, 2),
#'                                      nIndicator = c(4, 3, 5, 4),
#'                                      loadM = c(.5, .5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'                                      
#' # same as above, but ask for the required N to reject the hypothesis that 
#' # all three slopes are equal to each other
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = c(.2, .3, .4), corXX = corXX, 
#'                                      nullEffect = 'slopeX = slopeZ', 
#'                                      nullWhich = c(1, 2, 3),
#'                                      nIndicator = c(4, 3, 5, 4),
#'                                      loadM = c(.5, .5, .6, .7),
#'                                      alpha = .05, beta = .05)
#'     
#' # get required N to detect that 
#' # the slope for X2 group 1 (of b2 = .3) differs from the slope for X2 in group 2 (of b = .0). 
#' # The remaining slopes are equal in both groups (b1 = .2, b3 = .4).
#' # The measurement model is identical in both groups:
#' # The criterion (Y) is measured by 4 indicators loading by .5 each, 
#' # Predictors X1 and X3 are both measured by 5 indicators loading by .6 each,
#' # Predictor X2 is measured by 3 indicators loading by .7 each.
#' # Both groups are sized equally (N = list(1, 1)).
#' powerReg <- semPower.powerRegression(type = 'a-priori',
#'                                      slopes = list(c(.2, .3, .4), 
#'                                      c(.2, .0, .4)), 
#'                                      corXX = corXX, 
#'                                      nullEffect = 'slopeA = slopeB', 
#'                                      nullWhich = 2,
#'                                      nIndicator = c(4, 5, 3, 5),
#'                                      loadM = c(.5, .6, .7, .6),
#'                                      alpha = .05, beta = .05, 
#'                                      N = list(1, 1))
#'
#'# request a simulated post-hoc power analysis with 500 replications 
#'# to detect that the slope of X1 differs from zero.
#' set.seed(300121)
#' powerReg <- semPower.powerRegression(type = 'post-hoc',
#'                                      slopes = c(.2, .1), 
#'                                      nullWhich = 1,
#'                                      nIndicator = c(4, 3, 3), loadM = .5,
#'                                      alpha = .05, N = 500, 
#'                                      simulatedPower = TRUE, 
#'                                      simOptions = list(nReplications = 500))
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerRegression <- function(type, comparison = 'restricted',
                                     slopes = NULL, 
                                     corXX = NULL, 
                                     nullEffect = 'slope = 0',
                                     nullWhich = NULL,
                                     nullWhichGroups = NULL,
                                     standardized = TRUE,
                                     ...){
  
  args <- list(...)
  comparison <- checkComparisonModel(comparison)
  checkEllipsis(...)
  
  # we override Phi and Sigma later, so let's make sure it is not set in ellipsis argument
  if('Phi' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Phi, because the factor correlations depend on corXX and the slopes.')
  if('Sigma' %in% names(match.call(expand.dots = FALSE)$...)) stop('Cannot set Sigma, because Sigma is determined as function of corXX and the slopes.')
  
  # validate input
  nullEffect <- checkNullEffect(nullEffect, c('slope=0', 'slopex=slopez', 'slopea=slopeb'))
  if(is.null(slopes)) stop('slopes cannot be NULL.')
  if(nullEffect == 'slopea=slobeb' && !is.list(slopes)) stop('slopes must be a list when a multiple group analysis is requested.')
  if(is.list(slopes)) if(length(unique(unlist(lapply(slopes, length)))) > 1) stop('the same number of slopes must be provided for each group.')
  if(!is.list(slopes)) slopes <- list(slopes)
  if(!is.null(corXX) && !is.list(corXX)) corXX <- list(corXX)

  if(any(!unlist(lapply(slopes, is.vector))) && any(!unlist(lapply(slopes, is.matrix)))) stop('slopes must be a single number or a vector')
  lapply(slopes, function(y) invisible(lapply(y, function(x) checkBounded(x, 'All slopes ', bound = c(-1, 1), inclusive = TRUE))))
  if(any(unlist(lapply(slopes, function(x) sum(x^2) > 1)))) stop('slopes imply a negative residual variance for Y, make sure that the sum of the squared slopes is < 1')
  slopes <- lapply(slopes, function(x) if(!is.matrix(x)) matrix(x, nrow = length(x)) else x)
  isMultigroup <- length(slopes) > 1

  if(is.null(corXX)) corXX <- lapply(slopes, function(x) diag(nrow(x))) 
  if(any(unlist(lapply(corXX, is.vector))) && length(unique(unlist(lapply(corXX, length)))) > 1) stop('corXX must be a single number or a matrix') 
  if(isMultigroup && length(corXX) == 1) corXX <- rep(corXX, length(slopes)) # assume same corXX for all groups
  corXX <- lapply(corXX, function(x) {
    if(!is.matrix(x)){
      xx <- matrix(x, nrow = 2, ncol = 2) 
      diag(xx) <- 1
      xx
    }else{
      x
    }
  })
  invisible(lapply(corXX, checkPositiveDefinite))
  if(any(unlist(lapply(seq_along(corXX), function(x) ncol(corXX[[x]]) != nrow(slopes[[x]]))))) stop('Dimension of corXX does not match number of predictors.')
  
  if(is.null(nullWhich)) stop('nullWhich must be defined.')
  if(any(nullWhich < 1) || any(nullWhich > nrow(slopes[[1]]))) stop('nullWhich is invalid.')
  if(nullEffect == 'slopex=slopez'){
    if(length(nullWhich) < 2 || length(nullWhich) > nrow(slopes[[1]])) stop('nullWhich must contain at least two slopes when nullEffect is slopex=slopez, but not more slopes than available')
  }else{
    if(length(nullWhich) > 1) stop('nullWhich must be a single number when nullEffect is slope=0 or slopeA=slopeB')
  }
  nullWhich <- nullWhich + 1 # because first factor is criterion
  
  # get temporary Lambda so that we can check whether number of factors matches slopes + 1
  tLambda  <- args[['Lambda']]
  if(is.null(tLambda)){
      tLambda <- genLambda(if(is.list(args[['loadings']])) args[['loadings']][[1]] else args[['loadings']], 
                           if(is.list(args[['nIndicator']])) args[['nIndicator']][[1]] else args[['nIndicator']],
                           if(is.list(args[['loadM']])) args[['loadM']][[1]] else args[['loadM']], 
                           if(is.list(args[['loadSD']])) args[['loadSD']][[1]] else args[['loadSD']], 
                           if(is.list(args[['loadMinMax']])) args[['loadMinMax']][[1]] else args[['loadMinMax']])
  }
  if(ncol(tLambda) != (1 + nrow(slopes[[1]]))) stop('The number of factors does not match the number of slopes + 1. Remember to define a measurement model including both the criterion (Y) and all predictors (X).')
  
  ### calc implied sigma. 
  # standardized case: transform B and Psi to Phi
  if(standardized){
    B <- lapply(slopes, function(x){
      cB <- matrix(0, ncol = (length(x) + 1), nrow = (length(x) + 1))
      cB[nrow(cB), ] <- c(x, 0)
      cB
    })
    Psi <- lapply(seq_along(slopes), function(x){
      cPsi <- diag(ncol(B[[x]]))
      cPsi[1:(nrow(cPsi) - 1), 1:(ncol(cPsi) - 1)] <- corXX[[x]]
      cPsi
    })
    cPhi <- lapply(seq_along(B), function(x) getPhi.B(B[[x]], Psi[[x]]))
    # change order so that Y is first factor
    Phi <- lapply(cPhi, function(x){
      nf <- ncol(x)
      nPhi <- diag(nf)
      nPhi[2:nf, 2:nf] <- x[1:(nf - 1), 1:(nf - 1)]
      nPhi[1, 1:nf] <- c(1, x[nf, 1:(nf-1)])
      nPhi[1:nf, 1] <- c(1, x[1:(nf-1), nf])
      nPhi
    })
    
    generated <- semPower.genSigma(Phi = if(length(Phi) > 1) Phi else Phi[[1]], 
                                   useReferenceIndicator = TRUE, ...)  
    
  # unstandardized case
  }else{
    B <- lapply(slopes, function(x){
      cB <- matrix(0, ncol = (length(x) + 1), nrow = (length(x) + 1))
      cB[1, ] <- c(0, x)  # we want the first factor endogenous
      cB
    })
    Psi <- lapply(seq_along(slopes), function(x){
      cPsi <- diag(ncol(B[[x]]))
      cPsi[2:nrow(cPsi), 2:ncol(cPsi)] <- corXX[[x]]
      cPsi
    })
    generated <- semPower.genSigma(Beta = if(length(B) > 1) B else B[[1]], 
                                   Psi = if(length(Psi) > 1) Psi else Psi[[1]],
                                   useReferenceIndicator = TRUE, ...)
  }


  ### create ana model string
  # add regressions 
  if(isMultigroup) model <- generated[[1]][['modelTrueCFA']] else model <- generated[['modelTrueCFA']]
  np <- (1 + 1:ncol(corXX[[1]]))
  
  if(nullEffect == 'slopex=slopez'){
    tok <- ''
    for(i in 1:(length(nullWhich) - 1)){
      for(j in (i + 1):length(nullWhich)){
        tok <- paste(tok, paste0('pf', nullWhich[i], ' == ', 'pf', nullWhich[j]), sep = '\n')
      }
    }
    modelH0 <- paste(model, 
                     paste0('f1 ~ ', paste0(paste0('pf', np), '*f', np, collapse = ' + ')),
                     tok,
                     sep = '\n')
  }else if(nullEffect == 'slope=0'){
    modelH0 <- paste(model, 
                     paste0('f1 ~ ', paste0(paste0('pf', np), '*f', np, collapse = ' + ')),
                     paste0('pf', nullWhich,' == 0'),
                     sep = '\n')
  }else if(nullEffect == 'slopea=slopeb'){
    if(is.null(nullWhichGroups)) nullWhichGroups <- seq_along(slopes)
    lab <- paste0('ff', seq_along(slopes))
    lab[nullWhichGroups] <- 'pf1'
    lab <- paste0('c(', paste(lab, collapse = ','), ')*')
    modelH0 <- paste(model, 
                     paste0('f1 ~ ', paste0(lab, 'f', np[(nullWhich - 1)], ' + '), paste0('f',np[-(nullWhich - 1)], collapse = ' + ')),
                     sep = '\n')
  }else{
    stop('nullEffect not defined.')
  }

  # we always enforce metric invariance in the multigroup case
  if(isMultigroup){
    args[['lavOptions']] <- append(args[['lavOptions']], list(group.equal = c('loadings')))
  } 
  
  modelH1 <- NULL
  fitH1model <- FALSE
  if(comparison == 'restricted'){
    if(!isMultigroup){
      modelH1 <- paste(model, 
                       paste0('f1 ~ ', paste0(paste0('pf',(1 + 1:ncol(corXX[[1]]))), '*f',(1 + 1:ncol(corXX[[1]])), collapse = ' + ')),
                       sep = '\n')
    }else{
      # no slope labels in multigroup case
      modelH1 <- paste(model, 
                       paste0('f1 ~ ', paste0('f',(1 + 1:ncol(corXX[[1]])), collapse = ' + ')),
                       sep = '\n')
      
    }
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