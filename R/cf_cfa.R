
#' semPower.powerCFA
#'
#' Convenience function for performing power analyses for CFA models to reject one of the following hypotheses: 
#' (a) a zero correlation between two factors, (b) the equality of two correlations between factors,
#' or (c) the equality of a correlation between two factors across two or more groups. 
#' This requires the lavaan package.
#' 
#' @param type type of power analysis, one of `'a-priori'`, `'post-hoc'`, `'compromise'`.
#' @param comparison comparison model, one of `'saturated'` or `'restricted'` (the default). This determines the df for power analyses. `'saturated'` provides power to reject the model when compared to the saturated model, so the df equal the one of the hypothesized model. `'restricted'` provides power to reject the hypothesized model when compared to an otherwise identical model that just omits the restrictions defined in `nullEffect`, so the df equal the number of restrictions.
#' @param Phi either a single number defining the correlation between exactly two factors or the factor correlation matrix. A list for multiple group models.
#' @param nullEffect defines the hypothesis of interest, must be one of `'cor = 0'` (the default) to test whether a correlation is zero, `'corX = corZ'` to test for the equality of correlations, `'corA = corB'` to test for the equality of a correlation across groups, and `loading = 0` to test whether a loading is zero. Define the correlations to be set to equality in `nullWhich` and the groups in `nullWhichGroups`. 
#' @param nullWhich vector of size 2 indicating which element in Lambda should equal zero when `nullEffect = 'loading = 0'`, or which factor correlation in `Phi` is hypothesized to equal zero when `nullEffect = 'cor = 0'`, or to restrict to equality across groups when `nullEffect = 'corA = corB'`, or list of vectors defining which correlations to restrict to equality when `nullEffect = 'corX = corZ'`. Can also contain more than two correlations, e.g., `list(c(1, 2), c(1, 3), c(2, 3))` to set `Phi[1, 2] = Phi[1, 3] = Phi[2, 3]`. If omitted, the correlation between the first and the second factor is targeted, i. e., `nullWhich = c(1, 2)`.
#' @param nullWhichGroups for `nullEffect = 'corA = corB'`, vector indicating the groups for which equality constrains should be applied, e.g. `c(1, 3)` to constrain the relevant parameters of the first and the third group. If `NULL`, all groups are constrained to equality.
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
#' 
#' This function performs a power analysis to reject various hypotheses arising
#' in standard CFA models:
#' * `nullEffect = 'loading = 0'`: Tests the hypothesis that a loading is zero.
#' * `nullEffect = 'cor = 0'`: Tests the hypothesis that the correlation between two factors is zero. 
#' * `nullEffect = 'corX = corZ'`: Tests the hypothesis that two or more correlations between three or more factors are equal to each other.
#' * `nullEffect = 'corA = corB'`: Tests the hypothesis that the correlation between two factors is equal in two or more groups (always assuming metric invariance).
#' 
#' For hypotheses regarding regression relationships between factors, see [semPower.powerRegression()].
#' For hypotheses regarding mediation effects, see [semPower.powerMediation()].
#' For hypotheses regarding measurement invariance, see [semPower.powerMI()].
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
#' * `missingVars`: vector specifying the variables containing missing data (defaults to NULL).
#' * `missingVarProp`: can be used instead of `missingVars`: The proportion of variables containing missing data (defaults to zero).
#' * `missingProp`: The proportion of missingness for variables containing missing data (defaults to zero), either a single value or a vector giving the probabilities for each variable.
#' * `missingMechanism`: The missing data mechanism, one of `MCAR` (the default), `MAR`, or `NMAR`.
#' The approaches generating non-normal data require additional arguments detailed below.
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
#' # get required N to detect a correlation of >= .2 between two factors
#' # with a power of 95% on alpha = 5%, where the factors are  
#' # measured by 5 and 6 indicators, respectively, and all loadings are equal to .5
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               Phi = .2, 
#'                               nIndicator = c(5, 6), loadM = .5,
#'                               alpha = .05, beta = .05)
#' # show summary
#' summary(powercfa)
#' # optionally use lavaan to verify the model was set-up as intended
#' lavaan::sem(powercfa$modelH1, sample.cov = powercfa$Sigma, 
#' sample.nobs = powercfa$requiredN, sample.cov.rescale = FALSE)
#' lavaan::sem(powercfa$modelH0, sample.cov = powercfa$Sigma, 
#' sample.nobs = powercfa$requiredN, sample.cov.rescale = FALSE)
#' 
#' # same as above, but determine power with N = 500 on alpha = .05
#' powercfa <- semPower.powerCFA(type = 'post-hoc',
#'                               Phi = .2, 
#'                               nIndicator = c(5, 6), loadM = .5,
#'                               alpha = .05, N = 500)
#' 
#' # same as above, but determine the critical chi-square with N = 500 so that alpha = beta
#' powercfa <- semPower.powerCFA(type = 'compromise',
#'                               Phi = .2, 
#'                               nIndicator = c(5, 6), loadM = .5,
#'                               abratio = 1, N = 500)
#'                               
#' # same as above, but compare to the saturated model 
#' # (rather than to the less restricted model)
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               comparison = 'saturated',
#'                               Phi = .2, 
#'                               nIndicator = c(5, 6), loadM = .5,
#'                               alpha = .05, beta = .05)
#'                               
#' # same as above, but provide a reduced loading matrix defining
#' # three indicators with loadings of .7, .6, and .5 on the first factor and
#' # four indicators with loadings of .5, .6, .4, .8 on the second factor 
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               Phi = .2, 
#'                               loadings = list(c(.7, .6, .5), 
#'                                               c(.5, .6, .4, .8)),
#'                               alpha = .05, beta = .05)
#'
#' # detect that the loading of indicator 4 on the first factor differs from zero
#' Lambda <- matrix(c(
#'   c(.8, 0),
#'   c(.4, 0),
#'   c(.6, 0),
#'   c(.1, .5),
#'   c(0, .6),
#'   c(0, .7)
#' ), ncol = 2, byrow = TRUE)
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               Phi = .2,
#'                               nullEffect = 'loading = 0', 
#'                               nullWhich = c(4, 1), 
#'                               Lambda = Lambda,
#'                               alpha = .05, beta = .05)
#'
#'
#' # get required N to detect a correlation of >= .3 between factors 1 and 3  
#' # in a three factor model. Factors are measured by 3 indicators each, and all loadings 
#' # on the first, second, and third factor are .5, .6, and .7, respectively.
#' Phi <- matrix(c(
#'   c(1.00, 0.20, 0.30),
#'   c(0.20, 1.00, 0.10),
#'   c(0.30, 0.10, 1.00)
#' ), ncol = 3,byrow = TRUE)
#' 
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               Phi = Phi,
#'                               nullWhich = c(1, 3), 
#'                               nIndicator = c(3, 3, 3), loadM = c(.5, .6, .7),
#'                               alpha = .05, beta = .05)
#' 
#' # same as above, but ask for N to detect that 
#' # the correlation between factors 1 and 2 (of r = .2) differs from
#' # the correlation between factors 2 and 3 (of r = .3).
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               Phi = Phi,
#'                               nullEffect = 'corX = corZ',
#'                               nullWhich = list(c(1, 2), c(1, 3)), 
#'                               nIndicator = c(3, 3, 3), loadM = c(.5, .6, .7),
#'                               alpha = .05, beta = .05)
#'                               
#' # same as above, but ask for N to detect that all three correlations are unequal
#' powercfa <- semPower.powerCFA(type = 'a-priori',
#'                               Phi = Phi,
#'                               nullEffect = 'corX = corZ',
#'                               nullWhich = list(c(1, 2), c(1, 3), c(2, 3)), 
#'                               nIndicator = c(3, 3, 3), loadM = c(.5, .6, .7),
#'                               alpha = .05, beta = .05)
#'                               
#' # get required N to detect that the correlation between two factors
#' # in group 1 (of r = .2) differs from the one in group 2 (of r = .4). 
#' # The measurement model is identical for both groups:
#' # The first factor is measured by 3 indicators loading by .7 each, 
#' # the second factor is measured by 6 indicators loading by .5 each.
#' # Both groups are sized equally (N = list(1, 1)).
#' powercfa <- semPower.powerCFA(type = 'a-priori', 
#'                               nullEffect = 'corA = corB',
#'                               Phi = list(.2, .4), 
#'                               loadM = c(.7, .5), 
#'                               nIndicator = c(3, 6), 
#'                               alpha = .05, beta = .05, N = list(1, 1))
#'
#' # request a simulated post-hoc power analysis with 500 replications.
#' set.seed(300121)
#' powercfa <- semPower.powerCFA(type = 'post-hoc',
#'                               Phi = .2, 
#'                               nIndicator = c(5, 6), loadM = .5,
#'                               alpha = .05, N = 500, 
#'                               simulatedPower = TRUE, 
#'                               simOptions = list(nReplications = 500))
#' 
#' }
#' @seealso [semPower.genSigma()] [semPower.aPriori()] [semPower.postHoc()] [semPower.compromise()]
#' @export
semPower.powerCFA <- function(type, comparison = 'restricted', 
                              Phi = NULL,
                              nullEffect = 'cor = 0',
                              nullWhich = NULL, 
                              nullWhichGroups = NULL, 
                              ...){
  args <- list(...)
  
  # validate input
  checkEllipsis(...)
  comparison <- checkComparisonModel(comparison)
  nullEffect <- checkNullEffect(nullEffect, c('cor=0', 'corx=corz', 'cora=corb', 'loading=0'))
  if(nullEffect == 'loading=0' && is.null(Phi)) Phi <- diag(1)
  if(is.null(Phi)) stop('Phi must be defined')
  if(!is.null(nullWhichGroups) && !is.list(Phi)) stop('Phi must be provided for each group.')
  if(nullEffect == 'cora=corb' && !is.list(Phi)) stop('corA=corB refers to multigroup analysis, so Phi must be a list.')
  if(nullEffect == 'corx=corz' && (!is.matrix(Phi) || ncol(Phi) <= 2)) stop('corx=corz compares two correlations and thus requires at least three factors and a correlation matrix for Phi.')
  if(is.list(Phi) && !is.null(nullWhichGroups)) lapply(as.list(nullWhichGroups), function(x) checkBounded(x, 'Each element in nullWhichGroups', bound = c(1, length(Phi)), inclusive = TRUE)) 
    
  # generate sigma 
  generated <- semPower.genSigma(Phi = Phi, ...)
  
  ### now do validation of nullWhich, since we now know Phi
  isMultigroup <- is.list(Phi)
  if(isMultigroup) nfac <- ncol(generated[[1]][['Phi']]) else nfac <- ncol(generated[['Phi']])
  if(nullEffect != 'loading=0'){
    if(is.null(nullWhich) && nfac == 2) nullWhich <- c(1, 2)
    if(is.null(nullWhich)) stop('nullWhich must be defined.')
    if(!is.list(nullWhich)) nullWhich <- list(nullWhich)
    if(any(unlist(lapply(nullWhich, function(x) length(x) != 2)))) stop('nullWhich may only contain vectors of size two.')
    if(any(unlist(lapply(nullWhich, function(x) x[1] == x[2])))) stop('elements in nullWhich may not refer to variances.')
    if(any(unlist(lapply(nullWhich, function(x) (x[1] < 1 || x[2] < 1 || x[1] > nfac || x[2] > nfac))))) stop('At least one element in nullWhich is an out of bounds index concerning Phi.')
    if(length(nullWhich) > 1){
      for(i in 1:(length(nullWhich) - 1)){
        for(j in (i + 1):length(nullWhich)){
          if(nullWhich[[i]][1] == nullWhich[[j]][1] && nullWhich[[i]][2] == nullWhich[[j]][2]) stop('elements in nullWhich may not refer to the same correlation')
        }
      }
    }
  }else{
    if(is.null(nullWhich)) stop('nullWhich must be defined.')
    if(length(nullWhich) != 2) stop('nullWhich must be a vector with two elements.')
    if(nullWhich[1] > nrow(generated[['Lambda']]) || nullWhich[2] > ncol(generated[['Lambda']])  ) stop('nullWhich refers to an invalid element in Lambda. The first entry must be <= nIndicator, the second <= nFactors.')
    if(generated[['Lambda']][nullWhich[1], nullWhich[2]] == 0) stop('nullWhich refers to a loading with a population value of zero. The loading referred by nullWhich must differ from zero.')
  }

  # warn if loadings dont satisfy metric invariance
  if(isMultigroup){
    lambdas <- sapply(generated, '[[', 'Lambda')
    if(any(apply(lambdas, 1, function(x) length(unique(x)) != 1))) warning('At least one loading differs across groups, violating metric invariance. Verify that this is intended.')
  }
  
  ### H0 model
  if(isMultigroup) modelH0 <- generated[[1]][['modelTrueCFA']] else modelH0 <- generated[['modelTrueCFA']] 
  if(nullEffect == 'cor=0'){
    modelH0 <- paste(c(modelH0,
      paste0('f', nullWhich[[1]], collapse = ' ~~ 0*')),
      collapse = '\n')
  }else if(nullEffect == 'corx=corz'){
    labs <- list()
    tok <- ''
    for(i in 1:length(nullWhich)){
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
    modelH0 <- paste(c(modelH0, tok), collapse = '\n')
  }else if(nullEffect == 'cora=corb'){
    if(is.null(nullWhichGroups)) nullWhichGroups <- 1:length(Phi)
    lab <- paste0('ff', 1:length(Phi))
    lab[nullWhichGroups] <- 'pf1'
    lab <- paste0('c(', paste(lab, collapse = ','), ')*')
    modelH0 <- paste(c(modelH0,
                       paste0('f', nullWhich[[1]], collapse = paste0(' ~~ ', lab))),
                       collapse = '\n')
  }else if(nullEffect == 'loading=0'){
    tInd <- paste0('x', nullWhich[1])
    tFac <- paste0('f', nullWhich[2])
    tok <- strsplit(modelH0, '\n', fixed = TRUE)[[1]]
    tok <- lapply(tok, function(x){
      if(startsWith(x, tFac) && grepl('=~', x)){
        t <- lapply(strsplit(strsplit(x, '=~', fixed = TRUE)[[1]][2], '+', fixed = TRUE)[[1]], trimws)
        idx <- grep(tInd, unlist(t))
        if(grepl('NA', t[[idx]])){
          t[[idx]] <- sub('NA', '0', t[[idx]])
        }else{
          t[[idx]] <- paste0('0*', t[[idx]])
        }
        paste0(tFac, ' =~ ', paste0(unlist(t), collapse = ' + '))
      }else{
        x
      }
    })
    modelH0 <- paste(unlist(tok), collapse = '\n')
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
    if(isMultigroup) modelH1 <- generated[[1]][['modelTrueCFA']] else modelH1 <- generated[['modelTrueCFA']] 
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