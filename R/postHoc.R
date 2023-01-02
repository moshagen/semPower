#' semPower.postHoc
#'
#' Performs a post-hoc power analysis, i. e., determines power (= 1 - beta) given alpha, df, and and a measure of effect.
#'
#' @param effect effect size specifying the discrepancy between the null hypothesis (H0) and the alternative hypothesis (H1). A list for multiple group models; a vector of length 2 for effect-size differences.
#' @param effect.measure type of effect, one of `"F0"`, `"RMSEA"`, `"Mc"`, `"GFI"`, `"AGFI"`. Can be `NULL` if `Sigma` and `SigmaHat` are set.
#' @param alpha alpha error
#' @param N the number of observations (a list for multiple group models)
#' @param df the model degrees of freedom. See [semPower.getDf()] for a way to obtain the df of a specific model. 
#' @param p the number of observed variables, only required for `effect.measure = "GFI"` and `effect.measure = "AGFI"`.
#' @param SigmaHat can be used instead of `effect` and `effect.measure`: model implied covariance matrix (a list for multiple group models). Used in conjunction with `Sigma` to define the effect.
#' @param Sigma can be used instead of `effect` and `effect.measure`: population covariance matrix (a list for multiple group models). Used in conjunction with `SigmaHat` to define effect.
#' @param muHat can be used instead of `effect` and `effect.measure`: model implied mean vector. Used in conjunction with `mu`. If `NULL`, no meanstructure is involved.
#' @param mu can be used instead of `effect` and `effect.measure`: observed (or population) mean vector. Use in conjunction with `muHat`. If `NULL`, no meanstructure is involved.
#' @param simulatedPower whether to perform a simulated (`TRUE`, rather than analytical, `FALSE`) power analysis. Only available if `Sigma` and `SigmaHat` are defined.
#' @param modelH0 for simulated power: `lavaan` model string defining the (incorrect) analysis model.
#' @param modelH1 for simulated power: `lavaan` model string defining the comparison model. If omitted, the saturated model is the comparison model.
#' @param nReplications for simulated power: number of random samples drawn.
#' @param minConvergenceRate for simulated power: the minimum convergence rate required.
#' @param lavOptions a list of additional options passed to `lavaan`, e. g., `list(estimator = 'mlm')` to request robust ML estimation.
#' @param lavOptionsH1 alternative options passed to `lavaan` that are only used for the H1 model. If `NULL`, identical to `lavOptions`. Probably only useful for multigroup models.
#' @param seed for simulated power: the seed.
#' @return Returns a list. Use `summary()` to obtain formatted results.
#' @examples
#' \dontrun{
#' # achieved power with a sample of N = 250 to detect misspecifications corresponding
#' # to RMSEA >= .05 on 200 df on alpha = .05.
#' ph <- semPower.postHoc(effect = .05, effect.measure = "RMSEA", 
#'                        alpha = .05, N = 250, df = 200)
#' summary(ph)
#' 
#' # power analysis for to detect the difference between a model (with df = 200) exhibiting RMSEA = .05
#' # and a model (with df = 210) exhibiting RMSEA = .06.
#' ph <- semPower.postHoc(effect = c(.05, .06), effect.measure = "RMSEA", 
#'                        alpha = .05, N = 500, df = c(200, 210))
#' summary(ph)
#' 
#' # multigroup example
#' ph <- semPower.postHoc(effect = list(.02, .01), effect.measure = "F0", 
#'                         alpha = .05, N = list(250, 350), df = 200)
#' summary(ph)
#' 
#' # power analysis based on SigmaHat and Sigma (nonsense example)
#' ph <- semPower.postHoc(alpha = .05, N = 1000, df = 5,  
#'                        SigmaHat = diag(4), Sigma = cov(matrix(rnorm(4*1000),  ncol=4)))
#' summary(ph)
#' 
#' # simulated power analysis (nonsense example)
#' ph <- semPower.aPriori(alpha = .05, N = 500, df = 200,  
#'                        SigmaHat = list(diag(4), diag(4)), 
#'                        Sigma = list(cov(matrix(rnorm(4*1000), ncol=4)), cov(matrix(rnorm(4*1000), ncol=4))),
#'                        simulatedPower = TRUE, nReplications = 100)
#' summary(ph)
#' }
#' @seealso [semPower.aPriori()] [semPower.compromise()]
#' @importFrom stats qchisq pchisq 
#' @export
semPower.postHoc <- function(effect = NULL, effect.measure = NULL, alpha,
                             N, df = NULL, p = NULL,
                             SigmaHat = NULL, Sigma = NULL, muHat = NULL, mu = NULL,
                             simulatedPower = FALSE, 
                             modelH0 = NULL, modelH1 = NULL,
                             nReplications = 250, minConvergenceRate = .5, 
                             lavOptions = NULL, lavOptionsH1 = lavOptions,
                             seed = NULL,
                             ...){

  # validate input and do some preparations
  pp <- powerPrepare(type = 'post-hoc', 
                     effect = effect, effect.measure = effect.measure,
                     alpha = alpha, beta = NULL, power = NULL, abratio = NULL,
                     N = N, df = df, p = p,
                     SigmaHat = SigmaHat, Sigma = Sigma, muHat = muHat, mu = mu,
                     simulatedPower = simulatedPower, 
                     modelH0 = modelH0, modelH1 = modelH1,
                     nReplications = nReplications, minConvergenceRate = minConvergenceRate,
                     lavOptions = lavOptions)

  if(!simulatedPower){
    
    df <- pp[['df']]
    fmin <- pp[['fmin']]
    fmin.g <- pp[['fmin.g']]
    nrep <- NULL
    
    fit <- getIndices.F(fmin, df, pp[['p']], pp[['SigmaHat']], pp[['Sigma']], pp[['muHat']], pp[['mu']], pp[['N']])
    ncp <- getNCP(fmin.g, pp[['N']])
    
    beta <- pchisq(qchisq(alpha, df, lower.tail = FALSE), df, ncp = ncp)
    power <- pchisq(qchisq(alpha, df, lower.tail = FALSE), df, ncp = ncp, lower.tail = FALSE)
    
  }else{
    
    if(!is.null(seed)) set.seed(seed)
    sim <- simulate(modelH0 = modelH0, modelH1 = modelH1,
                    Sigma = pp[['Sigma']], mu = pp[['mu']], N = pp[['N']], alpha = alpha,
                    nReplications = nReplications, minConvergenceRate = minConvergenceRate,
                    lavOptions = lavOptions, lavOptionsH1 = lavOptionsH1)
    nrep <- sim[['nrep']]
    df <- sim[['df']]
    fmin <- fmin.g <- sim[['meanFmin']]
    if(!is.null(sim[['meanFminGroups']])) fmin.g <- sim[['meanFminGroups']]
    
    fit <- getIndices.F(fmin = fmin, df = df, p = pp[['p']], N = pp[['N']])
    ncp <- getNCP(fmin.g, pp[['N']])
    
    beta <- 1 - sim[['ePower']]
    power <- sim[['ePower']]
    
  }
  
  impliedAbratio <- alpha / beta

  result <- list(
    type = "post-hoc",
    alpha = alpha,
    beta = beta,
    power = power,
    impliedAbratio = impliedAbratio,
    ncp = ncp,
    fmin = fmin,
    effect = pp[['effect']],
    effect.measure = pp[['effect.measure']],
    N = pp[['N']],
    df = df,
    p = pp[['p']],
    chiCrit = qchisq(alpha, df, ncp = 0, lower.tail = FALSE),
    rmsea = fit[['rmsea']],
    mc = fit[['mc']],
    gfi = fit[['gfi']],
    agfi = fit[['agfi']],
    srmr = fit[['srmr']],
    cfi = fit[['cfi']],
    simulated = simulatedPower,
    nrep = nrep
  )

  class(result) <- "semPower.postHoc"
  result

}



#' semPower.postHoc.summary
#'
#' provide summary of post-hoc power analyses
#' @param object result object from semPower.posthoc
#' @param ... other
#' @export
summary.semPower.postHoc <- function(object, ...){

  out.table <- getFormattedResults('post-hoc', object)

  cat("\n semPower: Post-hoc power analysis\n")
  if(object[['simulated']]){
    cat(paste("\n Simulated power based on", object[['nrep']], "successful replications.\n"))
  }

  print(out.table, row.names = FALSE, right = FALSE)

  semPower.showPlot(chiCrit = object[['chiCrit']], ncp = object[['ncp']], df = object[['df']])
  
}

