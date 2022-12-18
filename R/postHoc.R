#' semPower.postHoc
#'
#' Determine power (1-beta) given alpha, df, and effect
#'
#' @param effect effect size specifying the discrepancy between H0 and H1 (a list for multiple group models; a vector of length 2 for effect-size differences)
#' @param effect.measure type of effect, one of "F0", "RMSEA", "Mc", "GFI", AGFI"
#' @param alpha alpha error
#' @param N the number of observations (a list for multiple group models)
#' @param df the model degrees of freedom 
#' @param p the number of observed variables, required for effect.measure = "GFI" and "AGFI"
#' @param SigmaHat model implied covariance matrix (a list for multiple group models). Use in conjunction with Sigma to define effect and effect.measure. 
#' @param Sigma population covariance matrix (a list for multiple group models). Use in conjunction with SigmaHat to define effect and effect.measure.
#' @param muHat model implied mean vector
#' @param mu observed (or population) mean vector
#' @param simulatedPower whether to perform a simulated (TRUE) (rather than analytical, FALSE) power analysis.
#' @param modelH0 for simulated power: lavaan model string defining the (incorrect) analysis model.
#' @param modelH1 for simulated power: lavaan model string defining the comparison model. If omitted, the saturated model is the comparison model.
#' @param fitH1model for simulated power: whether to fit the H1 model. If FALSE, the H1 model is assumed to show the same fit as the saturated model, and only the delta df are computed.
#' @param nReplications for simulated power: number of random samples drawn.
#' @param minConvergenceRate for simulated power: the minimum convergence rate required
#' @param lavOptions for simulated power: a list of additional options passed to lavaan, e.g., list(estimator = 'mlm') to request robust ML estimation
#' @param seed for simulated power: seed, by default 30012021 
#' @return list
#' @examples
#' \dontrun{
#' power <- semPower.postHoc(effect = .05, effect.measure = "RMSEA", alpha = .05, N = 250, df = 200)
#' summary(power)
#' power <- semPower.postHoc(effect = list(.02, .01), effect.measure = "F0", 
#'                           alpha = .05, N = list(250, 350), df = 200)
#' summary(power)
#' power <- semPower.postHoc(N = 1000, df = 5, alpha = .05,  
#'                           SigmaHat = diag(4), Sigma = cov(matrix(rnorm(4*1000),  ncol=4)))
#' summary(power)
#' }
#' @importFrom stats qchisq pchisq 
#' @export
semPower.postHoc <- function(effect = NULL, effect.measure = NULL, alpha,
                             N, df = NULL, p = NULL,
                             SigmaHat = NULL, Sigma = NULL, muHat = NULL, mu = NULL,
                             simulatedPower = FALSE, 
                             modelH0 = NULL, modelH1 = NULL, fitH1model = TRUE,
                             nReplications = 100, minConvergenceRate = .5, lavOptions = NULL, 
                             seed = 30012021,
                             ...){

  # validate input and do some preparations
  pp <- powerPrepare(type = 'post-hoc', 
                     effect = effect, effect.measure = effect.measure,
                     alpha = alpha, beta = NULL, power = NULL, abratio = NULL,
                     N = N, df = df, p = p,
                     SigmaHat = SigmaHat, Sigma = Sigma, muHat = muHat, mu = mu,
                     simulatedPower = simulatedPower, 
                     modelH0 = modelH0, modelH1 = modelH1, fitH1model = fitH1model,
                     nReplications = nReplications, minConvergenceRate = minConvergenceRate,
                     lavOptions = lavOptions)

  if(!simulatedPower){
    
    df <- pp$df
    fmin <- pp$fmin
    fmin.g <- pp$fmin.g
    nrep <- NULL
    
    fit <- getIndices.F(fmin, df, pp$p, SigmaHat, Sigma, muHat, mu, pp$N)
    ncp <- getNCP(fmin.g, pp$N)
    
    beta <- pchisq(qchisq(alpha, df, lower.tail = FALSE), df, ncp = ncp)
    power <- pchisq(qchisq(alpha, df, lower.tail = FALSE), df, ncp = ncp, lower.tail = FALSE)
    
  }else{
    
    set.seed(seed)
    sim <- simulate(modelH0 = modelH0, modelH1 = modelH1, fitH1model = fitH1model,
                    Sigma = Sigma, mu = mu, N = N, alpha = alpha,
                    nReplications = nReplications, minConvergenceRate = minConvergenceRate,
                    lavOptions = lavOptions)
    nrep <- sim$nrep
    df <- sim$df
    fmin <- sim$medianF
    fmin.g <- sim$medianF   ## TODO add multigroup support
    
    fit <- getIndices.F(fmin = fmin, df = df, p = pp$p, N = pp$N)
    ncp <- getNCP(fmin.g, pp$N)
    
    beta <- 1 - sim$ePower
    power <- sim$ePower
    
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
    effect = pp$effect,
    effect.measure = pp$effect.measure,
    N = pp$N,
    df = df,
    p = pp$p,
    chiCrit = qchisq(alpha, df, ncp = 0, lower.tail = FALSE),
    rmsea = fit$rmsea,
    mc = fit$mc,
    gfi = fit$gfi,
    agfi = fit$agfi,
    srmr = fit$srmr,
    cfi = fit$cfi,
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
  if(object$simulated){
    cat(paste("\n Simulated power based on", object$nrep, "successful replications.\n"))
  }

  print(out.table, row.names = FALSE, right = FALSE)

  semPower.showPlot(chiCrit = object$chiCrit, ncp = object$ncp, df = object$df)
  
}

