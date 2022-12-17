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
                             N, df, p = NULL,
                             SigmaHat = NULL, Sigma = NULL, muHat = NULL, mu = NULL,
                             ...){

  # validate input and do some preparations
  pp <- powerPrepare(type = 'post-hoc', 
                     effect = effect, effect.measure = effect.measure,
                     alpha = alpha, beta = NULL, power = NULL, abratio = NULL,
                     N = N, df = df, p = p,
                     SigmaHat = SigmaHat, Sigma = Sigma, muHat = muHat, mu = mu)
  
  fit <- getIndices.F(pp$fmin, pp$df, pp$p, SigmaHat, Sigma, muHat, mu, pp$N)
  ncp <- getNCP(pp$fmin.g, pp$N)

  beta <- pchisq(qchisq(alpha, pp$df, lower.tail = FALSE), pp$df, ncp = ncp)
  power <- pchisq(qchisq(alpha, pp$df, lower.tail = FALSE), pp$df, ncp = ncp, lower.tail = FALSE)
  impliedAbratio <- alpha / beta


  result <- list(
    type = "post-hoc",
    alpha = alpha,
    beta = beta,
    power = power,
    impliedAbratio = impliedAbratio,
    ncp = ncp,
    fmin = pp$fmin,
    effect = pp$effect,
    effect.measure = pp$effect.measure,
    N = pp$N,
    df = pp$df,
    p = pp$p,
    chiCrit = qchisq(alpha, df, ncp = 0, lower.tail = FALSE),
    rmsea = fit$rmsea,
    mc = fit$mc,
    gfi = fit$gfi,
    agfi = fit$agfi,
    srmr = fit$srmr,
    cfi = fit$cfi
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

  print(out.table, row.names = FALSE, right = FALSE)

  semPower.showPlot(chiCrit = object$chiCrit, ncp = object$ncp, df = object$df)
  
}

