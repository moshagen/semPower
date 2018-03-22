
#' semPower.postHoc
#'
#' Determine power (1-beta) given alpha, df, and effect
#'
#' @param effect effect size specifying the discrepancy between H0 and H1
#' @param effect.measure type of effect, one of "F0", "RMSEA", "Mc", "GFI", AGFI"
#' @param alpha alpha error
#' @param N the number of observations
#' @param df the model degrees of freedom
#' @param p the number of observed variables, required for effect.measure = "GammaHat", "GFI",  and "AGFI"
#' @param SigmaHat model implied covariance matrix. Use in conjuntion with Sigma to define effect and effect.measure. 
#' @param Sigma population covariance matrix. Use in conjuntion with SigmaHat to define effect and effect.measure.
#' @return list
#' @examples
#' \dontrun{
#' power <- semPower.postHoc(effect = .05, effect.measure = "RMSEA", alpha = .05, N = 250, df = 200)
#' power
#' power <- semPower.postHoc(N = 1000, df = 5, SigmaHat = diag(4), 
#'                           Sigma = cov(matrix(rnorm(4*1000),  ncol=4)))
#' power
#' }
#' @importFrom stats qchisq pchisq 
#' @export
semPower.postHoc <- function(effect = NULL, effect.measure = NULL, alpha,
                             N, df, p = NULL,
                             SigmaHat = NULL, Sigma = NULL){

  effect.measure <- toupper(effect.measure)
  
  validateInput('post-hoc', effect = effect, effect.measure = effect.measure,
                alpha = alpha, beta = NULL, power = NULL, abratio = NULL,
                N = N, df = df, p = p,
                SigmaHat = SigmaHat, Sigma = Sigma)

  if(!is.null(SigmaHat)){ # sufficient to check for on NULL matrix; primary validity check is in validateInput
    effect.measure <- 'F0'
    p <- ncol(SigmaHat)
  }

  fmin <- getF(effect, effect.measure, df, p, SigmaHat, Sigma)
  fit <- getIndices.F(fmin, df, p, SigmaHat, Sigma)
  ncp <- getNCP(fmin, N)

  beta <- pchisq(qchisq(alpha, df, lower.tail = F), df, ncp=ncp)
  power <- pchisq(qchisq(alpha, df, lower.tail = F), df, ncp=ncp, lower.tail = F)
  impliedAbratio <- alpha/beta


  result <- list(
    type = "post-hoc",
    alpha = alpha,
    beta = beta,
    power = power,
    impliedAbratio = impliedAbratio,
    ncp = ncp,
    fmin = fmin,
    effect = effect,
    effect.measure = effect.measure,
    N = N,
    df = df,
    p = p,
    chiCrit = qchisq(alpha, df,ncp = 0, lower.tail = F),
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

  print(out.table, row.names = F, right = F)

  semPower.showPlot(chiCrit = object$chiCrit, ncp = object$ncp, df = object$df)
  
}

