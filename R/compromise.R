
#' sempower.compromise
#'
#' Performs a compromise power analysis, i.e. determines the critical chi-square along with the implied alpha and beta, given a specified alpha/beta ratio, effect, N, and df
#'
#' @param effect effect size specifying the discrepancy between H0 and H1  (a list for multiple group models; a vector of length 2 for effect-size differences)
#' @param effect.measure type of effect, one of "F0","RMSEA", "Mc", "GFI", AGFI"
#' @param abratio the ratio of alpha to beta
#' @param N the number of observations  (a list for multiple group models)
#' @param df the model degrees of freedom
#' @param p the number of observed variables, required for effect.measure = "GFI" and "AGFI"
#' @param SigmaHat model implied covariance matrix (a list for multiple group models). Use in conjunction with Sigma to define effect and effect.measure.  
#' @param Sigma population covariance matrix (a list for multiple group models). Use in conjunction with SigmaHat to define effect and effect.measure.
#' @param muHat model implied mean vector
#' @param mu observed (or population) mean vector
#' @param simulatedPower whether to perform a simulated (TRUE) (rather than analytical, FALSE) power analysis. Not supported for compromise power analysis.
#' @return list
#' @examples
#' \dontrun{
#' cp.ph <- semPower.compromise(effect = .08, effect.measure = "RMSEA", abratio = 1, N = 250, df = 200)
#' summary(cp.ph)
#' }
#' @importFrom stats qchisq pchisq optim
#' @export
semPower.compromise  <- function(effect = NULL, effect.measure = NULL,
                                 abratio = 1,
                                 N, df = NULL, p = NULL,
                                 SigmaHat = NULL, Sigma = NULL, muHat = NULL, mu = NULL,
                                 simulatedPower = FALSE, ...){

  if(simulatedPower) stop('Simulated power is not available for compromise power analysis, because this would require a vast (infeasible) number of simulation runs to yield reliable results.')

  # validate input and do some preparations
  pp <- powerPrepare('compromise', effect = effect, effect.measure = effect.measure,
                     alpha = NULL, beta = NULL, power = NULL, abratio = abratio,
                     N = N, df = df, p = p,
                     SigmaHat = SigmaHat, Sigma = Sigma, muHat = muHat, mu = mu, 
                     simulatedPower = simulatedPower)

  df <- pp[['df']]
  
  fit <- getIndices.F(pp[['fmin']], df, pp[['p']], SigmaHat, Sigma, muHat, mu, pp[['N']])
  
  ncp <- getNCP(pp[['fmin.g']], pp[['N']])
  log.abratio <- log(abratio)

  if(ncp >= 1e5)
    warning('NCP is larger than 100000, this is going to cause trouble.')


  # determine max/min chi for valid alpha/beta prob
  max <- min <- NA
  # central chi always gives result up to 1e-320
  max <- qchisq(log(1e-320), df, lower.tail = FALSE, log.p = TRUE)

  # non-central chi accuracy is usually lower, depending on df and ncp
  pmin <- -Inf
  testp <- 1e-320
  while(is.infinite(pmin)){
    testp <- testp * 10
    testv <- max(log(1e-320), (log.abratio + log(testp)))
    min <- qchisq(testv, df, ncp, log.p = TRUE)
    pmin <- pchisq(min, df, ncp, log.p = TRUE) # beta
  }

  # cannot determine critChi when implied errors are too small
  bPrecisionWarning <- (min > max)

  if(!bPrecisionWarning){
    # rough estm
    start <- df + ncp / 3
    chiCritOptim <- optim(par = c(start), fn = getErrorDiff,
                          df = df, ncp = ncp, log.abratio = log.abratio,
                          method = 'L-BFGS-B', lower = min, upper = max)

    chiCrit <- chiCritOptim$par
    impliedAlpha <- pchisq(chiCrit, df, lower.tail = FALSE)
    impliedBeta <- pchisq(chiCrit, df, ncp)
    impliedAbratio <- impliedAlpha / impliedBeta
    impliedPower <- pchisq(chiCrit, df, ncp, lower.tail = FALSE)

  }else{
    # this is overriden later
    chiCrit <- 0
    impliedAlpha <- 0
    impliedBeta <- 0
    impliedAbratio <- 1
    impliedPower <- 1
  }


  result <- list(
    type = "post-hoc-compromise",
    desiredAbratio = abratio,
    chiCrit = chiCrit,
    impliedAlpha = impliedAlpha,
    impliedBeta = impliedBeta,
    impliedAbratio = impliedAbratio,
    impliedPower = impliedPower,
    ncp = ncp,
    fmin = pp[['fmin']],
    effect = pp[['effect']],
    effect.measure = pp[['effect.measure']],
    N = pp[['N']],
    df = df,
    p = pp[['p']],
    rmsea = fit[['rmsea']],
    mc = fit[['mc']],
    gfi = fit[['gfi']],
    agfi = fit[['agfi']],
    srmr = fit[['srmr']],
    cfi = fit[['cfi']],
    max = max,
    min = min,
    bPrecisionWarning = bPrecisionWarning
  )

  class(result) <- "semPower.compromise"

  result
}


#' getErrorDiff
#'
#' Determine the squared log-difference between alpha and beta error given a certain chi-square value from central chi-square(df) and a non-central chi-square(df, ncp) distribution.
#'
#' @param critChiSquare evaluated chi-squared value
#' @param df the model degrees of freedom
#' @param ncp the non-centrality parameter
#' @param log.abratio log(alpha/beta)
#' @return squared difference between alpha and beta on a log scale
#' @importFrom stats pchisq
getErrorDiff <- function(critChiSquare, df, ncp, log.abratio){

  alpha <- pchisq(critChiSquare, df, lower.tail = FALSE, log.p = TRUE)
  beta <- pchisq(critChiSquare, df, ncp, log.p = TRUE)

  if(is.infinite(beta) || is.infinite(alpha)){

    warning('Alpha or beta are too small.')
    diff <- 0

  }else{

    diff <- (alpha - (log.abratio + beta))^2     # note log scale

  }

  diff
}


#' summary.sempower.compromise
#'
#' provide summary of compromise post-hoc power analyses
#' @param object result object from semPower.compromise
#' @param ... other
#' @export
summary.semPower.compromise <- function(object, ...){

  out.table <- getFormattedResults('compromise', object)

  cat("\n semPower: Compromise power analysis\n")

  if(object[['bPrecisionWarning']])
    cat("\n\n WARNING: Alpha and/or Beta are smaller than 1e-240. Cannot determine critical Chi-Square exactly due to machine precision.")

  print(out.table, row.names = FALSE, right = FALSE)

  if(!object[['bPrecisionWarning']])
    semPower.showPlot(chiCrit = object[['chiCrit']], ncp = object[['ncp']], df = object[['df']])
  
}


