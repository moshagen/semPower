#' semPower.aPriori
#'
#' Determine required sample size given alpha, beta/power, df, and effect
#'
#' @param effect effect size specifying the discrepancy between H0 and H1
#' @param effect.measure type of effect, one of "F0", "RMSEA", "Mc", "GFI", AGFI"
#' @param alpha alpha error
#' @param beta beta error; set either beta or power
#' @param beta power (1-beta); set either beta or power
#' @param N the number of observations
#' @param df the model degrees of freedom
#' @param p the number of observed variables, required for effect.measure = "GFI" and "AGFI"
#' @param SigmaHat model implied covariance matrix. Use in conjuntion with Sigma to define effect and effect.measure. 
#' @param Sigma population covariance matrix. Use in conjuntion with SigmaHat to define effect and effect.measure.
#' @return list
#' @examples
#' \dontrun{
#' power <- semPower.aPriori(effect = .05, effect.measure = "RMSEA", alpha = .05, beta = .05, df = 200)
#' power
#' power <- semPower.aPriori(effect = .15, effect.measure = "F0", alpha = .05, power = .80, df = 100)
#' power
#' power <- semPower.aPriori(alpha = .01, beta = .05, df = 5, SigmaHat = diag(4), Sigma = cov(matrix(rnorm(4*1000),  ncol=4)))
#' power
#' }
#' @export
semPower.aPriori <- function(effect = NULL, effect.measure = NULL,
                             alpha = .05, beta = NULL, power = NULL,
                             df = NULL, p = NULL,
                             SigmaHat = NULL, Sigma = NULL){

  validateInput('a-priori', effect = effect, effect.measure = effect.measure,
                alpha = alpha, beta = beta, power = power, abratio = NULL,
                N = NULL, df = df, p = p,
                SigmaHat = SigmaHat, Sigma = Sigma)

  if(!is.null(SigmaHat)){ # sufficient to check for on NULL matrix; primary validity check is in validateInput
    effect.measure <- 'F0'
    p <- ncol(SigmaHat)
  }

  fmin <- getF(effect, effect.measure, df, p, SigmaHat, Sigma)
  fit <- getIndices.F(fmin, df, p, SigmaHat, Sigma)

  if(!is.null(beta)){
    desiredBeta <- beta
    desiredPower <- 1 - desiredBeta
  }else{
    desiredPower <- power
    desiredBeta <- 1 - power
  }
  logBetaTarget <- log(desiredBeta)
  critChi <- qchisq(alpha, df=df, ncp = 0, lower.tail = F)

  #TMP
  #  df <- 200
  #  alpha <- .05
  #  critChi <- qchisq(alpha, df=df, ncp = 0, lower.tail = F)
  #  logBetaTarget <- log(.05)
  #  fmin <- getF.RMSEA(.9, df)
  #TMP

  # make a reasonable guess about required sample size
  exponent <- -floor(log10(fmin))+1
  startN <- 5*10^(exponent)

  bPrecisionWarning <- (startN <= 10)

  if(!bPrecisionWarning){

    chiCritOptim <- suppressWarnings( # we dont want to hear that Nelder-Mead doesn't like unidimensional optimization
      optim(par = c(startN), fn = getBetadiff,
            critChi=critChi, logBetaTarget=logBetaTarget, fmin=fmin, df=df,
            method='Nelder-Mead')
    )

    requiredN <- ceiling(chiCritOptim$par)

    # even N = 10 achieves or exceeds desired power
    if(requiredN < 10){
      requiredN <- 10
      bPrecisionWarning = T
    }

  }else{
    # even N = 10 achieves or exceeds desired power
    requiredN <- 10
  }


  impliedNCP <- getNCP(fmin, requiredN)
  impliedBeta <- pchisq(critChi, df, impliedNCP)
  impliedPower <- pchisq(critChi, df, impliedNCP, lower.tail = F)
  impliedAbratio <- alpha / impliedBeta

  result <- list(
    type = "a-priori",
    alpha = alpha,
    desiredBeta = desiredBeta,
    desiredPower = desiredPower,
    impliedBeta = impliedBeta,
    impliedPower = impliedPower,
    impliedAbratio = impliedAbratio,
    impliedNCP = impliedNCP,
    fmin = fmin,
    effect = effect,
    effect.measure = effect.measure,
    requiredN = requiredN,
    df = df,
    p = p,
    chiCrit = critChi,
    rmsea = fit$rmsea,
    mc = fit$mc,
    gfi = fit$gfi,
    agfi = fit$agfi,
    srmr = fit$srmr,
    cfi = fit$cfi,
    bPrecisionWarning = bPrecisionWarning
  )

  class(result) <- "semPower.aPriori"
  result

}

#' getBetadiff
#'
#' get squared difference between requested and achieved beta on a logscale
#'
#' @param cN current N
#' @param critChi critical chi-square associated with chosen alpha error
#' @param logBetaTarget log(desired beta)
#' @param fmin minimum of the ML fit function
#' @param df the model degrees of freedom
#' @return squared difference requested and achieved beta on a log scale
#' @examples
#' \dontrun{
#' betaDiff <- getBetadiff(cN = 300, critChi = 350, logBetaTarget = log(.20), fmin = 0.50, df = 200)
#' betaDiff
#' }
#' @export
getBetadiff <- function(cN, critChi, logBetaTarget, fmin, df){
  diff <- .Machine$integer.max

  if(cN < 5){
    # avoid NA in pchisq; nelder-mead can handle NA diff
    diff <- NA

  }else{

    cNCP <- fmin * (cN - 1)
    cLogBeta <- pchisq(critChi, df, cNCP, log.p = T)

    diff <- (logBetaTarget - cLogBeta)^2

    # tmp
    #print(paste(cN, cLogBeta, logBetaTarget, diff))

  }

  diff
}


#' summary.semPower.aPriori
#'
#' provide summary of a-priori power analyses
#' @param result result object from semPower.aPriori
#' @export
summary.semPower.aPriori <- function(object, ...){

  out.table <- getFormattedResults('a-priori', object)

  cat("\n semPower: A-priori power analysis\n")

  if(object$bPrecisionWarning)
    cat("\n\n NOTE: Power is higher than requested even for a sample size < 10.\n\n")

  print(out.table, row.names = F, right = F)

}


############### UNIT TESTS

# pa.ap <- semPower.aPriori(effect = .05, effect.measure = "RMSEA", alpha = .05, beta = .05, df = 200, p=10)
# # summary(pa.ap)
# summary.semPower.aPriori(pa.ap)
#
# pa.ap <- semPower.aPriori(effect = .9, effect.measure = "RMSEA", alpha = .05, beta = .05, df = 200)
# # summary(pa.ap)
# summary.semPower.aPriori(pa.ap)
#
# pa.ap <- semPower.aPriori(effect = .00001, effect.measure = "RMSEA", alpha = .05, beta = .01, df = 200)
# # summary(pa.ap)
# summary.semPower.aPriori(pa.ap)




