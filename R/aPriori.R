#' semPower.aPriori
#'
#' Determine required sample size given alpha, beta/power, df, and effect
#'
#' @param effect effect size specifying the discrepancy between H0 and H1
#' @param effect.measure type of effect, one of "F0", "RMSEA", "Mc", "GFI", AGFI"
#' @param alpha alpha error
#' @param beta beta error; set either beta or power
#' @param power power (1-beta); set either beta or power
#' @param N a list of sample weights for multiple group power analyses, e.g. list(1,2) to make the second group twice as large as the first one 
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
#' power <- semPower.aPriori(alpha = .01, beta = .05, df = 5, 
#'                           SigmaHat = diag(4), Sigma = cov(matrix(rnorm(4*1000),  ncol=4)))
#' power
#' }
#' @importFrom stats qchisq pchisq optim
#' @export
semPower.aPriori <- function(effect = NULL, effect.measure = NULL,
                             alpha, beta = NULL, power = NULL,
                             N = NULL, df, p = NULL,
                             SigmaHat = NULL, Sigma = NULL){

  if(!is.null(effect.measure)) effect.measure <- toupper(effect.measure)

  # convert vectors to lists
  if(!is.list(effect) && length(effect) > 1) effect <- as.list(effect)
  if(!is.list(N) && length(N) > 1) N <- as.list(N)  # sample weights
  
  validateInput('a-priori', effect = effect, effect.measure = effect.measure,
                alpha = alpha, beta = beta, power = power, abratio = NULL,
                N = N, df = df, p = p,
                SigmaHat = SigmaHat, Sigma = Sigma)

  if(!is.null(SigmaHat)){ # sufficient to check for on NULL matrix; primary validity check is in validateInput
    effect.measure <- 'F0'
    p <- ifelse(is.list(SigmaHat), ncol(SigmaHat[[1]]), ncol(SigmaHat))
  }
  
  # make sure N/effects have the same length
  if((is.list(effect) || is.list(SigmaHat)) && length(N) == 1){
    N <- as.list(rep(N, ifelse(is.null(SigmaHat), length(effect), length(SigmaHat))))
  }
  if(is.null(SigmaHat) && is.list(N) && length(effect) == 1){
    effect <- as.list(rep(effect, length(N)))
  }
  ngroups <- ifelse(is.null(N), 1, length(N))
  
  # obsolete, single group case only
  # fmin <- getF(effect, effect.measure, df, p, SigmaHat, Sigma)
  # fit <- getIndices.F(fmin, df, p, SigmaHat, Sigma)
  
  if(!is.null(effect)){
    fmin.g <- sapply(effect, FUN = getF, effect.measure = effect.measure, df = df, p = p)
  }
  if(!is.null(SigmaHat)){
    if(is.list(Sigma)){
      fmin.g <- sapply(seq_along(SigmaHat), FUN = function(x) {getF.Sigma(SigmaHat = SigmaHat[[x]], S = Sigma[[x]]) })
    }else{
      fmin.g <- getF.Sigma(SigmaHat = SigmaHat, S = Sigma)
    }
  }
  
  fmin <- sum(fmin.g)
  
  if(!is.null(beta)){
    desiredBeta <- beta
    desiredPower <- 1 - desiredBeta
  }else{
    desiredPower <- power
    desiredBeta <- 1 - power
  }
  logBetaTarget <- log(desiredBeta)
  critChi <- qchisq(alpha, df=df, ncp = 0, lower.tail = F)

  # make a reasonable guess about required sample size
  exponent <- -floor(log10(fmin))+1
  startN <- 5*10^(exponent)

  bPrecisionWarning <- (startN <= 10)

  if(!bPrecisionWarning){

    weights <- 1
    if(!is.null(N) && length(N) > 1){
      weights <- unlist(N)/sum(unlist(N))
    }

    chiCritOptim <- optim(par = c(startN), fn = getBetadiff,
            critChi=critChi, logBetaTarget=logBetaTarget, fmin=unlist(fmin.g), df=df, weights=weights,
            method='Nelder-Mead', control = list(warn.1d.NelderMead=F))

    requiredN <- sum(ceiling(weights*chiCritOptim$par))  

    # even N = 10 achieves or exceeds desired power
    if(requiredN < 10){
      requiredN <- 10
      bPrecisionWarning = T
    }

  }else{
    # even N = 10 achieves or exceeds desired power
    requiredN <- 10
  }
  
  # N by group
  requiredN.g <- ceiling(weights*requiredN)
  
  # need to compute this after having determined Ns, because some indices rely on sample weights in multigroup case
  fit <- getIndices.F(fmin, df, p, SigmaHat, Sigma, requiredN.g)
  
  impliedNCP <- getNCP(fmin.g, requiredN.g)
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
    requiredN.g = requiredN.g,
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
#' @param weights sample weights for multiple group models 
#' @return squared difference requested and achieved beta on a log scale
#' @importFrom stats pchisq
getBetadiff <- function(cN, critChi, logBetaTarget, fmin, df, weights = NULL){
  diff <- .Machine$integer.max

  if(cN < 5){
    # avoid NA in pchisq; nelder-mead can handle NA diff
    diff <- NA
    
  }else{
    cNCP <- sum(fmin * ((weights *cN) - 1) )
    
    cLogBeta <- pchisq(critChi, df, cNCP, log.p = T)
    
    diff <- (logBetaTarget - cLogBeta)^2
    
  }
  
  diff
}


#' summary.semPower.aPriori
#'
#' provide summary of a-priori power analyses
#' @param object result object from semPower.aPriori
#' @param ... other
#' @export
summary.semPower.aPriori <- function(object, ...){

  out.table <- getFormattedResults('a-priori', object)

  cat("\n semPower: A-priori power analysis\n")

  if(object$bPrecisionWarning)
    cat("\n\n NOTE: Power is higher than requested even for a sample size < 10.\n\n")

  print(out.table, row.names = F, right = F)
  
  semPower.showPlot(chiCrit = object$chiCrit, ncp = object$impliedNCP, df = object$df)

}



