#' semPower.aPriori
#'
#' Determine required sample size given alpha, beta/power, df, and effect
#'
#' @param effect effect size specifying the discrepancy between H0 and H1 (a list for multiple group models; a vector of length 2 for effect-size differences)
#' @param effect.measure type of effect, one of "F0", "RMSEA", "Mc", "GFI", AGFI"
#' @param alpha alpha error
#' @param beta beta error; set either beta or power
#' @param power power (1-beta); set either beta or power
#' @param N a list of sample weights for multiple group power analyses, e.g. list(1,2) to make the second group twice as large as the first one
#' @param df the model degrees of freedom
#' @param p the number of observed variables, required for effect.measure = "GFI" and "AGFI"
#' @param SigmaHat model implied covariance matrix (a list for multiple group models). Use in conjunction with Sigma to define effect and effect.measure.
#' @param Sigma population covariance matrix (a list for multiple group models). Use in conjunction with SigmaHat to define effect and effect.measure.
#' @param muHat model implied mean vector
#' @param mu observed (or population) mean vector
#' @param simulatedPower whether to perform a simulated (TRUE) (rather than analytical, FALSE) power analysis.
#' @param modelH0 for simulated power: lavaan model string defining the (incorrect) analysis model.
#' @param modelH1 for simulated power: lavaan model string defining the comparison model. If omitted, the saturated model is the comparison model.
#' @param nReplications for simulated power: number of random samples drawn.
#' @param minConvergenceRate for simulated power: the minimum convergence rate required
#' @param lavOptions for simulated power: a list of additional options passed to lavaan, e.g., list(estimator = 'mlm') to request robust ML estimation
#' @param seed for simulated power
#' @return list
#' @examples
#' \dontrun{
#' power <- semPower.aPriori(effect = .05, effect.measure = "RMSEA", alpha = .05, beta = .05, df = 200)
#' summary(power)
#' power <- semPower.aPriori(effect = .15, effect.measure = "F0", alpha = .05, power = .80, df = 100)
#' summary(power)
#' power <- semPower.aPriori(effect = list(.05, .10), effect.measure = "F0", alpha = .05, 
#'                           power = .80, N = list(1, 1), df = 100)
#' summary(power)
#' power <- semPower.aPriori(alpha = .01, beta = .05, df = 5, 
#'                           SigmaHat = diag(4), Sigma = cov(matrix(rnorm(4*1000),  ncol=4)))
#' summary(power)
#' }
#' @importFrom stats qchisq pchisq optim
#' @export
semPower.aPriori <- function(effect = NULL, effect.measure = NULL,
                             alpha, beta = NULL, power = NULL,
                             N = NULL, df = NULL, p = NULL,
                             SigmaHat = NULL, Sigma = NULL, muHat = NULL, mu = NULL,
                             simulatedPower = FALSE, 
                             modelH0 = NULL, modelH1 = NULL,
                             nReplications = 250, minConvergenceRate = .5, lavOptions = NULL, 
                             seed = NULL,
                             ...){

  # validate input and do some preparations
  pp <- powerPrepare('a-priori', effect = effect, effect.measure = effect.measure,
                     alpha = alpha, beta = beta, power = power, abratio = NULL,
                     N = N, df = df, p = p,
                     SigmaHat = SigmaHat, Sigma = Sigma, muHat = muHat, mu = mu,
                     simulatedPower = simulatedPower, 
                     modelH0 = modelH0, modelH1 = modelH1,
                     nReplications = nReplications, minConvergenceRate = minConvergenceRate,
                     lavOptions = lavOptions)

  if(!is.null(beta)){
    desiredBeta <- beta
    desiredPower <- 1 - desiredBeta
  }else{
    desiredPower <- power
    desiredBeta <- 1 - power
  }
  logBetaTarget <- log(desiredBeta)

  weights <- 1
  if(!is.null(pp[['N']]) && length(pp[['N']]) > 1){
    weights <- unlist(pp[['N']]) / sum(unlist(pp[['N']]))
  }
  
  # analytical approach
  if(!simulatedPower){
    
    df <- pp[['df']]
    fmin <- pp[['fmin']]
    fmin.g <- pp[['fmin.g']]
    nrep <- NULL
    
    critChi <- qchisq(alpha, df = df, ncp = 0, lower.tail = FALSE)
    
    # make a reasonable guess about required sample size
    exponent <- -floor(log10(fmin)) + 1
    startN <- 5*10^(exponent)
    
    bPrecisionWarning <- (startN < 5)  # skip estm for N < 5, but take N = 10 as effective minimum 
    
    if(!bPrecisionWarning){
      
      chiCritOptim <- optim(par = c(startN), fn = getBetadiff,
                            critChi = critChi, logBetaTarget = logBetaTarget, fmin = unlist(fmin.g), df = df, weights = weights,
                            method = 'Nelder-Mead', control = list(warn.1d.NelderMead = FALSE))
      
      requiredN <- sum(ceiling(weights*chiCritOptim$par))  
      
      # even N = 10 achieves or exceeds desired power
      if(requiredN < 10){
        requiredN <- 10
        bPrecisionWarning <- TRUE
      }
      
    }else{
      # even N = 10 achieves or exceeds desired power
      requiredN <- 10
    }
    
    # N by group
    requiredN.g <- ceiling(weights * requiredN)

    impliedNCP <- getNCP(fmin.g, requiredN.g)
    impliedBeta <- pchisq(critChi, df, impliedNCP)
    impliedPower <- pchisq(critChi, df, impliedNCP, lower.tail = FALSE)
    
  # simulated power  
  }else{
    if(!is.null(seed)) set.seed(seed)
    iterationCounter <<- 1

    # determine starting N using analytical power
    ap <- semPower.powerLav(type = 'a-priori', 
                            alpha = alpha, beta = beta, power = power,
                            modelH0 = modelH0, modelH1 = modelH1,
                            Sigma = Sigma, mu = mu)
    startN <- ceiling(.95 * ap[['power']][['requiredN']]) # lets start a bit lower

    
    # for simulated power, we refuse to do anything
    bPrecisionWarning <- (ap[['power']][['requiredN']] < ncol(Sigma))  
    if(bPrecisionWarning) stop("Required N is smaller than the number of variables. Simulated power will not work well in this case because of very high nonconvergence rates.")

    # we need a pretty high tolerance because of sampling error: it doesn't make sense 
    # to suggest high accuracy when there it is in fact quite limited
    # one option would be to define tolerance in terms of expected sampling error, e.g., 
    # pchisq(chiCrit +/- sqrt(df), df, ncp, lower.tail = T), but this leads to wide margins.
    # the option below increases tolerance with decreasing nrep and beta, but is more restrictive.
    tolerance <- logBetaTarget^2 * 2 / nReplications
    chiCritOptim <- optim(par = c(N = startN), fn = getBetadiff,
                          logBetaTarget = logBetaTarget,
                          modelH0 = modelH0, modelH1 = modelH1,
                          Sigma = Sigma, mu = mu,  
                          alpha = alpha,
                          nReplications = nReplications, minConvergenceRate = minConvergenceRate,
                          lavOptions = lavOptions, simulatedPower = TRUE,
                          returnF = FALSE,
                          method = 'Nelder-Mead', 
                          control = list(warn.1d.NelderMead = FALSE, abstol = tolerance)
                          )
    
    if(chiCritOptim$convergence != 0) warning('A priori power analyses did not converge, results may be inaccurate.')    
    
    requiredN <- ceiling(chiCritOptim$par)  
    requiredN.g <- ceiling(chiCritOptim$par)  
    # TODO add multigroup support
    # requiredN <- sum(ceiling(weights*chiCritOptim$par))  
    # requiredN.g <- ceiling(weights * requiredN)
    
    # now call simulate with final N again to get all relevant parameters
    sim <- simulate(modelH0 = modelH0, modelH1 = modelH1,
                    Sigma = Sigma, mu = mu,
                    alpha = alpha, N = requiredN,
                    nReplications = nReplications, minConvergenceRate = minConvergenceRate,
                    lavOptions = lavOptions)
    
    nrep <- sim[['nrep']]
    df <- sim[['df']]
    fmin <- sim[['medianF']]
    fmin.g <- sim[['medianF']]   ## TODO add multigroup support
    
    critChi <- qchisq(alpha, df = df, ncp = 0, lower.tail = FALSE)
    
    impliedNCP <- getNCP(fmin.g, requiredN.g)
    impliedBeta <- 1 - sim[['ePower']]
    impliedPower <- sim[['ePower']]

  }
  
  # need to compute this after having determined Ns, because some indices rely on sample weights in multigroup case
  fit <- getIndices.F(fmin, df, p, SigmaHat, Sigma, muHat, mu, requiredN.g)
  
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
    effect = pp[['effect']],
    effect.measure = pp[['effect.measure']],
    requiredN = requiredN,
    requiredN.g = requiredN.g,
    df = df,
    p = pp[['p']],
    chiCrit = critChi,
    rmsea = fit[['rmsea']],
    mc = fit[['mc']],
    gfi = fit[['gfi']],
    agfi = fit[['agfi']],
    srmr = fit[['srmr']],
    cfi = fit[['cfi']],
    bPrecisionWarning = bPrecisionWarning,
    simulated = simulatedPower,
    nrep = nrep
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
#' @param simulatedPower whether to perform a simulated (TRUE) (rather than analytical, FALSE) power analysis.
#' @param ... other parameter passed to simulate()
#' @return squared difference requested and achieved beta on a log scale
#' @importFrom stats pchisq
getBetadiff <- function(cN, critChi, logBetaTarget, fmin, df, weights = NULL, 
                        simulatedPower = FALSE, ...){
  diff <- .Machine$integer.max
  
  if(!simulatedPower){

    if(cN < 5){
      # avoid NA in pchisq; nelder-mead can handle NA diff
      diff <- NA
    }else{
      cNCP <- sum(fmin * ((weights * cN) - 1) )
      cLogBeta <- pchisq(critChi, df, cNCP, log.p = TRUE)
      diff <- (logBetaTarget - cLogBeta)^2
    }
    
  }else{
    
    iterationCounter <<- iterationCounter + 1 
    ePower <- simulate(N = round(cN), ...)
    if(ePower == 1) ePower <- 1 - 1e-5
    diff <- (logBetaTarget - log(1 - ePower))^2
    print(paste("Iteration", iterationCounter, ":", cN, (1 - ePower), diff))
    
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

  if(object[['simulated']]){
    cat(paste("\n Simulated power based on", object[['nrep']], "successful replications.\n Note that simulated a-priori power analyses are only approximate,\n unless the number of replications is large.\n"))
  }
  
  if(object[['bPrecisionWarning']])
    cat("\n\n NOTE: Power is higher than requested even for a sample size < 10.\n\n")

  print(out.table, row.names = FALSE, right = FALSE)
  
  semPower.showPlot(chiCrit = object[['chiCrit']], ncp = object[['ncp']], df = object[['df']])
  
}



