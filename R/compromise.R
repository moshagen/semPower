
#' sempower.compromise
#'
#' Performs a compromise power analysis, i.e. determines the critical chi-square along with the implied alpha and beta, given a specified alpha/beta ratio, effect, N, and df
#'
#' @param effect effect size specifying the discrepancy between H0 and H1  (a list for multiple group models)
#' @param effect.measure type of effect, one of "F0","RMSEA", "Mc", "GFI", AGFI"
#' @param abratio the ratio of alpha to beta
#' @param N the number of observations  (a list for multiple group models)
#' @param df the model degrees of freedom
#' @param p the number of observed variables, required for effect.measure = "GammaHat", "GFI",  and "AGFI"
#' @param SigmaHat model implied covariance matrix (a list for multiple group models). Use in conjuntion with Sigma to define effect and effect.measure.  
#' @param Sigma population covariance matrix (a list for multiple group models). Use in conjuntion with SigmaHat to define effect and effect.measure.
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
                                 N, df, p = NULL,
                                 SigmaHat = NULL, Sigma = NULL){

  if(!is.null(effect.measure)) effect.measure <- toupper(effect.measure)
  
  # convert vectors to lists
  if(!is.list(effect) && length(effect) > 1) effect <- as.list(effect)
  if(!is.list(N) && length(N) > 1) N <- as.list(N) 
  
  validateInput('compromise', effect = effect, effect.measure = effect.measure,
                alpha = NULL, beta = NULL, power = NULL, abratio = abratio,
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
  
  fmin <- sum(unlist(fmin.g) * unlist(N) / sum(unlist(N)))
  fit <- getIndices.F(fmin, df, p, SigmaHat, Sigma, N)
  
  ncp <- getNCP(fmin.g, N)
  log.abratio <- log(abratio)

  if(ncp >= 1e5)
    warning('NCP is larger than 100000, this is going to cause trouble.')


  # determine max/min chi for valid alpha/beta prob
  max <- min <- NA
  # central chi always gives reusult up to 1e-320
  max <- qchisq(log(1e-320), df, lower.tail = F, log.p = T)

  # non-central chi accuracy is usually lower, depending on df and ncp
  pmin <- -Inf
  testp <- 1e-320
  while(is.infinite(pmin)){
    testp <- testp * 10
    testv <- max(log(1e-320), (log.abratio + log(testp)))
    min <- qchisq(testv, df, ncp, log.p = T)
    pmin <- pchisq(min, df, ncp, log.p = T) # beta
  }

  # cannot determine critchi when implied errors are too small
  bPrecisionWarning <- (min > max)

  if(!bPrecisionWarning){
    # rough estm
    start <- df + ncp/3
    chiCritOptim <- optim(par = c(start), fn = getErrorDiff,
                          df=df, ncp=ncp, log.abratio = log.abratio,
                          method='L-BFGS-B', lower=min, upper=max)

    chiCrit <- chiCritOptim$par
    impliedAlpha <- pchisq(chiCrit, df, lower.tail = F)
    impliedBeta <- pchisq(chiCrit, df, ncp)
    impliedAbratio <- impliedAlpha/impliedBeta
    impliedPower <- pchisq(chiCrit, df, ncp, lower.tail = F)

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
    fmin = fmin,
    effect = effect,
    effect.measure = effect.measure,
    N = N,
    df = df,
    p = p,
    rmsea = fit$rmsea,
    mc = fit$mc,
    gfi = fit$gfi,
    agfi = fit$agfi,
    srmr = fit$srmr,
    cfi = fit$cfi,
    max = max,
    min = min,
    bPrecisionWarning = bPrecisionWarning
  )

  class(result) <- "semPower.compromise"

  result
}


#' getErrorDiff
#'
#' determine the squared log-difference between alpha and beta error given a certain chi-square value from central chi-square(df) and a non-central chi-square(df, ncp) distribution.
#'
#' @param critChiSquare evaluated chi-squared value
#' @param df the model degrees of freedom
#' @param ncp the non-centrality parameter
#' @param log.abratio log(alpha/beta)
#' @return squared difference between alpha and beta on a log scale
#' @importFrom stats pchisq
getErrorDiff <- function(critChiSquare, df, ncp, log.abratio){

  alpha <- pchisq(critChiSquare, df, lower.tail = F, log.p = T)
  beta <- pchisq(critChiSquare, df, ncp, log.p = T)

  if(is.infinite(beta) || is.infinite(alpha)){

    warning('alpha or beta is too small')
    diff <- 0

  }else{

    diff <- (alpha - (log.abratio+beta))^2     # note log scale

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

  if(object$bPrecisionWarning)
    cat("\n\n WARNING: Alpha and/or Beta are smaller than 1e-240. Cannot determine critical Chi-Square exactly due to machine precision.")

  print(out.table, row.names = F, right = F)

  if(!object$bPrecisionWarning)
    semPower.showPlot(chiCrit = object$chiCrit, ncp = object$ncp, df = object$df)
  
}


