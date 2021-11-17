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
                             SigmaHat = NULL, Sigma = NULL){

  if(!is.null(effect.measure)) effect.measure <- toupper(effect.measure)
  
  # convert vectors to lists
  if(!is.list(N) && length(N) > 1) N <- as.list(N)
  
  validateInput('post-hoc', effect = effect, effect.measure = effect.measure,
                alpha = alpha, beta = NULL, power = NULL, abratio = NULL,
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

  if(!is.null(effect)){
    if(is.list(effect) || length(effect) == 1){
      fmin.g <- sapply(effect, FUN = getF, effect.measure = effect.measure, df = df, p = p)
    }else{
      # power for effect differences
      f1 <- getF(effect[1], effect.measure, df[1], p)
      f2 <- getF(effect[2], effect.measure, df[2], p)
      fdiff <- abs(f2 - f1)   # let's make order arbitrary  
      fmin.g <- rep(fdiff, length(N)) 
      df <- abs(df[2] - df[1]) # let's make order arbitrary  
    }
  }
  if(!is.null(SigmaHat)){
    if(is.list(Sigma)){
      fmin.g <- sapply(seq_along(SigmaHat), FUN = function(x) {getF.Sigma(SigmaHat = SigmaHat[[x]], S = Sigma[[x]])})
    }else{
      fmin.g <- getF.Sigma(SigmaHat = SigmaHat, S = Sigma)
    }
  }

  fmin <- sum(unlist(fmin.g) * unlist(N) / sum(unlist(N)))
  fit <- getIndices.F(fmin, df, p, SigmaHat, Sigma, N)
  ncp <- getNCP(fmin.g, N)

  beta <- pchisq(qchisq(alpha, df, lower.tail = FALSE), df, ncp = ncp)
  power <- pchisq(qchisq(alpha, df, lower.tail = FALSE), df, ncp = ncp, lower.tail = FALSE)
  impliedAbratio <- alpha / beta


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

