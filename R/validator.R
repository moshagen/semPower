
##########################  do input validation  #####################

#' validateInput
#'
#' Validates input for power calcuation function
#'
#'
#' @param power.type type of power analyses, one of "a-priori", post-hoc", "compromise", "powerplot.byN", "powerplot.byEffect"
#' @param effect effect size specifying the discrepancy between H0 and H1
#' @param effect.measure type of effect, one of "F0, "RMSEA", "Mc", "GFI", AGFI"
#' @param alpha alpha error
#' @param beta beta error
#' @param power power (1-beta)
#' @param abratio ratio alpha/beta
#' @param N the number of observations
#' @param df the model degrees of freedom
#' @param p the number of observed variables, required for effect.measure = "GFI" and "AGFI"
#' @param SigmaHat model implied covariance matrix
#' @param Sigma population covariance matrix
#' @param power.min for plotting: minimum power
#' @param power.max for plotting: maximum power
#' @param effect.min for plotting: minimum effect
#' @param effect.max for plotting: maximum effect
#' @param steps for plotting: number of sampled points
#' @param linewidth for plotting: linewidth
validateInput <- function(power.type = NULL, effect = NULL, effect.measure = NULL,
                          alpha = NULL, beta = NULL, power = NULL, abratio = NULL,
                          N = NULL, df = NULL, p = NULL,
                          SigmaHat = NULL, Sigma = NULL,
                          power.min = alpha, power.max = .999,
                          effect.min = NULL, effect.max = NULL,
                          steps = 50, linewidth = 1){

  known.effects.measures <- c("F0","RMSEA","Mc","GFI", "AGFI")

  # do input validation
  if(is.null(SigmaHat) && is.null(Sigma)){

    if(!effect.measure %in% known.effects.measures){
      stop(paste("effect measure is unknown, must be one of", paste(known.effects.measures, collapse = ", ")))
    }

    if(power.type != 'powerplot.byEffect')
      checkPositive(effect, effect.measure)

    if(effect.measure == "GFI" || effect.measure == "AGFI"){

      if(is.null(p)){
        stop("effect.measure GFI and AGFI require specification of p")
      }

      checkPositive(p, 'p')

    }
  }

  if((!is.null(SigmaHat) && is.null(Sigma)) || (is.null(Sigma) && !is.null(SigmaHat))){
    stop("Both SigmaHat and Sigma must be defined when effect is determined from covariance matrices")
  }
  if(!is.null(SigmaHat) && !is.null(Sigma)){
    if(ncol(SigmaHat) != ncol(Sigma) || nrow(SigmaHat) != nrow(Sigma))
      stop("Sigma and SigmaHat must be of equal size")
    if(!isSymmetric(Sigma) || !isSymmetric(SigmaHat))
      stop("Sigma and SigmaHat must be symmetric square matrices")
    if(ncol(SigmaHat) < 2)
      stop("Sigma and SigmaHat must be at least of dimension 2*2")
    if(sum(eigen(Sigma)$values < 0) > 0 || sum(eigen(SigmaHat)$values < 0) > 0)
      stop("Sigma and SigmaHat must be positive definite")
    if(!is.null(effect) || !is.null(effect.measure))
      warning("ignoring effect and effect.measure when Sigma and SigmaHat are set")
  }

  checkPositive(df, 'df')

  # specifics depending on type of power analyses

  if(power.type == "post-hoc" || power.type == "compromise"){
    checkPositive(N, 'N')
  }

  if(power.type == "a-priori" || power.type == "post-hoc"){
    checkBounded(alpha, 'alpha')
  }

  if(power.type == "a-priori"){
    if(is.null(beta) && is.null(power))
      stop("Need to define either beta or power in a-priori power analyis")
    if(!is.null(beta) && !is.null(power) && power != (1-beta))
      stop("Either set beta or set power, but not both.")
    if(!is.null(beta))
      checkBounded(beta, 'beta')
    if(!is.null(power))
      checkBounded(power, 'power')
  }

  if(power.type == "compromise"){
    checkPositive(abratio, 'abratio')
  }

  # specifics for power plots
  if(power.type == "powerplot.byN"){
    checkBounded(alpha, 'alpha')
    checkBounded(power.max, "power.max")
    checkBounded(power.min, "power.min")
    checkPositive(steps, "steps")
    checkPositive(linewidth, "linewidth")
  }
  
  if(power.type == "powerplot.byEffect"){
    checkPositive(N, 'N')
    checkBounded(alpha, 'alpha')
    checkPositive(effect.min, "effect.min")
    checkPositive(effect.max, "effect.max")
    checkPositive(steps, "steps")
    checkPositive(linewidth, "linewidth")
  }
  
  
}


#' checkPositive
#'
#' checks whether x is defined and a positive number, stop otherwise
#' @param x x
#' @param message identifier for x
checkPositive <- function(x, message){
  if(is.null(x) || is.na(x) || x <= 0){
    stop(paste(message," must be larger than zero"))
  }
}

#' checkBounded
#'
#' checks whether x is defined and lies within the specified bound
#' @param x x
#' @param message identifier for x
#' @param bound the boundaries, array of size two
checkBounded <- function(x, message, bound = c(0,1)){
  if(is.null(x) || is.na(x) || x <= bound[1] || x >= bound[2]){
    stop(paste(message," must must lie within",bound[1],'and', bound[2]))
  }
}

#' checkPositiveDefinite
#'
#' checks whether x is positive definite
#' @param x x
#' @param message identifier for x
checkPositiveDefinite <- function(x, message){
  if(is.null(x))
    stop(paste(message, " may not be NULL"))
  if(!is.numeric(x))
    stop(paste(message, " must contain numeric elements only"))
  if(!isSymmetric(x))
    stop(paste(message, " must be a symmetric square matrix"))
  if(sum(eigen(x)$values < 0) > 0)
    stop(paste(message, " must be positive definite"))
}
