#' powerPrepare
#'
#' Performs some preparations common to all types of power analyses.
#'
#' @param type type of power analysis
#' @param effect effect size specifying the discrepancy between H0 and H1 (a list for multiple group models; a vector of length 2 for effect-size differences)
#' @param effect.measure type of effect, one of `"F0"`, `"RMSEA"`, `"Mc"`, `"GFI"`, `"AGFI"`
#' @param alpha alpha error
#' @param beta beta error; set either beta or power
#' @param power power (=1 - beta); set either beta or power
#' @param abratio the ratio of alpha to beta
#' @param N the number of observations (a list for multiple group models)
#' @param df the model degrees of freedom 
#' @param p the number of observed variables, required for `effect.measure = "GFI"` and `effect.measure = "AGFI"`
#' @param SigmaHat model implied covariance matrix (a list for multiple group models). Use in conjunction with `Sigma` to define `effect` and `effect.measure`. 
#' @param Sigma observed (or population) covariance matrix (a list for multiple group models). Use in conjunction with `SigmaHat` to define `effect` and `effect.measure`.
#' @param muHat model implied mean vector
#' @param mu observed (or population) mean vector
#' @param simulatedPower whether to perform a simulated (`TRUE`) (rather than analytical, `FALSE`) power analysis.
#' @param modelH0 for simulated power: `lavaan` model string defining the (incorrect) analysis model.
#' @param nReplications for simulated power: number of random samples drawn.
#' @param minConvergenceRate for simulated power: the minimum convergence rate required
#' @param lavOptions for simulated power: a list of additional options passed to `lavaan`, e. g., `list(estimator = 'mlm')` to request robust ML estimation.
#' @return list
#' @importFrom utils installed.packages menu
powerPrepare <- function(type = NULL, 
                         effect = NULL, effect.measure = NULL,
                         alpha = NULL, beta = NULL, power = NULL, abratio = NULL,
                         N = NULL, df = NULL, p = NULL,
                         SigmaHat = NULL, Sigma = NULL, muHat = NULL, mu = NULL,
                         simulatedPower = FALSE, 
                         modelH0 = NULL, 
                         nReplications = NULL, minConvergenceRate = NULL, lavOptions = NULL){
  
  if(!is.null(effect.measure)) effect.measure <- toupper(effect.measure)
  
  if(!is.list(N) && length(N) > 1) N <- as.list(N)

  validateInput(type, effect = effect, effect.measure = effect.measure,
                alpha = alpha, beta = beta, power = power, abratio = abratio,
                N = N, df = df, p = p,
                SigmaHat = SigmaHat, Sigma = Sigma, muHat = muHat, mu = mu,
                simulatedPower = simulatedPower, 
                modelH0 = modelH0)

  # convert vectors to lists and vice versa
  if(is.list(N) && length(N) == 1) N <- N[[1]]
  if(is.list(effect) && length(effect) == 1) effect <- effect[[1]]
  if(is.list(SigmaHat) && length(SigmaHat) == 1) SigmaHat <- SigmaHat[[1]]
  if(is.list(Sigma) && length(Sigma) == 1) Sigma <- Sigma[[1]]
  if(is.list(muHat) && length(muHat) == 1) muHat <- muHat[[1]]
  if(is.list(mu) && length(mu) == 1) mu <- mu[[1]]
  
  
  if(!is.null(Sigma)){
    effect.measure <- 'F0'
    p <- ifelse(is.list(Sigma), ncol(Sigma[[1]]), ncol(Sigma))
  }
  
  # make sure N/effects have the same length
  if((is.list(effect) || is.list(Sigma)) && length(N) == 1){
    N <- as.list(rep(N, ifelse(is.null(Sigma), length(effect), length(Sigma))))
  }
  if(type == 'a-priori'){
    if(!is.null(Sigma) && !is.list(Sigma)){
      N <- 1 # single weight for single group model
    }else if(length(effect) == 1 || (length(effect) == 2 && is.null(N))){
      N <- 1 # single weight for single group model
    }
  }
  
  if(is.null(Sigma) && is.list(N) && length(effect) == 1){
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
      fmin.g <- sapply(seq_along(SigmaHat), FUN = function(x) {getF.Sigma(SigmaHat = SigmaHat[[x]], S = Sigma[[x]], muHat = muHat[[x]], mu = mu[[x]])})
    }else{
      fmin.g <- getF.Sigma(SigmaHat = SigmaHat, S = Sigma, muHat = muHat, mu = mu)
    }
  }
  
  if(!simulatedPower){
    fmin <- sum(unlist(fmin.g) * unlist(N) / sum(unlist(N)))
  }else{
    fmin <- fmin.g <- NULL
  }
  
  if(simulatedPower){
    if(!'lavaan' %in% rownames(installed.packages())) stop('This function depends on the lavaan package, so install lavaan first.')
  }
  
  list(
    effect.measure = effect.measure,
    SigmaHat = SigmaHat, 
    Sigma = Sigma, 
    muHat = muHat, 
    mu = mu,    
    N = N,
    p = p,
    effect = effect,
    fmin = fmin,
    fmin.g = fmin.g,
    df = df
  )
}

#' validateInput
#'
#' Validates input for power functions.
#'
#' @param power.type type of power analyses, one of `"a-priori"`, `"post-hoc"`, `"compromise"`, `"powerplot.byN"`, `"powerplot.byEffect"`
#' @param effect effect size specifying the discrepancy between H0 and H1
#' @param effect.measure type of effect, one of `"F0"`, `"RMSEA"`, `"Mc"`, `"GFI"`, `"AGFI"`
#' @param alpha alpha error
#' @param beta beta error
#' @param power power (= 1 - beta)
#' @param abratio ratio alpha/beta
#' @param N the number of observations
#' @param df the model degrees of freedom
#' @param p the number of observed variables, required for `effect.measure = "GFI"` and `effect.measure = "AGFI"`
#' @param SigmaHat model implied covariance matrix
#' @param Sigma observed (or population) covariance matrix
#' @param muHat model implied mean vector
#' @param mu observed (or population) mean vector
#' @param simulatedPower whether to perform a simulated (`TRUE`) (rather than analytical, `FALSE`) power analysis.
#' @param modelH0 for simulated power: `lavaan` model string defining the (incorrect) analysis model.
#' @param power.min for plotting: minimum power
#' @param power.max for plotting: maximum power
#' @param effect.min for plotting: minimum effect
#' @param effect.max for plotting: maximum effect
#' @param steps for plotting: number of sampled points
#' @param linewidth for plotting: linewidth
validateInput <- function(power.type = NULL, effect = NULL, effect.measure = NULL,
                          alpha = NULL, beta = NULL, power = NULL, abratio = NULL,
                          N = NULL, df = NULL, p = NULL,
                          SigmaHat = NULL, Sigma = NULL, muHat = NULL, mu = NULL,
                          simulatedPower = FALSE, modelH0 = NULL,
                          power.min = alpha, power.max = .999,
                          effect.min = NULL, effect.max = NULL,
                          steps = 50, linewidth = 1){

  known.effects.measures <- c("F0", "RMSEA", "MC", "GFI", "AGFI")
  
  if(!simulatedPower){
    
    if(is.null(df)) stop("df must be defined")
    lapply(df, checkPositive, message = 'df')
    
    # generic power analyses
    if(is.null(SigmaHat) && is.null(Sigma)){
      
      if(!effect.measure %in% known.effects.measures){
        stop(paste("Effect measure is unknown, must be one of", paste(known.effects.measures, collapse = ", ")))
      }
      
      # for effect-size differences, check matching length
      if(!is.null(effect) && !is.list(effect)){
        if(length(effect) > 2){
          stop("Power analyses with multiple groups requires specification of the effect size in each group provided as a list.")
        }else{
          if(length(effect) == 2 && length(df) != 2){
            stop("Power analyses for effect size differences requires specification of the df of the model pairs.")
          }
        }
      }
      
      # for multiple group analyses, check matching length
      if(is.list(N)){
        if(is.list(effect) && length(effect) == 1){
          warning("Only single effect size provided in multiple group power analyses, assuming equal effects in each group.")
        }else if(is.list(effect) && length(N) != length(effect)){
          stop("Power analyses with multiple groups requires specification of the effect size in each group.")
        }
      }
      if(is.list(effect) && power.type == "a-priori"){   # special messages for a priori, given weights need to be provided 
        if(is.null(N) || length(N) == 1){
          warning("No or only a single sample weight provided in multiple group power analyses, assuming equal weights.")        
        }else if(!is.null(N) && length(N) != length(effect)){
          stop("A priori power analyses with multiple groups requires specification of sample size weights for each group via the N argument")      
        } 
      }
      if(is.list(effect) && (power.type == "post-hoc" || power.type == "compromise")){
        if(!is.null(N) && length(N) == 1){
          warning("Only single sample size provided in multiple group power analyses, assuming equal sample sizes for each group.")
        }else if(!is.null(N) && length(N) != length(effect)){
          stop("Power analyses with multiple groups requires specification of sample sizes for each group")        
        }
      }
      
      if(power.type != 'powerplot.byEffect'){
        if(is.null(effect)) stop('Effect is not defined.')
        lapply(effect, checkPositive, message = effect.measure)
      }
      
      if(effect.measure == "GFI" || effect.measure == "AGFI"){
        if(is.null(p)){
          stop("effect.measure GFI and AGFI require specification of p")
        }
        checkPositive(p)
      }
    }
    
    # power analyses based on covariance matrices  
    if((!is.null(SigmaHat) && is.null(Sigma)) || (!is.null(Sigma) && is.null(SigmaHat))){
      stop("Both SigmaHat and Sigma must be defined when effect is determined from covariance matrices")
    }
    if((!is.null(unlist(muHat)) && is.null(unlist(mu))) || (!is.null(unlist(mu)) && is.null(unlist(muHat)))){
      stop("Both muHat and mu must be defined when effect is determined from covariance matrices and meanstructures.")
    }
    if((!is.null(unlist(muHat)) || !is.null(unlist(mu))) && (is.null(Sigma) || is.null(SigmaHat))){
      stop("When effect is to be determined from covariance matrices and meanstructures, Sigma, SigmaHat, mu, and muHat must be provided.")
    }
    
    if(!is.null(SigmaHat) && !is.null(Sigma)){
      
      if(!is.null(modelH0)) warning("ModelH0 is ignored when SigmaHat is defined.")
      if(!is.null(effect) || !is.null(effect.measure))
        warning("Ignoring effect and effect.measure when Sigma and SigmaHat are defined.")
      
      # convert to list in any case so validiation below always applies 
      if(!is.list(Sigma)) Sigma <- list(Sigma)
      if(!is.list(SigmaHat)) SigmaHat <- list(SigmaHat)
      
      if(any(sapply(Sigma, ncol) != sapply(SigmaHat, ncol)) || any(sapply(Sigma, nrow) != sapply(SigmaHat, nrow)))
        stop("Sigma and SigmaHat must be of equal size")
      if(any(!sapply(c(Sigma, SigmaHat), isSymmetric)))
        stop("Sigma and SigmaHat must be symmetric square matrices")
      if(any(sapply(c(Sigma, SigmaHat), ncol) < 2))
        stop("Sigma and SigmaHat must be at least of dimension 2*2")
      if(any(sapply(c(Sigma, SigmaHat), function(x){sum(eigen(x)$values < 0) > 0})))
        stop("Sigma and SigmaHat must be positive definite")
      
      # now remove list structure for length = 1
      if(is.list(SigmaHat) && length(SigmaHat) == 1) SigmaHat <- SigmaHat[[1]]
      if(is.list(Sigma) && length(Sigma) == 1) Sigma <- Sigma[[1]]
      
      # multigroup case
      if(is.list(SigmaHat) || is.list(Sigma)){
        if(length(SigmaHat) != length(Sigma)) stop("Multiple group power analyses require specification of SigmaHat and Sigma for each group.")
        if(is.null(N) && power.type != "a-priori") stop("Multiple group power analyses require specification of N for each group (representing weights in a priori power analysis).")
        if(is.null(N) && power.type == "a-priori") warning("No sample weights provided, assuming equal sample sizes in each group.")
        if(length(N) == 1){
          warning("Only single sample size provided in multiple group power analyses, assuming equal sample sizes (weights in a priori power analyses) for each group.")
        }else if(length(SigmaHat) != length(N)){
          stop("Multiple group power analyses require specification of N for each group.")
        } 
        
        if(!is.null(unlist(muHat)) && !is.null(unlist(mu))){
          if(!is.list(muHat)) muHat <- list(muHat)
          if(!is.list(mu)) mu <- list(mu)
          if(length(muHat) != length(mu)) stop("Multiple group power analyses require specification of mu and muHat for each group.")
          if(length(mu) != length(Sigma)) stop("Multiple group power analyses require specification of mu and muHat for each group.")
          if(any(sapply(c(mu, muHat), is.na)))
            stop("mu and muHat must not contain NA values.")
          if(any(sapply(mu, length) != sapply(muHat, length)))
            stop("mu and muHat must be of same size.")
          if(is.list(Sigma) && any(sapply(mu, length) != sapply(Sigma, ncol)))
            stop("mu and muHat must be of same length as the columns/rows of Sigma and SigmaHat.")
        }
      }
      
    }
    
    ### TODO remove code duplication by merging the following with the preceding block
    # simulated power specifics  
  }else{
    if(is.null(Sigma) || is.null(modelH0)) stop("Simulated power requires specification of Sigma and modelH0.")
    if(!is.null(effect) || !is.null(effect.measure))
      warning("Ignoring effect and effect.measure in simulated power.")
    
    if(!is.list(Sigma)) Sigma <- list(Sigma)
    if(any(!sapply(Sigma, isSymmetric)))
      stop("Sigma must be symmetric square matrices")
    if(any(sapply(Sigma, ncol) < 2))
      stop("Sigma must be at least of dimension 2*2")
    if(any(sapply(Sigma, function(x){sum(eigen(x)$values < 0) > 0})))
      stop("Sigma must be positive definite")
    
    if(length(Sigma) == 1) Sigma <- Sigma[[1]]
    if(!is.null(mu) && length(mu) == 1) mu <- mu[[1]]
    if(is.list(Sigma)){
      if(!is.null(mu) && length(mu) != length(Sigma))
        stop("Multiple group power analyses require specification of mu for each group.")
      if(is.null(N) && power.type != "a-priori") stop("Multiple group power analyses require specification of N for each group (representing weights in a priori power analysis).")
      if(is.null(N) && power.type == "a-priori") warning("No sample weights provided, assuming equal sample sizes in each group.")
      if(length(N) == 1){
        warning("Only single sample size provided in multiple group power analyses, assuming equal sample sizes (weights in a priori power analyses) for each group.")
      }else if(length(Sigma) != length(N)){
        stop("Multiple group power analyses require specification of N for each group.")
      } 
    }else{
      if(!is.null(mu) && length(mu) != ncol(Sigma)) stop("Mu must have the same length as ncol(Sigma).")
    }

  }
  
  # specifics depending on type of power analyses
  if(power.type == "post-hoc" || power.type == "compromise" || (power.type == "a-priori" && is.list(effect))){
    if(is.null(N)) stop('N is not defined.')
    lapply(N, checkPositive, message = 'N')
  }
  
  if(power.type == "a-priori" || power.type == "post-hoc"){
    checkBounded(alpha)
  }
  
  if(power.type == "a-priori"){
    if(is.null(beta) && is.null(power))
      stop("Need to define either beta or power in a-priori power analyis")
    if(!is.null(beta) && !is.null(power) && power != (1 - beta))
      stop("Either set beta or set power, but not both.")
    if(!is.null(beta))
      checkBounded(beta)
    if(!is.null(power))
      checkBounded(power)
  }
  
  if(power.type == "compromise"){
    checkPositive(abratio)
  }
  
  # specifics for power plots
  if(power.type == "powerplot.byN"){
    checkBounded(alpha)
    checkBounded(power.max)
    checkBounded(power.min)
    checkPositive(steps)
    checkPositive(linewidth)
  }
  
  if(power.type == "powerplot.byEffect"){
    checkPositive(N)
    checkBounded(alpha)
    checkPositive(effect.min)
    checkPositive(effect.max)
    checkPositive(steps)
    checkPositive(linewidth)
  }

}


#' checkPositive
#'
#' Checks whether `x` is defined and a positive number, stop otherwise.
#' @param x x
#' @param message identifier for `x`
checkPositive <- function(x, message = NULL){
  if(is.null(message)) message <- deparse(substitute(x))
  if(is.null(x) || is.na(x) || x <= 0){
    stop(paste(message, " must be larger than zero"))
  }
}

#' checkBounded
#'
#' Checks whether x is defined and lies within the specified bound, stop otherwise.
#' @param x x
#' @param message identifier for x
#' @param bound the boundaries, array of size two
#' @param inclusive whether x might lie on boundary
checkBounded <- function(x, message = NULL, bound = c(0, 1), inclusive = FALSE){
  if(is.null(message)) message <- deparse(substitute(x))
  inv <- is.null(x) || is.na(x) || !is.numeric(x)
  if(!inv && !inclusive && (x <= bound[1] || x >= bound[2])) inv <- TRUE
  if(!inv && inclusive && (x < bound[1] || x > bound[2])) inv <- TRUE
  if(inv) stop(paste(message, "must must lie within", bound[1], 'and', bound[2]))
}

#' checkPositiveDefinite
#'
#' Checks whether `x` is positive definite, stop otherwise.
#' 
#' @param x x
#' @param message identifier for `x`
#' @param stop whether to stop or to throw a warning
checkPositiveDefinite <- function(x, message = NULL, stop = TRUE){
  if(is.null(message)) message <- deparse(substitute(x))
  checkSymmetricSquare(x)
  if(sum(eigen(x)$values < 0) > 0){
    if(stop) stop(paste(message, " must be positive definite"))
    if(!stop) warning(paste(message, " must be positive definite"))
  }
}

#' checkSymmetricSquare
#'
#' Checks whether `x` is a symmetric square matrix, stop otherwise.
#' 
#' @param x x
#' @param message identifier for `x`
checkSymmetricSquare <- function(x, message = NULL){
  if(is.null(message)) message <- deparse(substitute(x))
  checkSquare(x)
  if(!isSymmetric(x))
    stop(paste(message, " must be a symmetric square matrix"))
}

#' checkSquare
#'
#' Checks whether `x` is a square matrix, stop otherwise.
#' 
#' @param x x
#' @param message identifier for `x`
checkSquare <- function(x, message = NULL){
  if(is.null(message)) message <- deparse(substitute(x))
  if(is.null(x))
    stop(paste(message, " may not be NULL"))
  if(!is.numeric(x))
    stop(paste(message, " must contain numeric elements only"))
  if(!is.matrix(x))
    stop(paste(message, " must be a matrix"))
  if(ncol(x) != nrow(x))
    stop(paste(message, " must be a square matrix"))
}

#' checkPowerTypes
#'
#' Checks whether type is one of `'a-priori'`, `'post-hoc'`, or `'compromise'` (or respective shortcuts), stop otherwise.
#' 
#' @param type type
#' @return Returns cleaned type
checkPowerTypes <- function(type){
  if(is.null(type) || length(type) != 1 || typeof(type) != 'character') stop('Type must be one of a-priori, post-hoc, or compromise.')
  type <- tolower(trimws(type))
  if(type == 'a priori' || type == 'apriori' || type == 'a_priori' || type == 'ap') type <- 'a-priori'
  if(type == 'post hoc' || type == 'posthoc' || type == 'post_hoc' || type == 'ph') type <- 'post-hoc'
  if(type == 'co' || type == 'comp') type <- 'compromise'
  if(!type %in% c('a-priori', 'post-hoc', 'compromise')) stop('Type must be one of a-priori, post-hoc, or compromise.')
  type
}

#' checkComparisonModel
#'
#' Checks whether comparison is one of `'restricted'` or `'saturated'` (or respective shortcuts), stop otherwise.
#' 
#' @param comparison comparison
#' @return Returns cleaned comparison
checkComparisonModel <- function(comparison){
  if(is.null(comparison) || length(comparison) != 1 || typeof(comparison) != 'character') stop('Comparison model must be one of "saturated" or "restricted".')
  comparison <- tolower(trimws(comparison))
  if(comparison == 'saturated' || comparison == 'sat') comparison <- 'saturated'
  if(comparison == 'restricted' || comparison == 'restr') comparison <- 'restricted'
  if(!comparison %in% c('saturated', 'restricted')) stop('Comparison model must be one of "saturated" or "restricted"')
  comparison
}

#' checkNullEffect
#'
#' Checks whether `nullEffect` is one of the valid effects, stop otherwise.
#' 
#' @param nullEffect nullEffect
#' @param valid vector of valid effects 
#' @param message message 
#' @return Returns cleaned nullEffect
checkNullEffect <- function(nullEffect, valid, message = NULL){
  if(is.null(message)) message <- deparse(substitute(nullEffect))
  if(is.null(nullEffect)) stop(paste(message, 'must be defined.'))
  if(length(nullEffect) > 1) stop(paste(message, 'must contain a single hypothesis'))
  nullEffect <- unlist(lapply(nullEffect, function(x) tolower(trimws(x))))
  nullEffect <- gsub(" ", "", nullEffect, fixed = TRUE)
  if(any(unlist(lapply(nullEffect, function(x) !x %in% valid)))) stop(paste(message, 'must be one of', paste(valid, collapse = ' ')))
  nullEffect
}

#' checkDataGenerationTypes
#'
#' Checks whether data generation type is one of `'normal'`, `'IG'`, `'mnonr'`, `'RK'`, or `'VM'`, stop otherwise.
#' 
#' @param type type
#' @return Returns cleaned data generation type
checkDataGenerationTypes <- function(type){
  if(is.null(type) || length(type) != 1 || typeof(type) != 'character') stop('Data generation type is invalid.')
  type <- tolower(trimws(type))
  if(type == 'norm') type <- 'normal'
  if(!type %in% c('normal', 'ig', 'mnonr', 'rk', 'vm')) stop('Data generation type must be one of "normal", "IG", "mnonr", "RK", or "VM"')
  type
}

#' checkMissingTypes
#'
#' Checks whether missing generation type is one of `'mcar'`, `'mar'`, or `'nmar'`, stop otherwise.
#' 
#' @param type type
#' @return Returns cleaned data generation type
checkMissingTypes <- function(type){
  if(is.null(type) || length(type) != 1 || typeof(type) != 'character') stop('Missing mechanism must be one of "mcar", "mar", or "mnar".')
  type <- tolower(trimws(type))
  if(!type %in% c('mcar', 'mar', 'nmar')) stop('Missing mechanism must be one of "mcar", "mar", or "mnar"')
  type
}

#' checkEllipsis
#'
#' Checks whether `...` contains arguments related to loadings and to power, stop otherwise.
#' 
#' @param ... the parameters to search.
checkEllipsis <- function(...){
  args <- names(list(...))
  if(!any(c('Lambda', 'loadings', 'nIndicator', 'loadM') %in% args)) stop('Missing arguments specifying the factor model. Remember to set either Lambda, loadings, or nIndicator and loadM.')
  if(!('alpha' %in% args || 'abratio' %in% args)) stop('Missing arguments related to power analysis. Remember to set alpha, power, N, and/or abratio.')
}
