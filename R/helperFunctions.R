#' semPower.genSigma
#'
#' Generate a covariance matrix (and mean vector) and associated lavaan model strings based on defined model features.
#' 
#' @param Lambda Factor loading matrix. 
#' @param Phi Factor correlation (or covariance) matrix or single number giving correlation between all factors or NULL for uncorrelated factors. 
#' @param Beta Regression slopes between latent variables (all-y notation). 
#' @param Psi Variance-covariance matrix of latent residuals when Beta is specified. If NULL, a diagonal matrix is assumed. 
#' @param Theta Variance-covariance matrix of manifest residuals. If NULL and Lambda is not a square matrix, Theta is diagonal so that the manifest variances are 1. If NULL and Lambda is square, Theta is 0.      
#' @param tau Intercepts. If NULL and Alpha is set, these are assumed to be zero. 
#' @param Alpha Factor means. If NULL and tau is set, these are assumed to be zero. 
#' @return A list containing the implied variance-covariance matrix (Sigma), the loading matrix (Lambda), the factor-covariance matrix (Phi) or the slopes (Beta) and the residual variances (Psi), the variance-covariance matrix between the manifest residuals (Theta), the implied indicator means (mu), intercepts (tau), and latent means (Alpha), as well as the associated lavaan model string defining the population (modelPop) and two lavaan models string defining a corresponding true (modelTrue) or pure cfa analysis model (modelTrueCFA) omitting any regression relationships.
#' @details
#' This function generates the variance-covariance matrix of the \eqn{p} observed variables \eqn{\Sigma} and their means \eqn{\mu} via a confirmatory factor (CFA) model or a more general structural equation model. 
#' 
#' In the CFA model, 
#' \deqn{\Sigma = \Lambda \Phi \Lambda' + \Theta}
#' where \eqn{\Lambda} is the \eqn{p \cdot m} loading matrix, \eqn{\Phi} is the variance-covariance matrix of the \eqn{m} factors, and \eqn{\Theta} is the residual variance-covariance matrix of the observed variables. The means are
#' \deqn{\mu = \tau + \Lambda \alpha}
#' with the \eqn{p} indicator intercepts \eqn{\tau} and the \eqn{m} factor means \eqn{\alpha}.
#' 
#' In the structural equation model, 
#' \deqn{\Sigma = \Lambda (I - \Beta)^{-1} \Psi [(I - \Beta)^{-1}]'  \Lambda' + \Theta } 
#' where \eqn{\Beta} is the \eqn{m \cdot m} matrix containing the regression slopes and \eqn{\Psi} is the (residual) variance-covariance matrix of the \eqn{m} factors. The means are
#' \deqn{\mu = \tau + \Lambda (I - \Beta)^{-1} \alpha}
#' 
#' In either model, the meanstructure can be omitted, leading to factors with zero means and zero intercepts. 
#' 
#' When \eqn{\Lambda = I}, the models above do not contain any factors and reduce to ordinary regression or path models.  
#' 
#' If `Phi` is defined, a CFA model is used, if `Beta` is defined, a structural equation model. 
#' When both `Phi` and `Beta` are NULL, a CFA model is used with \eqn{\Phi = I}, i.e., uncorrelated factors.
#' When `Phi` is a single number, all factor correlations are equal to this number.
#' 
#' When `Beta` is defined and `Psi` is NULL, \eqn{\Psi = I}.
#' 
#' When `Theta` is NULL, \eqn{\Theta} is a diagonal matrix with all elements such that the variances of the observed variables are 1. When there is only a single observed indicator for a factor, the corresponding element in \eqn{\Theta} is set to zero.
#'
#' Instead of providing the loading matrix \eqn{\Lambda} via `Lambda`, there are several shortcuts:
#' \itemize{
#' \item `loadings`: defines the primary loadings for each factor in a list structure, e. g. `loadings = list(c(.5, .4, .6), c(.8, .6, .6, .4))` defines a two factor model with three indicators loading on the first factor by .5, , 4., and .6, and four indicators loading in the second factor by .8, .6, .6, and .4.  
#' \item `nIndicator`: used in conjunction with `loadM` or `loadMinmax`, defines the number of indicators by factor, e. g., `nIndicator = c(3, 4)` defines a two factor model with three and four indicators for the first and second factor, respectively. `nIndicator` can also be a single number to define the same number of indicators for each factor. 
#' \item `loadM`: defines the mean loading either for all indicators (if a single number is provided) or separately for each factor (if a vector is provided), e. g. `loadM = c(.5, .6)` defines the mean loadings of the first factor to equal .5 and those of the second factor do equal .6
#' \item `loadSD`: used in conjunction with `loadM`, defines the standard deviations of the loadings. If omitted or NULL, the standard deviations are zero. Otherwise, the loadings are sampled from a normal distribution with N(loadM, loadSD) for each factor. 
#' \item `loadMinMax`: defines the minimum and maximum loading either for all factors or separately for each factor (as a list). The loadings are then sampled from a uniform distribution. For example, `loadMinMax = list(c(.4, .6), c(.4, .8))` defines the loadings for the first factor lying between .4 and .6, and those for the second factor between .4 and .8. 
#' }      
#'  
#' @examples
#' \dontrun{
#' # Provide the factor correlation for a two-factor model, the number of indicators by factor, 
#' # and a single loading which is equal for all indicators
#' genSigma <- semPower.genSigma(Phi = .2, nIndicator = c(5, 6), loadM = .5)
#' genSigma$Sigma # implied covariance matrix
#' 
#' # Provide Beta for a three factor model 
#' # with 3, 4, and 5 indicators loading by .6, 5, and .4, respectively. 
#' Beta <- matrix(c(
#'                 c(0.0, 0.0, 0.0),
#'                 c(0.3, 0.0, 0.0),  # f2 = .3*f1
#'                 c(0.2, 0.4, 0.0)   # f3 = .2*f1 + .4*f2
#'                ), byrow = TRUE, ncol = 3)
#' genSigma <- semPower.genSigma(Beta = Beta, nIndicator = c(3, 4, 5), loadM = c=(.6, .5, .4)
#' 
#' # Provide factor correlation matrix and loading matrix
#' Phi <- matrix(c(
#'                 c(1.0, 0.2, 0.5),
#'                 c(0.2, 1.0, 0.3),
#'                 c(0.5, 0.3, 1.0)
#'                ), byrow = TRUE, ncol = 3)
#' Lambda <- matrix(c(
#'                 c(0.5, 0.0, 0.0),
#'                 c(0.4, 0.0, 0.0),
#'                 c(0.3, 0.0, 0.0),
#'                 c(0.0, 0.7, 0.0),
#'                 c(0.0, 0.8, 0.0),
#'                 c(0.0, 0.5, 0.0),
#'                 c(0.0, 0.0, 0.5),
#'                 c(0.0, 0.0, 0.4),
#'                 c(0.0, 0.0, 0.6),
#'                ), byrow = TRUE, ncol = 3)
#'                
#' genSigma <- semPower.genSigma(Phi = Phi, Lambda = Lambda)
#' 
#' # same as above, but providing a reduced loading matrix, i.e.
#' # only define the primary loadings for each factor
#' loadings <- list(
#'                c(0.4, 0.5, 0.8),
#'                c(0.7, 0.6, 0.5, 0.4, 0.5),
#'                c(0.5, 0.5, 0.6, 0.8)
#'                )
#'                
#' genSigma <- semPower.genSigma(Phi = Phi, loadings = loadings)
#'   
#' # same as above, but providing
#' # the number of indicators by factor 
#' # and min-max loading for all factors (sampling from a uniform distribution)
#' genSigma <- semPower.genSigma(Phi = Phi, nIndicator = c(3, 5, 4), 
#'                               loadMinMax = c(.3, .8))
#'
#' # same as above, but providing
#' # mean and SD loading for all factors (sampling from a normal distribution)
#' genSigma <- semPower.genSigma(Phi = Phi, nIndicator = c(3, 5, 4), 
#'                               loadM = .5, loadSD = .1)
#'
#'
#' # same as above, but providing mean and sd of loadings for each factor
#' genSigma <- semPower.genSigma(Phi = Phi, nIndicator = c(3, 5, 4), 
#'                               loadM = c(.5, .6, .7), loadSD = c(0, .05, .01))
#'
#' # same as above, but using min-max loadings for each factor
#' loadMinMax <- list(
#'                    c(.4, .6),
#'                    c(.5, .8),
#'                    c(.3, .7)
#'                    )
#' genSigma <- semPower.genSigma(Phi = Phi, nIndicator = c(3, 5, 4), 
#'                               loadMinMax = loadMinMax)
#'                               
#' }
#' @export
semPower.genSigma <- function(Lambda = NULL,
                              Phi = NULL, 
                              Beta = NULL,  # capital Beta, to distinguish from beta error
                              Psi = NULL,
                              Theta = NULL,
                              tau = NULL,
                              Alpha = NULL,  # capital Alpha, to distinguish from alpha error
                              ...){
  
  args <- list(...)

  # multigroup case
  argsMG <- c(as.list(environment()), list(...))
  argsMG <- argsMG[!unlist(lapply(argsMG, is.null)) & names(argsMG) %in% c('Lambda', 'Phi', 'Beta', 'Psi', 'Theta', 'tau', 'Alpha', 'nIndicator', 'loadM', 'loadSD', 'loadings', 'loadMinMax')]
  nGroups <- unlist(lapply(seq_along(argsMG), function(x){
    len <- 1
    if(names(argsMG)[x] %in% c('loadings', 'loadMinMax')){
      if(is.list(argsMG[[x]][[1]])) len <- length(argsMG[[x]][[1]])
    }else{
      if(is.list(argsMG[[x]])) len <- length(argsMG[[x]])
    }
    len
  }))
  if(any(nGroups > 1)){
    if(sum(nGroups != 1) > 1 && var(nGroups[nGroups != 1]) != 0) stop('All list arguments in multiple group analysis must have the same length.')
    if(!is.null(argsMG[['loadM']]) && !is.list(argsMG[['loadM']])) argsMG[['loadM']] <- as.list(rep(list(argsMG[['loadM']]), max(nGroups)))
    if(!is.null(argsMG[['loadSD']]) && !is.list(argsMG[['loadSD']])) argsMG[['loadSD']] <- as.list(rep(list(argsMG[['loadSD']]), max(nGroups)))
    if(!is.null(argsMG[['loadMinMax']]) && !is.list(argsMG[['loadMinMax']])) argsMG[['loadMinMax']] <- as.list(rep(list(argsMG[['loadMinMax']]), max(nGroups)))
    params <- lapply(1:max(nGroups), function(x) lapply(argsMG, '[[', x))
    return(
      lapply(params, function(x) do.call(semPower.genSigma, x)) 
    )
  }

  # the following applies to single group cases 
  if(is.null(Lambda)){
    Lambda <- genLambda(args[['loadings']], args[['nIndicator']], 
                        args[['loadM']], args[['loadSD']], args[['loadMinMax']]) 
  }

  nfac <- ncol(Lambda)
  nIndicator <- apply(Lambda, 2, function(x) sum(x != 0))
  
  ### validate input
  if(ncol(Lambda) > 99) stop("Models with >= 100 factors are not supported.")
  if(!is.null(Beta) && !is.null(Phi)) stop('Either provide Phi or Beta, but not both. Did you mean to set Beta and Psi?')
  if(is.null(Beta) && is.null(Phi)) Phi <- diag(nfac)
  if(is.null(Beta)){
    if(length(Phi) == 1){
      Phi <- matrix(Phi, ncol = nfac, nrow = nfac)
      diag(Phi) <- 1
    } 
    if(ncol(Phi) != nfac) stop('Phi must have the same number of rows/columns as the number of factors.') 
    checkPositiveDefinite(Phi)
    if(!any(diag(Phi) != 1)){
      if(any(Phi > 1) || any(Phi > 1)) stop('Phi implies correlations outside [-1 - 1].')
    }
  }else{
    if(ncol(Beta) != nfac) stop('Beta must have the same number of rows/columns as the number of factors.')
    if(!is.null(Psi) && ncol(Psi) != nfac) stop('Psi must have the same number of rows/columns as the number of factors.')
    if(is.null(Psi)) Psi <- diag(ncol(Beta))
    if(any(diag(Psi) < 0)) stop('Model implies negative residual variances for latent variables (Psi)')
  }
  if(!is.null(tau)){
    if(length(tau) != sum(nIndicator)) stop('Intercepts (tau) must be of same length as the number of indicators')
  }
  if(!is.null(Alpha)){
    if(length(Alpha) != nfac) stop('Latent means (Alpha) must be of same length as the number of factors')
  }
  if(!is.null(tau) && is.null(Alpha)) Alpha <- rep(0, nfac)
  if(!is.null(Alpha) && is.null(tau)) tau <- rep(0, sum(nIndicator))

  ### compute Sigma
  if(is.null(Beta)){
    SigmaND <- Lambda %*% Phi %*% t(Lambda) 
  }else{
    invIB <- solve(diag(ncol(Beta)) - Beta)    
    SigmaND <- Lambda %*% invIB %*% Psi %*% t(invIB) %*% t(Lambda) 
  }
  # compute theta
  if(is.null(Theta)){
    # catch observed only models
    if(!any(nIndicator > 1)){
      Theta <- matrix(0, ncol = ncol(SigmaND), nrow = nrow(SigmaND))
    # catch latent+observed mixed models
    }else if(any(nIndicator <= 1)){
      fLat <- which(nIndicator > 1)
      indLat <- unlist(lapply(fLat, function(x) which(Lambda[, x] != 0)))
      indObs <- (1:nrow(Lambda))[-indLat]
      Theta <- diag(1 - diag(SigmaND))
      diag(Theta)[indObs] <- 0
    }else{
      Theta <- diag(1 - diag(SigmaND))
    }
  }
  # see whether there are negative variances
  if(any(diag(Theta) < 0)) stop('Model implies negative residual variances for observed variables (Theta)')
  
  Sigma <- SigmaND + Theta
  colnames(Sigma) <- rownames(Sigma) <- paste0('x', 1:ncol(Sigma)) # set row+colnames to make isSymmetric() work

  ### compute Mu
  mu <- NULL
  if(!is.null(tau)){
    if(is.null(Beta)){
      mu <- tau + Lambda %*% Alpha
    }else{
      mu <- tau + Lambda %*% invIB %*% Alpha
    }
    names(mu) <- paste0('x', 1:length(mu))
  }
  
  # also provide several lav model strings
  modelStrings <- genModelString(Lambda = Lambda,
                                 Phi = Phi, Beta = Beta, Psi = Psi, Theta = Theta,
                                 tau = tau, Alpha = Alpha,
                                 useReferenceIndicator = ifelse(is.null(args[['useReferenceIndicator']]), !is.null(Beta), args[['useReferenceIndicator']]), 
                                 metricInvariance = args[['metricInvariance']])

  append(
    list(Sigma = Sigma, 
       mu = mu,
       Lambda = Lambda, 
       Phi = Phi,
       Beta = Beta, 
       Psi = Psi,
       Theta = Theta,
       tau = tau,
       Alpha = Alpha),
    modelStrings
    )
}

#' genLambda
#'
#' Generate Lambda from various shortcuts
#' 
#' @param loadings A list providing the factor loadings by factor. Must not contain secondary loadings.   
#' @param nIndicator Vector indicating the number of indicators for each factor, e.g. c(4, 6) to define two factors with 4 and 6 indicators, respectively 
#' @param loadM Either a vector giving mean loadings for each factor or a single number to use for every loading
#' @param loadSD Either a vector giving the standard deviation of loadings for each factor or a single number, for use in conjunction with loadM. When NULL, SDs are set to zero.
#' @param loadMinMax A list giving the minimum and maximum loading for each factor or vector to apply to all factors 
#' @return A list containing the implied variance-covariance matrix (Sigma), the loading matrix (Lambda), the factor-covariance matrix (Phi) or the slopes (Beta) and the residual variances (Psi), the variance-covariance matrix between the manifest residuals (Theta), the implied indicator means (mu), intercepts (tau), and latent means (Alpha), as well as the associated lavaan model string defining the population (modelPop) and two lavaan models string defining a corresponding true (modelTrue) or pure cfa analysis model (modelTrueCFA) omitting any regression relationships.
#' @examples
#' \dontrun{
#' }
#' @importFrom stats rnorm runif 
genLambda <- function(loadings = NULL, 
                      nIndicator = NULL,
                      loadM = NULL,
                      loadSD = NULL,
                      loadMinMax = NULL){

  # validate input  
  if(!is.null(nIndicator) && !is.null(loadings)) stop('Either provide loadings or number of indicators, but not both.')
  if(is.null(nIndicator) && !(is.null(loadM) || is.null(loadMinMax))) stop('Provide loadings and number of indicators by factor.')
  if(is.null(nIndicator) && !is.list(loadings)) stop('loadings must be a list')
  nfac <- ifelse(is.null(loadings), length(nIndicator), length(loadings))
  if(is.null(loadings)){
    if(any(!sapply(nIndicator, function(x) x %% 1 == 0))) stop('Number of indicators must be a integer')
    invisible(sapply(nIndicator, function(x) checkBounded(x, 'Number of indicators ', bound = c(1, 10000), inclusive = TRUE)))
    if(is.null(loadM) && is.null(loadMinMax)) stop('Either mean loading or min-max loading need to be defined')
    if(is.null(loadMinMax) && length(loadM) == 1) loadM <- rep(loadM, nfac)
    if(is.null(loadMinMax) && length(loadM) != nfac) stop('Nindicator and mean loading must of same size')
    
    if(!is.null(loadMinMax) && (!is.null(loadM) || !is.null(loadSD))) stop('Either specify mean and SD of loadings or specify min-max loading, both not both.')
    if(length(loadMinMax) == 2) loadMinMax <- lapply(1:nfac, function(x) loadMinMax)
    if(is.null(loadMinMax) && is.null(loadSD)) loadSD <- rep(0, nfac)
    if(is.null(loadMinMax) && length(loadSD) == 1) loadSD <- rep(loadSD,nfac)
    if(is.null(loadMinMax) && !is.null(loadSD)) invisible(sapply(loadSD, function(x) checkBounded(x, 'Standard deviations', bound = c(0, .5), inclusive = TRUE)))
  }else{
    nIndicator <- unlist(lapply(loadings, length))  # crossloadings are disallowed
  }  

  Lambda <- matrix(0, ncol = nfac, nrow = sum(nIndicator))
  sidx <- 1
  for(f in 1:nfac){
    eidx <- sidx + (nIndicator[f] - 1)
    if(!is.null(loadM)){
      cload <- round(rnorm(nIndicator[f], loadM[f], loadSD[f]), 2)
      if(any(cload < -1) || any(cload > 1)) warning('Sampled loadings outside [-1, 1] were set to -1 or 1.')
      cload[cload < -1] <- -1
      cload[cload > 1] <- 1
    }else if(!is.null(loadMinMax)){
      cload <- round(runif(nIndicator[f], loadMinMax[[f]][1], loadMinMax[[f]][2]), 2)
    }else if(!is.null(loadings)){
      cload <- loadings[[f]]  # loadings only contains primary loadings
    }else{
      stop('loading not found')  
    }
    Lambda[sidx:eidx, f] <- cload
    sidx <- eidx + 1
  } 
  
  Lambda
}

#' genModelString
#'
#' creates lavaan model strings from model matrices.
#' @param Lambda Factor loading matrix. 
#' @param Phi Factor correlation (or covariance) matrix or single number giving correlation between all factors or NULL for uncorrelated factors. 
#' @param Beta Regression slopes between latent variables (all-y notation). 
#' @param Psi Variance-covariance matrix of latent residuals when Beta is specified. If NULL, a diagonal matrix is assumed. 
#' @param Theta Variance-covariance matrix of manifest residuals. If NULL and Lambda is not a square matrix, Theta is diagonal so that the manifest variances are 1. If NULL and Lambda is square, Theta is 0.      
#' @param tau Intercepts. If NULL and Alpha is set, these are assumed to be zero. 
#' @param Alpha Factor means. If NULL and tau is set, these are assumed to be zero. 
#' @param useReferenceIndicator Whether to identify factors in accompanying model strings by a reference indicator (TRUE) or by setting their variance to 1 (FALSE). When Beta is defined, a reference indicator is used by default, otherwise the variance approach. 
#' @param metricInvariance A list containing the factor indices for which the accompanying model strings should apply metric invariance labels, e.g. list(c(1,2), c(3,4)) to assume invariance for f1 and f2 as well as f3 and f4.  
#' @return A list containing the lavaan model string defining the population (modelPop) and two lavaan models string defining a corresponding true (modelTrue) or pure cfa analysis model (modelTrueCFA) omitting any regression relationships.
#' @importFrom stats var
genModelString <- function(Lambda = NULL,
                           Phi = NULL,
                           Beta = NULL,  # capital Beta, to distinguish from beta error
                           Psi = NULL,
                           Theta = NULL,
                           tau = NULL,
                           Alpha = NULL,  # capital Alpha, to distinguish from alpha error
                           useReferenceIndicator = !is.null(Beta),
                           metricInvariance = NULL){

  nfac <- ncol(Lambda)
  nIndicator <- apply(Lambda, 2, function(x) sum(x != 0))
  
  # validate input
  if(!is.null(metricInvariance)){
    if(!is.list(metricInvariance)) stop('metricInvariance must be a list')
    if(any(unlist(lapply(metricInvariance, function(x) length(x))) < 2)) stop('each list entry in metricInvariance must involve at least two factors')
    if(max(unlist(metricInvariance)) > nfac || min(unlist(metricInvariance)) <= 0) stop('factor index < 1 or > nfactors in metricInvariance')
    if(any(unlist(lapply(metricInvariance, function(x) var(nIndicator[x]))) != 0)) stop('factors in metricInvariance must have the same number of indicators')
    metricInvarianceLabels <- lapply(1:length(metricInvariance), function(x) paste0('l', formatC(x, width = 2, flag = 0), formatC(1:nIndicator[metricInvariance[[x]][1]], width = 2, flag = 0), '*')) 
  }
  
  # remove names if set. we currently don't support labels anyway
  # and this causes issues later in building the model strings
  rownames(Lambda) <- colnames(Lambda) <- NULL  
  rownames(Theta) <- colnames(Theta) <- NULL  
  rownames(Phi) <- colnames(Phi) <- NULL  
  rownames(Beta) <- colnames(Beta) <- NULL  
  rownames(Psi) <- colnames(Psi) <- NULL  
  names(tau) <- NULL
  names(Alpha) <- NULL
  
  # create lav population model string
  tok <- list()
  for(f in 1:ncol(Lambda)){
    iIdx <- which(Lambda[, f] != 0)
    if(!identical(iIdx, integer(0))){
      cload <- Lambda[iIdx, f]
      tok <- append(tok, paste0('f', f, ' =~ ', paste0(cload, '*', paste0('x', iIdx), collapse = ' + ')))
    }
    if(!is.null(Phi)){
      tok <- append(tok, paste0('f', f, ' ~~ ', Phi[f, f],'*', 'f', f))
    }else{
      tok <- append(tok, paste0('f', f, ' ~~ ', Psi[f, f],'*', 'f', f))
    }
  }
  # manifest residuals
  for(i in 1:nrow(Lambda)){
    tok <- append(tok, paste0(paste0('x', i), ' ~~ ', Theta[i, i], '*', paste0('x', i), collapse = '\n'))
  }
  # define factor cor / residual cor / regressions
  if(nfac > 1){
    if(!is.null(Phi)){
      for(f in 1:(nfac - 1)){
        for(ff in (f + 1):nfac){
          tok <- append(tok, paste0('f', f, ' ~~ ', Phi[f, ff], '*f', ff))
        }
      }
    }else{
      # regressions
      for(f in 1:ncol(Beta)){
        idx <- which(Beta[f, ] != 0)
        if(!identical(idx, integer(0)))
          tok <- append(tok, paste0('f', f, ' ~ ', paste0(Beta[f, idx], '*f', idx, collapse = '+'))) 
      }
      # (residual) correlations
      for(f in 1:(nfac - 1)){
        for(ff in (f + 1):nfac){
          if(Psi[f, ff] != 0)
            tok <- append(tok, paste0('f', f, ' ~~ ', Psi[f, ff], '*f', ff))
        }
      }
    }
  }
  # add means
  if(!is.null(tau)){
    tok <- append(tok, paste0(paste0('x', 1:nrow(Lambda)), ' ~ ', tau, '*1', collapse = '\n'))
  }
  if(!is.null(Alpha)){
    tok <- append(tok, paste0(paste0('f', 1:ncol(Lambda)), ' ~ ', Alpha, '*1', collapse = '\n'))
  }
  
  modelPop <- paste(unlist(tok), collapse = '\n')
  
  # build  analysis models strings, one pure cfa based, and one including regression relationships 
  tok <- list()
  for(f in 1:ncol(Lambda)){
    iIdx <- which(Lambda[, f] != 0)
    if(!identical(iIdx, integer(0))){
      # add invariance constraints
      if(any(unlist(lapply(metricInvariance, function(x) f %in% x)))){
        labelIdx <- which(unlist(lapply(metricInvariance, function(x) f %in% x)))
        clabel <- metricInvarianceLabels[[labelIdx]]
        if(useReferenceIndicator){
          # scale by first loading instead of 1
          tok <- append(tok, paste0('f', f, ' =~ ', Lambda[iIdx[1], f], '*x', iIdx[1], ' + ', paste0(clabel, 'x', iIdx, collapse = ' + ')))
        }else{
          tok <- append(tok, paste0('f', f, ' =~ NA*x', iIdx[1],' + ', paste0(clabel, 'x', iIdx, collapse = ' + ')))
        }
      }else{
        if(useReferenceIndicator){
          # scale by first loading instead of 1
          tok <- append(tok, paste0('f', f, ' =~ ', Lambda[iIdx[1], f], '*', paste0('x', iIdx, collapse = ' + ')))
        }else{
          tok <- append(tok, paste0('f', f, ' =~ NA*', paste0('x', iIdx, collapse = ' + ')))
        }
      }
    }
  }
  modelTrueCFA <- paste(c(unlist(tok)), collapse = '\n')
  if(!is.null(Beta)){
    # regressions
    for(f in 1:ncol(Beta)){
      idx <- which(Beta[f, ] != 0)
      if(!identical(idx, integer(0)))
        tok <- append(tok, paste0('f', f, ' ~ ', paste0('f', idx, collapse = '+'))) 
    }
    # (residual) correlations
    for(f in 1:(nfac - 1)){
      for(ff in (f + 1):nfac){
        if(Psi[f, ff] != 0)
          tok <- append(tok, paste0('f', f, ' ~~ f', ff))
      }
    }
  }
  modelTrue <- paste(c(unlist(tok)), collapse = '\n')
  if(!useReferenceIndicator){
    fvars <- sapply(1:nfac, function(f) paste0('f', f, ' ~~ 1*f', f))
    # factor variances are always 1, regardless of phi
    modelTrue <- paste(c(modelTrue, fvars), collapse = '\n')
    modelTrueCFA <- paste(c(modelTrueCFA, fvars), collapse = '\n')
  }
  
  list(modelPop = modelPop, modelTrue = modelTrue, modelTrueCFA = modelTrueCFA)
}


#' getPhi.B
#'
#' Computes implied correlations (completely standardized) from Beta matrix, disallowing recursive paths.
#' 
#' @param B matrix of regression coefficients (all-y notation). Must only contain non-zero lower-triangular elements, so the first row only includes zeros. 
#' @param lPsi (lesser) matrix of residual correlations. This is not the Psi matrix, but a lesser version ignoring all variances and containing correlations off the diagonal. Can be omitted for no correlations beyond those implied by B. 
#' @return the implied correlation matrix
#' @examples
#' \dontrun{
#' # mediation model
#' B <- matrix(c(
#'   c(.00, .00, .00),
#'   c(.10, .00, .00),
#'   c(.20, .30, .00)
#' ), byrow = TRUE, ncol = 3)
#' Phi <- getPhi.B(B)
#' 
#' # CLPM with residual correlations 
#' B <- matrix(c(
#'   c(.00, .00, .00, .00),
#'   c(.30, .00, .00, .00),
#'   c(.70, .10, .00, .00),
#'   c(.20, .70, .00, .00)
#' ), byrow = TRUE, ncol = 4)
#' lPsi <- matrix(c(
#'   c(.00, .00, .00, .00),
#'   c(.00, .00, .00, .00),
#'   c(.00, .00, .00, .30),
#'   c(.00, .00, .30, .00)
#' ), byrow = TRUE, ncol = 4)
#' Phi <- getPhi.B(B, lPsi)
#' }
getPhi.B <- function(B, lPsi = NULL){
  
  checkSquare(B)
  if(any(B[upper.tri(B, diag = TRUE)] != 0)) stop('B may not contain any non-zero values on or upper the diagonal.')
  if(any(rowSums(B^2) > 1)) stop('B implies negative residual variances.')
  invisible(lapply(B, function(x) lapply(x, function(x) checkBounded(x, 'All elements in B', bound = c(-1, 1)))))
  if(!is.null(lPsi)){
    checkSymmetricSquare(lPsi)
    if(ncol(lPsi) != ncol(B)) stop('lPsi must be of same dimension as B')
    invisible(lapply(lPsi, function(x) lapply(x, function(x) checkBounded(x, 'All elements in lPsi', bound = c(-1, 1)))))
    diag(lPsi) <- 0
  }

  ## to obtain phi in std metrc, exploit the structure of B to build Phi recursively
  ## there must be a simpler way to do this...
  exog <- apply(B, 1, function(x) !any(x != 0))
  Be <- B[!exog, ]
  if(!is.matrix(Be)) Be <- t(matrix(Be))
  
  Phi <- diag(ncol(B)) 
  for(i in 1:nrow(Be)){
    idx <- i + sum(exog) - 1
    cb <- matrix(Be[i, 1:idx])
    cr <- Phi[1:idx, 1:idx]
    predR <- (cr %*% cb)
    # add residual covariances
    if(!is.null(lPsi) && any(lPsi[idx, 1:(idx + 1)] != 0)){
      cR <- rbind(cr, t(predR))
      cR <- cbind(cR, t(t(c((predR),1))))
      cB <- B[1:(idx + 1), 1:(idx + 1)]
      rootResidualVar <-  sqrt(diag(1 - diag(cB %*% cR %*% t(cB)))) 
      cPsi <- lPsi[1:(idx + 1), 1:(idx + 1)]
      corR <- cR + rootResidualVar %*% cPsi %*% t(rootResidualVar) 
      predR <- corR[(idx + 1), 1:idx]
    }
    Phi[(idx + 1), 1:idx] <- t(predR)
    Phi[1:idx, (idx + 1)] <- predR
  }
  
  checkPositiveDefinite(Phi)
  
  Phi
}


#' semPower.getDf
#'
#' Determines the degrees of freedom of a given model provided as lavaan model string. 
#' This requires the lavaan package.
#' 
#' @param lavModel the lavaan model string 
#' @param nGroups for multigroup models: the number of groups 
#' @param group.equal for multigroup models: vector of strings defining the type(s) of cross-group equality constraints following the lavaan conventions (loadings, intercepts, means, residuals, residual.covariances, lv.variances, lv.covariances, regressions)
#' @return the df of the model
#' @examples
#' \dontrun{
#' lavModel <- '
#' f1 =~ x1 + x2 + x3 + x4
#' f2 =~ x5 + x6 + x7 + x8
#' f3 =~ y1 + y2 + y3
#' f3 ~ f1 + f2
#' '
#' semPower.getDf(lavModel)
#' 
#' # multigroup version
#' semPower.getDf(lavModel, nGroups = 3)  
#' semPower.getDf(lavModel, nGroups = 3, group.equal = c('loadings'))
#' semPower.getDf(lavModel, nGroups = 3, group.equal = c('loadings', 'intercepts'))
#' }
#' @importFrom utils installed.packages
#' @export
semPower.getDf <- function(lavModel, nGroups = NULL, group.equal = NULL){
  # check whether lavaan is available
  if(!'lavaan' %in% rownames(installed.packages())) stop('This function depends on the lavaan package, so install lavaan first.')
  
  # Fit model to dummy covariance matrix instead of counting parameters. 
  # This should also account for parameter restrictions and other intricacies
  # Model fitting will give warnings we just can ignore
  tryCatch({
    params <- lavaan::sem(lavModel)
    dummyS <- diag(params@Model@nvar)
    rownames(dummyS) <- params@Model@dimNames[[1]][[1]]
    if(is.null(nGroups) || nGroups == 1){
      dummyFit <- suppressWarnings(lavaan::sem(lavModel, sample.cov = dummyS, sample.nobs = 1000, warn = FALSE))
    }else{
      if(is.null(group.equal)){
        dummyFit <- suppressWarnings(lavaan::sem(lavModel, sample.cov = lapply(1:nGroups, function(x) dummyS), sample.nobs = rep(1000, nGroups), warn = FALSE))
      }else{
        dummyFit <- suppressWarnings(lavaan::sem(lavModel, sample.cov = lapply(1:nGroups, function(x) dummyS), sample.nobs = rep(1000, nGroups), group.equal = group.equal, warn = FALSE))
      }
    }
    dummyFit@test[['standard']][['df']]
  }, 
  warning = function(w){
    warning(w)
  }, 
  error = function(e){
    stop(e)
  }
  )
}

#' getLavOptions
#'
#' returns lavaan options including defaults as set in sem() as a list to be passed to lavaan() 
#' 
#' @param lavOptions additional options to be added to (or overwriting) the defaults  
#' @param isCovarianceMatrix if TRUE, also adds sample.nobs = 1000 and sample.cov.rescale = FALSE to lavoptions 
#' @param nGroupes the number of groups, 1 by default
#' @return a list of lav defaults
#' @examples
#' \dontrun{
#' }
getLavOptions <- function(lavOptions = NULL, isCovarianceMatrix = TRUE, nGroups = 1){
  # defaults as defined in lavaan::sem()
  lavOptionsDefaults <- list(int.ov.free = TRUE, int.lv.free = FALSE, auto.fix.first = TRUE,
                             auto.fix.single = TRUE, auto.var = TRUE, auto.cov.lv.x = TRUE,
                             auto.efa = TRUE, auto.th = TRUE, auto.delta = TRUE, auto.cov.y = TRUE) 
  
  # we also want N - 1 for ml based estm
  if(is.null(lavOptions[['estimator']]) || 
     (!is.null(lavOptions[['estimator']]) && toupper(lavOptions[['estimator']]) %in% c("ML", "MLM", "MLMV", "MLMVS", "MLF", "MLR")))
    lavOptionsDefaults <- append(list(likelihood = 'Wishart'), lavOptionsDefaults)
  
  if(isCovarianceMatrix){
    ns <- 1000
    if(nGroups > 1) ns <- as.list(rep(1000, nGroups))
    lavOptionsDefaults <- append(list(sample.nobs = ns, sample.cov.rescale = FALSE), 
                                 lavOptionsDefaults)
  }

  # append lavoptions to defaults, overwriting any duplicate key
  lavOptions <- append(lavOptions, lavOptionsDefaults)[!duplicated(c(names(lavOptions), names(lavOptionsDefaults)))]
  
  lavOptions
}
