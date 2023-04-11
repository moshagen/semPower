#' simulate
#'
#' Estimates empirical power using a simulation approach.
#' 
#' @param modelH0 `lavaan` model string defining the (incorrect) analysis model.
#' @param modelH1 `lavaan` model string defining the comparison model. If omitted, the saturated model is the comparison model.
#' @param Sigma population covariance matrix.
#' @param mu population means.
#' @param N sample size
#' @param alpha alpha error probability
#' @param simOptions a list of additional options specifying simulation details, see details. 
#' @param lavOptions a list of additional options passed to `lavaan`, e. g., `list(estimator = 'mlm')` to request robust ML estimation
#' @param lavOptionsH1 lavoptions when fitting `modelH1`. If `NULL`, the same as `lavOptions`.
#' @param returnFmin whether to return the mean unbiased Fmin over replications (i. e., `fmin_0 = fmin_hat - df/N`) 
#' @return Returns empirical power: `sum(p < alpha) / nReplications` or a list (if `returnFmin = TRUE`) with the following components:
#' \item{`ePower`}{the empirical power.}
#' \item{`meanFmin`}{the estimated mean unbiased Fmin over replications (i. e., `fmin_0 = fmin_hat - df/N`).}
#' \item{`meanFminGroups`}{the estimated mean unbiased Fmin by groups given as a vector, assuming the df spread equally over groups. Therefore, `meanFmin != sum(meanFminGroups)`}
#' \item{`df`}{the model df.}
#' \item{`nrep`}{the number of successful replications.}
#' \item{`convergenceRate`}{the convergence rate of the H0 model.}
#' \item{`bChiSq`}{median chi-square bias of the H1 model}
#' \item{`bLambda`}{average median bias in lambda in the H1 model}
#' \item{`bPhi`}{average median bias in phi in the H1 model}
#' \item{`bPsi`}{average median bias in psi in the H1 model}
#' \item{`bBeta`}{average median bias in beta in the H1 model}
#' @details 
#' 
#' The details of the simulation are specified in `simOptions`, which is a list that may have the following components:
#' * `nReplications`: The targeted number of valid simulation runs, defaults to 250.
#' * `minConvergenceRate`:  The minimum convergence rate required, defaults to .5. The maximum actual simulation runs are increased by a factor of 1/minConvergenceRate.
#' * `type`: specifies whether the data should be generated from a population assuming multivariate normality (`'normal'`; the default), or based on an approach generating non-normal data (`'IG'`, `'mnonr'`, or `'RK'`). 
#' The approaches generating non-normal data require additional arguments detailed below.
#' * `missingVars`: vector specifying the variables containing missing data (defaults to `NULL`).
#' * `missingVarProp`: can be used instead of `missingVars`: The proportion of variables containing missing data (defaults to zero).
#' * `missingProp`: The proportion of missingness for variables containing missing data (defaults to zero), either a single value or a vector giving the probabilities for each variable.
#' * `missingMechanism`: The missing data mechanism, one of `'MCAR'` (the default), `'MAR'`, or `'NMAR'`.
#' 
#' `type = 'IG'` implements the independent generator approach (IG, Foldnes & Olsson, 2016) approach 
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.
#' 
#' `type = 'mnonr'` implements the approach suggested by Qu, Liu, & Zhang (2020) and requires provision of  Mardia's multivariate skewness (`skewness`)  and kurtosis  (`kurtosis`), where 
#' skewness must be non-negative and kurtosis must be at least 1.641 skewness + p (p + 0.774), where p is the number of variables.
#' 
#' `type = 'RK'` implements the approach suggested by Ruscio & Kaczetow (2008) and requires provision of the population distributions
#'  of each variable (`distributions`). `distributions` must be a list (if all variables shall be based on the same population distribution) or a list of lists. 
#'  Each component must specify the population distribution (e.g. `rchisq`) and additional arguments (`list(df = 2)`).
#' 
#' 
#' Foldnes, N. & Olsson, U. H. (2016) A Simple Simulation Technique for Nonnormal Data with Prespecified Skewness, Kurtosis, and Covariance Matrix. *Multivariate Behavioral Research, 51*, 207-219. doi: 10.1080/00273171.2015.1133274
#'
#' Qu, W., Liu, H., & Zhang, Z. (2020). A method of generating multivariate non-normal random numbers with desired multivariate skewness and kurtosis. *Behavior Research Methods, 52*, 939-946. doi: 10.3758/s13428-019-01291-5
#' 
#' Ruscio, J., & Kaczetow, W. (2008). Simulating multivariate nonnormal data using an iterative algorithm. *Multivariate Behavioral Research, 43*, 355-381. doi: 10.1080/00273170802285693
#'  
#' @examples
#' \dontrun{
#' # create Sigma and modelH0 using powerCFA
#' powerCFA <- semPower.powerCFA(type = 'a-priori', alpha = .05, beta = .05,
#'                               comparison = 'saturated',
#'                               Phi = .2, loadings = list(rep(.5, 3), rep(.7, 3)))
#'                               
#' # perform simulated power analysis using defaults       
#' simulate(modelH0 = powerCFA$modelH0, 
#'          Sigma = powerCFA$Sigma,
#'          N = powerCFA$power$requiredN,
#'          alpha = .05,
#'          simulatedPower = TRUE)
#'          
#' 
#' # same with additional options       
#' simulate(modelH0 = powerCFA$modelH0, 
#'          Sigma = powerCFA$Sigma,
#'          N = powerCFA$power$requiredN,
#'          alpha = .05,
#'          simulatedPower = TRUE, 
#'          simOptions = list(nReplications = 500, minConvergenceRate = .80))
#' 
#' 
#' # same with IG as data generation routine
#' simulate(modelH0 = powerCFA$modelH0, 
#'          Sigma = powerCFA$Sigma,
#'          N = powerCFA$power$requiredN,
#'          alpha = .05,
#'          simulatedPower = TRUE, 
#'          simOptions = list(type = 'IG', 
#'                            skewness = c(0, 1, -2, 6, 5, 4), 
#'                            kurtosis = c(-3, 6, 9, 0, 2, -2)))
#'                            
#'                            
#' # same with mnonr as data generation routine
#' simulate(modelH0 = powerCFA$modelH0, 
#'          Sigma = powerCFA$Sigma,
#'          N = powerCFA$power$requiredN,
#'          alpha = .05,
#'          simulatedPower = TRUE, 
#'          simOptions = list(type = 'mnonr', 
#'                            skewness = 1, 
#'                            kurtosis = 50))
#'                            
#'                            
#' # same with RK as data generation routine
#' distributions <- list(
#'   list('rnorm', list(mean = 0, sd = 10)),
#'   list('runif', list(min = 0, max = 1)),
#'   list('rbeta', list(shape1 = 1, shape2 = 2)),
#'   list('rexp', list(rate = 1)),
#'   list('rpois', list(lambda = 4)),
#'   list('rbinom', list(size = 1, prob = .5))
#' )
#' simulate(modelH0 = powerCFA$modelH0, 
#'          Sigma = powerCFA$Sigma,
#'          N = powerCFA$power$requiredN,
#'          alpha = .05,
#'          simulatedPower = TRUE, 
#'          simOptions = list(type = 'RK', 
#'                            distributions = distributions))
#' }
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats median
simulate <- function(modelH0 = NULL, modelH1 = NULL,
                     Sigma = NULL, mu = NULL, 
                     N = NULL,
                     alpha = NULL,
                     simOptions = list(
                       nReplications = 250, 
                       minConvergenceRate = .5,
                       type = 'normal',
                       missingVars = NULL,
                       missingVarProp = 0,
                       missingProp = 0,
                       missingMechanism = 'MCAR'
                     ),
                     lavOptions = NULL, lavOptionsH1 = lavOptions,
                     returnFmin = TRUE){
  
  
  nReplications <- ifelse(!is.null(simOptions[['nReplications']]), simOptions[['nReplications']], 250)
  minConvergenceRate <- ifelse(!is.null(simOptions[['minConvergenceRate']]), simOptions[['minConvergenceRate']], .5)
  
  if(is.list(Sigma)) nvar <- ncol(Sigma[[1]]) else nvar <- ncol(Sigma) 
  if(nvar >= min(unlist(N))) stop('Simulated power is not possible when the number of variables exceeds N.')
  if(min(unlist(N)) < 50 || 2*nvar >= min(unlist(N))) warning('In small N situations, simulated power will likely be inaccurate and yield high non-convergence rates.')
  checkBounded(nReplications, bound = c(10, 100000), inclusive = TRUE)
  checkBounded(minConvergenceRate, bound = c(0.01, 1), inclusive = TRUE)
  if(nReplications < 100) warning("Empirical power estimate with < 100 replications will be unreliable.")
  if(minConvergenceRate < .25) warning("Empirical power estimate allowing a low convergence rate might be unreliable.")
  # remaining checks are done in corresponding power fnc
  # data gen checks are done in gendata
  
  # warn in case of long computation times
  lavEstimators <- c("MLM", "MLMV", "MLMVS", "MLF", "MLR", "WLS", "DWLS", "WLSM", "WLSMV", "ULSM", "ULSMV")
  costlyEstm <- (!is.null(lavOptions[['estimator']]) && toupper(lavOptions[['estimator']]) %in% lavEstimators)
  projectedLong <- (costlyEstm && ncol(Sigma) > 50 && nReplications > 100) || (ncol(Sigma) > 100 && nReplications > 500) || nReplications > 10000
  if(projectedLong){
    mResp <- menu(c("Yes", "No"), title = "Simulated power with the specified model, the number of replications, and the type of estimator will presumably take a long time.\nDo you really want to proceed?")
    if(mResp != 1){
      stop("Simulated power aborted.")
    }else{
      cat("You have been warned.\n")
    }
  }
  
  
  # we need to call lavaan() directly with defaults as defined in sem()
  lavOptions <- append(getLavOptions(lavOptions, isCovarianceMatrix = FALSE, nGroups = length(Sigma)),
                       list(check.gradient = FALSE, check.post = FALSE, check.vcov = FALSE))
  if(!is.null(modelH1)){
    lavOptionsH1 <- append(getLavOptions(lavOptionsH1, isCovarianceMatrix = FALSE, nGroups = length(Sigma)),
                           list(check.gradient = FALSE, check.post = FALSE, check.vcov = FALSE))
  }
  
  # generate data 
  if(!is.list(Sigma)){
    simData <- genData(N = N, Sigma = Sigma, mu = mu, gIdx = NULL, nSets = nReplications / minConvergenceRate, modelH0 = modelH0, simOptions = simOptions)
  }else{
    simData <- lapply(seq_along(Sigma), function(x) genData(N = N[[x]], Sigma = Sigma[[x]], mu = mu[[x]], gIdx = x, nSets = nReplications / minConvergenceRate, modelH0 = modelH0, simOptions = simOptions))
  }

  # if we have missings, notify lav and check estimator 
  if(sum(is.na(simData[[1]])) > 0){
    lavOptions <- append(lavOptions, list(missing = 'ml'))
    if(!is.null(lavOptions[['estimator']]) && tolower(lavOptions[['estimator']]) %in% c('mlm', 'mlmv', 'mlmvs')) stop('Estimators mlm, mlmv, and mlmvs are not supported with missing data. Use mlr or mlf instead.')
    if(!is.null(modelH1)) lavOptionsH1 <- append(lavOptionsH1, list(missing = 'ml'))
  }
  
  efmin <- efminGroups <- list()
  rChiSq <- rLambda <- rPhi <- rPsi <- rBeta <- list()
  ePower <- 0
  r <- rr <- 1
  progress <- txtProgressBar(min = 0, max = nReplications, initial = 0, style = 3)
  while(r <= nReplications && rr <= nReplications / minConvergenceRate){
    setTxtProgressBar(progress, r)
    tryCatch({
      
      if(!is.list(Sigma)){
        cdata <- simData[[r]]
      }else{
        cdata <- rbind(simData[[1]][[r]], simData[[2]][[r]])
      }

      lavArgs <- list(model = modelH0, data = cdata)
      if(is.list(Sigma)) lavArgs <- append(lavArgs, list(group = 'gIdx'))
      lavresH0 <- do.call(lavaan::lavaan, append(lavArgs, lavOptions))

      if(lavresH0@optim[["converged"]]){
        
        # fmin by group
        cfminGroups <- NULL
        if(is.list(Sigma)) cfminGroups <- lavresH0@Fit@test[[lavresH0@Options[['test']]]][['stat.group']] / (unlist(N) - 1)
        
        cfmin <- 2 * lavaan::fitMeasures(lavresH0, 'fmin') # lav reports .5*fmin
        p <- lavaan::fitMeasures(lavresH0, 'pvalue')
        df <- lavaan::fitMeasures(lavresH0, 'df')
        if(lavresH0@Options[['test']] != 'standard'){
          p <- lavaan::fitMeasures(lavresH0, 'pvalue.scaled')
          df <- lavaan::fitMeasures(lavresH0, 'df.scaled')
        }

        # handle restricted comparison model 
        # (modelH1 must always get fit because sampling error does not allow just using modelH0 estm with different df)
        if(!is.null(modelH1)){
          
          lavArgs <- list(model = modelH1, data = cdata)
          if(is.list(Sigma)) lavArgs <- append(lavArgs, list(group = 'gIdx'))
          
          lavresH1 <- do.call(lavaan::lavaan, append(lavArgs, lavOptionsH1))

          if(lavresH1@optim[["converged"]]){
            
            mcomp <- lavaan::anova(lavresH0, lavresH1) 
            p <- mcomp$`Pr(>Chisq)`[2]
            df <- mcomp$`Df diff`[2]
            cfmin <- 2 * (lavaan::fitMeasures(lavresH0, 'fmin') - lavaan::fitMeasures(lavresH1, 'fmin'))
            if(is.list(Sigma)) cfminGroups - lavresH1@Fit@test[[lavresH0@Options[['test']]]][['stat.group']] / (unlist(N) - 1)
            
            # store param est
            cChi <- lavaan::fitMeasures(lavresH1, 'chisq')
            if(lavresH1@Options[['test']] != 'standard'){
              cChi <- lavaan::fitMeasures(lavresH1, 'chisq.scaled')
            }
            rChiSq <- append(rChiSq, cChi)
            cLambda <- lavresH1@Model@GLIST[which(names(lavresH1@Model@GLIST) %in% 'lambda')]
            rLambda <- append(rLambda, list(unlist(lapply(cLambda, function(x) x[x != 0])))) # only non-zero loadings  
            cPsi <- lavresH1@Model@GLIST[which(names(lavresH1@Model@GLIST) %in% 'psi')]
            rPsi <- append(rPsi, list(unlist(cPsi)))
            cBeta <- lavresH1@Model@GLIST[which(names(lavresH1@Model@GLIST) %in% 'beta')]
            rBeta <- append(rBeta, list(unlist(cBeta)))
            cPhi <- lavresH1@Model@GLIST[which(names(lavresH1@Model@GLIST) %in% 'phi')]
            rPhi <- append(rPhi, list(unlist(cPhi)))
          }
        }
      }

      if(lavresH0@optim[["converged"]] && (is.null(modelH1)|| lavresH1@optim[["converged"]])){
        efmin <- append(efmin, cfmin)
        efminGroups <- append(efminGroups, list(cfminGroups)) 
        
        if(p < alpha)
          ePower <- ePower + 1
        
        r <- r + 1 
      }

    }, warning = function(w) {
      # print(paste('WARNING: ',w))
    }, error = function(e) {
       print(paste('ERROR: ',e))
    })
    rr <- rr + 1
  }
  close(progress)
  
  
  if((r - 1) == 0) stop("Something went wrong during model estimation, no replication converged.")
  ePower <- ePower / (r - 1)
  if((r - 1) / (rr - 1) < minConvergenceRate){ 
    warning(paste("Actual convergence rate of", round((r - 1) / (rr - 1), 2), "is below minConvergenceRate of", minConvergenceRate, ". Results are based on", (r - 1),"replications."))
  }
  
  if(returnFmin){
    # sample fmin is biased, we need unbiased pop fmin
    ubFmean <- mean(unlist(efmin) - df / sum(unlist(N)))
    if(ubFmean <= 0) warning('Simulated estimate of F0 is zero or lower. Try to increase the number of replications.')
    efminGroups <- do.call(rbind, efminGroups)
    # we assume that each group contributes proportional df
    ubFmeanGroups <- NULL
    if(!is.null(efminGroups)) ubFmeanGroups <- unlist(lapply(1:ncol(efminGroups), function(x) mean( efminGroups[ ,x] - (df / length(N)) / N[[x]] ))) / length(N)

    out <- list(
      ePower = ePower,
      meanFmin = ubFmean,
      meanFminGroups = ubFmeanGroups,
      df = df,
      nrep = (r - 1),
      convergenceRate = (r - 1) / (rr - 1)
    )
    
    # also store mean param estm and pop params
    if(!is.null(modelH1)){
      
      # assume that modelH1 is properly specified, so that fitting modelH1 to Sigma
      # yields population parameters
      if(is.list(Sigma)) sample.nobs <- list(1000, 1000) else sample.nobs <- 1000
      lavOptionsH1[['estimator']] <- 'ML'  # needs to be overwritten in case this is set, since we are not working with observed variables
      lavOptionsH1[['missing']] <- NULL
      lavresPop <- do.call(lavaan::lavaan, 
                           append(list(model = modelH1, sample.cov = Sigma, sample.mean = mu, sample.nobs = sample.nobs, 
                                       sample.cov.rescale = FALSE),
                                  lavOptionsH1))
      if(lavresPop@optim[['fx']] > 1e-6) warning('H1 model is not properly specified.')

      # calc bias
      bChiSq <- (median(unlist(rChiSq)) - lavaan::fitMeasures(lavresPop, 'df')) / lavaan::fitMeasures(lavresPop, 'df')
      
      cLambda <- lavresPop@Model@GLIST[which(names(lavresPop@Model@GLIST) %in% 'lambda')]
      pLambda <- unlist(lapply(cLambda, function(x) x[x != 0]))
      bLambda <- mean( (apply(do.call(rbind, rLambda), 2, median) -  pLambda) / pLambda )
      
      bPhi <- bPsi <- bBeta <- NULL
      # for phi/psi/beta, we only consider population parameters that are larger than
      # a small constant to avoid absurd relative biases for parameters with true values close to zero 
      cPhi <- lavresPop@Model@GLIST[which(names(lavresPop@Model@GLIST) %in% 'phi')]
      if(length(cPhi) > 0){
        nonzero <- unlist(lapply(cPhi, function(x) which(abs(x) > .01)))
        pPhi <- unlist(lapply(cPhi, function(x) x[nonzero]))
        if(!is.null(pPhi)) bPhi <- mean( (apply(do.call(rbind, rPhi), 2, median) -  pPhi) / pPhi ) else bPhi <- NULL
      }

      cPsi <- lavresPop@Model@GLIST[which(names(lavresPop@Model@GLIST) %in% 'psi')]
      if(length(cPsi) > 0){
        nonzero <- unlist(lapply(cPsi, function(x) which(abs(x) > .01)))
        pPsi <- unlist(lapply(cPsi, function(x) x[nonzero]))
        if(!is.null(pPsi)) bPsi <- mean( (apply(do.call(rbind, rPsi)[, nonzero], 2, median) -  pPsi) / pPsi ) else bPsi <- NULL
      }
      
      cBeta <- lavresPop@Model@GLIST[which(names(lavresPop@Model@GLIST) %in% 'beta')]
      if(length(cBeta) > 0){
        nonzero <- unlist(lapply(cBeta, function(x) which(abs(x) > .01)))
        pBeta <- unlist(lapply(cBeta, function(x) x[nonzero]))
        if(!is.null(pBeta)) bBeta <- mean( (apply(do.call(rbind, rBeta)[, nonzero], 2, median) -  pBeta) / pBeta ) else bBeta <- NULL
      }
      
      out <- append(out, list(
        bChiSq = bChiSq,
        bLambda = bLambda,
        bPhi = bPhi,
        bBeta = bBeta,
        bPsi = bPsi
      ))
    }
    
    out
    
  }else{
    ePower
  }
}


#' genData
#' 
#' Generates random data from population variance-covariance matrix and population means, either
#' from a multivariate normal distribution, or using one of various approaches to generate 
#' non-normal data.
#'
#' @param type one of 'normal' (the default), 'IG', 'mnonr', or 'RK'.
#' @param N sample size.
#' @param Sigma population covariance matrix.
#' @param mu population means.
#' @param nSets number of data sets to generate
#' @param gIdx if not `NULL`, add gIdx as numeric group index as additional variable to generated data
#' @param modelH0 a `lavaan` model string, only used to determine the number of factors when `type = 'RK'`
#' @param simOptions additional arguments passed to specific data generation routine
#' @return Returns the generated data
#' @examples
#' \dontrun{
#' gen <- semPower.genSigma(Phi = .2, loadings = list(rep(.5, 3), rep(.7, 3)))
#' data <- genData(N = 500, Sigma = gen$Sigma) 
#' }
#' @importFrom stats rbinom ecdf 
genData <- function(type = 'normal', 
                    N = NULL, Sigma = NULL, mu = NULL, 
                    nSets = 1, gIdx = NULL, modelH0 = NULL, simOptions = NULL){
  
  type <- checkDataGenerationTypes(type)
  
  missingMechanism <- ifelse(!is.null(simOptions[['missingMechanism']]), simOptions[['missingMechanism']], 'mcar')
  missingMechanism <- checkMissingTypes(missingMechanism)
  if(!is.null(simOptions[['missingVarsProp']]) && !is.null(simOptions[['missingVars']])) stop('Either set missingVarsProp or set missingVars, but not both.')
  missingVarsProp <- ifelse(!is.null(simOptions[['missingVarsProp']]), simOptions[['missingVarsProp']], 0)
  if(length(missingVarsProp) > 1) stop('missingVarsProp must be a single number.')
  checkBounded(missingVarsProp, inclusive = TRUE)
  if(!is.null(simOptions[['missingVars']])){
    missingVars <- simOptions[['missingVars']] 
  }else{
    missingVars <- sample(seq(ncol(Sigma)), round(missingVarsProp*(ncol(Sigma)))) 
  } 
  if(!is.null(simOptions[['missingProp']])) missingProp  <- simOptions[['missingProp']] else missingProp <- 0  
  if(length(missingProp) == 1) missingProp <- rep(missingProp, length(missingVars))
  if(length(missingProp) != length(missingVars)) stop('Either specify a single value for missingProp or define a missingProp for each missingvar.')
  lapply(missingProp, function(x) checkBounded(x, 'missingProp must lie', inclusive = TRUE))
  if((any(missingProp > 0) && length(missingVars) == 0) || (any(missingProp == 0) && length(missingVars) > 0)) stop('missingProp and either missingVarProp or missingVars must be larger than zero to produce missings.')

  switch(type,
         # normal by cholesky decomposition
         normal = {
           rdat <- genData.normal(N, Sigma, nSets)
         },
         # non-normal using IG
         ig = {
           rdat <- genData.IG(N, Sigma, nSets, simOptions[['skewness']], simOptions[['kurtosis']])
         },
         # non-normal using mnonr
         mnonr = {
           rdat <- genData.mnonr(N, Sigma, nSets, simOptions[['skewness']], simOptions[['kurtosis']])
         },
         # non-normal using ruscio
         rk = {
           rdat <- genData.RK(N, Sigma, nSets, simOptions[['distributions']], modelH0)
         }
  )  
  
  # add mu
  if(!is.null(mu)){
    rdat <- lapply(rdat, function(x) x + matrix(t(rep(mu, N)), ncol = ncol(Sigma), byrow = TRUE))
  }
  
  ### missings
  # implementation based on https://psyarxiv.com/rq6yb/
  # always produces maximum number of missing patterns
  if(any(missingProp > 0) && missingVarsProp > 0){

    rdat <- lapply(rdat, function(x){
      
      if(missingMechanism == 'mcar'){
        
        for(i in seq(missingVars)){
          NArows <- as.logical(rbinom(nrow(x), 1, missingProp[i]))
          x[NArows, missingVars[i]] <- NA
        }

      }else if(missingMechanism == 'mar'){
        
        # pick a conditioning variable
        if(length(missingVars) < ncol(Sigma)){
          conditioningVar <- sample(seq(ncol(Sigma))[!seq(ncol(Sigma)) %in% missingVars], 1)
        }else{
          warning('For MAR, one variable acts as conditioning variable and must not contain missings. Removing a single variable from missingVars.')
          conditioningVar <- sample(ncol(Sigma), 1)
          missingVars <- missingVars[-conditioningVar]
        }

        pecdf <- ecdf(x[ , conditioningVar])
        for(i in seq(missingVars)){
          NArows <- logical()
          for(j in 1:nrow(x)){
            percentile <- pecdf(x[j, conditioningVar])
            NArows[j] <- as.logical(rbinom(1, 1, prob = 2*missingProp[i]*percentile))
          }
          x[NArows, missingVars[i]] <- NA
        }
        
      }else if(missingMechanism == 'nmar'){
        
        for(i in seq(missingVars)){
          pecdf <- ecdf(x[ , missingVars[i]])
          NArows <- logical()
          for(j in 1:nrow(x)){
            percentile <- pecdf(x[j, missingVars[i]])
            NArows[j] <- as.logical(rbinom(1, 1, prob = 2*missingProp[i]*percentile))
          }
          x[NArows, missingVars[i]] <- NA
        }
        
      }
      
      x
      
    })
    
  }
  
  # add gidx
  if(!is.null(gIdx)){
    rdat <- lapply(rdat, function(x) cbind(x, matrix(rep(gIdx, N), ncol = 1, dimnames = list(NULL, c('gIdx')))) )
  }

  rdat 
}


#' genData.normal
#' 
#' Generates multivariate normal random data from population variance-covariance matrix.
#'
#' @param N sample size.
#' @param Sigma population covariance matrix.
#' @param nSets number of data sets to generate
#' @return Returns the generated data
#' @importFrom stats rnorm 
genData.normal <- function(N = NULL, Sigma = NULL, nSets = 1){
  lapply(seq(nSets), function(x){
    randomData <- matrix(rnorm(N * ncol(Sigma)), N, ncol(Sigma)) 
    t(t(chol(Sigma)) %*% t(randomData))
  })
}


#' genData.IG
#' 
#' Generates random data from population variance-covariance matrix using 
#' the independent generator approach (IG, Foldnes & Olsson, 2016) approach 
#' specifying third and fourth moments of the marginals.
#' 
#' This function is a wrapper for the respective function of the ´covsim´ package. 
#' 
#' For details, see 
#' Foldnes, N. & Olsson, U. H. (2016) A Simple Simulation Technique for Nonnormal Data with Prespecified Skewness, Kurtosis, and Covariance Matrix. *Multivariate Behavioral Research, 51*, 207-219. 10.1080/00273171.2015.1133274
#'
#' @param N sample size.
#' @param Sigma population covariance matrix.
#' @param nSets number of data sets to generate
#' @param skewness vector specifying skewness for each variable 
#' @param kurtosis vector specifying excess kurtosis for each variable
#' @return Returns the generated data
#' @importFrom utils installed.packages
genData.IG <- function(N = NULL, Sigma = NULL, nSets = 1,  
                       skewness = NULL, kurtosis = NULL){
  
  if(!'covsim' %in% rownames(installed.packages())) stop('Generation of non-normal random data using IG requires the covsim package, so install covsim first.')
  if(is.null(skewness) || is.null(kurtosis)) stop('skewness and kurtosis must not be NULL.')
  if(length(skewness) != ncol(Sigma)) stop('skewness must match ncol(Sigma), i.e., must be specified for each variable ')
  if(length(kurtosis) != ncol(Sigma)) stop('kurtosis must match ncol(Sigma), i.e., must be specified for each variable ')
  
  covsim::rIG(N = N, sigma = Sigma,
              skewness = skewness,
              excesskurt = kurtosis,
              reps = nSets)
  
}


#' genData.mnonr
#' 
#' Generates random data from population variance-covariance matrix using 
#' the approach by Qu, Liu, & Zhang (2020) specifying Mardia's multivariate skewness and kurtosis.
#' 
#' This function is a wrapper for the respective function of the `mnonr` package. 
#' 
#' For details, see 
#' Qu, W., Liu, H., & Zhang, Z. (2020). A method of generating multivariate non-normal random numbers with desired multivariate skewness and kurtosis. *Behavior Research Methods, 52*, 939-946. doi: 10.3758/s13428-019-01291-5
#'
#' @param N sample size.
#' @param Sigma population covariance matrix.
#' @param nSets number of data sets to generate
#' @param skewness multivariate skewness. May not be negative.
#' @param kurtosis multivariate kurtosis. Must be >= 1.641 skewness + p (p + 0.774), where p is the number of variables.
#' @return Returns the generated data
#' @importFrom utils installed.packages
genData.mnonr <- function(N = NULL, Sigma = NULL, nSets = 1,  
                       skewness = NULL, kurtosis = NULL){
  
  if(!'mnonr' %in% rownames(installed.packages())) stop('Generation of non-normal random data using mnonr requires the mnonr package, so install mnonr first.')
  if(is.null(skewness) || is.null(kurtosis)) stop('skewness and kurtosis must not be NULL.')
  if(length(skewness) != 1) stop('multivariate skewness must be a single number.')
  if(length(kurtosis) != 1) stop('multivariate kurtosis must be a single number.')
  
  lapply(seq(nSets), function(x){
    rd <- mnonr::mnonr(n = N, Sigma = Sigma, p = ncol(Sigma), 
                 ms = skewness, 
                 mk = kurtosis)
    colnames(rd) <- paste0('x', 1:ncol(Sigma))
    rd
  })
  
}

#' genData.RK
#' 
#' Generates random data from population variance-covariance matrix using 
#' the approach by Ruscio & Kaczetow (2008)
#' specifying distributions for the marginals.
#' 
#' This function is based on the implementation by Ruscio & Kaczetow (2008).
#' 
#' For details, see 
#' Ruscio, J., & Kaczetow, W. (2008). Simulating multivariate nonnormal data using an iterative algorithm. *Multivariate Behavioral Research, 43*, 355-381.
#'
#' @param N sample size.
#' @param Sigma population covariance matrix.
#' @param nSets number of data sets to generate
#' @param distributions a list specifying the population distribution and  additional arguments in a list either to apply to all variables (e.g. `list(rchisq, list(df = 2))`) or a list of lists specifying the distributions for each variable. See examples.
#' @param modelH0 a `lavaan` model string, only used to determine the number of factors.
#' @param maxIter maximum number of iterations, defaults to 10.
#' @return Returns the generated data
#' @examples
#' \dontrun{
#' distributions <- list(
#'   list('rchisq', list(df = 2)),
#'   list('runif', list(min = 0, max = 1)),
#'   list('rexp', list(rate = 1))
#' )
#' data <- genData.ruscio(N = 100, Sigma = diag(3),
#'                        distributions = distributions, 
#'                        modelH0 = 'f =~ x1 + x2 + x3')
#'                        
#' distributions <- list(
#'   list('rnorm', list(mean = 0, sd = 10)),
#'   list('runif', list(min = 0, max = 1)),
#'   list('rbeta', list(shape1 = 1, shape2 = 2)),
#'   list('rexp', list(rate = 1)),
#'   list('rpois', list(lambda = 4)),
#'   list('rbinom', list(size = 1, prob = .5))
#' )
#' data <- genData.ruscio(N = 100, Sigma = diag(6),
#'                        distributions = distributions, 
#'                        modelH0 = 'f1=~x1+x2+x3\nf2=~x4+x5+x6')
#' 
#'}
#' @importFrom stats factanal cov2cor cor
genData.RK <- function(N = NULL, Sigma = NULL, nSets = 1,
                           distributions = NULL, modelH0 = NULL, maxIter = 10){
  
  if(!is.list(distributions)) stop('distributions must be a list.')
  if(is.list(distributions) && !is.list(distributions[[1]])) distributions <- (rep(list(distributions), ncol(Sigma)))
  lapply(distributions, function(x){
    if(length(x) != 2 || !is.character(x[[1]]) || !is.list(x[[2]])) 
      stop('Each component of distributions must contain two parts, the first specifying the distribution (e.g. rnorm) and the second a (perhaps empty) list specifying the arguments, but omitting n (e.g. list(mean = 0, sd = 1))')
  }) 
   
  getLoadings <- function(Sigma, nFactors){
    fa <- factanal(covmat = Sigma, factors = nFactors)
    loadings <- matrix(fa$loadings, ncol = nFactors)
    loadings[loadings > 1] <- 1
    loadings[loadings < -1] <- -1
    if (loadings[1, 1] < 0) loadings <- loadings * -1
    loadings
  }
  
  # determine nFactors through model string; this could lead to incorrect results, 
  # but we want to avoid both parallel analysis and the need to supply nFactors directly 
  tok <- as.list(strsplit(modelH0, split = '\n', fixed = TRUE)[[1]])
  nFactors <- sum(unlist(lapply(tok, function(x) if(grepl('=~', x)) length(strsplit(x, split = '+', fixed = TRUE)[[1]]) > 1 )))
  nFactors <- max(1, nFactors)
  
  Sigma <- cov2cor(Sigma)
  k <- ncol(Sigma)
  rdat <- matrix(0, nrow = N, ncol = k) 
  
  nnDistributions <- lapply(distributions, function(x){
    sort(do.call(x[[1]], append(list(n = N), x[[2]])))
  })
  nnDistributions <- t(do.call(rbind, nnDistributions))

  # Generate random normal data for shared and unique components
  lapply(seq(nSets), function(x){
    Shared.Comp <- matrix(rnorm(N * nFactors, 0, 1), nrow = N, ncol = nFactors)
    Unique.Comp <- matrix(rnorm(N * k, 0, 1), nrow = N, ncol = k)
    residuals <- matrix(0, nrow = k, ncol = 1)
    
    # find best intermediate correlation matrix
    RMSRBest <- 1 
    SigmaInterm <- Sigma
    iter <- 0 
    while(iter < maxIter) {
      iter <- iter + 1
      # Calculate factor loadings and apply to reproduce desired correlations (steps 7, 8)
      loadings <- getLoadings(Sigma, nFactors)
      residuals <- 1 - rowSums(loadings^2)
      residuals[residuals < 0] <- 0
      residuals <- sqrt(residuals)
      for(i in 1:k){
        rdat[ ,i] <- (Shared.Comp %*% t(loadings))[ ,i] + Unique.Comp[ ,i] * residuals[i]
      }
      
      # Replace normal with nonnormal distributions (step 9)
      for(i in 1:k) {
        rdat <- rdat[sort.list(rdat[,i]),]
        rdat[,i] <- nnDistributions[,i]
      }
      # Calculate RMSR correlation, compare to lowest value, take appropriate action (steps 10, 11, 12)
      ResCor <- Sigma - cor(rdat)
      RMSR <- sqrt(sum(ResCor[lower.tri(ResCor)] * ResCor[lower.tri(ResCor)]) / (.5 * (k * k - k)))
      
      if(RMSR < RMSRBest){
        RMSRBest <- RMSR
        SigmaBest <- SigmaInterm
        ResCorBest <- ResCor
        SigmaInterm <- SigmaInterm + ResCor
        iter <- 0
      }else{
        iter <- iter + 1
        SigmaInterm <- SigmaBest + .5 ^ iter * ResCorBest
      }
    }
    
    # Construct the data set with the lowest RMSR
    loadings <- getLoadings(SigmaBest, nFactors)
    residuals <- 1 - rowSums(loadings^2)
    residuals[residuals < 0] <- 0
    residuals <- sqrt(residuals)
    for(i in 1:k){
      rdat[ ,i] <- (Shared.Comp %*% t(loadings))[ ,i] + Unique.Comp[ ,i] * residuals[i]
    }
    rdat <- apply(rdat, 2, scale) 
    for(i in 1:k) {
      rdat <- rdat[sort.list(rdat[ ,i]),]
      rdat[ ,i] <- nnDistributions[ ,i]
    }
    colnames(rdat) <- paste0('x', 1:ncol(Sigma))
    rdat
  })
  
}


