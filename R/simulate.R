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
#' * `nReplications`: The targeted number of valid simulation runs, defaults to 500.
#' * `minConvergenceRate`:  The minimum convergence rate required, defaults to .75. The maximum actual simulation runs are increased by a factor of 1/minConvergenceRate.
#' * `type`: specifies whether the data should be generated from a population assuming multivariate normality (`'normal'`; the default), or based on an approach generating non-normal data (`'IG'`, `'mnonr'`, `'RK'`, or `'VM'`). 
#' The approaches generating non-normal data require additional arguments detailed below.
#' * `missingVars`: vector specifying the variables containing missing data (defaults to `NULL`).
#' * `missingVarProp`: can be used instead of `missingVars`: The proportion of variables containing missing data (defaults to zero).
#' * `missingProp`: The proportion of missingness for variables containing missing data (defaults to zero), either a single value or a vector giving the probabilities for each variable.
#' * `missingMechanism`: The missing data mechanism, one of `'MCAR'` (the default), `'MAR'`, or `'NMAR'`.
#' * `nCores`: The number of cores to use for parallel processing. Defaults to 1 (= no parallel processing). This requires the `doFuture` package.
#' * `futureStrategy`: A string specifying the strategy how to resolve a future when `nCores` is larger than 1. Defaults to `'multisession'`. This is passed to the `plan` method of the `doFuture` package. See the `doFuture` package for valid strategies.
#' 
#' `type = 'IG'` implements the independent generator approach (IG, Foldnes & Olsson, 2016) approach 
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors. This requires the `covsim` package.
#' 
#' `type = 'mnonr'` implements the approach suggested by Qu, Liu, & Zhang (2020) and requires provision of  Mardia's multivariate skewness (`skewness`)  and kurtosis  (`kurtosis`), where 
#' skewness must be non-negative and kurtosis must be at least 1.641 skewness + p (p + 0.774), where p is the number of variables. This requires the `mnonr` package.
#' 
#' `type = 'RK'` implements the approach suggested by Ruscio & Kaczetow (2008) and requires provision of the population distributions
#'  of each variable (`distributions`). `distributions` must be a list (if all variables shall be based on the same population distribution) or a list of lists. 
#'  Each component must specify the population distribution (e.g. `rchisq`) and additional arguments (`list(df = 2)`).
#' 
#' `type = 'VM'` implements the third-order polynomial method (Vale & Maurelli, 1983) 
#' specifying third and fourth moments of the marginals, and thus requires that skewness (`skewness`) and excess kurtosis (`kurtosis`) for each variable are provided as vectors.
#' 
#' 
#' Foldnes, N. & Olsson, U. H. (2016) A Simple Simulation Technique for Nonnormal Data with Prespecified Skewness, Kurtosis, and Covariance Matrix. *Multivariate Behavioral Research, 51*, 207-219. doi: 10.1080/00273171.2015.1133274
#'
#' Qu, W., Liu, H., & Zhang, Z. (2020). A method of generating multivariate non-normal random numbers with desired multivariate skewness and kurtosis. *Behavior Research Methods, 52*, 939-946. doi: 10.3758/s13428-019-01291-5
#' 
#' Ruscio, J., & Kaczetow, W. (2008). Simulating multivariate nonnormal data using an iterative algorithm. *Multivariate Behavioral Research, 43*, 355-381. doi: 10.1080/00273170802285693
#' 
#' Vale, C. & Maurelli, V. (1983). Simulating multivariate nonnormal distributions. *Psychometrika, 48*, 465-471.
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
#'          N = powerCFA$requiredN,
#'          alpha = .05,
#'          simulatedPower = TRUE)
#'          
#' 
#' # same with additional options       
#' simulate(modelH0 = powerCFA$modelH0, 
#'          Sigma = powerCFA$Sigma,
#'          N = powerCFA$requiredN,
#'          alpha = .05,
#'          simulatedPower = TRUE, 
#'          simOptions = list(nReplications = 500, 
#'                            minConvergenceRate = .80, 
#'                            nCores = 8))
#' 
#' 
#' # same with IG as data generation routine
#' simulate(modelH0 = powerCFA$modelH0, 
#'          Sigma = powerCFA$Sigma,
#'          N = powerCFA$requiredN,
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
#'          N = powerCFA$requiredN,
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
#'          N = powerCFA$requiredN,
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
                       nReplications = 500, 
                       minConvergenceRate = .75,
                       type = 'normal',
                       missingVars = NULL,
                       missingVarProp = 0,
                       missingProp = 0,
                       missingMechanism = 'MCAR',
                       nCores = 1,
                       futureStrategy = 'multisession'
                     ),
                     lavOptions = NULL, lavOptionsH1 = lavOptions,
                     returnFmin = TRUE){

  nCores <- ifelse(!is.null(simOptions[['nCores']]), simOptions[['nCores']], 1)
  if(nCores > 1 && !'doFuture' %in% rownames(installed.packages())) stop('Parallel processing requires the doFuture package, so install doFuture first.')    
  if(nCores > 1 && !'progressr' %in% rownames(installed.packages())) stop('Parallel processing requires the progressr package, so install progressr first.')    

  nReplications <- ifelse(!is.null(simOptions[['nReplications']]), simOptions[['nReplications']], 500)
  minConvergenceRate <- ifelse(!is.null(simOptions[['minConvergenceRate']]), simOptions[['minConvergenceRate']], .75)
  maxReplications <- ceiling(nReplications / minConvergenceRate)
  
  if(is.list(Sigma)) nvar <- ncol(Sigma[[1]]) else nvar <- ncol(Sigma) 
  if(nvar >= min(unlist(N))) stop('Simulated power is not possible when the number of variables exceeds N.')
  if(min(unlist(N)) < 50 || 2*nvar >= min(unlist(N))) warning('In small N situations, simulated power will likely be inaccurate and yield high non-convergence rates.')
  checkBounded(nReplications, bound = c(10, 100000), inclusive = TRUE)
  checkBounded(minConvergenceRate, bound = c(0.01, 1), inclusive = TRUE)
  if(nReplications < 100) warning("Empirical power estimate with < 100 replications will be unreliable.")
  if(minConvergenceRate < .25) warning("Empirical power estimate allowing a low convergence rate might be unreliable.")
  # remaining checks are done in corresponding power fnc
  # data gen checks are done in gendata
  
  # don't allow ULS estimation (but ULSM, ULSMV are ok)
  if(!is.null(lavOptions[['estimator']]) && lavOptions[['estimator']] == 'ULS') stop('The ULS test-statistic has no known asymptotic distribution. You can try ULSM or ULSMV.')
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
  
  # generate data. we always generate maxreplications, though this incurs a bit overhead 
  if(!is.list(Sigma)){
    simData <- genData(N = N, Sigma = Sigma, mu = mu, gIdx = NULL, nSets = maxReplications, modelH0 = modelH0, simOptions = simOptions)
  }else{
    simData <- lapply(seq_along(Sigma), function(x) genData(N = N[[x]], Sigma = Sigma[[x]], mu = mu[[x]], gIdx = x, nSets = maxReplications, modelH0 = modelH0, simOptions = simOptions))
  }

  # if we have missings, notify lav and check estimator 
  if(sum(is.na(unlist(simData[[1]]))) > 0){
    lavOptions <- append(lavOptions, list(missing = 'ml'))
    if(!is.null(lavOptions[['estimator']]) && tolower(lavOptions[['estimator']]) %in% c('mlm', 'mlmv', 'mlmvs')) stop('Estimators mlm, mlmv, and mlmvs are not supported with missing data. Use mlr or mlf instead.')
    if(!is.null(modelH1)) lavOptionsH1 <- append(lavOptionsH1, list(missing = 'ml'))
  }
  
  # parallel
  if(nCores > 1){
    `%dofuture%` <- doFuture::`%dofuture%`
    future::plan(futureStrategy, workers = nCores)
    progressr::with_progress({
      p <- progressr::progressor(along = seq(nReplications)) # progressbar
      res <- foreach::foreach(r = seq(nReplications), .options.future = list(seed = TRUE)) %dofuture% {
        p()   # update progress bar
        doSim(r = r, 
              simData = simData,
              isMultigroup = is.list(Sigma),
              modelH0 = modelH0, modelH1 = modelH1, 
              lavOptions = lavOptions, lavOptionsH1 = lavOptionsH1)
      }
    })
    future::plan(future::sequential()) # explicitly close multisession workers

    fit <- lapply(res, '[[', 1)
    # replace non-converged by NA
    fit[which(unlist(lapply(fit, function(x) length(x) == 0)))] <- list(rep(list(rep(NA, 5)), 3))
    nConverged <- sum(!is.na(do.call(rbind, lapply(fit, '[[', 1))[, 1]))
    
    rr <- 1
  }else{
    # single core case
    res <- list()
    nConverged <- 0
  }
  
  
  ## check convergence and do additional non-parallel runs to reach target replications
  if(nConverged < nReplications){
    
    # fit models
    progressBar <- txtProgressBar(min = 0, max = (nReplications - nConverged), initial = 0, style = 3)
    r <- rr <- (nConverged + 1)
    while(r <= nReplications && rr <= maxReplications){
      setTxtProgressBar(progressBar, (r - nConverged))
      cr <- doSim(r = r, 
                  simData = simData,
                  isMultigroup = is.list(Sigma),
                  modelH0 = modelH0, modelH1 = modelH1,
                  lavOptions = lavOptions, lavOptionsH1 = lavOptionsH1)
      # [fit][fitH0][fmin] 
      if(!is.null(cr[[1]][[1]][1]) && !is.na(cr[[1]][[1]][1])){
        res <- append(res, list(cr))
        r <- r + 1
      }
      rr <- rr + 1
    }
    close(progressBar)
  }

  # check convergence
  fit <- lapply(res, '[[', 1)
  # replace non-converged by NA
  fit[which(unlist(lapply(fit, function(x) length(x) == 0)))] <- list(rep(list(rep(NA, 5)), 3))
  nConverged <- sum(!is.na(do.call(rbind, lapply(fit, '[[', 1))[, 1]))
  convergenceRate <- nConverged / length(res)
  if(nConverged == 0) stop("Something went wrong during model estimation, no replication converged.")
  if(convergenceRate < minConvergenceRate){ 
    warning(paste("Actual convergence rate of", round(convergenceRate, 2), "is below minConvergenceRate of", minConvergenceRate, ". Results are based on", nConverged,"replications."))
  }
  
  ## eval res
  fitH0 <- do.call(rbind, lapply(fit, '[[', 1))
  fitDiff <- do.call(rbind, lapply(fit, '[[', 3))
  colnames(fitH0) <- colnames(fitDiff) <- c('fmin', 'chisq', 'df', 'p', paste0('fminGroup', seq(length(N))))
  
  fitH0 <- fitH0[!is.na(fitH0[ ,'fmin']), ]
  fitDiff <- fitDiff[!is.na(fitDiff[ ,'fmin']), ]
  
  rrH0 <- sum(fitH0[, 'p'] < alpha) / nConverged # same as epower when no h1 model is provided
  ePower <- sum(fitDiff[, 'p'] < alpha) / nConverged
  
  if(!returnFmin){
    
    ePower
    
  }else{
    
    # sample fmin is biased, we need unbiased pop fmin
    df <- fitDiff[1, 'df']
    ubFmean <- mean(fitDiff[ , 'fmin'] - df / sum(unlist(N)))
    if(ubFmean <= 0) warning('Simulated estimate of F0 is zero or lower. Try to increase the number of replications.')
    
    efminGroups <- fitDiff[, paste0('fminGroup', seq(length(N)))]
    # we assume that each group contributes proportional df
    ubFmeanGroups <- NULL
    if(!is.na(efminGroups[1])) ubFmeanGroups <- unlist(lapply(1:ncol(efminGroups), function(x) mean( efminGroups[ , x] - (df / length(N)) / N[[x]] ))) / length(N)

    outSimData <- simData
    if(is.list(Sigma)) outSimData <- lapply(seq(nReplications), function(r) do.call(rbind, lapply(simData, '[[', r )))
    
    out <- list(
      ePower = ePower,
      meanFmin = ubFmean,
      meanFminGroups = ubFmeanGroups,
      df = df,
      dfH0 = fitH0[1, 'df'],
      nrep = nConverged,
      convergenceRate = convergenceRate,
      chiSqH0 = fitH0[, 'chisq'],
      chiSqDiff = fitDiff[, 'chisq'],
      rrH0 = rrH0,
      simRes = list(
        simData = outSimData,
        fitH0 = fitH0
      )
    )
    
    
    if(!is.null(modelH1)){
      fitH1 <- do.call(rbind, lapply(fit, '[[', 2))
      fitH1 <- fitH1[!is.na(fitH1[ ,'fmin']), ]
      
      colnames(fitH1) <- colnames(fitH0) 
      
      # eval param bias
      param <- lapply(res, '[[', 2)
      # only consider converged param est
      param <- param[unlist(lapply(param, function(x) length(x) > 0))]  
      rLambda <- do.call(rbind, lapply(param, '[[', 1))
      rPsi <- do.call(rbind, lapply(param, '[[', 2))
      rBeta <- do.call(rbind, lapply(param, '[[', 3))
      rPhi <- do.call(rbind, lapply(param, '[[', 4))
      
      # assume that modelH1 is properly specified, so that fitting modelH1 to Sigma
      # yields population parameters
      lavOptionsH1[['estimator']] <- 'ML'  # needs to be overwritten in case this is set, since we are not working with observed variables
      lavOptionsH1[['missing']] <- NULL
      lavresPop <- do.call(lavaan::lavaan, 
                           append(list(model = modelH1, sample.cov = Sigma, sample.mean = mu, sample.nobs = N, 
                                       sample.cov.rescale = FALSE),
                                  lavOptionsH1))
      if(lavresPop@optim[['fx']] > 1e-6) warning('H1 model is not properly specified.')
      
      # h1 chi bias. we compute this here, but cant do this for H0/diff because need the ncp
      bChiSqH1 <- ksChiSqH1 <- rrH1 <- NA 
      if(fitH1[1, 'df'] > 0 ){
        bChiSqH1 <- (mean(fitH1[, 'chisq']) - fitH1[1, 'df']) / fitH1[1, 'df']
        ksChiSqH1 <- getKSdistance(fitH1[, 'chisq'], fitH1[1, 'df'])  
        rrH1 <- sum(fitH1[, 'p'] < alpha) / nConverged 
      }

      # parameter bias
      cLambda <- lavresPop@Model@GLIST[which(names(lavresPop@Model@GLIST) %in% 'lambda')]
      pLambda <- unlist(lapply(cLambda, function(x) x[x != 0]))
      bLambda <- mean( (apply(rLambda, 2, median) -  pLambda) / pLambda )
      
      bPhi <- bPsi <- bBeta <- NULL
      # for phi/psi/beta, we only consider population parameters that are larger than
      # a small constant to avoid absurd relative biases for parameters with true values close to zero 
      cPhi <- lavresPop@Model@GLIST[which(names(lavresPop@Model@GLIST) %in% 'phi')]
      if(length(cPhi) > 0){
        nonzero <- unlist(lapply(cPhi, function(x) which(abs(x) > .01)))
        pPhi <- unlist(lapply(cPhi, function(x) x[nonzero]))
        if(!is.null(pPhi)) bPhi <- mean( (apply(rPhi, 2, median) - pPhi) / pPhi ) else bPhi <- NULL
      }
      
      cPsi <- lavresPop@Model@GLIST[which(names(lavresPop@Model@GLIST) %in% 'psi')]
      if(length(cPsi) > 0){
        nonzero <- unlist(lapply(cPsi, function(x) which(abs(x) > .01)))
        pPsi <- unlist(lapply(cPsi, function(x) x[nonzero]))
        if(!is.null(pPsi)) bPsi <- mean( (apply(rPsi[, nonzero], 2, median) - pPsi) / pPsi ) else bPsi <- NULL
      }
      
      cBeta <- lavresPop@Model@GLIST[which(names(lavresPop@Model@GLIST) %in% 'beta')]
      if(length(cBeta) > 0){
        nonzero <- unlist(lapply(cBeta, function(x) which(abs(x) > .01)))
        pBeta <- unlist(lapply(cBeta, function(x) x[nonzero]))
        if(!is.null(pBeta)) bBeta <- mean( (apply(rBeta[, nonzero], 2, median) - pBeta) / pBeta ) else bBeta <- NULL
      }    
      
      out <- append(out, list(
        bChiSqH1 = bChiSqH1,
        ksChiSqH1 = ksChiSqH1,
        rrH1 = rrH1,
        bLambda = bLambda,
        bPhi = bPhi,
        bBeta = bBeta,
        bPsi = bPsi
      ))
      out[['simRes']] <- append(out[['simRes']], 
                             list(
                               fitH1 = fitH1, 
                               fitDiff = fitDiff,
                               paramEst = param
                             ))
      
    }
    
    out
  }
  
}

#' doSim
#' 
#' Performs simulation on single data set
#'
#' @param r replication id
#' @param simData list of datafiles
#' @param isMultigroup multigroup flag
#' @param modelH0 `lavaan` model string defining the (incorrect) analysis model.
#' @param modelH1 `lavaan` model string defining the comparison model. If omitted, the saturated model is the comparison model.
#' @param lavOptions a list of additional options passed to `lavaan`, e. g., `list(estimator = 'mlm')` to request robust ML estimation
#' @param lavOptionsH1 lavoptions when fitting `modelH1`. If `NULL`, the same as `lavOptions`.
#' @return list 
doSim <- function(r, 
                  simData,
                  isMultigroup = FALSE,
                  modelH0, modelH1, 
                  lavOptions, lavOptionsH1){
  
  cres <- list(list(), list())
  
  tryCatch({
    
    if(!isMultigroup){
      cdata <- simData[[r]]
    }else{
      cdata <- do.call(rbind, lapply(simData, '[[', r ))
    }
    
    lavArgs <- list(model = modelH0, data = cdata)
    if(isMultigroup) lavArgs <- append(lavArgs, list(group = 'gIdx'))
    lavresH0 <- do.call(lavaan::lavaan, append(lavArgs, lavOptions))
    
    if(lavresH0@optim[["converged"]]){
      
      lavfitH0 <- lavaan::fitMeasures(lavresH0)
      fitH0 <- lavfitH0[c('fmin','chisq', 'df', 'pvalue')]
      if(lavresH0@Options[['test']] != 'standard'){
        fitH0 <- lavfitH0[c('fmin','chisq.scaled', 'df.scaled', 'pvalue.scaled')]
      }
      fitH0['fmin'] <- 2*fitH0['fmin'] # lav reports .5*fmin
      
      # fmin by group
      cfminGroupsH0 <- NA
      if(isMultigroup) cfminGroupsH0 <- lavresH0@Fit@test[[lavresH0@Options[['test']]]][['stat.group']] / (unlist(lavresH0@Data@nobs) - 1)
      
      fitH0 <- c(fitH0, cfminGroupsH0)
      fitDiff <- fitH0
      
      # handle restricted comparison model 
      # (modelH1 must always get fit because sampling error does not allow just using modelH0 estm with different df)
      if(!is.null(modelH1)){
        
        lavArgs <- list(model = modelH1, data = cdata)
        if(isMultigroup) lavArgs <- append(lavArgs, list(group = 'gIdx'))
        
        lavresH1 <- do.call(lavaan::lavaan, append(lavArgs, lavOptionsH1))
        
        if(lavresH1@optim[["converged"]]){
          
          lavfitH1 <- lavaan::fitMeasures(lavresH1)
          fitH1 <- lavfitH1[c('fmin','chisq', 'df', 'pvalue')]
          if(lavresH1@Options[['test']] != 'standard'){
            fitH1 <- lavfitH1[c('fmin','chisq.scaled', 'df.scaled', 'pvalue.scaled')]
          }
          fitH1['fmin'] <- 2*fitH1['fmin'] # lav reports .5*fmin
          
          cfminGroupsH1 <- NA
          if(isMultigroup) cfminGroupsH1 <- lavresH1@Fit@test[[lavresH1@Options[['test']]]][['stat.group']] / (unlist(lavresH1@Data@nobs) - 1)
          
          fitH1 <- c(fitH1, cfminGroupsH1)
          
          if(lavresH1@Options[['test']] != 'standard'){
            mcomp <- lavaan::anova(lavresH0, lavresH1, method = 'satorra.bentler.2010') 
          }else{
            mcomp <- lavaan::anova(lavresH0, lavresH1) 
          }
          
          cfminDiff <- fitH0['fmin'] - fitH1['fmin']
          cfminGroupsDiff <- NA
          if(isMultigroup) cfminGroupsDiff <- cfminGroupsH0 - cfminGroupsH1
          
          fitDiff <- c(cfminDiff, 
                       mcomp[['Chisq diff']][2], mcomp[['Df diff']][2], mcomp[['Pr(>Chisq)']][2], 
                       cfminGroupsDiff)
          
          
          # store param est
          cLambda <- lavresH1@Model@GLIST[which(names(lavresH1@Model@GLIST) %in% 'lambda')]
          cLambda <- unlist(lapply(cLambda, function(x) x[x != 0])) # only non-zero loadings  
          cPsi <- unlist(lavresH1@Model@GLIST[which(names(lavresH1@Model@GLIST) %in% 'psi')])
          cBeta <- unlist(lavresH1@Model@GLIST[which(names(lavresH1@Model@GLIST) %in% 'beta')])
          cPhi <- unlist(lavresH1@Model@GLIST[which(names(lavresH1@Model@GLIST) %in% 'phi')])
          
          cres <- list(
            list(fitH0, fitH1, fitDiff),
            list(cLambda, cPsi, cBeta, cPhi)
          )
          
        }
        
        # h0 model only  
      }else{
        
        cres <- list(
          list(fitH0, rep(NA, length(fitH0)), fitDiff),
          list()
        )
        
      }
    }
    
    
  }, warning = function(w) {
     print(paste('WARNING: ', w))
  }, error = function(e) {
    print(paste('ERROR: ', e))
  })
  
  cres
  
}

#' genData
#' 
#' Generates random data from population variance-covariance matrix and population means, either
#' from a multivariate normal distribution, or using one of various approaches to generate 
#' non-normal data.
#'
#' @param N sample size.
#' @param Sigma population covariance matrix.
#' @param mu population means.
#' @param nSets number of data sets to generate
#' @param gIdx if not `NULL`, add gIdx as numeric group index as additional variable to generated data
#' @param modelH0 a `lavaan` model string, only used to determine the number of factors when `type = 'RK'`
#' @param simOptions additional arguments specifying the data generation routine
#' @return Returns the generated data
#' @examples
#' \dontrun{
#' gen <- semPower.genSigma(Phi = .2, loadings = list(rep(.5, 3), rep(.7, 3)))
#' data <- genData(N = 500, Sigma = gen$Sigma) 
#' }
#' @importFrom stats rbinom ecdf 
genData <- function(N = NULL, Sigma = NULL, mu = NULL, 
                    nSets = 1, gIdx = NULL, modelH0 = NULL, simOptions = NULL){
  
  type <- ifelse(!is.null(simOptions[['type']]), simOptions[['type']], 'normal')
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
         # non-normal using vm
         vm = {
           rdat <- genData.VM(N, Sigma, nSets, simOptions[['skewness']], simOptions[['kurtosis']])
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
#' Generates multivariate normal random data conforming to a population variance-covariance matrix.
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

#' genData.VM
#' 
#' Generates random data conforming to a population variance-covariance matrix using 
#' the third-order polynomial method  (Vale & Maurelli, 1983) 
#' specifying third and fourth moments of the marginals.
#' 
#' This function is a slightly adapted copy of `lavaan`'s ValeMaurelli1983 implementation
#' that avoids computing the intermediate correlation for each data sets
#' and uses Sigma as input.
#' 
#' For details, see 
#' Vale, C. & Maurelli, V. (1983). Simulating multivariate nonnormal distributions. *Psychometrika, 48*, 465-471.
#'
#' @param N sample size.
#' @param Sigma population covariance matrix.
#' @param nSets number of data sets to generate
#' @param skewness vector specifying skewness for each variable 
#' @param kurtosis vector specifying excess kurtosis for each variable
#' @return Returns the generated data
#' @importFrom stats nlminb
genData.VM <- function(N = NULL, Sigma = NULL, nSets = 1,  
                       skewness = NULL, kurtosis = NULL){
  
  if(is.null(skewness) || is.null(kurtosis)) stop('skewness and kurtosis must not be NULL.')
  if(length(skewness) != ncol(Sigma)) stop('skewness must match ncol(Sigma), i.e., must be specified for each variable.')
  if(length(kurtosis) != ncol(Sigma)) stop('kurtosis must match ncol(Sigma), i.e., must be specified for each variable.')
  if(any(kurtosis < (skewness^2 - 2))) stop('For VM, each marginal kurtosis must be larger than its marginal skewness^2 - 2')
  
  tryCatch({
    
    # the following is essentially a copy of lavaan's implementation
    
    fleishman1978_abcd <- function(skewness, kurtosis) {
      system.function <- function(x, skewness, kurtosis) {
        b.=x[1L]; c.=x[2L]; d.=x[3L]
        eq1 <- b.*b. + 6*b.*d. + 2*c.*c. + 15*d.*d. - 1
        eq2 <- 2*c.*(b.*b. + 24*b.*d. + 105*d.*d. + 2) - skewness
        eq3 <- 24*(b.*d. + c.*c.*(1 + b.*b. + 28*b.*d.) +
                     d.*d.*(12 + 48*b.*d. + 141*c.*c. + 225*d.*d.)) - kurtosis
        eq <- c(eq1,eq2,eq3)
        sum(eq*eq) ## SS
      }
      
      out <- nlminb(start=c(1,0,0), objective=system.function,
                    scale=10,
                    control=list(trace=0),
                    skewness=skewness, kurtosis=kurtosis)
      if(out$convergence != 0 || out$objective > 1e-5) {
        warning("lavaan WARNING: ValeMaurelli1983 method did not convergence, or it did not find the roots")
      }
      b. <- out$par[1L]; c. <- out$par[2L]; d. <- out$par[3L]; a. <- -c.
      c(a.,b.,c.,d.)
    }
    
    getICOV <- function(b1, c1, d1, b2, c2, d2, R) {
      objectiveFunction <- function(x, b1, c1, d1, b2, c2, d2, R) {
        rho=x[1L]
        eq <- rho*(b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2) +
          rho*rho*(2*c1*c2) + rho*rho*rho*(6*d1*d2) - R
        eq*eq
      }
      
      out <- nlminb(start=R, objective=objectiveFunction,
                    scale=10, control=list(trace=0),
                    b1=b1, c1=c1, d1=d1, b2=b2, c2=c2, d2=d2, R=R)
      if(out$convergence != 0 || out$objective > 1e-5) warning("no convergence")
      rho <- out$par[1L]
      rho
    }
    
    # number of variables
    nvar <- ncol(Sigma)

    # create Fleishman table
    FTable <- matrix(0, nvar, 4L)
    for(i in 1:nvar) {
      FTable[i,] <- fleishman1978_abcd(skewness=skewness[i], kurtosis=kurtosis[i])
    }
    
    # compute intermediate correlations between all pairs
    ICOR <- diag(nvar)
    for(j in 1:(nvar-1L)) {
      for(i in (j+1):nvar) {
        if(Sigma[i,j] == 0) next
        ICOR[i,j] <- ICOR[j,i] <-
          getICOV(FTable[i,2], FTable[i,3], FTable[i,4],
                  FTable[j,2], FTable[j,3], FTable[j,4], R=Sigma[i,j])
      }
    }
    
    # generate all requested data sets rather than just a single one
    # generate Z 
    lX <- lZ <- genData.normal(N, ICOR, nSets)
    for(d in seq(lX)){
      Z <- lZ[[d]]
      # transform Z using Fleishman constants
      for(i in 1:nvar) {
        lX[[d]][,i] <- FTable[i,1L] + FTable[i,2L]*Z[,i] + FTable[i,3L]*Z[,i]^2 + FTable[i,4L]*Z[,i]^3
      }
      colnames(lX[[d]]) <- paste0('x', 1:ncol(Sigma))
    }
    
    lX
    
  }, warning = function(e) {
    stop('Data generation via VM and the supplied arguments did not succeed. Either change the values for skewness and kurtosis, or try a different data generating routine such as IG.')
  }, error = function(e) {
    stop('Data generation via VM and the supplied arguments did not succeed. Either change the values for skewness and kurtosis, or try a different data generating routine such as IG.')
  })

}

#' genData.IG
#' 
#' Generates random data conforming to a population variance-covariance matrix using 
#' the independent generator approach (IG, Foldnes & Olsson, 2016) approach 
#' specifying third and fourth moments of the marginals.
#' 
#' This function is a wrapper for the respective function of the `covsim` package. 
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
#' Generates random data conforming to a population variance-covariance matrix using 
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
#' Generates random data conforming to a population variance-covariance matrix using 
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
      rdat <- matrix(0, nrow = N, ncol = k) 
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
    scale(rdat)  # scale because lav tends to complain about diverging variances  
  })
  
}


