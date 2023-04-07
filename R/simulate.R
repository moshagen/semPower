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
#' @param simOptions a list of additional options specifying simulation details, see detais. 
#' @param lavOptions a list of additional options passed to `lavaan`, e. g., `list(estimator = 'mlm')` to request robust ML estimation
#' @param lavOptionsH1 lavoptions when fitting `modelH1`. If `NULL`, the same as `lavOptions`.
#' @param returnFmin whether to return the mean unbiased Fmin over replications (i. e., `fmin_0 = fmin_hat - df/N`) 
#' @return Returns empirical power: `sum(p < alpha) / nReplications` or a list (if `returnFmin = TRUE`) with the following components:
#' \item{`ePower`}{the empirical power.}
#' \item{`meanFmin`}{the estimated mean unbiased Fmin over replications (i. e., `fmin_0 = fmin_hat - df/N`).}
#' \item{`meanFminGroups`}{the estimated mean unbiased Fmin by groups given as a vector, assuming the df spread equally over groups. Therefore, `meanFmin != sum(meanFminGroups)`}
#' \item{`df`}{the model df.}
#' \item{`nrep`}{the successful number of replications.}
#' \item{`bChiSq`}{median chi-square bias of the H1 model}
#' \item{`bLambda`}{average median bias in lambda in the H1 model}
#' \item{`bPhi`}{average median bias in phi in the H1 model}
#' \item{`bPsi`}{average median bias in psi in the H1 model}
#' \item{`bBeta`}{average median bias in beta in the H1 model}
#' @details 
#' 
#' `simOptions` is a list that may have the following components:
#' * `nReplications`: The targeted number of valid simulation runs, defaults to 250.
#' * `minConvergenceRate`:  The minimum convergence rate required, defaults to .5. The maximum actual simulation runs are increased by a factor of 1/minConvergenceRate.
#'  
#' @examples
#' \dontrun{
#' # get Sigma and modelH0
#' powerCFA <- semPower.powerCFA(type = 'a-priori',
#'                               comparison = 'saturated',
#'                               Phi = .2, loadings = list(rep(.5, 3), rep(.7, 3)), 
#'                               alpha = .05, beta = .05)
#' # perform simulation       
#' simulate(modelH0 = powerCFA$modelH0, 
#'          powerCFA$Sigma,
#'          N = powerCFA$power$requiredN,
#'          alpha = powerCFA$power$alpha)
#' }
#' @importFrom utils txtProgressBar setTxtProgressBar
simulate <- function(modelH0 = NULL, modelH1 = NULL,
                     Sigma = NULL, mu = NULL, 
                     N = NULL,
                     alpha = NULL,
                     simOptions = list(
                       nReplications = 250, 
                       minConvergenceRate = .5
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
  
  # warn in case of long computation times
  lavEstimators <- c("MLM", "MLMV", "MLMVS", "MLF", "MLR", "WLS", "DWLS", "WLSM", "WLSMV", "ULSM", "ULSMV")
  costlyEstm <- (!is.null(lavOptions[['estimator']]) && toupper(lavOptions[['estimator']]) %in% lavEstimators)
  projectedLong <- (costlyEstm && ncol(Sigma) > 50 && nReplications > 100) || (ncol(Sigma) > 100 && nReplications > 500) || nReplications > 10000
  if(projectedLong){
    mResp <- menu(c("Yes", "No"), title = "Simulated power with the specified model, the number of replications, and the type of estimator will presumably take a long time.\nDo you really want to go on?")
    if(mResp != 1){
      stop("Simulated power aborted.")
    }else{
      cat("You have been warned.\n")
    }
  }
  
  
  # we need to call lavaan() directly with defaults as defined in sem()
  lavOptions <- getLavOptions(lavOptions, isCovarianceMatrix = FALSE, nGroups = length(Sigma))
  lavOptions <- append(lavOptions,
                       list(check.gradient = FALSE, check.post = FALSE, check.vcov = FALSE))
  if(!is.null(modelH1)){
    lavOptionsH1 <- getLavOptions(lavOptionsH1, isCovarianceMatrix = FALSE, nGroups = length(Sigma))
    lavOptionsH1 <- append(lavOptionsH1,
                           list(check.gradient = FALSE, check.post = FALSE, check.vcov = FALSE))
  }
  
  # generate data 
  if(!is.list(Sigma)){
    simData <- genData(N, Sigma, mu, gIdx = NULL, nSets = nReplications / minConvergenceRate)
  }else{
    simData <- lapply(seq_along(Sigma), function(x) genData(N[[x]], Sigma[[x]], mu[[x]], gIdx = x, nSets = nReplications / minConvergenceRate))
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
            rChiSq <- append(rChiSq, lavaan::fitMeasures(lavresH1, 'chisq'))
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
      lavresPop <- do.call(lavaan::lavaan, 
                           append(list(model = modelH1, sample.cov = Sigma, sample.mean = mu, sample.nobs = sample.nobs, 
                                       sample.cov.rescale = FALSE),
                                  lavOptionsH1))
      if(lavresPop@optim$fx > 1e-6) warning('H1 model is not properly specified.')

      # calc bias
      bChiSq <- (median(unlist(rChiSq)) - lavaan::fitMeasures(lavresPop, 'df')) / lavaan::fitMeasures(lavresPop, 'df')
      
      cLambda <- lavresPop@Model@GLIST[which(names(lavresPop@Model@GLIST) %in% 'lambda')]
      pLambda <- unlist(lapply(cLambda, function(x) x[x != 0]))
      bLambda <- mean( (apply(do.call(rbind, rLambda), 2, median) -  pLambda) / pLambda )
      
      bPhi <- bPsi <- bBeta <- NULL
      # for phi/psi/beta, we only consider population parameters that are larger than
      # a small constant to avoid absurd relative biases for parameter with true values close to zero 
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
#' Generates random data from population variance-covariance matrix and population means.
#' Currently only the multivariate normal case.
#'
#' @param N sample size.
#' @param Sigma population covariance matrix.
#' @param mu population means.
#' @param nSets number of data sets to generate
#' @param gIdx if not `NULL`, add gIdx as numeric group index as additional variable to generated data
#' @return Returns the generated data
#' @examples
#' \dontrun{
#' gen <- semPower.genSigma(Phi = .2, loadings = list(rep(.5, 3), rep(.7, 3)))
#' data <- genData(N = 500, Sigma = gen$Sigma) 
#' }
#' @importFrom stats rnorm 
genData <- function(N = NULL, Sigma = NULL, mu = NULL, nSets = 1, gIdx = NULL){
  lapply(seq(nSets), function(x){
    randomData <- matrix(rnorm(N * ncol(Sigma)), N, ncol(Sigma)) 
    rdat <- t(t(chol(Sigma)) %*% t(randomData))
    if(!is.null(mu)) rdat <- rdat + matrix(t(rep(mu, N)), ncol = ncol(Sigma), byrow = TRUE)
    if(!is.null(gIdx)) rdat <- cbind(rdat, matrix(rep(gIdx, N), ncol = 1, dimnames = list(NULL, c('gIdx')))) 
    rdat
  })
}

