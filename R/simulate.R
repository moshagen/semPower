#' simulate
#'
#' estimate empirical power using a simulation approach
#' 
#' @param modelH0 lavaan model string defining the (incorrect) analysis model.
#' @param modelH1 lavaan model string defining the comparison model. If omitted, the saturated model is the comparison model.
#' @param Sigma population covariance matrix.
#' @param mu population means.
#' @param N sample size
#' @param alpha alpha error probability
#' @param nReplications number of random samples drawn.
#' @param minConvergenceRate the minimum convergence rate required
#' @param lavOptions a list of additional options passed to lavaan, e.g., list(estimator = 'mlm') to request robust ML estimation
#' @param returnFmin whether return the median unbiased Fmin over replications (i.e, fmin_0 = fmin_hat - df/N) 
#' @return empirical power: sum(p < alpha) / nReplications
#' @examples
#' \dontrun{
#' }
#' @importFrom utils txtProgressBar 
simulate <- function(modelH0 = NULL, modelH1 = NULL,
                     Sigma = NULL, mu = NULL, 
                     N = NULL,
                     alpha = NULL,
                     nReplications = 250, minConvergenceRate = .5,
                     lavOptions = NULL,
                     returnFmin = TRUE){

  checkBounded(nReplications, bound = c(10, 100000), inclusive = TRUE)
  checkBounded(minConvergenceRate, bound = c(0.01, 1), inclusive = TRUE)
  # remaining checks are done in corresponding power fnc
  
  # we need to call lavaan() directly with defaults as defined in sem()
  lavOptions <- getLavOptions(lavOptions, isCovarianceMatrix = FALSE, nGroups = length(Sigma))
  
  efmin <- efminGroups <- list()
  ePower <- 0
  r <- rr <- 1
  progress <- txtProgressBar(min = 0, max = nReplications, initial = 0, style = 3)
  while(r <= nReplications && rr <= nReplications / minConvergenceRate){
    setTxtProgressBar(progress, r)
    tryCatch({
      if(!is.list(Sigma)){
        # single group case
        cdata <- genData(N, Sigma, mu)
        lavresH0 <- do.call(lavaan::lavaan, 
                            append(list(model = modelH0, data = cdata), lavOptions))
        cfminGroups <- NULL
      }else{
        # multigroup group case
        gdata <- lapply(1:length(Sigma), function(x) genData(N[[x]], Sigma[[x]], mu[[x]], gIdx = x))
        cdata <- unlist(do.call(rbind, gdata))
        lavresH0 <- do.call(lavaan::lavaan, 
                            append(list(model = modelH0, data = cdata), append(list(group = 'gIdx'), lavOptions)))
        # store fmin by group
        if(lavresH0@Options[['estimator']] %in% c("ML", "MLF", "WLS", "DWLS", "ULS")) testType <- 'standard' 
        if(lavresH0@Options[['estimator']] %in% c("MLM", "WLSM", "ULSM")) testType <- 'satorra.bentler' 
        if(lavresH0@Options[['estimator']] %in% c("MLR")) testType <- 'yuan.bentler.mplus' 
        if(lavresH0@Options[['estimator']] %in% c("MLMV", "WLSMV")) testType <- 'scaled.shifted' 
        if(lavresH0@Options[['estimator']] %in% c("MLMVS")) testType <- 'mean.var.adjusted' 
        cfminGroups <- lavresH0@Fit@test[[testType]][['stat.group']] / (unlist(N) - 1)
      }
      cfmin <- 2 * lavaan::fitMeasures(lavresH0, 'fmin') # lav reports .5*fmin
      p <- lavaan::fitMeasures(lavresH0, 'pvalue')
      df <- lavaan::fitMeasures(lavresH0, 'df')
      if(lavresH0@Options[['estimator']] %in% c("MLM", "MLMV", "MLMVS", "MLF", "MLR", "WLS", "DWLS", "WLSM", "WLSMV", "ULSM", "ULSMV")){
        p <- lavaan::fitMeasures(lavresH0, 'pvalue.scaled')
        df <- lavaan::fitMeasures(lavresH0, 'df.scaled')
      }

      # handle restricted comparison model 
      # (modelH1 must always get fit because sampling error does not allow just using modelH0 estm with different df)
      if(!is.null(modelH1)){
        if(!is.list(Sigma)){
          # single group case
          lavresH1 <- do.call(lavaan::lavaan, 
                              append(list(model = modelH1, data = cdata), lavOptions))
        }else{
          # multigroup group case
          lavresH1 <- do.call(lavaan::lavaan, 
                              append(list(model = modelH1, data = cdata), append(list(group = 'gIdx'), lavOptions)))
          cfminGroups <- cfminGroups - lavresH1@Fit@test[[testType]][['stat.group']] / (unlist(N) - 1)
        }
        mcomp <- lavaan::anova(lavresH0, lavresH1) 
        p <- mcomp$`Pr(>Chisq)`[2]
        df <- mcomp$`Df diff`[2]
        cfmin <- 2 * (lavaan::fitMeasures(lavresH0, 'fmin') - lavaan::fitMeasures(lavresH1, 'fmin'))
      }
      efmin <- append(efmin, cfmin)
      efminGroups <- append(efminGroups, list(cfminGroups)) 
                            
      if(p < alpha)
        ePower <- ePower + 1
      
      r <- r + 1
      
      # TODO 
      # handle relevant errors above, as this would also allow
      # to include non-properly converged models, as these are probably also relevant for
      # empirical power
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
  if(round((r - 1) / (rr - 1), 2) < minConvergenceRate){ 
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
    list(
      ePower = ePower,
      meanFmin = ubFmean,
      meanFminGroups = ubFmeanGroups,
      df = df,
      nrep = (r - 1)
    )
  }else{
    ePower
  }
}

#' genData
#' 
#' generate random data from population variance-covariance matrix and population means.
#' currently normal case only.
#'
#' @param N sample size.
#' @param Sigma population covariance matrix.
#' @param mu population means.
#' @param gIdx if not NULL, add gIdx as numeric group index as additional variable to generated data
#' @return data
#' @examples
#' \dontrun{
#' }
#' @importFrom stats rnorm 
genData <- function(N = NULL, Sigma = NULL, mu = NULL, gIdx = NULL){
  randomData <- matrix(rnorm(N * ncol(Sigma)), N, ncol(Sigma)) 
  rdat <- t(t(chol(Sigma)) %*% t(randomData))
  if(!is.null(mu)) rdat <- rdat + matrix(t(rep(mu, N)), ncol = ncol(Sigma), byrow = TRUE)
  if(!is.null(gIdx)) rdat <- cbind(rdat, matrix(rep(gIdx, N), ncol = 1, dimnames = list(NULL,c('gIdx')))) 
  rdat
}

