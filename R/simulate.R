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
  lavOptionsDefaults <- list(int.ov.free = TRUE, int.lv.free = FALSE, auto.fix.first = TRUE,
                             auto.fix.single = TRUE, auto.var = TRUE, auto.cov.lv.x = TRUE,
                             auto.efa = TRUE, auto.th = TRUE, auto.delta = TRUE, auto.cov.y = TRUE) 

  # we also want N - 1 for ml based estm
  if(is.null(lavOptions$estimator) || (!is.null(lavOptions$estimator) && toupper(lavOptions$estimator) %in% c("ML", "MLM", "MLMV", "MLMVS", "MLF", "MLR")))
    lavOptionsDefaults <- append(list(likelihood = 'Wishart'), lavOptionsDefaults)
  
  lavOptions <- append(lavOptions, lavOptionsDefaults)
  
  ef <- list()
  ePower <- 0
  r <- rr <- 1
  progress <- txtProgressBar(min = 0, max = nReplications, initial = 0, style = 3)
  while(r <= nReplications && rr <= nReplications / minConvergenceRate){
    setTxtProgressBar(progress, r)
    tryCatch({
      
      cdata <- genData(N, Sigma, mu)
      lavresH0 <- do.call(lavaan::lavaan, 
                          append(list(model = modelH0, data = cdata), lavOptions))
      ef <- append(ef, 2 * lavaan::fitMeasures(lavresH0, 'fmin')) # lav reports .5*fmin
      p <- lavaan::fitMeasures(lavresH0, 'pvalue')
      df <- lavaan::fitMeasures(lavresH0, 'df')
      if(!is.null(lavOptions$estimator) && toupper(lavOptions$estimator) %in% c("MLM", "MLMV", "MLMVS", "MLF", "MLR", "WLS", "DWLS", "WLSM", "WLSMV", "ULSM", "ULSMV")){
        p <- lavaan::fitMeasures(lavresH0, 'pvalue.scaled')
        df <- lavaan::fitMeasures(lavresH0, 'df.scaled')
      }

      # handle restricted comparison model 
      # (modelH1 must always get fit because sampling error does not allow just using modelH0 estm with differen df)
      if(!is.null(modelH1)){
        lavresH1 <- do.call(lavaan::lavaan, 
                            append(list(model = modelH1, data = cdata), lavOptions))
        mcomp <- lavaan::anova(lavresH0, lavresH1) 
        p <- mcomp$`Pr(>Chisq)`[2]
        df <- mcomp$`Df diff`[2]
        ef <- 2 * (lavaan::fitMeasures(lavresH0, 'fmin') - lavaan::fitMeasures(lavresH1, 'fmin'))
      }
      
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
      # print(paste('ERROR: ',e))
    })
    rr <- rr + 1
  }
  close(progress)
  
  ePower <- ePower / nReplications
  if((rr - 1) > (nReplications / minConvergenceRate)){
    warning(paste("Actual convergence rate of", round((r - 1) / (rr - 1), 2), "is below minConvergenceRate of", minConvergenceRate, ". Results are based on", nReplications,"replications."))
  }
  
  if(returnFmin){
    # sample fmin is biased, we report need unbiased pop fmin
    ubFmedian <- median(unlist(ef) - df / N)
    list(
      ePower = ePower,
      medianF = ubFmedian,
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
#' @return data
#' @examples
#' \dontrun{
#' }
#' @importFrom stats rnorm 
genData <- function(N = NULL, Sigma = NULL, mu = NULL){
  randomData <- matrix(rnorm(N * ncol(Sigma)), N, ncol(Sigma)) 
  rdat <- t(t(chol(Sigma)) %*% t(randomData))
  if(!is.null(mu)) rdat <- rdat + matrix(t(rep(mu, N)), ncol = ncol(Sigma), byrow = TRUE)
  rdat
}

