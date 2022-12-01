##
## unittests
##
### TODO
### multigroup
### cfaPower
### regressionPower
### mediationPower
### clpmPower
###
### + plug all test into testthat

helper_lav <- function(model, sigma, sample.nobs = 1000, ...){
  lavres <- lavaan::sem(model, sample.cov = sigma, sample.nobs = sample.nobs, sample.cov.rescale = FALSE, ...)
  list(res = lavres, 
       summary = lavaan::summary(lavres, stand = T),
       par = lavaan::parameterestimates(lavres, stand = T),
       fit = lavaan::fitmeasures(lavres))
}


test_powerConsistency <- function(){
  ap <- semPower(type='a-priori', effect = .05, 'RMSEA', 
                 alpha = .05, beta = .05, df = 100)
  ph <- semPower(type='post-hoc', ap$fmin, 'f0', 
                 alpha = .05, N = ap$requiredN, df = 100)
  cp <- semPower(type='compromise', ap$fmin, 'F0', abratio = 1, 
                 N = ap$requiredN, df = 100)
  
  valid <- round(ap$impliedPower - ph$power, 3) == 0 &&
    round(ap$chiCrit - cp$chiCrit) == 0
  
  if(valid){
    print('powerConsistency: OK')
  }else{
    warning('Invalid')
  }
}

test_effectSizeConsistency <- function(){
  ph1 <- semPower.postHoc(.05, 'RMSEA', alpha = .05, N = 200, df = 200, p = 22)
  ph2 <- semPower.postHoc(ph1$fmin, 'F0', alpha = .05, N = 200, df = 200, p = 22)
  ph3 <- semPower.postHoc(ph1$mc, 'mc', alpha = .05, N = 200, df = 200, p = 22)
  ph4 <- semPower.postHoc(ph1$gfi, 'gfi', alpha = .05, N = 200, df = 200, p = 22)
  ph5 <- semPower.postHoc(ph1$agfi, 'agfi', alpha = .05, N = 200, df = 200, p = 22)
  
  valid <- round(ph1$power - ph2$power, 3) == 0 &&
    round(ph2$power - ph3$power, 3) == 0 &&
    round(ph3$power - ph4$power, 3) == 0 &&
    round(ph4$power - ph5$power, 3) == 0 
  
  if(valid){
    print('effectSizeConsistency: OK')
  }else{
    warning('Invalid')
  }
}

test_powerEffectDifference <- function(){
  RMSEA1 <- 0.06
  RMSEA2 <- 0.04
  df1 <- 27
  df2 <- 24
  deltaF <- df1*RMSEA1^2 - df2*RMSEA2^2
  
  ph <- semPower.postHoc(effect = deltaF, effect.measure =
                           "F0", alpha = .05, N = 200, df = df1 - df2)
  ph2 <- semPower.postHoc(effect = c(RMSEA1, RMSEA2), effect.measure = "RMSEA",
                          alpha = .05, N = 200, df = c(df1, df2))
  
  
  valid <- round(ph$power - ph2$power, 3) == 0 
  
  if(valid){
    print('test_powerEffectDifference: OK')
  }else{
    warning('Invalid')
  }
}

test_df <- function(){
  lavModel1 <- '
  f1 =~ x1 + x2 + x3 + x4
  f2 =~ x5 + x6 + x7 + x8
  f3 =~ y1 + y2 + y3
  f3 ~ f1 + f2
  '
  lavModel2 <- '
  f1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + y1 + y2 + y3
  '
  lavModel3 <- '
  f1 =~ x1 + x2 + x3 + x4
  f2 =~ x5 + x6 + x7 + x8
  '
  
  valid <-  semPower.getDf(lavModel1) == 41 &&
    semPower.getDf(lavModel2) == 44 &&
    semPower.getDf(lavModel3) == 19 &&
    semPower.getDf(lavModel3, nGroups = 3) == 57 &&
    semPower.getDf(lavModel3, nGroups = 3, group.equal = c('loadings')) == 69 &&
    # TODO shall we estimate latent means? lav does not by default, see also https://lavaan.ugent.be/tutorial/groups.html 
    semPower.getDf(lavModel3, nGroups = 3, group.equal = c('loadings', 'intercepts')) == 81
  
  if(valid){
    print('test_df: OK')
  }else{
    warning('Invalid df')
  }
}

test_generateSigma <- function(){
  
  # single loading
  generated <- semPower.genSigma(Phi = .2, nIndicator = c(5, 6), loadM = .5)
  lavres <- helper_lav(generated$modelTrue, generated$Sigma)
  par <- lavres$par
  
  valid1 <- round(lavres$fit['fmin'], 4) == 0 &&
    sum(par$lhs == 'f1' & par$op == '=~') == 5 &&
    sum(par$lhs == 'f2' & par$op == '=~') == 6 &&
    round(sum((par[par$op == '=~', 'std.all'] - .5)^2), 4) == 0 &&
    round(par[par$lhs == 'f1' & par$rhs == 'f2', 'std.all'], 4) == .2    
  
  # different loadings
  generated2 <- semPower.genSigma(Phi = .5, nIndicator = c(4, 3), loadM = c(.7, .8))
  lavres2 <- helper_lav(generated2$modelTrue, generated2$Sigma)
  par <- lavres2$par
  
  valid2 <- valid1 && 
    round(lavres2$fit['fmin'], 4) == 0 &&
    sum(par$lhs == 'f1' & par$op == '=~') == 4 &&
    sum(par$lhs == 'f2' & par$op == '=~') == 3 &&
    round(sum(par[par$lhs == 'f1' & par$op == '=~', 'std.all'] - .7)^2, 4) == 0 &&
    round(sum(par[par$lhs == 'f2' & par$op == '=~', 'std.all'] - .8)^2, 4) == 0 &&
    round(par[par$lhs == 'f1' & par$rhs == 'f2', 'std.all'], 4) == .5    
  
  # phi + different loadings for each factor
  Phi <- matrix(c(
    c(1.0, 0.2, 0.5),
    c(0.2, 1.0, 0.3),
    c(0.5, 0.3, 1.0)
  ), byrow = T, ncol=3)
  loadings <- list(
    c(0.4, 0.5, 0.8),
    c(0.7, 0.6, 0.5, 0.6),
    c(0.8, 0.8, 0.5)
  )
  generated3 <- semPower.genSigma(Phi = Phi, loadings = loadings)
  lavres3 <- helper_lav(generated3$modelTrue, generated3$Sigma)
  par <- lavres3$par
  
  valid3 <- valid2 &&
    round(lavres3$fit['fmin'], 4) == 0 &&
    sum(par$lhs == 'f1' & par$op == '=~') == length(loadings[[1]]) &&
    sum(par$lhs == 'f2' & par$op == '=~') == length(loadings[[2]]) &&
    sum(par$lhs == 'f3' & par$op == '=~') == length(loadings[[3]]) &&
    round(sum((par[par$lhs == 'f1' & par$op == '=~', 'std.all'] - loadings[[1]])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$op == '=~', 'std.all'] - loadings[[2]])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '=~', 'std.all'] - loadings[[3]])^2), 4) == 0 &&
    round(par[par$lhs == 'f1' & par$rhs == 'f2', 'std.all'], 4) == .2 && 
    round(par[par$lhs == 'f1' & par$rhs == 'f3', 'std.all'], 4) == .5 && 
    round(par[par$lhs == 'f2' & par$rhs == 'f3', 'std.all'], 4) == .3
  
  # use reference indicator
  generated3a <- semPower.genSigma(Phi = Phi, loadings = loadings, useReferenceIndicator = TRUE)
  lavres3a <- helper_lav(generated3a$modelTrue, generated3a$Sigma)
  valid3a <- valid3 &&
    # all std equal, except for residual variances
    round(sum((lavres3$par[lavres3$par$lhs != lavres3$par$rhs, 'std.all'] - lavres3a$par[lavres3a$par$lhs != lavres3a$par$rhs, 'std.all'])^2), 4) == 0
  
  # implied covariance matrix
  gSigma <- generated3$Lambda %*% generated3$Phi %*% t(generated3$Lambda)
  diag(gSigma) <- 1
  
  valid4 <- valid3a &&
    round(sum((gSigma - generated3$Sigma)^2), 4) == 0    
  
  # intercepts 
  tau <- (1:10)/10
  generated4 <- semPower.genSigma(Phi = .2, nIndicator = c(5, 5), loadM = .5, 
                                  tau = tau, alpha = c(0, 0))
  lavres4 <- helper_lav(generated4$modelTrue, generated4$Sigma, sample.mean = generated4$mu)
  par <- lavres4$par
  
  valid5 <- valid4 &&
    round(sum((par[par$lhs %in% paste0('x', 1:10) & par$op == '~1', 'est'] - tau)^2), 4) == 0
  
  # latent means
  generated5a <- semPower.genSigma(nIndicator = c(3, 3), loadM = .5, 
                                   tau = rep(0, 6), Alpha = c(0, 0))
  generated5b <- semPower.genSigma(nIndicator = c(3, 3), loadM = .5, 
                                   tau = rep(0, 6), Alpha = c(0, 1))
  
  m <- paste(generated5a$modelTrue, 'f1 ~ 0*1', 'f2 ~ 1', sep='\n') 
  suppressWarnings(lavres <- lavaan::sem(m,
                                         sample.cov = list(generated5a$Sigma, generated5b$Sigma),
                                         sample.mean = list(generated5a$mu, generated5b$mu),
                                         sample.nobs = list(500, 500), sample.cov.rescale = FALSE,
                                         group.equal = c('loadings', 'intercepts')))
  par <- lavaan::parameterestimates(lavres)
  
  valid6 <- valid5 &&
    round((par[par$lhs == 'f2' & par$op == '~1' & par$group == 2, 'est'] -
             par[par$lhs == 'f2' & par$op == '~1' & par$group == 1, 'est']) - 1, 3) == 0
  
  # phi + Lambda
  Phi <- matrix(c(
    c(1.0, 0.5, 0.1),
    c(0.5, 1.0, 0.2),
    c(0.1, 0.2, 1.0)
  ), byrow = T, ncol=3)
  Lambda <- matrix(c(
    c(0.4, 0.0, 0.0),
    c(0.7, 0.0, 0.0),
    c(0.8, 0.0, 0.0),
    c(0.0, 0.6, 0.0),
    c(0.0, 0.7, 0.0),
    c(0.0, 0.4, 0.0),
    c(0.0, 0.0, 0.8),
    c(0.0, 0.0, 0.7),
    c(0.0, 0.0, 0.8)
  ), byrow = T, ncol = 3)
  
  generated6 <- semPower.genSigma(Phi = Phi, Lambda = Lambda)
  lavres6 <- helper_lav(generated6$modelTrue, generated6$Sigma)
  par <- lavres6$par
  
  valid7 <- valid6 &&
    round(lavres6$fit['fmin'], 4) == 0 && 
    round(sum((par[par$lhs == 'f1' & par$op == '=~', 'std.all'] - Lambda[, 1][Lambda[, 1] != 0])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$op == '=~', 'std.all'] - Lambda[, 2][Lambda[, 2] != 0])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '=~', 'std.all'] - Lambda[, 3][Lambda[, 3] != 0])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f1' & par$rhs == 'f2', 'std.all'] - Phi[1, 2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f1' & par$rhs == 'f3', 'std.all'] - Phi[1, 3])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$rhs == 'f3', 'std.all'] - Phi[2, 3])^2), 4) == 0
  
  ## test lambda with cross loadings
  Lambda <- matrix(c(
    c(0.4, 0.0, 0.0),
    c(0.7, 0.0, 0.0),
    c(0.8, 0.0, 0.0),
    c(0.1, 0.6, 0.1),
    c(0.0, 0.7, 0.0),
    c(0.0, 0.4, 0.1),
    c(0.0, 0.0, 0.8),
    c(0.0, 0.0, 0.7),
    c(0.0, 0.0, 0.8)
  ), byrow = T, ncol = 3)
  generated7 <- semPower.genSigma(Phi = Phi, Lambda = Lambda)
  lavres7 <- helper_lav(generated7$modelTrue, generated7$Sigma)
  par <- lavres7$par
  
  valid8 <- valid7 &&
    round(lavres7$fit['fmin'], 4) == 0 && 
    round(sum((par[par$lhs == 'f1' & par$op == '=~', 'std.all'] - Lambda[, 1][Lambda[, 1] != 0])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$op == '=~', 'std.all'] - Lambda[, 2][Lambda[, 2] != 0])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '=~', 'std.all'] - Lambda[, 3][Lambda[, 3] != 0])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f1' & par$rhs == 'f2', 'std.all'] - Phi[1, 2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f1' & par$rhs == 'f3', 'std.all'] - Phi[1, 3])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$rhs == 'f3', 'std.all'] - Phi[2, 3])^2), 4) == 0
  
  if(valid8){
    print('test_generateSigma: OK')
  }else{
    warning('Invalid')
  }
}

test_genPhi <- function(){
  # std, no psi
  m <- '
  f4 ~ f1 + f2 + f3
  f3 ~ f1 + f2
  f2 ~ f1'
  B <- matrix(c(
    c(.00, .00, .00, .00),
    c(.30, .00, .00, .00),
    c(.20, .40, .00, .00),
    c(.10, .05, .50, .00)
  ), byrow = TRUE, ncol = 4)
  
  Phi <- getPhi.B(B)
  colnames(Phi) <- paste0('f', 1:4)
  par <- helper_lav(m, Phi)$par
  
  valid <- round(sum((par[par$lhs == 'f4' & par$op == '~', 'std.all'] - B[4, 1:3])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '~', 'std.all'] - B[3, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$op == '~', 'std.all'] - B[2, 1])^2), 4) == 0
  
  # unstd, no psi
  Phi <- getPhi.B(B, standardized = FALSE)
  colnames(Phi) <- paste0('f', 1:4)
  par <- helper_lav(m, Phi)$par
  
  valid2 <- valid && 
    round(sum((par[par$lhs == 'f4' & par$op == '~', 'est'] - B[4, 1:3])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '~', 'est'] - B[3, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$op == '~', 'est'] - B[2, 1])^2), 4) == 0
  
  # std, psi
  m2 <- '
  f3 + f4 ~ f1 + f2
  f1 ~~ f2
  f3 ~~ f4
  '
  B <- matrix(c(
    c(.00, .00, .00, .00),
    c(.30, .00, .00, .00),
    c(.70, .10, .00, .00),
    c(.20, .70, .00, .00)
  ), byrow = TRUE, ncol = 4)
  lPsi <- matrix(c(
    c(.00, .00, .00, .00),
    c(.00, .00, .00, .00),
    c(.00, .00, .00, .30),
    c(.00, .00, .30, .00)
  ), byrow = TRUE, ncol = 4)
  
  Phi <- getPhi.B(B, lPsi)
  colnames(Phi) <- paste0('f', 1:4)
  par <- helper_lav(m2, Phi)$par
  
  valid3 <- valid2 && 
    round(sum((par[par$lhs == 'f4' & par$op == '~', 'std.all'] - B[4, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '~', 'std.all'] - B[3, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$op == '~', 'std.all'] - B[2, 1])^2), 4) == 0 &&
    round((par[par$lhs == 'f3' & par$rhs == 'f4', 'std.all'] - lPsi[3, 4])^2, 4) == 0
  
  # unstd, psi
  Phi <- getPhi.B(B, lPsi, standardized = FALSE)
  colnames(Phi) <- paste0('f', 1:4)
  par <- helper_lav(m2, Phi)$par
  
  valid4 <- valid3 && 
    round(sum((par[par$lhs == 'f4' & par$op == '~', 'est'] - B[4, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '~', 'est'] - B[3, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$op == '~', 'est'] - B[2, 1])^2), 4) == 0 &&
    round((par[par$lhs == 'f3' & par$rhs == 'f4', 'est'] - lPsi[3, 4])^2, 4) == 0    
  
  # gen sigma and unstd phi + different loadings for each factor
  B <- matrix(c(
    c(.00, .00, .00, .00),
    c(.30, .00, .00, .00),
    c(.70, .10, .00, .00), 
    c(.20, .60, .00, .00)
  ), byrow = TRUE, ncol = 4)
  lPsi <- matrix(c(
    c(.00, .00, .00, .00),
    c(.00, .00, .00, .00),
    c(.00, .00, .00, .10),
    c(.00, .00, .10, .00)
  ), byrow = TRUE, ncol = 4)
  loadings <- list(
    c(0.4, 0.5, 0.8),
    c(0.6, 0.6, 0.5),
    c(0.4, 0.5, 0.8),
    c(0.6, 0.6, 0.5)
  )
  Phi <- getPhi.B(B, lPsi, standardized = FALSE)
  generated6 <- semPower.genSigma(Phi = Phi, loadings = loadings, useReferenceIndicator = TRUE)
  m <- paste(generated6$modelTrue, 
             'f3 + f4 ~ f1 + f2
             f1 ~~ f2
             f3 ~~ f4'
             , sep = '\n')
  lavres6 <- helper_lav(m, generated6$Sigma)
  par <- lavres6$par
  
  valid5 <- valid4 &&
    round(lavres6$fit['fmin'], 4) == 0 &&
    sum(par$lhs == 'f1' & par$op == '=~') == length(loadings[[1]]) &&
    sum(par$lhs == 'f2' & par$op == '=~') == length(loadings[[2]]) &&
    sum(par$lhs == 'f3' & par$op == '=~') == length(loadings[[3]]) &&
    round(sum((par[par$lhs == 'f1' & par$op == '=~', 'est'] - loadings[[1]])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$op == '=~', 'est'] - loadings[[2]])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '=~', 'est'] - loadings[[3]])^2), 4) == 0 &&
    round(par[par$lhs == 'f3' & par$rhs == 'f1', 'est'] - B[3, 1], 4) == 0 && 
    round(par[par$lhs == 'f3' & par$rhs == 'f2', 'est'] - B[3, 2], 4) == 0 && 
    round(par[par$lhs == 'f4' & par$rhs == 'f1', 'est'] - B[4, 1], 4) == 0 && 
    round(par[par$lhs == 'f4' & par$rhs == 'f2', 'est'] - B[4, 2], 4) == 0 && 
    round(par[par$lhs == 'f1' & par$rhs == 'f2', 'est'] - B[2, 1], 4) == 0 && 
    round(par[par$lhs == 'f3' & par$rhs == 'f4', 'est'] - lPsi[3, 4], 4) == 0 
  
  
  # gen sigma and std phi + observed factors
  B <- matrix(c(
    c(.00, .00, .00),
    c(.30, .00, .00),
    c(.20, .40, .00)
  ), byrow = TRUE, ncol = 3)
  Phi <- getPhi.B(B)
  generated6 <- semPower.genSigma(Phi = Phi, nIndicator = rep(1, 3), loadM = 1)
  round(sum((generated6$Sigma - Phi)^2), 4) == 0
  
  # same with Lambda as diag matrix
  generated6b <- semPower.genSigma(Phi = Phi, Lambda = diag(3))
  
  
  valid6 <- valid5 &&
    round(sum((generated6$Sigma - Phi)^2), 4) == 0 &&
    round(sum((generated6b$Sigma - Phi)^2), 4) == 0
  
  if(valid5){
    print('test_genPhi: OK')
  }else{
    warning('Invalid')
  }
}

test_all <- function(){
  test_powerConsistency()  
  test_effectSizeConsistency()
  test_powerEffectDifference()
  test_df()
  test_generateSigma()
  test_genPhi()
}

test_all()



#####


