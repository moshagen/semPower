##
## unittests
##

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
    print('test_powerConsistency: OK')
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
    print('test_effectSizeConsistency: OK')
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

test_generateSigmaB <- function(){
  # manifest, unstd, no psi
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

  generated <- semPower.genSigma(Beta = B, Lambda = diag(ncol(B)))
  
  colnames(generated$Sigma) <- rownames(generated$Sigma) <- paste0('f', 1:4)
  par <- helper_lav(m, generated$Sigma)$par
  
  valid <- round(sum((par[par$lhs == 'f4' & par$op == '~', 'est'] - B[4, 1:3])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '~', 'est'] - B[3, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$op == '~', 'est'] - B[2, 1])^2), 4) == 0

  # manifest, unstd,  psi
  m <- '
  f4 ~ f1 + f2
  f3 ~ f1 + f2
  f2 ~~ f1
  f4 ~~ f3
  '
  B <- matrix(c(
    c(.00, .00, .00, .00),
    c(.00, .00, .00, .00),
    c(.40, .20, .00, .00),
    c(.10, .50, .00, .00)
  ), byrow = TRUE, ncol = 4)
  Psi <- matrix(c(
    c(1, .30, .00, .00),
    c(.30, 1, .00, .00),
    c(.00, .00, 1, .30),
    c(.00, .00, .30, 1)
  ), byrow = TRUE, ncol = 4)
  
  generated <- semPower.genSigma(Beta = B, Psi = Psi, Lambda = diag(ncol(B)))
  colnames(generated$Sigma) <- rownames(generated$Sigma)  <- paste0('f', 1:4)
  par <- helper_lav(m, generated$Sigma)$par

  valid2 <- valid &&
    round(sum((par[par$lhs == 'f4' & par$op == '~', 'est'] - B[4, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '~', 'est'] - B[3, 1:2])^2), 4) == 0 &&
    round((par[par$lhs == 'f1' & par$rhs == 'f2', 'est'] - .3), 4) == 0 &&
    round((par[par$lhs == 'f4' & par$rhs == 'f3', 'est'] - .3), 4) == 0 
  
  # manifest, clpm type model with residual correlations at wave 2 + 3
  m <- '
  x3 + x4 ~ x1 + x2
  x5 + x6 ~ x3 + x4
  x1 ~~ x2
  x3 ~~ x4
  x5 ~~ x6
  '
  B <- matrix(c(
    c(.0, .0, .0, .0, 0, 0),  # X1
    c(.0, .0, .0, .0, 0, 0),  # Y1
    c(.7, .1, .0, .0, 0, 0),  # X2
    c(.2, .8, .0, .0, 0, 0),  # Y2
    c(.0, .0, .7, .1, 0, 0),  # X3
    c(.0, .0, .2, .8, 0, 0)   # Y3
  ), byrow = TRUE, ncol = 6)
  Psi <- diag(ncol(B))
  Psi[3,4] <- Psi[4,3] <- .2
  Psi[5,6] <- Psi[6,5] <- .3
  
  generated <- semPower.genSigma(Beta = B, Psi = Psi, Lambda = diag(ncol(B)), useReferenceIndicator = T)
  par <- helper_lav(m, generated$Sigma)$par
  # test modelPop
  lSigma <- lavaan::fitted(lavaan::sem(generated$modelPop))$cov
  # test modelTrue and modelPop
  lavres <- helper_lav(generated$modelTrue, generated$Sigma)
  par2 <- lavres$par

  valid3 <- valid2 &&
    round(sum((par[par$lhs == 'x3' & par$op == '~', 'est'] - B[3, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'x4' & par$op == '~', 'est'] - B[4, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'x5' & par$op == '~', 'est'] - B[5, 3:4])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'x6' & par$op == '~', 'est'] - B[6, 3:4])^2), 4) == 0 &&
    round(par[par$lhs == 'x3' & par$rhs == 'x4', 'est'] - Psi[3, 4], 4) == 0 &&
    round(par[par$lhs == 'x5' & par$rhs == 'x6', 'est'] - Psi[5, 6], 4) == 0 &&
    round(sum((lSigma - generated$Sigma)^2), 4) == 0 &&
    round(lavres$fit['fmin'], 4) == 0 &&
    round(sum((par2[par2$lhs %in% paste0('f', 1:6) & par2$rhs %in% paste0('f', 1:6) & par2$op == '~', 'est'] - par[par$lhs %in% paste0('x', 1:6) & par$op == '~', 'est'])^2), 4) == 0 && 
    round(sum((par2[par2$lhs %in% c('f3', 'f5') & par2$rhs %in% c('f4', 'f6') & par2$op == '~~', 'est'] - par[par$lhs %in% c('x3', 'x5') & par$rhs %in% c('x4', 'x6') & par$op == '~~', 'est'])^2), 4) == 0

  # latent, unstd,  psi
  B <- matrix(c(
    c(.00, .00, .00, .00),
    c(.00, .00, .00, .00),
    c(.40, .20, .00, .00),
    c(.10, .50, .00, .00)
  ), byrow = TRUE, ncol = 4)
  Psi <- matrix(c(
    c(1, .30, .00, .00),
    c(.30, 1, .00, .00),
    c(.00, .00, 1, .30),
    c(.00, .00, .30, 1)
  ), byrow = TRUE, ncol = 4)
  
  generated <- semPower.genSigma(Beta = B, Psi = Psi, 
                                 loadM = c(.5, .6, .5, .6), nIndicator = rep(3, 4))
  par <- helper_lav(generated$modelTrue, generated$Sigma)$par
  
  # same with referent ind
  generated <- semPower.genSigma(Beta = B, Psi = Psi, 
                                 loadM = c(.5, .6, .5, .6), nIndicator = rep(3, 4),
                                 useReferenceIndicator = TRUE)
  par2 <- helper_lav(generated$modelTrue, generated$Sigma)$par

  valid4 <- valid3 &&
    round(sum((par[par$lhs %in% c('f1', 'f3') & par$op == '=~', 'est'] - .5)^2), 4) == 0 &&
    round(sum((par[par$lhs %in% c('f2', 'f4') & par$op == '=~', 'est'] - .6)^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '~', 'est'] - B[3, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f4' & par$op == '~', 'est'] - B[4, 1:2])^2), 4) == 0 &&
    round(par[par$lhs == 'f1' & par$rhs == 'f2', 'est'] - Psi[1, 2], 4) == 0 &&
    round(par[par$lhs == 'f3' & par$rhs == 'f4', 'est'] - Psi[3, 4], 4) == 0 &&
    round(sum((par2[par2$lhs != par2$rhs, 'est']  - par[par$lhs != par$rhs, 'est'])^2), 4) == 0

  # mixed observed and latent, unstd,  psi
  generated <- semPower.genSigma(Beta = B, Psi = Psi, 
                                 loadM = c(.5, .6, .5, .6), 
                                 nIndicator = c(3, 1, 1, 4))
  par <- helper_lav(generated$modelTrue, generated$Sigma)$par
  
  valid5 <- valid4 &&
    round(sum((par[par$lhs %in% c('f1', 'f3') & par$op == '=~', 'est'] - .5)^2), 4) == 0 &&
    round(sum((par[par$lhs %in% c('f2', 'f4') & par$op == '=~', 'est'] - .6)^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '~', 'est'] - B[3, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f4' & par$op == '~', 'est'] - B[4, 1:2])^2), 4) == 0 &&
    round(par[par$lhs == 'f1' & par$rhs == 'f2', 'est'] - Psi[1, 2], 4) == 0 &&
    round(par[par$lhs == 'f3' & par$rhs == 'f4', 'est'] - Psi[3, 4], 4) == 0 

  # test genLambda
  loadings <- list(
    c(0.4, 0.5, 0.8),
    c(0.7, 0.6, 0.5, 0.6),
    c(0.8, 0.8, 0.5)
  )
  Lambda <- genLambda(loadings)
  generated <- semPower.genSigma(Lambda = Lambda, Phi = .3)
  generated2 <- semPower.genSigma(loadings = loadings, Phi = .3)
  
  valid6 <- valid5 && 
    round(sum((generated$Sigma - generated2$Sigma)^2), 4) == 0
  
  if(valid6){
    print('test_generateSigmaB: OK')
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
  par <- helper_lav(m, Phi)$par

  valid <- round(sum((par[par$lhs == 'f4' & par$op == '~', 'std.all'] - B[4, 1:3])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '~', 'std.all'] - B[3, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$op == '~', 'std.all'] - B[2, 1])^2), 4) == 0

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

  valid2 <- valid &&
    round(sum((par[par$lhs == 'f4' & par$op == '~', 'std.all'] - B[4, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '~', 'std.all'] - B[3, 1:2])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$op == '~', 'std.all'] - B[2, 1])^2), 4) == 0 &&
    round((par[par$lhs == 'f3' & par$rhs == 'f4', 'std.all'] - lPsi[3, 4])^2, 4) == 0

  # gen sigma and phi + different loadings for each factor
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
  Phi <- getPhi.B(B, lPsi)
  generated6 <- semPower.genSigma(Phi = Phi, loadings = loadings, useReferenceIndicator = TRUE)
  m <- paste(generated6$modelTrue,
             'f3 + f4 ~ f1 + f2
             f1 ~~ f2
             f3 ~~ f4'
             , sep = '\n')
  lavres6 <- helper_lav(m, generated6$Sigma)
  par <- lavres6$par

  valid3 <- valid2 &&
    round(lavres6$fit['fmin'], 4) == 0 &&
    sum(par$lhs == 'f1' & par$op == '=~') == length(loadings[[1]]) &&
    sum(par$lhs == 'f2' & par$op == '=~') == length(loadings[[2]]) &&
    sum(par$lhs == 'f3' & par$op == '=~') == length(loadings[[3]]) &&
    round(sum((par[par$lhs == 'f1' & par$op == '=~', 'std.all'] - loadings[[1]])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f2' & par$op == '=~', 'std.all'] - loadings[[2]])^2), 4) == 0 &&
    round(sum((par[par$lhs == 'f3' & par$op == '=~', 'std.all'] - loadings[[3]])^2), 4) == 0 &&
    round(par[par$lhs == 'f3' & par$rhs == 'f1', 'std.all'] - B[3, 1], 4) == 0 &&
    round(par[par$lhs == 'f3' & par$rhs == 'f2', 'std.all'] - B[3, 2], 4) == 0 &&
    round(par[par$lhs == 'f4' & par$rhs == 'f1', 'std.all'] - B[4, 1], 4) == 0 &&
    round(par[par$lhs == 'f4' & par$rhs == 'f2', 'std.all'] - B[4, 2], 4) == 0 &&
    round(par[par$lhs == 'f1' & par$rhs == 'f2', 'std.all'] - B[2, 1], 4) == 0 &&
    round(par[par$lhs == 'f3' & par$rhs == 'f4', 'std.all'] - lPsi[3, 4], 4) == 0

  # gen sigma and std phi + observed factors
  B <- matrix(c(
    c(.00, .00, .00),
    c(.30, .00, .00),
    c(.20, .40, .00)
  ), byrow = TRUE, ncol = 3)
  Phi <- getPhi.B(B)
  generated6 <- semPower.genSigma(Phi = Phi, nIndicator = rep(1, 3), loadM = 1)

  # same with Lambda as diag matrix
  generated6b <- semPower.genSigma(Phi = Phi, Lambda = diag(3))

  valid4 <- valid3 &&
    round(sum((generated6$Sigma - Phi)^2), 4) == 0 &&
    round(sum((generated6b$Sigma - Phi)^2), 4) == 0

  if(valid4){
    print('test_genPhi: OK')
  }else{
    warning('Invalid')
  }
}

test_multigroup <- function(){
  
  # metric invariance model
  generated <- semPower.genSigma(loadings = list(c(.5, .6, .7)))
  generated2 <- semPower.genSigma(loadings = list(c(.5, .5, .7)))
  lavres <- helper_lav(generated$modelTrue, 
                       list(generated$Sigma, generated2$Sigma),
                       sample.nobs = list(500, 500),
                       group.equal = c('loadings'))
  sigmaHat1 <- lavaan::fitted(lavres$res)$`Group 1`$cov
  sigmaHat2 <- lavaan::fitted(lavres$res)$`Group 2`$cov
  
  ph <- semPower.postHoc(SigmaHat = list(sigmaHat1, sigmaHat2),
                         Sigma = list(generated$Sigma, generated2$Sigma),
                         alpha = .05, N = list(500, 500), df = 3)
  
  # f contrib by group
  getFgroup <- function(S, SigmaHat){
    sum(diag(S %*% solve(SigmaHat))) + log(det(SigmaHat)) - log(det(S)) - ncol(S)
  }
  f1 <- getFgroup(generated$Sigma, sigmaHat1)
  f2 <- getFgroup(generated2$Sigma, sigmaHat2)
  
  # power given effects by subgroups (fmin)
  ph2 <- semPower.postHoc(effect = list(f1, f2), effect.measure = 'F0', 
                          alpha = .05, N = list(500, 500), df = 3)
  
  valid <- round(ph$fmin - 2*lavres$fit['fmin'], 4) == 0 &&
    round(sum((c(f1, f2) - (lavres$summary$test$standard$stat.group / 500))^2), 4) == 0 &&
    round(ph2$power - ph$power, 4) == 0  
  
  # scalar invariance model
  generated <- semPower.genSigma(loadings = list(c(.5, .6, .7)), tau = c(0, 0, 0))
  generated2 <- semPower.genSigma(loadings = list(c(.5, .5, .7)), tau = c(0, .1, 0))
  lavres <- helper_lav(generated$modelTrue, 
                       list(generated$Sigma, generated2$Sigma),
                       sample.nobs = list(500, 500),
                       sample.mean = list(generated$mu, generated2$mu),
                       group.equal = c('loadings', 'intercepts'))
  sigmaHat1 <- lavaan::fitted(lavres$res)$`Group 1`$cov
  sigmaHat2 <- lavaan::fitted(lavres$res)$`Group 2`$cov
  muHat1 <- lavaan::fitted(lavres$res)$`Group 1`$mean
  muHat2 <- lavaan::fitted(lavres$res)$`Group 2`$mean
  
  ph <- semPower.postHoc(SigmaHat = list(sigmaHat1, sigmaHat2),
                         Sigma = list(generated$Sigma, generated2$Sigma),
                         muHat = list(muHat1, muHat2),
                         mu = list(generated$mu, generated2$mu),
                         alpha = .05, N = list(500, 500), df = 3)
  
  valid2 <- valid & round(ph$fmin - 2*lavres$fit['fmin'], 4) == 0
  
  if(valid2){
    print('test_multigroup: OK')
  }else{
    warning('Invalid')
  }
}

test_powerLav <- function(){
  mPop <- '
    f1 =~ .5*x1 + .6*x2 + .4*x3
    f2 =~ .7*x4 + .8*x5 + .3*x6
    f3 =~ .7*x7 + .8*x8 + .3*x9
    x1 ~~ .75*x1
    x2 ~~ .64*x2
    x3 ~~ .84*x3
    x4 ~~ .51*x4
    x5 ~~ .36*x5
    x6 ~~ .91*x6
    x7 ~~ .51*x7
    x8 ~~ .36*x8
    x9 ~~ .91*x9
    f1 ~~ 1*f1
    f2 ~~ 1*f2
    f3 ~~ 1*f3
    f1 ~~ .2*f2
    f1 ~~ .3*f3
    f2 ~~ .4*f3
  '
  mAna <- '
    f1 =~ x1 + x2 + x3
    f2 =~ x4 + x5 + x6
    f3 =~ x7 + x8 + x9
    f1 ~~ 0*f2
  '
  ph <- semPower.powerLav(type = 'post-hoc', comparison = 'saturated',
                          modelPop = mPop, modelH0 = mAna,
                          alpha = .05, N = 250)
  
  Sigma <- lavaan::fitted(lavaan::sem(mPop))$cov
  lavres <- helper_lav(mAna, Sigma)
  
  valid <- round(ph$power$fmin - 2*lavres$fit['fmin'], 4) == 0
  
  # use other comparison model
  mAna2 <- '
    f1 =~ x1 + x2 + x3
    f2 =~ x4 + x5 + x6
    f3 =~ x7 + x8 + x9
    f1 ~~ 0*f2
    f1 ~~ 0*f3
  '
  ph2 <- semPower.powerLav(type = 'post-hoc',
                           modelPop = mPop, modelH0 = mAna2, modelH1 = mAna,
                           alpha = .05, N = 250)
  
  lavres2 <- helper_lav(mAna2, Sigma)
  
  valid2 <- valid && 
    round(2*(lavres2$fit['fmin'] - lavres$fit['fmin']) - ph2$power$fmin, 4) == 0
  
  # consistency with f difference
  SigmaHat1 <- lavaan::fitted(lavres$res)$cov
  SigmaHat2 <- lavaan::fitted(lavres2$res)$cov
  f1 <- getF.Sigma(SigmaHat1, Sigma)
  f2 <- getF.Sigma(SigmaHat2, Sigma)
  deltaF <- f2 - f1
  ph3 <- semPower.postHoc(effect = deltaF, effect.measure = "F0",
                          alpha = .05, N = 250, df = 1)
  ph4 <- semPower.postHoc(effect = c(f1, f2), effect.measure = "F0",
                          alpha = .05, N = 250, 
                          df = c(lavres$res@test$standard$df, lavres2$res@test$standard$df))
  
  valid3 <- valid2 && 
    round(ph2$power$fmin - deltaF, 4) == 0 &&
    round(ph2$power$power - ph3$power, 4) == 0 &&
    round(ph3$power - ph4$power, 4) == 0
  
  if(valid3){
    print('test_powerLav: OK')
  }else{
    warning('Invalid')
  }
}

test_powerCFA <- function(){
  
  # simple model with restricted comparison model
  ph <- semPower.powerCFA(type = 'post-hoc', comparison = 'restricted',
                          nullEffect = 'cor=0',
                          nullWhich = c(1, 2),
                          Phi = .2, nIndicator = c(5, 6), loadM = .5,
                          alpha = .05, N = 250)
  
  lavres <- helper_lav(ph$modelH0, ph$Sigma)
  lavres2 <- helper_lav(ph$modelH1, ph$Sigma)
  
  ph2 <- semPower.powerLav('post-hoc', 
                           Sigma = ph$Sigma, modelH0 = ph$modelH0, 
                           modelH1 = ph$modelH1, fitH1model = FALSE,
                           alpha = .05, N = 250)
  
  ph3 <- semPower.postHoc(SigmaHat = ph$SigmaHat, Sigma = ph$Sigma,
                          df = 1, alpha = .05, N = 250)
  
  # saturated comparison model
  ph4 <- semPower.powerCFA(type = 'post-hoc', comparison = 'saturated',
                           nullEffect = 'cor=0',
                           nullWhich = c(1, 2),
                           Phi = .2, nIndicator = c(5, 6), loadM = .5,
                           alpha = .05, N = 250)
  
  valid <- round(lavres2$fit['fmin'], 4) == 0 &&
    round(2*lavres$fit['fmin'] - ph$power$fmin, 4) == 0 &&
    round(ph2$power$fmin - ph$power$fmin, 4) == 0 &&
    round(ph3$fmin - ph$power$fmin, 4) == 0 &&
    round(ph$power$power - ph2$power$power, 4) == 0 &&
    ph$power$df == 1 &&
    ph4$power$df == 44 &&
    round(ph$power$power - ph4$power$power, 4) != 0
  
  # phi matrix
  Phi <- matrix(c(
    c(1, .1, .2),
    c(.1, 1, .3),
    c(.2, .3, 1)
  ), byrow = T, ncol = 3)
  
  # nullWhich = c(1, 2)
  ph5 <- semPower.powerCFA(type = 'post-hoc', comparison = 'restricted', 
                           nullEffect = 'cor=0',
                           nullWhich = c(1, 2),
                           Phi = Phi, nIndicator = rep(3, 3), loadM = .5,
                           alpha = .05, N = 250)
  lavres3 <- helper_lav(ph5$modelH0, ph5$Sigma)
  par <- lavres3$par
  
  # nullWhich = c(2, 3)
  ph6 <- semPower.powerCFA(type = 'post-hoc', comparison = 'restricted', 
                           nullEffect = 'cor=0',
                           nullWhich = c(2, 3),
                           Phi = Phi, nIndicator = rep(3, 3), loadM = .5,
                           alpha = .05, N = 250)
  lavres6 <- helper_lav(ph6$modelH0, ph6$Sigma)
  par2 <- lavres6$par
  
  valid2 <- valid &&
    par[par$lhs == 'f1' & par$rhs == 'f2', 'est'] == 0 &&
    par[par$lhs == 'f2' & par$rhs == 'f3', 'est'] != 0 &&
    par2[par2$lhs == 'f1' & par2$rhs == 'f1', 'est'] != 0 &&
    par2[par2$lhs == 'f2' & par2$rhs == 'f3', 'est'] == 0 &&
    round(ph6$power$power, 4) > round(ph5$power$power, 4)
  
  # corx=cory, two equal
  ph7 <- semPower.powerCFA(type = 'post-hoc', comparison = 'restricted', 
                           nullEffect = 'corx=corz',
                           nullWhich = list(c(1,2), c(2, 3)),
                           Phi = Phi, nIndicator = rep(3, 3), loadM = .5,
                           alpha = .05, N = 250)
  lavres7 <- helper_lav(ph7$modelH0, ph7$Sigma)
  par3 <- lavres7$par

  # corx=cory, all equal
  ph8 <- semPower.powerCFA(type = 'post-hoc', comparison = 'restricted', 
                           nullEffect = 'corx=corz',
                           nullWhich = list(c(1,2), c(2, 3), c(1, 3)),
                           Phi = Phi, nIndicator = rep(3, 3), loadM = .5,
                           alpha = .05, N = 250)
  lavres8 <- helper_lav(ph8$modelH0, ph8$Sigma)
  par4 <- lavres8$par

  valid3 <- valid2 &&
    round(2*lavres7$fit['fmin'] - ph7$power$fmin, 4) == 0 &&
    round(par3[par3$lhs == 'f1' & par3$rhs == 'f2', 'std.all'] - par3[par3$lhs == 'f2' & par3$rhs == 'f3', 'std.all'], 4) == 0 &&
    round(var(
      c(par4[par4$lhs == 'f1' & par4$rhs == 'f2', 'std.all'], 
        par4[par4$lhs == 'f2' & par4$rhs == 'f3', 'std.all'], 
        par4[par4$lhs == 'f1' & par4$rhs == 'f3', 'std.all'])), 4) == 0 &&
    ph8$power$df == 2
  
  if(valid3){
    print('test_powerCFA: OK')
  }else{
    warning('Invalid')
  }
}

test_powerRegression <- function(){
  
  # restricted comparison model
  ph <- semPower.powerRegression(type = 'post-hoc', comparison = 'restricted',
                                 slopes = c(.2, .3), corXX = .4, nullWhich = 1,
                                 nIndicator = c(3, 3, 3), loadM = .5,
                                 alpha = .05, N = 250)
  lavres <- helper_lav(ph$modelH0, ph$Sigma)
  par <- lavres$par

  lavres2 <- helper_lav(ph$modelH1, ph$Sigma)
  par2 <- lavres2$par
  
  # power for second slope
  ph2 <- semPower.powerRegression(type = 'post-hoc', comparison = 'restricted',
                                  slopes = c(.2, .3), corXX = .4,  nullWhich = 2,
                                  nIndicator = c(3, 3, 3), loadM = .5,
                                  alpha = .05, N = 250)

  lavres3 <- helper_lav(ph2$modelH0, ph2$Sigma)
  
  # regression with 3 predictors, power for third slope, saturated comparison model
  corXX <- matrix(c(
    c(1.00, 0.20, 0.30),
    c(0.20, 1.00, 0.10),
    c(0.30, 0.10, 1.00)
  ), ncol = 3,byrow = TRUE)
  ph3 <- semPower.powerRegression(type = 'post-hoc', comparison = 'saturated',
                                  slopes = c(.2, .3, .4), corXX = corXX, nullWhich = 3,
                                  nIndicator = c(4, 3, 5, 4),
                                  loadM = c(.5, .5, .6, .7),
                                  alpha = .05, N = 250)
  lavres4 <- helper_lav(ph3$modelH0, ph3$Sigma)
  # same with restricted
  ph4 <- semPower.powerRegression(type = 'post-hoc', comparison = 'restricted',
                                  slopes = c(.2, .3, .4), corXX = corXX, nullWhich = 3,
                                  nIndicator = c(4, 3, 5, 4),
                                  loadM = c(.5, .5, .6, .7),
                                  alpha = .05, N = 250)
  lavres5 <- helper_lav(ph4$modelH1, ph4$Sigma)
  par4 <- lavres5$par

  # slope equality
  ph5 <- semPower.powerRegression(type = 'post-hoc', comparison = 'restricted',
                                  slopes = c(.2, .3, .4), corXX = corXX, 
                                  nullEffect = 'slopex=slopez', nullWhich = c(1, 2),
                                  nIndicator = c(4, 3, 5, 4),
                                  loadM = c(.5, .5, .6, .7),
                                  alpha = .05, N = 250)
  lavres6 <- helper_lav(ph5$modelH0, ph5$Sigma)
  par5 <- lavres6$par
  
  # slope equality of all slopes
  ph6 <- semPower.powerRegression(type = 'post-hoc', comparison = 'restricted',
                                  slopes = c(.2, .3, .4), corXX = corXX, 
                                  nullEffect = 'slopex=slopez', nullWhich = c(1, 2, 3),
                                  nIndicator = c(4, 3, 5, 4),
                                  loadM = c(.5, .5, .6, .7),
                                  alpha = .05, N = 250)
  
  lavres7 <- helper_lav(ph6$modelH0, ph6$Sigma)
  par6 <- lavres7$par

  # observed variables
  ph7 <- semPower.powerRegression(type = 'post-hoc', comparison = 'restricted',
                                 slopes = c(.2, .3), corXX = .4, nullWhich = 1,
                                 nIndicator = rep(1, 3), loadM = 1,
                                 alpha = .05, N = 250)
  
  
  lavres9 <- helper_lav(ph7$modelH1, ph7$Sigma)
  par7 <- lavres9$par
  lavres8 <- helper_lav(ph7$modelH0, ph7$Sigma)
  lavres10 <- helper_lav('x1 ~ p1*x2 + p2*x3 \n p1==0', ph7$Sigma)

  valid <- round(ph$power$fmin - 2*lavres$fit['fmin'], 4) == 0 &&
    round(ph2$power$fmin - 2*lavres3$fit['fmin'], 4) == 0 &&
    round(ph3$power$fmin - 2*lavres4$fit['fmin'], 4) == 0 &&
    round(par2[par2$lhs == 'f1' & par2$rhs == 'f2', 'std.all'] - .2, 4) == 0 &&
    round(par2[par2$lhs == 'f1' & par2$rhs == 'f3', 'std.all'] - .3, 4) == 0 &&
    round(par2[par2$lhs == 'f2' & par2$rhs == 'f3', 'std.all'] - .4, 4) == 0 &&
    lavres4$fit['df'] - ph3$power$df == 0 && 
    round(sum(par4[par4$op == '~~' & par4$lhs != par4$rhs, 'std.all'] - corXX[lower.tri(corXX)])^2, 4) == 0  &&
    round(sum(par4[par4$op == '~' & par4$lhs == 'f1', 'std.all'] - c(.2, .3, .4))^2, 4) == 0 && 
    round(var(par5[par5$lhs == 'f1' & par5$rhs %in% c('f2', 'f3'), 'est']), 4) == 0 &&
    round(2*lavres6$fit['fmin'] - ph5$power$fmin, 4) == 0 &&
    round(var(par6[par6$lhs == 'f1' & par6$rhs %in% c('f2', 'f3', 'f4'), 'est']), 4) == 0 &&
    round(2*lavres7$fit['fmin'] - ph6$power$fmin, 4) == 0 &&
    ph5$power$df == 1 && ph6$power$df == 2 &&
    round(ph7$power$fmin - 2*lavres8$fit['fmin'], 4) == 0 &&    
    round(ph7$power$fmin - 2*lavres10$fit['fmin'], 4) == 0 &&    
    round(par7[par7$lhs == 'f1' & par7$rhs == 'f2', 'std.all'] - .2, 4) == 0 &&
    round(par7[par7$lhs == 'f1' & par7$rhs == 'f3', 'std.all'] - .3, 4) == 0 &&
    round(par7[par7$lhs == 'f2' & par7$rhs == 'f3', 'std.all'] - .4, 4) == 0    

  if(valid){
    print('test_powerRegression: OK')
  }else{
    warning('Invalid')
  }
}

test_powerMediation <- function(){
  
  # simple mediation
  ph <- semPower.powerMediation(type = 'post-hoc', comparison = 'restricted',
                                bYX = .25, bMX = .3, bYM = .4,
                                nIndicator = c(3, 3, 3), loadM = .5,
                                alpha = .05, N = 250)

  lavres <- helper_lav(ph$modelH1, ph$Sigma)
  par <- lavres$par
  lavres2 <- helper_lav(ph$modelH0, ph$Sigma)
  par2 <- lavres2$par
  
  # same with saturated
  ph2 <- semPower.powerMediation(type = 'post-hoc', comparison = 'saturated',
                                bYX = .25, bMX = .3, bYM = .4,
                                nIndicator = c(3, 3, 3), loadM = .5,
                                alpha = .05, N = 250)

  # same with B
  B <- matrix(c(
                c(.00, .00, .00),
                c(.30, .00, .00),
                c(.25, .40, .00)
                ), byrow = TRUE, ncol = 3)
  ph3 <- semPower.powerMediation(type = 'post-hoc', comparison = 'restricted',
                                Beta = B, indirect = list(c(2,1), c(3,2)),
                                nIndicator = c(3, 3, 3), loadM = .5,
                                alpha = .05, N = 250)
  

  valid <- round(par[par$lhs == 'f2' & par$rhs == 'f1', 'std.all'] - .3, 4) == 0 &&
    round(par[par$lhs == 'f3' & par$rhs == 'f1', 'std.all'] - .25, 4) == 0 &&
    round(par[par$lhs == 'f3' & par$rhs == 'f2', 'std.all'] - .4, 4) == 0 &&
    round(par[par$lhs == 'ind', 'std.all'] - (.3*.4), 4) == 0  &&
    round(par2[par2$lhs == 'ind', 'std.all'], 4) == 0 &&
    round(2*lavres2$fit['fmin'] - ph$power$fmin, 4) == 0 &&
    round(ph2$power$fmin - ph$power$fmin, 4) == 0 &&
    ph$power$df == 1 && ph2$power$df == lavres2$res@Fit@test$standard$df &&
    round(ph3$power$fmin - ph$power$fmin, 4) == 0
  
  # more complex beta
  B <- matrix(c(
                c(.00, .00, .00, .00),
                c(.20, .00, .00, .00),
                c(.05, .30, .00, .00),
                c(.20, .10, .40, .00)
                ), byrow = TRUE, ncol = 4)
  ph4 <- semPower.powerMediation(type = 'post-hoc', comparison = 'restricted',
                                 Beta = B, indirect = list(c(2,1), c(3,2), c(4,3)),
                                 nIndicator = c(3, 3, 3, 3), loadM = .5,
                                 alpha = .05, N = 250)
  
  lavres3 <- helper_lav(ph4$modelH1, ph4$Sigma)
  par3 <- lavres3$par
  lavres4 <- helper_lav(ph4$modelH0, ph4$Sigma)
  par4 <- lavres4$par
  
  valid2 <- valid &&
    round(par3[par3$lhs == 'ind', 'std.all'] - .2*.3*.4, 4) == 0 &&
    round(2*lavres4$fit['fmin'] - ph4$power$fmin, 4) == 0 &&
    round(par4[par4$lhs == 'f2' & par4$rhs == 'f1', 'std.all'], 4) == 0

  # observed variables
  ph5 <- semPower.powerMediation(type = 'post-hoc', comparison = 'restricted',
                                bYX = .25, bMX = .3, bYM = .4,
                                nIndicator = rep(1, 3), loadM = 1,
                                alpha = .05, N = 250)
  
  lavres5 <- helper_lav(ph5$modelH0, ph5$Sigma)
  lavres5a <- helper_lav(ph5$modelH1, ph5$Sigma)
  par5 <- lavres5a$par  
  
  # mixed observed and latent, more complex beta
  B <- matrix(c(
    c(.00, .00, .00, .00),
    c(.20, .00, .00, .00),
    c(.05, .30, .00, .00),
    c(.20, .10, .40, .00)
  ), byrow = TRUE, ncol = 4)
  ph6 <- semPower.powerMediation(type = 'post-hoc', comparison = 'restricted',
                                 Beta = B, indirect = list(c(2,1), c(3,2), c(4,3)),
                                 nIndicator = c(3, 1, 3, 1), loadM = .5,
                                 alpha = .05, N = 250)
  
  lavres6 <- helper_lav(ph6$modelH1, ph6$Sigma)
  par6 <- lavres6$par  

  valid3 <- valid2 &&
    round(2*lavres5$fit['fmin'] - ph5$power$fmin, 4) == 0 && 
    round(par5[par5$lhs == 'ind', 'std.all'] - .3*.4, 4) == 0 && 
    round(par6[par6$lhs == 'ind', 'std.all'] - .2*.3*.4, 4) == 0 

  if(valid3){
    print('test_powerMediation: OK')
  }else{
    warning('Invalid')
  }
}

test_powerCLPM <- function(){
  
  # 2 waves, crossedX=0
  ph <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                           nWaves = 2,
                           stabilities = c(.8, .7), 
                           crossedEffects = c(.2, .1),
                           rXY = NULL, # diagonal
                           nullEffect = 'crossedx=0',
                           nIndicator = rep(3, 4), loadM = c(.5, .6, .5, .6),
                           metricInvariance = TRUE,
                           waveEqual = NULL,
                           alpha = .05, N = 250)
  
  lavres <- helper_lav(ph$modelH1, ph$Sigma)
  par <- lavres$par
  lavres2 <- helper_lav(ph$modelH0, ph$Sigma)
  par2 <- lavres2$par

  valid <- round(sum((par[par$op == '~', 'est'] - c(.8, .1, .2, .7))^2), 4) == 0 &&
    round(par[par$lhs == 'f1' & par$rhs == 'f2', 'est'], 4) == 0 &&
    round(par[par$lhs == 'f3' & par$rhs == 'f4', 'est'], 4) == 0 &&
    round(par2[par2$lhs == 'f4' & par$rhs == 'f1', 'est'], 4) == 0 &&
    round(2*lavres2$fit['fmin'] - ph$power$fmin, 4) == 0 &&
    !any(!(par2[par2$lhs == 'f1' & par2$op == '=~', 'label'] == par2[par2$lhs == 'f3' & par2$op == '=~', 'label'])) &&
    !any(!(par2[par2$lhs == 'f2' & par2$op == '=~', 'label'] == par2[par2$lhs == 'f4' & par2$op == '=~', 'label']))  

  # 2 waves, crossedY=0
  ph2 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                            nWaves = 2,
                            stabilities = c(.8, .7), 
                            crossedEffects = c(.2, .1),
                            rXY = NULL, # diagonal
                            nullEffect = 'crossedy=0',
                            nIndicator = rep(3, 4), loadM = c(.5, .6, .5, .6),
                            metricInvariance = TRUE,
                            waveEqual = NULL,
                            alpha = .05, N = 250)
  
  lavres3 <- helper_lav(ph2$modelH0, ph2$Sigma)
  par3 <- lavres3$par
  
  # 2 waves, stabX=0
  ph3 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                            nWaves = 2,
                            stabilities = c(.8, .7), 
                            crossedEffects = c(.2, .1),
                            rXY = NULL, # diagonal
                            nullEffect = 'stabX = 0',
                            nIndicator = rep(3, 4), loadM = c(.5, .6, .5, .6),
                            metricInvariance = TRUE,
                            waveEqual = NULL,
                            alpha = .05, N = 250)
  
  lavres4 <- helper_lav(ph3$modelH0, ph3$Sigma)
  par4 <- lavres4$par  

  # 2 waves, stabY=0
  ph4 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                            nWaves = 2,
                            stabilities = c(.8, .7), 
                            crossedEffects = c(.2, .1),
                            rXY = NULL, # diagonal
                            nullEffect = 'stabY = 0',
                            nIndicator = rep(3, 4), loadM = c(.5, .6, .5, .6),
                            metricInvariance = TRUE,
                            waveEqual = NULL,
                            alpha = .05, N = 250)
  
  lavres5 <- helper_lav(ph4$modelH0, ph4$Sigma)
  par5 <- lavres5$par    

  # 2 waves, stabX=stabY
  ph5 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                            nWaves = 2,
                            stabilities = c(.8, .7), 
                            crossedEffects = c(.2, .1),
                            rXY = NULL, # diagonal
                            nullEffect = 'stabX=stabY',
                            nIndicator = rep(3, 4), loadM = c(.5, .6, .5, .6),
                            metricInvariance = TRUE,
                            waveEqual = NULL,
                            alpha = .05, N = 250)
  
  lavres6 <- helper_lav(ph5$modelH0, ph5$Sigma)
  par6 <- lavres6$par    
  
  # 2 waves, crossedX=crossedY
  ph6 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                            nWaves = 2,
                            stabilities = c(.8, .7), 
                            crossedEffects = c(.2, .1),
                            rXY = NULL, # diagonal
                            nullEffect = 'crossedX=crossedY',
                            nIndicator = rep(3, 4), loadM = c(.5, .6, .5, .6),
                            metricInvariance = TRUE,
                            waveEqual = NULL,
                            alpha = .05, N = 250)
  
  lavres7 <- helper_lav(ph6$modelH0, ph6$Sigma)
  par7 <- lavres7$par     
  
  valid2 <- valid && 
    round(2*lavres3$fit['fmin'] - ph2$power$fmin, 4) == 0 &&
    round(par3[par3$lhs == 'f3' & par3$rhs == 'f2', 'est'], 4) == 0 &&
    round(2*lavres4$fit['fmin'] - ph3$power$fmin, 4) == 0 &&
    round(par4[par4$lhs == 'f3' & par4$rhs == 'f1', 'est'], 4) == 0 &&
    round(2*lavres5$fit['fmin'] - ph4$power$fmin, 4) == 0 &&
    round(par5[par5$lhs == 'f4' & par5$rhs == 'f2', 'est'], 4) == 0 &&
    round(2*lavres6$fit['fmin'] - ph5$power$fmin, 4) == 0 &&
    round(par6[par6$lhs == 'f4' & par6$rhs == 'f2', 'est'] - par6[par6$lhs == 'f3' & par6$rhs == 'f1', 'est'], 4) == 0 &&
    round(2*lavres7$fit['fmin'] - ph6$power$fmin, 4) == 0 &&
    round(par7[par7$lhs == 'f4' & par7$rhs == 'f1', 'est'] - par7[par7$lhs == 'f3' & par7$rhs == 'f2', 'est'], 4) == 0
    
  # 2 waves, crossedX=crossedY, cor XY
  ph7 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                            nWaves = 2,
                            stabilities = c(.8, .7), 
                            crossedEffects = c(.2, .1),
                            rXY = c(.3, .2), 
                            nullEffect = 'crossedX=crossedY',
                            nIndicator = rep(3, 4), loadM = c(.5, .6, .5, .6),
                            metricInvariance = TRUE,
                            waveEqual = NULL,
                            alpha = .05, N = 250)
  
  lavres8 <- helper_lav(ph7$modelH1, ph7$Sigma)
  par8 <- lavres8$par     

  # 2 waves, crossedX=crossedY, cor XY, no invariance
  ph8 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                            nWaves = 2,
                            stabilities = c(.8, .7), 
                            crossedEffects = c(.2, .1),
                            rXY = c(.3, .2), 
                            nullEffect = 'crossedX=crossedY',
                            nIndicator = rep(3, 4), loadM = c(.4, .5, .6, .7),
                            metricInvariance = FALSE,
                            waveEqual = NULL,
                            alpha = .05, N = 250)
  
  lavres9 <- helper_lav(ph8$modelH1, ph8$Sigma)
  par9 <- lavres9$par   

  valid3 <- valid2 &&
    round(par8[par8$lhs == 'f1' & par8$rhs == 'f2', 'est'] - .3, 4) == 0 && 
    round(par8[par8$lhs == 'f3' & par8$rhs == 'f4', 'est'] - .2, 4) == 0 &&
    round(par9[par9$lhs == 'f1' & par9$rhs == 'f2', 'est'] - .3, 4) == 0 && 
    round(par9[par9$lhs == 'f3' & par9$rhs == 'f4', 'est'] - .2, 4) == 0 &&
    round(sum((par9[par9$op == '~', 'est'] - c(.8, .1, .2, .7))^2), 4) == 0 &&
    !any(!(par9[par9$lhs == 'f1' & par9$op == '=~', 'est'] != par9[par9$lhs == 'f3' & par9$op == '=~', 'est'])) &&
    !any(!(par9[par9$lhs == 'f2' & par9$op == '=~', 'est'] != par9[par9$lhs == 'f4' & par9$op == '=~', 'est']))
    
  # 3 waves, wave invariant stabilities/crossed effects, crossedX=0 for first effect, 
  ph10 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                           nWaves = 3,
                           stabilities = c(.8, .7), 
                           crossedEffects = c(.2, .1),
                           rXY = c(.3, .2, .1),
                           nullEffect = 'crossedx=0',
                           nIndicator = rep(3, 6), loadM = c(.5, .6, .5, .6, .5, .6),
                           metricInvariance = TRUE,
                           nullWhich = NULL,
                           waveEqual = c('stabX','stabY','crossedX','crossedY'),
                           alpha = .05, N = 250)
  
  lavres10 <- helper_lav(ph10$modelH1, ph10$Sigma)
  par10 <- lavres10$par
  lavres11 <- helper_lav(ph10$modelH0, ph10$Sigma)
  par11 <- lavres11$par
  
  # 3 waves, no wave-equality constrains (=> less power)
  ph11 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                             nWaves = 3,
                             stabilities = c(.8, .7), 
                             crossedEffects = c(.2, .1),
                             rXY = c(.3, .2, .1),
                             nullEffect = 'crossedx=0',
                             nIndicator = rep(3, 6), loadM = c(.5, .6, .5, .6, .5, .6),
                             metricInvariance = TRUE,
                             nullWhich = 1,
                             waveEqual = NULL,
                             alpha = .05, N = 250)

  lavres12 <- helper_lav(ph11$modelH0, ph11$Sigma)

  # 3 waves, no wave invariant stabilities/crossed effects, crossedX=0 for second effect, 
  ph12 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                             nWaves = 3,
                             stabilities = list(c(.8, .7), c(.7, .6)), 
                             crossedEffects = list(c(.2, .1), c(.3, .1)),
                             rXY = c(.3, .2, .1),
                             nullEffect = 'crossedx=0',
                             nIndicator = rep(3, 6), loadM = c(.5, .6, .5, .6, .5, .6),
                             metricInvariance = TRUE,
                             nullWhich = 2,
                             waveEqual = NULL,
                             alpha = .05, N = 250)
  
  lavres13 <- helper_lav(ph12$modelH1, ph12$Sigma)
  par12 <- lavres13$par
  lavres14 <- helper_lav(ph12$modelH0, ph12$Sigma)
  par13 <- lavres14$par
  
  # 3 waves, no wave invariant stabilities/crossed effects, stabx equal, 
  ph13 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                             nWaves = 3,
                             stabilities = list(c(.8, .7), c(.7, .6)), 
                             crossedEffects = list(c(.2, .1), c(.3, .1)),
                             rXY = c(.3, .2, .1),
                             nullEffect = 'stabX',
                             nIndicator = rep(3, 6), loadM = c(.5, .6, .5, .6, .5, .6),
                             metricInvariance = TRUE,
                             nullWhich = NULL,
                             waveEqual = NULL,
                             alpha = .05, N = 250)
  
  lavres15 <- helper_lav(ph13$modelH0, ph13$Sigma)
  par14 <- lavres15$par
  
  # 3 waves, no wave invariant stabilities/crossed effects, staby equal, 
  ph14 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                             nWaves = 3,
                             stabilities = list(c(.8, .7), c(.7, .6)), 
                             crossedEffects = list(c(.2, .1), c(.3, .1)),
                             rXY = c(.3, .2, .1),
                             nullEffect = 'stabY',
                             nIndicator = rep(3, 6), loadM = c(.5, .6, .5, .6, .5, .6),
                             metricInvariance = TRUE,
                             nullWhich = NULL,
                             waveEqual = NULL,
                             alpha = .05, N = 250)
  
  lavres16 <- helper_lav(ph14$modelH0, ph14$Sigma)
  par15 <- lavres16$par  

  # 3 waves, no wave invariant stabilities/crossed effects, crossedx equal, 
  ph15 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                             nWaves = 3,
                             stabilities = list(c(.8, .7), c(.7, .6)), 
                             crossedEffects = list(c(.2, .1), c(.3, .1)),
                             rXY = c(.3, .2, .1),
                             nullEffect = 'crossedX',
                             nIndicator = rep(3, 6), loadM = c(.5, .6, .5, .6, .5, .6),
                             metricInvariance = TRUE,
                             nullWhich = NULL,
                             waveEqual = NULL,
                             alpha = .05, N = 250)
  
  lavres17 <- helper_lav(ph15$modelH0, ph15$Sigma)
  par16 <- lavres17$par

  # 3 waves, no wave invariant stabilities/crossed effects, crossedy equal, 
  ph16 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                             nWaves = 3,
                             stabilities = list(c(.8, .7), c(.7, .6)), 
                             crossedEffects = list(c(.2, .1), c(.3, .1)),
                             rXY = c(.3, .2, .1),
                             nullEffect = 'crossedY',
                             nIndicator = rep(3, 6), loadM = c(.5, .6, .5, .6, .5, .6),
                             metricInvariance = TRUE,
                             nullWhich = NULL,
                             waveEqual = NULL,
                             alpha = .05, N = 250)
  
  lavres18 <- helper_lav(ph16$modelH0, ph16$Sigma)
  par17 <- lavres18$par
  
  # 3 waves, wave-invariant stabilities/crossed effects, corxy equal, 
  ph17 <- semPower.powerCLPM(type = 'post-hoc', comparison = 'restricted',
                             nWaves = 3,
                             stabilities = list(c(.8, .7), c(.7, .6)), 
                             crossedEffects = list(c(.2, .1), c(.3, .1)),
                             rXY = c(.3, .2, .1),
                             nullEffect = 'corXY',
                             nIndicator = rep(3, 6), loadM = c(.5, .6, .5, .6, .5, .6),
                             metricInvariance = TRUE,
                             nullWhich = NULL,
                             waveEqual = c('stabX','stabY','crossedX','crossedY'),
                             alpha = .05, N = 250)
  
  lavres19 <- helper_lav(ph17$modelH0, ph17$Sigma)
  par18 <- lavres19$par

  valid4 <- valid3 && 
    round(sum((par10[par10$lhs %in% c('f3', 'f4') & par10$op == '~', 'est'] - par10[par10$lhs %in% c('f5', 'f6') & par10$op == '~', 'est'])^2), 4) == 0 &&  round(2*lavres11$fit['fmin'] - ph10$power$fmin, 4) == 0 &&
    round(par11[par11$lhs == 'f4' & par11$rhs == 'f1', 'est'], 4) ==  0 &&
    round(par11[par11$lhs == 'f6' & par11$rhs == 'f3', 'est'], 4) ==  0 &&
    ph10$power$fmin > ph11$power$fmin &&
    lavres12$res@Fit@test$standard$df < lavres11$res@Fit@test$standard$df &&
    round(sum((par12[par12$op == '~', 'est'] - c(.8, .3, .2, .7, .7, .1, .1, .6))^2), 4) == 0 && 
    round(2*lavres14$fit['fmin'] - ph12$power$fmin, 4) == 0 &&
    round(par13[par13$lhs == 'f4' & par13$rhs == 'f1', 'est'], 4) != 0 &&
    round(par13[par13$lhs == 'f6' & par13$rhs == 'f3', 'est'], 4) == 0 &&
    round(2*lavres15$fit['fmin'] - ph13$power$fmin, 4) == 0 &&
    round(par14[par14$lhs == 'f3' & par14$rhs == 'f1', 'est'] - par14[par14$lhs == 'f5' & par14$rhs == 'f3', 'est'], 4) == 0 &&
    round(par15[par15$lhs == 'f4' & par15$rhs == 'f2', 'est'] - par15[par15$lhs == 'f6' & par15$rhs == 'f4', 'est'], 4) == 0 &&
    round(par16[par16$lhs == 'f4' & par16$rhs == 'f1', 'est'] - par16[par16$lhs == 'f6' & par16$rhs == 'f3', 'est'], 4) == 0 &&
    round(2*lavres17$fit['fmin'] - ph15$power$fmin, 4) == 0 &&
    round(par17[par17$lhs == 'f3' & par17$rhs == 'f2', 'est'] - par17[par17$lhs == 'f5' & par17$rhs == 'f4', 'est'], 4) == 0 &&
    round(sum((ph16$Sigma - ph15$Sigma)^2), 4) == 0 &&   
    round(sum((ph14$Sigma - ph13$Sigma)^2), 4) == 0 &&   
    round(sum((ph13$Sigma - ph16$Sigma)^2), 4) == 0 &&
    round(2*lavres19$fit['fmin'] - ph17$power$fmin) == 0 &&
    round(sum((par18[par18$lhs == 'f3' & par18$rhs == 'f4', 'est'] - par18[par18$lhs == 'f5' & par18$rhs == 'f6', 'est'])^2), 4) == 0  
  
  if(valid4){
    print('test_powerCLPM: OK')
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
  test_generateSigmaB()
  test_genPhi()
  test_multigroup()
  test_powerLav()
  test_powerCFA()
  test_powerRegression()
  test_powerMediation()
  test_powerCLPM()
}

test_all()
