#' getFormattedResults
#'
#' Return data.frame containing formatted results.
#' 
#' @param type type of power analysis
#' @param result result object (list)
#' @param digits number of significant digits
#' @return data.frame
getFormattedResults <- function(type, result, digits = 6){
  
  isSimulated <- ifelse(!is.null(result[['simulated']]), result[['simulated']], FALSE)
  
  ########### common across power types
  
  if(!is.null(result[['srmr']]) && !is.null(result[['gfi']])){
    rows <- c('F0', 'RMSEA', 'SRMR', 'Mc', 'GFI', 'AGFI', 'CFI', '')
    head <- data.frame(rows)
    values <- result[c('fmin', 'rmsea', 'srmr', 'mc', 'gfi', 'agfi', 'cfi')]
    if(isSimulated) simValues <- result[c('meanFmin', 'simrmsea', 'simsrmr', 'simmc', 'simgfi', 'simagfi', 'simcfi')]
  }else if(!is.null(result[['gfi']])){
    rows <- c('F0', 'RMSEA', 'Mc', 'GFI', 'AGFI', '')
    head <- data.frame(rows)
    values <- result[c('fmin', 'rmsea', 'mc', 'gfi', 'agfi')]
    if(isSimulated) simValues <- result[c('meanFmin', 'simrmsea', 'simmc', 'simgfi', 'simagfi')]
  }else{
    rows <- c('F0', 'RMSEA', 'Mc', '')
    head <- data.frame(rows)
    values <- result[c('fmin', 'rmsea', 'mc')]
    if(isSimulated) simValues <- result[c('meanFmin', 'simrmsea', 'simmc')]
  }
  
  head$values <- c(formatC(unlist(values), format='f', digits = digits), '')
  if(isSimulated) head$simValues <- c(formatC(unlist(simValues), format='f', digits = digits), '')
  

  
  ########### a-priori
  
  if(type == 'a-priori'){
    
    ifelse(length(result[['requiredN.g']]) == 1, rows <- c('df', 'Required Num Observations', ''), rows <- c('df', 'Required Num Observations', ' ', ''))
    body <- data.frame(rows)
    if(length(result[['requiredN.g']]) == 1){
      body$values <- c(result[c('df', 'requiredN')], '')
      if(isSimulated) body$simValues <- c(result[c('simDf', 'simRequiredN')], '')
    }else{
      body$values <- c(result[c('df', 'requiredN')], paste0('(', paste(result[['requiredN.g']], collapse = ', '), ')'), '')
      if(isSimulated) body$simValues <- c(result[c('simDf', 'simRequiredN')], paste0('(', paste(result[['simRequiredN.g']], collapse = ', '), ')'), '')
    }

    rows <- c('Critical Chi-Square', 'NCP', 'Alpha', 'Beta', 'Power (1 - Beta)', 'Implied Alpha/Beta Ratio')
    foot <- data.frame(rows)
    
    v1 <- formatC(unlist(result[c('chiCrit', 'impliedNCP')]), format = 'f', digits = digits)
    v1 <- sapply(v1, substr, 1, digits + 2)
    if(isSimulated){
      simV1 <- formatC(unlist(result[c('simChiCrit', 'simImpliedNCP')]), format = 'f', digits = digits)      
      simV1 <- sapply(simV1, substr, 1, digits + 2)
    }
    v2 <- unlist(result[c('alpha', 'impliedBeta', 'impliedPower', 'impliedAbratio')])
    # determine whether to use float or scientific number format
    v2.f <- rep('f', length(v2))
    v2.f[v2 < 1e-5 | v2 > 1e5] <- 'e'
    v2 <- sapply(seq_along(v2), function(y, z, i){formatC(x = v2[i], format = v2.f[i], digits = digits)}, y = v2, z = v2.f)
    foot$values <- c(v1, v2)
    if(isSimulated){
      simV2 <- unlist(result[c('alpha', 'simImpliedBeta', 'simImpliedPower', 'simImpliedAbratio')])
      # determine whether to use float or scientific number format
      simV2.f <- rep('f', length(simV2))
      simV2.f[simV2 < 1e-5 | simV2 > 1e5] <- 'e'
      simV2 <- sapply(seq_along(simV2), function(y, z, i){formatC(x = simV2[i], format = simV2.f[i], digits = digits)}, y = simV2, z = simV2.f)
      foot$simValues <- c(simV1, simV2)
    }
    
    # manually correct some quirks
    if(result[['impliedBeta']] < 1e-5){
      foot$values[foot$rows == 'Power (1 - Beta)'] <- '> 0.9999'
    }
    if(isSimulated && result[['simImpliedPower']] == 1){
      foot$simValues[foot$rows == 'Beta'] <- '0'
      foot$simValues[foot$rows == 'Implied Alpha/Beta Ratio'] <- ''
    }
      
    out <- rbind(head, body, foot)
    rownames(out) <- colnames(out) <- NULL
    if(isSimulated) out <- rbind(c('', 'Analytical', 'Simulated'), c('', '', ''), out)
    
  }
  
  
  ########### post-hoc
  
  if(type == 'post-hoc'){
    
    ifelse(!is.list(result[['N']]), rows <- c('df', 'Num Observations'), rows <- c('df', 'Num Observations', ' '))
    body <- data.frame(rows)
    if(!is.list(result[['N']])){
      body$values <- result[c('df', 'N')]
      if(isSimulated) body$simValues <- result[c('simDf', 'N')]
    }else{
      body$values <- c(result[['df']], sum(unlist(result[['N']])), paste0('(', paste(result[['N']], collapse = ', '), ')'))
      if(isSimulated) body$simValues <- c(result[['simDf']], sum(unlist(result[['N']])), paste0('(', paste(result[['N']], collapse = ', '), ')'))
    }
    
    rows <- c('NCP', '', 'Critical Chi-Square')
    body2 <- data.frame(rows)
    v1 <- c(formatC(result[['ncp']], format = 'f', digits = digits), '',
            formatC(result[['chiCrit']], format = 'f', digits = digits)
    )
    body2$values <- sapply(v1, substr, 1, (digits + 2))
    if(isSimulated){
      simV1 <- c(formatC(result[['simNcp']], format = 'f', digits = digits), '',
              formatC(result[['simChiCrit']], format = 'f', digits = digits)
      )
      body2$simValues <- sapply(simV1, substr, 1, (digits + 2))
    }
    
    rows <- c('Alpha', 'Beta', 'Power (1 - Beta)', 'Implied Alpha/Beta Ratio')
    foot <- data.frame(rows)
    v <- unlist(result[c('alpha', 'beta', 'power', 'impliedAbratio')])
    # determine whether to use float or scientific number format
    v.f <- rep('f', length(v))
    v.f[v < 1e-5 | v > 1e5] <- 'e'
    foot$values <- sapply(seq_along(v),  function(y, z, i) {formatC(x = v[i], format = v.f[i], digits = digits)}, y = v, z = v.f)
    if(isSimulated){
      simV <- unlist(result[c('alpha', 'simBeta', 'simPower', 'simImpliedAbratio')])
      # determine whether to use float or scientific number format
      simV.f <- rep('f', length(v))
      simV.f[simV < 1e-5 | simV > 1e5] <- 'e'
      foot$simValues <- sapply(seq_along(simV),  function(y, z, i) {formatC(x = simV[i], format = simV.f[i], digits = digits)}, y = simV, z = simV.f)
    }
    
    # manually correct some quirks
    if(result[['beta']] < 1e-5){
      foot$values[foot$rows == 'Power (1 - Beta)'] <- '> 0.9999'
    }
    if(result[['beta']] == 0){
      foot$values[foot$rows == 'Beta'] <- '< 1.00e-320'
      foot$values[foot$rows == 'Implied Alpha/Beta Ratio'] <- '> 1.00e-320'
    }
    if(isSimulated && result[['simPower']] == 1){
      foot$simValues[foot$rows == 'Beta'] <- '0'
      foot$simValues[foot$rows == 'Implied Alpha/Beta Ratio'] <- ''
    }
    
    out <- rbind(head, body, body2, foot)
    rownames(out) <- colnames(out) <- NULL
    if(isSimulated) out <- rbind(c('', 'Analytical', 'Simulated'), c('', '', ''), out)
    
  }
  
  ########### compromise
  
  if(type == 'compromise'){
    
    ifelse(!is.list(result[['N']]), 
           rows <- c('df','Num Observations', 'Desired Alpha/Beta Ratio', '', 'Critical Chi-Square'), 
           rows <- c('df','Num Observations', ' ', 'Desired Alpha/Beta Ratio', '', 'Critical Chi-Square')
    )
    body <- data.frame(rows)
    
    if(!result[['bPrecisionWarning']]){
      sChiCrit <- substr(formatC(result[['chiCrit']], format = 'f', digits = digits), 1, digits + 2)
    }else{
      smax <- substr(formatC(result[['max']], format = 'f', digits = digits), 1, digits + 2)
      smin <- substr(formatC(result[['min']], format = 'f', digits = digits), 1, digits + 2)
      sChiCrit <- paste(smax,'< Chi-Square < ', smin)
    }
    if(!is.list(result[['N']])){
      body$values <- c(result[['df']], result[['N']], 
                       formatC(result[['desiredAbratio']], format = 'f', digits = digits), '', sChiCrit)
    }else{
      body$values <- c(result[['df']], sum(unlist(result[['N']])), paste0('(', paste(result[['N']], collapse = ', '), ')'), 
                       formatC(result[['desiredAbratio']], format = 'f', digits = digits), '', sChiCrit)
    }
    
    rows <- c('Implied Alpha', 'Implied Beta', 'Implied Power (1 - Beta)', 'Actual Alpha/Beta Ratio')
    foot <- data.frame(rows)
    v <- unlist(result[c('impliedAlpha', 'impliedBeta', 'impliedPower', 'impliedAbratio')])
    # determine whether to use float or scientific number format
    v.f <- rep('f', length(v))
    v.f[v < 1e-5 | v > 1e5] <- 'e'
    foot$values <- sapply(seq_along(v), function(y, z, i) {formatC(x = v[i], format = v.f[i], digits = digits)}, y = v, z = v.f)
    
    # manually correct some quirks
    if(result[['impliedBeta']] < 1e-5){
      foot$values[foot$rows == 'Implied Power (1 - Beta)'] <- '> 0.9999'
    }
    if(result[['bPrecisionWarning']]){
      foot$values[foot$rows == 'Implied Beta'] <- '< 1.00e-240'
      foot$values[foot$rows == 'Implied Alpha'] <- '< 1.00e-320'
      foot$values[foot$rows == 'Actual Alpha/Beta Ratio'] <- ' '
    }
    
    out <- rbind(head, body, foot)
    rownames(out) <- colnames(out) <- NULL

  }
  
  out
}


#' getFormattedResults
#'
#' Return data.frame containing formatted results.
#' 
#' @param object result object (list)
#' @param digits number of significant digits
#' @return data.frame
getFormattedSimulationResults <- function(object, digits = 2){

  simOut <- data.frame(NA, ncol = 2)
  simOut[1, ] <- c('Convergence Rate (%) of the H0 model', formatC(c(100*object[['convergenceRate']]), format = 'f', digits = digits))
  simOut[2, ] <- c('', '')

  if(!is.null(object[['bChiSqH0']])){
    simOut <- rbind(simOut, c('Chi-Square Bias (%)', ''))
    simOut <- rbind(simOut, c('H0 Model', formatC(c(100*object[['bChiSqH0']]), format = 'f', digits = digits)))
    if(!is.null(object[['bChiSqH1']])){
      simOut <- rbind(simOut, c('H1 Model', formatC(c(100*object[['bChiSqH1']]), format = 'f', digits = digits)))
      simOut <- rbind(simOut, c('H0-H1 Difference', formatC(c(100*object[['bChiSqDiff']]), format = 'f', digits = digits)))
    }
    simOut <- rbind(simOut, c('', ''))
    simOut <- rbind(simOut, c('Chi-Square KS-Distance', ''))
    simOut <- rbind(simOut, c('H0 Model', formatC(c(object[['ksChiSqH0']]), format = 'f', digits = 6)))
    if(!is.null(object[['bChiSqH1']])){
      simOut <- rbind(simOut, c('H1 Model', formatC(c(object[['ksChiSqH1']]), format = 'f', digits = 6)))
      simOut <- rbind(simOut, c('H0-H1 Difference', formatC(c(object[['ksChiSqDiff']]), format = 'f', digits = 6)))
      simOut <- rbind(simOut, c('', ''))
      simOut <- rbind(simOut, c('Rejection Rate (%)', ''))
      simOut <- rbind(simOut, c('H0 Model', formatC(c(100*object[['rrH0']]), format = 'f', digits = digits)))
      simOut <- rbind(simOut, c('H1 Model', formatC(c(100*object[['rrH1']]), format = 'f', digits = digits)))
    }
    simOut <- rbind(simOut, c('', ''))
  }

  if(!is.null(object[['bLambda']])){
    simOut <- rbind(simOut, c('Average Parameter Bias (%) in the H1 Model:', ''))
    simOut <- rbind(simOut, c('Loadings', formatC(c(100*object[['bLambda']]), format = 'f', digits = digits)))
    if(!is.null(object[['bPsi']])){
      simOut <- rbind(simOut, c('Variances/Covariances', formatC(c(100*object[['bPsi']]), format = 'f', digits = digits)))
    }
    if(!is.null(object[['bBeta']])){
      simOut <- rbind(simOut, c('Regression Coefficients', formatC(c(100*object[['bBeta']]), format = 'f', digits = digits)))
    }
    if(!is.null(object[['bPhi']])){
      simOut <- rbind(simOut, c('Variances/Covariances', formatC(c(100*object[['bPhi']]), format = 'f', digits = digits)))
    }
  }
  
  colnames(simOut) <- rownames(simOut) <- NULL
  
  simOut
}
