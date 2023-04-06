### TODO result should be accessed using [[x]] instead of $x 

#' getFormattedResults
#'
#' Return data.frame containing formatted results.
#' 
#' @param type type of power analysis
#' @param result result object (list)
#' @param digits number of significant digits
#' @return data.frame
getFormattedResults <- function(type, result, digits = 6){
  
  ########### common across power types
  
  if(!is.null(result$srmr) && !is.null(result$gfi)){
    rows <- c('F0', 'RMSEA', 'SRMR', 'Mc', 'GFI', 'AGFI', 'CFI', '')
    head <- data.frame(rows)
    values <- c(result$fmin, result$rmsea, result$srmr, result$mc, result$gfi, result$agfi, result$cfi)
  }else if(!is.null(result$gfi)){
    rows <- c('F0', 'RMSEA', 'Mc', 'GFI', 'AGFI', '')
    head <- data.frame(rows)
    values <- c(result$fmin, result$rmsea, result$mc, result$gfi, result$agfi)
  }else{
    rows <- c('F0', 'RMSEA', 'Mc', '')
    head <- data.frame(rows)
    values <- c(result$fmin, result$rmsea, result$mc)
  }
  
  head$values <- c(formatC(values, format='f', digits = digits), '')
  
  
  ########### a-priori
  
  if(type == 'a-priori'){
    
    ifelse(length(result$requiredN.g) == 1, rows <- c('df', 'Required Num Observations', ''), rows <- c('df', 'Required Num Observations', ' ', ''))
    body <- data.frame(rows)
    ifelse(length(result$requiredN.g) == 1, 
           body$values <- c(result$df, result$requiredN, ''),
           body$values <- c(result$df, result$requiredN, paste0('(', paste(result$requiredN.g, collapse = ', '), ')'), '')
    )    
    
    rows <- c('Critical Chi-Square', 'NCP', 'Alpha', 'Beta', 'Power (1 - Beta)', 'Implied Alpha/Beta Ratio')
    foot <- data.frame(rows)
    
    v1 <- formatC(c(result$chiCrit, result$impliedNCP), format = 'f', digits = digits)
    v1 <- sapply(v1, substr, 1, digits + 2)
    
    v2 <- c(result$alpha, result$impliedBeta, result$impliedPower, result$impliedAbratio)
    # determine whether to use float or scientific number format
    v2.f <- rep('f', length(v2))
    v2.f[v2 < 1e-5 | v2 > 1e5] <- 'e'
    v2 <- sapply(seq_along(v2), function(y, z, i){formatC(x = v2[i], format = v2.f[i], digits = digits)}, y = v2, z = v2.f)
    foot$values <- c(v1, v2)
    
    # manually correct some quirks
    if(result$impliedBeta < 1e-5){
      foot$values[foot$rows == 'Power (1 - Beta)'] <- '> 0.9999'
    }
    
    
    out <- rbind(head, body, foot)
    rownames(out) <- colnames(out) <- NULL
    
  }
  
  
  ########### post-hoc
  
  if(type == 'post-hoc'){
    
    ifelse(!is.list(result$N), rows <- c('df', 'Num Observations'), rows <- c('df', 'Num Observations', ' '))
    body <- data.frame(rows)
    ifelse(!is.list(result$N), 
           body$values <- c(result$df, result$N), 
           body$values <- c(result$df, sum(unlist(result$N)), paste0('(', paste(result$N, collapse = ', '), ')'))
    )
    
    rows <- c('NCP', '', 'Critical Chi-Square')
    body2 <- data.frame(rows)
    v1 <- c(formatC(result$ncp, format = 'f', digits = digits), '',
            formatC(result$chiCrit, format = 'f', digits = digits)
    )
    body2$values <- sapply(v1, substr, 1, (digits + 2))
    
    rows <- c('Alpha', 'Beta', 'Power (1 - Beta)', 'Implied Alpha/Beta Ratio')
    foot <- data.frame(rows)
    v <- c(result$alpha, result$beta, result$power, result$impliedAbratio)
    # determine whether to use float or scientific number format
    v.f <- rep('f', length(v))
    v.f[v < 1e-5 | v > 1e5] <- 'e'
    foot$values <- sapply(seq_along(v),  function(y, z, i) {formatC(x = v[i], format = v.f[i], digits = digits)}, y = v, z = v.f)
    
    # manually correct some quirks
    if(result$beta < 1e-5){
      foot$values[foot$rows == 'Power (1 - Beta)'] <- '> 0.9999'
    }
    if(result$beta == 0){
      foot$values[foot$rows == 'Beta'] <- '< 1.00e-320'
      foot$values[foot$rows == 'Implied Alpha/Beta Ratio'] <- '> 1.00e-320'
    }
    
    out <- rbind(head, body, body2, foot)
    rownames(out) <- colnames(out) <- NULL
    
  }
  
  ########### compromise
  
  if(type == 'compromise'){
    
    ifelse(!is.list(result$N), 
           rows <- c('df','Num Observations', 'Desired Alpha/Beta Ratio', '', 'Critical Chi-Square'), 
           rows <- c('df','Num Observations', ' ', 'Desired Alpha/Beta Ratio', '', 'Critical Chi-Square')
    )
    body <- data.frame(rows)
    
    if(!result$bPrecisionWarning){
      sChiCrit <- substr(formatC(result$chiCrit, format = 'f', digits = digits), 1, digits + 2)
    }else{
      smax <- substr(formatC(result$max, format = 'f', digits = digits), 1, digits + 2)
      smin <- substr(formatC(result$min, format = 'f', digits = digits), 1, digits + 2)
      sChiCrit <- paste(smax,'< Chi-Square < ', smin)
    }
    if(!is.list(result$N)){
      body$values <- c(result$df, result$N, 
                       formatC(result$desiredAbratio, format = 'f', digits = digits), '', sChiCrit)
    }else{
      body$values <- c(result$df, sum(unlist(result$N)), paste0('(', paste(result$N, collapse = ', '), ')'), 
                       formatC(result$desiredAbratio, format = 'f', digits = digits), '', sChiCrit)
    }
    
    rows <- c('Implied Alpha', 'Implied Beta', 'Implied Power (1 - Beta)', 'Actual Alpha/Beta Ratio')
    foot <- data.frame(rows)
    v <- c(result$impliedAlpha, result$impliedBeta, result$impliedPower, result$impliedAbratio)
    # determine whether to use float or scientific number format
    v.f <- rep('f', length(v))
    v.f[v < 1e-5 | v > 1e5] <- 'e'
    foot$values <- sapply(seq_along(v), function(y, z, i) {formatC(x = v[i], format = v.f[i], digits = digits)}, y = v, z = v.f)
    
    # manually correct some quirks
    if(result$impliedBeta < 1e-5){
      foot$values[foot$rows == 'Implied Power (1 - Beta)'] <- '> 0.9999'
    }
    if(result$bPrecisionWarning){
      foot$values[foot$rows == 'Implied Beta'] <- '< 1.00e-240'
      foot$values[foot$rows == 'Implied Alpha'] <- '< 1.00e-320'
      foot$values[foot$rows == 'Actual Alpha/Beta Ratio'] <- ' '
    }
    
    out <- rbind(head, body, foot)
    rownames(out) <- colnames(out) <- NULL
    
    
  }
  
  out
}
