#' semPower.showPlot
#'
#' Shows a plot showing central and non-central chi-square distribution.
#' 
#' @param chiCrit critical chi-square, e. g. `qchisq(alpha, df, ncp = 0, lower.tail = FALSE)`
#' @param ncp non-centrality parameter under H1
#' @param df degrees of freedom
#' @param linewidth linewidth
#' @param showLabels whether to add labels
#' @importFrom stats qchisq dchisq
#' @importFrom graphics plot abline lines polygon text
#' @importFrom grDevices rgb
#' @export
semPower.showPlot <- function(chiCrit, ncp, df, linewidth = 1, showLabels = TRUE){
  
  # define central and non-central chi
  maxvalue <- qchisq(.99999, df, ncp)
  minvalue <- max(.1, qchisq(.00001, df, 0))
  x <- seq(minvalue, maxvalue, length = 1000)
  
  xchi <- dchisq(x, df, ncp = 0)
  xncchi <- dchisq(x, df, ncp)
  
  plot(x, 
       xchi, 
       type = "l", 
       lty = 1,
       lwd = linewidth,
       col = "red",
       xlab = "Chi-Square",
       ylab = "Density",
       ylim = c(0, max(xchi)),
       xlim = c(minvalue, maxvalue)
  )

  # add noncentral chi
  lines(x, xncchi, lty = 2, lwd = linewidth, col = 'blue')
  
  # add critical value
  abline(v = chiCrit, col = 1, lty = 3, lwd = 1)

  # shade alpha error
  color <- rgb(255, 0, 0, alpha = 70, maxColorValue = 255)
  lb <- chiCrit
  ub <- maxvalue
  i <- x >= lb & x <= ub
  polygon(c(lb, x[i], ub), c(0, xchi[i], 0), col = color, border = NA)

  # shade beta error
  color <- rgb(0, 0, 255, alpha = 70, maxColorValue = 255)
  lb <- minvalue
  ub <- chiCrit
  i <- x >= lb & x <= ub
  polygon(c(lb, x[i], ub), c(0, xncchi[i], 0), col = color, border = NA)
  
  # add labels
  if(showLabels){
    text(expression(paste("Central ", chi^2)), x = df, y = max(xchi), pos = 4, cex = .8, col = 'red')
    text(expression(paste("Non-Central ", chi^2)), x = (df + ncp), y = .9*max(xncchi), pos = 4, cex = .8, col = 'blue')
    text(expression(paste("Critical ", chi^2)), x = chiCrit, y = .8*max(xchi), pos = 4, cex = .8)
    text(expression(paste(alpha, "-error")), x = chiCrit, y = dchisq(qchisq(.985, df), df), pos = 4, cex = .8, col = 'red')
    text(expression(paste(beta, "-error")), x = chiCrit, y = dchisq(qchisq(.985, df), df), pos = 2, cex = .8, col = 'blue')
  }
}



#' semPower.powerPlot
#'
#' Shows a plot showing power as function of N given the output of a power analysis. 
#' 
#' @param semPowerRes results of a semPower analysis 
#' @param ... other parameters passed to [semPower.powerPlot.byN()]
#' @return powerplot
#' @examples
#' \dontrun{
#' # perform a power analysis
#' powerCFA <- semPower.powerCFA(type = 'post-hoc', alpha = .05, N = 300,
#'                               Phi = .15, nIndicator = c(5, 4), loadM = c(.5, .6))
#' # show plot
#' semPower.powerPlot(powerCFA)
#' }
#' @export
semPower.powerPlot <- function(semPowerRes, ...){
  if(!is.null(semPowerRes[['simulated']]) && semPowerRes[['simulated']]){
    semPower.powerPlot.byN(effect.measure = 'f0', effect = semPowerRes[['simFmin']], 
                           alpha = semPowerRes[['alpha']], df = semPowerRes[['simDf']], 
                           ...)
  }else{
    semPower.powerPlot.byN(effect.measure = 'f0', effect = semPowerRes[['fmin']], 
                           alpha = semPowerRes[['alpha']], df = semPowerRes[['df']], 
                           ...)
  }
}


#' semPower.powerPlot.byN
#'
#' Shows a plot showing power as function of N for a given effect and alpha.
#' 
#' @param effect effect size specifying the discrepancy between H0 and H1
#' @param effect.measure type of effect, one of `"F0"`, `"RMSEA"`, `"Mc"`, `"GFI"`, `"AGFI"`
#' @param alpha alpha error
#' @param df the model degrees of freedom
#' @param p the number of observed variables, required for `effect.measure = "GFI"` and `effect.measure = "AGFI"`
#' @param SigmaHat model implied covariance matrix. Use in conjunction with `Sigma`` to define `effect` and `effect.measure`. 
#' @param Sigma population covariance matrix. Use in conjunction with `SigmaHat` to define `effect` and `effect.measure`.
#' @param power.min minimum power, must not be smaller than `alpha`.
#' @param power.max maximum power
#' @param steps number of steps 
#' @param linewidth linewidth
#' @return powerplot
#' @examples
#' \dontrun{
#' semPower.powerPlot.byN(effect = .05, effect.measure = "RMSEA", 
#'                        alpha = .05, power.min = .05, power.max = .99, df = 200)
#' }
#' @importFrom stats smooth.spline
#' @importFrom graphics plot
#' @export
semPower.powerPlot.byN <- function(effect = NULL, effect.measure = NULL,
                                   alpha, df, p = NULL,
                                   SigmaHat = NULL, Sigma = NULL, 
                                   power.min = alpha, power.max = .99,
                                   steps = 50, linewidth = 1){
  
  if(!is.null(effect.measure)) effect.measure <- toupper(effect.measure)
  
  validateInput('powerplot.byN', effect = effect, effect.measure = effect.measure,
                alpha = alpha, beta = NULL, power = NULL, abratio = NULL,
                N = NULL, df = df, p = p,
                SigmaHat = SigmaHat, Sigma = Sigma,
                power.min = power.min, power.max = power.max,
                steps = steps, linewidth = linewidth)
  
  # do this here instead of in validator
  if(power.min < alpha){
    power.min <- alpha
    warning("power cannot be lower than alpha, setting power.min=alpha")
  }
  
  # determine N
  power <- seq(power.min, power.max, length.out = steps)
  N <- power.act <- array()
  for(i in 1:length(power)){
    ap <- semPower.aPriori(effect, effect.measure, alpha = alpha, power = power[[i]], df = df, p = p)
    N[[i]] <- ap[['requiredN']]
    power.act[[i]] <- ap[['impliedPower']]
  }
  
  plot(
    smooth.spline(N, power.act), 
    type = "l", 
    lty = 1,
    lwd = linewidth,
    main = paste('Power for', effect.measure, '=', effect),
    xlab = "N",
    ylab = "Power",
    ylim = c(0, 1),
    xlim = c(max(0, min(N) - 10), max(N))
  )
  
} 




#' semPower.powerPlot.byEffect
#'
#' Shows a plot showing power as function of the effect for a given N and alpha.
#' 
#' @param effect.measure type of effect, one of `"F0"`, `"RMSEA"`, `"Mc"`, `"GFI"`, `"AGFI"`
#' @param alpha alpha error
#' @param N the number of observations
#' @param df the model degrees of freedom
#' @param p the number of observed variables, required for `effect.measure = "GFI"` and `effect.measure = "AGFI"`
#' @param effect.min minimum effect
#' @param effect.max maximum effect
#' @param steps number of steps 
#' @param linewidth linewidth
#' @return powerplot
#' @examples
#' \dontrun{
#' semPower.powerPlot.byEffect(effect.measure = "RMSEA", alpha = .05, 
#'                             N = 500, effect.min = .01, effect.max = .15, df = 200)
#' }
#' @importFrom stats smooth.spline
#' @importFrom graphics plot
#' @export
semPower.powerPlot.byEffect <- function(effect.measure = NULL,
                                        alpha, N, df, p = NULL,
                                        effect.min = NULL, effect.max = NULL,
                                        steps = 50, linewidth = 1){
  
  effect.measure <- toupper(effect.measure)
  
  validateInput('powerplot.byEffect', effect = NULL, effect.measure = effect.measure,
                alpha = alpha, beta = NULL, power = NULL, abratio = NULL,
                N = N, df = df, p = p,
                SigmaHat = NULL, Sigma = NULL,
                effect.min = effect.min, effect.max = effect.max,
                steps = steps, linewidth = linewidth)
  
  # determine effect
  effect <- seq(effect.min, effect.max, length.out = steps)
  power <- array()
  for(i in 1:length(effect)){
    ph <- semPower.postHoc(effect[[i]], effect.measure, alpha = alpha, N = N, df = df, p = p)
    power[[i]] <- ph[['power']]
  }
  
  plot(smooth.spline(effect, power), 
       type = "l", 
       lty = 1,
       lwd = linewidth,
       main = paste('Power for', effect.measure, 'with N = ', N),
       xlab = "Effect",
       ylab = "Power",
       ylim = c(0, 1),
       xlim = c(max(0, min(effect) - .1 * effect), max(effect))
  )
  
} 
