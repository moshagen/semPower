#' showPlot
#'
#' show a plot showing central and non-central chi-square distribution
#' 
#' @param chiCrit critical chi-square, e.g. qchisq(alpha, df, ncp=0, lower.tail = F)
#' @param ncp non-centrality parameter under H1
#' @param df degrees of freedom
#' @param linewidth linewidth
#' @export
showPlot <- function(chiCrit, ncp, df, linewidth = 1){
  
  # create random data
  maxvalue <- qchisq(.99999, df, ncp)
  minvalue <- qchisq(.00001, df, 0)
  x <- seq(minvalue, maxvalue, length=1000)
  
  # define central and non-central chi
  xchi <- dchisq(x, df, ncp = 0)
  xncchi <- dchisq(x, df, ncp)
  
  # draw plot
  plot(x, 
       xchi, 
       type="l", 
       lty=1,
       lwd=linewidth,
       col="red",
       xlab="Chi-Square",
       ylab="Density",
       ylim=c(0,max(xchi)),
       xlim=c(minvalue,maxvalue)
  )

  # add noncentral chi
  lines(x, xncchi, lty=2, lwd=linewidth, col='blue')
  
  # add critical value
  abline(v=chiCrit,col=1,lty=1,lwd=1)

  # shade alpha error
  color <- rgb(255, 0, 0, alpha=70, maxColorValue=255)
  lb <- chiCrit
  ub <- maxvalue
  i <- x >= lb & x <= ub
  polygon(c(lb,x[i],ub), c(0,xchi[i],0), col=color, border = NA)
  
  
  # shade beta error
  color <- rgb(0, 0, 255, alpha=70, maxColorValue=255)
  lb = minvalue
  ub = chiCrit
  i <- x >= lb & x <= ub
  polygon(c(lb,x[i],ub), c(0,xncchi[i],0), col=color, border = NA)
  
}


