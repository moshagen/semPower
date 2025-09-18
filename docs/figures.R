library(semPower)
library(semPlot)
library(lavaan)
library(semptools)

#lavres <- sem(powerARMA$modelH1, sample.cov = powerARMA$Sigma, sample.nobs = 1000, sample.mean = powerARMA$mu, sample.cov.rescale = F)
#summary(lavres)
prepParam <- function(modelH1, Sigma, mu){
  n <- 1000
  if(is.list(Sigma)) n <- list(1000, 1000)
  lavres <- sem(modelH1, sample.cov = Sigma, sample.nobs = n, sample.mean = mu, sample.cov.rescale = F)
  spm <- semPlotModel(lavres)
  # omit fixed to zero covariances
  spm@Pars <- spm@Pars[!(spm@Pars$fixed == TRUE & spm@Pars$est == 0 & spm@Pars$edge == '<->'), ]
  # round 
  spm@Pars$est <- round(spm@Pars$est, 2)
  spm@Pars$est[spm@Pars$est < 0 & spm@Pars$est > -0.005] <- 0 
  spm@Pars$std <- round(spm@Pars$std, 2)
  spm@Pars$std[spm@Pars$std < 0 & spm@Pars$std > -0.005] <- 0 
  spm
}

##################################

##### CFA

##################################

powerCFA <- semPower.powerCFA(
  type = 'ph', alpha = .05, N = 800,
  Phi = .25,
  nullEffect = 'cor = 0',
  nIndicator = c(4, 3), loadM = c(.5, .6))

spm <- prepParam(powerCFA$modelH1, powerCFA$Sigma, mu = powerCFA$mu)

splot <- semPaths(spm, 'mod', 'par', style = 'ram', layout = 'tree', 
                  intercepts =  T, residuals = T, fixedStyle = c("#666666", 1), 
                  rotation  = 1)
p_pa3 <- set_curve(splot, c("f2 ~~ f1" = 2))
p_pa3$Arguments$edge.color[!is.na(p_pa3$Arguments$edge.color)] <- "#666666" # make darker
p_pa3$graphAttributes$Edges$color <- rep("#666666", length(p_pa3$graphAttributes$Edges$color))
p_pa3 <- set_edge_attribute(p_pa3,
                            values = c("f1 ~~ f2" = "red"),
                            attribute_name = "color")

#pdf(file = "docs/fig/cfa.pdf", width = 14, height = 8)
png(filename = "docs/fig/cfa.png", width = 14000, height = 7000, units = "px", pointsize = 12, res = 1200)
plot(p_pa3)
dev.off()


##################################

##### Bifactor

##################################


bfLoadings <- rep(.5, 10)
bfWhichFactors <- c(1, 2, 3)
loadings <- list(
  c(.30, .20, .10),       # specific factor 1
  c(.05, .10, .15),       # specific factor 2
  c(.20, .05, .15),       # specific factor 3
  c(.70, .75, .80, .85)   # covariate
)
Phi <- matrix(c(
  c(1, .3),   # bifactor
  c(.3, 1)    # covariate
), ncol = 2, byrow = TRUE)
powerBF <- semPower.powerBifactor(
  type = 'ph', alpha = .05, N = 800,
  nullEffect = 'cor = 0',
  nullWhich = c(1, 2),
  bfLoadings = bfLoadings,
  bfWhichFactors = bfWhichFactors,
  Phi = Phi,
  loadings = loadings
)

spm <- prepParam(powerBF$modelH1, powerBF$Sigma, mu = powerBF$mu)


layout <- matrix(
  c(NA, NA, NA, NA, NA, "f1", NA, NA, NA, NA, NA, NA, NA, NA, NA,
    "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", NA, "x11", "x12", "x13", "x14",
    NA, NA, "f2", NA, NA, "f3", NA, NA, "f4", NA, NA, NA, NA, "f5", NA
  ), byrow = TRUE, 3, 15)

splot <- semPaths(spm, 'mod', 'par', style = 'ram', layout = layout, 
                  intercepts =  T, residuals = F, fixedStyle = c("#666666", 1), 
                  rotation  = 1)

p_pa3 <- set_curve(splot, c(
  "f5 ~~ f1" = .4,
  "f5 ~~ f2" = -1.2,
  "f5 ~~ f3" = -1,
  "f5 ~~ f4" = -.4
))
# make darker
p_pa3$Arguments$edge.color[!is.na(p_pa3$Arguments$edge.color)] <- "#666666"
p_pa3$graphAttributes$Edges$color <- rep("#666666", length(p_pa3$graphAttributes$Edges$color))
p_pa3 <- set_edge_attribute(p_pa3,
                            values = c("f5 ~~ f1" = "red"),
                            attribute_name = "color")

#pdf(file = "docs/fig/bf.pdf", width = 14, height = 7)
png(filename = "docs/fig/bf.png", width = 14000, height = 7000, units = "px", pointsize = 12, res = 1200)
plot(p_pa3)
dev.off()



##################################

##### Regression

##################################

powerReg <- semPower.powerRegression(
  type = 'a-priori', alpha = .05, power = .80,
  slopes = c(.2, .3), 
  corXX = .4, 
  nullEffect = 'slope = 0',
  nullWhich = 1,
  nIndicator = c(3, 5, 4), 
  loadM = c(.5, .6, .7))

spm <- prepParam(powerReg$modelH1, powerReg$Sigma, mu = powerReg$mu)

layout <- matrix(
  c("x4", NA, NA, NA,
    "x5", "f2", NA, NA,
    "x6", NA, NA, NA,
    "x7", NA, NA, NA,
    "x8", NA, NA, "x1",
    NA, NA, "f1", "x2",
    "x9", NA, NA, "x3",
    "x10", NA, NA, NA,
    "x11", "f3", NA, NA,
    "x12", NA, NA, NA
  ), byrow = TRUE, 10, 4)


splot <- semPaths(spm, 'mod', 'par', style = 'ram', layout = layout, 
                  intercepts =  T, residuals = T, fixedStyle = c("#666666", 1), 
                  rotation  = 2)
p_pa3 <- set_curve(splot, c("f2 ~~ f3" = -1))
p_pa3$Arguments$edge.color[!is.na(p_pa3$Arguments$edge.color)] <- "#666666" # make darker
p_pa3$graphAttributes$Edges$color <- rep("#666666", length(p_pa3$graphAttributes$Edges$color))
#p_pa3 <- safe_resid_position(p_pa3)
p_pa3 <- rotate_resid(p_pa3, c('f1'= 165, 'x12' = 270))
p_pa3 <- set_edge_attribute(p_pa3,
                   values = c("f1 ~ f2" = "red"),
                   attribute_name = "color")

#pdf(file = "docs/fig/reg.pdf", width = 7, height = 4)
png(filename = "docs/fig/reg.png", width = 14000, height = 8000, units = "px", pointsize = 12, res = 1200)
plot(p_pa3)
dev.off()



##################################

##### Mediation

##################################

powerMed <- semPower.powerMediation(
  type = 'ph', alpha = .05, N = 800,
  bYX = .25, 
  bMX = .3, 
  bYM = .4,
  nullEffect = 'ind = 0',
  nIndicator = c(3, 4, 5),
  loadM = c(.5, .6, .7)
)

spm <- prepParam(powerMed$modelH1, powerMed$Sigma, mu = powerMed$mu)

layout <- matrix(
  c(NA,NA,"x4","x5","x6","x7",NA,NA,
    NA,NA,NA,"f2",NA,NA,NA,"x8",
    "x1",NA,NA,NA,NA,NA,NA,"x9", 
    "x2","f1",NA,NA,NA,NA,"f3","x10",
    "x3",NA,NA,NA,NA,NA,NA,"x11",
    NA,NA,NA,NA,NA,NA,NA,"x12" 
  ), byrow = TRUE, 6, 8)
splot <- semPaths(spm, 'mod', 'par', style = 'ram', layout = layout, 
                  intercepts =  T, residuals = T, fixedStyle = c("#666666", 1), 
                  rotation  = 3)
p_pa3 <- rotate_resid(splot, c('f1'= 140, 'f2' = 180, 'f3' = -140, 'x4' = 360, 'x5' = 360, 'x7' = 360))
p_pa3$Arguments$edge.color[!is.na(p_pa3$Arguments$edge.color)] <- "#666666" # make darker
p_pa3$graphAttributes$Edges$color <- rep("#666666", length(p_pa3$graphAttributes$Edges$color))
p_pa3 <- set_edge_attribute(p_pa3,
                            values = c("f3 ~ f2" = "red", "f2 ~ pf0201*f1" = "red"),
                            attribute_name = "color")
p_pa3 <- change_node_label(p_pa3,
                           c(
                             f1 = "X",
                             f2 = "M",
                             f3 = "Y"
                           ))

#pdf(file = "docs/fig/med.pdf", width = 7, height = 4)
png(filename = "docs/fig/med.png", width = 14000, height = 8500, units = "px", pointsize = 12, res = 1200)
plot(p_pa3)
dev.off()



##################################

##### Path

##################################

### path power analysis example
Beta <- matrix(c(
  c(.00, .00, .00, .00),       # f1 = .00*f1 + .00*f2 + .00*f3 + .00*f4
  c(.20, .00, .00, .00),       # f2 = .20*f1 + .00*f2 + .00*f3 + .00*f4
  c(.00, .30, .00, .00),       # f3 = .00*f1 + .30*f2 + .00*f3 + .00*f4
  c(.10, .00, .40, .00)        # f4 = .10*f1 + .00*f2 + .40*f3 + .00*f4
), byrow = TRUE, ncol = 4)
powerPath <- semPower.powerPath(
  type = 'ph', alpha = .05, N = 800,
  Beta = Beta,
  nullWhich = c(4, 1),
  nIndicator = c(3, 4, 5, 6),
  loadM = c(.7, .5, .6, .8),
)

spm <- prepParam(powerPath$modelH1, powerPath$Sigma, mu = powerPath$mu)

layout <- matrix(
  c(NA,NA,NA,NA,NA,NA,"x1","x2","x3",NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,"f1",NA,NA,NA,NA,NA,NA,NA,
    NA,NA,NA,"f2",NA,NA,NA,NA,NA,NA,NA,"f4",NA,NA,NA,
    NA,NA,NA,NA,NA,NA,NA,"f3",NA,NA,NA,NA,NA,NA,NA,
    "x4","x5","x6","x7",NA,NA,NA,NA,NA,"x13","x14","x15","x16","x17","x18",
    NA,NA,NA,NA,NA,"x8","x9","x10","x11","x12",NA,NA,NA,NA,NA
  ), byrow = TRUE, 6, 15)

splot <- semPaths(spm, 'mod', 'par', style = 'ram', layout = layout, 
                  intercepts =  T, residuals = T, fixedStyle = c("#666666", 1), 
                  rotation  = 1)
p_pa3$Arguments$edge.color[!is.na(p_pa3$Arguments$edge.color)] <- "#666666" # make darker
p_pa3$graphAttributes$Edges$color <- rep("#666666", length(p_pa3$graphAttributes$Edges$color))
p_pa3 <- set_edge_attribute(splot,
                            values = c("f4 ~~ f1" = "red"),
                            attribute_name = "color")
p_pa3 <- rotate_resid(p_pa3, c('f1'= 180, 'x1' = 0, 'x2' = 0, 'x4' = 180, 'x5' = 180, 'x6' = 180, 'x7' = 180, 'x8' = 180,
                               'x13' = 180, 'x14' = 180, 'x15' = 180, 'x16' = 180, 'x17' = 180, 'x18' = 180))

#pdf(file = "docs/fig/path.pdf", width = 7, height = 4)
png(filename = "docs/fig/path.png", width = 14000, height = 8000, units = "px", pointsize = 12, res = 1200)
plot(p_pa3)
dev.off()



### generate cov example
Beta <- matrix(c(
  c(.00, .00, .00, .00),
  c(.00, .00, .00, .00),
  c(.20, .30, .00, .00),
  c(.40, .50, .00, .00)
), byrow = TRUE, ncol = 4)
Psi <- matrix(c(
  c(1, .60, .00, .00),
  c(.60, 1, .00, .00),
  c(.00, .00, 1, .10),
  c(.00, .00, .10, 1)
), byrow = TRUE, ncol = 4)

powerPath <- semPower.powerPath(
  type = 'ph', alpha = .05, N = 800,
  Beta = Beta,
  Psi = Psi,
  nullWhich = c(3, 1),
  nIndicator = rep(3, 4),
  loadM = c(.7, .5, .6, .8),
  standardized = F
)

spm <- prepParam(powerPath$modelH1, powerPath$Sigma, mu = powerPath$mu)

splot <- semPaths(spm, 'mod', 'par', style = 'ram', layout = 'tree', 
                  intercepts =  T, residuals = T, fixedStyle = c("#666666", 1), structural = T, 
                  rotation  = 2)
p_pa3 <- set_curve(splot, c("f3 ~~ f4" = -1))
p_pa3 <- set_edge_label_position(p_pa3, c('f4 ~ f1' = .25, 'f3 ~ f2' = .25))
p_pa3 <- change_node_label(p_pa3,
                           c(
                             f1 = "X1",
                             f2 = "X2",
                             f3 = "X3",
                             f4 = "X4"
                           ))

p_pa3$Arguments$edge.color[!is.na(p_pa3$Arguments$edge.color)] <- "#666666" # make darker
p_pa3$graphAttributes$Edges$color <- rep("#666666", length(p_pa3$graphAttributes$Edges$color))

png(filename = "docs/fig/pathx.png", width = 4000, height = 4000, units = "px", pointsize = 12, res = 1200)
plot(p_pa3)
dev.off()



##################################

##### MI

##################################


powerMI <- semPower.powerMI(
  type = 'ph', alpha = .05,  N = list(800, 800),
  comparison = 'configural', 
  nullEffect = 'metric',
  loadings= list(list(c(.5, .5, .5, .5, .5)), # group 1
                 list(c(.7, .4, .6, .4, .7))) # group 2
)

spm <- prepParam(powerMI$modelH1, powerMI$Sigma, mu = powerMI$mu)

splot <- semPaths(spm, 'mod', 'par', style = 'ram', layout = 'tree', 
                  intercepts =  T, residuals = T, fixedStyle = c("#666666", 1), 
                  rotation  = 1, panelGroups = T, ask = F)
splot[[1]]$Arguments$edge.color[!is.na(splot[[1]]$Arguments$edge.color)] <- "#666666" # make darker
splot[[1]]$graphAttributes$Edges$color <- rep("#666666", length(splot[[1]]$graphAttributes$Edges$color))
splot[[2]]$Arguments$edge.color[!is.na(splot[[2]]$Arguments$edge.color)] <- "#666666" # make darker
splot[[2]]$graphAttributes$Edges$color <- rep("#666666", length(splot[[2]]$graphAttributes$Edges$color))
p_pa3 <- set_edge_attribute(splot[[1]],
                            values = c("f1 =~ x1" = "red1",
                                       "f1 =~ x2" = "red2",
                                       "f1 =~ x3" = "red3",
                                       "f1 =~ x4" = "red4",
                                       "f1 =~ x5" = "#5B0000" 
                            ),
                            attribute_name = "color")
p_pa3 <- change_node_label(p_pa3,
                           c(f1 = "f1g1",
                             x1 = "x1g1",
                             x2 = "x2g1",
                             x3 = "x3g1",
                             x4 = "x4g1",
                             x5 = "x5g1"))
p_pa3b <- set_edge_attribute(splot[[2]],
                            values = c("f1 =~ x1" = "red1",
                                       "f1 =~ x2" = "red2",
                                       "f1 =~ x3" = "red3",
                                       "f1 =~ x4" = "red4",
                                       "f1 =~ x5" = "#5B0000" 
                            ),
                            attribute_name = "color")
p_pa3b <- change_node_label(p_pa3b,
                           c(f1 = "f1g2",
                             x1 = "x1g2",
                             x2 = "x2g2",
                             x3 = "x3g2",
                             x4 = "x4g2",
                             x5 = "x5g2"))

#pdf(file = "docs/fig/mi.pdf", width = 7, height = 4)
png(filename = "docs/fig/mi.png", width = 14000, height = 6000, units = "px", pointsize = 12, res = 1200)
par(mfrow=c(1,2))
plot(p_pa3)
plot(p_pa3b)
dev.off()
par(mfrow=c(1,1))





##################################

##### LI

##################################

powerLI <- semPower.powerLI(
  type = 'ph', alpha = .05, N = 800,
  comparison = 'configural', 
  nullEffect = 'metric',
  loadings = list(
    c(.7, .7, .7, .7),   # factor 1 at T1
    c(.4, .5, .6, .6)    # factor 1 at T2
  ),
  Phi = .3
)

spm <- prepParam(powerLI$modelH1, powerLI$Sigma, mu = powerLI$mu)

splot <- semPaths(spm, 'mod', 'par', style = 'ram', layout = 'tree', 
                  intercepts =  T, residuals = T, fixedStyle = c("#666666", 1), 
                  rotation  = 1)
splot$Arguments$edge.color[!is.na(splot$Arguments$edge.color)] <- "#666666" # make darker
splot$graphAttributes$Edges$color <- rep("#666666", length(splot$graphAttributes$Edges$color))
p_pa3 <- set_edge_attribute(splot,
                            values = c("f1 =~ x1" = "red1",
                                       "f2 =~ x5" = "red1",
                                       "f1 =~ x2" = "red2",
                                       "f2 =~ x6" = "red2",
                                       "f1 =~ x3" = "red3",
                                       "f2 =~ x7" = "red3",
                                       "f1 =~ x4" = "red4",
                                       "f2 =~ x8" = "red4" 
                                       ),
                            attribute_name = "color")
p_pa3 <- change_node_label(p_pa3,
                           c(
                             x1 = "x11",
                             x2 = "x21",
                             x3 = "x31",
                             x4 = "x41",
                             x5 = "x12",
                             x6 = "x22",
                             x7 = "x32",
                             x8 = "x42"
                           ))

#pdf(file = "docs/fig/li.pdf", width = 7, height = 4)
png(filename = "docs/fig/li.png", width = 14000, height = 8000, units = "px", pointsize = 12, res = 1200)
plot(p_pa3)
dev.off()



##################################

##### Autoreg

##################################

powerAutoreg <- semPower.powerAutoreg(
  'ph', alpha = .05, N = 800,
  nWaves = 3,
  autoregEffects = c(.5, .7),  # x1->x2, x2->x3
  nullEffect = 'autoreg = 0',
  nullWhich = 1,
  nIndicator = c(3, 3, 3), loadM = .5
)
spm <- prepParam(powerAutoreg$modelH1, powerAutoreg$Sigma, mu = powerAutoreg$mu)

layout <- matrix(
  c(NA,"f1",NA,NA,"f2",NA,NA,"f3",NA,
    "x1","x2","x3","x4","x5","x6","x7","x8","x9",
    NA,NA,NA,NA,NA,NA,NA,NA,NA
  ), byrow = TRUE, 3, 9)

splot <- semPaths(spm, 'mod', 'est', style = 'ram', layout = layout, 
                  intercepts =  T, residuals = T, fixedStyle = c("#666666", 1), 
                  rotation  = 1)
p_pa3 <- rotate_resid(splot, c('x5'= -90))
p_pa3$Arguments$edge.color[!is.na(p_pa3$Arguments$edge.color)] <- "#666666" # make darker
p_pa3$graphAttributes$Edges$color <- rep("#666666", length(p_pa3$graphAttributes$Edges$color))
p_pa3 <- set_edge_attribute(p_pa3,
                            values = c("f2 ~ pf0201*f1" = "red"),
                            attribute_name = "color")
p_pa3 <- set_curve(splot, c("x1 ~~ x4" = -2, "x2 ~~ x5" = -2, "x3 ~~ x6" = -2, 
                            "x4 ~~ x7" = -2, "x5 ~~ x8" = -2, "x6 ~~ x9" = -2,
                            "x1 ~~ x7" = -1.8, "x2 ~~ x8" = -1.8, "x3 ~~ x9" = -1.8))


#pdf(file = "docs/fig/autoreg.pdf", width = 7, height = 4)
png(filename = "docs/fig/autoreg.png", width = 14000, height = 8000, units = "px", pointsize = 12, res = 1200)
plot(p_pa3)
dev.off()




##################################

##### ARMA

##################################

powerARMA <- semPower.powerARMA(
  type = 'ph', alpha = .05, N = 800,
  nWaves = 5,
  autoregLag1 = c(.5, .4, .3, .2),  # x1->x2, x2->x3, x3->x4, x4->x5 
  mvAvgLag1 = c(.4, .4, .4, .4),    # n1->x2, n2->x3, n3->x4, n4->x5 
  variances = c(1, 1, 1, 1, 1),     # n1, n2, n3, n4, n5
  waveEqual = c('var','mvAvg'),
  nullEffect = 'autoreg',
  nIndicator = rep(3, 5), loadM = .5
)


spm <- prepParam(powerARMA$modelH1, powerARMA$Sigma, mu = powerARMA$mu)

splot <- semPaths(spm, 'mod', 'par', style = 'ram', layout = 'tree', 
                  intercepts =  T, residuals = F, fixedStyle = c("#666666", 1), 
                  rotation  = 1)
p_pa3$Arguments$edge.color[!is.na(p_pa3$Arguments$edge.color)] <- "#666666" # make darker
p_pa3$graphAttributes$Edges$color <- rep("#666666", length(p_pa3$graphAttributes$Edges$color))
p_pa3 <- set_edge_attribute(splot,
                            values = c("f2 ~~ f1" = "red", "f3 ~~ f2" = "red", "f4 ~~ f3" = "red", "f5 ~~ f4" = "red"),
                            attribute_name = "color")

#pdf(file = "docs/fig/arma.pdf", width = 14, height = 7)
png(filename = "docs/fig/arma.png", width = 18000, height = 8000, units = "px", pointsize = 12, res = 1200)
plot(p_pa3)
dev.off()




##################################

##### CLPM

##################################

### definition example
powerCLPM <- semPower.powerCLPM(
  type = 'ph', alpha = .05, N = 800,
  nullEffect = 'crossedX = 0',
  nullWhich = 1,
  nWaves = 3,
  autoregEffects = list(c(.8, .7), c(.6, .5)),
  crossedEffects = list(c(.1, .2), c(.3, .4)),
  rXY = c(.3, .1, .1),
  nIndicator = rep(3, 6),
  loadM = .5
)

spm <- prepParam(powerCLPM$modelH1, powerCLPM$Sigma, mu = powerCLPM$mu)
layout <- matrix(
  c("x1","x2","x3",NA,"x7","x8","x9", NA, "x13","x14","x15",
    NA,"f1",NA,NA,NA,"f3",NA,NA,NA,"f5",NA,
    NA,"f2",NA,NA,NA,"f4",NA,NA,NA,"f6",NA,
    "x4","x5","x6",NA,"x10","x11","x12",NA, "x16","x17","x18"
  ), byrow = TRUE, 4, 11)

splot <- semPaths(spm, 'mod', 'std', style = 'ram', layout = layout, 
                  intercepts =  T, residuals = T, fixedStyle = c("#666666", 1), structural = T,
                  rotation  = 1)
p_pa3 <- set_curve(splot, c("f1 ~~ f2" = -1, "f3 ~~ f4" = 1, "f5 ~~ f6" = 1))
p_pa3$Arguments$edge.color[!is.na(p_pa3$Arguments$edge.color)] <- "#666666" # make darker
p_pa3$graphAttributes$Edges$color <- rep("#666666", length(p_pa3$graphAttributes$Edges$color))
p_pa3 <- set_edge_label_position(p_pa3, c('f3 ~ f2' = .25, 'f4 ~ f1' = .25, 'f5 ~ f4' = .25, 'f6 ~ f3' = .25))
p_pa3 <- change_node_label(p_pa3,
                           c(f1 = "X1",
                             f2 = "Y1",
                             f3 = "X2",
                             f4 = "Y2",
                             f5 = "X3",
                             f6 = "Y3"
                           ),
                           label.cex = .8)


png(filename = "docs/fig/clpmi.png", width = 10000, height = 5000, units = "px", pointsize = 12, res = 1200)
plot(p_pa3)
dev.off()


### power example
powerCLPM <- semPower.powerCLPM(
  type = 'ph', alpha = .05, N = 800,
  nullEffect = 'crossedX = 0',
  nWaves = 2,
  autoregEffects = c(.60, .70),
  crossedEffects = c(.10, .15),
  rXY = c(.3, .1),
  nIndicator = c(5, 3, 5, 3),
  loadM = c(.5, .6, .5, .6)
)

spm <- prepParam(powerCLPM$modelH1, powerCLPM$Sigma, mu = powerCLPM$mu)
layout <- matrix(
  c("x1","x2","x3","x4","x5","x9","x10","x11","x12","x13",
    NA,NA,"f1",NA,NA,NA,NA,"f3",NA,NA,
    NA,NA,"f2",NA,NA,NA,NA,"f4",NA,NA,
    NA,"x6","x7","x8",NA,NA,"x14","x15","x16",NA
  ), byrow = TRUE, 4, 10)

splot <- semPaths(spm, 'mod', 'std', style = 'ram', layout = layout, 
                  intercepts =  T, residuals = T, fixedStyle = c("#666666", 1), 
                  rotation  = 1)
p_pa3 <- set_curve(splot, c("x1 ~~ x9" = 1.5, "x2 ~~ x10" = 1.5, "x3 ~~ x11" = 1.5, "x4 ~~ x12" = 1.5, "x5 ~~ x13" = 1.5,
                            "x6 ~~ x14" = -1.5, "x7 ~~ x15" = -1.5, "x8 ~~ x16" = -1.5,
                            "f1 ~~ f2" = -1, "f3 ~~ f4" = 1))
p_pa3$Arguments$edge.color[!is.na(p_pa3$Arguments$edge.color)] <- "#666666" # make darker
p_pa3$graphAttributes$Edges$color <- rep("#666666", length(p_pa3$graphAttributes$Edges$color))
p_pa3 <- set_edge_attribute(p_pa3,
                            values = c("f4 ~~ f1" = "red"),
                            attribute_name = "color")
p_pa3 <- rotate_resid(p_pa3, c('x6'= 180, 'x14'= 180, 'x15'= 180, 'x16'= 180, 
                               'x1'= 0, 'x2'= 0, 'x3'= 0, 'x4'= 0, 'x5'= 0,
                               'x11'= 0, 'x12'= 0, 'x13'= 0))
p_pa3 <- set_edge_label_position(p_pa3, c('f3 ~ f2' = .25, 'f4 ~ f1' = .25))
p_pa3 <- change_node_label(p_pa3,
                           c(f1 = "X1",
                             f2 = "Y1",
                             f3 = "X2",
                             f4 = "Y2"),
                           label.cex = .8)


#pdf(file = "docs/fig/clpm.pdf", width = 7, height = 7)
png(filename = "docs/fig/clpm.png", width = 14000, height = 14000, units = "px", pointsize = 12, res = 1200)
plot(p_pa3)
dev.off()



##################################

##### RICLPM

##################################

powerRICLPM <- semPower.powerRICLPM(
  type = 'ph', alpha = .05, N = 800,
  nullEffect = 'corBXBY = 0',
  nWaves = 3,
  nullWhich = 1,
  autoregEffects = c(.5, .4),   # (X, Y)
  crossedEffects = c(.2, .1),   # (X, Y)
  rXY = NULL,                   # diagonal
  rBXBY = .30,                  
  nIndicator = rep(3, 6),
  loadM = rep(c(.5, .6), 3)
)
writeLines(powerRICLPM$modelH1)
spm <- prepParam(powerRICLPM$modelH1, powerRICLPM$Sigma, mu = powerRICLPM$mu)
layout <- matrix(
  c("f1",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
    NA,NA,"x1","x2","x3","x7","x8","x9","x13","x14","x15",
    NA,NA,NA,"f9",NA,NA,"f11",NA,NA,"f13",NA,
    NA,NA,NA,"f3",NA,NA,"f5",NA,NA,"f7",NA,
    NA,NA,NA,"f4",NA,NA,"f6",NA,NA,"f8",NA,
    NA,NA,NA,"f10",NA,NA,"f12",NA,NA,"f14",NA,
    NA,NA,"x4","x5","x6","x10","x11","x12","x16","x17","x18",
    "f2",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA
  ), byrow = TRUE, 8, 11)

splot <- semPaths(spm, 'mod', 'par', style = 'ram', layout = layout, 
                  intercepts =  T, residuals = F, fixedStyle = c("#666666", 1), 
                  rotation  = 1, sizeMan = 3, sizeLat = 4)
p_pa3 <- set_curve(splot, c("f1 ~~ f2" = -.2, "f3 ~~ f4" = -1.5, "f5 ~~ f6" = -1.5, "f7 ~~ f8" = -1.5,
                            "x1 ~~ x7" = 1, "x7 ~~ x13" = 1, "x2 ~~ x8" = 1, "x8 ~~ x14" = 1,
                            "x3 ~~ x9" = 1, "x9 ~~ x15" = 1,
                            "x1 ~~ x13" = 1, "x2 ~~ x14" = 1, "x3 ~~ x15" = 1,
                            "x4 ~~ x10" = -1, "x10 ~~ x16" = -1, "x5 ~~ x11" = -1, "x11 ~~ x17" = -1,
                            "x6 ~~ x12" = -1, "x12 ~~ x18" = -1,
                            "x4 ~~ x16" = -1, "x5 ~~ x17" = -1, "x6 ~~ x18" = -1                            
                            ))
p_pa3 <- set_edge_label_position(p_pa3, c('f5 ~ f4' = .25, 'f6 ~ f3' = .25, 'f7 ~ f6' = .25, 'f8 ~ f5' = .25,
                                          'f1 =~ f9' = .25, 'f1 =~ f11' = .225, 'f1 =~ f13' = .2,
                                          'f2 =~ f10' = .25, 'f2 =~ f12' = .225, 'f2 =~ f14' = .2))
p_pa3$Arguments$edge.color[!is.na(p_pa3$Arguments$edge.color)] <- "#666666" # make darker
p_pa3$graphAttributes$Edges$color <- rep("#666666", length(p_pa3$graphAttributes$Edges$color))
p_pa3 <- set_edge_attribute(p_pa3,
                            values = c("f1 ~~ f2" = "red"),
                            attribute_name = "color")
p_pa3 <- change_node_label(p_pa3,
                           c(f1 = "Ix",
                             f2 = "Iy",
                             f9 = "FX1",
                             f11 = "FX2",
                             f13 = "FX3",
                             f10 = "FY1",
                             f12 = "YX2",
                             f14 = "FY3",
                             f3 = "cFX1",
                             f5 = "cFX2",
                             f7 = "cFX3",
                             f4 = "cFY1",
                             f6 = "cFY2",
                             f8 = "cFY3"
                           ),
                           label.cex = 1)


#pdf(file = "docs/fig/riclpm.pdf", width = 7, height = 7)
png(filename = "docs/fig/riclpm.png", width = 14000, height = 14000, units = "px", pointsize = 12, res = 1200)
plot(p_pa3)
dev.off()




##################################

##### LGCM

##################################

powerLGCM <- semPower.powerLGCM(
  type = 'ph', alpha = .05, N = 800,
  nWaves = 3,
  means = c(.5, .2),     # i, s
  variances = c(1, .5),  # i, s
  covariances = .25,
  nullEffect = 'sMean = 0',
  nIndicator = rep(3, 3), loadM = .5
)
spm <- prepParam(powerLGCM$modelH1, powerLGCM$Sigma, mu = powerLGCM$mu)
# layout <- matrix(
#   c(NA,NA,NA,NA,"i",NA,NA,"s",NA,NA,NA,
#     NA,"f1",NA,NA,NA,"f2",NA,NA,NA,"f3",NA,
#     "x1","x2","x3",NA,"x4","x5","x6",NA,"x7","x8","x9"
#   ), byrow = TRUE, 3, 11)

splot <- semPaths(spm, 'mod', 'par', style = 'ram', layout = 'tree', 
                  intercepts =  T, residuals = F, fixedStyle = c("#666666", 1), intAtSide=F,
                  rotation  = 1)
p_pa3 <- move_node(splot,list(f2 = c(x = .12)))
p_pa3 <- set_edge_label_position(p_pa3, c('i =~ f1' = .2, 'i =~ f2' = .2, 'i =~ f3' = .15,
                                          's =~ f1' = .15, 's =~ f2' = .2, 's =~ f3' = .2))
p_pa3$Arguments$edge.color[!is.na(p_pa3$Arguments$edge.color)] <- "#666666" # make darker
p_pa3$graphAttributes$Edges$color <- rep("#666666", length(p_pa3$graphAttributes$Edges$color))
p_pa3$graphAttributes$Edges$color[30] <- 'red'


#pdf(file = "docs/fig/lgcm.pdf", width = 7, height = 7)
png(filename = "docs/fig/lgcm.png", width = 14000, height = 10000, units = "px", pointsize = 12, res = 1200)
plot(p_pa3)
dev.off()




