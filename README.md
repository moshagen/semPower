[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/semPower)](https://www.r-pkg.org/badges/version/semPower)
[![Licence](https://img.shields.io/badge/licence-GPL--2-green.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
[![monthly downloads](https://cranlogs.r-pkg.org:443/badges/semPower)](https://cranlogs.r-pkg.org:443/badges/semPower)
[![total downloads](https://cranlogs.r-pkg.org:443/badges/grand-total/semPower)](https://cranlogs.r-pkg.org:443/badges/grand-total/semPower)

semPower
=====
  
semPower is an R-package that provides several functions to perform a-priori, compromise, and post-hoc power analyses for structural equation models (SEM). 

(Very) basic functionality is also provided as a shiny app, which you can use online at [https://sempower.shinyapps.io/sempower](https://sempower.shinyapps.io/sempower/).

### Installation

Install `semPower` via [CRAN](https://CRAN.R-project.org/package=semPower) or as follows:
  ```
install.packages("devtools")
devtools::install_github("moshagen/semPower")
```

### Manual

Read the manual by typing
```
vignette("semPower")
```
or [view the manual online](https://github.com/moshagen/semPower/blob/master/vignettes/semPower.pdf). I also warmly recommend this in-depth tutorial on power analyses in SEM using semPower:
  
  Jobst, L., Bader, M., & Moshagen, M. (in press). A Tutorial on Assessing Statistical Power and Determining Sample Size for Structural Equation Models. *Psychological Methods*. [https://doi.org/10.1037/met0000423](https://doi.org/10.1037/met0000423)  [preprint](https://github.com/moshagen/semPower/blob/master/vignettes/semPowerTutorial.pdf)



### Quick Examples for model-free power analyses

Determine the required sample size to detect misspecifications of a model (involving df = 100 degrees of freedom) corresponding to RMSEA = .05 with a power of 80% on an alpha error of .05:
  
```
ap <- semPower.aPriori(effect = .05, effect.measure = 'RMSEA', 
                       alpha = .05, power = .80, df = 100)
summary(ap)
```

Determine the achieved power with a sample size of N = 1000 to detect misspecifications of a model (involving df = 100 degrees of freedom) corresponding to RMSEA = .05 on an alpha error of .05:
  
```
ph <- semPower.postHoc(effect = .05, effect.measure = 'RMSEA', 
                       alpha = .05, N = 1000, df = 100)
summary(ph)
```

Determine the critical chi-square such that the associated alpha and beta errors are equal, assuming sample size of N = 1000, a model involving df = 100 degrees of freedom, and misspecifications corresponding to RMSEA = .05:
  
```
cp <- semPower.compromise(effect = .05, effect.measure = 'RMSEA', 
                          abratio = 1, N = 1000, df = 100)
summary(cp)
```

Plot power as function of the sample size to detect misspecifications corresponding to RMSEA = .05 (assuming df = 100) on alpha = .05:
  
```
semPower.powerPlot.byN(effect = .05, effect.measure = 'RMSEA', 
                       alpha = .05, df = 100, power.min = .05, power.max = .99)
```

Plot power as function of the magnitude of effect (measured through the RMSEA assuming df = 100) at N = 500 on alpha = .05:
  
```
semPower.powerPlot.byEffect(effect.measure = 'RMSEA', alpha = .05, N = 500, 
                            df = 100, effect.min = .001, effect.max = .10)
```

Obtain the df of a model provided as lavaan model string (this requires the lavaan package):
  
```
lavModel <- '
f1 =~ x1 + x2 + x3
f2 =~ x4 + x5 + x6
'
semPower.getDf(lavModel)
```

Determine the required sample size to discriminate a model exhibiting an RMSEA of .04 on 44 df from a model with RMSEA = .05 on 41 df with a power of 80% on an alpha error of .05:
  
```
ap <- semPower.aPriori(effect = c(.04, .05), effect.measure = 'RMSEA', 
                       alpha = .05, power = .80, df = c(44, 41))
summary(ap)
```

See the [manual](https://github.com/moshagen/semPower/blob/master/vignettes/semPower.pdf) for details.


### Quick Examples for model-based power analyses


All the following examples determine the required sample size to detect the specified effect  (a priori power analysis)
with a power of 80% on alpha .05 based on a factor model involving 4 factors that are measured by 3 indicators each 
with (primary) loading as specified in the `loadings` argument. See the [manual](https://github.com/moshagen/semPower/blob/master/vignettes/semPower.pdf) for details and for other ways to specify the 
measurement model.



#### CFA models

Determine sample size to detect that a correlation between the first and the second factor of at least .3 differs from zero:
```
Phi <- matrix(c(
  c(1, .3, .4, .5),
  c(.3, 1, .2, .6),
  c(.4, .2, 1, .1),
  c(.5, .6, .1, 1)
), ncol = 4, byrow = TRUE)

powerCFA <- semPower.powerCFA(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80,
  # define hypothesis
  Phi = Phi,
  nullEffect = 'cor = 0',
  nullWhich = c(1, 2),
  # define measurement model
  loadings = list(
    c(.7, .6, .5),
    c(.5, .8, .6),
    c(.7, .6, .5),
    c(.5, .8, .6))
  )
summary(powerCFA$power)
```

Determine sample size to detect that the correlations between factor 1 and 2 (of .3) as well as between 3 and 4  (of .1)  differ from each other:
```
Phi <- matrix(c(
  c(1, .3, .4, .5),
  c(.3, 1, .2, .6),
  c(.4, .2, 1, .1),
  c(.5, .6, .1, 1)
), ncol = 4, byrow = TRUE)

powerCFA <- semPower.powerCFA(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80,
  # define hypothesis
  Phi = Phi,
  nullEffect = 'corX = corZ',
  nullWhich = list(c(1, 2), c(3, 4)),
  # define measurement model
  loadings = list(
    c(.7, .6, .5),
    c(.5, .8, .6),
    c(.7, .6, .5),
    c(.5, .8, .6))
)
summary(powerCFA$power)
```

Determine sample size to detect that the correlations between factor 1 and 2 in group (of .3) differs from the one in group 2 (of .5):
```
Phi1 <- matrix(c(
  c(1, .3, .4, .5),
  c(.3, 1, .2, .6),
  c(.4, .2, 1, .1),
  c(.5, .6, .1, 1)
), ncol = 4, byrow = TRUE)
Phi2 <- matrix(c(
  c(1, .5, .4, .5),
  c(.5, 1, .2, .6),
  c(.4, .2, 1, .1),
  c(.5, .6, .1, 1)
), ncol = 4, byrow = TRUE)

powerCFA <- semPower.powerCFA(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80, N = list(1, 1),
  # define hypothesis
  Phi = list(Phi1, Phi2),
  nullEffect = 'corA = corB',
  nullWhich = c(1, 2),
  # define measurement model
  loadings = list(
    c(.7, .6, .5),
    c(.5, .8, .6),
    c(.7, .6, .5),
    c(.5, .8, .6))
)
summary(powerCFA$power)
```


#### Latent regression models
Determine sample size to detect that the first slope (of .2) differs from zero:

```
corXX <- matrix(c(
  c(1, .2, .6),
  c(.2, 1, .1),
  c(.6, .1, 1)
), ncol = 3, byrow = TRUE)

powerReg <- semPower.powerRegression(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80,
  # define hypothesis
  slopes = c(.2, .3, .4), 
  corXX = corXX, 
  nullEffect = 'slope = 0',
  nullWhich = 1,
  # define measurement model
  loadings = list(
    c(.7, .6, .5),
    c(.5, .8, .6),
    c(.7, .6, .5),
    c(.5, .8, .6))
)
summary(powerReg$power)
```

Determine sample size to detect that the first slope (of .2) differs from the third slope (of .4):

```
corXX <- matrix(c(
  c(1, .2, .6),
  c(.2, 1, .1),
  c(.6, .1, 1)
), ncol = 3, byrow = TRUE)

powerReg <- semPower.powerRegression(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80,
  # define hypothesis
  slopes = c(.2, .3, .4), 
  corXX = corXX, 
  nullEffect = 'slopeX = slopeZ',
  nullWhich = c(1, 3),
  # define measurement model
  loadings = list(
    c(.7, .6, .5),
    c(.5, .8, .6),
    c(.7, .6, .5),
    c(.5, .8, .6))
)
summary(powerReg$power)
```

Determine sample size to detect that the first slope in group 1 (of .2) differs from the first slope in group 2 (of .4):

```
corXX <- matrix(c(
  c(1, .2, .6),
  c(.2, 1, .1),
  c(.6, .1, 1)
), ncol = 3, byrow = TRUE)

powerReg <- semPower.powerRegression(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80, N = list(1, 1),
  # define hypothesis
  slopes = list(c(.2, .3, .4),
                c(.4, .3, .2)), 
  corXX = corXX, 
  nullEffect = 'slopeA = slopeB',
  nullWhich = 1,
  # define measurement model
  loadings = list(
    c(.7, .6, .5),
    c(.5, .8, .6),
    c(.7, .6, .5),
    c(.5, .8, .6))
)
summary(powerReg$power)
```

### Mediation models

Determine sample size to detect an indirect effect of at least .12  (= .3*.4) in a simple `X -> M -> Y` mediation based on an observed variable only model:

```
powerMed <- semPower.powerMediation(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80,
  # define hypothesis
  bYX = .25, 
  bMX = .3, 
  bYM = .4,
  nullEffect = 'ind = 0',
  # define observed only 
  Lambda = diag(3)
)
summary(powerMed$power)
```

Determine sample size to detect the indirect effect in group 1 (of .12) differs from the indirect effect in group 2 (of .25) in a simple `X -> M -> Y` mediation based on an observed variable only model:

```
powerMed <- semPower.powerMediation(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80, N = list(1, 1),
  # define hypothesis
  bYX = list(.25, .25), 
  bMX = list(.3, .5), 
  bYM = list(.4, .5),
  nullEffect = 'indA = indB',
  # define observed only 
  Lambda = diag(3)
)
summary(powerMed$power)
```

### Multigroup invariance 

Determine sample size to detect metric-noninvariance across two groups of magnitude as defined through the different loadings:

```
powerMI <- semPower.powerMI(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80, N = list(1, 1),
  # define hypothesis
  comparison = 'configural', 
  nullEffect = 'metric',
  # define measurement model
  loadings = list(
    # group 1
    list(
      c(.7, .6, .5),
      c(.5, .8, .6),
      c(.7, .6, .5),
      c(.5, .8, .6)),
    # group 2
    list(
      c(.6, .5, .4),
      c(.5, .8, .6),
      c(.6, .5, .4),
      c(.5, .8, .6))      
  )
)
summary(powerMI$power)
```

Determine sample size to detect that the latent means differ across groups:

```
powerMI <- semPower.powerMI(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80, N = list(1, 1),
  # define hypothesis
  comparison = c('loadings', 'intercepts'), 
  nullEffect = c('loadings', 'intercepts', 'means'),
  # define measurement model (same for all groups)
  loadings = list(
    c(.7, .6, .5),
    c(.5, .8, .6),
    c(.7, .6, .5),
    c(.5, .8, .6)),
  # define indicator intercepts
  tau = list(rep(0, 12), rep(0, 12)),
  # define latent means
  Alpha = list(
    # group 1
    rep(0, 4), 
    # group 2
    c(0.5, 0, 0.5, 0)
  )
)
summary(powerMI$power)
```

### Cross-lagged panel models (with or without random intercept)

Determine sample size to detect a cross-lagged effect of X on Y of at least .10 in a two-wave CLPM:
  
```
powerCLPM <- semPower.powerCLPM(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80,
  # define hypothesis
  nullEffect = 'crossedX = 0',
  nWaves = 2,
  autoregEffects = c(.60, .70),
  crossedEffects = c(.10, .15),
  rXY = c(.3, .1),
  # define measurement model
  loadings = list(
    c(.7, .6, .5),
    c(.5, .8, .6),
    c(.7, .6, .5),
    c(.5, .8, .6))
)
summary(powerCLPM$power)

```

Same as above, but in a random-intercept CLPM involving 3 waves with observed variables only:
  
```
powerRICLPM <- semPower.powerRICLPM(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80,
  # define hypothesis
  nullEffect = 'crossedX = 0',
  nWaves = 3,
  autoregEffects = c(.60, .70),
  crossedEffects = c(.10, .15),
  rXY = c(.3, .1, .1),
  waveEqual = c('autoregX', 'autoregY', 'crossedX', 'crossedY'),
  # define measurement model
  Lambda = diag(6)
)
summary(powerRICLPM$power)

```


Determine sample size to detect the cross-lagged effect of X on Y differs from the one of Y on X a two-wave CLPM:
  
```
powerCLPM <- semPower.powerCLPM(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80,
  # define hypothesis
  nullEffect = 'crossedX = crossedY',
  nWaves = 2,
  autoregEffects = c(.60, .70),
  crossedEffects = c(.10, .15),
  rXY = c(.3, .1),
  # define measurement model
  loadings = list(
    c(.7, .6, .5),
    c(.5, .8, .6),
    c(.7, .6, .5),
    c(.5, .8, .6))
)
summary(powerCLPM$power)
```


Determine sample size to detect the cross-lagged effect of X on Y in group 1 (of .10) differs from the one in group 2 (of .2):
  
```
powerCLPM <- semPower.powerCLPM(
  # define type of power analysis
  type = 'a-priori', alpha = .05, power = .80, N =list(1, 1),
  # define hypothesis
  nullEffect = 'crossedXA = crossedXB',
  nWaves = 2,
  autoregEffects = c(.60, .70),
  crossedEffects = list(
    # group 1
    list(.10, .15),
    # group 2
    list(.20, .15)
  ),
  rXY = c(.3, .1),
  # define measurement model
  loadings = list(
    c(.7, .6, .5),
    c(.5, .8, .6),
    c(.7, .6, .5),
    c(.5, .8, .6))
)
summary(powerCLPM$power)
```


For more details, see the [manual](https://github.com/moshagen/semPower/blob/master/vignettes/semPower.pdf).



### Citation

If you use `semPower` in publications, please cite the package as follows:
  
  Moshagen, M., & Erdfelder, E. (2016). A new strategy for testing structural equation models. *Structural Equation Modeling, 23*, 54-60. doi: 10.1080/10705511.2014.950896
