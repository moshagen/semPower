[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/semPower)](https://www.r-pkg.org/badges/version/semPower)
[![Licence](https://img.shields.io/badge/licence-GPL--2-green.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
[![monthly downloads](https://cranlogs.r-pkg.org:443/badges/semPower)](https://cranlogs.r-pkg.org:443/badges/semPower)
[![total downloads](https://cranlogs.r-pkg.org:443/badges/grand-total/semPower)](https://cranlogs.r-pkg.org:443/badges/grand-total/semPower)

semPower
=====

semPower is an R-package that provides several functions to perform a-priori, post-hoc, and compromise power analyses for structural equation models (SEM). 

### Installation

Install `semPower` via [CRAN](https://CRAN.R-project.org/package=semPower) or as follows:
```
install.packages("devtools")
library("devtools")
install_github("moshagen/semPower")
```

Basic functionality is also provided as a shiny app, which you can use online at [https://sempower.shinyapps.io/sempower](https://sempower.shinyapps.io/sempower/).

### Manual

Read the manual by typing
```
vignette("semPower")
```
or [view the manual online](https://github.com/moshagen/semPower/blob/master/vignettes/semPower.pdf).


### Quick Examples

Determine the required sample size to detect misspecifications of a model (involving df = 100 degrees of freedom) corresponding to RMSEA = .05 with a power of 80% on an alpha error of .05:

```
ap <- semPower.aPriori(effect = .05, effect.measure = 'RMSEA', alpha = .05, power = .80, df = 100)
summary(ap)
```

Determine the achieved power with a sample size of N = 1000 to detect misspecifications of a model (involving df = 100 degrees of freedom) corresponding to RMSEA = .05 on an alpha error of .05:

```
ph <- semPower.postHoc(effect = .05, effect.measure = 'RMSEA', alpha = .05, N = 1000, df = 100)
summary(ph)
```

Determine the critical chi-square such that the associated alpha and beta errors are equal, assuming sample size of N = 1000, a model involving df = 100 degrees of freedom, and misspecifications corresponding to RMSEA = .05:

```
cp <- semPower.compromise(effect = .05, effect.measure = 'RMSEA', abratio = 1, N = 1000, df = 100)
summary(cp)
```

Plot power as function of the sample size to detect misspecifications corresponding to RMSEA = .05 (assuming df = 100) on alpha = .05:

```
semPower.powerPlot.byN(effect = .05, effect.measure = 'RMSEA', alpha = .05, df = 100, power.min = .05, power.max = .99)
```

Plot power as function of the magnitude of effect (measured through the RMSEA assuming df = 100) at N = 500 on alpha = .05:

```
semPower.powerPlot.byEffect(effect.measure = 'RMSEA', alpha = .05, N = 500, df = 100, effect.min = .001, effect.max = .10)
```

For more details and for a description how to express the magnitude of effect in terms of model parameters, see the [manual](https://github.com/moshagen/semPower/blob/master/vignettes/semPower.pdf).



### Citation

If you use `semPower` in publications, please cite the package as follows:

Moshagen, M., & Erdfelder, E. (2016). A new strategy for testing structural equation models.*Structural Equation Modeling, 23*, 54-60. doi: 10.1080/10705511.2014.950896
 
