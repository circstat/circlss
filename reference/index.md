# Package index

## Model fitting

Fit circular-response GAM, the classical (parametric) circular
regressions

- [`circ_gam()`](https://circstat.github.io/circlss/reference/circ_gam.md)
  [`predict(`*`<circ_gam>`*`)`](https://circstat.github.io/circlss/reference/circ_gam.md)
  [`fitted(`*`<circ_gam>`*`)`](https://circstat.github.io/circlss/reference/circ_gam.md)
  [`print(`*`<circ_gam>`*`)`](https://circstat.github.io/circlss/reference/circ_gam.md)
  : Circular-response GAM
- [`circ_lm()`](https://circstat.github.io/circlss/reference/circ_lm.md)
  [`coef(`*`<circ_lm>`*`)`](https://circstat.github.io/circlss/reference/circ_lm.md)
  [`fitted(`*`<circ_lm>`*`)`](https://circstat.github.io/circlss/reference/circ_lm.md)
  [`residuals(`*`<circ_lm>`*`)`](https://circstat.github.io/circlss/reference/circ_lm.md)
  [`logLik(`*`<circ_lm>`*`)`](https://circstat.github.io/circlss/reference/circ_lm.md)
  [`predict(`*`<circ_lm>`*`)`](https://circstat.github.io/circlss/reference/circ_lm.md)
  [`print(`*`<circ_lm>`*`)`](https://circstat.github.io/circlss/reference/circ_lm.md)
  : Classical circular regression
- [`plot(`*`<circ_gam>`*`)`](https://circstat.github.io/circlss/reference/plot.circ_gam.md)
  : Plot a circular-response GAM fit
- [`plot(`*`<circ_lm>`*`)`](https://circstat.github.io/circlss/reference/plot.circ_lm.md)
  : Plot a classical circular regression fit

## Diagnostics

Circular residuals and the diagnostic panel grid for circ_gam() and
circ_lm() fits – the circular analogue of plot.lm() and gam.check(),
complementing the effect-display plot methods.

- [`circ_resid()`](https://circstat.github.io/circlss/reference/circ_resid.md)
  : Circular regression residuals
- [`circ_check()`](https://circstat.github.io/circlss/reference/circ_check.md)
  : Diagnostic panels for a circular regression fit

## Mixture models

An EM engine for finite mixtures of circular distributional GAMs –
density clustering and the circular-linear / circular-circular /
linear-circular regression trio.

- [`circ_mix()`](https://circstat.github.io/circlss/reference/circ_mix.md)
  [`circ_mix.control()`](https://circstat.github.io/circlss/reference/circ_mix.md)
  [`print(`*`<circ_mix>`*`)`](https://circstat.github.io/circlss/reference/circ_mix.md)
  [`summary(`*`<circ_mix>`*`)`](https://circstat.github.io/circlss/reference/circ_mix.md)
  [`logLik(`*`<circ_mix>`*`)`](https://circstat.github.io/circlss/reference/circ_mix.md)
  [`coef(`*`<circ_mix>`*`)`](https://circstat.github.io/circlss/reference/circ_mix.md)
  [`predict(`*`<circ_mix>`*`)`](https://circstat.github.io/circlss/reference/circ_mix.md)
  : Finite mixtures of circular distributional GAMs by EM
- [`plot(`*`<circ_mix>`*`)`](https://circstat.github.io/circlss/reference/plot.circ_mix.md)
  : Plot a fitted circular mixture
- [`circ_kmeans()`](https://circstat.github.io/circlss/reference/circ_kmeans.md)
  : Circular k-means clustering

## Response families

mgcv general families for angular responses: one linear predictor — and
so one formula, with smooths — per distribution parameter.

- [`vmlss()`](https://circstat.github.io/circlss/reference/vmlss.md) :
  von Mises location-scale family
- [`pnlss()`](https://circstat.github.io/circlss/reference/pnlss.md) :
  Projected normal location family
- [`wclss()`](https://circstat.github.io/circlss/reference/wclss.md) :
  Wrapped Cauchy location-scale family
- [`wnlss()`](https://circstat.github.io/circlss/reference/wnlss.md) :
  Wrapped normal location-scale family
- [`cardlss()`](https://circstat.github.io/circlss/reference/cardlss.md)
  : Cardioid location-scale family
- [`cartlss()`](https://circstat.github.io/circlss/reference/cartlss.md)
  : Cartwright power-of-cosine location-scale family
- [`jplss()`](https://circstat.github.io/circlss/reference/jplss.md) :
  Jones-Pewsey location-concentration-shape family
- [`ssjplss()`](https://circstat.github.io/circlss/reference/ssjplss.md)
  : Sine-skewed Jones-Pewsey location-concentration-shape-skewness
  family
- [`kjlss()`](https://circstat.github.io/circlss/reference/kjlss.md) :
  Kato-Jones Mobius four-parameter family
- [`vmftlss()`](https://circstat.github.io/circlss/reference/vmftlss.md)
  : Flat-topped von Mises location-concentration-shape family
- [`ibslss()`](https://circstat.github.io/circlss/reference/ibslss.md) :
  Inverse Batschelet location-concentration-skewness-peakedness family
- [`ajplss()`](https://circstat.github.io/circlss/reference/ajplss.md) :
  Asymmetric Jones-Pewsey location-concentration-shape-asymmetry family

## Linear-response families

Weight-aware location-scale families for a linear response with circular
covariates (l~c), drop-in compatible with circ_gam().

- [`gausslss()`](https://circstat.github.io/circlss/reference/gausslss.md)
  : Weight-aware Gaussian location-scale family
- [`gammalss()`](https://circstat.github.io/circlss/reference/gammalss.md)
  : Weight-aware gamma location-scale family

## Utilities

Angle conversion and wrapping helpers (radians/degrees, branch
wrapping).

- [`rad()`](https://circstat.github.io/circlss/reference/rad.md)
  [`deg()`](https://circstat.github.io/circlss/reference/rad.md)
  [`wrap()`](https://circstat.github.io/circlss/reference/rad.md) :
  Radian conversion helpers
