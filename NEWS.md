# circlss 0.1.0

First release.

## Circular families for mgcv

Twelve circular-response distributional families written against mgcv's
general-family interface, so `gam()` fits them directly — each distribution
parameter on its own (optionally smooth) predictor, for circular–linear and
circular–circular regression:

* von Mises (`vmlss()`), projected normal (`pnlss()`), wrapped Cauchy
  (`wclss()`), wrapped normal (`wnlss()`), cardioid (`cardlss()`),
  Cartwright (`cartlss()`)
* Jones–Pewsey (`jplss()`) with its sine-skewed (`ssjplss()`) and asymmetric
  (`ajplss()`) forms, Kato–Jones (`kjlss()`), flat-topped von Mises
  (`vmftlss()`), and inverse Batschelet (`ibslss()`)

Weight-aware Gaussian (`gausslss()`) and gamma (`gammalss()`) location-scale
families cover a linear response over a circular covariate.

## Front-end and helpers

* `circ_gam()` — thin `gam()` wrapper for these families: auto-pins
  cyclic-smooth period knots, fills unmodeled parameters with `~ 1`, defaults
  to `method = "REML"`, and adds geometry-aware `print` / `plot`.
* `circ_mix()` — finite mixtures of circlss families fitted by EM.
* `circ_lm()` — circular-response linear models.
* Helpers: `rad()` / `deg()` / `wrap()` for radian and branch handling;
  `circ_resid()` / `circ_check()` for diagnostics.
