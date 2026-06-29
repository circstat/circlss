# Changelog

## circlss 0.1.0

First release.

### Circular families for mgcv

Twelve circular-response distributional families written against mgcv’s
general-family interface, so
[`gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) fits them directly —
each distribution parameter on its own (optionally smooth) predictor,
for circular–linear and circular–circular regression:

- von Mises
  ([`vmlss()`](https://circstat.github.io/circlss/reference/vmlss.md)),
  projected normal
  ([`pnlss()`](https://circstat.github.io/circlss/reference/pnlss.md)),
  wrapped Cauchy
  ([`wclss()`](https://circstat.github.io/circlss/reference/wclss.md)),
  wrapped normal
  ([`wnlss()`](https://circstat.github.io/circlss/reference/wnlss.md)),
  cardioid
  ([`cardlss()`](https://circstat.github.io/circlss/reference/cardlss.md)),
  Cartwright
  ([`cartlss()`](https://circstat.github.io/circlss/reference/cartlss.md))
- Jones–Pewsey
  ([`jplss()`](https://circstat.github.io/circlss/reference/jplss.md))
  with its sine-skewed
  ([`ssjplss()`](https://circstat.github.io/circlss/reference/ssjplss.md))
  and asymmetric
  ([`ajplss()`](https://circstat.github.io/circlss/reference/ajplss.md))
  forms, Kato–Jones
  ([`kjlss()`](https://circstat.github.io/circlss/reference/kjlss.md)),
  flat-topped von Mises
  ([`vmftlss()`](https://circstat.github.io/circlss/reference/vmftlss.md)),
  and inverse Batschelet
  ([`ibslss()`](https://circstat.github.io/circlss/reference/ibslss.md))

Weight-aware Gaussian
([`gausslss()`](https://circstat.github.io/circlss/reference/gausslss.md))
and gamma
([`gammalss()`](https://circstat.github.io/circlss/reference/gammalss.md))
location-scale families cover a linear response over a circular
covariate.

### Front-end and helpers

- [`circ_gam()`](https://circstat.github.io/circlss/reference/circ_gam.md)
  — thin [`gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) wrapper for
  these families: auto-pins cyclic-smooth period knots, fills unmodeled
  parameters with `~ 1`, defaults to `method = "REML"`, and adds
  geometry-aware `print` / `plot`.
- [`circ_mix()`](https://circstat.github.io/circlss/reference/circ_mix.md)
  — finite mixtures of circlss families fitted by EM.
- [`circ_lm()`](https://circstat.github.io/circlss/reference/circ_lm.md)
  — circular-response linear models.
- Helpers:
  [`rad()`](https://circstat.github.io/circlss/reference/rad.md) /
  [`deg()`](https://circstat.github.io/circlss/reference/rad.md) /
  [`wrap()`](https://circstat.github.io/circlss/reference/rad.md) for
  radian and branch handling;
  [`circ_resid()`](https://circstat.github.io/circlss/reference/circ_resid.md)
  /
  [`circ_check()`](https://circstat.github.io/circlss/reference/circ_check.md)
  for diagnostics.
