# Plot a circular-response GAM fit

The default plot method for
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
fits: geometry-aware curves that
[`plot.gam`](https://rdrr.io/pkg/mgcv/man/plot.gam.html) does not
provide, covering the three regression geometries – circular–linear
(cylinder), circular–circular (torus) and linear–circular (the upright
can). For a fit with several covariates, no covariate, or an ordinary
(non-circlss) family it falls back to mgcv's per-smooth term plots.

## Usage

``` r
# S3 method for class 'circ_gam'
plot(
  x,
  view = c("flat", "geometry", "both"),
  n = 200,
  se = TRUE,
  pages = 1,
  rug = TRUE,
  ...
)
```

## Arguments

- x:

  A fitted
  [`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
  model.

- view:

  `"flat"` (default) draws one panel per distribution parameter on the
  response scale against the covariate. A circular location or direction
  carries a band of \\\pm\\ one circular standard deviation
  \\\sqrt{-2\log R}\\ of the fitted law (its predictive angular spread,
  with \\R\\ the mean resultant length), wrapped around the \\\pm\pi\\
  branch cut; a non-circular parameter (the concentration, ...) keeps a
  delta-method 2-SE band. The location panel is broken at the \\\pm\pi\\
  branch jump and the observed responses are overlaid. `"geometry"`
  draws the fitted curve on its natural surface – a **cylinder** for
  circular–linear (the response angle wraps the tube, the covariate runs
  along the axis), a **torus** for circular–circular (covariate around
  the ring), or an upright **can** for linear–circular (a linear
  response over a cyclic covariate: the covariate wraps the ring, the
  response is the height) – chosen from whether the response is circular
  and whether the covariate is cyclic. `"both"` places the geometry
  canvas beside the full set of flat parameter panels – exactly the
  panels `"flat"` draws, so the geometry and flat views never disagree.

- n:

  Number of grid points along the covariate.

- se:

  Draw the uncertainty band: \\\pm\\ the circular standard deviation for
  a circular location/direction, a 2-SE interval for a non-circular
  parameter.

- pages:

  If 1, lay the flat panels out on a single page.

- rug:

  Add a covariate rug to the location panel.

- ...:

  Passed to [`plot.gam`](https://rdrr.io/pkg/mgcv/man/plot.gam.html) in
  the fallback cases.

## Value

The fitted model, invisibly.

## Details

The geometry canvas is base-graphics only
([`persp`](https://rdrr.io/r/graphics/persp.html) +
[`trans3d`](https://rdrr.io/r/grDevices/trans3d.html)) and is chosen
descriptively from the covariate's basis; it is never a fitting input.
For [`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md) the
curve drawn is the derived mean direction \\\mathrm{atan2}(\mu_2,
\mu_1)\\.

## See also

[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md),
[`plot.gam`](https://rdrr.io/pkg/mgcv/man/plot.gam.html)
