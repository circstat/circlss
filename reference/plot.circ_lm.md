# Plot a classical circular regression fit

The default plot method for
[`circ_lm`](https://circstat.github.io/circlss/reference/circ_lm.md)
fits: the closed-form sibling of
[`plot.circ_gam`](https://circstat.github.io/circlss/reference/plot.circ_gam.md),
drawing the same three regression geometries – circular–linear
(cylinder), circular–circular (torus) and linear–circular (the upright
can) – from the fit's `type`. A multi-covariate `"cl"` fit has no single
covariate axis and defers with a message to
[`coef`](https://rdrr.io/r/stats/coef.html) /
[`predict`](https://rdrr.io/r/stats/predict.html) / `summary`.

## Usage

``` r
# S3 method for class 'circ_lm'
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
  [`circ_lm`](https://circstat.github.io/circlss/reference/circ_lm.md)
  model.

- view:

  `"flat"` (default) draws one panel per modelled parameter on the
  response scale against the covariate. A circular location carries a
  band of \\\pm\\ one circular standard deviation of the fitted law (the
  von Mises spread \\\sqrt{-2\log A_1(\kappa)}\\ of the fitted
  concentration); the concentration panel keeps a delta-method 2-SE
  band. A circular location is broken at the \\\pm\pi\\ branch jump and
  the observed responses are overlaid. `type = "cl"` draws the mean
  direction and the concentration \\\kappa\\; `type = "cc"` / `"lc"`
  draw the single fitted location. `"geometry"` draws the fitted
  location curve on its natural surface – a **cylinder** for
  circular–linear (`"cl"`: the response angle wraps the tube, the linear
  covariate runs along the axis), a **torus** for circular–circular
  (`"cc"`: the covariate around the ring), or an upright **can** for
  linear–circular (`"lc"`: the cyclic covariate wraps the ring, the
  linear response is the height). `"both"` places the geometry canvas
  beside the full set of flat panels (`"cl"`: the mean direction and
  \\\kappa\\; `"cc"`/`"lc"`: the single location) – exactly the panels
  `"flat"` draws, so the two views never disagree.

- n:

  Number of grid points along the covariate.

- se:

  Draw the uncertainty band (filled shadow on the flat panels,
  translucent ribbon on the geometry surface): \\\pm\\ the circular
  standard deviation for a circular location, a 2-SE interval for the
  concentration.

- pages:

  If 1, lay the flat panels out on a single page.

- rug:

  Add a covariate rug to the location panel.

- ...:

  Currently ignored.

## Value

The fitted model, invisibly.

## Details

The geometry canvas is base-graphics only
([`persp`](https://rdrr.io/r/graphics/persp.html) +
[`trans3d`](https://rdrr.io/r/grDevices/trans3d.html)) and shares the
surface, panel and band helpers with
[`plot.circ_gam`](https://circstat.github.io/circlss/reference/plot.circ_gam.md),
so a `circ_lm` leg renders in the same idiom as the matching `circ_gam`
leg. A circular location is banded by \\\pm\\ the circular standard
deviation of the fitted law – the von Mises spread \\\sqrt{-2\log
A_1(\kappa)}\\ for the per-point concentration of `"cl"` and the
residual concentration of `"cc"` – the predictive angular spread, not a
confidence interval of the mean. The remaining bands are the usual
intervals: the `"cl"` concentration on the log scale through \\Z\\V
Z'\\, the `"lc"` linear response the ordinary least-squares prediction
band.

## See also

[`circ_lm`](https://circstat.github.io/circlss/reference/circ_lm.md),
[`plot.circ_gam`](https://circstat.github.io/circlss/reference/plot.circ_gam.md)

## Examples

``` r
set.seed(1)
n <- 80
x <- rnorm(n)
theta <- (1 + 2 * atan(1.5 * x) + rnorm(n) / 4) %% (2 * pi)
dat <- data.frame(theta = theta, x = x)

## cl: the mean direction on the cylinder, with the kappa panel beside it
m <- circ_lm(theta ~ x, dat, type = "cl")
plot(m)                      # flat: mu (circular) and kappa
plot(m, view = "geometry")   # the fitted angle on the cylinder

## cc / lc: a single location on the torus / can
phi <- runif(n, 0, 2 * pi)
dcc <- data.frame(psi = (phi / 2 + rnorm(n) / 5) %% (2 * pi), phi = phi)
plot(circ_lm(psi ~ phi, dcc, type = "cc"), view = "both")

dlc <- data.frame(y = 5 + 2 * cos(phi) + rnorm(n) / 2, phi = phi)
plot(circ_lm(y ~ phi, dlc, type = "lc"), view = "geometry")
```
