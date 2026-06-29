# Plot a fitted circular mixture

Clustered views of a fitted mixture. The flat view colours the
observations by their MAP cluster, with – for a regression cell – each
component's fitted location curve over the single covariate, or – for a
density cell (`theta ~ 1`) – the per-cluster spread of the response with
each component's fitted mean direction. The geometry view draws those
per-component curves on the leg's natural 3D surface (cylinder / torus /
can).

## Usage

``` r
# S3 method for class 'circ_mix'
plot(x, view = c("flat", "geometry", "both"), n = 200, se = TRUE, ...)
```

## Arguments

- x:

  A fitted
  [`circ_mix`](https://circstat.github.io/circlss/reference/circ_mix.md)
  object.

- view:

  Which view to draw. `"flat"` (default) draws the flat panel.
  `"geometry"` draws every component's fitted location curve on the
  leg's natural surface – the cylinder (c~l), torus (c~c) or can (l~c),
  sharing
  [`plot.circ_gam`](https://circstat.github.io/circlss/reference/plot.circ_gam.md)'s
  surface canvas; `"both"` places the two side by side. A joint
  (product) component has only the flat torus-square projection, so it
  draws that whatever the `view`.

- n:

  Number of grid points for each component's fitted curve.

- se:

  If `TRUE`, band each component's fitted curve where the view supports
  it: \\\pm\\ the component law's circular standard deviation (its
  predictive angular spread) for a circular response, a pointwise 2-SE
  interval for a linear response.

- ...:

  Further graphical arguments (currently unused).

## Value

`x`, invisibly. Called for the plot it draws.

## See also

[`circ_mix`](https://circstat.github.io/circlss/reference/circ_mix.md),
[`plot.circ_gam`](https://circstat.github.io/circlss/reference/plot.circ_gam.md)

## Examples

``` r
library(mgcv)
set.seed(1); n <- 400
z  <- sample.int(2L, n, replace = TRUE)
x  <- runif(n, -1, 1)
mu <- 2 * atan(c(1, -1)[z] + c(2, -2)[z] * x)
y  <- vmlss()$rd(cbind(mu, rep(6, n)), rep(1, n), 1)
# \donttest{
m  <- circ_mix(y ~ x, data = data.frame(y, x), family = vmlss(), K = 2)
plot(m)
# }
```
