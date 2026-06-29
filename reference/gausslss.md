# Weight-aware Gaussian location-scale family

A Gaussian location-scale family for distributional regression,
modelling a real-valued response \\y\\ with \$\$y \sim N(\mu,
\sigma^2),\$\$ where the mean \\\mu\\ and the precision \\\tau =
1/\sigma\\ each get their own linear predictor, which may contain smooth
terms. It is a weight-aware, metadata-carrying adaptation of mgcv's
[`gaulss`](https://rdrr.io/pkg/mgcv/man/gaulss.html): unlike `gaulss` it
honours prior `weights` (needed for a weighted MLE and for EM mixtures),
and it carries the circlss parameter metadata so
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
treats it as a first-class location-scale family with named,
response-scale output.

## Usage

``` r
gausslss(link = list("identity", "logb"), b = 0.01)
```

## Arguments

- link:

  Two-element list of link names for the mean and the precision,
  following [`gaulss`](https://rdrr.io/pkg/mgcv/man/gaulss.html):
  `"identity"` (or `"log"`, `"inverse"`, `"sqrt"`) for the mean and
  `"logb"` for the precision.

- b:

  The `logb` link's offset, as in
  [`gaulss`](https://rdrr.io/pkg/mgcv/man/gaulss.html).

## Value

An object of class `c("general.family", "extended.family", "family")`
for use with [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (or its
front end
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)).

## Details

In the circlss regression trio this is the **linear-circular** (l~c)
family: a real-valued response over a *circular* covariate – a level
that varies around a cycle (time of day, season, phase) – fitted with a
cyclic smooth,

    circ_gam(list(y ~ s(phi, bs = "cc"), ~ s(phi, bs = "cc")),
             family = gausslss())

[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
then places the fitted mean on the "can" (an upright cylinder: the
circular covariate wraps the ring, the linear response is the height).

The parameterization follows `gaulss` exactly: the mean uses an identity
link and the second parameter is the *precision* \\\tau = 1/\sigma\\ on
the `logb` link, so the second fitted column is \\1/\sigma\\, not the
standard deviation. Log-likelihood derivatives up to fourth order are
implemented, so the family supports full Newton REML
(`method = "REML"`); `optimizer = "efs"` also works. At unit weights the
fit matches `gaulss`; integer prior weights reproduce a row-replicated
fit.

This family adapts GPL-licensed code from mgcv; see the package's
`inst/COPYRIGHTS`.

## See also

[`gammalss`](https://circstat.github.io/circlss/reference/gammalss.md),
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md),
[`gaulss`](https://rdrr.io/pkg/mgcv/man/gaulss.html)

## Examples

``` r
library(mgcv)
set.seed(1); n <- 300
phi <- runif(n, -pi, pi)                       # circular covariate (radians)
y <- 2 + 1.5 * sin(phi) + 0.8 * cos(2 * phi) + rnorm(n) * 0.3
b <- circ_gam(list(y ~ s(phi, bs = "cc"), ~ s(phi, bs = "cc")),
              data = data.frame(y, phi), family = gausslss())
head(predict(b, type = "response"))            # columns named mu, tau
plot(b, view = "both")                         # the l~c "can" + the mean panel
```
