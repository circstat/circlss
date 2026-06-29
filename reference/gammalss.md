# Weight-aware gamma location-scale family

A gamma location-scale family for distributional regression of a
*positive*, right-skewed response, with a (log) mean and a log scale
each modelled by its own linear predictor. It is a weight-aware,
metadata-carrying adaptation of mgcv's
[`gammals`](https://rdrr.io/pkg/mgcv/man/gammals.html): unlike `gammals`
it honours prior `weights` (needed for a weighted MLE and for EM
mixtures), and it carries the circlss parameter metadata so
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
treats it as a first-class location-scale family.

## Usage

``` r
gammalss(link = list("identity", "log"), b = -7)
```

## Arguments

- link:

  Two-element list of link names, following
  [`gammals`](https://rdrr.io/pkg/mgcv/man/gammals.html): `"identity"`
  for the (log) mean and `"log"` for the log scale.

- b:

  The log-scale link's offset, as in
  [`gammals`](https://rdrr.io/pkg/mgcv/man/gammals.html).

## Value

An object of class `c("general.family", "extended.family", "family")`
for use with [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (or its
front end
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)).

## Details

This is the positive-response member of the circlss linear-circular
(l~c) leg: a positive, skewed quantity that varies around a cycle –
rainfall by season, a concentration or rate by time of day, a speed by
direction – over a circular covariate fitted with a cyclic smooth. As
for
[`gausslss`](https://circstat.github.io/circlss/reference/gausslss.md),
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
places the fitted mean on the "can".

The parameterization follows `gammals`: a (log) mean and a log scale. At
unit weights the fit matches `gammals`; integer prior weights reproduce
a row-replicated fit. Log-likelihood derivatives up to fourth order are
implemented, so the family supports full Newton REML
(`method = "REML"`); `optimizer = "efs"` also works.

This family adapts GPL-licensed code from mgcv; see the package's
`inst/COPYRIGHTS`.

## See also

[`gausslss`](https://circstat.github.io/circlss/reference/gausslss.md),
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md),
[`gammals`](https://rdrr.io/pkg/mgcv/man/gammals.html)

## Examples

``` r
library(mgcv)
set.seed(1); n <- 300
phi <- runif(n, -pi, pi)                       # circular covariate (radians)
y <- rgamma(n, shape = 4, rate = 4 / exp(0.4 + 0.8 * sin(phi)))
b <- circ_gam(list(y ~ s(phi, bs = "cc"), ~ s(phi, bs = "cc")),
              data = data.frame(y, phi), family = gammalss())
head(predict(b, type = "response"))            # columns named mu, scale
```
