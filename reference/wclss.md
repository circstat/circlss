# Wrapped Cauchy location-scale family

A general family implementing distributional regression for a circular
(angular) response \\y\\ in radians under the wrapped Cauchy law
\$\$f(y) = \frac{1 - \rho^2}{2\pi\\(1 + \rho^2 - 2\rho\cos(y -
\mu))},\$\$ with both the mean direction \\\mu\\ and the mean resultant
length \\\rho\\ getting their own linear predictor, each of which may
contain smooth terms. Used with
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) and a list of two
formulas: the first names the response and models \\\mu\\, the second
models \\\mathrm{logit}(\rho)\\.

## Usage

``` r
wclss(link = list("tanhalf", "logit"))
```

## Arguments

- link:

  Two-element list of link names, for the mean direction and the mean
  resultant length. Currently only the defaults are available:
  `"tanhalf"` for the location and `"logit"` for the concentration
  parameter \\\rho \in (0, 1)\\.

## Value

An object of class `c("general.family", "extended.family", "family")`
for use with [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (or its
front end
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)).

## Details

The wrapped Cauchy is the heavy-tailed counterpart of the von Mises
([`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md)):
sharply peaked with fat circular tails, so it is the more robust choice
when the data contain angular outliers. Its trigonometric moments are
simply \\\rho^p\\, which gives clean residual conventions: Pearson
residuals standardize by \\\mathrm{Var}\\\sin(y-\mu)\\ = (1-\rho^2)/2\\.

The mean direction uses the Fisher-Lee tan-half link (\\\mu \in (-\pi,
\pi)\\, antipode unrepresentable, winding number zero – see
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md) when
the mean direction must wind). Log-likelihood derivatives up to fourth
order are implemented, so the family supports full Newton REML
(`method = "REML"`); `optimizer = "efs"` also works. Internally the
density denominator is computed in the cancellation-free form
\\(1-\rho)^2 + 4\rho\sin^2((y-\mu)/2)\\ so the log-likelihood stays
exact as \\\rho \to 1\\.

## References

Fisher, N. I. and Lee, A. J. (1992) Regression models for an angular
response. *Biometrics* 48, 665-677.

Wood, S. N., Pya, N. and Saefken, B. (2016) Smoothing parameter and
model selection for general smooth models. *Journal of the American
Statistical Association* 111, 1548-1575.

## See also

[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md),
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html)

## Examples

``` r
library(mgcv)
set.seed(1)
n <- 300
x <- runif(n)
mu <- 2 * atan(1.4 * sin(2 * pi * x))
rho <- plogis(0.5 + 0.8 * cos(2 * pi * x))
y <- atan2(sin(mu - log(rho) * rcauchy(n)), cos(mu - log(rho) * rcauchy(n)))
b <- gam(list(y ~ s(x), ~ s(x)), family = wclss(), method = "REML")
summary(b)
```
