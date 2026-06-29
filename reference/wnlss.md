# Wrapped normal location-scale family

A general family implementing distributional regression for a circular
(angular) response \\y\\ in radians under the wrapped normal law
\$\$f(y) = \frac{1}{2\pi}\left(1 + 2\sum\_{p=1}^{\infty}
\rho^{p^2}\cos\\p(y - \mu)\\\right),\$\$ the wrapping of \\N(\mu,
\sigma^2)\\ with \\\rho = e^{-\sigma^2/2}\\. Both the mean direction
\\\mu\\ and the mean resultant length \\\rho\\ get their own linear
predictor, each of which may contain smooth terms. Used with
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) and a list of two
formulas: the first names the response and models \\\mu\\, the second
models \\\mathrm{logit}(\rho)\\.

## Usage

``` r
wnlss(link = list("tanhalf", "logit"))
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

The wrapped normal is the bell-shaped circular law obtained by wrapping
a Gaussian onto the circle – close to the von Mises
([`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md)) in
shape but defined through its mean resultant length \\\rho\\ rather than
a concentration. Its trigonometric moments are \\\rho^{p^2}\\, so
Pearson residuals standardize by \\\mathrm{Var}\\\sin(y-\mu)\\ =
(1-\rho^4)/2\\.

Unlike the von Mises and wrapped Cauchy, the wrapped normal has no
closed-form normalizer. The log-density and its derivatives are
evaluated by a hybrid that switches per observation at \\\rho = 0.8\\:
the Fourier series above for \\\rho \le 0.8\\ (29 terms; truncation
\\\le 2\times10^{-82}\\), and a log-sum-exp over wrapped Gaussian images
for \\\rho \> 0.8\\, where the Fourier partial sums lose accuracy in the
tails.

**Optimizer.** Only first- and second-order log-likelihood derivatives
are provided, so `available.derivs = 0` and the family is fitted by the
extended Fellner-Schall optimizer rather than full Newton REML.
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) selects
`optimizer = "efs"` automatically; passing it explicitly is recommended
so the fitted object is labelled correctly (which
[`gam.check`](https://rdrr.io/pkg/mgcv/man/gam.check.html) relies on).

The mean direction uses the Fisher-Lee tan-half link (\\\mu \in (-\pi,
\pi)\\, antipode unrepresentable, winding number zero – see
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md) when
the mean direction must wind).

## References

Pewsey, A., Neuhaeuser, M. and Ruxton, G. D. (2013) *Circular Statistics
in R*. Oxford University Press.

Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method
for smoothing parameter optimization with application to Tweedie
location, scale and shape models. *Biometrics* 73, 1071-1081.

## See also

[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
[`wclss`](https://circstat.github.io/circlss/reference/wclss.md),
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md),
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html)

## Examples

``` r
library(mgcv)
set.seed(1)
n <- 300
x <- runif(n)
mu <- 2 * atan(1.2 * sin(2 * pi * x))
rho <- plogis(0.8 + 1.0 * cos(2 * pi * x))
sigma <- sqrt(-2 * log(rho))
y <- atan2(sin(mu + sigma * rnorm(n)), cos(mu + sigma * rnorm(n)))
b <- gam(list(y ~ s(x), ~ s(x)), family = wnlss(), optimizer = "efs")
summary(b)
```
