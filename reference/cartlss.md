# Cartwright power-of-cosine location-scale family

A general family implementing distributional regression for a circular
(angular) response \\y\\ in radians under Cartwright's power-of-cosine
law \$\$f(y) = \frac{2^{1/\zeta - 1}\\\Gamma(1 + 1/\zeta)^2}
{\pi\\\Gamma(1 + 2/\zeta)}\\\bigl(1 + \cos(y - \mu)\bigr)^{1/\zeta},
\qquad \zeta \> 0.\$\$ Both the mean direction \\\mu\\ and the
peakedness \\\zeta\\ get their own linear predictor, each of which may
contain smooth terms. Used with
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) and a list of two
formulas: the first names the response and models \\\mu\\, the second
models \\\log\zeta\\.

## Usage

``` r
cartlss(link = list("tanhalf", "log"))
```

## Arguments

- link:

  Two-element list of link names, for the mean direction and the
  peakedness parameter. Currently only the defaults are available:
  `"tanhalf"` for the location and `"log"` for the shape parameter
  \\\zeta \> 0\\.

## Value

An object of class `c("general.family", "extended.family", "family")`
for use with [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (or its
front end
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)).

## Details

Cartwright's distribution is a one-parameter **peakedness** family:
\\\zeta\\ raises \\1 + \cos(y-\mu)\\ to the power \\1/\zeta\\,
sharpening or flattening the single mode. Its mean resultant length is
\\\rho = 1/(\zeta + 1)\\, so \\\zeta \to 0\\ is sharply peaked (\\\rho
\to 1\\), \\\zeta \to \infty\\ is the circular uniform (\\\rho \to 0\\),
and \\\zeta = 1\\ is exactly the cardioid
([`cardlss`](https://circstat.github.io/circlss/reference/cardlss.md))
at its concentration ceiling, \\(1 + \cos)/2\pi\\ (\\\rho = 1/2\\). It
shares the power-of-cosine form with the Jones-Pewsey family only at the
latter's concentrated boundary; it is not an interior Jones-Pewsey
special case.

The density is evaluated in the half-angle form \\1 + \cos d =
2\cos^2(d/2)\\, which keeps the log-density exact near the antipode \\y
= \mu \pm \pi\\, where it has an *honest zero* for every \\\zeta\\. The
second trigonometric moment is \\\alpha_2 =
(1-\zeta)/\\(1+\zeta)(1+2\zeta)\\\\, so Pearson residuals standardize by
\\\mathrm{Var}\\\sin(y-\mu)\\ = (1 - \alpha_2)/2\\.

**Optimizer.** Only first- and second-order log-likelihood derivatives
are provided, so `available.derivs = 0` and the family is fitted by the
extended Fellner-Schall optimizer rather than full Newton REML.
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) selects
`optimizer = "efs"` automatically; passing it explicitly is recommended
so the fitted object is labelled correctly (which
[`gam.check`](https://rdrr.io/pkg/mgcv/man/gam.check.html) relies on).
The \\\zeta\\-derivatives involve the digamma and trigamma functions,
because the normalizer is built from \\\Gamma(1 + 1/\zeta)\\ and
\\\Gamma(1 + 2/\zeta)\\ – the first family whose normalizing constant is
not elementary.

The mean direction uses the Fisher-Lee tan-half link (\\\mu \in (-\pi,
\pi)\\, antipode unrepresentable, winding number zero – see
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md) when
the mean direction must wind).

## References

Cartwright, D. E. (1963) The use of directional spectra in studying the
output of a wave recorder on a moving ship. In *Ocean Wave Spectra*,
203-218. Prentice-Hall.

Jammalamadaka, S. R. and SenGupta, A. (2001) *Topics in Circular
Statistics*. World Scientific.

Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method
for smoothing parameter optimization with application to Tweedie
location, scale and shape models. *Biometrics* 73, 1071-1081.

## See also

[`cardlss`](https://circstat.github.io/circlss/reference/cardlss.md),
[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
[`wnlss`](https://circstat.github.io/circlss/reference/wnlss.md),
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md),
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html)

## Examples

``` r
library(mgcv)
set.seed(1)
n <- 400
x <- runif(n)
mu <- 2 * atan(sin(2 * pi * x))
zeta <- exp(-0.2 + 0.8 * cos(2 * pi * x))   # peakedness > 0
# Cartwright draws via the Beta-to-angle transform
ang <- 2 * asin(sqrt(rbeta(n, 0.5, 1 / zeta + 0.5)))
th <- mu + sample(c(-1, 1), n, replace = TRUE) * ang
y <- atan2(sin(th), cos(th))
b <- gam(list(y ~ s(x), ~ s(x)), family = cartlss(), optimizer = "efs")
summary(b)
```
