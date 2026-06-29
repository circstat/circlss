# Cardioid location-scale family

A general family implementing distributional regression for a circular
(angular) response \\y\\ in radians under the cardioid law \$\$f(y) =
\frac{1}{2\pi}\left(1 + 2\rho\cos(y - \mu)\right), \qquad 0 \le \rho \le
1/2,\$\$ a first-harmonic perturbation of the circular uniform. Both the
mean direction \\\mu\\ and the mean resultant length \\\rho\\ get their
own linear predictor, each of which may contain smooth terms. Used with
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) and a list of two
formulas: the first names the response and models \\\mu\\, the second
models \\\mathrm{logithalf}(\rho)\\.

## Usage

``` r
cardlss(link = list("tanhalf", "logithalf"))
```

## Arguments

- link:

  Two-element list of link names, for the mean direction and the mean
  resultant length. Currently only the defaults are available:
  `"tanhalf"` for the location and `"logithalf"` for the concentration
  parameter \\\rho \in (0, 1/2)\\.

## Value

An object of class `c("general.family", "extended.family", "family")`
for use with [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (or its
front end
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)).

## Details

The cardioid is the simplest departure from circular uniformity: a
single cosine ripple of amplitude \\2\rho\\ on the flat density. It is
therefore a **low-concentration / near-uniform** law – \\\rho = 0\\ is
exactly uniform, and concentration is capped at \\\rho = 1/2\\, beyond
which the density would go negative at the antimode. That hard upper
bound is why the concentration uses the **logit-half** link \\\eta =
\log\\\rho / (1/2 - \rho)\\\\, i.e. \\\rho =
\tfrac{1}{2}\\\mathrm{plogis}(\eta) \in (0, 1/2)\\ – the one new link
this family brings (the von Mises, wrapped Cauchy and wrapped normal all
use the ordinary logit on \\(0, 1)\\).

The cardioid is a pure first-harmonic distribution: its second and
higher trigonometric moments vanish, so the first moment is \\\rho\\
(hence the moment estimator \\\hat\rho = \bar R\\) and Pearson residuals
standardize by the constant \\\mathrm{Var}\\\sin(y-\mu)\\ = 1/2\\.

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

Jammalamadaka, S. R. and SenGupta, A. (2001) *Topics in Circular
Statistics*. World Scientific.

Pewsey, A., Neuhaeuser, M. and Ruxton, G. D. (2013) *Circular Statistics
in R*. Oxford University Press.

Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method
for smoothing parameter optimization with application to Tweedie
location, scale and shape models. *Biometrics* 73, 1071-1081.

## See also

[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
[`wclss`](https://circstat.github.io/circlss/reference/wclss.md),
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
rho <- 0.5 * plogis(0.3 + 1.2 * cos(2 * pi * x))   # mean resultant length < 1/2
# cardioid draws by rejection from a uniform envelope
phi <- runif(n, -pi, pi)
keep <- runif(n) <= (1 + 2 * rho * cos(phi)) / (1 + 2 * rho)
while (any(!keep)) {
  j <- which(!keep)
  phi[j] <- runif(length(j), -pi, pi)
  keep[j] <- runif(length(j)) <= (1 + 2 * rho[j] * cos(phi[j])) / (1 + 2 * rho[j])
}
y <- atan2(sin(mu + phi), cos(mu + phi))
b <- gam(list(y ~ s(x), ~ s(x)), family = cardlss(), optimizer = "efs")
summary(b)
```
