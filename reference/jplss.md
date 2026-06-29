# Jones-Pewsey location-concentration-shape family

A general family implementing distributional regression for a circular
(angular) response \\y\\ in radians under the Jones-Pewsey law \$\$f(y)
= c(\kappa, \psi)\\\bigl(\cosh(\kappa\psi) + \sinh(\kappa\psi)\cos(y -
\mu)\bigr)^{1/\psi}, \qquad \kappa \> 0,\\ \psi \in \mathbb{R},\$\$ with
\\c(\kappa,\psi)\\ a normalizing constant. This is the symmetric
**peakedness umbrella** that nests several of the package's other
families through the single shape parameter \\\psi\\.

## Usage

``` r
jplss(link = list("tanhalf", "log", "identity"))
```

## Arguments

- link:

  Three-element list of link names, for the mean direction, the
  concentration, and the shape. Currently only the defaults are
  available: `"tanhalf"` for the location \\\mu\\, `"log"` for the
  concentration \\\kappa \> 0\\, and `"identity"` for the shape \\\psi
  \in \mathbb{R}\\.

## Value

An object of class `c("general.family", "extended.family", "family")`
for use with [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (or its
front end
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)).

## Details

It is the first family with **three linear predictors**: the mean
direction \\\mu\\, the concentration \\\kappa\\, and the shape \\\psi\\
each get their own, any of which may contain smooth terms. Used with
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) and a list of
*three* formulas; the first names the response and models \\\mu\\, the
second models \\\log\kappa\\, and the third models \\\psi\\. The shape
is most often held global, i.e. fitted intercept-only with `~ 1`.

The shape \\\psi\\ indexes a family of symmetric circular laws and
recovers several special cases exactly:

- \\\psi \to 0\\: the von Mises
  ([`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md)),
  \\\propto \exp(\kappa\cos(y-\mu))\\;

- \\\psi = 1\\: the cardioid
  ([`cardlss`](https://circstat.github.io/circlss/reference/cardlss.md))
  with mean resultant length \\\rho = \tfrac12\tanh\kappa\\;

- \\\psi = -1\\: the wrapped Cauchy
  ([`wclss`](https://circstat.github.io/circlss/reference/wclss.md));

- \\\kappa \to 0\\: the circular uniform.

The Cartwright family
([`cartlss`](https://circstat.github.io/circlss/reference/cartlss.md))
sits on this family's concentrated boundary (\\\kappa \to \infty\\ at
\\\psi = \zeta\\), not in its interior.

**Normalizer.** Unlike the earlier families, \\c(\kappa, \psi)\\ has no
elementary form. The log-density value *and* the \\\kappa\\- and
\\\psi\\-score and Hessian carry kernel-weighted moments of the
normalizer, evaluated by composite 24-point Gauss-Legendre quadrature on
a feature-scale break-point ladder (so the concentrated \\\psi \< 0\\
spike and the antipodal near-kink are resolved at any concentration).
This is the first circlss family to integrate numerically inside the
likelihood.

**Optimizer.** Only first- and second-order log-likelihood derivatives
are provided, so `available.derivs = 0` and the family is fitted by the
extended Fellner-Schall optimizer rather than full Newton REML.
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) selects
`optimizer = "efs"` automatically; passing it explicitly is recommended
so the fitted object is labelled correctly (which
[`gam.check`](https://rdrr.io/pkg/mgcv/man/gam.check.html) relies on).
The shape \\\psi\\ is weakly identified when the response is diffuse or
its concentration smooth is over-flexible; holding \\\psi\\ global
(`~ 1`) and keeping the concentration well identified is the robust
default.

The mean direction uses the Fisher-Lee tan-half link (\\\mu \in (-\pi,
\pi)\\, antipode unrepresentable, winding number zero – see
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md) when
the mean direction must wind). The shape rides the identity link and may
be any real number.

## References

Jones, M. C. and Pewsey, A. (2005) A family of symmetric distributions
on the circle. *Journal of the American Statistical Association* 100,
1422-1428.

Pewsey, A., Neuhaeuser, M. and Ruxton, G. D. (2013) *Circular Statistics
in R*. Oxford University Press.

Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method
for smoothing parameter optimization with application to Tweedie
location, scale and shape models. *Biometrics* 73, 1071-1081.

## See also

[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
[`cardlss`](https://circstat.github.io/circlss/reference/cardlss.md),
[`wclss`](https://circstat.github.io/circlss/reference/wclss.md),
[`cartlss`](https://circstat.github.io/circlss/reference/cartlss.md),
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md),
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html)

## Examples

``` r
library(mgcv)
set.seed(1)
n <- 300
x <- runif(n)
mu <- 2 * atan(sin(2 * pi * x))
kappa <- exp(1.0 + 0.5 * cos(2 * pi * x))
psi <- 0.6
# Jones-Pewsey draws by inverse transform on the centered kernel, then shift
phi <- seq(-pi, pi, length.out = 4001)
dev <- vapply(seq_len(n), function(i) {
  g <- (cosh(kappa[i] * psi) + sinh(kappa[i] * psi) * cos(phi))^(1 / psi)
  cdf <- cumsum(g); cdf <- cdf / cdf[length(cdf)]
  approx(cdf, phi, runif(1), rule = 2)$y
}, numeric(1))
y <- atan2(sin(mu + dev), cos(mu + dev))
# smooth mean direction and concentration; global (intercept-only) shape
b <- gam(list(y ~ s(x), ~ s(x), ~ 1), family = jplss(), optimizer = "efs")
summary(b)
```
