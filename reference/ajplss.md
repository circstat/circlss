# Asymmetric Jones-Pewsey location-concentration-shape-asymmetry family

A general family implementing distributional regression for a circular
(angular) response \\y\\ in radians under the asymmetric Jones-Pewsey
law \$\$f(y) = c(\kappa, \psi, \nu)\\\bigl(\cosh(\kappa\psi) +
\sinh(\kappa\psi)\cos g\bigr)^{1/\psi}, \qquad g = \phi + \nu\cos\phi,\\
\\ \phi = y - \xi,\$\$ the Jones-Pewsey family
([`jplss`](https://circstat.github.io/circlss/reference/jplss.md)) with
its angle warped *forward* by \\g(\phi) = \phi + \nu\cos\phi\\ (Abe,
Pewsey & Shimizu 2013). It adds an asymmetry axis to the symmetric
Jones-Pewsey umbrella: \\\nu \> 0\\ and \\\nu \< 0\\ tilt the density to
opposite sides. (For \\\psi \to 0\\ the Jones-Pewsey kernel is the von
Mises \\\exp(\kappa\cos g)\\.)

## Usage

``` r
ajplss(link = list("tanhalf", "log", "identity", "tanh"))
```

## Arguments

- link:

  Four-element list of link names, for the location, concentration,
  shape and asymmetry. Currently only the defaults are available:
  `"tanhalf"` for the location \\\xi\\, `"log"` for the concentration
  \\\kappa \> 0\\, `"identity"` for the shape \\\psi \in \mathbb{R}\\,
  and `"tanh"` for the asymmetry \\\nu \in (-1, 1)\\.

## Value

An object of class `c("general.family", "extended.family", "family")`
for use with [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (or its
front end
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)).

## Details

It has **four linear predictors**: the location \\\xi\\, the
concentration \\\kappa\\, the shape \\\psi\\ and the asymmetry \\\nu\\
each get their own, any of which may contain smooth terms. Used with
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) and a list of
*four* formulas; the first names the response and models \\\xi\\, then
\\\log\kappa\\, then \\\psi\\, then \\\nu\\. The shape and asymmetry are
most often held global, i.e. fitted intercept-only with `~ 1`.

The shape \\\psi\\ and asymmetry \\\nu\\ index a family of symmetric and
skewed circular laws:

- \\\nu = 0\\: the symmetric Jones-Pewsey
  ([`jplss`](https://circstat.github.io/circlss/reference/jplss.md)),
  whose \\\psi\\ in turn nests the von Mises (\\\psi \to 0\\), the
  cardioid (\\\psi = 1\\) and the wrapped Cauchy (\\\psi \to -\infty\\);

- \\\nu \ne 0\\: an asymmetric (skewed) law;

- \\\kappa \to 0\\: the circular uniform.

**Mode anchor, not mean.** As with
[`ssjplss`](https://circstat.github.io/circlss/reference/ssjplss.md) and
[`ibslss`](https://circstat.github.io/circlss/reference/ibslss.md), once
\\\nu \ne 0\\ the location \\\xi\\ is the *mode anchor* of the
asymmetric density, not its mean direction; fitted-direction summaries
inherit that reading.

**Forward warp – no implicit differentiation.** The warp \\g(\phi) =
\phi + \nu\cos\phi\\ is a monotone reparameterization of the angle
(\\g'(\phi) = 1 - \nu\sin\phi \> 0\\ for \\\|\nu\| \< 1\\), so – unlike
[`ibslss`](https://circstat.github.io/circlss/reference/ibslss.md)'s
inverse warps – the score needs *no* implicit differentiation: it reuses
the Jones-Pewsey kernel terms evaluated at the warped angle \\g\\ and
chains them through \\g\\ (with \\g\_\phi = 1 - \nu\sin\phi\\ and
\\g\_\nu = \cos\phi\\).

**Normalizer and cross terms.** Unlike
[`ssjplss`](https://circstat.github.io/circlss/reference/ssjplss.md)'s
sine-skew factor – which integrates to 1 against the symmetric kernel
and so leaves the Jones-Pewsey normalizer *unchanged* – the forward warp
*moves* the normalizer: \\c(\kappa, \psi, \nu)\\ is the \\u = g(\phi)\\
substitution integrating the kernel against the warp Jacobian \\1/g'\\
over an asymmetry-aware Gauss-Legendre ladder. Consequently the
\\(\kappa,\nu)\\ and \\(\psi,\nu)\\ cross-derivatives do **not** vanish
(they do for `ssjplss`); the \\\nu\\-score and the full cross-parameter
Hessian carry the \\(\kappa, \psi, \nu)\\ log-normalizer moments
returned by that quadrature.

**Optimizer.** Only first- and second-order log-likelihood derivatives
are provided, so `available.derivs = 0` and the family is fitted by the
extended Fellner-Schall optimizer rather than full Newton REML.
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) selects
`optimizer = "efs"` automatically; passing it explicitly is recommended
so the fitted object is labelled correctly (which
[`gam.check`](https://rdrr.io/pkg/mgcv/man/gam.check.html) relies on).
The shape \\\psi\\ and asymmetry \\\nu\\ are weakly identified when the
response is diffuse; holding them global (`~ 1`) and keeping the
concentration well identified is the robust default.

The location uses the Fisher-Lee tan-half link (\\\xi \in (-\pi, \pi)\\,
antipode unrepresentable, winding number zero – see
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md) when
the location must wind). The asymmetry rides the tanh link (shared with
[`ssjplss`](https://circstat.github.io/circlss/reference/ssjplss.md),
[`vmftlss`](https://circstat.github.io/circlss/reference/vmftlss.md) and
[`ibslss`](https://circstat.github.io/circlss/reference/ibslss.md)),
bounded to \\(-1, 1)\\.

## References

Jones, M. C. and Pewsey, A. (2005) A family of symmetric distributions
on the circle. *Journal of the American Statistical Association* 100,
1422-1428.

Abe, T., Pewsey, A. and Shimizu, K. (2013) Extending circular
distributions through transformation of argument. *Annals of the
Institute of Statistical Mathematics* 65, 833-858.

Batschelet, E. (1981) *Circular Statistics in Biology*. Academic Press.

Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method
for smoothing parameter optimization with application to Tweedie
location, scale and shape models. *Biometrics* 73, 1071-1081.

## See also

[`jplss`](https://circstat.github.io/circlss/reference/jplss.md),
[`ssjplss`](https://circstat.github.io/circlss/reference/ssjplss.md),
[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
[`ibslss`](https://circstat.github.io/circlss/reference/ibslss.md),
[`vmftlss`](https://circstat.github.io/circlss/reference/vmftlss.md),
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md),
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html)

## Examples

``` r
library(mgcv)
set.seed(1)
n <- 300
x <- runif(n)
xi <- 2 * atan(sin(2 * pi * x))
kappa <- exp(1.0 + 0.4 * cos(2 * pi * x))
psi <- 0.5
nu <- 0.4       # asymmetry
# asymmetric Jones-Pewsey draws by inverse transform on the warped kernel:
# the forward warp g(phi) is closed-form, so gridding phi gives a
# (phi, density) grid to invert -- no root-finding needed for simulation.
phi <- seq(-pi, pi, length.out = 4001)
g <- phi + nu * cos(phi)
dev <- vapply(seq_len(n), function(i) {
  d <- (cosh(kappa[i] * psi) + sinh(kappa[i] * psi) * cos(g))^(1 / psi)
  d[d < 0] <- 0
  cdf <- cumsum((d[-1] + d[-length(d)]) / 2 * diff(phi))
  cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
  approx(cdf, phi, runif(1), rule = 2)$y
}, numeric(1))
y <- atan2(sin(xi + dev), cos(xi + dev))
# smooth location and concentration; global (intercept-only) shape and asymmetry
b <- gam(list(y ~ s(x), ~ s(x), ~ 1, ~ 1), family = ajplss(), optimizer = "efs")
summary(b)
```
