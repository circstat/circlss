# Inverse Batschelet location-concentration-skewness-peakedness family

A general family implementing distributional regression for a circular
(angular) response \\y\\ in radians under the inverse Batschelet law
\$\$f(y) = c(\kappa, \lambda)\\\exp\\\bigl(\kappa\cos A\bigr), \qquad
\phi = y - \xi,\\ \kappa \> 0,\\ \nu, \lambda \in (-1, 1),\$\$ where the
angle is warped forward into the von Mises cosine kernel by two
*inverse* maps: a skewness warp \\\phi^\star = t\_\nu^{-1}(\phi)\\
solving \\y - \nu(1 + \cos y) = \phi\\, then a peakedness warp \\u^\star
= s\_\lambda^{-1}(\phi^\star)\\ solving \\u - \tfrac{1+\lambda}{2}\sin u
= \phi^\star\\, with the kernel argument \\A = u^\star -
\tfrac{1-\lambda}{2}\sin u^\star\\. The constant \\c(\kappa, \lambda)\\
depends only on \\\kappa\\ and \\\lambda\\.

## Usage

``` r
ibslss(link = list("tanhalf", "log", "tanh", "tanh"))
```

## Arguments

- link:

  Four-element list of link names, for the location, the concentration,
  the skewness, and the peakedness. Currently only the defaults are
  available: `"tanhalf"` for the location \\\xi\\, `"log"` for the
  concentration \\\kappa \> 0\\, `"tanh"` for the skewness \\\nu \in
  (-1, 1)\\, and `"tanh"` for the peakedness \\\lambda \in (-1, 1)\\.

## Value

An object of class `c("general.family", "extended.family", "family")`
for use with [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (or its
front end
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)).

## Details

The skewness \\\nu\\ tilts the density (mode toward smaller or larger
angles) and the peakedness \\\lambda\\ flattens (\\\lambda \< 0\\) or
sharpens (\\\lambda \> 0\\) the peak; together they give a flexible
skew-and- peaked alternative to the von Mises. Because the warps are
*inverse* maps the transform has no closed form (a monotone solver
inverts it) and, unlike the flat-topped von Mises
([`vmftlss`](https://circstat.github.io/circlss/reference/vmftlss.md), a
forward warp), the score requires implicit differentiation.

It has **four linear predictors**: the location \\\xi\\, the
concentration \\\kappa\\, the skewness \\\nu\\, and the peakedness
\\\lambda\\ each get their own, any of which may contain smooth terms.
Used with [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) and a
list of *four* formulas; the first names the response and models
\\\xi\\, the second models \\\log\kappa\\, the third \\\nu\\, the fourth
\\\lambda\\. The skewness and peakedness are most often held global,
i.e. fitted intercept-only with `~ 1`.

The two shape parameters index a family of skew and/or peaked circular
laws and recover the von Mises exactly:

- \\\nu = \lambda = 0\\: the von Mises
  ([`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md)),
  \\\propto \exp(\kappa\cos(y-\xi))\\;

- \\\nu \ne 0\\: a skewed (asymmetric) law;

- \\\lambda \> 0\\: a sharper-than-von-Mises peak; \\\lambda \< 0\\: a
  flatter, flat-topped peak;

- \\\kappa \to 0\\: the circular uniform.

**Mode anchor.** As with
[`ssjplss`](https://circstat.github.io/circlss/reference/ssjplss.md),
once \\\nu \ne 0\\ the location \\\xi\\ is the *mode anchor*, not the
mean direction.

**Warps and normalizer.** The two inverse warps are inverted by a
vectorized monotone Newton iteration with a bisection mop-up near the
boundary (\\\nu, \lambda \to \pm 1\\); the per-observation score reuses
the solved roots \\\phi^\star, u^\star\\ and the warp slopes. The
normalizer \\c(\kappa, \lambda) = (1-\lambda) / \[(1+\lambda)\\2\pi
I_0(\kappa) - 2\lambda\int e^{\kappa\cos B}\]\\ has no elementary form:
the integral is an overflow-safe equispaced trapezoid (the warped kernel
is smooth and periodic, so the trapezoid is spectrally accurate) and the
\\\kappa\\- and \\\lambda\\-derivatives of \\\log c\\ are central finite
differences of it.

**Optimizer.** Only first- and second-order log-likelihood derivatives
are provided, so `available.derivs = 0` and the family is fitted by the
extended Fellner-Schall optimizer rather than full Newton REML.
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) selects
`optimizer = "efs"` automatically; passing it explicitly is recommended
so the fitted object is labelled correctly (which
[`gam.check`](https://rdrr.io/pkg/mgcv/man/gam.check.html) relies on).
The skewness and peakedness are weakly identified when the response is
diffuse; holding them global (`~ 1`) and keeping the concentration well
identified is the robust default.

The location uses the Fisher-Lee tan-half link (\\\xi \in (-\pi, \pi)\\,
antipode unrepresentable, winding number zero – see
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md) when
the location must wind). The skewness and peakedness both ride the tanh
link (shared with
[`ssjplss`](https://circstat.github.io/circlss/reference/ssjplss.md) and
[`vmftlss`](https://circstat.github.io/circlss/reference/vmftlss.md)),
bounded to \\(-1, 1)\\.

## References

Jones, M. C. and Pewsey, A. (2012) Inverse Batschelet distributions for
circular data. *Biometrics* 68, 183-193.

Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method
for smoothing parameter optimization with application to Tweedie
location, scale and shape models. *Biometrics* 73, 1071-1081.

## See also

[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
[`vmftlss`](https://circstat.github.io/circlss/reference/vmftlss.md),
[`ssjplss`](https://circstat.github.io/circlss/reference/ssjplss.md),
[`jplss`](https://circstat.github.io/circlss/reference/jplss.md),
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md),
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html)

## Examples

``` r
library(mgcv)
set.seed(1)
n <- 300
x <- runif(n)
xi <- 2 * atan(sin(2 * pi * x))
kappa <- exp(1.0 + 0.5 * cos(2 * pi * x))
nu <- 0.5      # skew
lmbd <- 0.3    # peakedness
# Inverse Batschelet draws: the forward warps phi(u) are closed-form and
# monotone, so gridding the base angle u gives a (phi, density) grid to invert
# -- no root-finding needed for simulation.
u <- seq(-pi, pi, length.out = 4001)
phistar <- u - 0.5 * (1 + lmbd) * sin(u)
phi <- phistar - nu * (1 + cos(phistar))
A <- u - 0.5 * (1 - lmbd) * sin(u)
dev <- vapply(seq_len(n), function(i) {
  g <- exp(kappa[i] * cos(A))
  cdf <- cumsum((g[-1] + g[-length(g)]) / 2 * diff(phi))
  cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
  approx(cdf, phi, runif(1), rule = 2)$y
}, numeric(1))
y <- atan2(sin(xi + dev), cos(xi + dev))
# smooth location and concentration; global (intercept-only) skew and peakedness
b <- gam(list(y ~ s(x), ~ s(x), ~ 1, ~ 1), family = ibslss(), optimizer = "efs")
summary(b)
```
