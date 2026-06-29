# Kato-Jones Mobius four-parameter family

A general family implementing distributional regression for a circular
(angular) response \\y\\ in radians under the Kato-Jones (2015) law
\$\$g(y) = \frac{1}{2\pi}\left\[1 + \frac{2\gamma\\ (\cos(y - \mu) -
\rho\cos\lambda)} {1 + \rho^2 - 2\rho\cos(y - \mu -
\lambda)}\right\],\$\$ a tractable four-parameter family obtained from a
Mobius transformation of the circle. The parameters control the first
two trigonometric moments: \\\mu\\ the mean direction, \\\gamma\\ the
mean resultant length, and \\(\rho, \lambda)\\ the magnitude and phase
of the second-order moment. It is the last and most general family in
the package, bringing both peakedness and skewness through a single
construction.

## Usage

``` r
kjlss(link = list("tanhalf", "logit", "identity", "identity"))
```

## Arguments

- link:

  Four-element list of link names, for the location, the mean resultant
  length, and the two disc-chart coordinates. Currently only the
  defaults are available: `"tanhalf"` for the location \\\mu\\,
  `"logit"` for the mean resultant length \\\gamma \in (0, 1)\\, and
  `"identity"` for the unconstrained chart coordinates \\u_1, u_2 \in
  \mathbb{R}\\.

## Value

An object of class `c("general.family", "extended.family", "family")`
for use with [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (or its
front end
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)).

## Details

It has **four linear predictors**. To keep the second-order pair
\\(\rho, \lambda)\\ inside its feasible region for *every* coefficient
vector, the family is parameterized by the **disc chart**: the Theorem-1
feasible set for the Cartesian shape pair \\(a, b) = (\rho\cos\lambda,
\rho\sin\lambda)\\ is the disc of centre \\(\gamma, 0)\\ and radius
\\1-\gamma\\, and the chart \$\$(a, b) = (\gamma, 0) +
(1-\gamma)\\u/\sqrt{1 + \lVert u\rVert^2}, \qquad u = (u_1, u_2) \in
\mathbb{R}^2,\$\$ maps an unconstrained \\u\\ onto its interior. So the
smooths ride the unconstrained chart coordinates \\u_1, u_2\\ (identity
links) and the coupled feasibility constraint can never be violated.
Used with [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) and a
list of *four* formulas: the first names the response and models
\\\mu\\, then \\\gamma\\, then \\u_1\\, then \\u_2\\. The chart
coordinates are most often held global, i.e. fitted intercept-only with
`~ 1`.

**Exact normalizer.** Unlike every other shape family here, the
Kato-Jones normalizer is exactly \\2\pi\\: the density is a
first-/second-trigonometric- moment perturbation of the circular uniform
that integrates to 1 with no special function, so the log-density and
all of its derivatives are elementary rational functions of \\(a, b)\\
and \\\cos(y-\mu)\\, \\\sin(y-\mu)\\ – there is no quadrature, Bessel or
Gamma term. The derivatives are taken with respect to the chart
coordinates by pushing the Cartesian scores through the chart Jacobian
and Hessian.

**Special and nested cases.** \\u = 0\\ is exactly the wrapped Cauchy
\\\mathrm{WC}(\mu, \gamma)\\
([`wclss`](https://circstat.github.io/circlss/reference/wclss.md)), the
natural reduced model, so intercept-only \\u_1, u_2\\ is the wrapped
Cauchy with covariate-free second-order shape. \\\gamma \to 0\\ gives
the circular uniform and \\\rho \to 0\\ (\\u\\ radial to 0) gives a
cardioid
([`cardlss`](https://circstat.github.io/circlss/reference/cardlss.md)).

**Mean direction, not mode.** Once \\\rho \neq 0\\ the density is
asymmetric and \\\mu\\ is the direction of the first trigonometric
moment, not the mode; fitted-direction summaries inherit that reading.

**Optimizer.** Only first- and second-order log-likelihood derivatives
are provided, so `available.derivs = 0` and the family is fitted by the
extended Fellner-Schall optimizer rather than full Newton REML.
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) selects
`optimizer = "efs"` automatically; passing it explicitly is recommended
so the fitted object is labelled correctly (which
[`gam.check`](https://rdrr.io/pkg/mgcv/man/gam.check.html) relies on).
The chart coordinates \\u_1, u_2\\ are weakly identified when the
response is diffuse; holding them global (`~ 1`) and keeping the mean
direction and \\\gamma\\ well identified is the robust default. With
four linear predictors and two flat shape directions, a cyclic
(`bs = "cc"`) model of all four can exceed `mgcv`'s fixed
extended-Fellner-Schall iteration cap, so prefer parametric or
thin-plate location terms.

The mean direction uses the Fisher-Lee tan-half link (\\\mu \in (-\pi,
\pi)\\, antipode unrepresentable, winding number zero – see
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md) when
the mean direction must wind).

## References

Kato, S. and Jones, M. C. (2015) A tractable and interpretable
four-parameter family of unimodal distributions on the circle.
*Biometrika* 102, 181-190.

Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method
for smoothing parameter optimization with application to Tweedie
location, scale and shape models. *Biometrics* 73, 1071-1081.

## See also

[`wclss`](https://circstat.github.io/circlss/reference/wclss.md),
[`cardlss`](https://circstat.github.io/circlss/reference/cardlss.md),
[`jplss`](https://circstat.github.io/circlss/reference/jplss.md),
[`ssjplss`](https://circstat.github.io/circlss/reference/ssjplss.md),
[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md),
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html)

## Examples

``` r
library(mgcv)
set.seed(1)
n <- 300
x <- runif(n)
mu <- 2 * atan(0.8 * sin(2 * pi * x))
gamma <- plogis(0.8)         # mean resultant length ~ 0.69
u1 <- 0.4; u2 <- -0.5        # a fixed second-order shape
# Kato-Jones draws by inverse transform on the centered kernel (the chart maps
# (gamma, u1, u2) to the Cartesian shape pair (a, b), guaranteed feasible)
r <- sqrt(1 + u1^2 + u2^2); om <- 1 - gamma
a <- gamma + om * u1 / r; b <- om * u2 / r
phi <- seq(-pi, pi, length.out = 4001)
D <- 1 + a^2 + b^2 - 2 * a * cos(phi) - 2 * b * sin(phi)
g <- pmax(1 + 2 * gamma * (cos(phi) - a) / D, 0)
cdf <- cumsum(g); cdf <- cdf / cdf[length(cdf)]
dev <- approx(cdf, phi, runif(n), rule = 2)$y
y <- atan2(sin(mu + dev), cos(mu + dev))
# smooth location; global (intercept-only) gamma and chart coordinates
b <- gam(list(y ~ s(x), ~ 1, ~ 1, ~ 1), family = kjlss(), optimizer = "efs")
summary(b)
```
