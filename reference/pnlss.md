# Projected normal location family

A general family implementing distributional regression for a circular
(angular) response \\y\\ in radians under the projected normal law:
\\y\\ is the angle of a bivariate normal with mean \\(\mu_1, \mu_2)\\
and identity covariance. Each mean component gets its own linear
predictor with an identity link, so
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) is called with a
list of two formulas: the first names the response and models \\\mu_1\\,
the second models \\\mu_2\\.

## Usage

``` r
pnlss(link = list("identity", "identity"))
```

## Arguments

- link:

  Two-element list of link names for the two Cartesian mean components.
  Only `"identity"` is available for both.

## Value

An object of class `c("general.family", "extended.family", "family")`
for use with [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (or its
front end
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)).

## Details

The fitted mean direction is \\\mathrm{atan2}(\mu_2, \mu_1)\\ and the
implied concentration grows with \\\\\mu\\\\. Because the direction is
assembled from two unconstrained components, there is no link branch
cut: unlike the tan-half parameterization of
[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md), the
fitted mean direction can cross any angle and can *wind* around the
circle (e.g. \\\mu(\varphi) = \varphi\\ with a cyclic covariate), which
makes `pnlss` the natural family for circular-circular regression with
rotation-type association. The trade-off is interpretability: location
and concentration are entangled in \\(\mu_1, \mu_2)\\ rather than
separated into distinct parameters.

Fitted values and `predict(..., type = "response")` return the two
*Cartesian components* \\(\mu_1, \mu_2)\\ as columns, matching the
identity links; compute the direction with `atan2(fit[, 2], fit[, 1])`
and the concentration scale with `sqrt(rowSums(fit^2))`.

Log-likelihood derivatives up to fourth order are implemented, so the
family supports full Newton REML (`method = "REML"`);
`optimizer = "efs"` also works. `"pearson"` residuals alias `"deviance"`
(both are the signed root of twice the log-likelihood gap to the
fitted-direction mode).

## References

Presnell, B., Morrison, S. P. and Littell, R. C. (1998) Projected
multivariate linear models for directional data. *Journal of the
American Statistical Association* 93, 1068-1077.

Wood, S. N., Pya, N. and Saefken, B. (2016) Smoothing parameter and
model selection for general smooth models. *Journal of the American
Statistical Association* 111, 1548-1575.

## See also

[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html)

## Examples

``` r
library(mgcv)
set.seed(1)
n <- 300
x <- runif(n)
m1 <- 1.5 * sin(2 * pi * x) + 0.5
m2 <- 1.2 * cos(2 * pi * x)
y <- atan2(m2 + rnorm(n), m1 + rnorm(n))  # exact projected normal draws
b <- gam(list(y ~ s(x), ~ s(x)), family = pnlss(), method = "REML")
summary(b)
fv <- fitted(b)
direction <- atan2(fv[, 2], fv[, 1])
```
