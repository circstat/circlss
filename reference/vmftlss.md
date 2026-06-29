# Flat-topped von Mises location-concentration-shape family

A general family implementing distributional regression for a circular
(angular) response \\y\\ in radians under the flat-topped von Mises law
\$\$f(y) = c(\kappa, \nu)\\\exp\\\bigl(\kappa\cos(\phi +
\nu\sin\phi)\bigr), \qquad \phi = y - \mu,\\ \kappa \> 0,\\ \nu \in (-1,
1),\$\$ with \\c(\kappa,\nu)\\ a normalizing constant. The peakedness
\\\nu\\ warps the angle forward inside the cosine: the warp \\B = \phi +
\nu\sin\phi\\ is odd in \\\phi\\, so the density stays symmetric about
\\\mu\\ while the peak is sharpened (\\\nu \> 0\\) or flattened (\\\nu
\< 0\\) – a flat-topped or sharply-peaked alternative to the von Mises
that keeps the same mean direction.

## Usage

``` r
vmftlss(link = list("tanhalf", "log", "tanh"))
```

## Arguments

- link:

  Three-element list of link names, for the mean direction, the
  concentration, and the peakedness. Currently only the defaults are
  available: `"tanhalf"` for the location \\\mu\\, `"log"` for the
  concentration \\\kappa \> 0\\, and `"tanh"` for the peakedness \\\nu
  \in (-1, 1)\\.

## Value

An object of class `c("general.family", "extended.family", "family")`
for use with [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (or its
front end
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)).

## Details

It has **three linear predictors**: the mean direction \\\mu\\, the
concentration \\\kappa\\, and the peakedness \\\nu\\ each get their own,
any of which may contain smooth terms. Used with
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) and a list of
*three* formulas; the first names the response and models \\\mu\\, the
second models \\\log\kappa\\, and the third models \\\nu\\. The
peakedness is most often held global, i.e. fitted intercept-only with
`~ 1`.

The peakedness \\\nu\\ indexes a family of symmetric circular laws and
recovers the von Mises exactly:

- \\\nu = 0\\: the von Mises
  ([`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md)),
  \\\propto \exp(\kappa\cos(y-\mu))\\;

- \\\nu \> 0\\: a sharper-than-von-Mises peak;

- \\\nu \< 0\\: a flatter, flat-topped peak;

- \\\kappa \to 0\\: the circular uniform.

**Normalizer.** The constant \\c(\kappa, \nu)\\ has no elementary form.
The log-density value *and* the \\\kappa\\- and \\\nu\\-score and
Hessian carry kernel-weighted moments of the normalizer, evaluated by an
adaptive equispaced trapezoid (the warped kernel is smooth and periodic,
so the trapezoid is spectrally accurate). The grid resolution grows with
\\\kappa\\ and \\\|\nu\|\\ so the sharpened peak is always resolved.

**Optimizer.** Only first- and second-order log-likelihood derivatives
are provided, so `available.derivs = 0` and the family is fitted by the
extended Fellner-Schall optimizer rather than full Newton REML.
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) selects
`optimizer = "efs"` automatically; passing it explicitly is recommended
so the fitted object is labelled correctly (which
[`gam.check`](https://rdrr.io/pkg/mgcv/man/gam.check.html) relies on).
The peakedness \\\nu\\ is weakly identified when the response is
diffuse; holding it global (`~ 1`) and keeping the concentration well
identified is the robust default.

The mean direction uses the Fisher-Lee tan-half link (\\\mu \in (-\pi,
\pi)\\, antipode unrepresentable, winding number zero – see
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md) when
the mean direction must wind). The peakedness rides the tanh link
(shared with
[`ssjplss`](https://circstat.github.io/circlss/reference/ssjplss.md)'s
skewness), bounded to \\(-1, 1)\\.

## References

Batschelet, E. (1981) *Circular Statistics in Biology*. Academic Press.

Abe, T., Pewsey, A. and Shimizu, K. (2013) Extending circular
distributions through transformation of argument. *Annals of the
Institute of Statistical Mathematics* 65, 833-858.

Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method
for smoothing parameter optimization with application to Tweedie
location, scale and shape models. *Biometrics* 73, 1071-1081.

## See also

[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
[`jplss`](https://circstat.github.io/circlss/reference/jplss.md),
[`ssjplss`](https://circstat.github.io/circlss/reference/ssjplss.md),
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
nu <- 0.5
# flat-topped vM draws by inverse transform on the centered kernel, then shift
phi <- seq(-pi, pi, length.out = 4001)
dev <- vapply(seq_len(n), function(i) {
  B <- phi + nu * sin(phi)
  g <- exp(kappa[i] * cos(B))
  cdf <- cumsum(g); cdf <- cdf / cdf[length(cdf)]
  approx(cdf, phi, runif(1), rule = 2)$y
}, numeric(1))
y <- atan2(sin(mu + dev), cos(mu + dev))
# smooth mean direction and concentration; global (intercept-only) peakedness
b <- gam(list(y ~ s(x), ~ s(x), ~ 1), family = vmftlss(), optimizer = "efs")
summary(b)
```
