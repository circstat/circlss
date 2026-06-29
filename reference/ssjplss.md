# Sine-skewed Jones-Pewsey location-concentration-shape-skewness family

A general family implementing distributional regression for a circular
(angular) response \\y\\ in radians under the sine-skewed Jones-Pewsey
law \$\$f(y) = c(\kappa, \psi)\\\bigl(\cosh(\kappa\psi) +
\sinh(\kappa\psi)\cos(y - \xi)\bigr)^{1/\psi}\\ \bigl(1 + \lambda
\sin(y - \xi)\bigr),\$\$ the Jones-Pewsey family
([`jplss`](https://circstat.github.io/circlss/reference/jplss.md))
multiplied by a sine-skew factor (Umbach-Jammalamadaka). It adds an
asymmetry axis to the symmetric Jones-Pewsey umbrella.

## Usage

``` r
ssjplss(link = list("tanhalf", "log", "identity", "tanh"))
```

## Arguments

- link:

  Four-element list of link names, for the location, concentration,
  shape and skewness. Currently only the defaults are available:
  `"tanhalf"` for the location \\\xi\\, `"log"` for the concentration
  \\\kappa \> 0\\, `"identity"` for the shape \\\psi \in \mathbb{R}\\,
  and `"tanh"` for the skewness \\\lambda \in (-1, 1)\\.

## Value

An object of class `c("general.family", "extended.family", "family")`
for use with [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (or its
front end
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)).

## Details

It is the first family with **four linear predictors**: the location
\\\xi\\, the concentration \\\kappa\\, the shape \\\psi\\ and the
skewness \\\lambda\\ each get their own, any of which may contain smooth
terms. Used with [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html)
and a list of *four* formulas; the first names the response and models
\\\xi\\, then \\\log\kappa\\, then \\\psi\\, then \\\lambda\\. The shape
and skewness are most often held global, i.e. fitted intercept-only with
`~ 1`.

The sine-skew factor \\1 + \lambda\sin(y-\xi)\\ integrates to 1 against
the symmetric Jones-Pewsey kernel (the first trigonometric moment of the
centered kernel is zero), so it leaves the Jones-Pewsey normalizer
\\c(\kappa,\psi)\\ *unchanged*. At \\\lambda = 0\\ the family reduces
exactly to
[`jplss`](https://circstat.github.io/circlss/reference/jplss.md). The
skewness rides the **tanh** link \\\lambda = \tanh(\eta)\\; as
\\\|\lambda\| \to 1\\ the density touches 0 where \\\lambda\sin(y-\xi) =
-1\\, so the link keeps \\\lambda\\ strictly interior.

**Mode anchor, not mean.** Once \\\lambda \neq 0\\, the location \\\xi\\
is the *mode anchor* of the asymmetric density, not its mean direction;
fitted-direction summaries inherit that reading.

**Normalizer and derivatives.** Because the normalizer is the
Jones-Pewsey one, the family reuses that family's Gauss-Legendre
quadrature machinery wholesale: the \\\kappa\\- and \\\psi\\-score and
Hessian are exactly the Jones-Pewsey terms, the \\(\kappa,\lambda)\\ and
\\(\psi,\lambda)\\ cross-derivatives vanish identically, and only the
\\\xi\\- and \\\lambda\\-directions carry the (elementary) skew terms.

**Optimizer.** Only first- and second-order log-likelihood derivatives
are provided, so `available.derivs = 0` and the family is fitted by the
extended Fellner-Schall optimizer rather than full Newton REML.
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) selects
`optimizer = "efs"` automatically; passing it explicitly is recommended
so the fitted object is labelled correctly (which
[`gam.check`](https://rdrr.io/pkg/mgcv/man/gam.check.html) relies on).
The shape \\\psi\\ and skew \\\lambda\\ are weakly identified when the
response is diffuse; holding them global (`~ 1`) and keeping the
concentration well identified is the robust default.

The location uses the Fisher-Lee tan-half link (\\\xi \in (-\pi, \pi)\\,
antipode unrepresentable, winding number zero – see
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md) when
the location must wind).

## References

Umbach, D. and Jammalamadaka, S. R. (2009) Building asymmetry into
circular distributions. *Statistics & Probability Letters* 79, 659-663.

Abe, T. and Pewsey, A. (2011) Sine-skewed circular distributions.
*Statistical Papers* 52, 683-707.

Jones, M. C. and Pewsey, A. (2005) A family of symmetric distributions
on the circle. *Journal of the American Statistical Association* 100,
1422-1428.

Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method
for smoothing parameter optimization with application to Tweedie
location, scale and shape models. *Biometrics* 73, 1071-1081.

## See also

[`jplss`](https://circstat.github.io/circlss/reference/jplss.md),
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
xi <- 2 * atan(sin(2 * pi * x))
kappa <- exp(1.0 + 0.4 * cos(2 * pi * x))
psi <- 0.5
lambda <- 0.5
# sine-skewed Jones-Pewsey draws by inverse transform on the centered kernel
phi <- seq(-pi, pi, length.out = 4001)
dev <- vapply(seq_len(n), function(i) {
  g <- (cosh(kappa[i] * psi) + sinh(kappa[i] * psi) * cos(phi))^(1 / psi) *
    (1 + lambda * sin(phi))
  g[g < 0] <- 0
  cdf <- cumsum(g); cdf <- cdf / cdf[length(cdf)]
  approx(cdf, phi, runif(1), rule = 2)$y
}, numeric(1))
y <- atan2(sin(xi + dev), cos(xi + dev))
# smooth location and concentration; global (intercept-only) shape and skew
b <- gam(list(y ~ s(x), ~ s(x), ~ 1, ~ 1), family = ssjplss(),
         optimizer = "efs")
summary(b)
```
