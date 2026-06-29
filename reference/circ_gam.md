# Circular-response GAM

A wrapper around [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) for
circular regression with the circlss families. It supplies the defaults
a circular fit needs: the cyclic-smooth knots, `method = "REML"`, a
trailing `~ 1` fill for any distribution parameters left without a
formula, and response-scale output columns named by the family's
parameters. Formulas and smoothing bases are passed to mgcv unchanged.

## Usage

``` r
circ_gam(
  formula,
  data,
  family = vmlss(),
  method = "REML",
  knots = NULL,
  weights = NULL,
  subset = NULL,
  na.action,
  offset = NULL,
  center = TRUE,
  ...
)

# S3 method for class 'circ_gam'
predict(object, ...)

# S3 method for class 'circ_gam'
fitted(object, ...)

# S3 method for class 'circ_gam'
print(x, ...)
```

## Arguments

- formula:

  A formula, or a list of formulas (one per distribution parameter),
  exactly as [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) expects for
  a general family. The first formula is two-sided and names the
  response. If fewer formulas than the family has parameters are
  supplied, the remainder are filled with `~ 1` (intercept-only), so
  `theta ~ s(x)` with
  [`vmlss()`](https://circstat.github.io/circlss/reference/vmlss.md)
  models the mean direction and holds the concentration constant.

- data:

  A data frame holding the response and covariates. Required: the
  covariate values are what set the cyclic-smooth knots.

- family:

  Any mgcv family. A circlss *circular* family
  ([`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
  [`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md), ...)
  selects a circular (angular) response; a circlss *linear*
  location-scale family
  ([`gausslss`](https://circstat.github.io/circlss/reference/gausslss.md),
  [`gammalss`](https://circstat.github.io/circlss/reference/gammalss.md))
  selects a real-valued response on a circular covariate, the
  linear–circular leg. Either switches on the named, response-scale
  output, the trailing `~ 1` fill, and the geometry-aware
  `print`/`plot`. An ordinary mgcv family such as
  [`gaussian`](https://rdrr.io/r/stats/family.html) is forwarded
  unchanged (it still gets the cyclic-knot default).

- method:

  Smoothing-parameter selection criterion; defaults to `"REML"` (mgcv's
  own default is `"GCV.Cp"`).

- knots:

  Optional knot list passed to `gam`. Knots you supply are always
  respected; for any cyclic smooth (`bs = "cc"` or `"cp"`) you do *not*
  supply, `circ_gam` fills the knots with one full radian period
  bracketing the data – `c(-pi, pi)` when the covariate takes negative
  values, else `c(0, 2*pi)`.

- weights:

  Optional prior weights on the observations, one per row of `data`, as
  in [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html).

- subset:

  Optional vector selecting the rows of `data` to fit, as in
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html).

- na.action:

  How missing values are handled, as in
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html); when omitted, mgcv's
  option-driven default
  ([`na.omit`](https://rdrr.io/r/stats/na.fail.html)) applies.

- offset:

  Optional model offset, as in
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html).

- center:

  For a circular response on the tan-half link (`vmlss` and the other
  tan-link families; never `pnlss`, which has no wall), rotate the
  response to a frame where the unreachable antipode \\\theta = \pi\\
  lies away from the data before fitting, then report all directions
  back in the original frame. `TRUE` (default) chooses the reference
  automatically (the centre of the occupied arc) and is an exact no-op
  when the data already clear the wall; `FALSE` disables it; a numeric
  value sets the reference angle directly. Only response-scale
  directions (`predict(type = "response")`, `fitted`, `plot`) are
  rotated back; link-scale output (`coef`, `predict(type = "link")`)
  stays in the centred frame. The applied rotation is stored in
  `fit$circ_center`. A mean that must *wind* across the wall cannot be
  centred away – use
  [`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md)
  there.

- ...:

  Further arguments passed verbatim to
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) (`optimizer`,
  `control`, `sp`, `select`, ...).

- object, x:

  A fitted `circ_gam` model.

## Value

A fitted model of class `c("circ_gam", "gam", "glm", "lm")`: the
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) object, with
circular-aware `predict`, `fitted`, `print` and `plot` methods. Every
other mgcv method (`summary`, `AIC`, `logLik`, `gam.check`, ...) is
inherited unchanged.

## Details

**Knots.** A cyclic smooth's period is its knot span. Because the data
are in radians the period is \\2\pi\\ and the wrap points are \\\pm\pi\\
(or \\0\\ and \\2\pi\\); without explicit knots,
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) wraps the cyclic basis at
the observed data range instead. `circ_gam` fills these knots
automatically, and stops if a cyclic covariate is not on a single radian
branch (convert it first with
[`rad`](https://circstat.github.io/circlss/reference/rad.md) or
[`wrap`](https://circstat.github.io/circlss/reference/rad.md)). Knots
supplied through `knots` override this for that covariate.

**Output columns.** For a circlss family, `predict(type = "response")`
and [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) label their
columns with the family's parameters – for example `mu, kappa` for
[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
`mu1, mu2` for
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md), and
`xi, kappa, psi, lambda` for
[`ssjplss`](https://circstat.github.io/circlss/reference/ssjplss.md).

**Optimizer.** [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) selects
the extended Fellner–Schall optimizer automatically for families that
supply only first- and second-order derivatives
(`available.derivs = 0`); pass `optimizer = "efs"` through `...` to
force it on a family that would otherwise use full Newton REML.

## See also

[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md),
[`rad`](https://circstat.github.io/circlss/reference/rad.md),
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html)

## Examples

``` r
library(mgcv)
set.seed(1); n <- 300

## circular-linear: von Mises, both parameters smooth in a covariate
x <- runif(n); y <- 2 * atan(1.2 * sin(2 * pi * x)) + rnorm(n) / 3
b <- circ_gam(list(y ~ s(x), ~ s(x)), data = data.frame(y, x), family = vmlss())
head(predict(b, type = "response"))   # columns named mu, kappa
plot(b)

## circular-circular: projected normal, cyclic covariate -- knots auto-pinned
phi <- runif(n, -pi, pi); yc <- 2 * atan(0.9 * sin(phi)) + rnorm(n) / 3
b2 <- circ_gam(list(yc ~ s(phi, bs = "cc"), ~ s(phi, bs = "cc")),
               data = data.frame(yc, phi), family = pnlss())

## linear-circular: a real response on a cyclic covariate (the "can")
yl <- 2 + 1.5 * sin(phi) + rnorm(n) / 3
b3 <- circ_gam(list(yl ~ s(phi, bs = "cc"), ~ s(phi, bs = "cc")),
               data = data.frame(yl, phi), family = gausslss())
plot(b3, view = "both")
```
