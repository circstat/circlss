# Classical circular regression

The closed-form, unpenalized counterpart of
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md):
the textbook circular regressions that the literature reports. `circ_lm`
is **parametric only** – a smooth term
([`s()`](https://rdrr.io/pkg/mgcv/man/s.html),
[`te()`](https://rdrr.io/pkg/mgcv/man/te.html), ...) is an error
pointing to
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md) –
and carries no `family` argument: `"cl"` is von Mises by construction,
while `"cc"` and `"lc"` are ordinary least squares with a residual
concentration reported as a summary. Reach for `circ_lm` for the
textbook fit or a fast unpenalized baseline; reach for
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
for penalized smooths, per-parameter modelling, or any family beyond von
Mises.

## Usage

``` r
circ_lm(
  formula,
  data,
  type = c("cl", "cc", "lc"),
  order = 1L,
  init = NULL,
  tol = 1e-08,
  maxit = 1000L,
  se = c("asymptotic", "bootstrap"),
  R = 999L,
  verbose = FALSE
)

# S3 method for class 'circ_lm'
coef(object, ...)

# S3 method for class 'circ_lm'
fitted(object, ...)

# S3 method for class 'circ_lm'
residuals(object, ...)

# S3 method for class 'circ_lm'
logLik(object, ...)

# S3 method for class 'circ_lm'
predict(object, newdata, type = c("direction", "kappa"), ...)

# S3 method for class 'circ_lm'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- formula:

  A formula, or (for `type = "cl"`) a list of one or two formulas. The
  first formula is two-sided and names the response.

- data:

  A data frame holding the response and covariates.

- type:

  For `circ_lm`, which classical fit. `"cl"` circular response on linear
  covariate(s) (Fisher–Lee von Mises); `"cc"` circular on circular
  (harmonic); `"lc"` linear on circular (harmonic). Hyphenated spellings
  (`"c-l"`, `"c-c"`, `"l-c"`) are accepted. For `predict` on a `"cl"`
  fit, the quantity returned: `"direction"` (the mean direction, the
  default) or `"kappa"` (the fitted concentration).

- order:

  Order of the trigonometric polynomial for `"cc"` and `"lc"` (number of
  harmonics of the angular predictor). Ignored for `"cl"`.

- init:

  Starting values for the `"cl"` iteration. `NULL` (default) starts
  cold: all coefficients zero, so \\\kappa \equiv 1\\. A named list
  `list(beta=, alpha=, gamma=)` sets explicit starts – any component
  omitted falls back to its cold value; `beta` and `gamma` take one
  value per covariate, `alpha` a single number. This lets you seed the
  joint (mixed) fit with estimates from separately fitted mean-only and
  kappa-only models, as Fisher (1993, Sec. 6.4.4) suggests. A bare
  numeric vector is taken, as before, as the mean-direction coefficients
  (`beta`).

- tol, maxit, verbose:

  IRLS convergence tolerance, iteration cap, and per-iteration logging
  (`"cl"` only).

- se:

  How standard errors are computed (`"cl"` only). `"asymptotic"`
  (default) uses the expected-information formulae; `"bootstrap"`
  replaces them with a parametric bootstrap – simulate from the fitted
  vM(\\\mu_i, \kappa_i\\), refit, and take the spread of the estimates –
  which Fisher (1993, Sec. 8.4) recommends below \\n \approx 25\\ –
  \\30\\, where the asymptotic SEs are unreliable. Stochastic, so set a
  seed for reproducibility.

- R:

  Number of bootstrap resamples when `se = "bootstrap"`.

- object, x:

  A fitted `circ_lm` model.

- ...:

  Unused.

- newdata:

  A data frame of new predictor values. For `predict`, omitting it
  returns the fitted values.

- digits:

  Number of significant digits for `print`.

## Value

An object of class `circ_lm`: a list whose contents depend on `type`.
`"cl"` carries `mu`, `kappa`, the coefficient vectors
(`beta`/`alpha`/`gamma`) with standard errors, and `loglik`/`aic`/`bic`;
`"cc"` carries the cos/sin `coefficients`, `rho`, residual `kappa`, and
higher-order-test `p_values`; `"lc"` carries the `coefficients`,
per-harmonic `amplitude`/`phase` (with delta-method SEs), and the usual
least-squares fit metrics. `predict`, `coef`, `fitted`, `residuals`, and
(for `"cl"`/`"lc"`) `logLik` methods are provided.

## Details

**cl – circular response, linear covariates.** The Fisher and Lee (1992)
von Mises model fitted by Green's (1984) IRLS, with the concentration
extensions of Fisher (1993) Sec. 6.4. The mean direction is \\\mu_i =
\mu_0 + 2\\\mathrm{atan}(x_i^\top\beta)\\ (the offset is *outside* the
link, as the textbooks and
[`circular::lm.circular`](https://rdrr.io/pkg/circular/man/lm.circular.html)
write it – distinct from `circ_gam`, which puts the intercept inside the
link). Because the reference direction \\\mu_0\\ is a free angle, this
model is exactly rotation-equivariant and never sits against the
tan-half wall, so `circ_lm` has no `center` argument (the `circ_gam`
counterpart to that property); likewise the `cc` harmonic model recovers
its mean by `atan2`, with no wall. A one- or two-formula list selects
the sub-model the same way `circ_gam` reads a formula list:
`list(theta ~ x, ~ 1)` (or just `theta ~ x`) models the mean with
constant \\\kappa\\; `list(theta ~ 1, ~ z)` models \\\log\kappa =
\alpha + z^\top\gamma\\ with constant \\\mu\\; `list(theta ~ x, ~ x)` is
the mixed model. The mixed model ties \\\mu\\ and \\\kappa\\ to one
shared design. Any of these may carry several covariates
(`theta ~ x + z`); only `cc`/`lc` are single-predictor. The mixed
iteration starts cold (\\\kappa \equiv 1\\) by default; pass
`init = list(beta=, alpha=, gamma=)` to seed it from your own starting
values – e.g. the estimates of separately fitted mean-only and
kappa-only models, the two-stage start Fisher (1993) Sec. 6.4.4
describes.

**cc – circular on circular.** The Sarma and Jammalamadaka (1993)
harmonic fit: `cos(theta)` and `sin(theta)` regressed by least squares
on a degree-`order` trigonometric polynomial of the angular predictor
and reassembled as \\\hat\mu = \mathrm{atan2}(\hat s, \hat c)\\, with
the circular correlation \\\rho\\, a residual concentration, and the
test for significance of the next harmonic order.

**lc – linear on circular.** Harmonic regression of a linear response on
a Fourier basis of the angular predictor, reporting each harmonic's
amplitude and phase with delta-method standard errors.

**Parity with circular::lm.circular.** The regression outputs
(\\\beta\\, \\\mu\\, the cos/sin coefficients, \\\rho\\, fitted values,
p-values) reproduce `circular` to machine precision. The reported
\\\kappa\\ can differ slightly: it is the only quantity passing through
the inverse Bessel ratio, and circlss returns the machine-precision
inverse (the exact \\\kappa\\ solving \\A_1(\kappa) = R\\) where
`circular` uses the classical piecewise approximation – a gap at that
approximation's error level (~1e-3), largest at high concentration. One
deliberate departure: the reported `logLik` (and `aic`/`bic`) is the
*full* von Mises log-likelihood with every estimated parameter counted
(\\\mu_0\\ and \\\kappa\\ included), so it exceeds `lm.circular`'s
printed `log.lik` by the \\n\log 2\pi\\ normalisation `circular` drops –
putting circlss's AIC on the standard scale, comparable to
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
or a `glm`.

## References

Fisher, N. I. and Lee, A. J. (1992) Regression models for an angular
response. *Biometrics* 48, 665-677.

Fisher, N. I. (1993) *Statistical Analysis of Circular Data*. Cambridge
University Press.

Sarma, Y. and Jammalamadaka, S. R. (1993) Circular regression. In
*Statistical Sciences and Data Analysis*, 109-128. VSP, Utrecht.

Pewsey, A., Neuhaeuser, M. and Ruxton, G. D. (2013) *Circular Statistics
in R*. Oxford University Press.

## See also

[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
for the penalized-spline distributional models;
[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md) for the
von Mises family used there.

## Examples

``` r
set.seed(1)
n <- 80
x <- rnorm(n)
theta <- (1 + 2 * atan(1.5 * x) + rnorm(n) / 4) %% (2 * pi)
dat <- data.frame(theta = theta, x = x)

## cl: Fisher-Lee mean direction (constant kappa)
m <- circ_lm(theta ~ x, dat, type = "cl")
m
predict(m, data.frame(x = c(-1, 0, 1)))

## cl: mixed model -- mean and log-kappa both linear in x
circ_lm(list(theta ~ x, ~ x), dat, type = "cl")

## cl: bootstrap SEs (Fisher 1993 Sec. 8.4) -- preferred at small n
set.seed(1)
circ_lm(theta ~ x, dat, type = "cl", se = "bootstrap", R = 199)

## cl: seed the mixed fit from separately fitted mean-only / kappa-only models
b0 <- circ_lm(theta ~ x, dat, type = "cl")
k0 <- circ_lm(list(theta ~ 1, ~ x), dat, type = "cl")
circ_lm(list(theta ~ x, ~ x), dat, type = "cl",
        init = list(beta = b0$beta, alpha = k0$alpha, gamma = k0$gamma))

## cl: several covariates (mean and concentration share the design)
dat$z <- rnorm(n)
circ_lm(list(theta ~ x + z, ~ x + z), dat, type = "cl")

## cc / lc: harmonic fits on an angular predictor
phi <- runif(n, 0, 2 * pi)
dcc <- data.frame(psi = (phi / 2 + rnorm(n) / 5) %% (2 * pi), phi = phi)
circ_lm(psi ~ phi, dcc, type = "cc", order = 1)

dlc <- data.frame(y = 5 + 2 * cos(phi) + rnorm(n) / 2, phi = phi)
circ_lm(y ~ phi, dlc, type = "lc", order = 1)
```
