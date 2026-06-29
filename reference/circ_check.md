# Diagnostic panels for a circular regression fit

The diagnostic display for circlss fits – the circular analogue of
`stats:::plot.lm` and
[`gam.check`](https://rdrr.io/pkg/mgcv/man/gam.check.html), and the
counterpart to the *effect*-display methods
[`plot.circ_gam`](https://circstat.github.io/circlss/reference/plot.circ_gam.md)
/
[`plot.circ_lm`](https://circstat.github.io/circlss/reference/plot.circ_lm.md)
(which answer “what did the model fit?”; `circ_check` answers “is the
fit any good?”). It lays out a panel grid of
[`circ_resid`](https://circstat.github.io/circlss/reference/circ_resid.md)-based
diagnostics, prints a goodness-of-fit table, and returns the statistics
invisibly. Dispatches over both front doors
([`circ_lm`](https://circstat.github.io/circlss/reference/circ_lm.md)
and
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)).

## Usage

``` r
circ_check(object, which = NULL, nsim = 1000L, rug = TRUE, ...)
```

## Arguments

- object:

  A fitted
  [`circ_lm`](https://circstat.github.io/circlss/reference/circ_lm.md)
  or
  [`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
  model.

- which:

  Character vector of panel keys to draw; `NULL` (default) picks a
  sensible set for the response type, and `"all"` draws every panel.
  Keys: `"rose"` (rose diagram of angular residuals), `"obsfit"`
  (observed vs fitted with the wrapped calibration diagonal),
  `"residcov"` (residual vs covariate, on a circular axis when the
  covariate is cyclic), `"qq.unif"` (quantile-residual uniform Q-Q with
  the Watson \\U^2\\ p-value), `"qq.norm"` (normal Q-Q of the deviance
  residuals), `"scaleloc"` (\\\sqrt{\|\\\mathrm{deviance\\
  residual}\\\|}\\ vs the fitted location, a concentration-adequacy
  check), `"hist"` (deviance-residual histogram with the standard-normal
  reference), and `"cook"` (residuals vs leverage with Cook's-distance
  contours;
  [`circ_lm`](https://circstat.github.io/circlss/reference/circ_lm.md)
  only – see Details).

- nsim:

  Number of simulation replicates for the quantile residual when the
  family has no closed-form distribution function; see
  [`circ_resid`](https://circstat.github.io/circlss/reference/circ_resid.md).

- rug:

  Add a covariate rug to the residual-vs-covariate panel.

- ...:

  Currently ignored.

## Value

The goodness-of-fit statistics, invisibly: a list with the sample size,
the residual mean direction and resultant length (circular response) or
mean and standard deviation (linear response), the Watson \\U^2\\
statistic and p-value, and the backend goodness-of-fit summary
([`k.check`](https://rdrr.io/pkg/mgcv/man/k.check.html) for
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md);
the higher-order harmonic test, convergence flag, or least-squares fit
metrics for
[`circ_lm`](https://circstat.github.io/circlss/reference/circ_lm.md)).

## Details

The default panel set is `c("rose", "obsfit", "residcov", "qq.unif")`
for a circular response and drops `"rose"` for a linear response (a rose
diagram needs an angular residual). The headline calibration check is
`"qq.unif"`: the probability-integral-transform residuals are uniform
under a correct fit regardless of how the concentration varies, and the
Watson \\U^2\\ test (rotation-invariant, unlike Kolmogorov-Smirnov)
quantifies the departure. The rose diagram should show one tight mode at
zero; an off-centre mode signals location bias, a multimodal one signals
missed structure (too few harmonics, or a smoothing basis with too small
a dimension). For a circular response `"obsfit"` carries the wrapped
diagonal and its \\\pm 2\pi\\ copies, so calibration is read across the
branch cut.

The deviance-residual panels (`"qq.norm"`, `"scaleloc"`, `"hist"`) are
opt-in via `which` (or `which = "all"`). They exploit the deviance
residual being constructed \\\approx N(0, 1)\\: `"scaleloc"` in
particular is the concentration-adequacy check with no linear-model
analogue – a trend in \\\sqrt{\|\\\mathrm{deviance\\ residual}\\\|}\\
against the fitted location means the dispersion model is wrong.

The influence panel `"cook"` is available for
[`circ_lm`](https://circstat.github.io/circlss/reference/circ_lm.md)
fits (the closed-form IRLS / least-squares leverage); for a
general-family
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
mgcv exposes no per-observation leverage, so the panel is dropped with a
message – read basis influence from the
[`k.check`](https://rdrr.io/pkg/mgcv/man/k.check.html) table and the
effective degrees of freedom instead.

For the full base set of linear-response panels on a `"lc"` fit, call
`plot(fit$lm)` directly.

## References

Watson, G. S. (1961) Goodness-of-fit tests on a circle. *Biometrika* 48,
109-114.

Wood, S. N. (2017) *Generalized Additive Models: An Introduction with
R*. Chapman and Hall/CRC, second edition.

## See also

[`circ_resid`](https://circstat.github.io/circlss/reference/circ_resid.md)
for the underlying residuals;
[`plot.circ_gam`](https://circstat.github.io/circlss/reference/plot.circ_gam.md),
[`plot.circ_lm`](https://circstat.github.io/circlss/reference/plot.circ_lm.md)
for the effect displays;
[`gam.check`](https://rdrr.io/pkg/mgcv/man/gam.check.html).

## Examples

``` r
set.seed(1)
n <- 120
x <- rnorm(n)
theta <- (1 + 2 * atan(1.5 * x) + rnorm(n) / 4) %% (2 * pi)
m <- circ_lm(theta ~ x, data.frame(theta, x), type = "cl")
circ_check(m)
```
