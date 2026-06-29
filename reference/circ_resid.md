# Circular regression residuals

The residual primitive for circlss fits – the circular analogue of
[`residuals`](https://rdrr.io/r/stats/residuals.html), and the quantity
every
[`circ_check`](https://circstat.github.io/circlss/reference/circ_check.md)
panel is a function of. “Observed minus fitted” is undefined when both
are angles, so `circ_resid` returns one of four residual definitions
that are well posed on the circle. It dispatches over both front doors
([`circ_lm`](https://circstat.github.io/circlss/reference/circ_lm.md)
and
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md))
so the residual vocabulary is identical across the parametric and
penalized fits.

## Usage

``` r
circ_resid(
  object,
  type = c("quantile", "deviance", "angular", "pearson"),
  nsim = 1000L,
  scale = c("uniform", "normal")
)
```

## Arguments

- object:

  A fitted
  [`circ_lm`](https://circstat.github.io/circlss/reference/circ_lm.md)
  or
  [`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
  model.

- type:

  The residual definition. `"quantile"` (default) is the
  probability-integral-transform residual, calibrated even when the
  concentration varies across observations; `"deviance"` is the signed
  root of the per-observation deviance, constructed to be approximately
  \\N(0, 1)\\ under a good fit; `"angular"` is the wrapped residual
  \\y - \hat\mu \in (-\pi, \pi\]\\ (the raw response residual, returned
  on the response scale for a linear-response fit); `"pearson"` is the
  score-standardized residual.

- nsim:

  Number of simulation replicates for the `"quantile"` residual when the
  family has no closed-form distribution function (every family except
  the von Mises and Gaussian cases, which are computed analytically).
  Set a seed for a reproducible simulated residual.

- scale:

  For `type = "quantile"` only: `"uniform"` returns the
  probability-integral transform on \\(0, 1)\\ (pairs with the Watson
  \\U^2\\ uniform Q-Q); `"normal"` maps it through
  [`qnorm`](https://rdrr.io/r/stats/Normal.html) to the Dunn–Smyth
  \\N(0, 1)\\ residual (pairs with a normal Q-Q).

## Value

A numeric vector of residuals, one per observation, carrying an
`attr(, "type")` tag (and, for `"quantile"`, an `attr(, "scale")` tag).

## Details

For a circular response the quantile residual is computed in the
*residual frame*: the wrapped residual \\y - \hat\mu\\ is ranked against
its fitted (von Mises/Gaussian) or simulated distribution, with the cut
placed at the antipode \\\pm\pi\\ – the least probable region of a
concentrated residual, which minimises the wrapping artifact. The von
Mises (`"cl"`, `"cc"`,
[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md)) and
Gaussian (`"lc"`,
[`gausslss`](https://circstat.github.io/circlss/reference/gausslss.md))
cases use a closed-form distribution function and are deterministic;
every other family is transformed by simulation from the family's
random-deviate generator, so a seed is needed for reproducibility.

The deviance, angular and Pearson residuals are read straight from each
family's own `residuals` method (for
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md))
or reconstructed from the stored von Mises / least-squares fit (for
[`circ_lm`](https://circstat.github.io/circlss/reference/circ_lm.md));
the centring is always taken from the fitted direction the model
reports, including the derived direction of
[`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md).

## References

Dunn, P. K. and Smyth, G. K. (1996) Randomized quantile residuals.
*Journal of Computational and Graphical Statistics* 5, 236-244.

Fisher, N. I. (1993) *Statistical Analysis of Circular Data*. Cambridge
University Press.

## See also

[`circ_check`](https://circstat.github.io/circlss/reference/circ_check.md)
for the diagnostic panel grid built on these residuals;
[`circ_lm`](https://circstat.github.io/circlss/reference/circ_lm.md),
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md).

## Examples

``` r
set.seed(1)
n <- 80
x <- rnorm(n)
theta <- (1 + 2 * atan(1.5 * x) + rnorm(n) / 4) %% (2 * pi)
m <- circ_lm(theta ~ x, data.frame(theta, x), type = "cl")

r <- circ_resid(m, type = "quantile")   # von Mises: closed-form, deterministic
head(r)
head(circ_resid(m, type = "angular"))    # wrapped y - mu_hat in (-pi, pi]
```
