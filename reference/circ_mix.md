# Finite mixtures of circular distributional GAMs by EM

Fits a \\K\\-component finite mixture of circular distributional GAMs by
the EM algorithm. It does not touch the families or mgcv internals: each
M-step is a weighted
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
fit and each E-step reads the family's per-observation density. Because
a component is reached only through that small interface, one engine
spans density clustering, the circular–linear / circular–circular /
linear–circular regression trio, and everything in between – the
response geometry is set entirely by `family`.

## Usage

``` r
circ_mix(
  formula,
  data,
  family = vmlss(),
  K = 2,
  search = c("fixed", "greedy", "grid"),
  assign = c("soft", "hard"),
  group = NULL,
  weights_model = ~1,
  control = circ_mix.control(),
  ...
)

circ_mix.control(
  lambda = NULL,
  kmin = 1L,
  kmax = 20L,
  penalty = c("auto", "fixed", "scheduled"),
  sp = NULL,
  optimizer = "efs",
  start = NULL,
  sp_every = 5L,
  restarts = 10L,
  init = "kmeans",
  moves = c("split", "merge", "death"),
  tol = 1e-06,
  max_iter = 200L,
  min_size = 5L,
  kappa_cap = Inf,
  wfloor = 1e-08,
  seed = NULL,
  cores = 1L,
  verbose = FALSE
)

# S3 method for class 'circ_mix'
print(x, ...)

# S3 method for class 'circ_mix'
summary(object, ...)

# S3 method for class 'circ_mix'
logLik(object, ...)

# S3 method for class 'circ_mix'
coef(object, ...)

# S3 method for class 'circ_mix'
predict(
  object,
  newdata,
  type = c("cluster", "density", "response"),
  log = FALSE,
  ...
)
```

## Arguments

- formula:

  A model formula for one component, exactly as
  [`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
  expects: `theta ~ 1` (density clustering), `theta ~ x1 + x2`
  (circular–linear regression), `theta ~ cos(phi) + sin(phi)` or
  `theta ~ s(phi, bs = "cc")` (circular–circular), or
  `y ~ s(phi, bs = "cc")` with a linear-response family
  (linear–circular). A *list* of formulas carrying two or more distinct
  responses fits a joint (torus) density by the chain rule – e.g.
  `list(psi ~ cos(phi) + sin(phi), phi ~ 1)` factorises \\f(\psi,\phi) =
  f(\psi \mid \phi)\\f(\phi)\\ into two factors; a response named in
  another formula's right-hand side is conditioned on it, and the list
  order is the chain-rule order. (A single response with two or more
  *location-scale* predictors is still one component, written
  `list(theta ~ s(x), ~ s(x))`; the joint reading needs two or more
  distinct left-hand sides. Joint densities over more than two responses
  are not yet supported.)

- data:

  A data frame holding the response and covariates.

- family:

  A circlss location-scale family object (one carrying `param_names`):
  any circular family
  ([`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
  [`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md), ...)
  for a circular response, or
  [`gausslss`](https://circstat.github.io/circlss/reference/gausslss.md)
  /
  [`gammalss`](https://circstat.github.io/circlss/reference/gammalss.md)
  for the linear-circular leg. The family is the only thing that sets
  the response geometry; the EM machinery is identical across all of
  them.

- K:

  The number of mixture components. Under `search = "fixed"` it is held;
  under `search = "greedy"` it is the *starting* count from which the
  data grows or shrinks \\K\\ (not a ceiling or floor); it is ignored
  under `search = "grid"` (which sweeps `control$kmin:control$kmax`).

- search:

  How the number of components is decided. `"fixed"` (default) holds
  `K`. Automatic-K search is opt-in: `"greedy"` runs bidirectional split
  / merge / death moves from the init `K`, accepting any move that
  lowers the penalised objective \\J = -2\\\mathrm{logLik} +
  \lambda\\\mathrm{df}\\ (\\= \\ BIC when `lambda = log(n)`) – the warm
  heuristic, which grows reliably from a small init `K`; `"grid"` fits
  every `K` in `control$kmin:control$kmax` with restarts and picks the
  minimum-\\J\\ `K` – the robust, embarrassingly-parallel selector and
  the cross-check on the moves.

- assign:

  The E-step assignment rule. `"soft"` (default) is EM with fractional
  responsibilities. `"hard"` is classification EM (CEM): each unit seats
  wholly at its argmax component, and the engine maximises the
  classification log-likelihood \\\sum_u \max_k (\log\pi_k + \log
  f_k(y_u))\\ rather than the mixture log-likelihood (its `logLik`/BIC
  are on that classification scale). Combined with `search = "greedy"`
  it is the circular DP-means / k-means-style hard clustering.

- group:

  The clustering unit. `NULL` (default) clusters *rows* – one
  responsibility per observation. A one-sided formula naming a grouping
  variable, `group = ~ id`, clusters *subjects / curves*: a subject's
  whole trajectory seats at one component (the longitudinal /
  latent-class-growth case). Under a group the responsibilities, MAP
  labels and the BIC sample size are all per subject.

- weights_model:

  Reserved for covariate-driven mixing weights (not yet implemented);
  must be `~ 1` (constant mixing proportions).

- control:

  A list of tuning parameters from `circ_mix.control()`.

- ...:

  Further arguments forwarded to the per-component
  [`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
  / [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) M-step (`knots`,
  `method`, ...).

- lambda:

  The penalty multiplier on the degrees of freedom in the model
  objective \\J\\; `NULL` uses `log(n_units)`, making \\J\\ the BIC.

- kmin, kmax:

  The lower and upper bounds on \\K\\ for automatic-K search (both
  `"greedy"` and `"grid"`).

- penalty:

  How the per-component M-step handles the smoothing parameters of
  penalised (smooth) terms. `"auto"` (default) lets REML select them
  every M-step, so each component gets its own automatically-chosen
  smoothness – the usual GAM behaviour; the trade-off is that the moving
  penalty makes the EM non-monotone for smooth components. `"fixed"`
  selects a single smoothness *once*, from a pooled single-component
  pilot fit on all the data, and holds it for every component and
  iteration – an opt-in for speed (no per-M-step REML search), a
  monotone EM, or robustness when `"auto"`'s per-component adaptivity
  lets a component over-flex and absorb a neighbouring cluster; the cost
  is one shared smoothness rather than a per-component one.
  `"scheduled"` starts from that same pooled value and re-selects every
  `sp_every` iterations *after* the first. For parametric (penalty-free)
  components the three modes coincide. Under `"auto"` / `"scheduled"`
  each per-M-step REML search is *warm-started* from the previous
  M-step's selected smoothing parameters (mgcv's `in.out`), so it
  converges in a step or two without changing the value it converges to
  – the per-component, per-iteration smoothness is unchanged, only
  reached faster.

- sp:

  Optional smoothing parameters to hold fixed (a numeric vector for a
  single-response component, a per-factor list for a joint one), in
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) order. When supplied it
  overrides the pilot fit under `penalty = "fixed"`; `NULL` (default)
  selects them automatically.

- optimizer:

  The mgcv outer optimiser for that REML smoothing-parameter search
  (ignored by `penalty = "fixed"` and by parametric components, which
  run no search). `"efs"` (default) is the extended Fellner–Schall
  method (Wood & Fasiolo, 2017): for the families fitted by outer Newton
  (`available.derivs = 2`, e.g.
  [`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
  [`pnlss`](https://circstat.github.io/circlss/reference/pnlss.md)) it
  is markedly faster here and avoids the flat-REML step-failure warnings
  while selecting the same smoothness; for the `available.derivs = 0`
  families it is already mgcv's default. Pass `"outer"` for outer
  Newton.

- start, kappa_cap:

  Accepted but not yet used (kappa_cap is reserved for concentration
  capping).

- sp_every:

  Under `penalty = "scheduled"`, the number of EM iterations between
  REML re-selections of the smoothing parameters (held fixed in
  between).

- restarts:

  Number of random-restart EM runs; the largest-log-likelihood run is
  kept.

- init:

  Initialisation of the responsibilities: `"kmeans"` (default; restart 1
  seeds by clustering the response – circular k-means
  ([`circ_kmeans`](https://circstat.github.io/circlss/reference/circ_kmeans.md))
  on the circle / torus for an angular response, ordinary
  [`kmeans`](https://rdrr.io/r/stats/kmeans.html) for the linear `l~c`
  leg – and later restarts are random) or `"random"`. `"emEM"` is not
  yet supported.

- moves:

  The structure moves the greedy search may attempt, any of `"split"`
  (grow), `"merge"` and `"death"` (shrink); `"birth"` (grow, redundant
  with split) is also accepted.

- tol, max_iter:

  EM convergence tolerance (relative change in the log-likelihood) and
  the maximum number of EM iterations per run.

- min_size:

  The soft-size floor \\n_k = \sum_i \gamma\_{ik}\\ below which a
  component is dropped by a death move, and the minimum members a
  component must have to be split.

- wfloor:

  Lower bound applied to the responsibilities used as M-step prior
  weights, keeping them strictly positive without perturbing the fit.

- seed:

  Optional integer seed set once before the restarts, for
  reproducibility. Each restart then seeds itself from a stream drawn
  here, so the result is identical whether the restarts run serially or
  in parallel.

- cores:

  Number of parallel workers for the embarrassingly-parallel restart
  runs (and each grid-\\K\\'s restarts). `1` (default) runs serially;
  `> 1` forks that many workers via
  [`mclapply`](https://rdrr.io/r/parallel/mclapply.html) (no speed-up on
  Windows, which cannot fork – it falls back to serial). Results do not
  depend on `cores`. A formula with penalised smooths always runs
  serially: a forked large BLAS call can crash a non-fork-safe threaded
  BLAS (e.g. macOS Accelerate) on some builds.

- verbose:

  If `TRUE`, report per-iteration progress.

- object, x:

  A fitted `circ_mix` model.

- newdata:

  A data frame of new data. For `predict`, omitting it uses the training
  frame.

- type:

  For `predict`, what to return: `"cluster"` (the \\n \times K\\
  responsibility matrix, the default), `"density"` (the mixture density
  per row) or `"response"` (per-component response-scale fitted curves,
  a list of length `K`).

- log:

  For `predict(type = "density")`, return the log-density.

## Value

An object of class `"circ_mix"` with, among others:

- components:

  the list of `K` fitted components – each wrapping a weighted
  [`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md),
  or, for a joint density, a *product* of several (one per chain-rule
  factor).

- gating:

  the mixing-weight object; `gating$pi` are the component proportions.

- gamma:

  the \\n \times K\\ matrix of soft responsibilities.

- cluster:

  the per-observation MAP cluster labels.

- loglik, df, bic:

  the mixture log-likelihood, degrees of freedom \\(K-1) + \sum_k
  \mathrm{edf}\_k\\, and BIC.

- converged, iter, ll_path:

  convergence flag, iteration count, and the recorded log-likelihood
  path of the selected run.

- K_init, search:

  the starting component count and the search mode used.

- move_trace:

  for `search = "greedy"`, a data frame of the accepted moves (`move`,
  `K_from`, `K_to`, `J_from`, `J_to`); for `"grid"`, the `K`-by-`K`
  sweep (`loglik`, `df`, `bic`, `J`); `NULL` for `"fixed"`.

- restarts:

  the per-restart log-likelihoods and the basin-hit count.

Methods provided: `print`, `summary`, `predict`, `coef`, `logLik` (so
`AIC`/`BIC` work), and
[`plot.circ_mix`](https://circstat.github.io/circlss/reference/plot.circ_mix.md).

## Details

**EM and restarts.** Each run alternates a weighted M-step (one
`circ_gam` per component, weighted by the responsibilities) with an
E-step that records the observed-data mixture log-likelihood and updates
the responsibilities. For parametric (penalty-free) components the EM is
monotone. The fit is repeated from `restarts` random responsibility
seeds and the largest-log-likelihood run is kept; `$restarts$basin_hits`
reports how many restarts reached it (a health signal). For parametric
components, and for smooth components under `penalty = "fixed"`, the EM
is monotone and a non-monotone step is warned about; under
`penalty = "auto"` a moving smoothing penalty makes small dips expected,
so they are not flagged.

**Model selection.** `df = (K-1) + sum_k edf_k` and
`BIC = -2 logLik + df * log(n)`.
[`logLik()`](https://rdrr.io/r/stats/logLik.html) carries `df` and
`nobs`, so [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) work generically. For a
joint component the per-component `edf` is summed over its factors, so
`df = (K-1) + sum_k sum_j edf_{kj}`.

**Joint (torus) density.** A multi-response `formula` (two distinct
left-hand sides) makes each component a *product* of weighted
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
fits – one per chain-rule factor – whose joint log-density is the sum of
the factor log-densities. The EM loop, restarts, MAP, \\J\\/BIC and the
automatic-K moves are unchanged: the joint case is a
component-implementation swap, not a different engine. The circular
k-means initialisation
([`circ_kmeans`](https://circstat.github.io/circlss/reference/circ_kmeans.md))
seeds on all angular responses jointly – one torus coordinate per
response, e.g. \\(\phi, \psi)\\ – and a greedy split divides the worst
component on its joint angular residuals, so it grows along whichever
response is the more over-dispersed.

**Automatic K.** With `search = "greedy"` (opt-in) `K` is only a
starting count: after EM converges the engine attempts structure moves –
a *split* of the worst-fit component (2-means on its angular residuals;
grows \\K\\), a *merge* of the two most similar components and a *death*
of any below-`min_size` component (both shrink \\K\\) – and accepts the
move that most lowers \\J = -2\\\mathrm{logLik} +
\lambda\\\mathrm{df}\\, repeating until no move improves \\J\\. Because
\\J\\ strictly decreases, the search is monotone and cannot cycle;
`$move_trace` records the accepted moves. The greedy search grows
reliably from a small init `K`; `search = "grid"` (fit every `K` in
`kmin:kmax` with restarts, pick the minimum-\\J\\) is the robust
selector and the recommended cross-check, especially when `K` may be
over-specified or the per-`K` optimum is hard to reach (raise `restarts`
there).

**Longitudinal / curve clustering.** With `group = ~ id` the unit is a
subject: the E-step sums each component's per-row log-densities within
subject (so the whole trajectory shares one responsibility), the M-step
broadcasts that responsibility back to the subject's rows in the
weighted
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md),
and [`logLik`](https://rdrr.io/r/stats/logLik.html) / BIC count \\n =\\
the number of subjects. The trajectory model is just the per-component
formula (e.g. `theta ~ cos(t) + sin(t)`); this is the circular
latent-class-growth (LCGA) case, and it composes with automatic K (moves
reassign whole curves).

**Penalised smooths and hard assignment.** Smooth component terms
(`s(x)`, `s(phi, bs = "cc")`, factor-smooth `s(t, id, bs = "fs")` random
curves) fit through the same weighted M-step; the only subtlety is the
smoothing parameters, handled by `control$penalty` – `"auto"` (default)
selects them by REML per component (automatic smoothness, the GAM way;
the per-iteration search is warm-started across iterations and driven by
the extended Fellner–Schall optimiser by default, so it stays cheap –
see `optimizer`), `"fixed"` freezes one pooled value (faster and
monotone, and a guard against a component over-flexing to absorb another
cluster). Their degrees of freedom enter the BIC as the *effective*
`sum(edf)`, which mgcv reports at the fixed or selected smoothing
parameters, so the criterion stays well defined under regularisation.
`assign = "hard"` switches the E-step to classification EM (CEM); with
`group = ~ id` it clusters whole curves by hard assignment, and with
`search = "greedy"` it is the circular DP-means.

**Unsupported arguments.** An argument value that is not yet implemented
(a non-trivial `weights_model`) is validated and raises an informative
error rather than being silently ignored.

## See also

[`plot.circ_mix`](https://circstat.github.io/circlss/reference/plot.circ_mix.md),
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md),
[`vmlss`](https://circstat.github.io/circlss/reference/vmlss.md),
[`gausslss`](https://circstat.github.io/circlss/reference/gausslss.md)

## Examples

``` r
library(mgcv)
set.seed(1); n <- 400
z  <- sample.int(2L, n, replace = TRUE, prob = c(0.4, 0.6))
x  <- runif(n, -1, 1)
mu <- 2 * atan(c(0.9, -0.9)[z] + c(2.2, -2.2)[z] * x)
y  <- vmlss()$rd(cbind(mu, rep(6, n)), rep(1, n), 1)
# \donttest{
## a two-component von Mises regression on a linear covariate (c~l); fixed K (default)
m  <- circ_mix(y ~ x, data = data.frame(y, x), family = vmlss(), K = 2)
m
table(truth = z, cluster = m$cluster)
plot(m)

## automatic K (opt-in): grow/shrink from the init K by greedy moves
auto <- circ_mix(y ~ x, data = data.frame(y, x), family = vmlss(), K = 2,
                 search = "greedy")
auto$K          # the selected number of components
auto$move_trace # the accepted split/merge/death moves

## the robust cross-check: a brute K = 1..6 BIC sweep
g <- circ_mix(y ~ x, data = data.frame(y, x), family = vmlss(),
              search = "grid", control = circ_mix.control(kmax = 6))

## a joint (torus) density over two angles: each component is a PRODUCT of a
## conditional f(psi | phi) and a marginal f(phi) -- compact 2-D blobs.
set.seed(2); n2 <- 300L
zz  <- sample.int(2L, n2, replace = TRUE)
phi <- vmlss()$rd(cbind(c(-2, 1)[zz], rep(4, n2)), rep(1, n2), 1)
psi <- vmlss()$rd(cbind(c(1, -1.5)[zz] + 0.8 * sin(phi), rep(5, n2)), rep(1, n2), 1)
j <- circ_mix(list(psi ~ cos(phi) + sin(phi), phi ~ 1),
              data = data.frame(psi, phi), family = vmlss(), K = 2)
j               # a "joint torus density" mixture; each component is a product
plot(j)         # the torus-square scatter, coloured by MAP cluster

## longitudinal: cluster whole CURVES, not rows (group = ~id). Two latent classes
## of circular growth trajectory; each subject seats at one component.
set.seed(3); nsub <- 40L; nt <- 8L
cl   <- sample.int(2L, nsub, replace = TRUE)
long <- do.call(rbind, lapply(seq_len(nsub), function(s) {
  tt <- sort(runif(nt)); ph <- 2 * pi * tt; mu <- c(1.2, -1.2)[cl[s]] * sin(ph)
  data.frame(id = s, phase = ph, a = vmlss()$rd(cbind(mu, rep(10, nt)), rep(1, nt), 1))
}))
lc <- circ_mix(a ~ cos(phase) + sin(phase), data = long,
               family = vmlss(), K = 2, group = ~ id)
lc                                       # K components over 40 subjects
table(truth = cl, cluster = lc$cluster)  # cluster is per subject
# }
```
