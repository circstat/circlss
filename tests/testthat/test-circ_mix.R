# Unit tests for circ_mix -- DETERMINISTIC by design.
# Following the package convention (circ_gam is tested with fit = FALSE, the
# families at the ll level), these exercise the engine's pure logic: the
# dormant-argument guards, formula dispatch, the E-step numerics (log-sum-exp /
# softmax), constant gating, the information-criterion arithmetic, and the S3
# methods on constructed objects. The stochastic end-to-end behaviour (EM
# convergence, label recovery, restarts) is validated separately on simulated
# data, outside the test suite -- never fitted here. The lone
# exception is the E-step density primitive .circ_logpdf (bottom of file): two
# single DETERMINISTIC (parametric, REML) fits whose per-observation $l0 must sum
# to the fitted logLik -- a self-consistency invariant, robust across BLAS.

## a tiny fixed frame; the guards below all fire BEFORE any fit
d <- data.frame(y = c(0.1, -0.2, 0.3, -0.4, 0.5),
                x = c(-1, -0.5, 0, 0.5, 1), id = c(1, 1, 2, 2, 3))

## ---- not-yet-supported / invalid arguments error (guards fire BEFORE fitting) -
## (search = "greedy"/"grid" and assign = "hard" are active and would fit, so
## their recovery is validated separately on simulated data, not here.)
test_that("not-yet-supported and invalid arguments are rejected clearly", {
  expect_error(circ_mix(y ~ x, d, vmlss(), weights_model = ~ x), "not yet supported")
  expect_error(circ_mix(y ~ x, d, family = gaussian()),        "param_names")
  expect_error(circ_mix(y ~ x, d, vmlss(), K = 0),             "K must")
  expect_error(circ_mix(y ~ x, d, vmlss(), K = 99),            "exceeds")
})

## ---- group = ~id: the clustering unit, and its guards (pre-fit) --------------
## (group = ~id is active and would fit, so subject recovery is validated
## separately on simulated data, not here.)
test_that("the group index maps rows to units; bad group errors pre-fit", {
  ## rows (group = NULL): each row its own unit
  ur <- circlss:::.circ_mix_group_index(NULL, d, nrow(d))
  expect_identical(ur$kind, "row")
  expect_identical(ur$grp, seq_len(nrow(d)))
  expect_equal(ur$n_units, nrow(d))
  ## subjects (group = ~id): d$id = c(1,1,2,2,3) -> 3 units, rows mapped 1,1,2,2,3
  us <- circlss:::.circ_mix_group_index(~ id, d, nrow(d))
  expect_identical(us$kind, "subject")
  expect_equal(us$n_units, 3L)
  expect_identical(us$grp, c(1L, 1L, 2L, 2L, 3L))
  expect_identical(us$labels, c("1", "2", "3"))
  ## a missing group variable and K > #subjects error before any fitting
  expect_error(circ_mix(y ~ x, d, vmlss(), group = ~ nope), "not a column")
  expect_error(circ_mix(y ~ x, d, vmlss(), K = 4, group = ~ id), "subjects")
})

## ---- per-row feature aggregates to per-unit means (identity for rows) --------
test_that("the unit feature averages rows within unit", {
  rf  <- matrix(c(1, 3,  2, 4,  10, 10), nrow = 3, byrow = TRUE)  # 3 rows x 2 cols
  ## two units: rows {1,2} -> unit 1 (mean of rows 1,2), row 3 -> unit 2
  uf <- circlss:::.circ_mix_unit_feature(rf, c(1L, 1L, 2L))
  expect_equal(dim(uf), c(2L, 2L))
  expect_equal(unname(uf[1, ]), c(1.5, 3.5))      # (1+2)/2, (3+4)/2
  expect_equal(unname(uf[2, ]), c(10, 10))
  ## grp = seq_len(n) is the identity (the row case)
  expect_equal(unname(circlss:::.circ_mix_unit_feature(rf, 1:3)), unname(rf))
})

## ---- joint guards: d > 2 rejected, missing response caught (pre-fit) ---------
## (a 2-response joint spec is active and would fit, so it is validated
## separately on simulated data, not here.)
test_that("a joint spec over > 2 responses and a missing response error pre-fit", {
  expect_error(circ_mix(list(y ~ 1, x ~ 1, id ~ 1), d, vmlss()), "d > 2")
  expect_error(circ_mix(list(y ~ phi, phi ~ 1), d, vmlss()), "not found in 'data'")
})

## ---- Phase-2 surface: default is "fixed" (auto-K is opt-in), all modes accepted -
test_that("search defaults to fixed and accepts greedy/grid; bad kmin/kmax errors", {
  expect_identical(eval(formals(circ_mix)$search)[1], "fixed")      # the default
  expect_identical(eval(formals(circ_mix)$search),
                   c("fixed", "greedy", "grid"))
  expect_error(circ_mix(y ~ x, d, vmlss(),
                        control = circ_mix.control(kmin = 5, kmax = 2)),
               "kmin")                                              # guard fires pre-fit
})

## ---- the pure responsibility reshapes behind the moves (no fitting) -----------
## Each move must hand the local EM a valid responsibility matrix: row sums stay 1
## and K changes by exactly +-1 (split/birth/merge) or to |keep| (death).
test_that("gamma reshapes preserve row sums and change K correctly", {
  g <- matrix(c(0.7, 0.3,
                0.2, 0.8,
                0.5, 0.5,
                0.9, 0.1), nrow = 4, byrow = TRUE)        # 4 x 2, rows sum to 1

  ## SPLIT column 1 by a per-row label in {1,2}: K -> 3, the two new columns sum
  ## to the old column 1, every other column untouched.
  gs <- circlss:::.circ_mix_gamma_split(g, 1L, c(1L, 2L, 1L, 2L))
  expect_equal(dim(gs), c(4L, 3L))
  expect_equal(unname(rowSums(gs)), rep(1, 4))
  expect_equal(gs[, 1] + gs[, 2], g[, 1])                 # split mass conserved
  expect_equal(gs[, 3], g[, 2])                           # other column carried over
  expect_equal(gs[2, 1], 0)                               # row 2 labelled 2 -> col a empty

  ## MERGE columns 1,2: K -> 1, the merged column is their sum (here all rows 1).
  gm <- circlss:::.circ_mix_gamma_merge(g, 1L, 2L)
  expect_equal(dim(gm), c(4L, 1L))
  expect_equal(unname(gm[, 1]), rep(1, 4))

  ## DEATH keeping only column 2: renormalised to 1.
  gd <- circlss:::.circ_mix_gamma_death(g, keep = 2L)
  expect_equal(dim(gd), c(4L, 1L))
  expect_equal(unname(gd[, 1]), rep(1, 4))

  ## BIRTH seeding rows {1,3}: K -> 3, seeded rows become one-hot on the new col.
  gb <- circlss:::.circ_mix_gamma_birth(g, idx = c(1L, 3L))
  expect_equal(dim(gb), c(4L, 3L))
  expect_equal(unname(rowSums(gb)), rep(1, 4))
  expect_equal(gb[c(1, 3), 3], c(1, 1))                   # seeded rows in the new comp
  expect_equal(gb[c(2, 4), 3], c(0, 0))                   # untouched rows are not
})

## ---- the J objective: -2 logLik + lambda * df, df = (K-1) + sum(edf) ----------
test_that("the search objective J matches its closed form", {
  st <- list(K = 2L, loglik = -100,
             components = list(list(fit = list(edf = c(1, 1))),
                               list(fit = list(edf = c(1, 1)))))
  ## df = (2-1) + (2 + 2) = 5
  expect_equal(circlss:::.circ_mix_objective(st, lambda = log(100)),
               -2 * -100 + log(100) * 5)
  expect_equal(circlss:::.circ_mix_objective(st, lambda = 2),     # AIC-like lambda
               200 + 2 * 5)
})

## ---- circ_mix.control() ------------------------------------------------------
test_that("circ_mix.control returns the committed fields and guards init", {
  ct <- circ_mix.control()
  expect_equal(ct$restarts, 10L)
  expect_equal(ct$wfloor, 1e-8)
  expect_equal(ct$tol, 1e-6)
  expect_true(is.null(ct$lambda))
  expect_identical(ct$init, "kmeans")                       # the default
  expect_identical(circ_mix.control(init = "random")$init, "random")
  expect_error(circ_mix.control(init = "emEM"), "not supported")
  expect_identical(ct$cores, 1L)                            # serial by default
  expect_identical(circ_mix.control(cores = 4)$cores, 4L)
  expect_identical(circ_mix.control(cores = 0)$cores, 1L)   # coerced to >= 1
  ## penalty handling -- "auto" default, validated, sp_every floored at 1
  expect_identical(ct$penalty, "auto")                      # the default
  expect_identical(circ_mix.control(penalty = "fixed")$penalty, "fixed")
  expect_error(circ_mix.control(penalty = "nope"), "should be one of")
  expect_equal(ct$sp_every, 5L)
  expect_identical(circ_mix.control(sp_every = 0)$sp_every, 1L)  # coerced to >= 1
})

## ---- formula dispatch: single-response vs joint (#distinct LHS) --------------
test_that("formula parsing finds the response and counts distinct responses", {
  expect_identical(circlss:::.circ_mix_response(y ~ x), "y")
  expect_identical(circlss:::.circ_mix_response(list(psi ~ phi, phi ~ 1)), "psi")
  expect_equal(circlss:::.circ_mix_n_responses(y ~ s(x)), 1L)
  expect_equal(circlss:::.circ_mix_n_responses(list(th ~ s(x), ~ s(x))), 1L)  # LSS of one resp
  expect_equal(circlss:::.circ_mix_n_responses(list(psi ~ s(phi), phi ~ 1)), 2L)  # joint
  ## all distinct responses, in chain order
  expect_identical(circlss:::.circ_mix_responses(y ~ x), "y")
  expect_identical(circlss:::.circ_mix_responses(list(psi ~ cos(phi), phi ~ 1)),
                   c("psi", "phi"))
  ## nested per-factor LSS list (a factor that is itself a circ_gam LSS spec)
  expect_equal(circlss:::.circ_mix_n_responses(
    list(list(psi ~ s(phi), ~ s(phi)), phi ~ 1)), 2L)
})

## ---- joint spec -> chain-rule factors (the product component shape) ----------
test_that(".circ_mix_factor_specs splits a spec into per-factor circ_gam specs", {
  ## a single formula is one factor (single-response cell)
  fs1 <- circlss:::.circ_mix_factor_specs(y ~ x)
  expect_length(fs1, 1L)
  expect_true(inherits(fs1[[1L]], "formula"))
  ## a list with one response stays one factor (circ_gam's own LSS grammar)
  fs2 <- circlss:::.circ_mix_factor_specs(list(th ~ s(x), ~ s(x)))
  expect_length(fs2, 1L)
  ## a two-response joint spec splits into two bare formulas, in chain order
  fs3 <- circlss:::.circ_mix_factor_specs(list(psi ~ cos(phi) + sin(phi), phi ~ 1))
  expect_length(fs3, 2L)
  expect_identical(circlss:::.circ_mix_response(fs3[[1L]]), "psi")
  expect_identical(circlss:::.circ_mix_response(fs3[[2L]]), "phi")
  ## the nested case: factor 1 keeps its own LSS list, factor 2 is the marginal
  fs4 <- circlss:::.circ_mix_factor_specs(list(list(psi ~ s(phi), ~ s(phi)), phi ~ 1))
  expect_length(fs4, 2L)
  expect_true(is.list(fs4[[1L]]) && !inherits(fs4[[1L]], "formula"))  # LSS list kept
  expect_true(inherits(fs4[[2L]], "formula"))
})

## ---- the product component accessors (constructed, no fitting) ---------------
## A product component holds d factor fits; the four cmp_*() accessors sum / list
## over factors. Build one from mock fits (just the fields the accessors read).
test_that("the product component accessors sum and list over factors", {
  mk <- function(co, edf, smooth = NULL)
    list(coefficients = co, edf = edf, smooth = smooth)
  prod <- circlss:::.circ_mix_product_component(list(
    psi = mk(c(`(Intercept)` = 0.7, b = 0.2), c(1, 1, 1)),   # 3 edf
    phi = mk(c(`(Intercept)` = -1.0), c(1, 1))))             # 2 edf
  expect_true(inherits(prod, "circ_mix_product"))
  expect_true(inherits(prod, "circ_mix_component"))          # IS-A component
  expect_equal(circlss:::cmp_edf(prod), 5)                   # 3 + 2 summed over factors
  cf <- circlss:::cmp_coef(prod)                             # per-factor list (warm-start shape)
  expect_type(cf, "list"); expect_named(cf, c("psi", "phi"))
  expect_equal(unname(cf$psi[["b"]]), 0.2)
  ## smooth detection is product-aware (drives em_core's monotonicity expectation)
  expect_false(circlss:::.circ_mix_has_smooth(prod))
  prod_s <- circlss:::.circ_mix_product_component(list(
    mk(0, 1.7, smooth = list("a cc smooth")), mk(0, c(1, 1))))
  expect_true(circlss:::.circ_mix_has_smooth(prod_s))
  ## a single (non-product) component is detected too
  cp1 <- structure(list(fit = mk(0, c(1, 1), smooth = list("s"))),
                   class = "circ_mix_component")
  expect_true(circlss:::.circ_mix_has_smooth(cp1))
})

## ---- hard (CEM) classification of a per-unit logmix matrix -------------------
## The hard analogue of the soft E-step helpers below: each unit seats wholly at
## its argmax, the loglik is the classification log-likelihood (sum of the maxima).
test_that("the hard E-step seats each unit at its argmax (classification loglik)", {
  logmix <- matrix(c(-1, -2,        # unit 1 -> component 1
                     -0.5, -0.3,    # unit 2 -> component 2
                     -10, -30),     # unit 3 -> component 1
                   nrow = 3, byrow = TRUE)
  cl <- circlss:::.circ_mix_classify(logmix)
  expect_identical(cl$z, c(1L, 2L, 1L))
  expect_equal(unname(rowSums(cl$gamma)), rep(1, 3))     # one-hot rows
  expect_equal(cl$gamma[2, ], c(0, 1))
  expect_equal(cl$loglik, -1 + -0.3 + -10)              # sum of the per-unit maxima
})

## ---- cmp_sp reads the fitted smoothing parameters ---------------------------
## A single component returns the sp vector; a product returns a per-factor list
## (the shape .circ_mix_fit_component()'s `sp` consumes for penalty = "fixed").
test_that("cmp_sp returns the sp vector (single) / per-factor list (product)", {
  cp <- structure(list(fit = list(sp = c(`s(x)` = 2.5))), class = "circ_mix_component")
  expect_equal(unname(circlss:::cmp_sp(cp)), 2.5)
  prod <- circlss:::.circ_mix_product_component(list(
    psi = list(sp = c(`s(phi)` = 1.1)), phi = list(sp = numeric(0))))
  sp <- circlss:::cmp_sp(prod)
  expect_type(sp, "list"); expect_named(sp, c("psi", "phi"))
  expect_equal(unname(sp$psi), 1.1)
  expect_length(sp$phi, 0L)                              # a parametric factor: no sp
})

## ---- E-step numerics: log-sum-exp mixture loglik and softmax responsibilities -
test_that("the E-step row reductions match closed forms and are stable", {
  ## logmix[i,k] = log pi_k + log f_k(y_i); hand-built, no fitting
  logmix <- matrix(c(-1, -2,
                     -0.5, -0.5,
                     -10, -30), nrow = 3, byrow = TRUE)
  rows <- circlss:::.circ_mix_loglik_rows(logmix)
  expect_equal(rows[1], log(exp(-1) + exp(-2)))
  expect_equal(rows[2], log(2 * exp(-0.5)))
  expect_equal(rows[3], -10 + log1p(exp(-20)))            # row-max stabilised
  g <- circlss:::.circ_mix_responsibilities(logmix)
  expect_equal(unname(rowSums(g)), rep(1, 3))
  expect_equal(g[1, ], c(exp(-1), exp(-2)) / (exp(-1) + exp(-2)))
  expect_equal(g[3, 1], 1, tolerance = 1e-8)              # no underflow to NaN
})

## ---- constant gating returns the pi matrix at any newdata --------------------
test_that("the constant gating broadcasts pi over rows", {
  pm <- circlss:::.circ_mix_gate(list(type = "constant", pi = c(0.3, 0.7)),
                                 data.frame(a = 1:4), K = 2L)
  expect_equal(dim(pm), c(4L, 2L))
  expect_equal(pm[1, ], c(0.3, 0.7))
  expect_true(all(pm[, 1] == 0.3))
})

## ---- component edf accessor (the df building block) --------------------------
test_that("cmp_edf sums the fit's edf", {
  cp <- structure(list(fit = list(edf = c(1, 1, 1))), class = "circ_mix_component")
  expect_equal(circlss:::cmp_edf(cp), 3)
})

## ---- IC arithmetic + methods on a constructed object (no fitting) ------------
## a mock component whose $fit carries just enough for coef()/edf
.mk_cp <- function(co) structure(list(fit = list(coefficients = co, edf = rep(1, length(co)))),
                                 class = "circ_mix_component")
.fake_mix <- function() structure(list(
  call = quote(circ_mix(y ~ x)), formula = y ~ x, family = list(family = "vmlss"),
  response = "y", K = 2L, K_init = 2L,
  components = list(.mk_cp(c(`(Intercept)` = 0.5, x = 1.2)),
                   .mk_cp(c(`(Intercept)` = -0.5, x = -1.2))),
  gating = list(type = "constant", pi = c(0.4, 0.6)),
  unit = list(kind = "row", index = 1:100, n_units = 100L),
  gamma = NULL, cluster = NULL, nk = c(40L, 60L), Gtilde = 2L,
  loglik = -210.5, df = 7, edf = c(3, 3), bic = -2 * -210.5 + 7 * log(100),
  objective = list(J = 0, loglik = -210.5, df = 7, bic = 0, lambda = log(100)),
  iter = 12L, converged = TRUE, ll_path = c(-260, -215, -210.5),
  monotone = TRUE, restarts = list(R = 3L, lls = c(-210.5, -211, -212), basin_hits = 2L),
  geometry = "cl", control = circ_mix.control()), class = "circ_mix")

test_that("logLik / AIC / BIC ride the stored loglik, df and unit count", {
  m  <- .fake_mix()
  ll <- logLik(m)
  expect_equal(as.numeric(ll), -210.5)
  expect_equal(attr(ll, "df"), 7)
  expect_equal(attr(ll, "nobs"), 100L)
  expect_equal(as.numeric(BIC(m)), -2 * -210.5 + 7 * log(100))   # generic via logLik
  expect_equal(as.numeric(AIC(m)), -2 * -210.5 + 2 * 7)
})

test_that("coef returns one named vector per component", {
  cf <- coef(.fake_mix())
  expect_named(cf, c("component1", "component2"))
  expect_named(cf[[1]], c("(Intercept)", "x"))
})

test_that("print and summary render the leg, components and fit summary", {
  m <- .fake_mix()
  expect_output(print(m), "Finite mixture")
  expect_output(print(m), "circular-linear")        # the c~l leg label
  expect_output(print(m), "BIC")
  expect_output(summary(m), "per-component coefficients")
})

## ---- the E-step density primitive .circ_logpdf -- the ONLY fits in this file.
## .circ_logpdf wraps a FITTED model (it reads coef + the lpmatrix), so two small
## DETERMINISTIC fits anchor its one contract: the per-observation $l0, summed over
## the training rows, IS the model log-likelihood -- on a circular and a linear
## (l~c) response.
test_that(".circ_logpdf is per-observation and sums to the fitted log-likelihood", {
  set.seed(1); n <- 150; x <- runif(n)
  y <- vmlss()$rd(cbind(2 * atan(1.1 * sin(2 * pi * x)), rep(6, n)), rep(1, n), 1)
  fit <- circ_gam(y ~ x, data = data.frame(y, x), family = vmlss())   # parametric MLE
  lp  <- circlss:::.circ_logpdf(fit)
  expect_length(lp, n)
  expect_true(all(is.finite(lp)))
  expect_equal(sum(lp), as.numeric(logLik(fit)), tolerance = 1e-4)
  expect_equal(lp, circlss:::.circ_logpdf(fit, fit$model))            # default frame
  expect_length(circlss:::.circ_logpdf(fit, data.frame(y = y[1:5], x = x[1:5])), 5)
  expect_error(circlss:::.circ_logpdf(fit, data.frame(x = 0)), "response")
})

test_that(".circ_logpdf works for the linear-response (l~c) leg via gausslss", {
  set.seed(2); n <- 150; phi <- runif(n, -pi, pi)
  y <- 1.5 + 1.2 * sin(phi) + rnorm(n, 0, 0.5)
  fit <- circ_gam(y ~ sin(phi) + cos(phi), data = data.frame(y, phi),
                  family = gausslss())
  expect_equal(sum(circlss:::.circ_logpdf(fit)), as.numeric(logLik(fit)),
               tolerance = 1e-4)
})
