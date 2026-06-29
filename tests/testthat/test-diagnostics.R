# Tests for the circ_resid()/circ_check() diagnostics. Following the package's
# deterministic-output policy (no model fitting in the suite -- fits are not
# bit-stable across BLAS), these exercise only the pure helpers on hand-built
# inputs: the Watson U2 statistic, the von Mises CDF, the PIT machinery on a
# synthetic "parts" object, the fitted-direction reader, and the rose binning.
# circ_resid()/circ_check() on real fits are smoke-tested offline, not here.

test_that(".circ_watson_u2 matches precomputed values and ranks (in)uniformity", {
  u1 <- c(0.05, 0.15, 0.4, 0.55, 0.6, 0.75, 0.85, 0.95)
  w1 <- circlss:::.circ_watson_u2(u1)
  # closed form for this u1: sum_sq - n*(ubar-1/2)^2 + 1/(12n) = 0.0175 + 1/96
  expect_equal(w1$stat, 0.0175 + 1 / 96)
  expect_equal(w1$p, 0.9457492029, tolerance = 1e-6)

  # a near-perfectly uniform grid -> tiny statistic, p ~ 1
  wg <- circlss:::.circ_watson_u2((seq_len(500) - 0.5) / 500)
  expect_lt(wg$stat, 1e-3)
  expect_gt(wg$p, 0.99)

  # a degenerate cluster -> large statistic, p ~ 0
  wc <- circlss:::.circ_watson_u2(rep(c(0.49, 0.5, 0.51), length.out = 120))
  expect_gt(wc$stat, 1)
  expect_lt(wc$p, 1e-6)

  # the statistic is non-negative and the p-value is a probability
  expect_gte(w1$stat, 0)
  expect_true(w1$p >= 0 && w1$p <= 1)
})

test_that(".circ_pvonmises is a valid von Mises CDF in the residual frame", {
  pvm <- circlss:::.circ_pvonmises
  # F(mu) = 1/2 (the mode bisects the mass), for every concentration
  for (k in c(0.3, 1, 5, 20))
    expect_equal(pvm(0.4, 0.4, k), 0.5, tolerance = 1e-12)
  # spans [0, 1] from the antipode -pi to +pi
  expect_equal(pvm(-pi, 0, 4), 0, tolerance = 1e-12)
  expect_equal(pvm(pi, 0, 4), 1, tolerance = 1e-12)
  # monotone increasing across the branch in the residual frame (mu = 0, so the
  # argument is the residual itself; with mu != 0 the residual wraps at mu +/- pi)
  th <- seq(-pi + 1e-3, pi - 1e-3, length.out = 50)
  expect_true(all(diff(pvm(th, 0, 3)) > 0))
  # kappa = 0 collapses to the circular uniform
  expect_equal(pvm(1.0, 0, 0), (1 + pi) / (2 * pi), tolerance = 1e-12)
  # a precomputed interior value (cross-checked vs circular::pvonmises offline)
  expect_equal(pvm(1.0, 0, 3), 0.9404336090, tolerance = 1e-9)
  # vectorises over theta, mu and kappa together
  expect_length(pvm(c(0, 1, 2), c(0, 0, 0), c(1, 2, 3)), 3L)
})

test_that(".circ_quantile_resid uses the analytic CDF, clamps, and maps to normal", {
  parts <- list(cdf = function() c(0, 0.3, 0.7, 1), simulate = NULL,
                n = 4L, response_circular = TRUE, family_label = "test")
  u <- circlss:::.circ_quantile_resid(parts, 10L, "uniform")
  expect_length(u, 4L)
  expect_true(all(u > 0 & u < 1))            # 0 and 1 clamped off the edges
  expect_equal(u[2:3], c(0.3, 0.7))
  z <- circlss:::.circ_quantile_resid(parts, 10L, "normal")
  expect_equal(z[2:3], stats::qnorm(c(0.3, 0.7)))
})

test_that(".circ_sim_pit ranks observed residuals by the Hazen rule", {
  # n = 2, a fixed simulated matrix; linear response (no wrapping)
  parts <- list(cdf = NULL, n = 2L, response_circular = FALSE,
                fitted_dir = c(0, 0), y = c(0, 1.5), family_label = "test",
                simulate = function(nsim)
                  matrix(c(-2, -1, 1, 2, -2, -1, 1, 2), nrow = 2, byrow = TRUE))
  u <- circlss:::.circ_sim_pit(parts, 4L)
  # row 1: dobs = 0, 2 of 4 sims below  -> (2 + 0.5)/(4 + 1) = 0.5
  # row 2: dobs = 1.5, 3 of 4 sims below -> (3 + 0.5)/(4 + 1) = 0.7
  expect_equal(u, c(0.5, 0.7))
})

test_that("the simulated quantile path is reproducible under a fixed seed", {
  parts <- list(cdf = NULL, n = 3L, response_circular = FALSE,
                fitted_dir = rep(0, 3), y = c(-1, 0, 1), family_label = "test",
                simulate = function(nsim) matrix(stats::rnorm(3 * nsim), 3, nsim))
  set.seed(11); a <- circlss:::.circ_quantile_resid(parts, 200L, "uniform")
  set.seed(11); b <- circlss:::.circ_quantile_resid(parts, 200L, "uniform")
  expect_identical(a, b)
  expect_length(a, 3L)
  expect_true(all(a > 0 & a < 1))
})

test_that(".circ_fitted_dir reads derived, then circular, then mean", {
  fm <- cbind(c(0.1, 0.2), c(2, 3))
  # first circular parameter column (vmlss-like)
  expect_equal(circlss:::.circ_fitted_dir(list(param_circular = c(TRUE, FALSE)), fm),
               c(0.1, 0.2))
  # derived direction wins even when no parameter is flagged circular (pnlss)
  pn <- list(param_circular = c(FALSE, FALSE),
             derived = list(direction = function(f) atan2(f[, 2], f[, 1])))
  expect_equal(circlss:::.circ_fitted_dir(pn, fm), atan2(fm[, 2], fm[, 1]))
  # linear families fall back to the mean column
  expect_equal(circlss:::.circ_fitted_dir(list(param_circular = c(FALSE, FALSE)), fm),
               fm[, 1])
})

test_that(".circ_rose_bins counts equal-width sectors over (-pi, pi]", {
  b <- circlss:::.circ_rose_bins(c(0, 0.01, -0.01, pi - 1e-6, -pi + 1e-6, pi / 2), 4L)
  expect_equal(b$breaks, seq(-pi, pi, length.out = 5L))
  expect_equal(b$counts, c(1L, 1L, 2L, 2L))
  expect_equal(sum(b$counts), 6L)
  # every angle lands in exactly one bin, even after wrapping past the cut
  bb <- circlss:::.circ_rose_bins(c(3 * pi, -3 * pi, 5 * pi / 2), 8L)
  expect_equal(sum(bb$counts), 3L)
})

test_that(".circ_check_keys resolves default/all/explicit and drops rose linearly", {
  ck <- circlss:::.circ_check_keys
  # default set depends on the response type
  expect_equal(ck(NULL, TRUE)$keys, c("rose", "obsfit", "residcov", "qq.unif"))
  expect_equal(ck(NULL, FALSE)$keys, c("obsfit", "residcov", "qq.unif"))
  expect_false(ck(NULL, FALSE)$warn_rose)         # silent drop for the default
  # "all" expands to every panel (rose dropped for a linear response)
  expect_length(ck("all", TRUE)$keys, 8L)
  expect_false("rose" %in% ck("all", FALSE)$keys)
  expect_length(ck("all", FALSE)$keys, 7L)
  expect_true("cook" %in% ck("all", TRUE)$keys)
  # an explicit vector is taken literally
  expect_equal(ck(c("qq.norm", "hist"), TRUE)$keys, c("qq.norm", "hist"))
  # explicitly asking for rose on a linear fit warns (flag) and drops it
  r <- ck(c("rose", "qq.unif"), FALSE)
  expect_true(r$warn_rose)
  expect_equal(r$keys, "qq.unif")
  # unknown keys error
  expect_error(ck("bogus", TRUE), "unknown panel key")
})

test_that("the angular residual wraps straddling pairs to (-pi, pi]", {
  # the diagnostic angular residual is wrap(theta - mu); a pair straddling +/-pi
  # must fold across the cut, not register as a ~2pi error
  expect_equal(wrap(3.1 - (-3.1)), 6.2 - 2 * pi, tolerance = 1e-12)
  expect_true(all(abs(wrap(c(3.1, -3.05) - c(-3.1, 3.05))) < pi))
})
