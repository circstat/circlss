# Tests for circ_gam(). circ_gam supplies
# defaults (cyclic-knot bracket, trailing ~ 1 fill, method = "REML") and forwards
# everything else to mgcv::gam() unchanged. These exercise that wiring at the
# input level (fit = FALSE, no model fit), plus the named response-scale output
# and the radian utilities.

test_that("cyclic-knot default brackets the radian period", {
  set.seed(1); n <- 200
  phi <- runif(n, -pi, pi)
  d <- data.frame(y = 2 * atan(sin(phi)), phi = phi)
  f <- list(y ~ s(phi, bs = "cc", k = 10), ~ s(phi, bs = "cc", k = 10))
  # covariate takes negative values -> [-pi, pi]
  expect_equal(circlss:::.circ_cc_knots(f, d, NULL), list(phi = c(-pi, pi)))
  # strictly non-negative covariate -> [0, 2*pi]
  d2 <- data.frame(y = d$y, phi = runif(n, 0.2, 2 * pi - 0.2))
  expect_equal(circlss:::.circ_cc_knots(f, d2, NULL), list(phi = c(0, 2 * pi)))
  # user-supplied knots win
  expect_equal(circlss:::.circ_cc_knots(f, d, list(phi = c(-3, 3)))$phi,
               c(-3, 3))
})

test_that("trailing ~ 1 fill builds the same model as an explicit intercept LP", {
  set.seed(1); n <- 150; x <- runif(n)
  d <- data.frame(y = 2 * atan(-0.3 + 1.1 * x) + rnorm(n) / 3, x = x)
  a <- circ_gam(y ~ s(x), data = d, family = vmlss(), fit = FALSE)            # filled
  b <- circ_gam(list(y ~ s(x), ~ 1), data = d, family = vmlss(), fit = FALSE) # explicit
  expect_equal(ncol(a$X), ncol(b$X))
  expect_equal(unname(a$X), unname(b$X))
})

test_that("predict(type = 'response') columns are named by the family", {
  set.seed(1); n <- 150; x <- runif(n)
  d <- data.frame(y = 2 * atan(-0.3 + 1.1 * x) + rnorm(n) / 3, x = x)
  b <- circ_gam(list(y ~ s(x), ~ s(x)), data = d, family = vmlss())
  pr <- predict(b, newdata = data.frame(x = c(0, .5, 1)), type = "response")
  expect_identical(colnames(pr), c("mu", "kappa"))
  prl <- predict(b, newdata = data.frame(x = 0.5), type = "link")  # passed through
  expect_true(is.null(colnames(prl)) || all(colnames(prl) == ""))
})

test_that("input checks: off-period covariate and over-long formula error", {
  set.seed(1); n <- 100; phi <- runif(n, -pi, pi)
  d <- data.frame(y = 2 * atan(sin(phi)) + rnorm(n) / 3, phi = phi + 10)
  expect_error(circ_gam(list(y ~ s(phi, bs = "cc"), ~ 1), data = d,
                        family = vmlss()), "radian")
  d2 <- data.frame(y = rnorm(50), x = runif(50))
  expect_error(circ_gam(list(y ~ x, ~ x, ~ 1), data = d2, family = vmlss()),
               "2 parameters")
})

test_that("radian utils convert and wrap correctly", {
  expect_equal(rad(c(0, 90, 180, 270)), c(0, pi / 2, pi, 3 * pi / 2))
  expect_equal(rad(c(0, 6, 12, 18), 24), c(0, pi / 2, pi, 3 * pi / 2))
  expect_equal(deg(pi), 180)
  expect_equal(wrap(3 * pi), pi)
  expect_equal(wrap(-0.5, "2pi"), 2 * pi - 0.5)
})

test_that("weights are forwarded to mgcv by name, not dropped", {
  set.seed(1); n <- 80; x <- runif(n)
  mu <- 2 * atan(1.2 * sin(2 * pi * x))
  y <- atan2(sin(mu + rnorm(n, sd = 0.5)), cos(mu + rnorm(n, sd = 0.5)))
  d <- data.frame(y = y, x = x)
  set.seed(11); w <- sample(1:3, n, replace = TRUE)
  # mgcv evaluates `weights` by NSE, so circ_gam takes it as an explicit formal
  # and forwards it by name; fit = FALSE resolves it into the model frame without
  # the (nondeterministic) fit.
  G  <- circ_gam(y ~ s(x, k = 6), data = d, family = vmlss(), weights = w,
                 fit = FALSE)
  G0 <- circ_gam(y ~ s(x, k = 6), data = d, family = vmlss(), fit = FALSE)
  expect_equal(as.numeric(G$w), as.numeric(w))   # reached mgcv, not dropped
  expect_true(all(G0$w == 1))                    # default unit weights
})

# ---- center = TRUE machinery (fit-free) -------------------------------------
# The reference chooser, wall-column gate and response rotate-back. .circ_center_ref
# rotates the (weighted) CIRCULAR MEAN to 0, but only when it sits near the wall.
# The hot wall-collapse recovery / round-trip / circ_mix behaviour are validated
# separately by model-fitting scripts, kept out of this fast suite.

test_that(".circ_center_ref is an exact no-op when the wall is already clear", {
  set.seed(1)
  # concentrated near 0: circular mean ~0, far from the wall -> ref 0
  expect_equal(circlss:::.circ_center_ref(wrap(rnorm(200, 0, 0.4))), 0)
  # a wide symmetric fan also has circular mean ~0 -> still a no-op
  expect_equal(circlss:::.circ_center_ref(wrap(seq(-2.0, 2.0, length.out = 200))), 0)
})

test_that(".circ_center_ref moves the wall off data that straddle pi", {
  set.seed(2)
  th  <- wrap(rnorm(200, pi, 0.4))             # mass sitting on the wall
  ref <- circlss:::.circ_center_ref(th)
  expect_gt(abs(ref), 0.5)
  # rotating by ref pulls the circular mean back to ~0 (off the wall)
  expect_lt(abs(atan2(mean(sin(th - ref)), mean(cos(th - ref)))), 0.3)
})

test_that(".circ_center_ref weighted branch centers a component's own mode", {
  set.seed(3)
  th <- wrap(c(rnorm(100, pi, 0.3), rnorm(100, 0, 0.3)))
  # weight only the near-pi mode -> reference ~ pi
  ref <- circlss:::.circ_center_ref(th, weights = c(rep(1, 100), rep(0, 100)))
  expect_lt(abs(wrap(ref - pi)), 0.3)
  # weight only the near-0 mode (clear of the wall) -> no-op
  expect_equal(circlss:::.circ_center_ref(th, weights = c(rep(0, 100), rep(1, 100))), 0)
})

test_that(".circ_wall_loc flags tan-half circular locations only", {
  expect_equal(circlss:::.circ_wall_loc(vmlss()), 1L)
  expect_equal(circlss:::.circ_wall_loc(wclss()), 1L)
  expect_true(is.na(circlss:::.circ_wall_loc(pnlss())))     # derived atan2, no wall
  expect_true(is.na(circlss:::.circ_wall_loc(gausslss())))  # linear response, no wall
})

test_that(".circ_rotate_response shifts only the circular column, wrapped", {
  p   <- cbind(mu = c(-3.0, 0, 3.0), kappa = c(2, 5, 9))
  out <- circlss:::.circ_rotate_response(p, vmlss(), 0.5)
  expect_equal(out[, "mu"], wrap(c(-3.0, 0, 3.0) + 0.5))
  expect_identical(out[, "kappa"], p[, "kappa"])            # scale untouched
  expect_identical(circlss:::.circ_rotate_response(p, vmlss(), 0), p)   # ref 0 = no-op
  expect_identical(circlss:::.circ_rotate_response(p, pnlss(), 0.5), p) # no wall = no-op
})
