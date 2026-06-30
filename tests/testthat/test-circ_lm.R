# Tests for circ_lm() -- the classical (parametric) circular-regression front
# door. Following the package's deterministic-output policy, these exercise the
# BLAS-independent surface only: the A1inv inverse (closed-form Bessel math), the
# input guards (all of which error BEFORE any model is fitted), and the design
# helpers. Convergence/parity against circular is checked offline, not here.

test_that("A1inv inverts A1 to machine precision and handles the edges", {
  ks <- c(0.05, 0.25, 0.7, 1.5, 4, 12, 40)
  expect_equal(circlss:::A1inv(circlss:::A1(ks)), ks, tolerance = 1e-8)
  # the defining property: A1(A1inv(R)) == R across the approximation's seams
  Rs <- c(0.1, 0.52, 0.53, 0.84, 0.85, 0.95)
  expect_equal(circlss:::A1(circlss:::A1inv(Rs)), Rs, tolerance = 1e-10)
  expect_identical(circlss:::A1inv(0), 0)
  expect_identical(circlss:::A1inv(1), Inf)
  expect_identical(circlss:::A1inv(-0.3), 0)   # clamped to [0, 1]
})

test_that("type accepts hyphenated and cased spellings, rejects unknown", {
  expect_identical(circlss:::.circ_lm_type("c-l"), "cl")
  expect_identical(circlss:::.circ_lm_type("C-C"), "cc")
  expect_identical(circlss:::.circ_lm_type(c("lc", "cc")), "lc")  # default vector
  expect_error(circlss:::.circ_lm_type("zz"), "type must be one of")
})

test_that("smooth terms are rejected with a pointer to circ_gam()", {
  d <- data.frame(theta = c(0.1, 1, 2, 3), x = c(0, 1, 2, 3))
  expect_error(circ_lm(theta ~ s(x), d, type = "cl"), "circ_gam")
  expect_error(circ_lm(list(theta ~ x, ~ te(x)), d, type = "cl"), "circ_gam")
  expect_error(circ_lm(y ~ s(x, bs = "cc"), d, type = "lc"), "circ_gam")
})

test_that("cc / lc require exactly one angular predictor", {
  d <- data.frame(y = c(0.1, 1, 2, 3), a = c(0, 1, 2, 3), b = c(1, 1, 2, 2))
  expect_error(circ_lm(y ~ a + b, d, type = "cc"), "exactly one")
  expect_error(circ_lm(y ~ a + b, d, type = "lc"), "exactly one")
  expect_error(circ_lm(y ~ 1, d, type = "lc"), "exactly one")
})

test_that("cl guards: nothing to regress, and mixed needs a shared design", {
  d <- data.frame(theta = c(0.1, 1, 2, 3), x = c(0, 1, 2, 3), z = c(3, 2, 1, 0))
  expect_error(circ_lm(theta ~ 1, d, type = "cl"), "nothing to regress")
  expect_error(circ_lm(list(theta ~ 1, ~ 1), d, type = "cl"), "nothing to regress")
  expect_error(circ_lm(list(theta ~ x, ~ z), d, type = "cl"), "shared design")
  expect_error(circ_lm(list(theta ~ x, ~ x, ~ x), d, type = "cl"), "at most two")
})

test_that("cl init resolves cold / named list / legacy numeric, with guards", {
  X <- matrix(c(1, 2, 3, 4), ncol = 1L)        # p = 1
  st <- function(model, init) circlss:::.circ_lm_cl_starts(X, model, init)

  expect_equal(st("mixed", NULL), list(beta = 0, gamma = 0, alpha = 0))   # cold default

  s <- st("mixed", list(beta = 0.5, alpha = 1.2, gamma = -0.3))           # explicit starts
  expect_equal(c(s$beta, s$gamma, s$alpha), c(0.5, -0.3, 1.2))
  expect_equal(st("mixed", list(gamma = 0.4))$beta, 0)                    # omitted -> cold
  expect_equal(st("mean", 0.7)$beta, 0.7)                                 # legacy numeric -> beta

  expect_error(st("mixed", list(beta = c(1, 2))), "length 1")             # length guards
  expect_error(st("mixed", list(alpha = c(1, 2))), "single number")
  expect_warning(st("mean", list(gamma = 0.1)), "unused")                 # misuse warnings
  expect_warning(st("kappa", list(beta = 0.1)), "unused")
  expect_warning(st("mixed", list(bad = 1)), "unknown init component")
})

test_that(".circ_lm_harmonic builds the cos/sin basis at the requested order", {
  x <- c(0, pi / 2, pi)
  H <- circlss:::.circ_lm_harmonic(x, 2L)
  expect_identical(names(H), c("cos1", "sin1", "cos2", "sin2"))
  expect_equal(H$cos1, cos(x))
  expect_equal(H$sin2, sin(2 * x))
})

test_that(".circ_lm_covariate finds the single plotting covariate, else NULL", {
  # cc / lc: the one angular predictor is stored directly
  expect_identical(circlss:::.circ_lm_covariate(list(type = "cc", var = "phi")), "phi")
  expect_identical(circlss:::.circ_lm_covariate(list(type = "lc", var = "phi")), "phi")
  # cl pools the mean and log-kappa formulas
  cl <- function(mu, ka) list(type = "cl", mu_formula = mu, kappa_formula = ka)
  expect_identical(circlss:::.circ_lm_covariate(cl(theta ~ x, ~ 1)), "x")
  expect_identical(circlss:::.circ_lm_covariate(cl(theta ~ 1, ~ x)), "x")
  expect_identical(circlss:::.circ_lm_covariate(cl(theta ~ x, ~ x)), "x")
  # several covariates have no single axis -> NULL (the geometry-view fallback)
  expect_null(circlss:::.circ_lm_covariate(cl(theta ~ x + z, ~ 1)))
  expect_null(circlss:::.circ_lm_covariate(cl(theta ~ x, ~ z)))
})

test_that(".circ_lm_design drops the intercept and names covariates", {
  d <- data.frame(x = c(1, 2, 3), z = c(4, 5, 6))
  X <- circlss:::.circ_lm_design(theta ~ x + z, d)
  expect_identical(colnames(X), c("x", "z"))
  expect_false("(Intercept)" %in% colnames(X))
  # an intercept-only RHS yields a zero-column design
  expect_identical(ncol(circlss:::.circ_lm_design(~ 1, d)), 0L)
})
