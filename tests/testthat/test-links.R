# FD checks of the tan-half link derivative chain. Each analytic derivative
# is checked against a central difference of the previous one, so every
# level gets one FD step (accurate to ~1e-9), not a noisy high-order stencil.

fd <- function(f, x, h = 1e-5) (f(x + h) - f(x - h)) / (2 * h)

test_that("tanhalf link derivatives are mutually consistent", {
  lk <- circlss:::tanhalf.link()
  mu <- seq(-2.8, 2.8, length.out = 41)

  # linkinv inverts linkfun
  expect_equal(lk$linkinv(lk$linkfun(mu)), mu, tolerance = 1e-12)

  # mu.eta(eta) = d linkinv / d eta
  eta <- seq(-8, 8, length.out = 41)
  expect_equal(lk$mu.eta(eta), fd(lk$linkinv, eta), tolerance = 1e-8)

  # g'(mu) implied by mu.eta: g' = 1 / mu.eta(g(mu))
  g1 <- function(m) 1 / lk$mu.eta(lk$linkfun(m))
  expect_equal(g1(mu), fd(lk$linkfun, mu), tolerance = 1e-7)

  # d2link = d g'/d mu, d3link = d d2link/d mu, d4link = d d3link/d mu
  expect_equal(lk$d2link(mu), fd(g1, mu), tolerance = 1e-6)
  expect_equal(lk$d3link(mu), fd(lk$d2link, mu), tolerance = 1e-6)
  expect_equal(lk$d4link(mu), fd(lk$d3link, mu), tolerance = 1e-5)
})

test_that("logithalf link derivatives are mutually consistent", {
  lk <- circlss:::logithalf.link()
  rho <- seq(0.05, 0.45, length.out = 41)   # interior of (0, 1/2)

  # linkinv inverts linkfun
  expect_equal(lk$linkinv(lk$linkfun(rho)), rho, tolerance = 1e-12)

  # mu.eta(eta) = d linkinv / d eta (eta range keeps rho strictly inside)
  eta <- seq(-6, 6, length.out = 41)
  expect_equal(lk$mu.eta(eta), fd(lk$linkinv, eta), tolerance = 1e-8)

  # g'(rho) implied by mu.eta: g' = 1 / mu.eta(g(rho))
  g1 <- function(r) 1 / lk$mu.eta(lk$linkfun(r))
  expect_equal(g1(rho), fd(lk$linkfun, rho), tolerance = 1e-7)

  # d2link = d g'/d rho, d3link = d d2link/d rho, d4link = d d3link/d rho
  expect_equal(lk$d2link(rho), fd(g1, rho), tolerance = 1e-6)
  expect_equal(lk$d3link(rho), fd(lk$d2link, rho), tolerance = 1e-6)
  expect_equal(lk$d4link(rho), fd(lk$d3link, rho), tolerance = 1e-5)
})

test_that("tanh link derivatives are mutually consistent", {
  lk <- circlss:::tanh_link()
  lam <- seq(-0.85, 0.85, length.out = 41)   # interior of (-1, 1)

  # linkinv inverts linkfun
  expect_equal(lk$linkinv(lk$linkfun(lam)), lam, tolerance = 1e-12)

  # mu.eta(eta) = d linkinv / d eta (eta range keeps lambda strictly inside)
  eta <- seq(-3, 3, length.out = 41)
  expect_equal(lk$mu.eta(eta), fd(lk$linkinv, eta), tolerance = 1e-8)

  # g'(lambda) implied by mu.eta: g' = 1 / mu.eta(g(lambda))
  g1 <- function(l) 1 / lk$mu.eta(lk$linkfun(l))
  expect_equal(g1(lam), fd(lk$linkfun, lam), tolerance = 1e-7)

  # d2link = d g'/d lambda, d3link = d d2link/d lambda, d4link = d d3link/d lambda
  expect_equal(lk$d2link(lam), fd(g1, lam), tolerance = 1e-6)
  expect_equal(lk$d3link(lam), fd(lk$d2link, lam), tolerance = 1e-5)
  expect_equal(lk$d4link(lam), fd(lk$d3link, lam), tolerance = 1e-4)
})

test_that("A1 derivative chain is FD-consistent, including the series branch", {
  a1 <- circlss:::A1
  d1 <- function(k) circlss:::A1prime(k)
  d2 <- function(k) circlss:::A1prime2(k)
  d3 <- function(k) circlss:::A1prime3(k)
  # d1 down to tiny kappa (step kept below kappa so FD stays in-domain)
  kappa <- c(1e-6, 0.005, 0.009, 0.011, 0.05, 0.5, 2, 10, 100, 700)
  expect_equal(d1(kappa), fd(a1, kappa, 1e-7), tolerance = 1e-5)
  # d2/d3 across the series/recurrence switch at kappa = 0.01
  kappa <- c(0.005, 0.009, 0.011, 0.05, 0.5, 2, 10, 100, 700)
  expect_equal(d2(kappa), fd(d1, kappa, 1e-5), tolerance = 1e-6)
  expect_equal(d3(kappa), fd(d2, kappa, 1e-4), tolerance = 1e-5)
})
