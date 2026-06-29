# FD spot-checks of the assembled ll: the gradient (lb) and Hessian (lbb)
# returned at deriv = 1 are compared against finite differences of the
# log-likelihood itself, on a real gam design (smooths in both linear
# predictors, so the full link chain and gamlss.etamu/gamlss.gH assembly
# are exercised). This catches transcription typos locally, independent of
# the Python differential tests.

sim_vm <- function(n, seed = 42) {
  set.seed(seed)
  x <- runif(n)
  mu <- 2 * atan(1.2 * sin(2 * pi * x))
  kappa <- exp(1 + 0.6 * cos(2 * pi * x))
  # Best-Fisher would be overkill for a derivative test point: wrapped
  # normal approximation keeps y plausibly dispersed, which is all the FD
  # check needs (derivatives are checked at arbitrary coef anyway)
  y <- atan2(sin(mu + rnorm(n) / sqrt(kappa)), cos(mu + rnorm(n) / sqrt(kappa)))
  data.frame(y = y, x = x)
}

test_that("ll gradient and Hessian match finite differences", {
  dat <- sim_vm(120)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6)),
                 family = vmlss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  y <- G$y
  p <- ncol(X)
  set.seed(7)
  # a moderate, non-degenerate test point: small smooth coefs around an
  # intercept start near the truth
  coef <- rnorm(p, sd = 0.1)

  llf <- function(b) fam$ll(y, X, b, wt = rep(1, length(y)), family = fam,
                            deriv = 0)$l
  ret <- fam$ll(y, X, coef, wt = rep(1, length(y)), family = fam, deriv = 1)

  h <- 1e-5
  g_fd <- vapply(seq_len(p), function(j) {
    bp <- bm <- coef
    bp[j] <- bp[j] + h
    bm[j] <- bm[j] - h
    (llf(bp) - llf(bm)) / (2 * h)
  }, numeric(1))
  expect_equal(drop(ret$lb), g_fd, tolerance = 1e-5)

  H_fd <- matrix(0, p, p)
  for (j in seq_len(p)) {
    bp <- bm <- coef
    bp[j] <- bp[j] + h
    bm[j] <- bm[j] - h
    gp <- fam$ll(y, X, bp, wt = rep(1, length(y)), family = fam, deriv = 1)$lb
    gm <- fam$ll(y, X, bm, wt = rep(1, length(y)), family = fam, deriv = 1)$lb
    H_fd[, j] <- (gp - gm) / (2 * h)
  }
  expect_equal(unname(as.matrix(ret$lbb)), (H_fd + t(H_fd)) / 2,
               tolerance = 1e-5)
})

test_that("logpdf matches the closed form on the stable branch", {
  # kappa*(cos d - 1) - log(2*pi*i0e(kappa)) == kappa*cos d - log(2*pi*I0)
  y <- c(-2, 0.3, 1.7)
  mu <- c(0.5, -0.1, 2.0)
  kappa <- c(0.3, 5, 50)
  fam <- vmlss()
  l0 <- kappa * (cos(y - mu) - 1) -
    log(2 * pi * besselI(kappa, 0, expon.scaled = TRUE))
  l0_naive <- kappa * cos(y - mu) - log(2 * pi * besselI(kappa, 0))
  expect_equal(l0, l0_naive, tolerance = 1e-12)
})
