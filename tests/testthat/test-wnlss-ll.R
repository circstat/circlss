# FD spot-checks of wnlss's assembled ll on a real gam design, with the
# concentration spanning BOTH density branches (Fourier rho <= 0.8 and
# wrapped-Gaussian-image rho > 0.8); plus a check that the per-observation
# log-density matches an independent many-image reference on both branches.

rwn <- function(n, mu, rho) {
  sigma <- sqrt(-2 * log(rho))
  th <- mu + sigma * rnorm(n)
  atan2(sin(th), cos(th))
}

test_that("wnlss ll gradient and Hessian match finite differences", {
  set.seed(31)
  n <- 160
  x <- runif(n)
  mu <- 2 * atan(1.2 * sin(2 * pi * x))
  rho <- plogis(0.6 + 1.6 * cos(2 * pi * x))   # ~[0.27, 0.90]: both branches
  expect_true(min(rho) < 0.8 && max(rho) > 0.8)
  dat <- data.frame(y = rwn(n, mu, rho), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6)),
                 family = wnlss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  p <- ncol(X)
  set.seed(9)
  coef <- rnorm(p, sd = 0.2)

  llf <- function(b) fam$ll(G$y, X, b, rep(1, n), fam, deriv = 0)$l
  ret <- fam$ll(G$y, X, coef, rep(1, n), fam, deriv = 1)

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
    gp <- fam$ll(G$y, X, bp, rep(1, n), fam, deriv = 1)$lb
    gm <- fam$ll(G$y, X, bm, rep(1, n), fam, deriv = 1)$lb
    H_fd[, j] <- (gp - gm) / (2 * h)
  }
  expect_equal(unname(as.matrix(ret$lbb)), (H_fd + t(H_fd)) / 2,
               tolerance = 1e-5)
})

test_that("wnlss log-density matches a many-image reference on both branches", {
  # independent reference: the wrapped normal is sum_k N(d + 2*pi*k; 0, sigma)
  wn_ref <- function(d, rho) {
    sigma <- sqrt(-2 * log(rho))
    k <- -60:60
    vapply(seq_along(d), function(i)
      log(sum(dnorm(d[i] + 2 * pi * k, 0, sigma[i]))), numeric(1))
  }
  d <- c(-2.5, -0.7, 0.0, 0.4, 1.9, 3.0)
  for (rho in c(0.3, 0.6, 0.79, 0.81, 0.9, 0.97)) {  # straddle the 0.8 split
    l0 <- circlss:::.wn_terms(d, rep(rho, length(d)))$l0
    expect_equal(l0, wn_ref(d, rep(rho, length(d))), tolerance = 1e-9)
  }
})
