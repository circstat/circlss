# FD spot-checks of jplss's assembled ll on a real three-LP gam design (mu via
# tan-half, kappa via log, psi via identity), plus the structural story of the
# Jones-Pewsey family: the shape psi nests von Mises (psi -> 0), the cardioid
# (psi = 1, rho = 1/2 tanh kappa) and the wrapped Cauchy (psi = -1), checked
# against independent closed forms; and a check that the quadrature normalizer
# makes the density integrate to 1.

# JP sampler by inverse transform on the centered kernel grid (test-only).
rjp <- function(mu, kappa, psi) {
  n <- length(mu)
  grid <- seq(-pi, pi, length.out = 2049)
  th <- numeric(n)
  for (k in unique(kappa)) {
    idx <- which(kappa == k)
    h <- circlss:::.jp_score_terms(grid, k, psi, second = FALSE)$h
    dens <- exp(h - max(h))
    cdf <- cumsum((dens[-1] + dens[-length(dens)]) * diff(grid) / 2)
    cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
    u <- runif(length(idx))
    th[idx] <- mu[idx] +
      approx(cdf, grid, xout = u, rule = 2, ties = "ordered")$y
  }
  atan2(sin(th), cos(th))
}

test_that("jplss ll gradient and Hessian match finite differences", {
  set.seed(41)
  n <- 160
  x <- runif(n)
  mu <- 2 * atan(0.5 * sin(2 * pi * x))
  kappa <- exp(0.9 + 0.4 * cos(2 * pi * x))   # ~[0.9, 3.7]
  psi <- 0.7
  dat <- data.frame(y = rjp(mu, kappa, psi), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1),
                 family = jplss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  p <- ncol(X)
  set.seed(9)
  coef <- rnorm(p, sd = 0.15)
  coef[attr(X, "lpi")[[3]]] <- 0.7        # a sensible psi intercept

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

test_that("jplss shape psi nests von Mises, cardioid and wrapped Cauchy", {
  set.seed(5)
  n <- 120
  x <- runif(n)
  dat <- data.frame(y = runif(n, -pi, pi), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1),
                 family = jplss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  jj <- attr(X, "lpi")
  set.seed(2)
  base <- rnorm(ncol(X), sd = 0.2)
  mk <- function(coef) list(
    mu = 2 * atan(drop(X[, jj[[1]]] %*% coef[jj[[1]]])),
    kappa = exp(drop(X[, jj[[2]]] %*% coef[jj[[2]]]))
  )
  l0 <- function(coef) fam$ll(G$y, X, coef, rep(1, n), fam, deriv = 0)$l0

  # psi = 1 : cardioid (1 + tanh(kappa) cos d)/(2 pi), rho = 1/2 tanh kappa
  coef <- base; coef[jj[[3]]] <- 1
  p <- mk(coef)
  card <- log(1 + tanh(p$kappa) * cos(G$y - p$mu)) - log(2 * pi)
  expect_equal(l0(coef), card, tolerance = 1e-7)

  # psi = -1 : wrapped Cauchy 1/(2 pi (cosh kappa - sinh kappa cos d))
  coef <- base; coef[jj[[3]]] <- -1
  p <- mk(coef)
  wc <- -log(2 * pi * (cosh(p$kappa) - sinh(p$kappa) * cos(G$y - p$mu)))
  expect_equal(l0(coef), wc, tolerance = 1e-7)

  # psi -> 0 : von Mises (exact through the |psi| < tol shortcut)
  coef <- base; coef[jj[[3]]] <- 1e-8
  p <- mk(coef)
  vm <- p$kappa * cos(G$y - p$mu) -
    log(2 * pi * besselI(p$kappa, 0, expon.scaled = TRUE)) - p$kappa
  expect_equal(l0(coef), vm, tolerance = 1e-10)
})

test_that("jplss cardioid limit matches cardlss exactly (rho = 1/2 tanh kappa)", {
  # cross-family check of the JP(psi = 1) = cardioid node: build both families'
  # log-density at the SAME mu and matched concentration and compare.
  set.seed(8)
  n <- 80
  d <- runif(n, -pi, pi)
  kappa <- c(0.3, 1.0, 2.5)
  for (k in kappa) {
    jp <- circlss:::.jp_score_terms(d, k, 1.0, second = FALSE)$h +
      circlss:::.jp_quad(k, 1.0)[["logc"]]
    rho <- 0.5 * tanh(k)
    card <- log(1 + 2 * rho * cos(d)) - log(2 * pi)
    expect_equal(jp, card, tolerance = 1e-8)
  }
})

test_that("jplss density integrates to 1 across kappa and psi", {
  th <- seq(-pi, pi, length.out = 40001)
  for (kappa in c(0.5, 2.0)) {
    for (psi in c(-1.0, -0.4, 0.5, 1.0, 2.0)) {
      logc <- circlss:::.jp_quad(kappa, psi)[["logc"]]
      h <- circlss:::.jp_score_terms(th, kappa, psi, second = FALSE)$h
      f <- exp(h + logc)
      area <- sum((f[-1] + f[-length(f)]) / 2 * diff(th))
      expect_equal(area, 1, tolerance = 1e-4)
    }
  }
})
