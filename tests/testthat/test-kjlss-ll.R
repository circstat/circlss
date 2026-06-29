# FD spot-checks of kjlss's assembled ll on a real four-LP gam design (mu via
# tan-half, gamma via logit, the disc-chart coordinates u1, u2 via identity),
# plus the structural facts: the disc chart reduces to the wrapped Cauchy
# exactly at u = 0, the density integrates to 1 (the normalizer is exactly 2*pi),
# and the chart keeps the second-order shape (rho, lambda) inside the feasible
# disc for any u.

# Kato-Jones sampler: inverse transform on the centered kernel (mu-invariant
# shape), one grid per unique (gamma, u1, u2).
rkj <- function(mu, gamma, u1, u2) {
  n <- length(mu)
  grid <- seq(-pi, pi, length.out = 4001)
  cg <- cos(grid); sg <- sin(grid)
  th <- numeric(n)
  key <- paste(gamma, u1, u2, sep = "\r")
  for (k in unique(key)) {
    idx <- which(key == k)
    ch <- circlss:::.kj_chart(gamma[idx[1]], u1[idx[1]], u2[idx[1]])
    a <- ch$a; b <- ch$b; g <- gamma[idx[1]]
    D <- pmax(1 + a * a + b * b - 2 * a * cg - 2 * b * sg, 1e-300)
    dens <- pmax(1 + 2 * g * (cg - a) / D, 0)
    cdf <- cumsum((dens[-1] + dens[-length(dens)]) * diff(grid) / 2)
    cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
    th[idx] <- mu[idx] +
      approx(cdf, grid, xout = runif(length(idx)), rule = 2,
             ties = "ordered")$y
  }
  atan2(sin(th), cos(th))
}

test_that("kjlss ll gradient and Hessian match finite differences", {
  set.seed(91)
  n <- 160
  x <- runif(n)
  mu <- 2 * atan(0.5 * sin(2 * pi * x))
  dat <- data.frame(y = rkj(mu, rep(0.6, n), rep(0.4, n), rep(-0.5, n)), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1, ~ 1),
                 family = kjlss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  p <- ncol(X)
  jj <- attr(X, "lpi")
  set.seed(9)
  coef <- rnorm(p, sd = 0.15)
  coef[jj[[3]]] <- 0.4       # u1 intercept
  coef[jj[[4]]] <- -0.5      # u2 intercept

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

test_that("kjlss reduces to wrapped Cauchy at u = 0 (l0 identical)", {
  set.seed(5)
  yy <- runif(200, -pi, pi)
  for (mu in c(-1.5, 0.3, 2.2)) {
    for (gamma in c(0.2, 0.6, 0.85)) {
      l0 <- circlss:::.kj_terms(yy, mu, gamma, 0, 0, deriv = FALSE)$l0
      # WC(mu, rho) with rho = gamma: (1/2pi)(1-rho^2)/(1+rho^2-2rho cos(y-mu))
      wc <- log((1 - gamma^2) /
                  (2 * pi * (1 + gamma^2 - 2 * gamma * cos(yy - mu))))
      expect_equal(l0, wc, tolerance = 1e-12)
    }
  }
})

test_that("kjlss density integrates to 1 across gamma, u1, u2", {
  th <- seq(-pi, pi, length.out = 40001)
  for (gamma in c(0.2, 0.6, 0.9)) {
    for (u1 in c(-1.0, 0.0, 0.8)) {
      for (u2 in c(-0.7, 0.5)) {
        f <- exp(circlss:::.kj_terms(th, 0.4, gamma, u1, u2, deriv = FALSE)$l0)
        area <- sum((f[-1] + f[-length(f)]) / 2 * diff(th))
        expect_equal(area, 1, tolerance = 1e-4)
      }
    }
  }
})

test_that("kjlss disc chart keeps (rho, lam) inside the feasible disc", {
  # the Theorem-1 constraint: (rho cos lam - gamma)^2 + (rho sin lam)^2
  #                            <= (1 - gamma)^2  for every chart point
  set.seed(8)
  for (i in 1:200) {
    gamma <- runif(1, 0.01, 0.99)
    u1 <- rnorm(1, 0, 3); u2 <- rnorm(1, 0, 3)
    ch <- circlss:::.kj_chart(gamma, u1, u2)
    lhs <- (ch$a - gamma)^2 + ch$b^2
    expect_lte(lhs, (1 - gamma)^2 + 1e-12)
  }
})
