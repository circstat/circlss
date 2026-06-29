# FD spot-checks of vmftlss's assembled ll on a real three-LP gam design (mu via
# tan-half, kappa via log, nu via tanh), plus the structural story of the
# flat-topped von Mises: the peakedness nu = 0 recovers the von Mises exactly,
# and the trapezoid normalizer makes the density integrate to 1. The peakedness
# warps the angle forward (B = phi + nu sin phi), keeping the law symmetric.

# flat-topped vM sampler by inverse transform on the centered kernel (test-only).
rvmft <- function(mu, kappa, nu) {
  n <- length(mu)
  grid <- seq(-pi, pi, length.out = 2049)
  sg <- sin(grid)
  th <- numeric(n)
  key <- paste(kappa, nu, sep = "\r")
  for (k in unique(key)) {
    idx <- which(key == k)
    ka <- kappa[idx[1]]; nv <- nu[idx[1]]
    B <- grid + nv * sg
    dens <- exp(ka * cos(B) - max(ka * cos(B)))
    cdf <- cumsum((dens[-1] + dens[-length(dens)]) * diff(grid) / 2)
    cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
    u <- runif(length(idx))
    th[idx] <- mu[idx] +
      approx(cdf, grid, xout = u, rule = 2, ties = "ordered")$y
  }
  atan2(sin(th), cos(th))
}

test_that("vmftlss ll gradient and Hessian match finite differences", {
  set.seed(41)
  n <- 160
  x <- runif(n)
  mu <- 2 * atan(0.5 * sin(2 * pi * x))
  kappa <- exp(0.9 + 0.4 * cos(2 * pi * x))   # ~[0.9, 3.7]
  nu <- 0.4
  dat <- data.frame(y = rvmft(mu, kappa, rep(nu, n)), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1),
                 family = vmftlss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  p <- ncol(X)
  set.seed(9)
  coef <- rnorm(p, sd = 0.15)
  coef[attr(X, "lpi")[[3]]] <- atanh(nu)   # a sensible nu intercept (tanh link)

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

test_that("vmftlss peakedness nu = 0 recovers the von Mises exactly", {
  set.seed(5)
  n <- 120
  x <- runif(n)
  dat <- data.frame(y = runif(n, -pi, pi), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1),
                 family = vmftlss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  jj <- attr(X, "lpi")
  set.seed(2)
  coef <- rnorm(ncol(X), sd = 0.2)
  coef[jj[[3]]] <- atanh(0)                 # nu = 0
  mu <- 2 * atan(drop(X[, jj[[1]]] %*% coef[jj[[1]]]))
  kappa <- exp(drop(X[, jj[[2]]] %*% coef[jj[[2]]]))
  vm <- kappa * cos(G$y - mu) -
    log(2 * pi * besselI(kappa, 0, expon.scaled = TRUE)) - kappa
  l0 <- fam$ll(G$y, X, coef, rep(1, n), fam, deriv = 0)$l0
  expect_equal(l0, vm, tolerance = 1e-8)
})

test_that("vmftlss l0 matches vmlss exactly at nu = 0", {
  # cross-family check of the nu = 0 reduction: the assembled vmftlss l0 equals
  # vmlss's at the same (mu, kappa).
  set.seed(8)
  n <- 100
  x <- runif(n)
  dat <- data.frame(y = runif(n, -pi, pi), x = x)
  Gv <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6)),
                  family = vmlss(), data = dat, fit = FALSE)
  Gf <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1),
                  family = vmftlss(), data = dat, fit = FALSE)
  set.seed(3)
  bv <- rnorm(ncol(Gv$X), sd = 0.2)
  jjf <- attr(Gf$X, "lpi")
  bf <- numeric(ncol(Gf$X))
  bf[jjf[[1]]] <- bv[attr(Gv$X, "lpi")[[1]]]
  bf[jjf[[2]]] <- bv[attr(Gv$X, "lpi")[[2]]]
  bf[jjf[[3]]] <- atanh(0)
  l0v <- Gv$family$ll(Gv$y, Gv$X, bv, rep(1, n), Gv$family, deriv = 0)$l0
  l0f <- Gf$family$ll(Gf$y, Gf$X, bf, rep(1, n), Gf$family, deriv = 0)$l0
  expect_equal(l0f, l0v, tolerance = 1e-9)
})

test_that("vmftlss density integrates to 1 across kappa and nu", {
  th <- seq(-pi, pi, length.out = 40001)
  for (kappa in c(0.5, 2.0, 8.0)) {
    for (nu in c(-0.7, -0.3, 0.3, 0.7)) {
      logc <- circlss:::.vmft_log_c_vec(kappa, nu)
      B <- th + nu * sin(th)
      f <- exp(kappa * cos(B) + logc)
      area <- sum((f[-1] + f[-length(f)]) / 2 * diff(th))
      expect_equal(area, 1, tolerance = 1e-4)
    }
  }
})
