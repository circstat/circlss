# FD spot-checks of ajplss's assembled ll on a real four-LP gam design (xi via
# tan-half, kappa via log, psi via identity, nu via the tanh link), plus the
# structural facts: the forward warp reduces to Jones-Pewsey exactly at nu = 0,
# the per-obs log-density is the JP log-kernel at the WARPED angle g = phi + nu
# cos phi plus the asymmetric normalizer, and the density still integrates to 1.
# Unlike ssjplss's sine-skew (which leaves the normalizer untouched), the warp
# MOVES the normalizer -- so the l0 reference is h(g) + .jp_log_c_asym_vec, with
# no log(1 + . sin) factor and no symmetric-JP shortcut.

# asymmetric JP sampler by inverse transform on the warped kernel (mirrors rd()).
rajp <- function(xi, kappa, psi, nu) {
  n <- length(xi)
  grid <- seq(-pi, pi, length.out = 2049)
  th <- numeric(n)
  key <- paste(kappa, psi, nu, sep = "\r")
  for (k in unique(key)) {
    idx <- which(key == k)
    g <- grid + nu[idx[1]] * cos(grid)
    h <- circlss:::.jp_score_terms(g, kappa[idx[1]], psi[idx[1]],
                                   second = FALSE)$h
    dens <- exp(h - max(h))
    dens[dens < 0] <- 0
    cdf <- cumsum((dens[-1] + dens[-length(dens)]) * diff(grid) / 2)
    cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
    th[idx] <- xi[idx] +
      approx(cdf, grid, xout = runif(length(idx)), rule = 2,
             ties = "ordered")$y
  }
  atan2(sin(th), cos(th))
}

test_that("ajplss ll gradient and Hessian match finite differences", {
  set.seed(51)
  n <- 160
  x <- runif(n)
  xi <- 2 * atan(0.5 * sin(2 * pi * x))
  kappa <- exp(0.9 + 0.4 * cos(2 * pi * x))
  psi <- 0.6
  nu <- 0.4
  dat <- data.frame(y = rajp(xi, kappa, rep(psi, n), rep(nu, n)), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1, ~ 1),
                 family = ajplss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  p <- ncol(X)
  jj <- attr(X, "lpi")
  set.seed(9)
  coef <- rnorm(p, sd = 0.15)
  coef[jj[[3]]] <- 0.6      # psi intercept
  coef[jj[[4]]] <- 0.42     # nu intercept (tanh link)

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

test_that("ajplss reduces to jplss at nu = 0 (l0 identical)", {
  set.seed(5)
  n <- 120
  x <- runif(n)
  dat <- data.frame(y = runif(n, -pi, pi), x = x)
  Ga <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1, ~ 1),
                  family = ajplss(), data = dat, fit = FALSE)
  Gj <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1),
                  family = jplss(), data = dat, fit = FALSE)
  ja <- attr(Ga$X, "lpi")
  jj <- attr(Gj$X, "lpi")
  set.seed(3)
  ca <- rnorm(ncol(Ga$X), sd = 0.2)
  ca[ja[[4]]] <- 0                       # nu intercept 0 -> nu = 0
  cj <- numeric(ncol(Gj$X))
  cj[jj[[1]]] <- ca[ja[[1]]]
  cj[jj[[2]]] <- ca[ja[[2]]]
  cj[jj[[3]]] <- ca[ja[[3]]]
  l0a <- Ga$family$ll(Ga$y, Ga$X, ca, rep(1, n), Ga$family, deriv = 0)$l0
  l0j <- Gj$family$ll(Gj$y, Gj$X, cj, rep(1, n), Gj$family, deriv = 0)$l0
  expect_equal(l0a, l0j, tolerance = 1e-12)
})

test_that("ajplss l0 is the JP log-kernel at the warped angle plus the asym normalizer", {
  set.seed(7)
  n <- 100
  x <- runif(n)
  dat <- data.frame(y = runif(n, -pi, pi), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1, ~ 1),
                 family = ajplss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  jj <- attr(X, "lpi")
  set.seed(2)
  coef <- rnorm(ncol(X), sd = 0.2)
  coef[jj[[4]]] <- 0.5
  l0 <- fam$ll(G$y, X, coef, rep(1, n), fam, deriv = 0)$l0
  xi <- 2 * atan(drop(X[, jj[[1]], drop = FALSE] %*% coef[jj[[1]]]))
  kappa <- exp(drop(X[, jj[[2]], drop = FALSE] %*% coef[jj[[2]]]))
  psi <- drop(X[, jj[[3]], drop = FALSE] %*% coef[jj[[3]]])
  nu <- tanh(drop(X[, jj[[4]], drop = FALSE] %*% coef[jj[[4]]]))
  phi <- G$y - xi
  g <- phi + nu * cos(phi)
  hh <- circlss:::.jp_score_terms(g, kappa, psi, second = FALSE)$h
  logc <- circlss:::.jp_log_c_asym_vec(kappa, psi, nu)
  expect_equal(l0, hh + logc, tolerance = 1e-12)
})

test_that("ajplss density integrates to 1 across kappa, psi, nu", {
  # f(y) = exp(h(g(y - xi)) + logc_asym), g(phi) = phi + nu cos phi. With xi = 0
  # the warped kernel times its asym normalizer integrates to 1 over the circle.
  th <- seq(-pi, pi, length.out = 40001)
  for (kappa in c(0.6, 2.0)) {
    for (psi in c(-0.5, 0.5, 1.0)) {
      for (nu in c(-0.6, 0.0, 0.7)) {
        logc <- circlss:::.jp_log_c_asym_vec(kappa, psi, nu)
        g <- th + nu * cos(th)
        h <- circlss:::.jp_score_terms(g, kappa, psi, second = FALSE)$h
        f <- exp(h + logc)
        area <- sum((f[-1] + f[-length(f)]) / 2 * diff(th))
        expect_equal(area, 1, tolerance = 1e-4)
      }
    }
  }
})
