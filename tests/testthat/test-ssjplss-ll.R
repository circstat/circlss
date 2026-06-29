# FD spot-checks of ssjplss's assembled ll on a real four-LP gam design (xi via
# tan-half, kappa via log, psi via identity, lambda via the tanh link), plus the
# structural facts: the sine-skew reduces to Jones-Pewsey exactly at lambda = 0,
# the per-obs log-density is the JP log-density plus log(1 + lambda sin(y - xi)),
# and the density still integrates to 1 (the skew factor leaves the normalizer
# unchanged).

# sine-skewed JP sampler by inverse transform on the centered skewed kernel.
rssjp <- function(xi, kappa, psi, lambda) {
  n <- length(xi)
  grid <- seq(-pi, pi, length.out = 2049)
  sg <- sin(grid)
  th <- numeric(n)
  key <- paste(kappa, psi, lambda, sep = "\r")
  for (k in unique(key)) {
    idx <- which(key == k)
    h <- circlss:::.jp_score_terms(grid, kappa[idx[1]], psi[idx[1]],
                                   second = FALSE)$h
    dens <- exp(h - max(h)) * (1 + lambda[idx[1]] * sg)
    dens[dens < 0] <- 0
    cdf <- cumsum((dens[-1] + dens[-length(dens)]) * diff(grid) / 2)
    cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
    th[idx] <- xi[idx] +
      approx(cdf, grid, xout = runif(length(idx)), rule = 2,
             ties = "ordered")$y
  }
  atan2(sin(th), cos(th))
}

test_that("ssjplss ll gradient and Hessian match finite differences", {
  set.seed(51)
  n <- 160
  x <- runif(n)
  xi <- 2 * atan(0.5 * sin(2 * pi * x))
  kappa <- exp(0.9 + 0.4 * cos(2 * pi * x))
  psi <- 0.6
  lambda <- 0.5
  dat <- data.frame(y = rssjp(xi, kappa, rep(psi, n), rep(lambda, n)), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1, ~ 1),
                 family = ssjplss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  p <- ncol(X)
  jj <- attr(X, "lpi")
  set.seed(9)
  coef <- rnorm(p, sd = 0.15)
  coef[jj[[3]]] <- 0.6      # psi intercept
  coef[jj[[4]]] <- 0.55     # lambda intercept (tanh link)

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

test_that("ssjplss reduces to jplss at lambda = 0 (l0 identical)", {
  set.seed(5)
  n <- 120
  x <- runif(n)
  dat <- data.frame(y = runif(n, -pi, pi), x = x)
  Gs <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1, ~ 1),
                  family = ssjplss(), data = dat, fit = FALSE)
  Gj <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1),
                  family = jplss(), data = dat, fit = FALSE)
  js <- attr(Gs$X, "lpi")
  jj <- attr(Gj$X, "lpi")
  set.seed(3)
  cs <- rnorm(ncol(Gs$X), sd = 0.2)
  cs[js[[4]]] <- 0                       # lambda intercept 0 -> lambda = 0
  cj <- numeric(ncol(Gj$X))
  cj[jj[[1]]] <- cs[js[[1]]]
  cj[jj[[2]]] <- cs[js[[2]]]
  cj[jj[[3]]] <- cs[js[[3]]]
  l0s <- Gs$family$ll(Gs$y, Gs$X, cs, rep(1, n), Gs$family, deriv = 0)$l0
  l0j <- Gj$family$ll(Gj$y, Gj$X, cj, rep(1, n), Gj$family, deriv = 0)$l0
  expect_equal(l0s, l0j, tolerance = 1e-12)
})

test_that("ssjplss l0 is the JP log-density plus log(1 + lambda sin)", {
  set.seed(7)
  n <- 100
  x <- runif(n)
  dat <- data.frame(y = runif(n, -pi, pi), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1, ~ 1),
                 family = ssjplss(), data = dat, fit = FALSE)
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
  lambda <- tanh(drop(X[, jj[[4]], drop = FALSE] %*% coef[jj[[4]]]))
  hh <- circlss:::.jp_score_terms(G$y - xi, kappa, psi, second = FALSE)$h
  logc <- circlss:::.jp_logc_vec(kappa, psi)
  ref <- hh + logc + log1p(lambda * sin(G$y - xi))
  expect_equal(l0, ref, tolerance = 1e-12)
})

test_that("ssjplss density integrates to 1 across kappa, psi, lambda", {
  th <- seq(-pi, pi, length.out = 40001)
  for (kappa in c(0.6, 2.0)) {
    for (psi in c(-0.5, 0.5, 1.0)) {
      logc <- circlss:::.jp_quad(kappa, psi)[["logc"]]
      h <- circlss:::.jp_score_terms(th, kappa, psi, second = FALSE)$h
      for (lambda in c(-0.6, 0.0, 0.7)) {
        f <- exp(h + logc) * (1 + lambda * sin(th))
        area <- sum((f[-1] + f[-length(f)]) / 2 * diff(th))
        expect_equal(area, 1, tolerance = 1e-4)
      }
    }
  }
})
