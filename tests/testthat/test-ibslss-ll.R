# FD spot-checks of ibslss's assembled ll on a real four-LP gam design (xi via
# tan-half, kappa via log, nu and lambda via tanh), plus the structural story of
# the inverse Batschelet: nu = lambda = 0 recovers the von Mises exactly, and the
# warp + trapezoid normalizer makes the density integrate to 1. The assembled
# gradient/Hessian wiring (gamlss.etamu / gamlss.gH over four predictors) is what
# the FD check guards. The normalizer kappa/lambda
# derivatives are finite-differenced inside .ibs_terms (as in pycircstat2), so the
# Hessian's kappa/lambda blocks carry an FD-of-FD noise floor and are checked at a
# looser tolerance than the analytic xi/nu blocks.

# Inverse Batschelet sampler (test-only): the FORWARD warps phi(u) are
# closed-form and monotone, so gridding the base angle u gives a (phi, density)
# grid to invert -- no root-finding needed for simulation.
ribs <- function(xi, kappa, nu, lmbd) {
  n <- length(xi)
  u <- seq(-pi, pi, length.out = 2049)
  th <- numeric(n)
  key <- paste(kappa, nu, lmbd, sep = "\r")
  for (k in unique(key)) {
    idx <- which(key == k)
    ka <- kappa[idx[1]]; nv <- nu[idx[1]]; lm <- lmbd[idx[1]]
    phistar <- u - 0.5 * (1 + lm) * sin(u)
    phi <- phistar - nv * (1 + cos(phistar))
    A <- u - 0.5 * (1 - lm) * sin(u)
    dens <- exp(ka * cos(A) - max(ka * cos(A)))
    cdf <- cumsum((dens[-1] + dens[-length(dens)]) / 2 * diff(phi))
    cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
    uu <- runif(length(idx))
    th[idx] <- xi[idx] + approx(cdf, phi, xout = uu, rule = 2,
                                ties = "ordered")$y
  }
  atan2(sin(th), cos(th))
}

test_that("ibslss ll gradient and Hessian match finite differences", {
  set.seed(41)
  n <- 160
  x <- runif(n)
  xi <- 2 * atan(0.5 * sin(2 * pi * x))
  kappa <- exp(0.9 + 0.4 * cos(2 * pi * x))   # ~[0.9, 3.7]
  nu <- 0.4
  lmbd <- 0.3
  dat <- data.frame(y = ribs(xi, kappa, rep(nu, n), rep(lmbd, n)), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1, ~ 1),
                 family = ibslss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  p <- ncol(X)
  jj <- attr(X, "lpi")
  set.seed(9)
  coef <- rnorm(p, sd = 0.15)
  coef[jj[[3]]] <- atanh(nu)      # sensible nu intercept (tanh link)
  coef[jj[[4]]] <- atanh(lmbd)    # sensible lambda intercept (tanh link)

  llf <- function(b) fam$ll(G$y, X, b, rep(1, n), fam, deriv = 0)$l
  ret <- fam$ll(G$y, X, coef, rep(1, n), fam, deriv = 1)

  h <- 1e-5
  g_fd <- vapply(seq_len(p), function(j) {
    bp <- bm <- coef
    bp[j] <- bp[j] + h
    bm[j] <- bm[j] - h
    (llf(bp) - llf(bm)) / (2 * h)
  }, numeric(1))
  # xi/nu blocks are analytic (match tightly); kappa/lambda carry the FD'd
  # normalizer gradient, so the aggregate (~3e-8) sits a touch above machine-FD.
  expect_equal(drop(ret$lb), g_fd, tolerance = 1e-6)

  H_fd <- matrix(0, p, p)
  for (j in seq_len(p)) {
    bp <- bm <- coef
    bp[j] <- bp[j] + h
    bm[j] <- bm[j] - h
    gp <- fam$ll(G$y, X, bp, rep(1, n), fam, deriv = 1)$lb
    gm <- fam$ll(G$y, X, bm, rep(1, n), fam, deriv = 1)$lb
    H_fd[, j] <- (gp - gm) / (2 * h)
  }
  # FD-of-lb of the FD'd-normalizer second derivatives (kappa/lambda blocks) is
  # the noisiest path (~3e-3 here); a wiring bug would be O(1), far above this.
  expect_equal(unname(as.matrix(ret$lbb)), (H_fd + t(H_fd)) / 2,
               tolerance = 1e-2)
})

test_that("ibslss nu = lambda = 0 recovers the von Mises exactly", {
  set.seed(5)
  n <- 120
  x <- runif(n)
  dat <- data.frame(y = runif(n, -pi, pi), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1, ~ 1),
                 family = ibslss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  jj <- attr(X, "lpi")
  set.seed(2)
  coef <- rnorm(ncol(X), sd = 0.2)
  coef[jj[[3]]] <- atanh(0)                 # nu = 0
  coef[jj[[4]]] <- atanh(0)                 # lambda = 0
  xi <- 2 * atan(drop(X[, jj[[1]]] %*% coef[jj[[1]]]))
  kappa <- exp(drop(X[, jj[[2]]] %*% coef[jj[[2]]]))
  vm <- kappa * cos(G$y - xi) -
    log(2 * pi * besselI(kappa, 0, expon.scaled = TRUE)) - kappa
  l0 <- fam$ll(G$y, X, coef, rep(1, n), fam, deriv = 0)$l0
  expect_equal(l0, vm, tolerance = 1e-8)
})

test_that("ibslss l0 matches vmlss exactly at nu = lambda = 0", {
  # cross-family check of the reduction: the assembled ibslss l0 equals vmlss's
  # at the same (xi, kappa).
  set.seed(8)
  n <- 100
  x <- runif(n)
  dat <- data.frame(y = runif(n, -pi, pi), x = x)
  Gv <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6)),
                  family = vmlss(), data = dat, fit = FALSE)
  Gf <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6), ~ 1, ~ 1),
                  family = ibslss(), data = dat, fit = FALSE)
  set.seed(3)
  bv <- rnorm(ncol(Gv$X), sd = 0.2)
  jjf <- attr(Gf$X, "lpi")
  bf <- numeric(ncol(Gf$X))
  bf[jjf[[1]]] <- bv[attr(Gv$X, "lpi")[[1]]]
  bf[jjf[[2]]] <- bv[attr(Gv$X, "lpi")[[2]]]
  bf[jjf[[3]]] <- atanh(0)
  bf[jjf[[4]]] <- atanh(0)
  l0v <- Gv$family$ll(Gv$y, Gv$X, bv, rep(1, n), Gv$family, deriv = 0)$l0
  l0f <- Gf$family$ll(Gf$y, Gf$X, bf, rep(1, n), Gf$family, deriv = 0)$l0
  expect_equal(l0f, l0v, tolerance = 1e-9)
})

test_that("ibslss density integrates to 1 across kappa, nu and lambda", {
  # uses the actual inverse-warp + normalizer code path (.ibs_warp_vec /
  # .ibs_log_c_vec): A(theta) from the two inverse warps, density exp(kappa cos A
  # + log c) integrated on a fine theta grid.
  th <- seq(-pi, pi, length.out = 20001)
  for (kappa in c(0.5, 2.0, 6.0)) {
    for (nu in c(-0.5, 0.4)) {
      for (lmbd in c(-0.4, 0.0, 0.5)) {
        w <- circlss:::.ibs_warp_vec(th, 0, nu, lmbd)
        A <- w$u_star - 0.5 * (1 - lmbd) * sin(w$u_star)
        logc <- circlss:::.ibs_log_c_vec(kappa, lmbd)
        f <- exp(kappa * cos(A) + logc)
        area <- sum((f[-1] + f[-length(f)]) / 2 * diff(th))
        expect_equal(area, 1, tolerance = 1e-3)
      }
    }
  }
})
