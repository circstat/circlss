# FD spot-checks of cartlss's assembled ll on a real gam design, with the
# peakedness zeta > 0 via the log link; plus a check that the per-observation
# log-density matches the closed-form Cartwright law (computed in the
# (1 + cos) form, independent of the family's half-angle evaluation) and that
# the Gamma-function normalizer makes the density integrate to 1.

rcart <- function(n, mu, zeta) {
  # exact: Beta-to-angle transform. T ~ Beta(1/2, 1/zeta + 1/2),
  # deviation 2*asin(sqrt(T)), reflected +/- around mu w.p. 1/2
  Tb <- rbeta(n, 0.5, 1 / zeta + 0.5)
  ang <- 2 * asin(sqrt(pmin(Tb, 1)))
  sgn <- sample(c(-1, 1), n, replace = TRUE)
  atan2(sin(mu + sgn * ang), cos(mu + sgn * ang))
}

test_that("cartlss ll gradient and Hessian match finite differences", {
  set.seed(31)
  n <- 160
  x <- runif(n)
  # tightly concentrated, gentle mu: keeps every observation well away from
  # the antipode of the (random) coef, where tan(d/2) blows up and a central
  # difference would itself be inaccurate -- so this is a crisp check of the
  # analytic derivatives, not of finite differencing across the density's zero
  mu <- 2 * atan(0.4 * sin(2 * pi * x))
  zeta <- exp(-1.2 + 0.3 * cos(2 * pi * x))   # ~[0.22, 0.41], rho ~ [0.71, 0.82]
  expect_true(min(zeta) > 0)
  dat <- data.frame(y = rcart(n, mu, zeta), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6)),
                 family = cartlss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  p <- ncol(X)
  set.seed(9)
  coef <- rnorm(p, sd = 0.15)

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

test_that("cartlss log-density is the closed-form Cartwright and normalizes", {
  set.seed(7)
  n <- 120
  x <- runif(n)
  mu <- 2 * atan(0.7 * sin(2 * pi * x))
  zeta <- exp(-0.3 + 0.6 * cos(2 * pi * x))
  dat <- data.frame(y = rcart(n, mu, zeta), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6)),
                 family = cartlss(), data = dat, fit = FALSE)
  fam <- G$family
  set.seed(3)
  coef <- rnorm(ncol(G$X), sd = 0.2)
  l0 <- fam$ll(G$y, G$X, coef, rep(1, n), fam, deriv = 0)$l0

  # independent reference: recover (mu, zeta) the family used through the
  # links and form log f in the (1 + cos) form (not the half-angle form ll uses)
  jj <- attr(G$X, "lpi")
  mm <- 2 * atan(drop(G$X[, jj[[1]]] %*% coef[jj[[1]]]))
  zz <- exp(drop(G$X[, jj[[2]]] %*% coef[jj[[2]]]))
  ref <- (1 / zz - 1) * log(2) + 2 * lgamma(1 + 1 / zz) - log(pi) -
    lgamma(1 + 2 / zz) + (1 / zz) * log(1 + cos(G$y - mm))
  expect_equal(l0, ref, tolerance = 1e-9)

  # the Cartwright density integrates to 1 over the circle (Gamma-normalizer
  # check), across the peakedness range
  cart_dens <- function(th, mu, zeta) {
    C <- 2^(1 / zeta - 1) * gamma(1 + 1 / zeta)^2 / (pi * gamma(1 + 2 / zeta))
    C * (1 + cos(th - mu))^(1 / zeta)
  }
  th <- seq(-pi, pi, length.out = 20001)
  for (z in c(0.3, 1.0, 3.0)) {
    f <- cart_dens(th, 0.4, z)
    area <- sum((f[-1] + f[-length(f)]) / 2 * diff(th))
    expect_equal(area, 1, tolerance = 1e-5)
  }
})
