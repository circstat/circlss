# FD spot-checks of cardlss's assembled ll on a real gam design, with the
# concentration spanning the cardioid's (0, 1/2) range via the logit-half
# link; plus a check that the per-observation log-density matches the
# closed-form cardioid law and that that law integrates to 1.

rcard <- function(n, mu, rho) {
  # exact: rejection from a uniform envelope (accept prob
  # (1 + 2 rho cos phi)/(1 + 2 rho))
  phi <- runif(n, -pi, pi)
  keep <- runif(n) <= (1 + 2 * rho * cos(phi)) / (1 + 2 * rho)
  while (any(!keep)) {
    j <- which(!keep)
    phi[j] <- runif(length(j), -pi, pi)
    keep[j] <- runif(length(j)) <= (1 + 2 * rho[j] * cos(phi[j])) / (1 + 2 * rho[j])
  }
  atan2(sin(mu + phi), cos(mu + phi))
}

test_that("cardlss ll gradient and Hessian match finite differences", {
  set.seed(31)
  n <- 160
  x <- runif(n)
  mu <- 2 * atan(sin(2 * pi * x))
  rho <- 0.5 * plogis(0.4 + 1.6 * cos(2 * pi * x))   # ~[0.06, 0.46], interior
  expect_true(min(rho) > 0 && max(rho) < 0.5)
  dat <- data.frame(y = rcard(n, mu, rho), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6)),
                 family = cardlss(), data = dat, fit = FALSE)
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

test_that("cardlss log-density is the closed-form cardioid and normalizes", {
  set.seed(7)
  n <- 120
  x <- runif(n)
  mu <- 2 * atan(0.8 * sin(2 * pi * x))
  rho <- 0.5 * plogis(0.2 + 1.5 * cos(2 * pi * x))
  dat <- data.frame(y = rcard(n, mu, rho), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6)),
                 family = cardlss(), data = dat, fit = FALSE)
  fam <- G$family
  set.seed(3)
  coef <- rnorm(ncol(G$X), sd = 0.3)
  l0 <- fam$ll(G$y, G$X, coef, rep(1, n), fam, deriv = 0)$l0

  # independent reference: recover (mu, rho) the family used through the
  # links (written out, not via the link objects) and form log f directly
  jj <- attr(G$X, "lpi")
  mm <- 2 * atan(drop(G$X[, jj[[1]]] %*% coef[jj[[1]]]))
  rr <- 0.5 * plogis(drop(G$X[, jj[[2]]] %*% coef[jj[[2]]]))
  ref <- log((1 + 2 * rr * cos(G$y - mm)) / (2 * pi))
  expect_equal(l0, ref, tolerance = 1e-12)

  # the cardioid density integrates to 1 over the circle (normalizer check),
  # including near the rho = 1/2 boundary
  th <- seq(-pi, pi, length.out = 20001)
  for (r in c(0.001, 0.25, 0.499)) {
    f <- (1 + 2 * r * cos(th - 0.4)) / (2 * pi)
    area <- sum((f[-1] + f[-length(f)]) / 2 * diff(th))
    expect_equal(area, 1, tolerance = 1e-6)
  }
})
