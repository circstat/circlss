# FD spot-checks of wclss's assembled ll on a real gam design, plus the
# stable-branch identity of the log-density.

test_that("wclss ll gradient and Hessian match finite differences", {
  rwc <- function(n, mu, rho) {
    th <- mu + (-log(rho)) * rcauchy(n)
    atan2(sin(th), cos(th))
  }
  set.seed(21)
  n <- 150
  x <- runif(n)
  mu <- 2 * atan(1.4 * sin(2 * pi * x))
  rho <- plogis(0.5 + 0.8 * cos(2 * pi * x))
  dat <- data.frame(y = rwc(n, mu, rho), x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6)),
                 family = wclss(), data = dat, fit = FALSE)
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

test_that("wclss log-density matches the textbook form away from rho = 1", {
  y <- c(-2, 0.3, 1.7)
  mu <- c(0.5, -0.1, 2.0)
  rho <- c(0.2, 0.6, 0.9)
  d <- y - mu
  stable <- log1p(-rho) + log1p(rho) - log(2 * pi) -
    log((1 - rho)^2 + 4 * rho * sin(0.5 * d)^2)
  naive <- log((1 - rho^2) / (2 * pi * (1 + rho^2 - 2 * rho * cos(d))))
  expect_equal(stable, naive, tolerance = 1e-12)
  # and the stable form survives rho within a few ulp of 1
  expect_true(is.finite(log1p(-(1 - 1e-14)) + log1p(1 - 1e-14) -
                          log(2 * pi) -
                          log((1e-14)^2 + 4 * (1 - 1e-14) * sin(0.5)^2)))
})
