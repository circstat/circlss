# FD spot-checks of pnlss's assembled ll on a real gam design, plus the
# stable-branch identities of the projected-normal density pieces.

test_that("pnlss ll gradient and Hessian match finite differences", {
  set.seed(11)
  n <- 150
  x <- runif(n)
  m1 <- 1.8 * sin(2 * pi * x) + 0.5
  m2 <- 1.2 * cos(2 * pi * x) - 0.3
  y <- atan2(m2 + rnorm(n), m1 + rnorm(n))
  dat <- data.frame(y = y, x = x)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6)),
                 family = pnlss(), data = dat, fit = FALSE)
  fam <- G$family
  X <- G$X
  p <- ncol(X)
  set.seed(8)
  coef <- rnorm(p, sd = 0.3)

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

test_that("pn_mills_inv is the inverse Mills ratio on the stable branch", {
  t <- c(-30, -8, -1, 0, 1, 8, 30)
  expect_equal(circlss:::pn_mills_inv(t), dnorm(t) / pnorm(t),
               tolerance = 1e-12)
  # extreme negative t: R(t) ~ |t| (naive dnorm/pnorm is 0/0 there)
  expect_equal(circlss:::pn_mills_inv(-1e4), 1e4, tolerance = 1e-4)
})
