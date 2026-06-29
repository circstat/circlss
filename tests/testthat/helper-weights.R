## Shared helpers for the prior-weights ("weighted likelihood") tests: a
## weight-aware *lss family must satisfy "weighting a row by w == duplicating
## that row w times", which is what lets circ_gam(..., weights = w) fit the
## weighted MLE. Both checks below work at the ll level (no model fitting).

## (1) EXACT, the core gate. At the ll level, for positive *integer* weights w,
## ll(wt = w) on n rows must equal ll(wt = 1) on the design with row i repeated
## w_i times -- bit-for-bit (to floating point), for the objective l, gradient
## lb and Hessian lbb. This is an identity, not an approximation: a weighted
## log-likelihood scales each per-observation derivative row by w_i, which is the
## same sum as stacking w_i identical rows. It isolates the family's derivative-
## block scaling from any optimiser/FD conditioning, so it is robust even for the
## families whose log-density has near-singular regions (cartlss antipode,
## ibslss Newton-inverted warps) where finite differences are stiff.
expect_weight_equals_duplication <- function(fam, X, y, coef, w, tol = 1e-8) {
  jj <- attr(X, "lpi")
  idx <- rep(seq_along(y), w)
  Xd <- X[idx, , drop = FALSE]
  attr(Xd, "lpi") <- jj                 # gamlss.gH needs the lpi attribute
  yd <- y[idx]
  rw <- fam$ll(y,  X,  coef, wt = w,                  family = fam, deriv = 1)
  rd <- fam$ll(yd, Xd, coef, wt = rep(1, length(yd)), family = fam, deriv = 1)
  expect_equal(rw$l, rd$l, tolerance = tol)
  expect_equal(drop(rw$lb), drop(rd$lb), tolerance = tol)
  expect_equal(unname(as.matrix(rw$lbb)), unname(as.matrix(rd$lbb)), tolerance = tol)
}

## (2) Weighted-gradient finite differences -- ll(deriv=1)$lb/$lbb match the FD of
## the weighted objective ll(deriv=0)$l, an independent cross-check (on a well-
## conditioned design) that the per-observation derivatives themselves are right
## under non-unit weights.
expect_weighted_gradient_fd <- function(fam, X, y, coef, wt, h = 1e-5,
                                        tol = 1e-5) {
  llf <- function(b) fam$ll(y, X, b, wt = wt, family = fam, deriv = 0)$l
  ret <- fam$ll(y, X, coef, wt = wt, family = fam, deriv = 1)
  p <- length(coef)
  g_fd <- vapply(seq_len(p), function(j) {
    bp <- bm <- coef; bp[j] <- bp[j] + h; bm[j] <- bm[j] - h
    (llf(bp) - llf(bm)) / (2 * h)
  }, numeric(1))
  expect_equal(drop(ret$lb), g_fd, tolerance = tol)

  H_fd <- matrix(0, p, p)
  for (j in seq_len(p)) {
    bp <- bm <- coef; bp[j] <- bp[j] + h; bm[j] <- bm[j] - h
    gp <- fam$ll(y, X, bp, wt = wt, family = fam, deriv = 1)$lb
    gm <- fam$ll(y, X, bm, wt = wt, family = fam, deriv = 1)$lb
    H_fd[, j] <- (gp - gm) / (2 * h)
  }
  expect_equal(unname(as.matrix(ret$lbb)), (H_fd + t(H_fd)) / 2, tolerance = tol)
}

## Family registry for the systematic across-family tests: constructor name and
## number of linear predictors. Single source of truth so a new family is one
## row. (vmlss also gets its own richer test-vmlss-weights.R.)
.WEIGHT_TEST_FAMILIES <- list(
  vmlss = 2, wclss = 2, pnlss = 2, cardlss = 2, cartlss = 2, wnlss = 2,
  jplss = 3, vmftlss = 3, ssjplss = 4, kjlss = 4, ibslss = 4, ajplss = 4
)

## Formula list for an nlp-parameter family: a smooth in the mean and the second
## parameter, the rest intercept-only (keeps the systematic fits well-posed).
.weight_test_formula <- function(nlp, smooth2 = TRUE) {
  f <- list(y ~ s(x, k = 6), if (smooth2) ~ s(x, k = 6) else ~ 1)
  if (nlp >= 3) f <- c(f, ~ 1)
  if (nlp >= 4) f <- c(f, ~ 1)
  f
}

## Generic circular data valid for every family: a smooth mean direction with
## moderate dispersion (kept tight enough to stay clear of antipodal cusps).
.weight_test_data <- function(n, seed) {
  set.seed(seed)
  x <- stats::runif(n)
  mu <- 2 * atan(1.2 * sin(2 * pi * x))
  y <- atan2(sin(mu + stats::rnorm(n, sd = 0.5)),
             cos(mu + stats::rnorm(n, sd = 0.5)))
  data.frame(y = y, x = x)
}
