# Prior-weights test specific to vmlss, the reference family: weighted-gradient
# finite differences, an FD cross-check of the derivative scaling that is
# independent of the exact duplication identity in test-weights.R. See
# helper-weights.R for the expect_* helper.

sim_vm_w <- function(n, seed = 42) {
  set.seed(seed)
  x <- runif(n)
  mu <- 2 * atan(1.2 * sin(2 * pi * x))
  kappa <- exp(1 + 0.6 * cos(2 * pi * x))
  y <- atan2(sin(mu + rnorm(n) / sqrt(kappa)), cos(mu + rnorm(n) / sqrt(kappa)))
  data.frame(y = y, x = x)
}

test_that("weighted ll gradient and Hessian match finite differences", {
  dat <- sim_vm_w(120)
  G <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(x, k = 6)),
                 family = vmlss(), data = dat, fit = FALSE)
  set.seed(7)
  coef <- rnorm(ncol(G$X), sd = 0.1)
  wt <- runif(length(G$y), 0.2, 3)     # random positive prior weights
  expect_weighted_gradient_fd(G$family, G$X, G$y, coef, wt)
})
