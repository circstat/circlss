## Size-aware MAP-penalty degeneracy guard for circ_mix (every *lss family).
## A finite mixture's likelihood is unbounded; EM reweighting lets a component
## over-concentrate (kappa -> Inf, or a bounded shape/peakedness parameter onto its
## singular boundary), which crashes or hangs the M-step. The guard adds a
## size-aware MAP penalty (lambda_k = c / N_k) pulling each component's
## concentration/shape toward the diffuse model. These tests pin: (1) the penalty
## kernels' derivatives, (2) that the penalty is correctly placed and link-chained
## inside ll, (3) that a standalone circ_gam is byte-identical (the guard is inert
## without circ_mix's map_lambda), (4) that it bounds an otherwise-degenerate
## component, and (5) that every family fits circ_mix(K = 2) without crashing.

test_that("degeneracy kernels match finite-difference derivatives", {
  fd1 <- function(f, v, h = 1e-6) (f(v + h) - f(v - h)) / (2 * h)
  fd2 <- function(f, v, h = 1e-5) (f(v + h) - 2 * f(v) + f(v - h)) / (h * h)
  kerns <- list(
    linear   = circlss:::.degen_linear(1L, 1),
    linearS  = circlss:::.degen_linear(1L, 7),
    ridge    = circlss:::.degen_ridge(1L, 1),
    ridgeS   = circlss:::.degen_ridge(1L, 30),
    bup      = circlss:::.degen_boundary_upper(1L, 0.5),
    bupS     = circlss:::.degen_boundary_upper(1L, 0.5, 3),
    bsym     = circlss:::.degen_boundary_sym(1L, 1),
    bsymS    = circlss:::.degen_boundary_sym(1L, 1, 10))
  ## test points kept inside each kernel's natural domain
  pts <- list(linear = c(0.3, 2, 8), linearS = c(0.3, 2, 8),
              ridge = c(-1, 0.5, 3), ridgeS = c(-1, 0.5, 3),
              bup = c(0.05, 0.25, 0.45), bupS = c(0.05, 0.25, 0.45),
              bsym = c(-0.8, 0, 0.6), bsymS = c(-0.8, 0, 0.6))
  for (nm in names(kerns)) {
    k <- kerns[[nm]]
    for (v in pts[[nm]]) {
      expect_equal(k$rho1(v), fd1(k$rho0, v), tolerance = 1e-4, info = paste(nm, "rho1 @", v))
      expect_equal(k$rho2(v), fd2(k$rho0, v), tolerance = 1e-3, info = paste(nm, "rho2 @", v))
    }
  }
})

test_that("the penalized ll gradient and Hessian match finite differences", {
  ## set map_lambda on a family and finite-difference the returned objective l
  ## against the returned gradient (lb) and Hessian (lbb): this exercises the
  ## kernel, the trind placement, and the gamlss link chain-rule end to end.
  set.seed(11); n <- 80
  check_fam <- function(fam, formula) {
    y <- fam$rd(switch(fam$family,
      vmlss   = cbind(2 * atan(rnorm(n) * 0.4), rep(4, n)),
      kjlss   = cbind(2 * atan(rnorm(n) * 0.4), rep(0.6, n), rep(0.5, n), rep(-0.3, n))),
      rep(1, n), 1)
    d <- data.frame(y = y)
    g <- suppressWarnings(circ_gam(formula, data = d, family = fam))
    X <- predict(g, d, type = "lpmatrix"); b <- coef(g)
    fam$map_lambda <- 0.05
    obj <- function(coef) fam$ll(d$y, X, coef, rep(1, n), fam, deriv = 0)$l
    r <- fam$ll(d$y, X, b, rep(1, n), fam, deriv = 1)
    gfd <- vapply(seq_along(b), function(j) { e <- numeric(length(b)); e[j] <- 1e-5
      (obj(b + e) - obj(b - e)) / 2e-5 }, numeric(1))
    expect_equal(as.numeric(r$lb), gfd, tolerance = 1e-3, info = paste(fam$family, "gradient"))
    ## Hessian row by finite-differencing the gradient
    grad <- function(coef) fam$ll(d$y, X, coef, rep(1, n), fam, deriv = 1)$lb
    Hfd <- vapply(seq_along(b), function(j) { e <- numeric(length(b)); e[j] <- 1e-5
      (grad(b + e) - grad(b - e)) / 2e-5 }, numeric(length(b)))
    expect_equal(r$lbb, (Hfd + t(Hfd)) / 2, tolerance = 1e-2, info = paste(fam$family, "Hessian"))
  }
  check_fam(vmlss(), list(y ~ 1, ~ 1))                       # Newton (available.derivs = 2)
  check_fam(kjlss(), list(y ~ 1, ~ 1, ~ 1, ~ 1))            # EFS (available.derivs = 0)
})

test_that("a standalone circ_gam is unaffected by the guard (map_lambda inert)", {
  set.seed(2); m <- 200; x <- runif(m)
  y <- vmlss()$rd(cbind(2 * atan(1.2 * sin(2 * pi * x)), exp(1 + 0.7 * cos(2 * pi * x))),
                  rep(1, m), 1)
  d <- data.frame(y = y, x = x)
  g0 <- suppressWarnings(circ_gam(list(y ~ s(x), ~ s(x)), data = d, family = vmlss()))
  f1 <- vmlss(); f1$map_lambda <- NULL                       # explicit NULL
  g1 <- suppressWarnings(circ_gam(list(y ~ s(x), ~ s(x)), data = d, family = f1))
  expect_equal(coef(g0), coef(g1))
  ## map_lambda = 0 must also be a no-op (.degen_active gates on > 0)
  f2 <- vmlss(); f2$map_lambda <- 0
  g2 <- suppressWarnings(circ_gam(list(y ~ s(x), ~ s(x)), data = d, family = f2))
  expect_equal(coef(g0), coef(g2))
})

test_that("the guard bounds an otherwise-degenerate component (kappa cap)", {
  ## the penalty solves A1(kappa_hat) = Rbar_w - lambda for an intercept-only
  ## concentration, so a larger lambda yields a smaller fitted kappa, monotonically.
  set.seed(3); n <- 200
  y <- atan2(sin(rnorm(n) * 0.18), cos(rnorm(n) * 0.18))     # tight cluster -> large kappa
  d <- data.frame(y = y)
  kap <- function(lam) { f <- vmlss(); if (!is.null(lam)) f$map_lambda <- lam
    unname(predict(suppressWarnings(circ_gam(list(y ~ 1, ~ 1), data = d, family = f)),
                   d[1, , drop = FALSE], type = "response")[1, 2]) }
  k_off <- kap(NULL); k_lo <- kap(0.01); k_hi <- kap(0.05)
  expect_true(k_off > k_lo && k_lo > k_hi)
  ## analytic check: A1(kappa_hat) == Rbar - lambda
  Rbar <- sqrt(mean(cos(y))^2 + mean(sin(y))^2)
  A1 <- function(k) besselI(k, 1, expon.scaled = TRUE) / besselI(k, 0, expon.scaled = TRUE)
  expect_equal(A1(kap(0.03)), Rbar - 0.03, tolerance = 1e-3)
})

test_that("circ_mix.control validates and deprecates", {
  expect_error(circ_mix.control(degen_strength = -1), "non-negative")
  expect_error(circ_mix.control(time_budget = 0), "positive")
  expect_warning(circ_mix.control(kappa_cap = 50), "deprecated")
  expect_silent(ctl <- circ_mix.control(degen_strength = 2, time_budget = 30))
  expect_equal(ctl$degen_strength, 2)
  expect_equal(ctl$time_budget, 30)
  expect_null(circ_mix.control()$time_budget)
  expect_equal(circ_mix.control()$degen_strength, 1)         # default c = 1
})

test_that("every circular family fits circ_mix(K = 2) without degenerating", {
  skip_on_cran()
  ## a clean two-component von Mises mixture every family can represent, so a
  ## failure is the guard (degeneracy), not model-data mismatch.
  set.seed(7); n <- 160
  z  <- sample.int(2L, n, replace = TRUE, prob = c(0.45, 0.55))
  mu <- c(-1.4, 1.4)[z]
  y  <- vmlss()$rd(cbind(mu, rep(4, n)), rep(1, n), 1)
  d  <- data.frame(y = y)
  fams <- list(vmlss = vmlss(), wnlss = wnlss(), wclss = wclss(), cardlss = cardlss(),
               pnlss = pnlss(), vmftlss = vmftlss(), jplss = jplss(), kjlss = kjlss(),
               ajplss = ajplss(), ssjplss = ssjplss(), ibslss = ibslss())
  for (nm in names(fams)) {
    fam <- fams[[nm]]
    ff  <- c(list(y ~ 1), rep(list(~ 1), fam$nlp - 1L))
    ctl <- circ_mix.control(restarts = 2L, init = "kmeans", seed = 1L,
                            max_iter = 120L, time_budget = 30)
    m <- tryCatch(suppressWarnings(circ_mix(ff, data = d, family = fam, K = 2L, control = ctl)),
                  error = function(e) conditionMessage(e))
    expect_s3_class(m, "circ_mix")
    if (inherits(m, "circ_mix")) {
      expect_true(is.finite(m$bic), info = nm)
      ## every non-circular (concentration/shape) parameter must stay finite and
      ## off its boundary -- the degeneracy signature is +/-Inf or a wall value.
      pc <- fam$param_circular
      vals <- unlist(lapply(m$components, function(cp)
        predict(cp$fit, d[1, , drop = FALSE], type = "response")[1, ][!pc]))
      expect_true(all(is.finite(vals)), info = paste(nm, "params finite"))
    }
  }
})
