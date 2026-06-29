# Tests for the weight-aware linear-response forks (gausslss, gammalss) that
# carry the l~c leg, and the geometry detection that makes circ_gam's print/plot
# first-class across the circular-linear / circular-circular / linear-circular
# trio. Parity is against the mgcv originals at unit weights; weight-awareness is
# the integer-weight == row-replication equivalence (stock gaulss/gammals ignore
# weights, which is exactly what these forks fix).

test_that("gausslss matches mgcv::gaulss at unit weights", {
  set.seed(1); n <- 250; x <- runif(n)
  d <- data.frame(y = 2 + 1.5 * sin(2 * pi * x) + rnorm(n) * 0.4, x = x)
  f <- list(y ~ s(x, k = 8), ~ s(x, k = 8))
  b0 <- mgcv::gam(f, family = mgcv::gaulss(), data = d, method = "REML")
  bx <- circ_gam(f, data = d, family = gausslss())
  expect_equal(unname(coef(bx)), unname(coef(b0)), tolerance = 1e-4)
})

test_that("gausslss honours prior weights (integer weights == replication)", {
  set.seed(1); n <- 200; x <- runif(n)
  d <- data.frame(y = 1 + 2 * x + rnorm(n) * 0.5, x = x)
  set.seed(2); w <- sample(1:3, n, TRUE)
  bw <- circ_gam(list(y ~ x, ~ x), data = d, family = gausslss(), weights = w)
  br <- circ_gam(list(y ~ x, ~ x), data = d[rep(seq_len(n), w), ], family = gausslss())
  expect_equal(unname(coef(bw)), unname(coef(br)), tolerance = 1e-5)
  b0 <- circ_gam(list(y ~ x, ~ x), data = d, family = gausslss())
  expect_gt(max(abs(coef(bw) - coef(b0))), 1e-3)   # weights actually bite
})

test_that("gammalss matches mgcv::gammals and honours weights", {
  set.seed(1); n <- 250; x <- runif(n)
  d <- data.frame(y = rgamma(n, 4, rate = 4 / exp(0.4 + 0.8 * x)), x = x)
  f <- list(y ~ s(x, k = 8), ~ s(x, k = 8))
  b0 <- mgcv::gam(f, family = mgcv::gammals(), data = d, method = "REML")
  bx <- circ_gam(f, data = d, family = gammalss())
  expect_equal(unname(coef(bx)), unname(coef(b0)), tolerance = 1e-4)
  set.seed(2); w <- sample(1:3, n, TRUE)
  bw <- circ_gam(list(y ~ x, ~ x), data = d, family = gammalss(), weights = w)
  br <- circ_gam(list(y ~ x, ~ x), data = d[rep(seq_len(n), w), ], family = gammalss())
  expect_equal(unname(coef(bw)), unname(coef(br)), tolerance = 1e-5)
})

test_that("the linear forks carry linear-response metadata", {
  for (fam in list(gausslss(), gammalss())) {
    expect_false(is.null(fam$param_names))
    expect_length(fam$param_names, fam$nlp)
    expect_false(fam$response_circular)        # the flag that selects the can
    expect_true(all(!fam$param_circular))
  }
})

test_that(".circ_geometry classifies the response x covariate 2x2", {
  set.seed(1); n <- 120; x <- runif(n); phi <- runif(n, -pi, pi)
  dcl <- data.frame(th = 2 * atan(sin(2 * pi * x)) + rnorm(n) / 3, x = x)
  dcc <- data.frame(th = 2 * atan(sin(phi)) + rnorm(n) / 3, phi = phi)
  dlc <- data.frame(y = 2 + sin(phi) + rnorm(n) / 3, phi = phi)
  dll <- data.frame(y = 2 + x + rnorm(n) / 3, x = x)
  cl <- circ_gam(th ~ s(x, k = 6), data = dcl, family = vmlss())
  cc <- circ_gam(th ~ s(phi, bs = "cc", k = 6), data = dcc, family = vmlss())
  lc <- circ_gam(y ~ s(phi, bs = "cc", k = 6), data = dlc, family = gausslss())
  ll <- circ_gam(y ~ s(x, k = 6), data = dll, family = gausslss())
  expect_identical(circlss:::.circ_geometry(cl)$kind, "cl")
  expect_identical(circlss:::.circ_geometry(cc)$kind, "cc")
  expect_identical(circlss:::.circ_geometry(lc)$kind, "lc")
  expect_identical(circlss:::.circ_geometry(ll)$kind, "ll")
  expect_false(circlss:::.circ_geometry(lc)$resp_circular)
  expect_true(circlss:::.circ_geometry(lc)$cov_circular)
})

test_that("circ_gam gives a linear fork named columns and the ~ 1 fill", {
  set.seed(1); n <- 150; phi <- runif(n, -pi, pi)
  d <- data.frame(y = 2 + sin(phi) + rnorm(n) / 3, phi = phi)
  b <- circ_gam(y ~ s(phi, bs = "cc", k = 6), data = d, family = gausslss())  # 1 formula
  pr <- predict(b, newdata = data.frame(phi = c(0, 1)), type = "response")
  expect_identical(colnames(pr), c("mu", "tau"))
  expect_equal(b$family$nlp, 2)            # filled to 2 LPs
})

test_that("plot.circ_gam draws the l~c can and flat view without error", {
  set.seed(1); n <- 150; phi <- runif(n, -pi, pi)
  d <- data.frame(y = 2 + 1.5 * sin(phi) + rnorm(n) / 3, phi = phi)
  b <- circ_gam(list(y ~ s(phi, bs = "cc", k = 8), ~ s(phi, bs = "cc", k = 8)),
                data = d, family = gausslss())
  tf <- tempfile(fileext = ".pdf"); grDevices::pdf(tf)
  on.exit({ grDevices::dev.off(); unlink(tf) }, add = TRUE)
  expect_error(plot(b, view = "flat"), NA)        # the loc = NA regression
  expect_error(plot(b, view = "geometry"), NA)    # the can
  expect_error(plot(b, view = "both"), NA)
})

test_that("plot.circ_gam bands the c~c torus with the fitted law's circular SD", {
  set.seed(1); n <- 120; phi <- runif(n, -pi, pi)
  dir <- phi + 0.5 * sin(phi)
  d <- data.frame(theta = atan2(2 * sin(dir) + rnorm(n), 2 * cos(dir) + rnorm(n)),
                  phi = phi)
  b <- circ_gam(list(theta ~ s(phi, bs = "cc", k = 6), ~ s(phi, bs = "cc", k = 6)),
                data = d, family = pnlss())
  expect_identical(circlss:::.circ_geometry(b)$kind, "cc")

  ## the direction band is +/- the projected-normal law's circular SD
  ## sqrt(-2 log R) -- finite, in [0, pi], and not identically 0
  nd  <- data.frame(phi = seq(-pi, pi, length.out = 25))
  rp  <- predict(b, newdata = nd, type = "response")
  csd <- circlss:::.circ_sd_quad(b$family, rp)
  expect_length(csd, 25)
  expect_true(all(is.finite(csd) & csd >= 0 & csd <= pi))
  expect_true(any(csd > 0))

  ## suppressWarnings: the base pdf() device substitutes the en-dash in the
  ## default title font (mbcsToSbcs) -- cosmetic, unrelated to the band
  tf <- tempfile(fileext = ".pdf"); grDevices::pdf(tf)
  on.exit({ grDevices::dev.off(); unlink(tf) }, add = TRUE)
  expect_error(suppressWarnings(plot(b, view = "geometry")), NA)             # torus + ribbon
  expect_error(suppressWarnings(plot(b, view = "both")), NA)
  expect_error(suppressWarnings(plot(b, view = "geometry", se = FALSE)), NA) # ribbon suppressed
})
