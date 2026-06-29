## circ_kmeans: circular k-means on the circle / torus. Deterministic (seeded,
## base-R arithmetic only -- no model fitting), so the partitions are stable.

test_that("circ_kmeans recovers two separated mean directions", {
  set.seed(42)
  theta <- c(stats::rnorm(60, 0, 0.25), stats::rnorm(60, pi, 0.25)) %% (2 * pi)
  km <- circ_kmeans(theta, 2)
  expect_length(km$cluster, 120L)
  expect_identical(sort(unique(km$cluster)), 1:2)
  expect_equal(sort(km$size), c(60L, 60L))
  ## one centre near 0, one near pi (robust to label swap and the +/-pi seam)
  c1     <- km$centers[, 1L]
  near0  <- abs(atan2(sin(c1),      cos(c1)))
  nearpi <- abs(atan2(sin(c1 - pi), cos(c1 - pi)))
  expect_true(any(near0 < 0.15))
  expect_true(any(nearpi < 0.15))
})

test_that("circ_kmeans clusters across the wrap-around seam", {
  set.seed(1)
  ## cluster A straddles 0 (points just above 0 and just below 2pi); B sits at pi.
  ## Euclidean k-means on the raw radians would tear A apart at the 0 / 2pi seam.
  A  <- stats::runif(40, -0.3, 0.3) %% (2 * pi)
  B  <- stats::runif(40, pi - 0.3, pi + 0.3)
  km <- circ_kmeans(c(A, B), 2)
  la <- km$cluster[seq_len(40)]
  lb <- km$cluster[41:80]
  expect_length(unique(la), 1L)            # A stays whole across the seam
  expect_length(unique(lb), 1L)
  expect_false(unique(la) == unique(lb))   # A and B are different clusters
})

test_that("circ_kmeans centre is the circular mean", {
  ## circular mean of {-0.2, -0.1, 0.1, 0.2} is exactly 0 (sin cancels)
  km <- circ_kmeans(c(-0.2, -0.1, 0.1, 0.2), 1)
  expect_equal(dim(km$centers), c(1L, 1L))
  expect_equal(as.numeric(km$centers), 0, tolerance = 1e-8)
  expect_identical(km$size, 4L)
})

test_that("circ_kmeans clusters on the torus (two angular coordinates)", {
  set.seed(7)
  n   <- 50L
  th1 <- c(stats::rnorm(n, 0, 0.2),    stats::rnorm(n, pi, 0.2))    %% (2 * pi)
  th2 <- c(stats::rnorm(n, pi / 2, 0.2), stats::rnorm(n, -pi / 2, 0.2)) %% (2 * pi)
  km  <- circ_kmeans(cbind(th1, th2), 2)
  expect_equal(dim(km$centers), c(2L, 2L))
  expect_equal(sort(km$size), c(n, n))
})

test_that("circ_kmeans guards its arguments", {
  expect_error(circ_kmeans(c(0, 1, 2), 5), "more cluster centres")
  expect_error(circ_kmeans(c(0, 1, 2), 0), "centers")
})

## ---- geometry dispatch used by circ_mix --------------------------------------
test_that(".circ_mix_kmeans routes on geometry and returns kmeans-shaped output", {
  set.seed(3)
  ang <- matrix(c(stats::rnorm(20, 0, 0.2), stats::rnorm(20, pi, 0.2)) %% (2 * pi),
                ncol = 1L)
  ck  <- circlss:::.circ_mix_kmeans(ang, 2L, circular = TRUE)
  expect_identical(names(ck),
                   c("cluster", "centers", "withinss", "tot.withinss", "size"))
  expect_length(ck$cluster, 40L)

  lin <- matrix(c(stats::rnorm(20, -4), stats::rnorm(20, 4)), ncol = 1L)
  ek  <- circlss:::.circ_mix_kmeans(lin, 2L, circular = FALSE)
  expect_true(all(c("cluster", "centers") %in% names(ek)))   # stats::kmeans shape
})

test_that(".circ_mix_assign labels by the geometry's nearest centre", {
  ## circular: points near 0 and near pi, centres at 0 and pi
  x  <- matrix(c(0.05, -0.05, pi + 0.05, pi - 0.05), ncol = 1L)
  mu <- matrix(c(0, pi), ncol = 1L)
  expect_identical(circlss:::.circ_mix_assign(x, mu, circular = TRUE),
                   c(1L, 1L, 2L, 2L))
  ## linear: squared-Euclidean nearest centre
  xl  <- matrix(c(-1, -0.9, 1, 1.1), ncol = 1L)
  mul <- matrix(c(-1, 1), ncol = 1L)
  expect_identical(circlss:::.circ_mix_assign(xl, mul, circular = FALSE),
                   c(1L, 1L, 2L, 2L))
})
