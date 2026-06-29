## Kato-Jones (2015) numerics for the kjlss general family: the Cartesian
## scores, the disc-chart pieces, the log-density and its first two derivatives,
## the moment fit, and the chart inverse.
##
## Kato-Jones is the only family whose normalizer is EXACTLY 2*pi: the density
##   g(theta) = (1/2pi)[1 + 2 gamma (cos(theta-mu) - rho cos lam)
##                          / (1 + rho^2 - 2 rho cos(theta-mu-lam))]
## is a first-/second-trigonometric-moment perturbation of the uniform that
## integrates to 1 with no special function. So there is no quadrature and no
## Bessel/Gamma term -- the entire family is elementary rational functions of
## the Cartesian shape pair (a, b) = (rho cos lam, rho sin lam) and c = cos d,
## s = sin d (d = theta - mu). In those coordinates
##   l = log N - log 2pi,   N = 1 + 2 gamma (c - a)/D,
##   D = 1 + a^2 + b^2 - 2 a c - 2 b s.
##
## The one genuinely novel trick is the DISC CHART: the
## Theorem-1 feasible region for (a, b) is the closed disc of centre (gamma, 0)
## and radius 1 - gamma, and the chart
##   (a, b) = (gamma, 0) + (1 - gamma) u / sqrt(1 + ||u||^2),   u in R^2,
## maps an unconstrained u = (u1, u2) onto its open interior. So kjlss smooths
## the unconstrained chart coordinates (u1, u2) -- identity links -- and the
## coupled feasibility constraint can never be violated by any coefficient
## vector. u = 0 is exactly the wrapped Cauchy WC(mu, gamma) (the nesting test
## point). All derivatives below are w.r.t. the CHART parameters (mu, gamma,
## u1, u2): the Cartesian scores pushed through the chart Jacobian/Hessian.

## Chart value plus the Jacobian/Hessian pieces the chain rule needs. The chart
## is linear in gamma, so only the mixed gamma-u and u-u second derivatives of
## (a, b) survive (carried as J = d v / d u and H1/H2 = d^2 v / d u^2).
.kj_chart <- function(gamma, u1, u2, second = FALSE) {
  r2 <- 1 + u1 * u1 + u2 * u2
  r <- sqrt(r2)
  r3 <- r2 * r
  v1 <- u1 / r
  v2 <- u2 / r
  om <- 1 - gamma
  J11 <- 1 / r - u1 * u1 / r3
  J12 <- -u1 * u2 / r3
  J22 <- 1 / r - u2 * u2 / r3
  out <- list(
    a = gamma + om * v1, b = om * v2,
    a_g = 1 - v1, b_g = -v2,
    au1 = om * J11, au2 = om * J12,
    bu1 = om * J12, bu2 = om * J22,
    J11 = J11, J12 = J12, J22 = J22, om = om
  )
  if (second) {
    r5 <- r3 * r2
    ## columns: (u1,u1), (u1,u2), (u2,u2)
    out$H1 <- cbind(-3 * u1 / r3 + 3 * u1^3 / r5,
                    -u2 / r3 + 3 * u1 * u1 * u2 / r5,
                    -u1 / r3 + 3 * u1 * u2 * u2 / r5)
    out$H2 <- cbind(-u2 / r3 + 3 * u1 * u1 * u2 / r5,
                    -u1 / r3 + 3 * u1 * u2 * u2 / r5,
                    -3 * u2 / r3 + 3 * u2^3 / r5)
  }
  out
}

## Per-observation log-density value AND (when deriv > 0) the chart-coordinate
## gradient l1 (n x 4: mu, gamma, u1, u2) and Hessian l2 (n x 10, the
## combinations_with_replacement order mm, mg, m1, m2, gg, g1, g2, 11, 12, 22).
## All inputs are length-n vectors (one parameter value per observation).
.kj_terms <- function(y, mu, gamma, u1, u2, deriv = TRUE) {
  ch <- .kj_chart(gamma, u1, u2, second = deriv)
  a <- ch$a
  b <- ch$b
  d <- y - mu
  cc <- cos(d)
  s <- sin(d)
  D <- pmax(1 + a * a + b * b - 2 * a * cc - 2 * b * s, 1e-300)
  E <- cc - a
  Fv <- E / D
  N <- pmax(1 + 2 * gamma * Fv, 1e-300)
  l0 <- log(N) - log(2 * pi)
  if (!deriv) return(list(l0 = l0, l1 = NULL, l2 = NULL))

  D2 <- D * D
  D3 <- D2 * D
  N2 <- N * N
  ## numerator E1 = d(c - a)/d.  : mu -> s, a -> -1, b -> 0
  ## denominator D1 : mu -> 2(b c - a s), a -> 2(a - c), b -> 2(b - s)
  D1_mu <- 2 * (b * cc - a * s)
  D1_a <- 2 * (a - cc)
  D1_b <- 2 * (b - s)
  F1_mu <- (s * D - E * D1_mu) / D2
  F1_a <- (-D - E * D1_a) / D2          # E1_a = -1
  F1_b <- (-E * D1_b) / D2              # E1_b = 0

  ## Cartesian gradient (mu, gamma, a, b)
  l1_mu <- 2 * gamma * F1_mu / N
  l1_g <- 2 * Fv / N
  la <- 2 * gamma * F1_a / N
  lb <- 2 * gamma * F1_b / N

  ## chart gradient
  l1 <- cbind(l1_mu,
              l1_g + la * ch$a_g + lb * ch$b_g,
              la * ch$au1 + lb * ch$bu1,
              la * ch$au2 + lb * ch$bu2)

  ## Cartesian second derivatives of F = E/D, then of l = log N
  ## Fxy = Exy/D - (E1x D1y + E1y D1x + E Dxy)/D2 + 2 E D1x D1y/D3
  Fmm <- (-cc) / D - (2 * s * D1_mu + E * (2 * (b * s + a * cc))) / D2 +
    2 * E * D1_mu * D1_mu / D3
  Fma <- -(s * D1_a - D1_mu + E * (-2 * s)) / D2 + 2 * E * D1_mu * D1_a / D3
  Fmb <- -(s * D1_b + E * (2 * cc)) / D2 + 2 * E * D1_mu * D1_b / D3
  Faa <- -(-2 * D1_a + E * 2) / D2 + 2 * E * D1_a * D1_a / D3
  Fab <- -(-D1_b) / D2 + 2 * E * D1_a * D1_b / D3
  Fbb <- -(E * 2) / D2 + 2 * E * D1_b * D1_b / D3

  lmm <- 2 * gamma * Fmm / N - 4 * gamma * gamma * F1_mu * F1_mu / N2
  lma <- 2 * gamma * Fma / N - 4 * gamma * gamma * F1_mu * F1_a / N2
  lmb <- 2 * gamma * Fmb / N - 4 * gamma * gamma * F1_mu * F1_b / N2
  laa <- 2 * gamma * Faa / N - 4 * gamma * gamma * F1_a * F1_a / N2
  lab <- 2 * gamma * Fab / N - 4 * gamma * gamma * F1_a * F1_b / N2
  lbb <- 2 * gamma * Fbb / N - 4 * gamma * gamma * F1_b * F1_b / N2
  lgg <- -4 * Fv * Fv / N2
  lgm <- 2 * F1_mu / N - 4 * gamma * Fv * F1_mu / N2
  lga <- 2 * F1_a / N - 4 * gamma * Fv * F1_a / N2
  lgb <- 2 * F1_b / N - 4 * gamma * Fv * F1_b / N2

  ag <- ch$a_g
  bg <- ch$b_g
  au1 <- ch$au1; au2 <- ch$au2; bu1 <- ch$bu1; bu2 <- ch$bu2
  om <- ch$om
  H1 <- ch$H1; H2 <- ch$H2

  ## chart Hessian, in combinations_with_replacement(mu, gamma, u1, u2) order
  l2 <- cbind(
    lmm,                                                  # (mu, mu)
    lgm + lma * ag + lmb * bg,                            # (mu, gamma)
    lma * au1 + lmb * bu1,                                # (mu, u1)
    lma * au2 + lmb * bu2,                                # (mu, u2)
    lgg + 2 * (lga * ag + lgb * bg) +
      laa * ag * ag + 2 * lab * ag * bg + lbb * bg * bg,  # (gamma, gamma)
    lga * au1 + lgb * bu1 + (laa * ag + lab * bg) * au1 +
      (lab * ag + lbb * bg) * bu1 - la * ch$J11 - lb * ch$J12,  # (gamma, u1)
    lga * au2 + lgb * bu2 + (laa * ag + lab * bg) * au2 +
      (lab * ag + lbb * bg) * bu2 - la * ch$J12 - lb * ch$J22,  # (gamma, u2)
    laa * au1 * au1 + 2 * lab * au1 * bu1 + lbb * bu1 * bu1 +
      la * om * H1[, 1] + lb * om * H2[, 1],             # (u1, u1)
    laa * au1 * au2 + lab * (au1 * bu2 + bu1 * au2) + lbb * bu1 * bu2 +
      la * om * H1[, 2] + lb * om * H2[, 2],             # (u1, u2)
    laa * au2 * au2 + 2 * lab * au2 * bu2 + lbb * bu2 * bu2 +
      la * om * H1[, 3] + lb * om * H2[, 3]              # (u2, u2)
  )
  list(l0 = l0, l1 = l1, l2 = l2)
}

## Closed-form chart inverse: book shape (rho, lam) -> chart (u1, u2). Boundary
## or infeasible (rho, lam) are pulled radially to vmax * disc radius so the
## inverse stays finite (a moments fit can land on the Theorem-1 circle).
.kj_disc_chart_inverse <- function(gamma, rho, lam, vmax = 1 - 1e-9) {
  om <- max(1 - gamma, 1e-12)
  v1 <- (rho * cos(lam) - gamma) / om
  v2 <- rho * sin(lam) / om
  nrm <- sqrt(v1 * v1 + v2 * v2)
  scale <- if (nrm > vmax) vmax / max(nrm, 1e-300) else 1
  v1 <- v1 * scale
  v2 <- v2 * scale
  den <- sqrt(max(1 - (v1 * v1 + v2 * v2), 1e-18))
  c(v1 / den, v2 / den)
}

## Method-of-moments fit (mu, gamma, rho, lam) used for the constant-target
## starts in initialize. gamma = mean resultant length R-bar; the second-order
## moment (alpha2, beta2) is projected into the feasible disc before rho/lam are
## read off. Transcribed from katojones_gen._fit_moments / _project_second_order.
.kj_moment_fit <- function(y) {
  sy <- mean(sin(y))
  cy <- mean(cos(y))
  mu0 <- atan2(sy, cy)
  Rbar <- sqrt(sy * sy + cy * cy)
  g0 <- min(max(Rbar, 0), 1 - 1e-9)
  centered <- y - mu0
  a2 <- mean(cos(2 * centered))
  b2 <- mean(sin(2 * centered))
  ## project (a2, b2) into the disc of centre (g0^2, 0), radius g0 (1 - g0)
  radius <- g0 * (1 - g0)
  ca <- g0 * g0
  va <- a2 - ca
  vb <- b2
  dist <- sqrt(va * va + vb * vb)
  if (radius <= 0) {
    ap <- ca; bp <- 0
  } else if (dist <= radius) {
    ap <- a2; bp <- b2
  } else if (dist == 0) {
    ap <- ca + radius; bp <- 0
  } else {
    sc <- radius / dist
    ap <- ca + va * sc
    bp <- vb * sc
  }
  if (g0 < 1e-12) {
    rho0 <- 0
    lam0 <- 0
  } else {
    r2 <- sqrt(ap * ap + bp * bp)
    rho0 <- min(max(r2 / max(g0, 1e-12), 0), 1 - 1e-9)
    lam0 <- if (rho0 < 1e-12) 0 else atan2(bp, ap) %% (2 * pi)
  }
  c(mu0, g0, rho0, lam0)
}

## Per-observation saturated log-density (max over theta) for the deviance.
## The density profile in d = theta - mu depends only on (gamma, u1, u2) -- mu
## is a pure shift -- so the max is grouped by the shape triple.
.kj_lsat <- function(mu, gamma, u1, u2) {
  grid <- seq(-pi, pi, length.out = 721L)
  cc <- cos(grid)
  s <- sin(grid)
  out <- numeric(length(mu))
  key <- paste(gamma, u1, u2, sep = "\r")
  for (k in unique(key)) {
    idx <- which(key == k)
    ch <- .kj_chart(gamma[idx[1]], u1[idx[1]], u2[idx[1]])
    a <- ch$a
    b <- ch$b
    g <- gamma[idx[1]]
    D <- pmax(1 + a * a + b * b - 2 * a * cc - 2 * b * s, 1e-300)
    N <- pmax(1 + 2 * g * (cc - a) / D, 1e-300)
    out[idx] <- max(log(N) - log(2 * pi))
  }
  out
}
