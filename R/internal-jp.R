## Jones-Pewsey internal numerics: the log-kernel, its parameter derivatives,
## and the quadrature normalizer/log-Z moments.
##
## The Jones-Pewsey density is
##   f(y) = c(kappa, psi) * (cosh(kappa psi) + sinh(kappa psi) cos(y - mu))^(1/psi),
## written through the overflow-safe decomposition g = e^A cos^2(d/2) +
## e^{-A} sin^2(d/2) with A = kappa psi and d = y - mu, so the log-kernel is
## h = (1/psi) log g. The normalizer c has no elementary form (it is an
## off-cut Legendre function); both it and the kappa/psi derivatives of
## log Z = -log c are kernel-weighted moments evaluated by composite 24-point
## Gauss-Legendre on a feature-scale break-point ladder -- the first circlss
## family to need numerical integration inside ll(). These functions are
## package-internal; the family closures (ll/residuals/rd) reach them by
## lexical scope, postproc by getFromNamespace (it is eval'd in mgcv's frame).

## Tolerances and caps.
.JP_KAPPA_TOL <- 1e-3   # below this kappa the density value is the uniform 1/2pi
.JP_PSI_TOL <- 1e-6     # below this |psi| the value reduces to von Mises
.JP_A_SMALL <- 1e-4     # |kappa psi| below this the psi-direction uses a series
.JP_LOG_CAP <- 250.0    # cap on log(Vtilde); tames the 0*inf antipode corner

## 24-point Gauss-Legendre nodes/weights on [-1, 1] by Golub-Welsch (eigen of
## the symmetric tridiagonal Jacobi matrix) -- the dependency-free twin of
## scipy.special.roots_legendre(24); agreement to ~1e-15 is far inside the
## quadrature error and the 2e-6 parity tolerance. Nodes sorted ascending so
## the composite-panel node order matches the Python side.
.gauss_legendre <- function(n) {
  k <- seq_len(n - 1)
  beta <- k / sqrt(4 * k * k - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- beta
  J[cbind(k + 1L, k)] <- beta
  ev <- eigen(J, symmetric = TRUE)
  nodes <- ev$values
  w <- 2 * (ev$vectors[1, ])^2
  ord <- order(nodes)
  list(nodes = nodes[ord], weights = w[ord])
}

.jp_gl24 <- .gauss_legendre(24L)

## Numerically stable log(exp(a) + exp(b)) = max + log1p(exp(-|a-b|)); both
## -Inf returns -Inf (the antipode corner, where one of a, b is finite anyway).
.jp_logaddexp <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  bad <- !is.finite(m)
  if (any(bad)) out[bad] <- m[bad]
  out
}

## The two length scales of the JP kernel in its own angle: the peak-curvature
## width at d = 0 and the antipodal near-kink scale 2 e^{-|A|} at d = pi, both
## clipped to [1e-320, 1] (a positivity guard at the denormal floor, NOT a
## resolution floor -- a floor above the true feature scale would truncate the
## ladder and lose the psi < 0 spike). Returns c(w_peak, w_anti).
.jp_feature_scales <- function(kappa, psi) {
  A <- kappa * psi
  if (abs(A) < 1e-12) {
    keff <- max(kappa, 1e-12)
    w_peak <- min(max(1 / sqrt(keff), 1e-320), 1.0)
  } else {
    log_keff <- if (A >= 0) {
      log(-expm1(-2 * A)) - log(2 * abs(psi))
    } else {
      log(expm1(-2 * A)) - log(2 * abs(psi))
    }
    w_peak <- min(max(exp(-0.5 * log_keff), 1e-320), 1.0)
  }
  w_anti <- min(max(2 * exp(-abs(A)), 1e-320), 1.0)
  c(w_peak, w_anti)
}

## Break-point ladder over the half circle [0, pi]: rungs grow geometrically
## (decade steps) from both ends at the peak (d = 0) and antipodal (d = pi)
## feature scales, so every panel spans at most one decade of its scale.
.jp_ladder_edges <- function(kappa, psi) {
  ws <- .jp_feature_scales(kappa, psi)
  edges <- c(0.0, pi)
  r <- ws[1]
  while (r < pi) { edges <- c(edges, r); r <- r * 10 }
  r <- ws[2]
  while (r < pi) { edges <- c(edges, pi - r); r <- r * 10 }
  sort.int(unique(edges))
}

## Composite 24-point Gauss-Legendre nodes/weights over arbitrary sorted panel
## edges. Built by one broadcast over the panels rather than a panel loop --
## bitwise the same as the per-panel fill it replaces (addition is commutative):
## panels in edge order, each in ascending node order, so the floating summation
## is order-deterministic. Shared by the symmetric half-circle sweep
## (.jp_gl_panels) and the asymmetric full-window sweep (.jp_asym_size_groups).
.jp_panels_from_edges <- function(edges) {
  xi <- .jp_gl24$nodes
  wgl <- .jp_gl24$weights
  L <- length(xi)
  lo <- edges[-length(edges)]
  hi <- edges[-1L]
  mid <- 0.5 * (lo + hi)
  hw <- 0.5 * (hi - lo)
  np <- length(lo)
  ## panel-major flatten (panel i's L entries contiguous):
  ##   node[(i-1)L + j] = mid[i] + hw[i] * xi[j],  weight = hw[i] * wgl[j].
  ## Built by rep() instead of as.vector(t(outer(hw, .) + mid)) -- no outer()/t()
  ## matrix allocations; bit-for-bit (the multiply is identical and add commutes).
  hw_rep <- rep(hw, each = L)
  list(nodes = rep(mid, each = L) + hw_rep * rep(xi, times = np),
       weights = hw_rep * rep(wgl, times = np))
}

## Composite GL panels over the ladder of [0, pi] (the symmetric kernel is even
## in d, so a half-circle sweep suffices).
.jp_gl_panels <- function(kappa, psi) {
  .jp_panels_from_edges(.jp_ladder_edges(kappa, psi))
}

## Unique (kappa, psi) GL grids grouped by node count for batched evaluation:
## each group stacks into one rectangular (g x L) node/weight matrix with NO
## padding -- every pair keeps its own adaptive .jp_ladder_edges grid, so a
## batched .jp_score_terms sweep is bit-for-bit the per-pair scalar, and a lone
## deep-psi spike (a longer ladder) forms its own one-row group instead of
## inflating every pair's node count. This is what lets
## .jp_logc_vec/.jp_quad_vec replace a per-pair loop with one sweep
## per distinct grid size. Returns a list of list(idx, nodes, wts), idx indexing
## the input (kappa, psi).
.jp_size_groups <- function(kappa, psi) {
  grids <- lapply(seq_along(kappa), function(i) .jp_gl_panels(kappa[i], psi[i]))
  sizes <- vapply(grids, function(g) length(g$nodes), integer(1))
  lapply(sort.int(unique(sizes)), function(L) {
    idx <- which(sizes == L)
    nodes <- matrix(unlist(lapply(grids[idx], `[[`, "nodes"), use.names = FALSE),
                    nrow = length(idx), ncol = L, byrow = TRUE)
    wts <- matrix(unlist(lapply(grids[idx], `[[`, "weights"), use.names = FALSE),
                  nrow = length(idx), ncol = L, byrow = TRUE)
    list(idx = idx, nodes = nodes, wts = wts)
  })
}

## The JP log-kernel h(d; kappa, psi) alone, for scalar kappa/psi and a vector
## of angles d (the quadrature integrand). psi = 0 is the von Mises member;
## |A| < A_SMALL the cumulant series; else the overflow-safe logaddexp form.
.jp_log_kernel <- function(d, kappa, psi) {
  if (psi == 0) return(kappa * cos(d))
  A <- kappa * psi
  if (abs(A) < .JP_A_SMALL) {
    s <- sin(d); cc <- cos(d); s2 <- s * s
    k3 <- -2 * cc * s2
    k4 <- s2 * (4 - 6 * s2)
    return(kappa * (cc + s2 * A / 2 + k3 * A * A / 6 + k4 * A^3 / 24))
  }
  half <- 0.5 * d
  lp <- 2 * log(abs(cos(half)))
  lq <- 2 * log(abs(sin(half)))
  .jp_logaddexp(A + lp, -A + lq) / psi
}

## Stable per-observation derivatives of the JP log-kernel. d, kappa, psi are
## recycled to a common length (vector d with scalar params for quadrature;
## all-vector at the data points). Returns the log-kernel h and first
## derivatives hphi/hk/hp; with second = TRUE adds the six second derivatives.
## Names and formulas:
##   T = tanh(arg), arg = A + log|cot(d/2)|;  hk = T,  hphi = -Vtilde sin d,
##   Vtilde = kappa sinhc(A)/g >= 0; the psi-direction uses the cumulant view
##   log g = K(A), with the |A| < A_SMALL series guarding the A -> 0 blow-up.
.jp_score_terms <- function(d, kappa, psi, second = TRUE) {
  n <- length(d)
  kappa <- rep_len(kappa, n)
  psi <- rep_len(psi, n)
  A <- kappa * psi
  absA <- abs(A)
  half <- 0.5 * d
  s <- sin(d); cc <- cos(d); s2 <- s * s
  lp <- 2 * log(abs(cos(half)))
  lq <- 2 * log(abs(sin(half)))
  logk <- log(kappa)
  logabss <- log(abs(s))
  lg <- .jp_logaddexp(A + lp, -A + lq)
  arg <- A + 0.5 * (lp - lq)
  T <- tanh(arg)
  a2 <- exp(-2 * abs(arg))
  sech2 <- 4 * a2 / (1 + a2)^2
  small <- absA < .JP_A_SMALL
  anysmall <- any(small)
  absA_s <- pmax(absA, .Machine$double.xmin)
  A2 <- A * A
  A2_s <- pmax(A2, .Machine$double.xmin)
  A3_s <- sign(A) * pmax(absA^3, .Machine$double.xmin)
  ## Each ifelse(small, A->0 series, main) below is rewritten as "main for all,
  ## then overwrite the small rows with the series" -- bit-for-bit, since the
  ## *_s pmax guards keep the main branch finite everywhere (the small rows are
  ## discarded). The series and its k3/k4/k5 intermediates are computed ONLY when
  ## some |A| is actually small (rare once psi has moved off 0), skipping ~8
  ## vector ops per call in the common case. (R has no block scope, so k3/k4/k5
  ## set in the guard below are visible to the second-derivative patches too.)
  if (anysmall) {
    k3 <- -2 * cc * s2
    k4 <- s2 * (4 - 6 * s2)
    k5 <- -8 * cc * s2 * (1 - 3 * s2)
  }
  lsinhc <- absA + log1p(-exp(-2 * absA_s)) - log(2 * absA_s)
  if (anysmall) lsinhc[small] <- (A * A / 6)[small]
  lVt <- pmin(logk + lsinhc - lg, .JP_LOG_CAP)
  Vt <- exp(lVt)
  sVt <- sign(s) * exp(logabss + lVt)   # Vtilde * sin d, 0*inf-safe
  h <- lg / psi
  if (anysmall)
    h[small] <- (kappa * (cc + s2 * A / 2 + k3 * A2 / 6 + k4 * A^3 / 24))[small]
  W <- (A * T - lg) / A2_s
  if (anysmall)
    W[small] <- (s2 / 2 + k3 * A / 3 + k4 * A2 / 8 + k5 * A^3 / 30)[small]
  out <- list(h = h, hphi = -sVt, hk = T, hp = kappa * kappa * W)
  if (!second) return(out)
  Tphi <- -sech2 / s
  s0 <- s == 0
  if (any(s0)) Tphi[s0] <- 0
  Wp <- (A2 * sech2 - 2 * A * T + 2 * lg) / A3_s
  if (anysmall) Wp[small] <- (k3 / 3 + k4 * A / 4 + k5 * A2 / 10)[small]
  Wphi <- (A * Tphi + psi * sVt) / A2_s
  if (anysmall)
    Wphi[small] <- (s * cc + 2 * s * (1 - 3 * cc * cc) * A / 3 +
                      s * cc * (1 - 3 * s2) * A2)[small]
  out$hphiphi <- -cc * Vt - psi * sVt * sVt
  out$hphik <- Tphi
  out$hphip <- kappa * kappa * Wphi
  out$hkk <- psi * sech2
  out$hkp <- kappa * sech2
  out$hpp <- kappa^3 * Wp
  out
}

## Scalar (kappa, psi) quadrature: the log normalizer log c = -log Z and the
## first/second (kappa, psi) derivatives of log Z, as kernel-weighted moments
## (one Gauss-Legendre node sweep serves all). The per-observation value
## shortcuts
## (uniform, von Mises) are applied by the callers, not here.
##   d log Z / d theta_a = E[h_a]
##   d2 log Z / d theta_a d theta_b = E[h_ab + h_a h_b] - E[h_a] E[h_b]
.jp_quad <- function(kappa, psi) {
  gp <- .jp_gl_panels(kappa, psi)
  nodes <- gp$nodes
  wts <- gp$weights
  t <- .jp_score_terms(nodes, kappa, psi, second = TRUE)
  hh <- t$h
  ## log c: peak value e^kappa factored out, full circle = 2 * half integral
  integral <- 2 * sum(wts * exp(hh - kappa))
  logc <- -(kappa + log(max(integral, .Machine$double.xmin)))
  ## log Z moments: e^{max h} factored out (cancels in the ratios)
  e <- exp(hh - max(hh)) * wts
  Z <- max(sum(e), .Machine$double.xmin)
  m <- function(v) sum(e * v) / Z
  dk <- m(t$hk)
  dp <- m(t$hp)
  dkk <- m(t$hkk + t$hk * t$hk) - dk * dk
  dkp <- m(t$hkp + t$hk * t$hp) - dk * dp
  dpp <- m(t$hpp + t$hp * t$hp) - dp * dp
  c(logc = logc, dk = dk, dp = dp, dkk = dkk, dkp = dkp, dpp = dpp)
}

## q-th cosine moment alpha_q = E[cos(q d)] of the centered JP law (the sine
## moments vanish, the kernel being even). Serves the Pearson residual scale
## Var{sin(d)} = (1 - alpha_2)/2.
.jp_cos_moment <- function(kappa, psi, q) {
  if (q == 0) return(1.0)
  if (kappa < .JP_KAPPA_TOL) return(0.0)
  if (abs(psi) < .JP_PSI_TOL)
    return(besselI(kappa, q, expon.scaled = TRUE) /
             besselI(kappa, 0, expon.scaled = TRUE))
  gp <- .jp_gl_panels(kappa, psi)
  nodes <- gp$nodes
  e <- gp$weights * exp(.jp_log_kernel(nodes, kappa, psi) - kappa)
  denom <- max(sum(e), .Machine$double.xmin)
  sum(cos(q * nodes) * e) / denom
}

## Per-observation log normalizer log c(kappa_i, psi_i) for the l0 VALUE, with
## value shortcuts: uniform (kappa < KAPPA_TOL -> -log 2pi), von Mises
## (|psi| < PSI_TOL -> -log(2pi I0(kappa))), else general quadrature. General
## rows are deduped to unique pairs, grouped by GL node count, and swept one
## size-group at a time -- an h-only (second = FALSE) sweep, since the value
## needs only the kernel. Bit-for-bit the per-pair .jp_quad logc.
.jp_logc_vec <- function(kappa, psi) {
  n <- length(kappa)
  out <- rep(-log(2 * pi), n)
  live <- kappa >= .JP_KAPPA_TOL
  vm <- live & (abs(psi) < .JP_PSI_TOL)
  if (any(vm))
    out[vm] <- -(log(2 * pi * besselI(kappa[vm], 0, expon.scaled = TRUE)) +
                   kappa[vm])
  gen <- which(live & !vm)
  if (length(gen)) {
    key <- paste(kappa[gen], psi[gen], sep = "\r")
    uk <- !duplicated(key)
    ukappa <- kappa[gen][uk]
    upsi <- psi[gen][uk]
    uvals <- numeric(length(ukappa))
    tiny <- .Machine$double.xmin
    for (grp in .jp_size_groups(ukappa, upsi)) {
      gi <- grp$idx
      k <- ukappa[gi]
      h <- .jp_score_terms(as.vector(grp$nodes), k, upsi[gi], second = FALSE)$h
      hmat <- matrix(h, nrow = length(gi))
      ## peak e^kappa factored out; full circle = 2 * half integral
      integral <- 2 * rowSums(grp$wts * exp(hmat - k))
      uvals[gi] <- -(k + log(pmax(integral, tiny)))
    }
    out[gen] <- uvals[match(key, key[uk])]
  }
  out
}

## Per-observation saturated reference for the sine-skewed family: the maximum
## over the circle of the normalizer-free centred log-density
## h(phi; kappa, psi) + log(1 + lambda sin phi). The skew tilts the mode off
## phi = 0, so unlike the symmetric families this peak is not h(0) = kappa and
## has no closed form; a fine grid suffices (deviance is a diagnostic, never
## parity-checked). Used by ssjplss residuals (deviance) and postproc.
.ssjp_lsat <- function(kappa, psi, lambda) {
  grid <- seq(-pi, pi, length.out = 721)
  sg <- sin(grid)
  vapply(seq_along(kappa), function(i) {
    h <- .jp_score_terms(grid, kappa[i], psi[i], second = FALSE)$h
    max(h + log1p(lambda[i] * sg))
  }, numeric(1))
}

## Per-observation general quadrature: an n x 6 matrix of (logc, dk, dp, dkk,
## dkp, dpp), evaluated once per unique (kappa, psi) pair. This is the single
## Gauss-Legendre sweep ll() needs: the value path takes column 1 (then
## overrides the uniform/von-Mises rows with closed forms), the derivative path
## takes columns 2:6 directly (the log-Z moments carry NO value shortcut -- the
## score stays live in the near-uniform corner). Unique pairs are grouped by GL
## node count and swept one size-group at a time (one .jp_score_terms call per
## distinct grid size, not per pair); bit-for-bit the per-pair .jp_quad.
.jp_quad_vec <- function(kappa, psi) {
  key <- paste(kappa, psi, sep = "\r")
  uk <- !duplicated(key)
  ukappa <- kappa[uk]
  upsi <- psi[uk]
  uvals <- matrix(0, nrow = length(ukappa), ncol = 6L)
  tiny <- .Machine$double.xmin
  for (grp in .jp_size_groups(ukappa, upsi)) {
    gi <- grp$idx
    g <- length(gi)
    L <- ncol(grp$nodes)
    k <- ukappa[gi]
    t <- .jp_score_terms(as.vector(grp$nodes), k, upsi[gi], second = TRUE)
    hmat <- matrix(t$h, g, L)
    ## logc: peak e^kappa factored out, full circle = 2 * half integral
    integral <- 2 * rowSums(grp$wts * exp(hmat - k))
    logc <- -(k + log(pmax(integral, tiny)))
    ## log Z moments: e^{row-wise max h} factored out (cancels in the ratios)
    rmax <- apply(hmat, 1L, max)
    e <- grp$wts * exp(hmat - rmax)
    Z <- pmax(rowSums(e), tiny)
    m <- function(v) rowSums(e * matrix(v, g, L)) / Z
    hk <- t$hk
    hp <- t$hp
    dk <- m(hk)
    dp <- m(hp)
    dkk <- m(t$hkk + hk * hk) - dk * dk
    dkp <- m(t$hkp + hk * hp) - dk * dp
    dpp <- m(t$hpp + hp * hp) - dp * dp
    uvals[gi, ] <- cbind(logc, dk, dp, dkk, dkp, dpp)
  }
  colnames(uvals) <- c("logc", "dk", "dp", "dkk", "dkp", "dpp")
  uvals[match(key, key[uk]), , drop = FALSE]
}

## ===========================================================================
## Asymmetric Jones-Pewsey extensions (ajplss). The asymmetry warps the JP
## kernel argument FORWARD, g(phi) = phi + nu cos phi (phi = y - xi), so the
## score reuses .jp_score_terms chained through g (g_phi = 1 - nu sin phi,
## g_nu = cos phi) -- NO implicit differentiation. Unlike ssjplss's sine-skew
## factor (which leaves the normalizer untouched), the warp MOVES Z, so c
## depends on (kappa, psi, nu) and the kappa/psi/nu log-Z moments come from
## .jp_logZ_moments_asym_vec. nu = 0 reduces to jplss exactly.
## ===========================================================================

## phi = g^{-1}(u) on the principal branch of g(phi) = phi + nu cos phi: 5
## bisection halvings warm a tight bracket, then 6 clamped-Newton steps (the
## closed-form g' = 1 - nu sin phi) polish to <=3e-14. The Newton step is clamped
## to the live bracket so it cannot run away where g' -> 0 (nu -> 1 near phi =
## pi/2). u and nu recycle to a common length (column-major per-row in the
## batched callers). Mirrors _jp_warp_inv.
.jp_warp_inv <- function(u, nu) {
  n <- length(u)
  nu <- rep_len(nu, n)
  a <- rep.int(-pi, n)
  b <- rep.int(pi, n)
  for (i in 1:5) {
    m <- 0.5 * (a + b)
    too_high <- (m + nu * cos(m)) > u
    b[too_high] <- m[too_high]
    a[!too_high] <- m[!too_high]
  }
  phi <- 0.5 * (a + b)
  for (i in 1:6)
    phi <- pmin.int(pmax.int(phi - (phi + nu * cos(phi) - u) / (1 - nu * sin(phi)), a), b)
  phi
}

## Break-point ladder in the kernel's own angle u = g(phi) over one period
## [-pi - nu, pi - nu]: decade rungs at the .jp_feature_scales of the peak
## (u = 0) and the antipodal near-kink (u ~ +-pi), plus rungs for the warp-weight
## 1/g''s bump at u = g(+-pi/2) = +-pi/2 of u-width ~(1-|nu|)^{3/2}. Candidates
## outside the window are dropped. Mirrors _jp_ladder_edges_asym.
.jp_ladder_edges_asym <- function(kappa, psi, nu) {
  two_pi <- 2 * pi
  lo <- -pi - nu
  hi <- pi - nu
  ws <- .jp_feature_scales(kappa, psi)
  cand <- numeric(0)
  r <- ws[1]
  while (r < two_pi) { cand <- c(cand, -r, r); r <- r * 10 }
  r <- ws[2]
  while (r < two_pi) { cand <- c(cand, -pi - r, -pi + r, pi - r, pi + r); r <- r * 10 }
  r <- max((1 - abs(nu))^1.5, 1e-3)
  while (r < two_pi) {
    cand <- c(cand, 0.5 * pi - r, 0.5 * pi + r, -0.5 * pi - r, -0.5 * pi + r)
    r <- r * 10
  }
  cand <- cand[cand > lo & cand < hi]
  sort.int(unique(c(lo, hi, 0.0, cand)))
}

## Unique (kappa, psi, nu) asymmetric GL grids grouped by node count for batched
## evaluation -- the asym twin of .jp_size_groups (full window, not half-circle;
## each tuple keeps its own adaptive _asym ladder, so a batched sweep is
## bit-for-bit the per-tuple scalar). Returns list(idx, nodes = g x L, wts = g x L).
.jp_asym_size_groups <- function(kappa, psi, nu) {
  grids <- lapply(seq_along(kappa), function(i)
    .jp_panels_from_edges(.jp_ladder_edges_asym(kappa[i], psi[i], nu[i])))
  sizes <- vapply(grids, function(g) length(g$nodes), integer(1))
  lapply(sort.int(unique(sizes)), function(L) {
    idx <- which(sizes == L)
    nodes <- matrix(unlist(lapply(grids[idx], `[[`, "nodes"), use.names = FALSE),
                    nrow = length(idx), ncol = L, byrow = TRUE)
    wts <- matrix(unlist(lapply(grids[idx], `[[`, "weights"), use.names = FALSE),
                  nrow = length(idx), ncol = L, byrow = TRUE)
    list(idx = idx, nodes = nodes, wts = wts)
  })
}

## Per-observation log normalizer log c(kappa_i, psi_i, nu_i) of the asymmetric
## kernel, by the substitution u = g(phi): integrate the JP kernel over the asym
## ladder times the warp Jacobian weight 1/g'(g^{-1}(u)) = 1/(1 - nu sin phi).
## kappa < KAPPA_TOL -> uniform. Unique triples grouped by node count, one h-only
## sweep per group; bit-for-bit the per-triple scalar. Mirrors _jp_log_c_asym_vec.
.jp_log_c_asym_vec <- function(kappa, psi, nu) {
  len <- max(length(kappa), length(psi), length(nu))
  kappa <- rep_len(kappa, len)
  psi <- rep_len(psi, len)
  nu <- rep_len(nu, len)
  out <- rep(-log(2 * pi), len)
  live <- which(kappa >= .JP_KAPPA_TOL)
  if (!length(live)) return(out)
  key <- paste(kappa[live], psi[live], nu[live], sep = "\r")
  uk <- !duplicated(key)
  ukappa <- kappa[live][uk]
  upsi <- psi[live][uk]
  unu <- nu[live][uk]
  uvals <- numeric(length(ukappa))
  tiny <- .Machine$double.xmin
  for (grp in .jp_asym_size_groups(ukappa, upsi, unu)) {
    gi <- grp$idx
    g <- length(gi)
    L <- ncol(grp$nodes)
    k <- ukappa[gi]
    nn <- unu[gi]
    phi <- matrix(.jp_warp_inv(as.vector(grp$nodes), nn), g, L)
    weight <- 1 / (1 - nn * sin(phi))
    h <- matrix(.jp_score_terms(as.vector(grp$nodes), k, upsi[gi], second = FALSE)$h, g, L)
    integral <- rowSums(grp$wts * exp(h - k) * weight)
    uvals[gi] <- -(k + log(pmax(integral, tiny)))
  }
  out[live] <- uvals[match(key, key[uk])]
  out
}

## Per-observation asymmetric JP log-density (the l0 VALUE): the stable log-kernel
## at the warped angle g = phi + nu cos phi plus the u-substituted normalizer,
## each datum its own (xi, kappa, psi, nu). kappa < KAPPA_TOL -> uniform.
## Mirrors _ajp_logpdf_vec.
.ajp_logpdf_vec <- function(x, xi, kappa, psi, nu) {
  len <- max(length(x), length(xi), length(kappa), length(psi), length(nu))
  x <- rep_len(x, len)
  xi <- rep_len(xi, len)
  kappa <- rep_len(kappa, len)
  psi <- rep_len(psi, len)
  nu <- rep_len(nu, len)
  phi <- x - xi
  g <- phi + nu * cos(phi)
  h <- .jp_score_terms(g, kappa, psi, second = FALSE)$h
  h[kappa < .JP_KAPPA_TOL] <- 0
  h + .jp_log_c_asym_vec(kappa, psi, nu)
}

## Per-observation FUSED asymmetric quadrature: an n x 10 matrix of (logc, dk,
## dp, dnu, dkk, dkp, dpp, dknu, dpnu, dnunu), evaluated in ONE Gauss-Legendre
## sweep per unique (kappa, psi, nu) triple -- the asym twin of .jp_quad_vec.
## The single u = g(phi) integral Z = sum(w e^{h-kappa}/g') yields BOTH the
## normalizer logc = -(kappa + log Z) (column 1) and the kernel-weighted log-Z
## moments (columns 2:10), so the derivative hot path no longer integrates the
## normalizer twice (.jp_log_c_asym_vec + a separate moment sweep). By the warp
## chain rule h_nu = h_phi cos phi, h_knu = h_phik cos phi, h_psinu = h_phip cos
## phi, h_nunu = h_phiphi cos^2 phi (with h_. the .jp_score_terms derivatives at
## the warped argument). kappa < KAPPA_TOL -> uniform: logc = -log 2pi, moments 0
## (the scalar contract; the score stays live in the near-uniform corner). The
## logc column is bit-for-bit .jp_log_c_asym_vec (same Z); the moment columns are
## bit-for-bit the old standalone sweep. This is the ajplss deriv = 1 hot path.
.jp_quad_asym_vec <- function(kappa, psi, nu) {
  len <- max(length(kappa), length(psi), length(nu))
  kappa <- rep_len(kappa, len)
  psi <- rep_len(psi, len)
  nu <- rep_len(nu, len)
  out <- cbind(rep(-log(2 * pi), len), matrix(0, len, 9L))
  live <- which(kappa >= .JP_KAPPA_TOL)
  if (!length(live)) return(out)
  key <- paste(kappa[live], psi[live], nu[live], sep = "\r")
  uk <- !duplicated(key)
  ukappa <- kappa[live][uk]
  upsi <- psi[live][uk]
  unu <- nu[live][uk]
  uvals <- matrix(0, length(ukappa), 10L)
  tiny <- .Machine$double.xmin
  for (grp in .jp_asym_size_groups(ukappa, upsi, unu)) {
    gi <- grp$idx
    g <- length(gi)
    L <- ncol(grp$nodes)
    k <- ukappa[gi]
    p <- upsi[gi]
    nn <- unu[gi]
    nodes <- as.vector(grp$nodes)
    phi <- matrix(.jp_warp_inv(nodes, nn), g, L)
    cphi <- cos(phi)
    weight <- 1 / (1 - nn * sin(phi))
    t <- .jp_score_terms(nodes, k, p, second = TRUE)
    hk <- matrix(t$hk, g, L)
    hp <- matrix(t$hp, g, L)
    hphi <- matrix(t$hphi, g, L)
    e <- grp$wts * exp(matrix(t$h, g, L) - k) * weight
    Z <- pmax(rowSums(e), tiny)
    logc <- -(k + log(Z))
    m <- function(v) rowSums(e * v) / Z
    hnu <- hphi * cphi
    dk <- m(hk)
    dp <- m(hp)
    dnu <- m(hnu)
    dkk <- m(matrix(t$hkk, g, L) + hk * hk) - dk * dk
    dkp <- m(matrix(t$hkp, g, L) + hk * hp) - dk * dp
    dpp <- m(matrix(t$hpp, g, L) + hp * hp) - dp * dp
    dknu <- m(matrix(t$hphik, g, L) * cphi + hk * hnu) - dk * dnu
    dpnu <- m(matrix(t$hphip, g, L) * cphi + hp * hnu) - dp * dnu
    dnunu <- m(matrix(t$hphiphi, g, L) * cphi * cphi + hnu * hnu) - dnu * dnu
    uvals[gi, ] <- cbind(logc, dk, dp, dnu, dkk, dkp, dpp, dknu, dpnu, dnunu)
  }
  out[live, ] <- uvals[match(key, key[uk]), , drop = FALSE]
  out
}

## Per-observation first/second (kappa, psi, nu) log-Z moments of the asymmetric
## kernel (an n x 9 matrix: dk, dp, dnu, dkk, dkp, dpp, dknu, dpnu, dnunu) -- the
## moment columns of the fused .jp_quad_asym_vec sweep. A named entry point for
## callers that want moments without the normalizer. Mirrors
## _jp_logZ_moments_asym_vec.
.jp_logZ_moments_asym_vec <- function(kappa, psi, nu)
  .jp_quad_asym_vec(kappa, psi, nu)[, -1L, drop = FALSE]

## Per-observation saturated reference for the asymmetric family: the maximum
## over the circle of the normalizer-free centred log-kernel h(g(phi)) at the
## warped angle g = phi + nu cos phi. The warp tilts the mode off phi = 0, so the
## peak has no closed form; a fine grid suffices (deviance is a diagnostic, never
## parity-checked). Used by ajplss residuals (deviance) and postproc.
.ajp_lsat <- function(kappa, psi, nu) {
  grid <- seq(-pi, pi, length.out = 721)
  vapply(seq_along(kappa), function(i) {
    g <- grid + nu[i] * cos(grid)
    max(.jp_score_terms(g, kappa[i], psi[i], second = FALSE)$h)
  }, numeric(1))
}
