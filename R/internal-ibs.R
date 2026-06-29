## Inverse Batschelet internal numerics: the two inverse angular warps, the grid
## normalizer log c(kappa, lambda), its finite-differenced (kappa, lambda)
## gradient/Hessian, the analytic log-kernel gradient, and the assembled
## (xi, kappa, nu, lambda) score/Hessian at the data.
##
## The inverse Batschelet density (Jones & Pewsey 2012) is
##   f(theta) = c(kappa, lambda) exp(kappa cos A),    phi = (theta - xi) wrapped,
## where the angle is warped FORWARD into the von Mises kernel by two inverse
## maps: a skewness warp phi* = t_nu^{-1}(phi) solving y - nu(1 + cos y) = phi,
## then a peakedness warp u* = s_lambda^{-1}(phi*) solving u - (1+lambda)/2 sin u
## = phi*, and A = u* - (1-lambda)/2 sin u* = a phi* + b u* with a=(1-lambda)/
## (1+lambda), b=2 lambda/(1+lambda) (the two forms agree because a+b=1 and
## a (1+lambda)/2 = (1-lambda)/2). Because the warps are INVERSE maps the score
## needs implicit differentiation (unlike vmftlss's forward warp): the
## per-observation kernel gradient reuses the solved roots phi*, u* and the warp
## slopes T_nu = 1 + nu sin phi*, S_lambda = 1 - (1+lambda)/2 cos u*,
## B' = 1 - (1-lambda)/2 cos u*. nu = lambda = 0 reduces to vmlss; kappa -> 0 to
## the circular uniform.
##
## The normalizer c(kappa, lambda) = (1-lambda)/[(1+lambda) 2pi I0(kappa)
## - 2 lambda J], J = int exp(kappa cos B(phi)) dphi with B(phi) = phi
## - (1-lambda)/2 sin phi, has no elementary form: J is an equispaced-trapezoid
## grid integral (the warped kernel is smooth and periodic, so the trapezoid is
## spectrally accurate) and the kappa/lambda derivatives of log c are CENTRAL
## FINITE DIFFERENCES of log c -- on the coarser .IBS_DERIV_GRID for the
## derivatives, the full .IBS_NUMERIC_GRID for the value, with grids and steps
## chosen so the score stays consistent with the logpdf the derivative tests
## difference. The kernel block of the Hessian is a central FD
## of the ANALYTIC kernel gradient (FD of a closed form, not FD-of-FD); the
## normalizer block adds the FD'd log c second derivatives only to the kappa,kappa
## / kappa,lambda / lambda,lambda pairs (c is free of xi and nu). These functions
## are package-internal; the family closures reach them by lexical scope,
## postproc by getFromNamespace.

## Tolerances, clamps and grids.
.IBS_KAPPA_TOL <- 1e-9       # below this kappa the density is the uniform 1/2pi
.IBS_KAPPA_UPPER <- 700.0    # kappa clamp (besselI / grid stability)
.IBS_LMBDA_TOL <- 1e-12      # |lambda| within this of 0/1 -> exact-limit edge
.IBS_NUMERIC_GRID <- 4096L   # trapezoid grid for the log c VALUE
.IBS_DERIV_GRID <- 256L      # coarser grid for the FD'd log c derivatives
.IBS_MOMENT_GRID <- 2048L    # grid for the Pearson cosine moment (diagnostic)
.IBS_NEWTON_TOL <- 1e-14     # monotone-solver Newton step tolerance
.IBS_NEWTON_MAXITER <- 20L   # Newton iterations before the bisection mop-up
.IBS_BISECT_ITER <- 60L      # bisection iterations for the near-boundary points

## Vectorized root of a smooth, monotone-increasing g(y) = rhs on [lo, hi] with
## g(lo) <= rhs <= g(hi). Newton from x0 (default rhs), each step clamped into
## the bracket, then a vectorized bisection mop-up for any point Newton leaves
## with residual above 1e-12 (the near-boundary cases where gprime -> 0 as
## nu, lambda -> +-1). Mirrors _solve_monotone_increasing; replaces a per-point
## brentq loop with a handful of array ops. Returns the root vector (the slope
## the Python returns alongside is unused on the regression path -- the kernel
## gradient recomputes T_nu, S_lambda directly).
.ibs_solve_monotone <- function(rhs, g, gprime, lo = -pi, hi = pi, x0 = NULL,
                                tol = .IBS_NEWTON_TOL,
                                max_iter = .IBS_NEWTON_MAXITER) {
  rhs <- as.numeric(rhs)
  y <- if (is.null(x0)) pmin.int(pmax.int(rhs, lo), hi)
       else pmin.int(pmax.int(rep_len(as.numeric(x0), length(rhs)), lo), hi)
  if (length(y)) {
    for (it in seq_len(max_iter)) {
      gp <- gprime(y)
      ## ifelse(abs(gp) > 1e-15, (g(y)-rhs)/gp, 0) as main-then-patch; gprime is
      ## the (finite, positive) monotone-warp slope, so the condition has no NA.
      step <- (g(y) - rhs) / gp
      step[abs(gp) <= 1e-15] <- 0.0
      y_new <- pmin.int(pmax.int(y - step, lo), hi)
      if (max(abs(y_new - y)) < tol) { y <- y_new; break }
      y <- y_new
    }
    bad <- abs(g(y) - rhs) > 1e-12
    if (any(bad)) {
      a <- rep(lo, length(y)); b <- rep(hi, length(y))
      for (it in seq_len(.IBS_BISECT_ITER)) {
        m <- 0.5 * (a + b)
        left <- g(m) <= rhs            # g increasing -> root at/above m
        a[left] <- m[left]
        b[!left] <- m[!left]
      }
      y[bad] <- (0.5 * (a + b))[bad]
    }
  }
  y
}

## Per-observation inverse warps for the regression path: phi* = t_nu^{-1}(phi)
## then u* = s_lambda^{-1}(phi*), with phi = (x - xi) wrapped and every parameter
## an array (each datum its own nu_i, lambda_i). Calls the shared monotone solver
## with array-valued closures. Returns (phi*, u*), both wrapped to [-pi, pi).
## Mirrors _invbat_warp_vec.
.ibs_warp_vec <- function(x, xi, nu, lmbd) {
  L <- max(length(x), length(xi), length(nu), length(lmbd))
  x <- rep_len(as.numeric(x), L); xi <- rep_len(as.numeric(xi), L)
  nu <- rep_len(as.numeric(nu), L); lmbd <- rep_len(as.numeric(lmbd), L)
  phi <- ((x - xi + pi) %% (2.0 * pi)) - pi
  root1 <- .ibs_solve_monotone(phi,
    function(yy) yy - nu * (1.0 + cos(yy)),
    function(yy) 1.0 + nu * sin(yy),
    x0 = phi)
  phi_star <- ((root1 + pi) %% (2.0 * pi)) - pi
  cc <- 0.5 * (1.0 + lmbd)
  root2 <- .ibs_solve_monotone(phi_star,
    function(uu) uu - cc * sin(uu),
    function(uu) 1.0 - cc * cos(uu),
    x0 = phi_star)
  u_star <- ((root2 + pi) %% (2.0 * pi)) - pi
  list(phi_star = phi_star, u_star = u_star)
}

## The four analytic log-KERNEL gradient components (no normalizer), as an
## n x 4 matrix in book order [xi, kappa, nu, lambda]. Reuses the two warps'
## solved roots phi*, u* and the recomputed warp slopes. Isolated so the Hessian
## can central-difference it for the kernel block (FD of an analytic gradient,
## not FD-of-FD). kappa is NOT clamped here -- the kernel gradient is exact in
## kappa.
.ibs_dlogkernel <- function(x, xi, kappa, nu, lmbd) {
  L <- max(length(x), length(xi), length(kappa), length(nu), length(lmbd))
  x <- rep_len(as.numeric(x), L); xi <- rep_len(as.numeric(xi), L)
  kappa <- rep_len(as.numeric(kappa), L); nu <- rep_len(as.numeric(nu), L)
  lmbd <- rep_len(as.numeric(lmbd), L)
  w <- .ibs_warp_vec(x, xi, nu, lmbd)
  phi_star <- w$phi_star; u_star <- w$u_star
  half <- 0.5 * (1.0 - lmbd)
  A <- u_star - half * sin(u_star)
  Bp <- 1.0 - half * cos(u_star)                 # B'(u*)
  Tnu <- 1.0 + nu * sin(phi_star)                # t_nu warp slope
  Slam <- 1.0 - 0.5 * (1.0 + lmbd) * cos(u_star) # s_lambda warp slope
  S <- kappa * sin(A)
  chain <- Bp / (Tnu * Slam)
  d_xi <- S * chain
  d_nu <- -S * chain * (1.0 + cos(phi_star))
  d_lmbd <- -S * 0.5 * sin(u_star) * (1.0 + Bp / Slam)
  d_kappa <- cos(A)
  cbind(xi = d_xi, kappa = d_kappa, nu = d_nu, lmbd = d_lmbd)
}

## log J = log int_{-pi}^{pi} exp(kappa cos B(phi)) dphi, B(phi) = phi
## - (1-lambda)/2 sin phi, by an overflow-safe (max-subtracted) trapezoid --
## the scalar form used only by the exact-limit edge fallback. Mirrors
## _log_invbatschelet_kernel_integral.
.ibs_kernel_log_integral <- function(kappa, lmbd, grid_size) {
  phi <- seq(-pi, pi, length.out = grid_size + 1L)
  log_kernel <- kappa * cos(phi - 0.5 * (1.0 - lmbd) * sin(phi))
  max_log <- max(log_kernel)
  w <- rep(1.0, grid_size + 1L); w[1L] <- 0.5; w[grid_size + 1L] <- 0.5
  log(2.0 * pi / grid_size) + max_log + log(sum(w * exp(log_kernel - max_log)))
}

## log c(kappa, lambda) by the overflow-safe denominator split D = (1+lambda)
## 2pi I0 - 2 lambda J, on grid_size nodes -- the scalar exact-limit form
## (kappa<=tol uniform, lambda=-1 -> -log J, lambda=+1 -> the I0/I1 closed form
## handling the 0/0). Mirrors _c_invbatschelet (returns log c). Only the rare
## edges of _invbat_log_c_array defer here; the log/tanh links keep fits in the
## interior where the vectorized form is used directly.
.ibs_log_c_numeric <- function(k, l, grid_size) {
  log_int <- .ibs_kernel_log_integral(k, l, grid_size)
  log_mult <- log(2.0 * pi) + log(besselI(k, 0, expon.scaled = TRUE)) + k
  log_t1 <- log1p(l) + log_mult
  log_t2 <- if (abs(l) <= .IBS_LMBDA_TOL) -Inf else log(2.0 * abs(l)) + log_int
  m <- max(log_t1, log_t2)
  t1 <- exp(log_t1 - m)
  t2 <- if (is.finite(log_t2)) exp(log_t2 - m) else 0.0
  denom <- if (l >= 0.0) t1 - t2 else t1 + t2
  if (denom <= 0.0 || !is.finite(denom)) return(NaN)
  log1p(-l) - (m + log(denom))
}

.ibs_log_c_scalar <- function(kappa, lmbd) {
  k <- min(max(kappa, 0.0), .IBS_KAPPA_UPPER)
  l <- lmbd
  if (k <= .IBS_KAPPA_TOL) return(-log(2.0 * pi))
  if (abs(l - 1.0) <= .IBS_LMBDA_TOL) {
    log_mult <- log(2.0 * pi) + log(besselI(k, 0, expon.scaled = TRUE)) + k
    A1 <- besselI(k, 1, expon.scaled = TRUE) / besselI(k, 0, expon.scaled = TRUE)
    K <- 1.0 - A1
    if (!is.finite(K) || K <= 0.0)
      return(.ibs_log_c_numeric(k, l, .IBS_NUMERIC_GRID * 2L))
    return(-log_mult - log(K))
  }
  if (abs(l + 1.0) <= .IBS_LMBDA_TOL)
    return(-.ibs_kernel_log_integral(k, l, .IBS_NUMERIC_GRID))
  lc <- .ibs_log_c_numeric(k, l, .IBS_NUMERIC_GRID)
  if (!is.finite(lc)) lc <- .ibs_log_c_numeric(k, l, .IBS_NUMERIC_GRID * 2L)
  lc
}

## Vectorized log c(kappa, lambda) for arrays on one shared grid -- a faithful,
## overflow-safe vectorization of the denominator split over the interior
## (kappa>0, |lambda|<1) the log/tanh links guarantee. Edge pairs (kappa~0,
## |lambda|~1, or any non-finite result) defer to the scalar form for
## exact-limit parity. Mirrors _invbat_log_c_array.
.ibs_log_c_array <- function(kappa, lmbd, grid_size) {
  L <- max(length(kappa), length(lmbd))
  kappa <- rep_len(as.numeric(kappa), L); lmbd <- rep_len(as.numeric(lmbd), L)
  kf <- pmin(pmax(kappa, 0.0), .IBS_KAPPA_UPPER)
  lf <- lmbd
  phi <- seq(-pi, pi, length.out = grid_size + 1L)
  sin_phi <- sin(phi)
  w <- rep(1.0, grid_size + 1L); w[1L] <- 0.5; w[grid_size + 1L] <- 0.5
  B <- sweep(outer(-0.5 * (1.0 - lf), sin_phi), 2L, phi, "+")  # L x (G+1)
  log_kernel <- kf * cos(B)                       # kf recycles per row (col-major)
  max_log <- apply(log_kernel, 1L, max)
  weighted_sum <- drop(exp(log_kernel - max_log) %*% w)
  log_int <- log(2.0 * pi / grid_size) + max_log + log(weighted_sum)
  log_mult <- log(2.0 * pi) + log(besselI(kf, 0, expon.scaled = TRUE)) + kf
  log_t1 <- log1p(lf) + log_mult
  log_t2 <- log(2.0 * abs(lf)) + log_int          # -Inf at lambda = 0 -> term 0
  m <- pmax(log_t1, log_t2)
  t1 <- exp(log_t1 - m)
  t2 <- ifelse(is.finite(log_t2), exp(log_t2 - m), 0.0)
  denom_scaled <- ifelse(lf >= 0.0, t1 - t2, t1 + t2)
  log_c <- log1p(-lf) - (m + log(denom_scaled))
  edge <- (kf <= .IBS_KAPPA_TOL) | (abs(abs(lf) - 1.0) <= .IBS_LMBDA_TOL) |
    !is.finite(log_c) | (denom_scaled <= 0.0)
  if (any(edge)) {
    for (idx in which(edge)) log_c[idx] <- .ibs_log_c_scalar(kf[idx], lf[idx])
  }
  log_c
}

## Shared unique-(kappa, lambda) collapse for the value/gradient/Hessian: an
## intercept-only or repeated-covariate fit collapses to a handful of grid
## passes. Returns the unique pair vectors and the map back to the input length.
.ibs_unique_pairs <- function(kappa, lmbd) {
  L <- max(length(kappa), length(lmbd))
  kappa <- rep_len(as.numeric(kappa), L); lmbd <- rep_len(as.numeric(lmbd), L)
  key <- paste(kappa, lmbd, sep = "\r")
  uk <- !duplicated(key); uidx <- which(uk)
  list(k = kappa[uidx], l = lmbd[uidx], idx = match(key, key[uk]))
}

## Vectorized log c over per-observation params, once per unique (kappa, lambda)
## at the full value grid -- the normalizer for the regression-path l0. Mirrors
## _invbat_log_c_vec.
.ibs_log_c_vec <- function(kappa, lmbd) {
  up <- .ibs_unique_pairs(kappa, lmbd)
  vals <- .ibs_log_c_array(up$k, up$l, grid_size = .IBS_NUMERIC_GRID)
  vals[up$idx]
}

## (d log c/d kappa, d log c/d lambda) by central finite-difference of
## .ibs_log_c_array on the coarse derivative grid, for all unique (kappa, lambda)
## pairs at once. The lambda step shrinks near the boundary.
.ibs_logc_grad_vec <- function(kappa, lmbd) {
  up <- .ibs_unique_pairs(kappa, lmbd)
  k <- up$k; lm <- up$l
  hk <- 1e-6 * pmax(1.0, abs(k))
  hl <- pmin(1e-6, 0.25 * (1.0 - abs(lm)))
  lc <- function(kk, ll) .ibs_log_c_array(kk, ll, grid_size = .IBS_DERIV_GRID)
  gk <- (lc(k + hk, lm) - lc(k - hk, lm)) / (2.0 * hk)
  gl <- (lc(k, lm + hl) - lc(k, lm - hl)) / (2.0 * hl)
  list(dk = gk[up$idx], dl = gl[up$idx])
}

## (d2 log c/d kappa2, d2 log c/d kappa d lambda, d2 log c/d lambda2) by direct
## second central differences of .ibs_log_c_array on the derivative grid -- the
## normalizer block of the Hessian. Direct second differences (one FD level on
## the smooth log c) avoid the FD-of-FD noise of differencing the FD'd gradient.
## Larger steps than the gradient. Mirrors _invbat_logc_hess_vec.
.ibs_logc_hess_vec <- function(kappa, lmbd) {
  up <- .ibs_unique_pairs(kappa, lmbd)
  k <- up$k; lm <- up$l
  hk <- 1e-4 * pmax(1.0, abs(k))
  hl <- pmin(1e-4, 0.25 * (1.0 - abs(lm)))
  lc <- function(kk, ll) .ibs_log_c_array(kk, ll, grid_size = .IBS_DERIV_GRID)
  f0 <- lc(k, lm)
  hkk <- (lc(k + hk, lm) - 2.0 * f0 + lc(k - hk, lm)) / (hk * hk)
  hll <- (lc(k, lm + hl) - 2.0 * f0 + lc(k, lm - hl)) / (hl * hl)
  hkl <- (lc(k + hk, lm + hl) - lc(k + hk, lm - hl) -
            lc(k - hk, lm + hl) + lc(k - hk, lm - hl)) / (4.0 * hk * hl)
  list(hkk = hkk[up$idx], hkl = hkl[up$idx], hll = hll[up$idx])
}

## q-th cosine moment alpha_q = E[cos(q phi)] of the centered inverse-Batschelet
## law, a grid moment under the warped kernel on .IBS_MOMENT_GRID. Serves the
## Pearson residual scale Var{sin phi} ~ (1 - alpha_2)/2 (a diagnostic
## standardization for the skewed law, as in ssjplss; xi is the mode anchor when
## nu != 0). Mirrors the .vmft_cos_moment role.
.ibs_cos_moment <- function(kappa, nu, lmbd, q) {
  if (q == 0) return(1.0)
  if (kappa < .IBS_KAPPA_TOL) return(0.0)
  grid <- .IBS_MOMENT_GRID
  phi <- seq(-pi, pi, length.out = grid + 1L)
  wt <- rep(1.0, grid + 1L); wt[1L] <- 0.5; wt[grid + 1L] <- 0.5
  w <- .ibs_warp_vec(phi, 0.0, nu, lmbd)          # xi = 0 -> phi unchanged
  A <- w$u_star - 0.5 * (1.0 - lmbd) * sin(w$u_star)
  lk <- kappa * cos(A)
  e <- wt * exp(lk - max(lk))
  denom <- max(sum(e), .Machine$double.xmin)
  sum(cos(q * phi) * e) / denom
}

## The inverse-Batschelet log-density value and its (xi, kappa, nu, lambda)
## score/Hessian at the data -- the single entry point the family ll() and the
## standalone validator share (the .vmft_terms convention). l0 reproduces the
## class _logpdf (clamped kappa, uniform shortcut); l1 the dlogpdf (analytic
## kernel gradient + FD'd normalizer gradient on kappa, lambda); l2 the d2logpdf
## (central FD of the analytic kernel gradient over all four params, symmetrized,
## plus the FD'd log c second derivatives on the kappa,kappa / kappa,lambda /
## lambda,lambda pairs). l2 columns are the 10 unique unordered pairs in
## combinations_with_replacement order (xi, kappa, nu, lambda):
## (xx, xk, xn, xl, kk, kn, kl, nn, nl, ll). The kernel gradient/Hessian take the
## RAW kappa (exact in kappa); only l0 and the normalizer clamp it. Returns l0,
## and (deriv) l1 [n x 4], l2 [n x 10].
.ibs_terms <- function(y, xi, kappa, nu, lmbd, deriv = TRUE) {
  n <- length(y)
  xi <- rep_len(as.numeric(xi), n); kappa <- rep_len(as.numeric(kappa), n)
  nu <- rep_len(as.numeric(nu), n); lmbd <- rep_len(as.numeric(lmbd), n)
  kc <- pmin(pmax(kappa, 0.0), .IBS_KAPPA_UPPER)
  w <- .ibs_warp_vec(y, xi, nu, lmbd)
  A <- w$u_star - 0.5 * (1.0 - lmbd) * sin(w$u_star)
  l0 <- kc * cos(A) + .ibs_log_c_vec(kc, lmbd)
  unif <- kc <= .IBS_KAPPA_TOL
  if (any(unif)) l0[unif] <- -log(2.0 * pi)
  out <- list(l0 = l0)
  if (deriv) {
    g0 <- .ibs_dlogkernel(y, xi, kappa, nu, lmbd)   # raw kappa
    gc <- .ibs_logc_grad_vec(kappa, lmbd)
    l1 <- cbind(g0[, 1L],
                g0[, 2L] + gc$dk,
                g0[, 3L],
                g0[, 4L] + gc$dl)
    ## kernel Hessian: central FD of the analytic kernel gradient over each
    ## param, H[, a, b] = d(grad_a)/d(param_b).
    base <- list(xi, kappa, nu, lmbd)
    steps <- list(1e-6, 1e-5 * pmax(1.0, abs(kappa)), 1e-6, 1e-6)
    H <- array(0.0, dim = c(n, 4L, 4L))
    for (b in 1:4) {
      hb <- steps[[b]]
      pp <- base; pm <- base
      pp[[b]] <- base[[b]] + hb
      pm[[b]] <- base[[b]] - hb
      gp <- .ibs_dlogkernel(y, pp[[1L]], pp[[2L]], pp[[3L]], pp[[4L]])
      gm <- .ibs_dlogkernel(y, pm[[1L]], pm[[2L]], pm[[3L]], pm[[4L]])
      H[, , b] <- (gp - gm) / (2.0 * hb)
    }
    H <- 0.5 * (H + aperm(H, c(1L, 3L, 2L)))         # symmetrize
    hc <- .ibs_logc_hess_vec(kappa, lmbd)
    l2 <- cbind(H[, 1L, 1L], H[, 1L, 2L], H[, 1L, 3L], H[, 1L, 4L],
                H[, 2L, 2L] + hc$hkk, H[, 2L, 3L], H[, 2L, 4L] + hc$hkl,
                H[, 3L, 3L], H[, 3L, 4L], H[, 4L, 4L] + hc$hll)
    colnames(l1) <- c("xi", "kappa", "nu", "lmbd")
    colnames(l2) <- c("xx", "xk", "xn", "xl", "kk", "kn", "kl", "nn", "nl", "ll")
    out$l1 <- l1
    out$l2 <- l2
  }
  out
}
