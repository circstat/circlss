## Flat-topped von Mises internal numerics: the grid normalizer log c(kappa,nu)
## and the first/second (kappa,nu)-derivatives of log Z as kernel-weighted grid
## moments.
##
## The flat-topped von Mises density is
##   f(theta) = c(kappa, nu) * exp(kappa cos(phi + nu sin phi)),   phi = theta - mu,
## the peakedness factor B = phi + nu sin phi warping the angle FORWARD (B is odd
## in phi, so the density stays symmetric -- nu is a peakedness/flat-top knob,
## nu>0 sharper, nu<0 flatter, NOT a skew one, and the score needs no implicit
## differentiation, unlike inverse_batschelet). The normalizer Z(kappa,nu) = int
## exp(kappa cos B) dphi has no elementary form; it and the kappa/nu derivatives
## of log Z = -log c are equispaced-trapezoid grid moments under the density (the
## jplss `ell_kappa = h_kappa - E[h_kappa]` pattern, but on a uniform grid rather
## than a Gauss-Legendre ladder -- the warped kernel is smooth and periodic, so
## the trapezoid is spectrally accurate). nu=0 reduces to plain vmlss; kappa->0
## to the circular uniform. These functions are package-internal; the family
## closures reach them by lexical scope, postproc by getFromNamespace.

## Tolerances, clamps and grid bounds.
.VMFT_KAPPA_TOL <- 1e-9     # below this kappa the density is the uniform 1/2pi
.VMFT_KAPPA_UPPER <- 1e3    # kappa clamp (numerical stability of the grid)
.VMFT_GRID_BASE <- 64.0     # additive base in the adaptive grid-size formula
.VMFT_GRID_SHARPNESS <- 12.0
.VMFT_MIN_GRID <- 512L
.VMFT_MAX_GRID <- 8192L

## Adaptive trapezoid grid size: ~ sqrt(kappa) (1 + |nu|), clamped to
## [512, 8192] and rounded UP to a power of two (even by construction). Higher
## concentration / |nu| -> a sharper peak -> more nodes. Mirrors _vmft_grid_size.
.vmft_grid_size <- function(kappa, nu) {
  sharpness <- (1.0 + abs(nu)) * sqrt(max(kappa, 0.0) + 1.0)
  target <- .VMFT_GRID_BASE + .VMFT_GRID_SHARPNESS * sharpness
  target <- min(max(target, .VMFT_MIN_GRID), .VMFT_MAX_GRID)
  power <- ceiling(log2(target))
  size <- bitwShiftL(1L, as.integer(power))
  size <- min(max(size, .VMFT_MIN_GRID), .VMFT_MAX_GRID)
  if (size %% 2L != 0L) size <- size + 1L
  as.integer(size)
}

## Vectorized log c(kappa, nu) = -log int exp(kappa cos(phi + nu sin phi)) dphi
## over per-observation params, once per unique (kappa, nu) on ONE shared grid
## sized at the max (kappa, |nu|) -- the overflow-safe max-subtracted trapezoid,
## matching _vmft_log_c_vec (and the cached scalar table's log_normalizer at the
## same grid size). kappa <= tol -> uniform -log 2pi.
.vmft_log_c_vec <- function(kappa, nu) {
  L <- max(length(kappa), length(nu))
  kappa <- rep_len(as.numeric(kappa), L)
  nu <- rep_len(as.numeric(nu), L)
  out <- rep(-log(2.0 * pi), L)
  live <- kappa > .VMFT_KAPPA_TOL
  if (!any(live)) return(out)
  kl <- pmin(pmax(kappa[live], 0.0), .VMFT_KAPPA_UPPER)
  nl <- nu[live]
  key <- paste(kl, nl, sep = "\r")
  uk <- !duplicated(key)
  uidx <- which(uk)
  k <- kl[uidx]
  nv <- nl[uidx]
  grid <- .vmft_grid_size(max(k), max(abs(nv)))
  phi <- seq(-pi, pi, length.out = grid + 1L)
  sin_phi <- sin(phi)
  w <- rep(1.0, grid + 1L)
  w[1L] <- 0.5
  w[grid + 1L] <- 0.5
  B <- sweep(outer(nv, sin_phi), 2L, phi, "+")    # U x (grid+1)
  log_ker <- k * cos(B)                            # k recycles per row (col-major)
  mx <- apply(log_ker, 1L, max)
  log_Z <- log(2.0 * pi / grid) + mx +
    log(drop(exp(log_ker - mx) %*% w))
  vals <- -log_Z
  out[live] <- vals[match(key, key[uk])]
  out
}

## First and second (kappa,nu)-derivatives of log Z(kappa,nu) as kernel-weighted
## moments under the density, vectorized over all unique (kappa,nu) pairs in one
## shared-grid pass -- the normalizer block of the flat-topped vM score/Hessian.
## With h_kappa = cos B and h_nu = -kappa sin B sin phi:
##   d logZ/d kappa = E[h_kappa],    d logZ/d nu = E[h_nu],
##   d2 logZ/da db  = E[h_ab] + Cov(h_a, h_b),
## using the kernel second derivatives h_kk = 0, h_kn = -sin B sin phi,
## h_nn = -kappa cos B sin^2 phi. Each E[.] is a normalizer-free ratio (the
## e^{kappa cosB - max} factor cancels), accurate at any concentration. Returns a
## list of (dk, dnu, dkk, dknu, dnunu) aligned to the input length. No uniform
## shortcut -- the score stays live in the near-uniform corner.
.vmft_logZ_moments_vec <- function(kappa, nu) {
  L <- max(length(kappa), length(nu))
  kappa <- rep_len(as.numeric(kappa), L)
  nu <- rep_len(as.numeric(nu), L)
  kf <- pmin(pmax(kappa, 0.0), .VMFT_KAPPA_UPPER)
  nf <- nu
  key <- paste(kf, nf, sep = "\r")
  uk <- !duplicated(key)
  uidx <- which(uk)
  k <- kf[uidx]
  nv <- nf[uidx]
  grid <- .vmft_grid_size(max(k), max(abs(nv)))
  phi <- seq(-pi, pi, length.out = grid + 1L)
  sin_phi <- sin(phi)
  sin2 <- sin_phi * sin_phi
  w <- rep(1.0, grid + 1L)
  w[1L] <- 0.5
  w[grid + 1L] <- 0.5
  B <- sweep(outer(nv, sin_phi), 2L, phi, "+")    # U x (grid+1)
  cosB <- cos(B)
  sinB <- sin(B)
  log_ker <- k * cosB
  rowmax <- apply(log_ker, 1L, max)
  e <- exp(log_ker - rowmax)
  e <- sweep(e, 2L, w, "*")
  Z <- pmax(rowSums(e), .Machine$double.xmin)
  m <- function(V) rowSums(e * V) / Z
  sBsp <- sweep(sinB, 2L, sin_phi, "*")           # sin B sin phi
  E_cosB <- m(cosB)
  E_sBsp <- m(sBsp)
  dk <- E_cosB                                    # d logZ / d kappa
  dnu <- -k * E_sBsp                              # d logZ / d nu = E[h_nu]
  dkk <- m(cosB * cosB) - E_cosB * E_cosB         # Var(cos B); h_kk = 0
  cov_k_nu <- -k * (m(cosB * sBsp) - E_cosB * E_sBsp)   # Cov(h_kappa, h_nu)
  dknu <- -E_sBsp + cov_k_nu                      # E[h_kn] = -E[sin B sin phi]
  var_nu <- k * k * (m(sBsp * sBsp) - E_sBsp * E_sBsp)  # Var(h_nu)
  dnunu <- -k * m(sweep(cosB, 2L, sin2, "*")) + var_nu  # E[h_nn] + Var(h_nu)
  idx <- match(key, key[uk])
  list(dk = dk[idx], dnu = dnu[idx], dkk = dkk[idx],
       dknu = dknu[idx], dnunu = dnunu[idx])
}

## q-th cosine moment alpha_q = E[cos(q phi)] of the centered flat-topped vM law
## (the sine moments vanish, the density being symmetric in phi). Serves the
## Pearson residual scale Var{sin phi} = (1 - alpha_2)/2. Grid moment on the same
## adaptive trapezoid as the normalizer. Mirrors the .jp_cos_moment role.
.vmft_cos_moment <- function(kappa, nu, q) {
  if (q == 0) return(1.0)
  if (kappa < .VMFT_KAPPA_TOL) return(0.0)
  grid <- .vmft_grid_size(kappa, nu)
  phi <- seq(-pi, pi, length.out = grid + 1L)
  w <- rep(1.0, grid + 1L)
  w[1L] <- 0.5
  w[grid + 1L] <- 0.5
  B <- phi + nu * sin(phi)
  lk <- kappa * cos(B)
  e <- w * exp(lk - max(lk))
  denom <- max(sum(e), .Machine$double.xmin)
  sum(cos(q * phi) * e) / denom
}

## The flat-topped vM log-density value and its (mu,kappa,nu) score/Hessian at
## the data, the single entry point the family ll() and the standalone validator
## share (the .kj_terms / .jp_score_terms convention). With phi = y - mu,
## B = phi + nu sin phi, B_phi = 1 + nu cos phi:
##   l0   = kappa cos B + log c(kappa, nu)           (uniform -log 2pi if k<=tol)
##   l_mu = kappa sin B B_phi                         (Z is mu-free: no norm term)
##   l_kappa = cos B - E[cos B]
##   l_nu = -kappa sin B sin phi - E[h_nu]
## l2 columns are the unique unordered pairs in combinations_with_replacement
## order (mu,kappa,nu): (mm, mk, mn, kk, kn, nn); the (mu,.) blocks are pure
## kernel (d^2 logZ/d mu d. = 0), the (k,k)/(k,n)/(n,n) blocks subtract the
## normalizer second derivatives. Returns l0, and (deriv) l1 [n x 3], l2 [n x 6].
.vmft_terms <- function(y, mu, kappa, nu, deriv = TRUE) {
  n <- length(y)
  mu <- rep_len(as.numeric(mu), n)
  kappa <- rep_len(as.numeric(kappa), n)
  nu <- rep_len(as.numeric(nu), n)
  phi <- y - mu
  s <- sin(phi)
  cc <- cos(phi)
  B <- phi + nu * s
  sinB <- sin(B)
  cosB <- cos(B)
  Bphi <- 1.0 + nu * cc
  l0 <- kappa * cosB + .vmft_log_c_vec(kappa, nu)
  unif <- kappa <= .VMFT_KAPPA_TOL
  if (any(unif)) l0[unif] <- -log(2.0 * pi)
  out <- list(l0 = l0)
  if (deriv) {
    M <- .vmft_logZ_moments_vec(kappa, nu)
    l1 <- cbind(kappa * sinB * Bphi,
                cosB - M$dk,
                -kappa * sinB * s - M$dnu)
    l2 <- cbind(-kappa * cosB * Bphi * Bphi + kappa * nu * sinB * s,
                sinB * Bphi,
                kappa * (cosB * s * Bphi + sinB * cc),
                -M$dkk,
                -sinB * s - M$dknu,
                -kappa * cosB * s * s - M$dnunu)
    colnames(l1) <- c("mu", "kappa", "nu")
    colnames(l2) <- c("mm", "mk", "mn", "kk", "kn", "nn")
    out$l1 <- l1
    out$l2 <- l2
  }
  out
}
