## Bessel-ratio helpers for the von Mises family: A1, A1inv, A1prime, A1prime2
## and A1prime3.
## Exponentially scaled Bessel functions keep
## every expression finite at large kappa, where I0/I1 themselves overflow.

A1 <- function(kappa) {
  besselI(kappa, 1, expon.scaled = TRUE) / besselI(kappa, 0, expon.scaled = TRUE)
}

## Inverse of A1: solve A1(kappa) = R for the von Mises concentration kappa >= 0.
## This is the MLE of kappa given a mean resultant length R, so accuracy is the
## reconstruction error |A1(A1inv(R)) - R|. The classical closed form (Fisher
## 1993, as in circular::A1inv) is a three-branch rational approximation whose
## error reaches ~2.5e-3 near the branch seams (R = 0.53, 0.85); that R-space
## error is amplified in kappa where A1 is flat (high concentration, A1' -> 0).
## Here that closed form only SEEDS a short Newton refinement against A1()/
## A1prime(), driving the reconstruction error to machine precision on (0, 1) --
## i.e. the returned kappa exactly inverts A1. (So circ_lm()'s kappa can differ
## from circular::lm.circular()'s at the approximation's error level, largest at
## high concentration; every non-kappa quantity, which does not pass through this
## inverse, matches to machine precision.) R is clamped to [0, 1]; R = 0 maps to
## 0, R = 1 to +Inf.
A1inv <- function(R) {
  R <- pmin(pmax(as.numeric(R), 0), 1)
  k <- ifelse(R < 0.53, 2 * R + R^3 + 5 * R^5 / 6,
       ifelse(R < 0.85, -0.4 + 1.39 * R + 0.43 / (1 - R),
              1 / pmax(R^3 - 4 * R^2 + 3 * R, .Machine$double.eps)))
  ok <- is.finite(k) & R > 0 & R < 1
  for (.i in seq_len(8L)) {
    kk <- k[ok]
    if (!length(kk)) break
    kk <- kk - (A1(kk) - R[ok]) / A1prime(kk)
    kk[!is.finite(kk) | kk < 0] <- 0
    k[ok] <- pmin(kk, 1e8)
  }
  k[R <= 0] <- 0
  k[R >= 1] <- Inf
  k
}

A1prime <- function(kappa, a1 = A1(kappa)) {
  ## A1'(k) = 1 - A1(k)/k - A1(k)^2; A1(k)/k -> 1/2 as k -> 0 is the
  ## removable singularity, filled in to avoid 0/0.
  r <- a1 / kappa
  r[kappa == 0] <- 0.5
  1 - r - a1 * a1
}

A1prime2 <- function(kappa, a1 = A1(kappa), d1 = A1prime(kappa, a1)) {
  ## A1''(k) = -A1'/k + A1/k^2 - 2 A1 A1'. The two 1/k terms cancel
  ## catastrophically as k -> 0 (true value ~ -3k/8), so below k < 0.01 a
  ## Maclaurin series is used; both branches agree to ~1e-12 at the switch.
  rec <- -d1 / kappa + a1 / kappa^2 - 2 * a1 * d1
  k2 <- kappa * kappa
  series <- kappa * (-3 / 8 + k2 * (5 / 24 - k2 * (77 / 1024)))
  out <- ifelse(kappa < 0.01, series, rec)
  out
}

A1prime3 <- function(kappa, a1 = A1(kappa), d1 = A1prime(kappa, a1),
                     d2 = A1prime2(kappa, a1, d1)) {
  ## A1'''(k) = -A1''/k + 2 A1'/k^2 - 2 A1/k^3 - 2 A1'^2 - 2 A1 A1'';
  ## same removable cancellation as A1prime2, series below k < 0.01.
  rec <- -d2 / kappa + 2 * d1 / kappa^2 - 2 * a1 / kappa^3 -
    2 * d1 * d1 - 2 * a1 * d2
  k2 <- kappa * kappa
  series <- -3 / 8 + k2 * (5 / 8 - k2 * (385 / 1024))
  out <- ifelse(kappa < 0.01, series, rec)
  out
}
