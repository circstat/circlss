#' Asymmetric Jones-Pewsey location-concentration-shape-asymmetry family
#'
#' A general family implementing distributional regression for a circular
#' (angular) response \eqn{y} in radians under the asymmetric Jones-Pewsey law
#' \deqn{f(y) = c(\kappa, \psi, \nu)\,\bigl(\cosh(\kappa\psi) +
#'   \sinh(\kappa\psi)\cos g\bigr)^{1/\psi}, \qquad
#'   g = \phi + \nu\cos\phi,\ \ \phi = y - \xi,}
#' the Jones-Pewsey family (\code{\link{jplss}}) with its angle warped
#' \emph{forward} by \eqn{g(\phi) = \phi + \nu\cos\phi} (Abe, Pewsey & Shimizu
#' 2013). It adds an asymmetry axis to the symmetric
#' Jones-Pewsey umbrella: \eqn{\nu > 0} and \eqn{\nu < 0} tilt the density to
#' opposite sides. (For \eqn{\psi \to 0} the Jones-Pewsey kernel is the von Mises
#' \eqn{\exp(\kappa\cos g)}.)
#'
#' It has \strong{four linear predictors}: the location \eqn{\xi}, the
#' concentration \eqn{\kappa}, the shape \eqn{\psi} and the asymmetry \eqn{\nu}
#' each get their own, any of which may contain smooth terms. Used with
#' \code{mgcv::gam} and a list of \emph{four} formulas; the first names the
#' response and models \eqn{\xi}, then \eqn{\log\kappa}, then \eqn{\psi}, then
#' \eqn{\nu}. The shape and asymmetry are most often held global, i.e. fitted
#' intercept-only with \code{~ 1}.
#'
#' @param link Four-element list of link names, for the location, concentration,
#'   shape and asymmetry. Currently only the defaults are available:
#'   \code{"tanhalf"} for the location \eqn{\xi}, \code{"log"} for the
#'   concentration \eqn{\kappa > 0}, \code{"identity"} for the shape
#'   \eqn{\psi \in \mathbb{R}}, and \code{"tanh"} for the asymmetry
#'   \eqn{\nu \in (-1, 1)}.
#' @return An object of class \code{c("general.family", "extended.family", "family")}
#'   for use with \code{\link[mgcv]{gam}} (or its front end \code{\link{circ_gam}}).
#' @details
#' The shape \eqn{\psi} and asymmetry \eqn{\nu} index a family of symmetric and
#' skewed circular laws:
#' \itemize{
#'   \item \eqn{\nu = 0}: the symmetric Jones-Pewsey (\code{\link{jplss}}), whose
#'     \eqn{\psi} in turn nests the von Mises (\eqn{\psi \to 0}), the cardioid
#'     (\eqn{\psi = 1}) and the wrapped Cauchy (\eqn{\psi \to -\infty});
#'   \item \eqn{\nu \ne 0}: an asymmetric (skewed) law;
#'   \item \eqn{\kappa \to 0}: the circular uniform.
#' }
#'
#' \strong{Mode anchor, not mean.} As with \code{\link{ssjplss}} and
#' \code{\link{ibslss}}, once \eqn{\nu \ne 0} the location \eqn{\xi} is the
#' \emph{mode anchor} of the asymmetric density, not its mean direction;
#' fitted-direction summaries inherit that reading.
#'
#' \strong{Forward warp -- no implicit differentiation.} The warp
#' \eqn{g(\phi) = \phi + \nu\cos\phi} is a monotone reparameterization of the
#' angle (\eqn{g'(\phi) = 1 - \nu\sin\phi > 0} for \eqn{|\nu| < 1}), so -- unlike
#' \code{\link{ibslss}}'s inverse warps -- the score needs \emph{no} implicit
#' differentiation: it reuses the Jones-Pewsey kernel terms evaluated at the
#' warped angle \eqn{g} and chains them through \eqn{g} (with
#' \eqn{g_\phi = 1 - \nu\sin\phi} and \eqn{g_\nu = \cos\phi}).
#'
#' \strong{Normalizer and cross terms.} Unlike \code{\link{ssjplss}}'s sine-skew
#' factor -- which integrates to 1 against the symmetric kernel and so leaves the
#' Jones-Pewsey normalizer \emph{unchanged} -- the forward warp \emph{moves} the
#' normalizer: \eqn{c(\kappa, \psi, \nu)} is the \eqn{u = g(\phi)} substitution
#' integrating the kernel against the warp Jacobian \eqn{1/g'} over an
#' asymmetry-aware Gauss-Legendre ladder. Consequently the \eqn{(\kappa,\nu)} and
#' \eqn{(\psi,\nu)} cross-derivatives do \strong{not} vanish (they do for
#' \code{ssjplss}); the \eqn{\nu}-score and the full cross-parameter Hessian carry
#' the \eqn{(\kappa, \psi, \nu)} log-normalizer moments returned by that
#' quadrature.
#'
#' \strong{Optimizer.} Only first- and second-order log-likelihood derivatives
#' are provided, so \code{available.derivs = 0} and the family is fitted by the
#' extended Fellner-Schall optimizer rather than full Newton REML.
#' \code{\link[mgcv]{gam}} selects \code{optimizer = "efs"} automatically; passing
#' it explicitly is recommended so the fitted object is labelled correctly (which
#' \code{\link[mgcv]{gam.check}} relies on). The shape \eqn{\psi} and asymmetry
#' \eqn{\nu} are weakly identified when the response is diffuse; holding them
#' global (\code{~ 1}) and keeping the concentration well identified is the robust
#' default.
#'
#' The location uses the Fisher-Lee tan-half link (\eqn{\xi \in (-\pi, \pi)},
#' antipode unrepresentable, winding number zero -- see \code{\link{pnlss}} when
#' the location must wind). The asymmetry rides the tanh link (shared with
#' \code{\link{ssjplss}}, \code{\link{vmftlss}} and \code{\link{ibslss}}), bounded
#' to \eqn{(-1, 1)}.
#' @examples
#' library(mgcv)
#' set.seed(1)
#' n <- 300
#' x <- runif(n)
#' xi <- 2 * atan(sin(2 * pi * x))
#' kappa <- exp(1.0 + 0.4 * cos(2 * pi * x))
#' psi <- 0.5
#' nu <- 0.4       # asymmetry
#' # asymmetric Jones-Pewsey draws by inverse transform on the warped kernel:
#' # the forward warp g(phi) is closed-form, so gridding phi gives a
#' # (phi, density) grid to invert -- no root-finding needed for simulation.
#' phi <- seq(-pi, pi, length.out = 4001)
#' g <- phi + nu * cos(phi)
#' dev <- vapply(seq_len(n), function(i) {
#'   d <- (cosh(kappa[i] * psi) + sinh(kappa[i] * psi) * cos(g))^(1 / psi)
#'   d[d < 0] <- 0
#'   cdf <- cumsum((d[-1] + d[-length(d)]) / 2 * diff(phi))
#'   cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
#'   approx(cdf, phi, runif(1), rule = 2)$y
#' }, numeric(1))
#' y <- atan2(sin(xi + dev), cos(xi + dev))
#' # smooth location and concentration; global (intercept-only) shape and asymmetry
#' b <- gam(list(y ~ s(x), ~ s(x), ~ 1, ~ 1), family = ajplss(), optimizer = "efs")
#' summary(b)
#' @references
#' Jones, M. C. and Pewsey, A. (2005) A family of symmetric distributions on the
#' circle. \emph{Journal of the American Statistical Association} 100, 1422-1428.
#'
#' Abe, T., Pewsey, A. and Shimizu, K. (2013) Extending circular distributions
#' through transformation of argument. \emph{Annals of the Institute of
#' Statistical Mathematics} 65, 833-858.
#'
#' Batschelet, E. (1981) \emph{Circular Statistics in Biology}. Academic Press.
#'
#' Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method for
#' smoothing parameter optimization with application to Tweedie location,
#' scale and shape models. \emph{Biometrics} 73, 1071-1081.
#' @seealso
#' \code{\link{jplss}}, \code{\link{ssjplss}}, \code{\link{vmlss}},
#' \code{\link{ibslss}}, \code{\link{vmftlss}}, \code{\link{pnlss}},
#' \code{\link[mgcv]{gam}}
#' @export
ajplss <- function(link = list("tanhalf", "log", "identity", "tanh")) {
  if (length(link) != 4)
    stop("ajplss requires 4 links specified as character strings")
  if (!(link[[1]] %in% "tanhalf"))
    stop(link[[1]], " link not available for the location parameter of ajplss")
  if (!(link[[2]] %in% "log"))
    stop(link[[2]], " link not available for the concentration parameter of ajplss")
  if (!(link[[3]] %in% "identity"))
    stop(link[[3]], " link not available for the shape parameter of ajplss")
  if (!(link[[4]] %in% "tanh"))
    stop(link[[4]], " link not available for the asymmetry parameter of ajplss")

  stats <- list()
  stats[[1]] <- tanhalf.link()
  stats[[2]] <- stats::make.link("log")
  fam <- structure(list(link = "log", canonical = "none",
                        linkfun = stats[[2]]$linkfun,
                        mu.eta = stats[[2]]$mu.eta),
                   class = "family")
  fam <- mgcv::fix.family.link(fam)
  stats[[2]]$d2link <- fam$d2link
  stats[[2]]$d3link <- fam$d3link
  stats[[2]]$d4link <- fam$d4link
  stats[[3]] <- stats::make.link("identity")
  stats[[3]]$d2link <- function(mu) numeric(length(mu))
  stats[[3]]$d3link <- function(mu) numeric(length(mu))
  stats[[3]]$d4link <- function(mu) numeric(length(mu))
  stats[[4]] <- tanh_link()

  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    type <- match.arg(type)
    xi <- object$fitted[, 1]
    kappa <- object$fitted[, 2]
    psi <- object$fitted[, 3]
    nu <- object$fitted[, 4]
    d <- object$y - xi
    if (type == "response") return(atan2(sin(d), cos(d)))
    sind <- sin(d)
    if (type == "pearson") {
      ## base (symmetric) JP scale; the warp tilts the mode, so this is a
      ## diagnostic standardization, not an exact residual variance
      a2 <- vapply(seq_along(kappa),
                   function(i) .jp_cos_moment(kappa[i], psi[i], 2), numeric(1))
      v <- (1 - a2) / 2
      v[v <= 0] <- .Machine$double.eps
      return(sind / sqrt(v))
    }
    ## deviance: 2*(l_sat - l_hat); the normalizer cancels, leaving the
    ## warp-tilted peak (grid max) minus the attained centered log-kernel
    g <- d + nu * cos(d)
    lhat <- .jp_score_terms(g, kappa, psi, second = FALSE)$h
    lsat <- .ajp_lsat(kappa, psi, nu)
    sign(sind) * sqrt(pmax(2 * (lsat - lhat), 0))
  }

  postproc <- expression({
    ## Deviance and null deviance against the common saturated reference (the
    ## warp-tilted per-obs peak). The deviance is normalizer-free; the null
    ## (constant von Mises MLE, psi -> 0, nu -> 0) needs the fitted normalizer,
    ## one quadrature pass via the internal helpers. local() so mgcv's frames
    ## stay intact.
    .ajplss.dev <- local({
      st <- getFromNamespace(".jp_score_terms", "circlss")
      lcv <- getFromNamespace(".jp_log_c_asym_vec", "circlss")
      lsatf <- getFromNamespace(".ajp_lsat", "circlss")
      xi <- object$fitted[, 1]
      kappa <- object$fitted[, 2]
      psi <- object$fitted[, 3]
      nu <- object$fitted[, 4]
      y <- object$y
      d <- y - xi
      g <- d + nu * cos(d)
      h <- st(g, kappa, psi, second = FALSE)$h
      lsat0 <- lsatf(kappa, psi, nu)                   # normalizer-free peak
      lhat0 <- h                                       # normalizer-free attained
      dev <- sum(2 * pmax(lsat0 - lhat0, 0))
      logc <- lcv(kappa, psi, nu)
      lsat <- sum(lsat0 + logc)
      sy <- mean(sin(y)); cy <- mean(cos(y))
      mu0 <- atan2(sy, cy)
      Rbar <- min(sqrt(sy * sy + cy * cy), 1 - 1e-08)
      kappa0 <- if (Rbar < 0.53) 2 * Rbar + Rbar^3 + 5 * Rbar^5 / 6 else
        if (Rbar < 0.85) -0.4 + 1.39 * Rbar + 0.43 / (1 - Rbar) else
          1 / (Rbar^3 - 4 * Rbar^2 + 3 * Rbar)
      kappa0 <- min(max(kappa0, 0.01), 500)
      lnull <- sum(kappa0 * (cos(y - mu0) - 1) -
                     log(2 * pi * besselI(kappa0, 0, expon.scaled = TRUE)))
      c(dev, 2 * (lsat - lnull))
    })
    object$deviance <- .ajplss.dev[1]
    object$null.deviance <- .ajplss.dev[2]
    rm(.ajplss.dev)
  })

  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                 D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
    if (is.null(offset)) offset <- vector("list", 5) else offset[[5]] <- 0
    if (is.null(wt)) wt <- 1            # prior weights; sandwich/direct callers may pass NULL
    jj <- attr(X, "lpi")
    if (is.null(eta)) {
      eta0 <- X[, jj[[1]], drop = FALSE] %*% coef[jj[[1]]]
      if (!is.null(offset[[1]])) eta0 <- eta0 + offset[[1]]
      eta1 <- X[, jj[[2]], drop = FALSE] %*% coef[jj[[2]]]
      if (!is.null(offset[[2]])) eta1 <- eta1 + offset[[2]]
      eta2 <- X[, jj[[3]], drop = FALSE] %*% coef[jj[[3]]]
      if (!is.null(offset[[3]])) eta2 <- eta2 + offset[[3]]
      eta3 <- X[, jj[[4]], drop = FALSE] %*% coef[jj[[4]]]
      if (!is.null(offset[[4]])) eta3 <- eta3 + offset[[4]]
    } else {
      eta0 <- eta[, 1]; eta1 <- eta[, 2]; eta2 <- eta[, 3]; eta3 <- eta[, 4]
    }
    eta0 <- drop(eta0); eta1 <- drop(eta1); eta2 <- drop(eta2); eta3 <- drop(eta3)
    xi <- family$linfo[[1]]$linkinv(eta0)
    kappa <- family$linfo[[2]]$linkinv(eta1)
    psi <- family$linfo[[3]]$linkinv(eta2)
    nu <- family$linfo[[4]]$linkinv(eta3)
    phi <- y - xi
    cphi <- cos(phi); sphi <- sin(phi)
    g <- phi + nu * cphi                    # forward asymmetry warp
    gphi <- 1 - nu * sphi                   # g'(phi); g_xi = -gphi

    ts <- .jp_score_terms(g, kappa, psi, second = deriv > 0)
    ## Normalizer + (when deriv) the (kappa, psi, nu) log-Z moments. deriv = 1
    ## takes ONE fused quadrature sweep -- .jp_quad_asym_vec returns logc in
    ## column 1 and the nine moments in columns 2:10 -- instead of integrating
    ## the normalizer twice (a separate .jp_log_c_asym_vec plus a moment sweep);
    ## deriv = 0 takes the cheaper value-only sweep.
    if (deriv) {
      Q <- .jp_quad_asym_vec(kappa, psi, nu)
      logc <- Q[, 1L]
      M <- Q[, -1L, drop = FALSE]   # n x 9: dk dp dnu dkk dkp dpp dknu dpnu dnunu
    } else {
      logc <- .jp_log_c_asym_vec(kappa, psi, nu)
    }
    ## l0 VALUE: log-kernel at the warped angle + u-substituted normalizer.
    ## kappa < KAPPA_TOL -> uniform (zero the kernel; logc is already -log 2pi).
    h0 <- ts$h
    unif <- kappa < .JP_KAPPA_TOL
    if (any(unif)) h0[unif] <- 0
    l0 <- h0 + logc
    l <- sum(wt * l0)

    if (deriv) {
      ## l1: d l / d(xi, kappa, psi, nu)
      l1 <- cbind(-ts$hphi * gphi,
                  ts$hk - M[, 1],
                  ts$hp - M[, 2],
                  ts$hphi * cphi - M[, 3])
      ## l2 columns in combinations_with_replacement order:
      ## (xx, xk, xp, xn, kk, kp, kn, pp, pn, nn)
      l2 <- cbind(ts$hphiphi * gphi * gphi - nu * cphi * ts$hphi,
                  -ts$hphik * gphi,
                  -ts$hphip * gphi,
                  -ts$hphiphi * cphi * gphi + ts$hphi * sphi,
                  ts$hkk - M[, 4],
                  ts$hkp - M[, 5],
                  ts$hphik * cphi - M[, 7],
                  ts$hpp - M[, 6],
                  ts$hphip * cphi - M[, 8],
                  ts$hphiphi * cphi * cphi - M[, 9])
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta0),
                   family$linfo[[2]]$mu.eta(eta1),
                   family$linfo[[3]]$mu.eta(eta2),
                   family$linfo[[4]]$mu.eta(eta3))
      g2 <- cbind(family$linfo[[1]]$d2link(xi),
                  family$linfo[[2]]$d2link(kappa),
                  family$linfo[[3]]$d2link(psi),
                  family$linfo[[4]]$d2link(nu))
      l3 <- l4 <- g3 <- g4 <- 0
      i2 <- family$tri$i2
      i3 <- family$tri$i3
      i4 <- family$tri$i4
      ## prior weights: a weighted log-likelihood scales every per-observation
      ## derivative row by wt. gamlss.etamu is linear per row in (l1,l2,l3,l4),
      ## so scaling these inputs == scaling the eta-derivatives gamlss.gH builds.
      ## scalar-0 placeholders for absent orders are left untouched.
      l1 <- wt * l1
      l2 <- wt * l2
      if (is.matrix(l3)) l3 <- wt * l3
      if (is.matrix(l4)) l4 <- wt * l4
      de <- mgcv::gamlss.etamu(l1, l2, l3, l4, ig1, g2, g3, g4,
                               i2, i3, i4, deriv - 1)
      ret <- mgcv::gamlss.gH(X, jj, de$l1, de$l2, i2, l3 = de$l3, i3 = i3,
                             l4 = de$l4, i4 = i4, d1b = d1b, d2b = d2b,
                             deriv = deriv - 1, fh = fh, D = D,
                             sandwich = sandwich)
      if (ncv) {
        ret$l1 <- de$l1
        ret$l2 <- de$l2
        ret$l3 <- de$l3
      }
    } else ret <- list()
    ret$l <- l
    ret$l0 <- l0
    ret
  }

  sandwich <- function(y, X, coef, wt, family, offset = NULL) {
    ll(y, X, coef, wt, family, offset = NULL, deriv = 1, sandwich = TRUE)$lbb
  }

  initialize <- expression({
    ## Location (mode anchor) starts from the data-following projected pilot
    ## (.tanhalf_pilot_coef). Concentration starts from A1-inverse kappa0, shape
    ## from psi0 = 0 (von Mises) and asymmetry from nu0 = 0 (symmetric), as
    ## constant targets through their links. local() so mgcv's frames are safe.
    n <- rep(1, nobs)
    use.unscaled <- if (!is.null(attr(E, "use.unscaled"))) TRUE else FALSE
    if (is.null(start)) {
      start <- local({
        jj <- attr(x, "lpi")
        pls <- function(Xd, Ed, yt) {
          b <- qr.coef(qr(rbind(Xd, Ed)), c(yt, rep(0, nrow(Ed))))
          b[!is.finite(b)] <- 0
          b
        }
        st <- rep(0, ncol(x))
        sy <- mean(sin(y)); cy <- mean(cos(y))
        Rbar <- min(sqrt(sy * sy + cy * cy), 1 - 1e-08)
        kappa0 <- if (Rbar < 0.53) 2 * Rbar + Rbar^3 + 5 * Rbar^5 / 6 else
          if (Rbar < 0.85) -0.4 + 1.39 * Rbar + 0.43 / (1 - Rbar) else
            1 / (Rbar^3 - 4 * Rbar^2 + 3 * Rbar)
        kappa0 <- min(max(kappa0, 0.01), 500)
        ## location: data-following projected pilot
        off1 <- if (!is.null(offset) && length(offset) >= 1 &&
                    !is.null(offset[[1]])) offset[[1]] else 0
        st[jj[[1]]] <- getFromNamespace(".tanhalf_pilot_coef", "circlss")(
          x[, jj[[1]], drop = FALSE], E[, jj[[1]], drop = FALSE], y, off1,
          family$linfo[[1]]$linkfun)
        ## concentration (log), shape (identity), asymmetry (tanh): constant targets
        targ <- c(NA, family$linfo[[2]]$linkfun(kappa0),
                  family$linfo[[3]]$linkfun(0), family$linfo[[4]]$linkfun(0))
        for (k in 2:4) {
          off <- if (!is.null(offset) && length(offset) >= k &&
                     !is.null(offset[[k]])) offset[[k]] else 0
          st[jj[[k]]] <- pls(x[, jj[[k]], drop = FALSE],
                             E[, jj[[k]], drop = FALSE],
                             rep(targ[k], nobs) - off)
        }
        st
      })
    }
  })

  rd <- function(mu, wt, scale) {
    ## Simulation by inverse transform on a fine centered grid of the WARPED
    ## kernel, one grid per unique (kappa, psi, nu).
    n <- nrow(mu)
    xiv <- mu[, 1]; kappa <- mu[, 2]; psi <- mu[, 3]; nu <- mu[, 4]
    grid <- seq(-pi, pi, length.out = 2049L)
    th <- numeric(n)
    key <- paste(kappa, psi, nu, sep = "\r")
    for (k in unique(key)) {
      idx <- which(key == k)
      g <- grid + nu[idx[1]] * cos(grid)
      h <- .jp_score_terms(g, kappa[idx[1]], psi[idx[1]], second = FALSE)$h
      dens <- exp(h - max(h))
      dens[dens < 0] <- 0
      cdf <- cumsum((dens[-1] + dens[-length(dens)]) * diff(grid) / 2)
      cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
      u <- stats::runif(length(idx))
      th[idx] <- xiv[idx] +
        stats::approx(cdf, grid, xout = u, rule = 2, ties = "ordered")$y
    }
    atan2(sin(th), cos(th))
  }

  structure(list(family = "ajplss", ll = ll, link = paste(link), nlp = 4,
                 param_names = c("xi", "kappa", "psi", "nu"),
                 param_circular = c(TRUE, FALSE, FALSE, FALSE),
                 sandwich = sandwich, tri = mgcv::trind.generator(4),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 0, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
