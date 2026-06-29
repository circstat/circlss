#' Inverse Batschelet location-concentration-skewness-peakedness family
#'
#' A general family implementing distributional regression for a circular
#' (angular) response \eqn{y} in radians under the inverse Batschelet law
#' \deqn{f(y) = c(\kappa, \lambda)\,\exp\!\bigl(\kappa\cos A\bigr), \qquad
#'   \phi = y - \xi,\ \kappa > 0,\ \nu, \lambda \in (-1, 1),}
#' where the angle is warped forward into the von Mises cosine kernel by two
#' \emph{inverse} maps: a skewness warp \eqn{\phi^\star = t_\nu^{-1}(\phi)}
#' solving \eqn{y - \nu(1 + \cos y) = \phi}, then a peakedness warp
#' \eqn{u^\star = s_\lambda^{-1}(\phi^\star)} solving
#' \eqn{u - \tfrac{1+\lambda}{2}\sin u = \phi^\star}, with the kernel argument
#' \eqn{A = u^\star - \tfrac{1-\lambda}{2}\sin u^\star}. The constant
#' \eqn{c(\kappa, \lambda)} depends only on \eqn{\kappa} and \eqn{\lambda}.
#'
#' The skewness \eqn{\nu} tilts the density (mode toward smaller or larger
#' angles) and the peakedness \eqn{\lambda} flattens (\eqn{\lambda < 0}) or
#' sharpens (\eqn{\lambda > 0}) the peak; together they give a flexible skew-and-
#' peaked alternative to the von Mises. Because the warps are \emph{inverse}
#' maps the transform has no closed form (a monotone solver inverts it) and,
#' unlike the flat-topped von Mises (\code{\link{vmftlss}}, a forward warp), the
#' score requires implicit differentiation.
#'
#' It has \strong{four linear predictors}: the location \eqn{\xi}, the
#' concentration \eqn{\kappa}, the skewness \eqn{\nu}, and the peakedness
#' \eqn{\lambda} each get their own, any of which may contain smooth terms. Used
#' with \code{mgcv::gam} and a list of \emph{four} formulas; the first names the
#' response and models \eqn{\xi}, the second models \eqn{\log\kappa}, the third
#' \eqn{\nu}, the fourth \eqn{\lambda}. The skewness and peakedness are most
#' often held global, i.e. fitted intercept-only with \code{~ 1}.
#'
#' @param link Four-element list of link names, for the location, the
#'   concentration, the skewness, and the peakedness. Currently only the defaults
#'   are available: \code{"tanhalf"} for the location \eqn{\xi}, \code{"log"} for
#'   the concentration \eqn{\kappa > 0}, \code{"tanh"} for the skewness
#'   \eqn{\nu \in (-1, 1)}, and \code{"tanh"} for the peakedness
#'   \eqn{\lambda \in (-1, 1)}.
#' @return An object of class \code{c("general.family", "extended.family", "family")}
#'   for use with \code{\link[mgcv]{gam}} (or its front end \code{\link{circ_gam}}).
#' @details
#' The two shape parameters index a family of skew and/or peaked circular laws
#' and recover the von Mises exactly:
#' \itemize{
#'   \item \eqn{\nu = \lambda = 0}: the von Mises (\code{\link{vmlss}}),
#'     \eqn{\propto \exp(\kappa\cos(y-\xi))};
#'   \item \eqn{\nu \ne 0}: a skewed (asymmetric) law;
#'   \item \eqn{\lambda > 0}: a sharper-than-von-Mises peak;
#'     \eqn{\lambda < 0}: a flatter, flat-topped peak;
#'   \item \eqn{\kappa \to 0}: the circular uniform.
#' }
#'
#' \strong{Mode anchor.} As with \code{\link{ssjplss}}, once \eqn{\nu \ne 0} the
#' location \eqn{\xi} is the \emph{mode anchor}, not the mean direction.
#'
#' \strong{Warps and normalizer.} The two inverse warps are inverted by a
#' vectorized monotone Newton iteration with a bisection mop-up near the
#' boundary (\eqn{\nu, \lambda \to \pm 1}); the per-observation score reuses the
#' solved roots \eqn{\phi^\star, u^\star} and the warp slopes. The normalizer
#' \eqn{c(\kappa, \lambda) = (1-\lambda) / [(1+\lambda)\,2\pi I_0(\kappa) -
#' 2\lambda\int e^{\kappa\cos B}]} has no elementary form: the integral is an
#' overflow-safe equispaced trapezoid (the warped kernel is smooth and periodic,
#' so the trapezoid is spectrally accurate) and the \eqn{\kappa}- and
#' \eqn{\lambda}-derivatives of \eqn{\log c} are central finite differences of
#' it.
#'
#' \strong{Optimizer.} Only first- and second-order log-likelihood derivatives
#' are provided, so \code{available.derivs = 0} and the family is fitted by the
#' extended Fellner-Schall optimizer rather than full Newton REML.
#' \code{\link[mgcv]{gam}} selects \code{optimizer = "efs"} automatically; passing
#' it explicitly is recommended so the fitted object is labelled correctly (which
#' \code{\link[mgcv]{gam.check}} relies on). The skewness and peakedness are
#' weakly identified when the response is diffuse; holding them global
#' (\code{~ 1}) and keeping the concentration well identified is the robust
#' default.
#'
#' The location uses the Fisher-Lee tan-half link (\eqn{\xi \in (-\pi, \pi)},
#' antipode unrepresentable, winding number zero -- see \code{\link{pnlss}} when
#' the location must wind). The skewness and peakedness both ride the tanh link
#' (shared with \code{\link{ssjplss}} and \code{\link{vmftlss}}), bounded to
#' \eqn{(-1, 1)}.
#' @examples
#' library(mgcv)
#' set.seed(1)
#' n <- 300
#' x <- runif(n)
#' xi <- 2 * atan(sin(2 * pi * x))
#' kappa <- exp(1.0 + 0.5 * cos(2 * pi * x))
#' nu <- 0.5      # skew
#' lmbd <- 0.3    # peakedness
#' # Inverse Batschelet draws: the forward warps phi(u) are closed-form and
#' # monotone, so gridding the base angle u gives a (phi, density) grid to invert
#' # -- no root-finding needed for simulation.
#' u <- seq(-pi, pi, length.out = 4001)
#' phistar <- u - 0.5 * (1 + lmbd) * sin(u)
#' phi <- phistar - nu * (1 + cos(phistar))
#' A <- u - 0.5 * (1 - lmbd) * sin(u)
#' dev <- vapply(seq_len(n), function(i) {
#'   g <- exp(kappa[i] * cos(A))
#'   cdf <- cumsum((g[-1] + g[-length(g)]) / 2 * diff(phi))
#'   cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
#'   approx(cdf, phi, runif(1), rule = 2)$y
#' }, numeric(1))
#' y <- atan2(sin(xi + dev), cos(xi + dev))
#' # smooth location and concentration; global (intercept-only) skew and peakedness
#' b <- gam(list(y ~ s(x), ~ s(x), ~ 1, ~ 1), family = ibslss(), optimizer = "efs")
#' summary(b)
#' @references
#' Jones, M. C. and Pewsey, A. (2012) Inverse Batschelet distributions for
#' circular data. \emph{Biometrics} 68, 183-193.
#'
#' Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method for
#' smoothing parameter optimization with application to Tweedie location,
#' scale and shape models. \emph{Biometrics} 73, 1071-1081.
#' @seealso
#' \code{\link{vmlss}}, \code{\link{vmftlss}}, \code{\link{ssjplss}},
#' \code{\link{jplss}}, \code{\link{pnlss}}, \code{\link[mgcv]{gam}}
#' @export
ibslss <- function(link = list("tanhalf", "log", "tanh", "tanh")) {
  if (length(link) != 4)
    stop("ibslss requires 4 links specified as character strings")
  if (!(link[[1]] %in% "tanhalf"))
    stop(link[[1]], " link not available for the location parameter of ibslss")
  if (!(link[[2]] %in% "log"))
    stop(link[[2]], " link not available for the concentration parameter of ibslss")
  if (!(link[[3]] %in% "tanh"))
    stop(link[[3]], " link not available for the skewness parameter of ibslss")
  if (!(link[[4]] %in% "tanh"))
    stop(link[[4]], " link not available for the peakedness parameter of ibslss")

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
  ## skewness nu and peakedness lambda both ride the tanh link, bounded (-1, 1)
  stats[[3]] <- tanh_link()
  stats[[4]] <- tanh_link()

  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    type <- match.arg(type)
    xi <- object$fitted[, 1]
    kappa <- object$fitted[, 2]
    nu <- object$fitted[, 3]
    lambda <- object$fitted[, 4]
    d <- object$y - xi
    if (type == "response") return(atan2(sin(d), cos(d)))
    sind <- sin(d)
    if (type == "pearson") {
      ## Var{sin(y - xi)} ~ (1 - alpha2)/2 from the 2nd cosine moment of the
      ## warped kernel (no elementary form). Diagnostic standardization; with
      ## nu != 0 the law is skewed and xi the mode anchor, so this is not an
      ## exact residual variance (as in ssjplss).
      a2 <- vapply(seq_along(kappa),
                   function(i) .ibs_cos_moment(kappa[i], nu[i], lambda[i], 2),
                   numeric(1))
      v <- (1 - a2) / 2
      v[v <= 0] <- .Machine$double.eps
      return(sind / sqrt(v))
    }
    ## deviance: 2*(l_sat - l_hat); the normalizer cancels (l_sat is the peak
    ## log-density, value kappa, at the warped mode A = 0), leaving
    ## 2*kappa*(1 - cos A), A the warped kernel angle. sign(A) reduces to
    ## sign(d) for the von Mises member (nu = lambda = 0, A = d).
    w <- .ibs_warp_vec(object$y, xi, nu, lambda)
    A <- w$u_star - 0.5 * (1 - lambda) * sin(w$u_star)
    sign(A) * sqrt(pmax(2 * kappa * (1 - cos(A)), 0))
  }

  postproc <- expression({
    ## Deviance and null deviance against the common saturated reference (the
    ## per-obs peak log-density kappa + log c at the fitted kappa, lambda). The
    ## deviance itself is normalizer-free, 2*kappa*(1 - cos A); the null
    ## (constant von Mises MLE, the nu -> 0, lambda -> 0 baseline) needs the
    ## fitted normalizer, one grid pass via the internal helpers (fetched by
    ## name, as ll runs in mgcv's frame). All work in local() so mgcv's frames
    ## stay intact.
    .ibslss.dev <- local({
      lcv <- getFromNamespace(".ibs_log_c_vec", "circlss")
      warpf <- getFromNamespace(".ibs_warp_vec", "circlss")
      xi <- object$fitted[, 1]
      kappa <- object$fitted[, 2]
      nu <- object$fitted[, 3]
      lambda <- object$fitted[, 4]
      y <- object$y
      w <- warpf(y, xi, nu, lambda)
      A <- w$u_star - 0.5 * (1 - lambda) * sin(w$u_star)
      dev <- sum(2 * pmax(kappa * (1 - cos(A)), 0))
      lsat <- sum(kappa + lcv(kappa, lambda))
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
    object$deviance <- .ibslss.dev[1]
    object$null.deviance <- .ibslss.dev[2]
    rm(.ibslss.dev)
  })

  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                 D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
    ## Tier-2: EFS calls this with deriv <= 1; the gradient/Hessian assembly
    ## runs at gamlss order deriv - 1 = 0. No deriv > 1 branches.
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
    nu <- family$linfo[[3]]$linkinv(eta2)
    lambda <- family$linfo[[4]]$linkinv(eta3)

    ts <- .ibs_terms(y, xi, kappa, nu, lambda, deriv = deriv > 0)
    l0 <- ts$l0
    l <- sum(wt * l0)

    if (deriv) {
      ## l1 columns (xi, kappa, nu, lambda) and l2 columns (xx, xk, xn, xl, kk,
      ## kn, kl, nn, nl, ll) come straight from .ibs_terms in mgcv's
      ## combinations_with_replacement order.
      l1 <- ts$l1
      l2 <- ts$l2
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta0),
                   family$linfo[[2]]$mu.eta(eta1),
                   family$linfo[[3]]$mu.eta(eta2),
                   family$linfo[[4]]$mu.eta(eta3))
      g2 <- cbind(family$linfo[[1]]$d2link(xi),
                  family$linfo[[2]]$d2link(kappa),
                  family$linfo[[3]]$d2link(nu),
                  family$linfo[[4]]$d2link(lambda))
      ## no third/fourth log-likelihood or link derivatives for this family
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
    ## (.tanhalf_pilot_coef). Concentration starts from A1-inverse kappa0, skew
    ## from nu0 = 0 and peakedness from lambda0 = 0 (the von Mises member), as
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
        ## concentration (log), skew (tanh), peakedness (tanh): constant targets
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
    ## Simulation by inverse transform on a fine centered grid of the warped
    ## kernel, one grid per unique (kappa, nu, lambda) -- the warps have no
    ## closed form, so the density is built from .ibs_warp_vec (mirrors the
    ## ssjplss/vmftlss rd; the reference inverts a cached spectral cdf, the same
    ## idea).
    n <- nrow(mu)
    xiv <- mu[, 1]; kappa <- mu[, 2]; nu <- mu[, 3]; lambda <- mu[, 4]
    grid <- seq(-pi, pi, length.out = 2049L)
    th <- numeric(n)
    key <- paste(kappa, nu, lambda, sep = "\r")
    for (k in unique(key)) {
      idx <- which(key == k)
      ka <- kappa[idx[1]]; nv <- nu[idx[1]]; lm <- lambda[idx[1]]
      w <- .ibs_warp_vec(grid, 0.0, nv, lm)
      A <- w$u_star - 0.5 * (1 - lm) * sin(w$u_star)
      dens <- exp(ka * cos(A) - max(ka * cos(A)))
      dens[dens < 0] <- 0
      cdf <- cumsum((dens[-1] + dens[-length(dens)]) * diff(grid) / 2)
      cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
      u <- stats::runif(length(idx))
      th[idx] <- xiv[idx] +
        stats::approx(cdf, grid, xout = u, rule = 2, ties = "ordered")$y
    }
    atan2(sin(th), cos(th))
  }

  structure(list(family = "ibslss", ll = ll, link = paste(link), nlp = 4,
                 param_names = c("xi", "kappa", "nu", "lambda"),
                 param_circular = c(TRUE, FALSE, FALSE, FALSE),
                 sandwich = sandwich, tri = mgcv::trind.generator(4),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 0, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
