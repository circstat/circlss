#' Flat-topped von Mises location-concentration-shape family
#'
#' A general family implementing distributional regression for a circular
#' (angular) response \eqn{y} in radians under the flat-topped von Mises law
#' \deqn{f(y) = c(\kappa, \nu)\,\exp\!\bigl(\kappa\cos(\phi + \nu\sin\phi)\bigr),
#'   \qquad \phi = y - \mu,\ \kappa > 0,\ \nu \in (-1, 1),}
#' with \eqn{c(\kappa,\nu)} a normalizing constant. The peakedness \eqn{\nu}
#' warps the angle forward inside the cosine: the warp \eqn{B = \phi +
#' \nu\sin\phi} is odd in \eqn{\phi}, so the density stays symmetric about
#' \eqn{\mu} while the peak is sharpened (\eqn{\nu > 0}) or flattened
#' (\eqn{\nu < 0}) -- a flat-topped or sharply-peaked alternative to the von
#' Mises that keeps the same mean direction.
#'
#' It has \strong{three linear predictors}: the mean direction \eqn{\mu}, the
#' concentration \eqn{\kappa}, and the peakedness \eqn{\nu} each get their own,
#' any of which may contain smooth terms. Used with \code{mgcv::gam} and a list
#' of \emph{three} formulas; the first names the response and models \eqn{\mu},
#' the second models \eqn{\log\kappa}, and the third models \eqn{\nu}. The
#' peakedness is most often held global, i.e. fitted intercept-only with
#' \code{~ 1}.
#'
#' @param link Three-element list of link names, for the mean direction, the
#'   concentration, and the peakedness. Currently only the defaults are available:
#'   \code{"tanhalf"} for the location \eqn{\mu}, \code{"log"} for the
#'   concentration \eqn{\kappa > 0}, and \code{"tanh"} for the peakedness
#'   \eqn{\nu \in (-1, 1)}.
#' @return An object of class \code{c("general.family", "extended.family", "family")}
#'   for use with \code{\link[mgcv]{gam}} (or its front end \code{\link{circ_gam}}).
#' @details
#' The peakedness \eqn{\nu} indexes a family of symmetric circular laws and
#' recovers the von Mises exactly:
#' \itemize{
#'   \item \eqn{\nu = 0}: the von Mises (\code{\link{vmlss}}),
#'     \eqn{\propto \exp(\kappa\cos(y-\mu))};
#'   \item \eqn{\nu > 0}: a sharper-than-von-Mises peak;
#'   \item \eqn{\nu < 0}: a flatter, flat-topped peak;
#'   \item \eqn{\kappa \to 0}: the circular uniform.
#' }
#'
#' \strong{Normalizer.} The constant \eqn{c(\kappa, \nu)} has no elementary
#' form. The log-density value \emph{and} the \eqn{\kappa}- and \eqn{\nu}-score
#' and Hessian carry kernel-weighted moments of the normalizer, evaluated by an
#' adaptive equispaced trapezoid (the warped kernel is smooth and periodic, so
#' the trapezoid is spectrally accurate). The grid resolution grows with
#' \eqn{\kappa} and \eqn{|\nu|} so the sharpened peak is always resolved.
#'
#' \strong{Optimizer.} Only first- and second-order log-likelihood derivatives
#' are provided, so \code{available.derivs = 0} and the family is fitted by the
#' extended Fellner-Schall optimizer rather than full Newton REML.
#' \code{\link[mgcv]{gam}} selects \code{optimizer = "efs"} automatically; passing
#' it explicitly is recommended so the fitted object is labelled correctly (which
#' \code{\link[mgcv]{gam.check}} relies on). The peakedness \eqn{\nu} is weakly
#' identified when the response is diffuse; holding it global (\code{~ 1}) and
#' keeping the concentration well identified is the robust default.
#'
#' The mean direction uses the Fisher-Lee tan-half link
#' (\eqn{\mu \in (-\pi, \pi)}, antipode unrepresentable, winding number
#' zero -- see \code{\link{pnlss}} when the mean direction must wind). The
#' peakedness rides the tanh link (shared with \code{\link{ssjplss}}'s
#' skewness), bounded to \eqn{(-1, 1)}.
#' @examples
#' library(mgcv)
#' set.seed(1)
#' n <- 300
#' x <- runif(n)
#' mu <- 2 * atan(sin(2 * pi * x))
#' kappa <- exp(1.0 + 0.5 * cos(2 * pi * x))
#' nu <- 0.5
#' # flat-topped vM draws by inverse transform on the centered kernel, then shift
#' phi <- seq(-pi, pi, length.out = 4001)
#' dev <- vapply(seq_len(n), function(i) {
#'   B <- phi + nu * sin(phi)
#'   g <- exp(kappa[i] * cos(B))
#'   cdf <- cumsum(g); cdf <- cdf / cdf[length(cdf)]
#'   approx(cdf, phi, runif(1), rule = 2)$y
#' }, numeric(1))
#' y <- atan2(sin(mu + dev), cos(mu + dev))
#' # smooth mean direction and concentration; global (intercept-only) peakedness
#' b <- gam(list(y ~ s(x), ~ s(x), ~ 1), family = vmftlss(), optimizer = "efs")
#' summary(b)
#' @references
#' Batschelet, E. (1981) \emph{Circular Statistics in Biology}. Academic Press.
#'
#' Abe, T., Pewsey, A. and Shimizu, K. (2013) Extending circular distributions
#' through transformation of argument. \emph{Annals of the Institute of
#' Statistical Mathematics} 65, 833-858.
#'
#' Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method for
#' smoothing parameter optimization with application to Tweedie location,
#' scale and shape models. \emph{Biometrics} 73, 1071-1081.
#' @seealso
#' \code{\link{vmlss}}, \code{\link{jplss}}, \code{\link{ssjplss}},
#' \code{\link{pnlss}}, \code{\link[mgcv]{gam}}
#' @export
vmftlss <- function(link = list("tanhalf", "log", "tanh")) {
  if (length(link) != 3)
    stop("vmftlss requires 3 links specified as character strings")
  if (!(link[[1]] %in% "tanhalf"))
    stop(link[[1]], " link not available for the location parameter of vmftlss")
  if (!(link[[2]] %in% "log"))
    stop(link[[2]], " link not available for the concentration parameter of vmftlss")
  if (!(link[[3]] %in% "tanh"))
    stop(link[[3]], " link not available for the shape parameter of vmftlss")

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
  ## peakedness nu rides the tanh link, bounded to (-1, 1)
  stats[[3]] <- tanh_link()

  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    type <- match.arg(type)
    mu <- object$fitted[, 1]
    kappa <- object$fitted[, 2]
    nu <- object$fitted[, 3]
    d <- object$y - mu
    if (type == "response") return(atan2(sin(d), cos(d)))
    sind <- sin(d)
    if (type == "pearson") {
      ## Var{sin(y - mu)} = (1 - alpha2)/2, with the 2nd cosine moment alpha2
      ## from the trapezoid grid (no elementary form for flat-topped vM)
      a2 <- vapply(seq_along(kappa),
                   function(i) .vmft_cos_moment(kappa[i], nu[i], 2), numeric(1))
      v <- (1 - a2) / 2
      v[v <= 0] <- .Machine$double.eps
      return(sind / sqrt(v))
    }
    ## deviance: 2*(l_sat - l_hat); the normalizer cancels (l_sat is the peak
    ## log-density at the fitted kappa, at the mode phi = 0 where B = 0, value
    ## kappa), leaving 2*kappa*(1 - cos B), B = d + nu sin d
    B <- d + nu * sind
    sign(sind) * sqrt(pmax(2 * kappa * (1 - cos(B)), 0))
  }

  postproc <- expression({
    ## Deviance and null deviance against the common saturated reference (the
    ## per-obs peak log-density kappa + log c at the fitted kappa, nu). The
    ## deviance itself is normalizer-free, 2*kappa*(1 - cos B); the null
    ## (constant von Mises MLE, the nu -> 0 baseline) needs the fitted
    ## normalizer, one trapezoid pass via the internal helper (fetched by name,
    ## as ll runs in mgcv's frame). All work in local() so mgcv's frames stay
    ## intact.
    .vmftlss.dev <- local({
      lcv <- getFromNamespace(".vmft_log_c_vec", "circlss")
      mu <- object$fitted[, 1]
      kappa <- object$fitted[, 2]
      nu <- object$fitted[, 3]
      y <- object$y
      d <- y - mu
      B <- d + nu * sin(d)
      dev <- sum(2 * pmax(kappa * (1 - cos(B)), 0))
      lsat <- sum(kappa + lcv(kappa, nu))
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
    object$deviance <- .vmftlss.dev[1]
    object$null.deviance <- .vmftlss.dev[2]
    rm(.vmftlss.dev)
  })

  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                 D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
    ## Tier-2: EFS calls this with deriv <= 1; the gradient/Hessian assembly
    ## runs at gamlss order deriv - 1 = 0. No deriv > 1 branches.
    if (is.null(offset)) offset <- vector("list", 4) else offset[[4]] <- 0
    if (is.null(wt)) wt <- 1            # prior weights; sandwich/direct callers may pass NULL
    jj <- attr(X, "lpi")
    if (is.null(eta)) {
      eta0 <- X[, jj[[1]], drop = FALSE] %*% coef[jj[[1]]]
      if (!is.null(offset[[1]])) eta0 <- eta0 + offset[[1]]
      eta1 <- X[, jj[[2]], drop = FALSE] %*% coef[jj[[2]]]
      if (!is.null(offset[[2]])) eta1 <- eta1 + offset[[2]]
      eta2 <- X[, jj[[3]], drop = FALSE] %*% coef[jj[[3]]]
      if (!is.null(offset[[3]])) eta2 <- eta2 + offset[[3]]
    } else {
      eta0 <- eta[, 1]; eta1 <- eta[, 2]; eta2 <- eta[, 3]
    }
    eta0 <- drop(eta0); eta1 <- drop(eta1); eta2 <- drop(eta2)
    mu <- family$linfo[[1]]$linkinv(eta0)
    kappa <- family$linfo[[2]]$linkinv(eta1)
    nu <- family$linfo[[3]]$linkinv(eta2)

    ts <- .vmft_terms(y, mu, kappa, nu, deriv = deriv > 0)
    l0 <- ts$l0
    l <- sum(wt * l0)

    if (deriv) {
      ## l1 columns (mu, kappa, nu) and l2 columns (mm, mk, mn, kk, kn, nn) come
      ## straight from .vmft_terms in mgcv's combinations_with_replacement order
      l1 <- ts$l1
      l2 <- ts$l2
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta0),
                   family$linfo[[2]]$mu.eta(eta1),
                   family$linfo[[3]]$mu.eta(eta2))
      g2 <- cbind(family$linfo[[1]]$d2link(mu),
                  family$linfo[[2]]$d2link(kappa),
                  family$linfo[[3]]$d2link(nu))
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
    ## Location starts from the data-following projected mean-direction pilot
    ## (.tanhalf_pilot_coef). Concentration
    ## starts from Fisher's A1-inverse kappa0 and the peakedness from nu0 = 0
    ## (the von Mises member), as constant targets through their links. local()
    ## so mgcv's frames (gam.fit5, initial.spg) are never clobbered.
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
        ## concentration (log) and peakedness (tanh): constant targets
        targ <- c(NA, family$linfo[[2]]$linkfun(kappa0),
                  family$linfo[[3]]$linkfun(0))
        for (k in 2:3) {
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
    ## Simulation by inverse transform on a fine centered grid of the kernel,
    ## one grid per unique (kappa, nu) -- adequate for predictive draws at
    ## fitted parameters (mirrors jplss rd; the reference inverts a cached
    ## spectral cdf table, the same idea).
    n <- nrow(mu)
    muv <- mu[, 1]
    kappa <- mu[, 2]
    nu <- mu[, 3]
    grid <- seq(-pi, pi, length.out = 2049L)
    sg <- sin(grid)
    th <- numeric(n)
    key <- paste(kappa, nu, sep = "\r")
    for (k in unique(key)) {
      idx <- which(key == k)
      ka <- kappa[idx[1]]
      nv <- nu[idx[1]]
      B <- grid + nv * sg
      dens <- exp(ka * cos(B) - max(ka * cos(B)))
      cdf <- cumsum((dens[-1] + dens[-length(dens)]) * diff(grid) / 2)
      cdf <- c(0, cdf)
      cdf <- cdf / cdf[length(cdf)]
      u <- stats::runif(length(idx))
      phi <- stats::approx(cdf, grid, xout = u, rule = 2, ties = "ordered")$y
      th[idx] <- muv[idx] + phi
    }
    atan2(sin(th), cos(th))
  }

  structure(list(family = "vmftlss", ll = ll, link = paste(link), nlp = 3,
                 param_names = c("mu", "kappa", "nu"),
                 param_circular = c(TRUE, FALSE, FALSE),
                 sandwich = sandwich, tri = mgcv::trind.generator(3),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 0, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
