#' Sine-skewed Jones-Pewsey location-concentration-shape-skewness family
#'
#' A general family implementing distributional regression for a circular
#' (angular) response \eqn{y} in radians under the sine-skewed Jones-Pewsey law
#' \deqn{f(y) = c(\kappa, \psi)\,\bigl(\cosh(\kappa\psi) +
#'   \sinh(\kappa\psi)\cos(y - \xi)\bigr)^{1/\psi}\,
#'   \bigl(1 + \lambda \sin(y - \xi)\bigr),}
#' the Jones-Pewsey family (\code{\link{jplss}}) multiplied by a sine-skew factor
#' (Umbach-Jammalamadaka). It adds an asymmetry axis to the symmetric Jones-Pewsey
#' umbrella.
#'
#' It is the first family with \strong{four linear predictors}: the location
#' \eqn{\xi}, the concentration \eqn{\kappa}, the shape \eqn{\psi} and the
#' skewness \eqn{\lambda} each get their own, any of which may contain smooth
#' terms. Used with \code{mgcv::gam} and a list of \emph{four} formulas; the first
#' names the response and models \eqn{\xi}, then \eqn{\log\kappa}, then
#' \eqn{\psi}, then \eqn{\lambda}. The shape and skewness are most often held
#' global, i.e. fitted intercept-only with \code{~ 1}.
#'
#' @param link Four-element list of link names, for the location, concentration,
#'   shape and skewness. Currently only the defaults are available:
#'   \code{"tanhalf"} for the location \eqn{\xi}, \code{"log"} for the
#'   concentration \eqn{\kappa > 0}, \code{"identity"} for the shape
#'   \eqn{\psi \in \mathbb{R}}, and \code{"tanh"} for the skewness
#'   \eqn{\lambda \in (-1, 1)}.
#' @return An object of class \code{c("general.family", "extended.family", "family")}
#'   for use with \code{\link[mgcv]{gam}} (or its front end \code{\link{circ_gam}}).
#' @details
#' The sine-skew factor \eqn{1 + \lambda\sin(y-\xi)} integrates to 1 against the
#' symmetric Jones-Pewsey kernel (the first trigonometric moment of the centered
#' kernel is zero), so it leaves the Jones-Pewsey normalizer \eqn{c(\kappa,\psi)}
#' \emph{unchanged}. At \eqn{\lambda = 0} the family reduces exactly to
#' \code{\link{jplss}}. The skewness rides the \strong{tanh} link
#' \eqn{\lambda = \tanh(\eta)}; as \eqn{|\lambda| \to 1} the density touches 0
#' where \eqn{\lambda\sin(y-\xi) = -1}, so the link keeps \eqn{\lambda} strictly
#' interior.
#'
#' \strong{Mode anchor, not mean.} Once \eqn{\lambda \neq 0}, the location
#' \eqn{\xi} is the \emph{mode anchor} of the asymmetric density, not its mean
#' direction; fitted-direction summaries inherit that reading.
#'
#' \strong{Normalizer and derivatives.} Because the normalizer is the
#' Jones-Pewsey one, the family reuses that family's Gauss-Legendre quadrature
#' machinery wholesale: the \eqn{\kappa}- and \eqn{\psi}-score and Hessian are
#' exactly the Jones-Pewsey terms, the \eqn{(\kappa,\lambda)} and
#' \eqn{(\psi,\lambda)} cross-derivatives vanish identically, and only the
#' \eqn{\xi}- and \eqn{\lambda}-directions carry the (elementary) skew terms.
#'
#' \strong{Optimizer.} Only first- and second-order log-likelihood derivatives are
#' provided, so \code{available.derivs = 0} and the family is fitted by the extended
#' Fellner-Schall optimizer rather than full Newton REML. \code{\link[mgcv]{gam}}
#' selects \code{optimizer = "efs"} automatically; passing it explicitly is
#' recommended so the fitted object is labelled correctly (which
#' \code{\link[mgcv]{gam.check}} relies on). The shape \eqn{\psi} and skew
#' \eqn{\lambda} are weakly identified when the response is diffuse; holding them
#' global (\code{~ 1}) and keeping the concentration well identified is the robust
#' default.
#'
#' The location uses the Fisher-Lee tan-half link (\eqn{\xi \in (-\pi, \pi)},
#' antipode unrepresentable, winding number zero -- see \code{\link{pnlss}} when
#' the location must wind).
#' @examples
#' library(mgcv)
#' set.seed(1)
#' n <- 300
#' x <- runif(n)
#' xi <- 2 * atan(sin(2 * pi * x))
#' kappa <- exp(1.0 + 0.4 * cos(2 * pi * x))
#' psi <- 0.5
#' lambda <- 0.5
#' # sine-skewed Jones-Pewsey draws by inverse transform on the centered kernel
#' phi <- seq(-pi, pi, length.out = 4001)
#' dev <- vapply(seq_len(n), function(i) {
#'   g <- (cosh(kappa[i] * psi) + sinh(kappa[i] * psi) * cos(phi))^(1 / psi) *
#'     (1 + lambda * sin(phi))
#'   g[g < 0] <- 0
#'   cdf <- cumsum(g); cdf <- cdf / cdf[length(cdf)]
#'   approx(cdf, phi, runif(1), rule = 2)$y
#' }, numeric(1))
#' y <- atan2(sin(xi + dev), cos(xi + dev))
#' # smooth location and concentration; global (intercept-only) shape and skew
#' b <- gam(list(y ~ s(x), ~ s(x), ~ 1, ~ 1), family = ssjplss(),
#'          optimizer = "efs")
#' summary(b)
#' @references
#' Umbach, D. and Jammalamadaka, S. R. (2009) Building asymmetry into circular
#' distributions. \emph{Statistics & Probability Letters} 79, 659-663.
#'
#' Abe, T. and Pewsey, A. (2011) Sine-skewed circular distributions.
#' \emph{Statistical Papers} 52, 683-707.
#'
#' Jones, M. C. and Pewsey, A. (2005) A family of symmetric distributions on the
#' circle. \emph{Journal of the American Statistical Association} 100, 1422-1428.
#'
#' Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method for
#' smoothing parameter optimization with application to Tweedie location,
#' scale and shape models. \emph{Biometrics} 73, 1071-1081.
#' @seealso
#' \code{\link{jplss}}, \code{\link{vmlss}}, \code{\link{cardlss}},
#' \code{\link{wclss}}, \code{\link{cartlss}}, \code{\link{pnlss}},
#' \code{\link[mgcv]{gam}}
#' @export
ssjplss <- function(link = list("tanhalf", "log", "identity", "tanh")) {
  if (length(link) != 4)
    stop("ssjplss requires 4 links specified as character strings")
  if (!(link[[1]] %in% "tanhalf"))
    stop(link[[1]], " link not available for the location parameter of ssjplss")
  if (!(link[[2]] %in% "log"))
    stop(link[[2]], " link not available for the concentration parameter of ssjplss")
  if (!(link[[3]] %in% "identity"))
    stop(link[[3]], " link not available for the shape parameter of ssjplss")
  if (!(link[[4]] %in% "tanh"))
    stop(link[[4]], " link not available for the skewness parameter of ssjplss")

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
    lambda <- object$fitted[, 4]
    d <- object$y - xi
    if (type == "response") return(atan2(sin(d), cos(d)))
    sind <- sin(d)
    if (type == "pearson") {
      ## base (symmetric) JP scale; xi is the mode anchor, so this is a
      ## diagnostic standardization, not an exact residual variance
      a2 <- vapply(seq_along(kappa),
                   function(i) .jp_cos_moment(kappa[i], psi[i], 2), numeric(1))
      v <- (1 - a2) / 2
      v[v <= 0] <- .Machine$double.eps
      return(sind / sqrt(v))
    }
    ## deviance: 2*(l_sat - l_hat); the normalizer cancels, leaving the
    ## skew-tilted peak (grid max) minus the attained centered log-density
    h <- .jp_score_terms(d, kappa, psi, second = FALSE)$h
    lhat <- h + log1p(lambda * sind)
    lsat <- .ssjp_lsat(kappa, psi, lambda)
    sign(sind) * sqrt(pmax(2 * (lsat - lhat), 0))
  }

  postproc <- expression({
    ## Deviance and null deviance against the common saturated reference (the
    ## skew-tilted per-obs peak). The deviance is normalizer-free; the null
    ## (constant von Mises MLE, psi -> 0, lambda -> 0) needs the fitted
    ## normalizer, one quadrature pass via the internal helpers. local() so
    ## mgcv's frames stay intact.
    .ssjplss.dev <- local({
      st <- getFromNamespace(".jp_score_terms", "circlss")
      lcv <- getFromNamespace(".jp_logc_vec", "circlss")
      lsatf <- getFromNamespace(".ssjp_lsat", "circlss")
      xi <- object$fitted[, 1]
      kappa <- object$fitted[, 2]
      psi <- object$fitted[, 3]
      lambda <- object$fitted[, 4]
      y <- object$y
      d <- y - xi
      h <- st(d, kappa, psi, second = FALSE)$h
      lsat0 <- lsatf(kappa, psi, lambda)               # normalizer-free peak
      lhat0 <- h + log1p(lambda * sin(d))              # normalizer-free attained
      dev <- sum(2 * pmax(lsat0 - lhat0, 0))
      logc <- lcv(kappa, psi)
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
    object$deviance <- .ssjplss.dev[1]
    object$null.deviance <- .ssjplss.dev[2]
    rm(.ssjplss.dev)
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
    lambda <- family$linfo[[4]]$linkinv(eta3)
    d <- y - xi
    s <- sin(d); cc <- cos(d)
    Q <- 1 + lambda * s

    ts <- .jp_score_terms(d, kappa, psi, second = deriv > 0)
    Qd <- .jp_quad_vec(kappa, psi)          # n x 6: logc, dk, dp, dkk, dkp, dpp
    ## l0 VALUE: jplss log-density (with uniform / von Mises shortcuts) + log Q
    logc <- Qd[, 1]
    vm <- (kappa >= .JP_KAPPA_TOL) & (abs(psi) < .JP_PSI_TOL)
    if (any(vm))
      logc[vm] <- -(log(2 * pi * besselI(kappa[vm], 0, expon.scaled = TRUE)) +
                      kappa[vm])
    h0 <- ts$h
    if (any(vm)) h0[vm] <- kappa[vm] * cos(d[vm])
    l0 <- h0 + logc
    unif <- kappa < .JP_KAPPA_TOL
    if (any(unif)) l0[unif] <- -log(2 * pi)
    l0 <- l0 + log1p(lambda * s)
    l <- sum(wt * l0)

    if (deriv) {
      Q2 <- Q * Q
      z <- numeric(length(d))
      ## l1: d l / d(xi, kappa, psi, lambda)
      l1 <- cbind(-lambda * cc / Q - ts$hphi,
                  ts$hk - Qd[, 2],
                  ts$hp - Qd[, 3],
                  s / Q)
      ## l2 columns in combinations_with_replacement order:
      ## (xx, xk, xp, xl, kk, kp, kl, pp, pl, ll)
      l2 <- cbind(-lambda * (s * Q + lambda * cc * cc) / Q2 + ts$hphiphi,
                  -ts$hphik,
                  -ts$hphip,
                  -cc / Q2,
                  ts$hkk - Qd[, 4],
                  ts$hkp - Qd[, 5],
                  z,
                  ts$hpp - Qd[, 6],
                  z,
                  -s * s / Q2)
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta0),
                   family$linfo[[2]]$mu.eta(eta1),
                   family$linfo[[3]]$mu.eta(eta2),
                   family$linfo[[4]]$mu.eta(eta3))
      g2 <- cbind(family$linfo[[1]]$d2link(xi),
                  family$linfo[[2]]$d2link(kappa),
                  family$linfo[[3]]$d2link(psi),
                  family$linfo[[4]]$d2link(lambda))
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
    ## from psi0 = 0 (von Mises) and skew from lambda0 = 0 (symmetric), as
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
        ## concentration (log), shape (identity), skew (tanh): constant targets
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
    ## Simulation by inverse transform on a fine centered grid of the SKEWED
    ## kernel, one grid per unique (kappa, psi, lambda).
    n <- nrow(mu)
    xiv <- mu[, 1]; kappa <- mu[, 2]; psi <- mu[, 3]; lambda <- mu[, 4]
    grid <- seq(-pi, pi, length.out = 2049L)
    sg <- sin(grid)
    th <- numeric(n)
    key <- paste(kappa, psi, lambda, sep = "\r")
    for (k in unique(key)) {
      idx <- which(key == k)
      h <- .jp_score_terms(grid, kappa[idx[1]], psi[idx[1]], second = FALSE)$h
      dens <- exp(h - max(h)) * (1 + lambda[idx[1]] * sg)
      dens[dens < 0] <- 0
      cdf <- cumsum((dens[-1] + dens[-length(dens)]) * diff(grid) / 2)
      cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
      u <- stats::runif(length(idx))
      th[idx] <- xiv[idx] +
        stats::approx(cdf, grid, xout = u, rule = 2, ties = "ordered")$y
    }
    atan2(sin(th), cos(th))
  }

  structure(list(family = "ssjplss", ll = ll, link = paste(link), nlp = 4,
                 param_names = c("xi", "kappa", "psi", "lambda"),
                 param_circular = c(TRUE, FALSE, FALSE, FALSE),
                 sandwich = sandwich, tri = mgcv::trind.generator(4),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 0, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
