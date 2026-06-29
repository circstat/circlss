#' Cartwright power-of-cosine location-scale family
#'
#' A general family implementing distributional regression for a circular
#' (angular) response \eqn{y} in radians under Cartwright's power-of-cosine law
#' \deqn{f(y) = \frac{2^{1/\zeta - 1}\,\Gamma(1 + 1/\zeta)^2}
#'   {\pi\,\Gamma(1 + 2/\zeta)}\,\bigl(1 + \cos(y - \mu)\bigr)^{1/\zeta},
#'   \qquad \zeta > 0.}
#' Both the mean direction \eqn{\mu} and the peakedness \eqn{\zeta} get their
#' own linear predictor, each of which may contain smooth terms. Used with
#' \code{mgcv::gam} and a list of two formulas: the first names the response
#' and models \eqn{\mu}, the second models \eqn{\log\zeta}.
#'
#' @param link Two-element list of link names, for the mean direction and
#'   the peakedness parameter. Currently only the defaults are available:
#'   \code{"tanhalf"} for the location and \code{"log"} for the shape parameter
#'   \eqn{\zeta > 0}.
#' @return An object of class \code{c("general.family", "extended.family", "family")}
#'   for use with \code{\link[mgcv]{gam}} (or its front end \code{\link{circ_gam}}).
#' @details
#' Cartwright's distribution is a one-parameter \strong{peakedness} family:
#' \eqn{\zeta} raises \eqn{1 + \cos(y-\mu)} to the power \eqn{1/\zeta}, sharpening
#' or flattening the single mode. Its mean resultant length is
#' \eqn{\rho = 1/(\zeta + 1)}, so \eqn{\zeta \to 0} is sharply peaked
#' (\eqn{\rho \to 1}), \eqn{\zeta \to \infty} is the circular uniform
#' (\eqn{\rho \to 0}), and \eqn{\zeta = 1} is exactly the cardioid
#' (\code{\link{cardlss}}) at its concentration ceiling, \eqn{(1 + \cos)/2\pi}
#' (\eqn{\rho = 1/2}). It shares the power-of-cosine form with the Jones-Pewsey
#' family only at the latter's concentrated boundary; it is not an interior
#' Jones-Pewsey special case.
#'
#' The density is evaluated in the half-angle form
#' \eqn{1 + \cos d = 2\cos^2(d/2)}, which keeps the log-density exact near the
#' antipode \eqn{y = \mu \pm \pi}, where it has an \emph{honest zero} for every
#' \eqn{\zeta}. The second trigonometric moment is
#' \eqn{\alpha_2 = (1-\zeta)/\{(1+\zeta)(1+2\zeta)\}}, so Pearson residuals
#' standardize by \eqn{\mathrm{Var}\{\sin(y-\mu)\} = (1 - \alpha_2)/2}.
#'
#' \strong{Optimizer.} Only first- and second-order log-likelihood derivatives are
#' provided, so \code{available.derivs = 0} and the family is fitted by the extended
#' Fellner-Schall optimizer rather than full Newton REML. \code{\link[mgcv]{gam}}
#' selects \code{optimizer = "efs"} automatically; passing it explicitly is
#' recommended so the fitted object is labelled correctly (which
#' \code{\link[mgcv]{gam.check}} relies on). The \eqn{\zeta}-derivatives involve
#' the digamma and trigamma functions, because the normalizer is built from
#' \eqn{\Gamma(1 + 1/\zeta)} and \eqn{\Gamma(1 + 2/\zeta)} -- the first family
#' whose normalizing constant is not elementary.
#'
#' The mean direction uses the Fisher-Lee tan-half link
#' (\eqn{\mu \in (-\pi, \pi)}, antipode unrepresentable, winding number
#' zero -- see \code{\link{pnlss}} when the mean direction must wind).
#' @examples
#' library(mgcv)
#' set.seed(1)
#' n <- 400
#' x <- runif(n)
#' mu <- 2 * atan(sin(2 * pi * x))
#' zeta <- exp(-0.2 + 0.8 * cos(2 * pi * x))   # peakedness > 0
#' # Cartwright draws via the Beta-to-angle transform
#' ang <- 2 * asin(sqrt(rbeta(n, 0.5, 1 / zeta + 0.5)))
#' th <- mu + sample(c(-1, 1), n, replace = TRUE) * ang
#' y <- atan2(sin(th), cos(th))
#' b <- gam(list(y ~ s(x), ~ s(x)), family = cartlss(), optimizer = "efs")
#' summary(b)
#' @references
#' Cartwright, D. E. (1963) The use of directional spectra in studying the
#' output of a wave recorder on a moving ship. In \emph{Ocean Wave Spectra},
#' 203-218. Prentice-Hall.
#'
#' Jammalamadaka, S. R. and SenGupta, A. (2001) \emph{Topics in Circular
#' Statistics}. World Scientific.
#'
#' Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method for
#' smoothing parameter optimization with application to Tweedie location,
#' scale and shape models. \emph{Biometrics} 73, 1071-1081.
#' @seealso
#' \code{\link{cardlss}}, \code{\link{vmlss}}, \code{\link{wnlss}},
#' \code{\link{pnlss}}, \code{\link[mgcv]{gam}}
#' @export
cartlss <- function(link = list("tanhalf", "log")) {
  if (length(link) != 2) stop("cartlss requires 2 links specified as character strings")
  if (!(link[[1]] %in% "tanhalf"))
    stop(link[[1]], " link not available for the location parameter of cartlss")
  if (!(link[[2]] %in% "log"))
    stop(link[[2]], " link not available for the concentration parameter of cartlss")

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

  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    type <- match.arg(type)
    mu <- object$fitted[, 1]
    zeta <- object$fitted[, 2]
    d <- object$y - mu
    if (type == "response") return(atan2(sin(d), cos(d)))
    sind <- sin(d)
    if (type == "pearson") {
      ## second trig moment alpha2 = (1 - zeta)/((1 + zeta)(1 + 2 zeta)), so
      ## Var(sin(y - mu)) = (1 - alpha2)/2 (constant 1/2 at zeta = 1, the
      ## cardioid point; -> 0 as zeta -> 0, -> 1/2 as zeta -> inf)
      a2 <- (1 - zeta) / ((1 + zeta) * (1 + 2 * zeta))
      return(sind / sqrt((1 - a2) / 2))
    }
    ## deviance: 2*(peak log-density at d = 0 minus attained); all zeta-only
    ## terms cancel, leaving -(4/zeta) log|cos(d/2)| >= 0, signed by sin d
    sign(sind) * sqrt(pmax(-(4 / zeta) * log(abs(cos(0.5 * d))), 0))
  }

  postproc <- expression({
    ## Deviance pair against the common saturated reference (per-obs peak
    ## log-density at the fitted zeta), as in cardlss; null = the moment fit
    ## (mean direction, zeta0 = (1 - Rbar)/Rbar from rho = 1/(zeta + 1)). The
    ## density is closed form, so no namespace helper is needed. local() so
    ## mgcv's frames are safe.
    .cartlss.dev <- local({
      mu <- object$fitted[, 1]
      zeta <- object$fitted[, 2]
      y <- object$y
      n <- length(y)
      lpeak <- function(z) (2 / z - 1) * log(2) + 2 * lgamma(1 + 1 / z) -
        lgamma(1 + 2 / z) - log(pi)
      lcart <- function(d, z) lpeak(z) + (2 / z) * log(abs(cos(0.5 * d)))
      dev <- sum(-(4 / zeta) * log(abs(cos(0.5 * (y - mu)))))
      lsat <- sum(lpeak(zeta))
      sy <- mean(sin(y)); cy <- mean(cos(y))
      mu0 <- atan2(sy, cy)
      rbar <- min(max(sqrt(sy * sy + cy * cy), 0.05), 0.95)
      zeta0 <- (1 - rbar) / rbar
      lnull <- sum(lcart(y - mu0, zeta0))
      c(dev, 2 * (lsat - lnull))
    })
    object$deviance <- .cartlss.dev[1]
    object$null.deviance <- .cartlss.dev[2]
    rm(.cartlss.dev)
  })

  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                 D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
    ## Tier-2: only l1/l2 are available, so mgcv (EFS) calls this with
    ## deriv <= 1; the gradient/Hessian assembly runs at gamlss order
    ## deriv - 1 = 0. There are no deriv > 1 / deriv > 3 branches.
    if (!is.null(offset)) offset[[3]] <- 0
    if (is.null(wt)) wt <- 1            # prior weights; sandwich/direct callers may pass NULL
    jj <- attr(X, "lpi")
    if (is.null(eta)) {
      eta <- X[, jj[[1]], drop = FALSE] %*% coef[jj[[1]]]
      if (!is.null(offset[[1]])) eta <- eta + offset[[1]]
      eta1 <- X[, jj[[2]], drop = FALSE] %*% coef[jj[[2]]]
      if (!is.null(offset[[2]])) eta1 <- eta1 + offset[[2]]
    } else {
      eta1 <- eta[, 2]
      eta <- eta[, 1]
    }
    eta <- drop(eta); eta1 <- drop(eta1)
    mu <- family$linfo[[1]]$linkinv(eta)
    zeta <- family$linfo[[2]]$linkinv(eta1)
    d <- y - mu
    inv <- 1 / zeta
    logcos <- log(abs(cos(0.5 * d)))
    ## log normalizer (gammaln form, overflow-free) + half-angle data term
    l0 <- (2 * inv - 1) * log(2) + 2 * lgamma(1 + inv) - lgamma(1 + 2 * inv) -
      log(pi) + 2 * inv * logcos
    l <- sum(wt * l0)

    if (deriv) {
      t <- tan(0.5 * d)
      L <- log(2) + 2 * logcos                 # log(1 + cos d)
      ## B = d/d(1/zeta) of the log normalizer's data-independent part, minus L
      B <- -log(2) - 2 * digamma(1 + inv) + 2 * digamma(1 + 2 * inv) - L
      l_zeta <- B * inv * inv
      ## l1: d l / d(mu, zeta); l2 columns ordered (mm, mr, rr)
      l1 <- cbind(t * inv, l_zeta)
      l2 <- cbind(-0.5 * (1 + t * t) * inv,
                  -t * inv * inv,
                  (2 * trigamma(1 + inv) - 4 * trigamma(1 + 2 * inv)) * inv^4 -
                    2 * inv * l_zeta)
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),
                   family$linfo[[2]]$mu.eta(eta1))
      g2 <- cbind(family$linfo[[1]]$d2link(mu),
                  family$linfo[[2]]$d2link(zeta))
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
    ## Cartwright's density has an honest zero at the antipode for every zeta,
    ## so a flat mean-direction start strands antipodal observations on the
    ## likelihood's cliffs (infinite tan(d/2) gradients) and the fit collapses
    ## to a constant mu. So the location gets a DATA-FOLLOWING projected pilot
    ## -- smooth sin(y) and cos(y) over the mu design, take mu-hat = atan2 of
    ## the two, and solve for the coefficients reproducing tan(mu-hat/2). (The
    ## earlier closed-form families have no density zero and start mu flat; the
    ## reference uses a wrapped-normal mu pilot here for the same reason.) The
    ## peakedness keeps the constant moment start zeta0 = (1 - Rbar)/Rbar from
    ## the relation Rbar = 1/(zeta + 1). Everything except the contract
    ## variables stays inside local().
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
        ## location: projected pilot through the tan-half link
        Xmu <- x[, jj[[1]], drop = FALSE]
        Emu <- E[, jj[[1]], drop = FALSE]
        off1 <- if (!is.null(offset) && length(offset) >= 1 &&
                    !is.null(offset[[1]])) offset[[1]] else 0
        shat <- drop(Xmu %*% pls(Xmu, Emu, sin(y)))
        chat <- drop(Xmu %*% pls(Xmu, Emu, cos(y)))
        eta_mu <- pmax(pmin(tan(0.5 * atan2(shat, chat)), 1e3), -1e3) - off1
        eta_mu[!is.finite(eta_mu)] <- 0
        st[jj[[1]]] <- pls(Xmu, Emu, eta_mu)
        ## concentration: constant moment start through the log link
        sy <- mean(sin(y)); cy <- mean(cos(y))
        rbar <- min(max(sqrt(sy * sy + cy * cy), 0.05), 0.95)
        tz <- family$linfo[[2]]$linkfun((1 - rbar) / rbar)
        off2 <- if (!is.null(offset) && length(offset) >= 2 &&
                    !is.null(offset[[2]])) offset[[2]] else 0
        st[jj[[2]]] <- pls(x[, jj[[2]], drop = FALSE],
                           E[, jj[[2]], drop = FALSE], rep(tz, nobs) - off2)
        st
      })
    }
  })

  rd <- function(mu, wt, scale) {
    ## exact: the Beta-to-angle transform.
    ## Draw T ~ Beta(1/2, 1/zeta + 1/2), set the deviation angle to
    ## 2*asin(sqrt(T)), and reflect it +/- around mu with probability 1/2.
    n <- nrow(mu)
    muv <- mu[, 1]
    zeta <- pmax(mu[, 2], 1e-8)
    Tb <- stats::rbeta(n, 0.5, 1 / zeta + 0.5)
    ang <- 2 * asin(sqrt(pmin(Tb, 1)))
    sgn <- ifelse(stats::runif(n) < 0.5, -1, 1)
    th <- muv + sgn * ang
    atan2(sin(th), cos(th))
  }

  structure(list(family = "cartlss", ll = ll, link = paste(link), nlp = 2,
                 param_names = c("mu", "zeta"),
                 param_circular = c(TRUE, FALSE),
                 sandwich = sandwich, tri = mgcv::trind.generator(2),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 0, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
