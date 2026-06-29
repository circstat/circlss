#' Cardioid location-scale family
#'
#' A general family implementing distributional regression for a circular
#' (angular) response \eqn{y} in radians under the cardioid law
#' \deqn{f(y) = \frac{1}{2\pi}\left(1 + 2\rho\cos(y - \mu)\right),
#'   \qquad 0 \le \rho \le 1/2,}
#' a first-harmonic perturbation of the circular uniform. Both the mean
#' direction \eqn{\mu} and the mean resultant length \eqn{\rho} get their own
#' linear predictor, each of which may contain smooth terms. Used with
#' \code{mgcv::gam} and a list of two formulas: the first names the response
#' and models \eqn{\mu}, the second models \eqn{\mathrm{logithalf}(\rho)}.
#'
#' @param link Two-element list of link names, for the mean direction and
#'   the mean resultant length. Currently only the defaults are available:
#'   \code{"tanhalf"} for the location and \code{"logithalf"} for the
#'   concentration parameter \eqn{\rho \in (0, 1/2)}.
#' @return An object of class \code{c("general.family", "extended.family", "family")}
#'   for use with \code{\link[mgcv]{gam}} (or its front end \code{\link{circ_gam}}).
#' @details
#' The cardioid is the simplest departure from circular uniformity: a single
#' cosine ripple of amplitude \eqn{2\rho} on the flat density. It is therefore
#' a \strong{low-concentration / near-uniform} law -- \eqn{\rho = 0} is exactly
#' uniform, and concentration is capped at \eqn{\rho = 1/2}, beyond which the
#' density would go negative at the antimode. That hard upper bound is why the
#' concentration uses the \strong{logit-half} link
#' \eqn{\eta = \log\{\rho / (1/2 - \rho)\}}, i.e.
#' \eqn{\rho = \tfrac{1}{2}\,\mathrm{plogis}(\eta) \in (0, 1/2)} -- the one new
#' link this family brings (the von Mises, wrapped Cauchy and wrapped normal
#' all use the ordinary logit on \eqn{(0, 1)}).
#'
#' The cardioid is a pure first-harmonic distribution: its second and higher
#' trigonometric moments vanish, so the first moment is \eqn{\rho} (hence the
#' moment estimator \eqn{\hat\rho = \bar R}) and Pearson residuals standardize
#' by the constant \eqn{\mathrm{Var}\{\sin(y-\mu)\} = 1/2}.
#'
#' \strong{Optimizer.} Only first- and second-order log-likelihood derivatives are
#' provided, so \code{available.derivs = 0} and the family is fitted by the extended
#' Fellner-Schall optimizer rather than full Newton REML. \code{\link[mgcv]{gam}}
#' selects \code{optimizer = "efs"} automatically; passing it explicitly is
#' recommended so the fitted object is labelled correctly (which
#' \code{\link[mgcv]{gam.check}} relies on).
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
#' rho <- 0.5 * plogis(0.3 + 1.2 * cos(2 * pi * x))   # mean resultant length < 1/2
#' # cardioid draws by rejection from a uniform envelope
#' phi <- runif(n, -pi, pi)
#' keep <- runif(n) <= (1 + 2 * rho * cos(phi)) / (1 + 2 * rho)
#' while (any(!keep)) {
#'   j <- which(!keep)
#'   phi[j] <- runif(length(j), -pi, pi)
#'   keep[j] <- runif(length(j)) <= (1 + 2 * rho[j] * cos(phi[j])) / (1 + 2 * rho[j])
#' }
#' y <- atan2(sin(mu + phi), cos(mu + phi))
#' b <- gam(list(y ~ s(x), ~ s(x)), family = cardlss(), optimizer = "efs")
#' summary(b)
#' @references
#' Jammalamadaka, S. R. and SenGupta, A. (2001) \emph{Topics in Circular
#' Statistics}. World Scientific.
#'
#' Pewsey, A., Neuhaeuser, M. and Ruxton, G. D. (2013) \emph{Circular
#' Statistics in R}. Oxford University Press.
#'
#' Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method for
#' smoothing parameter optimization with application to Tweedie location,
#' scale and shape models. \emph{Biometrics} 73, 1071-1081.
#' @seealso
#' \code{\link{vmlss}}, \code{\link{wclss}}, \code{\link{wnlss}},
#' \code{\link{pnlss}}, \code{\link[mgcv]{gam}}
#' @export
cardlss <- function(link = list("tanhalf", "logithalf")) {
  if (length(link) != 2) stop("cardlss requires 2 links specified as character strings")
  if (!(link[[1]] %in% "tanhalf"))
    stop(link[[1]], " link not available for the location parameter of cardlss")
  if (!(link[[2]] %in% "logithalf"))
    stop(link[[2]], " link not available for the concentration parameter of cardlss")

  stats <- list()
  stats[[1]] <- tanhalf.link()
  stats[[2]] <- logithalf.link()

  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    type <- match.arg(type)
    mu <- object$fitted[, 1]
    rho <- object$fitted[, 2]
    d <- object$y - mu
    if (type == "response") return(atan2(sin(d), cos(d)))
    sind <- sin(d)
    if (type == "pearson") {
      ## the cardioid is first-harmonic only, so its second trig moment is 0
      ## and Var(sin(y - mu)) = 1/2 exactly, for every rho
      return(sind / sqrt(0.5))
    }
    ## deviance: 2*(peak log-density at d = 0 minus attained), signed by the
    ## sine of the residual angle (the -log 2pi cancels in the difference)
    sign(sind) * sqrt(pmax(2 * (log1p(2 * rho) - log1p(2 * rho * cos(d))), 0))
  }

  postproc <- expression({
    ## Deviance pair against the common saturated reference (per-obs peak
    ## log-density at the fitted rho), as in vmlss/wclss/wnlss; null = the
    ## moment fit (mean direction, rho = Rbar -- the cardioid moment estimator,
    ## since the first trig moment is rho, clipped below 1/2). The cardioid
    ## density is closed form, so unlike wnlss postproc needs no helper from
    ## the namespace. local() so mgcv's frames are safe.
    .cardlss.dev <- local({
      mu <- object$fitted[, 1]
      rho <- object$fitted[, 2]
      y <- object$y
      lpk <- log1p(2 * rho) - log(2 * pi)
      lat <- log1p(2 * rho * cos(y - mu)) - log(2 * pi)
      dev <- sum(2 * (lpk - lat))
      lsat <- sum(lpk)
      sy <- mean(sin(y)); cy <- mean(cos(y))
      mu0 <- atan2(sy, cy)
      rho0 <- min(max(sqrt(sy * sy + cy * cy), 1e-3), 0.5 - 1e-3)
      lnull <- sum(log1p(2 * rho0 * cos(y - mu0)) - log(2 * pi))
      c(dev, 2 * (lsat - lnull))
    })
    object$deviance <- .cardlss.dev[1]
    object$null.deviance <- .cardlss.dev[2]
    rm(.cardlss.dev)
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
    rho <- family$linfo[[2]]$linkinv(eta1)
    d <- y - mu
    s <- sin(d)
    cd <- cos(d)
    P <- 1 + 2 * rho * cd            # density numerator, 1 + 2*rho*cos(d)
    l0 <- log1p(2 * rho * cd) - log(2 * pi)
    l <- sum(wt * l0)

    if (deriv) {
      w <- 1 / P
      w2 <- w * w
      ## l1: d l / d(mu, rho); l2 columns ordered (mm, mr, rr)
      l1 <- cbind(2 * rho * s * w,
                  2 * cd * w)
      l2 <- cbind(-2 * rho * (cd * P + 2 * rho * s * s) * w2,
                  2 * s * w2,
                  -4 * cd * cd * w2)
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),
                   family$linfo[[2]]$mu.eta(eta1))
      g2 <- cbind(family$linfo[[1]]$d2link(mu),
                  family$linfo[[2]]$d2link(rho))
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
    ## Start from the moment fit: mean direction mu0 through the tan-half link,
    ## and rho0 = Rbar (the cardioid moment estimator, E cos(y - mu) = rho)
    ## through the logit-half link -- clipped strictly below 1/2, the family's
    ## upper bound. Constant targets per linear predictor fitted by penalized
    ## LS. Everything except the contract variables stays inside local().
    n <- rep(1, nobs)
    use.unscaled <- if (!is.null(attr(E, "use.unscaled"))) TRUE else FALSE
    if (is.null(start)) {
      start <- local({
        jj <- attr(x, "lpi")
        sy <- mean(sin(y))
        cy <- mean(cos(y))
        rho0 <- min(max(sqrt(sy * sy + cy * cy), 0.01), 0.49)
        st <- rep(0, ncol(x))
        ## location: data-following projected pilot
        off1 <- if (!is.null(offset) && length(offset) >= 1 &&
                    !is.null(offset[[1]])) offset[[1]] else 0
        st[jj[[1]]] <- getFromNamespace(".tanhalf_pilot_coef", "circlss")(
          x[, jj[[1]], drop = FALSE], E[, jj[[1]], drop = FALSE], y, off1,
          family$linfo[[1]]$linkfun)
        ## concentration: constant moment target through the logit-half link
        ytj <- rep(family$linfo[[2]]$linkfun(rho0), nobs)
        if (!is.null(offset) && length(offset) >= 2 && !is.null(offset[[2]]))
          ytj <- ytj - offset[[2]]
        x1 <- x[, jj[[2]], drop = FALSE]
        e1 <- E[, jj[[2]], drop = FALSE]
        stj <- qr.coef(qr(rbind(x1, e1)), c(ytj, rep(0, nrow(e1))))
        stj[!is.finite(stj)] <- 0
        st[jj[[2]]] <- stj
        st
      })
    }
  })

  rd <- function(mu, wt, scale) {
    ## exact: rejection sampling from a uniform envelope. The cardioid density
    ## is bounded by (1 + 2*rho)/(2*pi), so a draw phi ~ U(-pi, pi) is accepted
    ## with probability (1 + 2*rho*cos phi)/(1 + 2*rho) (mean acceptance
    ## 1/(1 + 2*rho) >= 1/2 since rho <= 1/2). theta = mu + phi.
    n <- nrow(mu)
    muv <- mu[, 1]
    rho <- pmin(pmax(mu[, 2], 0), 0.5)
    out <- numeric(n)
    todo <- seq_len(n)
    while (length(todo) > 0) {
      m <- length(todo)
      phi <- stats::runif(m, -pi, pi)
      u <- stats::runif(m)
      acc <- u <= (1 + 2 * rho[todo] * cos(phi)) / (1 + 2 * rho[todo])
      ok <- todo[acc]
      out[ok] <- muv[ok] + phi[acc]
      todo <- todo[!acc]
    }
    atan2(sin(out), cos(out))
  }

  structure(list(family = "cardlss", ll = ll, link = paste(link), nlp = 2,
                 param_names = c("mu", "rho"),
                 param_circular = c(TRUE, FALSE),
                 sandwich = sandwich, tri = mgcv::trind.generator(2),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 0, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
