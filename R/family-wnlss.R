.WN_FOURIER_RHO_MAX <- 0.8

## Per-observation wrapped-normal log-density l0 and, when deriv >= 1, its
## first/second derivatives w.r.t. the distribution parameters mu and rho.
## `d` is the centred angle phi = y - mu; `lmu`/`lrho` are d l / d mu and
## d l / d rho (the phi -> mu sign flip is folded into the trig terms, as in
## the reference). l2 entries: lmm = d2/dmu2, lmr = d2/dmu drho, lrr =
## d2/drho2. Returns length-n vectors (NULL derivative slots when deriv = 0).
.wn_terms <- function(d, rho, deriv = 0) {
  n <- length(d)
  l0 <- numeric(n)
  lmu <- lrho <- lmm <- lmr <- lrr <- if (deriv) numeric(n) else NULL

  lo <- rho <= .WN_FOURIER_RHO_MAX

  if (any(lo)) {
    ii <- which(lo)
    dd <- d[ii]; rr <- rho[ii]
    p <- 1:29
    P2 <- p * p
    rp <- outer(rr, P2, `^`)              # rho^(p^2)
    cs <- cos(outer(dd, p))               # cos(p*phi)
    sn <- sin(outer(dd, p))               # sin(p*phi)
    seriesc <- rowSums(rp * cs)
    S0 <- 1 + 2 * seriesc
    l0[ii] <- log1p(2 * seriesc) - log(2 * pi)
    if (deriv) {
      rpm1 <- outer(rr, P2 - 1, `^`)      # rho^(p^2 - 1)
      Smu <- 2 * drop((rp * sn) %*% p)
      Srho <- 2 * drop((rpm1 * cs) %*% P2)
      lmu_lo <- Smu / S0
      lrho_lo <- Srho / S0
      lmu[ii] <- lmu_lo
      lrho[ii] <- lrho_lo
      Smm <- -2 * drop((rp * cs) %*% P2)
      Smr <- 2 * drop((rpm1 * sn) %*% (p^3))
      ## S_rho rho starts at p = 2: the p = 1 coefficient p^2(p^2-1) is 0,
      ## so rho^(p^2-2) never sees a negative exponent
      p2 <- p[-1]; P2b <- p2 * p2
      rpm2 <- outer(rr, P2b - 2, `^`)     # rho^(p^2 - 2), p >= 2
      Srr <- 2 * drop((rpm2 * cs[, -1, drop = FALSE]) %*% (P2b * (P2b - 1)))
      lmm[ii] <- Smm / S0 - lmu_lo * lmu_lo
      lmr[ii] <- Smr / S0 - lmu_lo * lrho_lo
      lrr[ii] <- Srr / S0 - lrho_lo * lrho_lo
    }
  }

  hi <- !lo
  if (any(hi)) {
    jj <- which(hi)
    dd <- d[jj]
    rh <- pmin(rho[jj], 1 - 1e-15)
    sigma <- sqrt(-2 * log(rh))           # wrapped Gaussian sd
    k <- -2:2
    Z <- outer(dd, 2 * pi * k, `+`) / sigma   # (phi + 2pi k)/sigma, row-wise
    Z2 <- Z * Z
    minz2 <- apply(Z2, 1, min)
    Wun <- exp(-0.5 * (Z2 - minz2))       # max-stabilized image weights
    Wsum <- rowSums(Wun)
    l0[jj] <- (-0.5 * minz2 + log(Wsum)) - log(sigma) - 0.5 * log(2 * pi)
    if (deriv) {
      W <- Wun / Wsum
      m1 <- rowSums(W * Z)
      m2 <- rowSums(W * Z2)
      l_sig <- (m2 - 1) / sigma           # d l / d sigma
      sp <- -1 / (rh * sigma)             # d sigma / d rho
      lmu[jj] <- m1 / sigma
      lrho[jj] <- l_sig * sp
      m3 <- rowSums(W * Z2 * Z)
      m4 <- rowSums(W * Z2 * Z2)
      s2 <- sigma * sigma
      l_mm <- (m2 - 1 - m1 * m1) / s2
      l_ms <- (m3 - 3 * m1 - m1 * (m2 - 1)) / s2
      l_ss <- (m4 - 5 * m2 + 2 - (m2 - 1)^2) / s2
      spp <- (s2 - 1) / (rh * rh * sigma * s2)   # d2 sigma / d rho2
      lmm[jj] <- l_mm
      lmr[jj] <- l_ms * sp
      lrr[jj] <- l_ss * sp * sp + l_sig * spp
    }
  }
  list(l0 = l0, lmu = lmu, lrho = lrho, lmm = lmm, lmr = lmr, lrr = lrr)
}

#' Wrapped normal location-scale family
#'
#' A general family implementing distributional regression for a circular
#' (angular) response \eqn{y} in radians under the wrapped normal law
#' \deqn{f(y) = \frac{1}{2\pi}\left(1 + 2\sum_{p=1}^{\infty}
#'   \rho^{p^2}\cos\{p(y - \mu)\}\right),}
#' the wrapping of \eqn{N(\mu, \sigma^2)} with \eqn{\rho = e^{-\sigma^2/2}}.
#' Both the mean direction \eqn{\mu} and the mean resultant length
#' \eqn{\rho} get their own linear predictor, each of which may contain
#' smooth terms. Used with \code{mgcv::gam} and a list of two formulas: the
#' first names the response and models \eqn{\mu}, the second models
#' \eqn{\mathrm{logit}(\rho)}.
#'
#' @param link Two-element list of link names, for the mean direction and
#'   the mean resultant length. Currently only the defaults are available:
#'   \code{"tanhalf"} for the location and \code{"logit"} for the
#'   concentration parameter \eqn{\rho \in (0, 1)}.
#' @return An object of class \code{c("general.family", "extended.family", "family")}
#'   for use with \code{\link[mgcv]{gam}} (or its front end \code{\link{circ_gam}}).
#' @details
#' The wrapped normal is the bell-shaped circular law obtained by wrapping a
#' Gaussian onto the circle -- close to the von Mises (\code{\link{vmlss}})
#' in shape but defined through its mean resultant length \eqn{\rho} rather
#' than a concentration. Its trigonometric moments are \eqn{\rho^{p^2}}, so
#' Pearson residuals standardize by
#' \eqn{\mathrm{Var}\{\sin(y-\mu)\} = (1-\rho^4)/2}.
#'
#' Unlike the von Mises and wrapped Cauchy, the wrapped normal has no
#' closed-form normalizer. The log-density and its derivatives are evaluated
#' by a hybrid that switches per observation at \eqn{\rho = 0.8}: the Fourier
#' series above for \eqn{\rho \le 0.8} (29 terms; truncation
#' \eqn{\le 2\times10^{-82}}), and a log-sum-exp over wrapped Gaussian images
#' for \eqn{\rho > 0.8}, where the Fourier partial sums lose accuracy in the
#' tails.
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
#' n <- 300
#' x <- runif(n)
#' mu <- 2 * atan(1.2 * sin(2 * pi * x))
#' rho <- plogis(0.8 + 1.0 * cos(2 * pi * x))
#' sigma <- sqrt(-2 * log(rho))
#' y <- atan2(sin(mu + sigma * rnorm(n)), cos(mu + sigma * rnorm(n)))
#' b <- gam(list(y ~ s(x), ~ s(x)), family = wnlss(), optimizer = "efs")
#' summary(b)
#' @references
#' Pewsey, A., Neuhaeuser, M. and Ruxton, G. D. (2013) \emph{Circular
#' Statistics in R}. Oxford University Press.
#'
#' Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method for
#' smoothing parameter optimization with application to Tweedie location,
#' scale and shape models. \emph{Biometrics} 73, 1071-1081.
#' @seealso
#' \code{\link{vmlss}}, \code{\link{wclss}}, \code{\link{pnlss}},
#' \code{\link[mgcv]{gam}}
#' @export
wnlss <- function(link = list("tanhalf", "logit")) {
  if (length(link) != 2) stop("wnlss requires 2 links specified as character strings")
  if (!(link[[1]] %in% "tanhalf"))
    stop(link[[1]], " link not available for the location parameter of wnlss")
  if (!(link[[2]] %in% "logit"))
    stop(link[[2]], " link not available for the concentration parameter of wnlss")

  stats <- list()
  stats[[1]] <- tanhalf.link()
  stats[[2]] <- stats::make.link("logit")
  fam <- structure(list(link = "logit", canonical = "none",
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
    rho <- object$fitted[, 2]
    d <- object$y - mu
    if (type == "response") return(atan2(sin(d), cos(d)))
    sind <- sin(d)
    if (type == "pearson") {
      ## WN trig moments are rho^(p^2), so E cos 2(y - mu) = rho^4 and
      ## Var(sin(y - mu)) = (1 - rho^4)/2
      return(sind / sqrt((1 - rho^4) / 2))
    }
    ## deviance: per-obs 2*(peak log-density at the fitted rho minus
    ## attained), signed by the sine of the residual angle
    lpk <- .wn_terms(numeric(length(d)), rho)$l0
    lat <- .wn_terms(d, rho)$l0
    sign(sind) * sqrt(pmax(2 * (lpk - lat), 0))
  }

  postproc <- expression({
    ## Deviance pair against the common saturated reference (per-obs peak
    ## log-density at the fitted rho), as in vmlss/wclss; null = the moment
    ## fit (mean direction, rho = Rbar -- the WN moment estimator, since the
    ## first trig moment is rho). local() so mgcv's frames are safe.
    .wnlss.dev <- local({
      ## postproc is eval()'d in mgcv's frame, off circlss's search path, so
      ## the unexported helper is fetched explicitly (closures like ll /
      ## residuals capture it lexically and need no such fetch)
      wn_terms <- getFromNamespace(".wn_terms", "circlss")
      mu <- object$fitted[, 1]
      rho <- object$fitted[, 2]
      y <- object$y
      n <- length(y)
      lpk <- wn_terms(numeric(n), rho)$l0
      lat <- wn_terms(y - mu, rho)$l0
      dev <- sum(2 * (lpk - lat))
      lsat <- sum(lpk)
      sy <- mean(sin(y)); cy <- mean(cos(y))
      mu0 <- atan2(sy, cy)
      rho0 <- min(max(sqrt(sy * sy + cy * cy), 1e-3), 1 - 1e-3)
      lnull <- sum(wn_terms(y - mu0, rep(rho0, n))$l0)
      c(dev, 2 * (lsat - lnull))
    })
    object$deviance <- .wnlss.dev[1]
    object$null.deviance <- .wnlss.dev[2]
    rm(.wnlss.dev)
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

    tw <- .wn_terms(d, rho, deriv = if (deriv) 1 else 0)
    l0 <- tw$l0
    l <- sum(wt * l0)

    if (deriv) {
      ## l1: d l / d(mu, rho); l2 columns ordered (mm, mr, rr)
      l1 <- cbind(tw$lmu, tw$lrho)
      l2 <- cbind(tw$lmm, tw$lmr, tw$lrr)
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
    ## Start from the moment fit: mean direction mu0 through the tan-half
    ## link, and rho0 = Rbar (the WN moment estimator, E cos(y - mu) = rho)
    ## through the logit link; constant targets per linear predictor fitted
    ## by penalized LS. Everything except the contract variables stays inside
    ## local().
    n <- rep(1, nobs)
    use.unscaled <- if (!is.null(attr(E, "use.unscaled"))) TRUE else FALSE
    if (is.null(start)) {
      start <- local({
        jj <- attr(x, "lpi")
        sy <- mean(sin(y))
        cy <- mean(cos(y))
        rho0 <- min(max(sqrt(sy * sy + cy * cy), 0.01), 0.95)
        st <- rep(0, ncol(x))
        ## location: data-following projected pilot
        off1 <- if (!is.null(offset) && length(offset) >= 1 &&
                    !is.null(offset[[1]])) offset[[1]] else 0
        st[jj[[1]]] <- getFromNamespace(".tanhalf_pilot_coef", "circlss")(
          x[, jj[[1]], drop = FALSE], E[, jj[[1]], drop = FALSE], y, off1,
          family$linfo[[1]]$linkfun)
        ## concentration: constant moment target through its link
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
    ## exact: the wrapped normal is the wrapping of N(mu, sigma^2) with
    ## sigma = sqrt(-2 log rho)
    n <- nrow(mu)
    rho <- pmin(pmax(mu[, 2], 1e-12), 1 - 1e-12)
    sigma <- sqrt(-2 * log(rho))
    th <- mu[, 1] + sigma * stats::rnorm(n)
    atan2(sin(th), cos(th))
  }

  structure(list(family = "wnlss", ll = ll, link = paste(link), nlp = 2,
                 param_names = c("mu", "rho"),
                 param_circular = c(TRUE, FALSE),
                 sandwich = sandwich, tri = mgcv::trind.generator(2),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 0, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
