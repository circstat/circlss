#' Wrapped Cauchy location-scale family
#'
#' A general family implementing distributional regression for a circular
#' (angular) response \eqn{y} in radians under the wrapped Cauchy law
#' \deqn{f(y) = \frac{1 - \rho^2}{2\pi\,(1 + \rho^2 - 2\rho\cos(y - \mu))},}
#' with both the mean direction \eqn{\mu} and the mean resultant length
#' \eqn{\rho} getting their own linear predictor, each of which may contain
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
#' The wrapped Cauchy is the heavy-tailed counterpart of the von Mises
#' (\code{\link{vmlss}}): sharply peaked with fat circular tails, so it is
#' the more robust choice when the data contain angular outliers. Its
#' trigonometric moments are simply \eqn{\rho^p}, which gives clean
#' residual conventions: Pearson residuals standardize by
#' \eqn{\mathrm{Var}\{\sin(y-\mu)\} = (1-\rho^2)/2}.
#'
#' The mean direction uses the Fisher-Lee tan-half link
#' (\eqn{\mu \in (-\pi, \pi)}, antipode unrepresentable, winding number
#' zero -- see \code{\link{pnlss}} when the mean direction must wind).
#' Log-likelihood derivatives up to fourth order are implemented, so the family
#' supports full Newton REML (\code{method = "REML"}); \code{optimizer = "efs"}
#' also works.
#' Internally the density denominator is computed in the
#' cancellation-free form \eqn{(1-\rho)^2 + 4\rho\sin^2((y-\mu)/2)} so the
#' log-likelihood stays exact as \eqn{\rho \to 1}.
#' @examples
#' library(mgcv)
#' set.seed(1)
#' n <- 300
#' x <- runif(n)
#' mu <- 2 * atan(1.4 * sin(2 * pi * x))
#' rho <- plogis(0.5 + 0.8 * cos(2 * pi * x))
#' y <- atan2(sin(mu - log(rho) * rcauchy(n)), cos(mu - log(rho) * rcauchy(n)))
#' b <- gam(list(y ~ s(x), ~ s(x)), family = wclss(), method = "REML")
#' summary(b)
#' @references
#' Fisher, N. I. and Lee, A. J. (1992) Regression models for an angular
#' response. \emph{Biometrics} 48, 665-677.
#'
#' Wood, S. N., Pya, N. and Saefken, B. (2016) Smoothing parameter and model
#' selection for general smooth models. \emph{Journal of the American Statistical
#' Association} 111, 1548-1575.
#' @seealso
#' \code{\link{vmlss}}, \code{\link{pnlss}}, \code{\link[mgcv]{gam}}
#' @export
wclss <- function(link = list("tanhalf", "logit")) {
  if (length(link) != 2) stop("wclss requires 2 links specified as character strings")
  if (!(link[[1]] %in% "tanhalf"))
    stop(link[[1]], " link not available for the location parameter of wclss")
  if (!(link[[2]] %in% "logit"))
    stop(link[[2]], " link not available for the concentration parameter of wclss")

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
      ## WC trig moments are rho^p, so Var(sin(y - mu)) = (1 - rho^2)/2
      return(sind / sqrt((1 - rho * rho) / 2))
    }
    ## deviance: 2*(l_mode - l) = 2*log(D / (1 - rho)^2), signed
    D <- (1 - rho)^2 + 4 * rho * sin(0.5 * d)^2
    sign(sind) * sqrt(pmax(2 * (log(D) - 2 * log1p(-rho)), 0))
  }

  postproc <- expression({
    ## Deviance pair against the common saturated reference (per-obs
    ## log-lik at the fitted mode), as in vmlss; null = moment fit
    ## (mean direction, rho = Rbar -- the WC moment estimator). local()
    ## so mgcv's frames are safe.
    .wclss.dev <- local({
      mu <- object$fitted[, 1]
      rho <- object$fitted[, 2]
      y <- object$y
      D <- (1 - rho)^2 + 4 * rho * sin(0.5 * (y - mu))^2
      dev <- sum(2 * (log(D) - 2 * log1p(-rho)))
      lsat <- sum(log1p(rho) - log1p(-rho) - log(2 * pi))
      sy <- mean(sin(y)); cy <- mean(cos(y))
      mu0 <- atan2(sy, cy)
      rho0 <- min(max(sqrt(sy * sy + cy * cy), 1e-3), 1 - 1e-3)
      D0 <- (1 - rho0)^2 + 4 * rho0 * sin(0.5 * (y - mu0))^2
      lnull <- sum(log1p(-rho0) + log1p(rho0) - log(2 * pi) - log(D0))
      c(dev, 2 * (lsat - lnull))
    })
    object$deviance <- .wclss.dev[1]
    object$null.deviance <- .wclss.dev[2]
    rm(.wclss.dev)
  })

  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                 D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
    ## deriv: 0 - eval; 1 - grad and Hess; 2 - diagonal of first deriv of
    ## Hess; 3 - first deriv of Hess; 4 - everything (gaulss convention).
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
    Dq <- (1 - rho)^2 + 4 * rho * sin(0.5 * d)^2
    l0 <- log1p(-rho) + log1p(rho) - log(2 * pi) - log(Dq)
    l <- sum(wt * l0)
    ## size-aware MAP degeneracy penalty (circ_mix M-step only; inert otherwise).
    pen <- if (.degen_active(family))
      .lss_map_penalty(family, cbind(mu, rho), family$map_lambda) else NULL
    if (!is.null(pen)) l <- l + sum(wt * pen$l0)

    if (deriv) {
      w <- 1 / Dq
      w2 <- w * w
      om <- 1 - rho * rho            # 1 - rho^2
      Dm <- -2 * rho * s
      Dr <- 2 * rho - 2 * cd
      ## l1: d l / d(mu, rho); l2 columns ordered (mm, mr, rr)
      l1 <- cbind(2 * rho * s * w,
                  -2 * rho / om - Dr * w)
      l2 <- cbind(-2 * rho * cd * w + Dm * Dm * w2,
                  2 * s * w + Dm * Dr * w2,
                  -2 * (1 + rho * rho) / om^2 - 2 * w + Dr * Dr * w2)
      if (!is.null(pen)) { l1 <- l1 + pen$l1; l2 <- l2 + pen$l2 }
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),
                   family$linfo[[2]]$mu.eta(eta1))
      g2 <- cbind(family$linfo[[1]]$d2link(mu),
                  family$linfo[[2]]$d2link(rho))
    }
    l3 <- l4 <- g3 <- g4 <- 0
    if (deriv > 1) {
      ## l3 columns (mmm, mmr, mrr, rrr)
      w3 <- w2 * w
      Dmm <- 2 * rho * cd
      Dmr <- -2 * s
      Drr <- 2
      Dmmm <- 2 * rho * s
      Dmmr <- 2 * cd
      l3 <- cbind(
        -Dmmm * w + 3 * Dmm * Dm * w2 - 2 * Dm^3 * w3,
        -Dmmr * w + (Dmm * Dr + 2 * Dmr * Dm) * w2 - 2 * Dm * Dm * Dr * w3,
        (2 * Dmr * Dr + Drr * Dm) * w2 - 2 * Dm * Dr * Dr * w3,
        -4 * rho * (3 + rho * rho) / om^3 + 3 * Drr * Dr * w2 -
          2 * Dr^3 * w3
      )
      g3 <- cbind(family$linfo[[1]]$d3link(mu),
                  family$linfo[[2]]$d3link(rho))
    }
    if (deriv > 3) {
      ## l4 columns (mmmm, mmmr, mmrr, mrrr, rrrr)
      w4 <- w3 * w
      Dmmmm <- -2 * rho * cd
      Dmmmr <- 2 * s
      l4 <- cbind(
        -Dmmmm * w + (4 * Dmmm * Dm + 3 * Dmm * Dmm) * w2 -
          12 * Dmm * Dm * Dm * w3 + 6 * Dm^4 * w4,
        -Dmmmr * w + (Dmmm * Dr + 3 * Dmmr * Dm + 3 * Dmm * Dmr) * w2 -
          6 * (Dmm * Dm * Dr + Dmr * Dm * Dm) * w3 + 6 * Dm^3 * Dr * w4,
        (2 * Dmmr * Dr + Dmm * Drr + 2 * Dmr * Dmr) * w2 -
          2 * (Dmm * Dr * Dr + Drr * Dm * Dm + 4 * Dmr * Dm * Dr) * w3 +
          6 * Dm * Dm * Dr * Dr * w4,
        3 * Dmr * Drr * w2 - 6 * (Drr * Dm * Dr + Dmr * Dr * Dr) * w3 +
          6 * Dm * Dr^3 * w4,
        -12 * (1 + 6 * rho * rho + rho^4) / om^4 + 3 * Drr * Drr * w2 -
          12 * Drr * Dr * Dr * w3 + 6 * Dr^4 * w4
      )
      g4 <- cbind(family$linfo[[1]]$d4link(mu),
                  family$linfo[[2]]$d4link(rho))
    }
    if (deriv) {
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
    ## link, and rho0 = Rbar (the WC moment estimator, E cos(y - mu) =
    ## rho) through the logit link; constant targets per linear predictor
    ## fitted by penalized LS. Everything except the contract variables
    ## stays inside local().
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
    ## exact: the wrapped Cauchy is the wrapping of a Cauchy with scale
    ## gamma = -log(rho)
    n <- nrow(mu)
    gam_ <- -log(pmin(pmax(mu[, 2], 1e-12), 1 - 1e-12))
    th <- mu[, 1] + gam_ * stats::rcauchy(n)
    atan2(sin(th), cos(th))
  }

  structure(list(family = "wclss", ll = ll, link = paste(link), nlp = 2,
                 param_names = c("mu", "rho"),
                 param_circular = c(TRUE, FALSE),
                 degen = list(.degen_linear(2L)),
                 sandwich = sandwich, tri = mgcv::trind.generator(2),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 2, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
