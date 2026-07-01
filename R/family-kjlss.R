#' Kato-Jones Mobius four-parameter family
#'
#' A general family implementing distributional regression for a circular
#' (angular) response \eqn{y} in radians under the Kato-Jones (2015) law
#' \deqn{g(y) = \frac{1}{2\pi}\left[1 + \frac{2\gamma\,
#'   (\cos(y - \mu) - \rho\cos\lambda)}
#'   {1 + \rho^2 - 2\rho\cos(y - \mu - \lambda)}\right],}
#' a tractable four-parameter family obtained from a Mobius transformation of the
#' circle. The parameters control the first two trigonometric moments: \eqn{\mu}
#' the mean direction, \eqn{\gamma} the mean resultant length, and \eqn{(\rho,
#' \lambda)} the magnitude and phase of the second-order moment. It is the last and
#' most general family in the package, bringing both peakedness and skewness
#' through a single construction.
#'
#' It has \strong{four linear predictors}. To keep the second-order pair
#' \eqn{(\rho, \lambda)} inside its feasible region for \emph{every} coefficient
#' vector, the family is parameterized by the \strong{disc chart}: the Theorem-1
#' feasible set for the Cartesian shape pair \eqn{(a, b) = (\rho\cos\lambda,
#' \rho\sin\lambda)} is the disc of centre \eqn{(\gamma, 0)} and radius
#' \eqn{1-\gamma}, and the chart
#' \deqn{(a, b) = (\gamma, 0) + (1-\gamma)\,u/\sqrt{1 + \lVert u\rVert^2},
#'   \qquad u = (u_1, u_2) \in \mathbb{R}^2,}
#' maps an unconstrained \eqn{u} onto its interior. So the smooths ride the
#' unconstrained chart coordinates \eqn{u_1, u_2} (identity links) and the coupled
#' feasibility constraint can never be violated. Used with \code{mgcv::gam} and a
#' list of \emph{four} formulas: the first names the response and models \eqn{\mu},
#' then \eqn{\gamma}, then \eqn{u_1}, then \eqn{u_2}. The chart coordinates are most
#' often held global, i.e. fitted intercept-only with \code{~ 1}.
#'
#' @param link Four-element list of link names, for the location, the mean
#'   resultant length, and the two disc-chart coordinates. Currently only the
#'   defaults are available: \code{"tanhalf"} for the location \eqn{\mu},
#'   \code{"logit"} for the mean resultant length \eqn{\gamma \in (0, 1)}, and
#'   \code{"identity"} for the unconstrained chart coordinates
#'   \eqn{u_1, u_2 \in \mathbb{R}}.
#' @return An object of class \code{c("general.family", "extended.family", "family")}
#'   for use with \code{\link[mgcv]{gam}} (or its front end \code{\link{circ_gam}}).
#' @details
#' \strong{Exact normalizer.} Unlike every other shape family here, the Kato-Jones
#' normalizer is exactly \eqn{2\pi}: the density is a first-/second-trigonometric-
#' moment perturbation of the circular uniform that integrates to 1 with no special
#' function, so the log-density and all of its derivatives are elementary rational
#' functions of \eqn{(a, b)} and \eqn{\cos(y-\mu)}, \eqn{\sin(y-\mu)} -- there is no
#' quadrature, Bessel or Gamma term. The derivatives are taken with respect to the
#' chart coordinates by pushing the Cartesian scores through the chart Jacobian and
#' Hessian.
#'
#' \strong{Special and nested cases.} \eqn{u = 0} is exactly the wrapped Cauchy
#' \eqn{\mathrm{WC}(\mu, \gamma)} (\code{\link{wclss}}), the natural reduced model,
#' so intercept-only \eqn{u_1, u_2} is the wrapped Cauchy with covariate-free
#' second-order shape. \eqn{\gamma \to 0} gives the circular uniform and
#' \eqn{\rho \to 0} (\eqn{u} radial to 0) gives a cardioid (\code{\link{cardlss}}).
#'
#' \strong{Mean direction, not mode.} Once \eqn{\rho \neq 0} the density is
#' asymmetric and \eqn{\mu} is the direction of the first trigonometric moment, not
#' the mode; fitted-direction summaries inherit that reading.
#'
#' \strong{Optimizer.} Only first- and second-order log-likelihood derivatives
#' are provided, so \code{available.derivs = 0} and the family is fitted by the
#' extended Fellner-Schall optimizer rather than full Newton REML.
#' \code{\link[mgcv]{gam}} selects \code{optimizer = "efs"} automatically; passing
#' it explicitly is recommended so the fitted object is labelled correctly (which
#' \code{\link[mgcv]{gam.check}} relies on). The chart coordinates \eqn{u_1, u_2}
#' are weakly identified when the response is diffuse; holding them global
#' (\code{~ 1}) and keeping the mean direction and \eqn{\gamma} well identified is
#' the robust default. With four
#' linear predictors and two flat shape directions, a cyclic (\code{bs = "cc"})
#' model of all four can exceed \code{mgcv}'s fixed extended-Fellner-Schall
#' iteration cap, so prefer parametric or thin-plate location terms.
#'
#' The mean direction uses the Fisher-Lee tan-half link (\eqn{\mu \in (-\pi, \pi)},
#' antipode unrepresentable, winding number zero -- see \code{\link{pnlss}} when
#' the mean direction must wind).
#' @examples
#' library(mgcv)
#' set.seed(1)
#' n <- 300
#' x <- runif(n)
#' mu <- 2 * atan(0.8 * sin(2 * pi * x))
#' gamma <- plogis(0.8)         # mean resultant length ~ 0.69
#' u1 <- 0.4; u2 <- -0.5        # a fixed second-order shape
#' # Kato-Jones draws by inverse transform on the centered kernel (the chart maps
#' # (gamma, u1, u2) to the Cartesian shape pair (a, b), guaranteed feasible)
#' r <- sqrt(1 + u1^2 + u2^2); om <- 1 - gamma
#' a <- gamma + om * u1 / r; b <- om * u2 / r
#' phi <- seq(-pi, pi, length.out = 4001)
#' D <- 1 + a^2 + b^2 - 2 * a * cos(phi) - 2 * b * sin(phi)
#' g <- pmax(1 + 2 * gamma * (cos(phi) - a) / D, 0)
#' cdf <- cumsum(g); cdf <- cdf / cdf[length(cdf)]
#' dev <- approx(cdf, phi, runif(n), rule = 2)$y
#' y <- atan2(sin(mu + dev), cos(mu + dev))
#' # smooth location; global (intercept-only) gamma and chart coordinates
#' b <- gam(list(y ~ s(x), ~ 1, ~ 1, ~ 1), family = kjlss(), optimizer = "efs")
#' summary(b)
#' @references
#' Kato, S. and Jones, M. C. (2015) A tractable and interpretable four-parameter
#' family of unimodal distributions on the circle. \emph{Biometrika} 102, 181-190.
#'
#' Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method for
#' smoothing parameter optimization with application to Tweedie location,
#' scale and shape models. \emph{Biometrics} 73, 1071-1081.
#' @seealso
#' \code{\link{wclss}}, \code{\link{cardlss}}, \code{\link{jplss}},
#' \code{\link{ssjplss}}, \code{\link{vmlss}}, \code{\link{pnlss}},
#' \code{\link[mgcv]{gam}}
#' @export
kjlss <- function(link = list("tanhalf", "logit", "identity", "identity")) {
  if (length(link) != 4)
    stop("kjlss requires 4 links specified as character strings")
  if (!(link[[1]] %in% "tanhalf"))
    stop(link[[1]], " link not available for the location parameter of kjlss")
  if (!(link[[2]] %in% "logit"))
    stop(link[[2]], " link not available for the gamma parameter of kjlss")
  if (!(link[[3]] %in% "identity"))
    stop(link[[3]], " link not available for the u1 parameter of kjlss")
  if (!(link[[4]] %in% "identity"))
    stop(link[[4]], " link not available for the u2 parameter of kjlss")

  stats <- list()
  stats[[1]] <- tanhalf.link()
  ## gamma in (0, 1): the standard logit. mgcv's fix.family.link supplies the
  ## d2link/d3link/d4link the chain rule needs (verified == analytic), so no
  ## custom link object is required -- same idiom as the log link below it.
  stats[[2]] <- stats::make.link("logit")
  fam <- structure(list(link = "logit", canonical = "none",
                        linkfun = stats[[2]]$linkfun,
                        mu.eta = stats[[2]]$mu.eta),
                   class = "family")
  fam <- mgcv::fix.family.link(fam)
  stats[[2]]$d2link <- fam$d2link
  stats[[2]]$d3link <- fam$d3link
  stats[[2]]$d4link <- fam$d4link
  for (k in 3:4) {
    stats[[k]] <- stats::make.link("identity")
    stats[[k]]$d2link <- function(mu) numeric(length(mu))
    stats[[k]]$d3link <- function(mu) numeric(length(mu))
    stats[[k]]$d4link <- function(mu) numeric(length(mu))
  }

  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    type <- match.arg(type)
    mu <- object$fitted[, 1]
    gamma <- object$fitted[, 2]
    u1 <- object$fitted[, 3]
    u2 <- object$fitted[, 4]
    d <- object$y - mu
    if (type == "response") return(atan2(sin(d), cos(d)))
    sind <- sin(d)
    if (type == "pearson") {
      ## diagnostic standardization by the first-moment scale 1 - gamma^2
      v <- (1 - gamma * gamma) / 2
      v[v <= 0] <- .Machine$double.eps
      return(sind / sqrt(v))
    }
    ## deviance: 2 (l_sat - l_hat), both proper log-densities (normalizer 2pi)
    lhat <- .kj_terms(object$y, mu, gamma, u1, u2, deriv = FALSE)$l0
    lsat <- .kj_lsat(mu, gamma, u1, u2)
    sign(sind) * sqrt(pmax(2 * (lsat - lhat), 0))
  }

  postproc <- expression({
    ## Deviance against the per-observation saturated reference, and a null
    ## deviance against the constant von Mises MLE (the same intercept-only
    ## reference the other shape families use for "deviance explained"). All
    ## terms are proper log-densities, so the normalizers are handled. local()
    ## keeps mgcv's frames intact.
    .kjlss.dev <- local({
      kt <- getFromNamespace(".kj_terms", "circlss")
      lsatf <- getFromNamespace(".kj_lsat", "circlss")
      mu <- object$fitted[, 1]
      gamma <- object$fitted[, 2]
      u1 <- object$fitted[, 3]
      u2 <- object$fitted[, 4]
      y <- object$y
      lhat <- kt(y, mu, gamma, u1, u2, deriv = FALSE)$l0
      lsat <- lsatf(mu, gamma, u1, u2)
      dev <- sum(2 * pmax(lsat - lhat, 0))
      sy <- mean(sin(y)); cy <- mean(cos(y))
      mu0 <- atan2(sy, cy)
      Rbar <- min(sqrt(sy * sy + cy * cy), 1 - 1e-08)
      kappa0 <- if (Rbar < 0.53) 2 * Rbar + Rbar^3 + 5 * Rbar^5 / 6 else
        if (Rbar < 0.85) -0.4 + 1.39 * Rbar + 0.43 / (1 - Rbar) else
          1 / (Rbar^3 - 4 * Rbar^2 + 3 * Rbar)
      kappa0 <- min(max(kappa0, 0.01), 500)
      lnull <- sum(kappa0 * (cos(y - mu0) - 1) -
                     log(2 * pi * besselI(kappa0, 0, expon.scaled = TRUE)))
      c(dev, 2 * (sum(lsat) - lnull))
    })
    object$deviance <- .kjlss.dev[1]
    object$null.deviance <- .kjlss.dev[2]
    rm(.kjlss.dev)
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
    mu <- family$linfo[[1]]$linkinv(eta0)
    gamma <- family$linfo[[2]]$linkinv(eta1)
    u1 <- family$linfo[[3]]$linkinv(eta2)
    u2 <- family$linfo[[4]]$linkinv(eta3)

    tt <- .kj_terms(y, mu, gamma, u1, u2, deriv = deriv > 0)
    l0 <- tt$l0
    l <- sum(wt * l0)
    ## size-aware MAP degeneracy penalty (circ_mix M-step only; inert otherwise):
    ## a ridge on the disc-chart coordinates (u1, u2) toward 0 = wrapped Cauchy,
    ## keeping the shape off the feasibility-disc boundary where the kernel
    ## denominator D -> 0 and the Hessian (~ 1/(1-rho)^6) blows up.
    pen <- if (.degen_active(family))
      .lss_map_penalty(family, cbind(mu, gamma, u1, u2), family$map_lambda) else NULL
    if (!is.null(pen)) l <- l + sum(wt * pen$l0)

    if (deriv) {
      l1 <- tt$l1            # n x 4 : mu, gamma, u1, u2
      l2 <- tt$l2            # n x 10 : combinations_with_replacement order
      if (!is.null(pen)) { l1 <- l1 + pen$l1; l2 <- l2 + pen$l2 }
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta0),
                   family$linfo[[2]]$mu.eta(eta1),
                   family$linfo[[3]]$mu.eta(eta2),
                   family$linfo[[4]]$mu.eta(eta3))
      g2 <- cbind(family$linfo[[1]]$d2link(mu),
                  family$linfo[[2]]$d2link(gamma),
                  family$linfo[[3]]$d2link(u1),
                  family$linfo[[4]]$d2link(u2))
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
    ## Location starts from the data-following projected pilot
    ## (.tanhalf_pilot_coef). gamma starts from the moment R-bar through the
    ## logit; the chart coordinates u1, u2 start from the moment second-order
    ## shape pushed through the closed-form chart inverse -- all as constant
    ## targets. local() so mgcv's frames stay intact.
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
        mf <- getFromNamespace(".kj_moment_fit", "circlss")(y)
        ui <- getFromNamespace(".kj_disc_chart_inverse", "circlss")(
          mf[2], mf[3], mf[4])
        ## A covariate-driven mu makes the GLOBAL second-order moment pool
        ## different mean directions, so (rho, lam) can land on the feasibility
        ## circle and the chart inverse diverges (u -> +/-1e4), which makes
        ## mgcv's gam.fit5 report an indefinite penalized likelihood. Bounding
        ## |u| to a strong-but-finite shape keeps the start well-conditioned, and
        ## only bites when the moment fit is already wild (|u| large).
        ui <- pmax(pmin(ui, 8), -8)
        ## location: data-following projected pilot
        off1 <- if (!is.null(offset) && length(offset) >= 1 &&
                    !is.null(offset[[1]])) offset[[1]] else 0
        st[jj[[1]]] <- getFromNamespace(".tanhalf_pilot_coef", "circlss")(
          x[, jj[[1]], drop = FALSE], E[, jj[[1]], drop = FALSE], y, off1,
          family$linfo[[1]]$linkfun)
        ## gamma (logit), u1 / u2 (identity): constant targets through links
        g0 <- min(max(mf[2], 1e-04), 1 - 1e-04)
        targ <- c(NA, family$linfo[[2]]$linkfun(g0), ui[1], ui[2])
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
    ## Simulation by inverse transform on a fine centered grid of the kernel,
    ## one grid per unique (gamma, u1, u2); the kernel shape is mu-invariant so
    ## the draw is shifted by mu afterwards.
    n <- nrow(mu)
    muv <- mu[, 1]; gamma <- mu[, 2]; u1 <- mu[, 3]; u2 <- mu[, 4]
    grid <- seq(-pi, pi, length.out = 2049L)
    cg <- cos(grid); sg <- sin(grid)
    th <- numeric(n)
    key <- paste(gamma, u1, u2, sep = "\r")
    for (k in unique(key)) {
      idx <- which(key == k)
      ch <- .kj_chart(gamma[idx[1]], u1[idx[1]], u2[idx[1]])
      a <- ch$a; b <- ch$b; g <- gamma[idx[1]]
      D <- pmax(1 + a * a + b * b - 2 * a * cg - 2 * b * sg, 1e-300)
      dens <- pmax(1 + 2 * g * (cg - a) / D, 0)
      cdf <- cumsum((dens[-1] + dens[-length(dens)]) * diff(grid) / 2)
      cdf <- c(0, cdf); cdf <- cdf / cdf[length(cdf)]
      u <- stats::runif(length(idx))
      th[idx] <- muv[idx] +
        stats::approx(cdf, grid, xout = u, rule = 2, ties = "ordered")$y
    }
    atan2(sin(th), cos(th))
  }

  structure(list(family = "kjlss", ll = ll, link = paste(link), nlp = 4,
                 param_names = c("mu", "gamma", "u1", "u2"),
                 param_circular = c(TRUE, FALSE, FALSE, FALSE),
                 degen = list(.degen_ridge(3L), .degen_ridge(4L)),
                 sandwich = sandwich, tri = mgcv::trind.generator(4),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 0, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
