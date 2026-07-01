#' Jones-Pewsey location-concentration-shape family
#'
#' A general family implementing distributional regression for a circular
#' (angular) response \eqn{y} in radians under the Jones-Pewsey law
#' \deqn{f(y) = c(\kappa, \psi)\,\bigl(\cosh(\kappa\psi) +
#'   \sinh(\kappa\psi)\cos(y - \mu)\bigr)^{1/\psi},
#'   \qquad \kappa > 0,\ \psi \in \mathbb{R},}
#' with \eqn{c(\kappa,\psi)} a normalizing constant. This is the symmetric
#' \strong{peakedness umbrella} that nests several of the package's other
#' families through the single shape parameter \eqn{\psi}.
#'
#' It is the first family with \strong{three linear predictors}: the mean
#' direction \eqn{\mu}, the concentration \eqn{\kappa}, and the shape \eqn{\psi}
#' each get their own, any of which may contain smooth terms. Used with
#' \code{mgcv::gam} and a list of \emph{three} formulas; the first names the
#' response and models \eqn{\mu}, the second models \eqn{\log\kappa}, and the
#' third models \eqn{\psi}. The shape is most often held global, i.e. fitted
#' intercept-only with \code{~ 1}.
#'
#' @param link Three-element list of link names, for the mean direction, the
#'   concentration, and the shape. Currently only the defaults are available:
#'   \code{"tanhalf"} for the location \eqn{\mu}, \code{"log"} for the
#'   concentration \eqn{\kappa > 0}, and \code{"identity"} for the shape
#'   \eqn{\psi \in \mathbb{R}}.
#' @return An object of class \code{c("general.family", "extended.family", "family")}
#'   for use with \code{\link[mgcv]{gam}} (or its front end \code{\link{circ_gam}}).
#' @details
#' The shape \eqn{\psi} indexes a family of symmetric circular laws and recovers
#' several special cases exactly:
#' \itemize{
#'   \item \eqn{\psi \to 0}: the von Mises (\code{\link{vmlss}}),
#'     \eqn{\propto \exp(\kappa\cos(y-\mu))};
#'   \item \eqn{\psi = 1}: the cardioid (\code{\link{cardlss}}) with mean
#'     resultant length \eqn{\rho = \tfrac12\tanh\kappa};
#'   \item \eqn{\psi = -1}: the wrapped Cauchy (\code{\link{wclss}});
#'   \item \eqn{\kappa \to 0}: the circular uniform.
#' }
#' The Cartwright family (\code{\link{cartlss}}) sits on this family's
#' concentrated boundary (\eqn{\kappa \to \infty} at \eqn{\psi = \zeta}), not in
#' its interior.
#'
#' \strong{Normalizer.} Unlike the earlier families, \eqn{c(\kappa, \psi)} has no
#' elementary form. The log-density value \emph{and} the \eqn{\kappa}- and
#' \eqn{\psi}-score and Hessian carry kernel-weighted moments of the normalizer,
#' evaluated by composite 24-point Gauss-Legendre quadrature on a feature-scale
#' break-point ladder (so the concentrated \eqn{\psi < 0} spike and the antipodal
#' near-kink are resolved at any concentration). This is the first circlss family
#' to integrate numerically inside the likelihood.
#'
#' \strong{Optimizer.} Only first- and second-order log-likelihood derivatives are
#' provided, so \code{available.derivs = 0} and the family is fitted by the extended
#' Fellner-Schall optimizer rather than full Newton REML. \code{\link[mgcv]{gam}}
#' selects \code{optimizer = "efs"} automatically; passing it explicitly is
#' recommended so the fitted object is labelled correctly (which
#' \code{\link[mgcv]{gam.check}} relies on). The shape \eqn{\psi} is weakly
#' identified when the response is diffuse or its concentration smooth is
#' over-flexible; holding \eqn{\psi} global (\code{~ 1}) and keeping the
#' concentration well identified is the robust default.
#'
#' The mean direction uses the Fisher-Lee tan-half link
#' (\eqn{\mu \in (-\pi, \pi)}, antipode unrepresentable, winding number
#' zero -- see \code{\link{pnlss}} when the mean direction must wind). The shape
#' rides the identity link and may be any real number.
#' @examples
#' library(mgcv)
#' set.seed(1)
#' n <- 300
#' x <- runif(n)
#' mu <- 2 * atan(sin(2 * pi * x))
#' kappa <- exp(1.0 + 0.5 * cos(2 * pi * x))
#' psi <- 0.6
#' # Jones-Pewsey draws by inverse transform on the centered kernel, then shift
#' phi <- seq(-pi, pi, length.out = 4001)
#' dev <- vapply(seq_len(n), function(i) {
#'   g <- (cosh(kappa[i] * psi) + sinh(kappa[i] * psi) * cos(phi))^(1 / psi)
#'   cdf <- cumsum(g); cdf <- cdf / cdf[length(cdf)]
#'   approx(cdf, phi, runif(1), rule = 2)$y
#' }, numeric(1))
#' y <- atan2(sin(mu + dev), cos(mu + dev))
#' # smooth mean direction and concentration; global (intercept-only) shape
#' b <- gam(list(y ~ s(x), ~ s(x), ~ 1), family = jplss(), optimizer = "efs")
#' summary(b)
#' @references
#' Jones, M. C. and Pewsey, A. (2005) A family of symmetric distributions on the
#' circle. \emph{Journal of the American Statistical Association} 100, 1422-1428.
#'
#' Pewsey, A., Neuhaeuser, M. and Ruxton, G. D. (2013) \emph{Circular Statistics
#' in R}. Oxford University Press.
#'
#' Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method for
#' smoothing parameter optimization with application to Tweedie location,
#' scale and shape models. \emph{Biometrics} 73, 1071-1081.
#' @seealso
#' \code{\link{vmlss}}, \code{\link{cardlss}}, \code{\link{wclss}},
#' \code{\link{cartlss}}, \code{\link{pnlss}}, \code{\link[mgcv]{gam}}
#' @export
jplss <- function(link = list("tanhalf", "log", "identity")) {
  if (length(link) != 3)
    stop("jplss requires 3 links specified as character strings")
  if (!(link[[1]] %in% "tanhalf"))
    stop(link[[1]], " link not available for the location parameter of jplss")
  if (!(link[[2]] %in% "log"))
    stop(link[[2]], " link not available for the concentration parameter of jplss")
  if (!(link[[3]] %in% "identity"))
    stop(link[[3]], " link not available for the shape parameter of jplss")

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
  ## shape psi rides the identity link; its 2nd-4th link derivatives are zero
  stats[[3]] <- stats::make.link("identity")
  stats[[3]]$d2link <- function(mu) numeric(length(mu))
  stats[[3]]$d3link <- function(mu) numeric(length(mu))
  stats[[3]]$d4link <- function(mu) numeric(length(mu))

  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    type <- match.arg(type)
    mu <- object$fitted[, 1]
    kappa <- object$fitted[, 2]
    psi <- object$fitted[, 3]
    d <- object$y - mu
    if (type == "response") return(atan2(sin(d), cos(d)))
    sind <- sin(d)
    if (type == "pearson") {
      ## Var{sin(y - mu)} = (1 - alpha2)/2, with the 2nd cosine moment alpha2
      ## from the quadrature ladder (no elementary form for Jones-Pewsey)
      a2 <- vapply(seq_along(kappa),
                   function(i) .jp_cos_moment(kappa[i], psi[i], 2), numeric(1))
      v <- (1 - a2) / 2
      v[v <= 0] <- .Machine$double.eps
      return(sind / sqrt(v))
    }
    ## deviance: 2*(l_sat - l_hat); the normalizer cancels (l_sat is the peak
    ## log-density at the fitted kappa, h(0) = kappa), leaving 2*(kappa - h)
    h <- .jp_score_terms(d, kappa, psi, second = FALSE)$h
    sign(sind) * sqrt(pmax(2 * (kappa - h), 0))
  }

  postproc <- expression({
    ## Deviance and null deviance against the common saturated reference (the
    ## per-obs peak log-density kappa + log c at the fitted kappa, psi). The
    ## deviance itself is normalizer-free, 2*(kappa - h); the null (constant
    ## von Mises MLE, the psi -> 0 baseline) needs the fitted normalizer, one
    ## quadrature pass via the internal helpers (fetched by name, as ll runs in
    ## mgcv's frame). All work in local() so mgcv's frames stay intact.
    .jplss.dev <- local({
      st <- getFromNamespace(".jp_score_terms", "circlss")
      lcv <- getFromNamespace(".jp_logc_vec", "circlss")
      mu <- object$fitted[, 1]
      kappa <- object$fitted[, 2]
      psi <- object$fitted[, 3]
      y <- object$y
      h <- st(y - mu, kappa, psi, second = FALSE)$h
      dev <- sum(2 * pmax(kappa - h, 0))
      lsat <- sum(kappa + lcv(kappa, psi))
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
    object$deviance <- .jplss.dev[1]
    object$null.deviance <- .jplss.dev[2]
    rm(.jplss.dev)
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
    psi <- family$linfo[[3]]$linkinv(eta2)
    d <- y - mu

    ts <- .jp_score_terms(d, kappa, psi, second = deriv > 0)
    Q <- .jp_quad_vec(kappa, psi)           # n x 6: logc, dk, dp, dkk, dkp, dpp
    ## l0 VALUE, with value shortcuts (uniform, von Mises)
    logc <- Q[, 1]
    vm <- (kappa >= .JP_KAPPA_TOL) & (abs(psi) < .JP_PSI_TOL)
    if (any(vm))
      logc[vm] <- -(log(2 * pi * besselI(kappa[vm], 0, expon.scaled = TRUE)) +
                      kappa[vm])
    h0 <- ts$h
    if (any(vm)) h0[vm] <- kappa[vm] * cos(d[vm])
    l0 <- h0 + logc
    unif <- kappa < .JP_KAPPA_TOL
    if (any(unif)) l0[unif] <- -log(2 * pi)
    l <- sum(wt * l0)
    ## size-aware MAP degeneracy penalty (circ_mix M-step only; inert otherwise):
    ## kappa toward 0 (linear) and the shape psi toward 0 = von Mises (ridge).
    pen <- if (.degen_active(family))
      .lss_map_penalty(family, cbind(mu, kappa, psi), family$map_lambda) else NULL
    if (!is.null(pen)) l <- l + sum(wt * pen$l0)

    if (deriv) {
      ## l1: d l / d(mu, kappa, psi); the derivatives use the quadrature
      ## moments at every kappa, psi (no value shortcut)
      l1 <- cbind(-ts$hphi, ts$hk - Q[, 2], ts$hp - Q[, 3])
      ## l2 columns ordered (mm, mk, mp, kk, kp, pp) = upper triangle, row-major
      l2 <- cbind(ts$hphiphi, -ts$hphik, -ts$hphip,
                  ts$hkk - Q[, 4], ts$hkp - Q[, 5], ts$hpp - Q[, 6])
      if (!is.null(pen)) { l1 <- l1 + pen$l1; l2 <- l2 + pen$l2 }
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta0),
                   family$linfo[[2]]$mu.eta(eta1),
                   family$linfo[[3]]$mu.eta(eta2))
      g2 <- cbind(family$linfo[[1]]$d2link(mu),
                  family$linfo[[2]]$d2link(kappa),
                  family$linfo[[3]]$d2link(psi))
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
    ## (.tanhalf_pilot_coef), which finds a well-conditioned basin even on
    ## near-flat shape surfaces. Concentration starts
    ## from Fisher's A1-inverse kappa0 and the shape from psi0 = 0 (the von
    ## Mises member), as constant targets through their links. local() so mgcv's
    ## frames (gam.fit5, initial.spg) are never clobbered.
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
        ## concentration (log) and shape (identity): constant targets
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
    ## one grid per unique (kappa, psi). Adequate for predictive draws at
    ## fitted parameters (the reference samples through a feature-scale-graded
    ## table, needed only for very narrow psi < 0 spikes, not here).
    n <- nrow(mu)
    muv <- mu[, 1]
    kappa <- mu[, 2]
    psi <- mu[, 3]
    grid <- seq(-pi, pi, length.out = 2049L)
    th <- numeric(n)
    key <- paste(kappa, psi, sep = "\r")
    for (k in unique(key)) {
      idx <- which(key == k)
      ka <- kappa[idx[1]]
      ps <- psi[idx[1]]
      h <- .jp_score_terms(grid, ka, ps, second = FALSE)$h
      dens <- exp(h - max(h))
      cdf <- cumsum((dens[-1] + dens[-length(dens)]) * diff(grid) / 2)
      cdf <- c(0, cdf)
      cdf <- cdf / cdf[length(cdf)]
      u <- stats::runif(length(idx))
      phi <- stats::approx(cdf, grid, xout = u, rule = 2, ties = "ordered")$y
      th[idx] <- muv[idx] + phi
    }
    atan2(sin(th), cos(th))
  }

  structure(list(family = "jplss", ll = ll, link = paste(link), nlp = 3,
                 param_names = c("mu", "kappa", "psi"),
                 param_circular = c(TRUE, FALSE, FALSE),
                 ## a firmer prior (scale > 1) on the Jones-Pewsey shape psi keeps
                 ## it in the range where the normalizer quadrature stays fast and
                 ## non-singular, so c = 1 holds as the default. kappa stays at the
                 ## default scale: a stronger pull would zero out a diffuse
                 ## component (kappa -> 0, where the quadrature is slowest), and the
                 ## size-aware lambda already caps the concentration. Standalone fits
                 ## are unaffected (the penalty is inert without circ_mix's map_lambda).
                 degen = list(.degen_linear(2L), .degen_ridge(3L, 30)),
                 sandwich = sandwich, tri = mgcv::trind.generator(3),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 0, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
