#' Weight-aware Gaussian location-scale family
#'
#' A Gaussian location-scale family for distributional regression, modelling a
#' real-valued response \eqn{y} with
#' \deqn{y \sim N(\mu, \sigma^2),}
#' where the mean \eqn{\mu} and the precision \eqn{\tau = 1/\sigma} each get their
#' own linear predictor, which may contain smooth terms. It is a weight-aware,
#' metadata-carrying adaptation of \pkg{mgcv}'s \code{\link[mgcv]{gaulss}}: unlike
#' \code{gaulss} it honours prior \code{weights} (needed for a weighted MLE and for
#' EM mixtures), and it carries the circlss parameter metadata so
#' \code{\link{circ_gam}} treats it as a first-class location-scale family with
#' named, response-scale output.
#'
#' @param link Two-element list of link names for the mean and the precision,
#'   following \code{\link[mgcv]{gaulss}}: \code{"identity"} (or \code{"log"},
#'   \code{"inverse"}, \code{"sqrt"}) for the mean and \code{"logb"} for the
#'   precision.
#' @param b The \code{logb} link's offset, as in \code{\link[mgcv]{gaulss}}.
#' @return An object of class \code{c("general.family", "extended.family", "family")}
#'   for use with \code{\link[mgcv]{gam}} (or its front end \code{\link{circ_gam}}).
#' @details
#' In the circlss regression trio this is the \strong{linear-circular} (l~c)
#' family: a real-valued response over a \emph{circular} covariate -- a level that
#' varies around a cycle (time of day, season, phase) -- fitted with a cyclic
#' smooth,
#' \preformatted{circ_gam(list(y ~ s(phi, bs = "cc"), ~ s(phi, bs = "cc")),
#'          family = gausslss())}
#' \code{\link{circ_gam}} then places the fitted mean on the "can" (an upright
#' cylinder: the circular covariate wraps the ring, the linear response is the
#' height).
#'
#' The parameterization follows \code{gaulss} exactly: the mean uses an identity
#' link and the second parameter is the \emph{precision} \eqn{\tau = 1/\sigma} on
#' the \code{logb} link, so the second fitted column is \eqn{1/\sigma}, not the
#' standard deviation. Log-likelihood derivatives up to fourth order are
#' implemented, so the family supports full Newton REML (\code{method = "REML"});
#' \code{optimizer = "efs"} also works. At unit weights the fit matches
#' \code{gaulss}; integer prior weights reproduce a row-replicated fit.
#'
#' This family adapts GPL-licensed code from \pkg{mgcv}; see the package's
#' \code{inst/COPYRIGHTS}.
#' @examples
#' library(mgcv)
#' set.seed(1); n <- 300
#' phi <- runif(n, -pi, pi)                       # circular covariate (radians)
#' y <- 2 + 1.5 * sin(phi) + 0.8 * cos(2 * phi) + rnorm(n) * 0.3
#' b <- circ_gam(list(y ~ s(phi, bs = "cc"), ~ s(phi, bs = "cc")),
#'               data = data.frame(y, phi), family = gausslss())
#' head(predict(b, type = "response"))            # columns named mu, tau
#' plot(b, view = "both")                         # the l~c "can" + the mean panel
#' @seealso
#' \code{\link{gammalss}}, \code{\link{circ_gam}}, \code{\link[mgcv]{gaulss}}
#' @export
gausslss <- function(link = list("identity", "logb"), b = 0.01) {
  if (length(link) != 2) stop("gaulss requires 2 links specified as character strings")
  okLinks <- list(c("inverse", "log", "identity", "sqrt"), "logb")
  stats <- list()
  if (link[[1]] %in% okLinks[[1]]) stats[[1]] <- stats::make.link(link[[1]]) else
    stop(link[[1]], " link not available for mu parameter of gausslss")
  fam <- structure(list(link = link[[1]], canonical = "none",
                        linkfun = stats[[1]]$linkfun, mu.eta = stats[[1]]$mu.eta),
                   class = "family")
  fam <- mgcv::fix.family.link(fam)
  stats[[1]]$d2link <- fam$d2link
  stats[[1]]$d3link <- fam$d3link
  stats[[1]]$d4link <- fam$d4link
  if (link[[2]] %in% okLinks[[2]]) {           # the logb link: tau = 1/(exp(eta)+b)
    stats[[2]] <- list()
    stats[[2]]$valideta <- function(eta) TRUE
    stats[[2]]$link <- link[[2]]
    stats[[2]]$linkfun <- eval(parse(text = paste("function(mu) log(1/mu -", b, ")")))
    stats[[2]]$linkinv <- eval(parse(text = paste("function(eta) 1/(exp(eta) +", b, ")")))
    stats[[2]]$mu.eta <- eval(parse(text =
      paste("function(eta) { ee <- exp(eta); -ee/(ee +", b, ")^2 }")))
    stats[[2]]$d2link <- eval(parse(text =
      paste("function(mu) { mub <- pmax(1 - mu *", b, ",.Machine$double.eps);(2*mub-1)/(mub*mu)^2}")))
    stats[[2]]$d3link <- eval(parse(text =
      paste("function(mu) { mub <- pmax(1 - mu *", b, ",.Machine$double.eps);((1-mub)*mub*6-2)/(mub*mu)^3}")))
    stats[[2]]$d4link <- eval(parse(text =
      paste("function(mu) { mub <- pmax(1 - mu *", b, ",.Machine$double.eps);(((24*mub-36)*mub+24)*mub-6)/(mub*mu)^4}")))
  } else stop(link[[2]], " link not available for precision parameter of gausslss")

  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    type <- match.arg(type)
    rsd <- object$y - object$fitted[, 1]
    if (type == "response") return(rsd) else
      return(rsd * object$fitted[, 2])         # (y - mu)/sigma
  }

  postproc <- expression({
    object$null.deviance <- sum(((object$y - mean(object$y)) * object$fitted[, 2])^2)
  })

  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                 D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
    ## deriv: 0 eval; 1 grad+Hess; 2 diag of dHess; 3 dHess; 4 everything.
    if (!is.null(offset)) offset[[3]] <- 0
    if (is.null(wt)) wt <- 1               # prior weights; NULL from sandwich/direct callers
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
    mu <- family$linfo[[1]]$linkinv(eta)
    tau <- family$linfo[[2]]$linkinv(eta1)     # tau = 1/sig
    n <- length(y)
    l1 <- matrix(0, n, 2)
    ymu <- y - mu; ymu2 <- ymu^2; tau2 <- tau^2
    l0 <- -.5 * ymu2 * tau2 - .5 * log(2 * pi) + log(tau)
    l <- sum(wt * l0)

    if (deriv) {
      l1[, 1] <- tau2 * ymu
      l1[, 2] <- 1 / tau - tau * ymu2
      ## second derivatives, order mm, ms, ss
      l2 <- cbind(-tau2, 2 * l1[, 1] / tau, -ymu2 - 1 / tau2)
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta), family$linfo[[2]]$mu.eta(eta1))
      g2 <- cbind(family$linfo[[1]]$d2link(mu), family$linfo[[2]]$d2link(tau))
    }

    l3 <- l4 <- g3 <- g4 <- 0
    if (deriv > 1) {
      ## third derivatives, order mmm, mms, mss, sss
      l3 <- cbind(0, -2 * tau, 2 * ymu, 2 / tau^3)
      g3 <- cbind(family$linfo[[1]]$d3link(mu), family$linfo[[2]]$d3link(tau))
    }
    if (deriv > 3) {
      ## fourth derivatives, order mmmm, mmms, mmss, msss, ssss
      l4 <- cbind(0, 0, -2, 0, -6 / tau2^2)
      g4 <- cbind(family$linfo[[1]]$d4link(mu), family$linfo[[2]]$d4link(tau))
    }

    if (deriv) {
      i2 <- family$tri$i2; i3 <- family$tri$i3; i4 <- family$tri$i4
      ## prior weights: a weighted log-likelihood scales every per-observation
      ## derivative row by wt. gamlss.etamu is linear per row, so scaling these
      ## inputs == scaling the eta-derivatives gamlss.gH builds.
      l1 <- wt * l1
      l2 <- wt * l2
      if (is.matrix(l3)) l3 <- wt * l3
      if (is.matrix(l4)) l4 <- wt * l4
      de <- mgcv::gamlss.etamu(l1, l2, l3, l4, ig1, g2, g3, g4, i2, i3, i4, deriv - 1)
      ret <- mgcv::gamlss.gH(X, jj, de$l1, de$l2, i2, l3 = de$l3, i3 = i3,
                             l4 = de$l4, i4 = i4, d1b = d1b, d2b = d2b,
                             deriv = deriv - 1, fh = fh, D = D, sandwich = sandwich)
      if (ncv) { ret$l1 <- de$l1; ret$l2 <- de$l2; ret$l3 <- de$l3 }
    } else ret <- list()
    ret$l <- l; ret$l0 <- l0; ret
  }

  sandwich <- function(y, X, coef, wt, family, offset = NULL) {
    ll(y, X, coef, wt, family, offset = NULL, deriv = 1, sandwich = TRUE)$lbb
  }

  initialize <- expression({
    ## regress g(y) on the mean model matrix, then the log absolute residuals on
    ## the log-sigma model matrix. E only regularizes the start: pen.reg() scales
    ## it to the model-matrix norm, unless attr(E, "use.unscaled") forces the raw
    ## stacked-penalty QR. gaulss's discrete (is.list(x)) branch is dropped --
    ## circlss is dense-only.
    n <- rep(1, nobs)
    use.unscaled <- if (!is.null(attr(E, "use.unscaled"))) TRUE else FALSE
    if (is.null(start)) {
      jj <- attr(x, "lpi")
      if (!is.null(offset)) offset[[3]] <- 0
      yt1 <- if (family$link[[1]] == "identity") y else
        family$linfo[[1]]$linkfun(abs(y) + max(y) * 1e-7)
      if (!is.null(offset[[1]])) yt1 <- yt1 - offset[[1]]
      start <- rep(0, ncol(x))
      x1 <- x[, jj[[1]], drop = FALSE]; e1 <- E[, jj[[1]], drop = FALSE]
      if (use.unscaled) {
        x1 <- rbind(x1, e1)
        startji <- qr.coef(qr(x1), c(yt1, rep(0, nrow(E))))
        startji[!is.finite(startji)] <- 0
      } else startji <- pen.reg(x1, e1, yt1)
      start[jj[[1]]] <- startji
      lres1 <- log(abs(y - family$linfo[[1]]$linkinv(x[, jj[[1]], drop = FALSE] %*% start[jj[[1]]])))
      if (!is.null(offset[[2]])) lres1 <- lres1 - offset[[2]]
      x1 <- x[, jj[[2]], drop = FALSE]; e1 <- E[, jj[[2]], drop = FALSE]
      if (use.unscaled) {
        x1 <- rbind(x1, e1)
        startji <- qr.coef(qr(x1), c(lres1, rep(0, nrow(E))))
        startji[!is.finite(startji)] <- 0
      } else startji <- pen.reg(x1, e1, lres1)
      start[jj[[2]]] <- startji
    }
  })

  rd <- function(mu, wt, scale) {
    ## prior weights scale the variance (the gaulss convention)
    stats::rnorm(nrow(mu), mu[, 1], sqrt(scale / wt) / mu[, 2])
  }

  structure(list(family = "gausslss", ll = ll, link = paste(link), nlp = 2,
                 param_names = c("mu", "tau"),
                 param_circular = c(FALSE, FALSE),
                 response_circular = FALSE,
                 sandwich = sandwich, tri = mgcv::trind.generator(2),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 2, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
