#' Weight-aware gamma location-scale family
#'
#' A gamma location-scale family for distributional regression of a
#' \emph{positive}, right-skewed response, with a (log) mean and a log scale each
#' modelled by its own linear predictor. It is a weight-aware, metadata-carrying
#' adaptation of \pkg{mgcv}'s \code{\link[mgcv]{gammals}}: unlike \code{gammals}
#' it honours prior \code{weights} (needed for a weighted MLE and for EM
#' mixtures), and it carries the circlss parameter metadata so
#' \code{\link{circ_gam}} treats it as a first-class location-scale family.
#'
#' @param link Two-element list of link names, following
#'   \code{\link[mgcv]{gammals}}: \code{"identity"} for the (log) mean and
#'   \code{"log"} for the log scale.
#' @param b The log-scale link's offset, as in \code{\link[mgcv]{gammals}}.
#' @return An object of class \code{c("general.family", "extended.family", "family")}
#'   for use with \code{\link[mgcv]{gam}} (or its front end \code{\link{circ_gam}}).
#' @details
#' This is the positive-response member of the circlss linear-circular (l~c)
#' leg: a positive, skewed quantity that varies around a cycle -- rainfall by
#' season, a concentration or rate by time of day, a speed by direction -- over a
#' circular covariate fitted with a cyclic smooth. As for \code{\link{gausslss}},
#' \code{\link{circ_gam}} places the fitted mean on the "can".
#'
#' The parameterization follows \code{gammals}: a (log) mean and a log scale. At
#' unit weights the fit matches \code{gammals}; integer prior weights reproduce a
#' row-replicated fit. Log-likelihood derivatives up to fourth order are
#' implemented, so the family supports full Newton REML (\code{method = "REML"});
#' \code{optimizer = "efs"} also works.
#'
#' This family adapts GPL-licensed code from \pkg{mgcv}; see the package's
#' \code{inst/COPYRIGHTS}.
#' @examples
#' library(mgcv)
#' set.seed(1); n <- 300
#' phi <- runif(n, -pi, pi)                       # circular covariate (radians)
#' y <- rgamma(n, shape = 4, rate = 4 / exp(0.4 + 0.8 * sin(phi)))
#' b <- circ_gam(list(y ~ s(phi, bs = "cc"), ~ s(phi, bs = "cc")),
#'               data = data.frame(y, phi), family = gammalss())
#' head(predict(b, type = "response"))            # columns named mu, scale
#' @seealso
#' \code{\link{gausslss}}, \code{\link{circ_gam}}, \code{\link[mgcv]{gammals}}
#' @export
gammalss <- function(link = list("identity", "log"), b = -7) {
  ## first deal with links and their derivatives...
  if (length(link) != 2) stop("gammalss requires 2 links specified as character strings")
  okLinks <- list(c("identity"), c("identity", "log"))
  stats <- list()
  for (i in 1:2) {
    if (link[[i]] %in% okLinks[[i]]) stats[[i]] <- stats::make.link(link[[i]]) else
      stop(link[[i]], " link not available for gammalss")
    fam <- structure(list(link = link[[i]], canonical = "none", linkfun = stats[[i]]$linkfun,
                          mu.eta = stats[[i]]$mu.eta),
                     class = "family")
    fam <- mgcv::fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }
  if (link[[2]] == "log") { ## g^{-1}(eta) = b + log(1+exp(eta-b)) link
    stats[[2]]$valideta <- function(eta) TRUE

    stats[[2]]$linkfun <- eval(parse(text = paste("function(mu,b=", b, ") {\n eta <- mub <- mu-b;\n",
      "ii <- mub < .Machine$double.eps;\n if (any(ii)) eta[ii] <- log(.Machine$double.eps)+b;\n",
      "jj <- mub > -log(.Machine$double.eps);if (any(jj)) eta[jj] <- mub[jj]+b;\n",
      "jj <- !jj & !ii;if (any(jj)) eta[jj] <- log(exp(mub[jj])-1)+b;eta }")))

    stats[[2]]$mu.eta <- eval(parse(text = paste("function(eta,b=", b, ") {\n",
      "eta <- eta - b; ii <- eta < 0;eta <- exp(-eta*sign(eta))\n",
      "if (any(ii)) { ei <- eta[ii];eta[ii] <- ei/(1+ei)}\n",
      "ii <- !ii;if (any(ii)) eta[ii] <- 1/(1+eta[ii])\n",
      "eta }\n")))

    stats[[2]]$linkinv <- eval(parse(text = paste("function(eta,b=", b, ") {\n",
      "mu <- eta;ii <- eta-b < -log(.Machine$double.eps)\n",
      "if (any(ii)) mu[ii] <- b + log1p(exp(eta[ii]-b))\n",
      "mu }\n")))

    stats[[2]]$d2link <- eval(parse(text = paste("function(mu,b=", b, ") {\n",
      "mub <- mu - b; mub <-  exp(-mub*sign(mub))\n",
      "-mub/(mub-1)^2 }\n")))

    stats[[2]]$d3link <- eval(parse(text = paste("function(mu,b=", b, ") {\n",
      "mub <- mu - b; sm <- -sign(mub);mub <- exp(mub*sm) \n",
      "sm*(mub+mub^2)/(mub-1)^3 }\n")))

    stats[[2]]$d4link <- eval(parse(text = paste("function(mu,b=", b, ") {\n",
      "mub <- mu - b; sm <- -sign(mub);mub <- exp(mub*sm)\n",
      "sm*(mub+4*mub^2+mub^3)/(mub-1)^4 }\n")))
  }

  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    mu <- object$fitted.values[, 1]
    rho <- object$fitted.values[, 2]
    y <- object$y
    type <- match.arg(type)
    if (type == "deviance") {
      rsd <- 2 * ((y - mu) / mu - log(y / mu)) * exp(-rho)
      rsd <- sqrt(pmax(0, rsd)) * sign(y - mu)
    } else if (type == "pearson") {
      rsd <- (y - mu) / (exp(rho * .5) * mu)
    } else {
      rsd <- y - mu
    }
    rsd
  } ## gammalss residuals

  postproc <- expression({
    ## code to evaluate in estimate.gam, to evaluate null deviance
    object$fitted.values[, 1] <- exp(object$fitted.values[, 1])
    .my <- mean(object$y)
    object$null.deviance <- sum(((object$y - .my) / .my - log(object$y / .my)) * exp(-object$fitted.values[, 2])) * 2
  })

  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                 D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
    ## function defining the gamlss gamma model log lik.
    ## deriv: 0 - eval
    ##        1 - grad and Hess
    ##        2 - diagonal of first deriv of Hess
    ##        3 - first deriv of Hess
    ##        4 - everything.
    if (!is.null(offset)) offset[[3]] <- 0
    if (is.null(wt)) wt <- 1               # prior weights; NULL from sandwich/direct callers
    jj <- attr(X, "lpi")                   # extract linear predictor index
    if (is.null(eta)) {
      eta <- X[, jj[[1]], drop = FALSE] %*% coef[jj[[1]]]
      if (!is.null(offset[[1]])) eta <- eta + offset[[1]]   # log mu
      etat <- X[, jj[[2]], drop = FALSE] %*% coef[jj[[2]]]
      if (!is.null(offset[[2]])) etat <- etat + offset[[2]]
    } else {
      etat <- eta[, 2]
      eta <- eta[, 1]
    }
    mu <- family$linfo[[1]]$linkinv(eta)   # mean
    th <- family$linfo[[2]]$linkinv(etat)  # log sigma

    eth <- exp(-th)                        # 1/exp1^th;
    logy <- log(y)
    ethmu <- exp(-th - mu)
    ethmuy <- ethmu * y
    etlymt <- eth * (logy - mu - th)
    n <- length(y)

    l0 <- etlymt - logy - ethmuy - lgamma(eth)   # l
    l <- sum(wt * l0)
    if (!is.finite(l)) return(list(l = l, l0 = l0))

    if (deriv) {
      l1 <- matrix(0, n, 2)

      l1[, 1] <- ethmuy - eth                # lm
      digeth <- digamma(eth)
      l1[, 2] <- -etlymt + ethmuy + eth * digeth - eth   # lt

      ## the second derivatives
      l2 <- matrix(0, n, 3)
      ## order mm,mt,tt
      l2[, 1] <- -ethmuy                     # lmm
      l2[, 2] <- eth - ethmuy                 # lmt
      eth2 <- eth^2; treth <- trigamma(eth)
      l2[, 3] <- etlymt - ethmuy - treth * eth2 - eth * digeth + 2 * eth   # ltt
      ## need some link derivatives for derivative transform
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta), family$linfo[[2]]$mu.eta(etat))
      g2 <- cbind(family$linfo[[1]]$d2link(mu), family$linfo[[2]]$d2link(th))
    }

    l3 <- l4 <- g3 <- g4 <- 0               # defaults

    if (deriv > 1) {
      ## the third derivatives
      ## order mmm,mmt,mtt,ttt
      l3 <- matrix(0, n, 4)
      l3[, 1] <- ethmuy                       # lmmm
      l3[, 2] <- ethmuy                       # lmmt
      l3[, 3] <- ethmuy - eth                 # lmtt
      eth3 <- eth2 * eth; g3eth <- psigamma(eth, deriv = 2)
      l3[, 4] <- -etlymt + ethmuy + g3eth * eth3 + 3 * treth * eth2 + eth * digeth - 3 * eth   # lttt
      g3 <- cbind(family$linfo[[1]]$d3link(mu), family$linfo[[2]]$d3link(th))
    }

    if (deriv > 3) {
      ## the fourth derivatives
      ## order mmmm,mmmt,mmtt,mttt,tttt
      l4 <- matrix(0, n, 5)
      l4[, 1] <- -ethmuy                      # lmmmm
      l4[, 2] <- -ethmuy                      # lmmmt
      l4[, 3] <- -ethmuy                      # lmmtt
      l4[, 4] <- eth - ethmuy                 # lmttt
      eth4 <- eth3 * eth
      l4[, 5] <- etlymt - ethmuy - psigamma(eth, deriv = 3) * eth4 - 6 * g3eth * eth3 -
                 7 * treth * eth2 - eth * digeth + 4 * eth   # ltttt
      g4 <- cbind(family$linfo[[1]]$d4link(mu), family$linfo[[2]]$d4link(th))
    }

    if (deriv) {
      i2 <- family$tri$i2; i3 <- family$tri$i3
      i4 <- family$tri$i4
      ## prior weights: a weighted log-likelihood scales every per-observation
      ## derivative row by wt. gamlss.etamu is linear per row, so scaling these
      ## inputs == scaling the eta-derivatives gamlss.gH builds.
      l1 <- wt * l1
      l2 <- wt * l2
      if (is.matrix(l3)) l3 <- wt * l3
      if (is.matrix(l4)) l4 <- wt * l4

      ## transform derivates w.r.t. mu to derivatives w.r.t. eta...
      de <- mgcv::gamlss.etamu(l1, l2, l3, l4, ig1, g2, g3, g4, i2, i3, i4, deriv - 1)

      ## get the gradient and Hessian...
      ret <- mgcv::gamlss.gH(X, jj, de$l1, de$l2, i2, l3 = de$l3, i3 = i3,
                             l4 = de$l4, i4 = i4, d1b = d1b, d2b = d2b,
                             deriv = deriv - 1, fh = fh, D = D, sandwich = sandwich)
      if (ncv) {
        ret$l1 <- de$l1; ret$l2 <- de$l2; ret$l3 <- de$l3
      }
    } else ret <- list()
    ret$l <- l; ret$l0 <- l0; ret
  } ## end ll gammalss

  sandwich <- function(y, X, coef, wt, family, offset = NULL) {
    ## compute filling for sandwich estimate of cov matrix
    ll(y, X, coef, wt, family, offset = NULL, deriv = 1, sandwich = TRUE)$lbb
  }

  initialize <- expression({
    ## regress log(y) on the mean model matrix, then the (scale-link of the) log
    ## absolute residuals on the scale model matrix. E only regularizes the start:
    ## pen.reg() scales it to the model-matrix norm, unless attr(E, "use.unscaled")
    ## forces the raw stacked-penalty QR. gammals's discrete (is.list(x)) branch is
    ## dropped -- circlss is dense-only.
    n <- rep(1, nobs)
    use.unscaled <- if (!is.null(attr(E, "use.unscaled"))) TRUE else FALSE
    if (is.null(start)) {
      jj <- attr(x, "lpi")
      if (!is.null(offset)) offset[[3]] <- 0
      yt1 <- log(y + max(y) * .Machine$double.eps^.75)
      if (!is.null(offset[[1]])) yt1 <- yt1 - offset[[1]]
      start <- rep(0, ncol(x))
      x1 <- x[, jj[[1]], drop = FALSE]; e1 <- E[, jj[[1]], drop = FALSE]
      if (use.unscaled) {
        x1 <- rbind(x1, e1)
        startji <- qr.coef(qr(x1), c(yt1, rep(0, nrow(E))))
        startji[!is.finite(startji)] <- 0
      } else startji <- pen.reg(x1, e1, yt1)
      start[jj[[1]]] <- startji
      lres1 <- family$linfo[[2]]$linkfun(log(abs(y - family$linfo[[1]]$linkinv(x[, jj[[1]], drop = FALSE] %*% start[jj[[1]]]))))
      if (!is.null(offset[[2]])) lres1 <- lres1 - offset[[2]]
      x1 <- x[, jj[[2]], drop = FALSE]; e1 <- E[, jj[[2]], drop = FALSE]
      if (use.unscaled) {
        x1 <- rbind(x1, e1)
        startji <- qr.coef(qr(x1), c(lres1, rep(0, nrow(E))))
        startji[!is.finite(startji)] <- 0
      } else startji <- pen.reg(x1, e1, lres1)
      start[jj[[2]]] <- startji
    }
  }) ## initialize gammalss

  rd <- function(mu, wt, scale) {
    ## simulate responses
    phi <- exp(mu[, 2])
    stats::rgamma(nrow(mu), shape = 1 / phi, scale = mu[, 1] * phi)
  } ## rd

  predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                      beta = NULL, off = NULL, Vb = NULL) {
    ## predict.gam(..., type = "response") uses this: the mean is exp(eta_1)
    ## (gammals carries the log on the mean under the identity link, so the
    ## default link inverse would return the log-mean), the scale is
    ## linfo[[2]]$linkinv(eta_2). gammals's discrete branch is dropped -- dense-only.
    if (is.null(eta)) {
      lpi <- attr(X, "lpi"); if (is.null(lpi)) lpi <- list(1:ncol(X))
      eta <- matrix(0, nrow(X), 2); ve <- matrix(0, nrow(X), 2)
      for (i in 1:2) {
        Xi <- X[, lpi[[i]], drop = FALSE]
        eta[, i] <- Xi %*% beta[lpi[[i]]]
        if (!is.null(off[[i]])) eta[, i] <- eta[, i] + off[[i]]
        if (se) ve[, i] <- drop(pmax(0, rowSums((Xi %*% Vb[lpi[[i]], lpi[[i]]]) * Xi)))
      }
    } else se <- FALSE
    gamma <- cbind(exp(eta[, 1]), family$linfo[[2]]$linkinv(eta[, 2]))
    if (se) {
      vp <- gamma
      vp[, 1] <- abs(gamma[, 1]) * sqrt(ve[, 1])
      vp[, 2] <- abs(family$linfo[[2]]$mu.eta(eta[, 2])) * sqrt(ve[, 2])
      return(list(fit = gamma, se.fit = vp))
    }
    list(fit = gamma)
  } ## gammals predict (dense)

  structure(list(family = "gammalss", ll = ll, link = paste(link), nlp = 2,
                 param_names = c("mu", "scale"),
                 param_circular = c(FALSE, FALSE),
                 response_circular = FALSE,
                 sandwich = sandwich, tri = mgcv::trind.generator(2),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd, predict = predict,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 2, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
} ## end gammalss
