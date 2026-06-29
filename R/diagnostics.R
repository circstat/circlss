#' Circular regression residuals
#'
#' The residual primitive for circlss fits -- the circular analogue of
#' \code{\link[stats]{residuals}}, and the quantity every \code{\link{circ_check}}
#' panel is a function of. \dQuote{Observed minus fitted} is undefined when both
#' are angles, so \code{circ_resid} returns one of four residual definitions that
#' are well posed on the circle. It dispatches over both front doors
#' (\code{\link{circ_lm}} and \code{\link{circ_gam}}) so the residual vocabulary
#' is identical across the parametric and penalized fits.
#'
#' @param object A fitted \code{\link{circ_lm}} or \code{\link{circ_gam}} model.
#' @param type The residual definition. \code{"quantile"} (default) is the
#'   probability-integral-transform residual, calibrated even when the
#'   concentration varies across observations; \code{"deviance"} is the signed
#'   root of the per-observation deviance, constructed to be approximately
#'   \eqn{N(0, 1)} under a good fit; \code{"angular"} is the wrapped residual
#'   \eqn{y - \hat\mu \in (-\pi, \pi]} (the raw response residual, returned on the
#'   response scale for a linear-response fit); \code{"pearson"} is the
#'   score-standardized residual.
#' @param nsim Number of simulation replicates for the \code{"quantile"} residual
#'   when the family has no closed-form distribution function (every family except
#'   the von Mises and Gaussian cases, which are computed analytically). Set a seed
#'   for a reproducible simulated residual.
#' @param scale For \code{type = "quantile"} only: \code{"uniform"} returns the
#'   probability-integral transform on \eqn{(0, 1)} (pairs with the Watson
#'   \eqn{U^2} uniform Q-Q); \code{"normal"} maps it through
#'   \code{\link[stats]{qnorm}} to the Dunn--Smyth \eqn{N(0, 1)} residual (pairs
#'   with a normal Q-Q).
#' @return A numeric vector of residuals, one per observation, carrying an
#'   \code{attr(, "type")} tag (and, for \code{"quantile"}, an
#'   \code{attr(, "scale")} tag).
#' @details
#' For a circular response the quantile residual is computed in the
#' \emph{residual frame}: the wrapped residual \eqn{y - \hat\mu} is ranked against
#' its fitted (von Mises/Gaussian) or simulated distribution, with the cut placed
#' at the antipode \eqn{\pm\pi} -- the least probable region of a concentrated
#' residual, which minimises the wrapping artifact. The von Mises (\code{"cl"},
#' \code{"cc"}, \code{\link{vmlss}}) and Gaussian (\code{"lc"},
#' \code{\link{gausslss}}) cases use a closed-form distribution function and are
#' deterministic; every other family is transformed by simulation from the
#' family's random-deviate generator, so a seed is needed for reproducibility.
#'
#' The deviance, angular and Pearson residuals are read straight from each
#' family's own \code{residuals} method (for \code{\link{circ_gam}}) or
#' reconstructed from the stored von Mises / least-squares fit (for
#' \code{\link{circ_lm}}); the centring is always taken from the fitted direction
#' the model reports, including the derived direction of \code{\link{pnlss}}.
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996) Randomized quantile residuals.
#' \emph{Journal of Computational and Graphical Statistics} 5, 236-244.
#'
#' Fisher, N. I. (1993) \emph{Statistical Analysis of Circular Data}. Cambridge
#' University Press.
#' @seealso
#' \code{\link{circ_check}} for the diagnostic panel grid built on these
#' residuals; \code{\link{circ_lm}}, \code{\link{circ_gam}}.
#' @examples
#' set.seed(1)
#' n <- 80
#' x <- rnorm(n)
#' theta <- (1 + 2 * atan(1.5 * x) + rnorm(n) / 4) %% (2 * pi)
#' m <- circ_lm(theta ~ x, data.frame(theta, x), type = "cl")
#'
#' r <- circ_resid(m, type = "quantile")   # von Mises: closed-form, deterministic
#' head(r)
#' head(circ_resid(m, type = "angular"))    # wrapped y - mu_hat in (-pi, pi]
#' @export
circ_resid <- function(object,
                       type = c("quantile", "deviance", "angular", "pearson"),
                       nsim = 1000L,
                       scale = c("uniform", "normal")) {
  type <- match.arg(type)
  scale <- match.arg(scale)
  parts <- .circ_resid_parts(object)
  r <- switch(type,
    angular  = parts$resid_response(),
    deviance = parts$resid_deviance(),
    pearson  = {
      if (is.null(parts$resid_pearson))
        stop("Pearson residuals are not available for this fit.", call. = FALSE)
      parts$resid_pearson()
    },
    quantile = .circ_quantile_resid(parts, as.integer(nsim), scale))
  r <- as.numeric(r)
  attr(r, "type") <- type
  if (type == "quantile") attr(r, "scale") <- scale
  r
}

## ===========================================================================
## Parts adapter: one place where the residual math lives, branching on the two
## front-door classes. Returns the per-observation primitives circ_resid() and
## circ_check() both consume: the response and fitted direction, closures for the
## three classical residual types, an analytic distribution function (or NULL),
## and a simulator (or NULL).
## ===========================================================================

.circ_resid_parts <- function(object) {
  if (inherits(object, "circ_gam")) return(.circ_parts_gam(object))
  if (inherits(object, "circ_lm"))  return(.circ_parts_lm(object))
  stop("circ_resid()/circ_check() support 'circ_gam' and 'circ_lm' fits only; ",
       "got an object of class ", paste(sQuote(class(object)), collapse = "/"),
       ".", call. = FALSE)
}

## The fitted direction the model reports: the derived direction when the family
## defines one (pnlss's atan2(mu2, mu1)), else the first circular parameter
## column (vmlss's mu, ...), else the mean column for a linear response. Mirrors
## plot.circ_gam()'s location logic, so it never assumes a column called "mu".
.circ_fitted_dir <- function(fam, fm) {
  if (!is.null(fam$derived) && !is.null(fam$derived$direction))
    return(as.numeric(fam$derived$direction(fm)))
  pc <- which(as.logical(fam$param_circular))
  if (length(pc)) return(as.numeric(fm[, pc[1]]))
  as.numeric(fm[, 1])
}

.circ_parts_gam <- function(object) {
  fam <- object$family
  ## RAW (fit-frame) fitted, consistent with object$y which is also fit-frame --
  ## using the rotating fitted.circ_gam here would double-rotate a centered fit.
  ## Residuals (y - dir) and the cdf are rotation-invariant; only the returned
  ## display y / fitted_dir are rotated back by the fit's center.
  fm  <- .circ_fitted_raw(object)              # named response-scale matrix
  if (is.null(dim(fm))) fm <- cbind(fm)
  y   <- as.numeric(object$y)
  n   <- length(y)
  ref <- if (!is.null(object$circ_center)) object$circ_center else 0
  resp_circ <- !is.null(fam$param_names) && !isFALSE(fam$response_circular)
  dir <- .circ_fitted_dir(fam, fm)

  ## analytic distribution function for the two closed-form cases only
  cdf <- NULL
  if (identical(fam$family, "vmlss")) {
    kappa <- as.numeric(fm[, 2])
    cdf <- function() .circ_pvonmises(y, dir, kappa)
  } else if (identical(fam$family, "gausslss")) {
    mu <- as.numeric(fm[, 1]); tau <- as.numeric(fm[, 2])   # tau = 1 / sigma
    cdf <- function() stats::pnorm((y - mu) * tau)
  }

  ## simulator: the family random-deviate generator, called nsim times on the
  ## whole fitted matrix (no column indexing -- the family decodes its columns)
  sim <- if (!is.null(fam$rd)) {
    w  <- object$prior.weights; if (is.null(w)) w <- rep(1, n)
    sc <- object$scale; if (is.null(sc) || !is.finite(sc)) sc <- 1
    function(nsim)
      matrix(unlist(lapply(seq_len(nsim), function(i) fam$rd(fm, w, sc))), n, nsim)
  } else NULL

  list(
    y = if (resp_circ) wrap(y + ref) else y,
    response_circular = resp_circ,
    fitted_dir = if (resp_circ) wrap(dir + ref) else dir,
    n = n, family_label = fam$family,
    resid_response = function() as.numeric(stats::residuals(object, type = "response")),
    resid_deviance = function() as.numeric(stats::residuals(object, type = "deviance")),
    resid_pearson  = function() as.numeric(stats::residuals(object, type = "pearson")),
    cdf = cdf, simulate = sim)
}

.circ_parts_lm <- function(object) {
  switch(object$type,
    cl = ,
    cc = .circ_parts_lm_vm(object),
    lc = .circ_parts_lm_lc(object),
    stop("unknown circ_lm type ", sQuote(object$type), ".", call. = FALSE))
}

## cl and cc are both von Mises around the fitted direction: cl carries a
## per-observation kappa (mean/kappa/mixed), cc a single residual concentration.
.circ_parts_lm_vm <- function(object) {
  y   <- wrap(as.numeric(object$frame[[object$response]]))
  dir <- wrap(as.numeric(object$fitted))
  n   <- length(y)
  kappa <- rep_len(as.numeric(object$kappa), n)
  ang <- function() wrap(y - dir)
  list(
    y = y, response_circular = TRUE, fitted_dir = dir,
    n = n, family_label = paste0("circ_lm:", object$type),
    resid_response = ang,
    resid_deviance = function() {
      d <- ang(); sign(sin(d)) * sqrt(pmax(2 * kappa * (1 - cos(d)), 0))
    },
    resid_pearson = function() {
      d <- ang(); v <- A1(kappa) / kappa; v[kappa == 0] <- 0.5; sin(d) / sqrt(v)
    },
    cdf = function() .circ_pvonmises(y, dir, kappa),
    simulate = function(nsim) {
      vm <- vmlss(); fmm <- cbind(dir, kappa); w <- rep(1, n)
      matrix(unlist(lapply(seq_len(nsim), function(i) vm$rd(fmm, w, 1))), n, nsim)
    })
}

## lc is an ordinary Gaussian least-squares fit on the harmonic basis; angular
## residuals do not apply, so the response residual is the raw y - fit and the
## quantile residual is the ordinary Gaussian probability-integral transform.
.circ_parts_lm_lc <- function(object) {
  fit   <- as.numeric(object$fitted)
  resid <- as.numeric(object$residuals)
  y     <- fit + resid
  sigma <- object$sigma
  n     <- length(y)
  list(
    y = y, response_circular = FALSE, fitted_dir = fit,
    n = n, family_label = "circ_lm:lc",
    resid_response = function() resid,
    resid_deviance = function() resid / sigma,
    resid_pearson  = function() resid / sigma,
    cdf = function() stats::pnorm(resid / sigma),
    simulate = function(nsim)
      matrix(stats::rnorm(n * nsim, rep(fit, nsim), sigma), n, nsim))
}

## ===========================================================================
## Quantile (PIT) residual
## ===========================================================================

.circ_quantile_resid <- function(parts, nsim, scale) {
  u <- if (!is.null(parts$cdf)) parts$cdf()
       else if (!is.null(parts$simulate)) .circ_sim_pit(parts, nsim)
       else stop("no distribution function or simulator for family ",
                 sQuote(parts$family_label), "; cannot form a quantile residual.",
                 call. = FALSE)
  eps <- 1e-6                                  # keep qnorm() and 0/1 ties finite
  u <- pmin(pmax(u, eps), 1 - eps)
  if (scale == "normal") stats::qnorm(u) else u
}

## Simulation PIT (the DHARMa route): rank each observed residual against nsim
## simulated residuals from the fitted model. For a circular response the ranking
## is done in the residual frame y - mu_hat wrapped to (-pi, pi], i.e. cut at the
## antipode. The Hazen rank ((# below) + 1/2) / (nsim + 1) is strictly inside
## (0, 1) -- no ties at the edges, no jitter, deterministic given the draws.
.circ_sim_pit <- function(parts, nsim) {
  sims <- parts$simulate(nsim)                 # n x nsim response draws
  dsim <- sweep(sims, 1, parts$fitted_dir, "-")
  dobs <- parts$y - parts$fitted_dir
  if (parts$response_circular) { dsim <- wrap(dsim); dobs <- wrap(dobs) }
  less <- rowSums(dsim < dobs)
  eq   <- rowSums(dsim == dobs)
  (less + 0.5 * eq + 0.5) / (nsim + 1)
}

## ===========================================================================
## Deterministic goodness-of-fit helpers (von Mises CDF, Watson U-squared)
## ===========================================================================

## von Mises distribution function with the origin at the antipode -pi, evaluated
## in the residual frame d = theta - mu wrapped to (-pi, pi]:
##   F(d) = (d + pi)/(2*pi) + (1/pi) * sum_{j>=1} (I_j(kappa)/I_0(kappa)) sin(j d)/j
## (the Fourier integral of the von Mises density). F(-pi) = 0, F(mu) = 1/2,
## F(pi) = 1, so F(d) ~ U(0, 1) under the model -- the analytic PIT. The Bessel
## ratio uses exponentially scaled besselI so it stays finite at large kappa.
## Vectorised over theta, mu and kappa; kappa = 0 collapses to the circular
## uniform. Matches circular::pvonmises (from = mu - pi) to machine precision.
.circ_pvonmises <- function(theta, mu = 0, kappa = 1, jmax = 60L, tol = 1e-12) {
  d <- wrap(theta - mu)
  out <- (d + pi) / (2 * pi)
  i0 <- besselI(kappa, 0, expon.scaled = TRUE)
  for (j in seq_len(jmax)) {
    rj <- besselI(kappa, j, expon.scaled = TRUE) / i0
    out <- out + rj * sin(j * d) / (j * pi)
    if (max(abs(rj)) < tol) break
  }
  pmin(pmax(out, 0), 1)
}

## Watson's U-squared test of uniformity for PIT values u in [0, 1) -- the
## rotation-invariant goodness-of-fit statistic (the Kolmogorov test is wrong on
## the circle, where the origin is arbitrary). For sorted u_(1) <= ... <= u_(n):
##   U2 = sum( (u_(i) - (2i-1)/(2n))^2 ) - n*(ubar - 1/2)^2 + 1/(12n).
## The asymptotic upper-tail p-value is the standard theta series
##   p = 2 * sum_{m>=1} (-1)^(m-1) exp(-2 m^2 pi^2 U2),
## truncated at convergence and clamped to [0, 1] (the series is an asymptotic
## approximation that can drift outside the unit interval for tiny U2). Returns
## list(stat, p); deterministic given u, so it is unit-testable.
.circ_watson_u2 <- function(u) {
  u <- sort(as.numeric(u))
  n <- length(u)
  i <- seq_len(n)
  ubar <- mean(u)
  stat <- sum((u - (2 * i - 1) / (2 * n))^2) - n * (ubar - 0.5)^2 + 1 / (12 * n)
  p <- 0
  for (m in seq_len(50L)) {
    term <- (-1)^(m - 1) * exp(-2 * m^2 * pi^2 * stat)
    p <- p + term
    if (abs(term) < 1e-12) break
  }
  list(stat = stat, p = min(max(2 * p, 0), 1))
}

#' Diagnostic panels for a circular regression fit
#'
#' The diagnostic display for circlss fits -- the circular analogue of
#' \code{stats:::plot.lm} and \code{\link[mgcv]{gam.check}}, and the counterpart
#' to the \emph{effect}-display methods \code{\link{plot.circ_gam}} /
#' \code{\link{plot.circ_lm}} (which answer \dQuote{what did the model fit?};
#' \code{circ_check} answers \dQuote{is the fit any good?}). It lays out a panel
#' grid of \code{\link{circ_resid}}-based diagnostics, prints a goodness-of-fit
#' table, and returns the statistics invisibly. Dispatches over both front doors
#' (\code{\link{circ_lm}} and \code{\link{circ_gam}}).
#'
#' @param object A fitted \code{\link{circ_lm}} or \code{\link{circ_gam}} model.
#' @param which Character vector of panel keys to draw; \code{NULL} (default)
#'   picks a sensible set for the response type, and \code{"all"} draws every
#'   panel. Keys: \code{"rose"} (rose diagram of angular residuals), \code{"obsfit"}
#'   (observed vs fitted with the wrapped calibration diagonal), \code{"residcov"}
#'   (residual vs covariate, on a circular axis when the covariate is cyclic),
#'   \code{"qq.unif"} (quantile-residual uniform Q-Q with the Watson \eqn{U^2}
#'   p-value), \code{"qq.norm"} (normal Q-Q of the deviance residuals),
#'   \code{"scaleloc"} (\eqn{\sqrt{|\,\mathrm{deviance\ residual}\,|}} vs the
#'   fitted location, a concentration-adequacy check), \code{"hist"}
#'   (deviance-residual histogram with the standard-normal reference), and
#'   \code{"cook"} (residuals vs leverage with Cook's-distance contours;
#'   \code{\link{circ_lm}} only -- see Details).
#' @param nsim Number of simulation replicates for the quantile residual when the
#'   family has no closed-form distribution function; see \code{\link{circ_resid}}.
#' @param rug Add a covariate rug to the residual-vs-covariate panel.
#' @param ... Currently ignored.
#' @return The goodness-of-fit statistics, invisibly: a list with the sample size,
#'   the residual mean direction and resultant length (circular response) or mean
#'   and standard deviation (linear response), the Watson \eqn{U^2} statistic and
#'   p-value, and the backend goodness-of-fit summary (\code{\link[mgcv]{k.check}}
#'   for \code{\link{circ_gam}}; the higher-order harmonic test, convergence flag,
#'   or least-squares fit metrics for \code{\link{circ_lm}}).
#' @details
#' The default panel set is \code{c("rose", "obsfit", "residcov", "qq.unif")} for
#' a circular response and drops \code{"rose"} for a linear response (a rose
#' diagram needs an angular residual). The headline calibration check is
#' \code{"qq.unif"}: the probability-integral-transform residuals are uniform
#' under a correct fit regardless of how the concentration varies, and the Watson
#' \eqn{U^2} test (rotation-invariant, unlike Kolmogorov-Smirnov) quantifies the
#' departure. The rose diagram should show one tight mode at zero; an off-centre
#' mode signals location bias, a multimodal one signals missed structure (too few
#' harmonics, or a smoothing basis with too small a dimension). For a circular
#' response \code{"obsfit"} carries the wrapped diagonal and its \eqn{\pm 2\pi}
#' copies, so calibration is read across the branch cut.
#'
#' The deviance-residual panels (\code{"qq.norm"}, \code{"scaleloc"},
#' \code{"hist"}) are opt-in via \code{which} (or \code{which = "all"}). They
#' exploit the deviance residual being constructed \eqn{\approx N(0, 1)}:
#' \code{"scaleloc"} in particular is the concentration-adequacy check with no
#' linear-model analogue -- a trend in \eqn{\sqrt{|\,\mathrm{deviance\
#' residual}\,|}} against the fitted location means the dispersion model is wrong.
#'
#' The influence panel \code{"cook"} is available for \code{\link{circ_lm}} fits
#' (the closed-form IRLS / least-squares leverage); for a general-family
#' \code{\link{circ_gam}} \pkg{mgcv} exposes no per-observation leverage, so the
#' panel is dropped with a message -- read basis influence from the
#' \code{\link[mgcv]{k.check}} table and the effective degrees of freedom instead.
#'
#' For the full base set of linear-response panels on a \code{"lc"} fit, call
#' \code{plot(fit$lm)} directly.
#' @references
#' Watson, G. S. (1961) Goodness-of-fit tests on a circle. \emph{Biometrika} 48,
#' 109-114.
#'
#' Wood, S. N. (2017) \emph{Generalized Additive Models: An Introduction with R}.
#' Chapman and Hall/CRC, second edition.
#' @seealso
#' \code{\link{circ_resid}} for the underlying residuals; \code{\link{plot.circ_gam}},
#' \code{\link{plot.circ_lm}} for the effect displays; \code{\link[mgcv]{gam.check}}.
#' @examples
#' set.seed(1)
#' n <- 120
#' x <- rnorm(n)
#' theta <- (1 + 2 * atan(1.5 * x) + rnorm(n) / 4) %% (2 * pi)
#' m <- circ_lm(theta ~ x, data.frame(theta, x), type = "cl")
#' circ_check(m)
#' @export
circ_check <- function(object, which = NULL, nsim = 1000L, rug = TRUE, ...) {
  parts <- .circ_resid_parts(object)
  sel   <- .circ_check_keys(which, parts$response_circular)
  keys  <- sel$keys
  if (sel$warn_rose)
    warning("panel 'rose' needs a circular response; dropping it.", call. = FALSE)

  ## leverage for the influence panel: cheap for circ_lm, unavailable for a
  ## general-family GAM (mgcv exposes no per-observation hat) -- drop it there.
  lev <- if ("cook" %in% keys) .circ_check_leverage(object) else NULL
  if ("cook" %in% keys && is.null(lev)) {
    message("circ_check: per-observation leverage is unavailable for this fit ",
            "(a general-family GAM); dropping the 'cook' panel. See the k.check ",
            "table and edf for basis influence.")
    keys <- setdiff(keys, "cook")
  }
  if (!length(keys)) stop("no panels to draw.", call. = FALSE)

  ## residuals computed once and shared. The angular and quantile residuals (and
  ## thus the Watson test) feed the table unconditionally; the deviance residual
  ## is built only when a panel that needs it is drawn.
  ang    <- parts$resid_response()
  u      <- .circ_quantile_resid(parts, as.integer(nsim), "uniform")
  watson <- .circ_watson_u2(u)
  dev    <- if (any(c("qq.norm", "scaleloc", "hist", "cook") %in% keys))
              parts$resid_deviance() else NULL
  cinfo  <- .circ_check_cov(object, parts)
  if (is.null(cinfo))                          # multi-covariate: vs the fitted LP
    cinfo <- list(name = "fitted direction", values = parts$fitted_dir,
                  circular = parts$response_circular)

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  nc <- ceiling(sqrt(length(keys))); nr <- ceiling(length(keys) / nc)
  graphics::par(mfrow = c(nr, nc))
  for (key in keys) switch(key,
    rose     = .circ_rose(ang),
    obsfit   = .circ_obsfit(parts),
    residcov = .circ_residcov(ang, cinfo, parts$response_circular, rug),
    qq.unif  = .circ_qqunif(u, watson),
    qq.norm  = .circ_qqnorm(dev),
    scaleloc = .circ_scaleloc(dev, parts$fitted_dir, parts$response_circular),
    hist     = .circ_residhist(dev),
    cook     = .circ_cook(dev, lev))

  invisible(.circ_check_table(object, parts, ang, watson))
}

## Resolve the panel keys from `which`: NULL -> the default four, "all" -> every
## panel, else the literal vector (validated). `rose` is dropped for a linear
## response (a rose needs an angular residual) -- silently for the default/all
## sets, with a warning flagged when the user listed it explicitly. Pulled out of
## circ_check() so the branching is unit-testable without a fitted model.
.circ_check_keys <- function(which, response_circular) {
  known   <- c("rose", "obsfit", "residcov", "qq.unif", "qq.norm", "scaleloc",
               "hist", "cook")
  default <- c("rose", "obsfit", "residcov", "qq.unif")
  explicit <- !is.null(which) && !identical(which, "all")
  keys <- if (is.null(which)) default
          else if (identical(which, "all")) known
          else which
  bad <- setdiff(keys, known)
  if (length(bad))
    stop("unknown panel key(s): ", paste(sQuote(bad), collapse = ", "),
         ". Available: ", paste(sQuote(known), collapse = ", "),
         " (or \"all\").", call. = FALSE)
  warn_rose <- !response_circular && explicit && "rose" %in% keys
  if (!response_circular) keys <- setdiff(keys, "rose")
  if (!length(keys)) stop("no panels to draw.", call. = FALSE)
  list(keys = keys, warn_rose = warn_rose)
}

## The single covariate of a fit (reusing the plot-method detectors), or NULL when
## there is not exactly one -- the residual-vs-covariate panel then falls back to
## the fitted direction. Carries whether that covariate is cyclic, for the axis.
.circ_check_cov <- function(object, parts) {
  if (inherits(object, "circ_gam")) {
    cov <- .circ_covariate(object)
    if (is.null(cov)) return(NULL)
    fl <- object$formula; if (inherits(fl, "formula")) fl <- list(fl)
    list(name = cov, values = as.numeric(object$model[[cov]]),
         circular = cov %in% .circ_cc_terms(fl))
  } else {
    cov <- .circ_lm_covariate(object)
    if (is.null(cov)) return(NULL)
    list(name = cov, values = as.numeric(object$frame[[cov]]),
         circular = object$type %in% c("cc", "lc"))
  }
}

## Per-observation leverage h_i and the model dimension p for the influence panel,
## or NULL when leverage is unavailable. circ_lm is cheap: lc/cc are ordinary
## least squares (hatvalues; the cc cos/sin share one design), cl reconstructs the
## converged IRLS hat. A general-family GAM returns NULL -- mgcv exposes no
## per-observation hat for them (object$hat is empty), and the principled GAM
## influence summary is the k.check edf table instead.
.circ_check_leverage <- function(object) {
  if (inherits(object, "circ_gam")) {
    if (length(object$hat)) return(list(h = as.numeric(object$hat), p = sum(object$edf)))
    return(NULL)
  }
  switch(object$type,
    lc = list(h = as.numeric(stats::hatvalues(object$lm)), p = object$lm$rank),
    cc = list(h = as.numeric(stats::hatvalues(object$cos_lm)), p = object$cos_lm$rank),
    cl = .circ_cl_leverage(object),
    NULL)
}

## The cl IRLS hat reconstructed from stored quantities (no fitter change): the
## working design and weights are rebuilt from the formula + frame, and the
## coefficient covariance (Vbeta for the mean, Vag for the concentration) is the
## inverse cross-product the fit already solved, so h = w * diag(G V G') and its
## trace recovers the design dimension p. mean/mixed use the mean design G; the
## concentration-only model uses the log-kappa design.
.circ_cl_leverage <- function(object) {
  fr <- object$frame
  if (object$model %in% c("mean", "mixed") && !is.null(object$Vbeta)) {
    X <- .circ_lm_design(object$mu_formula, fr)
    if (!ncol(X)) return(NULL)
    eta <- as.vector(X %*% object$beta)
    G   <- (2 / (1 + eta^2)) * X
    kap <- rep_len(object$kappa, nrow(X))
    h   <- (kap * A1(kap)) * rowSums((G %*% object$Vbeta) * G)
    return(list(h = as.numeric(h), p = ncol(G)))
  }
  if (object$model == "kappa" && !is.null(object$Vag)) {
    X1  <- cbind(1, .circ_lm_design(object$kappa_formula, fr))
    kap <- rep_len(object$kappa, nrow(X1))
    h   <- (kap^2 * A1prime(kap)) * rowSums((X1 %*% object$Vag) * X1)
    return(list(h = as.numeric(h), p = ncol(X1)))
  }
  NULL
}

## ---- panel drawers (base graphics only) -----------------------------------

## equal-width sector counts over (-pi, pi], the rose diagram's data. Pulled out
## of the drawer so the binning is unit-testable without a graphics device.
.circ_rose_bins <- function(ang, nbins = 24L) {
  ang  <- wrap(ang)
  brks <- seq(-pi, pi, length.out = nbins + 1L)
  list(breaks = brks,
       counts = tabulate(findInterval(ang, brks, rightmost.closed = TRUE), nbins))
}

## rose diagram of the angular residuals: equal-width sectors with radius
## proportional to sqrt(count), so sector AREA encodes frequency. One tight wedge
## at 0 (east) is a good fit; the firebrick arrow is the residual mean resultant.
.circ_rose <- function(ang, nbins = 24L, main = "angular residuals",
                       col = "steelblue") {
  b    <- .circ_rose_bins(ang, nbins)
  brks <- b$breaks; cnt <- b$counts
  rad  <- if (max(cnt) > 0) sqrt(cnt / max(cnt)) else cnt
  ang  <- wrap(ang)
  op <- graphics::par(mar = c(1, 1, 2.5, 1)); on.exit(graphics::par(op))
  graphics::plot.new(); graphics::plot.window(c(-1.2, 1.2), c(-1.2, 1.2), asp = 1)
  tt <- seq(0, 2 * pi, length.out = 100L)
  for (rr in c(0.5, 1)) graphics::lines(rr * cos(tt), rr * sin(tt),
                                        col = "gray85", lwd = 0.6)
  fill <- grDevices::adjustcolor(col, 0.5)
  for (k in seq_len(nbins)) {
    if (rad[k] <= 0) next
    a <- seq(brks[k], brks[k + 1L], length.out = 8L)
    graphics::polygon(c(0, rad[k] * cos(a), 0), c(0, rad[k] * sin(a), 0),
                      col = fill, border = col, lwd = 0.6)
  }
  labs <- expression(0, pi / 2, "" %+-% pi, -pi / 2)   # plotmath: device-safe Greek
  angs <- c(0, pi / 2, pi, -pi / 2)
  for (i in seq_along(angs))
    graphics::text(1.12 * cos(angs[i]), 1.12 * sin(angs[i]), labs[i],
                   cex = 0.8, col = "gray40")
  sb <- mean(sin(ang)); cb <- mean(cos(ang)); rbar <- sqrt(sb * sb + cb * cb)
  if (rbar > 1e-8)
    graphics::arrows(0, 0, rbar * cos(atan2(sb, cb)), rbar * sin(atan2(sb, cb)),
                     col = "firebrick", lwd = 2, length = 0.08)
  graphics::title(main = main, cex.main = 1)
}

## observed vs fitted. Circular response: the wrapped diagonal (with its +/-2pi
## copies) is perfect calibration, so off-diagonal mass shows WHERE on the circle
## the fit fails. Linear response: the ordinary obs-vs-fitted scatter with y = x.
.circ_obsfit <- function(parts) {
  y <- parts$y; f <- parts$fitted_dir
  dot <- grDevices::adjustcolor("black", 0.4)
  op <- graphics::par(mar = c(4, 4, 2.5, 1)); on.exit(graphics::par(op))
  if (parts$response_circular) {
    graphics::plot(f, y, xlim = c(-pi, pi), ylim = c(-pi, pi), asp = 1,
                   xlab = "fitted direction", ylab = "observed",
                   main = "observed vs fitted", pch = 16, cex = 0.5, col = dot)
    for (off in c(-2 * pi, 0, 2 * pi))
      graphics::abline(off, 1, col = "steelblue", lwd = 1,
                       lty = if (off == 0) 1 else 2)
  } else {
    graphics::plot(f, y, xlab = "fitted", ylab = "observed",
                   main = "observed vs fitted", pch = 16, cex = 0.5, col = dot)
    graphics::abline(0, 1, col = "steelblue", lwd = 1)
  }
}

## residual vs covariate. Leftover trend means missed structure (add a harmonic
## for cc/lc, raise the basis dimension for a gam). A cyclic covariate gets a
## circular x-axis and a wrap-aware lowess trend; a linear one an ordinary lowess.
.circ_residcov <- function(resid, cinfo, resid_circular, rug) {
  x <- cinfo$values; nm <- cinfo$name
  ylab <- if (resid_circular) "angular residual" else "residual"
  ## a cyclic covariate spans one radian period; pick the branch from the data,
  ## as plot.circ_gam() does -- c(-pi, pi) when it goes negative, else c(0, 2*pi)
  xlim <- if (cinfo$circular) {
            if (min(x, na.rm = TRUE) < 0) c(-pi, pi) else c(0, 2 * pi)
          } else range(x, na.rm = TRUE)
  op <- graphics::par(mar = c(4, 4, 2.5, 1)); on.exit(graphics::par(op))
  graphics::plot(x, resid, xlab = nm, ylab = ylab, main = "residual vs covariate",
                 pch = 16, cex = 0.5, col = grDevices::adjustcolor("black", 0.4),
                 xlim = xlim)
  graphics::abline(h = 0, col = "steelblue")
  ok <- is.finite(x) & is.finite(resid)
  if (sum(ok) > 5L) {
    if (cinfo$circular) {                       # tile +/-2pi so the smooth wraps
      xt <- c(x[ok] - 2 * pi, x[ok], x[ok] + 2 * pi); rt <- rep(resid[ok], 3L)
      o <- order(xt); lw <- stats::lowess(xt[o], rt[o], f = 1 / 3)
      keep <- lw$x >= -pi & lw$x <= pi
      graphics::lines(lw$x[keep], lw$y[keep], col = "firebrick", lwd = 1.5)
    } else {
      o <- order(x[ok])
      graphics::lines(stats::lowess(x[ok][o], resid[ok][o]),
                      col = "firebrick", lwd = 1.5)
    }
  }
  if (rug) graphics::rug(x, col = grDevices::adjustcolor("black", 0.3))
}

## quantile-residual uniform Q-Q: ordered PIT residuals against (i - 1/2)/n on the
## unit square. Points on the line == uniform == calibrated; the subtitle carries
## the Watson U^2 statistic and its uniformity p-value.
.circ_qqunif <- function(u, watson) {
  n <- length(u); us <- sort(u); pp <- (seq_len(n) - 0.5) / n
  op <- graphics::par(mar = c(4, 4, 2.5, 1)); on.exit(graphics::par(op))
  graphics::plot(pp, us, xlim = c(0, 1), ylim = c(0, 1), asp = 1,
                 xlab = "theoretical U(0, 1) quantile", ylab = "ordered PIT residual",
                 main = "quantile-residual Q-Q", pch = 16, cex = 0.5,
                 col = grDevices::adjustcolor("black", 0.5))
  graphics::abline(0, 1, col = "steelblue")
  ## the Watson U^2 stat sits inside the empty top-left corner (above the
  ## diagonal), not in the title margin where it would collide in a dense grid
  graphics::legend("topleft", bty = "n", cex = 0.75,
                   legend = sprintf("Watson U2 = %.3f\np = %.3f",
                                    watson$stat, watson$p))
}

## ---- Tier-2 panels: all keyed off the deviance residual (~ N(0,1)) ---------

## normal Q-Q of the deviance residuals -- the linear-model diagnostic that holds
## on the circle because the deviance residual is constructed ~ N(0, 1). Points on
## the qqline == the per-observation deviance is calibrated.
.circ_qqnorm <- function(dev) {
  op <- graphics::par(mar = c(4, 4, 2.5, 1)); on.exit(graphics::par(op))
  stats::qqnorm(dev, main = "deviance-residual normal Q-Q",
                xlab = "theoretical N(0, 1) quantile", ylab = "deviance residual",
                pch = 16, cex = 0.5, col = grDevices::adjustcolor("black", 0.5))
  stats::qqline(dev, col = "steelblue")
}

## scale-location: sqrt|deviance residual| vs the fitted location. A trend is the
## circular failure mode plot.lm never had -- the dispersion/concentration model is
## wrong (raise circ_lm "mean" -> "mixed", or add kappa ~ s(x) in circ_gam). The
## trend line is wrap-aware when the fitted location is circular.
.circ_scaleloc <- function(dev, fitted, resp_circular) {
  rs <- sqrt(abs(dev))
  xlab <- if (resp_circular) "fitted direction" else "fitted"
  xlim <- if (resp_circular) c(-pi, pi) else range(fitted, na.rm = TRUE)
  op <- graphics::par(mar = c(4, 4, 2.5, 1)); on.exit(graphics::par(op))
  graphics::plot(fitted, rs, xlab = xlab, xlim = xlim,
                 ylab = expression(sqrt(abs("deviance residual"))),
                 main = "scale-location", pch = 16, cex = 0.5,
                 col = grDevices::adjustcolor("black", 0.4))
  ok <- is.finite(fitted) & is.finite(rs)
  if (sum(ok) > 5L) {
    if (resp_circular) {                        # periodic smooth over the circle
      ft <- c(fitted[ok] - 2 * pi, fitted[ok], fitted[ok] + 2 * pi)
      rt <- rep(rs[ok], 3L); o <- order(ft); lw <- stats::lowess(ft[o], rt[o], f = 1 / 3)
      keep <- lw$x >= -pi & lw$x <= pi
      graphics::lines(lw$x[keep], lw$y[keep], col = "firebrick", lwd = 1.5)
    } else {
      o <- order(fitted[ok])
      graphics::lines(stats::lowess(fitted[ok][o], rs[ok][o]),
                      col = "firebrick", lwd = 1.5)
    }
  }
}

## histogram of the deviance residuals (the gam.check fourth panel) with the
## standard-normal reference density overlaid -- the shape it should match.
.circ_residhist <- function(dev) {
  op <- graphics::par(mar = c(4, 4, 2.5, 1)); on.exit(graphics::par(op))
  h <- graphics::hist(dev, breaks = "FD", plot = FALSE)
  graphics::hist(dev, breaks = "FD", freq = FALSE, border = "white",
                 col = grDevices::adjustcolor("steelblue", 0.4),
                 xlab = "deviance residual", main = "deviance-residual histogram",
                 ylim = c(0, max(h$density, stats::dnorm(0)) * 1.05))
  xs <- seq(min(dev), max(dev), length.out = 200L)
  graphics::lines(xs, stats::dnorm(xs), col = "firebrick", lwd = 1.5)
}

## ---- Tier-3 influence panel -----------------------------------------------

## residuals vs leverage with Cook's-distance contours -- a high-leverage point
## with a large residual swings the (closed-form) fit. The standardized deviance
## residual is dev / sqrt(1 - h); Cook's D = rstd^2 * h / (p (1 - h)); points past
## the 4/n rule of thumb are labelled by index.
.circ_cook <- function(dev, lev) {
  h <- lev$h; p <- lev$p
  rstd <- dev / sqrt(pmax(1 - h, 1e-8))
  D <- rstd^2 * h / (p * pmax(1 - h, 1e-8))
  op <- graphics::par(mar = c(4, 4, 2.5, 1)); on.exit(graphics::par(op))
  graphics::plot(h, rstd, xlab = "leverage", ylab = "std. deviance residual",
                 main = "residuals vs leverage", pch = 16, cex = 0.5,
                 col = grDevices::adjustcolor("black", 0.45),
                 xlim = c(0, max(h, na.rm = TRUE) * 1.05))
  graphics::abline(h = 0, col = "steelblue", lty = 2)
  pos <- h > 1e-6
  if (any(pos)) {                              # Cook's D in {0.5, 1}: rstd = +/- sqrt(D p (1-h)/h)
    hs <- seq(min(h[pos]), max(h), length.out = 100L)
    for (d in c(0.5, 1)) {
      rr <- sqrt(d * p * (1 - hs) / hs)
      graphics::lines(hs, rr, col = "firebrick", lty = 2, lwd = 0.8)
      graphics::lines(hs, -rr, col = "firebrick", lty = 2, lwd = 0.8)
    }
  }
  hi <- which(D > 4 / length(h))
  if (length(hi))
    graphics::text(h[hi], rstd[hi], labels = hi, cex = 0.6, pos = 3, col = "firebrick")
}

## ---- goodness-of-fit table (the gam.check analogue) -----------------------

## Print the residual summary + the backend goodness-of-fit, and return the
## statistics invisibly. Residual location: a mean direction and resultant length
## for a circular response (mean direction ~ 0 and small Rbar == unbiased,
## concentrated residuals), a mean and sd for a linear response.
.circ_check_table <- function(object, parts, ang, watson) {
  cat("\ncirc_check:", parts$family_label, "  n =", parts$n, "\n")
  out <- list(n = parts$n, family = parts$family_label, watson = watson)
  if (parts$response_circular) {
    sb <- mean(sin(ang)); cb <- mean(cos(ang))
    out$resid_mean_direction <- atan2(sb, cb); out$Rbar <- sqrt(sb * sb + cb * cb)
    cat(sprintf("  residual mean direction = %+.4f rad   resultant length = %.4f\n",
                out$resid_mean_direction, out$Rbar))
  } else {
    out$resid_mean <- mean(ang); out$resid_sd <- stats::sd(ang)
    cat(sprintf("  residual mean = %+.4f   sd = %.4f\n", out$resid_mean, out$resid_sd))
  }
  cat(sprintf("  Watson U2 (PIT uniformity) = %.4f   p = %.4f\n",
              watson$stat, watson$p))
  if (inherits(object, "circ_gam")) {
    kc <- tryCatch(mgcv::k.check(object), error = function(e) NULL)
    if (!is.null(kc)) {
      cat("\n  k.check (basis dimension adequacy):\n")
      print(round(kc, 4)); out$k.check <- kc
    }
  } else if (object$type == "cc") {
    cat(sprintf("\n  higher-order harmonic test: p(cos) = %.4f   p(sin) = %.4f\n",
                object$p_values[["cos"]], object$p_values[["sin"]]))
    out$p_values <- object$p_values
  } else if (object$type == "lc") {
    cat(sprintf("\n  R-squared = %.4f   residual sigma = %.4f\n",
                object$r_squared, object$sigma))
    out$r_squared <- object$r_squared; out$sigma <- object$sigma
  } else {                                       # cl
    krange <- if (length(object$kappa) == 1L) sprintf("= %.3f", object$kappa)
              else sprintf("in [%.3f, %.3f]", min(object$kappa), max(object$kappa))
    cat(sprintf("\n  converged = %s   kappa %s\n", isTRUE(object$converged), krange))
    out$converged <- isTRUE(object$converged)
  }
  invisible(out)
}
