#' Plot a classical circular regression fit
#'
#' The default plot method for \code{\link{circ_lm}} fits: the closed-form sibling
#' of \code{\link{plot.circ_gam}}, drawing the same three regression geometries --
#' circular--linear (cylinder), circular--circular (torus) and linear--circular
#' (the upright can) -- from the fit's \code{type}. A multi-covariate \code{"cl"}
#' fit has no single covariate axis and defers with a message to
#' \code{\link{coef}} / \code{\link{predict}} / \code{summary}.
#'
#' @param x A fitted \code{\link{circ_lm}} model.
#' @param view \code{"flat"} (default) draws one panel per modelled parameter on the
#'   response scale against the covariate. A circular location carries a band of
#'   \eqn{\pm} one circular standard deviation of the fitted law (the von Mises
#'   spread \eqn{\sqrt{-2\log A_1(\kappa)}} of the fitted concentration); the
#'   concentration panel keeps a delta-method 2-SE band. A circular location is
#'   broken at the \eqn{\pm\pi} branch jump and the observed responses are overlaid.
#'   \code{type = "cl"} draws the mean direction and the concentration \eqn{\kappa};
#'   \code{type = "cc"} / \code{"lc"} draw the single fitted location. \code{"geometry"}
#'   draws the fitted location curve on its natural surface -- a \strong{cylinder} for
#'   circular--linear (\code{"cl"}: the response angle wraps the tube, the linear
#'   covariate runs along the axis), a \strong{torus} for circular--circular
#'   (\code{"cc"}: the covariate around the ring), or an upright \strong{can} for
#'   linear--circular (\code{"lc"}: the cyclic covariate wraps the ring, the linear
#'   response is the height). \code{"both"} places the geometry canvas beside the
#'   full set of flat panels (\code{"cl"}: the mean direction and \eqn{\kappa};
#'   \code{"cc"}/\code{"lc"}: the single location) -- exactly the panels \code{"flat"}
#'   draws, so the two views never disagree.
#' @param n Number of grid points along the covariate.
#' @param se Draw the uncertainty band (filled shadow on the flat panels,
#'   translucent ribbon on the geometry surface): \eqn{\pm} the circular standard
#'   deviation for a circular location, a 2-SE interval for the concentration.
#' @param pages If 1, lay the flat panels out on a single page.
#' @param rug Add a covariate rug to the location panel.
#' @param ... Currently ignored.
#' @return The fitted model, invisibly.
#' @details
#' The geometry canvas is base-graphics only (\code{\link[graphics]{persp}} +
#' \code{\link[grDevices]{trans3d}}) and shares the surface, panel and band helpers
#' with \code{\link{plot.circ_gam}}, so a \code{circ_lm} leg renders in the same
#' idiom as the matching \code{circ_gam} leg. A circular location is banded by
#' \eqn{\pm} the circular standard deviation of the fitted law -- the von Mises
#' spread \eqn{\sqrt{-2\log A_1(\kappa)}} for the per-point concentration of
#' \code{"cl"} and the residual concentration of \code{"cc"} -- the predictive
#' angular spread, not a confidence interval of the mean. The remaining bands are
#' the usual intervals: the \code{"cl"} concentration on the log scale through
#' \eqn{Z\,V Z'}, the \code{"lc"} linear response the ordinary least-squares
#' prediction band.
#' @examples
#' set.seed(1)
#' n <- 80
#' x <- rnorm(n)
#' theta <- (1 + 2 * atan(1.5 * x) + rnorm(n) / 4) %% (2 * pi)
#' dat <- data.frame(theta = theta, x = x)
#'
#' ## cl: the mean direction on the cylinder, with the kappa panel beside it
#' m <- circ_lm(theta ~ x, dat, type = "cl")
#' plot(m)                      # flat: mu (circular) and kappa
#' plot(m, view = "geometry")   # the fitted angle on the cylinder
#'
#' ## cc / lc: a single location on the torus / can
#' phi <- runif(n, 0, 2 * pi)
#' dcc <- data.frame(psi = (phi / 2 + rnorm(n) / 5) %% (2 * pi), phi = phi)
#' plot(circ_lm(psi ~ phi, dcc, type = "cc"), view = "both")
#'
#' dlc <- data.frame(y = 5 + 2 * cos(phi) + rnorm(n) / 2, phi = phi)
#' plot(circ_lm(y ~ phi, dlc, type = "lc"), view = "geometry")
#' @seealso
#' \code{\link{circ_lm}}, \code{\link{plot.circ_gam}}
#' @export
plot.circ_lm <- function(x, view = c("flat", "geometry", "both"),
                         n = 200, se = TRUE, pages = 1, rug = TRUE, ...) {
  view <- match.arg(view)
  surface <- switch(x$type, cl = "cylinder", cc = "torus", lc = "can")
  resp_circular <- x$type %in% c("cl", "cc")
  cov_circular  <- x$type %in% c("cc", "lc")

  cov <- .circ_lm_covariate(x)
  if (is.null(cov)) {
    message("plot.circ_lm: the geometry/curve views need exactly one covariate; ",
            "this is a multi-covariate fit. Use coef()/predict()/summary().")
    return(invisible(x))
  }

  xv  <- as.numeric(x$frame[[cov]])
  rng <- range(xv, na.rm = TRUE)
  if (cov_circular) rng <- if (rng[1] < 0) c(-pi, pi) else c(0, 2 * pi)
  grid <- seq(rng[1], rng[2], length.out = n)
  nd   <- stats::setNames(data.frame(grid), cov)

  yobs <- as.numeric(x$frame[[x$response]])
  if (resp_circular) yobs <- atan2(sin(yobs), cos(yobs))

  panels   <- .circ_lm_panels(x, grid, nd)
  surf_idx <- 1L                      # the location is always the first panel

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))

  if (view == "flat") {
    if (isTRUE(pages == 1) && length(panels) > 1L) {
      nc <- ceiling(sqrt(length(panels))); nr <- ceiling(length(panels) / nc)
      graphics::par(mfrow = c(nr, nc))
    }
    for (i in seq_along(panels))
      .circ_panel(grid, panels[[i]], cov, se,
                  xobs = if (isTRUE(i == surf_idx)) xv else NULL,
                  yobs = if (isTRUE(i == surf_idx)) yobs else NULL, rug = rug)
    return(invisible(x))
  }

  sp  <- panels[[surf_idx]]
  glo <- if (!se) NULL else if (!is.null(sp$csd)) sp$mid - sp$csd else sp$lo
  ghi <- if (!se) NULL else if (!is.null(sp$csd)) sp$mid + sp$csd else sp$hi
  if (view == "both") .circ_both_layout(length(panels))
  .circ_geometry_panel(grid, sp$mid, xv, yobs, surface, lo = glo, hi = ghi)
  ## "both" is the geometry canvas plus the FULL flat panel set ("cl": the mean
  ## direction and kappa; "cc"/"lc": the single location) -- the panels "flat"
  ## draws, so the two views never disagree. (.circ_both_layout lives in plot-circ_gam.R.)
  if (view == "both") {
    graphics::par(mar = c(4.1, 4.1, 2.6, 1.1))      # undo the geometry's tight margins
    for (i in seq_along(panels))
      .circ_panel(grid, panels[[i]], cov, se,
                  xobs = if (isTRUE(i == surf_idx)) xv else NULL,
                  yobs = if (isTRUE(i == surf_idx)) yobs else NULL, rug = rug)
  }
  invisible(x)
}

## the single covariate of a circ_lm fit, or NULL if not exactly one. "cc"/"lc"
## always carry one angular predictor; "cl" pools the mean and log-kappa formulas.
.circ_lm_covariate <- function(x) {
  if (x$type %in% c("cc", "lc")) return(x$var)
  v <- unique(c(all.vars(x$mu_formula[[3L]]), all.vars(x$kappa_formula[[2L]])))
  if (length(v) == 1L) v else NULL
}

## the response-scale panels for a fit, each list(name, circular, mid, ...) in the
## shape .circ_panel()/.circ_geometry_panel() expect. A circular location carries
## `csd`, the half-width of its +/- circular-SD predictive band (the spread of a
## von Mises with the fitted concentration, sqrt(-2 log A1(kappa))); a linear
## quantity (kappa, the lc response) carries a delta-method/OLS lo--hi. The `se`
## flag gates only the display. Circular locations are centred to (-pi, pi]
## (atan2) so the flat break and the c(-pi, pi) ylim line up.
.circ_lm_panels <- function(x, grid, nd) {
  if (x$type == "lc") {
    p <- .circ_lm_lc_loc(x, grid)
    return(list(list(name = x$response, circular = FALSE,
                     mid = p$mid, lo = p$lo, hi = p$hi)))
  }
  if (x$type == "cc") {
    csd <- .circ_sd_from_R(A1(rep_len(as.numeric(x$kappa), length(grid))))
    return(list(list(name = x$response, circular = TRUE,
                     mid = .circ_lm_cc_loc(x, grid), csd = csd)))
  }
  ka <- .circ_lm_cl_kappa(x, nd)
  list(list(name = x$response, circular = TRUE, mid = .circ_lm_cl_mu(x, nd),
            csd = .circ_sd_from_R(A1(ka$mid))),
       list(name = "kappa", circular = FALSE, mid = ka$mid, lo = ka$lo, hi = ka$hi))
}

## ---- cl: Fisher-Lee mean direction and concentration ----------------------
## mu_i = mu0 + 2*atan(X beta), wrapped to (-pi, pi]. The band on it is +/- the
## circular SD of the fitted concentration kappa_i (.circ_lm_panels), so only the
## direction itself is returned here.
.circ_lm_cl_mu <- function(x, nd) {
  Xn  <- .circ_lm_design(x$mu_formula, nd)
  eta <- if (ncol(Xn)) as.vector(Xn %*% x$beta) else rep(0, nrow(nd))
  atan2(sin(x$mu + 2 * atan(eta)), cos(x$mu + 2 * atan(eta)))
}

## log(kappa_i) = alpha + Z gamma. Band on the log scale through Vag, then exp().
## Constant kappa (the mean-only model, no kappa design) is kappa +/- 2 se_kappa.
.circ_lm_cl_kappa <- function(x, nd) {
  Zn <- .circ_lm_design(x$kappa_formula, nd)
  if (!ncol(Zn)) {
    mid <- rep(x$kappa[1L], nrow(nd))
    if (is.null(x$se_kappa)) return(list(mid = mid, lo = NULL, hi = NULL))
    return(list(mid = mid, lo = mid - 2 * x$se_kappa, hi = mid + 2 * x$se_kappa))
  }
  eta <- x$alpha + as.vector(Zn %*% x$gamma)
  if (is.null(x$Vag)) return(list(mid = exp(eta), lo = NULL, hi = NULL))
  Z1 <- cbind(1, Zn); se_eta <- sqrt(pmax(rowSums((Z1 %*% x$Vag) * Z1), 0))
  list(mid = exp(eta), lo = exp(eta - 2 * se_eta), hi = exp(eta + 2 * se_eta))
}

## ---- cc: harmonic circular-circular location ------------------------------
## direction = atan2(sin_fit, cos_fit) of the two harmonic OLS predictions. The
## band on it is +/- the circular SD of the fitted residual concentration
## (.circ_lm_panels), so only the direction is returned here.
.circ_lm_cc_loc <- function(x, grid) {
  H  <- .circ_lm_harmonic(grid %% (2 * pi), x$order)
  cf <- as.numeric(stats::predict(x$cos_lm, H))
  sf <- as.numeric(stats::predict(x$sin_lm, H))
  atan2(sf, cf)
}

## ---- lc: harmonic linear-circular location -------------------------------
## the linear mean over the cyclic covariate; the OLS prediction band directly.
.circ_lm_lc_loc <- function(x, grid) {
  H  <- .circ_lm_harmonic(grid %% (2 * pi), x$order)
  pr <- stats::predict(x$lm, H, se.fit = TRUE)
  fit <- as.numeric(pr$fit); se <- as.numeric(pr$se.fit)
  list(mid = fit, lo = fit - 2 * se, hi = fit + 2 * se)
}
