#' Plot a fitted circular mixture
#'
#' Clustered views of a fitted mixture. The flat view colours the observations by
#' their MAP cluster, with -- for a regression cell -- each component's fitted
#' location curve over the single covariate, or -- for a density cell
#' (\code{theta ~ 1}) -- the per-cluster spread of the response with each
#' component's fitted mean direction. The geometry view draws those per-component
#' curves on the leg's natural 3D surface (cylinder / torus / can).
#'
#' @param x A fitted \code{\link{circ_mix}} object.
#' @param view Which view to draw. \code{"flat"} (default) draws the flat panel.
#'   \code{"geometry"} draws every component's fitted location curve on the leg's
#'   natural surface -- the cylinder (c~l), torus (c~c) or can (l~c), sharing
#'   \code{\link{plot.circ_gam}}'s surface canvas; \code{"both"} places the two
#'   side by side. A joint (product) component has only the flat torus-square
#'   projection, so it draws that whatever the \code{view}.
#' @param n Number of grid points for each component's fitted curve.
#' @param se If \code{TRUE}, band each component's fitted curve where the view
#'   supports it: \eqn{\pm} the component law's circular standard deviation (its
#'   predictive angular spread) for a circular response, a pointwise 2-SE interval
#'   for a linear response.
#' @param ... Further graphical arguments (currently unused).
#' @return \code{x}, invisibly. Called for the plot it draws.
#' @examples
#' library(mgcv)
#' set.seed(1); n <- 400
#' z  <- sample.int(2L, n, replace = TRUE)
#' x  <- runif(n, -1, 1)
#' mu <- 2 * atan(c(1, -1)[z] + c(2, -2)[z] * x)
#' y  <- vmlss()$rd(cbind(mu, rep(6, n)), rep(1, n), 1)
#' \donttest{
#' m  <- circ_mix(y ~ x, data = data.frame(y, x), family = vmlss(), K = 2)
#' plot(m)
#' }
#' @seealso
#' \code{\link{circ_mix}}, \code{\link{plot.circ_gam}}
#' @export
plot.circ_mix <- function(x, view = c("flat", "geometry", "both"), n = 200,
                          se = TRUE, ...) {
  view <- match.arg(view)

  ## joint (product) component: only the flat torus-square projection exists
  if (identical(x$geometry, "joint")) {
    if (view != "flat")
      message("plot.circ_mix: the joint geometry surface is not drawn yet; ",
              "showing the flat torus-square.")
    return(.circ_mix_plot_joint(x))
  }

  g       <- .circ_geometry(x$components[[1L]]$fit)
  surface <- switch(if (is.na(g$kind)) "" else g$kind,
                    cl = "cylinder", cc = "torus", lc = "can", NA_character_)
  if (view != "flat" && (is.null(g$cov) || is.na(surface))) {
    message("plot.circ_mix: no surface for this fit (geometry needs exactly one ",
            "covariate); drawing the flat view.")
    view <- "flat"
  }

  if (view == "flat") return(.circ_mix_flat(x, g, n, se))

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  if (view == "both") graphics::par(mfrow = c(1, 2))
  .circ_mix_geometry(x, g, surface, n, se)
  ## the flat half of "both" is exactly the view = "flat" panel (the clustered
  ## scatter/stripchart); reset the geometry's tight margins so its axes have room.
  if (view == "both") {
    graphics::par(mar = c(4.1, 4.1, 2.6, 1.1))
    .circ_mix_flat(x, g, n, se)
  }
  invisible(x)
}

## per-component band for the LOCATION on the response scale, matching
## plot.circ_gam. For a circular response it is the half-width `csd` of the +/-
## circular-SD predictive band -- the component law's own spread sqrt(-2 log R)
## by .circ_sd_quad, which covers every family including pnlss's derived
## direction. For a linear response it is a delta-method 2-SE lo--hi on the mean.
## NULL on a failed predict.
.circ_mix_loc_band <- function(cp, nd, resp_circular) {
  if (resp_circular) {
    fm <- tryCatch(cmp_predict(cp, nd, "response"), error = function(e) NULL)
    if (is.null(fm)) return(NULL)
    return(list(csd = .circ_sd_quad(cp$fit$family, fm)))
  }
  fam <- cp$fit$family
  pr <- tryCatch(stats::predict(cp$fit, newdata = nd, type = "link", se.fit = TRUE),
                 error = function(e) NULL)
  if (is.null(pr)) return(NULL)
  eta <- as.matrix(pr$fit); sef <- as.matrix(pr$se.fit)
  li <- fam$linfo[[1L]]$linkinv
  list(lo = li(eta[, 1L] - 2 * sef[, 1L]), hi = li(eta[, 1L] + 2 * sef[, 1L]))
}

## --- the flat panel: density stripchart or regression curves, MAP-coloured ------
.circ_mix_flat <- function(x, g, n, se = TRUE) {
  cov  <- g$cov
  cols <- .circ_mix_palette(x$K); z <- x$cluster[x$unit$index]   # per-row (subject -> its rows)
  dat  <- x$components[[1L]]$fit$model
  yobs <- as.numeric(dat[[x$response]])
  yw   <- if (g$resp_circular) atan2(sin(yobs), cos(yobs)) else yobs

  ## density cell (no covariate): clusters of the response alone
  if (is.null(cov)) {
    grp <- sort(unique(z))
    graphics::stripchart(split(yw, z), method = "jitter", pch = 16, cex = 0.6,
                         col = cols[grp], xlab = x$response, las = 1,
                         main = sprintf("circ_mix: %d clusters", x$Gtilde))
    if (g$resp_circular) {
      mu <- vapply(x$components, function(cp)
        atan2(sin(cmp_predict(cp, dat[1, , drop = FALSE], "response")[1, 1]),
              cos(cmp_predict(cp, dat[1, , drop = FALSE], "response")[1, 1])),
        numeric(1))
      graphics::abline(v = mu, col = cols, lwd = 2, lty = 2)
    }
    return(invisible(x))
  }

  ## regression cell: response vs covariate, per-component location
  xv   <- as.numeric(dat[[cov]])
  rng  <- range(xv, na.rm = TRUE)
  if (g$cov_circular) rng <- if (rng[1] < 0) c(-pi, pi) else c(0, 2 * pi)
  grid <- seq(rng[1], rng[2], length.out = n)
  nd   <- stats::setNames(data.frame(grid), cov)
  ylim <- if (g$resp_circular) c(-pi, pi) else range(yw, na.rm = TRUE)
  graphics::plot(xv, yw, col = grDevices::adjustcolor(cols[z], 0.55),
                 pch = 16, cex = 0.6, xlab = cov, ylab = x$response,
                 ylim = ylim, main = sprintf("circ_mix: %d components", x$K))
  yk_raw <- lapply(seq_len(x$K), function(k)
    cmp_predict(x$components[[k]], nd, "response")[, 1])
  ## per-component band: +/- circular-SD wrapped fill (circular response) or a
  ## 2-SE shadow (linear response), matching plot.circ_gam + the geometry ribbon
  if (se) for (k in seq_len(x$K)) {
    bd <- .circ_mix_loc_band(x$components[[k]], nd, g$resp_circular); if (is.null(bd)) next
    col <- grDevices::adjustcolor(cols[k], 0.25)
    if (g$resp_circular) .circ_band_fill_circular(grid, yk_raw[[k]], bd$csd, col)
    else                 .circ_band_fill(grid, bd$lo, bd$hi, col)
  }
  for (k in seq_len(x$K)) {                    # curve tiled to wrap with its band
    if (g$resp_circular) .circ_lines_circular(grid, yk_raw[[k]], cols[k], 2.5)
    else                 graphics::lines(grid, yk_raw[[k]], col = cols[k], lwd = 2.5)
  }
  invisible(x)
}

## --- the geometry panel: every component's location curve on the shared surface --
.circ_mix_geometry <- function(x, g, surface, n, se = TRUE) {
  cov  <- g$cov
  dat  <- x$components[[1L]]$fit$model
  cols <- .circ_mix_palette(x$K); z <- x$cluster[x$unit$index]   # per-row (subject -> its rows)
  xv   <- as.numeric(dat[[cov]]); yobs <- as.numeric(dat[[x$response]])
  rng  <- range(xv, na.rm = TRUE)
  if (g$cov_circular) rng <- if (rng[1] < 0) c(-pi, pi) else c(0, 2 * pi)
  grid <- seq(rng[1], rng[2], length.out = n)
  nd   <- stats::setNames(data.frame(grid), cov)
  wrapc <- function(v) atan2(sin(v), cos(v))

  ## per-component fitted location: raw (for the band centre + ribbon, which wrap
  ## via xyz) and wrapped (the curve reads continuously on the surface)
  yk_raw <- lapply(seq_len(x$K), function(k)
    cmp_predict(x$components[[k]], nd, "response")[, 1])
  yk   <- lapply(yk_raw, function(v) if (g$resp_circular) wrapc(v) else v)
  vobs  <- if (g$resp_circular) wrapc(yobs) else yobs
  bds   <- if (se) lapply(x$components, .circ_mix_loc_band, nd = nd,
                          resp_circular = g$resp_circular) else vector("list", x$K)
  bspan <- unlist(lapply(bds, function(b) c(b$lo, b$hi)))   # NULL for a circular csd band
  yspan <- if (surface == "can") c(unlist(yk), vobs, bspan) else NULL

  graphics::par(mar = c(1.0, 0.3, 2.2, 0.3))
  s <- .circ_surface(surface, grid, xv, yspan)
  ## per-component band ribbon, beneath the points and curves: +/- circular SD
  ## about the raw location (circular) or the 2-SE lo--hi (linear)
  for (k in seq_len(x$K)) {
    bd <- bds[[k]]; if (is.null(bd)) next
    lo <- if (g$resp_circular) yk_raw[[k]] - bd$csd else bd$lo
    hi <- if (g$resp_circular) yk_raw[[k]] + bd$csd else bd$hi
    .circ_fade_ribbon(s$to_u(grid), s$to_v(lo), s$to_v(hi), s$xyz, s$pm, cols[k])
  }
  ## observations, depth-faded, one .circ_fade_points() call per MAP colour
  for (k in seq_len(x$K)) {
    idx <- z == k; if (!any(idx)) next
    w <- s$xyz(s$to_u(xv[idx]), s$to_v(vobs[idx]))
    .circ_fade_points(w$x, w$y, w$z, s$pm, cols[k], cex = 0.5)
  }
  for (k in seq_len(x$K)) {
    w <- s$xyz(s$to_u(grid), s$to_v(yk[[k]]))
    .circ_fade_curve(w$x, w$y, w$z, s$pm, cols[k], lwd = 3)
  }
  graphics::title(main = sprintf("%s \u00b7 circ_mix (%d components)", surface, x$K),
                  line = 0.5, cex.main = 1.15)
  invisible(x)
}

## the joint flat view: each observation on the response square, coloured by MAP
## cluster, with each component's weighted (parent, conditional) centroid marked.
## The honest torus surface is a later view = "geometry"; this is the flat
## projection (cf. the rama 05-joint-K4 figure). The parent response (a covariate
## in another factor) goes on x, the conditional response on y.
.circ_mix_plot_joint <- function(x) {
  resps <- .circ_mix_responses(x$formula)            # c(conditional, parent) order
  circ  <- .circ_mix_resp_circular(x$family)
  fits  <- x$components[[1L]]$fits
  vals  <- lapply(seq_along(resps), function(j)
    as.numeric(fits[[j]]$model[[resps[j]]]))
  wrap  <- function(v) if (circ) atan2(sin(v), cos(v)) else v
  cols  <- .circ_mix_palette(x$K); z <- x$cluster[x$unit$index]   # per-row
  graphics::plot(wrap(vals[[2L]]), wrap(vals[[1L]]),
                 col = grDevices::adjustcolor(cols[z], 0.6), pch = 16, cex = 0.6,
                 xlab = resps[2L], ylab = resps[1L],
                 main = sprintf("circ_mix joint: %d torus clusters", x$Gtilde))
  for (k in seq_len(x$K)) {
    w  <- x$gamma[, k]
    cx <- if (circ) atan2(sum(w * sin(vals[[2L]])), sum(w * cos(vals[[2L]])))
          else stats::weighted.mean(vals[[2L]], w)
    cy <- if (circ) atan2(sum(w * sin(vals[[1L]])), sum(w * cos(vals[[1L]])))
          else stats::weighted.mean(vals[[1L]], w)
    graphics::points(wrap(cx), wrap(cy), pch = 21, bg = cols[k], col = "black",
                     cex = 2, lwd = 2)
  }
  invisible(x)
}

## a fixed, version-stable categorical palette (recycled past 10 components)
.circ_mix_palette <- function(K) {
  base <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
            "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f")
  rep(base, length.out = K)
}
