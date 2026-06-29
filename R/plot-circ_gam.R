#' Plot a circular-response GAM fit
#'
#' The default plot method for \code{\link{circ_gam}} fits: geometry-aware curves
#' that \code{\link[mgcv]{plot.gam}} does not provide, covering the three
#' regression geometries -- circular--linear (cylinder), circular--circular
#' (torus) and linear--circular (the upright can). For a fit with several
#' covariates, no covariate, or an ordinary (non-circlss) family it falls back to
#' mgcv's per-smooth term plots.
#'
#' @param x A fitted \code{\link{circ_gam}} model.
#' @param view \code{"flat"} (default) draws one panel per distribution parameter
#'   on the response scale against the covariate. A circular location or direction
#'   carries a band of \eqn{\pm} one circular standard deviation
#'   \eqn{\sqrt{-2\log R}} of the fitted law (its predictive angular spread, with
#'   \eqn{R} the mean resultant length), wrapped around the \eqn{\pm\pi} branch
#'   cut; a non-circular parameter (the concentration, ...) keeps a delta-method
#'   2-SE band. The location panel is broken at the \eqn{\pm\pi} branch jump and
#'   the observed responses are overlaid. \code{"geometry"} draws the fitted curve
#'   on its natural
#'   surface -- a \strong{cylinder} for circular--linear (the response angle wraps
#'   the tube, the covariate runs along the axis), a \strong{torus} for
#'   circular--circular (covariate around the ring), or an upright \strong{can} for
#'   linear--circular (a linear response over a cyclic covariate: the covariate wraps
#'   the ring, the response is the height) -- chosen from whether the response is
#'   circular and whether the covariate is cyclic. \code{"both"} places the geometry
#'   canvas beside the full set of flat parameter panels -- exactly the panels
#'   \code{"flat"} draws, so the geometry and flat views never disagree.
#' @param n Number of grid points along the covariate.
#' @param se Draw the uncertainty band: \eqn{\pm} the circular standard deviation
#'   for a circular location/direction, a 2-SE interval for a non-circular
#'   parameter.
#' @param pages If 1, lay the flat panels out on a single page.
#' @param rug Add a covariate rug to the location panel.
#' @param ... Passed to \code{\link[mgcv]{plot.gam}} in the fallback cases.
#' @return The fitted model, invisibly.
#' @details
#' The geometry canvas is base-graphics only (\code{\link[graphics]{persp}} +
#' \code{\link[grDevices]{trans3d}}) and is chosen descriptively from the
#' covariate's basis; it is never a fitting input. For \code{\link{pnlss}} the
#' curve drawn is the derived mean direction \eqn{\mathrm{atan2}(\mu_2, \mu_1)}.
#' @seealso
#' \code{\link{circ_gam}}, \code{\link[mgcv]{plot.gam}}
#' @export
plot.circ_gam <- function(x, view = c("flat", "geometry", "both"),
                          n = 200, se = TRUE, pages = 1, rug = TRUE, ...) {
  view <- match.arg(view)
  fam <- x$family
  pnm <- fam$param_names
  if (is.null(pnm)) return(.circ_plot_gam(x, pages, ...))

  g   <- .circ_geometry(x)
  cov <- g$cov
  if (is.null(cov)) {
    message("plot.circ_gam: not a single-covariate fit; drawing mgcv term ",
            "plots. Use mgcv::plot.gam(.) directly for full control.")
    return(.circ_plot_gam(x, pages, ...))
  }

  xv  <- as.numeric(x$model[[cov]])
  rng <- range(xv, na.rm = TRUE)
  if (g$cov_circular) rng <- if (rng[1] < 0) c(-pi, pi) else c(0, 2 * pi)
  grid <- seq(rng[1], rng[2], length.out = n)
  nd   <- stats::setNames(data.frame(grid), cov)

  pr   <- stats::predict(x, newdata = nd, type = "link", se.fit = TRUE)
  eta  <- as.matrix(pr$fit); seta <- as.matrix(pr$se.fit)
  nlp  <- fam$nlp

  resp_fit <- vapply(seq_len(nlp),
                     function(j) fam$linfo[[j]]$linkinv(eta[, j]), numeric(n))
  if (is.null(dim(resp_fit))) resp_fit <- matrix(resp_fit, n, nlp)

  ## rotate fitted directions (and observed angles) back to the original frame
  ## when the fit was centered (center = TRUE); 0 otherwise. Concentration/shape
  ## and the circular-SD band width are rotation-invariant, so only the circular
  ## location curve, the derived direction and the observed overlay shift.
  ref <- if (!is.null(x$circ_center)) x$circ_center else 0

  ## A circular location/direction band shows the fitted law's PREDICTIVE angular
  ## spread: +/- its circular standard deviation sqrt(-2 log R), R the mean
  ## resultant length, computed once from the whole response-scale fit. A
  ## non-circular parameter (the concentration, ...) keeps a delta-method 2-SE
  ## band; a linear response has no circular SD and keeps 2-SE throughout.
  csd <- if (g$resp_circular) .circ_sd_quad(fam, resp_fit) else NULL

  panels <- list()
  for (j in seq_len(nlp)) {
    if (isTRUE(fam$param_circular[j]))
      panels[[length(panels) + 1L]] <- list(
        name = pnm[j], circular = TRUE, mid = wrap(resp_fit[, j] + ref), csd = csd)
    else
      panels[[length(panels) + 1L]] <- list(
        name = pnm[j], circular = FALSE, mid = resp_fit[, j],
        lo = fam$linfo[[j]]$linkinv(eta[, j] - 2 * seta[, j]),
        hi = fam$linfo[[j]]$linkinv(eta[, j] + 2 * seta[, j]))
  }
  ## derived parameters (e.g. pnlss's mean direction): a circular direction, so
  ## the same +/- circular-SD predictive band as a native circular location.
  if (!is.null(fam$derived))
    for (nm in names(fam$derived))
      panels[[length(panels) + 1L]] <- list(
        name = nm, circular = TRUE, mid = wrap(fam$derived[[nm]](resp_fit) + ref),
        csd = if (identical(nm, "direction")) csd else NULL)

  yobs <- if (!is.null(x$y)) wrap(as.numeric(x$y) + ref) else NULL

  ## The panel placed on the geometry surface: the circular location (the angle
  ## wrapping the tube) for a circular response, else the mean (the value taken
  ## as the height on the can) for a linear response. overlay_idx is the panel
  ## the observed responses are drawn on -- never NA when a surface exists.
  if (g$resp_circular) {
    loc <- which(vapply(panels, `[[`, logical(1), "circular"))
    surf_idx <- if (length(loc)) loc[1] else NA_integer_
  } else surf_idx <- 1L
  overlay_idx <- surf_idx

  ## the natural surface for the leg: cylinder (c~l), torus (c~c), can (l~c)
  surface <- switch(g$kind, cl = "cylinder", cc = "torus", lc = "can",
                    NA_character_)

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))

  if (view != "flat" && (is.na(surface) || is.na(surf_idx))) {
    message("plot.circ_gam: no surface for this fit; drawing the flat view.")
    view <- "flat"
  }

  if (view == "flat") {
    if (isTRUE(pages == 1)) {
      nc <- ceiling(sqrt(length(panels))); nr <- ceiling(length(panels) / nc)
      graphics::par(mfrow = c(nr, nc))
    }
    for (i in seq_along(panels))
      .circ_panel(grid, panels[[i]], cov, se,
                  xobs = if (isTRUE(i == overlay_idx)) xv else NULL,
                  yobs = if (isTRUE(i == overlay_idx)) yobs else NULL, rug = rug)
    return(invisible(x))
  }

  sp  <- panels[[surf_idx]]
  glo <- if (!se) NULL else if (!is.null(sp$csd)) sp$mid - sp$csd else sp$lo
  ghi <- if (!se) NULL else if (!is.null(sp$csd)) sp$mid + sp$csd else sp$hi
  if (view == "both") .circ_both_layout(length(panels))
  .circ_geometry_panel(grid, sp$mid, xv, yobs, surface, lo = glo, hi = ghi)
  ## "both" is the geometry canvas plus the FULL flat panel set -- the same panels
  ## "flat" draws (one per distribution parameter), so the two views never disagree.
  if (view == "both") {
    graphics::par(mar = c(4.1, 4.1, 2.6, 1.1))      # undo the geometry's tight margins
    for (i in seq_along(panels))
      .circ_panel(grid, panels[[i]], cov, se,
                  xobs = if (isTRUE(i == overlay_idx)) xv else NULL,
                  yobs = if (isTRUE(i == overlay_idx)) yobs else NULL, rug = rug)
  }
  invisible(x)
}

## The view = "both" page: the geometry canvas and the P flat parameter panels as
## ONE near-square grid -- geometry is just an extra cell the same size as each
## panel (the flat-view grid with one slot added), so no panel is squeezed into a
## shared half. Drawn first, geometry takes cell 1; the panels fill cells 2..P+1
## row-major, any trailing cell left blank. Shared by plot.circ_gam and plot.circ_lm.
.circ_both_layout <- function(P) {
  m  <- P + 1L
  nc <- ceiling(sqrt(m)); nr <- ceiling(m / nc)
  graphics::par(mfrow = c(nr, nc))
}

## a filled (shadow) SE band that breaks cleanly at the +/-pi wrap: polygon each
## maximal run of finite lo/hi (an NA marks a branch break). Shared by the circ_gam
## and circ_mix flat panels so the fill matches the geometry ribbon -- a shadow
## everywhere, never vertical segments.
.circ_band_fill <- function(grid, lo, hi, col) {
  ok <- is.finite(lo) & is.finite(hi)
  if (!any(ok)) return(invisible())
  r <- rle(ok); end <- cumsum(r$lengths); start <- end - r$lengths + 1L
  for (i in which(r$values)) {
    ix <- start[i]:end[i]
    graphics::polygon(c(grid[ix], rev(grid[ix])), c(lo[ix], rev(hi[ix])),
                      border = NA, col = col)
  }
}

## A circular location curve/band is drawn on the flat (-pi, pi] panel by
## UNWRAPPING it to a continuous phase, then drawing every 2pi copy that can reach
## the panel and letting the plot region clip the rest. That keeps both the curve
## and its band continuous across the +/-pi cut -- they exit one edge and re-enter
## the opposite, and where the location winds the wrapped copies of the curve sit
## under the wrapped band -- instead of a principal branch broken at the cut.

## unwrap a wrapped angle to a continuous phase within each maximal finite run
## (cumulate per-step differences folded into (-pi, pi]); returns runs of (ix, phase)
.circ_unwrap_runs <- function(v) {
  ok <- is.finite(v); if (!any(ok)) return(list())
  r <- rle(ok); end <- cumsum(r$lengths); start <- end - r$lengths + 1L
  lapply(which(r$values), function(i) {
    ix <- start[i]:end[i]; m <- v[ix]
    if (length(m) > 1L) m <- cumsum(c(m[1L], atan2(sin(diff(m)), cos(diff(m)))))
    list(ix = ix, phase = m)
  })
}

## the 2pi shifts k for which [lo_min, hi_max] + 2pi k intersects (-pi, pi]
.circ_tile_k <- function(lo_min, hi_max)
  seq.int(ceiling((-pi - hi_max) / (2 * pi)), floor((pi - lo_min) / (2 * pi)))

## the +/- circular-SD band of `mid` with angular half-width `csd` (a scalar or a
## per-grid vector; capped at pi upstream, so at most the whole circle is shaded)
.circ_band_fill_circular <- function(grid, mid, csd, col) {
  csd <- rep_len(csd, length(mid))
  for (run in .circ_unwrap_runs(mid)) {
    g <- grid[run$ix]; lo <- run$phase - csd[run$ix]; hi <- run$phase + csd[run$ix]
    for (k in .circ_tile_k(min(lo), max(hi)))
      graphics::polygon(c(g, rev(g)), c(lo + 2 * pi * k, rev(hi + 2 * pi * k)),
                        border = NA, col = col)
  }
}

## the circular location curve, tiled to match its band (continuous across the cut)
.circ_lines_circular <- function(grid, mid, col, lwd) {
  for (run in .circ_unwrap_runs(mid)) {
    g <- grid[run$ix]; m <- run$phase
    for (k in .circ_tile_k(min(m), max(m)))
      graphics::lines(g, m + 2 * pi * k, col = col, lwd = lwd)
  }
}

## Circular standard deviation sqrt(-2 log R) of the fitted law at each grid
## point -- the PREDICTIVE angular spread a circular band shows, not a CI of the
## mean. R is the mean resultant length, obtained by deterministic quadrature of
## the family's OWN density: its ll() l0 on a fine angle grid at the fitted
## parameters, stacked into ONE call (X[, k] holds the per-row link-scale
## parameter so coef = 1 recovers each eta). One exact code path for every
## circlss family -- von Mises (R = A1(kappa)), wrapped Cauchy/normal/cardioid
## (rho), Cartwright (1/(zeta+1)), Kato-Jones (gamma), projected normal and the
## skew/flexible laws whose R has no closed form alike. `fitted` is the n x nlp
## response-scale parameter matrix.
.circ_sd_quad <- function(fam, fitted, M = 512L) {
  fitted <- as.matrix(fitted); nlp <- fam$nlp; n <- nrow(fitted)
  tg  <- seq(-pi, pi, length.out = M + 1L)[-(M + 1L)]
  eta <- vapply(seq_len(nlp),
                function(k) fam$linfo[[k]]$linkfun(fitted[, k]), numeric(n))
  if (is.null(dim(eta))) eta <- matrix(eta, n, nlp)
  X   <- eta[rep(seq_len(n), each = M), , drop = FALSE]
  attr(X, "lpi") <- as.list(seq_len(nlp))
  l0  <- fam$ll(y = rep(tg, times = n), X = X, coef = rep(1, nlp),
                wt = rep(1, n * M), family = fam, deriv = 0)$l0
  f   <- matrix(exp(l0), M, n)
  R   <- sqrt(colSums(cos(tg) * f)^2 + colSums(sin(tg) * f)^2) / colSums(f)
  .circ_sd_from_R(R)
}

## sqrt(-2 log R) with R held off 0 and 1 (finite band) and the result capped at
## pi (a half-circle each way == the full circle, the honest near-uniform limit,
## rather than an over-wrapping ribbon). Shared with the closed-form circ_lm path
## (R = A1(kappa)).
.circ_sd_from_R <- function(R) {
  R <- pmin(pmax(R, .Machine$double.eps), 1 - 1e-12)
  pmin(sqrt(-2 * log(R)), pi)
}

## ---- flat panel -----------------------------------------------------------
.circ_panel <- function(grid, p, xlab, se, xobs, yobs, rug) {
  circ <- isTRUE(p$circular)
  ylim <- if (circ) c(-pi, pi) else range(c(p$mid, p$lo, p$hi), na.rm = TRUE)
  if (!is.null(yobs)) ylim <- range(c(ylim, yobs), na.rm = TRUE)
  band <- grDevices::adjustcolor("steelblue", 0.25)
  graphics::plot(grid, p$mid, type = "n", xlab = xlab, ylab = p$name,
                 main = p$name, ylim = ylim)
  if (!is.null(xobs) && !is.null(yobs))
    graphics::points(xobs, yobs, pch = 16, cex = 0.4,
                     col = grDevices::adjustcolor("black", 0.25))
  if (se) {
    if (circ && !is.null(p$csd))
      .circ_band_fill_circular(grid, p$mid, p$csd, band)   # +/- circular SD, wrapped
    else if (!is.null(p$lo) && !is.null(p$hi))
      .circ_band_fill(grid, p$lo, p$hi, band)              # 2-SE shadow
  }
  if (circ) .circ_lines_circular(grid, p$mid, "steelblue", 2)   # wraps with its band
  else      graphics::lines(grid, p$mid, col = "steelblue", lwd = 2)
  if (rug && !is.null(xobs)) graphics::rug(xobs, col = grDevices::adjustcolor("black", 0.3))
}

## ---- geometry panel: cylinder (CL) or torus (CC) --------------------------
## persp() for the projection, trans3d() for overlays, depth-faded curve/points
## so the back of the surface reads lighter than the front.
.circ_depth3d <- function(x, y, z, pm) { p <- cbind(x, y, z, 1) %*% pm; p[, 3] / p[, 4] }

.circ_fade_curve <- function(x, y, z, pm, col, lwd = 3, lo = 0.22) {
  dp <- .circ_depth3d(x, y, z, pm); cv <- grDevices::trans3d(x, y, z, pm)
  a  <- lo + (1 - lo) * (dp - min(dp)) / diff(range(dp))
  o  <- order(utils::head(dp, -1))
  graphics::segments(utils::head(cv$x, -1)[o], utils::head(cv$y, -1)[o],
                     utils::tail(cv$x, -1)[o], utils::tail(cv$y, -1)[o],
                     col = vapply(utils::head(a, -1)[o],
                                  function(ai) grDevices::adjustcolor(col, ai), ""),
                     lwd = lwd, lend = 1)
}

.circ_fade_points <- function(x, y, z, pm, col, cex = 0.5, lo = 0.12) {
  dp <- .circ_depth3d(x, y, z, pm); cv <- grDevices::trans3d(x, y, z, pm)
  a  <- lo + (1 - lo) * (dp - min(dp)) / diff(range(dp)); o <- order(dp)
  graphics::points(cv$x[o], cv$y[o], pch = 16, cex = cex,
                   col = vapply(a[o], function(ai) grDevices::adjustcolor(col, ai), ""))
}

## a translucent uncertainty ribbon between two response curves vlo(u) and vhi(u)
## on the surface: u is the covariate coordinate the curve is drawn against and
## vlo/vhi are the band edges in the same units as the curve (the same xyz() maps
## both). Subdivided across the band into quads, depth-sorted and painted
## back-to-front so it sits on the surface like the curve does.
.circ_fade_ribbon <- function(u, vlo, vhi, xyz, pm, col, K = 4L,
                              a_lo = 0.05, a_hi = 0.18) {
  i0 <- seq_len(length(u) - 1L)
  quads <- vector("list", K * length(i0)); m <- 0L
  for (k in seq_len(K)) {
    t0 <- vlo + (vhi - vlo) * (k - 1L) / K
    t1 <- vlo + (vhi - vlo) *  k        / K
    for (i in i0) {
      w <- xyz(c(u[i], u[i + 1L], u[i + 1L], u[i]),
               c(t0[i], t0[i + 1L], t1[i + 1L], t1[i]))
      m <- m + 1L
      quads[[m]] <- list(p = grDevices::trans3d(w$x, w$y, w$z, pm),
                         d = mean(.circ_depth3d(w$x, w$y, w$z, pm)))
    }
  }
  dq <- vapply(quads, `[[`, numeric(1), "d")
  span <- diff(range(dq)); if (span <= 0) span <- 1
  aq <- a_lo + (a_hi - a_lo) * (dq - min(dq)) / span
  for (j in order(dq))
    graphics::polygon(quads[[j]]$p$x, quads[[j]]$p$y, border = NA,
                      col = grDevices::adjustcolor(col, aq[j]))
}

## the shared surface canvas: draw the persp() projection + wireframe for a leg,
## and return the projection matrix `pm`, the internal-coordinate mapper `xyz`, and
## the transforms `to_u` (RAW covariate -> surface u-coordinate) and `to_v` (RAW
## response -> surface v-coordinate). Drawn ONCE; a caller overlays curves/points/
## bands onto it via xyz(to_u(.), to_v(.)). plot.circ_gam overlays one fitted curve,
## plot.circ_mix one per component -- same canvas, so a leg renders identically
## either way. `yspan` sets the can's height range (the response range it must span).
.circ_surface <- function(surface, grid, xobs, yspan = NULL) {
  wire <- "gray88"
  if (surface == "torus") {
    R <- 2.0; r <- 0.82
    xyz  <- function(u, v) list(x = (R + r * cos(v)) * cos(u),    # u: ring angle
                                y = (R + r * cos(v)) * sin(u), z = r * sin(v))  # v: tube angle
    to_u <- function(v) v; to_v <- function(v) v
    lim <- R + r + 0.1
    pm <- graphics::persp(c(-lim, lim), c(-lim, lim), matrix(c(-r, -r, r, r), 2, 2),
                          zlim = c(-r - 0.4, r + 0.4), theta = 30, phi = 40, d = 4,
                          scale = FALSE, expand = 1, col = NA, border = NA,
                          box = FALSE, axes = FALSE)
    thd <- seq(-pi, pi, length.out = 56)
    for (p in seq(-pi, pi, length.out = 33)[-1]) {
      w <- xyz(rep(p, length(thd)), thd)
      graphics::lines(grDevices::trans3d(w$x, w$y, w$z, pm), col = wire, lwd = 0.4)
    }
    phd <- seq(-pi, pi, length.out = 150)
    for (t in seq(-pi, pi, length.out = 17)[-1]) {
      w <- xyz(phd, rep(t, length(phd)))
      graphics::lines(grDevices::trans3d(w$x, w$y, w$z, pm), col = wire, lwd = 0.4)
    }
  } else if (surface == "cylinder") {
    rho <- 0.95; ax <- 3.2
    xr   <- range(c(grid, xobs), na.rm = TRUE)
    to_u <- function(v) (v - xr[1]) / diff(xr) * 2 * ax - ax      # covariate -> axis
    to_v <- function(v) v
    xyz  <- function(u, v) list(x = u, y = rho * cos(v), z = rho * sin(v))  # u: axis, v: tube angle
    limx <- 3.7
    pm <- graphics::persp(c(-limx, limx), c(-1.6, 1.6), matrix(c(-rho, -rho, rho, rho), 2, 2),
                          zlim = c(-1.4, 1.4), theta = 24, phi = 16, d = 5,
                          scale = FALSE, expand = 1, col = NA, border = NA,
                          box = FALSE, axes = FALSE)
    thd <- seq(-pi, pi, length.out = 70)
    for (xx in seq(-ax, ax, length.out = 17)) {
      w <- xyz(rep(xx, length(thd)), thd)
      graphics::lines(grDevices::trans3d(w$x, w$y, w$z, pm), col = wire, lwd = 0.4)
    }
    xd <- seq(-ax, ax, length.out = 120)
    for (t in seq(-pi, pi, length.out = 13)[-1]) {
      w <- xyz(xd, rep(t, length(xd)))
      graphics::lines(grDevices::trans3d(w$x, w$y, w$z, pm), col = wire, lwd = 0.4)
    }
  } else {
    ## "can": the circular COVARIATE wraps the ring, the linear RESPONSE is the
    ## height -- the c~l cylinder with the roles swapped (the founding l~c view).
    r <- 1.0; H <- 1.5
    yl   <- if (is.null(yspan)) c(-1, 1) else range(yspan, na.rm = TRUE)
    to_u <- function(v) v
    to_v <- function(v) (v - yl[1]) / diff(yl) * 2 * H - H        # response -> height
    xyz  <- function(u, v) list(x = r * cos(u), y = r * sin(u), z = v)  # v: height (pre-mapped)
    pm <- graphics::persp(c(-r - 0.3, r + 0.3), c(-r - 0.3, r + 0.3),
                          matrix(c(-H, -H, H, H), 2, 2), zlim = c(-H - 0.3, H + 0.3),
                          theta = 30, phi = 18, d = 4, scale = FALSE, expand = 1,
                          col = NA, border = NA, box = FALSE, axes = FALSE)
    hd <- seq(-H, H, length.out = 40)
    for (p in seq(-pi, pi, length.out = 33)[-1]) {        # vertical generators
      w <- xyz(rep(p, length(hd)), hd)
      graphics::lines(grDevices::trans3d(w$x, w$y, w$z, pm), col = wire, lwd = 0.4)
    }
    phd <- seq(-pi, pi, length.out = 150)
    for (t in seq(-H, H, length.out = 13)) {              # horizontal rings
      w <- xyz(phd, rep(t, length(phd)))
      graphics::lines(grDevices::trans3d(w$x, w$y, w$z, pm), col = wire, lwd = 0.4)
    }
  }
  list(pm = pm, xyz = xyz, to_u = to_u, to_v = to_v)
}

## one fitted curve on its surface (plot.circ_gam): the canvas, then the band, the
## location curve and the observed points -- output identical to before, now via
## the shared .circ_surface() so plot.circ_mix can draw K curves on the same canvas.
.circ_geometry_panel <- function(grid, zv, xobs, yobs, surface,
                                 lo = NULL, hi = NULL, main = NULL) {
  col_curve <- "#c0392b"; col_data <- "#1f4e79"
  band <- !is.null(lo) && !is.null(hi)
  graphics::par(mar = c(1.0, 0.3, 2.2, 0.3))
  yspan <- if (surface == "can") c(zv, yobs, lo, hi) else NULL
  s <- .circ_surface(surface, grid, xobs, yspan)
  if (band) .circ_fade_ribbon(s$to_u(grid), s$to_v(lo), s$to_v(hi), s$xyz, s$pm, col_curve)
  w <- s$xyz(s$to_u(grid), s$to_v(zv)); .circ_fade_curve(w$x, w$y, w$z, s$pm, col_curve, lwd = 3)
  if (!is.null(xobs)) { w <- s$xyz(s$to_u(xobs), s$to_v(yobs))
    .circ_fade_points(w$x, w$y, w$z, s$pm, col_data) }
  if (is.null(main)) main <- switch(surface, torus = "torus \u00b7 circular\u2013circular",
                                    cylinder = "cylinder \u00b7 circular\u2013linear",
                                    "can (linear-circular)")
  graphics::title(main = main, line = 0.5, cex.main = 1.15)
}

## the single covariate of a fit, or NULL if not exactly one
.circ_covariate <- function(x) {
  fl <- x$formula; if (inherits(fl, "formula")) fl <- list(fl)
  resp <- all.vars(fl[[1]][[2L]])
  rhs  <- unique(unlist(lapply(fl, function(f) all.vars(f[[length(f)]]))))
  cov  <- setdiff(rhs, resp)
  if (length(cov) == 1L) cov else NULL
}

## defer to mgcv's plot.gam (term plots) on the un-circ_gam'd object
.circ_plot_gam <- function(x, pages, ...) {
  class(x) <- setdiff(class(x), "circ_gam")
  graphics::plot(x, pages = pages, ...)
}
