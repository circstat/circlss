#' Circular-response GAM
#'
#' A wrapper around \code{\link[mgcv]{gam}} for circular regression with the
#' circlss families. It supplies the defaults a circular fit needs: the
#' cyclic-smooth knots, \code{method = "REML"}, a trailing \code{~ 1} fill for any
#' distribution parameters left without a formula, and response-scale output
#' columns named by the family's parameters. Formulas and smoothing bases are
#' passed to \pkg{mgcv} unchanged.
#'
#' @param formula A formula, or a list of formulas (one per distribution
#'   parameter), exactly as \code{\link[mgcv]{gam}} expects for a general family.
#'   The first formula is two-sided and names the response. If fewer formulas than
#'   the family has parameters are supplied, the remainder are filled with
#'   \code{~ 1} (intercept-only), so \code{theta ~ s(x)} with \code{vmlss()} models
#'   the mean direction and holds the concentration constant.
#' @param data A data frame holding the response and covariates. Required: the
#'   covariate values are what set the cyclic-smooth knots.
#' @param family Any \pkg{mgcv} family. A circlss \emph{circular} family
#'   (\code{\link{vmlss}}, \code{\link{pnlss}}, \dots) selects a circular (angular)
#'   response; a circlss \emph{linear} location-scale family (\code{\link{gausslss}},
#'   \code{\link{gammalss}}) selects a real-valued response on a circular covariate,
#'   the linear--circular leg. Either switches on the named, response-scale output,
#'   the trailing \code{~ 1} fill, and the geometry-aware \code{print}/\code{plot}. An
#'   ordinary mgcv family such as \code{\link[stats]{gaussian}} is forwarded unchanged
#'   (it still gets the cyclic-knot default).
#' @param method Smoothing-parameter selection criterion; defaults to
#'   \code{"REML"} (mgcv's own default is \code{"GCV.Cp"}).
#' @param knots Optional knot list passed to \code{gam}. Knots you supply are
#'   always respected; for any cyclic smooth (\code{bs = "cc"} or \code{"cp"}) you
#'   do \emph{not} supply, \code{circ_gam} fills the knots with one full radian
#'   period bracketing the data -- \code{c(-pi, pi)} when the covariate takes
#'   negative values, else \code{c(0, 2*pi)}.
#' @param weights Optional prior weights on the observations, one per row of
#'   \code{data}, as in \code{\link[mgcv]{gam}}.
#' @param subset Optional vector selecting the rows of \code{data} to fit, as in
#'   \code{\link[mgcv]{gam}}.
#' @param na.action How missing values are handled, as in \code{\link[mgcv]{gam}};
#'   when omitted, mgcv's option-driven default (\code{\link[stats]{na.omit}})
#'   applies.
#' @param offset Optional model offset, as in \code{\link[mgcv]{gam}}.
#' @param center For a circular response on the tan-half link (\code{vmlss} and
#'   the other tan-link families; never \code{pnlss}, which has no wall), rotate
#'   the response to a frame where the unreachable antipode \eqn{\theta = \pi}
#'   lies away from the data before fitting, then report all directions back in
#'   the original frame. \code{TRUE} (default) chooses the reference
#'   automatically (the centre of the occupied arc) and is an exact no-op when
#'   the data already clear the wall; \code{FALSE} disables it; a numeric value
#'   sets the reference angle directly. Only response-scale directions
#'   (\code{predict(type = "response")}, \code{fitted}, \code{plot}) are rotated
#'   back; link-scale output (\code{coef}, \code{predict(type = "link")}) stays
#'   in the centred frame. The applied rotation is stored in
#'   \code{fit$circ_center}. A mean that must \emph{wind} across the wall cannot
#'   be centred away -- use \code{\link{pnlss}} there.
#' @param object,x A fitted \code{circ_gam} model.
#' @param ... Further arguments passed verbatim to \code{\link[mgcv]{gam}}
#'   (\code{optimizer}, \code{control}, \code{sp}, \code{select}, \dots).
#' @return A fitted model of class \code{c("circ_gam", "gam", "glm", "lm")}: the
#'   \code{\link[mgcv]{gam}} object, with circular-aware \code{predict},
#'   \code{fitted}, \code{print} and \code{plot} methods. Every other \pkg{mgcv}
#'   method (\code{summary}, \code{AIC}, \code{logLik}, \code{gam.check}, \dots) is
#'   inherited unchanged.
#' @details
#' \strong{Knots.} A cyclic smooth's period is its knot span. Because the data are
#' in radians the period is \eqn{2\pi} and the wrap points are \eqn{\pm\pi} (or
#' \eqn{0} and \eqn{2\pi}); without explicit knots, \code{\link[mgcv]{gam}} wraps
#' the cyclic basis at the observed data range instead. \code{circ_gam} fills these
#' knots automatically, and stops if a cyclic covariate is not on a single radian
#' branch (convert it first with \code{\link{rad}} or \code{\link{wrap}}). Knots
#' supplied through \code{knots} override this for that covariate.
#'
#' \strong{Output columns.} For a circlss family, \code{predict(type = "response")}
#' and \code{fitted()} label their columns with the family's parameters -- for
#' example \code{mu, kappa} for \code{\link{vmlss}}, \code{mu1, mu2} for
#' \code{\link{pnlss}}, and \code{xi, kappa, psi, lambda} for \code{\link{ssjplss}}.
#'
#' \strong{Optimizer.} \code{\link[mgcv]{gam}} selects the extended Fellner--Schall
#' optimizer automatically for families that supply only first- and second-order
#' derivatives (\code{available.derivs = 0}); pass \code{optimizer = "efs"} through
#' \code{...} to force it on a family that would otherwise use full Newton REML.
#' @examples
#' library(mgcv)
#' set.seed(1); n <- 300
#'
#' ## circular-linear: von Mises, both parameters smooth in a covariate
#' x <- runif(n); y <- 2 * atan(1.2 * sin(2 * pi * x)) + rnorm(n) / 3
#' b <- circ_gam(list(y ~ s(x), ~ s(x)), data = data.frame(y, x), family = vmlss())
#' head(predict(b, type = "response"))   # columns named mu, kappa
#' plot(b)
#'
#' ## circular-circular: projected normal, cyclic covariate -- knots auto-pinned
#' phi <- runif(n, -pi, pi); yc <- 2 * atan(0.9 * sin(phi)) + rnorm(n) / 3
#' b2 <- circ_gam(list(yc ~ s(phi, bs = "cc"), ~ s(phi, bs = "cc")),
#'                data = data.frame(yc, phi), family = pnlss())
#'
#' ## linear-circular: a real response on a cyclic covariate (the "can")
#' yl <- 2 + 1.5 * sin(phi) + rnorm(n) / 3
#' b3 <- circ_gam(list(yl ~ s(phi, bs = "cc"), ~ s(phi, bs = "cc")),
#'                data = data.frame(yl, phi), family = gausslss())
#' plot(b3, view = "both")
#' @seealso
#' \code{\link{vmlss}}, \code{\link{pnlss}}, \code{\link{rad}},
#' \code{\link[mgcv]{gam}}
#' @export
circ_gam <- function(formula, data, family = vmlss(), method = "REML",
                     knots = NULL, weights = NULL, subset = NULL,
                     na.action, offset = NULL, center = TRUE, ...) {
  if (missing(data))
    stop("circ_gam() requires 'data': the covariate values set the ",
         "cyclic-smooth knots.")

  ## a circlss family carries param metadata; an ordinary mgcv family does not
  is_circlss <- !is.null(family$param_names)

  flist <- if (inherits(formula, "formula")) list(formula) else formula
  if (!is.list(flist) || !length(flist) ||
      !all(vapply(flist, inherits, logical(1), "formula")))
    stop("'formula' must be a formula or a non-empty list of formulas.")

  if (is_circlss) {
    nlp <- family$nlp
    pnm <- family$param_names
    if (length(flist[[1]]) != 3L)
      stop("the first formula must be two-sided (name the response), ",
           "e.g. theta ~ s(x).")
    if (length(flist) > nlp)
      stop(sprintf("family '%s' has %d parameters (%s) but %d formulas ",
                   family$family, nlp, paste(pnm, collapse = ", "),
                   length(flist)),
           "were supplied.")
    if (length(flist) < nlp)            # trailing ~ 1 fill (intercept-only)
      flist <- c(flist, rep(list(stats::as.formula("~ 1")),
                            nlp - length(flist)))
  }

  ## center: rotate the circular response to a safe frame so the tan-half wall
  ## (the antipode of the link origin, theta = pi) lands away from the data, fit
  ## there, and let predict/fitted/plot rotate the reported directions back. A
  ## no-op for families with no wall (pnlss's derived direction, linear-response
  ## l~c) and for data already clear of the wall (`.circ_center_ref` snaps to 0).
  ## Only the response LHS of the first formula is touched; the cyclic-smooth
  ## RHS (and therefore `knots` below) is untouched.
  ref <- 0
  if (is_circlss && !isFALSE(center) && !is.na(.circ_wall_loc(family))) {
    resp_lhs <- flist[[1L]][[2L]]
    yc <- tryCatch(as.numeric(eval(resp_lhs, data, environment(flist[[1L]]))),
                   error = function(e) NULL)
    if (!is.null(yc) && length(yc)) {
      wts <- if (!missing(weights) && is.numeric(weights) &&
                 length(weights) == length(yc)) weights else NULL
      ref <- if (is.numeric(center)) as.numeric(center)[1L]
             else .circ_center_ref(yc, wts)
      if (isTRUE(is.finite(ref)) && ref != 0) {
        if (is.name(resp_lhs) && !is.null(data[[as.character(resp_lhs)]]))
          data[[as.character(resp_lhs)]] <- wrap(yc - ref)
        else {
          data[[".circ_centered_y"]] <- wrap(yc - ref)
          flist[[1L]][[2L]] <- as.name(".circ_centered_y")
        }
      } else ref <- 0
    }
  }

  knots <- .circ_cc_knots(flist, data, knots)

  ## single formula + ordinary family stays a single formula (mgcv treats a
  ## one-element list as a one-LP general setup, which gaussian() is not)
  form_arg <- if (is_circlss || is.list(formula)) flist else formula

  ## weights/subset/offset/na.action are forwarded BY NAME, never through `...`:
  ## mgcv::gam runs non-standard evaluation (match.call) on them, so routing them
  ## through `...` corrupts the captured model-frame call and errors with
  ## "..1 used in an incorrect context, no ... to look in". na.action is passed
  ## only when the caller supplies it, so mgcv keeps its default (na.omit) otherwise.
  gam_args <- list(form_arg, family = family, data = data, method = method,
                   knots = knots, weights = weights, subset = subset,
                   offset = offset, ...)
  if (!missing(na.action)) gam_args[["na.action"]] <- na.action
  fit <- do.call(mgcv::gam, gam_args)
  class(fit) <- c("circ_gam", class(fit))
  fit$circ_center <- ref                # rotation applied at fit time (0 if none)
  fit
}

## ---- centering off the tan-half wall (`center = TRUE`) -----------------------

## Column index of the circular location parameter that sits on the tan-half
## link -- the one carrying the antipode wall -- or NA when the family has none:
## pnlss reports a DERIVED atan2 direction (no wall, param_circular all FALSE) and
## linear-response families have no circular location. Uses the same
## `param_circular` flag plot.circ_gam() keys on, so the identical column is
## rotated everywhere. (Every native circular location in circlss is tan-half.)
.circ_wall_loc <- function(fam) {
  pc <- fam$param_circular
  if (is.null(pc) || !any(pc)) return(NA_integer_)
  which(pc)[1L]
}

## Reference direction to rotate the response by so the tan-half wall lands away
## from the data. The centre is the (weighted) CIRCULAR MEAN -- the canonical
## location for circular data, and the only rotation-equivariant one: it makes
## already-mean-centred data an exact fixed point (no second shift), unlike a
## largest-gap rule whose cut disagrees with the mean for skewed/wide responses
## and can land the wall on a regression's swept mean trajectory.
##
## We rotate the mean to 0 ONLY when it sits in the wall's half of the circle;
## data whose mean is well clear of the wall is left untouched (exact no-op).
## Unweighted = circ_gam (the whole response); weighted = a circ_mix component's
## responsibility-weighted mean. The two branches now use the same rule.
.circ_center_ref <- function(theta, weights = NULL, snap = 0.05) {
  ok <- is.finite(theta)
  theta <- theta[ok]
  if (!length(theta)) return(0)
  w  <- if (is.null(weights)) rep.int(1, length(theta)) else weights[ok]
  mu <- atan2(sum(w * sin(theta)), sum(w * cos(theta)))   # (weighted) circular mean
  ## recentre only when the mean is in the wall's half of the circle
  ref <- if (abs(wrap(mu - pi)) < pi / 2 + snap) mu else 0
  if (!is.finite(ref) || abs(ref) < snap) 0 else ref
}

## Rotate the circular-location column of a response-scale matrix back to the
## original frame. No-op when ref is 0/NULL or the family has no wall column.
.circ_rotate_response <- function(p, fam, ref) {
  if (is.null(ref) || !is.finite(ref) || ref == 0) return(p)
  loc <- .circ_wall_loc(fam)
  if (is.na(loc)) return(p)
  if (is.matrix(p) && ncol(p) >= loc) p[, loc] <- wrap(p[, loc] + ref)
  p
}

## Raw (fit-frame) response-scale fitted matrix, bypassing the rotate-back in
## fitted.circ_gam so internal consumers (diagnostics) stay self-consistent with
## object$y, which is also stored in the fit frame.
.circ_fitted_raw <- function(object)
  .circ_name_response(object$fitted.values, object$family, "response")

## Variables entering a cyclic smooth (bs %in% {"cc","cp"}) anywhere in the
## formula list. Walks the call tree for s()/te()/ti()/t2() with a literal
## bs = "cc"/"cp"; reads no mgcv internals and changes nothing.
.circ_cc_terms <- function(flist) {
  cyclic <- c("cc", "cp")
  vars <- character(0)
  walk <- function(e) {
    if (!is.call(e)) return(invisible())
    fn <- e[[1L]]
    if (is.symbol(fn) && as.character(fn) %in% c("s", "te", "ti", "t2")) {
      a <- as.list(e)[-1L]
      bs <- a[["bs"]]
      if (is.character(bs) && any(bs %in% cyclic)) {
        nm <- names(a); if (is.null(nm)) nm <- rep("", length(a))
        for (i in seq_along(a))
          if (!nzchar(nm[i]) && is.symbol(a[[i]]))
            vars <<- c(vars, as.character(a[[i]]))
      }
    } else for (i in seq_along(e)) walk(e[[i]])
  }
  for (f in flist) walk(f[[length(f)]])   # RHS of each formula
  unique(vars)
}

## Default the knots of every cyclic covariate to one full radian period that
## brackets the data -- c(-pi, pi) when the covariate goes negative, else
## c(0, 2*pi) -- so the cc basis wraps at the true period, not the observed
## data range. User-supplied knots always win.
.circ_cc_knots <- function(flist, data, knots) {
  vars <- .circ_cc_terms(flist)
  if (!length(vars)) return(knots)
  if (is.null(knots)) knots <- list()
  period <- 2 * pi
  for (v in vars) {
    if (!is.null(knots[[v]])) next                 # user knots win
    if (is.null(data[[v]])) next                   # not in data; mgcv will say so
    rng <- range(as.numeric(data[[v]]), na.rm = TRUE)
    if (diff(rng) > period + 1e-6)
      stop(sprintf(
        "cyclic covariate '%s' spans %.3f > 2*pi; it must be in radians on one branch.\n  Convert it (see ?rad / ?wrap) or pass knots = list(%s = ...) explicitly.",
        v, diff(rng), v))
    k <- if (rng[1] < 0) c(-pi, pi) else c(0, period)
    if (rng[1] < k[1] - 1e-6 || rng[2] > k[2] + 1e-6)
      stop(sprintf(
        "cyclic covariate '%s' lies outside its radian period [%.3f, %.3f].\n  Convert it (see ?rad / ?wrap) or pass knots = list(%s = ...) explicitly.",
        v, k[1], k[2], v))
    knots[[v]] <- k
  }
  knots
}

## The geometry of a fit, from two independent facts: is the RESPONSE an angle
## (a family fact -- every circular family carries param_names and omits the
## response_circular flag; the forked linear families set it FALSE), and does the
## single covariate enter a cyclic smooth. Their 2x2 names the leg and the
## surface its plot draws. cov is NULL (kind NA) when the fit is not exactly one
## covariate. Shared by print.circ_gam and plot.circ_gam.
.circ_geometry <- function(x) {
  fam <- x$family
  resp_circular <- !is.null(fam$param_names) && !isFALSE(fam$response_circular)
  cov <- .circ_covariate(x)
  fl  <- x$formula; if (inherits(fl, "formula")) fl <- list(fl)
  cov_circular <- !is.null(cov) && cov %in% .circ_cc_terms(fl)
  kind <- if (is.null(cov)) NA_character_
          else if ( resp_circular && !cov_circular) "cl"   # circular-linear   cylinder
          else if ( resp_circular &&  cov_circular) "cc"   # circular-circular torus
          else if (!resp_circular &&  cov_circular) "lc"   # linear-circular   can
          else                                      "ll"   # linear-linear     none
  list(kind = kind, resp_circular = resp_circular,
       cov_circular = cov_circular, cov = cov)
}

## predict()/fitted() with response-scale columns named by the family params.
#' @rdname circ_gam
#' @export
predict.circ_gam <- function(object, ...) {
  type <- list(...)$type; if (is.null(type)) type <- "link"
  p <- NextMethod()
  p <- .circ_name_response(p, object$family, type)
  if (identical(type, "response"))            # link/terms/lpmatrix stay in fit frame
    p <- .circ_rotate_response(p, object$family, object$circ_center)
  p
}

#' @rdname circ_gam
#' @export
fitted.circ_gam <- function(object, ...) {
  p <- NextMethod()
  p <- .circ_name_response(p, object$family, "response")
  .circ_rotate_response(p, object$family, object$circ_center)
}

## Name the columns of a response-scale fitted matrix; pass everything else
## (link/terms/lpmatrix predictions, vectors) through untouched.
.circ_name_response <- function(p, family, type) {
  if (is.null(type)) type <- "link"
  pnm <- family$param_names
  if (identical(type, "response") && !is.null(pnm) &&
      is.matrix(p) && ncol(p) == length(pnm))
    colnames(p) <- pnm
  p
}

#' @rdname circ_gam
#' @export
print.circ_gam <- function(x, ...) {
  fam <- x$family
  if (!is.null(fam$param_names)) {
    g    <- .circ_geometry(x)
    key  <- if (is.na(g$kind)) "" else g$kind
    head <- switch(key,
      cl = "Circular GAM (circular-linear)",
      cc = "Circular GAM (circular-circular)",
      lc = "Linear-circular GAM",
      ll = "Location-scale GAM",
      if (isTRUE(g$resp_circular)) "Circular GAM" else "Location-scale GAM")
    cat(sprintf("%s via circ_gam() -- family %s, parameters: %s\n",
                head, fam$family, paste(fam$param_names, collapse = ", ")))
    if (!is.null(x$circ_center) && isTRUE(x$circ_center != 0))
      cat(sprintf(
        "  centered at %+.4g rad for fitting; directions reported in original frame\n",
        x$circ_center))
  }
  NextMethod()
}
