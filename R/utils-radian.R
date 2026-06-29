#' Radian conversion helpers
#'
#' Small helpers for the \code{\link{circ_gam}} data contract: the response, and
#' any angular covariate entering a cyclic smooth, must be in radians. \code{rad}
#' converts to radians from a named or numeric period, \code{deg} converts radians
#' to degrees, and \code{wrap} folds angles to a single branch with
#' \code{atan2(sin, cos)} (exact and NA-safe).
#'
#' @param x A numeric vector of angles.
#' @param from For \code{rad}, the period of the input scale: one of
#'   \code{"degrees"} (360), \code{"gradians"} (400), \code{"turns"} (1),
#'   \code{"radians"} (a no-op), or any positive numeric period -- e.g. \code{24}
#'   for hour-of-day, \code{12} for month.
#' @param to For \code{wrap}, the target branch: \code{"pi"} for \eqn{(-\pi, \pi]}
#'   (the tan-half branch) or \code{"2pi"} for \eqn{[0, 2\pi)}.
#' @return A numeric vector the same length as \code{x}.
#' @examples
#' rad(c(0, 90, 180, 270))     # degrees -> radians
#' rad(c(0, 6, 12, 18), 24)    # hour-of-day -> radians
#' deg(pi)                     # 180
#' wrap(3 * pi)                # pi
#' wrap(-0.5, "2pi")           # 2*pi - 0.5
#' @seealso
#' \code{\link{circ_gam}}
#' @export
rad <- function(x, from = "degrees") {
  period <- if (is.numeric(from)) {
    if (length(from) != 1L || !is.finite(from) || from <= 0)
      stop("'from' must be a positive period or one of \"degrees\", ",
           "\"gradians\", \"turns\", \"radians\".")
    from
  } else switch(from,
                degrees = 360, gradians = 400, turns = 1, radians = 2 * pi,
                stop("unknown 'from' = \"", from, "\"; give a numeric period or ",
                     "one of \"degrees\", \"gradians\", \"turns\", \"radians\"."))
  x / period * 2 * pi
}

#' @rdname rad
#' @export
deg <- function(x) x * 180 / pi

#' @rdname rad
#' @export
wrap <- function(x, to = c("pi", "2pi")) {
  to <- match.arg(to)
  w <- atan2(sin(x), cos(x))            # (-pi, pi]
  if (to == "2pi") w <- w %% (2 * pi)   # [0, 2*pi)
  w
}
