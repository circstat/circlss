#' Movements of small blue periwinkles
#'
#' Direction and distance moved by 31 small blue periwinkles
#' (\emph{Nodilittorina unifasciata}) after they were transplanted downshore from
#' the height at which they normally live. This is Fisher's (1993) data set B.20 --
#' the worked example threaded through his Section 6.4 on circular-response
#' regression (Examples 6.3 and 6.9--6.11). It is bundled here with the radian and
#' mean-centred columns the \code{\link{circ_lm}} fits use, so no manual conversion
#' is needed.
#'
#' @format A data frame with 31 rows and 4 columns:
#' \describe{
#'   \item{theta.deg}{Direction moved, in degrees (compass bearing; the sea lay at
#'     roughly \eqn{275^\circ}).}
#'   \item{x}{Distance moved, in metres.}
#'   \item{theta.rad}{\code{theta.deg} in radians (\code{theta.deg * pi / 180}).}
#'   \item{xc}{\code{x} centred at its mean.}
#' }
#' @source Fisher, N. I. (1993) \emph{Statistical Analysis of Circular Data}.
#'   Cambridge University Press, Appendix B.20 (pp. 252--253). 
#' @examples
#' ## Fisher-Lee mixed model: mean and log-kappa both linear in distance (Ex. 6.11)
#' circ_lm(list(theta.rad ~ x, ~ x), data = periwinkles, type = "cl")
"periwinkles"

#' Flight orientation of nocturnal migrating songbirds
#'
#' Flight orientation of nocturnal passerine migrants, tracked by X-band radar at
#' Falsterbo, Sweden over the autumn migrations of 2009--2011 (Sjoberg & Nilsson 2015).
#' Each row is one radar track that carries all three of orientation, altitude, and wind
#' direction -- the tidy per-bird table the additive \code{\link{circ_gam}} fit
#' \code{orient ~ s(altitude) + s(wind)} needs. It is the joint covariate set
#' Ameijeiras-Alonso & Gijbels (2025) studied one covariate at a time; an additive GAM
#' fits both at once. Descending tracks are excluded (none carry all three covariates),
#' leaving 224 climbing and 910 level tracks.
#'
#' @format A data frame with 1134 rows and 5 columns:
#' \describe{
#'   \item{orient}{Flight orientation, the circular response, in radians on
#'     \eqn{[0, 2\pi)} (the compass direction the bird flew).}
#'   \item{altitude}{Altitude above ground, in metres.}
#'   \item{wind}{Wind direction encountered aloft (the direction the wind blew
#'     \emph{from}), in radians on \eqn{[-\pi, \pi]}.}
#'   \item{flight}{Flight phase, a factor with levels \code{"climbing"} and
#'     \code{"level"}, classified from \code{vspeed} (climbing if \eqn{> 0.75} m/s,
#'     level if within \eqn{\pm 0.75} m/s).}
#'   \item{vspeed}{Vertical speed, in m/s (the basis for \code{flight}).}
#' }
#' @source Sjoberg, S. & Nilsson, C. (2015) data on Dryad,
#'   \doi{10.5061/dryad.86020}.
#' @references
#' Sjoberg, S. & Nilsson, C. (2015) Nocturnal migratory songbirds adjust their
#' travelling direction aloft: evidence from a radiotelemetry and radar study.
#' \emph{Biology Letters} \strong{11}, 20150337.
#'
#' Ameijeiras-Alonso, J. & Gijbels, I. (2025) Semiparametric regression for circular
#' response with application in ecology. \emph{Scandinavian Journal of Statistics}
#' \strong{53}, 54--101.
#' @examples
#' ## Additive multi-covariate fit on the level-flight tracks: orientation on both
#' ## altitude and wind direction, for the mean and the concentration alike.
#' level <- subset(songbirds, flight == "level")
#' circ_gam(
#'   list(orient ~ s(altitude) + s(wind, bs = "cc"),
#'             ~ s(altitude) + s(wind, bs = "cc")),
#'   data = level, family = vmlss(), knots = list(wind = c(-pi, pi))
#' )
"songbirds"

#' Wind direction at a South African wind farm (ten-minute records)
#'
#' Wind measurements recorded every ten minutes at a wind farm in the Eastern Cape, South
#' Africa, throughout January 2019. The circular response -- wind direction -- is strongly
#' \strong{bimodal}, with an easterly and a westerly regime (a sea-breeze / land-breeze
#' signature whose mix shifts with time of day and wind speed), which makes it a natural
#' showcase for a \emph{mixture} of circular regressions (\code{\link{circ_mix}}) and for
#' the smooth cylinder/torus geometry of \code{\link{circ_gam}}. This is the raw series of
#' Skhosana & Nakhaei Rad (2026); they aggregate it to hourly means before fitting, which
#' you can reproduce with \code{aggregate()} on the hour of \code{ts} (circular mean of
#' \code{wd}, arithmetic mean of the rest).
#'
#' @format A data frame with 4464 rows and 6 columns:
#' \describe{
#'   \item{ts}{Timestamp of the measurement, a \code{POSIXct} (UTC), in 10-minute steps.}
#'   \item{wd}{Wind direction, the circular response, in radians on \eqn{[0, 2\pi)}
#'     (compass bearing the wind blew \emph{from}; \eqn{\pi/2} = E, \eqn{3\pi/2} = W).}
#'   \item{ws}{Wind speed, in m/s.}
#'   \item{tair}{Air temperature, in degrees Celsius.}
#'   \item{rh}{Relative humidity, in percent.}
#'   \item{tod}{Time of day, in radians on \eqn{[0, 2\pi)} (the diurnal circular
#'     covariate; \code{tod = 0} at midnight, \eqn{\pi} at noon).}
#' }
#' @source The authors' repository \url{https://github.com/Sphiwe-Skhosana/MixCircReg}
#'   (\code{wind_data.csv}), accompanying Skhosana & Nakhaei Rad (2026).
#' @references
#' Skhosana, S. & Nakhaei Rad, N. (2026) Model-based clustering using a new mixture of
#' circular regressions. \emph{arXiv:2601.05345}.
#' @examples
#' ## Smooth von Mises GAMs (the geometry the wind showcase draws):
#' ## a cylinder -- direction vs wind speed (circular ~ linear)
#' circ_gam(wd ~ s(ws), data = windfarm, family = vmlss())
#'
#' ## a torus -- direction vs time of day (circular ~ circular)
#' circ_gam(wd ~ s(tod, bs = "cc"), data = windfarm, family = vmlss(),
#'          knots = list(tod = c(0, 2 * pi)))
"windfarm"
