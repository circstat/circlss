# Flight orientation of nocturnal migrating songbirds

Flight orientation of nocturnal passerine migrants, tracked by X-band
radar at Falsterbo, Sweden over the autumn migrations of 2009–2011
(Sjoberg & Nilsson 2015). Each row is one radar track that carries all
three of orientation, altitude, and wind direction – the tidy per-bird
table the additive
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
fit `orient ~ s(altitude) + s(wind)` needs. It is the joint covariate
set Ameijeiras-Alonso & Gijbels (2025) studied one covariate at a time;
an additive GAM fits both at once. Descending tracks are excluded (none
carry all three covariates), leaving 224 climbing and 910 level tracks.

## Usage

``` r
songbirds
```

## Format

A data frame with 1134 rows and 5 columns:

- orient:

  Flight orientation, the circular response, in radians on \\\[0,
  2\pi)\\ (the compass direction the bird flew).

- altitude:

  Altitude above ground, in metres.

- wind:

  Wind direction encountered aloft (the direction the wind blew *from*),
  in radians on \\\[-\pi, \pi\]\\.

- flight:

  Flight phase, a factor with levels `"climbing"` and `"level"`,
  classified from `vspeed` (climbing if \\\> 0.75\\ m/s, level if within
  \\\pm 0.75\\ m/s).

- vspeed:

  Vertical speed, in m/s (the basis for `flight`).

## Source

Sjoberg, S. & Nilsson, C. (2015) data on Dryad,
[doi:10.5061/dryad.86020](https://doi.org/10.5061/dryad.86020) .

## References

Sjoberg, S. & Nilsson, C. (2015) Nocturnal migratory songbirds adjust
their travelling direction aloft: evidence from a radiotelemetry and
radar study. *Biology Letters* **11**, 20150337.

Ameijeiras-Alonso, J. & Gijbels, I. (2025) Semiparametric regression for
circular response with application in ecology. *Scandinavian Journal of
Statistics* **53**, 54–101.

## Examples

``` r
## Additive multi-covariate fit on the level-flight tracks: orientation on both
## altitude and wind direction, for the mean and the concentration alike.
level <- subset(songbirds, flight == "level")
circ_gam(
  list(orient ~ s(altitude) + s(wind, bs = "cc"),
            ~ s(altitude) + s(wind, bs = "cc")),
  data = level, family = vmlss(), knots = list(wind = c(-pi, pi))
)
```
