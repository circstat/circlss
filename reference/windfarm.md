# Wind direction at a South African wind farm (ten-minute records)

Wind measurements recorded every ten minutes at a wind farm in the
Eastern Cape, South Africa, throughout January 2019. The circular
response – wind direction – is strongly **bimodal**, with an easterly
and a westerly regime (a sea-breeze / land-breeze signature whose mix
shifts with time of day and wind speed), which makes it a natural
showcase for a *mixture* of circular regressions
([`circ_mix`](https://circstat.github.io/circlss/reference/circ_mix.md))
and for the smooth cylinder/torus geometry of
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md).
This is the raw series of Skhosana & Nakhaei Rad (2026); they aggregate
it to hourly means before fitting, which you can reproduce with
[`aggregate()`](https://rdrr.io/r/stats/aggregate.html) on the hour of
`ts` (circular mean of `wd`, arithmetic mean of the rest).

## Usage

``` r
windfarm
```

## Format

A data frame with 4464 rows and 6 columns:

- ts:

  Timestamp of the measurement, a `POSIXct` (UTC), in 10-minute steps.

- wd:

  Wind direction, the circular response, in radians on \\\[0, 2\pi)\\
  (compass bearing the wind blew *from*; \\\pi/2\\ = E, \\3\pi/2\\ = W).

- ws:

  Wind speed, in m/s.

- tair:

  Air temperature, in degrees Celsius.

- rh:

  Relative humidity, in percent.

- tod:

  Time of day, in radians on \\\[0, 2\pi)\\ (the diurnal circular
  covariate; `tod = 0` at midnight, \\\pi\\ at noon).

## Source

The authors' repository <https://github.com/Sphiwe-Skhosana/MixCircReg>
(`wind_data.csv`), accompanying Skhosana & Nakhaei Rad (2026).

## References

Skhosana, S. & Nakhaei Rad, N. (2026) Model-based clustering using a new
mixture of circular regressions. *arXiv:2601.05345*.

## Examples

``` r
## Smooth von Mises GAMs (the geometry the wind showcase draws):
## a cylinder -- direction vs wind speed (circular ~ linear)
circ_gam(wd ~ s(ws), data = windfarm, family = vmlss())

## a torus -- direction vs time of day (circular ~ circular)
circ_gam(wd ~ s(tod, bs = "cc"), data = windfarm, family = vmlss(),
         knots = list(tod = c(0, 2 * pi)))
```
