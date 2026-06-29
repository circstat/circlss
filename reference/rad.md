# Radian conversion helpers

Small helpers for the
[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)
data contract: the response, and any angular covariate entering a cyclic
smooth, must be in radians. `rad` converts to radians from a named or
numeric period, `deg` converts radians to degrees, and `wrap` folds
angles to a single branch with `atan2(sin, cos)` (exact and NA-safe).

## Usage

``` r
rad(x, from = "degrees")

deg(x)

wrap(x, to = c("pi", "2pi"))
```

## Arguments

- x:

  A numeric vector of angles.

- from:

  For `rad`, the period of the input scale: one of `"degrees"` (360),
  `"gradians"` (400), `"turns"` (1), `"radians"` (a no-op), or any
  positive numeric period – e.g. `24` for hour-of-day, `12` for month.

- to:

  For `wrap`, the target branch: `"pi"` for \\(-\pi, \pi\]\\ (the
  tan-half branch) or `"2pi"` for \\\[0, 2\pi)\\.

## Value

A numeric vector the same length as `x`.

## See also

[`circ_gam`](https://circstat.github.io/circlss/reference/circ_gam.md)

## Examples

``` r
rad(c(0, 90, 180, 270))     # degrees -> radians
rad(c(0, 6, 12, 18), 24)    # hour-of-day -> radians
deg(pi)                     # 180
wrap(3 * pi)                # pi
wrap(-0.5, "2pi")           # 2*pi - 0.5
```
