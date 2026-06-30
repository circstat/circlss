# Movements of small blue periwinkles

Direction and distance moved by 31 small blue periwinkles
(*Nodilittorina unifasciata*) after they were transplanted downshore
from the height at which they normally live. This is Fisher's (1993)
data set B.20 – the worked example threaded through his Section 6.4 on
circular-response regression (Examples 6.3 and 6.9–6.11). It is bundled
here with the radian and mean-centred columns the
[`circ_lm`](https://circstat.github.io/circlss/reference/circ_lm.md)
fits use, so no manual conversion is needed.

## Usage

``` r
periwinkles
```

## Format

A data frame with 31 rows and 4 columns:

- theta.deg:

  Direction moved, in degrees (compass bearing; the sea lay at roughly
  \\275^\circ\\).

- x:

  Distance moved, in metres.

- theta.rad:

  `theta.deg` in radians (`theta.deg * pi / 180`).

- xc:

  `x` centred at its mean.

## Source

Fisher, N. I. (1993) *Statistical Analysis of Circular Data*. Cambridge
University Press, Appendix B.20 (pp. 252–253).

## Examples

``` r
## Fisher-Lee mixed model: mean and log-kappa both linear in distance (Ex. 6.11)
circ_lm(list(theta.rad ~ x, ~ x), data = periwinkles, type = "cl")
```
