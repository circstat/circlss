# Bundled example datasets — sources and attribution

These CSV extracts back the `circlss` articles (real-data examples). They
are redistributed from the **NPCirc** R package (License: GPL ≥ 2;
`circlss` is GPL-3, so redistribution is license-compatible). Only the
columns used by the articles are kept; values are otherwise unchanged.

## `periwinkles.csv`
**Distance** (cm) and **direction** (degrees) of movement of small blue
periwinkles after relocation downshore — the classic circular–linear
regression dataset of Fisher & Lee (1992). Used as the clean
circular-response-on-a-covariate example (direction of movement vs
displacement distance), in well-concentrated data.

- Source: NPCirc package, dataset `periwinkles`.
- Original: Fisher, N. I. & Lee, A. J. (1992). Regression models for an
  angular response. *Biometrics* **48**, 665–677.

## `sandhoppers.csv`
Escape **orientation** `angle` (radians) of two sandhopper species
(*Talitrus saltator*, *Talorchestia brito*) measured on Zouara beach
(Tunisia) under natural conditions, with **temperature** (°C), sun
**azimuth** (deg), **humidity** (%), `species` and `sex`. Sandhoppers
orient toward the sea by a sun compass; the dataset shows the *precision*
(concentration) of that orientation varying with conditions — the
distributional payoff. 1828 observations.

- Source: NPCirc package, dataset `sandhoppers`.
- Original: Marchetti, G. M. & Scapini, F. (2003); Scapini, F. et al.
  (2002). Behavioural plasticity in sandhoppers.

## `pm10.csv`
Wind **direction** (degrees) with the **PM10** particulate level and wind
**speed** — used as the multi-covariate shape *contrast* (a circular law
flatter than von Mises, opposite the sandhopper case).

- Source: NPCirc package, dataset `pm10`.
- Original: Oliveira, M., Crujeiras, R. M. & Rodríguez-Casal, A. (2013).
  Nonparametric circular methods for exploring environmental data.
  *Environmental and Ecological Statistics* **20**, 1–17.
