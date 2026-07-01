## ============================================================================
## Size-aware MAP-penalty degeneracy guard for circ_mix (every *lss family)
## ============================================================================
##
## A finite mixture's likelihood is unbounded: a component raises it by
## concentrating onto a responsibility-weighted subset -- its concentration runs
## off (kappa -> Inf, soft) or a bounded shape/peakedness parameter is driven to
## its singular boundary (rho -> 1/2, nu -> +/-1, the Kato-Jones shape disc edge),
## where the Hessian blows up and gam.fit5 reports an indefinite penalized
## likelihood (a hard crash) or grinds (a hang). A single circ_gam on the full
## data never sees this (every weight is 1); EM reweighting is what unanchors the
## concentration -- the same model as a circ_mix component reaches kappa ~ 5600
## where the standalone fit tops out near 220.
##
## The guard adds to each reweighted M-step a penalty pulling the
## degeneracy-prone parameter(s) toward the family's diffuse / reduced model, with
## strength lambda_k = c / N_k that VANISHES as the component grows (N_k = the
## effective component size sum_i gamma_ik). It is the dual of c diffuse
## pseudo-observations, worth ~c out of N_k -- so a well-populated component is
## unaffected (data dominates, smooths fit normally under REML) and a collapsing
## one cannot drive its parameter to the boundary. One regularizer kills all three
## symptoms (soft kappa->Inf, the boundary crash, the slow-grind hang) because
## they are one root.
##
## The penalty is expressed in the family's NATURAL parameters; the existing
## gamlss.etamu link chain-rule then carries it to eta, so families never touch
## chain-rule bookkeeping. It is folded into the per-row derivative blocks BEFORE
## the prior-weight (responsibility) scaling, so it is automatically gamma-weighted
## -- giving A1(kappa_hat) = R_bar_w - c/N_k for the intercept-kappa M-step.
## ret$l0 (the per-observation log-density the E-step reads) is left untouched, so
## responsibilities, mixture logLik and BIC stay on the data scale.

## ---- penalty kernels rho(v) on a single natural parameter v -----------------
## Each returns a spec list(index, rho0, rho1, rho2): the kernel value and its
## first two derivatives wrt the natural parameter. Only orders 0-2 are ever
## needed: the available.derivs = 2 (Newton, l3/l4) families use linear/ridge
## kernels with rho3 = rho4 = 0, and every boundary kernel sits on an EFS
## (available.derivs = 0, deriv <= 1) family. A small epsilon floors the
## boundary-kernel denominators so a parameter sitting exactly on the wall yields a
## large-but-finite penalty rather than Inf/NaN.

## Each kernel takes an intrinsic `scale` multiplier (default 1). The user dials
## one global strength `c` (control$degen_strength) that becomes lambda_k = c/N_k;
## `scale` lets a family declare a parameter's *relative* sensitivity, so c = 1
## stays the universal default. It is needed because the strength a degeneracy
## demands is parameter-specific, not kernel-type-specific: the Jones-Pewsey shape
## psi and its expensive quadrature normalizer want a much firmer prior (to stay
## fast and non-singular) than the Kato-Jones disc coordinates, though both ride a
## ridge. A family that needs no special handling simply omits `scale`.

## linear kernel rho(v) = scale*v: an exponential prior toward 0. For an
## unbounded-above concentration (kappa, log link). Caps via the constant gradient
## pull -scale*lambda; no curvature (rho2 = 0), which is right where the data
## Hessian flattens (kappa -> Inf) rather than going singular.
.degen_linear <- function(index, scale = 1)
  list(index = as.integer(index),
       rho0 = function(v) scale * v,
       rho1 = function(v) rep_len(scale, length(v)),
       rho2 = function(v) numeric(length(v)))

## ridge kernel rho(v) = scale*v^2: a Gaussian prior toward 0, for an unbounded
## shape coordinate whose diffuse value is 0 (kjlss u1/u2 -> wrapped Cauchy, pnlss
## mu1/mu2 -> uniform along the radius, jplss psi -> von Mises). rho2 = 2*scale
## adds positive curvature that also stabilizes the Hessian.
.degen_ridge <- function(index, scale = 1)
  list(index = as.integer(index),
       rho0 = function(v) scale * v * v,
       rho1 = function(v) 2 * scale * v,
       rho2 = function(v) rep_len(2 * scale, length(v)))

## one-sided boundary kernel rho(v) = -scale*log(1 - v/vmax) on v in [0, vmax):
## diverges as v -> vmax, so a small component cannot push the parameter onto the
## wall (cardlss rho -> 1/2). Pulls toward 0 (uniform).
.degen_boundary_upper <- function(index, vmax, scale = 1) {
  vmax <- as.numeric(vmax)
  list(index = as.integer(index),
       rho0 = function(v) { z <- pmax(1 - v / vmax, .Machine$double.eps); -scale * log(z) },
       rho1 = function(v) { z <- pmax(1 - v / vmax, .Machine$double.eps); scale * (1 / vmax) / z },
       rho2 = function(v) { z <- pmax(1 - v / vmax, .Machine$double.eps); scale * (1 / vmax^2) / (z * z) })
}

## two-sided boundary kernel rho(v) = -scale*log(1 - (v/vmax)^2) on v in
## (-vmax, vmax): diverges as v -> +/-vmax, pulling toward 0 (the von Mises /
## symmetric member). For the tanh-link peakedness / skewness parameters (vmftlss
## nu, ajplss nu, ssjplss lambda, ibslss nu & lambda) whose crash is the +/-1 edge.
.degen_boundary_sym <- function(index, vmax, scale = 1) {
  vmax <- as.numeric(vmax)
  list(index = as.integer(index),
       rho0 = function(v) { z <- pmax(1 - (v / vmax)^2, .Machine$double.eps); -scale * log(z) },
       rho1 = function(v) { z <- pmax(1 - (v / vmax)^2, .Machine$double.eps); scale * (2 * v / vmax^2) / z },
       rho2 = function(v) { z <- pmax(1 - (v / vmax)^2, .Machine$double.eps)
                            scale * ((2 / vmax^2) / z + (4 * v^2 / vmax^4) / (z * z)) })
}

## ---- the shared penalty assembler -------------------------------------------
## Given the fitted natural parameters (n x nlp, in param_names order), the
## family's `degen` spec, and lambda = c / N_k, return the penalty's contribution
## to (l0, l1, l2) in NATURAL coordinates, ready to ADD to the family's per-row
## blocks before the prior-weight scaling. Negation and lambda-scaling are baked
## in; placement into the packed l2 column uses the family's own trind index, so a
## family supplies only `degen` and the cbind of its natural parameters. The
## penalty is separable across the listed parameters, so it touches the l1
## gradient entries and the l2 DIAGONAL (i2[j,j]) only -- no cross-parameter
## blocks.
.lss_map_penalty <- function(family, natparams, lambda) {
  spec <- family$degen
  nlp  <- family$nlp
  i2   <- family$tri$i2
  n    <- nrow(natparams)
  dl0  <- numeric(n)
  dl1  <- matrix(0, n, nlp)
  dl2  <- matrix(0, n, nlp * (nlp + 1L) / 2L)
  for (s in spec) {
    j <- s$index
    v <- natparams[, j]
    dl0      <- dl0      - lambda * s$rho0(v)
    dl1[, j] <- dl1[, j] - lambda * s$rho1(v)
    cjj      <- i2[j, j]
    dl2[, cjj] <- dl2[, cjj] - lambda * s$rho2(v)
  }
  list(l0 = dl0, l1 = dl1, l2 = dl2)
}

## Is the degeneracy penalty active for this fit? TRUE only when circ_mix's M-step
## has set a positive map_lambda AND the family declares a degen spec; NULL/absent
## for a standalone circ_gam, so a non-mixture fit is byte-for-byte unchanged.
.degen_active <- function(family)
  !is.null(family$degen) && !is.null(family$map_lambda) &&
    is.finite(family$map_lambda) && family$map_lambda > 0
