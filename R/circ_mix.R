#' Finite mixtures of circular distributional GAMs by EM
#'
#' Fits a \eqn{K}-component finite mixture of circular distributional GAMs by the
#' EM algorithm. It does not touch the families or \pkg{mgcv} internals: each M-step
#' is a weighted \code{\link{circ_gam}} fit and each E-step reads the family's
#' per-observation density. Because a component
#' is reached only through that small interface, one engine spans density
#' clustering, the circular--linear / circular--circular / linear--circular
#' regression trio, and everything in between -- the response geometry is set
#' entirely by \code{family}.
#'
#' @param formula A model formula for one component, exactly as
#'   \code{\link{circ_gam}} expects: \code{theta ~ 1} (density clustering),
#'   \code{theta ~ x1 + x2} (circular--linear regression), \code{theta ~ cos(phi) +
#'   sin(phi)} or \code{theta ~ s(phi, bs = "cc")} (circular--circular), or
#'   \code{y ~ s(phi, bs = "cc")} with a linear-response family (linear--circular).
#'   A \emph{list} of formulas carrying two or more distinct responses fits a joint
#'   (torus) density by the chain rule -- e.g.
#'   \code{list(psi ~ cos(phi) + sin(phi), phi ~ 1)} factorises
#'   \eqn{f(\psi,\phi) = f(\psi \mid \phi)\,f(\phi)} into two factors; a response
#'   named in another formula's right-hand side is conditioned on it, and the list
#'   order is the chain-rule order. (A single response with two or more
#'   \emph{location-scale} predictors is still one component, written
#'   \code{list(theta ~ s(x), ~ s(x))}; the joint reading needs two or more distinct
#'   left-hand sides. Joint densities over more than two responses are not yet
#'   supported.)
#' @param data A data frame holding the response and covariates.
#' @param family A circlss location-scale family object (one carrying
#'   \code{param_names}): any circular family (\code{\link{vmlss}},
#'   \code{\link{pnlss}}, \dots) for a circular response, or \code{\link{gausslss}}
#'   / \code{\link{gammalss}} for the linear-circular leg. The family is the only
#'   thing that sets the response geometry; the EM machinery is identical across
#'   all of them.
#' @param K The number of mixture components. Under \code{search = "fixed"} it is
#'   held; under \code{search = "greedy"} it is the \emph{starting} count from which
#'   the data grows or shrinks \eqn{K} (not a ceiling or floor); it is ignored under
#'   \code{search = "grid"} (which sweeps \code{control$kmin:control$kmax}).
#' @param search How the number of components is decided. \code{"fixed"} (default)
#'   holds \code{K}. Automatic-K search is opt-in: \code{"greedy"} runs bidirectional
#'   split / merge / death moves from the init \code{K}, accepting any move that lowers
#'   the penalised objective \eqn{J = -2\,\mathrm{logLik} + \lambda\,\mathrm{df}}
#'   (\eqn{= } BIC when \code{lambda = log(n)}) -- the warm heuristic, which grows
#'   reliably from a small init \code{K}; \code{"grid"} fits every \code{K} in
#'   \code{control$kmin:control$kmax} with restarts and picks the minimum-\eqn{J}
#'   \code{K} -- the robust, embarrassingly-parallel selector and the cross-check on the
#'   moves.
#' @param assign The E-step assignment rule. \code{"soft"} (default) is EM with
#'   fractional responsibilities. \code{"hard"} is classification EM (CEM): each unit
#'   seats wholly at its argmax component, and the engine maximises the classification
#'   log-likelihood \eqn{\sum_u \max_k (\log\pi_k + \log f_k(y_u))} rather than the
#'   mixture log-likelihood (its \code{logLik}/BIC are on that classification scale).
#'   Combined with \code{search = "greedy"} it is the circular DP-means / k-means-style
#'   hard clustering.
#' @param group The clustering unit. \code{NULL} (default) clusters \emph{rows} --
#'   one responsibility per observation. A one-sided formula naming a grouping variable,
#'   \code{group = ~ id}, clusters \emph{subjects / curves}: a subject's whole trajectory
#'   seats at one component (the longitudinal / latent-class-growth case). Under a group
#'   the responsibilities, MAP labels and the BIC sample size are all per subject.
#' @param weights_model Reserved for covariate-driven mixing weights (not yet
#'   implemented); must be \code{~ 1} (constant mixing proportions).
#' @param control A list of tuning parameters from \code{circ_mix.control()}.
#' @param object,x A fitted \code{circ_mix} model.
#' @param newdata A data frame of new data. For \code{predict}, omitting it uses
#'   the training frame.
#' @param type For \code{predict}, what to return: \code{"cluster"} (the
#'   \eqn{n \times K} responsibility matrix, the default), \code{"density"} (the
#'   mixture density per row) or \code{"response"} (per-component response-scale
#'   fitted curves, a list of length \code{K}).
#' @param log For \code{predict(type = "density")}, return the log-density.
#' @param ... Further arguments forwarded to the per-component
#'   \code{\link{circ_gam}} / \code{\link[mgcv]{gam}} M-step (\code{knots},
#'   \code{method}, \dots).
#' @return
#' An object of class \code{"circ_mix"} with, among others:
#' \item{components}{the list of \code{K} fitted components -- each wrapping a
#' weighted \code{\link{circ_gam}}, or, for a joint density, a \emph{product} of
#' several (one per chain-rule factor).}
#' \item{gating}{the mixing-weight object; \code{gating$pi} are the component
#' proportions.}
#' \item{gamma}{the \eqn{n \times K} matrix of soft responsibilities.}
#' \item{cluster}{the per-observation MAP cluster labels.}
#' \item{loglik, df, bic}{the mixture log-likelihood, degrees of freedom
#' \eqn{(K-1) + \sum_k \mathrm{edf}_k}, and BIC.}
#' \item{converged, iter, ll_path}{convergence flag, iteration count, and the
#' recorded log-likelihood path of the selected run.}
#' \item{K_init, search}{the starting component count and the search mode used.}
#' \item{move_trace}{for \code{search = "greedy"}, a data frame of the accepted
#' moves (\code{move}, \code{K_from}, \code{K_to}, \code{J_from}, \code{J_to}); for
#' \code{"grid"}, the \code{K}-by-\code{K} sweep (\code{loglik}, \code{df},
#' \code{bic}, \code{J}); \code{NULL} for \code{"fixed"}.}
#' \item{restarts}{the per-restart log-likelihoods and the basin-hit count.}
#' Methods provided: \code{print}, \code{summary}, \code{predict}, \code{coef},
#' \code{logLik} (so \code{AIC}/\code{BIC} work), and \code{\link{plot.circ_mix}}.
#' @details
#' \strong{EM and restarts.} Each run alternates a weighted M-step (one
#' \code{circ_gam} per component, weighted by the responsibilities) with an E-step
#' that records the observed-data mixture log-likelihood and updates the
#' responsibilities. For parametric (penalty-free) components the EM is monotone.
#' The fit is repeated from \code{restarts} random responsibility seeds and the
#' largest-log-likelihood run is kept; \code{$restarts$basin_hits} reports how many
#' restarts reached it (a health signal). For parametric components, and for smooth
#' components under \code{penalty = "fixed"}, the EM is monotone and a non-monotone
#' step is warned about; under \code{penalty = "auto"} a moving smoothing penalty
#' makes small dips expected, so they are not flagged.
#'
#' \strong{Model selection.} \code{df = (K-1) + sum_k edf_k} and
#' \code{BIC = -2 logLik + df * log(n)}. \code{logLik()} carries \code{df} and
#' \code{nobs}, so \code{AIC()} and \code{BIC()} work generically. For a joint
#' component the per-component \code{edf} is summed over its factors, so
#' \code{df = (K-1) + sum_k sum_j edf_{kj}}.
#'
#' \strong{Joint (torus) density.} A multi-response \code{formula} (two distinct
#' left-hand sides) makes each component a \emph{product} of weighted
#' \code{\link{circ_gam}} fits -- one per chain-rule factor -- whose joint
#' log-density is the sum of the factor log-densities. The EM loop, restarts, MAP,
#' \eqn{J}/BIC and the automatic-K moves are unchanged: the joint case is a
#' component-implementation swap, not a different engine. The circular k-means
#' initialisation (\code{\link{circ_kmeans}}) seeds on all angular responses
#' jointly -- one torus coordinate per response, e.g. \eqn{(\phi, \psi)} -- and a
#' greedy split divides the worst component on its joint angular residuals, so it
#' grows along whichever response is the more over-dispersed.
#'
#' \strong{Automatic K.} With \code{search = "greedy"} (opt-in) \code{K} is only
#' a starting count: after EM converges the engine attempts structure moves -- a
#' \emph{split} of the worst-fit component (2-means on its angular residuals; grows
#' \eqn{K}), a \emph{merge} of the two most similar components and a \emph{death} of
#' any below-\code{min_size} component (both shrink \eqn{K}) -- and accepts the move
#' that most lowers \eqn{J = -2\,\mathrm{logLik} + \lambda\,\mathrm{df}}, repeating
#' until no move improves \eqn{J}. Because \eqn{J} strictly decreases, the search is
#' monotone and cannot cycle; \code{$move_trace} records the accepted moves. The
#' greedy search grows reliably from a small init \code{K}; \code{search = "grid"}
#' (fit every \code{K} in \code{kmin:kmax} with restarts, pick the minimum-\eqn{J}) is
#' the robust selector and the recommended cross-check, especially when \code{K} may
#' be over-specified or the per-\code{K} optimum is hard to reach (raise
#' \code{restarts} there).
#'
#' \strong{Longitudinal / curve clustering.} With \code{group = ~ id} the unit is a
#' subject: the E-step sums each component's per-row log-densities within subject (so
#' the whole trajectory shares one responsibility), the M-step broadcasts that
#' responsibility back to the subject's rows in the weighted \code{\link{circ_gam}},
#' and \code{\link{logLik}} / BIC count \eqn{n =} the number of subjects. The trajectory
#' model is just the per-component formula (e.g. \code{theta ~ cos(t) + sin(t)}); this
#' is the circular latent-class-growth (LCGA) case, and it composes with automatic K
#' (moves reassign whole curves).
#'
#' \strong{Penalised smooths and hard assignment.} Smooth component terms
#' (\code{s(x)}, \code{s(phi, bs = "cc")}, factor-smooth \code{s(t, id, bs = "fs")}
#' random curves) fit through the same weighted M-step; the only subtlety is the
#' smoothing parameters, handled by \code{control$penalty} -- \code{"auto"} (default)
#' selects them by REML per component (automatic smoothness, the GAM way; the
#' per-iteration search is warm-started across iterations and driven by the extended
#' Fellner--Schall optimiser by default, so it stays cheap -- see \code{optimizer}),
#' \code{"fixed"} freezes one pooled value (faster and monotone, and a guard against a
#' component over-flexing to absorb another cluster). Their degrees of freedom enter the BIC as the \emph{effective}
#' \code{sum(edf)}, which mgcv reports at the fixed or selected smoothing parameters,
#' so the criterion stays well defined under regularisation. \code{assign = "hard"}
#' switches the E-step to classification EM (CEM); with \code{group = ~ id} it
#' clusters whole curves by hard assignment, and with \code{search = "greedy"} it is
#' the circular DP-means.
#'
#' \strong{Unsupported arguments.} An argument value that is not yet implemented
#' (a non-trivial \code{weights_model}) is validated and raises an informative
#' error rather than being silently ignored.
#' @examples
#' library(mgcv)
#' set.seed(1); n <- 400
#' z  <- sample.int(2L, n, replace = TRUE, prob = c(0.4, 0.6))
#' x  <- runif(n, -1, 1)
#' mu <- 2 * atan(c(0.9, -0.9)[z] + c(2.2, -2.2)[z] * x)
#' y  <- vmlss()$rd(cbind(mu, rep(6, n)), rep(1, n), 1)
#' \donttest{
#' ## a two-component von Mises regression on a linear covariate (c~l); fixed K (default)
#' m  <- circ_mix(y ~ x, data = data.frame(y, x), family = vmlss(), K = 2)
#' m
#' table(truth = z, cluster = m$cluster)
#' plot(m)
#'
#' ## automatic K (opt-in): grow/shrink from the init K by greedy moves
#' auto <- circ_mix(y ~ x, data = data.frame(y, x), family = vmlss(), K = 2,
#'                  search = "greedy")
#' auto$K          # the selected number of components
#' auto$move_trace # the accepted split/merge/death moves
#'
#' ## the robust cross-check: a brute K = 1..6 BIC sweep
#' g <- circ_mix(y ~ x, data = data.frame(y, x), family = vmlss(),
#'               search = "grid", control = circ_mix.control(kmax = 6))
#'
#' ## a joint (torus) density over two angles: each component is a PRODUCT of a
#' ## conditional f(psi | phi) and a marginal f(phi) -- compact 2-D blobs.
#' set.seed(2); n2 <- 300L
#' zz  <- sample.int(2L, n2, replace = TRUE)
#' phi <- vmlss()$rd(cbind(c(-2, 1)[zz], rep(4, n2)), rep(1, n2), 1)
#' psi <- vmlss()$rd(cbind(c(1, -1.5)[zz] + 0.8 * sin(phi), rep(5, n2)), rep(1, n2), 1)
#' j <- circ_mix(list(psi ~ cos(phi) + sin(phi), phi ~ 1),
#'               data = data.frame(psi, phi), family = vmlss(), K = 2)
#' j               # a "joint torus density" mixture; each component is a product
#' plot(j)         # the torus-square scatter, coloured by MAP cluster
#'
#' ## longitudinal: cluster whole CURVES, not rows (group = ~id). Two latent classes
#' ## of circular growth trajectory; each subject seats at one component.
#' set.seed(3); nsub <- 40L; nt <- 8L
#' cl   <- sample.int(2L, nsub, replace = TRUE)
#' long <- do.call(rbind, lapply(seq_len(nsub), function(s) {
#'   tt <- sort(runif(nt)); ph <- 2 * pi * tt; mu <- c(1.2, -1.2)[cl[s]] * sin(ph)
#'   data.frame(id = s, phase = ph, a = vmlss()$rd(cbind(mu, rep(10, nt)), rep(1, nt), 1))
#' }))
#' lc <- circ_mix(a ~ cos(phase) + sin(phase), data = long,
#'                family = vmlss(), K = 2, group = ~ id)
#' lc                                       # K components over 40 subjects
#' table(truth = cl, cluster = lc$cluster)  # cluster is per subject
#' }
#' @seealso
#' \code{\link{plot.circ_mix}}, \code{\link{circ_gam}},
#' \code{\link{vmlss}}, \code{\link{gausslss}}
#' @export
circ_mix <- function(formula, data, family = vmlss(),
                     K = 2,
                     search = c("fixed", "greedy", "grid"),
                     assign = c("soft", "hard"),
                     group = NULL,
                     weights_model = ~ 1,
                     control = circ_mix.control(), ...) {
  search <- match.arg(search)
  assign <- match.arg(assign)

  ## ---- scope guards: dormant arguments error until their phase lands --------
  if (length(all.vars(weights_model)))
    stop("weights_model (covariate-driven mixing weights) is not yet supported; ",
         "use weights_model = ~ 1.")
  if (missing(data))
    stop("circ_mix() requires 'data'.")
  if (is.null(family$param_names))
    stop("'family' must be a circlss location-scale family (one carrying ",
         "param_names), e.g. vmlss(), pnlss(), gausslss().")
  K <- as.integer(K)
  if (length(K) != 1L || is.na(K) || K < 1L)
    stop("K must be a single integer >= 1.")

  n     <- nrow(data)
  resps <- .circ_mix_responses(formula)         # >= 2 ==> a joint product component
  resp  <- resps[1L]                             # the primary (conditional) response
  if (length(resps) > 2L)
    stop("joint densities over more than two responses (d > 2) arrive in a later ",
         "phase; use a two-response (torus) spec for now.")
  miss <- resps[!resps %in% names(data)]
  if (length(miss))
    stop("response(s) not found in 'data': ", paste(miss, collapse = ", "), ".")

  ## ---- the clustering unit (contract 1): rows, or subjects via group = ~id ---
  u    <- .circ_mix_group_index(group, data, n)   # grp (per-row unit index), n_units
  grp  <- u$grp; n_units <- u$n_units
  if (K > n_units)
    stop("K = ", K, " exceeds the number of ", u$kind, "s (", n_units, ").")
  if (control$kmin > control$kmax)
    stop("control$kmin (", control$kmin, ") exceeds kmax (", control$kmax, ").")

  ## ---- fit: dispatch on the K-search strategy (contract 4) ------------------
  ## the penalty/IC count is the number of UNITS (rows, or subjects under group).
  ## Safety: a forked worker (cores > 1) that makes a LARGE BLAS call can segfault
  ## when R is linked against a threaded BLAS that is not fork-safe -- notably macOS
  ## Accelerate/vecLib, whose GCD-dispatched paths a fork() leaves holding stale
  ## dispatch handles. The trigger is BLAS operation SIZE crossing that threaded
  ## path, not "smooth vs parametric" per se; but penalised SMOOTH fits reliably
  ## cross it (mgcv's magic / gam.fit3 do heavy LAPACK), while parametric fits
  ## usually stay below -- a heuristic, not a guarantee. Whether a build is even
  ## exposed depends on R version / architecture / which BLAS it links (the default
  ## has been in flux), hence "some platforms". An unconditional serial fallback is
  ## the robust guard (VECLIB_MAXIMUM_THREADS doesn't bind Apple Silicon's AMX path).
  if (control$cores > 1L && .circ_mix_formula_smooth(formula)) {
    warning("circ_mix: not using parallel restarts (cores > 1) for smooth ",
            "(penalised) components -- a forked large BLAS call can crash a ",
            "non-fork-safe threaded BLAS (e.g. macOS Accelerate) on some builds; ",
            "running serially.", call. = FALSE)
    control$cores <- 1L
  }
  if (!is.null(control$seed)) set.seed(control$seed)

  ## ---- carry the E-step rule + penalty handling on the control --------------
  ## `assign` (soft EM vs hard CEM) rides the control so the EM core reads it
  ## without a new threaded argument. For penalty = "fixed"/"scheduled" the M-step
  ## smoothing parameters are seeded ONCE here, from a pooled single-component pilot
  ## fit on all the data -- "fixed" holds it for every component and iteration (so
  ## the penalty never moves and the EM is monotone in the penalised objective, the
  ## plan's "select lambda once" choice), "scheduled" starts from it and re-adapts
  ## periodically. It matters only for smooth components: a parametric formula has no
  ## penalty, so the three modes coincide and no pilot is needed. A user-supplied
  ## control$sp wins (the pilot is skipped).
  control$assign <- assign
  if (control$penalty %in% c("fixed", "scheduled") && is.null(control$sp) &&
      .circ_mix_formula_smooth(formula)) {
    control$sp <- tryCatch(.circ_mix_pooled_sp(formula, data, family, ...),
                           error = function(e) NULL)
    if (is.null(control$sp)) {
      warning("circ_mix: penalty = \"", control$penalty, "\" pilot fit failed; ",
              "falling back to penalty = \"auto\" (per-iteration REML).", call. = FALSE)
      control$penalty <- "auto"
    }
  }

  lam <- if (is.null(control$lambda)) log(n_units) else control$lambda
  res <- switch(search,
    fixed  = .circ_mix_search_fixed (formula, data, family, K, control, resp, grp, lam, ...),
    grid   = .circ_mix_search_grid  (formula, data, family, K, control, resp, grp, lam, ...),
    greedy = .circ_mix_search_greedy(formula, data, family, K, control, resp, grp, lam, ...))
  best <- res$state
  Kf   <- best$K
  if (isFALSE(best$monotone))
    warning("circ_mix: an EM run's log-likelihood was non-monotone (worst drop ",
            signif(best$worst_drop, 3L), "); parametric or fixed-penalty EM should ",
            "ascend -- inspect $ll_path (use penalty = \"fixed\" for monotone ",
            "smooth EM).", call. = FALSE)

  ## ---- assemble the (phase-stable) fitted object ----------------------------
  edf_k <- vapply(best$components, cmp_edf, numeric(1))
  df    <- (Kf - 1L) + sum(edf_k)
  bic   <- -2 * best$loglik + df * log(n_units)       # n = #units (subjects under group)
  J     <- -2 * best$loglik + df * lam
  nk    <- as.integer(tabulate(best$cluster, Kf))     # MAP cluster sizes (per unit)

  structure(list(
    call      = match.call(),
    formula   = formula,
    family    = family,
    response  = resp,
    K         = Kf,
    K_init    = K,
    search    = search,
    components = best$components,
    gating    = list(type = "constant", pi = best$pi),
    unit      = list(kind = u$kind, index = grp, n_units = n_units, labels = u$labels),
    gamma     = best$gamma,
    cluster   = best$cluster,
    nk        = nk,
    Gtilde    = sum(nk > 0L),
    loglik    = best$loglik,
    df        = df,
    edf       = edf_k,
    bic       = bic,
    objective = list(J = J, loglik = best$loglik, df = df, bic = bic, lambda = lam),
    iter      = best$iter,
    converged = best$converged,
    ll_path   = best$ll_path,
    monotone  = best$monotone,
    restarts  = list(R = best$R, lls = best$restart_lls, basin_hits = best$basin_hits),
    move_trace = res$trace,             # NULL (fixed) | accepted moves (greedy) | sweep (grid)
    geometry  = if (inherits(best$components[[1L]], "circ_mix_product")) "joint"
                else .circ_geometry(best$components[[1L]]$fit)$kind,
    control   = control
  ), class = "circ_mix")
}

#' @rdname circ_mix
#' @param lambda The penalty multiplier on the degrees of freedom in the model
#'   objective \eqn{J}; \code{NULL} uses \code{log(n_units)}, making \eqn{J} the BIC.
#' @param restarts Number of random-restart EM runs; the largest-log-likelihood
#'   run is kept.
#' @param init Initialisation of the responsibilities: \code{"kmeans"} (default;
#'   restart 1 seeds by clustering the response -- circular k-means
#'   (\code{\link{circ_kmeans}}) on the circle / torus for an angular response,
#'   ordinary \code{\link[stats]{kmeans}} for the linear \code{l~c} leg -- and later
#'   restarts are random) or \code{"random"}. \code{"emEM"} is not yet supported.
#' @param tol,max_iter EM convergence tolerance (relative change in the
#'   log-likelihood) and the maximum number of EM iterations per run.
#' @param wfloor Lower bound applied to the responsibilities used as M-step prior
#'   weights, keeping them strictly positive without perturbing the fit.
#' @param seed Optional integer seed set once before the restarts, for
#'   reproducibility. Each restart then seeds itself from a stream drawn here, so the
#'   result is identical whether the restarts run serially or in parallel.
#' @param cores Number of parallel workers for the embarrassingly-parallel restart
#'   runs (and each grid-\eqn{K}'s restarts). \code{1} (default) runs serially;
#'   \code{> 1} forks that many workers via \code{\link[parallel]{mclapply}} (no
#'   speed-up on Windows, which cannot fork -- it falls back to serial). Results do not
#'   depend on \code{cores}. A formula with penalised smooths always runs serially: a
#'   forked large BLAS call can crash a non-fork-safe threaded BLAS (e.g. macOS
#'   Accelerate) on some builds.
#' @param verbose If \code{TRUE}, report per-iteration progress.
#' @param kmin,kmax The lower and upper bounds on \eqn{K} for automatic-K search
#'   (both \code{"greedy"} and \code{"grid"}).
#' @param moves The structure moves the greedy search may attempt, any of
#'   \code{"split"} (grow), \code{"merge"} and \code{"death"} (shrink); \code{"birth"}
#'   (grow, redundant with split) is also accepted.
#' @param min_size The soft-size floor \eqn{n_k = \sum_i \gamma_{ik}} below which a
#'   component is dropped by a death move, and the minimum members a component must
#'   have to be split.
#' @param penalty How the per-component M-step handles the smoothing parameters of
#'   penalised (smooth) terms. \code{"auto"} (default) lets REML select them every M-step,
#'   so each component gets its own automatically-chosen smoothness -- the usual GAM
#'   behaviour; the trade-off is that the moving penalty makes the EM non-monotone for
#'   smooth components. \code{"fixed"} selects a single smoothness \emph{once}, from a
#'   pooled single-component pilot fit on all the data, and holds it for every component
#'   and iteration -- an opt-in for speed (no per-M-step REML search), a monotone EM, or
#'   robustness when \code{"auto"}'s per-component adaptivity lets a component over-flex
#'   and absorb a neighbouring cluster; the cost is one shared smoothness rather than a
#'   per-component one. \code{"scheduled"} starts from that same pooled value and
#'   re-selects every \code{sp_every} iterations \emph{after} the first. For parametric
#'   (penalty-free) components the three modes coincide. Under \code{"auto"} /
#'   \code{"scheduled"} each per-M-step REML search is \emph{warm-started} from the
#'   previous M-step's selected smoothing parameters (\pkg{mgcv}'s \code{in.out}), so
#'   it converges in a step or two without changing the value it converges to -- the
#'   per-component, per-iteration smoothness is unchanged, only reached faster.
#' @param optimizer The \pkg{mgcv} outer optimiser for that REML smoothing-parameter
#'   search (ignored by \code{penalty = "fixed"} and by parametric components, which
#'   run no search). \code{"efs"} (default) is the extended Fellner--Schall method
#'   (Wood & Fasiolo, 2017): for the families fitted by outer Newton
#'   (\code{available.derivs = 2}, e.g. \code{\link{vmlss}}, \code{\link{pnlss}}) it is
#'   markedly faster here and avoids the flat-REML step-failure warnings while selecting
#'   the same smoothness; for the \code{available.derivs = 0} families it is already
#'   \pkg{mgcv}'s default. Pass \code{"outer"} for outer Newton.
#' @param sp Optional smoothing parameters to hold fixed (a numeric vector for a
#'   single-response component, a per-factor list for a joint one), in
#'   \code{\link[mgcv]{gam}} order. When supplied it overrides the pilot fit under
#'   \code{penalty = "fixed"}; \code{NULL} (default) selects them automatically.
#' @param sp_every Under \code{penalty = "scheduled"}, the number of EM iterations
#'   between REML re-selections of the smoothing parameters (held fixed in between).
#' @param start,kappa_cap Accepted but not yet used (kappa_cap is reserved for
#'   concentration capping).
#' @export
circ_mix.control <- function(lambda = NULL,
                             kmin = 1L, kmax = 20L,
                             penalty = c("auto", "fixed", "scheduled"),
                             sp = NULL, optimizer = "efs",
                             start = NULL, sp_every = 5L,
                             restarts = 10L, init = "kmeans",
                             moves = c("split", "merge", "death"),
                             tol = 1e-6, max_iter = 200L,
                             min_size = 5L, kappa_cap = Inf,
                             wfloor = 1e-8, seed = NULL, cores = 1L,
                             verbose = FALSE) {
  init <- match.arg(init, c("kmeans", "random", "emEM"))
  if (init == "emEM")
    stop("init = \"emEM\" is not supported; use \"kmeans\" or \"random\".")
  penalty <- match.arg(penalty)
  if (!is.character(optimizer) || !length(optimizer))
    stop("optimizer must be an mgcv outer-optimiser string, e.g. \"efs\" or \"outer\".")
  moves <- match.arg(moves, c("split", "merge", "death", "birth"), several.ok = TRUE)
  list(lambda = lambda, kmin = as.integer(kmin), kmax = as.integer(kmax),
       penalty = penalty, sp = sp, optimizer = optimizer, start = start,
       sp_every = max(1L, as.integer(sp_every)),
       restarts = as.integer(restarts), init = init, moves = moves,
       tol = tol, max_iter = as.integer(max_iter),
       min_size = as.integer(min_size), kappa_cap = kappa_cap,
       wfloor = wfloor, seed = seed, cores = max(1L, as.integer(cores)),
       verbose = isTRUE(verbose))
}

## ============================================================================
## K-SEARCH STRATEGIES (contract 4: every comparison goes through .objective)
## ============================================================================

## the penalised objective J = -2 logLik + lambda * df, df = (K-1) + sum(edf).
## lambda = log(#units) makes J the BIC, so "let the data decide K" == "greedily
## minimise BIC by local moves". This is the single comparator the search uses.
.circ_mix_objective <- function(state, lambda) {
  df <- (state$K - 1L) + sum(vapply(state$components, cmp_edf, numeric(1)))
  -2 * state$loglik + lambda * df
}

## search = "fixed": the Phase-1 behaviour. K=1 is a single weighted-trivial
## circ_gam; K>1 is best-of-restarts. Returns {state, trace = NULL}.
.circ_mix_search_fixed <- function(formula, data, family, K, control, resp, grp, lam, ...) {
  if (K == 1L) {
    st <- .circ_mix_em_once(formula, data, family, 1L, control, resp,
                            init_method = "kmeans", grp = grp, ...)
    st$R <- 1L; st$restart_lls <- st$loglik; st$basin_hits <- 1L
  } else {
    st <- .circ_mix_restarts(formula, data, family, K, control, resp,
                             control$restarts, grp = grp, ...)
    st$R <- control$restarts
  }
  list(state = st, trace = NULL)
}

## search = "grid": the brute kmin:kmax sweep + IC -- a validation cross-check on
## the greedy moves, embarrassingly parallel in principle. Picks the min-J K.
## Returns {state, trace = data.frame(K, loglik, df, bic, J)}.
.circ_mix_search_grid <- function(formula, data, family, K, control, resp, grp, lam, ...) {
  nu <- max(grp)                                       # #units (rows, or subjects)
  Ks <- seq.int(control$kmin, control$kmax); Ks <- Ks[Ks <= nu]
  rows <- vector("list", length(Ks)); best <- NULL; bestJ <- Inf
  for (i in seq_along(Ks)) {
    k  <- Ks[i]
    fk <- tryCatch(.circ_mix_search_fixed(formula, data, family, k, control, resp, grp, lam, ...)$state,
                   error = function(e) { if (control$verbose)
                     cat(sprintf("  [grid] K=%d FAILED: %s\n", k, conditionMessage(e))); NULL })
    if (is.null(fk)) next
    J  <- .circ_mix_objective(fk, lam)
    df <- (fk$K - 1L) + sum(vapply(fk$components, cmp_edf, numeric(1)))
    rows[[i]] <- data.frame(K = k, loglik = fk$loglik, df = df,
                            bic = -2 * fk$loglik + df * log(nu), J = J)
    if (control$verbose)
      cat(sprintf("  [grid] K=%d  loglik=%.2f  df=%.1f  J=%.2f\n", k, fk$loglik, df, J))
    if (J < bestJ) { bestJ <- J; best <- fk }
  }
  if (is.null(best)) stop("grid search: all K in ", control$kmin, ":", control$kmax, " failed.")
  list(state = best, trace = do.call(rbind, Filter(Negate(is.null), rows)))
}

## search = "greedy": start at init K, run EM to convergence, then attempt
## structure moves; accept the move that drops J most, repeat until none improves
## J. Greedy + strict decrease => deterministic, monotone in J, and cannot cycle.
## Returns {state, trace = data.frame(step, move, K_from, K_to, J_from, J_to, dJ)}.
.circ_mix_search_greedy <- function(formula, data, family, K, control, resp, grp, lam, ...) {
  base  <- .circ_mix_search_fixed(formula, data, family, K, control, resp, grp, lam, ...)
  state <- base$state; J <- .circ_mix_objective(state, lam)
  cx <- list(formula = formula, data = data, family = family, control = control,
             resp = resp, grp = grp, lam = lam, resp_circ = .circ_mix_resp_circular(family))
  fns <- list(split = .circ_mix_move_split, merge = .circ_mix_move_merge,
              death = .circ_mix_move_death, birth = .circ_mix_move_birth)
  fns <- fns[intersect(control$moves, names(fns))]
  ## The candidate moves at a step are INDEPENDENT, so evaluate them concurrently
  ## (over min(#moves, cores) workers) -- this is what parallelises greedy, whose
  ## common "warm" path is otherwise a single-core sequence of move EM fits. Each
  ## move's own restart-fallback then runs serially (cx_inner), forking at one
  ## level only. Seeded per move (split's 2-means uses the RNG) so serial == parallel.
  mc <- min(length(fns), control$cores)
  cx_inner <- cx; if (mc > 1L) cx_inner$control$cores <- 1L
  trace <- list(); maxrounds <- 3L * control$kmax + 5L
  for (step in seq_len(maxrounds)) {
    cl <- .circ_mix_parallel_seeded(length(fns), function(i)
      tryCatch(fns[[i]](state, cx_inner, J, ...), error = function(e) NULL), mc)
    names(cl) <- names(fns)
    cands <- Filter(Negate(is.null), cl)
    if (!length(cands)) break
    Js  <- vapply(cands, `[[`, numeric(1), "J")
    win <- which.min(Js)
    if (Js[win] >= J - 1e-6 * (abs(J) + 1)) break          # no improving move
    cand <- cands[[win]]
    trace[[length(trace) + 1L]] <- data.frame(
      step = step, move = cand$move, K_from = state$K, K_to = cand$state$K,
      J_from = J, J_to = cand$J, dJ = cand$J - J, stringsAsFactors = FALSE)
    if (control$verbose)
      cat(sprintf("  [greedy] step %d: %-5s  K %d->%d  J %.2f->%.2f (dJ %.2f)\n",
                  step, cand$move, state$K, cand$state$K, J, cand$J, cand$J - J))
    state <- cand$state; J <- cand$J
  }
  ## carry the init-K restart health onto the final state (the basin signal)
  state$R <- base$state$R; state$restart_lls <- base$state$restart_lls
  state$basin_hits <- base$state$basin_hits
  list(state = state, trace = if (length(trace)) do.call(rbind, trace) else NULL)
}

## ---- the four structure moves (each: propose -> local EM -> {move, state, J}) -
## A move reshapes the per-unit responsibility matrix (contract 1), re-fits by a
## warm local EM, and reports its objective. The reshapes are pure (testable);
## the fitting is the ordinary engine. Each returns NULL when it cannot apply.

## Fit + score a proposed K, robustly but cheaply. The move's WARM responsibilities
## are tried first: if that single EM already lowers J below the current state, it
## is accepted immediately (the common, fast path -- an informative split/merge
## warm-starts almost perfectly). Only when the warm fit FAILS to improve does it pay
## for independent restarts at the new K (kmeans #1 + random) -- the load-bearing
## robustness net that (a) confirms a non-improving move really has no better
## optimum at that K before greedy gives up, and (b) recovers from a poor warm seed
## or component collapse. So growing is cheap, and the decision to stop (or to shrink
## from a poor fit) is still made against grid-quality optima. `J_cur` is the
## current objective; returns NULL only if every candidate fit failed.
.circ_mix_eval <- function(g_new, move, cx, J_cur, ...) {
  warm <- tryCatch(.circ_mix_em_core(cx$formula, cx$data, cx$family,
                     cx$control, cx$resp, g_new, grp = cx$grp, ...), error = function(e) NULL)
  tolJ <- 1e-6 * (abs(J_cur) + 1)
  if (!is.null(warm)) {
    Jw <- .circ_mix_objective(warm, cx$lam)
    if (Jw < J_cur - tolJ) return(list(move = move, state = warm, J = Jw))
  }
  ## warm did not improve (or failed): verify with restarts at the new K, run over
  ## control$cores (each seeded once in the parent -> serial/parallel identical).
  Kp <- ncol(g_new); R <- cx$control$restarts
  seeds <- sample.int(.Machine$integer.max, R)
  one <- function(r) { set.seed(seeds[r])
    im <- if (r == 1L) cx$control$init else "random"
    tryCatch(.circ_mix_em_once(cx$formula, cx$data, cx$family, Kp, cx$control,
                               cx$resp, init_method = im, grp = cx$grp, ...), error = function(e) NULL) }
  fits <- Filter(Negate(is.null), c(list(warm), .circ_mix_lapply(seq_len(R), one, cx$control$cores)))
  if (!length(fits)) return(NULL)
  best <- fits[[which.max(vapply(fits, `[[`, numeric(1), "loglik"))]]
  list(move = move, state = best, J = .circ_mix_objective(best, cx$lam))
}

## the split residual feature for a component: per response, the wrapped angular
## residual theta - mu_hat (circular) or the scaled residual (linear), column-
## bound over all factors. For a product component this is the JOINT residual
## feature, so the circular 2-means splits the worst blob along whichever response
## it is most over-dispersed in -- the per-product split rule, decided jointly
## rather than by picking a single driving factor. Returns the per-ROW feature plus
## its geometry flag; the caller reduces it to per-unit and clusters on it.
.circ_mix_resid_feature <- function(cp, cx) {
  resps <- .circ_mix_responses(cx$formula)
  preds <- if (inherits(cp, "circ_mix_product")) cmp_predict(cp, cx$data, "response")
           else list(cmp_predict(cp, cx$data, "response"))
  cols <- lapply(seq_along(resps), function(j) {
    mu <- preds[[j]][, 1L]                                 # location param of factor j
    v  <- as.numeric(cx$data[[resps[j]]])
    if (cx$resp_circ) atan2(sin(v - mu), cos(v - mu))      # wrapped angular residual
    else as.numeric(scale(v - mu))
  })
  list(x = do.call(cbind, cols), circular = cx$resp_circ)
}

## SPLIT (grow): split the worst-fit component (lowest mean per-obs density among
## its MAP members) by circular 2-means (circ_kmeans) on its angular residuals --
## the over-dispersed / bimodal angular-residual signal. Source: Ueda et al. (SMEM).
.circ_mix_move_split <- function(state, cx, J_cur, ...) {
  K <- state$K; if (K >= cx$control$kmax) return(NULL)
  grp <- cx$grp; z <- state$cluster                        # z is per UNIT
  Lr <- vapply(state$components, function(cp) cmp_logpdf(cp, cx$data), numeric(nrow(cx$data)))
  L  <- rowsum(Lr, grp)                                    # nu x K  per-unit log-density
  sz <- tabulate(z, K)
  md <- vapply(seq_len(K), function(k) if (sz[k] > 0L) mean(L[z == k, k]) else Inf, numeric(1))
  md[sz < 2L * cx$control$min_size] <- Inf                 # need room to split
  if (!any(is.finite(md))) return(NULL)
  j   <- which.min(md)
  rf  <- .circ_mix_resid_feature(state$components[[j]], cx)  # per-ROW residuals + flag
  ftx <- if (rf$circular) .circ_mix_unit_angle(rf$x, grp)   # per-UNIT angle(s)
         else .circ_mix_unit_feature(rf$x, grp)             # per-UNIT scaled mean
  mem <- which(z == j)                                      # the units in component j
  km  <- tryCatch(.circ_mix_kmeans(ftx[mem, , drop = FALSE], 2L, rf$circular),
                  error = function(e) NULL)
  if (is.null(km)) return(NULL)
  label <- .circ_mix_assign(ftx, km$centers, rf$circular)   # per-unit split label
  .circ_mix_eval(.circ_mix_gamma_split(state$gamma, j, label), "split", cx, J_cur, ...)
}

## MERGE (shrink): merge the two most similar components (smallest wrapped distance
## between their responsibility-weighted mean directions). Source: Ueda et al.
.circ_mix_move_merge <- function(state, cx, J_cur, ...) {
  K <- state$K; if (K <= cx$control$kmin || K < 2L) return(NULL)
  grp <- cx$grp
  ## per-component responsibility-weighted centre on each response; the pair to
  ## merge is the closest by distance summed over responses (joint on the torus).
  ## Responsibilities are per UNIT, so reduce each response to a per-unit value
  ## first (its circular/arithmetic mean within the unit) -- identity for rows.
  ctrs <- lapply(.circ_mix_responses(cx$formula), function(r) {
    vr <- as.numeric(cx$data[[r]])
    vu <- if (cx$resp_circ) atan2(rowsum(sin(vr), grp)[, 1L], rowsum(cos(vr), grp)[, 1L])
          else rowsum(vr, grp)[, 1L] / as.numeric(table(grp))
    vapply(seq_len(K), function(k) { w <- state$gamma[, k]
      if (cx$resp_circ) atan2(sum(w * sin(vu)), sum(w * cos(vu))) else sum(w * vu) / sum(w)
    }, numeric(1))
  })
  pair <- NULL; bd <- Inf
  for (a in seq_len(K - 1L)) for (b in (a + 1L):K) {
    d <- sum(vapply(ctrs, function(c0)
      if (cx$resp_circ) atan2(sin(c0[a] - c0[b]), cos(c0[a] - c0[b]))^2
      else (c0[a] - c0[b])^2, numeric(1)))
    if (d < bd) { bd <- d; pair <- c(a, b) }
  }
  .circ_mix_eval(.circ_mix_gamma_merge(state$gamma, pair[1L], pair[2L]), "merge", cx, J_cur, ...)
}

## DEATH (shrink): drop components whose soft size n_k = sum_i gamma_ik falls below
## min_size (folds in the collapse guards). Source: Figueiredo & Jain annihilation.
.circ_mix_move_death <- function(state, cx, J_cur, ...) {
  K <- state$K; nk <- colSums(state$gamma)
  keep <- which(nk >= cx$control$min_size)
  if (length(keep) >= K) return(NULL)                      # nothing undersized
  if (length(keep) < cx$control$kmin)                      # never fall below kmin
    keep <- sort(order(nk, decreasing = TRUE)[seq_len(cx$control$kmin)])
  if (length(keep) >= K) return(NULL)
  .circ_mix_eval(.circ_mix_gamma_death(state$gamma, sort(keep)), "death", cx, J_cur, ...)
}

## BIRTH (grow): seed a fresh component from the worst-explained points (the CRP
## "new table"). Source: the DP-means / Figueiredo-Jain birth move.
.circ_mix_move_birth <- function(state, cx, J_cur, ...) {
  K <- state$K; if (K >= cx$control$kmax) return(NULL)
  grp <- cx$grp; nu <- nrow(state$gamma)
  Lr <- vapply(state$components, function(cp) cmp_logpdf(cp, cx$data), numeric(nrow(cx$data)))
  L  <- rowsum(Lr, grp)                                    # nu x K  per-unit
  ld <- .circ_mix_loglik_rows(sweep(L, 2L, log(state$pi), "+"))
  m  <- min(2L * cx$control$min_size, nu %/% 2L)
  if (m < cx$control$min_size) return(NULL)
  idx <- order(ld)[seq_len(m)]                             # the worst-explained units
  .circ_mix_eval(.circ_mix_gamma_birth(state$gamma, idx), "birth", cx, J_cur, ...)
}

## ---- pure responsibility reshapes (no fitting; unit-tested) ------------------
## Each preserves row sums (every row still sums to 1), so the post-move gamma is
## a valid responsibility matrix the local EM can warm-start from.
.circ_mix_gamma_split <- function(g, j, label) {        # K -> K+1: split column j
  K <- ncol(g); a <- g[, j] * (label == 1L); b <- g[, j] * (label == 2L)
  left  <- if (j > 1L) g[, seq_len(j - 1L), drop = FALSE] else NULL
  right <- if (j < K)  g[, (j + 1L):K,      drop = FALSE] else NULL
  unname(cbind(left, a, b, right))
}
.circ_mix_gamma_merge <- function(g, j, l) {            # K -> K-1: merge cols j,l
  keep <- setdiff(seq_len(ncol(g)), c(j, l))
  unname(cbind(g[, keep, drop = FALSE], g[, j] + g[, l]))
}
.circ_mix_gamma_death <- function(g, keep) {            # K -> |keep|: drop & renorm
  sub <- g[, keep, drop = FALSE]; rs <- rowSums(sub)
  if (any(rs <= 0)) sub[rs <= 0, ] <- 1 / length(keep)
  unname(sub / rowSums(sub))
}
.circ_mix_gamma_birth <- function(g, idx) {             # K -> K+1: seed a new comp
  out <- cbind(g, 0); out[idx, ] <- 0; out[idx, ncol(out)] <- 1; unname(out)
}

## ============================================================================
## EM CORE + RESTARTS + INIT
## ============================================================================

## response geometry is a family fact -- circular unless the family declares
## response_circular = FALSE (the gausslss / gammalss l~c leg).
.circ_mix_resp_circular <- function(family)
  !is.null(family$param_names) && !isFALSE(family$response_circular)

## the pooled smoothing parameters for penalty = "fixed": fit the component spec
## ONCE on all the data with unit weights (a single circ_gam, or the product of
## them for a joint spec) and read its REML-selected sp. That value is shared
## across the K components and held fixed for the whole EM, so the penalty never
## moves -- the deterministic, monotone "select lambda once (pooled single-
## component fit)" choice. Returns a numeric vector (single) or a per-factor list
## (product), the shape .circ_mix_fit_component()'s `sp` consumes.
.circ_mix_pooled_sp <- function(formula, data, family, ...) {
  cp <- .circ_mix_fit_component(formula, data, family,
                               weights = rep(1, nrow(data)), ...)
  cmp_sp(cp)
}

## aggregate a per-ROW feature matrix to per-UNIT means (identity when each row is
## its own unit, i.e. grp = seq_len(n)). The unit average of a circular response's
## (cos, sin) is its mean resultant -- a sensible curve summary for subject init.
.circ_mix_unit_feature <- function(row_feat, grp)
  rowsum(row_feat, grp) / as.numeric(table(grp))

## aggregate per-ROW angle(s) to the per-UNIT circular mean (identity for rows: a
## row's unit angle is its own wrapped angle). A vector -> per-unit vector; a
## matrix, one column per response -> per-unit matrix, column-wise -- the angular
## counterpart of .circ_mix_unit_feature, used to seed circ_kmeans under group.
.circ_mix_unit_angle <- function(theta, grp) {
  if (is.matrix(theta))
    return(atan2(rowsum(sin(theta), grp), rowsum(cos(theta), grp)))
  atan2(rowsum(sin(theta), grp)[, 1L], rowsum(cos(theta), grp)[, 1L])
}

## ---- initial responsibilities: a hard label per UNIT -------------------------
## "kmeans" clusters the (per-unit) response so the components START separated --
## it breaks the symmetric saddle a weak-signal mixture would otherwise stall on.
## "random" assigns labels uniformly (the diversity seed for restarts beyond the
## first). The clustering feature is the per-unit value: the response ANGLE for a
## circular family (clustered by circ_kmeans on the circle / torus -- one
## coordinate per response, so a joint spec seeds on both angles jointly) or the
## scaled value for the linear l~c leg (ordinary stats::kmeans). Under group the
## per-unit value is the within-subject circular (or arithmetic) mean, so subjects
## -- not rows -- are clustered. Returns the feature matrix plus its geometry flag.
.circ_mix_init_features <- function(formula, data, family, grp) {
  circ <- .circ_mix_resp_circular(family)
  cols <- lapply(.circ_mix_responses(formula), function(r) {
    v <- as.numeric(data[[r]])
    if (circ) .circ_mix_unit_angle(v, grp)                          # per-unit angle
    else .circ_mix_unit_feature(matrix(as.numeric(scale(v)), ncol = 1L), grp)[, 1L]
  })
  list(x = do.call(cbind, cols), circular = circ)                  # one row per unit
}
.circ_mix_init_gamma <- function(method, feat, K, n) {
  lab <- NULL
  if (method == "kmeans" && nrow(unique(feat$x)) > K)
    lab <- tryCatch(.circ_mix_kmeans(feat$x, K, feat$circular)$cluster,
                    error = function(e) NULL)
  if (is.null(lab)) lab <- sample.int(K, n, replace = TRUE)
  g <- matrix(0, n, K)
  g[cbind(seq_len(n), lab)] <- 1
  g
}

## ---- parallel-or-serial apply over independent jobs --------------------------
## cores > 1 forks the jobs (mclapply); cores == 1 (the default) or Windows (no
## fork) runs serially.
.circ_mix_lapply <- function(X, FUN, cores) {
  if (is.null(cores) || cores <= 1L || .Platform$OS.type == "windows")
    return(lapply(X, FUN))
  parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE)
}

## Run R independent jobs, each seeded from a stream drawn ONCE here, then RESTORE
## the parent RNG to its post-draw state. The seeds make the jobs independent (so
## the result is the same for any core count); the restore is what keeps serial and
## parallel IDENTICAL: a serial run's in-process set.seed() calls would otherwise
## advance the global stream and change what a *later* call (e.g. the next grid-K's
## restarts) draws, whereas the parallel run seeds inside isolated forks that never
## touch the parent. Restoring undoes that pollution, so `cores` is a pure speed knob.
.circ_mix_parallel_seeded <- function(R, FUN, cores) {
  seeds <- sample.int(.Machine$integer.max, R)
  state <- get(".Random.seed", envir = .GlobalEnv)
  out   <- .circ_mix_lapply(seq_len(R), function(r) { set.seed(seeds[r]); FUN(r) }, cores)
  assign(".Random.seed", state, envir = .GlobalEnv)
  out
}

## ---- random-restart wrapper: keep the largest-loglik run (P1) ----------------
## Restart 1 uses control$init (kmeans by default -- a strong, separated start);
## later restarts are random, for basin diversity. The R runs are the
## embarrassingly-parallel work -- independent and identical whether run serially
## or across `control$cores` forked workers. The best (largest-loglik) run is kept.
.circ_mix_restarts <- function(formula, data, family, K, control, resp, R, ...) {
  fits <- .circ_mix_parallel_seeded(R, function(r) {
    init_r <- if (r == 1L) control$init else "random"
    tryCatch(
      .circ_mix_em_once(formula, data, family, K, control, resp,
                        init_method = init_r, ...),
      error = function(e) {
        if (control$verbose)
          cat(sprintf("   restart %d FAILED: %s\n", r, conditionMessage(e)))
        NULL })
  }, control$cores)
  lls  <- vapply(fits, function(f) if (is.null(f)) NA_real_ else f$loglik, numeric(1))
  ok   <- which(!vapply(fits, is.null, logical(1)))
  if (!length(ok)) stop("all ", R, " restarts failed for K = ", K, ".")
  best <- fits[[ ok[which.max(lls[ok])] ]]
  best$restart_lls <- lls
  best$basin_hits  <- sum(abs(lls - best$loglik) < 1e-3, na.rm = TRUE)
  best
}

## ---- one EM run from a chosen init method ------------------------------------
.circ_mix_em_once <- function(formula, data, family, K, control, resp,
                              init_method = "random", grp = NULL, ...) {
  if (is.null(grp)) grp <- seq_len(nrow(data))
  feat <- .circ_mix_init_features(formula, data, family, grp)   # list(x, circular); one row per unit
  g    <- .circ_mix_init_gamma(init_method, feat, K, nrow(feat$x))
  .circ_mix_em_core(formula, data, family, control, resp, g, grp = grp, ...)
}

## ---- one EM run from a GIVEN responsibility matrix ---------------------------
## The shared inner loop -- used by a fresh init (em_once) AND by every structure
## move (which warm-starts it from the post-move gamma). Ordering is M-then-E so
## the returned state is mutually consistent: `loglik` is the observed-data
## log-likelihood of the returned (components, pi), and `gamma` are that state's
## E-step responsibilities. Each full E->M cycle is the standard EM ascent step,
## so the recorded log-likelihoods are monotone for parametric (penalty-free)
## components -- recorded here (the driver warns ONCE for the returned run).
.circ_mix_em_core <- function(formula, data, family, control, resp, g, grp = NULL, ...) {
  n <- nrow(data); K <- ncol(g)
  if (is.null(grp)) grp <- seq_len(n)               # rows: each its own unit
  nu <- nrow(g)                                     # #units (rows, or subjects)
  penalty  <- if (is.null(control$penalty)) "auto" else control$penalty
  hard     <- identical(control$assign, "hard")     # hard CEM vs soft EM
  sp_every <- if (is.null(control$sp_every)) 5L else control$sp_every
  ## the outer optimiser for the REML smoothing-parameter SEARCH (penalty
  ## "auto"/"scheduled"); NULL for a parametric spec (no search) so those fits stay
  ## bit-identical to mgcv's default. "efs" (extended Fellner-Schall) is markedly
  ## faster than outer Newton for the available.derivs = 2 families and selects the
  ## same smoothness (see the `optimizer` control arg).
  opt <- if (.circ_mix_formula_smooth(formula)) control$optimizer else NULL
  if (K == 1L) {
    cp <- .circ_mix_fit_component(formula, data, family, weights = rep(1, n),
                                  sp = control$sp,
                                  optimizer = if (is.null(control$sp)) opt else NULL, ...)
    ll <- sum(cmp_logpdf(cp, data))                 # = sum over units of unit density
    return(list(components = list(cp), pi = 1, gamma = matrix(1, nu, 1),
                cluster = rep(1L, nu), loglik = ll, ll_path = ll,
                iter = 1L, converged = TRUE, monotone = TRUE, worst_drop = 0, K = 1L))
  }
  components <- vector("list", K)
  start_k <- vector("list", K)                      # coefficient warm starts (NULL => pilot)
  warm_k  <- vector("list", K)                      # last selected sp per component:
                                                    # warm-starts the next REML search
                                                    # (mgcv in.out), unchanged optimum
  ## per-component smoothing parameters held between re-selections:
  ##   "auto"      -> NULL every iter (REML re-selects each M-step; non-monotone)
  ##   "fixed"     -> control$sp (a pooled pilot's sp), never moves -> monotone
  ##   "scheduled" -> START at the pooled pilot's sp (like "fixed"), then re-select
  ##                  every sp_every iters AFTER the first -- so the early, diffuse
  ##                  responsibilities use the stable pooled sp (no over-flexing),
  ##                  and sp only re-adapts once the partition has settled.
  sp_k <- rep(list(if (penalty %in% c("fixed", "scheduled")) control$sp else NULL), K)
  pi <- colMeans(g)
  prev <- -Inf; ll_path <- numeric(0); conv <- FALSE
  parametric <- NA; worst_drop <- 0; ll <- NA_real_; prev_z <- NULL
  for (iter in seq_len(control$max_iter)) {
    reselect <- switch(penalty, auto = TRUE, fixed = FALSE,
                       scheduled = (iter > 1L && iter %% sp_every == 0L))
    ## M-step: weighted component refits + mixing weights, each warm-started from
    ## the previous iteration's coefficients (tracks the moving weighted MLE, so
    ## the parametric / fixed-sp EM stays monotone and converges in fewer
    ## iterations). A unit's responsibility is broadcast to all its rows (g[grp, k])
    ## -- the subject-level weight under group, identity for rows. The penalty is
    ## either re-selected by REML (sp = NULL) or held at sp_k[[k]] this iteration.
    pi <- colMeans(g)
    for (k in seq_len(K)) {
      components[[k]] <- .circ_mix_fit_component(
        formula, data, family, weights = pmax(g[grp, k], control$wfloor),
        start = start_k[[k]], sp = if (reselect) NULL else sp_k[[k]],
        optimizer = if (reselect) opt else NULL,
        warm = if (reselect) warm_k[[k]] else NULL, ...)
      if (reselect) {                               # remember the selected sp:
        spk <- cmp_sp(components[[k]])              #   warm_k warm-starts the next
        warm_k[[k]] <- spk                          #   search; sp_k holds it between
        if (!identical(penalty, "auto")) sp_k[[k]] <- spk   # scheduled re-selections
      }
    }
    start_k <- lapply(components, cmp_coef)
    if (is.na(parametric))                          # fixed across iterations
      parametric <- !any(vapply(components, .circ_mix_has_smooth, logical(1)))
    ## E-step: per-row densities summed within unit (rowsum) -> per-unit log-lik
    ## and responsibilities. Soft EM uses the log-sum-exp mixture log-lik and
    ## softmax gamma; hard CEM seats each unit (curve) wholly at its argmax and
    ## uses the classification log-lik (sum of per-unit maxima).
    Lr <- vapply(components, function(cp) cmp_logpdf(cp, data), numeric(n))
    L  <- rowsum(Lr, grp)                           # nu x K
    logmix <- sweep(L, 2L, log(pi), "+")
    if (hard) {
      cl <- .circ_mix_classify(logmix); ll <- cl$loglik; g <- cl$gamma; z <- cl$z
    } else {
      ll <- sum(.circ_mix_loglik_rows(logmix)); g <- .circ_mix_responsibilities(logmix)
    }
    ll_path <- c(ll_path, ll)
    if (control$verbose)
      cat(sprintf("    iter %3d  ll = %.5f  pi = (%s)\n", iter, ll,
                  paste(sprintf("%.3f", pi), collapse = ", ")))
    if (iter > 1L) {
      worst_drop <- max(worst_drop, prev - ll)      # how far ll ever fell
      if (abs(ll - prev) < control$tol * (abs(prev) + control$tol)) {
        conv <- TRUE; break
      }
      if (hard && identical(z, prev_z)) { conv <- TRUE; break }  # CEM: partition stable
    }
    prev <- ll; if (hard) prev_z <- z
  }
  ## monotone is EXPECTED only when the penalty does not move: parametric (no
  ## penalty) or fixed-sp. Under auto/scheduled smooths the REML penalty moves
  ## between iterations, so a dip is not a bug and is not flagged.
  monotone_expected <- isTRUE(parametric) || identical(penalty, "fixed")
  monotone <- !monotone_expected || worst_drop <= 1e-6 * (abs(ll) + 1)
  list(components = components, pi = pi, gamma = g, cluster = max.col(g),
       loglik = ll, ll_path = ll_path, iter = iter, converged = conv,
       monotone = monotone, worst_drop = worst_drop, K = K)
}

## ---- mixture density / responsibilities at arbitrary data --------------------
## logmix[i,k] = log pi_k + log f_k(y_i | x_i); shared by the E-step and predict.
.circ_mix_logmix <- function(object, newdata) {
  K <- object$K
  pmat <- .circ_mix_gate(object$gating, newdata, K)
  L <- vapply(object$components, function(cp) cmp_logpdf(cp, newdata),
              numeric(nrow(newdata)))
  if (!is.matrix(L)) L <- matrix(L, ncol = K)
  L + log(pmat)
}

.circ_mix_loglik_rows <- function(logmix) {
  m <- apply(logmix, 1L, max)
  m + log(rowSums(exp(logmix - m)))
}

.circ_mix_responsibilities <- function(logmix) {
  m <- apply(logmix, 1L, max)
  g <- exp(logmix - m)
  g / rowSums(g)
}

## hard (CEM / DP-means) assignment of a per-unit logmix matrix: each unit seats
## wholly at its argmax component. Returns the labels z, the 0/1 responsibility
## matrix, and the classification log-likelihood sum_u max_k (log pi_k + log f_k)
## -- the criterion hard EM ascends (the soft analogue of the two helpers above).
.circ_mix_classify <- function(logmix) {
  z   <- max.col(logmix, ties.method = "first")
  sel <- cbind(seq_along(z), z)
  g   <- matrix(0, nrow(logmix), ncol(logmix)); g[sel] <- 1
  list(z = z, gamma = g, loglik = sum(logmix[sel]))
}

## ---- gating (contract 3): constant mixing weights ---------------------------
.circ_mix_gate <- function(gating, newdata, K) {
  if (identical(gating$type, "constant"))
    return(matrix(gating$pi, nrow(newdata), K, byrow = TRUE))
  stop("unknown gating type '", gating$type, "'.")
}

## ---- formula helpers ---------------------------------------------------------
## the response owning one FACTOR -- a formula, or (the nested case) a list of LSS
## formulas whose first carries the LHS. NA when the factor has no LHS (an LHS-less
## LSS predictor of the preceding factor).
.circ_mix_factor_response <- function(f) {
  f1 <- if (inherits(f, "formula")) f else f[[1L]]
  if (length(f1) == 3L) all.vars(f1[[2L]])[1L] else NA_character_
}

## response name = the LHS of the first formula.
.circ_mix_response <- function(formula) {
  f1 <- if (inherits(formula, "formula")) formula else formula[[1L]]
  all.vars(f1[[2L]])[1L]
}

## all DISTINCT responses across a formula spec, in chain-rule order. One element
## => a single-response cell; two or more => a joint product component.
.circ_mix_responses <- function(formula) {
  fl <- if (inherits(formula, "formula")) list(formula) else formula
  unique(stats::na.omit(vapply(fl, .circ_mix_factor_response, character(1))))
}

## number of distinct responses: 1 = single-response cell (circ_gam's own LSS
## grammar applies); >= 2 = a joint product component. Formulas/factors after the
## first with no LHS share the preceding response.
.circ_mix_n_responses <- function(formula) length(.circ_mix_responses(formula))

## does the spec contain a penalised smooth term -- s() / te() / ti() / t2()? A
## deparse heuristic across all factors (a false positive only forces a safe serial
## run); the boundary class excludes names like cos( / my.s(. Used only to keep the
## fork-based parallelism off smooth fits (see the cores guard).
.circ_mix_formula_smooth <- function(formula) {
  fl  <- if (inherits(formula, "formula")) list(formula) else formula
  txt <- paste(unlist(lapply(fl, function(f)
    deparse(if (inherits(f, "formula")) f else if (is.list(f)) f[[1L]] else f))),
    collapse = " ")
  grepl("(^|[^[:alnum:]._])(s|te|ti|t2)[[:space:]]*\\(", txt)
}

## ---- the clustering UNIT (contract 1) ----------------------------------------
## group = NULL ⇒ each ROW is its own unit (grp = 1..n); group = ~id ⇒ the SUBJECT
## index (grp in 1..m), so a whole trajectory seats at one cluster. The integer grp
## maps each row to its unit; rowsum(., grp) collapses rows to units and g[grp, ]
## broadcasts a unit's weight back to its rows. Returns the index plus the unit
## kind/labels for the fitted object.
.circ_mix_group_index <- function(group, data, n) {
  if (is.null(group))
    return(list(grp = seq_len(n), n_units = n, kind = "row", labels = NULL))
  gv <- all.vars(group)
  if (length(gv) != 1L)
    stop("group must name a single variable, e.g. group = ~ id.")
  if (is.null(data[[gv]]))
    stop("the group variable '", gv, "' is not a column of 'data'.")
  f <- factor(data[[gv]])
  list(grp = as.integer(f), n_units = nlevels(f), kind = "subject",
       labels = levels(f))
}


## ---- component interface (contract 2), E-step density, and S3 methods --------

## The circ_mix COMPONENT interface (contract 2 of the engine).
##
## The EM loop never touches a circ_gam object directly; it goes through the
## four cmp_*() accessors below -- logpdf / edf / coef / predict. A single-response
## component wraps ONE weighted circ_gam fit. The joint
## component is a second shape -- a PRODUCT of d weighted circ_gam fits, one per
## chain-rule factor f(y) = prod_j f(y_j | parents_j) -- whose joint log-density is
## the SUM of the factor log-densities. The four accessors carry an inherits()
## branch for it (they are plain functions, not S3 generics: an unregistered
## internal generic is not found when dispatched through vapply), and the inner EM
## loop is unchanged -- it never learns which shape it holds. Keeping every
## per-component access behind these four is exactly what makes the joint case
## additive.

.circ_mix_component <- function(fit)
  structure(list(fit = fit), class = "circ_mix_component")

## a product (joint) component: d weighted circ_gam fits in chain-rule order,
## named by their response so coef()/print read cleanly.
.circ_mix_product_component <- function(fits)
  structure(list(fits = fits), class = c("circ_mix_product", "circ_mix_component"))

## Split a formula spec into chain-rule FACTORS, each an ordinary circ_gam spec:
##   - a single formula              -> one factor (single-response cell)
##   - a list with <= 1 response     -> one factor (circ_gam's own LSS grammar:
##                                      LSS predictors of one response)
##   - a list with >= 2 responses    -> one factor per distinct response; an
##                                      element with an LHS starts a factor, any
##                                      LHS-less element after it (or a nested list)
##                                      supplies that factor's further LSS predictors
.circ_mix_factor_specs <- function(formula) {
  if (inherits(formula, "formula")) return(list(formula))
  if (.circ_mix_n_responses(formula) <= 1L) return(list(formula))
  factors <- list(); cur <- NULL
  for (f in formula) {
    if (is.list(f) && !inherits(f, "formula")) {       # a nested per-factor LSS list
      if (!is.null(cur)) factors[[length(factors) + 1L]] <- cur
      cur <- NULL; factors[[length(factors) + 1L]] <- f
    } else if (length(f) == 3L) {                       # has an LHS -> a new factor
      if (!is.null(cur)) factors[[length(factors) + 1L]] <- cur
      cur <- list(f)
    } else {                                            # LHS-less -> an LSS predictor
      cur[[length(cur) + 1L]] <- f
    }
  }
  if (!is.null(cur)) factors[[length(factors) + 1L]] <- cur
  ## unwrap a single-formula factor to a bare formula (circ_gam's plain entry)
  lapply(factors, function(fl)
    if (is.list(fl) && length(fl) == 1L && inherits(fl[[1L]], "formula")) fl[[1L]] else fl)
}

## Fit one component (the M-step). A single-response spec is one weighted
## circ_gam; a multi-response (joint) spec is d weighted circ_gam fits sharing the
## SAME responsibilities `weights`, returned as a product component. `start` warm-
## starts each factor: a coefficient vector for the single case, a per-factor list
## (as cmp_coef returns for a product) for the joint case. `sp` FIXES the smoothing
## parameters (penalty = "fixed"/"scheduled"): a numeric vector for the
## single case, a per-factor list for the joint case; NULL (the default) lets REML
## re-select them, as for parametric components, which have none. An empty-length
## sp (a factor with no smooths) is passed as NULL -- nothing to hold.
## `optimizer` (an mgcv outer optimiser, e.g. "efs") and `warm` (a previous
## M-step's selected smoothing parameters) speed up the REML SEARCH without changing
## what it selects: they are injected only for a factor that is actually searching
## (its fixed `sp` is NULL) and carries smooths (a non-empty warm value), and only
## when the caller did not already pass `optimizer`/`in.out` through `...`. `warm`
## becomes mgcv's `in.out` (a search start, NOT a fixed value), so the selected
## smoothness is identical -- just reached in a step or two. Shapes mirror `sp`: a
## numeric vector for a single-response component, a per-factor list for a product.
.circ_mix_fit_component <- function(formula, data, family, weights, start = NULL,
                                    sp = NULL, optimizer = NULL, warm = NULL, ...) {
  norm_sp <- function(s) if (is.null(s) || !length(s)) NULL else s
  dots <- list(...)
  ## per-component centering off the tan-half wall: default on, weight-aware
  ## (each component centers on its own responsibility-weighted mode), overridable
  ## via circ_mix(..., center = FALSE). The E-step (.circ_logpdf) frame-aligns, so
  ## responsibilities/loglik are unchanged; only the M-step basin improves.
  ctr_flag <- if (is.null(dots$center)) TRUE else dots$center
  dots$center <- NULL
  ## extra gam args that warm-start / re-route the REML search for one factor whose
  ## fixed sp is `fixed_j` (NULL => searching) and warm start is `warm_j`.
  search_args <- function(fixed_j, warm_j) {
    if (!is.null(fixed_j)) return(list())                 # sp fixed: no search
    a <- list()
    if (!is.null(optimizer) && is.null(dots$optimizer)) a$optimizer <- optimizer
    if (!is.null(warm_j) && length(warm_j) && is.null(dots$in.out))
      a$in.out <- list(sp = warm_j, scale = 1)            # mgcv: a search START
    a
  }
  specs <- .circ_mix_factor_specs(formula)
  if (length(specs) == 1L) {
    fx <- norm_sp(sp)
    return(.circ_mix_component(do.call(circ_gam, c(
      list(specs[[1L]], data = data, family = family, weights = weights,
           start = start, sp = fx, center = ctr_flag), search_args(fx, warm), dots))))
  }
  fits <- lapply(seq_along(specs), function(j) {
    fx <- norm_sp(if (is.list(sp) && length(sp) >= j) sp[[j]] else NULL)
    wj <- if (is.list(warm) && length(warm) >= j) warm[[j]] else NULL
    do.call(circ_gam, c(
      list(specs[[j]], data = data, family = family, weights = weights,
           start = if (is.list(start) && length(start) >= j) start[[j]] else NULL,
           sp = fx, center = ctr_flag), search_args(fx, wj), dots))
  })
  names(fits) <- vapply(specs, .circ_mix_factor_response, character(1))
  .circ_mix_product_component(fits)
}

## --- the four-accessor interface ----------------------------------------------
## per-observation log-density log f_k(y_i | x_i) on `newdata` (the E-step
## density). For a product component this is the JOINT log-density: the sum of the
## factor log-densities (chain rule), each the same .circ_logpdf() primitive.
cmp_logpdf <- function(cp, newdata) {
  if (inherits(cp, "circ_mix_product"))
    return(Reduce(`+`, lapply(cp$fits, .circ_logpdf, newdata)))
  .circ_logpdf(cp$fit, newdata)
}

## effective degrees of freedom (sum of edf == #coef for a parametric fit); for a
## product, summed over factors -- df = (K-1) + sum_k sum_j edf_kj falls straight out.
cmp_edf <- function(cp) {
  if (inherits(cp, "circ_mix_product"))
    return(sum(vapply(cp$fits, function(f) sum(f$edf), numeric(1))))
  sum(cp$fit$edf)
}

## the coefficient vector; for a product, the per-factor list cmp_coef returns is
## also the warm-start structure .circ_mix_fit_component consumes (coef.circ_mix
## flattens it for display).
cmp_coef <- function(cp) {
  if (inherits(cp, "circ_mix_product"))
    return(lapply(cp$fits, stats::coef))
  stats::coef(cp$fit)
}

## response/link prediction on `newdata`; for a product, a per-factor list.
cmp_predict <- function(cp, newdata, type = "response") {
  if (inherits(cp, "circ_mix_product"))
    return(lapply(cp$fits, function(f)
      stats::predict(f, newdata = newdata, type = type)))
  stats::predict(cp$fit, newdata = newdata, type = type)
}

## the fitted smoothing parameters -- a numeric vector for a single component, a
## per-factor list for a product. The shape .circ_mix_fit_component()'s `sp`
## consumes: penalty = "fixed" reads it from a pooled pilot fit and holds it,
## penalty = "scheduled" captures it on each re-selection to hold between. A
## parametric factor has none (length-0), which fit_component passes as NULL.
cmp_sp <- function(cp) {
  if (inherits(cp, "circ_mix_product"))
    return(lapply(cp$fits, function(f) f$sp))
  cp$fit$sp
}

## does a component carry any penalised smooth? (a product is "smooth" if ANY
## factor is). Drives em_core's monotonicity expectation -- parametric / fixed-sp
## EM is monotone; a live REML penalty (smooth terms) need not be.
.circ_mix_has_smooth <- function(cp) {
  fits <- if (inherits(cp, "circ_mix_product")) cp$fits else list(cp$fit)
  any(vapply(fits, function(f) length(f$smooth) > 0L, logical(1)))
}

## the training frame of a (single or product) component -- the default newdata
## for predict()/plot(). A single component carries it as fit$model; a product
## component's raw variables are spread across its factor model frames (the torus's
## phi lives in the marginal factor, psi in the conditional), so reassemble the
## union of factor variables, taking each from the first factor that carries it raw.
.circ_mix_model_frame <- function(cp) {
  if (!inherits(cp, "circ_mix_product")) return(cp$fit$model)
  vars <- unique(unlist(lapply(cp$fits, function(f)
    all.vars(if (is.list(f$formula)) f$formula[[1L]] else f$formula))))
  cols <- list()
  for (v in vars) for (f in cp$fits)
    if (!is.null(f$model[[v]])) { cols[[v]] <- f$model[[v]]; break }
  as.data.frame(cols)
}

## .circ_logpdf(): the per-fit E-step density -- the primitive cmp_logpdf() calls.
##
## The unweighted per-observation log-density log f(y_i | x_i) of a fitted
## circ_gam, recovered by feeding the prediction lpmatrix (which carries the
## `lpi` attribute the family's ll() needs) back through that ll() at deriv = 0.
## A family's $l0 is the unweighted per-observation log-density by construction
## -- the weighted-likelihood contract scales the objective and derivative
## blocks by the prior weights but leaves $l0 untouched -- so this is exactly the
## quantity a finite-mixture E-step needs, and (summed over factors) the
## primitive a joint product component is built from.
##
## Response-agnostic: it reads only $l0, so it works for any circlss family --
## the circular ones (vmlss, pnlss, ...) and the linear-response forks
## (gausslss, gammalss) that carry the l~c leg -- without knowing which.

.circ_logpdf <- function(fit, newdata) {
  fam <- fit$family
  if (is.null(fam$ll))
    stop(".circ_logpdf() needs a circlss/general family carrying an ll() method; ",
         "family '", fam$family, "' has none.")
  f1   <- if (is.list(fit$formula)) fit$formula[[1L]] else fit$formula
  resp <- all.vars(f1[[2L]])[1L]               # the LHS (response) name
  if (missing(newdata) || is.null(newdata))
    newdata <- fit$model
  if (is.null(newdata[[resp]]))
    stop("newdata must carry the response '", resp,
         "' for the density to be evaluated.")
  Xlp <- stats::predict(fit, newdata = newdata, type = "lpmatrix")  # carries lpi
  y   <- as.numeric(newdata[[resp]])
  ## coef/Xlp live in the component's (possibly centered) fit frame; align the
  ## original-frame response to it so the density is evaluated consistently. The
  ## value is frame-independent (cos(y - mu) is unchanged), so responsibilities,
  ## loglik and BIC are identical to an uncentered fit -- centering only changes
  ## which basin the M-step reaches.
  ctr <- fit$circ_center
  if (!is.null(ctr) && is.finite(ctr) && ctr != 0) y <- wrap(y - ctr)
  drop(fam$ll(y, Xlp, stats::coef(fit), rep(1, nrow(newdata)),
              fam, deriv = 0)$l0)
}

## S3 methods for the "circ_mix" fitted object. These read fields the object
## always carries (gamma is always per-unit, components always expose the
## contract-2 generics, gating always answers .circ_mix_gate). AIC()/BIC() are
## inherited from stats via logLik.circ_mix.

## human-readable leg label, mirroring print.circ_gam's geometry switch
.circ_mix_leg <- function(kind) switch(as.character(kind),
  cl = "circular-linear (c~l)", cc = "circular-circular (c~c)",
  lc = "linear-circular (l~c)", ll = "location-scale (l~l)",
  joint = "joint torus density", "circular")

## deparse a formula spec (a formula, or a list of formulas / nested LSS lists)
## to a single compact string for printing.
.circ_mix_deparse <- function(formula) {
  one <- function(f) paste(deparse(if (inherits(f, "formula")) f else
    if (is.list(f)) f[[1L]] else f), collapse = " ")
  if (inherits(formula, "formula")) one(formula)
  else paste(vapply(formula, one, ""), collapse = " | ")
}

#' @rdname circ_mix
#' @export
print.circ_mix <- function(x, ...) {
  hard <- identical(x$control$assign, "hard")
  cat(sprintf("Finite mixture of circular GAMs (circ_mix) -- %s\n",
              .circ_mix_leg(x$geometry)))
  units <- if (identical(x$unit$kind, "subject"))
    sprintf("%d subjects (%d rows)", x$unit$n_units, length(x$unit$index))
    else sprintf("%d obs", x$unit$n_units)
  cat(sprintf("  family %s | K = %d component%s | %s%s\n",
              x$family$family, x$K, if (x$K > 1L) "s" else "", units,
              if (hard) " | hard (CEM)" else ""))
  cat(sprintf("  formula: %s\n", .circ_mix_deparse(x$formula)))
  cat(sprintf("  logLik%s = %.2f | df = %.2f | BIC = %.2f\n",
              if (hard) " (classification)" else "", x$loglik, x$df, x$bic))
  pi <- x$gating$pi
  cat("  components (MAP):\n")
  for (k in seq_len(x$K))
    cat(sprintf("    %2d:  pi = %.3f   n = %3d   edf = %.2f\n",
                k, pi[k], x$nk[k], x$edf[k]))
  if (x$Gtilde < x$K)
    cat(sprintf("  %d of %d components are empty under MAP.\n",
                x$K - x$Gtilde, x$K))
  ## auto-K provenance: how the component count was decided
  if (identical(x$search, "greedy")) {
    nm <- if (is.null(x$move_trace)) 0L else nrow(x$move_trace)
    cat(sprintf("  auto-K (greedy): K_init %d -> %d via %d accepted move%s",
                x$K_init, x$K, nm, if (nm == 1L) "" else "s"))
    if (nm) cat(sprintf(" [%s]", paste(x$move_trace$move, collapse = ", ")))
    cat(".\n")
  } else if (identical(x$search, "grid") && !is.null(x$move_trace)) {
    cat(sprintf("  auto-K (grid): swept K = %d..%d, selected K = %d by min %s.\n",
                min(x$move_trace$K), max(x$move_trace$K), x$K,
                if (is.null(x$control$lambda)) "BIC" else "J"))
  }
  cat(sprintf("  %s in %d iterations",
              if (x$converged) "converged" else "STOPPED (max_iter)", x$iter))
  if (x$K > 1L)
    cat(sprintf("; restart basin hits %d/%d", x$restarts$basin_hits, x$restarts$R))
  cat(".\n")
  if (x$K > 1L && !isTRUE(x$monotone))
    cat("  note: a non-monotone EM step was seen -- inspect $ll_path.\n")
  invisible(x)
}

#' @rdname circ_mix
#' @export
summary.circ_mix <- function(object, ...) {
  print(object)
  cat("\n  per-component coefficients:\n")
  cf <- stats::coef(object)
  for (k in seq_along(cf)) {
    cat(sprintf("  [component %d]\n", k))
    print(round(cf[[k]], 4))
  }
  invisible(object)
}

#' @rdname circ_mix
#' @export
logLik.circ_mix <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = object$unit$n_units,
            class = "logLik")

#' @rdname circ_mix
#' @export
coef.circ_mix <- function(object, ...) {
  cf <- lapply(object$components, function(cp) {
    co <- cmp_coef(cp)                          # product => named per-factor list
    if (is.list(co)) unlist(co) else co         # flatten for display (resp.term names)
  })
  names(cf) <- paste0("component", seq_along(cf))
  cf
}

#' @rdname circ_mix
#' @export
predict.circ_mix <- function(object, newdata,
                             type = c("cluster", "density", "response"),
                             log = FALSE, ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata))
    newdata <- .circ_mix_model_frame(object$components[[1L]])
  if (type == "response") {
    out <- lapply(object$components,
                  function(cp) cmp_predict(cp, newdata, type = "response"))
    names(out) <- paste0("component", seq_along(out))
    return(out)
  }
  if (is.null(newdata[[object$response]]))
    stop("predict(type = \"", type, "\") needs the response '",
         object$response, "' in 'newdata'.")
  logmix <- .circ_mix_logmix(object, newdata)
  if (type == "cluster")
    return(.circ_mix_responsibilities(logmix))
  ld <- .circ_mix_loglik_rows(logmix)            # log mixture density per row
  if (log) ld else exp(ld)
}
