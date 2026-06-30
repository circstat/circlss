#' Classical circular regression
#'
#' The closed-form, unpenalized counterpart of \code{\link{circ_gam}}: the
#' textbook circular regressions that the literature reports. \code{circ_lm} is
#' \strong{parametric only} -- a smooth term (\code{s()}, \code{te()}, \dots) is an
#' error pointing to \code{\link{circ_gam}} -- and carries no \code{family}
#' argument: \code{"cl"} is von Mises by construction, while \code{"cc"} and
#' \code{"lc"} are ordinary least squares with a residual concentration reported
#' as a summary. Reach for \code{circ_lm} for the textbook fit or a fast
#' unpenalized baseline; reach for \code{\link{circ_gam}} for penalized smooths,
#' per-parameter modelling, or any family beyond von Mises.
#'
#' @param formula A formula, or (for \code{type = "cl"}) a list of one or two
#'   formulas. The first formula is two-sided and names the response.
#' @param data A data frame holding the response and covariates.
#' @param type For \code{circ_lm}, which classical fit. \code{"cl"} circular
#'   response on linear covariate(s) (Fisher--Lee von Mises); \code{"cc"} circular
#'   on circular (harmonic); \code{"lc"} linear on circular (harmonic). Hyphenated
#'   spellings (\code{"c-l"}, \code{"c-c"}, \code{"l-c"}) are accepted. For
#'   \code{predict} on a \code{"cl"} fit, the quantity returned: \code{"direction"}
#'   (the mean direction, the default) or \code{"kappa"} (the fitted concentration).
#' @param order Order of the trigonometric polynomial for \code{"cc"} and
#'   \code{"lc"} (number of harmonics of the angular predictor). Ignored for
#'   \code{"cl"}.
#' @param init Starting values for the \code{"cl"} iteration. \code{NULL}
#'   (default) starts cold: all coefficients zero, so \eqn{\kappa \equiv 1}. A
#'   named list \code{list(beta=, alpha=, gamma=)} sets explicit starts -- any
#'   component omitted falls back to its cold value; \code{beta} and \code{gamma}
#'   take one value per covariate, \code{alpha} a single number. This lets you seed
#'   the joint (mixed) fit with estimates from separately fitted mean-only and
#'   kappa-only models, as Fisher (1993, Sec. 6.4.4) suggests. A bare numeric
#'   vector is taken, as before, as the mean-direction coefficients (\code{beta}).
#' @param tol,maxit,verbose IRLS convergence tolerance, iteration cap, and
#'   per-iteration logging (\code{"cl"} only).
#' @param se How standard errors are computed (\code{"cl"} only).
#'   \code{"asymptotic"} (default) uses the expected-information formulae;
#'   \code{"bootstrap"} replaces them with a parametric bootstrap -- simulate from
#'   the fitted vM(\eqn{\mu_i, \kappa_i}), refit, and take the spread of the
#'   estimates -- which Fisher (1993, Sec. 8.4) recommends below \eqn{n \approx 25}
#'   -- \eqn{30}, where the asymptotic SEs are unreliable. Stochastic, so set a
#'   seed for reproducibility.
#' @param R Number of bootstrap resamples when \code{se = "bootstrap"}.
#' @param object,x A fitted \code{circ_lm} model.
#' @param newdata A data frame of new predictor values. For \code{predict},
#'   omitting it returns the fitted values.
#' @param digits Number of significant digits for \code{print}.
#' @param ... Unused.
#' @return
#' An object of class \code{circ_lm}: a list whose contents depend on \code{type}.
#' \code{"cl"} carries \code{mu}, \code{kappa}, the coefficient vectors
#' (\code{beta}/\code{alpha}/\code{gamma}) with standard errors, and
#' \code{loglik}/\code{aic}/\code{bic}; \code{"cc"} carries the cos/sin
#' \code{coefficients}, \code{rho}, residual \code{kappa}, and higher-order-test
#' \code{p_values}; \code{"lc"} carries the \code{coefficients}, per-harmonic
#' \code{amplitude}/\code{phase} (with delta-method SEs), and the usual
#' least-squares fit metrics. \code{predict}, \code{coef}, \code{fitted},
#' \code{residuals}, and (for \code{"cl"}/\code{"lc"}) \code{logLik} methods are
#' provided.
#' @details
#' \strong{cl -- circular response, linear covariates.} The Fisher and Lee (1992)
#' von Mises model fitted by Green's (1984) IRLS, with the concentration
#' extensions of Fisher (1993) Sec. 6.4. The mean direction is
#' \eqn{\mu_i = \mu_0 + 2\,\mathrm{atan}(x_i^\top\beta)} (the offset is
#' \emph{outside} the link, as the textbooks and \code{circular::lm.circular}
#' write it -- distinct from \code{circ_gam}, which puts the intercept inside the
#' link). Because the reference direction \eqn{\mu_0} is a free angle, this model
#' is exactly rotation-equivariant and never sits against the tan-half wall, so
#' \code{circ_lm} has no \code{center} argument (the \code{circ_gam} counterpart
#' to that property); likewise the \code{cc} harmonic model recovers its mean by
#' \code{atan2}, with no wall. A one- or two-formula list selects the sub-model the same way
#' \code{circ_gam} reads a formula list:
#' \code{list(theta ~ x, ~ 1)} (or just \code{theta ~ x}) models the mean with
#' constant \eqn{\kappa}; \code{list(theta ~ 1, ~ z)} models
#' \eqn{\log\kappa = \alpha + z^\top\gamma} with constant \eqn{\mu};
#' \code{list(theta ~ x, ~ x)} is the mixed model. The mixed model ties
#' \eqn{\mu} and \eqn{\kappa} to one shared design. Any of these may carry several
#' covariates (\code{theta ~ x + z}); only \code{cc}/\code{lc} are single-predictor.
#' The mixed iteration starts cold (\eqn{\kappa \equiv 1}) by default; pass
#' \code{init = list(beta=, alpha=, gamma=)} to seed it from your own starting
#' values -- e.g. the estimates of separately fitted mean-only and kappa-only
#' models, the two-stage start Fisher (1993) Sec. 6.4.4 describes.
#'
#' \strong{cc -- circular on circular.} The Sarma and Jammalamadaka (1993) harmonic
#' fit: \code{cos(theta)} and \code{sin(theta)} regressed by least squares on a
#' degree-\code{order} trigonometric polynomial of the angular predictor and
#' reassembled as \eqn{\hat\mu = \mathrm{atan2}(\hat s, \hat c)}, with the circular
#' correlation \eqn{\rho}, a residual concentration, and the test for significance
#' of the next harmonic order.
#'
#' \strong{lc -- linear on circular.} Harmonic regression of a linear response on a
#' Fourier  basis of the angular predictor, reporting each harmonic's amplitude and
#' phase with delta-method standard errors.
#'
#' \strong{Parity with circular::lm.circular.} The regression outputs (\eqn{\beta},
#' \eqn{\mu}, the cos/sin coefficients, \eqn{\rho}, fitted values, p-values)
#' reproduce \code{circular} to machine precision. The reported \eqn{\kappa} can
#' differ slightly: it is the only quantity passing through the inverse Bessel
#' ratio, and circlss returns the machine-precision inverse (the exact
#' \eqn{\kappa} solving \eqn{A_1(\kappa) = R}) where \code{circular} uses the
#' classical piecewise approximation -- a gap at that approximation's error level
#' (~1e-3), largest at high concentration. One deliberate departure: the reported
#' \code{logLik} (and \code{aic}/\code{bic}) is the \emph{full} von Mises
#' log-likelihood with every estimated parameter counted (\eqn{\mu_0} and
#' \eqn{\kappa} included), so it exceeds \code{lm.circular}'s printed \code{log.lik}
#' by the \eqn{n\log 2\pi} normalisation \code{circular} drops -- putting circlss's
#' AIC on the standard scale, comparable to \code{\link{circ_gam}} or a \code{glm}.
#' @references
#' Fisher, N. I. and Lee, A. J. (1992) Regression models for an angular response.
#' \emph{Biometrics} 48, 665-677.
#'
#' Fisher, N. I. (1993) \emph{Statistical Analysis of Circular Data}. Cambridge
#' University Press.
#'
#' Sarma, Y. and Jammalamadaka, S. R. (1993) Circular regression. In
#' \emph{Statistical Sciences and Data Analysis}, 109-128. VSP, Utrecht.
#'
#' Pewsey, A., Neuhaeuser, M. and Ruxton, G. D. (2013) \emph{Circular Statistics
#' in R}. Oxford University Press.
#' @seealso
#' \code{\link{circ_gam}} for the penalized-spline distributional models;
#' \code{\link{vmlss}} for the von Mises family used there.
#' @examples
#' set.seed(1)
#' n <- 80
#' x <- rnorm(n)
#' theta <- (1 + 2 * atan(1.5 * x) + rnorm(n) / 4) %% (2 * pi)
#' dat <- data.frame(theta = theta, x = x)
#'
#' ## cl: Fisher-Lee mean direction (constant kappa)
#' m <- circ_lm(theta ~ x, dat, type = "cl")
#' m
#' predict(m, data.frame(x = c(-1, 0, 1)))
#'
#' ## cl: mixed model -- mean and log-kappa both linear in x
#' circ_lm(list(theta ~ x, ~ x), dat, type = "cl")
#'
#' ## cl: bootstrap SEs (Fisher 1993 Sec. 8.4) -- preferred at small n
#' set.seed(1)
#' circ_lm(theta ~ x, dat, type = "cl", se = "bootstrap", R = 199)
#'
#' ## cl: seed the mixed fit from separately fitted mean-only / kappa-only models
#' b0 <- circ_lm(theta ~ x, dat, type = "cl")
#' k0 <- circ_lm(list(theta ~ 1, ~ x), dat, type = "cl")
#' circ_lm(list(theta ~ x, ~ x), dat, type = "cl",
#'         init = list(beta = b0$beta, alpha = k0$alpha, gamma = k0$gamma))
#'
#' ## cl: several covariates (mean and concentration share the design)
#' dat$z <- rnorm(n)
#' circ_lm(list(theta ~ x + z, ~ x + z), dat, type = "cl")
#'
#' ## cc / lc: harmonic fits on an angular predictor
#' phi <- runif(n, 0, 2 * pi)
#' dcc <- data.frame(psi = (phi / 2 + rnorm(n) / 5) %% (2 * pi), phi = phi)
#' circ_lm(psi ~ phi, dcc, type = "cc", order = 1)
#'
#' dlc <- data.frame(y = 5 + 2 * cos(phi) + rnorm(n) / 2, phi = phi)
#' circ_lm(y ~ phi, dlc, type = "lc", order = 1)
#' @export
circ_lm <- function(formula, data, type = c("cl", "cc", "lc"),
                    order = 1L, init = NULL, tol = 1e-8, maxit = 1000L,
                    se = c("asymptotic", "bootstrap"), R = 999L,
                    verbose = FALSE) {
  if (missing(data))
    stop("circ_lm() requires 'data'.", call. = FALSE)
  type <- .circ_lm_type(type)
  se <- match.arg(se)
  if (se == "bootstrap" && type != "cl")
    stop("se = 'bootstrap' is implemented for type = 'cl' only; ",
         "cc/lc report least-squares standard errors.", call. = FALSE)

  flist <- if (inherits(formula, "formula")) list(formula) else formula
  if (!is.list(flist) || !length(flist) ||
      !all(vapply(flist, inherits, logical(1), "formula")))
    stop("'formula' must be a formula or a non-empty list of formulas.",
         call. = FALSE)
  .circ_lm_no_smooth(flist)

  fit <- switch(type,
    cl = .circ_lm_cl(flist, data, init, tol, as.integer(maxit), verbose, se, as.integer(R)),
    cc = .circ_lm_cc(flist, data, as.integer(order)),
    lc = .circ_lm_lc(flist, data, as.integer(order)))
  fit$type <- type
  fit$call <- match.call()
  class(fit) <- "circ_lm"
  fit
}

## --- shared front-door helpers ---------------------------------------------

## Normalise the leg selector: accept "cl"/"c-l"/"C-L" etc.
.circ_lm_type <- function(type) {
  key <- gsub("-", "", tolower(type[1L]))
  if (!key %in% c("cl", "cc", "lc"))
    stop("type must be one of 'cl'/'c-l' (circular ~ linear), ",
         "'cc'/'c-c' (circular ~ circular), or 'lc'/'l-c' (linear ~ circular); got ",
         sQuote(type[1L]), ".", call. = FALSE)
  key
}

## Parametric-only guard: reject mgcv smooth constructors anywhere on a RHS.
.circ_lm_no_smooth <- function(flist) {
  bad <- character(0)
  walk <- function(e) {
    if (!is.call(e)) return(invisible())
    fn <- e[[1L]]
    if (is.symbol(fn) && as.character(fn) %in% c("s", "te", "ti", "t2"))
      bad <<- c(bad, as.character(fn))
    else for (i in seq_along(e)) walk(e[[i]])
  }
  for (f in flist) walk(f[[length(f)]])
  if (length(bad))
    stop("circ_lm() fits parametric models only; smooth term(s) ",
         paste0(unique(bad), "()", collapse = ", "),
         " are not allowed. Use circ_gam() for penalized smooths.",
         call. = FALSE)
  invisible()
}

## RHS model matrix WITHOUT an intercept column (Fisher-Lee carries the
## intercept outside the per-covariate design). `~ 1` yields a 0-column matrix.
.circ_lm_design <- function(formula, data) {
  tl <- attr(stats::terms(formula), "term.labels")
  if (!length(tl)) return(matrix(numeric(0), nrow = nrow(data), ncol = 0L))
  stats::model.matrix(stats::reformulate(tl, intercept = FALSE), data)
}

## The single angular predictor of a one-sided/two-sided "y ~ x" formula.
.circ_lm_one_predictor <- function(formula, data) {
  if (length(formula) != 3L)
    stop("the formula must name the response, e.g. theta ~ phi.", call. = FALSE)
  resp <- deparse(formula[[2L]])
  vars <- all.vars(formula[[3L]])
  if (length(vars) != 1L)
    stop("circ_lm(type = 'cc'/'lc') takes exactly one angular predictor; got ",
         if (length(vars)) paste(sQuote(vars), collapse = ", ") else "none",
         ". For several covariates use circ_gam().", call. = FALSE)
  keep <- stats::complete.cases(data[c(resp, vars)])
  list(response = resp, var = vars, y = as.numeric(data[[resp]])[keep],
       x = as.numeric(data[[vars]])[keep], n = sum(keep))
}

## Harmonic design data frame: cos(k.x), sin(k.x) for k = 1..order, with
## syntactic names cos1, sin1, cos2, ... (mapped to readable labels at print).
.circ_lm_harmonic <- function(x, order) {
  x <- as.numeric(x)
  out <- vector("list", 2L * order)
  nm <- character(2L * order)
  j <- 1L
  for (k in seq_len(order)) {
    out[[j]] <- cos(k * x); nm[j] <- paste0("cos", k); j <- j + 1L
    out[[j]] <- sin(k * x); nm[j] <- paste0("sin", k); j <- j + 1L
  }
  d <- as.data.frame(out)
  names(d) <- nm
  d
}

## ===========================================================================
## cl -- Fisher-Lee von Mises regression (mean / kappa / mixed)
## ===========================================================================

## Free-parameter count for AIC/BIC: every estimated parameter, not just the
## regression slopes -- mean has (mu0, beta_1..p, kappa); kappa has (mu, alpha,
## gamma_1..p); mixed has (mu0, beta_1..p, alpha, gamma_1..p). Standard AIC, so it
## matches a circ_gam()/glm coefficient count (lm.circular reports no AIC).
.circ_lm_cl_npar <- function(model, p)
  switch(model, mean = p + 2L, kappa = p + 2L, mixed = 2L * p + 2L)

.circ_lm_cl <- function(flist, data, init, tol, maxit, verbose,
                        se = "asymptotic", R = 999L) {
  if (length(flist) > 2L)
    stop("circ_lm(type = 'cl') takes at most two formulas: ",
         "list(mu ~ ..., logkappa ~ ...).", call. = FALSE)
  mu_f <- flist[[1L]]
  if (length(mu_f) != 3L)
    stop("the first formula must name the response, e.g. theta ~ x.",
         call. = FALSE)
  kappa_f <- if (length(flist) == 2L) flist[[2L]] else stats::as.formula("~ 1")

  vars <- unique(c(all.vars(mu_f), all.vars(kappa_f)))
  keep <- stats::complete.cases(data[vars])
  data <- data[keep, , drop = FALSE]

  theta <- as.numeric(eval(mu_f[[2L]], data, parent.frame())) %% (2 * pi)
  Xmu <- .circ_lm_design(mu_f, data)
  Xka <- .circ_lm_design(kappa_f, data)
  has_mu <- ncol(Xmu) > 0L
  has_ka <- ncol(Xka) > 0L
  model <- if (has_mu && !has_ka) "mean"
           else if (!has_mu && has_ka) "kappa"
           else if (has_mu && has_ka) "mixed"
           else stop("no predictors in either formula; nothing to regress.",
                     call. = FALSE)
  if (model == "mixed" && !identical(colnames(Xmu), colnames(Xka)))
    stop("the Fisher-Lee fitter ties mu and kappa to one shared design, but ",
         "mu uses (", paste(colnames(Xmu), collapse = ", "), ") and kappa uses (",
         paste(colnames(Xka), collapse = ", "),
         "). Use circ_gam() for different covariates per predictor.",
         call. = FALSE)

  X <- if (model == "kappa") Xka else Xmu       # shared design (mean/mixed use Xmu)
  fit <- .circ_lm_cl_fit(theta, X, model, init, tol, maxit, verbose)

  npar <- .circ_lm_cl_npar(model, ncol(X))
  fit$n <- length(theta)
  fit$npar <- npar
  fit$aic <- -2 * fit$loglik + 2 * npar
  fit$bic <- -2 * fit$loglik + log(fit$n) * npar
  fit$model <- model
  fit$response <- deparse(mu_f[[2L]])
  fit$mu_formula <- mu_f
  fit$kappa_formula <- kappa_f
  fit$mu_terms <- colnames(Xmu)
  fit$kappa_terms <- colnames(Xka)
  fit$fitted <- if (model == "kappa") rep(fit$mu %% (2 * pi), length(theta))
                else (fit$mu + 2 * atan(as.vector(Xmu %*% fit$beta))) %% (2 * pi)
  fit$residuals <- (theta - fit$fitted) %% (2 * pi)
  ## plotting frame: the covariate column(s) and the (mod 2pi) response, so
  ## plot.circ_lm() can overlay the data and bound the covariate grid.
  keep <- intersect(names(data), unique(c(all.vars(mu_f[[3L]]), all.vars(kappa_f))))
  frame <- data[, keep, drop = FALSE]
  frame[[fit$response]] <- theta
  fit$frame <- frame

  fit$se_method <- se
  if (se == "bootstrap") {
    bse <- .circ_lm_cl_boot(fit, X, model, tol, maxit, R)
    for (.nm in names(bse)) fit[[.nm]] <- bse[[.nm]]   # override se_*/V* with bootstrap, add nboot
  }
  fit
}

## Resolve `init` into starting values (beta, alpha, gamma) of the right lengths.
## init: NULL -> cold (beta = gamma = 0, alpha = 0); list(beta=, alpha=, gamma=)
## -> explicit starts, any omitted component falling back to cold; bare numeric
## -> the mean-direction coefficients (beta). Lets a caller hand in estimates from
## separately fitted mean-only / kappa-only models as the joint iteration's start.
.circ_lm_cl_starts <- function(X, model, init) {
  p <- ncol(X)
  cold <- list(beta = rep(0, p), gamma = rep(0, p), alpha = 0)
  if (is.null(init)) return(cold)

  ## explicit named list: take what is given, validate lengths, warn on misuse.
  if (is.list(init)) {
    bad <- setdiff(names(init), c("beta", "alpha", "gamma"))
    if (length(bad))
      warning("ignoring unknown init component(s): ", paste(bad, collapse = ", "),
              ".", call. = FALSE)
    g <- cold
    if (!is.null(init$beta)) {
      if (model == "kappa")
        warning("init$beta is unused when only the concentration is modelled.",
                call. = FALSE)
      else if (length(init$beta) != p)
        stop(sprintf("init$beta must have length %d (one per mean covariate).", p),
             call. = FALSE)
      else g$beta <- as.numeric(init$beta)
    }
    if (!is.null(init$gamma)) {
      if (model == "mean")
        warning("init$gamma is unused when only the mean direction is modelled.",
                call. = FALSE)
      else if (length(init$gamma) != p)
        stop(sprintf("init$gamma must have length %d (one per concentration covariate).", p),
             call. = FALSE)
      else g$gamma <- as.numeric(init$gamma)
    }
    if (!is.null(init$alpha)) {
      if (model == "mean")
        warning("init$alpha is unused when only the mean direction is modelled.",
                call. = FALSE)
      else if (length(init$alpha) != 1L)
        stop("init$alpha must be a single number.", call. = FALSE)
      else g$alpha <- as.numeric(init$alpha)
    }
    return(g)
  }

  ## legacy: a bare numeric vector seeds the mean-direction coefficients.
  if (model == "kappa")
    warning("numeric 'init' seeds the mean-direction coefficients, unused for a ",
            "concentration-only model; use init = list(alpha=, gamma=).",
            call. = FALSE)
  g <- cold
  g$beta <- rep_len(as.numeric(init), p)
  g
}

## Green (1984) IRLS for the von Mises MLE. Mean: mu_i = mu0 + 2*atan(X beta),
## constant kappa. Kappa: log kappa_i = alpha + X gamma, constant mu. Mixed: both.
## Scores and expected-information weights follow Fisher (1993) Sec. 6.4 (the
## same recipe pycircstat2's CLRegression uses, validated against circular).
.circ_lm_cl_fit <- function(theta, X, model, init, tol, maxit, verbose) {
  n <- length(theta)
  p <- ncol(X)
  X1 <- cbind(1, X)
  linkinv <- function(eta) 2 * atan(eta)
  mu_eta  <- function(eta) 2 / (1 + eta^2)
  log_i0  <- function(k) k + log(besselI(k, 0, expon.scaled = TRUE))
  ## kappa = exp(eta), capped where the scaled Bessel ratio stays finite:
  ## besselI(., expon.scaled) underflows to 0 (=> A1, log_i0 -> NaN) near
  ## kappa ~ 2e5, so hold eta <= 11 (kappa <~ 6e4). Real fits stay far below this;
  ## the cap only fences a divergent IRLS step (e.g. from a poor user-supplied init).
  kappa_of <- function(eta) exp(pmin(pmax(eta, -50), 11))
  ## a singular weighted design (a divergent step from a poor init) must not
  ## crash the fit: signal it with NA so the loop bails and reports non-convergence.
  safe_solve <- function(A, b) tryCatch(as.vector(solve(A, b)),
                                        error = function(e) rep(NA_real_, ncol(A)))
  ridgeX  <- 1e-8 * diag(p)
  ridge1  <- 1e-8 * diag(p + 1L)

  st    <- .circ_lm_cl_starts(X, model, init)
  beta  <- st$beta
  gamma <- st$gamma
  alpha <- st$alpha
  mu <- 0; kappa <- 1
  ll_old <- -Inf; diff <- tol + 1; iter <- 0L; diverged <- FALSE

  for (iter in seq_len(maxit)) {
    if (model == "mean") {
      eta <- as.vector(X %*% beta)
      rdev <- theta - linkinv(eta)
      S <- mean(sin(rdev)); C <- mean(cos(rdev))
      mu <- atan2(S, C); kappa <- A1inv(sqrt(S^2 + C^2))
      G <- mu_eta(eta) * X
      w <- kappa * A1(kappa)
      u <- kappa * sin(rdev - mu)
      GtG <- crossprod(G)
      beta <- safe_solve(w * GtG + ridgeX, crossprod(G, u) + w * GtG %*% beta)
      ll <- -n * log_i0(kappa) + kappa * sum(cos(rdev - mu))

    } else if (model == "kappa") {
      kappa <- kappa_of(alpha + as.vector(X %*% gamma))
      mu <- atan2(sum(kappa * sin(theta)), sum(kappa * cos(theta)))
      a1p <- pmax(A1prime(kappa), 1e-12)
      y <- (cos(theta - mu) - A1(kappa)) / (a1p * kappa)
      w <- kappa^2 * a1p
      upd <- safe_solve(crossprod(X1, w * X1) + ridge1, crossprod(X1, w * y))
      alpha <- alpha + upd[1L]; gamma <- gamma + upd[-1L]
      ll <- -sum(log_i0(kappa)) + sum(kappa * cos(theta - mu))

    } else {                                    # mixed
      kappa <- kappa_of(alpha + as.vector(X %*% gamma))
      eta <- as.vector(X %*% beta)
      rdev <- theta - linkinv(eta)
      mu <- atan2(sum(kappa * sin(rdev)), sum(kappa * cos(rdev)))
      G <- mu_eta(eta) * X
      wb <- kappa * A1(kappa)
      GtWG <- crossprod(G, wb * G)
      beta <- safe_solve(GtWG + ridgeX, crossprod(G, kappa * sin(rdev - mu)) + GtWG %*% beta)
      a1p <- pmax(A1prime(kappa), 1e-12)
      y <- (cos(rdev - mu) - A1(kappa)) / (a1p * kappa)
      wg <- kappa^2 * a1p
      upd <- safe_solve(crossprod(X1, wg * X1) + ridge1, crossprod(X1, wg * y))
      alpha <- alpha + upd[1L]; gamma <- gamma + upd[-1L]
      ll <- -sum(log_i0(kappa)) + sum(kappa * cos(rdev - mu))
    }

    if (!is.finite(ll)) { diverged <- TRUE; break }   # singular step / overflow: bail cleanly
    diff <- abs(ll - ll_old)
    if (verbose)
      cat(sprintf("iter %d: logLik = %.6f, diff = %.2e\n", iter, ll, diff))
    if (diff < tol) break
    ll_old <- ll
  }
  converged <- !diverged && is.finite(diff) && diff < tol
  if (!converged)
    warning(if (diverged)
              "circ_lm(type = 'cl') diverged before converging; try init = NULL (cold) or other starts."
            else sprintf("circ_lm(type = 'cl') did not converge in %d iterations (last diff = %.2e, tol = %.2e).",
                         maxit, diff, tol),
            call. = FALSE)

  ll <- ll - n * log(2 * pi)   # the proper vM normalisation; the loop dropped this
                               # constant (as lm.circular does), restore it so logLik
                               # is a real log-likelihood and AIC/BIC are comparable.
  se <- tryCatch(.circ_lm_cl_se(theta, X, model, beta, alpha, gamma, kappa, mu),
                 error = function(e)
                   list(se_beta = rep(NA_real_, length(beta)), se_gamma = rep(NA_real_, length(gamma)),
                        se_alpha = NA_real_, se_mu = NA_real_, se_kappa = NA_real_,
                        Vbeta = NULL, Vag = NULL))
  list(mu = mu, kappa = kappa, beta = beta, alpha = alpha, gamma = gamma,
       loglik = ll, iter = iter, converged = converged,
       se_beta = se$se_beta, se_gamma = se$se_gamma, se_alpha = se$se_alpha,
       se_mu = se$se_mu, se_kappa = se$se_kappa, Vbeta = se$Vbeta, Vag = se$Vag)
}

## Large-sample SEs from the expected information (Fisher 1993, eq. 6.62-6.64,
## 6.82). Per-observation kappa_i SE (kappa/mixed) by the delta method on (a, g).
.circ_lm_cl_se <- function(theta, X, model, beta, alpha, gamma, kappa, mu) {
  n <- length(theta)
  p <- ncol(X)
  X1 <- cbind(1, X)
  out <- list()
  if (model == "mean") {
    eta <- as.vector(X %*% beta)
    G <- (2 / (1 + eta^2)) * X
    w <- kappa * A1(kappa)
    cov_b <- solve(w * crossprod(G))
    out$se_beta <- sqrt(diag(cov_b))
    out$Vbeta <- cov_b
    out$se_mu <- 1 / sqrt(max((n - p) * w, 1e-12))
    out$se_kappa <- sqrt(1 / max(n * (1 - A1(kappa)^2 - A1(kappa) / kappa), 1e-12))
  } else {
    if (model == "mixed") {
      eta <- as.vector(X %*% beta)
      G <- (2 / (1 + eta^2)) * X
      cov_b <- solve(crossprod(G, (kappa * A1(kappa)) * G))
      out$se_beta <- sqrt(diag(cov_b))
      out$Vbeta <- cov_b
    }
    w <- kappa^2 * A1prime(kappa)
    cov_ag <- solve(crossprod(X1, w * X1))
    out$Vag <- cov_ag
    out$se_alpha <- sqrt(cov_ag[1L, 1L])
    out$se_gamma <- sqrt(diag(cov_ag)[-1L])
    out$se_mu <- 1 / sqrt(max(sum(kappa * A1(kappa)) - 0.5, 1e-12))
    quad <- rowSums((X1 %*% cov_ag) * X1)
    out$se_kappa <- kappa * sqrt(pmax(quad, 0))
  }
  out
}

## Parametric bootstrap SEs (Fisher 1993, sec 8.4), the small-sample alternative
## to the expected-information SEs above: simulate responses from the fitted
## vM(mu_i, kappa_i), refit, and take the spread of the bootstrap estimates.
## Refits warm-start from the point estimate (the resample sits near it) so they
## converge fast; non-converged draws are dropped. Reproducibility is the caller's
## set.seed(). Returns the same se_*/V* fields the asymptotic path fills, by model.
.circ_lm_cl_boot <- function(fit, X, model, tol, maxit, R) {
  mu_i <- fit$fitted
  k_i  <- if (length(fit$kappa) == 1L) rep_len(fit$kappa, length(mu_i)) else fit$kappa
  p <- ncol(X)
  start <- switch(model,
    mean  = list(beta = fit$beta),
    kappa = list(alpha = fit$alpha, gamma = fit$gamma),
    mixed = list(beta = fit$beta, alpha = fit$alpha, gamma = fit$gamma))
  ## refits need only SE-grade precision, and the mixed IRLS converges slowly
  ## (kappa-leverage ill-conditioning), so loosen tol and lift the iteration cap
  ## off the user's point-fit defaults -- otherwise slow resamples hit maxit, get
  ## dropped, and the survivors bias the SE downward.
  boot_tol <- max(tol, 1e-6); boot_maxit <- max(maxit, 1000L)
  reps <- vector("list", R)
  for (b in seq_len(R)) {
    yb <- .circ_rvonmises(mu_i, k_i) %% (2 * pi)
    fb <- tryCatch(suppressWarnings(.circ_lm_cl_fit(yb, X, model, start, boot_tol, boot_maxit, FALSE)),
                   error = function(e) NULL)
    if (!is.null(fb) && isTRUE(fb$converged)) reps[[b]] <- fb
  }
  reps <- reps[!vapply(reps, is.null, logical(1))]
  nb <- length(reps)
  if (nb < 2L)
    stop("se = 'bootstrap': only ", nb, " of ", R, " resamples converged; ",
         "use se = 'asymptotic'.", call. = FALSE)
  if (nb < R %/% 2L)
    warning(sprintf("se = 'bootstrap': only %d of %d resamples converged.", nb, R),
            call. = FALSE)

  col <- function(nm) vapply(reps, `[[`, numeric(1), nm)                    # scalar per rep
  mat <- function(nm) do.call(rbind, lapply(reps, `[[`, nm))               # nb x p (robust at p = 1)
  Mu  <- col("mu")
  out <- list(nboot = nb,
              se_mu = stats::sd(atan2(sin(Mu - fit$mu), cos(Mu - fit$mu))))  # circular SE of mu
  if (model %in% c("mean", "mixed")) {
    Beta <- mat("beta")
    out$se_beta <- apply(Beta, 2L, stats::sd); out$Vbeta <- stats::cov(Beta)
  }
  if (model %in% c("kappa", "mixed")) {
    Alpha <- col("alpha"); Gamma <- mat("gamma")
    out$se_alpha <- stats::sd(Alpha); out$se_gamma <- apply(Gamma, 2L, stats::sd)
    out$Vag <- stats::cov(cbind(Alpha, Gamma))
    Kmat <- exp(sweep(Gamma %*% t(X), 1L, Alpha, "+"))   # nb x n: per-observation kappa
    out$se_kappa <- apply(Kmat, 2L, stats::sd)
  } else {
    out$se_kappa <- stats::sd(col("kappa"))              # mean model: one kappa
  }
  out
}

## Best (1979) / Fisher rejection sampler for the von Mises, vectorised over
## (mu, kappa) by recycling; one draw per element. kappa ~ 0 gives a uniform angle.
.circ_rvonmises <- function(mu, kappa) {
  n <- max(length(mu), length(kappa))
  mu <- rep_len(mu, n); kappa <- rep_len(kappa, n)
  out <- numeric(n)
  for (i in seq_len(n)) {
    k <- kappa[i]
    if (k < 1e-8) { out[i] <- runif(1, -pi, pi) + mu[i]; next }
    a <- 1 + sqrt(1 + 4 * k * k); b <- (a - sqrt(2 * a)) / (2 * k); r <- (1 + b * b) / (2 * b)
    repeat {
      z <- cos(pi * runif(1)); f <- (1 + r * z) / (r + z); cc <- k * (r - f); u <- runif(1)
      if (cc * (2 - cc) - u > 0 || log(cc / u) + 1 - cc >= 0) {
        out[i] <- sign(runif(1) - 0.5) * acos(max(min(f, 1), -1)) + mu[i]; break
      }
    }
  }
  out
}

## ===========================================================================
## cc -- Sarma & Jammalamadaka harmonic circular-circular regression
## ===========================================================================

.circ_lm_cc <- function(flist, data, order) {
  if (length(flist) != 1L)
    stop("circ_lm(type = 'cc') takes a single formula, e.g. theta ~ phi.",
         call. = FALSE)
  s <- .circ_lm_one_predictor(flist[[1L]], data)
  x <- s$x %% (2 * pi); y <- s$y %% (2 * pi); n <- s$n
  H <- .circ_lm_harmonic(x, order)
  cy <- cos(y); sy <- sin(y)
  d <- cbind(data.frame(.cos_y = cy, .sin_y = sy), H)
  rhs <- names(H)
  cos_lm <- stats::lm(stats::reformulate(rhs, ".cos_y"), d)
  sin_lm <- stats::lm(stats::reformulate(rhs, ".sin_y"), d)

  cos_fit <- stats::fitted(cos_lm); sin_fit <- stats::fitted(sin_lm)
  fitted <- atan2(sin_fit, cos_fit) %% (2 * pi)
  residuals <- (y - fitted) %% (2 * pi)
  rho <- as.numeric(sqrt((crossprod(cos_fit) + crossprod(sin_fit)) / n))
  A_k <- mean(cos(residuals))
  if (A_k < 0)
    warning("mean residual cosine is negative; residuals anti-align with the ",
            "fit (kappa clamped to 0). Check for misspecification.", call. = FALSE)
  kappa <- A1inv(A_k)

  ## Higher-order test (Sarma & Jammalamadaka 1993): do the (order + 1) harmonics
  ## add signal beyond the fitted design? One statistic each for cos(y), sin(y).
  ones <- matrix(1, n, 1L)
  Xm <- cbind(ones, as.matrix(H))
  W <- cbind(cos((order + 1) * x), sin((order + 1) * x))
  M <- Xm %*% solve(crossprod(Xm), t(Xm))
  IM <- diag(n) - M
  N <- W %*% solve(crossprod(W, IM %*% W), t(W))
  cc <- n - (2 * order + 1)
  num1 <- as.numeric(crossprod(cy, IM %*% N %*% IM %*% cy))
  den1 <- as.numeric(crossprod(cy, IM %*% cy))
  num2 <- as.numeric(crossprod(sy, IM %*% N %*% IM %*% sy))
  den2 <- as.numeric(crossprod(sy, IM %*% sy))
  T1 <- cc * num1 / den1; T2 <- cc * num2 / den2
  p_values <- c(cos = 1 - stats::pchisq(T1, 2), sin = 1 - stats::pchisq(T2, 2))

  list(var = s$var, response = s$response, order = order, n = n,
       cos_lm = cos_lm, sin_lm = sin_lm,
       coefficients = cbind(cos = stats::coef(cos_lm), sin = stats::coef(sin_lm)),
       rho = rho, A_k = A_k, kappa = kappa,
       fitted = fitted, residuals = residuals, p_values = p_values,
       frame = stats::setNames(data.frame(x, y), c(s$var, s$response)))
}

## ===========================================================================
## lc -- harmonic linear-circular regression
## ===========================================================================

.circ_lm_lc <- function(flist, data, order) {
  if (length(flist) != 1L)
    stop("circ_lm(type = 'lc') takes a single formula, e.g. y ~ phi.",
         call. = FALSE)
  s <- .circ_lm_one_predictor(flist[[1L]], data)
  x <- s$x %% (2 * pi); y <- s$y; n <- s$n
  H <- .circ_lm_harmonic(x, order)
  d <- cbind(data.frame(.y = y), H)
  fit <- stats::lm(stats::reformulate(names(H), ".y"), d)

  cf <- stats::coef(fit); V <- stats::vcov(fit)
  ## amplitude sqrt(a^2 + b^2) and phase atan2(b, a) per harmonic, delta-method SEs
  harm <- lapply(seq_len(order), function(k) {
    cn <- paste0("cos", k); sn <- paste0("sin", k)
    a <- cf[[cn]]; b <- cf[[sn]]
    amp <- sqrt(a^2 + b^2); phase <- atan2(b, a)
    vc <- V[cn, cn]; vs <- V[sn, sn]; cs <- V[cn, sn]; r2 <- amp^2
    se_amp <- if (amp > 0)
      sqrt(max((a^2 * vc + 2 * a * b * cs + b^2 * vs) / r2, 0)) else NA_real_
    se_phase <- if (amp > 0)
      sqrt(max((b^2 * vc - 2 * a * b * cs + a^2 * vs) / r2^2, 0)) else NA_real_
    c(k = k, cos = a, sin = b, amplitude = amp, phase = phase,
      se_amplitude = se_amp, se_phase = se_phase)
  })
  sm <- summary(fit)

  list(var = s$var, response = s$response, order = order, n = n, lm = fit,
       coefficients = cf, harmonics = do.call(rbind, harm),
       sigma = sm$sigma, r_squared = sm$r.squared,
       loglik = as.numeric(stats::logLik(fit)), npar = length(cf) + 1L,
       aic = stats::AIC(fit), bic = stats::BIC(fit),
       fitted = stats::fitted(fit), residuals = stats::residuals(fit),
       frame = stats::setNames(data.frame(x, y), c(s$var, s$response)))
}

## --- methods ---------------------------------------------------------------

#' @rdname circ_lm
#' @export
coef.circ_lm <- function(object, ...) object$coefficients

#' @rdname circ_lm
#' @export
fitted.circ_lm <- function(object, ...) object$fitted

#' @rdname circ_lm
#' @export
residuals.circ_lm <- function(object, ...) object$residuals

#' @rdname circ_lm
#' @export
logLik.circ_lm <- function(object, ...) {
  if (is.null(object$loglik))
    stop("type = 'cc' is two separate least-squares fits and has no single ",
         "log-likelihood.", call. = FALSE)
  val <- object$loglik
  attr(val, "df") <- if (!is.null(object$npar)) object$npar
                     else length(object$coefficients)
  attr(val, "nobs") <- object$n
  class(val) <- "logLik"
  val
}

#' @rdname circ_lm
#' @export
predict.circ_lm <- function(object, newdata, type = c("direction", "kappa"), ...) {
  if (missing(newdata))
    return(object$fitted)
  newdata <- as.data.frame(newdata)
  if (object$type == "cc") {
    H <- .circ_lm_harmonic(as.numeric(newdata[[object$var]]) %% (2 * pi), object$order)
    return(atan2(stats::predict(object$sin_lm, H),
                 stats::predict(object$cos_lm, H)) %% (2 * pi))
  }
  if (object$type == "lc") {
    H <- .circ_lm_harmonic(as.numeric(newdata[[object$var]]) %% (2 * pi), object$order)
    return(as.numeric(stats::predict(object$lm, H)))
  }
  ## cl
  type <- match.arg(type)
  if (type == "kappa") {
    if (object$model == "mean")
      return(rep(object$kappa, nrow(newdata)))
    Z <- .circ_lm_design(object$kappa_formula, newdata)
    return(exp(object$alpha + as.vector(Z %*% object$gamma)))
  }
  if (object$model == "kappa")
    return(rep(object$mu %% (2 * pi), nrow(newdata)))
  Xn <- .circ_lm_design(object$mu_formula, newdata)
  (object$mu + 2 * atan(as.vector(Xn %*% object$beta))) %% (2 * pi)
}

#' @rdname circ_lm
#' @export
print.circ_lm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n  ", paste(deparse(x$call), collapse = " "), "\n", sep = "")
  switch(x$type,
    cl = .circ_lm_print_cl(x, digits),
    cc = .circ_lm_print_cc(x, digits),
    lc = .circ_lm_print_lc(x, digits))
  invisible(x)
}

.circ_lm_coefmat <- function(est, se, names) {
  z <- est / se
  p <- 2 * stats::pnorm(-abs(z))
  m <- cbind(Estimate = est, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = p)
  rownames(m) <- names
  m
}

.circ_lm_print_cl <- function(x, digits) {
  head <- switch(x$model,
    mean  = "Circular-linear regression (Fisher-Lee), mean direction",
    kappa = "Circular-linear regression (Fisher-Lee), concentration",
    mixed = "Circular-linear regression (Fisher-Lee), mean and concentration")
  cat("\n", head, "\n", sep = "")
  if (x$model %in% c("mean", "mixed")) {
    cat("\nMean direction  mu = mu0 + 2*atan(X beta):\n")
    stats::printCoefmat(.circ_lm_coefmat(x$beta, x$se_beta, x$mu_terms),
                        digits = digits, has.Pvalue = TRUE)
  }
  if (x$model %in% c("kappa", "mixed")) {
    cat("\nConcentration  log(kappa) = alpha + Z gamma:\n")
    stats::printCoefmat(
      .circ_lm_coefmat(c(x$alpha, x$gamma), c(x$se_alpha, x$se_gamma),
                       c("(Intercept)", x$kappa_terms)),
      digits = digits, has.Pvalue = TRUE)
  }
  cat("\nmu0: ", format(x$mu, digits = digits),
      " (", format(x$se_mu, digits = digits), ")", sep = "")
  if (length(x$kappa) == 1L)
    cat("    kappa: ", format(x$kappa, digits = digits),
        " (", format(x$se_kappa, digits = digits), ")\n", sep = "")
  else
    cat(sprintf("    kappa: %.*f to %.*f (per observation)\n",
                digits, min(x$kappa), digits, max(x$kappa)))
  cat(sprintf("logLik: %s   AIC: %s   BIC: %s   n: %d\n",
              format(x$loglik, digits = digits), format(x$aic, digits = digits),
              format(x$bic, digits = digits), x$n))
  if (!isTRUE(x$converged)) cat("** did not converge **\n")
  if (identical(x$se_method, "bootstrap"))
    cat(sprintf("SEs: parametric bootstrap, %d resamples.   ", x$nboot))
  cat("p-values use the normal approximation.\n")
}

.circ_lm_print_cc <- function(x, digits) {
  cat("\nCircular-circular regression (harmonic, order ", x$order, ")\n",
      "  ", x$response, " ~ ", x$var, "   n = ", x$n, "\n\nCoefficients:\n",
      sep = "")
  print(round(x$coefficients, digits))
  cat("\nrho: ", format(x$rho, digits = digits),
      "    residual kappa: ", format(x$kappa, digits = digits),
      " (A_k = ", format(x$A_k, digits = digits), ")\n", sep = "")
  cat("Higher-order test p-values:  cos = ", format(x$p_values[["cos"]], digits = digits),
      ", sin = ", format(x$p_values[["sin"]], digits = digits), "\n", sep = "")
  cat(if (all(x$p_values > 0.05))
        "Higher-order terms not significant at the 0.05 level.\n"
      else "Higher-order terms significant at the 0.05 level.\n")
}

.circ_lm_print_lc <- function(x, digits) {
  cat("\nLinear-circular regression (harmonic, order ", x$order, ")\n",
      "  ", x$response, " ~ ", x$var, "   n = ", x$n, "\n\n", sep = "")
  print(summary(x$lm)$coefficients, digits = digits)
  H <- x$harmonics
  cat("\nHarmonic amplitude / phase:\n")
  amp <- H[, "amplitude", drop = TRUE]
  ph  <- H[, "phase", drop = TRUE]
  tab <- cbind(amplitude = amp, `se(amp)` = H[, "se_amplitude"],
               phase = ph, `se(phase)` = H[, "se_phase"])
  rownames(tab) <- paste0("k=", H[, "k"])
  print(round(tab, digits))
  cat(sprintf("\nsigma: %s   R-squared: %s   AIC: %s   BIC: %s   n: %d\n",
              format(x$sigma, digits = digits), format(x$r_squared, digits = digits),
              format(x$aic, digits = digits), format(x$bic, digits = digits), x$n))
}
