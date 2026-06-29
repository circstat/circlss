## Fisher-Lee tan-half link for a circular location parameter:
##   eta = g(mu) = tan(mu/2),  mu = g^{-1}(eta) = 2*atan(eta) in (-pi, pi).
## Derivative conventions are
## mgcv's: mu.eta is d mu / d eta as a function of eta; d2link, d3link and
## d4link are the 2nd-4th derivatives of the link g with respect to mu, as
## functions of mu. In t = tan(mu/2) (so dt/dmu = (1+t^2)/2):
##   g'  = (1+t^2)/2,            g'' = t(1+t^2)/2,
##   g''' = (1+t^2)(1+3t^2)/4,   g'''' = t(1+t^2)(2+3t^2)/2.
## The antipode mu = pi (mod 2pi) is the one singularity: unrepresentable
## under this link, the standard Fisher-Lee caveat.
tanhalf.link <- function() {
  list(
    link = "tanhalf",
    linkfun = function(mu) tan(0.5 * mu),
    linkinv = function(eta) 2 * atan(eta),
    mu.eta = function(eta) 2 / (1 + eta * eta),
    d2link = function(mu) {
      t <- tan(0.5 * mu)
      0.5 * t * (1 + t * t)
    },
    d3link = function(mu) {
      t2 <- tan(0.5 * mu)^2
      0.25 * (1 + t2) * (1 + 3 * t2)
    },
    d4link = function(mu) {
      t <- tan(0.5 * mu)
      t2 <- t * t
      0.5 * t * (1 + t2) * (2 + 3 * t2)
    },
    valideta = function(eta) TRUE
  )
}

## Logit-half link for a concentration parameter bounded to (0, 1/2) -- the
## cardioid's mean resultant length rho. With h = 1/2:
##   eta = g(rho) = log(rho / (h - rho)),  rho = g^{-1}(eta) = h * plogis(eta).
## Conventions are
## mgcv's, as in tanhalf.link(): mu.eta is d rho / d eta as a function of eta;
## d2link, d3link, d4link are the 2nd-4th derivatives of g with respect to rho,
## as functions of rho:
##   g'  = 1/rho + 1/(h-rho),        g'' = 1/(h-rho)^2 - 1/rho^2,
##   g''' = 2/(h-rho)^3 + 2/rho^3,   g'''' = 6/(h-rho)^4 - 6/rho^4.
## rho is clamped to the open interval (clip the logistic to [eps, 1-eps]);
## mu.eta is floored at the machine epsilon, like
## stats::make.link("logit"). The cardioid is the only family using this link.
logithalf.link <- function() {
  hi <- 0.5
  eps <- .Machine$double.eps
  list(
    link = "logithalf",
    linkfun = function(mu) log(mu / (hi - mu)),
    linkinv = function(eta) hi * pmin(pmax(stats::plogis(eta), eps), 1 - eps),
    mu.eta = function(eta) {
      s <- stats::plogis(eta)
      pmax(hi * s * (1 - s), eps)
    },
    d2link = function(mu) 1 / (hi - mu)^2 - 1 / mu^2,
    d3link = function(mu) 2 / (hi - mu)^3 + 2 / mu^3,
    d4link = function(mu) 6 / (hi - mu)^4 - 6 / mu^4,
    valideta = function(eta) TRUE
  )
}

## Tanh link for a shape parameter bounded to (-1, 1) -- the sine-skew lambda.
##   eta = g(lambda) = atanh(lambda),  lambda = g^{-1}(eta) = tanh(eta).
## Conventions are mgcv's, as in the
## other links: mu.eta is d lambda / d eta as a function of eta (= sech^2 eta,
## written overflow-safe and floored at machine epsilon); d2link, d3link,
## d4link are the 2nd-4th derivatives of g with respect to lambda:
##   g'  = 1/(1-lambda^2),              g'' = 2 lambda/(1-lambda^2)^2,
##   g''' = (2+6 lambda^2)/(1-lambda^2)^3,
##   g'''' = 24 lambda (1+lambda^2)/(1-lambda^2)^4.
## lambda is clamped inside the open interval (clip tanh to [-1+eps, 1-eps]) as
## the reference does -- at lambda = +/-1 the sine-skewed density touches 0 and
## its log has no finite limit. The sine-skewed Jones-Pewsey is the only family
## using this link.
tanh_link <- function() {
  eps <- .Machine$double.eps
  list(
    link = "tanh",
    linkfun = function(mu) atanh(mu),
    linkinv = function(eta) pmin(pmax(tanh(eta), -1 + eps), 1 - eps),
    mu.eta = function(eta) {
      a <- exp(-2 * abs(eta))
      pmax(4 * a / (1 + a)^2, eps)
    },
    d2link = function(mu) 2 * mu / (1 - mu * mu)^2,
    d3link = function(mu) (2 + 6 * mu * mu) / (1 - mu * mu)^3,
    d4link = function(mu) 24 * mu * (1 + mu * mu) / (1 - mu * mu)^4,
    valideta = function(eta) TRUE
  )
}
