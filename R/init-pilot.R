## Projected mean-direction warm start for a single circular (tan-half) location
## linear predictor: penalized-smooth cos(y) and sin(y) over the location LP's
## design columns (the penalty root E stacked as a regularizer), recombine the
## fitted mean direction mu-hat = atan2(s, c), and start the LP from link(mu-hat)
## (= tan(mu-hat/2) under the tan-half link).
##
## A flat (constant) mu start strands antipodal observations on the log-density
## cliffs of any density with a zero on the circle -- Cartwright's
## (1 + cos)^(1/zeta) is exactly 0 at the antipode -- and, even where there is
## no zero, is start-sensitive on near-flat shape surfaces: the optimum a
## first-order optimizer reaches there depends on where it began. The
## data-following pilot sidesteps the cliffs. Used by every tan-half location
## family; the projected normal's two Cartesian location components have no
## antipode zero, so pnlss keeps the constant start and does not call this.
.tanhalf_pilot_coef <- function(Xmu, Emu, y, off, linkfun) {
  pls <- function(yt) {
    b <- qr.coef(qr(rbind(Xmu, Emu)), c(yt, rep(0, nrow(Emu))))
    b[!is.finite(b)] <- 0
    b
  }
  shat <- drop(Xmu %*% pls(sin(y)))
  chat <- drop(Xmu %*% pls(cos(y)))
  eta_mu <- pmax(pmin(linkfun(atan2(shat, chat)), 1e6), -1e6) - off
  eta_mu[!is.finite(eta_mu)] <- 0
  pls(eta_mu)
}
