# Systematic prior-weights ("weighted likelihood") tests across all 12 *lss
# families. The contract: weighting a row by w == duplicating that row w times,
# so circ_gam(..., weights = w) fits the weighted MLE. See helper-weights.R for
# the registry (.WEIGHT_TEST_FAMILIES) and the expect_* helpers.

# ---- (1) EXACT ll-level identity: weighting == row duplication ---------------
# The core per-family gate. Machine-precision, no fitting, no finite differences,
# so it is unaffected by optimiser paths or stiff log-density regions.
for (nm in names(.WEIGHT_TEST_FAMILIES)) {
  local({
    fam_name <- nm
    nlp <- .WEIGHT_TEST_FAMILIES[[nm]]
    test_that(paste("weighted ll equals row duplication (exact):", fam_name), {
      dat <- .weight_test_data(60, seed = 42)
      G <- mgcv::gam(.weight_test_formula(nlp, smooth2 = TRUE),
                     family = get(fam_name)(), data = dat, fit = FALSE)
      set.seed(7)
      coef <- rnorm(ncol(G$X), sd = 0.12)
      set.seed(21)
      w <- sample(1:4, length(G$y), replace = TRUE)
      expect_weight_equals_duplication(G$family, G$X, G$y, coef, w)
    })
  })
}

# ---- (2) ret$l0 stays unweighted + NULL wt no-op ----------------------------
# l0 is the per-observation log-density; it must NEVER be scaled by wt (only the
# scalar objective l is). A NULL wt must behave like unit weights.
for (nm in names(.WEIGHT_TEST_FAMILIES)) {
  local({
    fam_name <- nm
    nlp <- .WEIGHT_TEST_FAMILIES[[nm]]
    test_that(paste("l0 unweighted and NULL wt is a no-op:", fam_name), {
      dat <- .weight_test_data(50, seed = 42)
      G <- mgcv::gam(.weight_test_formula(nlp, smooth2 = TRUE),
                     family = get(fam_name)(), data = dat, fit = FALSE)
      fam <- G$family; X <- G$X; y <- G$y; n <- length(y)
      set.seed(5)
      coef <- rnorm(ncol(X), sd = 0.12)
      set.seed(33)
      wt <- runif(n, 0.2, 3)

      r_w <- fam$ll(y, X, coef, wt = wt,         family = fam, deriv = 0)
      r_1 <- fam$ll(y, X, coef, wt = rep(1, n),  family = fam, deriv = 0)
      expect_equal(r_w$l0, r_1$l0)                 # l0 never scaled by wt
      expect_equal(r_w$l, sum(wt * r_w$l0))        # objective IS weighted
      expect_equal(r_1$l, sum(r_1$l0))

      # NULL wt == unit wt, at deriv = 1 (objective, gradient and Hessian)
      r_unit <- fam$ll(y, X, coef, wt = rep(1, n), family = fam, deriv = 1)
      r_null <- fam$ll(y, X, coef, wt = NULL,      family = fam, deriv = 1)
      expect_equal(r_unit$l, r_null$l)
      expect_equal(drop(r_unit$lb), drop(r_null$lb))
      expect_equal(unname(as.matrix(r_unit$lbb)), unname(as.matrix(r_null$lbb)))
    })
  })
}
