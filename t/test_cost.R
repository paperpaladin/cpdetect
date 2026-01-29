# The four segment costs against R's own stats.
#
# This is the weak oracle and it is here because it is cheap, not because it
# is convincing: dnorm and our prefix sums run in the same process, on the
# same libm, with the same rounding mode. It can catch an algebra mistake
# (a missing factor of two, the wrong variance estimator) and it structurally
# cannot catch a numerics mistake. t/test_scipy.R is the one that can.

source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "harness.R"))
cp_load()

set.seed(20261019)

# Tolerance. Both sides sum O(n) terms of magnitude O(1..10), so the error
# is O(n * 2^-52 * magnitude), about 1e-13 at n = 300. 1e-9 is four orders
# of margin and still three orders tighter than any defect we could make by
# algebra. The maximum observed difference is printed with every assertion
# so that a tolerance quietly doing all the work would be visible.
TOL <- 1e-9

# Reference: the same quantity from stats, computed per segment with no
# prefix sums anywhere.
ref_mean <- function(v, sigma2) -2 * sum(dnorm(v, mean(v), sqrt(sigma2), log = TRUE))
ref_var <- function(v, mu0, floorv) {
  s2 <- max(mean((v - mu0) ^ 2), floorv)
  -2 * sum(dnorm(v, mu0, sqrt(s2), log = TRUE))
}
ref_meanvar <- function(v, floorv) {
  s2 <- max(mean((v - mean(v)) ^ 2), floorv)
  -2 * sum(dnorm(v, mean(v), sqrt(s2), log = TRUE))
}
ref_pois <- function(v) -2 * sum(dpois(v, mean(v), log = TRUE))

# Segment list covering every length from 1 to n and a spread of offsets,
# not just a handful of hand-picked pairs. 1-element segments are included
# deliberately: they are where the variance costs degenerate.
pairs_for <- function(n) {
  out <- list()
  for (s in 0:(n - 1)) for (e in (s + 1):n) out[[length(out) + 1]] <- c(s, e)
  do.call(rbind, out)
}

check_all <- function(y, type, label, sigma2 = 1, mu = NULL) {
  co <- cp_cost(y, type, sigma2 = sigma2, mu = mu)
  P <- pairs_for(length(y))
  got <- cp_seg_many(co, P[, 1], P[, 2])
  want <- vapply(seq_len(nrow(P)), function(i) {
    v <- y[(P[i, 1] + 1):P[i, 2]]
    switch(type,
           mean = ref_mean(v, sigma2),
           var = ref_var(v, if (is.null(mu)) mean(y) else mu, co$varfloor),
           meanvar = ref_meanvar(v, co$varfloor),
           poisson = ref_pois(v))
  }, numeric(1))
  ok_near(got, want, TOL, paste0(label, "_", nrow(P), "_segments"))
}

y1 <- rnorm(40, 0, 1)
y2 <- c(rnorm(20, 0, 1), rnorm(20, 4, 1))
y3 <- c(rnorm(20, 0, 1), rnorm(20, 0, 3))
yc <- rep(2.5, 30)                               # the degenerate one
yp <- rpois(40, 3)
yz <- rep(0L, 25)                                # all-zero counts

check_all(y1, "mean", "mean_iid")
check_all(y2, "mean", "mean_jump")
check_all(y2, "mean", "mean_jump_sigma2_4", sigma2 = 4)
check_all(y3, "var", "var_scalechange")
check_all(y3, "var", "var_mu_given", mu = 0)
check_all(y1, "meanvar", "meanvar_iid")
check_all(y3, "meanvar", "meanvar_scalechange")
check_all(yc, "meanvar", "meanvar_constant_series")
check_all(yc, "var", "var_constant_series")
check_all(yp, "poisson", "poisson_counts")
check_all(yz, "poisson", "poisson_all_zero")

# The degenerate segment, called out on its own because it is the case the
# floor exists for and "it agrees with dnorm" is not the interesting claim.
cc <- cp_cost(yc, "meanvar")
ok(is.finite(cc$seg(0, 1)), "meanvar_single_point_is_finite")
ok(is.finite(cc$seg(0, 30)), "meanvar_constant_segment_is_finite")
ok(cc$fit(0, 30)[["var"]] == CP_VARFLOOR, "meanvar_constant_segment_hits_the_floor")

# Poisson, the 0 log 0 limit. An all-zero segment has log likelihood exactly
# zero, and the code reaches that by a branch rather than by arithmetic.
cz <- cp_cost(yz, "poisson")
ok_near(cz$seg(0, 25), 0, 0, "poisson_all_zero_segment_costs_exactly_zero")
ok(is.finite(cz$seg(0, 1)), "poisson_single_zero_is_finite")

# Additivity of total(): the number the optimisers minimise has to be the sum
# of the segment costs and nothing else.
cm <- cp_cost(y2, "mean")
ok_near(cp_total(cm, integer(0)), cm$seg(0, 40), 0, "total_no_changepoints_is_one_segment")
ok_near(cp_total(cm, 20), cm$seg(0, 20) + cm$seg(20, 40), 0, "total_one_changepoint_is_two_segments")
ok_near(cp_total(cm, c(10, 20, 30)),
        cm$seg(0, 10) + cm$seg(10, 20) + cm$seg(20, 30) + cm$seg(30, 40), 0,
        "total_three_changepoints_is_four_segments")

# Superadditivity: splitting a segment never costs more. This is PELT's
# precondition with K = 0, checked here on the cost functions themselves
# over every split of every segment, before any optimiser depends on it.
superadd <- function(y, type, label, ...) {
  co <- cp_cost(y, type, ...)
  n <- length(y)
  worst <- Inf
  for (s in 0:(n - 2)) for (e in (s + 2):n) for (t in (s + 1):(e - 1)) {
    d <- co$seg(s, e) - co$seg(s, t) - co$seg(t, e)
    if (d < worst) worst <- d
  }
  ok(worst >= -TOL, paste0("superadditive_", label),
     sprintf("min C(s,e) - C(s,t) - C(t,e) = %.3g over all splits", worst))
}
superadd(y2[1:18], "mean", "mean")
superadd(y3[1:18], "var", "var")
superadd(y3[1:18], "meanvar", "meanvar")
superadd(yc[1:18], "meanvar", "meanvar_constant")
superadd(yp[1:18], "poisson", "poisson")

# Rejections. A cost object must refuse what it cannot model rather than
# return a number for it.
refuses <- function(expr, label) {
  r <- tryCatch({ force(expr); "no error" }, error = function(e) conditionMessage(e))
  ok(r != "no error", paste0("refuses_", label), r)
}
refuses(cp_cost(numeric(0), "mean"), "empty_series")
refuses(cp_cost(c(1, NA, 3), "mean"), "na_in_series")
refuses(cp_cost(c(1, Inf), "mean"), "inf_in_series")
refuses(cp_cost(c(1, 2, 3), "nonsense"), "unknown_type")
refuses(cp_cost(c(1, 2, 3), "mean", sigma2 = 0), "zero_sigma2")
refuses(cp_cost(c(1, 2, 3), "mean", sigma2 = -1), "negative_sigma2")
refuses(cp_cost(c(1, -2, 3), "poisson"), "negative_counts")
refuses(cp_cost(c(1, 2.5, 3), "poisson"), "non_integer_counts")

done()
