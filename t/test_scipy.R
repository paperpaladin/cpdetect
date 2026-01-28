# The cross-language check.
#
# t/oracle.py computes the same four costs from scipy.stats.norm.logpdf and
# scipy.stats.poisson.logpmf, two passes over each segment, with math.fsum
# where the summation order could matter. It shares no code with R/cost.R and
# no arithmetic: ours is a difference of prefix sums computed once over the
# whole series, its is a fresh pass over the segment. A disagreement here is
# a numerics defect that no amount of R-against-R testing can reach, which is
# the entire reason this file exists and t/test_cost.R is not enough.

source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "harness.R"))
cp_load()

if (!oracle_available()) {
  # Loud, and red. scipy is oracle 2 of 6 in PLAN.md and a run without it has
  # checked something weaker than it claims to have checked. A quiet skip
  # here would print a smaller number and look like a shorter suite.
  ok(FALSE, "scipy_oracle_present",
     "python3 with scipy and numpy is required; see PLAN.md Oracle 2")
  done()
}
ok(TRUE, "scipy_oracle_present")

set.seed(31337)

# Tolerance, and why. Both sides sum O(n) terms of magnitude O(1..20), so
# the representable floor on the disagreement is about n * eps * magnitude,
# which is 1e-13 at the sizes used here. 1e-9 leaves four orders of margin.
# It is the same tolerance as t/test_cost.R deliberately: if the two oracles
# needed different tolerances for the same quantity, one of them would be
# measuring something else. Every assertion prints the observed maximum so
# that a tolerance silently doing the work would be visible in the log.
TOL <- 1e-9

all_pairs <- function(n) {
  s <- unlist(lapply(0:(n - 1), function(i) rep(i, n - i)))
  e <- unlist(lapply(0:(n - 1), function(i) (i + 1):n))
  list(s = s, e = e)
}

cross <- function(y, type, label, ...) {
  co <- cp_cost(y, type, ...)
  P <- all_pairs(length(y))
  got <- cp_seg_many(co, P$s, P$e)
  want <- as.numeric(oracle("cost", oracle_cost_req(co, P$s, P$e), tag = label))
  if (length(want) != length(got)) {
    ok(FALSE, paste0("scipy_", label), sprintf("oracle returned %d of %d", length(want), length(got)))
    return(invisible(NULL))
  }
  ok_near(got, want, TOL, paste0("scipy_", label, "_", length(got), "_segments"))
}

y1 <- rnorm(30)
y2 <- c(rnorm(15, 0, 1), rnorm(15, 3, 1))
y3 <- c(rnorm(15, 0, 1), rnorm(15, 0, 4))
y4 <- rnorm(30) + 1e6                      # the offset that broke the prefix sums
y5 <- rep(7, 20)                           # degenerate, the floor
yp <- rpois(30, 4)
yz <- rep(0L, 20)
yb <- rpois(30, 400)

cross(y1, "mean", "mean_iid")
cross(y2, "mean", "mean_jump")
cross(y2, "mean", "mean_sigma2_9", sigma2 = 9)
cross(y4, "mean", "mean_offset_1e6")
cross(y3, "var", "var_mu_zero", mu = 0)
cross(y3, "var", "var_mu_default")
cross(y1, "meanvar", "meanvar_iid")
cross(y3, "meanvar", "meanvar_scalechange")
cross(y4, "meanvar", "meanvar_offset_1e6")
cross(y5, "meanvar", "meanvar_constant")
cross(y5, "var", "var_constant")
cross(yp, "poisson", "poisson_small")
cross(yb, "poisson", "poisson_lambda_400")
cross(yz, "poisson", "poisson_all_zero")

# The oracle has to be capable of disagreeing. If oracle.py echoed our own
# numbers back, or if the request never reached it, every line above would
# pass. Send it a request whose answer we know is different and require the
# comparison to notice: the mean cost at sigma2 = 1 against sigma2 = 4 on
# the same data differ by a known amount, n * log(4).
co <- cp_cost(y1, "mean", sigma2 = 1)
P <- all_pairs(10)
req <- oracle_cost_req(co, P$s[1:5], P$e[1:5])
req[2] <- paste("sigma2", g17(4))
wrong <- as.numeric(oracle("cost", req, tag = "vacuity"))
right <- cp_seg_many(co, P$s[1:5], P$e[1:5])
ok(!near(wrong, right, TOL), "oracle_can_disagree",
   sprintf("maxdiff %.4g with sigma2 changed under it", max(abs(wrong - right))))
# ... and that the disagreement is the one the change implies, so that a
# broken oracle returning rubbish is not mistaken for a working one.
# -2 log L at sigma2 = 4 is m log(2 pi 4) + rss/4, and ours at sigma2 = 1 is
# m log(2 pi) + rss, so rss can be recovered from one and the other
# predicted exactly. An oracle returning noise fails this and passes the
# assertion above it.
m <- P$e[1:5] - P$s[1:5]
rss <- right - m * log(2 * pi)
ok_near(wrong, m * log(2 * pi * 4) + rss / 4, TOL, "oracle_disagrees_by_the_right_amount")

done()
