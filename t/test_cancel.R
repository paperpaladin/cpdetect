# How wrong the O(1) segment cost gets, measured.
#
# Every cost here is a difference of prefix sums, which is what makes the
# dynamic program O(n^2) instead of O(n^3) and is also the classic unstable
# way to compute a variance. The claim "it is fine in practice" is not
# checkable; the size of the error is. This file measures it against a
# two-pass computation over the segment itself, at four data offsets, with
# and without the two mitigations in R/cost.R.
#
# It also has to show the mitigations are load-bearing, so it asserts that
# the unmitigated version is measurably worse. A numerical safeguard whose
# removal changes nothing is not a safeguard.

source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "harness.R"))
cp_load()

set.seed(777001)
TOL <- 1e-9

# Two-pass reference. R's sum() accumulates doubles in a long double
# register, so this is the better-conditioned computation by a wide margin,
# which is the only reason it can serve as a reference for the other one.
ref_meanvar <- function(v, floorv) {
  s2 <- max(sum((v - mean(v)) ^ 2) / length(v), floorv)
  rss <- sum((v - mean(v)) ^ 2)
  length(v) * (log(2 * pi) + log(s2)) + rss / s2
}

worst_err <- function(y, ...) {
  co <- cp_cost(y, "meanvar", ...)
  n <- length(y)
  w <- 0
  for (s in 0:(n - 1)) for (e in (s + 1):n) {
    d <- abs(co$seg(s, e) - ref_meanvar(y[(s + 1):e], co$varfloor))
    if (is.finite(d) && d > w) w <- d
  }
  w
}

base <- rnorm(60)
offsets <- c(0, 1e3, 1e6, 1e9)
cent <- numeric(length(offsets))
raw <- numeric(length(offsets))
for (i in seq_along(offsets)) {
  y <- base + offsets[i]
  cent[i] <- worst_err(y, center = TRUE)
  raw[i] <- worst_err(y, center = FALSE)
}

for (i in seq_along(offsets))
  ok(cent[i] <= TOL, sprintf("centred_offset_%g", offsets[i]),
     sprintf("worst segment cost error %.3g", cent[i]))

# Vacuity guard. If the uncentred version were just as good, every assertion
# above would be true and would be saying nothing about centring.
ok(raw[3] > 1000 * max(cent[3], .Machine$double.xmin),
   "uncentred_is_measurably_worse_at_1e6",
   sprintf("uncentred %.3g vs centred %.3g, ratio %.3g",
           raw[3], cent[3], raw[3] / cent[3]))
ok(raw[4] > TOL, "uncentred_exceeds_the_suite_tolerance_at_1e9",
   sprintf("worst %.3g", raw[4]))
cat("#   offset      centred    uncentred\n")
for (i in seq_along(offsets))
  cat(sprintf("#   %-9g %-12.3g %.3g\n", offsets[i], cent[i], raw[i]))

# The variance floor on data whose scale is large. An absolute floor of
# 1e-8 sits far below the noise in the residual sum of squares there, and
# without a mitigation the cost of a near-constant segment becomes that
# noise divided by 1e-8. Series: a long run of a single large value with
# one short noisy stretch, so that some segments genuinely do hit the floor.
big <- c(rep(5e5, 25), 5e5 + rnorm(10, 0, 1e4), rep(5e5, 25))
dflt <- worst_err(big, center = TRUE)
bigger <- worst_err(big, center = TRUE, varfloor = 1e-3)
ok(dflt <= TOL, "floored_segments_agree_with_two_pass",
   sprintf("worst %.3g at the default floor of %g", dflt, CP_VARFLOOR))

# This pair is the measurement that decided which mitigation is doing the
# work. A floor proportional to the scale of the data was tried first and
# got this series from 2.0e+1 to 4.6e-2; the exact fallback in rss_of()
# gets it to 2.3e-13 with the plain absolute floor, and raising the floor
# by five orders of magnitude on top of that changes nothing. So the floor
# is not what makes the arithmetic right, and the relative floor was
# deleted after a separate sweep found it changed no segmentation either.
ok(bigger <= TOL, "the_accuracy_does_not_depend_on_the_size_of_the_floor",
   sprintf("floor 1e-3 gives %.3g against %.3g at the default", bigger, dflt))

# The exact cases. A one-observation segment has residual sum of squares
# exactly zero whatever the offset, and the code must return that rather
# than a cancellation residue, because the floor multiplies it by 1e8.
for (off in offsets) {
  co <- cp_cost(base + off, "meanvar")
  rs <- vapply(1:60, function(i) co$fit(i - 1, i)[["var"]], numeric(1))
  ok(all(rs == co$varfloor), sprintf("single_point_var_is_exactly_the_floor_at_%g", off),
     sprintf("distinct values %d", length(unique(rs))))
}

# Poisson is not centred, deliberately: there are no squares in its cost and
# shifting counts would make them non-integer. Check the claim rather than
# repeating it, with counts large enough that a cancellation would show.
yp <- rpois(50, 5e5)
cp <- cp_cost(yp, "poisson")
refp <- function(v) -2 * sum(dpois(v, mean(v), log = TRUE))
wp <- 0
for (s in 0:49) for (e in (s + 1):50) {
  d <- abs(cp$seg(s, e) - refp(yp[(s + 1):e]))
  if (d > wp) wp <- d
}
# This one cannot be fixed and the number is the point. -2 log L for a
# Poisson segment is 2 * sum(lgamma(y+1) + lambda - y log lambda), and at
# lambda = 5e5 each of those three terms is about 6e6 while their sum is
# about 13. Six orders of magnitude cancel, so an absolute error of
# eps * 6e6 = 1e-9 per observation is the best any implementation can do
# with lgamma, and the prefix sums over 50 observations carry about 1e-7 of
# it. It is not the prefix sums: computing the same thing one observation at
# a time has the same floor. Measured against the SIZE OF THE TERMS, which
# is the only comparison with any content.
scale_pois <- sum(abs(2 * lgamma(yp + 1)))
ok(wp / scale_pois < 1e-14, "poisson_large_counts_relative_to_term_size",
   sprintf("worst absolute %.3g; cost %.6g; terms sum to %.3g; relative %.3g",
           wp, cp$seg(0, 50), scale_pois, wp / scale_pois))

# At counts of the size this project is actually about -- failed logins per
# minute, process executions per interval -- the absolute error is at the
# suite tolerance, which is the claim that matters for telemetry.
yt <- rpois(50, 20)
ct <- cp_cost(yt, "poisson")
wt <- 0
for (s in 0:49) for (e in (s + 1):50) {
  d <- abs(ct$seg(s, e) - refp(yt[(s + 1):e]))
  if (d > wt) wt <- d
}
ok(wt <= TOL, "poisson_telemetry_scale_counts_are_exact_to_tolerance",
   sprintf("worst %.3g at lambda = 20", wt))

done()
