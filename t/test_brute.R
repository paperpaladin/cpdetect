# The hard invariant: the dynamic program finds the exhaustive minimum.
#
# For a short enough series every admissible segmentation can be scored and
# the best taken, and that is not a second opinion about the answer, it is
# the answer. cp_op has to match it on every case. The sweep covers all four
# costs, four penalties, three minimum segment lengths and a spread of true
# change point counts, because what this catches is off-by-one errors in the
# admissibility bounds, and those live at the edges of the grid.
#
# Two regimes, and the first run of this file taught me the difference.
#
#   Unrestricted. Every admissible segmentation, no cap on the number of
#   change points. Feasible up to n = 13 at minseg 1 (2^12), and up to
#   n = 25 at minseg 2 or 3, where requiring segments of length 2 cuts the
#   count from 16777216 to 46368 -- the compositions of n into parts of
#   size at least 2 are Fibonacci, not binary.
#
#   Capped at three change points, which is the shape the brief names for
#   n = 25. A capped enumeration is only comparable with an uncapped search
#   if the uncapped optimum is inside the cap, and on the first run five of
#   sixteen capped cases failed because it was not: with penalty "none" the
#   optimum of a 22-point series has eight change points, so brute force
#   reported a minimum five units ABOVE the dynamic program's. That is the
#   comparison being invalid, not cp_op being wrong, and it is asserted
#   rather than assumed from here on.

source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "harness.R"))
cp_load()

# Tolerance. Both totals are sums of at most 25 segment costs of magnitude
# up to about 100, so the representable disagreement is of order
# 25 * 2^-52 * 100 = 5e-13. The dynamic program accumulates through F and
# brute force sums the segments and adds k*beta at the end, so the grouping
# of the betas differs and the last bits with it. 1e-9 is three orders of
# margin; the worst observed value is reported with the assertion.
TOL <- 1e-9

types <- c("mean", "var", "meanvar", "poisson")

# Segmentations of n into exactly j segments each at least ms long are the
# compositions of n into j parts of size >= ms: choose(n - j*ms + j - 1, j - 1)
# of them. Summed over j up to the cap, this is how many cp_brute must have
# scored. Checked per case, because an enumeration that quietly skips most
# of the space agrees with everything.
n_expected <- function(n, ms, maxparts) {
  js <- seq_len(min(maxparts, n %/% ms))
  sum(choose(n - js * ms + js - 1, js - 1))
}

cases <- list()
add <- function(...) cases[[length(cases) + 1L]] <<- list(...)
sd_ct <- 0L

# Short and unrestricted: the full grid.
for (ty in types) for (pn in c("none", "aic", "bic", "mbic")) for (ms in 1:3) {
  sd_ct <- sd_ct + 1L
  add(type = ty, pen = pn, minseg = ms,
      n = c(10L, 12L, 13L)[1 + (sd_ct %% 3L)], k = sd_ct %% 3L,
      maxcp = Inf, seed = 410000L + sd_ct)
}
# Long and unrestricted, which minseg >= 2 makes possible.
for (ty in types) for (cfg in list(c(20L, 2L), c(22L, 2L), c(24L, 3L), c(25L, 3L))) {
  sd_ct <- sd_ct + 1L
  add(type = ty, pen = if (sd_ct %% 2L) "bic" else "mbic", minseg = cfg[2],
      n = cfg[1], k = 1L + (sd_ct %% 3L), maxcp = Inf, seed = 420000L + sd_ct)
}
# Long and capped at three, at minseg 1. Not meanvar: at minseg 1 the
# variance costs are degenerate and the optimum is one segment per
# observation, which is a real property of the cost and is measured in
# t/test_degenerate.R rather than worked around here.
for (ty in c("mean", "var", "poisson")) {
  sd_ct <- sd_ct + 1L
  add(type = ty, pen = "bic", minseg = 1L, n = 25L, k = 2L,
      maxcp = 3, seed = 430000L + sd_ct)
}

worst_cost <- 0
worst_self <- 0
n_mismatch <- 0L
mismatch_detail <- character(0)
total_enumerated <- 0
n_cap_binding <- 0L
n_count_wrong <- 0L
count_detail <- character(0)
n_cases <- 0L

for (cs in cases) {
  sim <- cp_sim_random(cs$n, cs$k, cs$type, effect = 3, gap = cs$minseg,
                       seed = cs$seed)
  co <- cp_cost(sim$y, cs$type)
  pr <- cp_prepare(co, cs$pen)
  a <- cp_op(pr$cost, pr$beta, cs$minseg)
  b <- cp_brute(pr$cost, pr$beta, cs$minseg, maxcp = cs$maxcp)
  n_cases <- n_cases + 1L
  total_enumerated <- total_enumerated + b$nconsidered

  if (is.finite(cs$maxcp) && length(b$cp) >= cs$maxcp) n_cap_binding <- n_cap_binding + 1L
  want_n <- n_expected(cs$n, cs$minseg,
                       if (is.finite(cs$maxcp)) cs$maxcp + 1 else cs$n)
  if (b$nconsidered != want_n) {
    n_count_wrong <- n_count_wrong + 1L
    count_detail <- c(count_detail,
      sprintf("%s/%s/ms%d/n%d: enumerated %d, closed form %d",
              cs$type, cs$pen, cs$minseg, cs$n, b$nconsidered, want_n))
  }

  d <- abs(a$cost - b$cost)
  if (d > worst_cost) worst_cost <- d
  if (!identical(as.integer(a$cp), as.integer(b$cp))) {
    n_mismatch <- n_mismatch + 1L
    mismatch_detail <- c(mismatch_detail,
      sprintf("%s/%s/ms%d/n%d: op [%s] %.10f vs brute [%s] %.10f",
              cs$type, cs$pen, cs$minseg, cs$n,
              paste(a$cp, collapse = ","), a$cost,
              paste(b$cp, collapse = ","), b$cost))
  }
  # Self-consistency: the change points cp_op reports must actually produce
  # the cost it reports. A backtrack that drifts from the F table is a
  # defect the cost comparison on its own cannot see.
  ds <- abs(cp_total(pr$cost, a$cp) + length(a$cp) * pr$beta - a$cost)
  if (ds > worst_self) worst_self <- ds
}

ok(n_cases == length(cases), "all_cases_ran", sprintf("%d cases", n_cases))
ok(worst_cost <= TOL, "op_cost_equals_brute_force_minimum",
   sprintf("worst |op - brute| = %.3g over %d cases, tol %.3g",
           worst_cost, n_cases, TOL))
ok(n_mismatch == 0L, "op_segmentation_equals_brute_force_segmentation",
   if (n_mismatch) paste(mismatch_detail, collapse = " | ")
   else sprintf("%d cases, identical change point vectors", n_cases))
ok(worst_self <= TOL, "op_reported_cost_matches_its_own_segmentation",
   sprintf("worst %.3g", worst_self))
ok(n_cap_binding == 0L, "the_three_changepoint_cap_never_binds",
   sprintf("%d of %d capped cases used the whole cap", n_cap_binding,
           sum(vapply(cases, function(c) is.finite(c$maxcp), logical(1)))))
ok(n_count_wrong == 0L, "enumeration_count_matches_the_closed_form",
   if (n_count_wrong) paste(count_detail, collapse = " | ")
   else sprintf("%d cases, %.0f segmentations scored in total",
                n_cases, total_enumerated))

# And the enumeration has to be capable of beating a wrong search. A
# deliberately crippled one -- considering only the no-change-point
# segmentation -- must lose to it by a wide margin on a series with an
# obvious change.
sim <- cp_sim(40, c(20), "mean", mean = c(0, 5), seed = 99001)
co <- cp_cost(sim$y, "mean"); pr <- cp_prepare(co, "bic")
crippled <- cp_total(pr$cost, integer(0))
b <- cp_brute(pr$cost, pr$beta, 1L, maxcp = 2)
ok(b$cost < crippled - 1,
   "brute_force_beats_a_segmentation_that_misses_a_real_change",
   sprintf("brute %.4f vs no-change %.4f", b$cost, crippled))

# minseg is a constraint brute force applies by filtering and cp_op by
# bounding its loop: two mechanisms that have to agree about which
# segmentations exist. Checked on the segment lengths directly as well.
for (ms in c(2L, 3L, 4L)) {
  sim <- cp_sim_random(24L, 2L, "meanvar", effect = 4, gap = ms, seed = 7700L + ms)
  co <- cp_cost(sim$y, "meanvar"); pr <- cp_prepare(co, "bic")
  a <- cp_op(pr$cost, pr$beta, ms)
  lens <- diff(c(0L, a$cp, 24L))
  ok(all(lens >= ms), sprintf("op_respects_minseg_%d", ms),
     sprintf("segment lengths %s", paste(lens, collapse = ",")))
}

# Degenerate shapes, where an off-by-one becomes a crash or an empty answer
# rather than a wrong number.
one <- cp_cost(c(1.5), "mean")
ok(length(cp_op(one, 1, 1L)$cp) == 0L, "op_single_observation_has_no_changepoint")
two <- cp_cost(c(1.5, 9.5), "mean")
ok(length(cp_op(two, 0, 2L)$cp) == 0L, "op_minseg_forbids_the_only_split")
ok(length(cp_op(two, 0, 1L)$cp) == 1L, "op_minseg_1_takes_the_split_at_zero_penalty")
ok(length(cp_op(cp_cost(rnorm(10), "mean"), 0, 20L)$cp) == 0L,
   "op_minseg_larger_than_n_gives_one_segment")

done()
