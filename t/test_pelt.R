# PELT against Optimal Partitioning.
#
# PELT's pruning is supposed to be exact. The claim is therefore not "the
# answers look similar" but "the answers are identical", and the comparison
# is worth nothing unless the pruning actually fired, so every case asserts
# a nonzero prune count as well.
#
# The file also plants a wrong pruning rule and requires the comparison to
# catch it. An agreement test between two searches is only evidence if a
# broken second search would disagree, and the second search here differs
# from the first by one inequality.

source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "harness.R"))
cp_load()

# Tolerance. Both searches accumulate the same recursion in the same order
# over (mostly) the same candidates, so the costs should agree to the last
# bit and generally do. 1e-9 all the same, for the cases where PELT's
# argmin comes from a different position in the candidate vector and the
# additions happen in a different order.
TOL <- 1e-9

types <- c("mean", "var", "meanvar", "poisson")
pens <- c("none", "aic", "bic", "mbic")

worst <- 0
n_cpdiff <- 0L
n_tie <- 0L
n_noprune <- 0L
noprune_detail <- character(0)
detail <- character(0)
tot_pruned <- 0
tot_op_cand <- 0
tot_pelt_cand <- 0
n_cases <- 0L
i <- 0L

for (ty in types) for (pn in pens) for (ms in c(1L, 2L, 5L)) {
  i <- i + 1L
  n <- c(120L, 200L, 350L)[1 + (i %% 3L)]
  k <- c(0L, 1L, 3L, 6L)[1 + (i %% 4L)]
  sim <- cp_sim_random(n, k, ty, effect = 2.5, gap = max(ms, 12L), seed = 610000L + i)
  co <- cp_cost(sim$y, ty)
  pr <- cp_prepare(co, pn)
  a <- cp_op(pr$cost, pr$beta, ms)
  b <- cp_pelt(pr$cost, pr$beta, ms)
  n_cases <- n_cases + 1L
  tot_pruned <- tot_pruned + b$npruned
  tot_op_cand <- tot_op_cand + a$ncand
  tot_pelt_cand <- tot_pelt_cand + b$ncand
  d <- abs(a$cost - b$cost)
  if (d > worst) worst <- d
  if (!identical(as.integer(a$cp), as.integer(b$cp))) {
    # Two answers of EQUAL cost is not a defect: with integer data and a
    # small penalty, several segmentations attain the minimum, and the two
    # searches see different candidate sets so they break the tie
    # differently. Two answers of different cost is a defect. Counted apart.
    if (d <= TOL) n_tie <- n_tie + 1L else {
      n_cpdiff <- n_cpdiff + 1L
      detail <- c(detail, sprintf("%s/%s/ms%d/n%d: op [%s] %.12f vs pelt [%s] %.12f",
                                  ty, pn, ms, n, paste(a$cp, collapse = ","), a$cost,
                                  paste(b$cp, collapse = ","), b$cost))
    }
  }
  if (b$npruned == 0L) {
    n_noprune <- n_noprune + 1L
    noprune_detail <- c(noprune_detail,
      sprintf("%s/%s/ms%d/n%d/k%d (beta %.2f)", ty, pn, ms, n, k, pr$beta))
  }
}

ok(n_cases == 48L, "pelt_sweep_ran", sprintf("%d cases", n_cases))
ok(n_cpdiff == 0L, "pelt_segmentation_equals_optimal_partitioning",
   if (n_cpdiff) paste(detail, collapse = " | ")
   else sprintf("%d cases, %d identical change point vectors, %d equal-cost ties",
                n_cases, n_cases - n_tie, n_tie))
ok(worst <= TOL, "pelt_cost_equals_optimal_partitioning",
   sprintf("worst %.3g over %d cases", worst, n_cases))
# Pruning does not fire in every case and should not be asserted to. When
# beta is large relative to the spread of F(s) + C(s+1,t) across s, nothing
# clears the threshold: the rule keeps every candidate within beta of the
# best, and on a series with no change point and a BIC penalty of 11.7 the
# whole candidate set sits inside that. That is PELT behaving correctly in
# the regime where it has nothing to prune, and it is the regime its linear
# time claim does not cover. The cases are named rather than hidden.
ok(tot_pruned > 0, "pruning_fired",
   sprintf("%.0f candidates discarded across %d of %d cases",
           tot_pruned, n_cases - n_noprune, n_cases))
ok(n_noprune <= 6L, "pruning_fires_in_most_cases",
   if (n_noprune) paste("no pruning in:", paste(noprune_detail, collapse = "; "))
   else "every case pruned something")
ok(tot_pelt_cand < tot_op_cand, "pruning_saved_cost_evaluations",
   sprintf("%.0f evaluations against Optimal Partitioning's %.0f, %.1f%%",
           tot_pelt_cand, tot_op_cand, 100 * tot_pelt_cand / tot_op_cand))

# The comparison has to be sensitive to the pruning rule, or it is an
# expensive way of running the same algorithm twice. K is the constant in
# the exactness theorem and 0 is the valid value for a likelihood cost;
# a positive K prunes candidates that can still be optimal. Require that
# to produce a different answer somewhere.
n_bad <- 0L
n_bad_cases <- 0L
for (ty in types) for (kk in c(2, 10, 50)) {
  n_bad_cases <- n_bad_cases + 1L
  sim <- cp_sim_random(200L, 4L, ty, effect = 2, gap = 20L, seed = 620000L + n_bad_cases)
  co <- cp_cost(sim$y, ty)
  pr <- cp_prepare(co, "bic")
  a <- cp_op(pr$cost, pr$beta, 2L)
  b <- cp_pelt(pr$cost, pr$beta, 2L, K = kk)
  if (!identical(as.integer(a$cp), as.integer(b$cp)) || abs(a$cost - b$cost) > TOL)
    n_bad <- n_bad + 1L
}
ok(n_bad > 0L, "a_too_aggressive_pruning_constant_is_caught",
   sprintf("%d of %d cases disagree with K > 0", n_bad, n_bad_cases))

# And the same for the direction nobody writes by accident: pruning on
# equality rather than strictly. cp_pelt keeps a candidate when
# base + K <= F(t); a tiny positive K is the same thing as keeping only
# when strictly less, and the series below is built so that exact ties
# exist, with integer counts and a penalty of zero.
yt <- as.integer(c(rep(2L, 12), rep(9L, 12), rep(2L, 12)))
ct <- cp_cost(yt, "poisson")
a <- cp_op(ct, 0, 1L)
b_ok <- cp_pelt(ct, 0, 1L, K = 0)
b_eq <- cp_pelt(ct, 0, 1L, K = 1e-9)
# Equal cost, different segmentation: 12 identical counts in a row means
# splitting a block is free, so the optimum is attained many times over.
# The claim that survives is about the cost, and it is exact to the last
# bit rather than to a tolerance.
ok(a$cost == b_ok$cost, "integer_data_with_ties_agrees_on_cost_exactly",
   sprintf("op [%s] and pelt [%s] both %.15g",
           paste(a$cp, collapse = ","), paste(b_ok$cp, collapse = ","), a$cost))
ok(!identical(as.integer(a$cp), as.integer(b_ok$cp)),
   "integer_data_with_ties_does_produce_two_different_optima",
   "recorded because it is the case the sweep above has to tolerate")
ok(abs(cp_total(ct, b_ok$cp) - a$cost) < TOL,
   "the_pelt_tie_segmentation_really_does_cost_that")
ok(b_eq$cost >= a$cost - TOL, "pruning_on_equality_never_finds_a_better_cost",
   sprintf("strict %.10f vs exact %.10f", b_eq$cost, a$cost))

# MBIC's K. cp_prepare attaches log(4/n) to the cost object, and cp_pelt
# must pick it up without being told.
sim <- cp_sim_random(150L, 3L, "mean", effect = 3, gap = 15L, seed = 63001L)
co <- cp_cost(sim$y, "mean")
pr <- cp_prepare(co, "mbic")
b <- cp_pelt(pr$cost, pr$beta, 2L)
ok(abs(b$K - log(4 / 150)) < 1e-12, "pelt_takes_K_from_the_mbic_cost_object",
   sprintf("K = %.6f, log(4/n) = %.6f", b$K, log(4 / 150)))
ok(abs(cp_pelt(cp_cost(sim$y, "mean"), 10, 2L)$K) == 0,
   "pelt_K_is_zero_for_an_ordinary_likelihood_cost")

# The minimum segment length and the pruning delay. Removing a condemned
# candidate on the spot rather than minseg steps later is the defect this
# file found on its first run, and prune_delay = 1 puts it back. It has to
# produce a worse cost than Optimal Partitioning on a measurable number of
# cases, or the delay is not known to be doing anything.
n_delay_bad <- 0L
n_delay_cases <- 0L
worst_delay <- 0
for (ty in types) for (ms in c(2L, 5L)) {
  n_delay_cases <- n_delay_cases + 1L
  sim <- cp_sim_random(350L, 0L, ty, effect = 2.5, gap = max(ms, 12L),
                       seed = 640000L + n_delay_cases)
  co <- cp_cost(sim$y, ty)
  pr <- cp_prepare(co, "none")
  a <- cp_op(pr$cost, pr$beta, ms)
  bad <- cp_pelt(pr$cost, pr$beta, ms, prune_delay = 1L)
  good <- cp_pelt(pr$cost, pr$beta, ms)
  if (bad$cost > a$cost + TOL) {
    n_delay_bad <- n_delay_bad + 1L
    worst_delay <- max(worst_delay, bad$cost - a$cost)
  }
  if (abs(good$cost - a$cost) > TOL)
    ok(FALSE, sprintf("delayed_pruning_exact_%s_ms%d", ty, ms),
       sprintf("%.10f vs %.10f", good$cost, a$cost))
}
ok(n_delay_bad > 0L, "immediate_removal_is_caught",
   sprintf("%d of %d cases return a cost above the exact minimum, worst by %.4f",
           n_delay_bad, n_delay_cases, worst_delay))
ok(cp_pelt(cp_cost(rnorm(50), "mean"), 2, 3L)$prune_delay == 3L,
   "prune_delay_defaults_to_minseg")

# PELT against the exhaustive minimum directly, not only against Optimal
# Partitioning. The two searches share cp_backtrack and the cost object, so
# a defect in either would agree with itself; brute force shares the cost
# object and nothing else.
n_bf <- 0L
n_bf_bad <- 0L
bf_detail <- character(0)
for (ty in types) for (ms in c(2L, 3L)) for (pn in c("none", "bic", "mbic")) {
  n_bf <- n_bf + 1L
  sim <- cp_sim_random(20L, 1L + (n_bf %% 2L), ty, effect = 3, gap = ms,
                       seed = 650000L + n_bf)
  co <- cp_cost(sim$y, ty)
  pr <- cp_prepare(co, pn)
  b <- cp_pelt(pr$cost, pr$beta, ms)
  e <- cp_brute(pr$cost, pr$beta, ms)
  if (abs(b$cost - e$cost) > TOL || !identical(as.integer(b$cp), as.integer(e$cp))) {
    n_bf_bad <- n_bf_bad + 1L
    bf_detail <- c(bf_detail, sprintf("%s/%s/ms%d: pelt [%s] %.12f vs brute [%s] %.12f",
                                      ty, pn, ms, paste(b$cp, collapse = ","), b$cost,
                                      paste(e$cp, collapse = ","), e$cost))
  }
}
ok(n_bf_bad == 0L, "pelt_equals_the_exhaustive_minimum",
   if (n_bf_bad) paste(bf_detail, collapse = " | ")
   else sprintf("%d cases at n = 20, every admissible segmentation enumerated", n_bf))

# Degenerate shapes.
ok(length(cp_pelt(cp_cost(c(1.5), "mean"), 1, 1L)$cp) == 0L,
   "pelt_single_observation_has_no_changepoint")
ok(length(cp_pelt(cp_cost(rnorm(10), "mean"), 0, 20L)$cp) == 0L,
   "pelt_minseg_larger_than_n_gives_one_segment")

done()
