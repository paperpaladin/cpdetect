# MBIC, and the pruning constant it costs.
#
# Zhang and Siegmund's modified BIC is the one penalty in R/penalty.R that
# is not a constant added per change point. Its criterion is
#
#     -2 log L + 3k log n + sum_{i=1}^{k+1} log(l_i / n)
#
# and the sum over segment lengths is what makes it worth having: a change
# point that creates a very short segment is penalised more than one that
# does not. It is also what breaks PELT, because the search needs the cost
# to be superadditive and log(a) + log(b) against log(a+b) is not.
#
# Three separate things are checked here. That our rearrangement of the
# criterion into "a cost plus a constant per change point" is the same
# function. That the superadditivity constant K is where the algebra says
# it is, by exhaustive search over triples rather than by argument. And
# that running PELT at K = 0 under MBIC is not merely unproven but actually
# wrong, measured as a rate over seeds.

source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "harness.R"))
cp_load()

TOL <- 1e-9

# ------------------------------------------------- 1. the rearrangement
#
# cp_mbic_cost adds log(segment length) to each segment and cp_penalty
# returns beta = 2 log n. The claim is that
#
#     sum (C_i + log l_i) + 2k log n
#
# equals the criterion above, up to a constant that no comparison between
# segmentations of the same series can see. Checked by computing both from
# the definition on arbitrary segmentations, not only on optima.
#
# The constant is +log n:
#
#   theirs = sum C_i + 3k log n + sum log l_i - (k+1) log n
#   ours   = sum C_i + sum log l_i + 2k log n
#   ours - theirs = (2k - 3k + k + 1) log n = log n
#
# and the first version of this file asserted -log n, which is the sign
# error you get by reading "sum log(l_i / n)" as contributing -(k+1) log n
# to OUR side instead of to theirs. The assertion caught it: the spread
# over 40 segmentations was 1.1e-13, so the offset was constant and
# therefore invisible to the optimiser and harmless, and it was 8.19 away
# from the value claimed. A constant difference that is the wrong constant
# is exactly the defect that only a direct comparison can find, because
# every segmentation this project ranks would come out in the same order.
set.seed(880011)
n <- 60L
sim <- cp_sim_random(n, 3L, "meanvar", effect = 3, gap = 6L, seed = 880012L)
plain <- cp_cost(sim$y, "meanvar")
pr <- cp_prepare(plain, "mbic")
diffs <- numeric(0)
for (rep in 1:40) {
  k <- sample(0:5, 1)
  taus <- if (k == 0L) integer(0) else sort(sample.int(n - 1L, k))
  if (any(diff(c(0L, taus, n)) < 1L)) next
  lens <- diff(c(0L, taus, n))
  ours <- cp_total(pr$cost, taus) + length(taus) * pr$beta
  theirs <- cp_total(plain, taus) + 3 * length(taus) * log(n) + sum(log(lens / n))
  diffs <- c(diffs, ours - theirs)
}
ok(length(diffs) >= 30, "mbic_rearrangement_cases",
   sprintf("%d segmentations compared", length(diffs)))
ok(max(abs(diffs - log(n))) <= TOL,
   "mbic_rearrangement_differs_from_the_criterion_by_exactly_log_n",
   sprintf("worst |diff - log n| = %.3g, and log n = %.6f",
           max(abs(diffs - log(n))), log(n)))
# A constant is only invisible if it really is constant. Two segmentations
# with different k must differ by the same amount, which the max above
# already covers, but the spread is the readable form of it.
ok(diff(range(diffs)) <= TOL, "the_offset_does_not_depend_on_k",
   sprintf("range over %d segmentations %.3g", length(diffs), diff(range(diffs))))

# ----------------------------------- 2. the superadditivity constant K
#
# For every s < t < T, how far short does splitting fall?
#
#     gap(s,t,T) = C(s,T) - C(s,t) - C(t,T)
#
# K must be at most the minimum of that over all triples. For a plain
# likelihood cost the minimum is >= 0 and K = 0 is valid. For the MBIC cost
# the log terms make it negative, and the algebra says the bound is
# log(4/n), attained when the two halves are equal and span the series.
min_gap <- function(cost, nn) {
  w <- Inf
  for (s in 0:(nn - 2)) for (e in (s + 2):nn) for (t in (s + 1):(e - 1)) {
    g <- cost$seg(s, e) - cost$seg(s, t) - cost$seg(t, e)
    if (g < w) w <- g
  }
  w
}
nn <- 24L
for (ty in c("mean", "var", "meanvar", "poisson")) {
  s2 <- cp_sim_random(nn, 1L, ty, effect = 3, gap = 4L, seed = 881000L + nchar(ty))
  base <- cp_cost(s2$y, ty)
  g0 <- min_gap(base, nn)
  ok(g0 >= -TOL, paste0("K_zero_is_valid_for_", ty),
     sprintf("min gap over %d triples = %.3g", choose(nn, 2) * 0 + sum(choose(2:nn, 2)), g0))
  gm <- min_gap(cp_mbic_cost(base), nn)
  bound <- log(4 / nn)
  ok(gm >= bound - TOL, paste0("K_log4_over_n_is_valid_for_mbic_", ty),
     sprintf("min gap %.6f, bound log(4/n) = %.6f, slack %.3g", gm, bound, gm - bound))
  ok(gm < -TOL, paste0("K_zero_is_NOT_valid_for_mbic_", ty),
     sprintf("min gap %.6f is below zero, so the usual K = 0 argument does not hold", gm))
}
# The bound is claimed to be attained in the limit by an equal split of the
# whole series. For the pure log part, with no data in it, it is exact.
logpart <- list(n = nn, seg = function(s, e) log(e - s))
ok(abs(min_gap(logpart, nn) - log(4 / nn)) <= TOL,
   "the_log_term_alone_attains_log_4_over_n",
   sprintf("min gap %.10f against %.10f", min_gap(logpart, nn), log(4 / nn)))

# ------------------------------- 3. PELT at K = 0 under MBIC, measured
#
# Not proven exact is a weaker statement than wrong. This measures which.
n_wrong <- 0L
n_right_bad <- 0L
n_runs <- 0L
worst_excess <- 0
right_detail <- character(0)
for (ty in c("mean", "var", "meanvar", "poisson")) for (sd in 1:12) {
  n_runs <- n_runs + 1L
  s3 <- cp_sim_random(140L, 1L + (sd %% 4L), ty, effect = 2, gap = 8L,
                      seed = 882000L + n_runs)
  co <- cp_cost(s3$y, ty)
  p3 <- cp_prepare(co, "mbic")
  exact <- cp_op(p3$cost, p3$beta, 3L)
  right <- cp_pelt(p3$cost, p3$beta, 3L)                  # K = log(4/n)
  wrong <- cp_pelt(p3$cost, p3$beta, 3L, K = 0)
  if (abs(right$cost - exact$cost) > TOL) {
    n_right_bad <- n_right_bad + 1L
    right_detail <- c(right_detail,
      sprintf("%s seed %d: %.10f vs exact %.10f", ty, sd, right$cost, exact$cost))
  }
  if (wrong$cost > exact$cost + TOL) {
    n_wrong <- n_wrong + 1L
    worst_excess <- max(worst_excess, wrong$cost - exact$cost)
  }
}
ok(n_runs == 48L, "mbic_pelt_runs", sprintf("%d runs", n_runs))
ok(n_right_bad == 0L, "mbic_pelt_at_K_log4_over_n_is_exact",
   if (n_right_bad) paste(right_detail, collapse = " | ")
   else sprintf("%d runs, every one equal to the unpruned answer", n_runs))
# "Not proven exact" and "wrong" are different claims. This is the second
# one, and it has to be nonzero or the whole discussion of K is theoretical:
# a K = 0 that never actually went wrong would mean the constant could be
# dropped with no measurable consequence.
ok(n_wrong > 0L, "mbic_pelt_at_K_zero_is_measurably_wrong",
   sprintf("%d of %d runs returned a cost above the exact minimum, worst by %.4f",
           n_wrong, n_runs, worst_excess))

# And the saving is not the reason to use the wrong constant: measure what
# the correct K costs in pruning.
s4 <- cp_sim_random(300L, 4L, "mean", effect = 3, gap = 20L, seed = 883001L)
co <- cp_cost(s4$y, "mean")
p4 <- cp_prepare(co, "mbic")
a0 <- cp_pelt(p4$cost, p4$beta, 2L, K = 0)
a1 <- cp_pelt(p4$cost, p4$beta, 2L)
# And the correct constant is not free. A negative K widens the keep
# window, so fewer candidates die and the search does more work. This is
# the trade MBIC actually makes and it is a tenfold one here, which is
# worth knowing before choosing the penalty.
ok(a1$ncand >= a0$ncand, "the_correct_K_prunes_no_more_than_K_zero",
   sprintf("K = log(4/n): %d evaluations, %d pruned; K = 0: %d evaluations, %d pruned; ratio %.1fx",
           a1$ncand, a1$npruned, a0$ncand, a0$npruned, a1$ncand / a0$ncand))

done()
