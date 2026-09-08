# Recovery AND false positives, in the suite, at a size that can run on
# every commit. The full table is `make recovery`; this is the part that
# has to stay true.
#
# Both numbers for every cell. The assertions are the shape of the curve
# rather than its exact height, because the heights move with n, the match
# window and the minimum segment length and pinning them would be pinning
# a configuration. The heights are printed with every assertion so a
# regression shows up as a number that moved.

source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "harness.R"))
cp_load()

N <- 220L
GAP <- 30L
W <- 8L
MINSEG <- 15L
SEEDS <- 12L
BASE <- 1200000L

cell <- function(type, effect, k, penalty = "bic", method = "pelt", tag) {
  res <- lapply(seq_len(SEEDS), function(i) {
    s <- cp_sim_random(N, k, type, effect = effect, gap = GAP,
                       seed = BASE + tag * 100L + i)
    list(true = s$cps,
         found = cp_detect(s$y, type, method = method, penalty = penalty,
                           minseg = MINSEG)$cp)
  })
  cp_score(res, W)
}

report <- function(sc, label) {
  per <- sc$recovery_per_series
  sprintf("recovery %s (sd %s), false positives %.2f/series (sd %.2f), any-FP %.2f, displacement %s",
          if (is.na(sc$recovery)) "-" else sprintf("%.3f", sc$recovery),
          if (all(is.na(per))) "-" else sprintf("%.3f", stats::sd(per, na.rm = TRUE)),
          sc$fp_per_series, sc$fp_sd, sc$any_fp_rate,
          if (is.na(sc$disp_mean)) "-" else sprintf("%.1f", sc$disp_mean))
}

tag <- 0L
for (type in c("mean", "var", "meanvar", "poisson")) {
  effs <- if (type == "poisson") c(1.3, 3) else c(0.6, 3)
  low <- cell(type, effs[1], 3L, tag = (tag <- tag + 1L))
  high <- cell(type, effs[2], 3L, tag = (tag <- tag + 1L))
  nul <- cell(type, effs[2], 0L, tag = (tag <- tag + 1L))

  ok(TRUE, paste0("recovery_", type, "_small_effect"), report(low, type))
  ok(TRUE, paste0("recovery_", type, "_large_effect"), report(high, type))
  ok(TRUE, paste0("recovery_", type, "_null_series"), report(nul, type))

  ok(high$recovery > low$recovery, paste0("recovery_rises_with_effect_", type),
     sprintf("%.3f at effect %.2f, %.3f at effect %.2f",
             low$recovery, effs[1], high$recovery, effs[2]))
  ok(high$recovery >= 0.9, paste0("recovery_is_high_at_a_large_effect_", type),
     sprintf("%.3f over %d series with 3 change points each", high$recovery, SEEDS))
  # The other half, which is the half a single accuracy number hides.
  ok(nul$fp_per_series <= 0.5, paste0("false_positives_on_null_series_", type),
     sprintf("%.3f per series, %.2f of series have any", nul$fp_per_series, nul$any_fp_rate))
  ok(high$precision >= 0.8, paste0("precision_at_a_large_effect_", type),
     sprintf("%.3f: %d of %d detections matched a true change point",
             high$precision, high$tp, high$tp + high$fp))
}

# The failure modes, each with its own number. These are the cases the
# README has to describe, and describing them needs them measured.
#
# A change point close to the end of the series has less data on one side
# of it, so it is harder to find than one in the middle at the same effect
# size.
edge <- function(pos, tag) {
  res <- lapply(seq_len(20L), function(i) {
    s <- cp_sim(N, pos, "mean", mean = c(0, 1.5),
                seed = BASE + tag * 100L + i)
    list(true = s$cps, found = cp_detect(s$y, "mean", penalty = "bic",
                                         minseg = 5L)$cp)
  })
  cp_score(res, W)
}
mid <- edge(125L, 50L)
near <- edge(20L, 51L)
ok(TRUE, "recovery_change_in_the_middle", report(mid, "mid"))
ok(TRUE, "recovery_change_20_from_the_start", report(near, "edge"))
ok(mid$recovery > near$recovery, "an_edge_change_is_harder",
   sprintf("%.3f at index 125, %.3f at index 20", mid$recovery, near$recovery))

# Two changes close together: a short segment between two boundaries is
# the shape a minimum segment length is supposed to protect, and it is
# also the shape that gets missed.
# The effect size has to be in the range where the answer is in doubt.
# Written first at 2.5 standard deviations, where 8 apart is recovered at
# 0.975 just like 60 apart and the assertion failed because there was
# nothing to see. At 1.2 the curve is visible.
adj <- function(gap2, tag) {
  res <- lapply(seq_len(25L), function(i) {
    s <- cp_sim(N, c(100L, 100L + gap2), "mean", mean = c(0, 1.2, 0),
                seed = BASE + tag * 100L + i)
    list(true = s$cps, found = cp_detect(s$y, "mean", penalty = "bic",
                                         minseg = 3L)$cp)
  })
  cp_score(res, 4L)
}
adjr <- vapply(c(60L, 16L, 8L, 4L),
               function(g) adj(g, 60L + g)$recovery, numeric(1))
for (i in seq_along(c(60L, 16L, 8L, 4L)))
  ok(TRUE, sprintf("recovery_two_changes_%d_apart", c(60L, 16L, 8L, 4L)[i]),
     sprintf("%.3f over 25 series, jump 1.2 standard deviations", adjr[i]))
ok(adjr[1] > adjr[3] && adjr[3] > adjr[4], "adjacent_changes_are_harder",
   sprintf("recovery %.3f, %.3f, %.3f, %.3f at separations 60, 16, 8, 4",
           adjr[1], adjr[2], adjr[3], adjr[4]))

done()
