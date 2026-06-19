# The generators, and whether the seed is real.
#
# Every detection rate in this project is measured on series from
# R/simulate.R, so a generator that quietly ignores an argument would make
# every one of those numbers a measurement of something else. The specific
# failure to guard against is the one recorded in the collection's notes: a
# seeded check that passes with the seed unwired, because it compares a
# report about the material rather than the material. The comparisons here
# are on the series itself, with identical(), not on a summary of it.
#
# Both directions are needed. Same seed giving the same series catches a
# seed that is ignored. Different seeds giving different series catches a
# series that is hardcoded. Either one alone passes under the other defect.

source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "harness.R"))
cp_load()

a1 <- cp_sim(50, c(20, 35), "mean", mean = c(0, 3, -2), seed = 123)
a2 <- cp_sim(50, c(20, 35), "mean", mean = c(0, 3, -2), seed = 123)
a3 <- cp_sim(50, c(20, 35), "mean", mean = c(0, 3, -2), seed = 124)
ok(identical(a1$y, a2$y), "same_seed_gives_the_identical_series",
   sprintf("%d values compared with identical()", length(a1$y)))
ok(!identical(a1$y, a3$y), "a_different_seed_gives_a_different_series",
   sprintf("%d of %d values differ", sum(a1$y != a3$y), length(a1$y)))

# The RNG state must be left moving when no seed is given, or a loop that
# seeds once and draws many series would get the same series every time.
set.seed(1)
b1 <- cp_sim(20, integer(0), "mean", seed = NULL)
b2 <- cp_sim(20, integer(0), "mean", seed = NULL)
ok(!identical(b1$y, b2$y), "no_seed_means_the_stream_advances")

# Across types, because each one draws from a different distribution and a
# seed wired into one branch says nothing about the others.
for (ty in c("mean", "var", "meanvar", "poisson")) {
  p <- cp_sim_random(80, 2L, ty, effect = 3, gap = 10L, seed = 5150L)
  q <- cp_sim_random(80, 2L, ty, effect = 3, gap = 10L, seed = 5150L)
  r <- cp_sim_random(80, 2L, ty, effect = 3, gap = 10L, seed = 5151L)
  ok(identical(p$y, q$y) && identical(p$cps, q$cps),
     paste0("cp_sim_random_is_reproducible_", ty))
  ok(!identical(p$y, r$y) || !identical(p$cps, r$cps),
     paste0("cp_sim_random_depends_on_the_seed_", ty))
}

# The change point positions must depend on the seed too, not only the
# noise. A generator that placed them at fixed positions and only varied
# the draws would pass every assertion above.
pos <- lapply(1:30, function(sd) cp_sim_random(100, 3L, "mean", gap = 10L,
                                               seed = 6000L + sd)$cps)
ok(length(unique(vapply(pos, paste, character(1), collapse = ","))) > 20,
   "change_point_positions_vary_with_the_seed",
   sprintf("%d distinct placements over 30 seeds",
           length(unique(vapply(pos, paste, character(1), collapse = ",")))))

# cp_place's constraint, checked on the output rather than trusted.
for (gap in c(1L, 5L, 12L)) {
  viol <- 0L
  for (sd in 1:40) {
    set.seed(7000L + sd)
    cps <- cp_place(100L, 4L, gap)
    if (length(cps) != 4L || any(diff(c(0L, cps, 100L)) < gap)) viol <- viol + 1L
  }
  ok(viol == 0L, sprintf("cp_place_respects_gap_%d", gap),
     sprintf("%d violations over 40 placements", viol))
}
ok(tryCatch({ cp_place(10L, 5L, 4L); FALSE }, error = function(e) TRUE),
   "cp_place_refuses_an_impossible_request")

# The parameters have to reach the data. effect is the one every
# measurement sweeps over, so it gets a direct check: the difference
# between consecutive segment means must scale with it.
jump <- function(ef) {
  s <- cp_sim(400, 200, "mean", mean = c(0, ef), sd = 1, seed = 8080L)
  abs(mean(s$y[201:400]) - mean(s$y[1:200]))
}
j1 <- jump(1); j4 <- jump(4)
ok(j4 > 3 * j1, "effect_size_reaches_the_data",
   sprintf("observed jump %.3f at effect 1, %.3f at effect 4", j1, j4))

sp <- cp_sim(2000, 1000, "poisson", lambda = c(2, 9), seed = 8081L)
ok(abs(mean(sp$y[1:1000]) - 2) < 0.3 && abs(mean(sp$y[1001:2000]) - 9) < 0.5,
   "poisson_rates_reach_the_data",
   sprintf("observed rates %.3f and %.3f against 2 and 9",
           mean(sp$y[1:1000]), mean(sp$y[1001:2000])))
ok(all(sp$y == floor(sp$y)) && all(sp$y >= 0), "poisson_draws_are_non_negative_integers")

sv <- cp_sim(4000, 2000, "var", mean = 0, sd = c(1, 4), seed = 8082L)
ok(abs(sd(sv$y[2001:4000]) / sd(sv$y[1:2000]) - 4) < 0.4,
   "variance_ratio_reaches_the_data",
   sprintf("observed ratio %.3f against 4", sd(sv$y[2001:4000]) / sd(sv$y[1:2000])))

ok(tryCatch({ cp_sim(10, c(5, 5), "mean"); FALSE }, error = function(e) TRUE),
   "cp_sim_refuses_repeated_change_points")
ok(tryCatch({ cp_sim(10, 10, "mean"); FALSE }, error = function(e) TRUE),
   "cp_sim_refuses_a_change_point_at_n")
ok(tryCatch({ cp_sim(10, 3, "nonsense"); FALSE }, error = function(e) TRUE),
   "cp_sim_refuses_an_unknown_type")

done()
