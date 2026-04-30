# How often the greedy method is wrong, measured.
#
# Binary segmentation is what sits under most of the change point detection
# anyone actually runs, and it is not exact. The point of having written
# two exact searches is that the question "how often does greedy cost you
# the optimum" has an answer here rather than a citation.
#
# Two things are asserted and the rest is measured. Binary segmentation can
# never beat the exact minimum, which is a fact about the objective and
# would indicate a defect in cp_op if it were violated. And it must
# disagree somewhere, because a greedy search that matched an exhaustive
# one on every input in the sweep would mean the sweep is too easy.

source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "harness.R"))
cp_load()

TOL <- 1e-9

types <- c("mean", "var", "meanvar", "poisson")
effects <- c(0.75, 1.5, 3)
ks <- c(1L, 3L, 6L)

n_cases <- 0L
n_cpdiff <- 0L
n_worse <- 0L
n_better <- 0L
excess <- numeric(0)
by_effect <- list()

for (ty in types) for (ef in effects) for (k in ks) for (rep in 1:4) {
  n_cases <- n_cases + 1L
  sim <- cp_sim_random(240L, k, ty, effect = ef, gap = 15L, seed = 700000L + n_cases)
  co <- cp_cost(sim$y, ty)
  pr <- cp_prepare(co, "bic")
  a <- cp_op(pr$cost, pr$beta, 5L)
  g <- cp_binseg(pr$cost, pr$beta, 5L)
  if (!identical(as.integer(a$cp), as.integer(g$cp))) n_cpdiff <- n_cpdiff + 1L
  d <- g$cost - a$cost
  if (d > TOL) { n_worse <- n_worse + 1L; excess <- c(excess, d) }
  if (d < -TOL) n_better <- n_better + 1L
  key <- sprintf("effect_%.2f", ef)
  by_effect[[key]] <- c(by_effect[[key]], as.numeric(d > TOL))
}

ok(n_cases == 144L, "binseg_sweep_ran", sprintf("%d cases", n_cases))
ok(n_better == 0L, "binseg_never_beats_the_exact_minimum",
   sprintf("%d cases found a lower cost than cp_op, which would be a defect in cp_op",
           n_better))
ok(n_worse > 0L, "binseg_is_measurably_worse_somewhere",
   sprintf("worse in %d of %d cases (%.1f%%), differs in segmentation in %d (%.1f%%)",
           n_worse, n_cases, 100 * n_worse / n_cases,
           n_cpdiff, 100 * n_cpdiff / n_cases))
ok(length(excess) == n_worse, "excess_recorded_for_every_worse_case",
   if (n_worse) sprintf("excess penalised cost: median %.4f, max %.4f",
                        median(excess), max(excess)) else "none")
for (key in sort(names(by_effect)))
  ok(TRUE, paste0("binseg_failure_rate_", key),
     sprintf("%d of %d cases worse than exact (%.1f%%)",
             sum(by_effect[[key]]), length(by_effect[[key]]),
             100 * mean(by_effect[[key]])))

# The structural reason, isolated. Binary segmentation chooses its first
# split while modelling the rest of the series as one piece, so a pattern
# whose two halves cancel in the overall mean hides itself from the first
# cut. A plateau that goes up and comes back down symmetrically is the
# standard example: the best single split of the whole series is a
# compromise in the middle of the plateau, and the two real boundaries are
# then looked for inside stretches that no longer contain both.
set.seed(702001)
hide <- c(rnorm(60, 0, 1), rnorm(30, 2.2, 1), rnorm(60, 0, 1))
co <- cp_cost(hide, "mean")
pr <- cp_prepare(co, "bic")
a <- cp_op(pr$cost, pr$beta, 5L)
g <- cp_binseg(pr$cost, pr$beta, 5L)
ok(g$cost >= a$cost - TOL, "plateau_binseg_not_better",
   sprintf("exact [%s] %.4f, greedy [%s] %.4f",
           paste(a$cp, collapse = ","), a$cost, paste(g$cp, collapse = ","), g$cost))
ok(TRUE, "plateau_recorded",
   sprintf("true change points 60 and 90; exact finds [%s], greedy finds [%s]",
           paste(a$cp, collapse = ","), paste(g$cp, collapse = ",")))

# Cost, which is the reason anyone uses it.
sim <- cp_sim_random(1000L, 6L, "mean", effect = 2.5, gap = 40L, seed = 703001L)
co <- cp_cost(sim$y, "mean")
pr <- cp_prepare(co, "bic")
tm <- function(f) as.numeric(system.time(f())[["elapsed"]])
a <- cp_op(pr$cost, pr$beta, 5L)
b <- cp_pelt(pr$cost, pr$beta, 5L)
g <- cp_binseg(pr$cost, pr$beta, 5L)
ok(g$ncand < b$ncand && b$ncand < a$ncand, "evaluation_counts_are_ordered",
   sprintf("n = 1000: binseg %d, pelt %d, op %d", g$ncand, b$ncand, a$ncand))

# The easy case, written expecting agreement and not getting it, which is
# the more useful result. 1000 points, six change points, separated by at
# least 40, a jump of 2.5 standard deviations: greedy still places two of
# the six two points away from where the exact search does. Not a large
# error and not zero, and the point of the exact searches is that the
# difference is a measured number rather than an assumption.
ok(g$cost >= a$cost - TOL, "binseg_on_the_easy_case_is_not_better",
   sprintf("greedy [%s] %.5f, exact [%s] %.5f, excess %.5f",
           paste(g$cp, collapse = ","), g$cost,
           paste(a$cp, collapse = ","), a$cost, g$cost - a$cost))
ok(sum(abs(sort(g$cp) - sort(a$cp))) <= 6L, "binseg_on_the_easy_case_is_close",
   sprintf("total absolute displacement %d over %d change points",
           sum(abs(sort(g$cp) - sort(a$cp))), length(a$cp)))

done()
