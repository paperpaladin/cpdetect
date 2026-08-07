# Calibrating the permutation test.
#
# A test whose p-values are not uniform under the null is not a test, and
# the only way to find out is to run it on a few hundred series that have
# no change point in them and look at the distribution. A test with the
# right size and no power is a random number generator with good manners,
# so the power curve is measured in the same file.
#
# The third measurement is the one that matters for the application. The
# permutation null assumes exchangeability, which for a time series means
# independence, and authentication counts per minute are not independent.
# The size of the resulting error is measured rather than warned about.

source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "harness.R"))
cp_load()

TOL <- 1e-9
B <- 199L
NSIM <- 200L
MS <- 5L
N <- 100L

# ------------------------------------------- the statistic, across languages
if (oracle_available()) {
  worst <- 0
  nbad <- 0L
  for (ty in c("mean", "var", "meanvar", "poisson")) for (sd in 1:3) {
    s <- cp_sim_random(60L, 1L, ty, effect = 2.5, gap = 10L, seed = 950000L + sd)
    co <- cp_cost(s$y, ty)
    got <- cp_amoc(co, 4L)
    req <- c(paste("type", ty), paste("sigma2", g17(co$sigma2)),
             paste("mu", g17(if (is.null(co$mu)) 0 else co$mu)),
             paste("varfloor", g17(co$varfloor)), "minseg 4",
             paste("y", paste(g17(co$y), collapse = " ")))
    out <- oracle("amoc", req, tag = paste0("amoc", ty, sd))
    pstat <- as.numeric(sub("^stat ", "", grep("^stat ", out, value = TRUE)))
    pat <- as.integer(sub("^at ", "", grep("^at ", out, value = TRUE)))
    worst <- max(worst, abs(got$stat - pstat))
    if (got$at != pat) nbad <- nbad + 1L
  }
  ok(worst <= TOL, "amoc_statistic_agrees_with_scipy",
     sprintf("worst %.3g over 12 series", worst))
  ok(nbad == 0L, "amoc_location_agrees_with_scipy", sprintf("%d disagreements", nbad))
} else {
  ok(FALSE, "scipy_oracle_present_for_amoc", "python3 with scipy and numpy is required")
}

# The statistic cannot be negative: the best split is at worst as good as
# no split, by the superadditivity already established in t/test_cost.R.
neg <- 0L
for (sd in 1:40) {
  s <- cp_sim(40L, integer(0), "mean", seed = 951000L + sd)
  if (cp_amoc(cp_cost(s$y, "mean"), 3L)$stat < -TOL) neg <- neg + 1L
}
ok(neg == 0L, "amoc_statistic_is_never_negative", sprintf("%d of 40 negative", neg))

# ----------------------------------------------------- size under the null
pv <- numeric(NSIM)
for (i in seq_len(NSIM)) {
  s <- cp_sim(N, integer(0), "mean", seed = 960000L + i)
  pv[i] <- cp_permtest(s$y, "mean", B = B, minseg = MS, seed = 970000L + i)$p
}
rej05 <- mean(pv <= 0.05)
rej10 <- mean(pv <= 0.10)
# Binomial standard error at alpha = 0.05 over 200 draws is 0.0154, so a
# true size of 0.05 lands inside 0.019..0.081 about 95 per cent of the
# time. The window below is three standard errors, chosen before the run.
ok(abs(rej05 - 0.05) <= 3 * sqrt(0.05 * 0.95 / NSIM), "size_at_alpha_0.05",
   sprintf("%.3f over %d null series, B = %d, 3 s.e. window is +-%.3f",
           rej05, NSIM, B, 3 * sqrt(0.05 * 0.95 / NSIM)))
ok(abs(rej10 - 0.10) <= 3 * sqrt(0.10 * 0.90 / NSIM), "size_at_alpha_0.10",
   sprintf("%.3f over %d null series", rej10, NSIM))
# Uniformity of the whole distribution, not only two of its quantiles.
# ks.test warns about ties because a permutation p-value takes only B+1
# distinct values; the statistic is still the right measurement and the
# warning is suppressed rather than the test replaced.
kst <- suppressWarnings(ks.test(pv, "punif"))
ok(kst$p.value > 0.01, "p_values_are_uniform_under_the_null",
   sprintf("KS D = %.4f, p = %.4f over %d null series", kst$statistic, kst$p.value, NSIM))
ok(all(pv >= 1 / (B + 1) - TOL) && all(pv <= 1 + TOL), "p_values_are_in_range",
   sprintf("min %.5f (floor is 1/(B+1) = %.5f), max %.5f",
           min(pv), 1 / (B + 1), max(pv)))

# ------------------------------------------------------------- power curve
pw <- list()
for (ef in c(0.25, 0.5, 1, 2)) {
  r <- numeric(60L)
  for (i in seq_len(60L)) {
    s <- cp_sim(N, N %/% 2L, "mean", mean = c(0, ef), seed = 980000L + i * 10L + ef * 4)
    r[i] <- cp_permtest(s$y, "mean", B = B, minseg = MS, seed = 981000L + i)$p <= 0.05
  }
  pw[[sprintf("%.2f", ef)]] <- mean(r)
  ok(TRUE, sprintf("power_at_effect_%.2f", ef),
     sprintf("%.3f rejected at alpha 0.05 over 60 series, n = %d, change at %d",
             mean(r), N, N %/% 2L))
}
ok(pw[["2.00"]] > pw[["0.25"]], "power_increases_with_the_effect_size",
   sprintf("%.3f at 0.25 standard deviations, %.3f at 2", pw[["0.25"]], pw[["2.00"]]))
ok(pw[["2.00"]] > 0.9, "the_test_has_power_at_a_two_sigma_jump",
   sprintf("%.3f", pw[["2.00"]]))

# ------------------------------------------- where the assumption breaks
#
# AR(1) noise with no change point at all. The permutation null breaks the
# serial dependence, so the permuted statistics are drawn from a
# distribution with less structure than the real series has, the reference
# set is too narrow, and the test rejects far too often. This is the
# single biggest reason a detector like this misbehaves on real telemetry,
# and the number is the point.
ar1 <- function(n, rho, seed) {
  set.seed(seed)
  e <- rnorm(n, 0, sqrt(1 - rho ^ 2))
  y <- numeric(n); y[1] <- rnorm(1)
  for (i in 2:n) y[i] <- rho * y[i - 1] + e[i]
  y
}
for (rho in c(0.3, 0.7)) {
  r <- numeric(80L)
  for (i in seq_len(80L))
    r[i] <- cp_permtest(ar1(N, rho, 990000L + i), "mean", B = B, minseg = MS,
                        seed = 991000L + i)$p <= 0.05
  ok(TRUE, sprintf("size_under_AR1_rho_%.1f", rho),
     sprintf("%.3f rejected at alpha 0.05 over 80 series with no change point",
             mean(r)))
  if (rho == 0.7)
    ok(mean(r) > 0.15, "autocorrelation_inflates_the_size_a_lot",
       sprintf("%.3f at rho 0.7 against a nominal 0.05", mean(r)))
}

# ------------------------------------------------------------ mechanics
a <- cp_permtest(cp_sim(60, 30, "mean", mean = c(0, 2), seed = 12)$y,
                 "mean", B = 99L, minseg = 4L, seed = 777)
b <- cp_permtest(cp_sim(60, 30, "mean", mean = c(0, 2), seed = 12)$y,
                 "mean", B = 99L, minseg = 4L, seed = 777)
c3 <- cp_permtest(cp_sim(60, 30, "mean", mean = c(0, 2), seed = 12)$y,
                  "mean", B = 99L, minseg = 4L, seed = 778)
ok(identical(a$nulls, b$nulls) && a$p == b$p, "permtest_is_reproducible_from_its_seed",
   sprintf("%d null statistics compared with identical()", length(a$nulls)))
# Written in exactly the shape that failed. cp_sim is called inline, so
# its own set.seed is evaluated lazily inside cp_permtest; before the
# force() calls were added that happened AFTER cp_permtest set its seed
# and both runs produced bit-identical permutations. Assigning the series
# to a local variable first hides the defect entirely, so the regression
# has to keep the inline form.
ok(!identical(a$nulls, c3$nulls), "permtest_depends_on_its_seed",
   sprintf("%d of %d null statistics differ between seeds 777 and 778",
           sum(a$nulls != c3$nulls), length(a$nulls)))
ok(a$hits == sum(a$nulls >= a$stat), "hits_is_the_count_it_claims_to_be")
ok(abs(a$p - (1 + a$hits) / (a$B + 1)) < TOL, "p_is_the_formula_it_claims_to_be")
ok(tryCatch({ cp_permtest(rnorm(10), "mean", B = 0L); FALSE }, error = function(e) TRUE),
   "permtest_refuses_zero_resamples")
ok(tryCatch({ cp_permtest(rnorm(6), "mean", B = 9L, minseg = 5L); FALSE },
            error = function(e) TRUE),
   "permtest_refuses_when_no_split_is_admissible")

done()
