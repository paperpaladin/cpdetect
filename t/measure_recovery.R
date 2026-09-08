#!/usr/bin/env Rscript
# make recovery -- the detection quality numbers quoted in the README.
#
# Both rates, always, and the spread across seeds. A detector that finds
# every real change and twenty spurious ones is not a good detector, and a
# single accuracy number cannot tell you which one you have.
#
# Everything here is seeded from a fixed base so the table is reproducible;
# the seed of every cell is printed with it.

suppressWarnings(source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."),
                                  "t", "harness.R")))
cp_load()

BASE <- 1100000L
N <- 400L
GAP <- 40L
W <- 10L          # a detection within 10 of a true change point counts
MINSEG <- 20L
SEEDS <- as.integer(Sys.getenv("CPDETECT_SEEDS", unset = "30"))

cell <- function(type, effect, k, method, penalty, nseeds = SEEDS, tag = 0L) {
  res <- vector("list", nseeds)
  for (i in seq_len(nseeds)) {
    sd <- BASE + tag * 1000L + i
    s <- cp_sim_random(N, k, type, effect = effect, gap = GAP, seed = sd)
    f <- cp_detect(s$y, type, method = method, penalty = penalty, minseg = MINSEG)
    res[[i]] <- list(true = s$cps, found = f$cp)
  }
  sc <- cp_score(res, W)
  per <- sc$recovery_per_series
  list(type = type, effect = effect, k = k, method = method, penalty = penalty,
       nseeds = nseeds, seed_from = BASE + tag * 1000L + 1L,
       recovery = sc$recovery,
       recovery_sd = if (all(is.na(per))) NA_real_ else stats::sd(per, na.rm = TRUE),
       fp_per_series = sc$fp_per_series, fp_sd = sc$fp_sd,
       any_fp = sc$any_fp_rate, precision = sc$precision,
       disp = sc$disp_mean)
}

cat(sprintf("cpdetect recovery and false positives\n"))
cat(sprintf("n = %d, minimum segment %d, match window %d, %d seeds per cell\n\n",
            N, MINSEG, W, SEEDS))
hdr <- sprintf("%-8s %-7s %-6s %3s %7s %7s %7s %7s %7s %7s %6s",
               "cost", "method", "pen", "k", "effect", "recov", "rec.sd",
               "fp/ser", "fp.sd", "anyfp", "disp")
cat(hdr, "\n", strrep("-", nchar(hdr)), "\n", sep = "")

tag <- 0L
rows <- list()
for (type in c("mean", "var", "meanvar", "poisson")) {
  effs <- if (type == "poisson") c(1.25, 1.5, 2, 4) else c(0.5, 1, 2, 4)
  for (method in c("pelt", "binseg")) for (pen in c("bic", "mbic")) {
    if (method == "binseg" && pen == "mbic") next
    for (k in c(0L, 2L, 5L)) for (ef in (if (k == 0L) effs[1] else effs)) {
      tag <- tag + 1L
      r <- cell(type, ef, k, method, pen, tag = tag)
      rows[[length(rows) + 1L]] <- r
      cat(sprintf("%-8s %-7s %-6s %3d %7.2f %7s %7s %7.3f %7.3f %7.3f %6s\n",
                  r$type, r$method, r$penalty, r$k, r$effect,
                  if (is.na(r$recovery)) "-" else sprintf("%.3f", r$recovery),
                  if (is.na(r$recovery_sd)) "-" else sprintf("%.3f", r$recovery_sd),
                  r$fp_per_series, r$fp_sd, r$any_fp,
                  if (is.na(r$disp)) "-" else sprintf("%.1f", r$disp)))
    }
  }
}
cat("\nrecov   pooled recovery rate: recovered true change points / true change points\n")
cat("rec.sd  standard deviation of the per-series recovery rate across seeds\n")
cat("fp/ser  spurious change points per series (detections matching no true one)\n")
cat("anyfp   fraction of series with at least one spurious change point\n")
cat("disp    mean |detected - true| over recovered change points, in samples\n")
