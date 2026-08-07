# Is the change real?
#
# A detected change point is the argmax of a statistic over n - 1 candidate
# positions, so its apparent size is inflated by the maximisation: on pure
# noise, the best of 199 splits looks impressive. The null distribution of
# the MAXIMUM is what a p-value here has to be computed against, and under
# an independence null the permutation distribution of the series gives it
# without assuming anything about the shape of the data.
#
# The statistic is the cost reduction from the best admissible single
# split, which is the likelihood ratio for at-most-one-change against none:
#
#     D = C(1..n) - min over t of [ C(1..t) + C(t+1..n) ]
#
# Under the null, permuting the series destroys any change in distribution
# over time while leaving the marginal distribution exactly as it was, so
# the permuted statistics are draws from the null distribution of D with
# the real marginal rather than an assumed one. p = (1 + hits) / (B + 1),
# the standard conservative form: the observed value is counted as one of
# its own reference set, so p is never 0 and the test is valid at any B.
#
# What this does NOT test. It assumes exchangeability, which is
# independence here. Authentication failures per minute are autocorrelated,
# and against autocorrelated data a permutation null is too narrow and the
# test rejects far too often. That is measured in t/test_permtest.R rather
# than mentioned, because it is the single largest reason a detector like
# this misbehaves on real telemetry.

# The at-most-one-change statistic and where it is attained.
cp_amoc <- function(cost, minseg = 1L) {
  n <- cost$n
  minseg <- as.integer(minseg)
  lo <- minseg
  hi <- n - minseg
  if (lo > hi) return(list(stat = -Inf, at = NA_integer_, ncand = 0L))
  ts <- lo:hi
  whole <- cost$seg(0L, n)
  gain <- whole - (cost$seg(rep(0L, length(ts)), ts) +
                   cost$seg(ts, rep(n, length(ts))))
  i <- which.max(gain)
  list(stat = gain[i], at = ts[i], ncand = length(ts))
}

cp_permtest <- function(y, type = "mean", B = 499L, minseg = 1L,
                        seed = NULL, ...) {
  # Force every argument before touching the random number generator.
  #
  # R evaluates arguments lazily, at the point of first use inside the
  # function. set.seed() used to be the first statement here, and the
  # natural way to call this -- cp_permtest(cp_sim(...)$y, seed = 778) --
  # therefore ran cp_sim's own set.seed AFTER ours, wiping it. Two calls
  # differing only in the seed produced bit-identical permutations and the
  # test reported a reproducible answer, which is exactly what it would
  # report if the seed were unwired. It is the lazy-evaluation form of the
  # defect the collection's notes describe: the seed looks honoured
  # because the check cannot see the difference.
  #
  # The assertion that caught it compares the null statistics themselves
  # across two seeds, so the symptom was visible. Calling the same code
  # with the series in a local variable first hid it completely, which is
  # why the regression in t/test_permtest.R passes the generator call
  # inline, in the shape that fails.
  force(y); force(type); force(B); force(minseg)
  B <- as.integer(B)
  if (B < 1L) stop("cp_permtest: B must be at least 1")
  n <- length(y)
  obs <- cp_amoc(cp_cost(y, type, ...), minseg)
  if (!is.finite(obs$stat))
    stop("cp_permtest: no admissible split at minseg ", minseg, " for n = ", n)
  if (!is.null(seed)) set.seed(seed)
  nulls <- numeric(B)
  for (b in seq_len(B))
    nulls[b] <- cp_amoc(cp_cost(sample(y), type, ...), minseg)$stat
  # >= rather than >. A permuted series that ties the observed statistic is
  # evidence against the alternative, not for it, and with integer counts
  # ties are common enough to change the answer.
  hits <- sum(nulls >= obs$stat)
  list(stat = obs$stat, at = obs$at, B = B, hits = hits,
       p = (1 + hits) / (B + 1), nulls = nulls, type = type,
       minseg = minseg, seed = seed)
}
