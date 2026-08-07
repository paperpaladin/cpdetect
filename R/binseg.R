# Binary segmentation. Scott and Knott (1974); the method most people are
# actually running when they say "change point detection".
#
# Find the single split of a stretch that reduces the cost most, take it if
# the reduction beats beta, recurse into both halves. O(n log n) when it
# finds a balanced tree of change points and O(n^2) when it does not, with a
# constant small enough that it is faster than PELT on anything.
#
# It is here as a measurement, not as a feature. Binary segmentation is a
# greedy search and is NOT exact: the first split is chosen while the rest
# of the series is still being modelled as one piece, so two changes that
# partly cancel can hide each other and a first split in the wrong place is
# never revisited. How often that costs it the optimum is a fact about the
# method, and t/test_binseg.R measures it against cp_op rather than
# asserting it. The comparison exists because the exact searches do, which
# is the whole argument for having written them.

cp_binseg <- function(cost, beta, minseg = 1L) {
  n <- cost$n
  minseg <- as.integer(minseg)
  if (minseg < 1L) stop("cp_binseg: minseg must be at least 1")

  cps <- integer(0)
  ntested <- 0L
  nsplits <- 0L

  # Best split of y[(s+1)..e], or NULL if none is admissible or none helps.
  best_split <- function(s, e) {
    lo <- s + minseg
    hi <- e - minseg
    if (lo > hi) return(NULL)
    ts <- lo:hi
    gain <- cost$seg(s, e) -
      (cost$seg(rep(s, length(ts)), ts) + cost$seg(ts, rep(e, length(ts))))
    ntested <<- ntested + length(ts)
    i <- which.max(gain)
    list(t = ts[i], gain = gain[i])
  }

  recurse <- function(s, e) {
    b <- best_split(s, e)
    if (is.null(b) || b$gain <= beta) return(invisible(NULL))
    cps <<- c(cps, b$t)
    nsplits <<- nsplits + 1L
    recurse(s, b$t)
    recurse(b$t, e)
    invisible(NULL)
  }
  recurse(0L, n)
  cps <- sort(cps)

  structure(list(cp = cps, cost = cp_total(cost, cps) + length(cps) * beta,
                 F = NULL, ncand = ntested, npruned = 0L, method = "binseg",
                 beta = beta, minseg = minseg, nsplits = nsplits),
            class = "cpfit")
}
