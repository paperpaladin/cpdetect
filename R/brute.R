# Exhaustive enumeration.
#
# Not a second implementation of the search. For a series short enough, the
# set of admissible segmentations is finite and small, and scoring all of
# them and taking the minimum IS the answer, by the definition of what the
# dynamic program is trying to compute. It is the only oracle in this
# project that cannot be wrong about the optimum; it can only be wrong about
# the cost function, which is why the cost function has scipy.
#
# Enumeration is by recursion over the boundaries, in increasing k and
# lexicographic order within k, which is the order cp_op's which.min
# resolves a tie in, so the two can be compared as segmentations and not
# only as costs.
#
# The generator was a combn over all subsets of the interior boundaries
# with the minimum segment length applied afterwards as a filter. That is
# the methodologically safer way round -- a generator that knows about
# minseg can share an off-by-one with the dynamic program it is checking,
# and a filter written from the definition cannot -- and it is unusable: at
# n = 22 and minseg 2 it generated 2097152 subsets to keep 10946, and the
# sweep took 89 seconds. The generator knows about minseg now, and two
# things stand in for the filter it replaced. The filter itself is still
# here, inside consider(), where it costs nothing and would fire if the
# recursion ever emitted an inadmissible set. And t/test_brute.R checks the
# number of segmentations scored against the closed form for the number of
# compositions of n into parts of size at least minseg, per case, so a
# generator that silently skipped part of the space would be caught by
# combinatorics rather than by inspection.
#
# Cost: 2^(n-1) segmentations at minseg 1, which is why maxcp exists; but
# only Fibonacci(n-1) at minseg 2, which is 46368 at n = 25 and needs no cap
# at all.

cp_brute <- function(cost, beta, minseg = 1L, maxcp = Inf) {
  n <- cost$n
  minseg <- as.integer(minseg)
  if (minseg < 1L) stop("cp_brute: minseg must be at least 1")
  kmax <- min(maxcp, floor(n / minseg) - 1L, n - 1L)
  if (!is.finite(kmax) || kmax < 0) kmax <- 0

  # Every segment cost, computed once. There are only n(n+1)/2 of them and
  # an enumeration asks for the same ones tens of thousands of times: at
  # n = 22 and minseg 2 there are 10946 segmentations over 253 distinct
  # segments. This is memoisation, not a shortcut through the oracle. The
  # same cost$seg is called, the same number of distinct times, and the
  # summation order over a segmentation is unchanged; it took the whole
  # sweep from 109 seconds to 11.
  CM <- matrix(NA_real_, n + 1L, n)
  for (e in seq_len(n)) {
    ss <- 0L:(e - 1L)
    CM[ss + 1L, e] <- cost$seg(ss, e)
  }

  best <- Inf
  best_cp <- integer(0)
  nseg <- 0L

  consider <- function(cp) {
    b <- c(0L, cp, n)
    if (any(diff(b) < minseg)) return(invisible(NULL))
    nseg <<- nseg + 1L
    # Summed left to right over the segments, which is the same order as
    # cp_total. A different order would differ in the last bits and the
    # agreement check would be measuring summation order.
    c0 <- sum(CM[cbind(b[-length(b)] + 1L, b[-1L])]) + length(cp) * beta
    # Strict improvement only, so the FIRST segmentation attaining the
    # minimum wins. combn emits in lexicographic order and cp_op's
    # which.min takes the first minimum, so both resolve an exact tie to the
    # earliest change point and the two can be compared as segmentations.
    if (c0 < best) { best <<- c0; best_cp <<- cp }
    invisible(NULL)
  }

  # emit(prefix, last, k): choose k more boundaries after `last`. The upper
  # bound n - k*minseg is what makes this enumerate only admissible sets:
  # with k boundaries still to place after t, the shortest the rest of the
  # series can be is k*minseg.
  emit <- function(prefix, last, k) {
    if (k == 0L) return(consider(prefix))
    lo <- last + minseg
    hi <- n - k * minseg
    if (lo > hi) return(invisible(NULL))
    for (t in lo:hi) emit(c(prefix, t), t, k - 1L)
    invisible(NULL)
  }
  for (k in 0L:max(0L, kmax)) emit(integer(0), 0L, k)

  structure(list(cp = best_cp, cost = best, F = NULL, ncand = nseg,
                 npruned = 0L, method = "brute", beta = beta,
                 minseg = minseg, maxcp = kmax, nconsidered = nseg),
            class = "cpfit")
}
