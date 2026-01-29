# Optimal Partitioning. Jackson et al. (2005), "An algorithm for optimal
# partitioning of data on an interval".
#
# F(t) is the best penalised cost of segmenting y[1..t]:
#
#     F(0) = -beta
#     F(t) = min over s of  F(s) + C(s+1, t) + beta
#
# with s ranging over admissible previous boundaries. F(0) is -beta rather
# than 0 so that beta is paid per CHANGE POINT rather than per segment: k
# change points make k+1 segments and k+1 betas, and the first one is
# cancelled. Nothing downstream can see the difference, but the reported
# total then matches the textbook objective and the penalties in
# R/penalty.R are the textbook constants.
#
# O(n^2) cost evaluations, every one of them O(1) from R/cost.R's prefix
# sums. This is the base case: PELT prunes the same recursion and must
# return the same answer, and brute force enumerates the same objective and
# must find the same minimum. A defect here is invisible to both unless they
# exist, which is why they do.

cp_op <- function(cost, beta, minseg = 1L) {
  n <- cost$n
  minseg <- as.integer(minseg)
  if (minseg < 1L) stop("cp_op: minseg must be at least 1")
  if (!(is.numeric(beta) && length(beta) == 1L && is.finite(beta)))
    stop("cp_op: beta must be one finite number")
  if (2L * minseg > n) {
    # No split is admissible. One segment, and say so rather than returning
    # an empty answer that looks like "no change found".
    return(structure(list(cp = integer(0), cost = cost$seg(0L, n),
                          F = NULL, ncand = 0L, npruned = 0L,
                          method = "op", beta = beta, minseg = minseg,
                          nsplit_possible = FALSE),
                     class = "cpfit"))
  }

  Fv <- rep(Inf, n + 1L)          # Fv[t + 1] is F(t)
  Fv[1L] <- -beta
  prev <- rep(NA_integer_, n + 1L)
  ncand <- 0L

  for (t in minseg:n) {
    s <- 0L:(t - minseg)
    vals <- Fv[s + 1L] + cost$seg(s, t) + beta
    ncand <- ncand + length(s)
    i <- which.min(vals)
    # which.min takes the FIRST minimum, so exact ties resolve to the
    # earliest change point. Said out loud because brute force has to break
    # ties the same way for the two to be compared by segmentation rather
    # than only by cost, and because with integer data exact ties happen.
    Fv[t + 1L] <- vals[i]
    prev[t + 1L] <- s[i]
  }

  structure(list(cp = cp_backtrack(prev, n), cost = Fv[n + 1L], F = Fv,
                 ncand = ncand, npruned = 0L, method = "op",
                 beta = beta, minseg = minseg, nsplit_possible = TRUE),
            class = "cpfit")
}

# Walk the predecessor array back from n. Shared with PELT, which builds the
# same array over a smaller candidate set.
cp_backtrack <- function(prev, n) {
  out <- integer(0)
  t <- n
  guard <- 0L
  while (t > 0L) {
    s <- prev[t + 1L]
    if (is.na(s)) stop("cp_backtrack: no predecessor for t = ", t,
                       "; the recursion left a gap")
    if (s > 0L) out <- c(s, out)
    t <- s
    guard <- guard + 1L
    if (guard > n + 1L) stop("cp_backtrack: predecessor chain does not decrease")
  }
  out
}
