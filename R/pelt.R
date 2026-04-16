# PELT: Pruned Exact Linear Time. Killick, Fearnhead and Eckley (2012),
# "Optimal detection of changepoints with a linear computational cost".
#
# The same recursion as Optimal Partitioning, over a shrinking candidate
# set. The pruning rule is
#
#     if  F(s) + C(s+1, t) + K  >  F(t)   then s can be discarded for good,
#
# and the theorem that makes it exact rather than approximate needs a
# constant K with
#
#     C(s+1, t) + C(t+1, T) + K  <=  C(s+1, T)      for all s < t < T.
#
# For a cost that is minus twice a maximised log likelihood, K = 0: splitting
# a segment gives the fit more freedom, so two segments never cost more than
# one. That argument is not a proof about THIS code, so t/test_cost.R checks
# superadditivity over every split of every segment of several series, and
# t/test_pelt.R checks the pruning the only way that matters, by requiring
# PELT and Optimal Partitioning to return the same segmentation.
#
# Three things that are easy to get wrong here and each has its own check.
#
# The inequality is strict. s is discarded only when the left side is
# STRICTLY greater; discarding on equality drops candidates that are still
# optimal. This is the classic bug, and the only thing that can see it is an
# exact comparison against an unpruned search on data where the tie occurs.
#
# The count has to be reported. A PELT whose pruning never fires is Optimal
# Partitioning with a longer name, and a test of it proves nothing about the
# rule. Three counts come back: npruned, how many candidates were discarded
# for good; nadded, how many ever entered the live set, which is the
# denominator npruned is out of; and ncand, the number of cost evaluations,
# which is what PELT actually saves and is directly comparable with
# cp_op's. A candidate is pruned at most once, so npruned can never exceed
# n, and the saving shows up in ncand, not in npruned.
#
# K is not always zero. MBIC's log(segment length) term is not
# superadditive, so under MBIC the valid constant is log(4/n), and R/
# penalty.R attaches it to the cost object so that forgetting is not
# possible from the caller's side.
#
# ------------------------------------------- the theorem needs minseg = 1
#
# The proof that a discarded s is never optimal again runs
#
#   F(s) + C(s+1,T) + beta
#     >= F(s) + C(s+1,t) + C(t+1,T) + K + beta        (superadditivity)
#     >= F(t) + C(t+1,T) + beta                       (s failed the test)
#     >= F(T)                                         (t is a candidate for T)
#
# and the last step is false when T - t < minseg, because then t is not an
# admissible last change point for T and F(T) was never offered it. So a
# candidate that fails the test at step t is dominated only from step
# t + minseg onwards, and removing it immediately throws away answers that
# are still reachable.
#
# This is not a hypothetical. The first version removed on the spot and the
# comparison against Optimal Partitioning caught it at once: with minseg 2
# on a 350-point series it returned a cost of 754.0822 where the exact
# answer is 753.8801, and it was wrong on six of the forty-eight cases in
# t/test_pelt.R, every one of them at minseg above 1 and none at minseg 1.
# Killick et al. state the algorithm with no minimum segment length and the
# gap is invisible until something enumerates.
#
# Candidates therefore carry a death step rather than being deleted. The
# saving is unchanged in the regime that matters, since minseg is small
# next to n.

# prune_delay exists so that the defect described above can be put back.
# Setting it to 1 is the immediate removal that was wrong; the default,
# NULL, means minseg, which is the correct value. t/test_pelt.R runs both
# and requires the first to disagree with Optimal Partitioning, because a
# fix nobody has watched fail is a fix nobody knows is load-bearing.
cp_pelt <- function(cost, beta, minseg = 1L, K = NULL, prune_delay = NULL) {
  n <- cost$n
  minseg <- as.integer(minseg)
  if (minseg < 1L) stop("cp_pelt: minseg must be at least 1")
  if (!(is.numeric(beta) && length(beta) == 1L && is.finite(beta)))
    stop("cp_pelt: beta must be one finite number")
  if (is.null(K)) K <- if (is.null(cost$pelt_K)) 0 else cost$pelt_K
  if (!(is.numeric(K) && length(K) == 1L && is.finite(K)))
    stop("cp_pelt: K must be one finite number")
  if (is.null(prune_delay)) prune_delay <- minseg
  prune_delay <- as.integer(prune_delay)
  if (prune_delay < 1L) stop("cp_pelt: prune_delay must be at least 1")

  if (2L * minseg > n)
    return(structure(list(cp = integer(0), cost = cost$seg(0L, n),
                          F = NULL, ncand = 0L, npruned = 0L, nadded = 0L,
                          method = "pelt", beta = beta, minseg = minseg,
                          K = K, nsplit_possible = FALSE),
                     class = "cpfit"))

  Fv <- rep(Inf, n + 1L)
  Fv[1L] <- -beta
  prev <- rep(NA_integer_, n + 1L)

  live <- 0L                 # candidate set R_t, always sorted increasing
  dead <- Inf                # live[i] stops being usable at step dead[i]
  npruned <- 0L
  nadded <- 1L               # candidates ever admitted to the live set
  ncand <- 0L                # cost evaluations: the thing PELT saves

  for (t in minseg:n) {
    alive <- dead > t
    live <- live[alive]
    dead <- dead[alive]
    far <- t - live >= minseg
    cands <- live[far]
    held <- live[!far]       # too recent to end a segment at t, but alive
    held_dead <- dead[!far]
    if (length(cands) == 0L) {
      # Reachable only when minseg > 1 and every live candidate is still
      # too recent. F(t) stays Inf, which is correct: no admissible
      # segmentation of y[1..t] exists. The candidate t is still added.
      live <- c(held, t)
      dead <- c(held_dead, Inf)
      nadded <- nadded + 1L
      next
    }
    base <- Fv[cands + 1L] + cost$seg(cands, t)
    vals <- base + beta
    ncand <- ncand + length(cands)
    i <- which.min(vals)
    Fv[t + 1L] <- vals[i]
    prev[t + 1L] <- cands[i]

    # Strictly greater is discarded; the <= here is the whole theorem.
    # Strictly greater is discarded, and discarded means "from step
    # t + minseg", not "now". A candidate already condemned keeps its
    # earlier death step: the test only ever brings the date forward.
    cd <- dead[far]
    fail <- (base + K) > Fv[t + 1L]
    cd[fail] <- pmin(cd[fail], t + prune_delay)
    npruned <- npruned + sum(fail & is.infinite(dead[far]))
    live <- c(cands, held, t)
    dead <- c(cd, held_dead, Inf)
    o <- order(live)
    live <- live[o]
    dead <- dead[o]
    nadded <- nadded + 1L
  }

  if (!is.finite(Fv[n + 1L]))
    stop("cp_pelt: no admissible segmentation of ", n, " points at minseg ", minseg)

  structure(list(cp = cp_backtrack(prev, n), cost = Fv[n + 1L], F = Fv,
                 ncand = ncand, npruned = npruned, nadded = nadded,
                 method = "pelt", beta = beta, minseg = minseg, K = K,
                 prune_delay = prune_delay,
                 nsplit_possible = TRUE),
            class = "cpfit")
}
