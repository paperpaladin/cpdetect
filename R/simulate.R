# Known-answer series.
#
# Everything measured about detection quality in this project is measured
# against change points placed here, so this file is as load-bearing as any
# of the algorithms and has the same obligation: the seed has to be real.
# A seeded generator that ignores its seed passes every "same seed, same
# output" check that compares a report rather than the material, so
# t/test_simulate.R compares the series itself.
#
# set.seed() is called when a seed is given and the caller's RNG state is
# NOT restored, deliberately. Restoring it would make a loop over seeds look
# reproducible while the data it produced depended on where the loop
# started. Every caller here passes an explicit seed for every draw.

# Place k change points in 1..n-1 with at least `gap` between them and at
# least `gap` at each end. Returns them sorted and strictly increasing.
#
# Rejection sampling rather than a closed form: the closed form (choose k
# from n - (k+1)*gap + k positions and spread) is correct and is the kind of
# correct that is wrong by one somewhere, and this is test-corpus code where
# a retry costs nothing. It gives up rather than looping forever, because a
# generator that cannot meet its own constraints must say so.
cp_place <- function(n, k, gap = 1L, tries = 10000L) {
  if (k == 0L) return(integer(0))
  if ((k + 1L) * gap > n)
    stop("cp_place: ", k, " change points with gap ", gap, " do not fit in ", n)
  for (i in seq_len(tries)) {
    cps <- sort(sample.int(n - 1L, k))
    b <- c(0L, cps, n)
    if (all(diff(b) >= gap)) return(as.integer(cps))
  }
  stop("cp_place: could not place ", k, " change points with gap ", gap,
       " in ", n, " after ", tries, " tries")
}

# A series with the given change points and per-segment parameters. Each
# parameter vector is recycled to length(cps) + 1.
cp_sim <- function(n, cps = integer(0), type = "mean",
                   mean = 0, sd = 1, lambda = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cps <- as.integer(cps)
  b <- c(0L, cps, as.integer(n))
  if (any(diff(b) < 1L)) stop("cp_sim: change points must be strictly increasing within 1..n-1")
  k <- length(cps)
  lens <- diff(b)
  rep_par <- function(v) rep_len(as.numeric(v), k + 1L)
  mu <- rep_par(mean); s <- rep_par(sd); lam <- rep_par(lambda)
  y <- switch(type,
    mean = ,
    meanvar = unlist(lapply(seq_along(lens), function(i) rnorm(lens[i], mu[i], s[i]))),
    var = unlist(lapply(seq_along(lens), function(i) rnorm(lens[i], mu[i], s[i]))),
    poisson = unlist(lapply(seq_along(lens), function(i) rpois(lens[i], lam[i]))),
    stop("cp_sim: unknown type '", type, "'"))
  list(y = y, cps = cps, type = type, n = as.integer(n), seed = seed,
       par = list(mean = mu, sd = s, lambda = lam))
}

# A random instance: k change points at random positions with random
# per-segment parameters. `effect` is the size of the jump between
# consecutive segments, in units of the within-segment standard deviation
# for the Gaussian types and as a multiplicative factor for Poisson. It is
# an argument because the whole of the failure-mode measurement is a sweep
# over it.
cp_sim_random <- function(n, k, type = "mean", effect = 3, gap = 1L,
                          base_sd = 1, base_lambda = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cps <- cp_place(n, k, gap)
  m <- k + 1L
  signs <- sample(c(-1, 1), m, replace = TRUE)
  if (type == "poisson") {
    # Alternating rather than a random walk: a walk in +-effect steps runs
    # away to 10^9 over a dozen change points and the later segments stop
    # being counts of anything. Alternating keeps every consecutive pair a
    # factor of `effect` apart, which is the thing being varied.
    lam <- base_lambda * ifelse(seq_len(m) %% 2 == 1, 1, effect)
    out <- cp_sim(n, cps, "poisson", lambda = lam, seed = NULL)
  } else if (type == "var") {
    sds <- base_sd * ifelse(seq_len(m) %% 2 == 1, 1, effect)
    out <- cp_sim(n, cps, "var", mean = 0, sd = sds, seed = NULL)
  } else if (type == "meanvar") {
    mus <- cumsum(c(0, signs[-1])) * effect * base_sd
    sds <- base_sd * ifelse(seq_len(m) %% 2 == 1, 1, max(effect, 1.5))
    out <- cp_sim(n, cps, "meanvar", mean = mus, sd = sds, seed = NULL)
  } else {
    mus <- cumsum(c(0, signs[-1])) * effect * base_sd
    out <- cp_sim(n, cps, "mean", mean = mus, sd = base_sd, seed = NULL)
  }
  out$seed <- seed
  out$effect <- effect
  out
}
