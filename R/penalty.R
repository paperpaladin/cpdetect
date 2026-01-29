# Penalties, and the one that breaks PELT.
#
# Every search here minimises
#
#     sum over segments of C(segment) + k * beta
#
# where k is the number of change points. beta is the entire defence against
# a detector that returns a change point per observation, since the cost of
# any segmentation is non-increasing in k. Which beta is a modelling choice
# and the names below are the usual ones, written out rather than cited so
# that a reader can check the constant rather than take it.
#
# p is the number of free parameters in one segment: 1 for a change in mean,
# 1 for a change in variance, 2 for a change in both, 1 for Poisson. Adding a
# change point costs p new parameters plus the location itself, hence p + 1.

cp_penalties <- c("none", "aic", "bic", "sic", "hq", "mbic")

# Returns beta, and whether the cost has to be modified. MBIC is not a
# penalty in the sense the others are: see below.
cp_penalty <- function(name, cost, n = cost$n) {
  p <- cost$nparam
  if (is.numeric(name)) {
    if (length(name) != 1L || !is.finite(name) || name < 0)
      stop("cp_penalty: a numeric penalty must be one finite non-negative number")
    return(list(beta = as.numeric(name), mbic = FALSE, name = "manual"))
  }
  name <- tolower(as.character(name))
  beta <- switch(name,
    none = 0,
    aic  = 2 * (p + 1),
    bic  = (p + 1) * log(n),
    sic  = (p + 1) * log(n),          # the same thing under its other name
    hq   = 2 * (p + 1) * log(log(n)),
    # Zhang and Siegmund (2007). Their criterion, as a thing to minimise, is
    #
    #   -2 log L + 3k log n + sum_{i=1}^{k+1} log(l_i / n)
    #
    # with l_i the length of segment i. Expand the last sum:
    #
    #   sum log l_i - (k+1) log n
    #
    # so the whole thing is sum over segments of (C_i + log l_i), plus
    # k * (3 log n - log n), plus a constant -log n that no comparison can
    # see. That is exactly Optimal Partitioning with each segment's cost
    # raised by log(its length) and beta = 2 log n. It is segment-additive,
    # which is why it can be searched at all.
    mbic = 2 * log(n),
    stop("cp_penalty: unknown penalty '", name, "'"))
  if (name == "hq" && n < 16)
    stop("cp_penalty: Hannan-Quinn needs log(log(n)) > 0, so n must exceed e^1 = 3; ",
         "below n = 16 it is smaller than AIC and the name is misleading")
  list(beta = beta, mbic = identical(name, "mbic"), name = name)
}

# MBIC's log(segment length) term, as a wrapper round a cost object.
#
# This is where MBIC stops being interchangeable with the others. PELT's
# pruning is exact only if there is a constant K with
#
#     C(s,t) + C(t,T) + K <= C(s,T)    for all s < t < T
#
# and for a maximised log likelihood K = 0 works, because splitting a
# segment can only improve the fit. Adding log(length) destroys that:
# log(a) + log(b) against log(a + b) goes the wrong way as soon as
# a*b > a+b, which is almost always. The smallest value the left side can
# fall short by is log((a+b)/(a*b)) minimised at a = b = L/2 over L <= n,
# which is log(4/n). So MBIC needs K = log(4/n), a negative constant, and
# PELT run at K = 0 under MBIC is not exact.
#
# cp_pelt() picks this up from the cost object rather than being told, so
# that forgetting is not an option. t/test_mbic.R measures how often K = 0
# actually produces a different answer, because "not provably exact" and
# "wrong in practice" are different claims.
cp_mbic_cost <- function(cost) {
  base <- cost$seg
  out <- cost
  out$seg <- function(s, e) base(s, e) + log(e - s)
  out$pelt_K <- log(4 / cost$n)
  out$mbic <- TRUE
  out
}

# Beta and the (possibly wrapped) cost in one call, which is what every
# entry point above actually wants.
cp_prepare <- function(cost, penalty) {
  pen <- cp_penalty(penalty, cost)
  if (pen$mbic) cost <- cp_mbic_cost(cost)
  list(cost = cost, beta = pen$beta, name = pen$name)
}
