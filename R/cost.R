# Segment cost functions.
#
# A cost function here is minus twice the maximised log likelihood of one
# segment under its model. Minus twice, rather than minus once, so that the
# penalties below are the textbook BIC and AIC constants rather than those
# constants halved; and the full log likelihood including the terms that do
# not depend on the parameters, rather than the usual residual sum of
# squares, for one specific reason: the constants are what makes the number
# comparable against scipy.stats.norm.logpdf in another language. A cost that
# has had its additive constants dropped can only be checked against another
# copy of the same convention, which is to say against ourselves.
#
# Every cost is O(1) in the segment length, from prefix sums laid down once.
# That is what makes Optimal Partitioning O(n^2) rather than O(n^3), and it
# is also where catastrophic cancellation lives: sum(y^2) - sum(y)^2/m is the
# difference of two large numbers when the data has a large offset.
# t/test_cancel.R measures how wrong it gets rather than assuming it is fine.
#
# ------------------------------------------------------------------ indexing
#
# A segmentation of y[1..n] is a set of boundaries
#
#     0 = tau_0 < tau_1 < ... < tau_k < tau_{k+1} = n
#
# and segment i is y[(tau_{i-1} + 1) .. tau_i]. So a change point is reported
# as the index of the LAST observation before the change. seg(s, e) is the
# cost of y[(s+1)..e]; s is exclusive, e is inclusive, both 0-based in the
# sense that s = 0 means "from the start". Prefix vectors are stored with a
# leading zero so that the sum over (s+1)..e is pre[e + 1] - pre[s + 1] with
# no special case at s = 0.
#
# seg() is vectorised over s with e scalar, because that is the shape of the
# inner loop in both Optimal Partitioning and PELT, and an R loop over
# candidate start points costs more than the arithmetic it performs.

# The smallest variance any segment is allowed to be fitted with.
#
# Without a floor, a segment whose observations are all equal has maximum
# likelihood variance 0 and cost -Inf, so the optimiser's answer to any
# series containing two equal neighbours is "a change point on either side of
# them", which is both the global minimum and nonsense. The floor is a
# constraint on the parameter space, and every variance cost below is
# evaluated as the true -2 log likelihood AT the constrained maximum rather
# than as the unconstrained formula with a clamped argument. That distinction
# is not cosmetic: it is what keeps the cost superadditive, which is the
# precondition PELT's pruning rests on. See R/pelt.R and t/test_pelt.R.
# The absolute part of the floor. A second, relative part is added below.
CP_VARFLOOR <- 1e-8

# The relative part, in units of the mean square of the (centred) data.
#
# An absolute floor alone is not enough, and the reason is the same
# amplification: the residual sum of squares computed from prefix sums
# carries an absolute error of order eps * sum(y^2), so on a series of bytes
# per second with values around 10^6 that error is around 10^-4. Dividing it
# by a floor of 10^-8 puts 10^4 into a cost, which is larger than any real
# term in it. The floor therefore has to sit above the noise in the quantity
# it is flooring. 1024 * eps = 2.3e-13 of the mean square is about three
# orders of magnitude above that noise and, for data of order 1, is far
# below the 10^-8 absolute floor and so changes nothing.
#
# This is a modelling choice and it is visible in the answers: it says a
# segment may not be fitted with a variance smaller than the resolution at
# which a variance can be computed at all. t/test_cancel.R measures what it
# is worth.
CP_VARFLOOR_REL <- 1024 * .Machine$double.eps

cp_types <- c("mean", "var", "meanvar", "poisson")

# Free parameters per segment, which is what a change point buys. Used by the
# penalties in R/penalty.R and nowhere else.
cp_nparam <- c(mean = 1L, var = 1L, meanvar = 2L, poisson = 1L)

# center: build the Gaussian prefix sums from y - mean(y) rather than from y.
# The residual sum of squares is invariant to a shift in exact arithmetic and
# is not in doubles: sum(y^2) - sum(y)^2/m with y around 10^4 loses eight
# significant digits before the subtraction happens. Default TRUE. The
# argument exists so that t/test_cancel.R can measure what it buys instead of
# the README asserting that it buys something.
cp_cost <- function(y, type = "mean", sigma2 = 1, mu = NULL,
                    varfloor = NULL, center = TRUE) {
  if (!type %in% cp_types)
    stop("cp_cost: unknown cost type '", type, "'")
  y <- as.numeric(y)
  n <- length(y)
  if (n < 1L) stop("cp_cost: empty series")
  if (any(!is.finite(y))) stop("cp_cost: series contains non-finite values")
  if (!(is.numeric(sigma2) && length(sigma2) == 1L && is.finite(sigma2) && sigma2 > 0))
    stop("cp_cost: sigma2 must be a single positive finite number")

  # Poisson keeps the raw counts: there are no squares in it, so there is no
  # cancellation to avoid, and shifting the data would make it non-integer.
  shift <- if (center && type != "poisson") mean(y) else 0
  yc <- y - shift
  z <- c(0, cumsum(yc))                      # prefix sum of the shifted data
  zz <- c(0, cumsum(yc * yc))                # prefix sum of its squares
  log2pi <- log(2 * pi)
  if (is.null(varfloor))
    varfloor <- max(CP_VARFLOOR, CP_VARFLOOR_REL * (if (n) sum(yc * yc) / n else 0))
  if (!(is.numeric(varfloor) && length(varfloor) == 1L && is.finite(varfloor) && varfloor > 0))
    stop("cp_cost: varfloor must be a single positive finite number")

  self <- list(n = n, y = y, type = type, sigma2 = sigma2,
               varfloor = varfloor, center = center, shift = shift,
               nparam = unname(cp_nparam[[type]]))

  # Residual sum of squares of y[(s+1)..e] about its own mean, from the
  # prefix sums, with the one case that has an exact answer taken out of the
  # arithmetic entirely.
  #
  # A single observation has residual sum of squares exactly zero. The prefix
  # sums do not know that: zz[e+1] - zz[s+1] is the difference of two running
  # totals and comes back as the true y^2 plus a few units in the last place
  # of the TOTAL, which for y[33] = 0.0409 in a series of 40 was 5.4e-15
  # instead of 0. Harmless on its own, except that the variance floor then
  # divides it by 1e-8 and the cost is wrong in the seventh decimal place.
  # That is the whole story of this three-line function: the floor is a 10^8
  # amplifier sitting directly downstream of a cancelling subtraction.
  # The second half is the fallback. Where the prefix difference comes back
  # at or below the size of its own rounding error, its digits are noise, and
  # the variance floor downstream divides it by 1e-8 or smaller. Those few
  # segments are recomputed in two passes over the data, which is exact
  # enough and is reached only by segments that are numerically constant.
  #
  # shortcut: the fallback is O(m) per segment, so a series that is constant
  # over long stretches (bytes-out from an idle host, a counter that does not
  # move) can drive the dynamic program from O(n^2) to O(n^3). An exactly
  # constant run does not trigger it, because the prefix difference is then
  # exactly zero and zero is not above the noise threshold; a run that is
  # constant after centring but not before is the case that pays. The upgrade
  # path is a run-length prescan that answers those segments in O(1).
  rss_noise <- 8 * .Machine$double.eps * max(abs(zz))
  rss_of <- function(s, e) {
    m <- e - s
    S <- z[e + 1] - z[s + 1]
    r <- ifelse(m == 1, 0, pmax((zz[e + 1] - zz[s + 1]) - S * S / m, 0))
    if (rss_noise > 0) {
      bad <- which(r > 0 & r <= rss_noise)
      if (length(bad)) {
        sv <- rep_len(s, length(r))
        r[bad] <- vapply(bad, function(k) {
          v <- yc[(sv[k] + 1):e]
          sum((v - mean(v)) ^ 2)
        }, numeric(1))
      }
    }
    r
  }

  if (type == "mean") {
    # Known (or assumed) variance sigma2; only the mean may change.
    # -2 log L = m log(2 pi sigma2) + sum (y - ybar)^2 / sigma2
    self$seg <- function(s, e) {
      (e - s) * (log2pi + log(sigma2)) + rss_of(s, e) / sigma2
    }
    self$fit <- function(s, e) c(mean = shift + (z[e + 1] - z[s + 1]) / (e - s))

  } else if (type == "var") {
    # Known mean, variance may change. The mean defaults to the mean of the
    # whole series, which is the usual choice and is also a choice: a series
    # whose mean changes will be read as a variance change if you tell this
    # cost the mean is constant. Stated in the README under what does not work.
    mu0 <- if (is.null(mu)) mean(y) else as.numeric(mu)
    if (length(mu0) != 1L || !is.finite(mu0)) stop("cp_cost: mu must be one finite number")
    self$mu <- mu0
    d2 <- (y - mu0) ^ 2
    q <- c(0, cumsum(d2))
    # Same exact-at-m-equals-1 case as rss_of, for the same reason: a one
    # observation segment has sum of squares (y - mu)^2 exactly, and taking
    # it out of the prefix difference keeps the floor from amplifying the
    # cancellation residue.
    q_noise <- 8 * .Machine$double.eps * max(abs(q))
    q_of <- function(s, e) {
      m <- e - s
      r <- ifelse(m == 1, d2[e], pmax(q[e + 1] - q[s + 1], 0))
      if (q_noise > 0) {
        bad <- which(r > 0 & r <= q_noise)
        if (length(bad)) {
          sv <- rep_len(s, length(r))
          r[bad] <- vapply(bad, function(k) sum(d2[(sv[k] + 1):e]), numeric(1))
        }
      }
      r
    }
    self$seg <- function(s, e) {
      m <- e - s
      Q <- q_of(s, e)
      s2 <- pmax(Q / m, varfloor)
      m * (log2pi + log(s2)) + Q / s2
    }
    self$fit <- function(s, e) {
      m <- e - s
      c(mu = mu0, var = pmax(q_of(s, e) / m, varfloor))
    }

  } else if (type == "meanvar") {
    self$seg <- function(s, e) {
      m <- e - s
      rss <- rss_of(s, e)
      s2 <- pmax(rss / m, varfloor)
      m * (log2pi + log(s2)) + rss / s2
    }
    self$fit <- function(s, e) {
      m <- e - s
      rss <- rss_of(s, e)
      c(mean = shift + (z[e + 1] - z[s + 1]) / m, var = pmax(rss / m, varfloor))
    }

  } else {
    # Poisson. For count telemetry: failed logins per minute, exec events per
    # host per interval. The lgamma prefix is the log(y!) term of the pmf; it
    # cancels out of every comparison between segmentations of the same
    # series, and it is kept because without it the number is not a log
    # likelihood and scipy cannot be asked whether it is right.
    if (any(y < 0)) stop("cp_cost: the Poisson cost needs non-negative counts")
    if (any(y != floor(y))) stop("cp_cost: the Poisson cost needs integer counts")
    g <- c(0, cumsum(lgamma(y + 1)))
    self$seg <- function(s, e) {
      m <- e - s
      S <- z[e + 1] - z[s + 1]
      G <- g[e + 1] - g[s + 1]
      # 0 log 0 = 0 at the limit. An all-zero segment has lambda_hat 0 and
      # log likelihood exactly 0, which is right and is not what log(0 / m)
      # evaluates to.
      Slog <- ifelse(S > 0, S * log(S / m), 0)
      -2 * (Slog - S) + 2 * G
    }
    self$fit <- function(s, e) c(lambda = (z[e + 1] - z[s + 1]) / (e - s))
  }

  self
}

# Total cost of a segmentation, given as interior change points: the index of
# the last observation of every segment but the last. A free function rather
# than a member, because R/penalty.R wraps a cost object to add MBIC's
# log(segment length) term and a member closure would still be calling the
# unwrapped seg().
cp_total <- function(cost, taus) {
  b <- c(0L, as.integer(taus), cost$n)
  if (any(diff(b) < 1L)) stop("cp_total: change points must be strictly increasing within 1..n-1")
  sum(vapply(seq_len(length(b) - 1L),
             function(i) cost$seg(b[i], b[i + 1L]), numeric(1)))
}

# Vectorised over both arguments, for the oracle dumps where s and e are
# arbitrary pairs rather than a column of a DP table. Kept separate because
# the DP inner loop only ever needs the scalar-e form and paying for the
# general case there would be the most-executed line in the program.
cp_seg_many <- function(cost, s, e) {
  stopifnot(length(s) == length(e))
  vapply(seq_along(s), function(i) cost$seg(s[i], e[i]), numeric(1))
}
