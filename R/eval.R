# Scoring a detection against the truth.
#
# Two numbers, always both. A detector that finds every real change and
# twenty spurious ones is not a good detector, and a single accuracy figure
# hides which of the two it is. So:
#
#   recovery rate   recovered true change points / true change points
#   false positives detected change points matched to no true one
#
# A detected change point almost never lands on the exact index, because
# the maximum likelihood boundary of a segment is a random variable even
# when the change is real and large. So matching is within a tolerance
# window w, and w is an argument with no default: it is the single choice
# that moves these numbers most, and burying it in a default would make
# every rate here a statement about a number nobody saw.
#
# The matching is one-to-one and greedy by increasing distance. One-to-one
# matters: without it, two detections either side of one true change count
# as two recoveries and zero false positives, which flatters a detector
# that has split a change in half. Greedy-by-distance is not guaranteed
# optimal in the assignment sense, and the alternative is the Hungarian
# algorithm over a problem with at most a handful of points on each side.
# shortcut: greedy matching. It differs from optimal only when two true
# change points are within w of each other, which the generators here
# avoid by construction (gap is always larger than w in every measurement)
# and which t/test_eval.R has a case for. The upgrade path is a full
# assignment, about thirty lines.

cp_match <- function(true, found, w) {
  true <- sort(as.integer(true))
  found <- sort(as.integer(found))
  if (!(is.numeric(w) && length(w) == 1L && is.finite(w) && w >= 0))
    stop("cp_match: w must be one finite non-negative number")
  pairs <- matrix(integer(0), ncol = 2)
  tfree <- rep(TRUE, length(true))
  ffree <- rep(TRUE, length(found))
  if (length(true) && length(found)) {
    d <- outer(true, found, function(a, b) abs(a - b))
    repeat {
      cand <- d
      cand[!tfree, ] <- NA
      cand[, !ffree] <- NA
      cand[cand > w] <- NA
      if (all(is.na(cand))) break
      i <- which(cand == min(cand, na.rm = TRUE), arr.ind = TRUE)[1, ]
      pairs <- rbind(pairs, c(true[i[1]], found[i[2]]))
      tfree[i[1]] <- FALSE
      ffree[i[2]] <- FALSE
    }
  }
  list(tp = nrow(pairs), fn = sum(tfree), fp = sum(ffree),
       pairs = pairs, missed = true[tfree], spurious = found[ffree],
       w = w,
       # Absolute displacement of the matched pairs: how far off a
       # recovered change point is. Zero true change points makes this
       # NA rather than 0, because "perfectly located" and "nothing to
       # locate" are different facts.
       disp = if (nrow(pairs)) abs(pairs[, 1] - pairs[, 2]) else integer(0))
}

# Aggregate over a list of (true, found) pairs. Returns the two rates, the
# per-series counts, and the spread, because a single run is noise and a
# mean with no spread beside it is a single run with extra steps.
cp_score <- function(results, w) {
  m <- lapply(results, function(r) cp_match(r$true, r$found, w))
  tp <- vapply(m, function(x) x$tp, integer(1))
  fn <- vapply(m, function(x) x$fn, integer(1))
  fp <- vapply(m, function(x) x$fp, integer(1))
  ntrue <- tp + fn
  nfound <- tp + fp
  disp <- unlist(lapply(m, function(x) x$disp))
  list(
    nseries = length(m),
    # Pooled over all series: the number most people mean by "recall".
    recovery = if (sum(ntrue)) sum(tp) / sum(ntrue) else NA_real_,
    # Per-series recall, so that a spread can be quoted. NA for series
    # with no true change point, which are scored by false positives only.
    recovery_per_series = ifelse(ntrue > 0, tp / pmax(ntrue, 1), NA_real_),
    precision = if (sum(nfound)) sum(tp) / sum(nfound) else NA_real_,
    fp_total = sum(fp),
    fp_per_series = mean(fp),
    fp_sd = stats::sd(fp),
    # The number that matters on a null series: did it say anything at all.
    any_fp_rate = mean(fp > 0),
    disp_mean = if (length(disp)) mean(disp) else NA_real_,
    disp_max = if (length(disp)) max(disp) else NA_integer_,
    tp = sum(tp), fn = sum(fn), fp = sum(fp), w = w)
}
