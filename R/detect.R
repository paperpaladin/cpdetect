# The one entry point everything above the algorithms goes through.
#
# The command line in ./cpdetect, the recovery measurements in
# t/measure_recovery.R and the in-suite checks all want the same four
# steps: build the cost, resolve the penalty, wrap the cost if the penalty
# is MBIC, run the search. Doing that in three places is three chances to
# resolve the penalty against the wrong cost object, which is a defect that
# produces plausible answers.

cp_methods <- c("op", "pelt", "binseg", "brute")

cp_detect <- function(y, type = "mean", method = "pelt", penalty = "bic",
                      minseg = 1L, maxcp = Inf, sigma2 = 1, mu = NULL,
                      varfloor = NULL) {
  if (!method %in% cp_methods)
    stop("cp_detect: unknown method '", method, "'")
  cost <- cp_cost(y, type, sigma2 = sigma2, mu = mu, varfloor = varfloor)
  pr <- cp_prepare(cost, penalty)
  fit <- switch(method,
    op = cp_op(pr$cost, pr$beta, minseg),
    pelt = cp_pelt(pr$cost, pr$beta, minseg),
    binseg = cp_binseg(pr$cost, pr$beta, minseg),
    brute = cp_brute(pr$cost, pr$beta, minseg, maxcp))
  fit$penalty <- pr$name
  fit$type <- type
  fit$cost_obj <- pr$cost
  fit$plain_cost <- cost
  fit$n <- cost$n
  fit
}

# Fitted parameters per segment, for printing and for the tests that check
# a detector recovered the right rates rather than only the right indices.
cp_segments <- function(fit) {
  b <- c(0L, as.integer(fit$cp), fit$n)
  do.call(rbind, lapply(seq_len(length(b) - 1L), function(i) {
    p <- fit$plain_cost$fit(b[i], b[i + 1L])
    data.frame(start = b[i] + 1L, end = b[i + 1L], n = b[i + 1L] - b[i],
               param = paste(names(p), signif(p, 7), sep = "=", collapse = " "),
               cost = fit$plain_cost$seg(b[i], b[i + 1L]),
               stringsAsFactors = FALSE)
  }))
}
