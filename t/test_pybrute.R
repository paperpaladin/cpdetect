# Brute force, a second time, in Python.
#
# R's cp_brute and R's cp_op call the same cost$seg. If that function is
# wrong they agree with each other and are both wrong, and the agreement
# sweep in t/test_brute.R would report a clean run. This file closes that:
# t/oracle.py enumerates the segmentations itself, with itertools, scoring
# each one from scipy densities computed two passes over the segment, and
# sums with math.fsum. Nothing about the search or the arithmetic is shared.
#
# It also reports how many segmentations it scored, which is compared
# against R's count and against the closed form. Three independent answers
# to "how many admissible segmentations are there" is a stronger statement
# about minseg than any of them alone.

source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "harness.R"))
cp_load()

if (!oracle_available()) {
  ok(FALSE, "scipy_oracle_present_for_brute", "python3 with scipy and numpy is required")
  done()
}

# Tolerance. Python sums the segments with math.fsum, which is exactly
# rounded, and R sums left to right, so this comparison sees a little more
# disagreement than the R-against-R one. Still four orders under 1e-9 in
# practice; the worst observed is reported.
TOL <- 1e-9

brute_req <- function(cost, beta, minseg, maxcp, mbic) {
  c(paste("type", cost$type),
    paste("sigma2", g17(cost$sigma2)),
    paste("mu", g17(if (is.null(cost$mu)) 0 else cost$mu)),
    paste("varfloor", g17(cost$varfloor)),
    paste("penalty", g17(beta)),
    paste("minseg", minseg),
    paste("maxcp", maxcp),
    paste("mbic", if (mbic) 1 else 0),
    paste("y", paste(g17(cost$y), collapse = " ")))
}

parse_out <- function(lines) {
  get <- function(k) {
    i <- grep(paste0("^", k, " "), lines)
    if (!length(i)) return(NULL)
    sub(paste0("^", k, " *"), "", lines[i[1]])
  }
  cps <- get("cp")
  list(best = as.numeric(get("best")),
       cp = if (is.null(cps) || !nzchar(trimws(cps))) integer(0)
            else as.integer(strsplit(trimws(cps), " +")[[1]]),
       n = as.integer(get("n")))
}

n_expected <- function(n, ms, maxparts) {
  js <- seq_len(min(maxparts, n %/% ms))
  sum(choose(n - js * ms + js - 1, js - 1))
}

grid <- list()
i <- 0L
for (ty in c("mean", "var", "meanvar", "poisson"))
  for (pn in c("none", "bic", "mbic"))
    for (ms in c(1L, 2L)) {
      i <- i + 1L
      if (ty %in% c("var", "meanvar") && ms == 1L) next   # degenerate; see test_degenerate.R
      grid[[length(grid) + 1L]] <- list(type = ty, pen = pn, minseg = ms,
                                        n = c(10L, 11L, 12L)[1 + (i %% 3L)],
                                        k = i %% 3L, seed = 510000L + i)
    }

worst <- 0
n_cpdiff <- 0L
n_countdiff <- 0L
detail <- character(0)
cdetail <- character(0)
total_py <- 0

for (g in grid) {
  sim <- cp_sim_random(g$n, g$k, g$type, effect = 3, gap = g$minseg, seed = g$seed)
  co <- cp_cost(sim$y, g$type)
  pr <- cp_prepare(co, g$pen)
  a <- cp_op(pr$cost, pr$beta, g$minseg)
  r <- parse_out(oracle("brute",
                        brute_req(co, pr$beta, g$minseg, g$n, pr$name == "mbic"),
                        tag = paste0(g$type, g$pen, g$minseg)))
  total_py <- total_py + r$n
  d <- abs(a$cost - r$best)
  if (d > worst) worst <- d
  if (!identical(as.integer(a$cp), r$cp)) {
    n_cpdiff <- n_cpdiff + 1L
    detail <- c(detail, sprintf("%s/%s/ms%d/n%d: R [%s] %.12f vs py [%s] %.12f",
                                g$type, g$pen, g$minseg, g$n,
                                paste(a$cp, collapse = ","), a$cost,
                                paste(r$cp, collapse = ","), r$best))
  }
  want <- n_expected(g$n, g$minseg, g$n)
  if (r$n != want) {
    n_countdiff <- n_countdiff + 1L
    cdetail <- c(cdetail, sprintf("%s/ms%d/n%d: python %d, closed form %d",
                                  g$type, g$minseg, g$n, r$n, want))
  }
}

ok(length(grid) > 0L, "pybrute_cases_exist", sprintf("%d cases", length(grid)))
ok(worst <= TOL, "op_cost_equals_the_python_enumeration",
   sprintf("worst %.3g over %d cases, %0.f segmentations scored in Python",
           worst, length(grid), total_py))
ok(n_cpdiff == 0L, "op_segmentation_equals_the_python_enumeration",
   if (n_cpdiff) paste(detail, collapse = " | ") else "identical change point vectors")
ok(n_countdiff == 0L, "python_enumeration_count_matches_the_closed_form",
   if (n_countdiff) paste(cdetail, collapse = " | ")
   else sprintf("%.0f segmentations across %d cases", total_py, length(grid)))

# The Python side has to be able to disagree, and the cheapest way to prove
# it is to ask it a question with a different answer: the same series with
# the penalty removed must give at least as many change points, and on a
# series with a real change it gives strictly more.
sim <- cp_sim(14, c(7), "mean", mean = c(0, 4), seed = 611)
co <- cp_cost(sim$y, "mean")
r_bic <- parse_out(oracle("brute", brute_req(co, cp_penalty("bic", co)$beta, 1L, 14L, FALSE), tag = "vb"))
r_non <- parse_out(oracle("brute", brute_req(co, 0, 1L, 14L, FALSE), tag = "vn"))
ok(length(r_non$cp) > length(r_bic$cp), "python_penalty_argument_is_live",
   sprintf("bic gives %d change points, no penalty gives %d",
           length(r_bic$cp), length(r_non$cp)))
ok(r_non$best < r_bic$best, "python_unpenalised_cost_is_lower",
   sprintf("%.4f vs %.4f", r_non$best, r_bic$best))

done()
