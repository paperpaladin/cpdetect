# The matcher, against answers worked out by hand.
#
# Every recovery and false positive rate this project reports comes out of
# cp_match, so it gets cases whose answers are obvious by inspection rather
# than computed by the thing under test. The interesting ones are the
# degenerate shapes: nothing detected, nothing true, two detections around
# one true change, one detection between two true changes.

source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "harness.R"))
cp_load()

chk <- function(true, found, w, tp, fn, fp, label) {
  m <- cp_match(true, found, w)
  ok(m$tp == tp && m$fn == fn && m$fp == fp, label,
     sprintf("got tp=%d fn=%d fp=%d, want tp=%d fn=%d fp=%d",
             m$tp, m$fn, m$fp, tp, fn, fp))
}

chk(c(50, 100), c(50, 100), 5, 2, 0, 0, "exact_match")
chk(c(50, 100), c(52, 97), 5, 2, 0, 0, "within_window")
chk(c(50, 100), c(56, 94), 5, 0, 2, 2, "outside_window_counts_both_ways")
chk(c(50, 100), integer(0), 5, 0, 2, 0, "nothing_detected")
chk(integer(0), c(50, 100), 5, 0, 0, 2, "nothing_true")
chk(integer(0), integer(0), 5, 0, 0, 0, "nothing_at_all")
chk(c(50), c(48, 52), 5, 1, 0, 1, "two_detections_one_change_is_one_false_positive")
chk(c(48, 52), c(50), 5, 1, 1, 0, "one_detection_two_changes_is_one_miss")
chk(c(50), c(50), 0, 1, 0, 0, "zero_window_accepts_an_exact_hit")
chk(c(50), c(51), 0, 0, 1, 1, "zero_window_rejects_one_off")
chk(c(10, 200), c(205, 12), 5, 2, 0, 0, "unsorted_input_is_sorted_first")

# The greedy matching's known limit, kept as a case rather than as a
# sentence. Two true change points 2 apart with one detection between
# them: greedy takes the nearer, which is the same answer an optimal
# assignment gives here. The case where they differ needs both detections
# and both true points inside one window, and is constructed below.
m <- cp_match(c(50, 52), c(51), 5)
ok(m$tp == 1 && m$fn == 1, "greedy_takes_the_nearer_of_two_crowded_truths",
   sprintf("matched %d, missed %s", m$tp, paste(m$missed, collapse = ",")))
m2 <- cp_match(c(50, 54), c(51, 55), 5)
ok(m2$tp == 2 && m2$fp == 0, "crowded_but_resolvable_matches_both",
   sprintf("pairs %s", paste(apply(m2$pairs, 1, paste, collapse = "->"), collapse = " ")))

# Displacement, which is the quality of a recovery as opposed to the fact
# of it.
m3 <- cp_match(c(50, 100, 150), c(53, 100, 146), 6)
ok(identical(sort(m3$disp), c(0L, 3L, 4L)), "displacement_is_reported",
   sprintf("%s", paste(sort(m3$disp), collapse = ",")))
ok(is.integer(cp_match(integer(0), integer(0), 3)$disp) &&
   length(cp_match(integer(0), integer(0), 3)$disp) == 0L,
   "no_pairs_means_no_displacements")

# Aggregation. Hand-built results so the arithmetic is checkable.
res <- list(list(true = c(50, 100), found = c(51, 99)),      # 2 tp
            list(true = c(50, 100), found = c(51)),          # 1 tp 1 fn
            list(true = integer(0), found = c(30, 70)),      # 2 fp
            list(true = c(50), found = integer(0)))          # 1 fn
sc <- cp_score(res, 5)
ok(sc$tp == 3 && sc$fn == 2 && sc$fp == 2, "score_totals",
   sprintf("tp %d fn %d fp %d", sc$tp, sc$fn, sc$fp))
ok(abs(sc$recovery - 3 / 5) < 1e-12, "pooled_recovery_is_tp_over_true",
   sprintf("%.4f against 0.6", sc$recovery))
ok(abs(sc$precision - 3 / 5) < 1e-12, "precision_is_tp_over_found",
   sprintf("%.4f against 0.6", sc$precision))
ok(abs(sc$fp_per_series - 0.5) < 1e-12, "false_positives_per_series",
   sprintf("%.4f against 0.5", sc$fp_per_series))
ok(abs(sc$any_fp_rate - 0.25) < 1e-12, "fraction_of_series_with_any_false_positive",
   sprintf("%.4f against 0.25", sc$any_fp_rate))
ok(is.na(sc$recovery_per_series[3]), "a_null_series_has_no_recovery_rate",
   "NA rather than 0 or 1, because there is nothing to recover")

ok(tryCatch({ cp_match(1, 1, -1); FALSE }, error = function(e) TRUE),
   "match_refuses_a_negative_window")
ok(tryCatch({ cp_match(1, 1, c(1, 2)); FALSE }, error = function(e) TRUE),
   "match_refuses_a_vector_window")

done()
