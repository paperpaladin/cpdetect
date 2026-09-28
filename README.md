# cpdetect

Offline change point detection for security and operational telemetry,
written from the papers rather than from a package. Given a series of
numbers, it finds the instants where the generating distribution changed:
the mean, the variance, both at once, or the rate of a count process such
as failed logins per minute.

Base R and MASS only. No `changepoint`, no `strucchange`, no CRAN at all.
Writing the dynamic programme is the point of the exercise.

```sh
cpdetect [OPTIONS] FILE            one number per line; - reads stdin
  --cost mean|var|meanvar|poisson   which parameter may change (mean)
  --method op|pelt|binseg|brute     exact DP, pruned exact DP, greedy,
                                    exhaustive (pelt)
  --penalty bic|mbic|aic|hq|none|N  per-changepoint penalty (bic)
  --minseg N                        shortest admissible segment (1)
  --maxcp K                         cap on change points, brute only
  --sigma2 S                        known variance, --cost mean only (1)
  --mu M                            known mean, --cost var only
  --perm B                          permutation test with B resamples
  --seed S                          seed for the permutation draws
  --quiet                           print only the change point indices
```

`make test` runs the suite. It needs `Rscript` and `python3` with scipy,
and refuses to run without either rather than reporting a smaller number
having checked less.

## What it is checked against

Four oracles, and the first of them is not another implementation.

**Brute force.** For a series short enough, every possible segmentation can
be enumerated and the true cost minimum found by definition. That is not a
second opinion, it is the answer. Optimal partitioning agrees with it on
all 67 cases, worst absolute cost difference 2.84e-14 against a tolerance
of 1e-9, and the change point vectors are identical rather than merely
equally priced.

**PELT against optimal partitioning.** PELT's pruning is supposed to be
exact, not approximate, so the two must return the same segmentation and
not just a similar one. They do, on all 48 cases, with 0 equal-cost ties,
which matters because a tie would let them agree by luck.

The pruning is shown to be doing something: 8,013 candidates discarded
across 44 of the 48 cases, and at n = 1000 the three methods evaluate
3,892 (binseg), 117,379 (pelt) and 496,506 (op) candidates. A PELT whose
pruning never fired would be optimal partitioning wearing a different
name, and the counts are asserted rather than assumed.

Setting the pruning constant K to zero, which is the classic way to get
this wrong, makes 2 of 48 runs return a cost above the exact minimum, the
worst by 0.7881. So the constant is load-bearing and the suite can tell.

**Python and scipy, as a different language.** `t/oracle.py` shares nothing
with `R/cost.R` but the definitions: it computes segment costs from
`scipy.stats.norm.logpdf` and `scipy.stats.poisson.logpmf` in two passes
rather than from prefix sums, enumerates with `itertools.combinations`
rather than by dynamic programming, and accumulates with `math.fsum` where
summation order could matter. Agreement is then about the definitions
rather than about one person's arithmetic. Worst cost difference 7.54e-12
over 18 cases and 7,880 segmentations scored, identical segmentations, and
the enumeration count matches the closed form for the number of ways to
place k change points in n points.

**Known-answer series.** Data generated with change points placed
deliberately, so recovery can be measured instead of asserted.

## Recovery and false positives, separately

Both, always. A detector that finds every real change and twenty spurious
ones is not a good detector, and a single accuracy number hides that. The
evaluation reports pooled recovery (true positives over true changes),
false positives per series, and the fraction of series with any false
positive at all. A null series with no change points reports NA rather
than 0 or 1, because there is nothing there to recover.

## What does not work

The permutation test is implemented and its draws are seeded, but no claim
is made here about its calibration under the null. Measuring the actual
false positive rate of the test at a stated alpha needs more series than
the suite runs.

`binseg` is greedy and is not exact. It is in the tool because it is what
most people reach for and the evaluation counts show what it costs you,
not because its answers are trusted. Where this README quotes an exact
result it is `op` or `pelt`.

The costs assume independent observations. Telemetry is usually
autocorrelated, which inflates the number of change points found, and
nothing here corrects for it.
