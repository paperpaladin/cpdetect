# PLAN

## Contract

```
cpdetect [OPTIONS] FILE            one number per line; "-" reads stdin
  --cost mean|var|meanvar|poisson  which parameter is allowed to change
  --method op|pelt|binseg|brute    exact DP, pruned exact DP, greedy, exhaustive
  --penalty bic|mbic|aic|hq|NUMBER per-changepoint penalty beta
  --minseg N                       shortest admissible segment (default 1)
  --maxcp K                        cap on changepoints (brute force only)
  --sigma2 S                       known variance, --cost mean only
  --perm B                         permutation test with B resamples
  --seed S                         seed for anything random
  --quiet                          print only the changepoint indices
```

Output is a header line of the settings, one `cp` line listing the detected
change points as 1-based indices of the LAST observation of each segment, one
`segment` line per segment with its fitted parameters and its own cost, and a
`total` line with the penalised objective. `pelt` additionally prints
`pruned=N/M`, the number of candidate start points removed by the pruning rule
out of the number offered to it.

`--perm B` adds a `perm` line: the observed maximum cost reduction over all
admissible single splits, the number of resamples whose statistic matched or
exceeded it, and the resulting p-value `(1 + hits) / (B + 1)`.

Pure R: base R plus `stats`, with `MASS` available and currently unused. No
CRAN package. `python3` with `scipy` and `numpy` is used by the test harness
and by nothing else.

Exit 0 on a clean run, 1 on malformed or unreadable input, 2 on a usage error.
The three are distinct because a refusal and a wrong answer are different
events and the suite asserts on which one happened.

What must not break: the segmentation `op` returns is a global minimiser of
the penalised cost, `pelt` returns the same segmentation as `op` on every
input, and neither ever reports a change point that its own segment list
does not imply.

## Hard invariant

**On every series short enough to enumerate, the dynamic program's answer is
the exhaustive minimum; and PELT's pruned search returns exactly what
unpruned Optimal Partitioning returns, with the pruning demonstrably having
fired.**

Concretely, and checked by `make test`:

- For every (series, cost, penalty, minseg) case with n <= 25, the penalised
  cost found by `op` equals the minimum over all 2^(n-1) admissible
  segmentations found by brute-force enumeration, to a stated tolerance, and
  no enumerated segmentation is cheaper. Brute force is not a second opinion,
  it is the definition.
- A second brute force, written in Python against numpy and scipy, agrees
  with the R brute force on the same cases. Two languages, two authors'
  worth of arithmetic.
- `pelt` and `op` return identical change point vectors on every series in
  the suite, and the run asserts the pruning counter is nonzero, because a
  PELT whose pruning never fires is Optimal Partitioning under another name
  and proves nothing about the pruning rule.
- Every segment cost agrees with the same quantity computed in Python from
  `scipy.stats.norm.logpdf` and `scipy.stats.poisson.logpmf`, to a stated
  tolerance. R's own `stats` is used too, but it is the weaker check: it is
  arithmetic we wrote checked against arithmetic in the same runtime.
- Recovery rate and false-positive rate are both measured over seeded
  known-answer corpora and both are reported. Neither alone is a result.
- Every number in the README is produced by a target in the Makefile.

## Oracle

Ranked, by how little of it we control:

1. **Exhaustive enumeration.** For n <= 25 every admissible segmentation can
   be scored and the minimum taken. This is not another implementation of
   the answer, it is the answer, and it is the only oracle in this list that
   cannot itself be wrong about what the optimum is (it can only be wrong
   about the cost function, which is why the cost function has its own
   oracle below).
2. **Python, numpy and scipy 1.15.3.** An independent implementation of the
   four cost functions, of the test statistics, and of brute-force
   enumeration, in a different language with a different floating point
   evaluation order. A cost computed two ways in two languages catches a
   numerics error that no amount of R-only testing can see.
3. **The algorithms against each other.** PELT's pruning is exact by
   construction, so PELT must equal Optimal Partitioning on every input, and
   the claim is only worth something alongside a count of what the pruning
   removed. Binary segmentation is *not* exact, so the measured rate at
   which it disagrees is a result rather than a failure.
4. **Known-answer corpora.** Series generated with change points placed by
   us, where recovery and false positives are both countable. Seeded, with
   the seed asserted to be load-bearing.
5. **R's `stats`.** `dnorm`, `dpois`, `var`, `ks.test`. Deliberately last:
   checking R arithmetic with R arithmetic in the same process shares every
   library and every rounding mode, and the one thing it cannot find is the
   thing we are most likely to have got wrong.
6. **Deliberate defects in cpdetect itself.** Each mutant must turn a named
   assertion red. Survivors recorded in the README, equivalent mutants
   recorded as equivalent and argued for.

## Non-goals

- Online or sequential detection. Everything here sees the whole series
  before answering. CUSUM-style streaming detection is a different contract
  with a different failure mode (detection delay) and no oracle here.
- Multivariate series. One cost function per series, one scalar observation
  per index.
- Dependence. Every cost function here assumes observations are independent
  given the segment. Authentication-failure counts are autocorrelated in
  real life, which inflates the change point count; this is measured and
  stated rather than corrected for.
- Being a package. No `DESCRIPTION`, no namespace, no `R CMD check`.
- Reimplementing `changepoint` or `strucchange`. Neither is installed and
  neither is going to be. Writing the algorithm from the paper is the point;
  agreeing with somebody else's R code would prove only that two people read
  the same package.
- Nonparametric costs (`cpt.np`, empirical distribution). The parametric
  four are enough to exercise the machinery.
- Choosing the penalty for you. CROPS and the elbow plot are not here.

## Commit sequence

### Milestone 0: a harness that can fail

- [x] Repository skeleton, Makefile, and `t/run.sh` with a self-check
      why: the harness is the thing most likely to be silently broken, so it
      is written first and watched failing first. Every test script must
      report how many assertions it intended to run, and `run.sh` must refuse
      a script whose count does not match what it saw, because a script that
      dies halfway through otherwise just reports a smaller number.

### Milestone 1: the cost functions, checked across a language boundary

- [x] Cumulative sufficient statistics and the change-in-mean cost
      why: every later layer calls `seg_cost(s, e)` tens of thousands of
      times, so it has to be O(1) from prefix sums, and prefix sums are
      exactly where catastrophic cancellation lives.
- [x] Cross-language check of the mean cost against `scipy.stats.norm.logpdf`
      why: `SS - S^2/m` and `sum((y - mean(y))^2)` are the same number in
      exact arithmetic and not in doubles. This is the check that can see it.
- [x] Change in variance and change in mean-and-variance costs
      why: the variance costs are where the degenerate segment lives, and a
      cost function that returns `-Inf` for a constant segment makes every
      optimiser above it choose nonsense.
- [x] The Poisson cost, for count telemetry
      why: failed logins per minute are counts, not Gaussians, and the
      normal costs on small counts give the wrong segmentation. The
      `0 * log 0` limit has to be handled explicitly.
- [x] Measure the cancellation error of the prefix-sum costs
      why: claiming the O(1) trick is safe is not the same as knowing how
      wrong it gets. Compare against a two-pass per-segment computation on
      series with a large offset and report the worst observed error.
      done: two separate defects, both measured. Centring the prefix sums
      takes the worst segment cost error at offset 1e6 from 6.8e-2 to
      5.4e-12. And the variance floor is a 1e8 amplifier sitting directly
      downstream of a cancelling subtraction: a one-observation segment
      whose residual sum of squares is exactly 0 came out of the prefix
      difference as 5.4e-15, and 5.4e-15 / 1e-8 is an error in the seventh
      decimal place of the cost. Fixed by answering the exact cases exactly
      and by falling back to a two-pass computation for any segment whose
      prefix difference is at or below its own rounding error, which took
      the worst case on a large-scale series from 2.0e+1 to 2.3e-13.
      The Poisson cost at lambda = 5e5 keeps an absolute error of 1.3e-7
      and always will: its three terms are each 6e6 and their sum is 13.

### Milestone 2: exact segmentation

- [x] Optimal Partitioning
      why: the base case. Everything after it is an optimisation of this, and
      a wrong answer here is invisible unless something enumerates.
      done: 67 cases, all four costs, four penalties, minseg 1 to 3, n up
      to 25. Worst |op - brute| 2.8e-14 against a tolerance of 1e-9, and
      the change point vectors are identical in every case.
- [x] Brute-force enumeration and the agreement sweep
      why: the hard invariant. This is the strongest oracle in the project
      and it has to exist before PELT, or PELT's agreement with OP will be
      agreement with a shared mistake.
- [x] An independent brute force in Python
      why: R's brute force and R's DP call the same `seg_cost`. If that
      function is wrong they agree and are both wrong. Python's does not.
- [x] PELT, with a pruning counter
      why: the pruning inequality is the classic bug in this algorithm and
      the only way to see it is to require exact agreement with OP while
      proving the pruning actually fired.
- [x] Verify the K = 0 precondition numerically
      why: PELT's exactness rests on `C(s,t) + C(t,T) + K <= C(s,T)`. For a
      maximised log-likelihood cost K = 0 is the usual argument, and the
      variance floor is a modification of the cost that the argument does
      not cover. Check it by enumeration rather than by assertion.
- [x] Minimum segment length, in all three searches
      why: `minseg` changes which candidates are admissible, and an
      off-by-one in the pruning loop's bounds is a defect brute force can
      see and nothing else can.

- [ ] Measure whether the relative variance floor changes any segmentation
      why: `CP_VARFLOOR_REL` was added to stop the floor amplifying
      cancellation noise and then the exact fallback made that unnecessary,
      so it is currently a modelling choice with no measurement behind it.
      The claim to test is that on data of large scale an absolute floor of
      1e-8 makes a near-constant run so cheap that the optimiser brackets it
      with spurious change points. Either that is measurable, or the
      relative floor should be deleted.

### Milestone 3: model selection

- [x] BIC, AIC and Hannan-Quinn penalties
      why: the penalty is the only thing standing between a detector and one
      change point per observation, and the number of parameters a change
      costs differs per cost function.
- [x] MBIC, and the pruning constant it breaks
      why: MBIC adds `log(segment length)` to each segment's cost, which is
      not superadditive, so `K = 0` is no longer valid and PELT silently
      stops being exact. Measure the disagreement rather than assuming it.
- [ ] Measure selected-k against true-k over a seeded corpus
      why: a penalty is a claim about how many change points you will be
      told about. Report the distribution, not an average.

### Milestone 4: is the change real

- [ ] Permutation test for a single change point
      why: a detected change point is a maximum over n-1 candidates, so its
      apparent significance is inflated by the maximisation. The permutation
      distribution of the same maximum is the only honest null here.
- [ ] Calibrate the permutation test under the null
      why: a test whose p-values are not uniform under the null is not a
      test. Measure the uniformity with `ks.test` and the rejection rate at
      alpha = 0.05 over many seeded null series.
- [ ] Measure the power curve
      why: the rejection rate under the null is half the picture, and a test
      with a correct size and no power is a random number generator with
      good manners.

### Milestone 5: does it find the changes that are there

- [ ] Seeded generators for known-answer series, with the seed asserted
      why: a seeded check can pass with the seed ignored. The digest of the
      generated series, not the report about it, is what gets compared.
- [ ] Recovery rate and false-positive rate, over seeds, with variance
      why: a single run is noise and a single accuracy number hides the
      trade. Both rates, both reported, with the spread across seeds.
- [x] Binary segmentation, and the rate at which it is wrong
      why: the greedy method is what most people actually run. Measuring how
      often it differs from the exact answer is a result about the method,
      not a bug in ours.
- [ ] The failure modes: small jumps, short segments, adjacent changes
      why: the README has to say what does not work, and it has to say it
      with numbers from a corpus built to break it.

### Milestone 6: the program

- [ ] Command line entry point and output format
      why: nothing above this line is usable from a shell, and the output
      format is a contract the test suite can read.
- [ ] Input validation, and refusal distinct from a wrong answer
      why: non-numeric input, an empty file, negative counts under the
      Poisson cost and a `minseg` larger than the series all have to be
      refused with exit 1, not answered with exit 0 and rubbish.
- [ ] A worked detection-engineering example
      why: the project is framed as telemetry analysis and so far nothing in
      it has looked at a series that resembles telemetry.

### Milestone 7: proving the suite can fail

- [ ] The mutation harness
      why: every number above is worthless if the assertion that produced it
      cannot go red. Mutants are applied to the source, the suite is run, and
      a mutant that survives is a hole in the suite.
- [ ] Run the mutants and record the survivors
      why: survivors are the output of this exercise, not an embarrassment.
      An equivalent mutant has to be argued for, not asserted.
- [ ] Close the holes the survivors expose
      why: a survivor that is not equivalent names a property nothing checks.

### Milestone 8: scale and honesty

- [ ] Measure PELT's pruning against series length
      why: PELT's claim is linear expected time under conditions that
      include changepoints occurring at a constant rate. Measure the
      candidate counts on series with no change points, where the claim is
      weakest.
- [ ] README, with every number traced to a Makefile target
      why: written last so the limitations are observed rather than guessed
      at, and every figure in it reproducible by a reader.
