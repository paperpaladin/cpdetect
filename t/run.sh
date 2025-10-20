#!/bin/sh
# cpdetect test suite.
#
# Each t/*.R prints "ok NAME" / "FAIL NAME" lines and a final "done N". This
# script runs them, tallies, and prints "N passed, M failed" as its last line,
# which is the form the collection's tooling parses.
#
# Three ways a script can be wrong that are not an assertion going red, and
# all three are checked rather than assumed:
#
#   exit status nonzero      - R threw, or done() saw a failure
#   missing "done" line      - the script died before reaching the end
#   done N != lines seen     - the script skipped assertions it meant to run
#
# The last two overlap, and the overlap was measured rather than assumed:
# with the done-count cross-check disabled the "liar" self-check goes red,
# but with the done-line-present check disabled nothing goes red, because a
# script that dies mid-way is also caught by the count comparison (done 0
# against one assertion seen). The present-check earns its place on exactly
# one case, a script that calls done() twice, so that case is planted too.
#
# The last two are the ones that matter. A suite that quietly gets shorter
# reports a smaller number and nothing else. SELFCHECK below plants one of
# each and requires this script to catch it, because a check nobody has
# watched fail is not known to be a check.

set -u

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(dirname "$HERE")
cd "$ROOT" || exit 1

SCRATCH=$(mktemp -d "${TMPDIR:-/tmp}/cpdetect.XXXXXX") || exit 1
trap 'rm -rf "$SCRATCH"' EXIT INT TERM HUP
CPDETECT_SCRATCH="$SCRATCH"
CPDETECT_ROOT="$ROOT"
export CPDETECT_SCRATCH CPDETECT_ROOT

command -v Rscript >/dev/null 2>&1 || {
  echo "cpdetect: Rscript is not installed, and it is the only runtime there is" >&2
  echo "0 passed, 0 failed" >&2
  exit 2
}

pass=0
fail=0

# Sets RF_PASS, RF_FAIL, RF_PROBLEM. RF_PROBLEM is empty when the script
# behaved; the assertion tallies are still returned so a partial run is
# visible rather than discarded.
run_file() {
  rf_src=$1
  rf_out="$SCRATCH/$(basename "$rf_src").out"
  RF_PROBLEM=""
  if Rscript "$rf_src" >"$rf_out" 2>&1; then rf_st=0; else rf_st=$?; fi

  set -- $(awk '
    /^ok   /   { p++ }
    /^FAIL /   { f++ }
    /^done /   { d = $2; nd++ }
    END        { printf "%d %d %d %d", p+0, f+0, d+0, nd+0 }' "$rf_out")
  RF_PASS=$1
  RF_FAIL=$2
  rf_claim=$3
  rf_ndone=$4

  if [ "$rf_ndone" -ne 1 ]; then
    if [ "$rf_ndone" -eq 0 ]; then
      RF_PROBLEM="no done line: the script did not reach its end"
    else
      RF_PROBLEM="$rf_ndone done lines: the script ended more than once"
    fi
  elif [ "$rf_claim" -ne $((RF_PASS + RF_FAIL)) ]; then
    RF_PROBLEM="done says $rf_claim, saw $((RF_PASS + RF_FAIL)) assertion lines"
  elif [ "$rf_st" -ne 0 ] && [ "$RF_FAIL" -eq 0 ]; then
    RF_PROBLEM="exit status $rf_st with no FAIL line"
  elif [ $((RF_PASS + RF_FAIL)) -eq 0 ]; then
    RF_PROBLEM="ran and asserted nothing"
  fi
}

# ---------------------------------------------------------------- self-check
#
# Four planted scripts, each broken in one of the ways run_file claims to
# notice, plus one that is fine. If run_file waves any of the four through,
# every number below it is decoration.
echo "--- harness self-check"
SELF="$SCRATCH/self"
mkdir -p "$SELF"

cat > "$SELF/clean.R" <<'RB'
source(file.path(Sys.getenv("CPDETECT_ROOT"), "t", "harness.R"))
ok(TRUE, "planted_true")
done()
RB
cat > "$SELF/redassert.R" <<'RB'
source(file.path(Sys.getenv("CPDETECT_ROOT"), "t", "harness.R"))
ok(FALSE, "planted_false")
done()
RB
cat > "$SELF/midexit.R" <<'RB'
source(file.path(Sys.getenv("CPDETECT_ROOT"), "t", "harness.R"))
ok(TRUE, "planted_before")
quit(status = 0, save = "no")
ok(TRUE, "never_reached")
done()
RB
cat > "$SELF/throws.R" <<'RB'
source(file.path(Sys.getenv("CPDETECT_ROOT"), "t", "harness.R"))
ok(TRUE, "planted_before_throw")
stop("planted error")
done()
RB
cat > "$SELF/twodone.R" <<'RB'
source(file.path(Sys.getenv("CPDETECT_ROOT"), "t", "harness.R"))
ok(TRUE, "planted_one")
done()
ok(TRUE, "planted_two")
done()
RB
cat > "$SELF/liar.R" <<'RB'
source(file.path(Sys.getenv("CPDETECT_ROOT"), "t", "harness.R"))
ok(TRUE, "planted_one")
cat("done 7\n")
RB

self_expect() {  # name, want-problem(yes/no)
  run_file "$SELF/$1.R"
  if [ "$2" = yes ]; then
    if [ -n "$RF_PROBLEM" ] || [ "$RF_FAIL" -gt 0 ]; then
      echo "ok   selfcheck_$1 (caught: ${RF_PROBLEM:-$RF_FAIL FAIL lines})"
      pass=$((pass + 1))
    else
      echo "FAIL selfcheck_$1: harness accepted a script it must reject"
      fail=$((fail + 1))
    fi
  else
    if [ -z "$RF_PROBLEM" ] && [ "$RF_FAIL" -eq 0 ] && [ "$RF_PASS" -eq 1 ]; then
      echo "ok   selfcheck_$1"
      pass=$((pass + 1))
    else
      echo "FAIL selfcheck_$1: harness rejected a clean script: ${RF_PROBLEM:-ok=$RF_PASS fail=$RF_FAIL}"
      fail=$((fail + 1))
    fi
  fi
}

self_expect clean     no
self_expect redassert yes
self_expect midexit   yes
self_expect throws    yes
self_expect twodone   yes
self_expect liar      yes

# ------------------------------------------------------------------- the suite
problems=""
nfiles=0
for f in "$HERE"/test_*.R; do
  [ -e "$f" ] || continue
  nfiles=$((nfiles + 1))
  base=$(basename "$f" .R)
  echo "--- $base"
  run_file "$f"
  sed -n 's/^\(ok   \|FAIL \)/&/p' "$SCRATCH/$(basename "$f").out"
  pass=$((pass + RF_PASS))
  fail=$((fail + RF_FAIL))
  if [ -n "$RF_PROBLEM" ]; then
    echo "FAIL $base harness: $RF_PROBLEM"
    sed -n '1,40p' "$SCRATCH/$(basename "$f").out" | sed 's/^/    | /'
    fail=$((fail + 1))
    problems="$problems $base"
  fi
done

if [ "$nfiles" -eq 0 ]; then
  echo "FAIL no test_*.R found: the suite would otherwise report a clean zero"
  fail=$((fail + 1))
fi

echo
[ -n "$problems" ] && echo "scripts with harness problems:$problems"
echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ] || exit 1
exit 0
