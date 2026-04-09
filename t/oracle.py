#!/usr/bin/env python3
"""Independent implementation of cpdetect's arithmetic, in Python.

This exists to be a different language, so it deliberately shares nothing
with R/cost.R but the definitions. It computes every segment cost in two
passes over the segment, from scipy.stats.norm.logpdf and
scipy.stats.poisson.logpmf, rather than from prefix sums; it enumerates
segmentations with itertools.combinations rather than by dynamic
programming; and it accumulates with math.fsum where the order of summation
could matter. If the R side and this side agree, the agreement is about the
definitions rather than about one person's arithmetic.

scipy and numpy are used HERE and nowhere in the program under test.
cpdetect itself is base R plus stats.

Protocol: a request file on argv[1], a response on stdout. The request is
one "key value..." per line; the response is one value per line, printed at
full double precision with %.17g so that nothing is lost in the text.

Requests:

  cost      type/params/y[/mbic], then one "seg s e" per line
            -> one cost per seg line
  brute     type/params/y/penalty/minseg/maxcp
            -> "best <penalised cost>", "cp <t1> <t2> ...", "n <count>"
  amoc      type/params/y/minseg
            -> "stat <max cost reduction>", "at <argmax>"
"""

import sys
import math
from itertools import combinations

import numpy as np
from scipy import stats


def fmt(x):
    return "%.17g" % x


def read_request(path):
    req = {"seg": []}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            key, _, rest = line.partition(" ")
            if key == "seg":
                a, b = rest.split()
                req["seg"].append((int(a), int(b)))
            elif key == "y":
                req["y"] = np.array([float(v) for v in rest.split()])
            elif key in ("type",):
                req[key] = rest.strip()
            else:
                req[key] = float(rest)
    return req


# ------------------------------------------------------------------- costs
#
# Each takes the raw observations of one segment and returns -2 log L at the
# maximum likelihood parameters, with the variance constrained from below by
# `varfloor`. Written from the density functions rather than from a
# simplified formula on purpose: a simplification is where a factor of two
# goes missing, and simplifying the same way in both languages would hide it.

def cost_mean(v, p):
    mu = math.fsum(v) / len(v)
    return -2.0 * math.fsum(stats.norm.logpdf(v, loc=mu, scale=math.sqrt(p["sigma2"])))


def cost_var(v, p):
    mu = p["mu"]
    s2 = max(math.fsum((v - mu) ** 2) / len(v), p["varfloor"])
    return -2.0 * math.fsum(stats.norm.logpdf(v, loc=mu, scale=math.sqrt(s2)))


def cost_meanvar(v, p):
    mu = math.fsum(v) / len(v)
    s2 = max(math.fsum((v - mu) ** 2) / len(v), p["varfloor"])
    return -2.0 * math.fsum(stats.norm.logpdf(v, loc=mu, scale=math.sqrt(s2)))


def cost_poisson(v, p):
    lam = math.fsum(v) / len(v)
    if lam == 0.0:
        # poisson.logpmf(0, 0) is 0 in scipy, but logpmf(k>0, 0) is -inf and
        # an all-zero segment is the only way lam can be 0, so the branch is
        # about being explicit rather than about avoiding a nan.
        return 0.0
    return -2.0 * math.fsum(stats.poisson.logpmf(v, mu=lam))


COSTS = {"mean": cost_mean, "var": cost_var,
         "meanvar": cost_meanvar, "poisson": cost_poisson}


def seg_cost(req, s, e):
    c = COSTS[req["type"]](req["y"][s:e], req)
    if req.get("mbic", 0.0):
        # Zhang and Siegmund's log(segment length) term, which is what makes
        # their criterion segment-additive. Written here independently from
        # the R side rather than taken from it: the whole point of a second
        # implementation is that it is not a transcription.
        c += math.log(e - s)
    return c


# ------------------------------------------------------------------ actions

def do_cost(req, out):
    for s, e in req["seg"]:
        out.append(fmt(seg_cost(req, s, e)))


def do_brute(req, out):
    """Every admissible segmentation, scored. The definition of the answer.

    combinations() over the interior boundaries rather than a recursion, so
    that the count of segmentations considered is a thing this can report and
    the R side can be required to match. A disagreement about how many
    segmentations exist is a disagreement about minseg, and that is worth
    separating from a disagreement about which is cheapest.
    """
    n = len(req["y"])
    minseg = int(req.get("minseg", 1))
    maxcp = int(req.get("maxcp", n))
    beta = req["penalty"]
    best, best_cp, count = math.inf, None, 0
    interior = range(minseg, n - minseg + 1)
    kmax = min(maxcp, max(0, n // minseg - 1))
    for k in range(0, kmax + 1):
        for cps in combinations(interior, k):
            bounds = (0,) + cps + (n,)
            if any(bounds[i + 1] - bounds[i] < minseg for i in range(len(bounds) - 1)):
                continue
            count += 1
            c = math.fsum(seg_cost(req, bounds[i], bounds[i + 1])
                          for i in range(len(bounds) - 1)) + k * beta
            if c < best:
                best, best_cp = c, cps
    out.append("best " + fmt(best))
    out.append("cp " + " ".join(str(t) for t in best_cp))
    out.append("n %d" % count)


def do_amoc(req, out):
    """The at-most-one-change statistic: the largest cost reduction over all
    admissible single splits. The permutation test's observed value."""
    n = len(req["y"])
    minseg = int(req.get("minseg", 1))
    whole = seg_cost(req, 0, n)
    best, at = -math.inf, 0
    for t in range(minseg, n - minseg + 1):
        d = whole - (seg_cost(req, 0, t) + seg_cost(req, t, n))
        if d > best:
            best, at = d, t
    out.append("stat " + fmt(best))
    out.append("at %d" % at)


ACTIONS = {"cost": do_cost, "brute": do_brute, "amoc": do_amoc}


def main(argv):
    if len(argv) != 3 or argv[1] not in ACTIONS:
        sys.stderr.write("usage: oracle.py {%s} REQUEST\n" % "|".join(ACTIONS))
        return 2
    req = read_request(argv[2])
    req.setdefault("sigma2", 1.0)
    req.setdefault("varfloor", 1e-8)
    req.setdefault("mu", 0.0)
    out = []
    ACTIONS[argv[1]](req, out)
    sys.stdout.write("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
