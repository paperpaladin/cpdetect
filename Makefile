# cpdetect -- no compilation step, so `all` is a syntax gate.
#
# `cq buildall` builds the DEFAULT target in a fresh worktree at every commit,
# so the default target has to be something that can fail. R parses a script
# lazily at call time, which means a syntax error three functions down is
# invisible until that function runs; `parse()` on every source is the
# cheapest thing that catches it at commit time rather than at test time.

RSRC := $(wildcard R/*.R) $(wildcard cpdetect) $(wildcard t/*.R)

.PHONY: all test check-syntax recovery mutants bench clean

all: check-syntax

check-syntax:
	@Rscript -e 'bad <- 0; for (f in commandArgs(TRUE)) { r <- tryCatch({parse(f); ""}, error=function(e) conditionMessage(e)); if (nzchar(r)) { cat("parse error:", f, "\n", r, "\n"); bad <- bad + 1 } }; if (bad) quit(status=1); cat("syntax ok:", length(commandArgs(TRUE)), "files\n")' $(RSRC)

test: all
	@sh t/run.sh

# Measured numbers quoted in the README. Separate targets because they are
# slow and because a reader should be able to reproduce one without the rest.
recovery: all
	@Rscript t/measure_recovery.R

mutants: all
	@sh t/mutants.sh

bench: all
	@Rscript t/measure_pruning.R

clean:
	@rm -rf scratch
