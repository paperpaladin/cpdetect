# Assertion counting for the R-side test scripts.
#
# The contract with t/run.sh is two-sided on purpose. A script prints one
# "ok NAME" or "FAIL NAME: detail" line per assertion, and at the end prints
# "done N" where N is how many assertion lines it believes it printed. run.sh
# checks both the exit status and that N equals the number of lines it
# actually saw. Without the second half a script that dies two thirds of the
# way through reports a smaller number of passes and nothing else, which is
# indistinguishable from a suite that got shorter.

.H <- new.env(parent = emptyenv())
.H$pass <- 0L
.H$fail <- 0L

ok <- function(cond, name, detail = "") {
  if (!isTRUE(cond)) {
    .H$fail <- .H$fail + 1L
    cat("FAIL ", name, if (nzchar(detail)) paste0(": ", detail) else "", "\n", sep = "")
  } else {
    .H$pass <- .H$pass + 1L
    cat("ok   ", name, if (nzchar(detail)) paste0(" (", detail, ")") else "", "\n", sep = "")
  }
  invisible(isTRUE(cond))
}

# Floating point comparison with the tolerance named at the call site. There
# is no default: a tolerance you did not choose is a tolerance you cannot
# defend, and the two in use here differ by four orders of magnitude for
# reasons that are written down where they are used.
near <- function(a, b, tol) {
  if (length(a) != length(b)) return(FALSE)
  if (!all(is.finite(a)) || !all(is.finite(b))) return(identical(a, b))
  all(abs(a - b) <= tol)
}

ok_near <- function(a, b, tol, name) {
  d <- if (length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)))
    max(abs(a - b)) else NA_real_
  ok(near(a, b, tol), name,
     sprintf("maxdiff %.3g <= tol %.3g", d, tol))
}

done <- function() {
  cat("done ", .H$pass + .H$fail, "\n", sep = "")
  if (.H$fail > 0L) quit(status = 1, save = "no")
  invisible(NULL)
}

# Sourcing the library under test. Called from every t/*.R so that the set of
# files is in one place.
cp_load <- function() {
  root <- Sys.getenv("CPDETECT_ROOT", unset = ".")
  for (f in sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)))
    source(f)
  invisible(NULL)
}

scratch <- function(name) {
  d <- Sys.getenv("CPDETECT_SCRATCH", unset = "")
  if (!nzchar(d)) stop("CPDETECT_SCRATCH is not set; run through t/run.sh")
  file.path(d, name)
}

# ---------------------------------------------------------------- the oracle
#
# Talking to t/oracle.py. The request goes through a file rather than a pipe
# so that a failing case can be re-run by hand with the exact bytes R sent,
# and every number crosses as %.17g, which round-trips an IEEE double
# exactly. A comparison that loses two digits in the text is a comparison
# against a slightly different number, and would need a tolerance that
# describes the printf format rather than the arithmetic.

g17 <- function(x) sprintf("%.17g", x)

oracle_available <- function() {
  r <- suppressWarnings(system2("python3", c("-c", shQuote("import scipy, numpy")),
                                stdout = FALSE, stderr = FALSE))
  identical(r, 0L)
}

oracle <- function(action, lines, tag = action) {
  req <- scratch(paste0("req_", tag, "_", .H$pass + .H$fail, ".txt"))
  writeLines(lines, req)
  out <- suppressWarnings(
    system2("python3",
            c(shQuote(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "oracle.py")),
              action, shQuote(req)),
            stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status")
  if (!is.null(st) && st != 0)
    stop("oracle.py ", action, " failed (status ", st, "):\n", paste(out, collapse = "\n"))
  out
}

# The common shape: describe a cost object and a list of (s, e) pairs.
oracle_cost_req <- function(cost, S, E) {
  c(paste("type", cost$type),
    paste("sigma2", g17(cost$sigma2)),
    paste("mu", g17(if (is.null(cost$mu)) 0 else cost$mu)),
    paste("varfloor", g17(cost$varfloor)),
    paste("y", paste(g17(cost$y), collapse = " ")),
    sprintf("seg %d %d", S, E))
}
