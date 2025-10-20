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
