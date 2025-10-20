# The harness checking itself. t/run.sh plants broken scripts and requires
# them to be rejected; this file covers the pieces run.sh cannot see, which
# is the tolerance comparison every later test is built on.

source(file.path(Sys.getenv("CPDETECT_ROOT", unset = "."), "t", "harness.R"))

ok(near(1, 1 + 1e-12, 1e-9), "near_accepts_within_tolerance")
ok(!near(1, 1 + 1e-6, 1e-9), "near_rejects_outside_tolerance")
# Written first as near(1, 1 + 1e-9, 1e-9) and it failed, which is the point
# of writing the boundary case at all: 1 + 1e-9 is not representable, the
# nearest double is 1.000000082740371e-09 away, and the comparison is
# correctly outside the tolerance. The boundary has to be built out of
# numbers whose difference is exact.
ok(near(0, 1e-9, 1e-9), "near_boundary_is_inclusive")
ok(!near(1, 1 + 1e-9, 1e-9), "near_boundary_1e9_above_one_is_not_exact")
ok(!near(c(1, 2), c(1, 2, 3), 1e-9), "near_rejects_length_mismatch")

# Non-finite values are where a tolerance comparison stops meaning anything:
# Inf - Inf is NaN and NaN <= tol is NA, which isTRUE() turns into FALSE
# rather than into an error. near() takes the identical() branch instead, so
# Inf == Inf passes and Inf vs a finite number does not.
ok(near(Inf, Inf, 1e-9), "near_inf_equals_inf")
ok(!near(Inf, 1e300, 1e-9), "near_inf_is_not_a_large_number")
ok(!near(NaN, NaN, 1e-9) == FALSE, "near_nan_is_not_equal_to_itself_by_identical")
ok(!near(0, NA_real_, 1e-9), "near_na_is_not_near_anything")

# A tolerance of zero is exact equality, and the suite uses it in one place
# (integer-valued counts), so it has to behave.
ok(near(2, 2, 0), "near_zero_tolerance_exact")
# Also written wrong first: .Machine$double.eps is the spacing at 1, and at
# 2 the spacing is twice that, so 2 + eps rounds back to 2 and the
# comparison was true. One ulp at 2 is 2^-51.
ok(near(2, 2 + .Machine$double.eps, 0), "near_zero_tolerance_half_ulp_rounds_away")
ok(!near(2, 2 + 2^-51, 0), "near_zero_tolerance_rejects_one_ulp")

ok(nzchar(Sys.getenv("CPDETECT_SCRATCH")), "scratch_dir_is_exported")
ok(dir.exists(Sys.getenv("CPDETECT_SCRATCH")), "scratch_dir_exists")
ok(!file.exists(scratch("nothing_here")), "scratch_path_is_empty_to_start")

done()
