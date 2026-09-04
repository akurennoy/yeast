# Runs every formula-level test.
#
# Run from the repository root:  Rscript tests/run_all.R
#
# These are fast, deterministic checks of the boundary and betting formulas.
# They are intended to be run before any full experiment, and they must pass
# from a clean checkout.

TEST_FILES = c(
  "tests/test_aggregator.R",
  "tests/test_check_schedule.R",
  "tests/test_boundaries.R",
  "tests/test_msprt_boundary.R",
  "tests/test_ons_update.R"
)

for (test_file in TEST_FILES) {
  cat(sprintf("\n--- %s ---\n", test_file))
  source(test_file, local = new.env())
}

cat("\nAll tests passed.\n")
