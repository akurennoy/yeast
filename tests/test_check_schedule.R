# Pins the discrete check schedule shared by the scheduled-look methods and by
# the discrete-monitoring mode of run_simulations.R.
#
# The property that matters is that k checks over N observations land on whole
# observation indices, are strictly increasing, and end exactly at N. The
# floating-point form N * seq(1/k, 1, 1/k) fails the last of these at
# k = 14, 28, 56 for N = 500, dropping the terminal check.

source("utils.R")

check_that = function(condition, description) {
  if (!isTRUE(condition)) {
    stop(sprintf("FAILED: %s", description), call. = FALSE)
  }
  cat(sprintf("  ok: %s\n", description))
}

HORIZONS = c(500, 503, 1000)
CHECK_COUNTS = c(1, 2, 7, 14, 28, 42, 56, 100)

for (num_observations in HORIZONS) {
  for (num_checks in CHECK_COUNTS) {
    times = discrete_check_times(num_observations, num_checks)
    label = sprintf("N = %d, k = %d", num_observations, num_checks)
    check_that(length(times) == num_checks,
               sprintf("%s: k check times", label))
    check_that(all(times == as.integer(times)),
               sprintf("%s: whole observation indices", label))
    check_that(all(diff(times) > 0) || num_checks == 1,
               sprintf("%s: strictly increasing", label))
    check_that(times[num_checks] == num_observations,
               sprintf("%s: last check at N", label))
    check_that(times[1] >= 1, sprintf("%s: first check at or after 1", label))
  }
}

# A check count equal to the horizon must reduce to a check after every
# observation. The float form returns 0 as its first element here, and index 0
# silently drops an element in R rather than erroring.
check_that(all(discrete_check_times(503, 503) == 1:503),
           "k = N reduces to every observation")

cat("test_check_schedule.R: all checks passed\n")
