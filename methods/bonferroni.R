source("methods/sequential_test.R")
source("utils.R")


Bonferroni = R6Class(
  "Bonferroni",
  inherit = SequentialTest,
  public = list(
    significance_level = NULL,
    increment_std = NULL,
    num_checks = NULL,
    initialize = function(name,
                          significance_level,
                          increment_std,
                          num_checks) {
      super$initialize(name)
      self$significance_level = significance_level
      self$increment_std = increment_std
      self$num_checks = num_checks
    },
    boundary = function(num_observations) {
      # A check happens after a whole number of observations, so the schedule is
      # floored before it is used for the index and for the scale alike. Integer
      # division avoids the binary-rounding artefact of trunc(N * j / k), which
      # lands a check one observation early whenever N * j / k is an exact
      # integer that is not exactly representable.
      check_times = (num_observations * seq_len(self$num_checks)) %/% self$num_checks
      return(list(
        index = check_times,
        value = (qnorm(1 - self$significance_level / self$num_checks)
                 * self$increment_std * sqrt(check_times))
      ))
    }
  )
)


# set.seed(2024)
# print(round(measure_detection_rate(
#   Bonferroni$new("Bonferroni14", 0.05, 10, 14),
#   10,
#   500,
#   1000
# ), 2))
