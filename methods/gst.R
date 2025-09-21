if(!require("ldbounds", character.only=TRUE)) {
  install.packages("ldbounds")
  library("ldbounds", character.only=TRUE)
}

source("methods/sequential_test.R")
source("utils.R")


# Over sampling - Calculate bounds using the over-sample trick suggested in literature
# Help-function GST over sample
iter_gstb = function(actual_num_checks,
                     expected_num_checks,
                     significance_level,
                     phi) {
  counter = 0
  bound = numeric()
  if (actual_num_checks > expected_num_checks) {
    for (s in (actual_num_checks - expected_num_checks + 1):actual_num_checks) {
      counter = counter + 1
      bound[counter] = ldBounds(
        t = seq(1 / s, 1, 1 / s),
        iuse = 3,
        phi = phi,
        alpha = significance_level,
        sides = 1
      )$upper.bounds[s]
    }
  }
  return(bound)
}


GST = R6Class(
  "GST",
  inherit = SequentialTest,
  public = list(
    significance_level = NULL,
    phi = NULL,
    expected_num_checks = NULL,
    actual_num_checks = NULL,
    increment_std = NULL,
    thresholds = NULL,
    normalised_actual_check_times = NULL,
    initialize = function(name,
                          significance_level,
                          phi,
                          expected_num_checks,
                          actual_num_checks,
                          increment_std) {
      super$initialize(name)
      self$significance_level = significance_level
      self$phi = phi
      self$expected_num_checks = expected_num_checks
      self$actual_num_checks = actual_num_checks
      self$increment_std = increment_std
      self$thresholds = c(
        ldBounds(
          t = seq(1 / expected_num_checks, 1, 1 / expected_num_checks),
          iuse = 3,
          phi = phi,
          alpha = significance_level,
          sides = 1
        )$upper.bounds,
        iter_gstb(
          actual_num_checks,
          expected_num_checks,
          significance_level,
          phi
        )
      )[1:actual_num_checks]
      self$normalised_actual_check_times = seq(1 / actual_num_checks, 1, 1 / actual_num_checks)
    },
    monitor = function(trajectory, assignment_indicators=NULL) {
      N = length(trajectory)
      standardised_trajectory = trajectory / self$increment_std / sqrt(1:N)
      return(standardised_trajectory[N * self$normalised_actual_check_times]
             > self$thresholds)
    }
  )
)


# set.seed(2024)
# print(round(measure_detection_rate(
#   GST$new("GST14", 0.05, 1, 14, 14, 10), 10, 500, 1000
# ), 2))
#
# set.seed(2024)
# print(round(measure_detection_rate(
#   GST$new("GST14o", 0.05, 1, 14, 28, 10), 10, 500, 1000
# ), 2))
#
# set.seed(2024)
# print(round(measure_detection_rate(
#   GST$new("GST14u", 0.05, 1, 14, 7, 10), 10, 500, 1000
# ), 2))
#
# set.seed(2024)
# print(round(measure_detection_rate(
#   GST$new("GST14phi3o", 0.05, 3, 14, 28, 10), 10, 500, 1000
# ), 2))
