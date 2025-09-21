source("methods/sequential_test.R")
source("utils.R")


# Netflix's version of Always Valid F-test (mSPRT)
# Effectively implements Robbin's confidence sequence,
# see https://arxiv.org/abs/2210.08589v2, formula (3).

mSPRT = R6Class(
  "mSPRT",
  inherit = SequentialTest,
  public = list(
    significance_level = NULL,
    # significance level
    phi = NULL,
    increment_std = NULL,
    # the standard deviation of the increment
    initialize = function(name,
                          significance_level,
                          increment_std,
                          phi) {
      super$initialize(name)
      self$significance_level = significance_level
      self$increment_std = increment_std
      self$phi = phi
    },
    monitor = function(trajectory, assignment_indicators=NULL) {
      N = length(trajectory)
      two_sided_significance_level = 2 * self$significance_level
      n = 1:N
      boundary = self$increment_std / sqrt(n) * sqrt(log((self$phi + n) / (
        self$phi * two_sided_significance_level ^ 2
      ))
      * (self$phi + n) / n)
      return(trajectory / n > boundary)
    }
  )
)


# set.seed(2024)
# print(round(measure_detection_rate(
#   mSPRT$new("mSPRT25", 0.05, 10, 25), 10, 500, 1000
# ), 2))
