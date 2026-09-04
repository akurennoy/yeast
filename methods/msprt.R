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
    # Radius of the Robbins normal-mixture confidence sequence for the mean, so
    # that the test rejects when the running mean exceeds it. The factor
    # (phi + n) / n sits inside the square root.
    mean_scale_boundary = function(num_observations) {
      two_sided_significance_level = 2 * self$significance_level
      n = 1:num_observations
      return(self$increment_std / sqrt(n) * sqrt(log((self$phi + n) / (
        self$phi * two_sided_significance_level ^ 2
      ))
      * (self$phi + n) / n))
    },
    boundary = function(num_observations) {
      n = 1:num_observations
      return(list(index = n, value = n * self$mean_scale_boundary(num_observations)))
    }
  )
)


# set.seed(2024)
# print(round(measure_detection_rate(
#   mSPRT$new("mSPRT25", 0.05, 10, 25), 10, 500, 1000
# ), 2))
