source("methods/sequential_test.R")
source("utils.R")


YEAST = R6Class(
  "YEAST",
  inherit = SequentialTest,
  public = list(
    significance_level = NULL,
    # significance level
    expected_num_observations = NULL,
    increment_std = NULL,
    # the standard deviation of the increment
    initialize = function(name,
                          significance_level,
                          expected_num_observations,
                          increment_std) {
      super$initialize(name)
      self$significance_level = significance_level
      self$expected_num_observations = expected_num_observations
      self$increment_std = increment_std
    },
    monitor = function(trajectory) {
      boundary = (
        qnorm(1 - self$significance_level / 2)
        * sqrt(self$expected_num_observations) * self$increment_std
      )
      return(trajectory > boundary)
    }
  )
)


# set.seed(2024)
# print(round(measure_fdr(YEAST$new("YEAST", 0.05, 500, 10), 10, 500, 1000), 2))
