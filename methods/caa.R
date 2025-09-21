source("methods/sequential_test.R")
source("utils.R")


CAA = R6Class(
  "CAA",
  inherit = SequentialTest,
  public = list(
    significance_level = NULL,
    increment_std = NULL,
    initialize = function(name, significance_level, increment_std) {
      super$initialize(name)
      self$significance_level = significance_level
      self$increment_std = increment_std
    },
    monitor = function(trajectory, assignment_indicators=NULL) {
      N = length(trajectory)
      n = 1:N
      standardised_trajectory = trajectory / self$increment_std / sqrt(n)
      return(standardised_trajectory > qnorm(1 - self$significance_level) / (n / N))
    }
  )
)


# set.seed(2024)
# print(round(measure_detection_rate(CAA$new("CAA", 0.05, 10), 10, 500, 1000), 2))
