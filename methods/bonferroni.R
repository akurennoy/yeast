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
    monitor = function(trajectory, assignment_indicators=NULL) {
      N = length(trajectory)
      n = 1:N
      standardised_trajectory = trajectory / self$increment_std / sqrt(n)
      return(standardised_trajectory[N * seq(1 / self$num_checks, 1, 1 / self$num_checks)] > qnorm(1 - self$significance_level / self$num_checks))
    }
  )
)


# set.seed(2024)
# print(round(measure_fdr(
#   Bonferroni$new("Bonferroni14", 0.05, 10, 14),
#   10,
#   500,
#   1000
# ), 2))
