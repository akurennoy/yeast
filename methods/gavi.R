source("methods/sequential_test.R")
source("utils.R")


# Eppo's version of the Generalized Always Valid Inference (GAVI)
# See https://arxiv.org/pdf/1810.08240, Proposition 5 (formula (51)).

GAVI = R6Class(
  "GAVI",
  inherit = SequentialTest,
  public = list(
    significance_level = NULL,
    # significance level
    n_tune = NULL,
    increment_std = NULL,
    # the standard deviation of the increment
    initialize = function(name,
                          significance_level,
                          increment_std,
                          n_tune) {
      super$initialize(name)
      self$significance_level = significance_level
      self$increment_std = increment_std
      self$n_tune = n_tune
    },
    rho = function() {
      two_sided_significance_level = 2 * self$significance_level
      return(self$n_tune / (log(log(
        exp(1) * two_sided_significance_level ^ (-2)
      ))
      - 2 * log(two_sided_significance_level)))
    },
    boundary = function(num_observations) {
      two_sided_significance_level = 2 * self$significance_level
      rho = self$rho()
      n = 1:num_observations
      return(list(
        index = n,
        value = self$increment_std * sqrt((n + rho) * log((n + rho) / (rho * two_sided_significance_level ^
                                                                         2)))
      ))
    }
  )
)


# set.seed(2024)
# print(round(measure_detection_rate(
#   GAVI$new("GAVI250", 0.05, 10, 250), 10, 500, 1000
# ), 2))
