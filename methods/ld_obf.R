source("methods/sequential_test.R")
source("utils.R")

library(stats)


LanDeMetsOBF <- R6Class(
  "LanDeMetsOBF",
  inherit = SequentialTest,
  public = list(
    significance_level = NULL,
    increment_std = NULL,
    initialize = function(name,
                          significance_level,
                          increment_std) {
      super$initialize(name)
      self$significance_level <- significance_level
      self$increment_std <- increment_std
    },
    
    # O'Brien-Fleming alpha-spending function
    obf_spending = function(information_fraction) {
      if (information_fraction <= 0)
        return(0)
      if (information_fraction >= 1)
        return(self$significance_level)
      z_alpha = qnorm(1 - self$significance_level / 2)
      spent_alpha = 2 * (1 - pnorm(z_alpha / sqrt(information_fraction)))
      return(min(spent_alpha, self$significance_level))
    },
    
    monitor = function(trajectory, assignment_indicators=NULL) {
      n <- length(trajectory)
      stopping <- logical(n)
      for (i in seq_len(n)) {
        info_frac <- i / n
        spent_alpha <- self$obf_spending(info_frac)
        z <- trajectory[i] / sqrt(i) / self$increment_std
        critical_value <- qnorm(1 - spent_alpha / 2)
        stopping[i] <- z > critical_value
      }
      return(stopping)
    }
  )
)


# set.seed(2024)
# print(round(measure_detection_rate(
#   LanDeMetsOBF$new("Lan-DeMets OBF", 0.05, 10), 10, 500, 1000
# ), 2))
