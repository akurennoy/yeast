if(!require("ldbounds", character.only=TRUE)) {
  install.packages("ldbounds")
  library("ldbounds", character.only=TRUE)
}

source("methods/sequential_test.R")
source("utils.R")

library(stats)


LanDeMetsOBF <- R6Class(
  "LanDeMetsOBF",
  inherit = SequentialTest,
  public = list(
    significance_level = NULL,
    increment_std = NULL,
    critical_values = NULL,
    initialize = function(name,
                          significance_level,
                          increment_std) {
      super$initialize(name)
      self$significance_level <- significance_level
      self$increment_std <- increment_std
      # The 100 scheduled-look OBF critical values are cached on disk. Generate
      # them on first use so that the pipeline runs from a clean checkout.
      bounds_file <- "obf_bounds_100.csv"
      if (!file.exists(bounds_file)) {
        cb <- commonbounds(100, iuse = "OF", alpha = significance_level, sides = 1)
        write.csv(cb$upper.bounds, file = bounds_file, row.names = FALSE)
      }
      self$critical_values <- read.csv(bounds_file)$x
      stopifnot(length(self$critical_values) == 100)
    },
    
    # The 100 scheduled-look critical values are held piecewise constant between
    # looks, so the test is applied after every event rather than only at the
    # 100 scheduled looks.
    boundary = function(num_observations) {
      stopifnot(num_observations >= 100)
      n = 1:num_observations
      critical_values = rep(
        self$critical_values,
        each = ceiling(num_observations / length(self$critical_values))
      )[n]
      return(list(
        index = n,
        value = critical_values * self$increment_std * sqrt(n)
      ))
    }
  )
)

# cb <- commonbounds(100, iuse = "OF", alpha = 0.05, sides = 1)
# upper_bounds <- cb$upper.bounds
# write.csv(upper_bounds, file = "obf_bounds_100.csv", row.names = FALSE)


# set.seed(2024)
# print(round(measure_detection_rate(
#   LanDeMetsOBF$new("Lan-DeMets OBF", 0.05, 10), 10, 500, 1000
# ), 3))
