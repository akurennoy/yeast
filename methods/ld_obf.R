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
    boundary = NULL,
    initialize = function(name,
                          significance_level,
                          increment_std) {
      super$initialize(name)
      self$significance_level <- significance_level
      self$increment_std <- increment_std
      self$boundary <- read.csv("obf_bounds_100.csv")$x
    },
    
    monitor = function(trajectory, assignment_indicators=NULL) {
      n = length(trajectory)
      stopifnot(n >= 100)
      return(
        (trajectory / self$increment_std / sqrt(1:n))#[seq(1, n, by=n %/% 100)]
        > rep(self$boundary, each = ceiling(n / length(self$boundary)))[1:n] 
      )
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
