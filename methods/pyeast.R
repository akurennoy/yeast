library(mvtnorm)
library(stats)

source("methods/sequential_test.R")
source("utils.R")


ell = function(h, k, rho) {
  return (1 - (pnorm(h) + pnorm(k))
          + pmvnorm(
            upper = c(h, k),
            mean = c(0, 0),
            sigma = matrix(c(1, rho, rho, 1), nrow = 2)
          ))
}


compute_convolution = function(s, sigma, b_1, b_2) {
  # Computes the convolution of a centered normal and a truncated normal at b_2.
  #
  # s - the std of the centered normal
  # sigma - the std of the untruncated "version" of the truncated normal
  # b_1 - the (upper) truncation (cutoff) point
  #
  zeta = sqrt(sigma ^ 2 + s ^ 2)
  rho = sigma / sqrt(sigma ^ 2 + s ^ 2)
  h = b_2 / zeta
  k = b_1 / sigma
  return (1 - (1 - pnorm(h) - ell(h, k, rho)) / pnorm(k))
}


estimate_false_detection_rate_bound = function(cumulative_num_observations_by_period,
                                               thresholds,
                                               increment_std) {
  stopifnot(length(cumulative_num_observations_by_period) > 0)
  stopifnot(length(cumulative_num_observations_by_period) == length(thresholds))
  
  N_1 = cumulative_num_observations_by_period[1]
  false_detection_rate = (2 * (1 - pnorm(thresholds[1] / (
    sqrt(N_1) * increment_std
  ))))
  for (k in 2:length(cumulative_num_observations_by_period)) {
    Z_k = pnorm(thresholds[k - 1]
                / (
                  sqrt(cumulative_num_observations_by_period[k - 1]) * increment_std
                ))
    N_k = (
      cumulative_num_observations_by_period[k]
      - cumulative_num_observations_by_period[k - 1]
    )
    J_k = compute_convolution(
      sqrt(N_k) * increment_std,
      sqrt(cumulative_num_observations_by_period[k - 1]) * increment_std,
      thresholds[k - 1],
      thresholds[k]
    )
    false_detection_rate = false_detection_rate + 2 * Z_k * (1 - J_k)
  }
  return(false_detection_rate)
}


compute_pYEAST_thresholds = function(expected_cumulative_num_observations_by_period,
                                   increment_std,
                                   significance_level) {
  thresholds = (
    qnorm(1 - significance_level / 2)
    * sqrt(expected_cumulative_num_observations_by_period)
    * increment_std
  )
  while (estimate_false_detection_rate_bound(expected_cumulative_num_observations_by_period,
                                             thresholds,
                                             increment_std) > significance_level) {
    thresholds = thresholds * 1.005
  }
  # cat("FPR bound =", estimate_false_detection_rate_bound(
  #   expected_cumulative_num_observations_by_period, thresholds, increment_std)
  # )
  # cat("\nThresholds:", thresholds)
  return(thresholds)
}


pYEAST = R6Class(
  "pYEAST",
  inherit = SequentialTest,
  public = list(
    significance_level = NULL,
    actual_cumulative_num_observations_by_period = NULL,
    increment_std = NULL,
    thresholds = NULL,
    initialize = function(name,
                          significance_level,
                          expected_cumulative_num_observations_by_period,
                          # an array with K elements where K is the number of periods;
                          # contains the cumulative sum of the expected number of observations
                          # for each period
                          increment_std,
                          actual_cumulative_num_observations_by_period =
                            NULL) {
      super$initialize(name)
      self$significance_level = significance_level
      self$increment_std = increment_std
      self$thresholds = compute_pYEAST_thresholds(
        expected_cumulative_num_observations_by_period,
        increment_std,
        significance_level
      )  # obtaining the thresholds - one for each period
      if (is.null(actual_cumulative_num_observations_by_period)) {
        self$actual_cumulative_num_observations_by_period = expected_cumulative_num_observations_by_period
      } else {
        self$actual_cumulative_num_observations_by_period = actual_cumulative_num_observations_by_period
      }
    },
    monitor_ = function(trajectory,
                        actual_cumulative_num_observations_by_period) {
      N = length(trajectory)
      K = length(actual_cumulative_num_observations_by_period)
      stopifnot(N == actual_cumulative_num_observations_by_period[K])
      stopifnot(K == length(self$thresholds))
      
      # producing an alerting boundary by repeating the corresponding threshold
      # value for each observation within a period
      boundary = rep(NA, N)
      boundary[1:actual_cumulative_num_observations_by_period[1]] = self$thresholds[1]
      for (k in 2:K) {
        boundary[(actual_cumulative_num_observations_by_period[k - 1] + 1):actual_cumulative_num_observations_by_period[k]] = self$thresholds[k]
      }
      
      return(trajectory > boundary)
    },
    monitor = function(trajectory, assignment_indicators=NULL) {
      return(self$monitor_(
        trajectory,
        self$actual_cumulative_num_observations_by_period
      ))
    }
  )
)


# set.seed(2024)
# print(round(measure_detection_rate(
#   pYEAST$new("pYEAST7", 0.05, round((1:7) * (500 / 7)), 10), 10, 500, 1000
# ), 2))
# 
# set.seed(2024)
# print(round(measure_detection_rate(
#   pYEAST$new("pYEAST14", 0.05, round((1:14) * (500 / 14)), 10), 10, 500, 1000
# ), 2))
