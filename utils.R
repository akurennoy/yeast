library(data.table)
library(R6)

source("methods/sequential_test.R")


measure_fdr = function(sequential_test,
                       increment_std,
                       num_observations,
                       num_replications) {
  num_detections = 0
  for (r in 1:num_replications) {
    increments = rnorm(num_observations, 0, increment_std)
    if (any(sequential_test$monitor(cumsum(increments)))) {
      num_detections = num_detections + 1
    }
  }
  return(as.numeric(num_detections) / num_replications)
}


get_savings = function(detection_indicators,
                       check_times,
                       num_observations) {
  index_of_first_detection = which(detection_indicators == TRUE)[1]
  return(ifelse(
    is.na(index_of_first_detection),
    0.0,
    1.0 - ceiling(check_times[index_of_first_detection]) / num_observations
  ))
}


Aggregator = R6Class(
  "Aggregator",
  public = list(
    EFFECT_COL = "effect",
    MODE_COL = "mode",
    METHOD_COL = "method",
    NUM_TRIALS_COL = "num_trials",
    NUM_DETECTIONS_COL = "num_detections",
    TOTAL_SAVINGS_COL = "total_savings",
    DETECTION_RATE_COL = "detection_rate",
    AVERAGE_SAVINGS_COL = "average_savings",
    df = NULL,
    initialize = function() {
      
    },
    update = function(relative_effect,
                      monitoring_mode,
                      method_name,
                      detection_indicator,
                      savings) {
      row_index <- which(self$df[self$EFFECT_COL] == relative_effect
                         & self$df[self$MODE_COL] == monitoring_mode
                         & self$df[self$METHOD_COL] == method_name)
      if (length(row_index) == 0) {
        new_row = data.frame(relative_effect,
                             monitoring_mode,
                             method_name,
                             1,
                             detection_indicator,
                             savings)
        names(new_row) = list(
          self$EFFECT_COL,
          self$MODE_COL,
          self$METHOD_COL,
          self$NUM_TRIALS_COL,
          self$NUM_DETECTIONS_COL,
          self$TOTAL_SAVINGS_COL
        )
        self$df = rbind(self$df, new_row)
      } else {
        stopifnot(length(row_index) == 1)
        self$df[row_index, self$NUM_TRIALS_COL] = self$df[row_index, self$NUM_TRIALS_COL] + 1
        self$df[row_index, self$NUM_DETECTIONS_COL] = self$df[row_index, self$NUM_DETECTIONS_COL] + detection_indicator
        self$df[row_index, self$TOTAL_SAVINGS_COL] = self$df[row_index, self$TOTAL_SAVINGS_COL] + savings
      }
    },
    get_result = function() {
      result = data.table(self$df)
      result[, self$DETECTION_RATE_COL := as.numeric(.SD[[self$NUM_DETECTIONS_COL]]) / .SD[[self$NUM_TRIALS_COL]]]
      result[, self$AVERAGE_SAVINGS_COL := .SD[[self$TOTAL_SAVINGS_COL]] / .SD[[self$NUM_TRIALS_COL]]]
      output_columns = c(
        self$EFFECT_COL,
        self$MODE_COL,
        self$METHOD_COL,
        self$DETECTION_RATE_COL,
        self$AVERAGE_SAVINGS_COL
      )
      return(result)#[, ..output_columns])
    }
  )
)


ResultReporter = R6Class(
  "ResultReporter",
  public = list(
    DETECTION_RATE_COL = "detection_rate",
    AVERAGE_SAVINGS_COL = "average_savings",
    result = NULL,
    continuous_methods = NULL,
    initialize = function(result, continuous_methods) {
      self$result = result
      self$continuous_methods = continuous_methods
    },
    print_main_table = function(metric) {
      continuous_method_names = unlist(lapply(self$continuous_methods, function(st)
        st$name))
      print(dcast(self$result[mode == "stream"], method ~ effect, value.var = metric)[J(continuous_method_names)])
    },
    print_detection_rate = function() {
      print("False Detection Rate and Power")
      self$print_main_table(self$DETECTION_RATE_COL)
    },
    print_savings = function() {
      print("Average Sample/Time Savings")
      self$print_main_table(self$AVERAGE_SAVINGS_COL)
    },
    print_additional_results = function(relative_effect=0.2) {
      row_keys = c(unlist(lapply(self$continuous_methods, function(st)
        st$name)),
        "GST",
        "GSTu",
        "GSTo",
        "GSTphi3",
        "GSTuphi3",
        "GSTophi3")
      print("Additional Results: False Detection Rate (relative effect = 0.0)")
      print(dcast(self$result[effect == 0.0], method ~ mode, value.var = self$DETECTION_RATE_COL)[J(row_keys)])
      print(sprintf("Additional Results: Power (relative effect = %.1f)", relative_effect))
      print(dcast(self$result[effect == relative_effect], method ~ mode, value.var = self$DETECTION_RATE_COL)[J(row_keys)])
    }
  )
)


report = function(result, continuous_methods) {
  reporter = ResultReporter$new(result, continuous_methods)
  reporter$print_detection_rate()
  reporter$print_savings()
  reporter$print_additional_results()
}
