if(!require("data.table", character.only = TRUE)) {
  install.packages("data.table")
  library("data.table", character.only = TRUE)
}
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
      return(result)
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
    print_additional_results = function(relative_effect = 0.2) {
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
      print(sprintf(
        "Additional Results: Power (relative effect = %.1f)",
        relative_effect
      ))
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


DataCleaner = R6Class(
  "DataCleaner",
  public = list(
    metric_col = NULL,
    user_id_col = NULL,
    dttm_col = NULL,
    event_value_cutoff = NULL,
    initialize = function(historical_data,
                          metric_col,
                          user_id_col,
                          dttm_col,
                          q = 0.999) {
      dt = setorderv(data.table(historical_data), dttm_col)
      
      cols = c(user_id_col, metric_col)
      event_value_dt = dt[, ..cols][, .(metric_sum = sum(get(metric_col))), by = user_id_col]
      event_value_cutoff = quantile(event_value_dt[, .SD[["metric_sum"]]], q)
      
      self$metric_col = metric_col
      self$user_id_col = user_id_col
      self$dttm_col = dttm_col
      self$event_value_cutoff = event_value_cutoff
    },
    clean = function(data) {
      dt = setorderv(data.table(data), self$dttm_col)
      dt[, metric_cumsum := cumsum(get(self$metric_col)), by = eval(self$user_id_col)]
      return(dt[metric_cumsum <= self$event_value_cutoff])
    }
  )
)


test = function() {
  dt = data.table(
    user_id = c("u1", "u2", "u1", "u1", "u2", "u2"),
    metric = c(99.99, 99.99, 164.00, 100.00, 256.50, 256.50),
    dttm = c(
      "2023-08-01T10:00:00.146Z",
      "2023-08-01T10:00:00.146Z",
      "2023-08-01T11:00:00.146Z",
      "2023-08-01T10:30:00.146Z",
      "2023-08-02T09:00:00.146Z",
      "2023-08-02T09:01:00.146Z"
    )
  )
  
  dc = DataCleaner$new(dt, "metric", "user_id", "dttm")
  cleaned_dt = dc$clean(dt)
  stopifnot(nrow(cleaned_dt) == 5)
  stopifnot(max(cleaned_dt[user_id == "u1", "metric"]) == 164.00)
}


test()
