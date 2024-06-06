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


# compute_recursive_std = function(observations) {
#   num_observations = length(observations)
#   return(sqrt(
#     (cumsum(observations^2) / 1:num_observations
#      - (cumsum(observations) / 1:num_observations)^2)
#     * 1:num_observations / c(NA, 1:(num_observations - 1))
#   ))
# }


# test_compute_recursive_std = function() {
#   observations = c(1, 0.5, -1.2)
#
#   std = numeric(length(observations))
#
#   for (i in 1:length(observations)) {
#     std[i] = sqrt(var(observations[1:i]))
#   }
#   stopifnot(all.equal(std, compute_recursive_std(observations)))
# }


# test_compute_recursive_std()


apply_winsorisation = function(data, col, cutoff) {
  data[[col]] = pmin(data[[col]], cutoff)
  return(data)
}


test_apply_winsorisation = function() {
  df = apply_winsorisation(data.table(a = c(1, 5, 2), b = c(100, 100, 100)), "a", 4)
  stopifnot(max(df[, "a"]) == 4)
  stopifnot(max(df[, "b"]) == 100)
}


DataCleaner = R6Class(
  "DataCleaner",
  public = list(
    metric_col = NULL,
    user_id_col = NULL,
    NUM_EVENTS_COL = "num_events",
    num_events_cutoff = NULL,
    event_value_cap = NULL,
    initialize = function(historical_data,
                          metric_col,
                          user_id_col,
                          q = 0.99) {
      dt = data.table(historical_data)
      num_events_dt = dt[, .(num_events = .N), by = user_id_col]
      setnames(num_events_dt, "num_events", self$NUM_EVENTS_COL)
      num_events_cutoff = quantile(num_events_dt[, .SD[[self$NUM_EVENTS_COL]]], q)
      
      self$metric_col = metric_col
      self$user_id_col = user_id_col
      self$num_events_cutoff = 1 # num_events_cutoff
      self$event_value_cap = quantile(dt[num_events_dt[get(self$NUM_EVENTS_COL) <= num_events_cutoff], on =
                                           user_id_col][[metric_col]], q)
    },
    clean = function(data) {
      dt = data.table(data, key = self$user_id_col)
      num_events_dt = dt[, .(num_events = .N), by = eval(self$user_id_col)]
      setnames(num_events_dt, "num_events", self$NUM_EVENTS_COL)
      
      dt = dt[num_events_dt[get(self$NUM_EVENTS_COL) <= self$num_events_cutoff], on =
                self$user_id_col]
      return(dt)
      # return(apply_winsorisation(dt, self$metric_col, self$event_value_cap))
    }
  )
)


test = function() {
  test_apply_winsorisation()
  
  dt = data.table(
    user = c("u1", "u2", "u1", "u1", "u2"),
    metric = c(99.99, 15.99, 154.00, 100.00, 256.50)
  )
  dc = DataCleaner$new(dt, "metric", "user")
  cleaned_dt = dc$clean(dt)
  stopifnot(nrow(cleaned_dt) == 2)
  stopifnot(max(cleaned_dt[, "metric"]) < 256.50)
}


# test()
