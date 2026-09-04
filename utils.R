if(!require("data.table", character.only = TRUE)) {
  install.packages("data.table")
  library("data.table", character.only = TRUE)
}
library(R6)

source("methods/sequential_test.R")


measure_detection_rate = function(sequential_test,
                                  increment_std,
                                  num_observations,
                                  num_replications,
                                  effect_size = 0.0) {
  num_detections = 0
  for (r in 1:num_replications) {
    assignments = rbinom(num_observations, 1, 0.5)
    increments = rnorm(num_observations, effect_size * assignments, increment_std)
    if (any(sequential_test$monitor(cumsum(increments), assignments))) {
      num_detections = num_detections + 1
    }
  }
  return(as.numeric(num_detections) / num_replications)
}


# Wilson score interval. Preferred over the Wald interval because the Monte
# Carlo detection rates reach 0 and 1, where the Wald interval degenerates to a
# zero-width interval.
wilson_interval = function(num_successes, num_trials, confidence_level = 0.95) {
  z = qnorm(1 - (1 - confidence_level) / 2)
  p = num_successes / num_trials
  denominator = 1 + z ^ 2 / num_trials
  centre = (p + z ^ 2 / (2 * num_trials)) / denominator
  half_width = (z / denominator) * sqrt(p * (1 - p) / num_trials
                                        + z ^ 2 / (4 * num_trials ^ 2))
  return(list(lower = pmax(0, centre - half_width),
              upper = pmin(1, centre + half_width)))
}


# Whether a realised Type-I error is inflated beyond the nominal level by more
# than Monte Carlo error, used to decide which rows may be presented as powerful.
# Comparing the point estimate against the nominal level directly would turn a
# fraction of a standard error into an exclusion: with 100,000 replications at
# alpha = 0.05 the standard error is about 0.0007, so a realised 0.0501 is noise
# and not evidence of an inflated size.
is_size_inflated = function(num_successes,
                            num_trials,
                            significance_level,
                            confidence_level = 0.95) {
  lower = wilson_interval(num_successes, num_trials, confidence_level)$lower
  return(lower > significance_level)
}


# The observation indices at which the j-th of k equally spaced checks falls.
# Exact integer division, for the same reason as in the boundary methods:
# num_observations * seq(1/k, 1, 1/k) truncates to one observation early
# whenever the product is an integer that is not exactly representable, which
# at k = 14, 28, 56 drops the terminal check entirely.
discrete_check_times = function(num_observations, num_checks) {
  return((num_observations * seq_len(num_checks)) %/% num_checks)
}


# Fraction of the planned event budget left unused when the test stops. A
# replication in which no alert is raised contributes exactly 0, and the average
# reported in the tables is taken over ALL replications, detecting or not. This
# is a saving in events, not in calendar time: the two coincide only if the
# arrival rate is constant over the monitoring window.
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
    # Counters are held in a hashed environment keyed by (effect, mode, method)
    # and the table is assembled only in get_result(). The previous
    # implementation rescanned a growing data frame on every update, which is
    # linear in the number of cells and dominates the runtime of the full
    # experiment (tens of millions of updates).
    counts = NULL,
    meta = NULL,
    insertion_order = NULL,
    initialize = function() {
      self$counts = new.env(hash = TRUE, parent = emptyenv())
      self$meta = new.env(hash = TRUE, parent = emptyenv())
      self$insertion_order = character()
    },
    update = function(relative_effect,
                      monitoring_mode,
                      method_name,
                      detection_indicator,
                      savings) {
      key = paste(relative_effect, monitoring_mode, method_name, sep = "\r")
      cell = self$counts[[key]]
      increment = c(1, as.numeric(detection_indicator), savings)
      if (is.null(cell)) {
        self$counts[[key]] = increment
        self$meta[[key]] = list(effect = relative_effect,
                                mode = monitoring_mode,
                                method = method_name)
        self$insertion_order = c(self$insertion_order, key)
      } else {
        self$counts[[key]] = cell + increment
      }
    },
    get_result = function() {
      keys = self$insertion_order
      cells = lapply(keys, function(key) self$counts[[key]])
      metas = lapply(keys, function(key) self$meta[[key]])
      result = data.table(
        effect = vapply(metas, function(m) as.numeric(m$effect), numeric(1)),
        mode = vapply(metas, function(m) as.character(m$mode), character(1)),
        method = vapply(metas, function(m) as.character(m$method), character(1)),
        num_trials = vapply(cells, function(x) x[1], numeric(1)),
        num_detections = vapply(cells, function(x) x[2], numeric(1)),
        total_savings = vapply(cells, function(x) x[3], numeric(1))
      )
      setnames(result, c(self$EFFECT_COL, self$MODE_COL, self$METHOD_COL,
                         self$NUM_TRIALS_COL, self$NUM_DETECTIONS_COL,
                         self$TOTAL_SAVINGS_COL))
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


# test()
