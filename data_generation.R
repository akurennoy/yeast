library(R6)


EventValueGenerator = R6Class(
  "EventValueGenerator",
  public = list(
    generate_control_and_treatment_event_values = function(num_events, relative_effect) {
      stop("This method should be overridden by subclasses")
    }
  )
)


NormalEventValueGenerator = R6Class(
  "normal",
  inherit = EventValueGenerator,
  public = list(
    mean = NULL,
    sd = NULL,
    initialize = function(mean, sd) {
      self$mean = mean
      self$sd = sd
    },
    generate_control_and_treatment_event_values = function(num_events, relative_effect) {
      control_events = rnorm(num_events, self$mean, self$sd)
      treatment_events = rnorm(num_events, self$mean * (1 + relative_effect), self$sd)
      return(list(control = control_events, treatment = treatment_events))
    }
  )
)


ShiftedStudentEventValueGenerator = R6Class(
  "student",
  inherit = EventValueGenerator,
  public = list(
    mean = NULL,
    df = NULL,
    initialize = function(mean, df) {
      self$mean = mean
      self$df = df
    },
    generate_control_and_treatment_event_values = function(num_events, relative_effect) {
      control_events = rt(num_events, self$df) + self$mean
      treatment_events = (rt(num_events, self$df) + self$mean * (1 + relative_effect))
      return(list(control = control_events, treatment = treatment_events))
    }
  )
)


GammaEventValueGenerator = R6Class(
  "gamma",
  inherit = EventValueGenerator,
  public = list(
    k = NULL,
    theta = NULL,
    initialize = function(k, theta) {
      self$k = k
      self$theta = theta
    },
    generate_control_and_treatment_event_values = function(num_events, relative_effect) {
      control_events = rgamma(num_events, shape = self$k, rate = 1 / self$theta)
      treatment_events = rgamma(num_events,
                                shape = self$k,
                                rate = 1 / (self$theta * (1 + relative_effect)))
      return(list(control = control_events, treatment = treatment_events))
    }
  )
)


DataGenerator = R6Class(
  "DataGenerator",
  public = list(
    event_value_generator = NULL,
    initialize = function(event_value_generator) {
      self$event_value_generator = event_value_generator
    },
    generate_cumulative_difference_trajectory = function(num_observations, relative_effect) {
      event_values = self$event_value_generator$generate_control_and_treatment_event_values(num_observations, relative_effect)
      return(cumsum(event_values$treatment)
             - cumsum(event_values$control))
    }
  )
)


generate_assignments = function(
    unique_user_ids, batch_size, user_id_col, assignment_col_names
) {
  assignments = data.table(cbind(
    unique_user_ids,
    matrix(
      sample(
        c(0, 1),
        # 0 - control, 1 - treatment
        nrow(unique_user_ids) * batch_size,
        replace = TRUE
      ),
      nrow = nrow(unique_user_ids),
      ncol = batch_size
    )
  ))
  setnames(assignments, c(
    user_id_col,
    assignment_col_names
  ))
  setkeyv(assignments, user_id_col)
  return(assignments)
}

DataGeneratorFromRealEvents = R6Class(
  "DataGeneratorFromRealEvents",
  public = list(
    TREATMENT_INDICATOR_COL = "treatment_indicator",
    DTTM_COL = "occurred_at",
    metric_col = NULL,
    user_id_col = NULL,
    
    ASSIGNMENT_REPLICATION_COL_PREFIX = "a",
    assignment_batch_size = NULL,
    assignment_col_names = NULL,

    unique_user_ids = NULL,    
    real_events_data_table = NULL,
    
    next_replication_num = 1,
    
    initialize = function(
      real_events, metric_col, user_id_col, assignment_batch_size=100
    ) {
      real_events_data_table = data.table(real_events, key=user_id_col)
      unique_user_ids = unique(real_events_data_table[, ..user_id_col])
      assignment_col_names = paste0(
        self$ASSIGNMENT_REPLICATION_COL_PREFIX,
        1:assignment_batch_size
      )
      assignments = generate_assignments(
        unique_user_ids, assignment_batch_size, user_id_col, assignment_col_names
      )
      self$assignment_batch_size = assignment_batch_size
      self$assignment_col_names = assignment_col_names
      self$real_events_data_table = setorderv(real_events_data_table[assignments], self$DTTM_COL)
      self$unique_user_ids = unique_user_ids
      self$metric_col = metric_col
      self$user_id_col = user_id_col
    },
    
    replicate_treatment_indicators = function() {
      assignments = generate_assignments(
        self$unique_user_ids, self$assignment_batch_size, self$user_id_col, self$assignment_col_names
      )
      self$real_events_data_table[, (self$assignment_col_names) := NULL]
      setkeyv(self$real_events_data_table, self$user_id_col)
      self$real_events_data_table = setorderv(self$real_events_data_table[assignments], self$DTTM_COL)
    },
    
    generate_cumulative_difference_trajectory = function() {
      if(self$next_replication_num > self$assignment_batch_size) {
        self$replicate_treatment_indicators()
        self$next_replication_num = 1
      }
      assignment_col = paste0(self$ASSIGNMENT_REPLICATION_COL_PREFIX,
                              self$next_replication_num)
      w = self$real_events_data_table[, .SD[[assignment_col]]]
      y = self$real_events_data_table[, .SD[[self$metric_col]]]
      x = y * (1 - 2 * w)
      s = cumsum(x)
      
      self$next_replication_num = self$next_replication_num + 1
      return(s)
    }
  )
)


test = function() {
  df = data.table(
    user_id = c("u1", "u2", "u1", "u1", "u3"),
    occurred_at = c(
      "2023-06-17T15:10:40.146Z",
      "2023-06-17T15:11:40.146Z",
      "2023-06-17T16:10:41.146Z",
      "2023-06-17T16:10:40.146Z",
      "2023-06-16T00:01:00.146Z"
    ),
    gmv_euro = c(1, 2, 1, 3, 1.5)
  )
  set.seed(0)
  generator = DataGeneratorFromRealEvents$new(df, "gmv_euro", "user_id", 2)
  times = generator$real_events_data_table[, "occurred_at"]
  setorder(times, "occurred_at")
  stopifnot(all(generator$real_events_data_table[, "occurred_at"] == times))
  stopifnot(all(
    generator$real_events_data_table[, "user_id"] == c("u3", "u1", "u2", "u1", "u1", "u1")
  ))
  stopifnot(all(
    generator$generate_cumulative_difference_trajectory() == c(-1.5, -2.5, -0.5, -3.5, -4.5)
  ))
  stopifnot(all(
    generator$generate_cumulative_difference_trajectory() == c(-1.5, -0.5, 1.5, 4.5, 5.5)
  ))
  stopifnot(all(
    generator$generate_cumulative_difference_trajectory() == c(1.5, 2.5, 4.5, 7.5, 8.5)
  ))
  stopifnot(all(
    generator$generate_cumulative_difference_trajectory() == c(1.5, 0.5, -1.5, -4.5, -5.5)
  ))
}


# test()
