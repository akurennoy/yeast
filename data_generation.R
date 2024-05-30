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
