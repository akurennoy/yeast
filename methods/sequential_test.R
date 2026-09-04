library(R6)


# Contract for boundary(): returns list(index, value), where index holds the
# positions in the event stream at which the test is evaluated and value holds
# the alerting boundary at those positions, both on the cumulative-sum scale of
# the monitored trajectory. The default monitor() below is then exactly
# trajectory[index] > value, which keeps each boundary a single source of truth
# that tests/test_boundaries.R can pin to its closed form. Tests that are not
# boundary-based (SeqC2ST) override monitor() and leave boundary() undefined.
SequentialTest = R6Class("SequentialTest",
                         public = list(
                           name = NULL,
                           initialize = function(name) {
                             self$name = name
                           },
                           boundary = function(num_observations) {
                             stop("This method should be overridden by subclasses")
                           },
                           monitor = function(trajectory, assignment_indicators = NULL) {
                             b = self$boundary(length(trajectory))
                             return(trajectory[b$index] > b$value)
                           }
                         ))
