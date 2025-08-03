library(R6)


SequentialTest = R6Class("SequentialTest",
                         public = list(
                           name = NULL,
                           initialize = function(name) {
                             self$name = name
                           },
                           monitor = function(trajectory, assignment_indicators) {
                             stop("This method should be overridden by subclasses")
                           }
                         ))
