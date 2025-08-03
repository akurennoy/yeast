source("methods/sequential_test.R")
source("utils.R")


QDAStats <- R6Class("QDAStats", public = list(
  n = 0,
  mean = 0.0,
  var = 1.0,
  
  update = function(x) {
    self$n <- self$n + 1
    if (self$n == 1) {
      self$mean = x
      self$var = 1.0
    } else {
      delta = x - self$mean
      self$mean = self$mean + delta / self$n
      self$var = max(1e-6, ((self$n - 2) * self$var + delta * (x - self$mean)) / (self$n - 1))
    }
  }
))


SeqC2ST_QDA <- R6Class(
  "SeqC2ST_QDA",
  inherit = SequentialTest,
  public = list(
    significance_level = NULL,
    initial_wealth = NULL,
    initialize = function(name,
                          significance_level,
                          initial_wealth = 1.0) {
      super$initialize(name)
      self$significance_level = significance_level
      self$initial_wealth = initial_wealth
    },
    monitor = function(trajectory, assignment_indicators) {
      n = length(trajectory)
      stopping = logical(n)
      
      stats_test = QDAStats$new()
      stats_control = QDAStats$new()
      
      lam = 0.0
      a = 1.0
      K = self$initial_wealth
      
      for (i in 1:n) {
        # -- extracting the next observation (x) and the associated assignment (w)
        if (i == 1) {
          x = trajectory[1]
        } else {
          x = trajectory[i] - trajectory[i - 1]
        }
        w = 2 * assignment_indicators[i] - 1
        
        # -- computing the classifier score
        if (stats_test$n < 2 || stats_control$n < 2) {
          g = 0
        } else {
          g = 2 * (
            -0.5 * log(stats_test$var) - 0.5 * (x - stats_test$mean)^2 / stats_test$var
            > -0.5 * log(stats_control$var) - 0.5 * (x - stats_control$mean)^2 / stats_control$var
          ) - 1
        }
        
        # -- computing the payoff
        f = w * g
        
        # -- updating the wealth
        K = K * (1 + f * lam)
        
        # -- checking the detection condition
        stopping[i] = (K >= 1 / self$significance_level)
        
        # -- updating the stats
        if (w == 1) {
          stats_test$update(x)
        } else {
          stats_control$update(x)
        }
        
        # -- updating the bet
        denom = 1 - lam
        if (abs(denom) < 1e-12) {
          denom = sign(denom) * 1e-12
        }
        z = f / denom
        a = a + z^2
        lam <- max(min(lam - 2 / (2 - log(3)) * z / a, 0.5), -0.5)
      }

      return(stopping)
    }
    
  )
)

set.seed(2024)
print(round(measure_fdr(SeqC2ST_QDA$new("SeqC2ST_QDA", 0.05), 10, 500, 1000), 2))
