library(progress)

source("methods/sequential_test.R")
source("utils.R")
source("methods/yeast.R")


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


OnlineQDA <- R6Class("OnlineQDA",
                     public = list(
                       stats_pos = NULL,
                       stats_neg = NULL,

                       initialize = function(lr = 0.1) {
                         self$stats_pos = QDAStats$new()
                         self$stats_neg = QDAStats$new()
                       },

                       predict = function(x) {
                         (
                           -0.5 * log(self$stats_pos$var) - 0.5 * (x - self$stats_pos$mean)^2 / self$stats_pos$var
                           > -0.5 * log(self$stats_neg$var) - 0.5 * (x - self$stats_neg$mean)^2 / self$stats_neg$var
                         )
                       },

                       score = function(x) {
                         p = self$predict(x)
                         2 * p - 1
                       },

                       update = function(x, y) {
                         if (y == 1) {
                           self$stats_pos$update(x)
                         }
                         else {
                           self$stats_neg$update(x)
                         }
                       }
                     )
)


SeqC2ST <- R6Class("SeqC2ST",
                   inherit = SequentialTest,
                   
                   public = list(
                     alpha = NULL,
                     lr = 0.1,
                     K = NULL,
                     model = NULL,
                     z_sumsq = NULL,
                     v = NULL,
                     
                     
                     initialize = function(name, alpha = 0.05, lr = 0.1) {
                       super$initialize(name)
                       self$alpha <- alpha
                       self$lr <- lr
                       self$K <- 1.0
                       self$model <- OnlineQDA$new()
                       self$z_sumsq <- 0.0
                       self$v <- 0.0
                     },
                     
                     step = function(z_t, w_t) {
                       g_t <- self$model$score(z_t)
                       f_t <- w_t * g_t
                       self$K <- min(1e6, self$K + f_t * self$v * self$K)
                       
                       if (self$K >= 1 / self$alpha) {
                         return(TRUE)
                       }
                       
                       y01 <- (w_t + 1) / 2
                       self$model$update(z_t, y01)
                       
                       # z_ons <- g_t / max(1e-12, 1 - g_t * self$v)
                       z_ons <- g_t / max(1e-12, 1 - self$v)
                       self$z_sumsq <- self$z_sumsq + z_ons^2
                       A <- 1 + self$z_sumsq
                       step_size <- 2 / (2 - log(3))# / 2
                       v_new <- self$v - step_size * z_ons / A
                       self$v <- min(0.5, max(-0.0, v_new))
                       
                       return(FALSE)
                     },
                     
                     monitor = function(trajectory, assignment_indicators) {
                       self$model <- OnlineQDA$new()
                       self$K = 1.0
                       self$z_sumsq = 0.0
                       self$v = 0.0
                       
                       n = length(trajectory)
                       stopping = logical(n)
                       
                       for (i in 1:n) {
                         # -- extracting the next observation (x) and the associated assignment (w)
                         w = 2 * assignment_indicators[i] - 1
                         if (i == 1) {
                           x = trajectory[1] / (-w)
                         } else {
                           x = (trajectory[i] - trajectory[i - 1]) / (-w)
                         }
                         
                         stopping[i] = self$step(x, w)
                       }
                       return(stopping)
                     }
                   )
)


# set.seed(2024)
# print(round(measure_fdr(SeqC2ST$new("SeqC2ST"), 1, 500, 1000), 2))
