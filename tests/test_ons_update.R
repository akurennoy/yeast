# Formula-level tests for the ONS betting-fraction update used by SeqC2ST.
#
# Run from the repository root:  Rscript tests/test_ons_update.R
#
# The ONS step minimises the loss l_t(v) = -log(1 + v f_t), where f_t is the
# SIGNED payoff that enters the wealth update K_t = K_{t-1} (1 + v_t f_t).
# Its gradient at the current bet is z_t = -f_t / (1 + v_t f_t), the running
# sum of squares is A_t = 1 + sum_{s <= t} z_s^2, and the update is
#   v_{t+1} = Proj_[0, 1/2] ( v_t - (2 / (2 - log 3)) z_t / A_t ).
# The test compares hand-computed steps with the implementation.

source("methods/sec_c_2st_ons_qda.R")

TOL <- 1e-12
STEP <- 2 / (2 - log(3))

hand_step <- function(v, z_sumsq, f) {
  z <- -f / (1 + v * f)
  z_sumsq <- z_sumsq + z^2
  A <- 1 + z_sumsq
  v_new <- min(0.5, max(0.0, v - STEP * z / A))
  list(v = v_new, z_sumsq = z_sumsq, z = z, A = A)
}

# --- Step 1: from the initial state v = 0, a winning bet (f = +1) must
# increase the bet, and a losing bet (f = -1) must leave it at the lower end of
# the projection interval.
s <- hand_step(0, 0, +1)
stopifnot(abs(s$z + 1) < TOL)          # z = -1
stopifnot(abs(s$A - 2) < TOL)          # A = 1 + 1
stopifnot(abs(s$v - min(0.5, STEP / 2)) < TOL)
stopifnot(s$v > 0)

s_lose <- hand_step(0, 0, -1)
stopifnot(abs(s_lose$z - 1) < TOL)     # z = +1
stopifnot(abs(s_lose$v) < TOL)         # projected back to 0

# --- Step 2: a second winning bet from a positive bet level.
v1 <- s$v
s2 <- hand_step(v1, s$z_sumsq, +1)
stopifnot(abs(s2$z + 1 / (1 + v1)) < TOL)
stopifnot(s2$v >= v1 - TOL)

# --- The implementation must reproduce these steps. The QDA score is forced to
# a known value so that the payoff f_t is deterministic, and the wealth check is
# kept away from the stopping threshold.
make_test_object <- function(score_value) {
  obj <- SeqC2ST$new("test", alpha = 1e-12)
  obj$model <- local({
    fixed <- score_value
    list(score = function(x) fixed, update = function(x, y) invisible(NULL))
  })
  obj
}

for (f in c(+1, -1)) {
  for (w in c(+1, -1)) {
    g <- f * w                     # so that f_t = w_t * g_t = f
    obj <- make_test_object(g)
    v_before <- obj$v
    z_before <- obj$z_sumsq
    obj$step(0.0, w)
    expected <- hand_step(v_before, z_before, f)
    stopifnot(abs(obj$v - expected$v) < 1e-10)
    stopifnot(abs(obj$z_sumsq - expected$z_sumsq) < 1e-10)
  }
}

# --- The wealth process must use the signed payoff, not the raw score.
obj <- make_test_object(1)
obj$v <- 0.25
obj$step(0.0, +1)                  # f_t = +1, so wealth grows by a factor 1.25
stopifnot(abs(obj$K - 1.25) < 1e-12)

obj <- make_test_object(1)
obj$v <- 0.25
obj$step(0.0, -1)                  # f_t = -1, so wealth shrinks by a factor 0.75
stopifnot(abs(obj$K - 0.75) < 1e-12)

cat("test_ons_update.R: all checks passed\n")
