# Checks that the boundary() refactor is behaviour-preserving.
#
# Run from the repository root:  Rscript tests/test_boundaries.R
#
# Each method exposes boundary() and inherits monitor() from SequentialTest.
# The expected values below are the formulas as they stood before that refactor,
# so a discrepancy means the refactor changed a test's decisions.
#
# Two things matter here. First, the boundary itself is asserted, not only the
# resulting detections: an earlier version of this file compared detections on a
# single trajectory, which let a wrong boundary through because a small error
# only rarely flips a decision. Second, every check runs at a horizon that the
# scheduled-look methods divide evenly and at one they do not, because the
# flooring of a fractional check time is exactly where they differ.

source("methods/bonferroni.R")
source("methods/gavi.R")
source("methods/gst.R")
source("methods/ld_obf.R")
source("methods/msprt.R")
source("methods/pyeast.R")
source("methods/yeast.R")

TOL <- 1e-10
ALPHA <- 0.05
SD <- sqrt(2)
HORIZONS <- c(500, 503)


check_boundary <- function(label, test, N, expected_index, expected_value) {
  b <- test$boundary(N)
  stopifnot(length(b$index) == length(expected_index))
  stopifnot(all(b$index == expected_index))
  stopifnot(length(b$value) == length(expected_value))
  stopifnot(max(abs(b$value - expected_value)) < TOL)
}


check_detections <- function(label, test, trajectory, expected_indicators) {
  actual <- test$monitor(trajectory)
  stopifnot(length(actual) == length(expected_indicators))
  stopifnot(all(actual == expected_indicators))
  cat(sprintf("  %s: %d monitored points, %d crossings\n",
              label, length(actual), sum(actual)))
}


for (N in HORIZONS) {
  cat(sprintf("horizon N = %d\n", N))
  set.seed(20260903)
  trajectory <- cumsum(rnorm(N, 0, SD))
  n <- 1:N
  standardised <- trajectory / SD / sqrt(n)

  # --- YEAST: a constant boundary at the forecast horizon.
  for (expected_n in c(400, 450, 500, 550, 600)) {
    test <- YEAST$new("YEAST", ALPHA, expected_n, SD)
    threshold <- qnorm(1 - ALPHA / 2) * sqrt(expected_n) * SD
    check_boundary("YEAST", test, N, n, rep(threshold, N))
    check_detections(sprintf("YEAST(n=%d)", expected_n), test, trajectory,
                     trajectory > threshold)
  }

  # --- pYEAST: piecewise constant thresholds.
  for (num_periods in c(7, 14)) {
    test <- pYEAST$new("pYEAST", ALPHA, round((1:num_periods) * (N / num_periods)), SD)
    cutoffs <- test$actual_cumulative_num_observations_by_period
    expected_boundary <- rep(NA_real_, N)
    expected_boundary[1:cutoffs[1]] <- test$thresholds[1]
    for (k in 2:num_periods) {
      expected_boundary[(cutoffs[k - 1] + 1):cutoffs[k]] <- test$thresholds[k]
    }
    check_boundary("pYEAST", test, N, n, expected_boundary)
    check_detections(sprintf("pYEAST%02d", num_periods), test, trajectory,
                     trajectory > expected_boundary)
  }

  # --- mSPRT: the pre-refactor comparison was on the mean scale.
  for (phi in c(100, 25, 1 / 0.3 ^ 2)) {
    alpha_star <- 2 * ALPHA
    radius <- SD / sqrt(n) * sqrt(log((phi + n) / (phi * alpha_star ^ 2)) * (phi + n) / n)
    test <- mSPRT$new("mSPRT", ALPHA, SD, phi)
    stopifnot(max(abs(test$mean_scale_boundary(N) - radius)) < TOL)
    check_boundary("mSPRT", test, N, n, n * radius)
    check_detections(sprintf("mSPRT(phi=%.2f)", phi), test, trajectory,
                     trajectory / n > radius)
  }

  # --- GAVI: the cumulative-sum scale, unchanged.
  for (n_tune in c(250, 500, 750, 10000)) {
    alpha_star <- 2 * ALPHA
    rho <- n_tune / (log(log(exp(1) * alpha_star ^ (-2))) - 2 * log(alpha_star))
    expected_boundary <- SD * sqrt((n + rho) * log((n + rho) / (rho * alpha_star ^ 2)))
    test <- GAVI$new("GAVI", ALPHA, SD, n_tune)
    stopifnot(abs(test$rho() - rho) < TOL)
    check_boundary("GAVI", test, N, n, expected_boundary)
    check_detections(sprintf("GAVI(%d)", n_tune), test, trajectory,
                     trajectory > expected_boundary)
  }

  # --- Bonferroni: standardised comparison at equally spaced checks. A check
  # time that is not a whole number of observations must be floored for the
  # index and for the sqrt scale alike.
  for (num_checks in c(14, 28, 42, 56, N)) {
    test <- Bonferroni$new("Bonferroni", ALPHA, SD, num_checks)
    check_times <- (N * seq_len(num_checks)) %/% num_checks
    critical_value <- qnorm(1 - ALPHA / num_checks)
    check_boundary("Bonferroni", test, N, check_times,
                   critical_value * SD * sqrt(check_times))
    check_detections(sprintf("Bonferroni(%d)", num_checks), test, trajectory,
                     standardised[check_times] > critical_value)
  }

  # --- GST: standardised comparison against the alpha-spending thresholds.
  for (params in list(c(1, 14, 14), c(3, 14, 14), c(1, 28, 14), c(1, 7, 14))) {
    test <- GST$new("GST", ALPHA, params[1], params[2], params[3], SD)
    check_times <- (N * seq_len(params[3])) %/% params[3]
    check_boundary("GST", test, N, check_times,
                   test$thresholds * SD * sqrt(check_times))
    check_detections(sprintf("GST(phi=%d, expected=%d, actual=%d)",
                             params[1], params[2], params[3]), test, trajectory,
                     standardised[check_times] > test$thresholds)
  }

  # --- Lan-DeMets OBF: piecewise constant critical values, applied after every
  # event.
  test <- LanDeMetsOBF$new("LanDeMetsOBF", ALPHA, SD)
  piecewise <- rep(test$critical_values,
                   each = ceiling(N / length(test$critical_values)))[n]
  check_boundary("LanDeMetsOBF", test, N, n, piecewise * SD * sqrt(n))
  check_detections("LanDeMetsOBF", test, trajectory, standardised > piecewise)
}

cat("test_boundaries.R: all checks passed\n")
