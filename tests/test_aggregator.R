# Checks the hashed Aggregator against a direct reimplementation of the
# original data-frame accumulator.
#
# Run from the repository root:  Rscript tests/test_aggregator.R

source("utils.R")

set.seed(4242)

# The original implementation, kept here only as the reference to compare with.
reference_aggregate = function(records) {
  df = NULL
  for (r in records) {
    row_index = which(df$effect == r$effect & df$mode == r$mode
                      & df$method == r$method)
    if (length(row_index) == 0) {
      df = rbind(df, data.frame(effect = r$effect, mode = r$mode,
                                method = r$method, num_trials = 1,
                                num_detections = as.numeric(r$detected),
                                total_savings = r$savings))
    } else {
      stopifnot(length(row_index) == 1)
      df[row_index, "num_trials"] = df[row_index, "num_trials"] + 1
      df[row_index, "num_detections"] = df[row_index, "num_detections"] + r$detected
      df[row_index, "total_savings"] = df[row_index, "total_savings"] + r$savings
    }
  }
  df$detection_rate = df$num_detections / df$num_trials
  df$average_savings = df$total_savings / df$num_trials
  return(df)
}

effects = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.05, 0.5)
modes = c("stream", "14", "28", "42", "56", "non-sequential")
methods = c("YEAST", "mSPRTphi11", "GST", "SeqC2ST_QDA", "Bonferroni")

records = lapply(1:4000, function(i) list(
  effect = sample(effects, 1),
  mode = sample(modes, 1),
  method = sample(methods, 1),
  detected = sample(c(TRUE, FALSE), 1),
  savings = round(runif(1), 5)
))

aggregator = Aggregator$new()
for (r in records) {
  aggregator$update(r$effect, r$mode, r$method, r$detected, r$savings)
}
actual = as.data.frame(aggregator$get_result())
expected = reference_aggregate(records)

stopifnot(identical(dim(actual), dim(expected)))
stopifnot(identical(names(actual), names(expected)))
# Row order must also match: cells appear in first-encounter order.
for (column in names(expected)) {
  if (is.numeric(expected[[column]])) {
    stopifnot(max(abs(actual[[column]] - expected[[column]])) < 1e-12)
  } else {
    stopifnot(identical(as.character(actual[[column]]),
                        as.character(expected[[column]])))
  }
}

stopifnot(sum(actual$num_trials) == length(records))

cat(sprintf("test_aggregator.R: %d cells, %d updates, all checks passed\n",
            nrow(actual), length(records)))
