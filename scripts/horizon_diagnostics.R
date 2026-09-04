# Secondary diagnostics isolating the effect of the planned-versus-realised
# horizon mismatch in the Online Retail experiment.
#
#   Rscript scripts/horizon_diagnostics.R
#
# The primary experiment (evaluate_on_online_retail.R) builds the YEAST boundary
# from the event-count forecast Ntilde taken from the preceding six months and
# then monitors every event of the following six months, of which there are
# M_T. Because M_T > Ntilde on this series, the primary run sits outside the
# setting of Theorem 1, which assumes the horizon used to build the boundary is
# the horizon over which the maximum is taken.
#
# Three configurations are run on identical replications so that the mismatch is
# the only thing that varies:
#
#   primary    boundary from Ntilde, monitor all M_T events   (the reported design)
#   truncated  boundary from Ntilde, monitor min(M_T, Ntilde) (theorem-aligned)
#   oracle     boundary from M_T,    monitor all M_T events   (unattainable in
#                                                              practice: M_T is
#                                                              not known when the
#                                                              boundary is fixed)
#
# Only YEAST is run. The theorem being probed is YEAST's, and the competitors
# would add hours of runtime (almost all of it SeqC2ST-QDA) without bearing on
# the question.

if (!require("data.table", character.only = TRUE)) {
  install.packages("data.table")
  library("data.table", character.only = TRUE)
}
library(parallel)
library(sandwich)

source("data_generation.R")
source("methods/yeast.R")
source("utils.R")


SIGNIFICANCE_LEVEL = 0.05
NUM_ASSIGNMENT_REPLICATIONS = as.integer(
  Sys.getenv("YEAST_NUM_REPLICATIONS", unset = "100000"))
PROCESSED_DATA_DIRECTORY = "./online-retail-orders/"
OUTPUT_FILE = "horizon_diagnostics.csv"
EFFECTS = c(0.00, 0.05, 0.10, 0.20, 0.50)

USER_ID_COL = "user_id"
METRIC_COL = "order_value"
DTTM_COL = "occurred_at"

input_files = sort(list.files(PROCESSED_DATA_DIRECTORY))
stopifnot(length(input_files) >= 2)

# Reproduced from evaluate_on_online_retail.R so that the scale input, the
# capping and the forecast are identical to the primary run.
estimate_increment_std = function(data, metric_col, user_id_col) {
  dg = DataGeneratorFromRealEvents$new(data, metric_col, user_id_col, 1)
  dt = dg$real_events_data_table[, x := get(metric_col) * (1 - 2 * a1)]
  dt = dt[, year_month := format(get(DTTM_COL), "%Y-%m")]
  grouped = dt[, .(x = sum(x)), by = c(user_id_col, "year_month")]
  model = lm(x ~ 1, data = grouped)
  v = vcovCL(model, cluster = grouped[[user_id_col]], type = "HC3")[[1]]
  return(sqrt(v * nrow(grouped) ^ 2 / nrow(dt)))
}

num_cores = min(13, detectCores() - 1)
num_assignment_replications = rep(NUM_ASSIGNMENT_REPLICATIONS %/% num_cores,
                                  num_cores)
num_assignment_replications[num_cores] = (
  NUM_ASSIGNMENT_REPLICATIONS
  - sum(num_assignment_replications[-num_cores]))

cl = makeCluster(num_cores)
clusterExport(cl, varlist = c(
  "PROCESSED_DATA_DIRECTORY", "METRIC_COL", "USER_ID_COL", "DTTM_COL",
  "input_files", "DataCleaner", "estimate_increment_std",
  "generate_assignments", "DataFromRealEvents", "DataGeneratorFromRealEvents",
  "Aggregator", "SequentialTest", "YEAST", "get_savings",
  "num_assignment_replications", "SIGNIFICANCE_LEVEL", "EFFECTS"))

process_file = function(i) {
  library(data.table)
  library(sandwich)
  set.seed(i)

  raw_preceeding = readRDS(sprintf("%s/%s", PROCESSED_DATA_DIRECTORY, input_files[1]))
  data_cleaner = DataCleaner$new(raw_preceeding, METRIC_COL, USER_ID_COL,
                                 DTTM_COL, q = 0.999)
  preceeding = data_cleaner$clean(raw_preceeding)
  increment_std = estimate_increment_std(preceeding, METRIC_COL, USER_ID_COL)

  n_tilde = nrow(preceeding)
  data = data_cleaner$clean(readRDS(sprintf("%s/%s", PROCESSED_DATA_DIRECTORY,
                                            input_files[2])))
  m_t = nrow(data)

  generator = DataGeneratorFromRealEvents$new(data, METRIC_COL, USER_ID_COL)
  aggregator = Aggregator$new()

  # The truncated arm monitors a prefix of the same trajectory; the oracle arm
  # differs only in the horizon fed to the boundary.
  primary = YEAST$new("YEAST", SIGNIFICANCE_LEVEL, n_tilde, increment_std)
  oracle = YEAST$new("YEAST", SIGNIFICANCE_LEVEL, m_t, increment_std)
  truncation_point = min(m_t, n_tilde)

  for (relative_effect in EFFECTS) {
    for (r in 1:num_assignment_replications[i]) {
      trajectory = generator$generate_cumulative_difference_trajectory(
        -relative_effect)$trajectory

      for (arm in c("primary", "truncated", "oracle")) {
        if (arm == "primary") {
          path = trajectory
          test = primary
        } else if (arm == "truncated") {
          path = trajectory[1:truncation_point]
          test = primary
        } else {
          path = trajectory
          test = oracle
        }
        detections = test$monitor(path, NULL)
        aggregator$update(relative_effect, arm, "YEAST", any(detections),
                          get_savings(detections, 1:length(path), length(path)))
      }
    }
  }

  result = aggregator$get_result()
  result$n_tilde = n_tilde
  result$m_t = m_t
  return(result)
}

results = parLapply(cl, 1:num_cores, process_file)
stopCluster(cl)

combined = rbindlist(results)
totals = combined[, .(num_trials = sum(num_trials),
                      num_detections = sum(num_detections),
                      total_savings = sum(total_savings),
                      n_tilde = max(n_tilde),
                      m_t = max(m_t)),
                  by = .(effect, mode, method)]
totals[, `:=`(detection_rate = num_detections / num_trials,
              average_savings = total_savings / num_trials)]
setorder(totals, effect, mode)
write.csv(totals, OUTPUT_FILE, row.names = FALSE)
cat(sprintf("wrote %s\n", OUTPUT_FILE))
print(dcast(totals, mode ~ effect, value.var = "detection_rate"), digits = 5)
