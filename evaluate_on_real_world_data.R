# In this script, we will evaluate the sequential tests on real data from
# a major e-commerce platform. That data consists of orders (checkout events).
# We will randomly assign users that placed those orders to either control and
# treatment and then apply the sequential tests to the resulting trajectory of
# the cumulative revenue difference between the two groups. We will then repeat
# this multiple times and measure the fraction of times where a detection
# occurs. Since the users are assigned randomly post-factum and no real
# treatment gets applied to the treatment group, this fraction should be close
# to the nominal significance level.

# The data is split in 2-week periods. This duration corresponds to a typical
# experiment duration on the e-commerce platform. We will apply the described
# procedure to those periods separately and then pool the detection outcomes
# to compute the detection fraction (false detection rate) across both the
# 2-week periods and the generated assignments. When applying the sequential
# testing methods we will set the increment standard deviation to its estimate
# obtained from the previous 2-week period. GAVI and mSPRT will be additionally
# evaluated in the mode in which the increment standard deviation is estimated
# during the monitoring process (using the data arrived so far). The expected
# number of observations required by SST and pSST methods will also be
# estimated on the preceeding 2-week period. There are 25 2-week periods in
# total and we will use 23 of them for evaluation (the first one is used to
# compute the 99.9% percentile of the total order value which serves as as a cutoff
# for data filtering that gets applied to all subsequent data; the remaining 24
# fortnight periods are used in such a way that the input parameters mentioned
# above are estimated on the preceeding period, hence we leave the first of the
# 24 periods out, which leaves us with 23 periods for the evaluation).
# For each of the 2-week periods we will generate the random assignments 100
# times and will have 2400 cumulative revenue trajectories to "monitor" and
# measure the false detection rate.

library(arrow)
library(sandwich)

source("data_generation.R")
source("methods/bonferroni.R")
source("methods/caa.R")
source("methods/gavi.R")
source("methods/gst.R")
source("methods/msprt.R")
source("methods/psst.R")
source("methods/sst.R")
source("utils.R")


# Global Settings


SIGNIFICANCE_LEVEL = 0.05
NUM_ASSIGNMENT_REPLICATIONS = 1000 # for each 2-week period
DATA_DIRECTORY = "./checkouts"
INCREMENT_STD_NUM_BURN_IN_STEPS = 100 # when estimating the increment standard
# deviation recursively from the in-period
# this number of initial estimates will
# be substituted with the estimate
# computed from the previous period data


# Definitions


initialise_continuous_methods = function(robust_increment_std,
                                         non_robust_increment_std,
                                         expected_num_observations) {
  return(
    list(
      # -- SST
      SST = SST$new(
        "SST",
        SIGNIFICANCE_LEVEL,
        sum(expected_num_observations),
        robust_increment_std
      ),
      # SSTnr = SST$new(
      #   "SST-non-robust",
      #   SIGNIFICANCE_LEVEL,
      #   expected_num_observations,
      #   non_robust_increment_std
      # ),
      # -- pSST
      pSST14 = pSST$new(
        "pSST14",
        SIGNIFICANCE_LEVEL,
        cumsum(expected_num_observations),
        increment_std
      ),
      # -- mSPRT
      mSPRTphi100 = mSPRT$new("mSPRT100", SIGNIFICANCE_LEVEL, robust_increment_std, 100),
      mSPRTphi025 = mSPRT$new("mSPRT025", SIGNIFICANCE_LEVEL, robust_increment_std, 25),
      mSPRTphi011 = mSPRT$new(
        "mSPRT011",
        SIGNIFICANCE_LEVEL,
        robust_increment_std,
        1 / 0.3 ^
          2
      ),
      # mSPRTphi100nr = mSPRT$new(
      #   "mSPRT100-non-robust",
      #   SIGNIFICANCE_LEVEL,
      #   non_robust_increment_std,
      #   100
      # ),
      # mSPRTphi025nr = mSPRT$new(
      #   "mSPRT025-non-robust",
      #   SIGNIFICANCE_LEVEL,
      #   non_robust_increment_std,
      #   25
      # ),
      # mSPRTphi011nr = mSPRT$new(
      #   "mSPRT011-non-robust",
      #   SIGNIFICANCE_LEVEL,
      #   non_robust_increment_std,
      #   1 / 0.3 ^ 2
      # ),
      # -- GAVI
      GAVI250 = GAVI$new("GAVI250", SIGNIFICANCE_LEVEL, robust_increment_std, 250),
      GAVI500 = GAVI$new("GAVI500", SIGNIFICANCE_LEVEL, robust_increment_std, 500),
      GAVI750 = GAVI$new("GAVI750", SIGNIFICANCE_LEVEL, robust_increment_std, 750),
      # GAVI250nr = GAVI$new("GAVI250-non-robust", SIGNIFICANCE_LEVEL, non_robust_increment_std, 250),
      # GAVI500nr = GAVI$new("GAVI500-non-robust", SIGNIFICANCE_LEVEL, non_robust_increment_std, 500),
      # GAVI750nr = GAVI$new("GAVI750-non-robuts", SIGNIFICANCE_LEVEL, non_robust_increment_std, 750),
      # -- CAA (Statsig)
      CAA = CAA$new("CAA", SIGNIFICANCE_LEVEL, robust_increment_std),
      # CAAnr = CAA$new("CAA-nr", SIGNIFICANCE_LEVEL, non_robust_increment_std),
      Classical = Bonferroni$new("Classical", SIGNIFICANCE_LEVEL, robust_increment_std, 1)
      # Classicalnr = Bonferroni$new("Classical-non-robust", SIGNIFICANCE_LEVEL, non_robust_increment_std, 1)
    )
  )
}


estimate_increment_std = function(data,
                                  metric_col,
                                  user_id_col,
                                  sample_size = 10000) {
  dg = DataGeneratorFromRealEvents$new(data, metric_col, user_id_col, 1)
  dt = dg$real_events_data_table[, x := gmv_euro * (1 - 2 * a1)]
  dt = dt[, occurred_at := as.POSIXct(occurred_at, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")]
  dt = dt[, date := as.Date(occurred_at)][, hour := hour(occurred_at)]
  grouped = dt[, .(x = sum(x)), by = .(user_id, date, hour)]
  # TODO: take the user id column name from the parameter instead of using
  #.      a hardcoded value in the by parameter
  users = data.table(user_id = unique(grouped[[user_id_col]]))[sample(.N, sample_size)]
  subsampled_grouped = grouped[users, on = "user_id"]
  model <- lm(x ~ 1, data = subsampled_grouped)
  v = (
    vcovCL(model, cluster = subsampled_grouped[[user_id_col]], type = "HC3")[[1]]
    * nrow(subsampled_grouped) / nrow(grouped)
  )
  
  return(sqrt(v * nrow(grouped) ** 2 / nrow(dt)))
}


# Conducting the Evaluation


# Each input file contains orders placed within a 2-week period. The order
# of the file names corresponds to the chronological order of the 2-week period.
# There are 25 files in total covering the time frame between 2023-05-23 and
# 2024-05-12 (without any gaps or overlaps).
input_files = sort(list.files(DATA_DIRECTORY))

USER_ID_COL = "user_id"
METRIC_COL = "gmv_euro"
DTTM_COL = "occurred_at"
set.seed(2025)
result = NULL
data_cleaner = DataCleaner$new(read_parquet(sprintf("%s/%s", DATA_DIRECTORY, input_files[1])), METRIC_COL, USER_ID_COL, DTTM_COL, q =
                                 0.999)
print(
  sprintf(
    "Will remove users with more than %.2f of total order value.",
    data_cleaner$event_value_cutoff
  )
)

raw_preceeding_data = read_parquet(sprintf("%s/%s", DATA_DIRECTORY, input_files[2]))
preceeding_data = data_cleaner$clean(raw_preceeding_data)
num_unique_users_before = length(unique(raw_preceeding_data[[USER_ID_COL]]))
num_unique_users_after = length(unique(preceeding_data[[USER_ID_COL]]))
print(
  sprintf(
    "Removed %i users out of %i",
    num_unique_users_before - num_unique_users_after,
    num_unique_users_before
  )
)

for (i in 3:length(input_files)) {
  increment_std = estimate_increment_std(preceeding_data, METRIC_COL, USER_ID_COL)
  print(sprintf(
    "Increment std estimate using preceeding data = %.2f",
    increment_std
  ))
  expected_num_observations = setorder(preceeding_data[, date := as.Date(as.POSIXct(occurred_at, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"))][, .(num_obs = .N), by = date], "date")[["num_obs"]]
  expected_num_observations = nrow(preceeding_data)
  
  data = data_cleaner$clean(read_parquet(sprintf("%s/%s", DATA_DIRECTORY, input_files[i])))
  data_generator = DataGeneratorFromRealEvents$new(data, METRIC_COL, USER_ID_COL)
  aggregator = Aggregator$new()
  for (r in 1:NUM_ASSIGNMENT_REPLICATIONS) {
    if (r %% 100 == 0) {
      print(sprintf("Replication # %.03d", r))
    }
    trajectory = data_generator$generate_cumulative_difference_trajectory()
    # the trajectory of the cumulative difference in the revenue between
    # control and treatment
    
    continuous_methods = initialise_continuous_methods(increment_std,
                                                       in_period_increment_std,
                                                       expected_num_observations)
    
    for (statistical_test in continuous_methods) {
      detection_indicators = statistical_test$monitor(trajectory)
      # -- continuous monitoring mode
      aggregator$update(0.0,
                        "stream",
                        statistical_test$name,
                        any(detection_indicators),
                        0)
    }
  }
  
  result = aggregator$get_result()
  
  preceeding_data = data
  
  break
}

print(result)
print(result)