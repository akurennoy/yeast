# In this script, we will evaluate the sequential tests on real data from
# a major e-commerce platform. That data consists of orders (checkout events).
# We will randomly assign users that placed those orders to either control and
# treatment and then apply the sequential tests to the resulting trajectory of
# the cumulative revenue difference between the two groups. We will then repeat
# this multiple times and measure the number of detections. 
# Since the users are assigned randomly post-factum and no real
# treatment gets applied to the treatment group, the relative frequency of
# detection should be close to the nominal significance level.


library(arrow)
library(parallel)
library(sandwich)
library(stringr)

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
NUM_ASSIGNMENT_REPLICATIONS = 10000 # for each evaluation period
DATA_DIRECTORY = "./checkouts"
OUTPUT_DIRECTORY = "./results"


# Definitions


initialise_continuous_methods = function(robust_increment_std,
                                         non_robust_increment_std,
                                         expected_num_observations,
                                         actual_num_observations) {
  return(list(
    # -- SST
    SST = SST$new(
      "SST",
      SIGNIFICANCE_LEVEL,
      sum(expected_num_observations),
      robust_increment_std
    ),
    SSTnr = SST$new(
      "SST-non-robust",
      SIGNIFICANCE_LEVEL,
      sum(expected_num_observations),
      non_robust_increment_std
    ),
    # -- pSST
    pSST14 = pSST$new(
      "pSST14",
      SIGNIFICANCE_LEVEL,
      cumsum(expected_num_observations),
      robust_increment_std,
      cumsum(actual_num_observations)
    ),
    pSST14nr = pSST$new(
      "pSST14-non-robust",
      SIGNIFICANCE_LEVEL,
      cumsum(expected_num_observations),
      non_robust_increment_std,
      cumsum(actual_num_observations)
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
    mSPRTphi100nr = mSPRT$new(
      "mSPRT100-non-robust",
      SIGNIFICANCE_LEVEL,
      non_robust_increment_std,
      100
    ),
    mSPRTphi025nr = mSPRT$new(
      "mSPRT025-non-robust",
      SIGNIFICANCE_LEVEL,
      non_robust_increment_std,
      25
    ),
    mSPRTphi011nr = mSPRT$new(
      "mSPRT011-non-robust",
      SIGNIFICANCE_LEVEL,
      non_robust_increment_std,
      1 / 0.3 ^ 2
    ),
    # -- GAVI
    GAVI250 = GAVI$new("GAVI250", SIGNIFICANCE_LEVEL, robust_increment_std, 250),
    GAVI500 = GAVI$new("GAVI500", SIGNIFICANCE_LEVEL, robust_increment_std, 500),
    GAVI750 = GAVI$new("GAVI750", SIGNIFICANCE_LEVEL, robust_increment_std, 750),
    GAVI250nr = GAVI$new(
      "GAVI250-non-robust",
      SIGNIFICANCE_LEVEL,
      non_robust_increment_std,
      250
    ),
    GAVI500nr = GAVI$new(
      "GAVI500-non-robust",
      SIGNIFICANCE_LEVEL,
      non_robust_increment_std,
      500
    ),
    GAVI750nr = GAVI$new(
      "GAVI750-non-robust",
      SIGNIFICANCE_LEVEL,
      non_robust_increment_std,
      750
    ),
    # -- CAA (Statsig)
    CAA = CAA$new("CAA", SIGNIFICANCE_LEVEL, robust_increment_std),
    CAAnr = CAA$new("CAA-non-robust", SIGNIFICANCE_LEVEL, non_robust_increment_std),
    # -- Classical (a z-test conducted once at the end of the experiment)
    Classical = Bonferroni$new("Classical", SIGNIFICANCE_LEVEL, robust_increment_std, 1),
    Classicalnr = Bonferroni$new(
      "Classical-non-robust",
      SIGNIFICANCE_LEVEL,
      non_robust_increment_std,
      1
    )
  ))
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


estimate_expected_number_of_orders = function(data) {
  return(setorder(data[, date := as.Date(as.POSIXct(occurred_at, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"))][, .(num_obs = .N), by = date], "date")[["num_obs"]])
}


# Conducting the Evaluation

input_files = sort(list.files(DATA_DIRECTORY))

stopifnot(length(input_files) >= 2)

USER_ID_COL = "user_id"
METRIC_COL = "gmv_euro"
DTTM_COL = "occurred_at"

system(sprintf("mkdir -p %s", OUTPUT_DIRECTORY))

num_cores = detectCores() - 1
cl = makeCluster(num_cores)

print(num_cores)

clusterExport(
  cl,
  varlist=c(
    "DATA_DIRECTORY",
    "METRIC_COL",
    "USER_ID_COL",
    "DTTM_COL",
    "input_files",
    "DataCleaner",
    "estimate_increment_std",
    "estimate_expected_number_of_orders",
    "generate_assignments",
    "DataGeneratorFromRealEvents",
    "Aggregator",
    "initialise_continuous_methods",
    "SequentialTest",
    "Bonferroni",
    "SST",
    "ell",
    "compute_convolution",
    "estimate_false_detection_rate_bound",
    "compute_pSST_thresholds",
    "pSST",
    "GAVI",
    "mSPRT",
    "CAA",
    "NUM_ASSIGNMENT_REPLICATIONS",
    "SIGNIFICANCE_LEVEL",
    "OUTPUT_DIRECTORY"
  )
)


process_file = function(i) {
  library(arrow)
  library(data.table)
  library(sandwich)
  library(mvtnorm)
  library(stats)
  
  set.seed(2024 + i - 1)
  
  raw_preceeding_data = read_parquet(sprintf("%s/%s", DATA_DIRECTORY, input_files[i - 1]))
  data_cleaner = DataCleaner$new(raw_preceeding_data, METRIC_COL, USER_ID_COL, DTTM_COL, q =
                                   0.999)
  # print(sprintf(
  #   "Total order value cutoff = %.2f",
  #   data_cleaner$event_value_cutoff
  # ))
  preceeding_data = data_cleaner$clean(raw_preceeding_data)
  # num_unique_users_before = length(unique(raw_preceeding_data[[USER_ID_COL]]))
  # num_unique_users_after = length(unique(preceeding_data[[USER_ID_COL]]))
  
  robust_increment_std = estimate_increment_std(preceeding_data, METRIC_COL, USER_ID_COL)
  # print(
  #   sprintf(
  #     "Robust increment std estimate (using preceeding data) = %.2f",
  #     robust_increment_std
  #   )
  # )
  non_robust_increment_std = sqrt(mean(preceeding_data[[METRIC_COL]] ^ 2))
  # print(
  #   sprintf(
  #     "Non-robust increment std estimate (using preceeding data) = %.2f",
  #     non_robust_increment_std
  #   )
  # )
  expected_num_observations = estimate_expected_number_of_orders(preceeding_data)
  
  data = data_cleaner$clean(read_parquet(sprintf("%s/%s", DATA_DIRECTORY, input_files[i])))
  actual_num_observations = estimate_expected_number_of_orders(data)
  
  data_generator = DataGeneratorFromRealEvents$new(data, METRIC_COL, USER_ID_COL)
  aggregator = Aggregator$new()
  
  continuous_methods = initialise_continuous_methods(
    robust_increment_std,
    non_robust_increment_std,
    expected_num_observations,
    actual_num_observations
  )
  
  for (r in 1:NUM_ASSIGNMENT_REPLICATIONS) {
    #   if (r %% 100 == 0) {
    #     print(sprintf("Replication # %.04d", r))
    #   }
    trajectory = data_generator$generate_cumulative_difference_trajectory()
    # the trajectory of the cumulative difference in the revenue between
    # control and treatment
    
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
  output_file_name = paste(strsplit(input_files[i], ".parquet"), ".csv", sep =
                             "")
  write.csv(result,
            paste(OUTPUT_DIRECTORY, output_file_name, sep = "/"),
            row.names = FALSE)
  
  return(result)
}

results = parLapply(cl, 2:length(input_files), process_file)

# Stop the cluster
stopCluster(cl)

# Aggregating the results

# List all files in the directory
files = list.files(OUTPUT_DIRECTORY, full.names = TRUE)

# Read and combine all CSV files into a single data.table
df = rbindlist(lapply(files, function(file) fread(file) %>% .[, file := basename(file)]))

# Group by 'method' and calculate the mean detection rate
fdr_dt = df[, .(detection_rate = mean(detection_rate, na.rm = TRUE)), by = method]

# Add a 'variance_estimate' column based on the method name
fdr_dt[, `:=`(
  variance_estimate = ifelse(grepl("non-robust$", method), "non-robust", "robust"),
  method = str_split_fixed(method, "-", 2)[, 1]
)]

result_dt = dcast(fdr_dt, method ~ variance_estimate, value.var = "detection_rate")

# Round the values to 2 decimal places
result_dt[, `:=`(
  robust = round(robust, 2),
  `non-robust` = round(`non-robust`, 2)
)]

methods = c(
  "SST", "pSST14", "GAVI250", "GAVI500", "GAVI750",
  "mSPRT100", "mSPRT011", "mSPRT025",
  "GAVI250", "GAVI500", "GAVI750",
  "CAA", "Classical"
)

print(result_dt[methods, c("method", "robust", "non-robust")])
