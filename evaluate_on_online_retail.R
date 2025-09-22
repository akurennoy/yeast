# In this script, we evaluate the sequential tests on a standard, publicly
# available dataset. It contains transactions from a UK-based and
# registered non-store online retail. The dataset is called Online Retail and
# can be downloaded from https://archive.ics.uci.edu/dataset/352/online+retail.
# We randomly assign users that placed those orders to either control or
# treatment and then apply the sequential tests to the resulting trajectory of
# the cumulative revenue difference between the two groups. We will then repeat
# this multiple times and measure the number of detections. 
# Since the users are assigned randomly post-factum and no real
# treatment gets applied to the treatment group, the relative frequency of
# detection should be close to the nominal significance level.


library(arrow)
library(dplyr)
library(parallel)
library(readxl)
library(sandwich)
library(stringr)
library(xtable)

source("data_generation.R")
source("methods/bonferroni.R")
source("methods/caa.R")
source("methods/gavi.R")
source("methods/gst.R")
source("methods/msprt.R")
source("methods/yeast.R")
source("methods/pyeast.R")
source("methods/ld_obf.R")
source("methods/sec_c_2st_ons_qda.R")

source("utils.R")


# Global Settings


SIGNIFICANCE_LEVEL = 0.05
NUM_ASSIGNMENT_REPLICATIONS = 100000
PROCESSED_DATA_DIRECTORY = "./online-retail-orders/"
OUTPUT_DIRECTORY = "./online-retail-results"


# We being by downloading the data to the current directory, unpacking it, and
# loading into a data frame.


url = "https://archive.ics.uci.edu/static/public/352/online+retail.zip"
zip_file = "online_retail.zip"
extracted_file = "Online Retail.xlsx"

download.file(url, destfile = zip_file, mode = "wb")
unzip(zip_file, files = extracted_file)

data <- read_excel(extracted_file)


# Aggregating the records by order


orders = data %>%
  filter(
    !str_detect(InvoiceNo, "^c") 
    & !is.na(CustomerID)
  ) %>%
  group_by(InvoiceNo, CustomerID) %>%
  summarise(
    order_value = sum(Quantity * UnitPrice),
    occurred_at = min(InvoiceDate)
  ) %>%
  arrange(occurred_at) %>%
  # mutate(date = as.Date(occurred_at)) %>%
  mutate(year_month = format(occurred_at, "%Y-%m")) %>%
  filter(
    year_month != "2011-12"  # excluding the last month since the data for it is incomplete
  ) %>%
  rename(
    user_id = CustomerID
  )

# # ------------
# mu <- mean(orders$order_value, na.rm = TRUE)
# sigma <- sd(orders$order_value, na.rm = TRUE)
# 
# meanlog <- log(mu^2 / sqrt(sigma^2 + mu^2))
# sdlog <- sqrt(log(1 + (sigma^2 / mu^2)))
# 
# orders <- orders %>%
#   mutate(order_value = rlnorm(n(), meanlog, sdlog))
# # -------------

# We split the data into two parts (6 months each)
# The first part is used for estimating the input parameters.
# The second part is used for the validation.


# dir.create(PROCESSED_DATA_DIRECTORY)
write_parquet(
  orders %>% filter(year_month %in% c("2010-12", "2011-01", "2011-02", "2011-03", "2011-04", "2011-05")),
  paste0(
    PROCESSED_DATA_DIRECTORY,
    "00_201012_201105.parquet"
  )
)
write_parquet(
  orders %>% filter(year_month %in% c("2011-06", "2011-07", "2011-08", "2011-09", "2011-10", "2011-11")),
  paste0(
    PROCESSED_DATA_DIRECTORY,
    "01_201106_201111.parquet"
  )
)


# Definitions


initialise_continuous_methods = function(robust_increment_std,
                                         non_robust_increment_std,
                                         expected_num_observations,
                                         actual_num_observations) {
  return(list(
    # -- YEAST
    YEAST = YEAST$new(
      "YEAST",
      SIGNIFICANCE_LEVEL,
      sum(expected_num_observations),
      robust_increment_std
    ),
    YEASTnr = YEAST$new(
      "YEAST-non-robust",
      SIGNIFICANCE_LEVEL,
      sum(expected_num_observations),
      non_robust_increment_std
    ),
    # -- pYEAST
    pYEAST6 = pYEAST$new(
      "pYEAST6",
      SIGNIFICANCE_LEVEL,
      cumsum(expected_num_observations),
      robust_increment_std,
      cumsum(actual_num_observations)
    ),
    pYEAST6nr = pYEAST$new(
      "pYEAST6-non-robust",
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
    GAVI10K = GAVI$new("GAVI10K", SIGNIFICANCE_LEVEL, robust_increment_std, 10000),
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
    GAVI10Knr = GAVI$new(
      "GAVI10K-non-robust",
      SIGNIFICANCE_LEVEL,
      non_robust_increment_std,
      10000
    ),
    LanDeMetsOBF = LanDeMetsOBF$new(
      "LanDeMetsOBF", SIGNIFICANCE_LEVEL, robust_increment_std
    ),
    LanDeMetsOBFnr = LanDeMetsOBF$new(
      "LanDeMetsOBF-non-robust", SIGNIFICANCE_LEVEL, non_robust_increment_std
    ),
    SeqC2ST_QDA = SeqC2ST$new(
      "SeqC2ST_QDA", SIGNIFICANCE_LEVEL
    ),
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
                                  user_id_col) {
  dg = DataGeneratorFromRealEvents$new(data, metric_col, user_id_col, 1)
  dt = dg$real_events_data_table[, x := get(metric_col) * (1 - 2 * a1)]
  grouped = dt[, .(x = sum(x)), by = .(user_id, year_month)]
  # TODO: take the user id column name from the parameter instead of using
  #.      a hardcoded value in the by parameter
  model <- lm(x ~ 1, data = grouped)
  v = (
    vcovCL(model, cluster = grouped[[user_id_col]], type = "HC3")[[1]]
  )
  return(sqrt(v * nrow(grouped) ** 2 / nrow(dt)))
}

# estimate_increment_std = function(data,
#                                   metric_col,
#                                   user_id_col) {
#   dg = DataGeneratorFromRealEvents$new(data, metric_col, user_id_col, 1)
#   dt = dg$real_events_data_table[, x := get(metric_col) * (1 - 2 * a1)]
#   grouped = dt[, .(x = sum(x)), by = .(user_id)]
#   # TODO: take the user id column name from the parameter instead of using
#   #.      a hardcoded value in the by parameter

#   return(sqrt(var(grouped$x) * nrow(grouped) / nrow(dt)))
# }


estimate_expected_number_of_orders = function(data) {
  return(setorder(data[, .(num_obs = .N), by = year_month], "year_month")[["num_obs"]])
}


# Conducting the Evaluation

input_files = sort(list.files(PROCESSED_DATA_DIRECTORY))

USER_ID_COL = "user_id"
METRIC_COL = "order_value"
DTTM_COL = "occurred_at"

system(sprintf("mkdir -p %s", OUTPUT_DIRECTORY))

num_cores = detectCores() - 1
cl = makeCluster(num_cores)  # We parallelise the simulations over the cores

print(num_cores)

num_assignment_replications = rep(NUM_ASSIGNMENT_REPLICATIONS %/% num_cores, num_cores)
num_assignment_replications[num_cores] = num_assignment_replications[num_cores] + NUM_ASSIGNMENT_REPLICATIONS - sum(num_assignment_replications)
stopifnot(sum(num_assignment_replications) == NUM_ASSIGNMENT_REPLICATIONS)


clusterExport(
  cl,
  varlist=c(
    "PROCESSED_DATA_DIRECTORY",
    "METRIC_COL",
    "USER_ID_COL",
    "DTTM_COL",
    "input_files",
    "num_assignment_replications",
    "DataCleaner",
    "estimate_increment_std",
    "estimate_expected_number_of_orders",
    "generate_assignments",
    "DataFromRealEvents",
    "DataGeneratorFromRealEvents",
    "Aggregator",
    "initialise_continuous_methods",
    "SequentialTest",
    "Bonferroni",
    "YEAST",
    "ell",
    "compute_convolution",
    "estimate_false_detection_rate_bound",
    "compute_pYEAST_thresholds",
    "pYEAST",
    "GAVI",
    "mSPRT",
    "CAA",
    "LanDeMetsOBF",
    "QDAStats",
    "OnlineQDA",
    "SeqC2ST",
    "SIGNIFICANCE_LEVEL",
    "OUTPUT_DIRECTORY",
    "NUM_OBSERVATIONS",
    "get_savings"
  )
)


process_file = function(i) {
  library(arrow)
  library(data.table)
  library(sandwich)
  library(mvtnorm)
  library(stats)
  
  set.seed(i)
  
  file_index = 2
  raw_preceeding_data = read_parquet(sprintf("%s/%s", PROCESSED_DATA_DIRECTORY, input_files[file_index - 1]))
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
  
  data = data_cleaner$clean(read_parquet(sprintf("%s/%s", PROCESSED_DATA_DIRECTORY, input_files[file_index])))
  actual_num_observations = estimate_expected_number_of_orders(data)
  
  data_generator = DataGeneratorFromRealEvents$new(data, METRIC_COL, USER_ID_COL)
  aggregator = Aggregator$new()
  
  continuous_methods = initialise_continuous_methods(
    robust_increment_std,
    non_robust_increment_std,
    expected_num_observations,
    actual_num_observations
  )
  
  for (relative_effect in c(0.00, 0.05, 0.10, 0.20, 0.50)) { #-0.01, -0.02, -0.05, -0.10)) {
    for (r in 1:num_assignment_replications[i]) {
      #   if (r %% 100 == 0) {
      #     print(sprintf("Replication # %.04d", r))
      #   }
      generation_result = data_generator$generate_cumulative_difference_trajectory(-relative_effect)
      trajectory = generation_result$trajectory
      # the trajectory of the cumulative difference in the revenue between
      # control and treatment
      
      # computing the p-value of the standard non-sequential test
      # (applied at the end of the experiment)
      grouped_by_user = data.table(
        x=generation_result$signed_event_values, w=generation_result$assignment_indicators, user_id=generation_result$user_ids
      )[, .(x = sum(x), w =max(w)), by = .(user_id)]
      ttest_pvalue = t.test(
        -grouped_by_user[grouped_by_user$w == 1, x],
        grouped_by_user[grouped_by_user$w == 0, x],
        var.equal = TRUE
      )$p.value
      aggregator$update(relative_effect,
                        "non-sequential",
                        "ttest",
                        ttest_pvalue < SIGNIFICANCE_LEVEL,
                        0)
      
      for (statistical_test in continuous_methods) {
        detection_indicators = statistical_test$monitor(
          trajectory, generation_result$assignment_indicators
        )
        # -- continuous monitoring mode
        aggregator$update(relative_effect,
                          "stream",
                          statistical_test$name,
                          any(detection_indicators),
                          get_savings(
                            detection_indicators,
                            1:length(trajectory),
                            length(trajectory)
                          )
        )
      }
    }
  }
  
  result = aggregator$get_result()
  output_file_name = paste(strsplit(input_files[file_index], ".parquet"), sprintf("_%i", i), ".csv", sep = "")
  write.csv(result,
            paste(OUTPUT_DIRECTORY, output_file_name, sep = "/"),
            row.names = FALSE)
  
  return(result)
}

results = parLapply(cl, 1:num_cores, process_file)

stopCluster(cl)


# Aggregating the results


files = list.files(OUTPUT_DIRECTORY, full.names = TRUE)

df = rbindlist(lapply(files, function(file)
  fread(file) %>% .[, file := basename(file)]))

dr_dt = df[, .(num_detections = sum(num_detections),
                num_trials = sum(num_trials)), by = .(method, effect)]
num_methods = nrow(dr_dt)

dr_dt[, `:=`(detection_rate = num_detections / num_trials)]
dr_dt[, `:=`(
  variance_estimate = ifelse(grepl("non-robust", method), "non-robust", "robust"),
  method = str_split_fixed(method, "-", 2)[, 1],
  ci_pm = qnorm(1 - 0.05 / 2 / num_methods) * sqrt(detection_rate * (1 - detection_rate) / num_trials)  # this is the CI half-length
)]
dr_dt[, `:=`(
  detection_rate_with_ci = paste0(
    sprintf("%.4f", round(detection_rate, 4)),
    " ± ",
    sprintf("%.4f", round(ci_pm, 4))
  )
)]

# -- Computing the False Detection Rate

fdr_dt = dcast(dr_dt[dr_dt$effect == 0.0],
                  method ~ variance_estimate,
                  value.var = "detection_rate")

# -- Computing the Power

pow_dt = dcast(dr_dt[dr_dt$variance_estimate == 'robust'],
               method ~ effect,
               value.var = "detection_rate_with_ci")
# pow_dt = dcast(dr_dt[dr_dt$variance_estimate == 'robust'],
#                method ~ effect,
#                value.var = "detection_rate")
# ttest_values = pow_dt[method == "Classical",]
# pow_dt_numeric_cols = names(pow_dt)[sapply(pow_dt, is.numeric)]
# relative_pow_dt = copy(pow_dt)
# for (col in pow_dt_numeric_cols) {
#   relative_pow_dt[, (col) := round(get(col) / ttest_values[[col]], 4)]
# }


as_dt = df[, .(total_savings = sum(total_savings),
               num_detections = sum(num_detections)), by = .(method, effect)]
as_dt[, `:=`(average_savings = total_savings / num_detections)]
savings_dt = dcast(as_dt[as_dt$variance_estimate == 'robust'],
                   method ~ effect,
                   value.var = "average_savings")

# -- Printing Results

methods = c(
  "Classical",
  "YEAST",
  "mSPRT100",
  "mSPRT011",
  "mSPRT025",
  "GAVI250",
  "GAVI500",
  "GAVI750",
  "GAVI10K",
  "LanDeMetsOBF",
  "SeqC2ST_QDA"
)

print(fdr_dt[methods, c("method", "robust", "non-robust")])

print(xtable(fdr_dt[methods, c("method", "robust", "non-robust")]))

# pow_cols = names(pow_dt)
# pow_cols = c("method", sort(pow_dt_numeric_cols, decreasing = TRUE)[2:length(pow_dt_numeric_cols)])
pow_cols = c("method", names(pow_dt)[2:length(names(pow_dt))])

print(pow_dt[methods, ..pow_cols])

print(xtable(pow_dt[methods, ..pow_cols]))

# print(relative_pow_dt[methods, ..pow_cols])
# 
# print(xtable(relative_pow_dt[methods, ..pow_cols]))

savings_cols = names(savings_dt)
