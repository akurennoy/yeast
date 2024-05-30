# library(tidyverse)
# library(xtable)
# library(scales)

source("methods/bonferroni.R")
source("methods/caa.R")
source("methods/gavi.R")
source("methods/gst.R")
source("methods/msprt.R")
source("methods/psst.R")
source("methods/sst.R")
source("utils.R")


# References
# [1] https://engineering.atspotify.com/2023/03/choosing-sequential-testing-framework-comparisons-and-discussions/

# Global Simulation Settings

NUM_OBSERVATIONS = 500 # number of observations per experiment group,
                       # observations come in pairs
NUM_REPLICATIONS = 10
EVENT_VALUE_STD = 1
SIGNIFICANCE_LEVEL = 0.05
EFFECT_SIZES = c(0.0, 0.1, 0.2, 0.3, 0.4)

# Sequential Test Setup

increment_std = sqrt(2) * EVENT_VALUE_STD
continuous_methods = list(
  # -- SST
  SST=SST$new("SST", 0.05, NUM_OBSERVATIONS, increment_std),
  # -- pSST
  pSST07=pSST$new(
    "pSST07", 0.05, round((1: 7) * (NUM_OBSERVATIONS /  7)), increment_std
  ),
  pSST14=pSST$new(
    "pSST14", 0.05, round((1:14) * (NUM_OBSERVATIONS / 14)), increment_std
  ),
  # -- mSPRT
  mSPRTphi100=mSPRT$new("mSPRT100", 0.05, increment_std, 100),
  mSPRTphi025=mSPRT$new("mSPRT025", 0.05, increment_std,  25),
  mSPRTphi011=mSPRT$new("mSPRT011", 0.05, increment_std,  1 / 0.3^2),
  # -- GAVI
  GAVI250=GAVI$new("GAVI250", 0.05, increment_std, 750),
  GAVI500=GAVI$new("GAVI500", 0.05, increment_std, 500),
  GAVI750=GAVI$new("GAVI750", 0.05, increment_std, 750),
  Bonferroni_stream=Bonferroni$new(
    "Bonferroni", 0.05, increment_std, NUM_OBSERVATIONS
  ),
  # -- CAA (Statsig)
  CAA=CAA$new("CAA", 0.05, increment_std)
)
discrete_methods = list( # keys are the numbers of discrete checks
  "14"=list(
    # -- GST with accurate sampling
    GST14=GST$new("GST", 0.05, 1, 14, 14, increment_std),
    GST14phi3=GST$new("GSTphi3", 0.05, 3, 14, 14, increment_std),
    # -- GST with undersampling (e.g. fewer checks than expected)
    GST14u=GST$new("GSTu", 0.05, 1,  28, 14, increment_std),
    GST14uphi3=GST$new("GSTuphi3", 0.05, 3,  28, 14, increment_std),
    # -- GST with oversampling (e.g. more checks than expected)
    GST14o=GST$new("GSTo", 0.05, 1,  7, 14, increment_std),
    GST14ophi3=GST$new("GSTophi3", 0.05, 3,  7, 14, increment_std),
    # -- Bonferroni
    Bonferroni14=Bonferroni$new("Bonferroni", 0.05, increment_std, 14)
  ),
  "28"=list(
    # -- GST with accurate sampling
    GST28=GST$new("GST", 0.05, 1, 28, 28, increment_std),
    GST28phi3=GST$new("GSTphi3", 0.05, 3, 28, 28, increment_std),
    # -- GST with undersampling (e.g. fewer checks than expected)
    GST28u=GST$new("GSTu", 0.05, 1,  56, 28, increment_std),
    GST28uphi3=GST$new("GSTuphi3", 0.05, 3,  56, 28, increment_std),
    # -- GST with oversampling (e.g. more checks than expected)
    GST28o=GST$new("GSTo", 0.05, 1, 14, 28, increment_std),
    GST28ophi3=GST$new("GSTophi3", 0.05, 3, 14, 28, increment_std),
    # -- Bonferroni
    Bonferroni28=Bonferroni$new("Bonferroni", 0.05, increment_std, 28)
  ),
  "42"=list(
    # -- GST with accurate sampling
    GST42=GST$new("GST", 0.05, 1, 42, 42, increment_std),
    GST42phi3=GST$new("GSTphi3", 0.05, 3, 42, 42, increment_std),
    # -- GST with undersampling (e.g. fewer checks than expected)
    GST42u=GST$new("GSTu", 0.05, 1,  84, 42, increment_std),
    GST42uphi3=GST$new("GSTuphi3", 0.05, 3,  84, 42, increment_std),
    # -- GST with oversampling (e.g. more checks than expected)
    GST42o=GST$new("GSTo", 0.05, 1, 21, 42, increment_std),
    GST42ophi3=GST$new("GSTophi3", 0.05, 3, 21, 42, increment_std),
    # -- Bonferroni
    Bonferroni42=Bonferroni$new("Bonferroni", 0.05, increment_std, 42)
  ),
  "56"=list(
    # -- GST with accurate sampling
    GST56=GST$new("GST", 0.05, 1, 56, 56, increment_std),
    GST56phi3=GST$new("GSTphi3", 0.05, 3, 56, 56, increment_std),
    # -- GST with undersampling (e.g. fewer checks than expected)
    GST56u=GST$new("GSTu", 0.05, 1, 112, 56, increment_std),
    GST56uphi3=GST$new("GSTuphi3", 0.05, 3, 112, 56, increment_std),
    # -- GST with oversampling (e.g. more checks than expected)
    GST56o=GST$new("GSTo", 0.05, 1, 28, 56, increment_std),
    GST56ophi3=GST$new("GSTophi3", 0.05, 3, 28, 56, increment_std),
    # -- Bonferroni
    Bonferroni56=Bonferroni$new("Bonferroni", 0.05, increment_std, 56)
  )
)

# Simulations

run_experiment = function(
    seed,
    event_value_generator,
    continuous_methods,
    discrete_methods
) {
  set.seed(seed)
  data_generator = DataGenerator$new(event_value_generator)
  
  aggregator = Aggregator$new()
  
  for(r in 1:NUM_REPLICATIONS) {
    for(relative_effect in EFFECT_SIZES) {
      trajectory = data_generator$generate_cumulative_difference_trajectory(
        NUM_OBSERVATIONS, relative_effect
      ) # the trajectory of the cumulative difference in the metric of interest
        # between control and treatment
      
      # -- continuous monitoring methods
      
      for(statistical_test in continuous_methods) {
        detection_indicators = statistical_test$monitor(trajectory)
        # -- continuous monitoring mode
        aggregator$update(
          relative_effect,
          "stream",
          statistical_test$name,
          any(detection_indicators),
          get_savings(
            detection_indicators, 1:NUM_OBSERVATIONS, NUM_OBSERVATIONS
          )
        )
        # -- discrete monitoring mode
        for(num_checks_str in names(discrete_methods)) {
          num_checks = as.numeric(num_checks_str)
          check_times = NUM_OBSERVATIONS * seq(1 / num_checks, 1, 1 / num_checks)
          aggregator$update(
            relative_effect,
            num_checks_str,
            statistical_test$name,
            any(detection_indicators[check_times]),
            get_savings(
              detection_indicators[check_times], check_times, NUM_OBSERVATIONS
            )
          )
        }
      }
      
      # -- discrete monitoring methods
      
      for(num_checks_str in names(discrete_methods)) {
        num_checks = as.numeric(num_checks_str)
        for(statistical_test in discrete_methods[[num_checks_str]]) {
          detection_indicators = statistical_test$monitor(trajectory)
          check_times = NUM_OBSERVATIONS * seq(1 / num_checks, 1, 1 / num_checks)
          aggregator$update(
            relative_effect,
            num_checks_str,
            statistical_test$name,
            any(detection_indicators),
            get_savings(
              detection_indicators, check_times, NUM_OBSERVATIONS
            )
          )
        }
      }
    }
  }
  return(aggregator$get_result())
}


print("MAIN EXPERIMENT")
result = run_experiment(
  8163, # This is the seed value used in ref. [1].
        # We will use it for the normal event value distribution
        # to reproduce the results of ref. [1].
  NormalEventValueGenerator$new(1, EVENT_VALUE_STD),
  continuous_methods,
  discrete_methods
)
write.csv(result, "normal.csv", row.names=FALSE)
report(result, continuous_methods)

print("EXPERIMENTS WITH NON-NORMAL DATA")

# We are setting the seed differently depending on the chosen event value
# distribution so that the simulation results obtained for the different
# distribution types are not correlated.

# The parameters of the non-normal distributions were set to match the
# coefficient of variation of the normal distribution used in ref [1].

run_experiment(
  2023,
  ShiftedStudentEventValueGenerator$new(sqrt(3), 3),
  continuous_methods,
  discrete_methods
)
write.csv(result, "student.csv", row.names=FALSE)
report(result, continuous_methods)

run_experiment(
  2024,
  GammaEventValueGenerator$new(1, 2),
  continuous_methods,
  discrete_methods
)
write.csv(result, "gamma.csv", row.names=FALSE)
report(result, continuous_methods)
