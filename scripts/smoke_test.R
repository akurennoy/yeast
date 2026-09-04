# A short null-only simulation, to be run before the full experiment.
#
# Run from the repository root:  Rscript scripts/smoke_test.R
#
# With 2,000 replications the Monte Carlo standard error on a detection rate
# near 0.05 is about 0.005, which is enough to catch a boundary that is wrong by
# an order of magnitude but not enough to certify calibration. The assertions
# below are therefore deliberately loose: they are a sanity gate, not a result.

source("methods/bonferroni.R")
source("methods/gavi.R")
source("methods/ld_obf.R")
source("methods/msprt.R")
source("methods/pyeast.R")
source("methods/yeast.R")
source("methods/sec_c_2st_ons_qda.R")
source("data_generation.R")
source("utils.R")

NUM_OBSERVATIONS = 500
NUM_REPLICATIONS = as.integer(Sys.getenv("YEAST_SMOKE_REPLICATIONS", "2000"))
SIGNIFICANCE_LEVEL = 0.05
INCREMENT_STD = sqrt(2)

set.seed(31337)

methods = list(
  YEAST = YEAST$new("YEAST", SIGNIFICANCE_LEVEL, NUM_OBSERVATIONS, INCREMENT_STD),
  pYEAST14 = pYEAST$new("pYEAST14", SIGNIFICANCE_LEVEL,
                        round((1:14) * (NUM_OBSERVATIONS / 14)), INCREMENT_STD),
  mSPRTphi100 = mSPRT$new("mSPRTphi100", SIGNIFICANCE_LEVEL, INCREMENT_STD, 100),
  GAVI250 = GAVI$new("GAVI250", SIGNIFICANCE_LEVEL, INCREMENT_STD, 250),
  Bonferroni = Bonferroni$new("Bonferroni", SIGNIFICANCE_LEVEL, INCREMENT_STD,
                              NUM_OBSERVATIONS),
  LanDeMetsOBF = LanDeMetsOBF$new("LanDeMetsOBF", SIGNIFICANCE_LEVEL, INCREMENT_STD)
)

data_generator = DataGenerator$new(NormalEventValueGenerator$new(1, 1))
seq_c_2st = SeqC2ST$new("SeqC2ST_QDA", SIGNIFICANCE_LEVEL)

detections = setNames(integer(length(methods) + 1),
                      c(names(methods), "SeqC2ST_QDA"))
for (r in 1:NUM_REPLICATIONS) {
  generation_result = data_generator$generate_cumulative_difference_trajectory(
    NUM_OBSERVATIONS, 0.0
  )
  for (method_name in names(methods)) {
    if (any(methods[[method_name]]$monitor(generation_result$trajectory))) {
      detections[[method_name]] = detections[[method_name]] + 1
    }
  }
  event_values = generation_result$event_values
  crossed = seq_c_2st$monitor(
    trajectory = cumsum(as.vector(rbind(event_values$treatment, -event_values$control))),
    assignment_indicators = rep(c(1, 0), length.out = 2 * NUM_OBSERVATIONS)
  )
  if (any(crossed)) {
    detections[["SeqC2ST_QDA"]] = detections[["SeqC2ST_QDA"]] + 1
  }
}

rates = detections / NUM_REPLICATIONS
ci = wilson_interval(detections, NUM_REPLICATIONS)
for (i in seq_along(rates)) {
  cat(sprintf("%-16s %.4f [%.4f, %.4f]\n", names(rates)[i], rates[[i]],
              ci$lower[[i]], ci$upper[[i]]))
}

# Loose gate: a one-sided nominal 5% test must not be crossing wildly often, and
# must not be inert.
stopifnot(all(rates < 0.20))
stopifnot(rates[["YEAST"]] > 0.005)

cat("\nsmoke_test.R: null crossing rates are in the expected range\n")
