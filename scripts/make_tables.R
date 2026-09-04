# Regenerates every LaTeX table in the manuscript from the CSV output of
# run_simulations.R.
#
# Run from the repository root, after the experiments:
#   Rscript scripts/make_tables.R
#
# Tables are written to tables/*.tex and are meant to be \input into the
# manuscript rather than copied by hand, so that no reported figure can drift
# from the run that produced it.

library(data.table)

source("utils.R")

OUTPUT_DIR = "tables"
SIGNIFICANCE_LEVEL = 0.05

CONTINUOUS_METHOD_ORDER = c(
  "YEAST", "YEASTn110", "YEASTn120", "YEASTn80", "YEASTn90",
  "mSPRTphi100", "mSPRTphi11", "mSPRTphi25",
  "GAVI250", "GAVI500", "GAVI750",
  "Bonferroni", "LanDeMetsOBF", "SeqC2ST_QDA"
)

DISCRETE_METHOD_ORDER = c(
  "YEAST", "GST", "GSTphi3", "GSTo", "GSTophi3", "GSTu", "GSTuphi3"
)

DISCRETE_METHOD_LABELS = c(
  GSTo = "GSToversampled", GSTophi3 = "GSToversampledphi3",
  GSTu = "GSTundersampled", GSTuphi3 = "GSTundersampledphi3"
)


escape_latex = function(x) {
  return(gsub("_", "\\\\_", x, fixed = FALSE))
}

# The manuscript writes this benchmark with a hyphen; the code uses the
# underscore form as an R identifier.
DISPLAY_NAMES = c(SeqC2ST_QDA = "SeqC2ST-QDA")

display_name = function(x) {
  mapped = unname(DISPLAY_NAMES[x])
  ifelse(is.na(mapped), x, mapped)
}



# A cell is emboldened when it attains the highest rate in its column among the
# eligible methods, or when its interval overlaps the interval of the method
# that does. Eligibility is decided by the realised null rate, so that a method
# whose empirical Type-I error exceeds the nominal level is never presented as
# the most powerful.
bold_mask = function(rates, lower, upper, eligible) {
  mask = rep(FALSE, length(rates))
  if (!any(eligible)) {
    return(mask)
  }
  best = which(eligible)[which.max(rates[eligible])]
  mask[eligible] = lower[eligible] <= upper[best] & upper[eligible] >= lower[best]
  return(mask)
}


# Three decimals normally, but a value strictly inside (0, 1) must never print
# as 0.000 or 1.000: that degenerate presentation is what the switch from Wald
# to Wilson intervals was made to avoid.
fmt = function(...) {
  x = c(...)
  collapses = any(x > 0 & x < 1 & (round(x, 3) == 0 | round(x, 3) == 1))
  sprintf(if (collapses) "%.4f" else "%.3f", x)
}


rate_cell = function(rate, lower, upper, bold = FALSE) {
  parts = fmt(rate, lower, upper)
  text = sprintf("%s [%s,%s]", parts[1], parts[2], parts[3])
  if (isTRUE(bold)) sprintf("\\textbf{%s}", text) else text
}


write_table = function(lines, filename, caption, label, column_spec, header,
                       size = "small", colsep = NULL) {
  dir.create(OUTPUT_DIR, showWarnings = FALSE)
  body = c(
    "\\begin{table}[ht]",
    sprintf("\\caption{%s}\\label{%s}", caption, label),
    "\\centering",
    sprintf("\\begin{%s}", size),
    if (!is.null(colsep)) sprintf("\\setlength{\\tabcolsep}{%s}", colsep),
    sprintf("\\begin{tabular}{%s}", column_spec),
    "  \\hline",
    header,
    "  \\hline",
    lines,
    "   \\hline",
    "\\end{tabular}",
    sprintf("\\end{%s}", size),
    "\\end{table}"
  )
  writeLines(body, file.path(OUTPUT_DIR, filename))
  cat(sprintf("wrote %s\n", file.path(OUTPUT_DIR, filename)))
}


make_detection_rate_table = function(result, effects, filename, caption, label) {
  stream = result[mode == "stream" & method %in% CONTINUOUS_METHOD_ORDER
                  & effect %in% effects]
  stream[, method := factor(method, levels = CONTINUOUS_METHOD_ORDER)]
  setorder(stream, method, effect)

  num_cells = nrow(stream)
  confidence_level = 1 - 0.05 / num_cells
  ci = wilson_interval(stream$num_detections, stream$num_trials, confidence_level)
  stream[, `:=`(lower = ci$lower, upper = ci$upper)]

  null_rate = stream[effect == 0.0, .(method, inflated = is_size_inflated(
    num_detections, num_trials, SIGNIFICANCE_LEVEL, confidence_level))]
  stream = merge(stream, null_rate, by = "method", sort = FALSE)

  lines = character()
  methods = levels(droplevels(stream$method))
  for (this_effect in effects) {
    block = stream[effect == this_effect]
    setorder(block, method)
    # Only the alternatives are emboldened; the null column reports size.
    block[, is_bold := if (this_effect == 0.0) FALSE else bold_mask(
      detection_rate, lower, upper, !inflated
    )]
    stream[effect == this_effect, is_bold := block$is_bold]
  }

  for (i in seq_along(methods)) {
    this_method = methods[i]
    cells = character()
    for (this_effect in effects) {
      row = stream[method == this_method & effect == this_effect]
      cells = c(cells, rate_cell(row$detection_rate, row$lower, row$upper,
                                 row$is_bold))
    }
    lines = c(lines, sprintf("  %d & %s & %s \\\\", i,
                             escape_latex(display_name(this_method)), paste(cells, collapse = " & ")))
  }

  header = c(
    sprintf("  & effect size & %s\\\\", paste(sprintf("%.1f", effects), collapse = " & ")),
    sprintf("  & method %s \\\\", paste(rep("&", length(effects)), collapse = " "))
  )
  write_table(lines, filename, caption, label,
              paste0("rl", paste(rep("l", length(effects)), collapse = "")), header,
              size = "footnotesize", colsep = "2pt")
}


make_savings_table = function(result, effects, filename, caption, label) {
  savings_methods = setdiff(CONTINUOUS_METHOD_ORDER, grep("^YEASTn", CONTINUOUS_METHOD_ORDER, value = TRUE))
  stream = result[mode == "stream" & method %in% savings_methods
                  & effect %in% effects]
  stream[, method := factor(method, levels = savings_methods)]
  setorder(stream, method, effect)

  lines = character()
  methods = levels(droplevels(stream$method))
  for (i in seq_along(methods)) {
    cells = vapply(effects, function(e) sprintf(
      "%.0f", 100 * stream[method == methods[i] & effect == e]$average_savings
    ), character(1))
    lines = c(lines, sprintf("  %d & %s & %s \\\\", i,
                             escape_latex(display_name(methods[i])), paste(cells, collapse = " & ")))
  }
  header = c(
    sprintf("  & effect size & %s\\\\", paste(sprintf("%.1f", effects), collapse = " & ")),
    sprintf("  & method %s \\\\", paste(rep("&", length(effects)), collapse = " "))
  )
  write_table(lines, filename, caption, label,
              paste0("rl", paste(rep("r", length(effects)), collapse = "")), header)
}


make_discrete_table = function(result, this_effect, filename, caption, label,
                               with_stream = FALSE, drop_inflated = FALSE) {
  modes = c("14", "28", "42", "56", if (with_stream) "stream")
  block = result[effect == this_effect & mode %in% modes
                 & method %in% DISCRETE_METHOD_ORDER]
  eligible = DISCRETE_METHOD_ORDER
  if (drop_inflated) {
    # A method is excluded from the power table when its realised size is
    # inflated beyond Monte Carlo error at any of the discrete check counts.
    null_block = result[effect == 0.0 & mode %in% setdiff(modes, "stream")
                        & method %in% DISCRETE_METHOD_ORDER]
    confidence_level = 1 - 0.05 / nrow(null_block)
    inflated = unique(null_block[is_size_inflated(
      num_detections, num_trials, SIGNIFICANCE_LEVEL, confidence_level)]$method)
    eligible = setdiff(eligible, inflated)
  }
  block = block[method %in% eligible]
  lines = character()
  present = DISCRETE_METHOD_ORDER[DISCRETE_METHOD_ORDER %in% unique(block$method)]
  for (i in seq_along(present)) {
    this_method = present[i]
    cells = vapply(modes, function(m) {
      row = block[method == this_method & mode == m]
      if (nrow(row) == 0) "-" else sprintf("%.2f", row$detection_rate)
    }, character(1))
    label_text = DISCRETE_METHOD_LABELS[this_method]
    if (is.na(label_text)) label_text = this_method
    lines = c(lines, sprintf("  %d & %s & %s \\\\", i,
                             escape_latex(display_name(label_text)), paste(cells, collapse = " & ")))
  }
  header = sprintf("  & type & %s \\\\", paste(modes, collapse = " & "))
  write_table(lines, filename, caption, label,
              paste0("rl", strrep("r", length(modes))), header)
}


main = function() {
  normal = fread("normal.csv")
  make_detection_rate_table(
    normal, c(0.0, 0.1, 0.2, 0.3),
    "main_experiment.tex",
    "Simulation Experiment: \\\\ False Detection Rate and Power",
    "tbl:mainExperiment"
  )
  make_savings_table(
    normal, c(0.1, 0.2, 0.3, 0.4),
    "main_experiment_savings.tex",
    "Simulation Experiment: Sample/Time Savings, \\%",
    "tbl:mainExperimentSavings"
  )
  make_discrete_table(
    normal, 0.0, "appendix_fpr.tex",
    "False Positive Rate", "tbl:appendixFPR"
  )
  make_discrete_table(
    normal, 0.2, "appendix_power.tex",
    "Power (under a treatment effect of 0.2 standard deviations)",
    "tbl:appendixPower",
    with_stream = TRUE, drop_inflated = TRUE
  )
}


main()
