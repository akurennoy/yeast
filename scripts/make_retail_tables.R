# Regenerates the three Online Retail LaTeX tables from the per-worker CSVs in
# online-retail-results/, so that the tables can be rebuilt in seconds without
# repeating the two-hour evaluation.
#
#   Rscript scripts/make_retail_tables.R

library(data.table)

source("utils.R")

INPUT_DIR = "online-retail-results"
OUTPUT_DIR = "tables"
SIGNIFICANCE_LEVEL = 0.05
CALIBRATION_TOLERANCE = 0.01

METHOD_ORDER = c(
  "YEAST", "mSPRTphi100", "mSPRTphi11", "mSPRTphi25",
  "GAVI250", "GAVI500", "GAVI750", "GAVI10K",
  "LanDeMetsOBF", "SeqC2ST_QDA"
)
POWER_EFFECTS = c(0.05, 0.1, 0.2, 0.5)


escape_latex = function(x) gsub("_", "\\\\_", x)

# The manuscript writes this benchmark with a hyphen; the code uses the
# underscore form as an R identifier.
DISPLAY_NAMES = c(SeqC2ST_QDA = "SeqC2ST-QDA")

display_name = function(x) {
  mapped = unname(DISPLAY_NAMES[x])
  ifelse(is.na(mapped), x, mapped)
}



write_table = function(lines, filename, caption, label, column_spec, header,
                       size = "small", colsep = "4pt") {
  dir.create(OUTPUT_DIR, showWarnings = FALSE)
  body = c(
    "\\begin{table}[ht]",
    sprintf("    \\caption{%s}\\label{%s}", caption, label),
    "    \\centering",
    sprintf("    \\begin{%s}", size),
    sprintf("    \\setlength{\\tabcolsep}{%s}", colsep),
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


files = list.files(INPUT_DIR, full.names = TRUE, pattern = "[.]csv$")
stopifnot(length(files) > 0)
df = rbindlist(lapply(files, fread))

totals = df[, .(num_detections = sum(num_detections),
                num_trials = sum(num_trials),
                total_savings = sum(total_savings)), by = .(method, effect)]
totals[, `:=`(
  variance_estimate = ifelse(grepl("non-robust", method), "non-robust", "robust"),
  base_method = sub("-.*$", "", method)
)]
totals[, `:=`(detection_rate = num_detections / num_trials,
              average_savings = total_savings / num_trials)]

# One Bonferroni correction across every reported detection-rate cell.
num_cells = nrow(totals)
ci = wilson_interval(totals$num_detections, totals$num_trials,
                     confidence_level = 1 - 0.05 / num_cells)
totals[, `:=`(lower = ci$lower, upper = ci$upper)]

# Three decimals normally, but a value strictly inside (0, 1) must never be
# printed as 0.000 or 1.000: that is exactly the degenerate presentation the
# switch from Wald to Wilson intervals was made to avoid.
fmt = function(...) {
  x = c(...)
  collapses = any(x > 0 & x < 1 & (round(x, 3) == 0 | round(x, 3) == 1))
  sprintf(if (collapses) "%.4f" else "%.3f", x)
}

# compact drops the space after the comma; the power table needs the ~12pt this
# saves to fit inside \textwidth.
cell = function(rate, lower, upper, bold = FALSE, compact = FALSE) {
  parts = fmt(rate, lower, upper)
  text = sprintf(if (compact) "%s [%s,%s]" else "%s [%s, %s]",
                 parts[1], parts[2], parts[3])
  if (isTRUE(bold)) sprintf("\\textbf{%s}", text) else text
}

robust = totals[variance_estimate == "robust"]
null_rate = robust[effect == 0.0, .(base_method, inflated = is_size_inflated(
  num_detections, num_trials, SIGNIFICANCE_LEVEL, 1 - 0.05 / num_cells))]


# -- False detection rate: robust against non-robust variance estimate

fdr_lines = character()
for (i in seq_along(METHOD_ORDER)) {
  this_method = METHOD_ORDER[i]
  r = totals[base_method == this_method & effect == 0.0 & variance_estimate == "robust"]
  nr = totals[base_method == this_method & effect == 0.0 & variance_estimate == "non-robust"]
  # SeqC2ST-QDA takes no variance input, so its single result spans both columns.
  if (nrow(nr) == 0) {
    cells = sprintf("\\multicolumn{2}{c}{%s}",
                    cell(r$detection_rate, r$lower, r$upper,
                         abs(r$detection_rate - SIGNIFICANCE_LEVEL) <= CALIBRATION_TOLERANCE))
  } else {
    cells = paste(
      cell(r$detection_rate, r$lower, r$upper,
           abs(r$detection_rate - SIGNIFICANCE_LEVEL) <= CALIBRATION_TOLERANCE),
      cell(nr$detection_rate, nr$lower, nr$upper),
      sep = " & "
    )
  }
  fdr_lines = c(fdr_lines, sprintf("  %d & %s & %s \\\\", i,
                                   escape_latex(display_name(this_method)), cells))
}
write_table(fdr_lines, "retail_fdr.tex", paste0("Semi-synthetic Experiment: Empirical Type-I Error (Online Retail).",
            " Several entries are run outside the setting of their published guarantee; see Table~\\ref{tbl:applicability} for what was run and how each row should be read."),
            "tbl:real-world-experiment-online-retail", "llcc",
            " & method & robust var. est. & non-robust var. est.\\\\")


# -- Power, robust variance estimate only

power = robust[effect %in% POWER_EFFECTS]
power = merge(power, null_rate, by = "base_method", sort = FALSE)
power[, eligible := !inflated]
power[, is_bold := FALSE]
for (this_effect in POWER_EFFECTS) {
  block = power[effect == this_effect & base_method %in% METHOD_ORDER]
  if (any(block$eligible)) {
    best = which(block$eligible)[which.max(block[eligible == TRUE]$detection_rate)]
    flag = block$eligible &
      block$lower <= block$upper[best] & block$upper >= block$lower[best]
    power[effect == this_effect & base_method %in% METHOD_ORDER,
          is_bold := flag]
  }
}

# The non-sequential t-test is a descriptive reference, so it is set in italics
# and left out of the bold ranking.
t_test = robust[base_method == "Classical" & effect %in% POWER_EFFECTS]
setorder(t_test, effect)
power_lines = sprintf(
  " & \\textit{Non-seq. $t$-test} & %s \\\\",
  paste(vapply(seq_len(nrow(t_test)), function(k) sprintf(
    "\\textit{%s}", cell(t_test$detection_rate[k], t_test$lower[k], t_test$upper[k],
                          compact = TRUE)
  ), character(1)), collapse = " & ")
)
for (i in seq_along(METHOD_ORDER)) {
  this_method = METHOD_ORDER[i]
  cells = vapply(POWER_EFFECTS, function(e) {
    row = power[base_method == this_method & effect == e]
    cell(row$detection_rate, row$lower, row$upper, row$is_bold, compact = TRUE)
  }, character(1))
  # A row whose realised size exceeds the nominal level is set in italics: its
  # detection rates are not power figures comparable with the rest.
  if (!power[base_method == this_method][1]$eligible) {
    cells = sprintf("\\textit{%s}", cells)
  }
  power_lines = c(power_lines, sprintf("  %d & %s & %s \\\\", i,
                                       escape_latex(display_name(this_method)),
                                       paste(cells, collapse = " & ")))
}
write_table(
  power_lines, "retail_power.tex",
  paste0("Semi-synthetic Experiment: Power (Online Retail).", " Several entries are run outside the setting of their published guarantee; see Table~\\ref{tbl:applicability} for what was run and how each row should be read."),
  "tbl:semi-synthetic-power", "rlllll",
  paste(
    " & method & \\multicolumn{4}{c}{relative decrease in the order value}\\\\",
    sprintf(" &  & %s\\\\", paste(sprintf("%.2f", POWER_EFFECTS), collapse = " & ")),
    sep = "\n"
  ),
  size = "footnotesize", colsep = "2pt"
)


# -- Average sample savings, robust variance estimate only

savings_lines = character()
for (i in seq_along(METHOD_ORDER)) {
  this_method = METHOD_ORDER[i]
  cells = vapply(POWER_EFFECTS, function(e) {
    row = robust[base_method == this_method & effect == e]
    sprintf("%.0f", 100 * row$average_savings)
  }, character(1))
  savings_lines = c(savings_lines, sprintf("  %d & %s & %s \\\\", i,
                                           escape_latex(display_name(this_method)),
                                           paste(cells, collapse = " & ")))
}
write_table(
  savings_lines, "retail_savings.tex",
  "Semi-synthetic Experiment: Sample Savings, \\% (Online Retail)",
  "tbl:retail-savings", "rlrrrr",
  paste(
    " & method & \\multicolumn{4}{c}{relative decrease in the order value}\\\\",
    sprintf(" &  & %s\\\\", paste(sprintf("%.2f", POWER_EFFECTS), collapse = " & ")),
    sep = "\n"
  )
)
