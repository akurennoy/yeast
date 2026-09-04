# Regenerates the horizon-mismatch diagnostic table from horizon_diagnostics.csv.
#
#   Rscript scripts/make_horizon_table.R
#
# Run scripts/horizon_diagnostics.R first.

if (!require("data.table", character.only = TRUE)) {
  install.packages("data.table")
  library("data.table", character.only = TRUE)
}
source("utils.R")

INPUT_FILE = "horizon_diagnostics.csv"
OUTPUT_DIR = "tables"
SIGNIFICANCE_LEVEL = 0.05
EFFECTS = c(0.00, 0.05, 0.10, 0.20, 0.50)

ARM_ORDER = c("primary", "truncated", "oracle")
ARM_LABEL = c(
  primary = "Primary: $\\tilde N$ boundary, all $M_T$ events",
  truncated = "Theorem-aligned: stop at $\\min(M_T,\\,\\tilde N)$",
  oracle = "Oracle: $M_T$ boundary, all $M_T$ events")

result = fread(INPUT_FILE)
n_tilde = max(result$n_tilde)
m_t = max(result$m_t)

num_cells = nrow(result)
confidence_level = 1 - SIGNIFICANCE_LEVEL / num_cells
ci = wilson_interval(result$num_detections, result$num_trials, confidence_level)
result[, `:=`(lower = ci$lower, upper = ci$upper)]

cell = function(row) {
  sprintf("%.3f [%.3f, %.3f]", row$detection_rate, row$lower, row$upper)
}

# The table is laid out with the three configurations as columns: with the
# effects as columns instead, five "rate [lo, hi]" cells overflow the text
# width by about 150pt.
lines = character()
for (this_effect in EFFECTS) {
  cells = vapply(ARM_ORDER,
                 function(arm) cell(result[mode == arm & effect == this_effect]),
                 character(1))
  label = if (this_effect == 0) {
    "$\\xi = 0.00$ (size)"
  } else {
    sprintf("$\\xi = %.2f$", this_effect)
  }
  lines = c(lines, sprintf("  %s & %s \\\\", label,
                           paste(cells, collapse = " & ")))
}

caption = sprintf(
  paste("Horizon-Mismatch Diagnostics (Online Retail, YEAST only).",
        "The planned event budget is $\\tilde N = %s$ and the realised count is",
        "$M_T = %s$, a ratio of $%.3f$. The first row is the empirical Type-I",
        "error; the remaining rows are power at the stated relative decrease in",
        "the order value. Intervals are Bonferroni-corrected Wilson intervals",
        "over the %d cells of the table. All three configurations run on",
        "identical replications."),
  format(n_tilde, big.mark = ","), format(m_t, big.mark = ","),
  m_t / n_tilde, num_cells)

body = c(
  "\\begin{table}[ht]",
  sprintf("    \\caption{%s}\\label{tbl:horizon-diagnostics}", caption),
  "    \\centering",
  "    \\begin{footnotesize}",
  "    \\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{lccc}",
  "  \\hline",
  paste0(" effect & primary & theorem-aligned & oracle \\\\"),
  paste0("  & $\\tilde N$ bd., all $M_T$ & stop at $\\min(M_T,\\tilde N)$",
         " & $M_T$ bd., all $M_T$ \\\\"),
  "  \\hline",
  lines,
  "   \\hline",
  "\\end{tabular}",
  "    \\end{footnotesize}",
  "\\end{table}")

dir.create(OUTPUT_DIR, showWarnings = FALSE)
writeLines(body, file.path(OUTPUT_DIR, "horizon_diagnostics.tex"))
cat(sprintf("wrote %s\n", file.path(OUTPUT_DIR, "horizon_diagnostics.tex")))
