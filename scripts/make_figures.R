# Regenerates power-4.png (the non-normal power curves) from student.csv and
# gamma.csv.
#
# Run from the repository root, after run_simulations.R:
#   Rscript scripts/make_figures.R

library(data.table)

OUTPUT_FILE = "power-4.png"
SIGNIFICANCE_LEVEL = 0.05

METHODS = c(
  "YEAST", "YEASTn110", "YEASTn120",
  "mSPRTphi100", "mSPRTphi11", "mSPRTphi25",
  "GAVI250", "GAVI500", "GAVI750",
  "Bonferroni", "LanDeMetsOBF", "SeqC2ST_QDA"
)

COLOURS = c("#1f77b4", "#17becf", "#9edae5", "#9467bd", "#c5b0d5", "#e377c2",
            "#2ca02c", "#98df8a", "#bcbd22", "#ff7f0e", "#8c564b", "#7f7f7f")
SYMBOLS = c(15, 16, 18, 17, 2, 4, 19, 6, 17, 17, 8, 16)


panel = function(result, title, show_ylab) {
  effects = sort(unique(result$effect))
  plot(NA, xlim = range(effects), ylim = c(0, 1.02),
       xlab = "true relative effect",
       ylab = if (show_ylab) "detection rate" else "",
       main = title, yaxt = "n")
  axis(2, at = c(0, SIGNIFICANCE_LEVEL, 0.5, 1.0),
       labels = c("0.00", "0.05", "0.50", "1.00"), las = 1)
  abline(h = SIGNIFICANCE_LEVEL, col = "red")
  for (i in seq_along(METHODS)) {
    curve = result[method == METHODS[i]]
    if (nrow(curve) == 0) next
    setorder(curve, effect)
    lines(curve$effect, curve$detection_rate, col = COLOURS[i], lty = i %% 5 + 1)
    points(curve$effect, curve$detection_rate, col = COLOURS[i], pch = SYMBOLS[i])
  }
}


student = fread("student.csv")[mode == "stream"]
gamma = fread("gamma.csv")[mode == "stream"]

png(OUTPUT_FILE, width = 1600, height = 640, res = 150)
layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.55))
par(mar = c(4, 4, 3, 1))
panel(student, "student", TRUE)
panel(gamma, "gamma", FALSE)
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = METHODS, col = COLOURS, pch = SYMBOLS,
       lty = seq_along(METHODS) %% 5 + 1, bty = "n", cex = 0.8)
dev.off()

cat(sprintf("wrote %s\n", OUTPUT_FILE))
