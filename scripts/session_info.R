# Prints the exact R and package versions used for a run, in a form suitable
# for committing as ENVIRONMENT.md.
#
#   Rscript scripts/session_info.R > ENVIRONMENT.md

REQUIRED_PACKAGES = c(
  "data.table", "dplyr", "ldbounds", "mvtnorm", "progress",
  "R6", "readxl", "sandwich", "xtable"
)

cat("# Environment\n\n")
cat(sprintf("- Recorded: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
cat(sprintf("- R: %s\n", R.version.string))
cat(sprintf("- Platform: %s\n", R.version$platform))
commit = tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA)
dirty = tryCatch(length(system("git status --porcelain", intern = TRUE)) > 0,
                 error = function(e) NA)
cat(sprintf("- Commit: %s%s\n\n",
            if (length(commit) == 1) commit else "unknown",
            if (isTRUE(dirty)) " (working tree has uncommitted changes)" else ""))

cat("| package | version |\n|---|---|\n")
for (package_name in REQUIRED_PACKAGES) {
  version = tryCatch(as.character(packageVersion(package_name)),
                     error = function(e) "NOT INSTALLED")
  cat(sprintf("| %s | %s |\n", package_name, version))
}

# The results files are the provenance record for every number in the paper, so
# their modification times and realised replication counts are pinned here too.
# A count other than 100,000 means a partial run was left in place.
RESULTS_FILES = c("normal.csv", "student.csv", "gamma.csv")
cat("\n| results file | produced | replications |\n|---|---|---|\n")
for (results_file in RESULTS_FILES) {
  if (!file.exists(results_file)) {
    cat(sprintf("| %s | ABSENT | - |\n", results_file))
    next
  }
  replications = tryCatch(
    paste(format(unique(data.table::fread(results_file)$num_trials),
                 scientific = FALSE, big.mark = ","), collapse = ", "),
    error = function(e) "-")
  cat(sprintf("| %s | %s | %s |\n", results_file,
              format(file.mtime(results_file), "%Y-%m-%d %H:%M:%S %Z"),
              replications))
}
