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
# ENVIRONMENT.md is this script's own redirect target, so it is already dirty by
# the time git runs and must be excluded or the flag is always set.
dirty = tryCatch(
  length(system("git status --porcelain -- . ':(exclude)ENVIRONMENT.md'",
                intern = TRUE)) > 0,
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
# the commit that last changed each one, and its realised replication count, are
# pinned here too. A count other than 100,000 means a partial run was left in
# place. We deliberately report the commit date rather than the file mtime: a
# checkout rewrites mtimes, so for a tracked file the mtime records when the
# working tree was last switched, not when the run was produced.
RESULTS_FILES = c("normal.csv", "student.csv", "gamma.csv",
                  "horizon_diagnostics.csv")
cat("\n| results file | last changed in | date | replications |\n",
    "|---|---|---|---|\n", sep = "")
for (results_file in RESULTS_FILES) {
  if (!file.exists(results_file)) {
    cat(sprintf("| %s | ABSENT | - | - |\n", results_file))
    next
  }
  replications = tryCatch(
    paste(format(unique(data.table::fread(results_file)$num_trials),
                 scientific = FALSE, big.mark = ","), collapse = ", "),
    error = function(e) "-")
  # The format string must be quoted: its field separator would otherwise be
  # read by the shell as a pipe.
  log_line = tryCatch(
    system(sprintf("git log -1 --format='%%h@%%cd' --date=format:'%%Y-%%m-%%d %%H:%%M %%Z' -- %s",
                   shQuote(results_file)), intern = TRUE),
    error = function(e) character(0))
  if (length(log_line) == 1 && nzchar(log_line)) {
    parts = strsplit(log_line, "@", fixed = TRUE)[[1]]
    commit = parts[1]
    changed = trimws(gsub("'", "", parts[2]))
  } else {
    commit = "uncommitted"
    changed = format(file.mtime(results_file), "%Y-%m-%d %H:%M %Z")
  }
  cat(sprintf("| %s | `%s` | %s | %s |\n", results_file, commit, changed,
              replications))
}
