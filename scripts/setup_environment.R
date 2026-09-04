# Installs the packages the experiments need.
#
# Run from the repository root, once, in a fresh R installation:
#   Rscript scripts/setup_environment.R
#
# On macOS the Posit snapshot repository serves source only, and building arrow
# and stringi from source takes far longer than the experiments themselves, so
# the default here is the CRAN binary repository. Reproducibility is then
# recorded rather than pinned: run scripts/session_info.R after the install and
# commit the resulting ENVIRONMENT.md, which lists the exact resolved versions.
# Set YEAST_CRAN_SNAPSHOT to a date to use a frozen source snapshot instead.

CRAN_SNAPSHOT_DATE = Sys.getenv("YEAST_CRAN_SNAPSHOT", "")
REPOSITORY = if (nzchar(CRAN_SNAPSHOT_DATE)) {
  sprintf("https://packagemanager.posit.co/cran/%s", CRAN_SNAPSHOT_DATE)
} else {
  "https://cloud.r-project.org"
}

REQUIRED_PACKAGES = c(
  "data.table",
  "dplyr",
  "ldbounds",    # alpha-spending boundaries for GST and Lan-DeMets OBF
  "mvtnorm",     # bivariate normal probabilities in the pYEAST bound
  "progress",
  "R6",
  "readxl",      # reading the Online Retail workbook
  "sandwich",    # cluster-robust variance estimation
  "xtable"
)

missing_packages = setdiff(REQUIRED_PACKAGES, rownames(installed.packages()))
if (length(missing_packages) == 0) {
  cat("All required packages are already installed.\n")
} else {
  cat(sprintf("Installing from %s:\n  %s\n",
              REPOSITORY, paste(missing_packages, collapse = ", ")))
  install_type = if (nzchar(CRAN_SNAPSHOT_DATE)) "source" else getOption("pkgType")
  install.packages(missing_packages, repos = REPOSITORY, type = install_type)
}

still_missing = setdiff(REQUIRED_PACKAGES, rownames(installed.packages()))
if (length(still_missing) > 0) {
  stop(sprintf("Failed to install: %s", paste(still_missing, collapse = ", ")))
}
cat("Environment ready. Record the resolved versions with:\n")
cat("  Rscript scripts/session_info.R > ENVIRONMENT.md\n")
