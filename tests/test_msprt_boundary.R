# Formula-level tests for the mSPRT / normal-mixture boundaries.
#
# Run from the repository root:  Rscript tests/test_msprt_boundary.R
#
# What is checked:
#  1. the implemented mSPRT boundary equals the closed form obtained by
#     inverting Robbins' normal mixture, in particular that the factor
#     (phi + n) / n sits INSIDE the square root;
#  2. the cumulative-sum form of the mSPRT boundary at tuning value phi is
#     identical to the GAVI normal-mixture boundary at rho = phi, i.e. the two
#     families in the paper's comparison are one family under two
#     parameterisations;
#  3. the GAVI tuning conversion maps a target time to rho as documented in the
#     manuscript (target time divided by approximately 6.33 at alpha = 0.05).

source("methods/msprt.R")
source("methods/gavi.R")

TOL <- 1e-9

# Closed form: with mixing variance tau^2 = sigma^2 / phi, the mixture
# likelihood ratio exceeds 1 / alpha_star exactly when
#   mean^2 > sigma^2 (phi + n) n^{-2} log((phi + n) / (phi alpha_star^2)).
robbins_radius <- function(n, phi, alpha, sigma) {
  alpha_star <- 2 * alpha
  sqrt(sigma ^ 2 * (phi + n) / n ^ 2 * log((phi + n) / (phi * alpha_star ^ 2)))
}

gavi_cumulative_boundary <- function(n, rho, alpha, sigma) {
  alpha_star <- 2 * alpha
  sigma * sqrt((n + rho) * log((n + rho) / (rho * alpha_star ^ 2)))
}

N <- 5000
n <- 1:N
grid <- expand.grid(
  phi = c(1, 11, 25, 100, 1580),
  alpha = c(0.01, 0.05, 0.1),
  sigma = c(0.5, 1, 10, 137.4)
)

max_impl_err <- 0
max_family_err <- 0

for (i in seq_len(nrow(grid))) {
  phi <- grid$phi[i]; alpha <- grid$alpha[i]; sigma <- grid$sigma[i]

  test <- mSPRT$new("test", alpha, sigma, phi)
  implemented <- test$mean_scale_boundary(N)
  stopifnot(length(implemented) == N)
  max_impl_err <- max(max_impl_err,
                      max(abs(implemented - robbins_radius(n, phi, alpha, sigma))))

  # The cumulative-sum boundary of the mSPRT at phi against GAVI at rho = phi.
  cumulative <- test$boundary(N)$value
  max_family_err <- max(
    max_family_err,
    max(abs(cumulative - gavi_cumulative_boundary(n, phi, alpha, sigma)))
  )
}

cat(sprintf("max |implemented - Robbins closed form| = %.3e\n", max_impl_err))
cat(sprintf("max |mSPRT(phi) - GAVI(rho = phi)|      = %.3e\n", max_family_err))
stopifnot(max_impl_err < TOL)
stopifnot(max_family_err < TOL * 1000)

# A guard against the parenthesisation error of putting (phi + n) / n outside
# the square root: that variant must differ from the implementation.
wrong <- (grid$phi[1] + n) / n * robbins_radius(n, grid$phi[1], 0.05, 1)
stopifnot(max(abs(wrong - mSPRT$new("t", 0.05, 1, grid$phi[1])$mean_scale_boundary(N))) > 1e-3)

# Documented tuning conversion used by the GAVI wrapper.
alpha <- 0.05
alpha_star <- 2 * alpha
divisor <- log(log(exp(1) * alpha_star ^ (-2))) - 2 * log(alpha_star)
cat(sprintf("alpha = %.3f: target time divided by %.2f\n", alpha, divisor))
stopifnot(abs(divisor - 6.33) < 0.01)
for (target in c(250, 500, 750, 10000)) {
  implemented_rho <- GAVI$new("t", alpha, 1, target)$rho()
  stopifnot(abs(implemented_rho - target / divisor) < TOL)
}
stopifnot(all(abs(
  vapply(c(250, 500, 750, 10000),
         function(target) GAVI$new("t", alpha, 1, target)$rho(), numeric(1))
  - c(39.5, 79.0, 118.5, 1580.1)
) < 0.15))

cat("test_msprt_boundary.R: all checks passed\n")
