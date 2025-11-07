# helpers.R
# Utility helpers (no plotting here). Comments are in English.

# Evaluate fitted Bernstein copula density on an [0,1]^2 grid (kdecopula fit)
eval_cop_density_grid <- function(fit, N = 150) {
  # fit: object from kdecopula::kdecop(method = "bern")
  # N:   grid resolution per axis (N x N)
  u <- v <- seq(0, 1, length.out = N)
  grid <- expand.grid(U = u, V = v)
  # Namespace call to avoid masking issues
  dens <- kdecopula::dkdecop(as.matrix(grid), fit)
  list(u = u, v = v, df = transform(grid, density = dens))
}

# Bivariate normal pdf with correlation rho (standard normals, zero means)
dnorm2 <- function(x, y, rho) {
  z <- (x^2 - 2*rho*x*y + y^2) / (1 - rho^2)
  (1 / (2*pi*sqrt(1 - rho^2))) * exp(-0.5 * z)
}

