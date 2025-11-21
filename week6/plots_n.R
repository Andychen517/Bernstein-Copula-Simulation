# plots.R
# All plotting functions return ggplot objects.

# 1) Copula scatter in [0,1]^2
plot_cop_scatter <- function(UV) {
  df <- data.frame(U = UV[,1], V = UV[,2])
  ggplot2::ggplot(df, ggplot2::aes(U, V)) +
    ggplot2::geom_point(alpha = 0.25, size = 0.6) +
    ggplot2::coord_fixed() + ggplot2::xlim(0,1) + ggplot2::ylim(0,1) +
    ggplot2::labs(title = "Copula samples (U,V) in [0,1]^2", x = "U", y = "V") +
    ggplot2::theme_minimal()
}

# 2) 2D binned heatmap on [0,1]^2
plot_cop_bin2d <- function(UV, bins = 60) {
  df <- data.frame(U = UV[,1], V = UV[,2])
  ggplot2::ggplot(df, ggplot2::aes(U, V)) +
    ggplot2::geom_bin2d(bins = bins) +
    ggplot2::coord_fixed() + ggplot2::xlim(0,1) + ggplot2::ylim(0,1) +
    ggplot2::labs(title = "Copula heatmap (2D bin counts)", x = "U", y = "V") +
    ggplot2::theme_minimal()
}

# 3) KDE contours (visual shape; not the Bernstein fit)
plot_cop_kde <- function(UV, bins = 20) {
  df <- data.frame(U = UV[,1], V = UV[,2])
  ggplot2::ggplot(df, ggplot2::aes(U, V)) +
    ggplot2::stat_density_2d(geom = "contour",
                             contour_var = "density",
                             bins = bins) +
    ggplot2::coord_fixed() + ggplot2::xlim(0,1) + ggplot2::ylim(0,1) +
    ggplot2::labs(title = "Copula KDE contours", x = "U", y = "V") +
    ggplot2::theme_minimal()
}

# 4) Bernstein copula pointwise density on [0,1]^2 (from kdecopula fit)
plot_bern_density <- function(fit, N = 150, title_extra = NULL) {
  gr <- eval_cop_density_grid(fit, N = N)
  df <- gr$df
  ttl <- "Bernstein copula density (method='bern')"
  if (!is.null(title_extra)) ttl <- paste0(ttl, " — ", title_extra)
  ggplot2::ggplot(df, ggplot2::aes(U, V, z = density)) +
    ggplot2::geom_contour_filled(bins = 12) +
    ggplot2::coord_fixed() + ggplot2::xlim(0,1) + ggplot2::ylim(0,1) +
    ggplot2::labs(title = ttl, x = "U", y = "V", fill = "dens") +
    ggplot2::theme_minimal()
}

# 5) Bivariate normal contours in (X,Y) = (Φ^{-1}(U), Φ^{-1}(V))
plot_bivar_normal <- function(rho, xlim = c(-3,3), ylim = c(-3,3), N = 200) {
  xs <- seq(xlim[1], xlim[2], length.out = N)
  ys <- seq(ylim[1], ylim[2], length.out = N)
  g  <- expand.grid(X = xs, Y = ys)
  g$z <- dnorm2(g$X, g$Y, rho = rho)
  ggplot2::ggplot(g, ggplot2::aes(X, Y, z = z)) +
    ggplot2::geom_contour(bins = 15) +
    ggplot2::labs(title = "Bivariate normal contours (standard margins)",
                  x = "X = Φ^{-1}(U)", y = "Y = Φ^{-1}(V)") +
    ggplot2::theme_minimal()
}

# 6) Mean RMSE curve with SE ribbon (from replicate study)
plot_mean_rmse_curve <- function(df_curve) {
  # df_curve: columns m, mean_rmse, se
  ggplot2::ggplot(df_curve, ggplot2::aes(x = m, y = mean_rmse)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_rmse - 1.96*se,
                                      ymax = mean_rmse + 1.96*se),
                         alpha = 0.2) +
    ggplot2::labs(title = "Mean RMSE vs m (with 95% ribbon from replicates)",
                  x = "knots m", y = "Mean RMSE") +
    ggplot2::theme_minimal()
}
