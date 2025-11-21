# app.R
# Gaussian Copula → Bernstein Copula (bivariate) via kdecopula
library(shiny)
library(copula)
library(kdecopula)
library(ggplot2)

source("helpers_n.R")
source("plots_n.R")

cat("R version: ", R.version.string, "\n")
suppressWarnings({
  cat("kdecopula version: ",
      tryCatch(as.character(utils::packageVersion("kdecopula")),
               error = function(e) "NOT FOUND"), "\n")
  cat(".libPaths():\n"); print(.libPaths())
})

ui <- fluidPage(
  titlePanel("Gaussian Copula → Bernstein Copula (via kdecopula)"),
  sidebarLayout(
    sidebarPanel(
      # ---- Simulation parameters ----
      sliderInput("rho", "Gaussian copula correlation (rho):",
                  min = -0.95, max = 0.95, value = 0.5, step = 0.05),
      numericInput("n", "Sample size (n):", value = 2000, min = 200, step = 500),
      numericInput("dimd", "Dimension d (Gaussian copula):", value = 2, min = 2, max = 10, step = 1),
      hr(),
      # ---- Bernstein parameters ----
      checkboxInput("auto_knots", "Auto-select Bernstein knots", value = TRUE),
      numericInput("knots", "Knots (if manual):", value = 12, min = 4, max = 60, step = 1),
      numericInput("gridN", "Density grid points per axis:", value = 150, min = 50, max = 300, step = 10),
      numericInput("bins", "bin2d bins (heatmap):", value = 60, min = 20, max = 150, step = 5),
      numericInput("kde_bins", "KDE contour levels:", value = 20, min = 5, max = 60, step = 1),
      actionButton("go", "Run / Refresh", class = "btn-primary"),
      hr(),
      h4("RMSE / Grid comparison (single dataset)"),
      numericInput("gridK", "Grid K (cells per axis):", value = 16, min = 6, max = 40, step = 2),
      numericInput("m_min", "Knots m: min", value = 4, min = 2, max = 80, step = 1),
      numericInput("m_max", "Knots m: max", value = 40, min = 4, max = 120, step = 1),
      numericInput("m_by",  "Knots m: step", value = 2, min = 1, max = 20, step = 1),
      actionButton("go_rmse", "Compute RMSE sweep", class = "btn-secondary"),
      hr(),
      h4("Replicates & Hypothesis Test (NEW)"),
      numericInput("R_reps", "Replicates R:", value = 100, min = 10, step = 10),
      radioButtons("eval_scheme", "RMSE evaluation points:",
                   c("midpoint grid (K_eval x K_eval)" = "mid",
                     "random N_eval points"            = "rnd")),
      conditionalPanel(
        "input.eval_scheme == 'mid'",
        numericInput("K_eval", "K_eval (cells each axis):", value = 25, min = 6, max = 100, step = 1)
      ),
      conditionalPanel(
        "input.eval_scheme == 'rnd'",
        numericInput("N_eval", "N_eval (random points):", value = 2000, min = 200, step = 200)
      ),
      numericInput("seed_rep", "Random seed:", value = 1, min = 1, step = 1),
      actionButton("go_reps", "Run R replicates", class = "btn-warning"),
      hr(),
      h5("Hypothesis test (m* vs alternative m)"),
      numericInput("m_alt", "Alternative m (compare to m*):", value = 12, min = 2, step = 1),
      checkboxInput("use_boot", "Bootstrap mean diff CI", value = TRUE),
      numericInput("B_boot", "Bootstrap B:", value = 1000, min = 100, step = 100)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Copula Scatter",
                 plotOutput("cop_scatter", height = 420),
                 verbatimTextOutput("cor_stats")),
        tabPanel("Copula Heat (bin2d)",
                 plotOutput("cop_bin2d", height = 420)),
        tabPanel("Copula KDE Contours",
                 plotOutput("cop_kde", height = 420)),
        tabPanel("Bernstein Density (kdecopula, [0,1]^2)",
                 plotOutput("bern_density", height = 460)),
        tabPanel("Bivariate Normal (Contours)",
                 plotOutput("bi_norm", height = 460)),
        tabPanel("RMSE vs m (single dataset, grid integration)",
                 plotOutput("rmse_curve", height = 420),
                 verbatimTextOutput("rmse_best")),
        tabPanel("Grid Error Heatmap",
                 plotOutput("cell_error", height = 420),
                 verbatimTextOutput("grid_stats")),
        tabPanel("Replicate Mean RMSE (NEW)",
                 plotOutput("rmse_mean_curve", height = 420),
                 verbatimTextOutput("rep_summary")),
        tabPanel("Hypothesis Test (NEW)",
                 verbatimTextOutput("htest_out"))
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    U = NULL, UV = NULL, fit = NULL, knots_used = NA_integer_,
    rmse_df = NULL, m_best = NULL, grid_eval = NULL, grid_stats = NULL,
    reps_curve = NULL, reps_mstar = NULL, rmse_mat = NULL, reps_mstar_mean = NULL
  )
  
  # ====== Single dataset RMSE sweep (grid integration, your original logic) ======
  observeEvent(input$go_rmse, {
    req(rv$UV)
    withProgress(message = "RMSE sweep over knots m ...", value = 0, {
      rho  <- input$rho
      K    <- as.integer(input$gridK)
      mseq <- seq(as.integer(input$m_min), as.integer(input$m_max), by = as.integer(input$m_by))
      incProgress(0.2, detail = sprintf("Evaluating %d m values", length(mseq)))
      df_rmse <- rmse_sweep_m(rv$UV, rho, mseq, K = K)
      rv$rmse_df <- df_rmse
      best_row <- df_rmse[which.min(df_rmse$rmse), , drop = FALSE]
      rv$m_best <- as.integer(best_row$m)
      
      # Fit at m* and grid stats + cell χ²
      fit_best <- kdecopula::kdecop(rv$UV, method = "bern", knots = rv$m_best)
      g_best   <- grid_eval_copula(fit_best, rho = rho, K = K)
      stats    <- grid_rmse(g_best)
      p_emp <- cell_probs_from_samples(rv$UV, K = K)
      p_hat <- cell_probs_from_fit(fit_best, K = K)
      chi2  <- chisq_cells(p_hat, p_emp)
      rv$grid_eval  <- g_best
      rv$grid_stats <- list(rmse = stats$rmse, ise = stats$ise, chisq = chi2, K = K, m = rv$m_best)
      incProgress(0.8, detail = "Done")
    })
  })
  
  # ====== Simulate & fit once ======
  observeEvent(input$go, {
    withProgress(message = "Simulating & fitting...", value = 0, {
      set.seed(1)
      incProgress(0.25, detail = "Generating Gaussian copula data")
      d   <- as.integer(input$dimd)
      rho <- input$rho
      n   <- as.integer(input$n)
      nc  <- copula::normalCopula(param = rho, dim = d, dispstr = "ex")
      U   <- copula::rCopula(n, nc)
      UV  <- U[, 1:2, drop = FALSE]
      incProgress(0.5, detail = "Selecting knots & fitting Bernstein (kdecopula)")
      kts <- if (isTRUE(input$auto_knots)) {
        ns <- getNamespace("kdecopula")
        if ("bw_bern" %in% ls(ns, all.names = TRUE)) {
          getFromNamespace("bw_bern", "kdecopula")(UV)
        } else {
          max(8L, min(60L, as.integer(round(nrow(UV)^(1/3)))))
        }
      } else {
        as.integer(input$knots)
      }
      rv$knots_used <- kts
      fit <- kdecopula::kdecop(UV, method = "bern", knots = kts)
      rv$U <- U; rv$UV <- UV; rv$fit <- fit
      incProgress(0.25, detail = "Done")
    })
  }, ignoreInit = FALSE)
  
  # ====== Plots (original) ======
  output$cop_scatter <- renderPlot({ req(rv$UV); plot_cop_scatter(rv$UV) })
  output$cor_stats <- renderPrint({
    req(rv$UV)
    cat("Kendall's tau and Spearman's rho (rank correlations):\n")
    print(cor(rv$UV, method = "kendall"))
    print(cor(rv$UV, method = "spearman"))
  })
  output$cop_bin2d   <- renderPlot({ req(rv$UV); plot_cop_bin2d(rv$UV, bins = input$bins) })
  output$cop_kde     <- renderPlot({ req(rv$UV); plot_cop_kde(rv$UV, bins = input$kde_bins) })
  output$bern_density<- renderPlot({ req(rv$fit); plot_bern_density(rv$fit, N = input$gridN,
                                                                    title_extra = paste0("knots=", rv$knots_used)) })
  output$bi_norm     <- renderPlot({ plot_bivar_normal(rho = input$rho) })
  
  output$rmse_curve <- renderPlot({
    req(rv$rmse_df)
    df <- rv$rmse_df
    ggplot(df, aes(m, rmse)) +
      geom_line() + geom_point() +
      { if (!is.null(rv$m_best)) geom_vline(xintercept = rv$m_best, linetype = 2) } +
      labs(title = sprintf("Grid-RMSE vs Bernstein knots (K=%d)", input$gridK),
           x = "knots m", y = "RMSE") +
      theme_minimal()
  })
  output$rmse_best <- renderPrint({
    req(rv$rmse_df, rv$m_best)
    best <- rv$rmse_df[which.min(rv$rmse_df$rmse), , drop = FALSE]
    cat("Best m (min RMSE):\n"); print(best)
  })
  output$cell_error <- renderPlot({
    req(rv$grid_eval)
    g <- rv$grid_eval
    ggplot(g, aes(U, V, fill = diff2)) +
      geom_tile() +
      coord_fixed() + xlim(0,1) + ylim(0,1) +
      labs(title = sprintf("Cell-wise squared error (m=%d, K=%d)", rv$grid_stats$m, rv$grid_stats$K),
           x = "U", y = "V", fill = "(ĉ − c_true)^2") +
      theme_minimal()
  })
  output$grid_stats <- renderPrint({
    req(rv$grid_stats)
    with(rv$grid_stats, {
      cat(sprintf("RMSE  = %.6f\nISE   = %.6f\nChi^2  = %.6f\nK      = %d\nm*     = %d\n",
                  rmse, ise, chisq, K, m))
    })
  })
  
  # ====== NEW: R-replicate experiment & hypothesis test ======
  observeEvent(input$go_reps, {
    withProgress(message = "Running replicate experiment ...", value = 0, {
      rho  <- input$rho
      n    <- as.integer(input$n)
      mseq <- seq(as.integer(input$m_min), as.integer(input$m_max), by = as.integer(input$m_by))
      R    <- as.integer(input$R_reps)
      
      # build evaluation points as requested
      eval_scheme <- input$eval_scheme
      if (eval_scheme == "mid") {
        UV_eval <- make_eval_points(K_eval = as.integer(input$K_eval), scheme = "midpoint")
      } else {
        UV_eval <- make_eval_points(N_eval = as.integer(input$N_eval), scheme = "random")
      }
      
      incProgress(0.2, detail = "Simulating & computing RMSE curves across replicates")
      out <- replicate_rmse_experiment(
        R = R, n = n, rho = rho, m_vec = mseq, UV_eval = UV_eval,
        dimd = as.integer(input$dimd), seed = as.integer(input$seed_rep)
      )
      rv$reps_curve <- out$curve
      rv$reps_mstar <- out$m_star
      rv$rmse_mat   <- out$rmse_mat
      rv$reps_mstar_mean <- mean(out$m_star)
      
      incProgress(0.7, detail = "Hypothesis test prep")
      # default hypo test: m0 = argmin mean RMSE; m1 = input$m_alt
      m0 <- out$curve$m[ which.min(out$curve$mean_rmse) ]
      rv$htest_res <- compare_m_within_reps(
        rmse_mat = out$rmse_mat, m_vec = mseq, m0 = m0, m1 = as.integer(input$m_alt),
        B = if (isTRUE(input$use_boot)) as.integer(input$B_boot) else 0L,
        conf = 0.95, seed = as.integer(input$seed_rep)
      )
      rv$m0_best <- m0
      incProgress(1.0, detail = "Done")
    })
  })
  
  output$rmse_mean_curve <- renderPlot({
    req(rv$reps_curve)
    plot_mean_rmse_curve(rv$reps_curve) +
      { if (!is.null(rv$m0_best)) geom_vline(xintercept = rv$m0_best, linetype = 2) }
  })
  
  output$rep_summary <- renderPrint({
    req(rv$reps_curve, rv$reps_mstar)
    cat("Replicate summary (point-sampled RMSE as in professor's definition):\n")
    m0 <- rv$m0_best
    cat(sprintf("Best m by mean RMSE (across R replicates): m* = %d\n", m0))
    cat(sprintf("Mean of per-replicate argmin m*: %.3f\n", rv$reps_mstar_mean))
    cat("Distribution of m* (table):\n")
    print(table(rv$reps_mstar))
    cat("\nHead of mean RMSE curve:\n")
    print(head(rv$reps_curve))
  })
  
  output$htest_out <- renderPrint({
    req(rv$rmse_mat, rv$m0_best, rv$htest_res)
    cat("Hypothesis test: compare alternative m to best m* (from mean RMSE)\n")
    cat(sprintf("m* (by mean RMSE) = %d;  alternative m = %d\n", rv$m0_best, as.integer(input$m_alt)))
    h <- rv$htest_res
    cat(sprintf("\nPaired differences d = RMSE(m_alt) - RMSE(m*):\n"))
    cat(sprintf("Mean diff = %.6f, SE = %.6f, t = %.3f (df=%d), p = %.4g\n",
                h$mean_diff, h$se, h$t, h$df, h$p))
    if (!is.null(h$boot_mean_ci)) {
      cat(sprintf("Bootstrap %.0f%% CI for mean diff: [%.6f, %.6f]\n",
                  100*h$boot_conf, h$boot_mean_ci[1], h$boot_mean_ci[2]))
    }
    cat("\nInterpretation: if mean diff > 0 and CI excludes 0 (or p < 0.05),\n")
    cat("then m* is significantly better than the alternative m (lower RMSE).\n")
  })
}

shinyApp(ui, server)
