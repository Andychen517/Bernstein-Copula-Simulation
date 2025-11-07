# app.R
# Shiny app: Gaussian copula -> Bernstein copula (bivariate) via kdecopula
# This file wires UI + Server and sources helpers/plots.

# ---- packages ----
library(shiny)
library(copula)          # for Gaussian copula simulation
library(kdecopula)

# Source project modules
source("helpers.R")       # eval_cop_density_grid(), dnorm2()
source("plots.R")         # plot_* functions (use ggplot2:: internally)

# Optional sanity prints (helpful when teammates run on different machines)
cat("R version: ", R.version.string, "\n")
suppressWarnings({
  cat("kdecopula version: ",
      tryCatch(as.character(utils::packageVersion("kdecopula")),
               error = function(e) "NOT FOUND"), "\n")
  cat(".libPaths():\n"); print(.libPaths())
})

# ---- UI ----
ui <- fluidPage(
  titlePanel("Gaussian Copula → Bernstein Copula (via kdecopula)"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("rho", "Gaussian copula correlation (rho):",
                  min = -0.95, max = 0.95, value = 0.5, step = 0.05),
      numericInput("n", "Sample size (n):", value = 2000, min = 200, step = 500),
      numericInput("dimd", "Dimension d (Gaussian copula):", value = 2, min = 2, max = 10, step = 1),
      hr(),
      checkboxInput("auto_knots", "Auto-select Bernstein knots", value = TRUE),
      numericInput("knots", "Knots (if manual):", value = 12, min = 4, max = 60, step = 1),
      numericInput("gridN", "Density grid points per axis:", value = 150, min = 50, max = 300, step = 10),
      numericInput("bins", "bin2d bins (heatmap):", value = 60, min = 20, max = 150, step = 5),
      numericInput("kde_bins", "KDE contour levels:", value = 20, min = 5, max = 60, step = 1),
      actionButton("go", "Run / Refresh", class = "btn-primary")
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
                 plotOutput("bi_norm", height = 460))
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  rv <- reactiveValues(U = NULL, UV = NULL, fit = NULL, knots_used = NA_integer_)
  
  observeEvent(input$go, {
    withProgress(message = "Simulating & fitting...", value = 0, {
      set.seed(1)
      
      incProgress(0.25, detail = "Generating Gaussian copula data")
      # 1) Generate U ~ Gaussian copula (exchangeable)
      d   <- as.integer(input$dimd)
      rho <- input$rho
      n   <- as.integer(input$n)
      nc  <- copula::normalCopula(param = rho, dim = d, dispstr = "ex")
      U   <- copula::rCopula(n, nc)           # U in [0,1]^d
      UV  <- U[, 1:2, drop = FALSE]           # take first two margins (bivariate)
      
      incProgress(0.5, detail = "Selecting knots & fitting Bernstein (kdecopula)")
      # 2) Select Bernstein knots
      kts <- if (isTRUE(input$auto_knots)) {
        # Try internal (non-exported) bw_bern; if absent, use heuristic m ~ n^(1/3)
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
      cat("Knots used:", kts, "\n")
      
      # 3) Fit Bernstein copula via kdecopula
      fit <- kdecopula::kdecop(UV, method = "bern", knots = kts)
      
      rv$U <- U; rv$UV <- UV; rv$fit <- fit
      
      incProgress(0.25, detail = "Done")
    })
  }, ignoreInit = FALSE)
  
  # ---- Plots (use functions from plots.R) ----
  output$cop_scatter <- renderPlot({
    req(rv$UV)
    plot_cop_scatter(rv$UV)
  })
  
  output$cor_stats <- renderPrint({
    req(rv$UV)
    cat("Kendall's tau and Spearman's rho (rank correlations):\n")
    print(cor(rv$UV, method = "kendall"))  ## measures how consistently two variables move together.
    print(cor(rv$UV, method = "spearman"))  ## measures how well one variable can be described as a monotonic function of the other.
  })
  
  output$cop_bin2d <- renderPlot({
    req(rv$UV)
    plot_cop_bin2d(rv$UV, bins = input$bins)
  })
  
  output$cop_kde <- renderPlot({
    req(rv$UV)
    plot_cop_kde(rv$UV, bins = input$kde_bins)
  })
  
  output$bern_density <- renderPlot({
    req(rv$fit)
    plot_bern_density(rv$fit, N = as.integer(input$gridN),
                      title_extra = paste0("knots=", rv$knots_used))
  })
  
  output$bi_norm <- renderPlot({
    plot_bivar_normal(rho = input$rho)
  })
}

shinyApp(ui, server)
