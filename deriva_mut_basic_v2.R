################################################################################
#
#               SIMULADOR DE DERIVA Y MUTACION
#               Autor: Angel Criollo Rayo
#               Version 7: 23 - mayo - 2025
#               Grupo: Citogenetica, Filogenia y Evlucion de Poblaciones
#               Doctorado en Ciencias Biomedicas
#               Asignatura: Genetica Poblacional y Evolutiva (electiva)
#
################################################################################
library(shiny)
library(ggplot2)
library(reshape2)

simulate_drift <- function(num_pops, num_loci, generations, N_vec, p0_vec, mu_fwd_vec, mu_rev_vec) {
  results <- array(NA, dim = c(generations + 1, num_loci, num_pops),
                   dimnames = list(NULL, paste0("Locus", 1:num_loci), paste0("Pop", 1:num_pops)))
  
  for (pop in 1:num_pops) {
    N <- N_vec[pop]
    for (locus in 1:num_loci) {
      p <- numeric(generations + 1)
      p[1] <- p0_vec[locus]
      mu_fwd <- mu_fwd_vec[locus]
      mu_rev <- mu_rev_vec[locus]
      for (gen in 2:(generations + 1)) {
        # Drift step
        p_drifted <- rbinom(1, 2 * N, p[gen - 1]) / (2 * N)
        # Asymmetric mutation step
        p[gen] <- p_drifted * (1 - mu_fwd) + (1 - p_drifted) * mu_rev
      }
      results[, locus, pop] <- p
    }
  }
  
  return(results)
}

ui <- fluidPage(
  titlePanel("Genetic Drift Simulator with Asymmetric Mutation"),
  sidebarLayout(
    sidebarPanel(
      textInput("pop_sizes", "Population Sizes (comma-separated):", "100,100,100"),
      textInput("p0_values", "Initial Frequencies per Locus (comma-separated):", "0.5,0.5,0.5"),
      textInput("mu_fwd", "Forward Mutation Rates A→a (comma-separated):", "0.001,0.001,0.001"),
      textInput("mu_rev", "Reverse Mutation Rates a→A (comma-separated):", "0.001,0.001,0.001"),
      numericInput("generations", "Number of Generations:", 100, min = 1),
      actionButton("simulate", "Run Simulation")
    ),
    mainPanel(
      plotOutput("driftPlot")
    )
  )
)

server <- function(input, output, session) {
  sim_data <- reactiveVal(NULL)
  
  observeEvent(input$simulate, {
    tryCatch({
      N_vec <- as.numeric(strsplit(input$pop_sizes, ",")[[1]])
      p0_vec <- as.numeric(strsplit(input$p0_values, ",")[[1]])
      mu_fwd_vec <- as.numeric(strsplit(input$mu_fwd, ",")[[1]])
      mu_rev_vec <- as.numeric(strsplit(input$mu_rev, ",")[[1]])
      
      num_pops <- length(N_vec)
      num_loci <- length(p0_vec)
      generations <- input$generations
      
      if (any(is.na(N_vec)) || any(is.na(p0_vec)) || any(is.na(mu_fwd_vec)) || any(is.na(mu_rev_vec))) 
        stop("Inputs must be numeric and comma-separated.")
      if (length(mu_fwd_vec) != num_loci || length(mu_rev_vec) != num_loci)
        stop("Number of mutation rates must match number of loci.")
      if (any(mu_fwd_vec < 0 | mu_fwd_vec > 1) || any(mu_rev_vec < 0 | mu_rev_vec > 1))
        stop("Mutation rates must be between 0 and 1.")
      if (any(N_vec <= 0)) stop("Population sizes must be positive.")
      if (any(p0_vec < 0 | p0_vec > 1)) stop("Initial frequencies must be between 0 and 1.")
      
      sim <- simulate_drift(num_pops, num_loci, generations, N_vec, p0_vec, mu_fwd_vec, mu_rev_vec)
      sim_data(sim)
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  output$driftPlot <- renderPlot({
    req(sim_data())
    sim <- sim_data()
    df <- melt(sim)
    colnames(df) <- c("Generation", "Locus", "Population", "AlleleFreq")
    
    ggplot(df, aes(x = Generation, y = AlleleFreq, color = Locus)) +
      geom_line(alpha = 0.8) +
      facet_wrap(~ Population) +
      theme_minimal() +
      labs(
        title = "Genetic Drift with Mutation",
        y = "Allele Frequency",
        caption = "Mutation model: p' = p*(1 - mu_fwd) + (1 - p)*mu_rev"
      )
  })
}

shinyApp(ui = ui, server = server)
