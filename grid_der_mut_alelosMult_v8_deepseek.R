################################################################################
#
#               DERIVA Y MUTACION: ALELOS MULTIPLES
#               Autor: Angel Criollo Rayo & DeepSeek
#               Version 2: 24 - Mayo - 2025
#               Grupo: Citogenetica, Filogenia y Evlucion de Poblaciones
#               Doctorado en Ciencias Biomedicas
#               Asignatura: Genetica Poblacional y Evolutiva (electiva)
#
################################################################################
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

ui <- fluidPage(
  titlePanel("Genetic Drift Simulation with Coalescence Analysis"),
  sidebarLayout(
    sidebarPanel(
      h4("Simulation Parameters"),
      numericInput("popSize", "Population size (Nₑ)", value = 50, min = 2, max = 200),
      numericInput("generations", "Number of generations (t)", value = 20, min = 1, max = 100),
      numericInput("nAlleles", "Number of initial alleles", value = 2, min = 1, max = 10),
      numericInput("mutationRate", "Mutation rate (µ)", value = 0.01, min = 0, max = 1, step = 0.01),
      actionButton("runSim", "Run Simulation", class = "btn-primary"),
      hr(),
      h4("Coalescence Calculator"),
      numericInput("calcGen", "Calculate up to generation", value = 10, min = 1, max = 100),
      actionButton("calcCoal", "Calculate Coalescence"),
      hr(),
      h4("Display Options"),
      checkboxInput("showPedigree", "Show pedigree connections", value = TRUE),
      sliderInput("pointSize", "Point size", min = 1, max = 5, value = 3)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Pedigree & Diversity",
                 div(style = "overflow-y: auto; height: 500px;",
                     plotOutput("plotSim", height = "800px")),
                 plotOutput("plotDiversity", height = "300px")),
        tabPanel("Allele Frequencies",
                 plotOutput("plotFrequencies", height = "600px")),
        tabPanel("Coalescence Analysis",
                 h4("Coalescence and Heterozygosity"),
                 withMathJax(uiOutput("coalFormulas")),
                 verbatimTextOutput("coalResults"),
                 plotOutput("coalPlot", height = "400px")),
        tabPanel("Data Summary",
                 downloadButton("downloadData", "Download Data"),
                 dataTableOutput("summaryTable"))
      )
    )
  )
)

server <- function(input, output) {
  simulation_data <- reactiveValues(
    pedigree = NULL,
    diversity = NULL,
    frequencies = NULL,
    new_id_count = 1,
    allele_colors = NULL,
    coal_results = NULL
  )
  
  # Coalescence formulas with heterozygosity
  output$coalFormulas <- renderUI({
    withMathJax(
      helpText("Key formulas:"),
      helpText("1. Coalescence probability:"),
      helpText("$$P_{\\text{coal}}(t) = \\left(1 - \\frac{1}{2N_e}\\right)^{t-1} \\cdot \\frac{1}{2N_e}$$"),
      helpText("2. Mutation probability:"),
      helpText("$$P_{\\text{mut}}(t) = 1 - (1 - \\mu)^{2t}$$"),
      helpText("3. Mean heterozygosity (equilibrium):"),
      helpText("$$H = \\frac{4N_e\\mu}{1 + 4N_e\\mu} = \\frac{\\theta}{1 + \\theta}$$"),
      helpText("Where:"),
      helpText("- \\(N_e\\) = effective population size"),
      helpText("- \\(\\mu\\) = mutation rate"),
      helpText("- \\(\\theta = 4N_e\\mu\\) = population mutation rate")
    )
  })
  
  # Coalescence calculator with heterozygosity
  observeEvent(input$calcCoal, {
    req(input$popSize, input$mutationRate, input$calcGen)
    
    Ne <- input$popSize
    μ <- input$mutationRate
    max_t <- input$calcGen
    
    # Calculate coalescence probabilities
    t <- 1:max_t
    coal_prob <- (1 - 1/(2*Ne))^(t-1) * (1/(2*Ne))
    
    # Calculate mutation probabilities
    no_mut_prob <- (1 - μ)^(2*t)
    mut_prob <- 1 - no_mut_prob
    
    # Calculate mean heterozygosity (equilibrium)
    θ <- 4 * Ne * μ
    mean_het <- θ / (1 + θ)
    
    simulation_data$coal_results <- data.frame(
      Generation = t,
      P_Coalesce = coal_prob,
      P_NoMutation = no_mut_prob,
      P_Mutation = mut_prob,
      P_Identical = coal_prob * no_mut_prob,
      Cumulative_Coal = cumsum(coal_prob),
      Theta = θ,
      Mean_Heterozygosity = mean_het
    )
  })
  
  # Coalescence results output
  output$coalResults <- renderPrint({
    req(simulation_data$coal_results)
    
    θ <- unique(simulation_data$coal_results$Theta)
    mean_het <- unique(simulation_data$coal_results$Mean_Heterozygosity)
    
    cat("Coalescence Analysis Results:\n")
    cat("----------------------------\n")
    cat("Effective population size (Nₑ):", input$popSize, "\n")
    cat("Mutation rate (µ):", input$mutationRate, "\n")
    cat("Population mutation rate (θ = 4Nₑµ):", round(θ, 4), "\n")
    cat("Mean heterozygosity (H = θ/(1+θ)):", round(mean_het, 4), "\n\n")
    cat("By generation (first 5):\n")
    print(head(simulation_data$coal_results, 5))
  })
  
  # Coalescence plot
  output$coalPlot <- renderPlot({
    req(simulation_data$coal_results)
    
    df <- simulation_data$coal_results %>%
      select(Generation, P_Coalesce, P_Mutation, P_Identical) %>%
      pivot_longer(-Generation, names_to = "Probability_Type", values_to = "Value")
    
    ggplot(df, aes(x = Generation, y = Value, color = Probability_Type)) +
      geom_line(size = 1.2) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_minimal() +
      labs(title = "Coalescence and Mutation Probabilities",
           y = "Probability",
           color = "Probability Type") +
      scale_color_manual(
        values = c("P_Coalesce" = "red", 
                   "P_Mutation" = "blue", 
                   "P_Identical" = "green"),
        labels = c("Coalescence", "Mutation", "Identical")
      ) +
      theme(legend.position = "bottom")
  })
  
  # Generate consistent color mapping for alleles
  generate_allele_colors <- function(alleles) {
    all_alleles <- unique(alleles)
    n <- length(all_alleles)
    
    # Create a consistent color mapping
    colors <- rainbow(n)
    names(colors) <- all_alleles
    
    # Sort colors by allele name for consistency
    colors <- colors[order(names(colors))]
    
    return(colors)
  }
  
  observeEvent(input$runSim, {
    # Reset simulation data
    simulation_data$new_id_count <- 1
    
    # Run simulation
    withProgress(message = 'Running simulation...', value = 0, {
      popSize <- input$popSize
      generations <- input$generations
      mutation_rate <- input$mutationRate
      allele_pool <- LETTERS[1:input$nAlleles]
      
      pop <- list()
      gen0 <- replicate(popSize, {
        alleles <- sample(allele_pool, 2, replace = TRUE)
        list(alleles = alleles, parents = NA)
      }, simplify = FALSE)
      pop[[1]] <- gen0
      
      # Initialize data storage
      pedigree_df <- data.frame()
      connections_df <- data.frame()
      diversity_df <- data.frame()
      all_alleles <- character(0)
      
      # Track all possible alleles across generations
      allele_tracker <- unique(allele_pool)
      
      for (g in 1:(generations + 1)) {
        incProgress(1/(generations + 1), detail = paste("Generation", g))
        
        if (g > 1) {
          prev_gen <- pop[[g - 1]]
          new_gen <- list()
          for (i in 1:popSize) {
            parents <- sample(seq_along(prev_gen), 2)
            p1 <- prev_gen[[parents[1]]]$alleles
            p2 <- prev_gen[[parents[2]]]$alleles
            offspring <- get_offspring(sample(p1, 1), sample(p2, 1),
                                       mutation_rate, 
                                       simulation_data)
            new_gen[[i]] <- list(alleles = offspring, parents = parents)
            
            # Track new alleles
            allele_tracker <- unique(c(allele_tracker, offspring))
          }
          pop[[g]] <- new_gen
        }
        
        # Collect data for current generation
        alleles_in_gen <- c()
        for (i in seq_along(pop[[g]])) {
          indiv <- pop[[g]][[i]]
          pedigree_df <- rbind(pedigree_df, data.frame(
            gen = g,
            id = paste(g, i, sep = "_"),
            x = i,
            allele1 = indiv$alleles[1],
            allele2 = indiv$alleles[2],
            stringsAsFactors = FALSE
          ))
          alleles_in_gen <- c(alleles_in_gen, indiv$alleles)
          
          if (g > 1 && !is.null(indiv$parents) && !is.na(indiv$parents[1])) {
            for (p in indiv$parents) {
              connections_df <- rbind(connections_df, data.frame(
                x = i,
                y = g,
                xend = p,
                yend = g - 1,
                stringsAsFactors = FALSE
              ))
            }
          }
        }
        
        # Calculate frequencies for ALL tracked alleles (including those at 0)
        freq_table <- table(factor(alleles_in_gen, levels = allele_tracker))/length(alleles_in_gen)
        
        diversity_df <- rbind(diversity_df, data.frame(
          generation = g,
          allele_count = length(unique(alleles_in_gen)),
          stringsAsFactors = FALSE
        ))
        
        # Store frequencies for this generation
        gen_freq <- data.frame(
          generation = g,
          allele = names(freq_table),
          frequency = as.numeric(freq_table),
          stringsAsFactors = FALSE
        )
        
        if (g == 1) {
          simulation_data$frequencies <- gen_freq
        } else {
          simulation_data$frequencies <- rbind(simulation_data$frequencies, gen_freq)
        }
        
        all_alleles <- unique(c(all_alleles, alleles_in_gen))
      }
      
      # Generate consistent color mapping
      simulation_data$allele_colors <- generate_allele_colors(c(allele_pool, 
                                                                paste0("N", 1:simulation_data$new_id_count)))
      
      # Store other results
      simulation_data$pedigree <- pedigree_df
      simulation_data$connections <- connections_df
      simulation_data$diversity <- diversity_df
    })
  })
  
  get_offspring <- function(p1, p2, mutation_rate, sim_data) {
    alleles <- c(p1, p2)
    mutate <- function(allele) {
      if (runif(1) < mutation_rate) {
        new_id <- paste0("N", sim_data$new_id_count)
        sim_data$new_id_count <- sim_data$new_id_count + 1
        return(new_id)
      }
      return(allele)
    }
    gametes <- sapply(sample(alleles, 2, replace = TRUE), mutate)
    return(gametes)
  }
  
  # Pedigree plot with consistent colors
  output$plotSim <- renderPlot({
    req(simulation_data$pedigree, simulation_data$allele_colors)
    
    df <- simulation_data$pedigree
    
    p <- ggplot()
    
    if (input$showPedigree && nrow(simulation_data$connections) > 0) {
      p <- p + geom_segment(
        data = simulation_data$connections,
        aes(x = xend, xend = x, y = yend, yend = y),
        color = "grey", alpha = 0.3
      )
    }
    
    p + geom_point(
      data = df,
      aes(x = x - 0.1, y = gen, color = allele1),
      size = input$pointSize
    ) +
      geom_point(
        data = df,
        aes(x = x + 0.1, y = gen, color = allele2),
        size = input$pointSize
      ) +
      scale_color_manual(values = simulation_data$allele_colors) +
      scale_y_reverse(expand = expansion(mult = c(0.05, 0.05))) +
      theme_minimal() +
      labs(x = "Individual", y = "Generation") +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        axis.text.y = element_text(size = 10)
      ) +
      guides(color = guide_legend(nrow = 2, override.aes = list(size = 4)))
  })
  
  # Frequency plot with consistent colors
  output$plotFrequencies <- renderPlot({
    req(simulation_data$frequencies, simulation_data$allele_colors)
    
    # Filter to show only alleles that appear at least once
    active_alleles <- simulation_data$frequencies %>%
      group_by(allele) %>%
      filter(any(frequency > 0)) %>%
      pull(allele) %>%
      unique()
    
    plot_data <- simulation_data$frequencies %>%
      filter(allele %in% active_alleles)
    
    ggplot(plot_data, aes(x = generation, y = frequency, color = allele)) +
      geom_line(size = 0.5) +
      geom_point(size = 1) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      scale_color_manual(values = simulation_data$allele_colors) +
      theme_minimal() +
      labs(x = "Generation", y = "Allele Frequency",
           title = "Allele Frequency Dynamics") +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 8),  # Smaller legend text
            legend.key.size = unit(0.5, "lines"),  # Smaller legend items
            legend.margin = margin(t = 0, unit = "cm")) +  # Tighter margins
      guides(color = guide_legend(nrow = 2,         # 2 rows for wider layout
                                  override.aes = list(size = 2)))
  })
  
  # Diversity plot
  output$plotDiversity <- renderPlot({
    req(simulation_data$diversity)
    
    ggplot(simulation_data$diversity, aes(x = generation, y = allele_count)) +
      geom_line(color = "darkgreen", size = 1.2) +
      geom_point(color = "darkgreen", size = 2) +
      theme_minimal() +
      labs(x = "Generation", y = "Number of Unique Alleles",
           title = "Allele Diversity Over Time")
  })
  
  # Data summary table
  output$summaryTable <- renderDataTable({
    req(simulation_data$frequencies)
    
    simulation_data$frequencies %>%
      spread(key = allele, value = frequency, fill = 0) %>%
      arrange(generation)
  })
  
  # Download handler
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("genetic-drift-data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(simulation_data$frequencies, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
