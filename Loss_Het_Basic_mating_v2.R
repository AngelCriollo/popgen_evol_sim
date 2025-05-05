library(shiny)
library(ggplot2)
library(dplyr)
library(plotly)
library(tidyr)

ui <- fluidPage(
  titlePanel("Genetic Drift Simulation: Loss of Heterozygosity"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n_populations", "Number of Populations", value = 2, min = 1),
      numericInput("n_generations", "Number of Generations", value = 8, min = 2),
      numericInput("n_individuals", "Individuals per Generation", value = 20, min = 2),
      checkboxInput("selfing", "Allow Selfing", value = TRUE),
      selectInput("init_genotype", "Initial Genotype Configuration",
                  choices = c("Heterozygous (A/a)", "Random", "Custom Upload")),
      conditionalPanel(
        condition = "input.init_genotype == 'Custom Upload'",
        fileInput("custom_genotype_file", "Upload Custom Genotype CSV")
      ),
      actionButton("simulate", "Run Simulation"),
      br(),
      downloadButton("download_csv", "Download CSV"),
      downloadButton("download_plot", "Download Plot")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Allele Plot", uiOutput("plotUI")),
        tabPanel("Heterozygosity", plotOutput("heterozygosity_plot")),
        tabPanel("Animation", plotlyOutput("animation_plot"))
      )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$simulate, {
    
    n_populations <- input$n_populations
    n_generations <- input$n_generations
    n_individuals <- input$n_individuals
    allow_selfing <- input$selfing
    allele_colors <- c("A" = "blue", "a" = "red")
    
    all_data <- list()
    heterozygosity_summary <- data.frame()
    
    for (pop in 1:n_populations) {
      generations <- list()
      
      if (input$init_genotype == "Heterozygous (A/a)") {
        gen0 <- data.frame(id = 1:n_individuals, allele1 = "A", allele2 = "a", gen = 0)
      } else if (input$init_genotype == "Random") {
        genotypes <- sample(c("AA", "Aa", "aa"), n_individuals, replace = TRUE)
        allele1 <- substr(genotypes, 1, 1)
        allele2 <- substr(genotypes, 2, 2)
        gen0 <- data.frame(id = 1:n_individuals, allele1 = allele1, allele2 = allele2, gen = 0)
      } else {
        req(input$custom_genotype_file)
        gen0 <- read.csv(input$custom_genotype_file$datapath)
        stopifnot(all(c("id", "allele1", "allele2", "gen") %in% names(gen0)))
      }
      
      generations[[1]] <- gen0
      heterozygosity_summary <- rbind(heterozygosity_summary,
                                      data.frame(pop = pop, generation = 0,
                                                 heterozygosity = mean(gen0$allele1 != gen0$allele2)))
      
      for (g in 2:n_generations) {
        prev_gen <- generations[[g - 1]]
        new_gen <- data.frame(id = 1:n_individuals, allele1 = NA, allele2 = NA,
                              parent1 = NA, parent2 = NA, gen = g - 1)
        
        for (i in 1:n_individuals) {
          if (allow_selfing) {
            parents <- sample(prev_gen$id, 2, replace = TRUE)
          } else {
            parents <- sample(prev_gen$id, 2, replace = FALSE)
          }
          
          p1 <- prev_gen[prev_gen$id == parents[1], ]
          p2 <- prev_gen[prev_gen$id == parents[2], ]
          
          new_gen[i, "allele1"] <- sample(c(p1$allele1, p1$allele2), 1)
          new_gen[i, "allele2"] <- sample(c(p2$allele1, p2$allele2), 1)
          new_gen[i, "parent1"] <- p1$id
          new_gen[i, "parent2"] <- p2$id
        }
        
        heterozygosity_summary <- rbind(heterozygosity_summary,
                                        data.frame(pop = pop, generation = g - 1,
                                                   heterozygosity = mean(new_gen$allele1 != new_gen$allele2)))
        
        generations[[g]] <- new_gen
      }
      all_data[[pop]] <- generations
    }
    
    # Plot data prep
    plot_data <- data.frame()
    for (pop_idx in seq_along(all_data)) {
      gens <- all_data[[pop_idx]]
      for (g_idx in seq_along(gens)) {
        gen_data <- gens[[g_idx]]
        pd <- data.frame(pop = pop_idx, generation = gen_data$gen,
                         individual = gen_data$id,
                         allele1 = gen_data$allele1,
                         allele2 = gen_data$allele2,
                         parent1 = ifelse("parent1" %in% names(gen_data), gen_data$parent1, NA),
                         parent2 = ifelse("parent2" %in% names(gen_data), gen_data$parent2, NA))
        plot_data <- rbind(plot_data, pd)
      }
    }
    
    # Output heterozygosity plot
    output$heterozygosity_plot <- renderPlot({
      ggplot(heterozygosity_summary, aes(x = generation, y = heterozygosity, color = factor(pop))) +
        geom_line() + geom_point() +
        labs(title = "Heterozygosity Over Time", x = "Generation", y = "Heterozygosity") +
        theme_minimal()
    })
    
    # Animation plot
    output$animation_plot <- renderPlotly({
      long_data <- plot_data %>%
        mutate(a1x = individual - 0.2, a2x = individual + 0.2) %>%
        pivot_longer(cols = c("allele1", "allele2"), names_to = "allele", values_to = "value") %>%
        mutate(x = ifelse(allele == "allele1", a1x, a2x))
      
      plot_ly(long_data, x = ~x, y = ~generation, frame = ~generation, color = ~value,
              type = 'scatter', mode = 'markers', marker = list(size = 6)) %>%
        layout(yaxis = list(autorange = "reversed"), title = "Allele Transitions")
    })
    
    # Individual plots
    output$plotUI <- renderUI({
      lapply(1:n_populations, function(i) {
        plotname <- paste0("plot", i)
        plotOutput(plotname)
      })
    })
    
    for (i in 1:n_populations) {
      local({
        my_i <- i
        output[[paste0("plot", my_i)]] <- renderPlot({
          pdata <- plot_data %>% filter(pop == my_i)
          pdata <- pdata %>% mutate(x1 = individual - 0.2, x2 = individual + 0.2)
          
          ggplot() +
            geom_point(data = pdata, aes(x = x1, y = generation, color = allele1), size = 2.5) +
            geom_point(data = pdata, aes(x = x2, y = generation, color = allele2), size = 2.5) +
            scale_color_manual(values = allele_colors) +
            scale_y_reverse() +
            theme_minimal() +
            labs(title = paste("Population", my_i))
        })
      })
    }
    
    output$download_csv <- downloadHandler(
      filename = function() {
        paste("simulation_data_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(plot_data, file, row.names = FALSE)
      }
    )
    
    output$download_plot <- downloadHandler(
      filename = function() {
        paste("allele_plot_", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        ggsave(file, plot = last_plot(), device = "png")
      }
    )
  })
}

shinyApp(ui, server)

