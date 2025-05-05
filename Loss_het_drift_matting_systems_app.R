# chat gpt con shyni y mejoras: autofecundacion o cruzada

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- Helper Function for Mendelian Inheritance ---
sample_gamete <- function(allele1, allele2) {
  sample(c(allele1, allele2), size = 1)
}

# --- Core Simulation Function ---
simulate_population <- function(n_ind, n_gen,
                                initial_proportions = c(het = 1, homAA = 0, homBB = 0),
                                selfing_rate = 0.0) {
  pop <- vector("list", n_gen)
  
  if (sum(initial_proportions) != 1) stop("initial_proportions must sum to 1")
  if (selfing_rate < 0 || selfing_rate > 1) stop("selfing_rate must be between 0 and 1")
  
  types <- sample(c("het", "homAA", "homBB"), n_ind, replace = TRUE, prob = initial_proportions)
  allele1 <- allele2 <- character(n_ind)
  allele1[types == "het"] <- "A"; allele2[types == "het"] <- "B"
  allele1[types == "homAA"] <- "A"; allele2[types == "homAA"] <- "A"
  allele1[types == "homBB"] <- "B"; allele2[types == "homBB"] <- "B"
  
  pop[[1]] <- data.frame(id = 1:n_ind, allele1, allele2, parent1 = NA, parent2 = NA)
  
  for (g in 2:n_gen) {
    parent1 <- parent2 <- integer(n_ind)
    
    for (i in 1:n_ind) {
      if (runif(1) < selfing_rate) {
        p <- sample(1:n_ind, 1)
        parent1[i] <- p; parent2[i] <- p
      } else {
        parent1[i] <- sample(1:n_ind, 1)
        parent2[i] <- sample(1:n_ind, 1)
      }
    }
    
    allele1 <- mapply(function(p) {
      sample_gamete(pop[[g - 1]]$allele1[p], pop[[g - 1]]$allele2[p])
    }, parent1)
    
    allele2 <- mapply(function(p) {
      sample_gamete(pop[[g - 1]]$allele1[p], pop[[g - 1]]$allele2[p])
    }, parent2)
    
    pop[[g]] <- data.frame(id = 1:n_ind, allele1, allele2, parent1, parent2)
  }
  
  return(pop)
}

get_plot_data <- function(pop) {
  do.call(rbind, lapply(seq_along(pop), function(g) {
    df <- pop[[g]]
    data.frame(
      gen = g,
      id = df$id,
      allele1 = df$allele1,
      allele2 = df$allele2,
      parent1 = df$parent1,
      parent2 = df$parent2
    )
  }))
}

get_connection_data <- function(plot_data) {
  plot_data %>%
    filter(gen > 1) %>%
    mutate(
      x = id, y = -gen,
      x_parent1 = parent1, y_parent1 = -(gen - 1),
      x_parent2 = parent2, y_parent2 = -(gen - 1)
    ) %>%
    tidyr::pivot_longer(cols = c(x_parent1, x_parent2), names_to = "ptype", values_to = "x_parent") %>%
    mutate(y_parent = ifelse(ptype == "x_parent1", y_parent1, y_parent2)) %>%
    select(x, y, x_parent, y_parent)
}

get_allele_freqs <- function(pop) {
  do.call(rbind, lapply(seq_along(pop), function(g) {
    alleles <- unlist(pop[[g]][, c("allele1", "allele2")])
    freq_A <- mean(alleles == "A")
    data.frame(generation = g, A = freq_A, B = 1 - freq_A)
  }))
}

get_heterozygosity <- function(pop) {
  do.call(rbind, lapply(seq_along(pop), function(g) {
    df <- pop[[g]]
    het <- mean(df$allele1 != df$allele2)
    data.frame(generation = g, heterozygosity = het)
  }))
}

# --- UI ---
ui <- fluidPage(
  titlePanel("Genetic Drift Simulator"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("n_ind", "Individuals", min = 5, max = 100, value = 20),
      sliderInput("n_gen", "Generations", min = 5, max = 50, value = 15),
      sliderInput("selfing", "Selfing Rate", min = 0, max = 1, value = 0.5, step = 0.05),
      numericInput("prop_het", "Proportion Heterozygous", value = 0.5, min = 0, max = 1),
      numericInput("prop_homAA", "Proportion AA", value = 0.25, min = 0, max = 1),
      numericInput("prop_homBB", "Proportion BB", value = 0.25, min = 0, max = 1),
      actionButton("run", "Run Simulation")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Genotype Grid", plotOutput("genotype_plot", height = "600px")),
        tabPanel("Allele Frequencies", plotOutput("freq_plot")),
        tabPanel("Heterozygosity", plotOutput("het_plot"))
      )
    )
  )
)

# --- Server ---
server <- function(input, output) {
  sim_data <- eventReactive(input$run, {
    props <- c(het = input$prop_het, homAA = input$prop_homAA, homBB = input$prop_homBB)
    sim <- simulate_population(input$n_ind, input$n_gen, props, selfing_rate = input$selfing)
    list(
      pop = sim,
      plot_data = get_plot_data(sim),
      conn_data = get_connection_data(get_plot_data(sim)),
      freq_data = get_allele_freqs(sim),
      het_data = get_heterozygosity(sim)
    )
  })
  
  output$genotype_plot <- renderPlot({
    req(sim_data())
    dat <- sim_data()$plot_data
    conn <- sim_data()$conn_data
    
    dat <- dat %>%
      mutate(
        x = id, y = -gen,
        a1_x = id - 0.2, a2_x = id + 0.2,
        a1_col = ifelse(allele1 == "A", "blue", "red"),
        a2_col = ifelse(allele2 == "A", "blue", "red")
      )
    
    ggplot() +
      geom_segment(data = conn, aes(x = x_parent, y = y_parent, xend = x, yend = y),
                   color = "gray70", linewidth = 0.4) +
      geom_point(data = dat, aes(x = a1_x, y = y, color = a1_col), size = 3) +
      geom_point(data = dat, aes(x = a2_x, y = y, color = a2_col), size = 3) +
      scale_color_identity() +
      labs(title = "Genotype Grid", x = "Individual", y = "Generation") +
      theme_minimal()
  })
  
  output$freq_plot <- renderPlot({
    req(sim_data())
    df <- sim_data()$freq_data %>%
      pivot_longer(cols = c("A", "B"), names_to = "allele", values_to = "freq")
    
    ggplot(df, aes(x = generation, y = freq, color = allele)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("A" = "blue", "B" = "red")) +
      labs(title = "Allele Frequencies Over Time", y = "Frequency") +
      theme_minimal()
  })
  
  output$het_plot <- renderPlot({
    req(sim_data())
    df <- sim_data()$het_data
    
    ggplot(df, aes(x = generation, y = heterozygosity)) +
      geom_line(color = "darkgreen", size = 1.2) +
      labs(title = "Heterozygosity Over Time", y = "Proportion Heterozygous") +
      theme_minimal()
  })
  
}

shinyApp(ui = ui, server = server)

