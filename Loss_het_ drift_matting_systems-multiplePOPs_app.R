library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- Helper Function for Mendelian Inheritance ---
sample_gamete <- function(allele1, allele2) {
  sample(c(allele1, allele2), size = 1)
}

# --- Core Simulation Function for One Population ---
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

# --- Function to Simulate Multiple Populations with Custom Parameters ---
simulate_multiple_populations <- function(n_pops, n_ind, n_gen, list_props, list_selfing) {
  lapply(1:n_pops, function(i) simulate_population(n_ind, n_gen, list_props[[i]], list_selfing[i]))
}

get_plot_data <- function(pop, pop_id) {
  do.call(rbind, lapply(seq_along(pop), function(g) {
    df <- pop[[g]]
    data.frame(
      gen = g,
      id = df$id,
      allele1 = df$allele1,
      allele2 = df$allele2,
      parent1 = df$parent1,
      parent2 = df$parent2,
      population = paste0("Pop", pop_id)
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
    select(x, y, x_parent, y_parent, population)
}

get_allele_freqs <- function(pop, pop_id) {
  do.call(rbind, lapply(seq_along(pop), function(g) {
    alleles <- unlist(pop[[g]][, c("allele1", "allele2")])
    freq_A <- mean(alleles == "A")
    data.frame(generation = g, A = freq_A, B = 1 - freq_A, population = paste0("Pop", pop_id))
  }))
}

get_heterozygosity <- function(pop, pop_id) {
  do.call(rbind, lapply(seq_along(pop), function(g) {
    df <- pop[[g]]
    het <- mean(df$allele1 != df$allele2)
    data.frame(generation = g, heterozygosity = het, population = paste0("Pop", pop_id))
  }))
}

# --- UI ---
ui <- fluidPage(
  titlePanel("Genetic Drift Simulator - Multiple Populations with Custom Parameters"),
  sidebarLayout(
    sidebarPanel(
      numericInput("n_pops", "Number of Populations", min = 1, value = 3),
      sliderInput("n_ind", "Individuals per Population", min = 5, max = 100, value = 20),
      sliderInput("n_gen", "Generations", min = 5, max = 50, value = 15),
      textAreaInput("prop_matrix", "Initial Genotype Proportions (het, homAA, homBB) per population (comma-separated rows)",
                    value = "0.5,0.25,0.25\n0.2,0.4,0.4\n1,0,0", rows = 4),
      textInput("selfing_vec", "Selfing Rates (comma-separated)", value = "0.1,0.5,0.9"),
      numericInput("ncol", "Number of Columns in Grid Layout", value = 2, min = 1),
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
    prop_lines <- strsplit(input$prop_matrix, "\n")[[1]]
    list_props <- lapply(prop_lines, function(line) {
      parts <- as.numeric(strsplit(line, ",")[[1]])
      names(parts) <- c("het", "homAA", "homBB")
      parts / sum(parts)
    })
    
    selfing_rates <- as.numeric(strsplit(input$selfing_vec, ",")[[1]])
    
    if (length(list_props) != input$n_pops || length(selfing_rates) != input$n_pops) {
      showNotification("Number of rows and selfing rates must match number of populations", type = "error")
      return(NULL)
    }
    
    all_pops <- simulate_multiple_populations(input$n_pops, input$n_ind, input$n_gen, list_props, selfing_rates)
    
    plot_data <- do.call(rbind, lapply(seq_along(all_pops), function(i) get_plot_data(all_pops[[i]], i)))
    conn_data <- get_connection_data(plot_data)
    freq_data <- do.call(rbind, lapply(seq_along(all_pops), function(i) get_allele_freqs(all_pops[[i]], i)))
    het_data <- do.call(rbind, lapply(seq_along(all_pops), function(i) get_heterozygosity(all_pops[[i]], i)))
    
    list(
      plot_data = plot_data,
      conn_data = conn_data,
      freq_data = freq_data,
      het_data = het_data
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
      geom_segment(data = conn, aes(x = x_parent, y = y_parent, xend = x, yend = y, group = interaction(population, x)),
                   color = "gray70", linewidth = 0.4) +
      geom_point(data = dat, aes(x = a1_x, y = y, color = a1_col), size = 3) +
      geom_point(data = dat, aes(x = a2_x, y = y, color = a2_col), size = 3) +
      scale_color_identity() +
      facet_wrap(~population, scales = "free_x", ncol = input$ncol) +
      labs(title = "Genotype Grid", x = "Individual", y = "Generation") +
      theme_minimal()
  })
  
  output$freq_plot <- renderPlot({
    req(sim_data())
    df <- sim_data()$freq_data %>%
      pivot_longer(cols = c("A", "B"), names_to = "allele", values_to = "freq")
    
    ggplot(df, aes(x = generation, y = freq, color = allele)) +
      geom_line(size = 1.2) +
      facet_wrap(~population, ncol = input$ncol) +
      scale_color_manual(values = c("A" = "blue", "B" = "red")) +
      labs(title = "Allele Frequencies Over Time", y = "Frequency") +
      theme_minimal()
  })
  
  output$het_plot <- renderPlot({
    req(sim_data())
    df <- sim_data()$het_data
    
    ggplot(df, aes(x = generation, y = heterozygosity)) +
      geom_line(color = "darkgreen", size = 1.2) +
      facet_wrap(~population, ncol = input$ncol) +
      labs(title = "Heterozygosity Over Time", y = "Proportion Heterozygous") +
      theme_minimal()
  })
}

shinyApp(ui = ui, server = server)
