################################################################################
#
#               SIMULADOR MENDELIANO PERDIDA HETEROCIGOSIS Y EFECTOS DE DERIVA 
#               Autor: Angel Criollo Rayo
#               Version 7: 01 - Mayo - 2025, asistido por ChatGPT
#               Grupo: Citogenetica, Filogenia y Evlucion de Poblaciones
#               Doctorado en Ciencias Biomedicas
#               Asignatura: Genetica Poblacional y Evolutiva (electiva)
#
################################################################################
# Este script esta diseñado para simular la transmision mendeliana de un ge con 
# dos alelos y analizar los efectos de la deriva en la perdida de la heterocigosidad
#  a lo largo de varias generaciones en varias poblaciones.
# Simula uniones o apareamientos por autofecundacion en distintos grados o
# fecundación cruzada. Tambien pueden simular cuellos de botella y efecto fundador
# Se puede variar la frecuencia inicial de los genotipos y grados de "selfing".
################################################################################
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- Helper Function for Mendelian Inheritance ---
sample_gamete <- function(allele1, allele2) {
  sample(c(allele1, allele2), size = 1)
}

# --- Core Simulation Function for One Population ---
simulate_population <- function(gen_sizes, n_gen,
                                initial_proportions = c(het = 1, homA1A1 = 0, homA2A2 = 0),
                                selfing_rate = 0.0) {
  pop <- vector("list", n_gen)
  
  if (sum(initial_proportions) != 1) stop("initial_proportions must sum to 1")
  if (selfing_rate < 0 || selfing_rate > 1) stop("selfing_rate must be between 0 and 1")
  
  n_ind_1 <- gen_sizes[1]
  types <- sample(c("het", "homA1A1", "homA2A2"), n_ind_1, replace = TRUE, prob = initial_proportions)
  allele1 <- allele2 <- character(n_ind_1)
  allele1[types == "het"] <- "A1"; allele2[types == "het"] <- "A2"
  allele1[types == "homA1A1"] <- "A1"; allele2[types == "homA1A1"] <- "A1"
  allele1[types == "homA2A2"] <- "A2"; allele2[types == "homA2A2"] <- "A2"
  
  pop[[1]] <- data.frame(id = 1:n_ind_1, allele1, allele2, parent1 = NA, parent2 = NA)
  
  for (g in 2:n_gen) {
    n_ind_g <- gen_sizes[g]
    prev_gen <- pop[[g - 1]]
    parent1 <- parent2 <- integer(n_ind_g)
    
    for (i in 1:n_ind_g) {
      if (runif(1) < selfing_rate) {
        p <- sample(1:nrow(prev_gen), 1)
        parent1[i] <- p; parent2[i] <- p
      } else {
        repeat {
          p1 <- sample(1:nrow(prev_gen), 1)
          p2 <- sample(1:nrow(prev_gen), 1)
          if (p1 != p2) break
        }
        parent1[i] <- p1
        parent2[i] <- p2
      }
    }

    allele1 <- mapply(function(p) sample_gamete(prev_gen$allele1[p], prev_gen$allele2[p]), parent1)
    allele2 <- mapply(function(p) sample_gamete(prev_gen$allele1[p], prev_gen$allele2[p]), parent2)
    
    pop[[g]] <- data.frame(id = 1:n_ind_g, allele1, allele2, parent1, parent2)
  }
  
  return(pop)
}

# --- Function to Simulate Multiple Populations ---
simulate_multiple_populations <- function(n_pops, gen_sizes_list, n_gen, list_props, list_selfing) {
  lapply(1:n_pops, function(i) {
    simulate_population(gen_sizes_list[[i]], n_gen, list_props[[i]], list_selfing[i])
  })
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
    freq_A <- mean(alleles == "A1")
    data.frame(generation = g, A1 = freq_A, A2 = 1 - freq_A, population = paste0("Pop", pop_id))
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
  titlePanel("Genetic Drift Simulator - Variable Gen Sizes"),
  sidebarLayout(
    sidebarPanel(
      numericInput("n_pops", "Number of Populations", min = 1, value = 3),
      sliderInput("n_gen", "Generations", min = 2, max = 50, value = 10),
      textAreaInput("prop_matrix", "Initial Genotype Proportions (het, homA1A1, homA2A2) per population (comma-separated rows)",
                    value = "0.5,0.25,0.25\n0.2,0.4,0.4\n1,0,0", rows = 4),
      textInput("selfing_vec", "Selfing Rates (comma-separated)", value = "0.1,0.5,0.9"),
      textAreaInput("gen_size_matrix", "Generation Sizes per Population (comma-separated rows)", 
                    value = "20,20,20,20,20,20,20,20,20,20\n25,25,25,25,25,25,25,25,25,25\n30,30,30,30,30,30,30,30,30,30", 
                    rows = 4),
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
      names(parts) <- c("het", "homA1A1", "homA2A2")
      parts / sum(parts)
    })
    
    selfing_rates <- as.numeric(strsplit(input$selfing_vec, ",")[[1]])
    
    gen_size_lines <- strsplit(input$gen_size_matrix, "\n")[[1]]
    gen_sizes_list <- lapply(gen_size_lines, function(line) as.numeric(strsplit(line, ",")[[1]]))
    
    # Validate inputs
    if (length(list_props) != input$n_pops || 
        length(selfing_rates) != input$n_pops || 
        length(gen_sizes_list) != input$n_pops) {
      showNotification("Number of rows in each input must match number of populations", type = "error")
      return(NULL)
    }
    
    if (any(sapply(gen_sizes_list, length) != input$n_gen)) {
      showNotification("Each generation size row must have length equal to number of generations", type = "error")
      return(NULL)
    }
    
    all_pops <- simulate_multiple_populations(input$n_pops, gen_sizes_list, input$n_gen, list_props, selfing_rates)
    
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
    
    # Set up positions and allele-based fill colors
    dat <- dat %>%
      mutate(
        x = id,
        y = -gen,
        fill1 = ifelse(allele1 == "A1", "blue", "red"),
        fill2 = ifelse(allele2 == "A1", "blue", "red")
      )
    
    ggplot() +
      # Add lines showing parent-offspring connections
      geom_segment(data = conn, aes(x = x_parent, y = y_parent, xend = x, yend = y), 
                   color = "gray70", linewidth = 0.3) +
      
      # Draw two small tiles per individual (one for each allele)
      geom_tile(data = dat, aes(x = x - 0.15, y = y, fill = fill1), width = 0.25, height = 0.4) +
      geom_tile(data = dat, aes(x = x + 0.15, y = y, fill = fill2), width = 0.25, height = 0.4) +
      
      scale_fill_identity() +
      facet_wrap(~population, scales = "free_x", ncol = input$ncol) +
      labs(title = "Genotype Visualization", x = "Individual", y = "Generation") +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank())
  })
  
  output$freq_plot <- renderPlot({
    req(sim_data())
    df <- sim_data()$freq_data %>%
      pivot_longer(cols = c("A1", "A2"), names_to = "allele", values_to = "freq")
    
    ggplot(df, aes(x = generation, y = freq, color = allele)) +
      geom_line(size = 1.2) +
      facet_wrap(~population, ncol = input$ncol) +
      scale_color_manual(values = c("A1" = "blue", "A2" = "red")) +
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

