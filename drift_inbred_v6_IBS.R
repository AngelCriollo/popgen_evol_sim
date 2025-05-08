################################################################################
#
#               APP TRANSMISION MENDELIANA Y ENDOGAMIA
#               Autor: Angel Criollo Rayo 
#               Version 6: 07 - Mayo - 2025, asistido por ChatGPT
#               Grupo: Citogenetica, Filogenia y Evlucion de Poblaciones
#               Doctorado en Ciencias Biomedicas
#               Asignatura: Genetica Poblacional y Evolutiva (electiva)
#
################################################################################
#
# Este script esta diseñado para analizar la transmisión mendeliana de genes a lo
# largo de varias generaciones en una poblacion, bajo diferentes contextos de endo-
# gamia o apareamientos por autofecundacion en distintos grados o
# fecundación cruzada.
# Se puede variar la frecuencia inicial de los genotipos y grados de "selfing".
################################################################################
library(shiny)
library(DT)
library(ggplot2)
library(scales)

# Simulate multilocus Mendelian inheritance
generate_offspring_genotype <- function(mother, father, gene_names) {
  offspring <- ""
  for (i in seq_along(gene_names)) {
    m_alleles <- strsplit(substr(mother, 2 * i - 1, 2 * i), "")[[1]]
    f_alleles <- strsplit(substr(father, 2 * i - 1, 2 * i), "")[[1]]
    allele1 <- sample(m_alleles, 1)
    allele2 <- sample(f_alleles, 1)
    geno <- paste0(sort(c(allele1, allele2)), collapse = "")
    offspring <- paste0(offspring, geno)
  }
  offspring
}

get_phenotype <- function(genotype, gene_names) {
  phenotype <- c()
  for (i in seq_along(gene_names)) {
    alleles <- strsplit(substr(genotype, (2 * i - 1), 2 * i), "")[[1]]
    gene <- gene_names[i]
    if (any(alleles == toupper(gene))) {
      phenotype <- c(phenotype, paste0(toupper(gene), "_"))
    } else {
      phenotype <- c(phenotype, paste0(tolower(gene), tolower(gene)))
    }
  }
  paste(phenotype, collapse = " ")
}

summarize_phenotypes <- function(data, gene_names) {
  data$Phenotype <- sapply(data$Genotype, get_phenotype, gene_names = gene_names)
  summary <- as.data.frame(table(data$Generation, data$Phenotype))
  colnames(summary) <- c("Generation", "Phenotype", "Count")
  
  # Calculate proportions
  total_per_gen <- aggregate(Count ~ Generation, data = summary, sum)
  summary <- merge(summary, total_per_gen, by = "Generation", suffixes = c("", "_Total"))
  summary$Proportion <- summary$Count / summary$Count_Total
  summary <- summary[order(summary$Generation, summary$Phenotype), ]
  return(summary)
}

compute_multilocus_heterozygosity <- function(data, gene_names) {
  result <- data.frame()
  for (gen in unique(data$Generation)) {
    gen_data <- subset(data, Generation == gen)
    heterozygosities <- sapply(gen_data$Genotype, function(gt) {
      sum(sapply(seq_along(gene_names), function(i) {
        alleles <- strsplit(substr(gt, 2 * i - 1, 2 * i), "")[[1]]
        alleles[1] != alleles[2]
      }))
    })
    mlh <- heterozygosities / length(gene_names)
    result <- rbind(result, data.frame(
      Generation = gen,
      Mean_MLH = mean(mlh),
      Var_MLH = var(mlh),
      N = length(mlh)
    ))
  }
  return(result)
}

compute_ibs_matrix <- function(data, gene_names, generation = NULL) {
  if (!is.null(generation)) {
    data <- subset(data, Generation == generation)
  }
  
  n <- nrow(data)
  ibs_matrix <- matrix(NA, nrow = n, ncol = n)
  rownames(ibs_matrix) <- data$ID
  colnames(ibs_matrix) <- data$ID
  
  get_alleles <- function(geno, i) {
    strsplit(substr(geno, 2 * i - 1, 2 * i), "")[[1]]
  }
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      geno_i <- data$Genotype[i]
      geno_j <- data$Genotype[j]
      total_score <- 0
      max_score <- 2 * length(gene_names)
      
      for (g in seq_along(gene_names)) {
        a1 <- get_alleles(geno_i, g)
        a2 <- get_alleles(geno_j, g)
        shared <- length(intersect(a1, a2))
        total_score <- total_score + shared
      }
      
      score <- total_score / max_score
      ibs_matrix[i, j] <- score
      ibs_matrix[j, i] <- score
    }
    ibs_matrix[i, i] <- 1
  }
  ibs_matrix[n, n] <- 1
  return(ibs_matrix)
}

genotype_grid_plot <- function(data, gene_names, show_labels = TRUE) {
  offset <- 0.2
  data_long <- data.frame()
  
  # Order generations
  gen_levels <- unique(data$Generation)
  data$Generation <- factor(data$Generation, levels = gen_levels)
  
  # Long format
  for (i in seq_len(nrow(data))) {
    genotype <- data$Genotype[i]
    id <- data$ID[i]
    gen <- data$Generation[i]
    gen_index <- which(levels(data$Generation) == gen)
    
    for (j in seq_along(gene_names)) {
      alleles <- strsplit(substr(genotype, (2 * j - 1), 2 * j), "")[[1]]
      gene <- gene_names[j]
      
      y_base <- (length(gen_levels) - gen_index) * (length(gene_names) + 1)
      y_pos <- y_base + (length(gene_names) - j + 1)
      
      data_long <- rbind(data_long,
                         data.frame(ID = id, Generation = gen, Gene = gene,
                                    Allele = alleles[1], x_offset = -offset, y = y_pos))
      data_long <- rbind(data_long,
                         data.frame(ID = id, Generation = gen, Gene = gene,
                                    Allele = alleles[2], x_offset = offset, y = y_pos))
    }
  }
  
  # Assign x positions
  data$GenIndex <- as.numeric(data$Generation)
  data <- data[order(data$GenIndex, data$ID), ]
  # Scale individuals horizontally to avoid crowding
  max_individuals <- max(table(data$Generation))
  scaling_factor <- 2  # spacing multiplier
  
  data$x_base <- ave(seq_along(data$ID), data$Generation, FUN = function(x) {
    n <- length(x)
    start <- (max_individuals - n) / 2 + 1  # center individuals
    return(seq(from = start, by = scaling_factor, length.out = n))
  })
  
  data_long <- merge(data_long, data[, c("ID", "Generation", "x_base")], by = c("ID", "Generation"))
  data_long$x <- data_long$x_base + data_long$x_offset
  
  # Palette
  alleles <- sort(unique(data_long$Allele))
  palette <- setNames(brewer.pal(n = min(10, length(alleles)), name = "Set3"), alleles)
  
  # Axis labels
  y_breaks <- c()
  y_labels <- c()
  for (i in seq_along(gen_levels)) {
    for (j in seq_along(gene_names)) {
      y_breaks <- c(y_breaks, (length(gen_levels) - i) * (length(gene_names) + 1) + (length(gene_names) - j + 1))
      y_labels <- c(y_labels, paste0(gene_names[j], " (", gen_levels[i], ")"))
    }
  }
  
  # Plot
  p <- ggplot(data_long, aes(x = x, y = y, fill = Allele)) +
    geom_tile(width = 0.35, height = 0.9, color = "white", linewidth = 0.3) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      expand = expansion(add = 0.5)
    ) +
    scale_fill_manual(values = palette) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_text(angle = 45, hjust = 1, size = 9)
    ) +
    labs(fill = "Allele")
  
  # Labels below
  if (show_labels) {
    label_df <- unique(data_long[, c("ID", "x", "Generation")])
    label_df <- label_df[!duplicated(label_df[c("ID", "Generation")]), ]
    label_df <- label_df[seq(1, nrow(label_df), by = 2), ]  # Show every 2nd label
    
    p <- p + geom_text(data = label_df, aes(x = x, y = 0, label = ID),
                       angle = 90, vjust = 0.5, size = 2, inherit.aes = FALSE)
  }
  
  # Draw parent-child lines if applicable
  if (all(c("Parent1", "Parent2") %in% colnames(data))) {
    id_coords <- unique(data_long[, c("ID", "Generation", "x")])
    id_coords$GenIndex <- as.numeric(id_coords$Generation)
    
    offspring <- data[, c("ID", "Generation", "Parent1", "Parent2")]
    offspring$GenIndex <- as.numeric(offspring$Generation)
    offspring$x_offspring <- data[data$ID %in% offspring$ID, "x_base"]
    
    parent_coords <- id_coords
    parent_coords$GenIndex <- parent_coords$GenIndex + 1
    
    segs1 <- merge(offspring, parent_coords, by.x = c("Parent1", "GenIndex"), by.y = c("ID", "GenIndex"))
    segs1 <- segs1[, c("x", "x_offspring", "GenIndex")]
    segs1$y <- (length(gen_levels) - segs1$GenIndex + 1) * (length(gene_names) + 1) + 0.5
    segs1$yend <- segs1$y - (length(gene_names) + 1)
    
    segs2 <- merge(offspring, parent_coords, by.x = c("Parent2", "GenIndex"), by.y = c("ID", "GenIndex"))
    segs2 <- segs2[, c("x", "x_offspring", "GenIndex")]
    segs2$y <- (length(gen_levels) - segs2$GenIndex + 1) * (length(gene_names) + 1) + 0.5
    segs2$yend <- segs2$y - (length(gene_names) + 1)
    
    segments <- rbind(
      data.frame(x = segs1$x, xend = segs1$x_offspring, y = segs1$y, yend = segs1$yend),
      data.frame(x = segs2$x, xend = segs2$x_offspring, y = segs2$y, yend = segs2$yend)
    )
    
    p <- p + geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend),
                          color = "gray40", linewidth = 0.4, inherit.aes = FALSE)
  }
  
  return(p)
}

create_cousin_pairs <- function(g2, g1) {
  g2$ParentPair <- paste(g2$MotherID, g2$FatherID, sep = "-")
  family_groups <- split(g2, g2$ParentPair)
  
  cousin_pairs <- list()
  family_keys <- names(family_groups)
  for (i in seq_along(family_keys)) {
    for (j in (i+1):length(family_keys)) {
      group1 <- family_groups[[i]]
      group2 <- family_groups[[j]]
      for (ind1 in 1:nrow(group1)) {
        for (ind2 in 1:nrow(group2)) {
          if (group1$Sex[ind1] != group2$Sex[ind2]) {
            if (group1$Sex[ind1] == "M") {
              cousin_pairs[[length(cousin_pairs) + 1]] <- list(
                father = group1$Genotype[ind1],
                mother = group2$Genotype[ind2]
              )
            } else {
              cousin_pairs[[length(cousin_pairs) + 1]] <- list(
                father = group2$Genotype[ind2],
                mother = group1$Genotype[ind1]
              )
            }
          }
        }
      }
    }
  }
  cousin_pairs
}

simulate_population <- function(generation_sizes, gene_names,
                                male_genotypes, male_props,
                                female_genotypes, female_props,
                                mating_type = "random",
                                selfing_rate = 0, seed = 123) {
  set.seed(seed)
  all_generations <- list()
  n_generations <- length(generation_sizes)
  
  create_initial_generation <- function(n_individuals) {
    males <- sample(male_genotypes, size = n_individuals / 2, replace = TRUE, prob = male_props)
    females <- sample(female_genotypes, size = n_individuals / 2, replace = TRUE, prob = female_props)
    data.frame(ID = 1:n_individuals,
               Genotype = c(males, females),
               Sex = c(rep("M", n_individuals / 2), rep("F", n_individuals / 2)),
               Generation = 1,
               stringsAsFactors = FALSE)
  }
  
  all_generations[[1]] <- create_initial_generation(generation_sizes[1])
  last_id <- generation_sizes[1]
  
  for (gen in 2:n_generations) {
    n_individuals <- generation_sizes[gen]
    parents <- all_generations[[gen - 1]]
    males <- parents[parents$Sex == "M", ]
    females <- parents[parents$Sex == "F", ]
    
    if (nrow(males) == 0 || nrow(females) == 0) {
      warning(paste0("Generation ", gen - 1, " has no males or females. Switching to random mating."))
      mating_type <- "random"
    }
    
    if (mating_type == "random") {
      fathers <- sample(males$Genotype, size = n_individuals / 2, replace = TRUE)
      mothers <- sample(females$Genotype, size = n_individuals / 2, replace = TRUE)
    } else if (mating_type == "selfing") {
      fathers <- character(n_individuals / 2)
      mothers <- character(n_individuals / 2)
      for (i in 1:(n_individuals / 2)) {
        if (runif(1) < selfing_rate) {
          parent <- parents[sample(nrow(parents), 1), ]
          fathers[i] <- parent$Genotype
          mothers[i] <- parent$Genotype
        } else {
          father <- males[sample(nrow(males), 1), ]
          mother <- females[sample(nrow(females), 1), ]
          fathers[i] <- father$Genotype
          mothers[i] <- mother$Genotype
        }
      }
    } else if (mating_type == "cousin") {
      if (gen < 3) {
        fathers <- sample(males$Genotype, size = n_individuals / 2, replace = TRUE)
        mothers <- sample(females$Genotype, size = n_individuals / 2, replace = TRUE)
      } else {
        g2 <- all_generations[[gen - 1]]
        g1 <- all_generations[[gen - 2]]
        
        cousin_pairs <- create_cousin_pairs(g2, g1)
        num_pairs_needed <- n_individuals / 2
        
        if (length(cousin_pairs) >= num_pairs_needed) {
          selected_pairs <- cousin_pairs[1:num_pairs_needed]
        } else {
          selected_pairs <- cousin_pairs
          extra_needed <- num_pairs_needed - length(cousin_pairs)
          extra_fathers <- sample(males$Genotype, size = extra_needed, replace = TRUE)
          extra_mothers <- sample(females$Genotype, size = extra_needed, replace = TRUE)
          selected_pairs <- c(selected_pairs, mapply(function(f, m) list(father = f, mother = m),
                                                     extra_fathers, extra_mothers, SIMPLIFY = FALSE))
        }
        
        fathers <- sapply(selected_pairs, function(p) p$father)
        mothers <- sapply(selected_pairs, function(p) p$mother)
      }
    } else {
      fathers <- sample(males$Genotype, size = n_individuals / 2, replace = TRUE)
      mothers <- sample(females$Genotype, size = n_individuals / 2, replace = TRUE)
    }
    
    offspring <- data.frame(ID = (last_id + 1):(last_id + n_individuals),
                            Genotype = mapply(generate_offspring_genotype, mothers, fathers,
                                              MoreArgs = list(gene_names = gene_names)),
                            Sex = sample(c("M", "F"), size = n_individuals, replace = TRUE),
                            Generation = gen,
                            stringsAsFactors = FALSE)
    
    all_generations[[gen]] <- offspring
    last_id <- last_id + n_individuals
  }
  
  do.call(rbind, all_generations)
}

# Define UI
ui <- fluidPage(
  titlePanel("Multilocus Mendelian Inheritance Simulation"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n_genes", "Number of genes:", value = 2, min = 1, max = 10),
      textInput("gene_names", "Gene names (comma-separated):", value = "A,B"),
      
      textInput("male_genotypes", "Male genotypes (comma-separated):", value = "AABB,AaBB,AaBb"),
      textInput("male_props", "Male genotype proportions (comma-separated):", value = "0.3,0.4,0.3"),
      
      textInput("female_genotypes", "Female genotypes (comma-separated):", value = "aaBB,aaBb,aabb"),
      textInput("female_props", "Female genotype proportions (comma-separated):", value = "0.2,0.5,0.3"),
      
      textInput("generation_sizes", "Generation sizes (comma-separated):", value = "10,20,40"),
      
      selectInput("mating_type", "Mating type:",
                  choices = c("random", "selfing", "cousin"), selected = "random"),
      
      conditionalPanel(
        condition = "input.mating_type == 'selfing'",
        sliderInput("selfing_rate", "Selfing rate:", min = 0, max = 1, value = 0.5)
      ),
      uiOutput("generationSelect"),
      
      actionButton("simulate", "Run Simulation")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Genotypes", 
                 div(style = "overflow-x: scroll; overflow-y: scroll; height: 700px;", 
                     plotOutput("genotypePlot", height = "1000px"))),
        tabPanel("Allele Frequency", plotOutput("allelePlot")),
        tabPanel("H_obs, H_exp, F", plotOutput("combinedPlot")),
        #tabPanel("Heterozygosity", plotOutput("heteroPlot")),
        #tabPanel("Inbreeding Coefficient (F)", plotOutput("fPlot")),
        tabPanel("Phenotype Summary", dataTableOutput("phenotypeTable")),
        tabPanel("Multilocus Heterozygosity", plotOutput("mlhPlot")),
        tabPanel("IBS Matrix", plotOutput("ibsPlot", height = "600px"))
        
      )
    )
    
  )
)

compute_population_stats <- function(data, gene_names) {
  stats <- data.frame()
  
  generations <- unique(data$Generation)
  for (gen in generations) {
    gen_data <- subset(data, Generation == gen)
    n <- nrow(gen_data)
    
    for (g in seq_along(gene_names)) {
      alleles <- substr(gen_data$Genotype, 2 * g - 1, 2 * g)
      allele_counts <- table(unlist(strsplit(paste(alleles, collapse = ""), "")))
      total_alleles <- sum(allele_counts)
      freqs <- allele_counts / total_alleles
      
      # Get heterozygosity
      heterozygotes <- sum(sapply(alleles, function(x) {
        a <- strsplit(x, "")[[1]]
        return(a[1] != a[2])
      }))
      H_obs <- heterozygotes / n
      
      if (length(freqs) == 2) {
        p <- freqs[1]
        q <- freqs[2]
        H_exp <- 2 * p * q
        F <- 1 - H_obs / H_exp
      } else {
        H_exp <- NA
        F <- NA
      }
      
      stats <- rbind(stats, data.frame(
        Generation = gen,
        Gene = gene_names[g],
        p = freqs[1],
        q = ifelse(length(freqs) > 1, freqs[2], 0),
        H_obs = H_obs,
        H_exp = H_exp,
        F = F
      ))
    }
  }
  return(stats)
}

# Define Server
server <- function(input, output, session) {
  data <- eventReactive(input$simulate, {
    gene_names <- unlist(strsplit(input$gene_names, ","))
    male_genotypes <- unlist(strsplit(gsub("\\s+", "", input$male_genotypes), ","))
    male_props <- as.numeric(unlist(strsplit(input$male_props, ",")))
    female_genotypes <- unlist(strsplit(gsub("\\s+", "", input$female_genotypes), ","))
    female_props <- as.numeric(unlist(strsplit(input$female_props, ",")))
    generation_sizes <- as.numeric(unlist(strsplit(input$generation_sizes, ",")))
    
    simulate_population(generation_sizes = generation_sizes,
                        gene_names = gene_names,
                        male_genotypes = male_genotypes,
                        male_props = male_props,
                        female_genotypes = female_genotypes,
                        female_props = female_props,
                        mating_type = input$mating_type,
                        selfing_rate = ifelse(input$mating_type == "selfing", input$selfing_rate, 0))
  })
  
  output$genotypePlot <- renderPlot({
    req(data())
    gene_names <- unlist(strsplit(input$gene_names, ","))
    genotype_grid_plot(data(), gene_names)
  })
  
  stats <- reactive({
    req(data())
    gene_names <- unlist(strsplit(input$gene_names, ","))
    compute_population_stats(data(), gene_names)
  })
  
  output$allelePlot <- renderPlot({
    df <- stats()
    ggplot(df, aes(x = Generation, y = p, color = Gene)) +
      ylim(0,1) +
      geom_line() + geom_point() + labs(title = "Allele Frequency (p)")
  })
  
  # output$heteroPlot <- renderPlot({
  # df <- stats()
  #ggplot(df, aes(x = Generation, y = H_obs, color = Gene)) +
  #  ylim(0,1) +
  # geom_line() + geom_point() + labs(title = "Observed Heterozygosity")
  #})
  
  #output$fPlot <- renderPlot({
  #  df <- stats()
  # ggplot(df, aes(x = Generation, y = F, color = Gene)) +
  #  geom_line() + geom_point() + labs(title = "Inbreeding Coefficient (F)")
  #})
  
  output$combinedPlot <- renderPlot({
    df <- stats()
    # Gather H_obs, H_exp, and F into long format
    df_long <- reshape2::melt(df, id.vars = c("Generation", "Gene"), 
                              measure.vars = c("H_obs", "H_exp", "F"),
                              variable.name = "Metric", value.name = "Value")
    
    ggplot(df_long, aes(x = Generation, y = Value, color = Gene, group = Gene)) +
      geom_line() + geom_point() +
      facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
      labs(title = "Heterozygosity and Inbreeding Coefficient Over Time",
           y = "", x = "Generation",
           subtitle = expression(F == 1 - frac(H[obs], H[exp]))) +
      theme_minimal(base_size = 14)
  })
  
  phenotype_summary <- reactive({
    req(data())
    gene_names <- unlist(strsplit(input$gene_names, ","))
    summarize_phenotypes(data(), gene_names)
  })
  
  output$phenotypeTable <- renderDataTable({
    phenotype_summary()
  }, options = list(pageLength = 10, scrollX = TRUE))
  
  mlh_stats <- reactive({
    req(data())
    gene_names <- unlist(strsplit(input$gene_names, ","))
    compute_multilocus_heterozygosity(data(), gene_names)
  })
  
  output$mlhPlot <- renderPlot({
    df <- mlh_stats()
    ggplot(df, aes(x = Generation, y = Mean_MLH)) +
      geom_line() + geom_point() +
      ylim(0, 1) +
      labs(title = "Mean Multilocus Heterozygosity", y = "Mean MLH", x = "Generation") +
      theme_minimal()
  })
  
  output$ibsPlot <- renderPlot({
    req(data(), input$selectedGen)
    gene_names <- unlist(strsplit(input$gene_names, ","))
    ibs <- compute_ibs_matrix(data(), gene_names, generation = input$selectedGen)

    df <- as.data.frame(as.table(ibs))
    colnames(df) <- c("ID1", "ID2", "IBS")
    
    ggplot(df, aes(x = ID1, y = ID2, fill = IBS)) +
      geom_tile() +
      scale_fill_viridis_c(limits = c(0, 1)) +
      labs(title = paste("IBS Similarity (Generation", input$selectedGen, ")"),
           fill = "IBS", x = "", y = "",
           subtitle = expression(IBS[ij] == frac("shared alleles", 2 * L))) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 90, size = 7),
            axis.text.y = element_text(size = 7))
  })
  
  output$generationSelect <- renderUI({
    req(data())
    selectInput("selectedGen", "Select Generation for IBS", 
                choices = sort(unique(data()$Generation)),
                selected = max(data()$Generation))
  })
  
  
}

shinyApp(ui = ui, server = server)
