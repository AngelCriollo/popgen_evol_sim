################################################################################
#
#               APP SIMULADOR MENDELIANO: DISTRIBUCION INDEPENDIENTE
#               Autor: Angel Criollo Rayo 
#               Version 2: 01 - Mayo - 2025, asistido por ChatGPT
#               Grupo: Citogenetica, Filogenia y Evlucion de Poblaciones
#               Doctorado en Ciencias Biomedicas
#               Asignatura: Genetica Poblacional y Evolutiva (electiva)
#
################################################################################
#
# Este script esta diseñado para simular la distribución independiente
# Simula uniones o apareamientos por autofecundacion en distintos grados o
# fecundación cruzada.
# Se puede modular: la frecuencia inicial de los genotipos, grados de "selfing"
# numero de individuos/generacion, numero de generaciones, numero de genes.
################################################################################

library(shiny)
library(dplyr)
library(digest)
library(ggplot2)
library(tidyr)

# Define full gene pool
default_gene_pool <- LETTERS[1:10]

allele_colors <- setNames(lapply(default_gene_pool, function(g) {
  c(setNames(paste0("#", substr(digest::digest(g), 1, 6)), g),
    setNames(paste0("#", substr(digest::digest(tolower(g)), 1, 6)), tolower(g)))
}), default_gene_pool)

parse_genotype_string <- function(genostr, genes) {
  alleles <- strsplit(gsub("\\s", "", genostr), "")[[1]]
  if (length(alleles) != length(genes) * 2) return(NULL)
  setNames(lapply(seq(1, length(alleles), by = 2), function(i) alleles[i:(i+1)]), genes)
}

phenotype <- function(geno) {
  sapply(names(geno), function(g) {
    if (any(grepl("[A-Z]", geno[[g]]))) paste0(g, "_") else paste0(tolower(g), tolower(g))
  }) %>% paste0(collapse = "")
}

sample_parents <- function(geno_strs, props, genes, sex, count, id_start) {
  genos <- strsplit(geno_strs, ",")[[1]]
  genos <- trimws(genos)
  props <- as.numeric(strsplit(props, ",")[[1]])
  
  if (length(genos) != length(props) || sum(props) <= 0 || any(is.na(props))) return(NULL)
  props <- props / sum(props)
  
  sampled <- sample(genos, size = count, replace = TRUE, prob = props)
  ids <- seq(id_start, id_start + count - 1)
  lapply(seq_along(sampled), function(i) {
    geno <- parse_genotype_string(sampled[i], genes)
    list(
      id = ids[i],
      sex = sex,
      generation = 0,
      genotype = geno,
      phenotype = phenotype(geno),
      parents = NA
    )
  })
}

make_generation <- function(prev_gen, gen_num, selfing_ratio, n, id_start, genes) {
  offspring <- list()
  for (i in 1:n) {
    sex <- sample(c("M", "F"), 1)
    is_self <- runif(1) < selfing_ratio || length(prev_gen) == 0
    
    if (is_self) {
      parent <- sample(prev_gen, 1)[[1]]
      geno <- lapply(parent$genotype, function(a) sample(a, 2, replace = TRUE))
      parents <- rep(parent$id, 2)
    } else {
      mom <- sample(Filter(function(p) p$sex == "F", prev_gen), 1)[[1]]
      dad <- sample(Filter(function(p) p$sex == "M", prev_gen), 1)[[1]]
      geno <- mapply(
        function(m, d) c(sample(m, 1), sample(d, 1)),
        mom$genotype, dad$genotype, SIMPLIFY = FALSE
      )
      parents <- c(mom$id, dad$id)
    }
    
    offspring[[i]] <- list(
      id = id_start + i,
      sex = sex,
      generation = gen_num,
      genotype = geno,
      phenotype = phenotype(geno),
      parents = parents
    )
  }
  offspring
}

ui <- fluidPage(
  titlePanel("Mendelian Simulation (Configurable Genotypes, Validation, ggplot)"),
  sidebarLayout(
    sidebarPanel(
      numericInput("n_gen", "Generations (excluding parents):", 3, min = 1),
      numericInput("n_ind", "Individuals per generation:", 20, min = 2),
      sliderInput("selfing", "Selfing ratio:", 0, 1, 0, step = 0.1),
      numericInput("n_genes", "Number of genes:", 5, min = 1, max = 10),
      textInput("male_genos", "Male genotypes (comma-separated):", "AaBbCcDdEe"),
      textInput("male_props", "Male proportions:", "1.0"),
      textInput("female_genos", "Female genotypes (comma-separated):", "AaBbCcDdEe"),
      textInput("female_props", "Female proportions:", "1.0"),
      actionButton("simulate", "Run Simulation"),
      textOutput("validation")
    ),
    mainPanel(
      plotOutput("plot", height = "800px"),
      verbatimTextOutput("summary")
    )
  )
)

server <- function(input, output) {
  output$validation <- renderText({
    genes <- default_gene_pool[1:input$n_genes]
    m <- strsplit(input$male_genos, ",")[[1]]
    f <- strsplit(input$female_genos, ",")[[1]]
    len_check <- all(nchar(trimws(c(m, f))) == input$n_genes * 2)
    if (!len_check) return("❌ Each genotype must have 2 alleles per gene.")
    
    mp <- as.numeric(strsplit(input$male_props, ",")[[1]])
    fp <- as.numeric(strsplit(input$female_props, ",")[[1]])
    if (length(mp) != length(m) || length(fp) != length(f)) return("❌ Mismatched number of genotypes and proportions.")
    if (abs(sum(mp) - 1) > 1e-6 || abs(sum(fp) - 1) > 1e-6) return("❌ Proportions must sum to 1.")
    "✅ Input valid."
  })
  
  sim <- eventReactive(input$simulate, {
    genes <- default_gene_pool[1:input$n_genes]
    id_counter <- 1
    
    males <- sample_parents(input$male_genos, input$male_props, genes, "M", 5, id_counter)
    id_counter <- id_counter + 5
    females <- sample_parents(input$female_genos, input$female_props, genes, "F", 5, id_counter)
    id_counter <- id_counter + 5
    
    if (is.null(males) || is.null(females)) return(NULL)
    
    parents <- c(males, females)
    all_inds <- parents
    prev <- parents
    
    for (g in 1:input$n_gen) {
      next_gen <- make_generation(prev, g, input$selfing, input$n_ind, id_counter, genes)
      all_inds <- c(all_inds, next_gen)
      prev <- next_gen
      id_counter <- id_counter + input$n_ind
    }
    list(individuals = all_inds, genes = genes)
  })
  
  output$summary <- renderPrint({
    res <- sim()
    if (is.null(res)) return("Invalid input.")
    inds <- res$individuals
    gens <- split(inds, sapply(inds, function(i) i$generation))
    for (g in names(gens)) {
      cat(paste("Generation", g, "-", length(gens[[g]]), "individuals\n"))
      tab <- table(sapply(gens[[g]], function(i) i$phenotype))
      print(tab)
      cat("\n")
    }
  })
  
  output$plot <- renderPlot({
    res <- sim()
    if (is.null(res)) return(NULL)
    inds <- res$individuals
    genes <- res$genes
    
    positions <- data.frame(
      id = sapply(inds, `[[`, "id"),
      generation = sapply(inds, `[[`, "generation"),
      sex = sapply(inds, `[[`, "sex")
    ) %>%
      group_by(generation) %>%
      mutate(ind_index = row_number()) %>%
      ungroup()
    
    df <- do.call(rbind, lapply(inds, function(ind) {
      data.frame(
        id = ind$id,
        generation = ind$generation,
        sex = ind$sex,
        gene = rep(genes, each = 2),
        allele = unlist(ind$genotype),
        allele_num = rep(1:2, times = length(genes)),
        parents = paste(ind$parents, collapse = ",")
      )
    })) %>%
      left_join(positions, by = c("id", "generation", "sex")) %>%
      mutate(
        x = ind_index,
        y = -generation * 3 + as.numeric(factor(gene)) * 0.4 + (allele_num - 1) * 0.1,
        shape = ifelse(sex == "M", 22, 21),
        fill = mapply(function(g, a) allele_colors[[g]][[a]], gene, allele)
      )
    
    parent_edges <- do.call(rbind, lapply(inds, function(ind) {
      if (!is.null(ind$parents) && all(!is.na(ind$parents))) {
        child_row <- positions[positions$id == ind$id, ]
        p1_row <- positions[positions$id == ind$parents[1], ]
        p2_row <- positions[positions$id == ind$parents[2], ]
        rbind(
          data.frame(x = p1_row$ind_index, xend = child_row$ind_index,
                     y = -p1_row$generation * 3, yend = -child_row$generation * 3),
          data.frame(x = p2_row$ind_index, xend = child_row$ind_index,
                     y = -p2_row$generation * 3, yend = -child_row$generation * 3)
        )
      }
    }))
    
    ggplot(df, aes(x = x, y = y)) +
      geom_segment(data = parent_edges, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                   color = "gray60", linetype = "dashed", linewidth = 0.4) +
      geom_point(aes(fill = fill, shape = sex), size = 4, stroke = 0.7, color = "black") +
      scale_shape_manual(values = c("M" = 22, "F" = 21)) +
      scale_fill_identity() +
      theme_minimal() +
      labs(title = "Mendelian Simulation: Genotype Visualization with Ancestry",
           x = "Individuals", y = "Generation (stacked genes)") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  })
}

shinyApp(ui = ui, server = server)

