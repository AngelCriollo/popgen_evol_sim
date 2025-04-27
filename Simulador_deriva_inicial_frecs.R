################################################################################
#
#               SIMULADOR DE DERIVA
#               Autor: Angel Criollo Rayo & DeepSeek
#               Version 7: 26 - Abril - 2025
#               Grupo: Citogenetica, Filogenia y Evlucion de Poblaciones
#               Doctorado en Ciencias Biomedicas
#               Asignatura: Genetica Poblacional y Evolutiva (electiva)
#
################################################################################
#
# Este script esta diseñado para analizar los efectos del tamaño poblacional
# en las frecuencias alelicas a lo largo de multiples generaciones. 
# en varios marcadores a la vez, sin tener en cuenta SN, mutación, ligamiento o endogamia
#
################################################################################


library(tidyverse)
library(gridExtra)

set.seed(123) # For reproducibility

# Parameters
n_populations <- 1000       # Number of populations
pop_size <- 500           # Size of each population
n_markers <- 3           # Number of SNP markers
n_generations <- 50      # Number of generations to simulate
initial_freqs <- rep(0.5,n_markers) # runif(n_markers) Random initial frequencies for each marker

# Function to generate random biallelic SNPs
generate_snps <- function(n) {
  nucleotides <- c("A", "C", "T", "G")
  snps <- list()
  for (i in 1:n) {
    alleles <- sample(nucleotides, 2, replace = FALSE)
    snps[[i]] <- list(allele1 = alleles[1], allele2 = alleles[2])
  }
  return(snps)
}

# Function to initialize genotypes for a population
initialize_genotypes <- function(pop_size, markers, initial_freqs) {
  genotypes <- list()
  for (m in 1:length(markers)) {
    freq <- initial_freqs[m]
    marker <- markers[[m]]
    alleles <- c(marker$allele1, marker$allele2)
    
    # Generate genotypes according to initial frequency
    n_allele1 <- round(pop_size * 2 * freq)
    n_allele2 <- pop_size * 2 - n_allele1
    
    # Ensure we have exactly the right frequencies
    all_alleles <- c(rep(alleles[1], n_allele1), rep(alleles[2], n_allele2))
    all_alleles <- sample(all_alleles) # Shuffle
    
    # Create genotypes
    pop_genotypes <- character(pop_size)
    for (i in 1:pop_size) {
      a1 <- all_alleles[2*i - 1]
      a2 <- all_alleles[2*i]
      pop_genotypes[i] <- paste(sort(c(a1, a2)), collapse = "")
    }
    
    genotypes[[m]] <- pop_genotypes
  }
  return(genotypes)
}

# Function to calculate allele frequencies
calculate_allele_freqs <- function(genotypes, markers) {
  freqs <- numeric(length(markers))
  for (m in 1:length(markers)) {
    marker <- markers[[m]]
    alleles <- c(marker$allele1, marker$allele2)
    
    # Count alleles
    geno_string <- paste(genotypes[[m]], collapse = "")
    count1 <- str_count(geno_string, alleles[1])
    count2 <- str_count(geno_string, alleles[2])
    
    freqs[m] <- count1 / (count1 + count2)
  }
  return(freqs)
}

# Function to calculate heterozygosity
calculate_heterozygosity <- function(genotypes, markers) {
  obs_het <- numeric(length(markers))
  exp_het <- numeric(length(markers))
  
  for (m in 1:length(markers)) {
    marker <- markers[[m]]
    geno <- genotypes[[m]]
    
    # Observed heterozygosity
    het_count <- sum(str_sub(geno, 1, 1) != str_sub(geno, 2, 2))
    obs_het[m] <- het_count / length(geno)
    
    # Expected heterozygosity under HWE
    alleles <- c(marker$allele1, marker$allele2)
    geno_string <- paste(geno, collapse = "")
    count1 <- str_count(geno_string, alleles[1])
    count2 <- str_count(geno_string, alleles[2])
    p <- count1 / (count1 + count2)
    q <- 1 - p
    exp_het[m] <- 2 * p * q
  }
  
  return(list(obs_het = obs_het, exp_het = exp_het))
}

# Function to perform chi-square test for HWE
test_hwe <- function(genotypes, markers) {
  p_values <- numeric(length(markers))
  
  for (m in 1:length(markers)) {
    marker <- markers[[m]]
    geno <- genotypes[[m]]
    
    # Count genotypes
    alleles <- c(marker$allele1, marker$allele2)
    hom1 <- sum(geno == paste(rep(alleles[1], 2), collapse = ""))
    hom2 <- sum(geno == paste(rep(alleles[2], 2), collapse = ""))
    het <- length(geno) - hom1 - hom2
    
    # Calculate expected counts under HWE
    total <- hom1 + hom2 + het
    p <- (2 * hom1 + het) / (2 * total)
    q <- 1 - p
    
    exp_hom1 <- p^2 * total
    exp_hom2 <- q^2 * total
    exp_het <- 2 * p * q * total
    
    # Chi-square test
    observed <- c(hom1, het, hom2)
    expected <- c(exp_hom1, exp_het, exp_hom2)
    
    # Handle cases where expected counts are too small
    if (any(expected < 5)) {
      p_values[m] <- NA
    } else {
      chi_sq <- sum((observed - expected)^2 / expected)
      p_values[m] <- pchisq(chi_sq, df = 1, lower.tail = FALSE)
    }
  }
  
  return(p_values)
}

# Function to simulate next generation
next_generation <- function(genotypes, markers, pop_size) {
  new_genotypes <- list()
  
  for (m in 1:length(markers)) {
    marker <- markers[[m]]
    geno <- genotypes[[m]]
    
    # Create allele pool
    alleles <- unlist(str_split(geno, ""))
    
    # Sample new genotypes
    new_geno <- character(pop_size)
    for (i in 1:pop_size) {
      a1 <- sample(alleles, 1)
      a2 <- sample(alleles, 1)
      new_geno[i] <- paste(sort(c(a1, a2)), collapse = "")
    }
    
    new_genotypes[[m]] <- new_geno
  }
  
  return(new_genotypes)
}

# Generate random SNPs
markers <- generate_snps(n_markers)

# Initialize populations
populations <- list()
for (p in 1:n_populations) {
  populations[[p]] <- list(
    genotypes = initialize_genotypes(pop_size, markers, initial_freqs),
    marker_info = markers
  )
}

# Data structures to store results
history <- list()
freq_history <- list()
het_history <- list()
hwe_history <- list()
genotype_counts_history <- list()

# Main simulation loop
for (gen in 0:n_generations) {
  for (p in 1:n_populations) {
    # Calculate statistics for current generation
    freqs <- calculate_allele_freqs(populations[[p]]$genotypes, markers)
    het <- calculate_heterozygosity(populations[[p]]$genotypes, markers)
    hwe_p <- test_hwe(populations[[p]]$genotypes, markers)
    
    # Store results
    freq_history[[length(freq_history) + 1]] <- tibble(
      generation = gen,
      population = p,
      marker = 1:n_markers,
      allele_freq = freqs
    )
    
    het_history[[length(het_history) + 1]] <- tibble(
      generation = gen,
      population = p,
      marker = 1:n_markers,
      obs_het = het$obs_het,
      exp_het = het$exp_het
    )
    
    hwe_history[[length(hwe_history) + 1]] <- tibble(
      generation = gen,
      population = p,
      marker = 1:n_markers,
      hwe_p_value = hwe_p
    )
    
    # Store genotype counts
    for (m in 1:n_markers) {
      marker <- markers[[m]]
      geno <- populations[[p]]$genotypes[[m]]
      alleles <- c(marker$allele1, marker$allele2)
      
      hom1 <- sum(geno == paste(rep(alleles[1], 2), collapse = ""))
      hom2 <- sum(geno == paste(rep(alleles[2], 2), collapse = ""))
      het <- length(geno) - hom1 - hom2
      
      genotype_counts_history[[length(genotype_counts_history) + 1]] <- tibble(
        generation = gen,
        population = p,
        marker = m,
        hom1 = hom1,
        het = het,
        hom2 = hom2
      )
    }
    
    # Store complete genotype history (only for first few generations to save memory)
    if (gen <= 5) {
      for (m in 1:n_markers) {
        for (i in 1:pop_size) {
          history[[length(history) + 1]] <- tibble(
            generation = gen,
            population = p,
            marker = m,
            individual = i,
            genotype = populations[[p]]$genotypes[[m]][i]
          )
        }
      }
    }
    
    # Evolve to next generation (except for last generation)
    if (gen < n_generations) {
      populations[[p]]$genotypes <- next_generation(
        populations[[p]]$genotypes,
        markers,
        pop_size
      )
    }
  }
}

# Combine results into data frames
freq_df <- bind_rows(freq_history)
het_df <- bind_rows(het_history) %>% 
  pivot_longer(cols = c(obs_het, exp_het), names_to = "het_type", values_to = "heterozygosity")
hwe_df <- bind_rows(hwe_history)
genotype_counts_df <- bind_rows(genotype_counts_history)
genotype_history_df <- if (length(history) > 0) bind_rows(history) else NULL

# Create plots
# 1. Allele frequency dynamics
p1 <- ggplot(freq_df, aes(x = generation, y = allele_freq, color = factor(population))) +
  geom_line() +
  facet_wrap(~marker, labeller = labeller(marker = function(x) paste("Marker", x))) +
  labs(title = paste0("Allele Frequency Dynamics"," - ","N=",pop_size), 
       y = "Allele Frequency", 
       color = "Population") +
  ylim(0,1) +
  theme_minimal()

# 2. Heterozygosity dynamics
p2 <- ggplot(het_df, aes(x = generation, y = heterozygosity, 
                         color = factor(population), linetype = het_type)) +
  geom_line() +
  facet_wrap(~marker, labeller = labeller(marker = function(x) paste("Marker", x))) +
  labs(title = paste0("Heterozygosity Dynamics"," - ","N=",pop_size), 
       y = "Heterozygosity", 
       color = "Population",
       linetype = "Type") +
  theme_minimal()

# 3. HWE test p-values
p3 <- ggplot(hwe_df, aes(x = generation, y = -log10(hwe_p_value), color = factor(population))) +
  geom_line() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  facet_wrap(~marker, labeller = labeller(marker = function(x) paste("Marker", x))) +
  labs(title = paste0("HWE Test p-values (-log10 scale)", " - ","N=",pop_size), 
       y = "-log10(p-value)", 
       color = "Population") +
  theme_minimal()

# Display plots
grid.arrange(p1, p2, p3, ncol = 1)

# Output data frames
list(
  genotype_history = genotype_history_df,
  genotype_counts = genotype_counts_df,
  allele_frequencies = freq_df,
  heterozygosity = het_df,
  hwe_tests = hwe_df
)

#===============================================================================

# histograma para mas de 100 muestras o poblaciones
# este paso es para analizar: si p=q=0.5, las frecuencias de la siguiente
# generación, tendran una distribucion normal, probar con n=5,50 y 500 (para 1000 muestras o poblaciones)
library(fitdistrplus)

pop_sd <- function( input_list )
{ n <- length( input_list)
sd( input_list )*sqrt((n-1)/n)
}

filtered_data <- freq_df %>%
  filter(generation ==1 & marker==1) # 

media_f <- mean(filtered_data$allele_freq)
sigma_f <- pop_sd(filtered_data$allele_freq)

breaks <- pretty(range(filtered_data$allele_freq), n = 20, min.n = 1) 
bwidth <- breaks[2]-breaks[1]

y <- dnorm(filtered_data$allele_freq, mean = media_f, sd = sigma_f) 

ggplot(filtered_data, aes(x = allele_freq)) + 
  geom_histogram(binwidth = bwidth, aes(y = ..density..),
                 fill="lightgrey", colour="black")  +
  stat_function(fun = dnorm, args = list(mean= media_f, sd = sigma_f)) + 
  xlim(0, 1) +
  labs(title = paste0("Allele Freq. Distribution"," n_pop=",n_populations,",N=",pop_size), 
       y = "Allele Frequency", 
       color = "Population") +
  geom_area(data = filtered_data, 
            aes(x = allele_freq,
                y = ifelse(allele_freq < media_f + 1 * sigma_f & allele_freq > media_f - 1 * sigma_f,
                           y, 0)),
            fill = "red", alpha = 0.2)

descdist(filtered_data$allele_freq, boot = 1000)



