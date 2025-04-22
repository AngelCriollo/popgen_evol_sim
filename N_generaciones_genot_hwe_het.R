################################################################################
#
#         SIMULADOR: EQUILIBRIO DE HW A LO LARGO DE VARIAS GENERACIONES
#.        Angel Alexandro Criollo Rayo
#.        Grupo de Citogenetica, Filogenia y Evolucion de Poblaciones
#         Doctorado en Ciencias Biomedicas -Universidad del Tolima-
#         Asignatura: Genetica Poblacional y Evolutiva (electiva)
#
################################################################################
# Simulador para observar el efecto de la deriva y el tamano de la poblacional
# en las frecuencias alelicas, genotipicas, revisar el equilibrio HW y las 
# heterocigosidades en cada generacion.
#-------------------------------------------------------------------------------

library(data.table)
library(ggplot2)

################################################################################
# 1. Parametros iniciales -EDITABLES-
n_populations <- 3
n_individuals <- 100
n_snps <- 5
nucleotides <- c("A", "C", "T", "G")
n_generations <- 10

################################################################################

# 2. Alelos escogidos por cada SNP
set.seed(123)
snp_alleles <- data.table(
  SNP = paste0("SNP", 1:n_snps),
  Allele1 = sample(nucleotides, n_snps, replace = TRUE)
)
snp_alleles[, Allele2 := sapply(Allele1, function(x) sample(setdiff(nucleotides, x), 1))]

# 3. Poblaciones iniciales
initialize_population <- function(pop_name, n_ind) {
  pop <- data.table(
    Population = pop_name,
    Individual = paste0(pop_name, "_", 1:n_ind),
    Generation = 0
  )
  
  for (i in 1:n_snps) {
    a1 <- snp_alleles[i, Allele1]
    a2 <- snp_alleles[i, Allele2]
    p <- runif(1, 0.2, 0.8) # Random starting frequency
    
    pop[, paste0("SNP", i) := sample(
      c(paste0(a1, a1), paste0(a1, a2), paste0(a2, a2)),
      size = n_ind,
      replace = TRUE,
      prob = c(p^2, 2*p*(1-p), (1-p)^2)
    )]
  }
  return(pop)
}

# 4. Creacion poblacion inicial
pop_list <- lapply(paste0("Pop", 1:n_populations), initialize_population, n_ind = n_individuals)
genotypes_df <- rbindlist(pop_list)

# 5. Funcion para crear la siguiente generacion
generate_next_generation <- function(current_gen_df, gen_num) {
  # Melt to get all alleles
  allele_pool <- melt(current_gen_df,
                      id.vars = c("Population", "Individual", "Generation"),
                      variable.name = "SNP",
                      value.name = "Genotype")[,
                                               {
                                                 alleles <- c(substr(Genotype, 1, 1), substr(Genotype, 2, 2))
                                                 .(Allele = alleles)
                                               }, by = .(Population, SNP)
                      ]
  
  # Crear nuevos genotipos a partir de muestreo aleatorio
  new_genotypes <- allele_pool[,
                               {
                                 sampled_alleles <- sample(Allele, 2 * n_individuals, replace = TRUE)
                                 genotypes <- paste0(
                                   sampled_alleles[seq(1, 2*n_individuals, by = 2)],
                                   sampled_alleles[seq(2, 2*n_individuals, by = 2)]
                                 )
                                 .(Individual = paste0("Gen", gen_num, "_", 1:n_individuals),
                                   Generation = gen_num,
                                   Genotype = genotypes)
                               }, by = .(Population, SNP)
  ]
  
  # Convert back to wide format
  dcast(new_genotypes, Population + Individual + Generation ~ SNP, value.var = "Genotype")
}

# 6. Run simulation
evolution_df <- copy(genotypes_df)
for (gen in 1:n_generations) {
  new_gen <- generate_next_generation(evolution_df[Generation == gen-1], gen)
  evolution_df <- rbindlist(list(evolution_df, new_gen))
}

# 7. Calculate allele frequencies over time (CORRECTED)
melted_evolution <- melt(evolution_df,
                         id.vars = c("Generation", "Population", "Individual"),
                         variable.name = "SNP",
                         value.name = "Genotype")

allele_freqs_over_time <- melted_evolution[,
                                           {
                                             alleles <- snp_alleles[SNP == .BY$SNP]
                                             a1_count <- sum(substr(Genotype, 1, 1) == alleles$Allele1) + 
                                               sum(substr(Genotype, 2, 2) == alleles$Allele1)
                                             total <- 2 * .N
                                             .(Freq1 = a1_count/total, Freq2 = (total - a1_count)/total)
                                           }, by = .(Generation, Population, SNP)
]

# 8. Calculate heterozygosity over time (using melted data)
het_over_time <- melted_evolution[,
                                  {
                                    alleles <- snp_alleles[SNP == .BY$SNP]
                                    # Count heterozygotes (alleles differ)
                                    het_count <- sum(substr(Genotype, 1, 1) != substr(Genotype, 2, 2))
                                    .(Het = het_count/.N, 
                                      N = .N) # Include sample size for reference
                                  }, 
                                  by = .(Generation, Population, SNP)
]

# 9. Function to calculate HWE tests for a generation
calculate_hwe <- function(gen_data, gen_num) {
  melted_gen <- melt(gen_data,
                     id.vars = c("Population", "Individual", "Generation"),
                     variable.name = "SNP",
                     value.name = "Genotype")
  
  hwe_results <- melted_gen[,
                            {
                              alleles <- snp_alleles[SNP == .BY$SNP]
                              
                              # Count observed genotypes
                              obs_aa <- sum(Genotype == paste0(alleles$Allele1, alleles$Allele1))
                              obs_ab <- sum(Genotype == paste0(alleles$Allele1, alleles$Allele2))
                              obs_bb <- sum(Genotype == paste0(alleles$Allele2, alleles$Allele2))
                              total <- obs_aa + obs_ab + obs_bb
                              
                              # Calculate allele frequencies
                              p <- (2 * obs_aa + obs_ab) / (2 * total)
                              q <- 1 - p
                              
                              # Calculate expected counts under HWE
                              exp_aa <- p^2 * total
                              exp_ab <- 2 * p * q * total
                              exp_bb <- q^2 * total
                              
                              # Calculate chi-square (with continuity correction)
                              chi_sq <- sum(
                                (abs(obs_aa - exp_aa) - 0.5)^2 / exp_aa,
                                (abs(obs_ab - exp_ab) - 0.5)^2 / exp_ab,
                                (abs(obs_bb - exp_bb) - 0.5)^2 / exp_bb
                              )
                              
                              # Calculate p-value
                              p_val <- pchisq(chi_sq, df = 1, lower.tail = FALSE)
                              
                              # Return results
                              .(Obs_A1A1 = obs_aa,
                                Obs_A1A2 = obs_ab,
                                Obs_A2A2 = obs_bb,
                                Exp_A1A1 = exp_aa,
                                Exp_A1A2 = exp_ab,
                                Exp_A2A2 = exp_bb,
                                ChiSq = chi_sq,
                                Pval = p_val,
                                N = total)
                            },
                            by = .(Population, SNP)
  ]
  
  cbind(Generation = gen_num, hwe_results)
}

# Calculate HWE for all generations
all_hwe_results <- rbindlist(
  lapply(0:n_generations, function(gen) {
    calculate_hwe(evolution_df[Generation == gen], gen)
  })
)

# Format results for better readability
formatted_hwe <- all_hwe_results[, .(
  Generation,
  Population,
  SNP,
  `Obs A1A1` = Obs_A1A1,
  `Obs A1A2` = Obs_A1A2,
  `Obs A2A2` = Obs_A2A2,
  `Exp A1A1` = round(Exp_A1A1, 2),
  `Exp A1A2` = round(Exp_A1A2, 2),
  `Exp A2A2` = round(Exp_A2A2, 2),
  `Chi-Sq` = round(ChiSq, 4),
  `P-value` = ifelse(Pval < 0.0001, "<0.0001", format(round(Pval, 4), nsmall = 4)),
  `Significant` = ifelse(Pval < 0.05, "*", "")
)]

# View formatted results
print(formatted_hwe)

all_hwe_results <- rbindlist(
  lapply(0:n_generations, function(gen) {
    calculate_hwe(evolution_df[Generation == gen], gen)
  })
)


################################################################################
#
# Analisis visual
#
################################################################################

# 8. Plot results allele freq.
ggplot(allele_freqs_over_time, 
       aes(x = Generation, y = Freq1, color = Population, group = interaction(Population, SNP))) +
  geom_line() +
  facet_wrap(~ SNP) +
  labs(title = "Allele Frequency Change Over Generations",
       y = "Frequency of Allele 1") +
  theme_minimal()

# View heterozygosity results
print(het_over_time)

# Plot heterozygosity changes
ggplot(het_over_time, 
       aes(x = Generation, y = Het, color = Population, group = Population)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ SNP, ncol = 2) +
  labs(title = "Heterozygosity Change Over Generations",
       subtitle = "Expected decrease due to genetic drift",
       y = "Observed Heterozygosity",
       x = "Generation") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom")

# Combined plot with allele frequencies and heterozygosity
combined_results <- merge(allele_freqs_over_time, het_over_time, 
                          by = c("Generation", "Population", "SNP"))

ggplot(combined_results, aes(x = Generation)) +
  geom_line(aes(y = Freq1, color = Population, linetype = "Allele Freq")) +
  geom_line(aes(y = Het, color = Population, linetype = "Heterozygosity")) +
  facet_grid(Population ~ SNP) +
  labs(title = "Genetic Diversity Dynamics",
       y = "Value",
       color = "Population",
       linetype = "Metric") +
  scale_linetype_manual(values = c("Allele Freq" = 1, "Heterozygosity" = 2)) +
  theme_minimal()

# HWE p-val dynamics
hwe_plot <- ggplot(all_hwe_results, aes(x = Generation, y = -log10(Pval))) +
  geom_point(aes(color = Population), size = 1, alpha = 0.8) +
  geom_line(aes(group = interaction(Population, SNP), color = Population), alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  facet_wrap(~ SNP, ncol = 2) +
  labs(title = "Hardy-Weinberg Equilibrium Test Results",
       subtitle = "Dashed line indicates p = 0.05 significance threshold",
       y = "-log10(p-value)",
       x = "Generation") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom")

print(hwe_plot)




















