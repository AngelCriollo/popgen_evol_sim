################################################################################
#
#               EQUILIBRIO DE HW EN DOS O MAS POBLACIONES
#               Angel Criollo Rayo
#               Version 6: 26 - Marzo - 2025
#               Grupo: Citogenetica, Filogenia y Evlucion de Poblaciones
#               Doctorado en Ciencias Biomedicas -UT
#               Asignatura: Genetica Poblacional y Evolutiva (electiva)
#
################################################################################
#
# Simulador de varias poblaciones con datos genotipicos individuales.  
# para diferentes marcadores. Se puede controlar el n/genotipo en las 3 pob.
# El producto final de este codigo es un set de datos simulado: genotypes_df
# a partir de estos datos, elaboren un script y calculen: el valor p del chi cua-
# drado para cada SNP y revisar si estan o no en EHW. Diseñar una forma grafica de visualizar.
################################################################################
library(data.table)
library(ggplot2)
library(gridExtra) # Create a visual table using gridExtra
library(grid)

### PARAMETROS ###
n_populations <- 3       # Number of populations
n_individuals <- 100     # Individuals per population
n_snps <- 5              # Number of genetic markers
nucleotides <- c("A", "C", "T", "G") # Possible alleles

### Define Initial Genotype Counts ###
# Structure: list(Population = list(SNP = c(AA, Aa, aa)))
initial_genotypes <- list(
  Pop1 = list(
    SNP1 = c(25, 50, 25),  # 25 AA, 50 Aa, 25 aa
    SNP2 = c(10, 80, 10),
    SNP3 = c(40, 20, 40),
    SNP4 = c(30, 40, 30),
    SNP5 = c(5, 90, 5)
  ),
  Pop2 = list(
    SNP1 = c(45, 10, 45),
    SNP2 = c(20, 60, 20),
    SNP3 = c(10, 80, 10),
    SNP4 = c(25, 50, 25),
    SNP5 = c(30, 40, 30)
  ),
  Pop3 = list(
    SNP1 = c(10, 80, 10),
    SNP2 = c(40, 20, 40),
    SNP3 = c(25, 50, 25),
    SNP4 = c(15, 70, 15),
    SNP5 = c(50, 0, 50)    # Fixed for one allele
  )
)
################################################################################
# Primera generacion -g1
################################################################################
### 1. Initialize SNP Alleles ###
set.seed(123)
snp_alleles <- data.table(
  SNP = paste0("SNP", 1:n_snps),
  Allele1 = sample(nucleotides, n_snps, replace = TRUE)
)
snp_alleles[, Allele2 := sapply(Allele1, function(x) sample(setdiff(nucleotides, x), 1))]

### 2. Create Populations with Defined Genotype Counts ###
create_population <- function(pop_name, genotype_counts) {
  pop <- data.table(
    Population = pop_name,
    Individual = paste0(pop_name, "_", 1:n_individuals)
  )
  
  for (snp in paste0("SNP", 1:n_snps)) {
    snp_num <- as.numeric(gsub("SNP", "", snp))
    alleles <- snp_alleles[SNP == snp]
    
    # Create genotypes according to specified counts
    aa <- paste0(alleles$Allele1, alleles$Allele1)
    ab <- paste0(alleles$Allele1, alleles$Allele2)
    bb <- paste0(alleles$Allele2, alleles$Allele2)
    
    counts <- genotype_counts[[snp_num]]
    genotypes <- c(
      rep(aa, counts[1]),
      rep(ab, counts[2]),
      rep(bb, counts[3])
    )
    
    # If counts don't sum to n_individuals, fill remaining randomly
    if (length(genotypes) < n_individuals) {
      remaining <- n_individuals - length(genotypes)
      genotypes <- c(genotypes, sample(c(aa, ab, bb), remaining, replace = TRUE))
    }
    
    pop[, (snp) := sample(genotypes)]  # Randomize order
  }
  return(pop)
}

### 3. Generate All Populations ###
genotypes_df <- rbindlist(lapply(1:n_populations, function(i) {
  create_population(paste0("Pop", i), initial_genotypes[[i]])
}))

### 4. Verify Initial Genotype Counts ###
verify_counts <- melt(genotypes_df, 
                      id.vars = c("Population", "Individual"),
                      variable.name = "SNP",
                      value.name = "Genotype")[,
                                               {
                                                 alleles <- snp_alleles[SNP == .BY$SNP]
                                                 aa <- sum(Genotype == paste0(alleles$Allele1, alleles$Allele1))
                                                 ab <- sum(Genotype == paste0(alleles$Allele1, alleles$Allele2))
                                                 bb <- sum(Genotype == paste0(alleles$Allele2, alleles$Allele2))
                                                 .(AA = aa, Aa = ab, aa = bb)
                                               }, by = .(Population, SNP)
                      ]

print(verify_counts)

################################################################################
# Grafico de resumen primera generacion -g1
################################################################################

### 5. Visualize Initial Genotype Frequencies ###
# First reshape the data for plotting
genotype_long <- melt(verify_counts, 
                      id.vars = c("Population", "SNP"),
                      variable.name = "Genotype",
                      value.name = "Count")

# Define consistent genotype colors
genotype_colors <- c("AA" = "#1b9e77",  # Green
                     "Aa" = "#d95f02",   # Orange
                     "aa" = "#7570b3")   # Purple

# Create the plot
ggplot(genotype_long, aes(x = SNP, y = Count, fill = Genotype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Population, nrow = 1) +
  scale_fill_manual(values = genotype_colors,
                    labels = c("Homozygous1 (AA)", 
                               "Heterozygous (Aa)", 
                               "Homozygous2 (aa)")) +
  labs(title = "Initial Genotype Counts by Population",
       subtitle = "Showing all three genotype categories for each SNP",
       y = "Number of Individuals",
       x = "SNP Markers",
       fill = "Genotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3)
################################################################################
#
# Segunda generacion -g2
#
################################################################################
# Function to generate next generation based on HWE expectations
generate_next_generation <- function(current_gen_df) {
  # Calculate allele frequencies for current generation
  allele_freqs <- melt(current_gen_df, 
                       id.vars = c("Population", "Individual"),
                       variable.name = "SNP",
                       value.name = "Genotype")[,
                                                {
                                                  alleles <- snp_alleles[SNP == .BY$SNP]
                                                  a1_count <- sum(substr(Genotype, 1, 1) == alleles$Allele1) + 
                                                    sum(substr(Genotype, 2, 2) == alleles$Allele1)
                                                  total <- 2 * .N
                                                  .(p = a1_count/total, q = 1 - a1_count/total)
                                                }, by = .(Population, SNP)
                       ]
  
  # Create new generation based on HWE expectations
  new_gen <- allele_freqs[,
                          {
                            # Calculate HWE genotype probabilities
                            probs <- c(p^2, 2*p*q, q^2)
                            
                            # Sample new genotypes
                            genotypes <- sample(
                              c(paste0(snp_alleles[SNP == .BY$SNP, Allele1], 
                                       snp_alleles[SNP == .BY$SNP, Allele1]),
                                paste0(snp_alleles[SNP == .BY$SNP, Allele1], 
                                       snp_alleles[SNP == .BY$SNP, Allele2]),
                                paste0(snp_alleles[SNP == .BY$SNP, Allele2], 
                                       snp_alleles[SNP == .BY$SNP, Allele2])),
                              size = n_individuals,
                              replace = TRUE,
                              prob = probs
                            )
                            
                            .(Individual = paste0("Ind_", 1:n_individuals),
                              Genotype = genotypes)
                          }, by = .(Population, SNP)
  ]
  
  # Convert back to wide format
  dcast(new_gen, Population + Individual ~ SNP, value.var = "Genotype")
}

# Generate next generation
next_gen_df <- generate_next_generation(genotypes_df)

# Verify genotype counts in new generation
new_gen_counts <- melt(next_gen_df, 
                       id.vars = c("Population", "Individual"),
                       variable.name = "SNP",
                       value.name = "Genotype")[,
                                                {
                                                  alleles <- snp_alleles[SNP == .BY$SNP]
                                                  aa <- sum(Genotype == paste0(alleles$Allele1, alleles$Allele1))
                                                  ab <- sum(Genotype == paste0(alleles$Allele1, alleles$Allele2))
                                                  bb <- sum(Genotype == paste0(alleles$Allele2, alleles$Allele2))
                                                  .(AA = aa, Aa = ab, aa = bb)
                                                }, by = .(Population, SNP)
                       ]

# Combine original and new generation for comparison
combined_counts <- rbind(
  cbind(Generation = "Initial", verify_counts),
  cbind(Generation = "Next", new_gen_counts)
)

# Prepare data for plotting - melt to long format
plot_data <- melt(combined_counts, 
                  id.vars = c("Generation", "Population", "SNP"),
                  variable.name = "Genotype",
                  value.name = "Count")

#-------------------------------------------------------------------------------
# HWE comparado entre ambas generaciones: g1 y g2
#-------------------------------------------------------------------------------
# Function to calculate HWE statistics for a generation
calculate_hwe_stats <- function(genotype_data) {
  melted_data <- melt(genotype_data, 
                      id.vars = c("Population", "Individual"),
                      variable.name = "SNP",
                      value.name = "Genotype")
  
  hwe_results <- melted_data[,
                             {
                               alleles <- snp_alleles[SNP == .BY$SNP]
                               
                               # Count observed genotypes
                               obs_AA <- sum(Genotype == paste0(alleles$Allele1, alleles$Allele1))
                               obs_Aa <- sum(Genotype == paste0(alleles$Allele1, alleles$Allele2))
                               obs_aa <- sum(Genotype == paste0(alleles$Allele2, alleles$Allele2))
                               total <- obs_AA + obs_Aa + obs_aa
                               
                               # Calculate allele frequencies
                               p <- (2 * obs_AA + obs_Aa) / (2 * total)
                               q <- 1 - p
                               
                               # Calculate expected counts under HWE
                               exp_AA <- p^2 * total
                               exp_Aa <- 2 * p * q * total
                               exp_aa <- q^2 * total
                               
                               # Chi-square test with Yates' continuity correction
                               chi_sq <- sum((abs(c(obs_AA, obs_Aa, obs_aa) - c(exp_AA, exp_Aa, exp_aa)) - 0.5)^2 / 
                                               c(exp_AA, exp_Aa, exp_aa))
                               df <- 1  # degrees of freedom
                               p_value <- pchisq(chi_sq, df, lower.tail = FALSE)
                               
                               # Return results
                               .(Obs_AA = obs_AA, Obs_Aa = obs_Aa, Obs_aa = obs_aa,
                                 Exp_AA = exp_AA, Exp_Aa = exp_Aa, Exp_aa = exp_aa,
                                 ChiSq = chi_sq, P_value = p_value,
                                 Allele1_Freq = p, Allele2_Freq = q)
                             },
                             by = .(Population, SNP)
  ]
  
  return(hwe_results)
}

# Calculate HWE stats for initial generation
initial_hwe <- calculate_hwe_stats(genotypes_df)

# Generate next generation
next_gen_df <- generate_next_generation(genotypes_df)

# Calculate HWE stats for next generation
next_gen_hwe <- calculate_hwe_stats(next_gen_df)

# Combine results with generation indicator
all_hwe_results <- rbind(
  cbind(Generation = "Initial", initial_hwe),
  cbind(Generation = "Next", next_gen_hwe)
)

# Format the table for better display
formatted_hwe_table <- all_hwe_results[, .(
  Generation,
  Population,
  SNP,
  `Obs AA` = Obs_AA,
  `Obs Aa` = Obs_Aa,
  `Obs aa` = Obs_aa,
  `Exp AA` = round(Exp_AA, 2),
  `Exp Aa` = round(Exp_Aa, 2),
  `Exp aa` = round(Exp_aa, 2),
  `Chi-Sq` = round(ChiSq, 4),
  `P-value` = ifelse(P_value < 0.0001, "<0.0001", 
                     format(round(P_value, 4), nsmall = 4)),
  `Signif` = ifelse(P_value < 0.05, "*", "")
)]

# Print the formatted table
print(formatted_hwe_table)

# Function to create clean table grob
# Function to create clean table grob (corrected version)
# Function to create clean table grob (corrected version)
create_table_grob <- function(data) {
  # Define table theme
  tt <- ttheme_minimal(
    core = list(
      fg_params = list(hjust = 0, x = 0.05),
      bg_params = list(fill = c("white", "#f7f7f7"))
    ),
    colhead = list(
      fg_params = list(fontface = "bold"),
      bg_params = list(fill = "lightgray")
    )
  )
  
  # Create and arrange the table
  grid.arrange(
    tableGrob(data, theme = tt, rows = NULL),
    top = textGrob("HWE Test Results", 
                   gp = gpar(fontsize = 14, fontface = "bold")),
    bottom = textGrob("* = Significant deviation from HWE (p < 0.05)", 
                      gp = gpar(fontsize = 10)),
    nrow = 1
  )
}
# Display the table
create_table_grob(formatted_hwe_table)

# Save to CSV if needed
fwrite(formatted_hwe_table, "hwe_test_results.csv")

################################################################################
# Grafico de resumen genotipos: g1 y g2
#-------------------------------------------------------------------------------
# Create comparative genotype plot with vertical facets: primera y segunda generacion
ggplot(plot_data, aes(x = SNP, y = Count, fill = Genotype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(Generation ~ Population, 
             scales = "free_x", 
             space = "free_x",
             switch = "y") +
  scale_fill_manual(values = genotype_colors,
                    labels = c("Homozygous1 (AA)", 
                               "Heterozygous (Aa)", 
                               "Homozygous2 (aa)")) +
  labs(title = "Genotype Counts: Initial vs Next Generation",
       subtitle = "Top row: Initial generation | Bottom row: Next generation under HWE",
       y = "Number of Individuals",
       x = "SNP Markers",
       fill = "Genotype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y = element_text(face = "bold"),
        panel.spacing = unit(1, "lines")) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) # Make room for labels

################################################################################
# Grafico de resumen HWE g1 y g2
#-------------------------------------------------------------------------------

# Plot the HWE deviations
ggplot(all_hwe_results, aes(x = SNP, y = -log10(P_value), color = Generation)) +
  geom_point(aes(color = Generation, shape = Population), size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  facet_wrap(~ Population, nrow = 1) +
  labs(title = "HWE Test Results Across Generations",
       subtitle = "Dashed line indicates p = 0.05 significance threshold",
       y = "-log10(p-value)",
       x = "SNP Markers") +
  scale_color_manual(values = c("Initial" = "blue", "Next" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")





