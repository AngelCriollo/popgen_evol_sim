################################################################################
#
#               SIMULADOR DE DERIVA
#               Autor: Angel Criollo Rayo & DeepSeek
#               Version 7: 26 - Marzo - 2025
#               Grupo: Citogenetica, Filogenia y Evolucion de Poblaciones
#               Doctorado en Ciencias Biomedicas
#               Asignatura: Genetica Poblacional y Evolutiva (electiva)
#
################################################################################
#
# Este script esta diseñado para analizar los efectos del tamaño poblacional
# en las frecuencias alelicas y genotipicas a lo largo de multiples generaciones. 
# en varios marcadores a la vez, sin tener en cuenta SN, mutación, ligamiento o endogamia
#
################################################################################

library(data.table)
library(ggplot2)

rm(list = setdiff(ls(), lsf.str())) 

################################################################################
# Establezca Parametros 
n_populations <- 5       # EDITAR numero poblaciones
n_individuals <- 30     # EDITAR numero individuos t0
n_new_individuals <- 50 # EDITAR numero individuos cada poblacion/generacion
n_snps <- 3              # EDITAR numero marcadores
nucleotides <- c("A", "C", "T", "G") # EDITAR alelos posibles
n_iterations <- 50       # EDITAR numero generaciones
prob.0 <- c(0.25, 0.5, 0.25) #EDITAR prob. A1A1/A1A2/A2A2 en t0. Defecto en HWE

################################################################################
########################  NO EDITAR DESDE ACA  #################################
################################################################################

# 1. Initialize SNP alleles and tracked alleles
set.seed(123)
snp_alleles <- data.table(SNP = paste0("SNP", 1:n_snps))[
  , c("Allele1", "Allele2") := {
    a1 <- sample(nucleotides, .N, replace = TRUE)
    a2 <- sapply(a1, function(x) sample(setdiff(nucleotides, x), 1))
    .(a1, a2)
  }
][
  , TrackedAllele := ifelse(runif(.N) > 0.5, Allele1, Allele2)
]

# 2. Initialize genotype data
genotypes_df <- data.table(
  Population = rep(paste0("Pop", 1:n_populations), each = n_individuals)
)

for (snp in snp_alleles$SNP) {
  alleles <- unlist(snp_alleles[SNP == snp, .(Allele1, Allele2)])
  genotypes_df[, (snp) := sample(
    c(
      paste0(alleles[1], alleles[1]),
      paste0(alleles[1], alleles[2]),
      paste0(alleles[2], alleles[2])
    ),
    .N,
    replace = TRUE,
    prob = prob.0 # por defecto en HWE
  )]
}

# 3. Function to calculate frequencies
calculate_frequencies <- function(genotypes_df, iteration) {
  melted <- melt(
    genotypes_df,
    id.vars = "Population",
    variable.name = "SNP",
    value.name = "Genotype"
  )
  
  # Standardize heterozygotes (AC = CA)
  melted[, c("Allele1", "Allele2") := {
    alleles <- snp_alleles[SNP == .BY$SNP, .(Allele1, Allele2)]
    list(
      substr(Genotype, 1, 1),
      substr(Genotype, 2, 2)
    )
  }, by = SNP]
  
  # Allele frequencies
  allele_freq <- melted[
    , .(Freq = mean(c(Allele1, Allele2) == snp_alleles[SNP == .BY$SNP, TrackedAllele])),
    by = .(Population, SNP)
  ][, Iteration := iteration]
  
  # Genotype frequencies
  genotype_counts <- melted[
    , .(Count = .N), by = .(Population, SNP, Genotype)
  ][
    , Frequency := Count / sum(Count), by = .(Population, SNP)
  ][, Iteration := iteration]
  
  list(allele_freq = allele_freq, genotype_counts = genotype_counts)
}

# 4. Initialize history with iteration 0
initial_freqs <- calculate_frequencies(genotypes_df, 0)
allele_history <- initial_freqs$allele_freq
genotype_history <- initial_freqs$genotype_counts

# 5. Simulation function
generate_new_genotypes <- function(genotypes_df) {
  melted <- melt(
    genotypes_df,
    id.vars = "Population",
    variable.name = "SNP",
    value.name = "Genotype"
  )
  
  # Sample new genotypes
  allele_pools <- melted[
    , c("Allele1", "Allele2") := {
      alleles <- snp_alleles[SNP == .BY$SNP, .(Allele1, Allele2)]
      list(
        substr(Genotype, 1, 1),
        substr(Genotype, 2, 2)
      )
    }, by = SNP
  ][
    , .(Allele = c(Allele1, Allele2)), by = .(Population, SNP)
  ]
  
  new_individuals <- paste0("NewInd_", 1:n_new_individuals)
  new_genotypes <- allele_pools[
    , {
      sampled_alleles <- sample(Allele, 2 * n_new_individuals, replace = TRUE)
      genotypes <- paste0(
        sampled_alleles[seq(1, length(sampled_alleles), by = 2)],
        sampled_alleles[seq(2, length(sampled_alleles), by = 2)]
      )
      # Standardize heterozygotes
      genotypes <- sapply(genotypes, function(gt) {
        sorted <- sort(unlist(strsplit(gt, "")))
        paste0(sorted[1], sorted[2])
      })
      .(Genotype = genotypes, Individual = new_individuals)
    },
    by = .(Population, SNP)
  ]
  
  dcast(new_genotypes, Population + Individual ~ SNP, value.var = "Genotype")
}

# 6. Run simulation
for (i in 1:n_iterations) {  # 5 iterations
  genotypes_df <- generate_new_genotypes(genotypes_df)
  current_freqs <- calculate_frequencies(genotypes_df, i)
  allele_history <- rbind(allele_history, current_freqs$allele_freq)
  genotype_history <- rbind(genotype_history, current_freqs$genotype_counts)
}

################################################################################
#
# RESULTADOS 
#
################################################################################
allele_res <- subset(allele_history, !SNP %in% "Individual")
genot_res <- subset(genotype_history, !SNP %in% "Individual")

# 7. Dinamica de las frecuencias alelicas
ggplot(allele_res, aes(x = Iteration, y = Freq, color = Population)) +
  geom_line() +
  facet_wrap(~ SNP) +
  labs(title = "Tracked Allele Frequency Dynamics") +
  theme_minimal()

# 8. Dinamica de las frecuencias genotipicas
# 1. Define consistent color scheme
genotype_colors <- c(
  "Homozygous1" = "#1b9e77",  # Green
  "Heterozygous" = "#d95f02",  # Orange
  "Homozygous2" = "#7570b3"   # Purple
)

# 2. Classify genotypes for each SNP
genot_res[, GenotypeClass := fcase(
  Genotype == paste0(snp_alleles[SNP == .BY$SNP, Allele1], snp_alleles[SNP == .BY$SNP, Allele1]),
  "Homozygous1",
  Genotype == paste0(snp_alleles[SNP == .BY$SNP, Allele2], snp_alleles[SNP == .BY$SNP, Allele2]),
  "Homozygous2",
  default = "Heterozygous"
), by = SNP]

# 3. Create annotation data for facet legends
facet_legends <- unique(snp_alleles[, .(SNP, Allele1, Allele2)])
facet_legends[, Label := paste0(SNP, ": ", Allele1, Allele1, " (H1) | ", 
                                Allele1, Allele2, " (Het) | ",
                                Allele2, Allele2, " (H2)")]

# 4. Create the plot
ggplot(genot_res, aes(x = Iteration, y = Frequency, color = GenotypeClass)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 0.1) +
  
  # Faceting and labels
  facet_grid(Population ~ SNP) +
  labs(
    title = "Genotype Frequency Dynamics",
    x = "Generation",
    y = "Frequency",
    color = "Genotype Class"
  ) +
  
  # Color scheme
  scale_color_manual(values = genotype_colors) +
  
  # Theme adjustments
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  
  # Add facet-specific legends
  geom_text(
    data = facet_legends,
    aes(x = -Inf, y = Inf, label = Label),
    hjust = -0.05, vjust = 1.5,
    size = 2.8, color = "black",
    inherit.aes = FALSE
  ) +
  
  # Add universal color guide
  annotate(
    "text", x = -Inf, y = -Inf,
    label = "Color Key: H1=Green | Het=Orange | H2=Purple",
    hjust = -0.05, vjust = -1,
    size = 3, color = "black"
  )





