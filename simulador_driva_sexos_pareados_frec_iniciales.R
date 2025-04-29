# Load libraries
library(tidyverse)

# Set parameters
set.seed(123)
n_populations <- 3
n_individuals <- 50
n_snps <- 5
n_generations <- 20
initial_freqs <- rep(0.5, n_snps)# runif(n_snps, 0.3, 0.7) 
nucleotides <- c("A", "C", "G", "T")

# Generate SNP markers: random allele pairs
snp_alleles <- replicate(n_snps, sample(nucleotides, 2, replace = FALSE), simplify = FALSE)
names(snp_alleles) <- paste0("SNP", 1:n_snps)

# Functions
calculate_allele_freq <- function(genotypes, alleles) {
  allele_counts <- table(unlist(strsplit(paste(genotypes, collapse=""), "")))
  total_alleles <- sum(allele_counts)
  freq <- ifelse(alleles[1] %in% names(allele_counts), allele_counts[[alleles[1]]], 0) / total_alleles
  return(freq)
}

calculate_ho <- function(genotypes) {
  heterozygotes <- sum(substr(genotypes, 1, 1) != substr(genotypes, 2, 2))
  ho <- heterozygotes / length(genotypes)
  return(ho)
}

calculate_he <- function(freq) {
  he <- 2 * freq * (1 - freq)
  return(he)
}

hwe_test <- function(genotypes, alleles) {
  counts <- table(genotypes)
  homo1 <- counts[paste0(alleles[1], alleles[1])]
  homo2 <- counts[paste0(alleles[2], alleles[2])]
  het1 <- counts[paste0(alleles[1], alleles[2])]
  het2 <- counts[paste0(alleles[2], alleles[1])]
  hets <- sum(het1, het2, na.rm=TRUE)
  
  homo1 <- ifelse(is.na(homo1), 0, homo1)
  homo2 <- ifelse(is.na(homo2), 0, homo2)
  
  n <- homo1 + homo2 + hets
  if (n == 0) return(NA)
  
  p <- (2*homo1 + hets) / (2*n)
  q <- 1 - p
  
  expected <- c(p^2, 2*p*q, q^2) * n
  observed <- c(homo1, hets, homo2)
  
  if (any(expected < 5)) return(NA)
  
  chisq <- sum((observed - expected)^2 / expected)
  pvalue <- pchisq(chisq, df=1, lower.tail=FALSE)
  return(pvalue)
}

initialize_generation <- function(n_individuals, snp_alleles, initial_freqs) {
  genotypes <- map2(snp_alleles, initial_freqs, function(alleles, freq) {
    total_alleles <- n_individuals * 2
    
    n_allele1 <- round(freq * total_alleles)
    n_allele2 <- total_alleles - n_allele1
    
    allele_pool <- c(rep(alleles[1], n_allele1), rep(alleles[2], n_allele2))
    
    if (length(allele_pool) != total_alleles) {
      stop("Allele pool size mismatch")
    }
    
    shuffled_alleles <- sample(allele_pool)
    genotypes <- sapply(seq(1, total_alleles, by=2), function(i) {
      paste0(sort(c(shuffled_alleles[i], shuffled_alleles[i+1])), collapse="")
    })
    
    return(genotypes)
  }) %>%
    set_names(names(snp_alleles)) %>%
    as_tibble()
  
  sexes <- sample(c("M", "F"), n_individuals, replace=TRUE)
  if (abs(sum(sexes=="M") - sum(sexes=="F")) > 1) {
    n_m <- floor(n_individuals/2)
    n_f <- n_individuals - n_m
    sexes <- sample(c(rep("M", n_m), rep("F", n_f)))
  }
  
  tibble(ID=1:n_individuals, Sex=sexes) %>% bind_cols(genotypes)
}

next_generation <- function(prev_gen, snp_alleles) {
  males <- prev_gen %>% filter(Sex == "M")
  females <- prev_gen %>% filter(Sex == "F")
  
  n_individuals <- nrow(prev_gen)
  new_gen <- tibble(ID=1:n_individuals, Sex=sample(c("M", "F"), n_individuals, replace=TRUE))
  
  for (snp in names(snp_alleles)) {
    male_pool <- unlist(strsplit(paste(males[[snp]], collapse=""), ""))
    female_pool <- unlist(strsplit(paste(females[[snp]], collapse=""), ""))
    
    male_alleles <- sample(male_pool, n_individuals, replace=TRUE)
    female_alleles <- sample(female_pool, n_individuals, replace=TRUE)
    
    new_gen[[snp]] <- mapply(function(m,f) paste0(sort(c(m,f)), collapse=""), male_alleles, female_alleles)
  }
  
  return(new_gen)
}

# Initialize populations
populations <- map(1:n_populations, ~initialize_generation(n_individuals, snp_alleles, initial_freqs))
names(populations) <- paste0("Pop", 1:n_populations)

# Store results
genotype_history <- list()
genotype_counts <- list()

# Main loop
for (pop_name in names(populations)) {
  pop_hist <- list()
  pop_counts <- list()
  pop <- populations[[pop_name]]
  
  for (gen in 0:n_generations) {
    pop$Generation <- gen
    pop$Population <- pop_name
    
    pop_hist[[as.character(gen)]] <- pop
    
    for (snp in names(snp_alleles)) {
      alleles <- snp_alleles[[snp]]
      freq <- calculate_allele_freq(pop[[snp]], alleles)
      ho <- calculate_ho(pop[[snp]])
      he <- calculate_he(freq)
      pvalue <- hwe_test(pop[[snp]], alleles)
      
      pop_counts[[length(pop_counts)+1]] <- tibble(
        Population = pop_name,
        Generation = gen,
        SNP = snp,
        Allele = alleles[1],
        Frequency = freq,
        Ho = ho,
        He = he,
        HWE_pvalue = pvalue
      )
    }
    
    if (gen != n_generations) {
      pop <- next_generation(pop, snp_alleles)
    }
  }
  
  genotype_history[[pop_name]] <- bind_rows(pop_hist)
  genotype_counts[[pop_name]] <- bind_rows(pop_counts)
}

# Final dataframes
genotype_history_df <- bind_rows(genotype_history)
genotype_counts_df <- bind_rows(genotype_counts)

################################################################################
# Allele frequency dynamics
ggplot(genotype_counts_df, aes(x=Generation, y=Frequency, color=Population)) +
  geom_line() +
  facet_wrap(~SNP, scales="free_y") +
  theme_minimal() +
  ggtitle("Allele frequency dynamics over generations")

# Observed vs expected heterozygosity
genotype_counts_long <- genotype_counts_df %>%
  pivot_longer(cols=c(Ho, He), names_to="Type", values_to="Value")

ggplot(genotype_counts_long, aes(x=Generation, y=Value, color=Population, linetype=Type)) +
  geom_line() +
  facet_wrap(~SNP, scales="free_y") +
  theme_minimal() +
  ggtitle("Observed vs Expected Heterozygosity")

# HWE p-values
ggplot(genotype_counts_df, aes(x=Generation, y=HWE_pvalue, color=Population)) +
  geom_line() +
  facet_wrap(~SNP, scales="free_y") +
  geom_hline(yintercept=0.05, linetype="dashed", color="red") +
  theme_minimal() +
  ggtitle("P-values of HWE Chi-square tests")














