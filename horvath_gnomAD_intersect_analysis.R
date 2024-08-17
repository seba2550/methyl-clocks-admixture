## This is a script to analyze the intersection of Horvath clock CpGs and gnomAD (3.0) variants (may include SNPs and indels)
# Load some libraries
library(tidyverse)
library(ggpubr)

# Set the working directory
setwd("/Users/sgonzalez/Desktop/Capra Lab/Thesis Project/Aim_1")

# Load the bed file holding the intersect of the Horvath clock CpGs and the gnomAD variants
horvath_gnomAD <- read.table("separate_files/horvath_clock_cpgs_gnomAD_variants.bed")

# Get the allele frequencies for each SNP by doing a split on the INFO column
horvath_gnomAD$AF <- as.numeric(str_split_i(str_split_i(horvath_gnomAD$V14, ";", i = 3), "=", i = 2))

# Get additional information from the INFO space
horvath_gnomAD$type_of_variant <- str_split_i(str_split_i(horvath_gnomAD$V14, ";", -1), "\\|", i = 2)

horvath_gnomAD$gene_symbol <- str_split_i(str_split_i(horvath_gnomAD$V14, ";", -1), "\\|", i = 4)

horvath_gnomAD$gene_id <- str_split_i(str_split_i(horvath_gnomAD$V14, ";", -1), "\\|", i = 5)

# Count the number of distinct CpGs affected by gnomAD SNPs
length(unique(horvath_gnomAD$V4))

# Count the number of distinct SNPs in gnomAD, using the chromosome and position columns
length(unique(paste(horvath_gnomAD$V7, horvath_gnomAD$V8, sep = ":")))


##horvath_gnomAD$variant_effect_location_2 <- str_split_i(str_split_i(horvath_gnomAD$V14, ";", -1), "\\|", i = -10)

# Generate a plot showing the number of SNPs at all allele frequencies in gnomAD (using log10 scale)
ggplot(horvath_gnomAD, aes(x = AF)) +
  geom_histogram(fill = "blue", color = "black") +
  labs(x = expression(paste("Allele Frequency (", log[10], " scale)")), 
       y = "Number of Clock CpG Disrupting SNPs", 
       title = "Clock CpG-disrupting SNP Distribution by Allele Frequency") +
  scale_x_continuous(trans = "log10") +
  theme_minimal()

# ggpubr version of the above plot
gghistogram(
  data = horvath_gnomAD,
  x = "AF",
  fill = "blue",
  color = "black",
  ylab = "Number of Clock CpG Disrupting SNPs",
) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Allele Frequency (log10 scale)")

# Generate a plot showing the number of SNPs at positive CpGs vs negative CpGs
ggplot(horvath_gnomAD, aes(x = V5, fill = V5)) +
  geom_bar() +
  labs(x = "Marginal Age Relationship", 
       y = "Number of SNPs", 
       title = "Marginal Age Relationship of CpGs Affected by gnomAD SNPs") +
  theme_minimal() +
  theme(legend.position = "none")

# Generate a plot showing the type of variant distribution.
ggplot(horvath_gnomAD, aes(x = type_of_variant, fill = type_of_variant)) +
  geom_bar() +
  labs(x = "Type of Variant", 
       y = "Number of Variants", 
       title = "Type of Variant") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Generate a plot showing the variant effect location distribution.
ggplot(horvath_gnomAD, aes(x = variant_effect_location_2, fill = variant_effect_location_2)) +
  geom_bar() +
  labs(x = "Variant Effect Location", 
       y = "Number of Variants", 
       title = "Variant Effect Location") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Generate a plot showing the allele frequency distribution by type of variant
ggplot(horvath_gnomAD, aes(x = AF, fill = type_of_variant)) +
  geom_histogram() +
  labs(x = expression(paste("Allele Frequency (", log[10], " scale)")), 
       y = "Number of SNPs", 
       title = "SNP Distribution by Allele Frequency and Type of Variant") +
  scale_x_continuous()

########### Start the 450K Methylation Array analysis
# Load the result of intersecting the 450K array CpGs with gnomAD variants (3.0)
array_gnomAD <- read.table("methylation_450k_gnomAD_intersect.bed")

# Gather the allele frequencies for the 450K CpG-disrupting SNPs from gnomAD
array_gnomAD$AF <- as.numeric(str_split_i(str_split_i(array_gnomAD$V12, ";", i = 3), "=", i = 2))

# Get some additional data from the INFO column, just as we did for the Horvath clock CpGs
array_gnomAD$type_of_variant <- str_split_i(str_split_i(array_gnomAD$V12, ";", -1), "\\|", i = 2)

# Generate a plot showing the number of SNPs at all allele frequencies in gnomAD (using log10 scale)
ggplot(array_gnomAD, aes(x = AF)) +
  geom_histogram(fill = "blue", color = "black") +
  labs(x = expression(paste("Allele Frequency (", log[10], " scale)")), 
       y = "Number of SNPs", 
       title = "SNP Distribution by Allele Frequency") +
  scale_x_continuous(trans = "log10") +
     theme_minimal()

# Generate a plot showing the type of variant distribution.
ggplot(array_gnomAD, aes(x = type_of_variant, fill = type_of_variant)) +
  geom_bar() +
  labs(x = "Type of Variant", 
       y = "Number of Variants", 
       title = "Type of Variant") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Generate a plot showing the allele frequency distribution by type of variant
ggplot(array_gnomAD, aes(x = AF, fill = type_of_variant)) +
  geom_histogram() +
  labs(x = expression(paste("Allele Frequency (", log[10], " scale)")), 
       y = "Number of SNPs", 
       title = "SNP Distribution by Allele Frequency and Type of Variant") +
  scale_x_continuous(trans = "log10") +
  theme(legend.position = "none")

# Generate a plot comparing the count of SNPs at clock CpGs and 450K CpGs
ggplot() +
  geom_bar(data = horvath_gnomAD, aes(x = "Horvath Clock CpGs", fill = "Horvath Clock CpGs")) +
  geom_bar(data = array_gnomAD, aes(x = "450K CpGs", fill = "450K CpGs")) +
  labs(x = "CpG Type", 
       y = "Number of SNPs", 
       title = "Comparison of SNPs at Horvath Clock CpGs and 450K CpGs") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_y_continuous(trans = "log10")


##### Some stat questions
# How many unique clock CpGs are affected by gnomAD SNPs?
length(unique(horvath_gnomAD$V4))

# What is the proportion of Horvath Clock CpGs affected?
(length(unique(horvath_gnomAD$V4))/353) * 100

# How many unique 450K CpGs are affected by gnomAD SNPs?
length(unique(array_gnomAD$V4))

# What proportion of 450K CpGs are affected?
(length(unique(array_gnomAD$V4))/485512) * 100
