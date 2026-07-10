### This script loads metadata from Grady, GENOA, and GSE87571 (Swedish) datasets
### and generates a boxplot comparing the age distributions across all three cohorts.

# Set the working directory
setwd("~/Desktop/Capra Lab/Thesis Project/Aim_1/")

# Load libraries
library(tidyverse)
library(GEOquery)

# ============================================================================
# Load and process Grady metadata (GSE72680)
# African Americans from the Grady Trauma Project
# ============================================================================
message("Loading Grady (GSE72680) metadata...")
gse_grady <- getGEO("GSE72680", GSEMatrix = TRUE)
gse_grady_data <- gse_grady[[1]]
phenos_grady <- pData(gse_grady_data)

# Extract age
phenos_grady$age <- as.numeric(gsub(".*age: ([0-9]+).*", "\\1", phenos_grady$characteristics_ch1.1))

# Clean up memory
rm(gse_grady, gse_grady_data)

# ============================================================================
# Load and process GENOA metadata (GSE210254)
# African Americans from the GENOA study (450K array)
# ============================================================================
message("Loading GENOA (GSE210254) metadata...")
gse_genoa <- getGEO("GSE210254", GSEMatrix = TRUE)
gse_genoa_data <- gse_genoa[[1]]
phenos_genoa <- pData(gse_genoa_data)

# Extract age
phenos_genoa$age <- as.numeric(gsub(".+:\\s*([0-9.]+).*", "\\1", phenos_genoa$characteristics_ch1.2))

# Clean up memory
rm(gse_genoa, gse_genoa_data)

# ============================================================================
# Load and process GSE87571 metadata (Swedish cohort)
# White individuals from Sweden
# ============================================================================
message("Loading GSE87571 (Swedish) metadata...")
gse_swedish <- getGEO("GSE87571", GSEMatrix = TRUE)
gse_swedish_data <- gse_swedish[[1]]
phenos_swedish <- pData(gse_swedish_data)

# Extract age
phenos_swedish$age <- as.numeric(gsub(".*age: ([0-9]+).*", "\\1", phenos_swedish$characteristics_ch1.1))

# Clean up memory
rm(gse_swedish, gse_swedish_data)

# ============================================================================
# Combine data for plotting
# ============================================================================
message("Preparing data for visualization...")

# Create a combined dataframe with dataset labels
combined_ages <- rbind(
    data.frame(
        dataset = "Grady (AA)",
        age = phenos_grady$age
    ),
    data.frame(
        dataset = "GENOA (AA)",
        age = phenos_genoa$age
    ),
    data.frame(
        dataset = "GSE87571 (Swedish)",
        age = phenos_swedish$age
    )
)

# Remove any NA values
combined_ages <- combined_ages %>% filter(!is.na(age))

# ============================================================================
# Generate boxplot
# ============================================================================
message("Generating boxplot...")

# Define colors for each dataset
dataset_colors <- c(
    "Grady (AA)" = "#1f77b4",
    "GENOA (AA)" = "#ff7f0e",
    "GSE87571 (Swedish)" = "#2ca02c"
)

# Create the boxplot
age_boxplot <- ggplot(combined_ages, aes(x = dataset, y = age, fill = dataset)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2, outlier.alpha = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.2, size = 1) +
    scale_fill_manual(values = dataset_colors) +
    labs(
        title = "Age Distribution Across Methylation Cohorts",
        subtitle = "Comparison of Grady, GENOA, and Swedish (GSE87571) datasets",
        x = "Dataset",
        y = "Age (years)"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
    )

# Display the plot
print(age_boxplot)

# ============================================================================
# Print summary statistics
# ============================================================================
message("\n=== Summary Statistics ===")
summary_stats <- combined_ages %>%
    group_by(dataset) %>%
    summarise(
        n = n(),
        mean_age = mean(age, na.rm = TRUE),
        median_age = median(age, na.rm = TRUE),
        sd_age = sd(age, na.rm = TRUE),
        min_age = min(age, na.rm = TRUE),
        max_age = max(age, na.rm = TRUE),
        .groups = "drop"
    )

print(summary_stats)

# ============================================================================
# Subset to older individuals (age >= 55) to match MAGENTA cohort
# ============================================================================
message("\n=== Subsetting to Age >= 55 ===")

# Filter to older individuals
combined_ages_older <- combined_ages %>% filter(age >= 55)

# Generate boxplot for older subset
age_boxplot_older <- ggplot(combined_ages_older, aes(x = dataset, y = age, fill = dataset)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2, outlier.alpha = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.2, size = 1) +
    scale_fill_manual(values = dataset_colors) +
    labs(
        title = "Age Distribution Across Methylation Cohorts (Age >= 55)",
        subtitle = "Comparison of Grady, GENOA, and Swedish (GSE87571) datasets",
        x = "Dataset",
        y = "Age (years)"
    ) +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
    )

# Display the plot
print(age_boxplot_older)

# Print summary statistics for older subset
message("\n=== Summary Statistics (Age >= 55) ===")
summary_stats_older <- combined_ages_older %>%
    group_by(dataset) %>%
    summarise(
        n = n(),
        mean_age = mean(age, na.rm = TRUE),
        median_age = median(age, na.rm = TRUE),
        sd_age = sd(age, na.rm = TRUE),
        min_age = min(age, na.rm = TRUE),
        max_age = max(age, na.rm = TRUE),
        .groups = "drop"
    )

print(summary_stats_older)
