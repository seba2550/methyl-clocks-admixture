# Load the Functions you need to calculate the PCClocks (you need to change the path to the directory
# where you installed the code)
library(tidyverse)
library(readxl)
library(GEOquery)
library(broom)
library(data.table)

clocksDir <- "/Users/sgonzalez/Desktop/Capra Lab/Thesis Project/Aim_1/PC-Clocks-main/" # where the clocks directory was downloaded to
# must end with '/'

source(paste(clocksDir, "run_calcPCClocks.R", sep = ""))
source(paste(clocksDir, "run_calcPCClocks_Accel.R", sep = ""))

# Load the file with your pheno data and methylation data in it (Here we used the example data included)
load(paste(clocksDir, "Example_PCClock_Data_final.RData", sep = ""))
PCClock_DNAmAge <- calcPCClocks(
  path_to_PCClocks_directory = clocksDir,
  datMeth = datMeth,
  datPheno = datPheno
)

# IMPORTANT FORMATTING NOTE: If you are not using the example methylation and Pheno data, you will need to have specific
#     formatting. Please ensure that your Methylation dataframe/ matrix is of the methylation beta values and row names
#     are sample names, and column names are CpGs.
#     For the pheno data, ensure that the data frame/ matrix has rows as samples, and columns as whatever phenotype
#     variables you have/ wish. This can also include the original CpG clocks if you used the online Horvath calculator
#     as well. HOWEVER, the pheno matrix MUST have a column named "Age", and a column named "Female" (capital required),
#     especially if you want to calculate GrimAge and its components. Female = 1 is F and Female = 0 is M.
#
#     If you don't have this information, you can also just set it so Females is all 1 (all samples labeled female) and
#     all the same Age. Just know that this won't be accurate for PCGrimAge or components, and that you can't run
#     the acceleration calculations with calcPCClock_Accel.
#
#     The code below is going to ask if you would like to check the order of your datPheno and datMeth samples to ensure
#     they line up. For this to work, you will need to type the column name of datPheno with the names of the samples or
#     'skip'.

# Load MAGENTA betas
magenta_betas <- readRDS("/Users/sgonzalez/Desktop/Capra Lab/Thesis Project/Aim_1/betaMatrices/normalizedBetas/beta_QGCDPB_combined.rds")

# Load the metadata
sample_metadata <- read_xlsx("/Users/sgonzalez/Desktop/Capra Lab/Thesis Project/Aim_1/ADmethy_pheno.xlsx")

# We have some leftover samples for which there are no methylation data. Let's filter them out
sample_metadata <- subset(sample_metadata, Beta_ID %in% colnames(magenta_betas))

# Change the cohort names in the metadata to human-readable.
sample_metadata <- sample_metadata %>%
  mutate(COHORT = case_when(
    COHORT == "CuADI" ~ "Cuban",
    COHORT == "NHW" ~ "Non-Hispanic White",
    COHORT == "PERUVIAN" ~ "Peruvian",
    COHORT == "PRADI" ~ "Puerto Rican",
    COHORT == "REAAADI" ~ "African American",
    TRUE ~ COHORT # Keep the original value if no match is found
  ))

# Modify the rownames to remove everything after the underscore
rownames(magenta_betas) <- gsub("_.*", "", rownames(magenta_betas))

# Re-format some columns so that they work in the PC clocks function
sample_metadata <- sample_metadata %>% rename(Age = AGE_OF_EXAM)
sample_metadata <- sample_metadata %>% mutate(Female = case_when(
  SEX == "F" ~ 1,
  SEX == "M" ~ 0
))

# Transpose the betas into PC clock format
magenta_betas_transposed <- as.data.frame(t(magenta_betas))


# Get the PC Clocks values and the PC Clock Acceleration values
PCClock_DNAmAge <- calcPCClocks(
  path_to_PCClocks_directory = clocksDir,
  datMeth = magenta_betas_transposed,
  datPheno = sample_metadata
)
# in order to calculate Acceleration below, you will need a column called "Age", just as was needed for PCGrimAge.
PCClock_DNAmAge <- calcPCClocks_Accel(PCClock_DNAmAge)


######### Analysis of the above results
# Function to calculate MAE and MSE stratified by cohort
calculate_metrics_toy_dataset <- function(data) {
  metrics <- data %>%
    group_by(group) %>%
    summarise(
      n = n(),
      MAE = mean(abs(PCHorvath1 - Age)),
      MSE = mean((PCHorvath1 - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCHorvath1 - Age))
    ) %>%
    ungroup()

  return(metrics)
}
calculate_metrics_all <- function(data) {
  # Calculate metrics by cohort
  metrics <- data %>%
    group_by(COHORT) %>%
    summarise(
      n = n(),
      MAE = median(abs(PCHorvath1 - Age)),
      MSE = mean((PCHorvath1 - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCHorvath1 - Age))
    ) %>%
    ungroup()

  return(metrics)
}

calculate_metrics_controls <- function(data) {
  # Filter for controls only
  control_data <- data[data$STATUS == "CONTROL", ]

  # Calculate metrics by cohort
  metrics <- control_data %>%
    group_by(COHORT) %>%
    summarise(
      n = n(),
      MAE = median(abs(PCHorvath1 - Age)),
      MSE = mean((PCHorvath1 - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCHorvath1 - Age))
    ) %>%
    ungroup()

  return(metrics)
}
calculate_metrics_controls_skin_blood <- function(data) {
  # Filter for controls only
  control_data <- data[data$STATUS == "CONTROL", ]

  # Calculate metrics by cohort
  metrics <- control_data %>%
    group_by(COHORT) %>%
    summarise(
      n = n(),
      MAE = median(abs(PCHorvath2 - Age)),
      MSE = mean((PCHorvath2 - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCHorvath2 - Age))
    ) %>%
    ungroup()

  return(metrics)
}

calculate_metrics_controls_hannum <- function(data) {
  # Filter for controls only
  control_data <- data[data$STATUS == "CONTROL", ]

  # Calculate metrics by cohort
  metrics <- control_data %>%
    group_by(COHORT) %>%
    summarise(
      n = n(),
      MAE = median(abs(PCHannum - Age)),
      MSE = mean((PCHannum - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCHannum - Age))
    ) %>%
    ungroup()

  return(metrics)
}

calculate_metrics_controls_phenoage <- function(data) {
  # Filter for controls only
  control_data <- data[data$STATUS == "CONTROL", ]

  # Calculate metrics by cohort
  metrics <- control_data %>%
    group_by(COHORT) %>%
    summarise(
      n = n(),
      MAE = median(abs(PCPhenoAge - Age)),
      MSE = mean((PCPhenoAge - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCPhenoAge - Age))
    ) %>%
    ungroup()

  return(metrics)
}


# Example usage:
results_all <- calculate_metrics_all(PCClock_DNAmAge)
results_controls <- calculate_metrics_controls(PCClock_DNAmAge)
results_horvath_skin_blood_controls <- calculate_metrics_controls_skin_blood(PCClock_DNAmAge)
results_controls_hannum <- calculate_metrics_controls_hannum(PCClock_DNAmAge)
results_controls_phenoage <- calculate_metrics_controls_phenoage(PCClock_DNAmAge)
print(results_all)
print(results_controls)
print(results_horvath_skin_blood_controls)
print(results_controls_hannum)
print(results_controls_phenoage)

horvath_all_correlations <- PCClock_DNAmAge %>%
  group_by(COHORT) %>%
  summarise(correlation = cor(PCHorvath1, Age))

horvath_controls_correlations <- PCClock_DNAmAge %>%
  filter(STATUS == "CONTROL") %>%
  group_by(COHORT) %>%
  summarise(correlation = cor(PCHorvath1, Age))

horvath_skin_blood_controls_correlations <- PCClock_DNAmAge %>%
  filter(STATUS == "CONTROL") %>%
  group_by(COHORT) %>%
  summarise(correlation = cor(PCHorvath2, Age))
hannum_controls_correlations <- PCClock_DNAmAge %>%
  filter(STATUS == "CONTROL") %>%
  group_by(COHORT) %>%
  summarise(correlation = cor(PCHannum, Age))

# Prepare data for the plot: rename cohort and set factor order
plot_data <- PCClock_DNAmAge %>%
  mutate(
    COHORT = ifelse(COHORT == "Non-Hispanic White", "White", COHORT),
    COHORT = factor(COHORT, levels = c("White", "African American", "Puerto Rican", "Cuban", "Peruvian"))
  )
plot_data_controls <- PCClock_DNAmAge %>%
  filter(STATUS == "CONTROL") %>%
  mutate(
    COHORT = ifelse(COHORT == "Non-Hispanic White", "White", COHORT),
    COHORT = factor(COHORT, levels = c("White", "African American", "Puerto Rican", "Cuban", "Peruvian"))
  )

# Calculate correlation and MAE for each cohort to display in the plot
cor_labels <- plot_data %>%
  group_by(COHORT) %>%
  summarise(
    r = cor(PCHorvath1, Age, use = "complete.obs"),
    mae = mean(abs(PCHorvath1 - Age), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("R = ", round(r, 3), "\nMAE = ", round(mae, 2)))
ggplot(data = plot_data, mapping = aes(x = PCHorvath1, y = Age)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = F) +
  geom_text(
    data = cor_labels, aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.1, size = 4, fontface = "italic", inherit.aes = FALSE
  ) +
  theme_classic() +
  ggtitle("MAGENTA Cases and Controls Correlations (PC Horvath Clock)") +
  xlab("PC Horvath DNAmAge Estimate") +
  ylab("Chronological Age") +
  facet_wrap(~COHORT)

cor_labels_controls <- plot_data_controls %>%
  group_by(COHORT) %>%
  summarise(
    r = cor(PCHorvath1, Age, use = "complete.obs"),
    mae = mean(abs(PCHorvath1 - Age), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("R = ", round(r, 3), "\nMAE = ", round(mae, 2)))

ggplot(data = plot_data_controls, mapping = aes(x = PCHorvath1, y = Age)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = F) +
  geom_text(
    data = cor_labels_controls, aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.1, size = 4, fontface = "italic", inherit.aes = FALSE
  ) +
  theme_classic() +
  ggtitle("MAGENTA Controls Correlations (PC Horvath Clock)") +
  xlab("PC Horvath DNAmAge Estimate") +
  ylab("Chronological Age") +
  facet_wrap(~COHORT)
cor_labels_horvath2_controls <- plot_data_controls %>%
  group_by(COHORT) %>%
  summarise(
    r = cor(PCHorvath2, Age, use = "complete.obs"),
    mae = mean(abs(PCHorvath2 - Age), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("R = ", round(r, 3), "\nMAE = ", round(mae, 2)))

ggplot(data = plot_data_controls, mapping = aes(x = PCHorvath2, y = Age)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "blue", se = F) +
  geom_text(
    data = cor_labels_horvath2_controls, aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.1, size = 4, fontface = "italic", inherit.aes = FALSE
  ) +
  ggtitle("MAGENTA Controls Correlations (PC Horvath Clock Skin and Blood)") +
  xlab("PC Horvath DNAmAge Estimate") +
  ylab("Chronological Age") +
  facet_wrap(~COHORT)
cor_labels_hannum_controls <- plot_data_controls %>%
  group_by(COHORT) %>%
  summarise(
    r = cor(PCHannum, Age, use = "complete.obs"),
    mae = mean(abs(PCHannum - Age), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("R = ", round(r, 3), "\nMAE = ", round(mae, 2)))

ggplot(data = plot_data_controls, mapping = aes(x = PCHannum, y = Age)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "blue", se = F) +
  geom_text(
    data = cor_labels_hannum_controls, aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.1, size = 4, fontface = "italic", inherit.aes = FALSE
  ) +
  ggtitle("MAGENTA Controls Correlations (PC Hannum Clock Skin and Blood)") +
  xlab("PC Hannum DNAmAge Estimate") +
  ylab("Chronological Age") +
  facet_wrap(~COHORT)

# Plot correlations for MAGENTA Cases by cohort
plot_data_cases <- plot_data %>% filter(STATUS == "AD")

cor_labels_cases <- plot_data_cases %>%
  group_by(COHORT) %>%
  summarise(
    r = cor(PCHorvath1, Age, use = "complete.obs"),
    mae = mean(abs(PCHorvath1 - Age), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("R = ", round(r, 3), "\nMAE = ", round(mae, 2)))

ggplot(data = plot_data_cases, mapping = aes(x = PCHorvath1, y = Age)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", color = "blue", se = F) +
  geom_text(
    data = cor_labels_cases, aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.1, size = 4, fontface = "italic", inherit.aes = FALSE
  ) +
  ggtitle("MAGENTA Cases Correlations (PC Horvath Clock)") +
  xlab("PC Horvath DNAmAge Estimate") +
  ylab("Chronological Age") +
  facet_wrap(~COHORT)

horvath_cases_correlations <- PCClock_DNAmAge %>%
  filter(STATUS == "AD") %>%
  group_by(COHORT) %>%
  summarise(correlation = cor(PCHorvath1, Age))


t.test(PCClock_DNAmAge$PCHorvath1Resid ~ PCClock_DNAmAge$STATUS)
t.test(PCClock_DNAmAge$PCHorvath2Resid ~ PCClock_DNAmAge$STATUS)
ggplot(data = PCClock_DNAmAge) +
  geom_boxplot(mapping = aes(x = STATUS, y = PCHorvath1Resid, fill = STATUS), outlier.shape = NA) +
  geom_jitter(mapping = aes(x = STATUS, y = PCHorvath1Resid), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Cohort") +
  ylab("Age Acceleration (PC Horvath DNAm Age - Chronological Age)") +
  theme_classic()
ggplot(data = PCClock_DNAmAge) +
  geom_boxplot(mapping = aes(x = STATUS, y = PCHorvath2Resid, fill = STATUS), outlier.shape = NA) +
  geom_jitter(mapping = aes(x = STATUS, y = PCHorvath2Resid), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Cohort") +
  ylab("Age Acceleration (PC Horvath Skin and Blood DNAm Age - Chronological Age)") +
  theme_classic()

ggplot(data = PCClock_DNAmAge) +
  geom_boxplot(mapping = aes(x = STATUS, y = PCHorvath1Resid, fill = STATUS), outlier.shape = NA) +
  geom_jitter(mapping = aes(x = STATUS, y = PCHorvath1Resid), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Cohort") +
  ylab("Age Acceleration (PC Horvath DNAm Age - Chronological Age)") +
  theme_classic() +
  facet_wrap(~ factor(COHORT, c("Non-Hispanic White", "African American", "Cuban", "Puerto Rican", "Peruvian")))
t_test_cohort_status_pc_horvath <- PCClock_DNAmAge %>%
  group_by(COHORT) %>%
  nest() %>%
  mutate(
    t_test = map(data, ~ t.test(PCHorvath1Resid ~ STATUS, data = .x) %>% broom::tidy())
  ) %>%
  unnest(t_test)

ggplot(data = PCClock_DNAmAge) +
  geom_boxplot(mapping = aes(x = STATUS, y = PCHorvath2Resid, fill = STATUS), outlier.shape = NA) +
  geom_jitter(mapping = aes(x = STATUS, y = PCHorvath2Resid), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Cohort") +
  ylab("Age Acceleration (PC Horvath Skin and Blood DNAm Age - Chronological Age)") +
  theme_classic() +
  facet_wrap(~ factor(COHORT, c("Non-Hispanic White", "African American", "Cuban", "Puerto Rican", "Peruvian")))
t_test_cohort_status_pc_horvath_skin_blood <- PCClock_DNAmAge %>%
  group_by(COHORT) %>%
  nest() %>%
  mutate(
    t_test = map(data, ~ t.test(PCHorvath2Resid ~ STATUS, data = .x) %>% broom::tidy())
  ) %>%
  unnest(t_test)

###### AA Grady Trauma Project Comparison
# Load in the metadata using GEOquery package
# This will download the series matrix file directly from GEO
gse <- getGEO("GSE72680", GSEMatrix = TRUE)
gse_data <- gse[[1]] # Access the first data matrix
phenos <- pData(gse_data) # Access the phenotyping data

# Remove gse and gse_data for memory saving
rm(gse)
rm(gse_data)

# Get the ages and sex of the individuals sampled
phenos$age <- as.numeric(gsub(".*age: ([0-9]+).*", "\\1", phenos$characteristics_ch1.1))
phenos$sex <- gsub(".*Sex: ([a-zA-Z]+).*", "\\1", phenos$characteristics_ch1, ignore.case = TRUE)

# Now load the betas
# Path to your downloaded file
beta_file <- "~/Desktop/Capra Lab/Thesis Project/Aim_1/AA_Grady_Trauma_Project/GSE72680_beta_values.txt"

# Read the beta values file using fread
grady_beta_values <- fread(beta_file, header = TRUE, sep = "\t", check.names = FALSE)

# Filter out all the p-value columns, we only want the betas
cols_to_keep <- names(grady_beta_values)[!grepl("Detection PVal", names(grady_beta_values))]
grady_beta_values <- grady_beta_values[, ..cols_to_keep]

# Move V1 to row names, then transpose the betas so that they're in proper formatting for the PC clocks
grady_beta_values <- column_to_rownames(grady_beta_values, var = "V1")
grady_beta_values_transposed <- as.data.frame(t(grady_beta_values))

# Re-format the phenos dataframe
phenos <- phenos %>% rename(Age = age)
phenos <- phenos %>% mutate(Female = case_when(
  sex == "Female" ~ 1,
  sex == "Male" ~ 0
))

# Get the PC Clocks values and the PC Clock Acceleration values
grady_PCClock_DNAmAge <- calcPCClocks(
  path_to_PCClocks_directory = clocksDir,
  datMeth = grady_beta_values_transposed,
  datPheno = phenos
)
grady_PCClock_DNAmAge <- calcPCClocks_Accel(grady_PCClock_DNAmAge)

calculate_metrics <- function(data) {
  metrics <- data %>%
    summarise(
      n = n(),
      MAE = median(abs(PCHorvath1 - Age)),
      MSE = mean((PCHorvath1 - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCHorvath1 - Age))
    ) %>%
    return(metrics)
}
calculate_metrics_hannum <- function(data) {
  metrics <- data %>%
    summarise(
      n = n(),
      MAE = median(abs(PCHannum - Age)),
      MSE = mean((PCHannum - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCHannum - Age))
    ) %>%
    return(metrics)
}

grady_results <- calculate_metrics(grady_PCClock_DNAmAge)
print(grady_results)
grady_hannum_results <- calculate_metrics(grady_PCClock_DNAmAge)
print(grady_hannum_results)

ggplot(data = grady_PCClock_DNAmAge) +
  geom_point(mapping = aes(x = PCHorvath1, y = Age)) +
  theme_classic() +
  geom_smooth(mapping = aes(x = PCHorvath1, y = Age), method = "lm", color = "blue") +
  ggtitle("Grady Project African Americans Correlations (PC Horvath Clock)") +
  xlab("PC Horvath DNAmAge Estimate") +
  ylab("Chronological Age")

grady_PCClock_DNAmAge_older_subset <- grady_PCClock_DNAmAge[grady_PCClock_DNAmAge$Age >= 55, ]
grady_older_results <- calculate_metrics(grady_PCClock_DNAmAge_older_subset)
print(grady_older_results)

grady_hannum_older_results <- calculate_metrics(grady_PCClock_DNAmAge_older_subset)
print(grady_hannum_older_results)
# Calculate the correlation between Horvath and age for the subsetted data
older_correlation <- cor(grady_PCClock_DNAmAge_older_subset$PCHorvath1, grady_PCClock_DNAmAge_older_subset$Age, use = "complete.obs")
hannum_older_correlation <- cor(grady_PCClock_DNAmAge_older_subset$PCHannum, grady_PCClock_DNAmAge_older_subset$Age, use = "complete.obs")

# Print the correlation
print(older_correlation)
print(hannum_older_correlation)

# Calculate MAE for the subset
older_mae <- mean(abs(grady_PCClock_DNAmAge_older_subset$PCHorvath1 - grady_PCClock_DNAmAge_older_subset$Age), na.rm = TRUE)

# Scatter plot with regression line for subsetted data
p_grady <- ggplot(grady_PCClock_DNAmAge_older_subset, aes(x = PCHorvath1, y = Age)) +
  geom_point() + # Scatter plot points
  geom_smooth(method = "lm", color = "black", se = FALSE) + # Regression line
  annotate("text",
    x = -Inf, y = Inf,
    label = paste0("R = ", round(older_correlation, 3)),
    hjust = -0.1, vjust = 1.5, size = 4, fontface = "italic"
  ) +
  theme_classic() +
  ggtitle("Grady Project African Americans Correlations (Age >= 55) (PC Horvath Clock)") +
  xlab("PC Horvath DNAmAge Estimate") +
  ylab("Chronological Age")

# Display and save plot
print(p_grady)
saveRDS(p_grady, "~/Desktop/Capra Lab/Thesis Project/Aim_1/p_grady_older.rds")


older_absolute_errors <- (abs(grady_PCClock_DNAmAge_older_subset$PCHorvath1 - grady_PCClock_DNAmAge_older_subset$Age))
print(median(older_absolute_errors))

# Create a density plot of absolute errors
ggplot(data = data.frame(older_absolute_errors), aes(x = older_absolute_errors)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(
    title = "Density Plot of Absolute Errors",
    x = "Absolute Error",
    y = "Density"
  ) +
  theme_minimal()

###### White Swedish Comparison
gse <- getGEO("GSE87571", GSEMatrix = TRUE)
gse_data <- gse[[1]] # Access the first data matrix
phenos <- pData(gse_data) # Access the phenotyping data

# Remove gse and gse_data for memory saving
rm(gse)
rm(gse_data)

# Get the ages and sex of the individuals sampled
phenos$age <- as.numeric(gsub(".*age: ([0-9]+).*", "\\1", phenos$characteristics_ch1.1))
phenos$sex <- phenos$`gender:ch1`


# Now load the betas
# Path to your downloaded files
beta_file <- "~/Desktop/Capra Lab/Thesis Project/Aim_1/GSE87571/GSE87571_matrix1of2.txt"
beta_file2 <- "~/Desktop/Capra Lab/Thesis Project/Aim_1/GSE87571/GSE87571_matrix2of2.txt"

# Read the beta values files using fread
beta_values <- fread(beta_file, header = TRUE, sep = "\t", check.names = FALSE)
beta_values2 <- fread(beta_file2, header = TRUE, sep = "\t", check.names = FALSE)

# Clean them up by removing the p-val columns
cols_to_keep <- names(beta_values)[!grepl("\\.1", names(beta_values))]
cols_to_keep2 <- names(beta_values2)[!grepl("\\.1", names(beta_values2))]

beta_values <- beta_values[, ..cols_to_keep]
beta_values <- beta_values %>% rename(CpG = ID_REF)

beta_values2 <- beta_values2[, ..cols_to_keep2]
beta_values2 <- beta_values2 %>% rename(CpG = ID_REF)

# Merge the two beta matrices together, and remove the standalones to alleviate memory constraints
merged_beta_values <- merge(beta_values, beta_values2, by = "CpG", all = TRUE)
rm(beta_values)
rm(beta_values2)

#### Re-order the samples so that they match in the betas matrix and metadata (allows the DNAmAge function to calculate age acceleration)
phenos$sample_ID <- gsub("^([^ ]+).*", "\\1", phenos$title)



# Remove columns with NAs, and remove those samples from beta values dataframe
samples_to_remove <- phenos %>%
  filter(is.na(age)) %>%
  pull(sample_ID)

good_columns <- setdiff(colnames(merged_beta_values), samples_to_remove)

merged_beta_values <- subset(merged_beta_values, select = good_columns)

phenos <- phenos %>%
  filter(!is.na(age))


# Step 1: Extract the column names from beta_values (sample IDs)
sample_ids_in_beta <- colnames(merged_beta_values)

# Step 2: Match the sample IDs in phenos$sample_ID to the column order in beta_values
matching_indices <- match(sample_ids_in_beta, phenos$sample_ID)
if (any(is.na(matching_indices))) {
  warning("Some sample IDs in beta_values do not exist in phenos.")
} # We can ignore it cause its just the CpG marker column


# Step 3: Reorder the phenos dataframe based on the matching indices
phenos_reordered <- phenos[na.omit(matching_indices), ]


# Check the result
print(phenos_reordered)

# Rename columns
phenos_reordered <- phenos_reordered %>% rename(Age = age)
phenos_reordered <- phenos_reordered %>% mutate(Female = case_when(
  sex == "Female" ~ 1,
  sex == "Male" ~ 0
))

# Transpose beta values
merged_beta_values <- as.data.frame(merged_beta_values)
merged_beta_values <- merged_beta_values %>% column_to_rownames(var = "CpG")
merged_beta_values_transposed <- as.data.frame(t(merged_beta_values))

# Get the PC Clocks values and the PC Clock Acceleration values
swedish_PCClock_DNAmAge <- calcPCClocks(
  path_to_PCClocks_directory = clocksDir,
  datMeth = merged_beta_values_transposed,
  datPheno = phenos_reordered
)
swedish_PCClock_DNAmAge <- calcPCClocks_Accel(swedish_PCClock_DNAmAge)

calculate_metrics <- function(data) {
  metrics <- data %>%
    summarise(
      n = n(),
      MAE = median(abs(PCHorvath1 - Age)),
      MSE = mean((PCHorvath1 - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCHorvath1 - Age))
    ) %>%
    return(metrics)
}
calculate_metrics_skin_and_blood <- function(data) {
  metrics <- data %>%
    summarise(
      n = n(),
      MAE = median(abs(PCHorvath2 - Age)),
      MSE = mean((PCHorvath2 - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCHorvath2 - Age))
    ) %>%
    return(metrics)
}
calculate_metrics_hannum <- function(data) {
  metrics <- data %>%
    summarise(
      n = n(),
      MAE = median(abs(PCHannum - Age)),
      MSE = mean((PCHannum - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCHannum - Age))
    ) %>%
    return(metrics)
}

swedish_results <- calculate_metrics(swedish_PCClock_DNAmAge)
swedish_results_skin_and_blood <- calculate_metrics_skin_and_blood(swedish_PCClock_DNAmAge)
print(swedish_results)
print(swedish_results_skin_and_blood)

ggplot(data = swedish_PCClock_DNAmAge) +
  geom_point(mapping = aes(x = PCHorvath1, y = Age)) +
  theme_classic() +
  geom_smooth(mapping = aes(x = PCHorvath1, y = Age), method = "lm", color = "blue") +
  ggtitle("Swedish Whites Correlations (PC Horvath Clock)") +
  xlab("PC Horvath DNAmAge Estimate") +
  ylab("Chronological Age")
ggplot(data = swedish_PCClock_DNAmAge) +
  geom_point(mapping = aes(x = PCHorvath2, y = Age)) +
  theme_classic() +
  geom_smooth(mapping = aes(x = PCHorvath2, y = Age), method = "lm", color = "blue") +
  ggtitle("Swedish Whites Correlations (PC Horvath Clock (Skin and Blood))") +
  xlab("PC Horvath DNAmAge Estimate") +
  ylab("Chronological Age")

swedish_PCClock_DNAmAge_older_subset <- swedish_PCClock_DNAmAge[swedish_PCClock_DNAmAge$Age >= 55, ]

# Calculate the correlation between Horvath and age for the subsetted data
older_correlation <- cor(swedish_PCClock_DNAmAge_older_subset$PCHorvath1, swedish_PCClock_DNAmAge_older_subset$Age, use = "complete.obs")
older_correlation_skin_and_blood <- cor(swedish_PCClock_DNAmAge_older_subset$PCHorvath2, swedish_PCClock_DNAmAge_older_subset$Age, use = "complete.obs")
older_correlation_hannum <- cor(swedish_PCClock_DNAmAge_older_subset$PCHannum, swedish_PCClock_DNAmAge_older_subset$Age, use = "complete.obs")

# Print the correlation
print(older_correlation)
print(older_correlation_skin_and_blood)
print(older_correlation_hannum)

# Print results for older Swedish individuals
print(calculate_metrics(swedish_PCClock_DNAmAge_older_subset))
print(calculate_metrics_hannum(swedish_PCClock_DNAmAge_older_subset))

# Calculate MAE for the subset
older_mae <- mean(abs(swedish_PCClock_DNAmAge_older_subset$PCHorvath1 - swedish_PCClock_DNAmAge_older_subset$Age), na.rm = TRUE)

# Scatter plot with regression line for subsetted data
p_swedish <- ggplot(swedish_PCClock_DNAmAge_older_subset, aes(x = PCHorvath1, y = Age)) +
  geom_point() + # Scatter plot points
  geom_smooth(method = "lm", color = "black", se = FALSE) + # Regression line
  annotate("text",
    x = -Inf, y = Inf,
    label = paste0("R = ", round(older_correlation, 3)),
    hjust = -0.1, vjust = 1.5, size = 4, fontface = "italic"
  ) +
  theme_classic() +
  ggtitle("Swedish Whites Correlations (Age >= 55) (PC Horvath Clock)") +
  xlab("PC Horvath DNAmAge Estimate") +
  ylab("Chronological Age")

# Display and save plot
print(p_swedish)
saveRDS(p_swedish, "~/Desktop/Capra Lab/Thesis Project/Aim_1/p_swedish_older.rds")
# Calculate MAE for the subset (Skin and Blood)
older_mae_skin_and_blood <- mean(abs(swedish_PCClock_DNAmAge_older_subset$PCHorvath2 - swedish_PCClock_DNAmAge_older_subset$Age), na.rm = TRUE)

ggplot(swedish_PCClock_DNAmAge_older_subset, aes(x = Age, y = PCHorvath2)) +
  geom_point(color = "blue", alpha = 0.6) + # Scatter plot points
  geom_smooth(method = "lm", color = "red", se = FALSE) + # Regression line
  labs(
    title = "Correlation between PC Horvath DNAmAge (Skin and Blood) and Chronological Age (Age >= 55)",
    x = "Age",
    y = "PC Horvath (Skin and Blood)"
  ) +
  theme_classic() +
  annotate("text",
    x = min(swedish_PCClock_DNAmAge_older_subset$Age, na.rm = TRUE),
    y = max(swedish_PCClock_DNAmAge_older_subset$PCHorvath2, na.rm = TRUE),
    label = paste0("R = ", round(older_correlation_skin_and_blood, 2), "\nMAE = ", round(older_mae_skin_and_blood, 2)),
    hjust = 0, vjust = 1.1, color = "black"
  )

older_absolute_errors <- (abs(swedish_PCClock_DNAmAge_older_subset$PCHorvath1 - swedish_PCClock_DNAmAge_older_subset$Age))
print(median(older_absolute_errors))

older_absolute_errors_skin_and_blood <- (abs(swedish_PCClock_DNAmAge_older_subset$PCHorvath2 - swedish_PCClock_DNAmAge_older_subset$Age))
print(median(older_absolute_errors_skin_and_blood))
# Create a density plot of absolute errors
ggplot(data = data.frame(older_absolute_errors), aes(x = older_absolute_errors)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(
    title = "Density Plot of Absolute Errors",
    x = "Absolute Error",
    y = "Density"
  ) +
  theme_classic()


############# GENOA results
gse <- getGEO("GSE210255", GSEMatrix = TRUE)
gse_data <- gse[[1]] # Access the first data matrix
phenos <- pData(gse_data) # Access the phenotyping data

# Remove gse and gse_data for memory saving
rm(gse)
rm(gse_data)

# Get the ages and sex of the individuals sampled
phenos$age <- as.numeric(gsub(".+:\\s*([0-9.]+).*", "\\1", phenos$characteristics_ch1.2))

# Now load the betas
# Path to your downloaded file
beta_file <- "~/Desktop/Capra Lab/Thesis Project/Aim_1/AA_GENOA_Methylation_Data/GSE210255_Beta_value_EPIC.txt"

# Read the beta values file using fread
beta_values <- fread(beta_file, header = TRUE, sep = " ", check.names = FALSE, quote = "\"")

# Get the betas into the proper format
beta_values <- as.data.frame(beta_values)
beta_values <- beta_values %>% column_to_rownames(var = "CpGnames")
beta_values_transposed <- as.data.frame(t(beta_values))

# Clean up the sample names in phenos
phenos$title <- gsub("sample", "", phenos$title)

# Create Female column (1 for F, 0 for M)
phenos$Female <- ifelse(grepl("F", phenos$characteristics_ch1.1), 1, 0)

# Extract age as integer
phenos$Age <- as.integer(sub(".*: ", "", phenos$characteristics_ch1.2))


# Get the PC Clocks values and the PC Clock Acceleration values
genoa_PCClock_DNAmAge <- calcPCClocks(
  path_to_PCClocks_directory = clocksDir,
  datMeth = beta_values_transposed,
  datPheno = phenos
)
genoa_PCClock_DNAmAge <- calcPCClocks_Accel(genoa_PCClock_DNAmAge)

calculate_metrics <- function(data) {
  metrics <- data %>%
    summarise(
      n = n(),
      MAE = median(abs(PCHorvath1 - Age)),
      MSE = mean((PCHorvath1 - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCHorvath1 - Age))
    ) %>%
    return(metrics)
}
calculate_metrics_hannum <- function(data) {
  metrics <- data %>%
    summarise(
      n = n(),
      MAE = median(abs(PCHannum - Age)),
      MSE = mean((PCHannum - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCHannum - Age))
    ) %>%
    return(metrics)
}
calculate_metrics_skin_and_blood <- function(data) {
  metrics <- data %>%
    summarise(
      n = n(),
      MAE = median(abs(PCHorvath2 - Age)),
      MSE = mean((PCHorvath2 - Age)^2),
      RMSE = sqrt(MSE),
      MAX = max(abs(PCHorvath2 - Age))
    ) %>%
    return(metrics)
}
genoa_results <- calculate_metrics(genoa_PCClock_DNAmAge)
genoa_results_skin_and_blood <- calculate_metrics_skin_and_blood(genoa_PCClock_DNAmAge)
print(genoa_results)
print(genoa_results_skin_and_blood)

cor(genoa_PCClock_DNAmAge$PCHorvath1, genoa_PCClock_DNAmAge$Age)

# Subset data to just the older individuals
genoa_PCClock_DNAmAge_older_subset <- genoa_PCClock_DNAmAge[genoa_PCClock_DNAmAge$Age >= 55, ]
# Calculate the correlation between Horvath and age for the subsetted data
older_correlation <- cor(genoa_PCClock_DNAmAge_older_subset$PCHorvath1, genoa_PCClock_DNAmAge_older_subset$Age, use = "complete.obs")
older_correlation_hannum <- cor(genoa_PCClock_DNAmAge_older_subset$PCHannum, genoa_PCClock_DNAmAge_older_subset$Age, use = "complete.obs")

older_absolute_errors <- (abs(genoa_PCClock_DNAmAge_older_subset$PCHorvath1 - genoa_PCClock_DNAmAge_older_subset$Age))
print(median(older_absolute_errors))

print(calculate_metrics(genoa_PCClock_DNAmAge_older_subset))
print(calculate_metrics_hannum(genoa_PCClock_DNAmAge_older_subset))

# Calculate MAE for the subset
older_mae <- mean(abs(genoa_PCClock_DNAmAge_older_subset$PCHorvath1 - genoa_PCClock_DNAmAge_older_subset$Age), na.rm = TRUE)

# Scatter plot with regression line for older GENOA subset
p_genoa <- ggplot(genoa_PCClock_DNAmAge_older_subset, aes(x = PCHorvath1, y = Age)) +
  geom_point() + # Scatter plot points
  geom_smooth(method = "lm", color = "black", se = FALSE) + # Regression line
  annotate("text",
    x = -Inf, y = Inf,
    label = paste0("R = ", round(older_correlation, 3)),
    hjust = -0.1, vjust = 1.5, size = 4, fontface = "italic"
  ) +
  theme_classic() +
  ggtitle("GENOA African Americans Correlations (Age >= 55) (PC Horvath Clock)") +
  xlab("PC Horvath DNAmAge Estimate") +
  ylab("Chronological Age")

# Display and save plot
print(p_genoa)
saveRDS(p_genoa, "~/Desktop/Capra Lab/Thesis Project/Aim_1/p_genoa_older.rds")

# Code to load the plots individually if needed:
p_grady <- readRDS("~/Desktop/Capra Lab/Thesis Project/Aim_1/p_grady_older.rds")
p_swedish <- readRDS("~/Desktop/Capra Lab/Thesis Project/Aim_1/p_swedish_older.rds")
p_genoa <- readRDS("~/Desktop/Capra Lab/Thesis Project/Aim_1/p_genoa_older.rds")

# Calculate global scales and metrics from the plot data to ensure consistency across all three
all_x <- c(p_grady$data$PCHorvath1, p_swedish$data$PCHorvath1, p_genoa$data$PCHorvath1)
all_y <- c(p_grady$data$Age, p_swedish$data$Age, p_genoa$data$Age)
x_lims <- range(all_x, na.rm = TRUE)
y_lims <- range(all_y, na.rm = TRUE)

# Function to update plot with MAE label (R is already on the plot)
add_labels <- function(p) {
  mae <- mean(abs(p$data$PCHorvath1 - p$data$Age), na.rm = TRUE)
  p + annotate("text",
    x = -Inf, y = Inf,
    label = paste0("MAE = ", round(mae, 2)),
    hjust = -0.1, vjust = 3.0, size = 4, fontface = "italic"
  )
}

p_combined_older <- gridExtra::grid.arrange(
  add_labels(p_swedish) + ggtitle("Swedish Whites") + lims(x = x_lims, y = y_lims),
  add_labels(p_grady) + ggtitle("Grady African Americans") + lims(x = x_lims, y = y_lims),
  add_labels(p_genoa) + ggtitle("GENOA African Americans") + lims(x = x_lims, y = y_lims),
  ncol = 3,
  top = "Older Individuals (Age >= 55) PC Horvath Clock Correlations"
)

# Save combined plot
ggsave("combined_older_correlations.png", p_combined_older, width = 15, height = 5)
saveRDS(p_combined_older, "p_combined_older.rds")
