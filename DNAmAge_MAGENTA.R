# This script applies a series of methylation clocks (implemented within the methylclock R library) to methylation data from the MAGENTA study

# Establish working directory
setwd("/Users/sgonzalez/Desktop/Capra Lab/Thesis Project/Aim_1")

# Load our libraries as needed
library(methylclock)
library(tidyverse)
library(cowplot)
library(readxl)
library(broom)
library(table1)
library(ggpubr)
library(boot)
library(cocor)
library(DunedinPACE)

# Load the beta matrix
normalized_combined <- readRDS("betaMatrices/normalizedBetas/beta_QGCDPB_combined.rds")

# Load the complete sample metadata
sample_metadata <- read_xlsx("ADmethy_pheno.xlsx")

# We have some leftover samples for which there are no methylation data. Let's filter them out
sample_metadata <- subset(sample_metadata, Beta_ID %in% colnames(normalized_combined))

# For poster-making: change the cohort names in the metadata to human-readable. Then make the "Table 1".
sample_metadata_v2 <- sample_metadata %>%
  mutate(COHORT = case_when(
    COHORT == "CuADI" ~ "Cuban",
    COHORT == "NHW" ~ "Non-Hispanic White",
    COHORT == "PERUVIAN" ~ "Peruvian",
    COHORT == "PRADI" ~ "Puerto Rican",
    COHORT == "REAAADI" ~ "African American",
    TRUE ~ COHORT  # Keep the original value if no match is found
  ))
table1(~ ETHNICITY + COHORT + SEX | STATUS, data = sample_metadata)
table1(~ COHORT + SEX | STATUS, data = sample_metadata_v2)
# Modify the rownames to remove everything after the underscore
rownames(normalized_combined) <- gsub("_.*", "", rownames(normalized_combined))

# Get biological age estimates by using the DNAm Age clocks
bio_age_estimates_magenta <- DNAmAge(normalized_combined, cell.count = F, normalize = F)
bio_age_estimates_magenta_age_diff <- DNAmAge(normalized_combined, cell.count = F, age = sample_metadata$AGE_OF_EXAM, normalize = F)

# Join the metadata to the age estimates, using sample IDs as keys
bio_age_estimates_magenta_metadata <- bio_age_estimates_magenta %>% left_join(sample_metadata, join_by("id" == "Beta_ID"))
bio_age_estimates_magenta_age_diff_metadata <- bio_age_estimates_magenta_age_diff %>% left_join(sample_metadata, join_by("id" == "Beta_ID"))

dunedin_preds <- PACEProjector(normalized_combined)
# Get the DunedinPACE age estimates in rows of a data frame
dunedin_preds_df <- as.data.frame(dunedin_preds$DunedinPACE)
# Convert the rownames to a column (sample IDs)
dunedin_preds_df <- rownames_to_column(dunedin_preds_df, var = "ID")

# Join the DunedinPACE predictions to the bigger dataframe
bio_age_estimates_magenta_age_diff_metadata <- bio_age_estimates_magenta_age_diff_metadata %>% left_join(dunedin_preds_df, join_by("id" == "ID"))

# Write the data to a file
#write.csv(bio_age_estimates_magenta_metadata, "bio_age_estimates_magenta_metadata.csv")
#write.csv(bio_age_estimates_magenta_age_diff_metadata, "bio_age_estimates_magenta_age_diff_metadata.csv")

# Plot biological age estimates versus chronological age for the individuals sampled
# First up, the Horvath Clock for control individuals
ggplot(data = subset(bio_age_estimates_magenta_metadata, STATUS == "CONTROL")) +
  geom_point(mapping = aes(x = Horvath, y = AGE_OF_EXAM)) +
  theme_minimal() +
  geom_smooth(mapping = aes(x = Horvath, y = AGE_OF_EXAM), method = "lm", se = FALSE, color = 'blue') +
  xlab('Horvath DNAmAge Estimate') +
  ylab('Chronological Age') +
  ggtitle("Correlation between Horvath DNAmAge and Chrono Age for MAGENTA Controls") +
  annotate("text", x = max(subset(bio_age_estimates_magenta_metadata, STATUS == "CONTROL")$Horvath), y = min(subset(bio_age_estimates_magenta_metadata, STATUS == "CONTROL")$AGE_OF_EXAM), 
  label = paste("R2 =", round(summary(lm(AGE_OF_EXAM ~ Horvath, data = subset(bio_age_estimates_magenta_metadata, STATUS == "CONTROL")))$r.squared, 2)), 
  hjust = 1, vjust = 0)

# And for AD individuals
ggplot(data = subset(bio_age_estimates_magenta_metadata, STATUS == "AD")) +
  geom_point(mapping = aes(x = Horvath, y = AGE_OF_EXAM)) +
  theme_minimal() +
  geom_smooth(mapping = aes(x = Horvath, y = AGE_OF_EXAM), method = "lm", se = FALSE, color = 'blue') +
  xlab('Horvath DNAmAge Estimate') +
  ylab('Chronological Age') +
  ggtitle("Correlation between Horvath DNAmAge and Chrono Age for MAGENTA Cases") +
  annotate("text", x = max(subset(bio_age_estimates_magenta_metadata, STATUS == "AD")$Horvath), y = min(subset(bio_age_estimates_magenta_metadata, STATUS == "AD")$AGE_OF_EXAM), 
           label = paste("R2 =", round(summary(lm(AGE_OF_EXAM ~ Horvath, data = subset(bio_age_estimates_magenta_metadata, STATUS == "AD")))$r.squared, 2)), 
           hjust = 1, vjust = 0)

# Plot correlations for MAGENTA Controls by cohort
ggplot(data = subset(bio_age_estimates_magenta_metadata, STATUS == "CONTROL")) +
  geom_point(mapping = aes(x = Horvath, y = AGE_OF_EXAM)) +
  theme_minimal() +
  geom_smooth(mapping = aes(x = Horvath, y = AGE_OF_EXAM), method = "lm", color = 'blue') +
  ggtitle("MAGENTA Controls Correlations (Horvath Clock)") +
  xlab('Horvath DNAmAge Estimate') +
  ylab('Chronological Age') +
  facet_wrap(~COHORT)

ggplot(data = subset(bio_age_estimates_magenta_metadata, STATUS == "CONTROL" & COHORT == "NHW")) +
  geom_point(mapping = aes(x = Horvath, y = AGE_OF_EXAM)) +
  theme_minimal() +
  geom_smooth(mapping = aes(x = Horvath, y = AGE_OF_EXAM), method = "lm", color = 'blue') +
  #ggtitle("MAGENTA NHW Controls Correlations (Horvath Clock)") +
  xlab('Horvath DNAmAge Estimate') +
  ylab('Chronological Age') +
  xlim(60, 100) +
  ylim(60, 100)
ggplot(data = subset(bio_age_estimates_magenta_metadata, STATUS == "CONTROL" & COHORT == "REAAADI")) +
  geom_point(mapping = aes(x = Horvath, y = AGE_OF_EXAM)) +
  theme_minimal() +
  geom_smooth(mapping = aes(x = Horvath, y = AGE_OF_EXAM), method = "lm", color = 'blue') +
  #ggtitle("MAGENTA AA Controls Correlations (Horvath Clock)") +
  xlab('Horvath DNAmAge Estimate') +
  ylab('Chronological Age')+
  xlim(60, 100) +
  ylim(60, 100)
ggplot(data = subset(bio_age_estimates_magenta_metadata, STATUS == "CONTROL" & COHORT == "PRADI")) +
  geom_point(mapping = aes(x = Horvath, y = AGE_OF_EXAM)) +
  theme_minimal() +
  geom_smooth(mapping = aes(x = Horvath, y = AGE_OF_EXAM), method = "lm", color = 'blue') +
  #ggtitle("MAGENTA PR Controls Correlations (Horvath Clock)") +
  xlab('Horvath DNAmAge Estimate') +
  ylab('Chronological Age')+
  xlim(60, 100) +
  ylim(60, 100)
ggplot(data = subset(bio_age_estimates_magenta_metadata, STATUS == "CONTROL" & COHORT == "CuADI")) +
  geom_point(mapping = aes(x = Horvath, y = AGE_OF_EXAM)) +
  theme_minimal() +
  geom_smooth(mapping = aes(x = Horvath, y = AGE_OF_EXAM), method = "lm", color = 'blue') +
  #ggtitle("MAGENTA CUBAN Controls Correlations (Horvath Clock)") +
  xlab('Horvath DNAmAge Estimate') +
  ylab('Chronological Age')+
  xlim(60, 100) +
  ylim(60, 100)
ggplot(data = subset(bio_age_estimates_magenta_metadata, STATUS == "CONTROL" & COHORT == "PERUVIAN")) +
  geom_point(mapping = aes(x = Horvath, y = AGE_OF_EXAM)) +
  theme_minimal() +
  geom_smooth(mapping = aes(x = Horvath, y = AGE_OF_EXAM), method = "lm", color = 'blue') +
  #ggtitle("MAGENTA PERUVIAN Controls Correlations (Horvath Clock)") +
  xlab('Horvath DNAmAge Estimate') +
  ylab('Chronological Age')+
  xlim(60, 100) +
  ylim(60, 100)

horvath_controls_correlations <- bio_age_estimates_magenta_metadata %>%
  filter(STATUS == "CONTROL") %>%
  group_by(COHORT) %>%
  summarise(correlation = cor(Horvath, AGE_OF_EXAM))

# Plot correlations for MAGENTA Cases by cohort
ggplot(data = subset(bio_age_estimates_magenta_metadata, STATUS == "AD")) +
  geom_point(mapping = aes(x = Horvath, y = AGE_OF_EXAM)) +
  theme_minimal() +
  geom_smooth(mapping = aes(x = Horvath, y = AGE_OF_EXAM), method = "lm", color = 'blue') +
  ggtitle("MAGENTA Cases Correlations (Horvath Clock)") +
  xlab('Horvath DNAmAge Estimate') +
  ylab('Chronological Age') +
  facet_wrap(~COHORT)

horvath_cases_correlations <- bio_age_estimates_magenta_metadata %>%
  filter(STATUS == "AD") %>%
  group_by(COHORT) %>%
  summarise(correlation = cor(Horvath, AGE_OF_EXAM))


# Now Horvath for individuals with AGE_OF_EXAM > 70
ggplot(data = subset(bio_age_estimates_magenta_metadata, AGE_OF_EXAM > 70)) +
  geom_point(mapping = aes(x = Horvath, y = AGE_OF_EXAM)) +
  theme_minimal() +
  geom_smooth(mapping = aes(x = Horvath, y = AGE_OF_EXAM, color = STATUS), method = "lm", se = FALSE, color = 'blue') +
  xlab('Horvath DNAmAge Estimate') +
  ylab('Chronological Age') +
  annotate("text", x = max(subset(bio_age_estimates_magenta_metadata, AGE_OF_EXAM > 70)$Horvath), y = min(subset(bio_age_estimates_magenta_metadata, AGE_OF_EXAM > 70)$AGE_OF_EXAM), 
           label = paste("R2 =", round(summary(lm(AGE_OF_EXAM ~ Horvath, data = subset(bio_age_estimates_magenta_metadata, AGE_OF_EXAM > 70)))$r.squared, 2)), 
           hjust = 1, vjust = 0)
ggplot(data = subset(bio_age_estimates_magenta_metadata, AGE_OF_EXAM > 70)) +
  geom_point(mapping = aes(x = Horvath, y = AGE_OF_EXAM, color = STATUS)) +
  theme_minimal() +
  geom_smooth(mapping = aes(x = Horvath, y = AGE_OF_EXAM, color = STATUS), method = "lm", se = FALSE, color = 'blue') +
  xlab('Horvath DNAmAge Estimate') +
  ylab('Chronological Age')
ggplot(data = subset(bio_age_estimates_magenta_metadata, AGE_OF_EXAM > 70)) +
  geom_point(mapping = aes(x = Horvath, y = AGE_OF_EXAM, shape = COHORT)) +
  theme_minimal() +
  xlab('Horvath DNAmAge Estimate') +
  ylab('Chronological Age')
ggplot(data = subset(bio_age_estimates_magenta_metadata, AGE_OF_EXAM > 70)) +
  geom_point(mapping = aes(x = Horvath, y = AGE_OF_EXAM, shape = COHORT, color = STATUS)) +
  theme_minimal() +
  xlab('Horvath DNAmAge Estimate') +
  ylab('Chronological Age')



# Now PhenoAge
ggplot(data = bio_age_estimates_magenta_metadata) +
  geom_point(mapping = aes(x = Levine, y = AGE_OF_EXAM)) +
  theme_minimal() +
  geom_smooth(mapping = aes(x = Levine, y = AGE_OF_EXAM, color = STATUS), method = "lm", se = FALSE, color = 'blue') +
  xlab('PhenoAge Estimate') +
  ylab('Chronological Age') +
  annotate("text", x = max(bio_age_estimates_magenta_metadata$Levine), y = min(bio_age_estimates_magenta_metadata$AGE_OF_EXAM), 
           label = paste("R2 =", round(summary(lm(AGE_OF_EXAM ~ Levine, data = bio_age_estimates_magenta_metadata))$r.squared, 2)), 
           hjust = 1, vjust = 0)
ggplot(data = bio_age_estimates_magenta_metadata) +
  geom_point(mapping = aes(x = Levine, y = AGE_OF_EXAM, color = STATUS)) +
  theme_minimal() +
  geom_smooth(mapping = aes(x = Levine, y = AGE_OF_EXAM, color = STATUS), method = "lm", se = FALSE, color = 'blue') +
  xlab('PhenoAge Estimate') +
  ylab('Chronological Age')
ggplot(data = bio_age_estimates_magenta_metadata) +
  geom_point(mapping = aes(x = Levine, y = AGE_OF_EXAM, shape = COHORT)) +
  theme_minimal() +
  xlab('PhenoAge Estimate') +
  ylab('Chronological Age')
ggplot(data = bio_age_estimates_magenta_metadata) +
  geom_point(mapping = aes(x = Levine, y = AGE_OF_EXAM, shape = COHORT, color = STATUS)) +
  theme_minimal() +
  xlab('PhenoAge Estimate') +
  ylab('Chronological Age')
# Use the built-in functions from methylclock to make some plots and add stats on there
plotDNAmAge(bio_age_estimates_magenta_metadata$Horvath, bio_age_estimates_magenta_metadata$AGE_OF_EXAM, "Horvath Clock")
plotDNAmAge(bio_age_estimates_magenta_metadata$Levine, bio_age_estimates_magenta_metadata$AGE_OF_EXAM, "PhenoAge Clock")

#### Age Acceleration Analysis ####
# Plot the age acceleration for the Horvath clock, separating by status
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_violin(mapping = aes(x = STATUS, y = ageAcc.Horvath, fill = STATUS), trim = FALSE) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Horvath), width = 0.2, alpha = 0.7) +
  stat_summary(mapping = aes(x = STATUS, y = ageAcc.Horvath), fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  xlab("Disease Status") +
  ylab("Age Acceleration (Horvath DNAm Age - Chronological Age)") + 
  theme_minimal()
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc.Horvath, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Horvath), width = 0.2, alpha = 0.7) +
  xlab("Disease Status") +
  ylab("Age Acceleration (Horvath DNAm Age - Chronological Age)") + 
  theme_minimal() 

ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_violin(mapping = aes(x = STATUS, y = ageAcc.Horvath, fill = STATUS), trim = FALSE) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Age Acceleration (Horvath DNAm Age - Chronological Age)") + 
  theme_minimal()
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_violin(mapping = aes(x = STATUS, y = ageAcc2.Horvath, fill = STATUS), trim = FALSE) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Intrinsic Age Acceleration (Horvath DNAm Age - Chronological Age)") + 
  theme_minimal()

ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc.Horvath, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Age Acceleration (Horvath DNAm Age - Chronological Age)") + 
  theme_minimal()

ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Horvath, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Intrinsic Age Acceleration (Horvath DNAm Age - Chronological Age)") + 
  theme_minimal()

ggplot(data = subset(bio_age_estimates_magenta_age_diff_metadata, COHORT == "NHW")) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Horvath, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Intrinsic Age Acceleration (Horvath DNAm Age - Chronological Age)") +
  ggtitle("Intrinsic Age Acceleration for MAGENTA NHW Individuals (Horvath Clock)") +
  ylim(-10, 20) +
  theme_minimal()
ggplot(data = subset(bio_age_estimates_magenta_age_diff_metadata, COHORT == "REAAADI")) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Horvath, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Intrinsic Age Acceleration (Horvath DNAm Age - Chronological Age)") +
  ggtitle("Intrinsic Age Acceleration for MAGENTA AA Individuals (Horvath Clock)") +
  ylim(-10, 20) +
  theme_minimal()
ggplot(data = subset(bio_age_estimates_magenta_age_diff_metadata, COHORT == "PRADI")) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Horvath, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Intrinsic Age Acceleration (Horvath DNAm Age - Chronological Age)") +
  ggtitle("Intrinsic Age Acceleration for MAGENTA PUR Individuals (Horvath Clock)") +
  ylim(-10, 20) +
  theme_minimal()
ggplot(data = subset(bio_age_estimates_magenta_age_diff_metadata, COHORT == "CuADI")) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Horvath, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Intrinsic Age Acceleration (Horvath DNAm Age - Chronological Age)") +
  ggtitle("Intrinsic Age Acceleration for MAGENTA CUB Individuals (Horvath Clock)") +
  ylim(-10, 20) +
  theme_minimal()
ggplot(data = subset(bio_age_estimates_magenta_age_diff_metadata, COHORT == "PERUVIAN")) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Horvath, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Intrinsic Age Acceleration (Horvath DNAm Age - Chronological Age)") +
  ggtitle("Intrinsic Age Acceleration for MAGENTA PER Individuals (Horvath Clock)") +
  ylim(-10, 20) +
  theme_minimal()

# Perform a t-test to see if the age acceleration is different between the two groups
t.test(bio_age_estimates_magenta_age_diff_metadata$ageAcc.Horvath ~ bio_age_estimates_magenta_age_diff_metadata$STATUS)
t.test(bio_age_estimates_magenta_age_diff_metadata$ageAcc2.Horvath ~ bio_age_estimates_magenta_age_diff_metadata$STATUS)


# Plot the t-test results
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_violin(mapping = aes(x = STATUS, y = ageAcc.Horvath, fill = STATUS), trim = FALSE) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_minimal() +
  geom_segment(aes(x = 1, xend = 2, y = 65, yend = 65), color = "black") +
  geom_text(aes(x = 1.5, y = 67, label = "*"), size = 10) +
  xlab("Disease Status") +
  ylab("Age Acceleration (Horvath DNAm Age - Chronological Age)")
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc.Horvath, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_minimal() +
  geom_segment(aes(x = 1, xend = 2, y = 65, yend = 65), color = "black") +
  geom_text(aes(x = 1.5, y = 67, label = "*"), size = 10) +
  xlab("Disease Status") +
  ylab("Age Acceleration (Horvath DNAm Age - Chronological Age)")


# Now do all of the above but for the PhenoAge clock
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_violin(mapping = aes(x = STATUS, y = ageAcc.Levine, fill = STATUS), trim = FALSE) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Levine), width = 0.2, alpha = 0.7) +
  stat_summary(mapping = aes(x = STATUS, y = ageAcc.Levine), fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  xlab("Disease Status") +
  ylab("Age Acceleration (Levine DNAm Age - Chronological Age)") + 
  theme_minimal()
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc.Levine, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Levine), width = 0.2, alpha = 0.7) +
  xlab("Disease Status") +
  ylab("Age Acceleration (Levine DNAm Age - Chronological Age)") + 
  theme_minimal() 

ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_violin(mapping = aes(x = STATUS, y = ageAcc.Levine, fill = STATUS), trim = FALSE) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Levine), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Age Acceleration (Levine DNAm Age - Chronological Age)") + 
  theme_minimal()

ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Levine, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Levine), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Age Acceleration (Levine DNAm Age - Chronological Age)") + 
  theme_minimal() 

# Perform a t-test to see if the age acceleration is different between the two groups
t.test(bio_age_estimates_magenta_age_diff_metadata$ageAcc.Levine ~ bio_age_estimates_magenta_age_diff_metadata$STATUS)
t.test(bio_age_estimates_magenta_age_diff_metadata$ageAcc2.Levine ~ bio_age_estimates_magenta_age_diff_metadata$STATUS)


# Plot the t-test results
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_violin(mapping = aes(x = STATUS, y = ageAcc.Levine, fill = STATUS), trim = FALSE) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Levine), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_minimal() +
  geom_segment(aes(x = 1, xend = 2, y = 50, yend = 50), color = "black") +
  geom_text(aes(x = 1.5, y = 55, label = "NS"), size = 10) +
  xlab("Disease Status") +
  ylab("Age Acceleration (Levine DNAm Age - Chronological Age)")

ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc.Levine, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Levine), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_minimal() +
  geom_segment(aes(x = 1, xend = 2, y = 50, yend = 50), color = "black") +
  geom_text(aes(x = 1.5, y = 55, label = "NS"), size = 10) +
  xlab("Disease Status") +
  ylab("Age Acceleration (Levine DNAm Age - Chronological Age)")



## Cohort Analysis ##
# Plot the age acceleration for the Horvath clock, separating by cohort
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_violin(mapping = aes(x = COHORT, y = ageAcc.Horvath, fill = COHORT), trim = FALSE) +
  geom_jitter(mapping = aes(x = COHORT, y = ageAcc.Horvath), width = 0.2, alpha = 0.7) +
  stat_summary(mapping = aes(x = COHORT, y = ageAcc.Horvath), fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  xlab("Cohort") +
  ylab("Age Acceleration (Horvath DNAm Age - Chronological Age)") + 
  theme_minimal()

# Plot the age acceleration for the Horvath clock, separating by cohort and status
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc.Horvath, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Cohort") +
  ylab("Age Acceleration (Horvath DNAm Age - Chronological Age)") + 
  theme_minimal() +
  facet_wrap(COHORT ~ .)
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Horvath, fill = STATUS), outlier.shape = NA) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Cohort") +
  ylab("Age Acceleration (Horvath DNAm Age - Chronological Age)") + 
  theme_minimal() +
  facet_wrap(~factor(COHORT, c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")))
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Horvath, fill = STATUS), outlier.shape = NA) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Horvath), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Cohort") +
  ylim(NA, 20) +
  ylab("Intrinsic Age Acceleration (Horvath DNAm Age - Chronological Age)") + 
  theme_minimal() +
  facet_wrap(~factor(COHORT, c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")))


# Plot the age acceleration for the Levine clock, separating by cohort
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_violin(mapping = aes(x = COHORT, y = ageAcc.Levine, fill = COHORT), trim = FALSE) +
  geom_jitter(mapping = aes(x = COHORT, y = ageAcc.Levine), width = 0.2, alpha = 0.7) +
  stat_summary(mapping = aes(x = COHORT, y = ageAcc.Levine), fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  xlab("Cohort") +
  ylab("Age Acceleration (Levine DNAm Age - Chronological Age)") + 
  theme_minimal()
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_violin(mapping = aes(x = COHORT, y = ageAcc2.Levine, fill = COHORT), trim = FALSE) +
  geom_jitter(mapping = aes(x = COHORT, y = ageAcc2.Levine), width = 0.2, alpha = 0.7) +
  stat_summary(mapping = aes(x = COHORT, y = ageAcc2.Levine), fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  xlab("Cohort") +
  ylab("Age Acceleration (Adjusting for Cell Type Composition Effect") + 
  theme_minimal()

# Plot the age acceleration for the Levine clock, separating by cohort and status
cohort_labels_human_readable <- list(
  'NHW' = 'Non-Hispanic Whites',
  'REAAADI' = 'African Americans',
  'CuADI' = 'Cubans',
  'PRADI' = 'Puerto Ricans',
  'PERUVIAN' = 'Peruvians'
)
cohort_labeler <- function(variable,value){
  return(cohort_labels_human_readable[value])
}

ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc.Levine, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Levine), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Age Acceleration (Levine DNAm Age - Chronological Age)") + 
  theme_minimal() +
  facet_wrap(COHORT ~ .)
ggplot(data = bio_age_estimates_magenta_age_diff_metadata) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Levine, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Levine), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylim(NA, 40) +
  ylab("Intrinsic Age Acceleration (Adjusting for Cell Type Composition Effect)") + 
  theme_minimal() +
  facet_wrap(~factor(COHORT, c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")), nrow = 2, labeller = cohort_labeler)

ggplot(data = subset(bio_age_estimates_magenta_age_diff_metadata, COHORT == "NHW")) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Levine, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Levine), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Intrinsic Age Acceleration (PhenoAge DNAm Age - Chronological Age)") +
  ggtitle("Intrinsic Age Acceleration for MAGENTA NHW Individuals (PhenoAge Clock)") +
  theme_minimal()
ggplot(data = subset(bio_age_estimates_magenta_age_diff_metadata, COHORT == "REAAADI")) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Levine, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Levine), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Intrinsic Age Acceleration (PhenoAge DNAm Age - Chronological Age)") +
  ggtitle("Intrinsic Age Acceleration for MAGENTA AA Individuals (PhenoAge Clock)") +
  ylim(-20, 30) +
  theme_minimal()
ggplot(data = subset(bio_age_estimates_magenta_age_diff_metadata, COHORT == "PRADI")) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Levine, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Levine), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Intrinsic Age Acceleration (PhenoAge DNAm Age - Chronological Age)") +
  ggtitle("Intrinsic Age Acceleration for MAGENTA PUR Individuals (PhenoAge Clock)") +
  ylim(-20, 30) +
  theme_minimal()
ggplot(data = subset(bio_age_estimates_magenta_age_diff_metadata, COHORT == "CuADI")) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Levine, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Levine), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Intrinsic Age Acceleration (PhenoAge DNAm Age - Chronological Age)") +
  ggtitle("Intrinsic Age Acceleration for MAGENTA CUB Individuals (PhenoAge Clock)") +
  ylim(-20, 30) +
  theme_minimal()
ggplot(data = subset(bio_age_estimates_magenta_age_diff_metadata, COHORT == "PERUVIAN")) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Levine, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Levine), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Intrinsic Age Acceleration (PhenoAge DNAm Age - Chronological Age)") +
  ggtitle("Intrinsic Age Acceleration for MAGENTA PER Individuals (PhenoAge Clock)") +
  ylim(-20, 30) +
  theme_minimal()



# Perform a t-test to compare age acceleration between AD and Control groups within each cohort
t_test_cohort_status_horvath <- bio_age_estimates_magenta_age_diff_metadata %>%
  group_by(COHORT) %>%
  nest() %>%
  mutate(
    t_test = map(data, ~ t.test(ageAcc.Horvath ~ STATUS, data = .x) %>% tidy())
  ) %>%
  unnest(t_test)
t_test_cohort_status_horvath_2 <- bio_age_estimates_magenta_age_diff_metadata %>%
  group_by(COHORT) %>%
  nest() %>%
  mutate(
    t_test = map(data, ~ t.test(ageAcc2.Horvath ~ STATUS, data = .x) %>% tidy())
  ) %>%
  unnest(t_test)

t_test_cohort_status_phenoage <- bio_age_estimates_magenta_age_diff_metadata %>%
  group_by(COHORT) %>%
  nest() %>%
  mutate(
    t_test = map(data, ~ t.test(ageAcc.Levine ~ STATUS, data = .x) %>% tidy())
  ) %>%
  unnest(t_test)
t_test_cohort_status_phenoage_2 <- bio_age_estimates_magenta_age_diff_metadata %>%
  group_by(COHORT) %>%
  nest() %>%
  mutate(
    t_test = map(data, ~ t.test(ageAcc2.Levine ~ STATUS, data = .x) %>% tidy())
  ) %>%
  unnest(t_test)

### Coefficients Analysis ###
# Plot the coefficients for the Horvath clock
ggplot(data = coefHorvath, aes(x = CpGmarker, y = CoefficientTraining)) +
  geom_point() +
  theme_minimal() +
  xlab("Horvath Clock CpG Marker") +
  ylab("Horvath Clock Coefficient") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

# Plot the coefficients for the PhenoAge clock
ggplot(data = coefLevine, aes(x = CpGmarker, y = CoefficientTraining)) +
  geom_point() +
  theme_minimal() +
  xlab("PhenoAge Clock CpG Marker") +
  ylab("PhenoAge Clock Coefficient") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())  


### Bootstrap Analysis ###
# First define our generalized function to calculate correlations between two variables.
foo <- function(data, indices) {
  dt <- data[indices,]
  c(
    cor(dt[, 2], dt[, 28], method = 'p')
  )
}

# Create some filtered dataframes. We want to keep only the control individuals for each population group.
nhw_controls <- bio_age_estimates_magenta_age_diff_metadata %>% filter(COHORT == "NHW" & STATUS == "CONTROL")
aa_controls <- bio_age_estimates_magenta_age_diff_metadata %>% filter(COHORT == "REAAADI" & STATUS == "CONTROL")
pr_controls <- bio_age_estimates_magenta_age_diff_metadata %>% filter(COHORT == "PRADI" & STATUS == "CONTROL")
cub_controls <- bio_age_estimates_magenta_age_diff_metadata %>% filter(COHORT == "CuADI" & STATUS == "CONTROL")
per_controls <- bio_age_estimates_magenta_age_diff_metadata %>% filter(COHORT == "PERUVIAN" & STATUS == "CONTROL")

# Now leverage the "boot" function to bootstrap  these correlations 1,000 times.
set.seed(42)
nhw_cor_bootstrap <- boot(nhw_controls, foo, R = 1000)
aa_cor_bootstrap <- boot(aa_controls, foo, R = 1000)
pr_cor_bootstrap <- boot(pr_controls, foo, R = 1000)
cub_cor_bootstrap <- boot(cub_controls, foo, R = 1000)
per_cor_bootstrap <- boot(per_controls, foo, R = 1000)


# We'll use the cocor library to compare the correlations between independent groups
nhw_aa_compare <- cocor.indep.groups(r1.jk=nhw_cor_bootstrap$t0, r2.hm=aa_cor_bootstrap$t0,
                   n1=nrow(nhw_controls), n2=nrow(aa_controls),
                   alternative="greater", alpha=0.05, conf.level=0.95, null.value=0)
nhw_pr_compare <- cocor.indep.groups(r1.jk=nhw_cor_bootstrap$t0, r2.hm=pr_cor_bootstrap$t0,
                                     n1=nrow(nhw_controls), n2=nrow(pr_controls),
                                     alternative="greater", alpha=0.05, conf.level=0.95, null.value=0)
nhw_cub_compare <- cocor.indep.groups(r1.jk=nhw_cor_bootstrap$t0, r2.hm=cub_cor_bootstrap$t0,
                                     n1=nrow(nhw_controls), n2=nrow(cub_controls),
                                     alternative="greater", alpha=0.05, conf.level=0.95, null.value=0)
nhw_per_compare <- cocor.indep.groups(r1.jk=nhw_cor_bootstrap$t0, r2.hm=per_cor_bootstrap$t0,
                                     n1=nrow(nhw_controls), n2=nrow(per_controls),
                                     alternative="greater", alpha=0.05, conf.level=0.95, null.value=0)



# Now combine these bootstrapped correlations for pairwise comparisons
nhw_aa <- as.data.frame(cbind(nhw_cor_bootstrap$t, aa_cor_bootstrap$t))
nhw_pr <- as.data.frame(cbind(nhw_cor_bootstrap$t, pr_cor_bootstrap$t))
nhw_cub <- as.data.frame(cbind(nhw_cor_bootstrap$t, cub_cor_bootstrap$t))
nhw_per <- as.data.frame(cbind(nhw_cor_bootstrap$t, per_cor_bootstrap$t))
aa_pr <- as.data.frame(cbind(aa_cor_bootstrap$t, pr_cor_bootstrap$t))
aa_cub <- as.data.frame(cbind(aa_cor_bootstrap$t, cub_cor_bootstrap$t))
aa_per <- as.data.frame(cbind(aa_cor_bootstrap$t, per_cor_bootstrap$t))
pr_cub <- as.data.frame(cbind(pr_cor_bootstrap$t, cub_cor_bootstrap$t))
pr_per <- as.data.frame(cbind(pr_cor_bootstrap$t, per_cor_bootstrap$t))
cub_per <- as.data.frame(cbind(cub_cor_bootstrap$t, per_cor_bootstrap$t))

# Calculate the deltas between correlations for each population group
nhw_aa$delta <- nhw_aa$V1 - nhw_aa$V2
nhw_pr$delta <- nhw_pr$V1 - nhw_pr$V2
nhw_cub$delta <- nhw_cub$V1 - nhw_cub$V2
nhw_per$delta <- nhw_per$V1 - nhw_per$V2
aa_pr$delta <- aa_pr$V1 - aa_pr$V2
aa_cub$delta <- aa_cub$V1 - aa_cub$V2
aa_per$delta <- aa_per$V1 - aa_per$V2
pr_cub$delta <- pr_cub$V1 - pr_cub$V2
pr_per$delta <- pr_per$V1 - pr_per$V2
cub_per$delta <- cub_per$V1 - cub_per$V2

# Plot the density of deltas for all dataframes.
ggplot(data = nhw_aa, aes(x = delta)) +
  geom_density(fill = "blue", alpha = 0.1) +
  theme_minimal() +
  xlab("Delta Correlation (NHW - AA)") +
  ylab("Density")

ggplot(data = nhw_pr, aes(x = delta)) +
  geom_density(fill = "blue", alpha = 0.1) +
  theme_minimal() +
  xlab("Delta Correlation (NHW - PR)") +
  ylab("Density")

ggplot(data = nhw_cub, aes(x = delta)) +
  geom_density(fill = "blue", alpha = 0.1) +
  theme_minimal() +
  xlab("Delta Correlation (NHW - CuB)") +
  ylab("Density")

ggplot(data = nhw_per, aes(x = delta)) +
  geom_density(fill = "blue", alpha = 0.1) +
  theme_minimal() +
  xlab("Delta Correlation (NHW - PERUVIAN)") +
  ylab("Density")

ggplot(data = aa_pr, aes(x = delta)) +
  geom_density(fill = "blue", alpha = 0.1) +
  theme_minimal() +
  xlab("Delta Correlation (AA - PR)") +
  ylab("Density")

ggplot(data = aa_cub, aes(x = delta)) +
  geom_density(fill = "blue", alpha = 0.1) +
  theme_minimal() +
  xlab("Delta Correlation (AA - CuB)") +
  ylab("Density")

ggplot(data = aa_per, aes(x = delta)) +
  geom_density(fill = "blue", alpha = 0.1) +
  theme_minimal() +
  xlab("Delta Correlation (AA - PERUVIAN)") +
  ylab("Density")

ggplot(data = pr_cub, aes(x = delta)) +
  geom_density(fill = "blue", alpha = 0.1) +
  theme_minimal() +
  xlab("Delta Correlation (PR - CuB)") +
  ylab("Density")

ggplot(data = pr_per, aes(x = delta)) +
  geom_density(fill = "blue", alpha = 0.1) +
  theme_minimal() +
  xlab("Delta Correlation (PR - PERUVIAN)") +
  ylab("Density")

ggplot(data = cub_per, aes(x = delta)) +
  geom_density(fill = "blue", alpha = 0.1) +
  theme_minimal() +
  xlab("Delta Correlation (CuB - PERUVIAN)") +
  ylab("Density")



