# This is a script to analyze meQTLs that influence clock CpGs in gnomAD admixed (Latino) populations

# Set the working directory
setwd("/Users/sgonzalez/Desktop/Capra Lab/Thesis Project/Aim_1/")

# Load some libraries
library(tidyverse)

# Read the intersect file in: first columns refer to variant
horvath_meqtls_gnomAD_admixed <- read.table("horvath_cpg_meQTL_gnomAD_ancestry_intersected.bed")

# Get the AFs for the 3 local ancestry blocks
horvath_meqtls_gnomAD_admixed$AFR_AF <- as.numeric(str_split_i(str_split_i(horvath_meqtls_gnomAD_admixed$V22, ";", i = 3), "=", i = 2))
horvath_meqtls_gnomAD_admixed$AMR_AF <- as.numeric(str_split_i(str_split_i(horvath_meqtls_gnomAD_admixed$V22, ";", i = 6), "=", i = 2))
horvath_meqtls_gnomAD_admixed$EUR_AF <- as.numeric(str_split_i(str_split_i(horvath_meqtls_gnomAD_admixed$V22, ";", i = 9), "=", i = 2))

# Make a plot of the allele frequencies for the 3 local ancestry blocks
horvath_meqtls_gnomAD_admixed %>%
  gather(Ancestry, AF, AFR_AF:EUR_AF) %>%
  ggplot(aes(x = AF, fill = Ancestry)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(x = "Allele Frequency", y = "Number of meQTLs") +
  theme_minimal()


# Make a violin plot of the allele frequencies for the 3 local ancestry blocks
horvath_meqtls_gnomAD_admixed %>%
  gather(Ancestry, AF, AFR_AF:EUR_AF) %>%
  ggplot(aes(x = Ancestry, y = AF, fill = Ancestry)) +
  geom_boxplot() +
  geom_point(position = position_jitter(seed = 1, width = 0.2)) +
  labs(x = "Ancestry", y = "Allele Frequency") +
  theme_minimal()

# Make a violin plot of the allele frequencies for the 3 local ancestry blocks, but connect the points going across the ancestry blocks
# This is messy, but we'll refine this idea shortly
horvath_meqtls_gnomAD_admixed %>%
  gather(Ancestry, AF, AFR_AF:EUR_AF) %>%
  ggplot(aes(x = Ancestry, y = AF, fill = Ancestry)) +
  geom_violin() +
  geom_point() +
  geom_line(mapping = aes(group = V22),
            position = position_dodge(0.1)) +
  labs(x = "Ancestry", y = "Allele Frequency") +
  theme_minimal()

# Calculate pairwise differences for each variant per ancestry block
pairwise_diff <- horvath_meqtls_gnomAD_admixed %>%
  mutate(AFR_AMR_diff = AFR_AF - AMR_AF,
         AFR_EUR_diff = AFR_AF - EUR_AF,
         AMR_EUR_diff = AMR_AF - EUR_AF)

# Create a barplot to represent the overall difference in allele frequencies for each variant
pairwise_diff %>%
  gather(Difference, Value, AFR_AMR_diff:AMR_EUR_diff) %>%
  ggplot(aes(x = Difference, y = Value, fill = Difference)) +
  geom_boxplot() +
  geom_jitter() +
  labs(x = "Pairwise Difference", y = "Overall Difference in Allele Frequencies") +
  theme_minimal()

pairwise_diff %>%
  gather(Difference, Value, AFR_AMR_diff:AMR_EUR_diff) %>%
  ggplot(aes(x = Difference, y = Value, fill = Difference)) +
  geom_boxplot() +
  geom_line(mapping = aes(group = V22),
            position = position_dodge(0.1)) +
  labs(x = "Pairwise Difference", y = "Overall Difference in Allele Frequencies") +
  theme_minimal()

# Plot the number of variants that, for each ancestry, are not fixed (AF = 1) or not "ghosts" entirely (AF = 0)
horvath_meqtls_gnomAD_admixed %>%
  mutate(AFR_not_fixed = ifelse(AFR_AF != 1 & AFR_AF != 0, 1, 0),
         AMR_not_fixed = ifelse(AMR_AF != 1 & AMR_AF != 0, 1, 0),
         EUR_not_fixed = ifelse(EUR_AF != 1 & EUR_AF != 0, 1, 0)) %>%
  gather(Ancestry, Not_Fixed, AFR_not_fixed:EUR_not_fixed) %>%
  group_by(Ancestry) %>%
  summarize(Proportion = sum(Not_Fixed) / n()) %>%
  ggplot(aes(x = Ancestry, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity") +
  labs(x = "Ancestry", y = "Proportion of Variants") +
  theme_minimal()

# Make a density plot to show the frequencies of each meQTL in each ancestry
horvath_meqtls_gnomAD_admixed %>%
  gather(Ancestry, AF, AFR_AF:EUR_AF) %>%
  ggplot(aes(x = AF, fill = Ancestry)) +
  geom_density(alpha = 0.5) +
  labs(x = "Allele Frequency", y = "Density") +
  theme_minimal()

######## Start of the analysis of the three sets of meQTLs that influence clock CpGs ########
# The sets of meQTLs are from Hawe et al, GENOA, and EPIGEN (same as the above set for this last one)

# Read the file in
clock_cpg_three_sets_meqtls_gnomAD_local_ancestry <- read.table("clock_cpg_meqtls_hg38_gnomAD_local_ancestry.txt")

# Get the AFs for the 3 local ancestry blocks
clock_cpg_three_sets_meqtls_gnomAD_local_ancestry$AFR_AF <- as.numeric(str_split_i(str_split_i(clock_cpg_three_sets_meqtls_gnomAD_local_ancestry$V8, ";", i = 3), "=", i = 2))
clock_cpg_three_sets_meqtls_gnomAD_local_ancestry$AMR_AF <- as.numeric(str_split_i(str_split_i(clock_cpg_three_sets_meqtls_gnomAD_local_ancestry$V8, ";", i = 6), "=", i = 2))
clock_cpg_three_sets_meqtls_gnomAD_local_ancestry$EUR_AF <- as.numeric(str_split_i(str_split_i(clock_cpg_three_sets_meqtls_gnomAD_local_ancestry$V8, ";", i = 9), "=", i = 2))

# Make a plot of the allele frequencies for the 3 local ancestry blocks
clock_cpg_three_sets_meqtls_gnomAD_local_ancestry %>%
  gather(Ancestry, AF, AFR_AF:EUR_AF) %>%
  ggplot(aes(x = AF, fill = Ancestry)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(x = "Allele Frequency", y = "Number of meQTLs") +
  theme_minimal()

# For ggpubr version
# Transform the data
clock_cpg_three_sets_meqtls_gnomAD_local_ancestry_long <- clock_cpg_three_sets_meqtls_gnomAD_local_ancestry %>%
  gather(Ancestry, AF, AFR_AF:EUR_AF)

# Create the plot
p9 <- gghistogram(
  data = clock_cpg_three_sets_meqtls_gnomAD_local_ancestry_long,
  x = "AF",
  fill = "Ancestry",
  color = "black",
  alpha = 0.5,
  bins = 50,
  xlab = "Allele Frequency",
  ylab = "Number of meQTLs",
)

# Make a box plot of the allele frequencies for the 3 local ancestry blocks
clock_cpg_three_sets_meqtls_gnomAD_local_ancestry %>%
  gather(Ancestry, AF, AFR_AF:EUR_AF) %>%
  ggplot(aes(x = Ancestry, y = AF, fill = Ancestry)) +
  geom_boxplot() +
  labs(x = "Ancestry", y = "Allele Frequency") +
  theme_minimal()

# ggpubr version
p10 <- ggboxplot(
  data = clock_cpg_three_sets_meqtls_gnomAD_local_ancestry_long,
  x = "Ancestry",
  y = "AF",
  fill = "Ancestry",
  xlab = "Ancestry",
  ylab = "Allele Frequency",
)

# Plot the pairwise differences for each variant per ancestry block
pairwise_diff_clock_cpg <- clock_cpg_three_sets_meqtls_gnomAD_local_ancestry %>%
  mutate(AFR_AMR_diff = AFR_AF - AMR_AF,
         AFR_EUR_diff = AFR_AF - EUR_AF,
         AMR_EUR_diff = AMR_AF - EUR_AF)

pairwise_diff_clock_cpg %>%
  gather(Difference, Value, AFR_AMR_diff:AMR_EUR_diff) %>%
  ggplot(aes(x = Difference, y = Value, fill = Difference)) +
  geom_boxplot() +
  labs(x = "Pairwise Difference", y = "Overall Difference in Allele Frequencies") +
  theme_minimal()

# For ggpubr version
# Transform the data
pairwise_diff_clock_cpg_long <- pairwise_diff_clock_cpg %>%
  gather(Difference, Value, AFR_AMR_diff:AMR_EUR_diff)

# Create the plot
p11 <- ggboxplot(
  data = pairwise_diff_clock_cpg_long,
  x = "Difference",
  y = "Value",
  fill = "Difference",
  xlab = "Pairwise Difference",
  ylab = "Overall Difference in Allele Frequencies",
  outlier.shape = NA
) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") 

# Plot the number of variants that, for each ancestry, are not fixed (AF = 1) or not "ghosts" entirely (AF = 0)
clock_cpg_three_sets_meqtls_gnomAD_local_ancestry %>%
  mutate(AFR_not_fixed = ifelse(AFR_AF != 1 & AFR_AF != 0, 1, 0),
         AMR_not_fixed = ifelse(AMR_AF != 1 & AMR_AF != 0, 1, 0),
         EUR_not_fixed = ifelse(EUR_AF != 1 & EUR_AF != 0, 1, 0)) %>%
  gather(Ancestry, Not_Fixed, AFR_not_fixed:EUR_not_fixed) %>%
  group_by(Ancestry) %>%
  summarize(Proportion = sum(Not_Fixed) / n()) %>%
  ggplot(aes(x = Ancestry, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity") +
  labs(x = "Ancestry", y = "Proportion of Variants") +
  theme_minimal()

