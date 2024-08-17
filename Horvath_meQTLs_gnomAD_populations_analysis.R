# This is a script to analyze the meQTLs (from the three aforementioned sets of meQTLs) that influence clock CpGs in gnomAD (4.1) populations

# Load the necessary libraries
library(tidyverse)
library(ggpubr)

# Read the file in
horvath_meqtls_gnomAD <- read.table("clock_cpg_meqtls_hg38_gnomAD.txt")

# Get the AFs for each population. We'll match on the string "AF_pop" and then extract the value
# We'll use a function to do this. The function will take in the info string and return a named vector of allele frequencies
extract_AF <- function(info_string) {
  # Modify the regular expression to capture numbers in scientific notation
  af_matches <- str_match_all(info_string, "AF_(\\w+)=([0-9]+(?:\\.[0-9]+)?(?:e[+-]?[0-9]+)?)")[[1]]
  setNames(as.numeric(af_matches[, 3]), af_matches[, 2])
}
# Write another function that gets the global AF value for gnomAD
extract_global_AF <- function(info_string) {
  pattern <- "AF=(-?\\d+(\\.\\d+)?(e[+-]?\\d+)?)"
  match <- regmatches(info_string, regexpr(pattern, info_string, perl = TRUE))
  af_value <- sub("AF=", "", match)
  return(as.numeric(af_value))
}

# Apply the function to a dataframe column
extract_AF_column <- function(df, column_name) {
  df$new_AF <- sapply(df[[column_name]], extract_global_AF)
  return(df)
}
horvath_meqtls_gnomAD <- extract_AF_column(horvath_meqtls_gnomAD, "V8")

# We apply the function to the V8 column. We'll use rowwise to apply the function to each row.
horvath_meqtls_gnomAD <- horvath_meqtls_gnomAD %>%
  rowwise() %>%
  mutate(af_data = list(extract_AF(V8))) %>%
  unnest_wider(af_data)

# Make a filtered version of the dataframe, selecting columns V1-V5, and the columns for specific populations
horvath_meqtls_gnomAD_filtered <- horvath_meqtls_gnomAD %>%
  select(V1:V5, afr, ami, amr, asj, eas, fin, mid, nfe, sas, remaining)

# Plot the allele frequencies for the major global populations
horvath_meqtls_gnomAD_filtered %>%
  gather(Ancestry, AF, afr:sas) %>%
  ggplot(aes(x = AF, fill = Ancestry)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  labs(x = "Allele Frequency", y = "Number of meQTLs") +
  theme_minimal()

# for ggpubr version
# Transform the data
horvath_meqtls_gnomAD_filtered_long <- horvath_meqtls_gnomAD_filtered %>%
  gather(Ancestry, AF, afr:sas)

# Create the plot
p7 <- gghistogram(
  data = horvath_meqtls_gnomAD_filtered_long,
  x = "AF",
  fill = "Ancestry",
  color = "black",
  position = "identity",
  alpha = 0.5,
  bins = 50,
  xlab = "Allele Frequency",
  ylab = "Number of meQTLs",
)

# Plot the allele frequencies for the major global populations as boxplots
horvath_meqtls_gnomAD_filtered %>%
  gather(Ancestry, AF, afr:sas) %>%
  ggplot(aes(x = Ancestry, y = AF, fill = Ancestry)) +
  geom_boxplot() +
  labs(x = "Ancestry", y = "Allele Frequency") +
  theme_minimal()

# ggpubr version
p8 <- ggboxplot(
  data = horvath_meqtls_gnomAD_filtered_long,
  x = "Ancestry",
  y = "AF",
  xlab = "Ancestry",
  ylab = "Allele Frequency",
  outlier.shape = NA,
  legend = "none"
)

p8_violin <- ggviolin(
  data = horvath_meqtls_gnomAD_filtered_long,
  x = "Ancestry",
  y = "AF",
  xlab = "Ancestry",
  ylab = "Allele Frequency",
  legend = "none" 
)


# Extract the AF column as a vector
af_values <- horvath_meqtls_gnomAD$new_AF

# Create a data frame with AF values and an index for plotting
df_af <- data.frame(Index = seq_along(af_values), AF = af_values)

df_af <- df_af %>% arrange(desc(AF))


# Create a barplot of the AF values
gghistogram(df_af, x = "AF", 
                     fill = "lightblue", color = "black", 
                     x.text.angle = 0, ylab = "Number of meQTLs",
                     xlab = "Allele Frequency")

###### Mann Whitney U test to compare meQTL AF in AFR vs all other populations
comparison_columns <- colnames(horvath_meqtls_gnomAD_filtered)
comparison_columns <- comparison_columns[-c(1,2,3,4,5)] # Remove all columns that aren't a population AF

# Create a list to store the results
test_results <- list()

# Loop through each of the other columns and perform the Mann-Whitney U test
for (col in comparison_columns) {
  # Perform the Mann-Whitney U test
  test <- wilcox.test(horvath_meqtls_gnomAD_filtered$afr, horvath_meqtls_gnomAD_filtered[[col]], alternative = "two.sided")
  
  # Store the results
  test_results[[col]] <- list(
    statistic = test$statistic,
    p_value = test$p.value,
    alternative = test$alternative,
    method = test$method
  )
}

# Convert results to a data frame for easier viewing
test_results_df <- do.call(rbind, lapply(names(test_results), function(x) {
  cbind(Column = x, test_results[[x]])
}))

# View the results
test_results_df

# How many variants are found in just AFR pops?
afr_index <- which(names(horvath_meqtls_gnomAD_filtered) == "afr")

# Find rows where "afr" is non-zero and all subsequent columns are zero
result <- horvath_meqtls_gnomAD_filtered[ horvath_meqtls_gnomAD_filtered$afr != 0 & rowSums(horvath_meqtls_gnomAD_filtered[ , (afr_index + 1):ncol(horvath_meqtls_gnomAD_filtered)]) == 0, ]

# Display the result
print(result)
View(result)
