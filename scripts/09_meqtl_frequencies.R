## This is a script to analyze the clock CpG meQTL from three cohorts that are found in MAGENTA individuals
# Set the working directory
setwd("/Users/sgonzalez/Desktop/Capra Lab/Thesis Project/Aim_1/")

# Load some libraries
library(tidyverse)
library(vcfR)
library(readxl)
library(methylclock) # For clock coefficients


# Load the files
aa_meqtl <- read.vcfR("AA_clock_meqtl_ancestry.vcf.gz")
hisp_meqtl <- read.vcfR("HISP_clock_meqtl_ancestry.vcf.gz")

# Read sample IDs from VCF file
hisp_samples <- read.table("HISPANIC_sample_IDs.txt", header = FALSE, stringsAsFactors = FALSE)$V1
aa_samples <- read.table("AA_sample_IDs.txt", header = FALSE, stringsAsFactors = FALSE)$V1

# Convert the VCF FIXED fields into a dataframe (gets us the variants)
aa_clock_variants <- as.data.frame(getFIX(aa_meqtl))
hisp_clock_variants <- as.data.frame(getFIX(hisp_meqtl))

# Extract the INFO column (gets us the affected clock CpG and the beta for the meQTL)
aa_info_df <- aa_meqtl@fix[, "INFO"] %>%
    str_split_fixed(";", 2) %>%
    as_tibble() %>%
    separate(V1, into = c("X", "PROBE"), sep = "=") %>%
    separate(V2, into = c("X2", "BETA"), sep = "=") %>%
    select(PROBE, BETA)
hisp_info_df <- hisp_meqtl@fix[, "INFO"] %>%
    str_split_fixed(";", 2) %>%
    as_tibble() %>%
    separate(V1, into = c("X", "PROBE"), sep = "=") %>%
    separate(V2, into = c("X2", "BETA"), sep = "=") %>%
    select(PROBE, BETA)

aa_info <- aa_meqtl@fix[, "INFO"]
hisp_info <- hisp_meqtl@fix[, "INFO"]

# Convert genotype fields into a data frame
aa_geno_df <- extract.gt(aa_meqtl, element = "GT") # Extract GT (genotype)
aa_an1_df <- extract.gt(aa_meqtl, element = "AN1") # Extract AN1 (ancestry haplotype 1)
aa_an2_df <- extract.gt(aa_meqtl, element = "AN2") # Extract AN2 (ancestry haplotype 2)

hisp_geno_df <- extract.gt(hisp_meqtl, element = "GT")
hisp_an1_df <- extract.gt(hisp_meqtl, element = "AN1")
hisp_an2_df <- extract.gt(hisp_meqtl, element = "AN2")

# Combine genotype and ancestry
aa_geno_long <- pivot_longer(as_tibble(aa_geno_df, rownames = "Variant"),
    cols = -Variant,
    names_to = "Sample",
    values_to = "Genotype"
)
hisp_geno_long <- pivot_longer(as_tibble(hisp_geno_df, rownames = "Variant"),
    cols = -Variant,
    names_to = "Sample",
    values_to = "Genotype"
)

aa_an1_long <- pivot_longer(as_tibble(aa_an1_df, rownames = "Variant"),
    cols = -Variant,
    names_to = "Sample",
    values_to = "Ancestry1"
)
aa_an2_long <- pivot_longer(as_tibble(aa_an2_df, rownames = "Variant"),
    cols = -Variant,
    names_to = "Sample",
    values_to = "Ancestry2"
)

hisp_an1_long <- pivot_longer(as_tibble(hisp_an1_df, rownames = "Variant"),
    cols = -Variant,
    names_to = "Sample",
    values_to = "Ancestry1"
)

hisp_an2_long <- pivot_longer(as_tibble(hisp_an2_df, rownames = "Variant"),
    cols = -Variant,
    names_to = "Sample",
    values_to = "Ancestry2"
)

# Merge genotype and ancestry
aa_geno_full <- reduce(list(aa_geno_long, aa_an1_long, aa_an2_long), left_join, by = c("Variant", "Sample"))
hisp_geno_full <- reduce(list(hisp_geno_long, hisp_an1_long, hisp_an2_long), left_join, by = c("Variant", "Sample"))

# Filter to keep just the individuals in MAGENTA for which we have methylation data
aa_geno_full <- aa_geno_full %>% filter(Sample %in% aa_samples)
hisp_geno_full <- hisp_geno_full %>% filter(Sample %in% hisp_samples)

# Get the Puerto Ricans
magenta_metadata <- read_xlsx("ADmethy_pheno.xlsx")
magenta_puerto_ricans <- magenta_metadata %>% filter(COHORT == "PRADI")

pr_geno_full <- hisp_geno_full %>% filter(Sample %in% magenta_puerto_ricans$CGI)

# Extract variant info from aa_meqtl@fix
aa_variant_info <- data.frame(
    Variant = aa_meqtl@fix[, "ID"],
    INFO = as.character(aa_meqtl@fix[, "INFO"])
)

# Extract PROBE and BETA from INFO
aa_variant_info <- aa_variant_info %>%
    mutate(
        PROBE = str_extract(INFO, "PROBE=[^;]+") %>% str_remove("PROBE="),
        BETA = str_extract(INFO, "BETA=[-0-9.]+") %>% str_remove("BETA=") %>% as.numeric(),
        Variant = paste0(Variant, "_", row_number()) # Append row index
    ) %>%
    select(Variant, PROBE, BETA) # Keep only relevant columns

# Merge using ID as the key
aa_geno_full <- aa_geno_full %>%
    left_join(aa_variant_info, by = "Variant")

# Do all the same for the Hispanics
hisp_variant_info <- data.frame(
    Variant = hisp_meqtl@fix[, "ID"],
    INFO = as.character(hisp_meqtl@fix[, "INFO"])
)

# Extract PROBE and BETA from INFO
hisp_variant_info <- hisp_variant_info %>%
    mutate(
        PROBE = str_extract(INFO, "PROBE=[^;]+") %>% str_remove("PROBE="),
        BETA = str_extract(INFO, "BETA=[-0-9.]+") %>% str_remove("BETA=") %>% as.numeric(),
        Variant = paste0(Variant, "_", row_number()) # Append row index
    ) %>%
    select(Variant, PROBE, BETA) # Keep only relevant columns

# Merge  using ID as the key
hisp_geno_full <- hisp_geno_full %>%
    left_join(hisp_variant_info, by = "Variant")

pr_geno_full <- pr_geno_full %>%
    left_join(hisp_variant_info, by = "Variant")


aa_variant_counts <- aa_geno_full %>%
    group_by(PROBE) %>%
    summarise(unique_variants = n_distinct(Variant), .groups = "drop") %>%
    arrange(desc(unique_variants))
hisp_variant_counts <- hisp_geno_full %>%
    group_by(PROBE) %>%
    summarise(unique_variants = n_distinct(Variant), .groups = "drop") %>%
    arrange(desc(unique_variants))

ggplot(aa_variant_counts, aes(x = unique_variants)) +
    geom_histogram(bins = 50, fill = "blue", alpha = 0.7, color = "black") +
    labs(
        x = "Number of Unique Variants",
        y = "Clock CpGs Affected"
    ) +
    theme_classic()

ggplot(aa_variant_counts, aes(x = reorder(PROBE, unique_variants), y = unique_variants)) +
    coord_flip() +
    geom_col(fill = "steelblue", alpha = 0.7) +
    labs(
        x = "Horvath Clock CpG Sites (240/353)",
        y = "Number of Unique Variants Affecting Horvath Clock CpGs"
    ) +
    theme_classic() +
    theme(axis.text.y = element_blank())



########### PLOTS
###### African Americans
# Plot the number of homozygous and heterozygous individuals
aa_geno_full <- aa_geno_full %>%
    mutate(ALT_Copies = case_when(
        Genotype == "0|0" ~ 0,
        Genotype == "0|1" | Genotype == "1|0" ~ 1,
        Genotype == "1|1" ~ 2,
        TRUE ~ NA_real_
    )) %>%
    filter(!is.na(ALT_Copies)) # Remove NA values

# Count individuals per ALT copy number
aa_genotype_counts <- aa_geno_full %>%
    group_by(ALT_Copies) %>%
    summarise(Count = n())

# Plot as a bar chart
ggplot(aa_genotype_counts, aes(x = as.factor(ALT_Copies), y = Count, fill = as.factor(ALT_Copies))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("blue", "orange", "red")) +
    labs(x = "Number of ALT Alleles", y = "Number of Variants", fill = "ALT Copies") +
    ggtitle("Zygosity of clock meQTL in African Americans") +
    theme_minimal()

# Reshape the dataframe to long format for ancestries
aa_geno_full_long <- aa_geno_full %>%
    pivot_longer(cols = c(Ancestry1, Ancestry2), names_to = "Ancestry_Type", values_to = "Ancestry") %>%
    mutate(Ancestry = factor(Ancestry, levels = c(0, 1, 2), labels = c("AMR", "EUR", "AFR")))
# Stacked bar plot
ggplot(aa_geno_full_long, aes(x = ALT_Copies, fill = Ancestry)) +
    geom_bar(position = "stack") +
    labs(x = "Number of ALT Alleles", y = "Number of Variants", fill = "Ancestry") +
    scale_fill_manual(values = c("AMR" = "blue", "EUR" = "red", "AFR" = "green")) +
    ggtitle("Ancestry Distribution for clock meQTL (African Americans)") +
    theme_minimal()





###### Hispanics
hisp_geno_full <- hisp_geno_full %>%
    mutate(ALT_Copies = case_when(
        Genotype == "0|0" ~ 0,
        Genotype == "0|1" | Genotype == "1|0" ~ 1,
        Genotype == "1|1" ~ 2,
        TRUE ~ NA_real_
    )) %>%
    filter(!is.na(ALT_Copies)) # Remove NA values

# Count individuals per ALT copy number
hisp_genotype_counts <- hisp_geno_full %>%
    group_by(ALT_Copies) %>%
    summarise(Count = n())

# Plot as a bar chart
ggplot(hisp_genotype_counts, aes(x = as.factor(ALT_Copies), y = Count, fill = as.factor(ALT_Copies))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("blue", "orange", "red")) +
    labs(x = "Number of ALT Alleles", y = "Number of Variants", fill = "ALT Copies") +
    ggtitle("Zygosity of clock meQTL in Hispanics") +
    theme_minimal()

# Reshape the dataframe to long format for ancestries
hisp_geno_full_long <- hisp_geno_full %>%
    pivot_longer(cols = c(Ancestry1, Ancestry2), names_to = "Ancestry_Type", values_to = "Ancestry") %>%
    mutate(Ancestry = factor(Ancestry, levels = c(0, 1, 2), labels = c("AMR", "EUR", "AFR")))
# Stacked bar plot
ggplot(hisp_geno_full_long, aes(x = ALT_Copies, fill = Ancestry)) +
    geom_bar(position = "stack") +
    labs(x = "Number of ALT Alleles", y = "Number of Variants", fill = "Ancestry") +
    scale_fill_manual(values = c("AMR" = "blue", "EUR" = "red", "AFR" = "green")) +
    ggtitle("Ancestry Distribution for clock meQTL (Hispanics)") +
    theme_minimal()


####### Puerto Ricans
pr_geno_full <- pr_geno_full %>%
    mutate(ALT_Copies = case_when(
        Genotype == "0|0" ~ 0,
        Genotype == "0|1" | Genotype == "1|0" ~ 1,
        Genotype == "1|1" ~ 2,
        TRUE ~ NA_real_
    )) %>%
    filter(!is.na(ALT_Copies)) # Remove NA values

# Count individuals per ALT copy number
pr_genotype_counts <- pr_geno_full %>%
    group_by(ALT_Copies) %>%
    summarise(Count = n())

# Plot as a bar chart
ggplot(pr_genotype_counts, aes(x = as.factor(ALT_Copies), y = Count, fill = as.factor(ALT_Copies))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("blue", "orange", "red")) +
    labs(x = "Number of ALT Alleles", y = "Number of Variants", fill = "ALT Copies") +
    ggtitle("Zygosity of clock meQTL in Puerto Ricans") +
    theme_minimal()

# Reshape the dataframe to long format for ancestries
pr_geno_full_long <- pr_geno_full %>%
    pivot_longer(cols = c(Ancestry1, Ancestry2), names_to = "Ancestry_Type", values_to = "Ancestry") %>%
    mutate(Ancestry = factor(Ancestry, levels = c(0, 1, 2), labels = c("AMR", "EUR", "AFR")))
# Stacked bar plot
ggplot(pr_geno_full_long, aes(x = ALT_Copies, fill = Ancestry)) +
    geom_bar(position = "stack") +
    labs(x = "Number of ALT Alleles", y = "Number of Variants", fill = "Ancestry") +
    scale_fill_manual(values = c("AMR" = "blue", "EUR" = "red", "AFR" = "green")) +
    ggtitle("Ancestry Distribution for clock meQTL (Puerto Ricans)") +
    theme_minimal()



######### Beta matrix analysis
# Load the combined beta matrices
betas <- readRDS("betaMatrices/normalizedBetas/beta_QGCDPB_combined.rds")
# Modify the rownames to remove everything after the underscore
rownames(betas) <- gsub("_.*", "", rownames(betas))


# Extract CpGs affected by variants
aa_affected_cpgs <- unique(aa_geno_full_long$PROBE)

# Subset betas dataframe to keep only affected CpGs
aa_betas_filtered <- betas[rownames(betas) %in% aa_affected_cpgs, ]

# Convert betas to long format
aa_betas_long <- aa_betas_filtered %>%
    as.data.frame() %>%
    rownames_to_column(var = "CpG") %>%
    pivot_longer(-CpG, names_to = "Beta_ID", values_to = "Methylation_Level") # Use Beta_ID to match magenta_metadata

# Merge final_df with magenta_metadata to get correct sample IDs
aa_geno_full_long_corrected <- aa_geno_full_long %>%
    left_join(magenta_metadata, by = c("Sample" = "CGI")) %>% # Match CGI to get Beta_ID
    drop_na(Beta_ID) # Ensure only valid matches

# Merge on Beta_ID and CpG (PROBE)
aa_merged_df <- aa_betas_long %>%
    inner_join(aa_geno_full_long_corrected, by = c("Beta_ID" = "Beta_ID", "CpG" = "PROBE"))

# Standardize genotype notation: Convert "0|1" to "1|0"
aa_merged_df <- aa_merged_df %>%
    mutate(Genotype = ifelse(Genotype == "0|1", "1|0", Genotype))

# Summarize data: Compute mean BETA for each Sample-Genotype combination
summary_df <- aa_merged_df %>%
    group_by(Sample, Genotype) %>%
    summarise(mean_BETA = mean(BETA), .groups = "drop")

# Create boxplot
ggplot(summary_df, aes(x = Genotype, y = mean_BETA, fill = Genotype)) +
    geom_boxplot() +
    theme_minimal() +
    labs(
        title = "BETA Distributions by Genotype",
        x = "Genotype",
        y = "Mean BETA"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for readability


aa_merged_df_meqtl_individuals <- aa_merged_df %>% filter(ALT_Copies >= 1)



####### African differentiated meQTL (with assistance from gnomAD)
# Read in the file that has Horvath clock meQTL that were found only in AFR individuals, and that have an AF = 0 in all other gnomAD pops
afr_specific_meqtls <- read.table("horvath_gnomad_afr_specific_meqtl.txt", header = TRUE, sep = "\t")

# And read in the files that have three different approaches for finding African differentiated meQTL
afr_diff_meqtls_fc <- read_csv("african_differentiated_fc.csv")
afr_diff_meqtls_5_fc <- read_csv("african_differentiated_5_fc.csv")
afr_diff_meqtls_abs_diff <- read_csv("african_differentiated_abs_diff.csv")
afr_diff_meqtls_zscore <- read_csv("african_differentiated_zscore.csv")

# Create a new column id by joining V1, V2, V4, and V5, but removing the chr string from V1. For example: chr15 74171315 A G leads to 15:74171315:A:G
afr_specific_meqtls <- afr_specific_meqtls %>%
    mutate(id = paste0(gsub("chr", "", V1), ":", V2, ":", V4, ":", V5))

afr_diff_meqtls_fc <- afr_diff_meqtls_fc %>%
    mutate(id = paste0(gsub("chr", "", V1), ":", V2, ":", V4, ":", V5))

afr_diff_meqtls_5_fc <- afr_diff_meqtls_5_fc %>%
    mutate(id = paste0(gsub("chr", "", V1), ":", V2, ":", V4, ":", V5))

afr_diff_meqtls_abs_diff <- afr_diff_meqtls_abs_diff %>%
    mutate(id = paste0(gsub("chr", "", V1), ":", V2, ":", V4, ":", V5))

afr_diff_meqtls_zscore <- afr_diff_meqtls_zscore %>%
    mutate(id = paste0(gsub("chr", "", V1), ":", V2, ":", V4, ":", V5))


# Clean up the trailing underscore and number in the Variant column for aa_merged_df
aa_merged_df <- aa_merged_df %>%
    mutate(Variant = gsub("_\\d+$", "", Variant))

# Now we'll filter the aa_merged dataframe to only keep the variants that are in the afr-specific, or afr-differentiated dataframes
aa_geno_full_afr_specific <- aa_merged_df %>%
    filter(Variant %in% afr_specific_meqtls$id)

aa_geno_full_afr_diff_fc <- aa_merged_df %>%
    filter(Variant %in% afr_diff_meqtls_fc$id)

aa_geno_full_afr_diff_5_fc <- aa_merged_df %>%
    filter(Variant %in% afr_diff_meqtls_5_fc$id)

aa_geno_full_afr_diff_abs_diff <- aa_merged_df %>%
    filter(Variant %in% afr_diff_meqtls_abs_diff$id)

aa_geno_full_afr_diff_zscore <- aa_merged_df %>%
    filter(Variant %in% afr_diff_meqtls_zscore$id)

# Plot the number of variants for each of the four approaches above
# Combine the dataframes into a list
list_of_dataframes <- list(
    "AFR-Specific" = aa_geno_full_afr_specific,
    "AFR-Differentiated (Fold Change)" = aa_geno_full_afr_diff_fc,
    "AFR-Differentiated (5 Fold Change)" = aa_geno_full_afr_diff_5_fc,
    "AFR-Differentiated (Absolute Difference)" = aa_geno_full_afr_diff_abs_diff,
    "AFR-Differentiated (Z-score)" = aa_geno_full_afr_diff_zscore
)

# Count the number of unique variants in each dataframe
variant_counts <- sapply(list_of_dataframes, function(df) length(unique(df$Variant)))

# Create a dataframe for plotting
plot_data <- data.frame(
    Category = names(variant_counts),
    Count = unname(variant_counts)
)

# Order the categories by count in ascending order
plot_data$Category <- factor(plot_data$Category, levels = plot_data$Category[order(plot_data$Count)])

# Create the bar plot
ggplot(plot_data, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity") +
    labs(
        title = "Number of Unique Variants in Each Category",
        x = "Category",
        y = "Number of Unique Variants"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    ggtitle("") +
    scale_fill_brewer(palette = "Set2")


# Make a boxplot showing the methylation levels for a random CpG site, grouping by number of ALT_copies
# Select a random CpG site
random_cpg <- sample(unique(aa_merged_df$CpG), 1)

# Filter the data for the selected CpG
cpg_data <- aa_merged_df %>%
    filter(CpG == random_cpg)

# Create the boxplot
ggplot(cpg_data, aes(x = as.factor(ALT_Copies), y = Methylation_Level, fill = as.factor(ALT_Copies))) +
    geom_boxplot() +
    labs(
        title = paste("Methylation Levels for", random_cpg),
        x = "Number of ALT Copies",
        y = "Methylation Level",
        fill = "ALT Copies"
    ) +
    theme_minimal()

# Generate the above plot for all CpG sites
# Create a list to store the plots
plot_list <- list()

# Loop through each CpG site
for (cpg in unique(aa_merged_df$CpG)) {
    # Filter the data for the current CpG
    cpg_data <- aa_merged_df %>%
        filter(CpG == cpg)

    # Create the boxplot
    p <- ggplot(cpg_data, aes(x = as.factor(ALT_Copies), y = Methylation_Level, fill = as.factor(ALT_Copies))) +
        geom_boxplot() +
        labs(
            title = paste("Methylation Levels for", cpg),
            x = "Number of ALT Copies",
            y = "Methylation Level",
            fill = "ALT Copies"
        ) +
        theme_minimal()

    # Store the plot in the list
    plot_list[[cpg]] <- p
}



library(ggpubr) # For the boxplot

# Assuming aa_merged_df is your dataframe containing the methylation levels and ALT_Copies information

# Corrected Summarization with Handling for Insufficient Observations
methylation_diff_summary <- aa_merged_df %>%
    group_by(CpG) %>%
    summarise(
        mean_methylation_alt = mean(Methylation_Level[ALT_Copies >= 1], na.rm = TRUE),
        mean_methylation_ref = mean(Methylation_Level[ALT_Copies == 0], na.rm = TRUE),
        diff_methylation = mean_methylation_alt - mean_methylation_ref,
        p_value = if (sum(ALT_Copies >= 1, na.rm = TRUE) >= 2 & sum(ALT_Copies == 0, na.rm = TRUE) >= 2) {
            t.test(Methylation_Level[ALT_Copies >= 1], Methylation_Level[ALT_Copies == 0])$p.value
        } else {
            NA_real_ # Return NA if there are less than 2 observations
        },
        .groups = "drop"
    )

methylation_diff_summary_zscore_afr_diff <- aa_geno_full_afr_diff_zscore %>%
    group_by(CpG) %>%
    summarise(
        mean_methylation_alt = mean(Methylation_Level[ALT_Copies >= 1], na.rm = TRUE),
        mean_methylation_ref = mean(Methylation_Level[ALT_Copies == 0], na.rm = TRUE),
        diff_methylation = mean_methylation_alt - mean_methylation_ref,
        p_value = if (sum(ALT_Copies >= 1, na.rm = TRUE) >= 2 & sum(ALT_Copies == 0, na.rm = TRUE) >= 2) {
            t.test(Methylation_Level[ALT_Copies >= 1], Methylation_Level[ALT_Copies == 0])$p.value
        } else {
            NA_real_ # Return NA if there are less than 2 observations
        },
        .groups = "drop"
    )

# Adjust for multiple testing, only for the values that have a p-value
methylation_diff_summary <- methylation_diff_summary %>%
    mutate(
        p_value_adjusted = ifelse(!is.na(p_value), p.adjust(p_value, method = "BH"), NA_real_)
    )

methylation_diff_summary_zscore_afr_diff <- methylation_diff_summary_zscore_afr_diff %>%
    mutate(
        p_value_adjusted = ifelse(!is.na(p_value), p.adjust(p_value, method = "BH"), NA_real_)
    )

# Remove rows with NA p_value_adjusted
# methylation_diff_summary <- methylation_diff_summary %>%
#  drop_na(p_value_adjusted)

# Check for NA values in diff_methylation
if (any(is.na(methylation_diff_summary$diff_methylation))) {
    print("There are NA values in the diff_methylation column. Removing them now")
    methylation_diff_summary <- methylation_diff_summary %>%
        drop_na(diff_methylation)
}

if (any(is.na(methylation_diff_summary_zscore_afr_diff$diff_methylation))) {
    print("There are NA values in the diff_methylation column. Removing them now")
    methylation_diff_summary_zscore_afr_diff <- methylation_diff_summary_zscore_afr_diff %>%
        drop_na(diff_methylation)
}
# check if the dataframe has rows
if (nrow(methylation_diff_summary) == 0) {
    stop("The methylation_diff_summary dataframe is empty. Cannot create boxplot.")
}

if (nrow(methylation_diff_summary_zscore_afr_diff) == 0) {
    stop("The methylation_diff_summary_zscore_afr_diff dataframe is empty. Cannot create boxplot.")
}
# Create a histogram of the methylation difference
ggplot(methylation_diff_summary, aes(x = diff_methylation)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    labs(
        x = "Methylation Difference (ALT - REF)",
        y = "Number of CpG Sites"
    ) +
    theme_classic()

ggplot(methylation_diff_summary_zscore_afr_diff, aes(x = diff_methylation)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    labs(
        x = "Methylation Difference (ALT - REF)",
        y = "Number of CpG Sites"
    ) +
    theme_classic()

# Create a density plot of the methylation difference
ggplot(methylation_diff_summary, aes(x = diff_methylation)) +
    geom_density(fill = "skyblue", color = "black", alpha = 0.7) +
    labs(
        x = "Methylation Difference (ALT - REF)",
        y = "Density"
    ) +
    theme_classic()

ggplot(methylation_diff_summary_zscore_afr_diff, aes(x = diff_methylation)) +
    geom_density(fill = "skyblue", color = "black", alpha = 0.7) +
    labs(
        x = "Methylation Difference (ALT - REF)",
        y = "Density"
    ) +
    theme_classic()

# Remove rows with NA p_value_adjusted, because it is now used
methylation_diff_summary <- methylation_diff_summary %>%
    drop_na(p_value_adjusted)

methylation_diff_summary_zscore_afr_diff <- methylation_diff_summary_zscore_afr_diff %>%
    drop_na(p_value_adjusted)

# Create the volcano plot
# Handle p-values of 0 by replacing them with a very small number
# This prevents -log10(0) which would be infinity
methylation_diff_summary$p_value_adjusted_fix <- ifelse(
    methylation_diff_summary$p_value_adjusted == 0,
    1e-300, # A very small value instead of 0
    methylation_diff_summary$p_value_adjusted
)

methylation_diff_summary_zscore_afr_diff$p_value_adjusted_fix <- ifelse(
    methylation_diff_summary_zscore_afr_diff$p_value_adjusted == 0,
    1e-300, # A very small value instead of 0
    methylation_diff_summary_zscore_afr_diff$p_value_adjusted
)

# First, check the range of your methylation differences
# print(range(methylation_diff_summary$diff_methylation, na.rm = TRUE))

ggplot(methylation_diff_summary, aes(x = diff_methylation, y = -log10(p_value_adjusted_fix))) +
    geom_point(aes(color = ifelse(p_value_adjusted < 0.05, "Significant", "Not Significant")), size = 3) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    # Set both x and y limits to ensure we see all data points
    coord_cartesian(
        ylim = c(0, 325),
        # Expand x-axis limits by 5% on each side to ensure we see all points
        xlim = c(
            min(methylation_diff_summary$diff_methylation, na.rm = TRUE) * 1.05,
            max(methylation_diff_summary$diff_methylation, na.rm = TRUE) * 1.05
        )
    ) +
    labs(
        x = "Methylation Difference (ALT - REF)",
        y = "-log10(Adjusted p-value)",
        color = "Significance"
    ) +
    theme_classic() +
    theme(legend.position = "bottom")

ggplot(methylation_diff_summary_zscore_afr_diff, aes(x = diff_methylation, y = -log10(p_value_adjusted_fix))) +
    geom_point(aes(color = ifelse(p_value_adjusted < 0.05, "Significant", "Not Significant")), size = 3) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    # Set both x and y limits to ensure we see all data points
    coord_cartesian(
        ylim = c(0, 325),
        # Expand x-axis limits by 5% on each side to ensure we see all points
        xlim = c(
            min(methylation_diff_summary_zscore_afr_diff$diff_methylation, na.rm = TRUE) * 1.05,
            max(methylation_diff_summary_zscore_afr_diff$diff_methylation, na.rm = TRUE) * 1.05
        )
    ) +
    labs(
        x = "Methylation Difference (ALT - REF)",
        y = "-log10(Adjusted p-value)",
        color = "Significance"
    ) +
    theme_classic() +
    theme(legend.position = "bottom")

# Create a boxplot showing distribution of diff_methylation
if (nrow(methylation_diff_summary) > 0) {
    ggplot(methylation_diff_summary, aes(y = diff_methylation)) +
        geom_boxplot(fill = "lightblue") +
        labs(
            y = "Methylation Difference (ALT - REF)",
            title = "Boxplot of Methylation Difference (ALT vs. REF)"
        ) +
        theme_classic() +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )
} else {
    print("methylation_diff_summary has no rows. Skipping boxplot creation.")
}

if (nrow(methylation_diff_summary_zscore_afr_diff) > 0) {
    ggplot(methylation_diff_summary_zscore_afr_diff, aes(y = diff_methylation)) +
        geom_boxplot(fill = "lightblue") +
        labs(
            y = "Methylation Difference (ALT - REF)",
            title = "Boxplot of Methylation Difference (ALT vs. REF)"
        ) +
        theme_classic() +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )
} else {
    print("methylation_diff_summary_zscore_afr_diff has no rows. Skipping boxplot creation.")
}


###### Age prediction analysis
# Read in the age predictions for MAGENTA
age_preds <- read_csv("bio_age_estimates_magenta_metadata.csv")
age_preds <- age_preds %>% select(-...1) # Remove useless column

# Join the age predictions to the genotyping and methylation data
aa_merged_df <- aa_merged_df %>%
    left_join(age_preds, by = c("Beta_ID" = "id"))

aa_geno_full_afr_diff_zscore <- aa_geno_full_afr_diff_zscore %>%
    left_join(age_preds, by = c("Beta_ID" = "id"))

aa_geno_full_afr_diff_abs_diff <- aa_geno_full_afr_diff_abs_diff %>%
    left_join(age_preds, by = c("Beta_ID" = "id"))

aa_geno_full_afr_diff_5_fc <- aa_geno_full_afr_diff_5_fc %>%
    left_join(age_preds, by = c("Beta_ID" = "id"))


# Calculate the absolute age acceleration for the Horvath clock, which is our "error" measure
aa_merged_df$horvath_error <- abs(aa_merged_df$Horvath - aa_merged_df$AGE_OF_EXAM.x)

aa_geno_full_afr_diff_zscore$horvath_error <- abs(aa_geno_full_afr_diff_zscore$Horvath - aa_geno_full_afr_diff_zscore$AGE_OF_EXAM.x)
aa_geno_full_afr_diff_abs_diff$horvath_error <- abs(aa_geno_full_afr_diff_abs_diff$Horvath - aa_geno_full_afr_diff_abs_diff$AGE_OF_EXAM.x)
aa_geno_full_afr_diff_5_fc$horvath_error <- abs(aa_geno_full_afr_diff_5_fc$Horvath - aa_geno_full_afr_diff_5_fc$AGE_OF_EXAM.x)


variant_counts <- aa_merged_df %>%
    group_by(Sample) %>%
    summarize(
        Variant_Count = sum(ALT_Copies, na.rm = TRUE),
        horvath_abs_error = abs(first(horvath_error))
    )

afr_diff_zscore_variant_counts <- aa_geno_full_afr_diff_zscore %>%
    group_by(Sample) %>%
    summarize(
        Variant_Count = sum(ALT_Copies, na.rm = TRUE),
        horvath_abs_error = abs(first(horvath_error))
    )

afr_diff_abs_diff_variant_counts <- aa_geno_full_afr_diff_abs_diff %>%
    group_by(Sample) %>%
    summarize(
        Variant_Count = sum(ALT_Copies, na.rm = TRUE),
        horvath_abs_error = abs(first(horvath_error))
    )

afr_diff_5_fc_variant_counts <- aa_geno_full_afr_diff_5_fc %>%
    group_by(Sample) %>%
    summarize(
        Variant_Count = sum(ALT_Copies, na.rm = TRUE),
        horvath_abs_error = abs(first(horvath_error))
    )

ggplot(variant_counts, aes(x = Variant_Count, y = horvath_abs_error)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    labs(
        title = "Relationship Between Variant Count and Absolute Horvath Clock Error",
        x = "Number of Variants (ALT copies)",
        y = "Absolute Horvath Clock Error",
        caption = "Each point represents an individual"
    ) +
    theme_classic()

ggplot(afr_diff_zscore_variant_counts, aes(x = Variant_Count, y = horvath_abs_error)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    labs(
        title = "Relationship Between Variant Count (Z-Score AFR Diff) and Absolute Horvath Clock Error",
        x = "Number of Variants (ALT copies)",
        y = "Absolute Horvath Clock Error",
        caption = "Each point represents an individual"
    ) +
    theme_classic()

ggplot(afr_diff_abs_diff_variant_counts, aes(x = Variant_Count, y = horvath_abs_error)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    labs(
        title = "Relationship Between Variant Count (Abs Diff AFR) and Absolute Horvath Clock Error",
        x = "Number of Variants (ALT copies)",
        y = "Absolute Horvath Clock Error",
        caption = "Each point represents an individual"
    ) +
    theme_classic()

ggplot(afr_diff_5_fc_variant_counts, aes(x = Variant_Count, y = horvath_abs_error)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    labs(
        title = "Relationship Between Variant Count (5 FC AFR Diff) and Absolute Horvath Clock Error",
        x = "Number of Variants (ALT copies)",
        y = "Absolute Horvath Clock Error",
        caption = "Each point represents an individual"
    ) +
    theme_classic()



# Step 1: Calculate the number of variants per individual where ALT_copies >= 1
variant_counts <- aa_merged_df %>%
    group_by(Sample) %>%
    summarize(
        # Count sites where ALT_copies >= 1
        Variant_Count = sum(ALT_Copies >= 1, na.rm = TRUE),
        # Count sites where ALT_copies = 0
        Reference_Count = sum(ALT_Copies == 0, na.rm = TRUE),
        # Get the horvath_error for each individual
        horvath_error = first(horvath_error)
    )

afr_diff_zscore_variant_counts <- aa_geno_full_afr_diff_zscore %>%
    group_by(Sample) %>%
    summarize(
        # Count sites where ALT_copies >= 1
        Variant_Count = sum(ALT_Copies >= 1, na.rm = TRUE),
        # Count sites where ALT_copies = 0
        Reference_Count = sum(ALT_Copies == 0, na.rm = TRUE),
        # Get the horvath_error for each individual
        horvath_error = first(horvath_error)
    )

afr_diff_zscore_variant_counts <- aa_geno_full_afr_diff_zscore %>%
    group_by(Sample) %>%
    summarize(
        # Count sites where ALT_copies >= 1
        Variant_Count = sum(ALT_Copies >= 1, na.rm = TRUE),
        # Count sites where ALT_copies = 0
        Reference_Count = sum(ALT_Copies == 0, na.rm = TRUE),
        # Get the horvath_error for each individual
        horvath_error = first(horvath_error)
    )
# Step 2: Create a scatter plot showing the relationship
ggplot(variant_counts, aes(x = Variant_Count, y = horvath_error)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) + # Add a linear trend line
    labs(
        title = "Relationship Between Variant Sites and Horvath Clock Error",
        x = "Number of Sites with ALT_copies >= 1",
        y = "Horvath Clock Error",
        caption = "Each point represents an individual"
    ) +
    theme_classic()

ggplot(afr_diff_zscore_variant_counts, aes(x = Variant_Count, y = horvath_error)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) + # Add a linear trend line
    labs(
        title = "Relationship Between Variant Sites (Z-Score AFR Diff) and Horvath Clock Error",
        x = "Number of Sites with ALT_copies >= 1",
        y = "Horvath Clock Error",
        caption = "Each point represents an individual"
    ) +
    theme_classic()

ggplot(afr_diff_abs_diff_variant_counts, aes(x = Variant_Count, y = horvath_abs_error)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) + # Add a linear trend line
    labs(
        title = "Relationship Between Variant Sites (Abs Diff AFR-Diff) and Horvath Clock Error",
        x = "Number of Sites with ALT_copies >= 1",
        y = "Horvath Clock Error",
        caption = "Each point represents an individual"
    ) +
    theme_classic()

ggplot(afr_diff_5_fc_variant_counts, aes(x = Variant_Count, y = horvath_abs_error)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) + # Add a linear trend line
    labs(
        title = "Relationship Between Variant Sites (5 FC AFR-Diff) and Horvath Clock Error",
        x = "Number of Sites with ALT_copies >= 1",
        y = "Horvath Clock Error",
        caption = "Each point represents an individual"
    ) +
    theme_classic()

# Optional: Calculate and display correlation statistics
cor_result <- cor.test(afr_diff_zscore_variant_counts$Variant_Count, afr_diff_zscore_variant_counts$horvath_error)
print(paste("Correlation between variant count and horvath error:", round(cor_result$estimate, 3)))
print(paste("p-value:", round(cor_result$p.value, 4)))

cor_result <- cor.test(afr_diff_abs_diff_variant_counts$Variant_Count, afr_diff_abs_diff_variant_counts$horvath_abs_error)
print(paste("Correlation between variant count and horvath error:", round(cor_result$estimate, 3)))
print(paste("p-value:", round(cor_result$p.value, 4)))

cor_result <- cor.test(afr_diff_5_fc_variant_counts$Variant_Count, afr_diff_5_fc_variant_counts$horvath_abs_error)
print(paste("Correlation between variant count and horvath error:", round(cor_result$estimate, 3)))
print(paste("p-value:", round(cor_result$p.value, 4)))


###### meQTL burden analysis
# Step 1: Summarize the data by summing ALT_Copies for each sample
variant_counts <- aa_geno_full_afr_diff_zscore %>%
    group_by(Sample) %>%
    summarize(
        alt_copies_sum = sum(ALT_Copies, na.rm = TRUE),
        horvath_error = first(horvath_error) # Assuming horvath_error is constant per sample
    )
variant_counts_abs_diff <- aa_geno_full_afr_diff_abs_diff %>%
    group_by(Sample) %>%
    summarize(
        alt_copies_sum = sum(ALT_Copies, na.rm = TRUE),
        horvath_error = first(horvath_error) # Assuming horvath_error is constant per sample
    )
variant_counts_5_fc_diff <- aa_geno_full_afr_diff_5_fc %>%
    group_by(Sample) %>%
    summarize(
        alt_copies_sum = sum(ALT_Copies, na.rm = TRUE),
        horvath_error = first(horvath_error) # Assuming horvath_error is constant per sample
    )

# Step 2: Examine the relationship
print(variant_counts)
print(variant_counts_abs_diff)
print(variant_counts_5_fc_diff)


# Additional visualization: Create bins of ALT_Copies sums to check for patterns
variant_counts <- variant_counts %>%
    mutate(alt_copies_bin = cut(alt_copies_sum,
        breaks = quantile(alt_copies_sum, probs = seq(0, 1, 0.25), na.rm = TRUE),
        include.lowest = TRUE,
        labels = c("Q1", "Q2", "Q3", "Q4")
    ))
variant_counts_abs_diff <- variant_counts_abs_diff %>%
    mutate(alt_copies_bin = cut(alt_copies_sum,
        breaks = quantile(alt_copies_sum, probs = seq(0, 1, 0.25), na.rm = TRUE),
        include.lowest = TRUE,
        labels = c("Q1", "Q2", "Q3", "Q4")
    ))
variant_counts_5_fc_diff <- variant_counts_5_fc_diff %>%
    mutate(alt_copies_bin = cut(alt_copies_sum,
        breaks = quantile(alt_copies_sum, probs = seq(0, 1, 0.25), na.rm = TRUE),
        include.lowest = TRUE,
        labels = c("Q1", "Q2", "Q3", "Q4")
    ))

# Box plot of horvath error by ALT_Copies sum quartile
ggplot(variant_counts, aes(x = alt_copies_bin, y = horvath_error, fill = alt_copies_bin)) +
    geom_boxplot(alpha = 0.7) +
    labs(
        x = "ALT_Copies Sum Quartile",
        y = "Horvath Error"
    ) +
    theme_classic() +
    theme(legend.position = "none")
ggplot(variant_counts_abs_diff, aes(x = alt_copies_bin, y = horvath_error, fill = alt_copies_bin)) +
    geom_boxplot(alpha = 0.7) +
    labs(
        x = "ALT_Copies Sum Quartile",
        y = "Horvath Error"
    ) +
    theme_classic() +
    theme(legend.position = "none")
ggplot(variant_counts_5_fc_diff, aes(x = alt_copies_bin, y = horvath_error, fill = alt_copies_bin)) +
    geom_boxplot(alpha = 0.7) +
    labs(
        x = "ALT_Copies Sum Quartile",
        y = "Horvath Error"
    ) +
    theme_classic() +
    theme(legend.position = "none")

###### Distribution of meQTL betas
meqtl_effects <- aa_geno_full_afr_diff_zscore %>%
    select(Variant, CpG, BETA) %>%
    distinct() %>%
    mutate(abs_effect_beta = abs(BETA)) %>%
    arrange(desc(abs_effect_beta))

meqtl_effects %>%
    ggplot(aes(x = BETA)) +
    geom_density() +
    theme_classic()

meqtl_effects %>%
    ggplot(aes(x = abs_effect_beta)) +
    geom_density() +
    theme_classic()

meqtl_effects <- aa_geno_full_afr_diff_abs_diff %>%
    select(Variant, CpG, BETA) %>%
    distinct() %>%
    mutate(abs_effect_beta = abs(BETA)) %>%
    arrange(desc(abs_effect_beta))

meqtl_effects %>%
    ggplot(aes(x = BETA)) +
    geom_density() +
    theme_classic()

meqtl_effects %>%
    ggplot(aes(x = abs_effect_beta)) +
    geom_density() +
    theme_classic()

meqtl_effects <- aa_geno_full_afr_diff_5_fc %>%
    select(Variant, CpG, BETA) %>%
    distinct() %>%
    mutate(abs_effect_beta = abs(BETA)) %>%
    arrange(desc(abs_effect_beta))

meqtl_effects %>%
    ggplot(aes(x = BETA)) +
    geom_density() +
    theme_classic()

meqtl_effects %>%
    ggplot(aes(x = abs_effect_beta)) +
    geom_density() +
    theme_classic()

# Extract the top 10 strongest meQTLs
top_meqtls <- meqtl_effects %>%
    head(10)

# Create an empty list to store results for each variant
all_meqtl_data <- list()

# Loop through all top 10 variants
for (i in 1:nrow(top_meqtls)) {
    current_variant <- top_meqtls$Variant[i]

    # Extract data for the current variant
    current_data <- aa_geno_full_afr_diff_zscore %>%
        filter(Variant == current_variant)

    # Store in the list with the variant name as identifier
    all_meqtl_data[[current_variant]] <- current_data
}

# If you want to combine all results into a single dataframe with a variant identifier
all_meqtl_df <- bind_rows(all_meqtl_data, .id = "Variant_ID")

# Select data for the strongest meQTL
strongest_meqtl <- top_meqtls$Variant[1]
strongest_data <- aa_geno_full_afr_diff_zscore %>%
    filter(Variant == strongest_meqtl)

# Create individual profile plot for the strongest meQTL
individual_plot <- ggplot(strongest_data, aes(x = Sample, y = Methylation_Level, fill = factor(ALT_Copies))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
        title = paste("Individual Methylation Profiles for", strongest_meqtl),
        x = "Individual ID",
        y = "Methylation Value",
        fill = "Genotype"
    ) +
    scale_fill_manual(
        values = c("0" = "#E69F00", "1" = "#56B4E9", "2" = "#009E73"),
        labels = c("0" = "Homozygous Ref", "1" = "Heterozygous", "2" = "Homozygous Alt")
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

print(individual_plot)

# First, check if there are duplicates in the data
duplicates_check <- strongest_data %>%
    group_by(Sample, ALT_Copies) %>%
    summarise(count = n(), .groups = "drop") %>%
    filter(count > 1)

# If you have duplicates, you'll need to deduplicate before plotting
strongest_data_unique <- strongest_data %>%
    distinct(Sample, ALT_Copies, Methylation_Level, horvath_error)

all_meqtl_df_unique <- all_meqtl_df %>% distinct(Sample, ALT_Copies, Methylation_Level, horvath_error, Variant_ID)

# Now create the plot with the deduplicated data
# individual_plot <- ggplot(strongest_data_unique, aes(x = factor(ALT_Copies), y = Methylation_Level, fill = factor(ALT_Copies))) +
#   geom_boxplot() +
#   geom_jitter(width = 0.2, alpha = 0.5) +
#   labs(
#     title = paste("Methylation Levels by Genotype for", strongest_meqtl),
#     x = "Number of ALT Alleles",
#     y = "Methylation Level",
#     fill = "Genotype"
#   ) +
#   scale_fill_manual(values = c("0" = "#E69F00", "1" = "#56B4E9", "2" = "#009E73"),
#                     labels = c("0" = "Homozygous Ref", "1" = "Heterozygous", "2" = "Homozygous Alt")) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )
#
# print(individual_plot)

# ggplot(strongest_data_unique, aes(x = factor(ALT_Copies), y = horvath_error, fill = factor(ALT_Copies))) +
#   geom_boxplot() +
#   geom_jitter(width = 0.2, alpha = 0.5) +
#   labs(
#     title = paste("Horvath Error by Genotype for", strongest_meqtl),
#     x = "Number of ALT Alleles",
#     y = "Horvath Error",
#     fill = "Genotype"
#   ) +
#   scale_fill_manual(values = c("0" = "#E69F00", "1" = "#56B4E9", "2" = "#009E73"),
#                     labels = c("0" = "Homozygous Ref", "1" = "Heterozygous", "2" = "Homozygous Alt")) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )

# Function to remove outliers using Tukey's rule
remove_outliers <- function(x) {
    q1 <- quantile(x, 0.25, na.rm = TRUE)
    q3 <- quantile(x, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower <- q1 - 1.5 * iqr
    upper <- q3 + 1.5 * iqr
    x[x < lower] <- lower
    x[x > upper] <- upper
    return(x)
}

# Function to create both plots for a given variant
create_variant_plots <- function(variant_data, variant_name) {
    # Deduplicate and collapse genotypes
    variant_data_unique <- variant_data %>%
        distinct() %>%
        mutate(
            ALT_Carrier = ifelse(ALT_Copies == 0, "0", "1+"),
            Methylation_Level = remove_outliers(Methylation_Level),
            horvath_error = remove_outliers(horvath_error)
        )

    # Define fill colors and labels
    fill_colors <- c("0" = "#E69F00", "1+" = "#0072B2")
    fill_labels <- c("0" = "Homozygous Ref", "1+" = "ALT Carrier (1 or 2 copies)")

    # Plot 1: Methylation Levels by Collapsed Genotype
    methylation_plot <- ggplot(variant_data_unique, aes(x = ALT_Carrier, y = Methylation_Level, fill = ALT_Carrier)) +
        geom_boxplot() +
        geom_jitter(width = 0.2, alpha = 0.5) +
        labs(
            x = "ALT Allele Status",
            y = "Methylation Level",
            fill = "Genotype"
        ) +
        scale_fill_manual(values = fill_colors, labels = fill_labels) +
        scale_x_discrete(labels = fill_labels) +
        theme_pubr(base_size = 16) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none",
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14)
        )

    # Plot 2: Horvath Error by Collapsed Genotype
    horvath_plot <- ggplot(variant_data_unique, aes(x = ALT_Carrier, y = horvath_error, fill = ALT_Carrier)) +
        geom_boxplot() +
        geom_jitter(width = 0.2, alpha = 0.5) +
        labs(
            x = "ALT Allele Status",
            y = "Horvath Error",
            fill = "Genotype"
        ) +
        scale_fill_manual(values = fill_colors, labels = fill_labels) +
        scale_x_discrete(labels = fill_labels) +
        theme_pubr(base_size = 16) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none",
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14)
        )

    # Return both plots
    return(list(methylation_plot = methylation_plot, horvath_plot = horvath_plot))
}



strong_meqtl_effects <- meqtl_effects %>% filter(abs_effect_beta > 0.5)

# Create the function to find the top variant for each CpG site
filter_top_variants <- function(data) {
    # Group by CpG site and find the row with max effect size
    result <- data %>%
        group_by(CpG) %>%
        slice_max(abs_effect_beta, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        arrange(desc(abs_effect_beta))

    return(result)
}

strong_meqtl_effects_filtered <- filter_top_variants(strong_meqtl_effects)

# View the result
print(strong_meqtl_effects_filtered)

# Extract the top 10 strongest meQTLs
top_meqtls <- strong_meqtl_effects_filtered %>%
    head(10)

# Create an empty list to store all plots
all_plots <- list()

# Loop through all top 10 variants
for (i in 1:nrow(top_meqtls)) {
    current_variant <- top_meqtls$Variant[i]

    # Extract data for the current variant
    current_data <- aa_geno_full_afr_diff_zscore %>%
        filter(Variant == current_variant)

    # Generate plots for this variant
    variant_plots <- create_variant_plots(current_data, current_variant)

    # Store plots in the list
    all_plots[[current_variant]] <- variant_plots
}

# Now you can access plots for each variant
# For example, to print both plots for the first variant:
first_variant <- top_meqtls$Variant[1]
print(all_plots[[first_variant]]$methylation_plot)
print(all_plots[[first_variant]]$horvath_plot)

# To print all plots for all variants
for (variant in names(all_plots)) {
    print(all_plots[[variant]]$methylation_plot)
    print(all_plots[[variant]]$horvath_plot)
}

# Display 2x2 grid plots for all variants to RStudio viewer
library(gridExtra)
library(ggplot2)

# Loop through all variants and display 2x2 grid plots
for (i in 1:nrow(top_meqtls)) {
    variant <- top_meqtls$Variant[i]

    # Create and print grid plot for this variant
    cat("Displaying plots for variant", i, ":", variant, "\n")

    grid.arrange(
        all_plots[[variant]]$methylation_plot,
        all_plots[[variant]]$horvath_plot,
        ncol = 2,
        top = grid::textGrob(paste("Variant:", variant), gp = grid::gpar(fontsize = 14, fontface = "bold"))
    )
}
# Pick the first variant
variant <- top_meqtls$Variant[1]

# Create the figure as a ggarrange object
fig5d <- ggarrange(
    all_plots[[variant]]$methylation_plot,
    all_plots[[variant]]$horvath_plot,
    ncol = 2, nrow = 1,
    common.legend = FALSE
) %>%
    annotate_figure(
        top = text_grob(
            paste("Variant:", variant),
            face = "bold", size = 16
        )
    )

# Save the object to RDS
saveRDS(fig5d, "fig5d.rds")


# First, get the top variant per CpG site
get_top_variants <- function(data) {
    data %>%
        group_by(CpG) %>%
        slice_max(abs_effect_beta, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        arrange(desc(abs_effect_beta))
}

# Function that handles both the analysis and visualization in one workflow
analyze_and_visualize_variants <- function(meQTL_data, individual_data) {
    cat("Step 1: Identifying top variants for each CpG...\n")
    top_variants <- get_top_variants(meQTL_data)
    cat("Found", nrow(top_variants), "unique CpG-variant pairs.\n")

    cat("Step 2: Preparing data for analysis...\n")
    # Pre-filter and label individual data once
    # FIX: Ensure we drop NAs in ALT_Copies to avoid errors
    clean_data <- individual_data %>%
        filter(!is.na(ALT_Copies), !is.na(Methylation_Level), !is.na(horvath_error)) %>%
        mutate(
            ALT_Group = if_else(ALT_Copies == 0, "Reference (0)", "Alternative (1 or 2)"),
            ALT_Group = factor(ALT_Group, levels = c("Reference (0)", "Alternative (1 or 2)"))
        )

    # Initialize results containers
    analysis_results <- list()
    plots <- list()

    cat("Step 3: Running Mann-Whitney U tests and generating plots...\n")
    # Iterate through each variant-CpG pair
    for (i in 1:nrow(top_variants)) {
        current_variant <- top_variants$Variant[i]
        current_cpg <- top_variants$CpG[i]

        # Subset data for this variant
        variant_data <- clean_data %>%
            filter(Variant == current_variant)

        # Verify group counts
        ref_count <- sum(variant_data$ALT_Group == "Reference (0)")
        alt_count <- sum(variant_data$ALT_Group == "Alternative (1 or 2)")

        if (ref_count < 3 || alt_count < 3) {
            cat(sprintf("  Skipping %s: Insufficient n (Ref=%d, Alt=%d)\n", current_variant, ref_count, alt_count))
            next
        }

        # helper for safe testing
        safe_test <- function(formula, data) {
            tryCatch(
                {
                    # conf.int = TRUE is needed to get median difference estimate
                    broom::tidy(wilcox.test(formula, data = data, conf.int = TRUE))
                },
                error = function(e) {
                    return(NULL)
                }
            )
        }

        # Run tests
        meth_test <- safe_test(Methylation_Level ~ ALT_Group, variant_data)
        clock_test <- safe_test(horvath_error ~ ALT_Group, variant_data)

        if (is.null(meth_test) || is.null(clock_test)) {
            cat("  Skipping", current_variant, ": Test execution failed.\n")
            next
        }

        # Store results
        # Mann-Whitney U estimates magnitude of difference (location shift)
        analysis_results[[length(analysis_results) + 1]] <- data.frame(
            Variant = current_variant,
            CpG = current_cpg,
            abs_effect_beta = top_variants$abs_effect_beta[i],
            Methylation_p_value = meth_test$p.value,
            Methylation_diff = meth_test$estimate, # Median difference
            Horvath_p_value = clock_test$p.value,
            Horvath_diff = clock_test$estimate # Median difference
        )

        # Generate Plot
        # Combine measurments for faceting
        plot_data <- variant_data %>%
            select(ALT_Group, Methylation_Level, horvath_error) %>%
            pivot_longer(cols = c(Methylation_Level, horvath_error), names_to = "Variable", values_to = "Value")

        p <- ggplot(plot_data, aes(x = ALT_Group, y = Value, fill = ALT_Group)) +
            geom_boxplot(alpha = 0.7, outlier.shape = NA) +
            geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
            facet_wrap(~Variable, scales = "free_y") +
            labs(
                title = paste("Variant:", current_variant, "| CpG:", current_cpg),
                subtitle = paste(
                    "Effect size:", round(top_variants$abs_effect_beta[i], 4),
                    "| Horvath p:", round(clock_test$p.value, 4),
                    "| Meth p:", round(meth_test$p.value, 4)
                ),
                y = "Value", x = "ALT Copy Group"
            ) +
            theme_classic() +
            theme(legend.position = "none")

        plots[[length(plots) + 1]] <- p
    }

    cat("Step 4: Compiling summary results...\n")
    summary_df <- bind_rows(analysis_results) %>%
        arrange(Methylation_p_value)

    cat("Analysis complete. Processed", nrow(summary_df), "valid associations.\n")

    return(list(
        top_variants = top_variants,
        summary = summary_df,
        plots = plots
    ))
}


# Usage
results <- analyze_and_visualize_variants(top_meqtls, aa_geno_full_afr_diff_zscore)

# View summary table
print(results$summary)

# View the plots
results$plots

# --- Custom Analysis for Variant 3:51639064:A:G with Outlier Removal ---
cat("\n--- Custom Analysis for Variant 3:51639064:A:G with Outlier Removal ---\n")

# 1. Select the specific variant data
target_variant <- "3:51639064:A:G"
# Ensure we use the same dataset as the main analysis
variant_data_custom <- aa_geno_full_afr_diff_zscore %>%
    filter(Variant == target_variant) %>%
    filter(!is.na(ALT_Copies), !is.na(horvath_error)) %>%
    mutate(
        ALT_Group = if_else(ALT_Copies == 0, "Reference (0)", "Alternative (1 or 2)"),
        ALT_Group = factor(ALT_Group, levels = c("Reference (0)", "Alternative (1 or 2)"))
    )

# 2. Key step: Remove top 2 outliers from Reference Clean
ref_sorted <- variant_data_custom %>%
    filter(ALT_Group == "Reference (0)") %>%
    arrange(desc(horvath_error))

# Check we have enough data to remove 2
if (nrow(ref_sorted) > 2) {
    # Remove top 2
    ref_clean <- ref_sorted %>% slice(-(1:2))

    # Get Alt data
    alt_data <- variant_data_custom %>%
        filter(ALT_Group == "Alternative (1 or 2)")

    # Combine
    custom_data_clean <- bind_rows(ref_clean, alt_data)

    cat("Removed 2 highest outliers from Reference group.\n")
    cat("Top 2 removed values:", paste(round(ref_sorted$horvath_error[1:2], 3), collapse = ", "), "\n")

    # 3. Run Mann-Whitney U Test
    custom_test <- broom::tidy(t.test(horvath_error ~ ALT_Group, data = custom_data_clean, conf.int = TRUE))
    print(custom_test)

    # 4. Plot
    p_custom <- ggplot(custom_data_clean, aes(x = ALT_Group, y = horvath_error, fill = ALT_Group)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
        labs(
            title = paste("Variant:", target_variant, "(Top 2 Ref Outliers Removed)"),
            subtitle = paste("Mann-Whitney U p-value:", round(custom_test$p.value, 5)),
            y = "Horvath Error", x = "ALT Copy Group"
        ) +
        theme_classic() +
        theme(legend.position = "none")

    print(p_custom)
} else {
    cat("Not enough reference samples to remove 2 outliers safely.\n")
}

# --- Analysis of Dual Significant Hits (Methylation & Horvath Error) ---
cat("\n--- Dual Significance Analysis (Lenient Filter: beta > 0.1) ---\n")

# Create lenient set (beta > 0.1 instead of 0.5)
lenient_meqtl_effects <- meqtl_effects %>% filter(abs_effect_beta > 0.1)
lenient_meqtl_effects_filtered <- filter_top_variants(lenient_meqtl_effects)

cat("Analyzing", nrow(lenient_meqtl_effects_filtered), "variants (beta > 0.1)...\n")

# Run analysis on lenient top variants
full_results <- analyze_and_visualize_variants(lenient_meqtl_effects_filtered, aa_geno_full_afr_diff_zscore)

# Filter for dual hits (p < 0.05 for both)
dual_hits <- full_results$summary %>%
    filter(Methylation_p_value < 0.05 & Horvath_p_value < 0.05) %>%
    arrange(Horvath_p_value)

cat("Found", nrow(dual_hits), "variants with significant differences in BOTH Methylation and Horvath Error (p < 0.05).\n")

# Display results
if (nrow(dual_hits) > 0) {
    print(dual_hits)

    # Plotting loop
    for (i in 1:nrow(dual_hits)) {
        var_name <- dual_hits$Variant[i]
        cat("Displaying dual hit:", var_name, "\n")

        # Match index in full results to retrieve the pre-generated plot
        idx <- which(full_results$top_variants$Variant == var_name)
        if (length(idx) > 0) {
            print(full_results$plots[[idx]])
        }
    }
} else {
    cat("No dual hits found.\n")
}

# --- Analysis of Horvath Coefficient Weights ---
cat("\n--- Horvath Coefficient Weight Analysis ---\n")

# hypothesis: strong meQTLs might not affect the clock because they are on low-weight CpGs

# Join meQTL effects with Horvath coefficients
# 'coefHorvath' comes from methylclock library
load_DNAm_Clocks_data()
horvath_weights <- lenient_meqtl_effects_filtered %>%
    left_join(coefHorvath, by = c("CpG" = "CpGmarker")) %>%
    rename(Horvath_Coefficient = CoefficientTraining) %>%
    filter(!is.na(Horvath_Coefficient))

cat("Matched", nrow(horvath_weights), "variants to Horvath coefficients.\n")

# Summary of weights for top hits
cat("\nSummary of Horvath Coefficients for Top meQTLs:\n")
summary(horvath_weights$Horvath_Coefficient)

# Plot meQTL Beta vs Horvath Coefficient
p_weights <- ggplot(horvath_weights, aes(x = abs_effect_beta, y = abs(Horvath_Coefficient))) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
    labs(
        title = "meQTL Strength vs Horvath Weight",
        subtitle = "Are strong meQTLs on important clock CpGs?",
        x = "meQTL Effect Size (abs_beta)",
        y = "Absolute Horvath Coefficient (Weight)"
    ) +
    theme_classic()

print(p_weights)

# Top 10 strongest meQTLs and their weights
cat("\nTop 10 Strongest meQTLs and their Clock Weights:\n")
top_10_weights <- horvath_weights %>%
    arrange(desc(abs_effect_beta)) %>%
    select(Variant, CpG, abs_effect_beta, Horvath_Coefficient) %>%
    head(10)
print(top_10_weights)


library(ggrepel) # For non-overlapping labels if needed

# Function to create a scatterplot of meQTL count vs horvath error
plot_meqtl_count_vs_horvath_error <- function(individual_data, meQTL_data = NULL, top_variants_only = TRUE) {
    # If using only top variants per CpG and meQTL_data is provided
    if (top_variants_only && !is.null(meQTL_data)) {
        # Get the top variant for each CpG
        top_variants <- meQTL_data %>%
            group_by(CpG) %>%
            slice_max(abs_effect_beta, n = 1, with_ties = FALSE) %>%
            ungroup() %>%
            pull(Variant)

        # Filter individual data to only include top variants
        filtered_data <- individual_data %>%
            filter(Variant %in% top_variants)
    } else {
        # Use all variants
        filtered_data <- individual_data
    }

    # Calculate meQTL count per person (variants where ALT_Copies > 0)
    person_meQTL_counts <- filtered_data %>%
        filter(ALT_Copies > 0) %>%
        group_by(CGI) %>%
        summarize(meQTL_count = n_distinct(Variant))

    # Get horvath error for each person
    # Assuming horvath error is the same for a person across all variants
    person_horvath_errors <- individual_data %>%
        select(CGI, horvath_error) %>%
        distinct()

    # Join the datasets
    plot_data <- person_meQTL_counts %>%
        left_join(person_horvath_errors, by = "CGI")

    # Calculate correlation
    correlation <- cor.test(plot_data$meQTL_count, plot_data$horvath_error)
    cor_value <- round(correlation$estimate, 3)
    p_value <- round(correlation$p.value, 4)

    # Create the scatter plot
    p <- ggplot(plot_data, aes(x = meQTL_count, y = horvath_error)) +
        geom_point(alpha = 0.7, size = 2, aes(color = meQTL_count)) +
        geom_smooth(method = "lm", se = F, color = "red", linetype = "dashed") +
        scale_color_viridis_c() +
        labs(
            title = "Relationship Between Number of meQTL Variants and Horvath Error",
            subtitle = paste0("Correlation = ", cor_value, " (p = ", p_value, ")"),
            x = "Number of meQTL Variants (ALT_Copies > 0)",
            y = "Horvath Error",
            color = "Number of\nmeQTL Variants"
        ) +
        theme_classic() +
        theme(
            legend.position = "right",
            plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold")
        )

    # # Identify outliers (optional)
    # # Using IQR method to find potential outliers
    # q1 <- quantile(plot_data$horvath_error, 0.25)
    # q3 <- quantile(plot_data$horvath_error, 0.75)
    # iqr <- q3 - q1
    # upper_bound <- q3 + 1.5 * iqr
    # lower_bound <- q1 - 1.5 * iqr
    #
    # outliers <- plot_data %>%
    #   filter(horvath_error > upper_bound | horvath_error < lower_bound |
    #            meQTL_count > quantile(meQTL_count, 0.95))
    #
    # # If there are outliers, label them
    # if(nrow(outliers) > 0 && nrow(outliers) <= 10) {
    #   p <- p +
    #     geom_text_repel(
    #       data = outliers,
    #       aes(label = CGI),
    #       box.padding = 0.5,
    #       point.padding = 0.3,
    #       segment.color = "grey50"
    #     )


    # Return both the plot and the data
    return(list(
        plot = p,
        data = plot_data,
        correlation = list(
            coefficient = cor_value,
            p_value = p_value
        )
    ))
}

# Example usage:
results <- plot_meqtl_count_vs_horvath_error(aa_geno_full_afr_diff_zscore, strong_meqtl_effects_filtered, top_variants_only = F)
results_abs_diff <- plot_meqtl_count_vs_horvath_error(aa_geno_full_afr_diff_abs_diff, strong_meqtl_effects_filtered, top_variants_only = F)
results_5_fc <- plot_meqtl_count_vs_horvath_error(aa_geno_full_afr_diff_5_fc, meqtl_effects, top_variants_only = F)

# View the plot
results$plot
results_abs_diff$plot
results_5_fc$plot

# Check the correlation results
results$correlation
results_abs_diff$correlation
results_5_fc$correlation

# Get summary statistics for the data
summary(results$data)
summary(results_abs_diff$data)
summary(results_5_fc$data)


library(ggrepel) # For non-overlapping labels if needed

# Function to create a scatterplot of weighted meQTL score vs horvath error
plot_weighted_meqtl_vs_horvath_error <- function(individual_data, meQTL_data, top_variants_only = TRUE) {
    # If using only top variants per CpG and meQTL_data is provided
    if (top_variants_only) {
        # Get the top variant for each CpG with their effect sizes
        top_variant_data <- meQTL_data %>%
            group_by(CpG) %>%
            slice_max(abs_effect_beta, n = 1, with_ties = FALSE) %>%
            ungroup() %>%
            select(Variant, CpG, abs_effect_beta)

        # Filter individual data to only include top variants
        analysis_data <- individual_data %>%
            filter(Variant %in% top_variant_data$Variant) %>%
            # Join with the variant effect data
            left_join(top_variant_data, by = "Variant")
    } else {
        # Use all variants, but make sure we have the effect sizes
        analysis_data <- individual_data %>%
            left_join(meQTL_data %>% select(Variant, CpG, abs_effect_beta), by = "Variant")
    }

    # Calculate weighted meQTL score per person
    # For each person:
    # 1. Keep only variants where they have alternative alleles (ALT_Copies > 0)
    # 2. Weight each variant by its absolute effect size
    # 3. Sum these weighted values to get a total weighted score
    person_weighted_scores <- analysis_data %>%
        filter(ALT_Copies > 0) %>%
        # Weight by both ALT_Copies and effect size
        mutate(weighted_effect = ALT_Copies * abs_effect_beta) %>%
        group_by(CGI) %>%
        summarize(
            # Sum of weighted effects
            weighted_meQTL_score = sum(weighted_effect),
            # Regular count for comparison
            meQTL_count = n_distinct(Variant)
        )

    # Get horvath error for each person
    # Assuming horvath error is the same for a person across all variants
    person_horvath_errors <- individual_data %>%
        select(CGI, horvath_error) %>%
        distinct()

    # Join the datasets
    plot_data <- person_weighted_scores %>%
        left_join(person_horvath_errors, by = "CGI")

    # Calculate correlations for both weighted score and raw count
    # weighted_correlation <- cor.test(plot_data$weighted_meQTL_score, plot_data$horvath_error)
    # weighted_cor_value <- round(weighted_correlation$estimate, 3)
    # weighted_p_value <- round(weighted_correlation$p.value, 4)

    # count_correlation <- cor.test(plot_data$meQTL_count, plot_data$horvath_error)
    # count_cor_value <- round(count_correlation$estimate, 3)
    # count_p_value <- round(count_correlation$p.value, 4)

    # Create the scatter plot with weighted score
    p_weighted <- ggplot(plot_data, aes(x = weighted_meQTL_score, y = horvath_error)) +
        geom_point(alpha = 0.7, size = 2, aes(color = weighted_meQTL_score)) +
        geom_smooth(method = "lm", se = F, color = "red", linetype = "dashed") +
        scale_color_viridis_c() +
        labs(
            title = "Relationship Between Weighted meQTL Score and Horvath Error",
            # subtitle = paste0("Correlation = ", weighted_cor_value, " (p = ", weighted_p_value, ")"),
            x = "Weighted meQTL Score (sum of ALT_Copies × abs_effect_beta)",
            y = "Horvath Error",
            color = "Weighted\nmeQTL Score"
        ) +
        theme_classic() +
        theme(
            legend.position = "right",
            plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold")
        )

    # Create scatter plot with standard count for comparison
    p_count <- ggplot(plot_data, aes(x = meQTL_count, y = horvath_error)) +
        geom_point(alpha = 0.7, size = 2, aes(color = meQTL_count)) +
        geom_smooth(method = "lm", se = F, color = "red", linetype = "dashed") +
        scale_color_viridis_c() +
        labs(
            title = "Relationship Between Number of meQTL Variants and Horvath Error",
            # subtitle = paste0("Correlation = ", count_cor_value, " (p = ", count_p_value, ")"),
            x = "Number of meQTL Variants (ALT_Copies > 0)",
            y = "Horvath Error",
            color = "Number of\nmeQTL Variants"
        ) +
        theme_classic() +
        theme(
            legend.position = "right",
            plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold")
        )

    # # Identify outliers for the weighted score plot
    # q1 <- quantile(plot_data$horvath_error, 0.25)
    # q3 <- quantile(plot_data$horvath_error, 0.75)
    # iqr <- q3 - q1
    # upper_bound <- q3 + 1.5 * iqr
    # lower_bound <- q1 - 1.5 * iqr
    #
    # outliers <- plot_data %>%
    #   filter(horvath_error > upper_bound | horvath_error < lower_bound |
    #            weighted_meQTL_score > quantile(weighted_meQTL_score, 0.95))
    #
    # # If there are outliers, label them
    # if(nrow(outliers) > 0 && nrow(outliers) <= 10) {
    #   p_weighted <- p_weighted +
    #     geom_text_repel(
    #       data = outliers,
    #       aes(label = CGI),
    #       box.padding = 0.5,
    #       point.padding = 0.3,
    #       segment.color = "grey50"
    #     )
    # }
    #
    # Return both plots and the data
    return(list(
        weighted_plot = p_weighted,
        count_plot = p_count,
        data = plot_data
        # correlations = list(
        # weighted = list(coefficient = weighted_cor_value, p_value = weighted_p_value),
        # count = list(coefficient = count_cor_value, p_value = count_p_value)
    ))
}


results <- plot_weighted_meqtl_vs_horvath_error(aa_geno_full_afr_diff_zscore, meqtl_effects, top_variants_only = F)
results_abs_diff <- plot_weighted_meqtl_vs_horvath_error(aa_geno_full_afr_diff_abs_diff, meqtl_effects, top_variants_only = F)
results_5_fc <- plot_weighted_meqtl_vs_horvath_error(aa_geno_full_afr_diff_5_fc, meqtl_effects, top_variants_only = F)


# View the weighted plot
results$weighted_plot
results_abs_diff$weighted_plot
results_5_fc$weighted_plot

# View the count plot for comparison
results$count_plot
results_abs_diff$count_plot
results_5_fc$count_plot



library(broom) # For tidy statistical output

# Function to create plot comparing horvath error differences by variant effect size
plot_horvath_diff_by_beta <- function(individual_data, meQTL_data, top_variants_only = TRUE) {
    # Get the variants to analyze
    if (top_variants_only) {
        # Select only top variant per CpG
        variants_to_analyze <- meQTL_data %>%
            group_by(CpG) %>%
            slice_max(abs_effect_beta, n = 1, with_ties = FALSE) %>%
            ungroup()
    } else {
        # Use all variants
        variants_to_analyze <- meQTL_data
    }

    # Initialize empty list to store results
    variant_results <- list()

    # For each variant, calculate the difference in horvath error
    variant_stats <- data.frame()

    for (i in 1:nrow(variants_to_analyze)) {
        current_variant <- variants_to_analyze$Variant[i]
        current_beta <- variants_to_analyze$abs_effect_beta[i]
        current_cpg <- variants_to_analyze$CpG[i]

        # Get data for this variant
        variant_data <- individual_data %>%
            filter(Variant == current_variant) %>%
            mutate(ALT_Group = ifelse(ALT_Copies == 0, "Reference (0)", "Alternative (1+)"))

        # Skip if we don't have enough data in both groups
        if (sum(variant_data$ALT_Group == "Reference (0)") < 3 ||
            sum(variant_data$ALT_Group == "Alternative (1+)") < 3) {
            next
        }

        # Calculate mean horvath error for each group
        group_means <- variant_data %>%
            group_by(ALT_Group) %>%
            summarize(
                mean_horvath_error = mean(horvath_error, na.rm = TRUE),
                sd_horvath_error = sd(horvath_error, na.rm = TRUE),
                n = n()
            )

        # Calculate the difference (Alternative - Reference)
        if (nrow(group_means) == 2) {
            horvath_diff <- group_means$mean_horvath_error[group_means$ALT_Group == "Alternative (1+)"] -
                group_means$mean_horvath_error[group_means$ALT_Group == "Reference (0)"]

            # Perform Mann-Whitney U test (Wilcoxon rank-sum test)
            ref_values <- variant_data$horvath_error[variant_data$ALT_Group == "Reference (0)"]
            alt_values <- variant_data$horvath_error[variant_data$ALT_Group == "Alternative (1+)"]

            mann_whitney_result <- wilcox.test(alt_values, ref_values, exact = FALSE)

            # Add to results dataframe
            variant_stats <- rbind(variant_stats, data.frame(
                Variant = current_variant,
                CpG = current_cpg,
                beta = current_beta,
                horvath_diff = horvath_diff,
                p_value = mann_whitney_result$p.value,
                significant = mann_whitney_result$p.value < 0.05,
                ref_mean = group_means$mean_horvath_error[group_means$ALT_Group == "Reference (0)"],
                alt_mean = group_means$mean_horvath_error[group_means$ALT_Group == "Alternative (1+)"],
                ref_n = group_means$n[group_means$ALT_Group == "Reference (0)"],
                alt_n = group_means$n[group_means$ALT_Group == "Alternative (1+)"]
            ))
        }
    }

    # # Calculate correlation between beta and horvath difference
    # correlation <- cor.test(variant_stats$beta, variant_stats$horvath_diff)
    # cor_value <- round(correlation$estimate, 3)
    # p_value <- round(correlation$p.value, 4)

    # Create the scatter plot
    p <- ggplot(variant_stats, aes(x = beta, y = horvath_diff)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
        geom_point(aes(color = significant), alpha = 0.7) +
        scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
        scale_size_continuous(range = c(2, 6)) +
        labs(
            title = "Difference in Horvath Error by meQTL Variant Effect Size",
            subtitle = "Statistical significance based on Mann-Whitney U test",
            x = "Variant Effect Size (abs_effect_beta)",
            y = "Horvath Error Difference (Alternative - Reference)",
            color = "Significant\n(p < 0.05)"
        ) +
        theme_classic() +
        theme(
            legend.position = "right",
            plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(face = "italic"),
            axis.title = element_text(face = "bold")
        )

    # Label significant points or points with large differences
    # Find top variants by absolute difference or significance
    to_label <- variant_stats %>%
        filter(significant == TRUE | abs(horvath_diff) > quantile(abs(horvath_diff), 0.9)) %>%
        arrange(p_value) %>%
        head(10) # Label at most 10 points to avoid crowding

    if (nrow(to_label) > 0) {
        p <- p +
            geom_text_repel(
                data = to_label,
                aes(label = Variant),
                box.padding = 0.5,
                point.padding = 0.3,
                segment.color = "grey50",
                size = 3
            )
    }

    # Return the plot and the data
    return(list(
        plot = p,
        data = variant_stats
    ))
}

results <- plot_horvath_diff_by_beta(aa_geno_full_afr_diff_zscore, meqtl_effects, top_variants_only = F)
results_abs_diff <- plot_horvath_diff_by_beta(aa_geno_full_afr_diff_abs_diff, meqtl_effects, top_variants_only = F)
results_5_fc <- plot_horvath_diff_by_beta(aa_geno_full_afr_diff_5_fc, meqtl_effects, top_variants_only = F)


# View the plot
results$plot
results_abs_diff$plot
results_5_fc$plot

# Get the top variants by absolute difference
top_diff <- results$data %>%
    arrange(desc(abs(horvath_diff))) %>%
    head(10)
print(top_diff)

# Get the most significant variants
top_sig <- results$data %>%
    arrange(p_value) %>%
    head(10)
print(top_sig)



plot_horvath_diff_by_beta <- function(individual_data, meQTL_data, top_variants_only = F) {
    if (top_variants_only) {
        variants_to_analyze <- meQTL_data %>%
            group_by(CpG) %>%
            slice_max(abs_effect_beta, n = 1, with_ties = FALSE) %>%
            ungroup()
    } else {
        variants_to_analyze <- meQTL_data
    }

    variant_stats <- data.frame()

    for (i in 1:nrow(variants_to_analyze)) {
        current_variant <- variants_to_analyze$Variant[i]
        current_beta <- variants_to_analyze$abs_effect_beta[i]
        current_cpg <- variants_to_analyze$CpG[i]

        variant_data <- individual_data %>%
            filter(Variant == current_variant) %>%
            mutate(ALT_Group = ifelse(ALT_Copies == 0, "Reference (0)", "Alternative (1+)"))

        group_counts <- table(variant_data$ALT_Group)
        if (any(group_counts < 3)) next

        group_stats <- variant_data %>%
            group_by(ALT_Group) %>%
            summarize(
                mean_horvath_error = mean(horvath_error, na.rm = TRUE),
                median_horvath_error = median(horvath_error, na.rm = TRUE),
                sd_horvath_error = sd(horvath_error, na.rm = TRUE),
                n = n(),
                .groups = "drop"
            )

        if (nrow(group_stats) == 2) {
            ref_median <- group_stats$median_horvath_error[group_stats$ALT_Group == "Reference (0)"]
            alt_median <- group_stats$median_horvath_error[group_stats$ALT_Group == "Alternative (1+)"]
            horvath_diff_median <- alt_median - ref_median

            wilcox_test <- wilcox.test(horvath_error ~ ALT_Group, data = variant_data, exact = FALSE)

            total_n <- sum(group_stats$n)
            ref_n <- group_stats$n[group_stats$ALT_Group == "Reference (0)"]
            alt_n <- group_stats$n[group_stats$ALT_Group == "Alternative (1+)"]
            ref_prop <- ref_n / total_n
            alt_prop <- alt_n / total_n

            is_imbalanced <- ref_prop > 0.9 | alt_prop > 0.9

            variant_stats <- rbind(variant_stats, data.frame(
                Variant = current_variant,
                CpG = current_cpg,
                beta = current_beta,
                horvath_diff_median = horvath_diff_median,
                p_value = wilcox_test$p.value,
                significant = wilcox_test$p.value < 0.05,
                ref_median = ref_median,
                alt_median = alt_median,
                ref_n = ref_n,
                alt_n = alt_n,
                ref_prop = ref_prop,
                alt_prop = alt_prop,
                imbalanced = is_imbalanced
            ))
        }
    }

    # Plot
    p <- ggplot(variant_stats, aes(x = beta, y = horvath_diff_median)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
        geom_point(aes(
            color = significant,
            size = ref_n + alt_n,
            shape = imbalanced
        ), alpha = 0.7) +
        scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
        scale_size_continuous(range = c(2, 6)) +
        scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) + # Circle vs triangle
        labs(
            title = "Median Horvath Error Difference by meQTL Variant Effect Size",
            x = "Variant Effect Size (abs_effect_beta)",
            y = "Horvath Error Median Difference (Alt - Ref)",
            color = "Significant\n(p < 0.05)",
            size = "Sample Size\n(total n)",
            shape = "Imbalanced\nAllele Frequency"
        ) +
        theme_classic() +
        theme(
            legend.position = "right",
            plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold")
        )

    # Label top variants
    to_label <- variant_stats %>%
        filter(significant == TRUE | abs(horvath_diff_median) > quantile(abs(horvath_diff_median), 0.9)) %>%
        arrange(p_value) %>%
        head(10)

    if (nrow(to_label) > 0) {
        p <- p +
            geom_text_repel(
                data = to_label,
                aes(label = Variant),
                box.padding = 0.5,
                point.padding = 0.3,
                segment.color = "grey50",
                size = 3
            )
    }

    return(list(
        plot = p,
        data = variant_stats
    ))
}


results <- plot_horvath_diff_by_beta(aa_geno_full_afr_diff_zscore, strong_meqtl_effects_filtered)
results_5_fc <- plot_horvath_diff_by_beta(aa_geno_full_afr_diff_5_fc, strong_meqtl_effects_filtered)

results$plot
results_5_fc$plot

plot_horvath_error_by_meqtl_beta <- function(individual_data, meQTL_data, top_variants_only = TRUE) {
    # Get variants to use
    if (top_variants_only) {
        variants_to_analyze <- meQTL_data %>%
            group_by(CpG) %>%
            slice_max(abs_effect_beta, n = 1, with_ties = FALSE) %>%
            ungroup()
    } else {
        variants_to_analyze <- meQTL_data
    }

    # Join with individual-level data to get horvath error per variant
    merged_data <- individual_data %>%
        filter(ALT_Copies > 0) %>%
        inner_join(variants_to_analyze, by = "Variant")

    # Summarize: mean Horvath error per variant
    variant_summary <- merged_data %>%
        group_by(Variant, CGI) %>%
        summarize(
            mean_horvath_error = mean(horvath_error, na.rm = TRUE),
            sd_horvath_error = sd(horvath_error, na.rm = TRUE),
            n = n(),
            beta = first(abs_effect_beta),
            CpG = meQTL_data$CpG[match(Variant, meQTL_data$Variant)],
            .groups = "drop"
        )

    # Plot
    p <- ggplot(variant_summary, aes(x = beta, y = mean_horvath_error)) +
        geom_point(aes(size = n), alpha = 0.7, color = "steelblue") +
        scale_size_continuous(range = c(2, 6)) +
        labs(
            title = "Mean Horvath Error by meQTL Variant Effect Size",
            x = "Variant Effect on Methylation (abs_effect_beta)",
            y = "Mean Horvath Error",
            size = "Sample Size\n(n)"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold")
        )

    return(list(
        plot = p,
        data = variant_summary
    ))
}

res <- plot_horvath_error_by_meqtl_beta(aa_geno_full_afr_diff_zscore, strong_meqtl_effects_filtered, top_variants_only = F)
res
res_5_fc <- plot_horvath_error_by_meqtl_beta(aa_geno_full_afr_diff_5_fc, strong_meqtl_effects_filtered, top_variants_only = F)
res_5_fc




plot_horvath_error_by_meqtl_beta <- function(individual_data, meQTL_data, top_variants_only = TRUE) {
    # Get variants to use
    if (top_variants_only) {
        variants_to_analyze <- meQTL_data %>%
            group_by(CpG) %>%
            slice_max(abs_effect_beta, n = 1, with_ties = FALSE) %>%
            ungroup()
    } else {
        variants_to_analyze <- meQTL_data
    }

    # Join with individual-level data to get horvath error per variant
    # The key issue: we need to keep track of individual effects
    merged_data <- individual_data %>%
        filter(ALT_Copies > 0) %>%
        inner_join(variants_to_analyze, by = "Variant")

    # Create a new column that represents the actual methylation effect for each individual
    # This will be the variant effect size multiplied by the number of alternate alleles
    merged_data <- merged_data %>%
        mutate(individual_methylation_effect = abs_effect_beta * ALT_Copies)

    # Now plot the individual-level data
    p <- ggplot(merged_data, aes(x = individual_methylation_effect, y = horvath_error)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", color = "red") +
        facet_wrap(~Variant, scales = "free_x", ncol = 3) +
        labs(
            title = "Horvath Error by Individual Methylation Effect",
            x = "Individual Methylation Effect (variant effect × ALT copies)",
            y = "Horvath Error"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold"),
            strip.text = element_text(size = 8)
        )

    # Create an alternative plot showing all variants together
    p_all <- ggplot(merged_data, aes(x = individual_methylation_effect, y = horvath_error, color = Variant)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = FALSE, size = 0.7) +
        labs(
            title = "Horvath Error by Individual Methylation Effect",
            subtitle = "All variants",
            x = "Individual Methylation Effect (variant effect × ALT copies)",
            y = "Horvath Error"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(face = "bold"),
            axis.title = element_text(face = "bold")
        )

    # If there are too many variants, simplify the legend in the all-variants plot
    if (length(unique(merged_data$Variant)) > 10) {
        p_all <- p_all + guides(color = guide_legend(ncol = 2))
    }

    return(list(
        plot_faceted = p,
        plot_all = p_all,
        data = merged_data
    ))
}

res <- plot_horvath_error_by_meqtl_beta(aa_geno_full_afr_diff_zscore, strong_meqtl_effects_filtered, top_variants_only = F)
res
res_5_fc <- plot_horvath_error_by_meqtl_beta(aa_geno_full_afr_diff_5_fc, strong_meqtl_effects_filtered, top_variants_only = F)
res_5_fc

# meQTL Horvath Clock Accuracy Analysis
# Comparing African ancestry differentiated vs non-differentiated meQTLs

library(MatchIt)

# Function to prepare meQTL data from two separate dataframes - updated for integrated structure
prepare_meqtl_data <- function(aa_geno_full, aa_geno_full_afr_diff) {
    # Standardize column names if needed
    if ("Variant" %in% colnames(aa_geno_full) && !"variant_id" %in% colnames(aa_geno_full)) {
        aa_geno_full$variant_id <- aa_geno_full$Variant
    }
    if ("Sample" %in% colnames(aa_geno_full) && !"sample_id" %in% colnames(aa_geno_full)) {
        aa_geno_full$sample_id <- aa_geno_full$Sample
    }

    # Do the same for the African differentiated dataset
    if ("Variant" %in% colnames(aa_geno_full_afr_diff) && !"variant_id" %in% colnames(aa_geno_full_afr_diff)) {
        aa_geno_full_afr_diff$variant_id <- aa_geno_full_afr_diff$Variant
    }
    if ("Sample" %in% colnames(aa_geno_full_afr_diff) && !"sample_id" %in% colnames(aa_geno_full_afr_diff)) {
        aa_geno_full_afr_diff$sample_id <- aa_geno_full_afr_diff$Sample
    }

    # Add african_differentiated flag to full dataset
    aa_geno_full$african_differentiated <- FALSE

    # Mark African differentiated variants
    if (nrow(aa_geno_full_afr_diff) > 0) {
        # Find matching variants
        common_variants <- intersect(aa_geno_full$variant_id, aa_geno_full_afr_diff$variant_id)
        aa_geno_full$african_differentiated[aa_geno_full$variant_id %in% common_variants] <- TRUE

        cat("Total variant-sample combinations in full dataset:", nrow(aa_geno_full), "\n")
        cat(
            "African differentiated variant-sample combinations:",
            sum(aa_geno_full$african_differentiated), "\n"
        )
        cat("Unique African differentiated variants:", length(common_variants), "\n")
        cat("Total unique variants:", length(unique(aa_geno_full$variant_id)), "\n")
    }

    return(aa_geno_full)
}

# Function to perform propensity score matching on beta values - updated for integrated data
match_meqtls_by_beta <- function(meqtl_data, caliper = 0.1) {
    # Get unique variants with their properties for matching
    variant_summary <- meqtl_data %>%
        group_by(variant_id, african_differentiated, BETA) %>%
        summarise(n_samples = n(), .groups = "drop")

    # Create a binary treatment variable (1 = African differentiated, 0 = not)
    variant_summary$african_diff_binary <- ifelse(variant_summary$african_differentiated, 1, 0)

    # Check if we have both groups
    if (sum(variant_summary$african_diff_binary == 1) == 0) {
        stop("No African differentiated variants found in the data")
    }
    if (sum(variant_summary$african_diff_binary == 0) == 0) {
        stop("No non-differentiated variants found in the data")
    }

    cat("Pre-matching summary:\n")
    cat("African differentiated variants:", sum(variant_summary$african_diff_binary == 1), "\n")
    cat("Non-differentiated variants:", sum(variant_summary$african_diff_binary == 0), "\n")

    # Perform matching based on absolute beta values
    match_obj <- matchit(african_diff_binary ~ abs(BETA),
        data = variant_summary,
        method = "nearest",
        caliper = caliper,
        std.caliper = TRUE
    )

    # Extract matched variant IDs
    matched_variants <- match.data(match_obj)

    # Filter original data to include only matched variants
    matched_data <- meqtl_data[meqtl_data$variant_id %in% matched_variants$variant_id, ]

    # Print matching summary
    cat("\nMatching Summary:\n")
    print(summary(match_obj))

    return(matched_data)
}

# Function to calculate clock error effects by genotype - updated for integrated data
calculate_clock_effects <- function(meqtl_data_combined) {
    results <- data.frame()

    # Get unique variants from matched data
    unique_variants <- unique(meqtl_data_combined[, c("variant_id", "african_differentiated", "BETA")])

    for (i in 1:nrow(unique_variants)) {
        variant_id <- unique_variants$variant_id[i]
        african_diff <- unique_variants$african_differentiated[i]
        beta <- unique_variants$BETA[i]

        # Get data for this specific variant
        variant_data <- meqtl_data_combined[meqtl_data_combined$variant_id == variant_id, ]

        if (nrow(variant_data) == 0) next

        # Recode genotype: 0 = reference (0 copies), 1 = variant (1+ copies)
        variant_data$genotype_binary <- ifelse(variant_data$ALT_Copies >= 1, 1, 0)

        # Calculate median clock errors by genotype
        summary_stats <- variant_data %>%
            group_by(genotype_binary) %>%
            summarise(
                median_error = median(horvath_error, na.rm = TRUE),
                mean_error = mean(horvath_error, na.rm = TRUE),
                n = n(),
                .groups = "drop"
            )

        # Add variant information
        summary_stats$variant_id <- variant_id
        summary_stats$african_differentiated <- african_diff
        summary_stats$BETA <- beta

        results <- rbind(results, summary_stats)
    }

    return(results)
}

# Function to prepare data for plotting
prepare_plot_data <- function(clock_effects) {
    # Calculate total N per variant
    variant_totals <- clock_effects %>%
        group_by(variant_id) %>%
        summarise(total_n = sum(n), .groups = "drop")

    # Calculate difference in median error (variant - reference)
    effect_diff <- clock_effects %>%
        group_by(variant_id, african_differentiated, BETA) %>%
        summarise(
            error_difference = median_error[genotype_binary == 1] - median_error[genotype_binary == 0],
            ref_error = median_error[genotype_binary == 0],
            var_error = median_error[genotype_binary == 1],
            .groups = "drop"
        )

    # Create long format for boxplot
    plot_data <- clock_effects %>%
        filter(genotype_binary == 1) %>% # Only include variant carriers (1+ copies)
        left_join(variant_totals, by = "variant_id") %>%
        mutate(
            frequency = n / total_n,
            genotype_label = "Variant Carriers (1+)",
            ancestry_group = ifelse(african_differentiated,
                "African Ancestry Differentiated",
                "Non-Differentiated"
            )
        )

    return(list(effect_diff = effect_diff, plot_data = plot_data))
}

# Main analysis function - updated for integrated data structure
run_meqtl_clock_analysis <- function(aa_geno_full, aa_geno_full_afr_diff) {
    cat("Step 1: Preparing meQTL data from two dataframes...\n")
    meqtl_data <- prepare_meqtl_data(aa_geno_full, aa_geno_full_afr_diff)

    cat("\nStep 2: Matching meQTLs by beta values...\n")
    matched_meqtls <- match_meqtls_by_beta(meqtl_data)

    cat("\nStep 3: Calculating clock effects for matched meQTLs...\n")
    clock_effects <- calculate_clock_effects(matched_meqtls)

    cat("\nStep 4: Preparing data for visualization...\n")
    plot_data_list <- prepare_plot_data(clock_effects)

    # Statistical comparison
    cat("\nStep 5: Statistical comparison...\n")
    african_diff_effects <- plot_data_list$effect_diff$error_difference[
        plot_data_list$effect_diff$african_differentiated == TRUE
    ]
    non_diff_effects <- plot_data_list$effect_diff$error_difference[
        plot_data_list$effect_diff$african_differentiated == FALSE
    ]

    # Wilcoxon test for difference in median effects
    if (length(african_diff_effects) > 0 & length(non_diff_effects) > 0) {
        wilcox_test <- wilcox.test(african_diff_effects, non_diff_effects)

        cat(
            "Median effect difference (African differentiated):",
            median(african_diff_effects, na.rm = TRUE), "\n"
        )
        cat(
            "Median effect difference (Non-differentiated):",
            median(non_diff_effects, na.rm = TRUE), "\n"
        )
        cat("Wilcoxon test p-value:", wilcox_test$p.value, "\n")
    } else {
        wilcox_test <- NULL
        cat("Warning: Insufficient data for statistical comparison\n")
    }

    return(list(
        meqtl_data = meqtl_data,
        matched_meqtls = matched_meqtls,
        clock_effects = clock_effects,
        plot_data = plot_data_list$plot_data,
        effect_differences = plot_data_list$effect_diff,
        statistical_test = wilcox_test
    ))
}

# Create comparative boxplot
create_comparative_boxplot <- function(results) {
    # Main boxplot comparing clock errors by genotype and ancestry group
    p1 <- ggplot(
        results$plot_data,
        aes(x = genotype_label, y = median_error, fill = ancestry_group)
    ) +
        geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        scale_fill_manual(values = c(
            "African Ancestry Differentiated" = "#E74C3C",
            "Non-Differentiated" = "#3498DB"
        )) +
        ylim(0, 12) +
        labs(
            title = "Horvath Clock Error by meQTL Genotype and Ancestry Differentiation",
            x = "Genotype",
            y = "Horvath Clock Error",
            fill = "meQTL Type",
            subtitle = paste(
                "Matched meQTLs (n =",
                length(unique(results$plot_data$variant_id)), "pairs)"
            )
        ) +
        theme_classic() +
        theme(
            legend.position = "bottom",
            plot.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12)
        ) +
        facet_wrap(~ancestry_group, scales = "free_y")

    # Effect size comparison plot
    p2 <- ggplot(
        results$effect_differences,
        aes(
            x = african_differentiated, y = error_difference,
            fill = african_differentiated
        )
    ) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.5) +
        scale_fill_manual(values = c("FALSE" = "#3498DB", "TRUE" = "#E74C3C")) +
        scale_x_discrete(labels = c(
            "FALSE" = "Non-Differentiated",
            "TRUE" = "African Differentiated"
        )) +
        ylim(0, 12) +
        labs(
            title = "Effect Size Comparison: Change in Clock Error",
            x = "meQTL Type",
            y = "Difference in Median Clock Error\n(Variant - Reference)",
            subtitle = "Positive values = increased error with variant"
        ) +
        theme_classic() +
        theme(
            legend.position = "none",
            plot.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12)
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")

    # Combine plots
    combined_plot <- grid.arrange(p1, p2, nrow = 2, heights = c(2, 1))

    return(list(main_plot = p1, effect_plot = p2, combined = combined_plot))
}

# Clean up the OG dataframe with ALL variants so that we can use it with the function
# aa_geno_full <- aa_geno_full %>%
#   mutate(Variant = gsub("_\\d+$", "", Variant))
# # Remove NA rows (specifically present in the BETA column), otherwise the matching function doesn't work
# aa_geno_full <- aa_geno_full %>% drop_na()
# # Add the horvath error to the OG dataframe
# aa_geno_full <- aa_geno_full %>%
#   left_join(
#     aa_geno_full_afr_diff_zscore %>%
#       select(Sample, horvath_error) %>%
#       distinct(Sample, .keep_all = TRUE),
#     by = "Sample"
#   )

# To run with your actual data, use:
results <- run_meqtl_clock_analysis(
    aa_geno_full = aa_geno_full_afr_diff_zscore, # Your full meQTL dataframe
    aa_geno_full_afr_diff = aa_geno_full_afr_diff_5_fc # Your African differentiated subset
)



# Create plots
plots <- create_comparative_boxplot(results)

# Display the main comparative plot
print(plots$main_plot)


######### Run the above pipeline for Hispanics and then Puerto Ricans
# Modify the variant IDs so that we can match them to the AFR-diff meQTL variant IDs
hisp_geno_full <- hisp_geno_full %>%
    mutate(Variant = gsub("_\\d+$", "", Variant))
pr_geno_full <- pr_geno_full %>%
    mutate(Variant = gsub("_\\d+$", "", Variant))

# Join the age predictions to the genotyping data
hisp_geno_full <- hisp_geno_full %>%
    inner_join(age_preds, by = c("Sample" = "CGI"))
pr_geno_full <- pr_geno_full %>%
    inner_join(age_preds, by = c("Sample" = "CGI"))

# Calculate the Horvath clock's error
hisp_geno_full$horvath_error <- abs(hisp_geno_full$Horvath - hisp_geno_full$AGE_OF_EXAM)
pr_geno_full$horvath_error <- abs(pr_geno_full$Horvath - pr_geno_full$AGE_OF_EXAM)

# Filter the data to create dataframes with just the AFR-differentiated variants in Hispanics and Puerto Ricans
hisp_geno_full_afr_5_fc <- hisp_geno_full[hisp_geno_full$Variant %in% aa_geno_full_afr_diff_5_fc$Variant, ]
pr_geno_full_afr_5_fc <- pr_geno_full[pr_geno_full$Variant %in% aa_geno_full_afr_diff_5_fc$Variant, ]

# Drop NAs
hisp_geno_full <- hisp_geno_full[!is.na(hisp_geno_full$BETA), ]
hisp_geno_full_afr_5_fc <- hisp_geno_full_afr_5_fc[!is.na(hisp_geno_full_afr_5_fc$BETA), ]

pr_geno_full <- pr_geno_full[!is.na(pr_geno_full$BETA), ]
pr_geno_full_afr_5_fc <- pr_geno_full_afr_5_fc[!is.na(pr_geno_full_afr_5_fc$BETA), ]

# Filter out the Puerto Ricans from the Hispanics dataframes
hisp_geno_full <- hisp_geno_full %>%
    filter(COHORT != "PRADI")
hisp_geno_full_afr_5_fc <- hisp_geno_full_afr_5_fc %>%
    filter(COHORT != "PRADI")

# Generate dataframes for Cubans and Peruvians separately
cuban_geno_full <- hisp_geno_full %>%
    filter(COHORT == "CuADI")
cuban_geno_full_afr_5_fc <- hisp_geno_full_afr_5_fc %>%
    filter(COHORT == "CuADI")

peruvian_geno_full <- hisp_geno_full %>%
    filter(COHORT == "PERUVIAN")
peruvian_geno_full_afr_5_fc <- hisp_geno_full_afr_5_fc %>%
    filter(COHORT == "PERUVIAN")

# Run the pipeline
results_hispanics <- run_meqtl_clock_analysis(
    aa_geno_full = hisp_geno_full, # Your full meQTL dataframe
    aa_geno_full_afr_diff = hisp_geno_full_afr_5_fc # Your African differentiated subset
)

results_pr <- run_meqtl_clock_analysis(
    aa_geno_full = pr_geno_full, # Your full meQTL dataframe
    aa_geno_full_afr_diff = pr_geno_full_afr_5_fc # Your African differentiated subset
)

results_cub <- run_meqtl_clock_analysis(
    aa_geno_full = cuban_geno_full,
    aa_geno_full_afr_diff = cuban_geno_full_afr_5_fc
)

results_per <- run_meqtl_clock_analysis(
    aa_geno_full = peruvian_geno_full,
    aa_geno_full_afr_diff = peruvian_geno_full_afr_5_fc
)

# Create plots
plots_hispanics <- create_comparative_boxplot(results_hispanics)
plots_pr <- create_comparative_boxplot(results_pr)
plots_cub <- create_comparative_boxplot(results_cub)
plots_per <- create_comparative_boxplot(results_per)


# Display the main comparative plot
print(plots_hispanics$main_plot)
print(plots_pr$main_plot)
print(plots_cub$main_plot)
print(plots_per$main_plot)

cohort_data <- list(
    "AFR" = list(aa_geno_full = aa_geno_full_afr_diff_zscore, aa_geno_full_afr_diff = aa_geno_full_afr_diff_5_fc),
    "PR" = list(aa_geno_full = pr_geno_full, aa_geno_full_afr_diff = pr_geno_full_afr_5_fc),
    "CUB" = list(aa_geno_full = cuban_geno_full, aa_geno_full_afr_diff = cuban_geno_full_afr_5_fc),
    "PER" = list(aa_geno_full = peruvian_geno_full, aa_geno_full_afr_diff = peruvian_geno_full_afr_5_fc)
)


# Multi-cohort analysis function with error handling
run_multi_cohort_meqtl_analysis <- function(cohort_data_list) {
    # cohort_data_list should be a named list with elements like:
    # list("Cohort1" = list(aa_geno_full = df1, aa_geno_full_afr_diff = df2),
    #      "Cohort2" = list(aa_geno_full = df3, aa_geno_full_afr_diff = df4), ...)

    all_results <- list()
    successful_cohorts <- character()
    failed_cohorts <- character()

    cat("Running analysis for", length(cohort_data_list), "cohorts...\n")

    # Run analysis for each cohort with error handling
    for (cohort_name in names(cohort_data_list)) {
        cat("\nAnalyzing cohort:", cohort_name, "\n")

        cohort_data <- cohort_data_list[[cohort_name]]

        # Check data availability before running analysis
        tryCatch(
            {
                # Quick check for AFR differentiated variants
                cat("  Checking data structure...\n")
                cat("  aa_geno_full rows:", nrow(cohort_data$aa_geno_full), "\n")
                cat("  aa_geno_full_afr_diff rows:", nrow(cohort_data$aa_geno_full_afr_diff), "\n")

                # Run your existing analysis function for this cohort
                cohort_results <- run_meqtl_clock_analysis(
                    cohort_data$aa_geno_full,
                    cohort_data$aa_geno_full_afr_diff
                )

                # Add cohort identifier to the results
                if (!is.null(cohort_results$plot_data)) {
                    cohort_results$plot_data$cohort <- cohort_name
                }
                if (!is.null(cohort_results$effect_differences)) {
                    cohort_results$effect_differences$cohort <- cohort_name
                }

                all_results[[cohort_name]] <- cohort_results
                successful_cohorts <- c(successful_cohorts, cohort_name)
                cat("  ✓ Analysis completed successfully\n")
            },
            error = function(e) {
                cat("  ✗ Error in cohort", cohort_name, ":", e$message, "\n")
                failed_cohorts <- c(failed_cohorts, cohort_name)

                # Create empty results structure for failed cohort
                all_results[[cohort_name]] <<- list(
                    meqtl_data = NULL,
                    matched_meqtls = NULL,
                    clock_effects = NULL,
                    plot_data = NULL,
                    effect_differences = NULL,
                    statistical_test = NULL,
                    error = e$message
                )
            }
        )
    }

    cat("\n=== ANALYSIS SUMMARY ===\n")
    cat("Successful cohorts:", length(successful_cohorts), "-", paste(successful_cohorts, collapse = ", "), "\n")
    if (length(failed_cohorts) > 0) {
        cat("Failed cohorts:", length(failed_cohorts), "-", paste(failed_cohorts, collapse = ", "), "\n")
    }

    return(list(
        results = all_results,
        successful_cohorts = successful_cohorts,
        failed_cohorts = failed_cohorts
    ))
}

# Function to combine data from all cohorts and calculate statistics
prepare_multi_cohort_plot_data <- function(analysis_results) {
    all_results <- analysis_results$results
    successful_cohorts <- analysis_results$successful_cohorts

    # Only use successful cohorts
    successful_results <- all_results[successful_cohorts]

    if (length(successful_results) == 0) {
        stop("No cohorts completed analysis successfully")
    }

    # Combine effect differences from successful cohorts only
    effect_dfs <- list()
    for (cohort_name in names(successful_results)) {
        if (!is.null(successful_results[[cohort_name]]$effect_differences)) {
            effect_dfs[[cohort_name]] <- successful_results[[cohort_name]]$effect_differences
        }
    }

    if (length(effect_dfs) == 0) {
        stop("No effect difference data available from any cohort")
    }

    combined_effects <- do.call(rbind, effect_dfs)

    # Calculate summary statistics for each cohort and differentiation status
    summary_stats <- combined_effects %>%
        group_by(cohort, african_differentiated) %>%
        summarise(
            median_effect = median(error_difference, na.rm = TRUE),
            q25 = quantile(error_difference, 0.25, na.rm = TRUE),
            q75 = quantile(error_difference, 0.75, na.rm = TRUE),
            n = n(),
            .groups = "drop"
        )

    # Calculate statistical tests for each successful cohort
    stat_tests <- list()
    for (cohort_name in names(successful_results)) {
        cohort_effects <- combined_effects[combined_effects$cohort == cohort_name, ]

        african_diff_effects <- cohort_effects$error_difference[
            cohort_effects$african_differentiated == TRUE
        ]
        non_diff_effects <- cohort_effects$error_difference[
            cohort_effects$african_differentiated == FALSE
        ]

        if (length(african_diff_effects) > 2 & length(non_diff_effects) > 2) {
            wilcox_result <- wilcox.test(african_diff_effects, non_diff_effects)
            stat_tests[[cohort_name]] <- list(
                p_value = wilcox_result$p.value,
                significant = wilcox_result$p.value < 0.05,
                median_diff = median(african_diff_effects, na.rm = TRUE),
                median_non_diff = median(non_diff_effects, na.rm = TRUE),
                n_diff = length(african_diff_effects),
                n_non_diff = length(non_diff_effects)
            )
        } else {
            stat_tests[[cohort_name]] <- list(
                p_value = NA,
                significant = FALSE,
                median_diff = ifelse(length(african_diff_effects) > 0,
                    median(african_diff_effects, na.rm = TRUE), NA
                ),
                median_non_diff = ifelse(length(non_diff_effects) > 0,
                    median(non_diff_effects, na.rm = TRUE), NA
                ),
                n_diff = length(african_diff_effects),
                n_non_diff = length(non_diff_effects)
            )
        }
    }

    return(list(
        combined_effects = combined_effects,
        summary_stats = summary_stats,
        statistical_tests = stat_tests,
        successful_cohorts = successful_cohorts
    ))
}

# Function to create the multi-cohort comparison plot
create_multi_cohort_plot <- function(analysis_results) {
    # Prepare combined data
    plot_data_prep <- prepare_multi_cohort_plot_data(analysis_results)
    combined_effects <- plot_data_prep$combined_effects
    stat_tests <- plot_data_prep$statistical_tests
    successful_cohorts <- plot_data_prep$successful_cohorts

    if (nrow(combined_effects) == 0) {
        stop("No data available for plotting")
    }

    # Create labels for differentiation status
    combined_effects$diff_label <- ifelse(
        combined_effects$african_differentiated,
        "AFR Differentiated",
        "Non-Differentiated"
    )

    # Reorder cohorts to only include successful ones
    combined_effects$cohort <- factor(combined_effects$cohort,
        levels = successful_cohorts
    )

    # Create the main plot
    p <- ggplot(
        combined_effects,
        aes(x = cohort, y = error_difference, fill = diff_label)
    ) +
        geom_boxplot(
            position = position_dodge(width = 0.5),
            alpha = 0.7, outlier.alpha = 0.5
        ) +
        geom_jitter(
            position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
            alpha = 0.4, size = 0.8
        ) +
        scale_fill_manual(values = c(
            "AFR Differentiated" = "#E74C3C",
            "Non-Differentiated" = "#3498DB"
        )) +
        labs(
            x = "Cohort",
            y = "Difference in Median Clock Error\n(Variant - Reference)",
            fill = "meQTL Type"
        ) +
        theme_classic() +
        theme(
            legend.position = "bottom",
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10)
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7)

    # Add significance annotations
    sig_annotations <- data.frame()

    # Calculate a consistent y-position for all annotations
    data_range <- max(combined_effects$error_difference, na.rm = TRUE) - min(combined_effects$error_difference, na.rm = TRUE)
    annotation_y <- max(combined_effects$error_difference, na.rm = TRUE) + data_range * 0.2

    for (i in seq_along(stat_tests)) {
        cohort_name <- names(stat_tests)[i]
        test_result <- stat_tests[[cohort_name]]

        if (!is.na(test_result$p_value)) {
            # Determine significance symbol
            if (test_result$p_value < 0.001) {
                sig_symbol <- "***"
            } else if (test_result$p_value < 0.01) {
                sig_symbol <- "**"
            } else if (test_result$p_value < 0.05) {
                sig_symbol <- "*"
            } else {
                sig_symbol <- "ns"
            }

            sig_annotations <- rbind(sig_annotations, data.frame(
                cohort = cohort_name,
                y_pos = annotation_y,
                label = paste0(
                    "p = ", format(test_result$p_value, digits = 3),
                    " (", sig_symbol, ")"
                ),
                stringsAsFactors = FALSE
            ))
        }
    }

    # Add significance annotations to plot
    if (nrow(sig_annotations) > 0) {
        p <- p +
            geom_text(
                data = sig_annotations,
                aes(x = cohort, y = y_pos, label = label),
                inherit.aes = FALSE, size = 4, hjust = 0.5
            ) +
            # Add brackets for significance
            geom_segment(
                data = sig_annotations,
                aes(
                    x = as.numeric(factor(cohort)) - 0.25,
                    xend = as.numeric(factor(cohort)) + 0.25,
                    y = y_pos - 0.4, yend = y_pos - 0.4
                ),
                inherit.aes = FALSE, color = "black"
            )
    }

    return(list(
        plot = p,
        statistical_tests = stat_tests,
        plot_data = combined_effects,
        successful_cohorts = successful_cohorts
    ))
}
# Diagnostic function to help troubleshoot data issues
diagnose_cohort_data <- function(cohort_data_list) {
    cat("=== DATA DIAGNOSTIC ===\n")

    for (cohort_name in names(cohort_data_list)) {
        cat("\nCohort:", cohort_name, "\n")
        cohort_data <- cohort_data_list[[cohort_name]]

        # Check basic structure
        cat("  aa_geno_full: ")
        if (is.null(cohort_data$aa_geno_full)) {
            cat("NULL\n")
        } else {
            cat(nrow(cohort_data$aa_geno_full), "rows,", ncol(cohort_data$aa_geno_full), "cols\n")
            cat("    Column names:", paste(head(colnames(cohort_data$aa_geno_full)), collapse = ", "), "...\n")
        }

        cat("  aa_geno_full_afr_diff: ")
        if (is.null(cohort_data$aa_geno_full_afr_diff)) {
            cat("NULL\n")
        } else {
            cat(nrow(cohort_data$aa_geno_full_afr_diff), "rows,", ncol(cohort_data$aa_geno_full_afr_diff), "cols\n")
            if (nrow(cohort_data$aa_geno_full_afr_diff) == 0) {
                cat("    *** WARNING: No AFR differentiated variants in this cohort! ***\n")
            }
        }

        # Check if required functions exist
        cat("  Functions available: ")
        required_funcs <- c(
            "prepare_meqtl_data", "match_meqtls_by_beta",
            "calculate_clock_effects", "prepare_plot_data"
        )
        available_funcs <- sapply(required_funcs, exists)
        cat(paste(names(available_funcs)[available_funcs], collapse = ", "), "\n")
        if (any(!available_funcs)) {
            cat("    *** MISSING:", paste(names(available_funcs)[!available_funcs], collapse = ", "), "***\n")
        }
    }
}


diagnose_cohort_data(cohort_data)
all_results <- run_multi_cohort_meqtl_analysis(cohort_data)
# saveRDS(all_results, file = "all_results.rds")
plot_results <- create_multi_cohort_plot(all_results)
print(plot_results$plot)

pr_eda <- all_results$results$PR$plot_data
cub_eda <- all_results$results$CUB$plot_data
per_eda <- all_results$results$PER$plot_data
aa_eda <- all_results$results$AFR$plot_data

ggplot(pr_eda, aes(x = n)) +
    geom_density(fill = "steelblue", alpha = 0.7) +
    xlim(0, 300) +
    labs(
        title = "Density Distribution of Number of Individuals with a Variant (PR)",
        x = "n",
        y = "Density"
    ) +
    theme_classic()
ggplot(aa_eda, aes(x = n)) +
    geom_density(fill = "steelblue", alpha = 0.7) +
    xlim(0, 300) +
    labs(
        title = "Density Distribution of Number of Individuals with a Variant (AA)",
        x = "n",
        y = "Density"
    ) +
    theme_classic()
ggplot(cub_eda, aes(x = n)) +
    geom_density(fill = "steelblue", alpha = 0.7) +
    xlim(0, 300) +
    labs(
        title = "Density Distribution of Number of Individuals with a Variant (CUB)",
        x = "n",
        y = "Density"
    ) +
    theme_classic()
ggplot(per_eda, aes(x = n)) +
    geom_density(fill = "steelblue", alpha = 0.7) +
    xlim(0, 300) +
    labs(
        title = "Density Distribution of Number of Individuals with a Variant (PER)",
        x = "n",
        y = "Density"
    ) +
    theme_classic()


ggplot(cub_eda, aes(x = african_differentiated, y = n)) +
    geom_jitter(width = 0.2, alpha = 0.6, color = "steelblue") +
    geom_boxplot(alpha = 0.3, outlier.shape = NA) +
    labs(
        x = "African Differentiated",
        y = "Number of Individuals with Variant"
    ) +
    theme_classic()
ggplot(per_eda, aes(x = african_differentiated, y = n)) +
    geom_jitter(width = 0.2, alpha = 0.6, color = "steelblue") +
    geom_boxplot(alpha = 0.3, outlier.shape = NA) +
    labs(
        x = "African Differentiated",
        y = "Number of Individuals with Variant"
    ) +
    theme_classic()

# Combine all dataframes with a dataset identifier
combined_data <- bind_rows(
    aa_eda %>% mutate(dataset = "AA"),
    pr_eda %>% mutate(dataset = "PR"),
    cub_eda %>% mutate(dataset = "CUB"),
    per_eda %>% mutate(dataset = "PER")
)

# Set factor levels to control the order
combined_data$dataset <- factor(combined_data$dataset, levels = c("AA", "PR", "CUB", "PER"))

# Create the combined plot
ggplot(combined_data, aes(x = african_differentiated, y = n)) +
    geom_jitter(width = 0.2, alpha = 0.6, aes(color = dataset)) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA) +
    facet_wrap(~dataset, nrow = 1) +
    scale_color_manual(values = c(
        "AA" = "coral", "PR" = "purple",
        "CUB" = "steelblue", "PER" = "forestgreen"
    )) +
    labs(
        x = "African Differentiated",
        y = "Number of Individuals with Variant",
        color = "Dataset"
    ) +
    theme_classic() +
    theme(
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(face = "bold"),
        legend.position = "none"
    )

# Subset the data
combined_true <- combined_data %>% filter(african_differentiated == TRUE)
combined_false <- combined_data %>% filter(african_differentiated == FALSE)

# Define a plotting function to avoid duplication
make_plot <- function(data, title) {
    ggplot(data, aes(x = dataset, y = frequency)) +
        geom_jitter(width = 0.2, alpha = 0.6, aes(color = dataset)) +
        geom_boxplot(alpha = 0.3, outlier.shape = NA) +
        scale_color_manual(values = c(
            "AA" = "coral", "PR" = "purple",
            "CUB" = "steelblue", "PER" = "forestgreen"
        )) +
        labs(
            x = "Dataset",
            y = "Variant Frequency",
            color = "Dataset",
            title = title
        ) +
        theme_pubr(base_size = 16) +
        theme(
            strip.background = element_rect(fill = "lightgray"),
            strip.text = element_text(face = "bold"),
            legend.position = "none"
        )
}

# Generate the two plots
plot_true <- make_plot(combined_true, "African Differentiated Clock CpG meQTLs")
plot_false <- make_plot(combined_false, "Non-Differentiated Clock CpG meQTLs")

# Display them side by side if desired
library(patchwork)
fig5e <- plot_false + plot_true
saveRDS(fig5e, "plots_rds/fig5e.rds")

# Prepare datasets
datasets <- list(
    "AA_EDA" = aa_eda,
    "PR_EDA" = pr_eda,
    "CUB_EDA" = cub_eda,
    "PER_EDA" = per_eda
)

# Add dataset identifier to each dataframe
datasets_combined <- datasets %>%
    imap(~ mutate(.x, dataset = .y)) %>%
    bind_rows()

# Calculate common y-limit across all datasets
# max_count <- datasets_combined %>%
#   group_by(dataset, african_differentiated) %>%
#   summarise(count = n(), .groups = 'drop') %>%
#   pull(count) %>%
#   max()
#
# y_limit <- ceiling(max_count * 1.1)  # Add 10% padding

# OPTION 1: Individual plots with same y-limits
create_individual_plots <- function() {
    plots <- list()

    for (name in names(datasets)) {
        cat("\n=== Analysis for", name, "===\n")

        p <- ggplot(datasets[[name]], aes(x = n, fill = african_differentiated)) +
            geom_histogram(bins = 30, alpha = 0.8) +
            facet_wrap(~african_differentiated,
                labeller = labeller(african_differentiated = c(
                    "FALSE" = "Non-Differentiated",
                    "TRUE" = "African Differentiated"
                ))
            ) +
            scale_fill_manual(values = c("FALSE" = "skyblue", "TRUE" = "coral")) +
            scale_x_continuous(limits = c(0, 300)) +
            scale_y_continuous(limits = c(0, 200)) +
            labs(
                title = paste("Distribution of n -", name),
                x = "Sample Size (n)",
                y = "Count"
            ) +
            theme_classic() +
            theme(legend.position = "none")

        plots[[name]] <- p
        print(p)

        # Summary statistics
        summary_stats <- datasets[[name]] %>%
            group_by(african_differentiated) %>%
            summarise(
                count = n(),
                mean_n = mean(n, na.rm = TRUE),
                median_n = median(n, na.rm = TRUE),
                sd_n = sd(n, na.rm = TRUE),
                min_n = min(n, na.rm = TRUE),
                max_n = max(n, na.rm = TRUE),
                .groups = "drop"
            ) %>%
            mutate(african_differentiated = ifelse(african_differentiated,
                "African Differentiated",
                "Non-Differentiated"
            ))

        print(summary_stats)
        cat("\n")
    }

    return(plots)
}

# OPTION 2: All datasets plotted together
create_combined_plot <- function() {
    cat("\n=== Combined Plot for All Datasets ===\n")

    # Set cohort factor levels in desired order
    datasets_combined$cohort <- factor(datasets_combined$cohort, levels = c("AFR", "PR", "CUB", "PER"))

    combined_plot <- ggplot(datasets_combined, aes(x = n, fill = african_differentiated)) +
        geom_histogram(bins = 30, alpha = 0.8) +
        facet_grid(african_differentiated ~ cohort,
            labeller = labeller(african_differentiated = c(
                "FALSE" = "Non-Differentiated",
                "TRUE" = "African Differentiated"
            ))
        ) +
        scale_fill_manual(values = c("FALSE" = "skyblue", "TRUE" = "coral")) +
        scale_x_continuous(limits = c(0, 300)) +
        scale_y_continuous(limits = c(0, 75)) +
        labs(
            x = "Number of Individuals With Variant",
            y = "Variant Count"
        ) +
        theme_classic() +
        theme(legend.position = "none")

    print(combined_plot)

    return(combined_plot)
}


# Run individual plots with consistent y-limits
individual_plots <- create_individual_plots()

# Run combined plot
combined_plot <- create_combined_plot()
