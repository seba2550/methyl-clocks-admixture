# This is a script to investigate the Horvath Clock CpGs and sets of meQTLs
# (From Hawe et al., 2022: https://doi.org/10.1038/s41588-021-00969-x EUR and SAS individuals
# Shang et al., 2023: https://doi.org/10.1038/s41467-023-37961-4 African American individuals
# and Vilicaña et al., 2023: https://doi.org/10.1186/s13059-023-03011-x) EUR individuals from the UK

setwd("~/Desktop/Capra Lab/Thesis Project/Aim_1")

# Load some libraries
library(tidyverse)
library(data.table)
library(UpSetR)
library(cowplot)
library(ggpubr)
library(DunedinPACE)
library(methylclock)

# Load the Horvath Clock CpGs
horvath <- read_csv("Horvath_clock_coefficients.csv")

# Load the COSMO pairs data (meQTLs)
# load("cosmopairs_combined_151216.RData")

# Find overlaps between CpGs with meQTLs and the Horvath Clock CpGs
# clock_cpg_cosmo_meqtls <- cosmo %>% inner_join(horvath, join_by(cpg == CpGmarker))

# Save the data and load it in R (for ease of analysis, given that the COSMO pairs file is huge)
# write_csv(clock_cpg_cosmo_meqtls, "clock_cpg_cosmo_meqtls.csv")
clock_cpg_cosmo_meqtls <- read_csv("clock_cpg_cosmo_meqtls.csv")

# Remove the COSMO pairs from memory
# rm(cosmo)

# Load the GENOA meQTLs
# genoa_meqtls <- fread("./GENOA_meQTL_summary_stat_allchr/merged_data.csv")
# genoa_meqtls <- genoa_meqtls %>% filter(p_wald < 0.05)

# Find overlaps between CpGs with GENOA meQTLs and the Horvath Clock CpGs
# clock_cpg_genoa_meqtls <- genoa_meqtls %>% inner_join(horvath, join_by(CpG == CpGmarker))

# Save the data and load it in R (for ease of analysis, given that the GENOA meQTLs file is huge)
# write_csv(clock_cpg_genoa_meqtls, "clock_cpg_genoa_meqtls.csv")
clock_cpg_genoa_meqtls <- read_csv("clock_cpg_genoa_meqtls.csv")

# Correct for multiple testing, and then filter the GENOA meQTLs
clock_cpg_genoa_meqtls <- clock_cpg_genoa_meqtls %>% filter(p.adjust(p_wald, method = "BH") < 0.05)

# Remove the genoa meQTLs from memory
# rm(genoa_meqtls)

# Load the set of Horvath clock meQTLs from the 3 UK cohorts
clock_cpg_uk_meqtls <- read_csv("horvath_cpgs_uk_meqtls.csv")
# Split the Top.SNP column to get the REF and ALT allele
clock_cpg_uk_meqtls <- clock_cpg_uk_meqtls %>% separate(Top.SNP, c("SNP.ID", "REF", "ALT"), sep = "_", remove = FALSE)

# Join all three sets of clock-influencing meQTLs
# Rename the columns to match, and then bind by rows
clock_cpg_meqtls <- bind_rows(
  clock_cpg_cosmo_meqtls %>%
    select(snp.chr, snp.pos, A1, A2, cpg, beta.comb) %>%
    rename(chr = snp.chr, pos = snp.pos, REF = A1, ALT = A2, CpG = cpg, beta = beta.comb),
  clock_cpg_genoa_meqtls %>%
    select(chr, ps, allele1, allele0, CpG, beta) %>%
    rename(chr = chr, pos = ps, REF = allele0, ALT = allele1, CpG = CpG, beta = beta),
  clock_cpg_uk_meqtls %>%
    select(SNP.chr, SNP.pos, REF, ALT, CpGmarker, Beta) %>%
    rename(chr = SNP.chr, pos = SNP.pos, REF = REF, ALT = ALT, CpG = CpGmarker, beta = Beta)
)

# Write the combined data to a tab-delimited txt file
# write_tsv(clock_cpg_meqtls, "clock_cpg_meqtls.txt")

# Get the number of unique clock CpGs with meQTLs in each dataset
length(unique(clock_cpg_cosmo_meqtls$cpg))
length(unique(clock_cpg_genoa_meqtls$CpG))
length(unique(clock_cpg_uk_meqtls$CpGmarker))

# Get the number of unique SNPs in each dataset
length(unique(clock_cpg_cosmo_meqtls$snp))
length(unique(clock_cpg_genoa_meqtls$rs))
length(unique(clock_cpg_uk_meqtls$Top.SNP))

# Find overlap of clock CpGs with meQTLs in the 3 datasets
length(intersect(clock_cpg_cosmo_meqtls$cpg, clock_cpg_genoa_meqtls$CpG))
length(intersect(clock_cpg_cosmo_meqtls$cpg, clock_cpg_uk_meqtls$CpGmarker))
length(intersect(clock_cpg_genoa_meqtls$CpG, clock_cpg_uk_meqtls$CpGmarker))

# Get the number of unique pairs of clock CpGs and SNPs in each dataset
clock_cpg_cosmo_meqtls %>% summarise(n = n_distinct(cpg, snp))
clock_cpg_genoa_meqtls %>% summarise(n = n_distinct(CpG, rs))
clock_cpg_uk_meqtls %>% summarise(n = n_distinct(CpGmarker, Top.SNP))

# Create new IDs to tie in all three datasets. This works because all meQTLs here are in hg19 format. We'll use the format "chrN:posY
clock_cpg_cosmo_meqtls <- clock_cpg_cosmo_meqtls %>% mutate(new_id = paste0("chr", snp.chr, ":", snp.pos, "_", cpg))
clock_cpg_genoa_meqtls <- clock_cpg_genoa_meqtls %>% mutate(new_id = paste0("chr", chr, ":", ps, "_", CpG))
clock_cpg_uk_meqtls <- clock_cpg_uk_meqtls %>% mutate(new_id = paste0("chr", SNP.chr, ":", SNP.pos, "_", CpGmarker))

# Make an upset plot to show the overlaps between the datasets
# Combine the data frames and create a list of sets
list_of_sets <- list(
  Cosmo = clock_cpg_cosmo_meqtls$new_id,
  Genoa = clock_cpg_genoa_meqtls$new_id,
  UK = clock_cpg_uk_meqtls$new_id
)

# Convert the list to a data frame suitable for UpSetR
sets_df <- fromList(list_of_sets)

# Create a custom color palette
main_bar_color <- "#56B4E9"
sets_bar_color <- "#E69F00"
point_color <- "#0072B2"
line_color <- "#009E73"

# Create the UpSet plot with enhanced aesthetics
upset(
  sets_df,
  sets = c("Cosmo", "Genoa", "UK"),
  main.bar.color = main_bar_color,
  sets.bar.color = sets_bar_color,
  order.by = "freq",
  keep.order = TRUE,
  point.size = 3.5, # Size of points in the intersections
  line.size = 1.5, # Thickness of lines in the intersections
  mb.ratio = c(0.6, 0.4), # Ratio of main bar plot to set size plot
  text.scale = c(1.5, 1.5, 1, 1, 1.2, 1.2), # Scale text for different elements
  set_size.show = F, # Show set size numbers
)

######## Plot allele frequencies for the meQTLs affecting Horvath clock CpGs ########

# Plot the allele frequencies for the cosmo meQTLs affecting Horvath clock CpGs
plot1 <- ggplot(clock_cpg_cosmo_meqtls, aes(x = eaf)) +
  geom_histogram(fill = "red", color = "black") +
  xlim(0, 0.5) +
  labs(
    x = expression(paste("Allele Frequency")),
    y = "Number of SNPs"
  ) +
  theme_minimal()

# ggpubr version
plot1_ggpubr <- gghistogram(
  data = clock_cpg_cosmo_meqtls,
  x = "eaf",
  fill = "red",
  color = "black",
  xlim = c(0, 0.5),
  ylab = "Number of SNPs",
  title = NULL, # Add title if needed
) +
  xlab("Allele Frequency")

# Plot the allele frequencies for the GENOA meQTLs affecting Horvath clock CpGs
plot2 <- ggplot(clock_cpg_genoa_meqtls, aes(x = af)) +
  geom_histogram(fill = "green", color = "black") +
  xlim(0, 0.5) +
  labs(
    x = expression(paste("Allele Frequency")),
    y = "Number of SNPs"
  ) +
  theme_minimal()

# ggpubr version
plot2_ggpubr <- gghistogram(
  data = clock_cpg_genoa_meqtls,
  x = "af",
  fill = "green",
  color = "black",
  xlim = c(0, 0.5),
  ylab = "Number of SNPs",
  title = NULL, # Add title if needed
) +
  xlab("Allele Frequency")

# Plot the allele frequencies for the UK meQTLs affecting Horvath clock CpGs
plot3 <- ggplot(clock_cpg_uk_meqtls, aes(x = MAF)) +
  geom_histogram(fill = "blue", color = "black") +
  xlim(0, 0.5) +
  labs(
    x = expression(paste("Allele Frequency")),
    y = "Number of SNPs"
  ) +
  theme_minimal()

# ggpubr version
plot3_ggpubr <- gghistogram(
  data = clock_cpg_uk_meqtls,
  x = "MAF",
  fill = "blue",
  color = "black",
  xlim = c(0, 0.5),
  ylab = "Number of SNPs",
  title = NULL, # Add title if needed
) +
  xlab("Allele Frequency")

# Combine all the plots using cowplot, and add a legend
plot_grid(plot1, plot2, plot3, ncol = 3, labels = c("A", "B", "C"))
plot_grid(plot1_ggpubr, plot2_ggpubr, plot3_ggpubr, ncol = 3, labels = c("A", "B", "C"))

###### All clock meQTL analysis ######
# Get the DunedinPACE CpGs
dunedin_cpgs <- as.data.frame(getRequiredProbes())

# Find overlaps between DunedinPACE CpGs and the meQTL datasets
dunedin_cosmo_meqtls <- dunedin_cpgs %>% filter(DunedinPACE %in% clock_cpg_cosmo_meqtls$cpg)
dunedin_genoa_meqtls <- dunedin_cpgs %>% filter(DunedinPACE %in% clock_cpg_genoa_meqtls$CpG)
dunedin_uk_meqtls <- dunedin_cpgs %>% filter(DunedinPACE %in% clock_cpg_uk_meqtls$CpGmarker)

# Get the CpGs for all clocks
load_DNAm_Clocks_data()
caus_age_cpgs <- read_csv("YingCausAge.csv")
cortical_clock_cpgs <- read_table("CorticalClockCoefs.txt")


# Find overlaps for the Horvath, Hannum, BLUP, and PhenoAge clocks with the meQTL datasets
horvath_cosmo_meqtls <- coefHorvath %>% filter(CpGmarker %in% clock_cpg_cosmo_meqtls$cpg)
horvath_genoa_meqtls <- coefHorvath %>% filter(CpGmarker %in% clock_cpg_genoa_meqtls$CpG)
horvath_uk_meqtls <- coefHorvath %>% filter(CpGmarker %in% clock_cpg_uk_meqtls$CpGmarker)

hannum_cosmo_meqtls <- coefHannum %>% filter(CpGmarker %in% clock_cpg_cosmo_meqtls$cpg)
hannum_genoa_meqtls <- coefHannum %>% filter(CpGmarker %in% clock_cpg_genoa_meqtls$CpG)
hannum_uk_meqtls <- coefHannum %>% filter(CpGmarker %in% clock_cpg_uk_meqtls$CpGmarker)

blup_cosmo_meqtls <- coefBLUP %>% filter(CpGmarker %in% clock_cpg_cosmo_meqtls$cpg)
blup_genoa_meqtls <- coefBLUP %>% filter(CpGmarker %in% clock_cpg_genoa_meqtls$CpG)
blup_uk_meqtls <- coefBLUP %>% filter(CpGmarker %in% clock_cpg_uk_meqtls$CpGmarker)

en_cosmo_meqtls <- coefEN %>% filter(CpGmarker %in% clock_cpg_cosmo_meqtls$cpg)
en_genoa_meqtls <- coefEN %>% filter(CpGmarker %in% clock_cpg_genoa_meqtls$CpG)
en_uk_meqtls <- coefEN %>% filter(CpGmarker %in% clock_cpg_uk_meqtls$CpGmarker)

pheno_cosmo_meqtls <- coefLevine %>% filter(CpGmarker %in% clock_cpg_cosmo_meqtls$cpg)
pheno_genoa_meqtls <- coefLevine %>% filter(CpGmarker %in% clock_cpg_genoa_meqtls$CpG)
pheno_uk_meqtls <- coefLevine %>% filter(CpGmarker %in% clock_cpg_uk_meqtls$CpGmarker)

caus_age_cosmo_meqtls <- caus_age_cpgs %>% filter(CpGmarker %in% clock_cpg_cosmo_meqtls$cpg)
caus_age_genoa_meqtls <- caus_age_cpgs %>% filter(CpGmarker %in% clock_cpg_genoa_meqtls$CpG)
caus_age_uk_meqtls <- caus_age_cpgs %>% filter(CpGmarker %in% clock_cpg_uk_meqtls$CpGmarker)

cortical_clock_cosmo_meqtls <- cortical_clock_cpgs %>% filter(probe %in% clock_cpg_cosmo_meqtls$cpg)
cortical_clock_genoa_meqtls <- cortical_clock_cpgs %>% filter(probe %in% clock_cpg_genoa_meqtls$CpG)
cortical_clock_uk_meqtls <- cortical_clock_cpgs %>% filter(probe %in% clock_cpg_uk_meqtls$CpGmarker)


# Plot the number of meQTLs for each clock
plot4 <- ggplot(
  data.frame(
    Clock = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge"),
    Count = c(nrow(horvath_cosmo_meqtls), nrow(hannum_cosmo_meqtls), nrow(blup_cosmo_meqtls), nrow(en_cosmo_meqtls), nrow(pheno_cosmo_meqtls))
  ),
  aes(x = Clock, y = Count, fill = Clock)
) +
  geom_bar(stat = "identity") +
  labs(x = "Clock", y = "Number of meQTLs") +
  # Order the x-axis as Horvath, Hannum, BLUP, PhenoAge
  scale_x_discrete(limits = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge")) +
  theme_minimal() +
  theme(legend.position = "none")

plot5 <- ggplot(
  data.frame(
    Clock = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge"),
    Count = c(nrow(horvath_genoa_meqtls), nrow(hannum_genoa_meqtls), nrow(blup_genoa_meqtls), nrow(en_genoa_meqtls), nrow(pheno_genoa_meqtls))
  ),
  aes(x = Clock, y = Count, fill = Clock)
) +
  geom_bar(stat = "identity") +
  labs(x = "Clock", y = "Number of meQTLs") +
  # Order the x-axis as Horvath, Hannum, BLUP, PhenoAge
  scale_x_discrete(limits = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge")) +
  theme_minimal() +
  theme(legend.position = "none")

plot6 <- ggplot(
  data.frame(
    Clock = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge"),
    Count = c(nrow(horvath_uk_meqtls), nrow(hannum_uk_meqtls), nrow(blup_uk_meqtls), nrow(en_uk_meqtls), nrow(pheno_uk_meqtls))
  ),
  aes(x = Clock, y = Count, fill = Clock)
) +
  geom_bar(stat = "identity") +
  labs(x = "Clock", y = "Number of meQTLs") +
  # Order the x-axis as Horvath, Hannum, BLUP, PhenoAge
  scale_x_discrete(limits = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge")) +
  theme_minimal() +
  theme(legend.position = "none")

# Combine all the plots using cowplot, and add a legend
plot_grid(plot4, plot5, plot6, ncol = 3, labels = c("A", "B", "C"))

# Do the same as above, but with the proportion of clock CpGs with meQTLs
plot7 <- ggplot(
  data.frame(
    Clock = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge"),
    Proportion = c(
      nrow(horvath_cosmo_meqtls) / nrow(coefHorvath),
      nrow(hannum_cosmo_meqtls) / nrow(coefHannum),
      nrow(blup_cosmo_meqtls) / nrow(coefBLUP),
      nrow(en_cosmo_meqtls) / nrow(coefEN),
      nrow(pheno_cosmo_meqtls) / nrow(coefLevine)
    )
  ),
  aes(x = Clock, y = Proportion, fill = Clock)
) +
  geom_bar(stat = "identity") +
  labs(x = "Clock", y = "Proportion of CpGs with meQTLs") +
  # Order the x-axis as Horvath, Hannum, BLUP, PhenoAge
  scale_x_discrete(limits = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.75)) +
  theme_minimal() +
  theme(legend.position = "none")

plot8 <- ggplot(
  data.frame(
    Clock = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge"),
    Proportion = c(
      nrow(horvath_genoa_meqtls) / nrow(coefHorvath),
      nrow(hannum_genoa_meqtls) / nrow(coefHannum),
      nrow(blup_genoa_meqtls) / nrow(coefBLUP),
      nrow(en_genoa_meqtls) / nrow(coefEN),
      nrow(pheno_genoa_meqtls) / nrow(coefLevine)
    )
  ),
  aes(x = Clock, y = Proportion, fill = Clock)
) +
  geom_bar(stat = "identity") +
  labs(x = "Clock", y = "Proportion of CpGs with meQTLs") +
  # Order the x-axis as Horvath, Hannum, BLUP, PhenoAge
  scale_x_discrete(limits = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.75)) +
  theme_minimal() +
  theme(legend.position = "none")

plot9 <- ggplot(
  data.frame(
    Clock = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge"),
    Proportion = c(
      nrow(horvath_uk_meqtls) / nrow(coefHorvath),
      nrow(hannum_uk_meqtls) / nrow(coefHannum),
      nrow(blup_uk_meqtls) / nrow(coefBLUP),
      nrow(en_uk_meqtls) / nrow(coefEN),
      nrow(pheno_uk_meqtls) / nrow(coefLevine)
    )
  ),
  aes(x = Clock, y = Proportion, fill = Clock)
) +
  geom_bar(stat = "identity") +
  labs(x = "Clock", y = "Proportion of CpGs with meQTLs") +
  # Order the x-axis as Horvath, Hannum, BLUP, PhenoAge
  scale_x_discrete(limits = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.75)) +
  theme_minimal() +
  theme(legend.position = "none")

# Combine all the plots using cowplot, and add a legend
plot_grid(plot7, plot8, plot9, ncol = 3, labels = c("A", "B", "C"))

# For each clock, join the three sets of meQTLs. Combine based on the CpGmarker column. Do not use bind rows.
horvath_meqtls <- full_join(horvath_cosmo_meqtls, horvath_genoa_meqtls, by = "CpGmarker") %>%
  full_join(horvath_uk_meqtls, by = "CpGmarker")

hannum_meqtls <- full_join(hannum_cosmo_meqtls, hannum_genoa_meqtls, by = "CpGmarker") %>%
  full_join(hannum_uk_meqtls, by = "CpGmarker")

blup_meqtls <- full_join(blup_cosmo_meqtls, blup_genoa_meqtls, by = "CpGmarker") %>%
  full_join(blup_uk_meqtls, by = "CpGmarker")

en_meqtls <- full_join(en_cosmo_meqtls, en_genoa_meqtls, by = "CpGmarker") %>%
  full_join(en_uk_meqtls, by = "CpGmarker")

pheno_meqtls <- full_join(pheno_cosmo_meqtls, pheno_genoa_meqtls, by = "CpGmarker") %>%
  full_join(pheno_uk_meqtls, by = "CpGmarker")

caus_age_meqtls <- full_join(caus_age_cosmo_meqtls, caus_age_genoa_meqtls, by = "CpGmarker") %>%
  full_join(caus_age_uk_meqtls, by = "CpGmarker")

cortical_clock_meqtls <- full_join(cortical_clock_cosmo_meqtls, cortical_clock_genoa_meqtls, by = "probe") %>%
  full_join(cortical_clock_uk_meqtls, by = "probe")

# Now plot the number of meQTLs for each clock
plot10 <- ggplot(
  data.frame(
    Clock = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge", "DunedinPACE"),
    Count = c(nrow(horvath_meqtls), nrow(hannum_meqtls), nrow(blup_meqtls), nrow(en_meqtls), nrow(pheno_meqtls), 0)
  ),
  aes(x = Clock, y = Count, fill = Clock)
) +
  geom_bar(stat = "identity") +
  labs(x = "Clock", y = "Number of meQTLs") +
  # Order the x-axis as Horvath, Hannum, BLUP, PhenoAge
  scale_x_discrete(limits = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge", "DunedinPACE")) +
  theme_minimal() +
  theme(legend.position = "none")

# ggpubr version
plot10_ggpubr <- ggbarplot(
  data = data.frame(
    Clock = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge", "DunedinPACE"),
    Count = c(nrow(horvath_meqtls), nrow(hannum_meqtls), nrow(blup_meqtls), nrow(en_meqtls), nrow(pheno_meqtls), 0)
  ),
  x = "Clock",
  y = "Count",
  fill = "Clock",
  xlab = "Clock",
  ylab = "Number of meQTLs",
  order = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge", "DunedinPACE") # Orders the x-axis categories
) +
  theme(legend.position = "none")

# Now plot the proportion of CpGs with meQTLs for each clock
plot11 <- ggplot(
  data.frame(
    Clock = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge", "DunedinPACE"),
    Proportion = c(
      nrow(horvath_meqtls) / nrow(coefHorvath),
      nrow(hannum_meqtls) / nrow(coefHannum),
      nrow(blup_meqtls) / nrow(coefBLUP),
      nrow(en_meqtls) / nrow(coefEN),
      nrow(pheno_meqtls) / nrow(coefLevine),
      0
    )
  ),
  aes(x = Clock, y = Proportion, fill = Clock)
) +
  geom_bar(stat = "identity") +
  labs(x = "Clock", y = "Proportion of CpGs with meQTLs") +
  # Order the x-axis as Horvath, Hannum, BLUP, PhenoAge
  scale_x_discrete(limits = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge", "DunedinPACE")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.80)) +
  theme_minimal() +
  theme(legend.position = "none")

# ggpubr version
plot11_ggpubr <- ggbarplot(
  data = data.frame(
    Clock = c("Horvath", "Hannum", "EN", "PhenoAge", "DunedinPACE"),
    Proportion = c(
      nrow(horvath_meqtls) / nrow(coefHorvath),
      nrow(hannum_meqtls) / nrow(coefHannum),
      nrow(en_meqtls) / nrow(coefEN),
      nrow(pheno_meqtls) / nrow(coefLevine),
      0
    ),
    Count = c(
      nrow(horvath_meqtls),
      nrow(hannum_meqtls),
      nrow(en_meqtls),
      nrow(pheno_meqtls),
      0
    )
  ),
  x = "Clock",
  y = "Proportion",
  # fill = "Clock",
  xlab = "Clock",
  ylab = "Proportion of CpGs with meQTLs",
  order = c("Horvath", "Hannum", "BLUP", "EN", "PhenoAge", "DunedinPACE") # Orders the x-axis categories
) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.80)) +
  geom_text(
    aes(label = Count),
    vjust = -0.5, # Adjust vertical position of labels
    size = 5, # Size of the labels
    color = "black"
  ) +
  theme(legend.position = "none")

plot11_ggpubr <- plot11_ggpubr +
  theme(text = element_text(size = 16))
print(plot11_ggpubr)
saveRDS(plot11_ggpubr, "plots_RDS/fig5c.rds")

# Combine all the plots using cowplot, and add a legend
plot_grid(plot10, plot11, ncol = 2, labels = c("A", "B"))
plot_grid(plot10_ggpubr, plot11_ggpubr, ncol = 2, labels = c("A", "B"))

# Make a plot for Shea Andrews with the meQTLs for CausAge and the cortical methylation clock
plot12_ggpubr <- ggbarplot(
  data = data.frame(
    Clock = c("Horvath", "CausAge", "Cortical Clock"),
    Proportion = c(
      nrow(horvath_meqtls) / nrow(coefHorvath),
      nrow(caus_age_meqtls) / nrow(caus_age_cpgs),
      nrow(cortical_clock_meqtls) / nrow(cortical_clock_cpgs)
    ),
    Count = c(
      nrow(horvath_meqtls),
      nrow(caus_age_meqtls),
      nrow(cortical_clock_meqtls)
    )
  ),
  x = "Clock",
  y = "Proportion",
  fill = "Clock",
  xlab = "Clock",
  ylab = "Proportion of CpGs with meQTLs",
  order = c("Horvath", "CausAge", "Cortical Clock") # Orders the x-axis categories
) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.80)) +
  geom_text(
    aes(label = Count),
    vjust = -0.5, # Adjust vertical position of labels
    size = 5, # Size of the labels
    color = "black"
  ) +
  theme(legend.position = "none")

print(plot12_ggpubr)

plot13_ggpubr <- ggbarplot(
  data = data.frame(
    Clock = c("CausAge", "Cortical Clock"),
    Proportion = c(
      nrow(caus_age_meqtls) / nrow(caus_age_cpgs),
      nrow(cortical_clock_meqtls) / nrow(cortical_clock_cpgs)
    ),
    Count = c(
      nrow(caus_age_meqtls),
      nrow(cortical_clock_meqtls)
    )
  ),
  x = "Clock",
  y = "Proportion",
  fill = "Clock",
  xlab = "Clock",
  ylab = "Proportion of CpGs with meQTLs",
  order = c("CausAge", "Cortical Clock") # Orders the x-axis categories
) +
  scale_y_continuous(labels = scales::percent) +
  geom_text(
    aes(label = Count),
    vjust = -0.5, # Adjust vertical position of labels
    size = 5, # Size of the labels
    color = "black"
  ) +
  theme(legend.position = "none")

print(plot13_ggpubr)



# Create a combined dataframe with CpG sites and their clock sources
combined_meqtls <- rbind(
  data.frame(
    CpG = horvath_meqtls$CpGmarker,
    Clock = "Horvath"
  ),
  data.frame(
    CpG = hannum_meqtls$CpGmarker,
    Clock = "Hannum"
  ),
  data.frame(
    CpG = en_meqtls$CpGmarker,
    Clock = "EN"
  ),
  data.frame(
    CpG = pheno_meqtls$CpGmarker,
    Clock = "Pheno"
  )
)

# View the first few rows
head(combined_meqtls)

# Optional: Remove duplicate CpG-Clock combinations if any exist
combined_meqtls <- unique(combined_meqtls)

cpg_clock_counts <- aggregate(Clock ~ CpG, data = combined_meqtls, FUN = length)
names(cpg_clock_counts)[2] <- "n_clocks"

# Get CpGs that appear in multiple clocks (2 or more)
multi_clock_cpgs <- cpg_clock_counts[cpg_clock_counts$n_clocks > 1, ]
multi_clock_cpgs <- multi_clock_cpgs[order(-multi_clock_cpgs$n_clocks), ]

# Show summary
cat("CpGs appearing in multiple clocks:\n")
print(table(multi_clock_cpgs$n_clocks))

# Get detailed information about which clocks each multi-clock CpG appears in
multi_clock_details <- combined_meqtls[combined_meqtls$CpG %in% multi_clock_cpgs$CpG, ]
multi_clock_details <- multi_clock_details[order(multi_clock_details$CpG), ]

# Create a summary showing which clocks each CpG appears in
cpg_clock_list <- aggregate(Clock ~ CpG,
  data = multi_clock_details,
  FUN = function(x) paste(sort(x), collapse = ", ")
)
names(cpg_clock_list)[2] <- "Clocks"

# Merge with count information
cpg_summary <- merge(cpg_clock_list, multi_clock_cpgs, by = "CpG")
cpg_summary <- cpg_summary[order(-cpg_summary$n_clocks), ]

# Display results
cat("\nCpG sites affected by meQTL in multiple clocks:\n")
print(cpg_summary)


saveRDS(combined_meqtls, file = "~/Desktop/Capra Lab/Thesis Project/Aim_1/all_meqtl_affected_cpgs.rds")

# Count the number of unique variants per Horvath CpG
horvath_variants_per_cpg <- clock_cpg_meqtls %>%
  mutate(variant_id = paste(chr, pos, REF, ALT, sep = ":")) %>%
  group_by(CpG) %>%
  summarise(n_variants = n_distinct(variant_id)) %>%
  arrange(desc(n_variants))

# Calculate statistics for the lines
mean_variants <- mean(horvath_variants_per_cpg$n_variants)
median_variants <- median(horvath_variants_per_cpg$n_variants)
max_variants <- max(horvath_variants_per_cpg$n_variants)

# Visualize the number of unique variants per CpG
plot_horvath_variants <- ggbarplot(
  horvath_variants_per_cpg,
  x = "CpG",
  y = "n_variants",
  fill = "#56B4E9",
  color = "#56B4E9",
  sort.val = "desc", # Sort CpGs by number of variants
  title = "Unique Variants per Horvath Clock CpG",
  xlab = "CpG Site (Ordered by variant count)",
  ylab = "Number of Unique Variants"
) +
  geom_hline(yintercept = mean_variants, linetype = "dashed", color = "#D55E00", size = 1) +
  geom_hline(yintercept = median_variants, linetype = "dashed", color = "#009E73", size = 1) +
  geom_hline(yintercept = max_variants, linetype = "dashed", color = "#0072B2", size = 1) +
  annotate("text", x = 10, y = mean_variants, label = paste("Mean:", round(mean_variants, 1)), color = "#D55E00", vjust = -1, fontface = "bold") +
  annotate("text", x = 10, y = median_variants, label = paste("Median:", median_variants), color = "#009E73", vjust = 1.5, fontface = "bold") +
  annotate("text", x = 10, y = max_variants, label = paste("Max:", max_variants), color = "#0072B2", vjust = -1, fontface = "bold") +
  theme_pubclean() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_blank(), # Hidden due to large number of CpGs
    axis.ticks.x = element_blank()
  )

# Display the plot
print(plot_horvath_variants)
