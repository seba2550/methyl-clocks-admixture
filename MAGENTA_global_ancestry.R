# This script analyzes the global ancestry proportions for the individuals in the MAGENTA study, including a look at the DNAmAge predictions.

# Establish working directory
setwd("/Users/sgonzalez/Desktop/Capra Lab/Thesis Project/Aim_1")

# Load our libraries as needed
library(tidyverse)
library(readxl)
library(moderndive)
library(cowplot)
library(ggpubr)
library(pheatmap)

# Read in the results of the DNAmAge calculations along with the corresponding metadata
bio_age <- read.csv("bio_age_estimates_magenta_age_diff_metadata.csv", row.names = 1)

# Read in the global ancestry proportions
AA_global_anc <- read.table("genotyping_data/AA_ancestry_proportions.txt", header = T)
NHW_global_anc <- read.table("genotyping_data/NHW_ancestry_proportions.txt", header = T)
HISP_global_anc <- read.table("genotyping_data/HISPANIC_ancestry_proportions.txt", header = T)

#tmp <- read_xlsx("AANHW_472_id4genotyping.xlsx")

# Merge the three ancestry proportion tables together
global_anc <- rbind(AA_global_anc, NHW_global_anc, HISP_global_anc)

#sample_metadata <- read_xlsx("ADmethy_pheno.xlsx")
#methylation_meta <- subset(sample_metadata, sample_metadata$Beta_ID %in% colnames(normalized_combined))
#tmp <- subset(methylation_meta, !CGI %in% global_anc$sample_ID)

#write_csv(tmp, "MAGENTA_methylation_samples_missing_genetics.csv")

# Now add the global ancestry proportions to the biological age table, joining by CGI and sample_id
bio_age_global_anc <- bio_age %>%
  left_join(global_anc, by = c("CGI" = "sample_ID"))

########## Let's make some plots ########## 
ggplot(bio_age_global_anc, aes(x = ETHNICITY, y = CEU)) +
  geom_jitter(aes(color = COHORT)) +
  theme_minimal()
ggplot(bio_age_global_anc, aes(x = ETHNICITY, y = YRI)) +
  geom_jitter(aes(color = COHORT)) +
  theme_minimal()
ggplot(bio_age_global_anc, aes(x = ETHNICITY, y = PEL)) +
  geom_jitter(aes(color = COHORT)) +
  theme_minimal()

# Melt the data. I want to have the global ancestry proportions in a single column, and to plot those proportions for each cohort separately.
bio_age_global_anc_melted <- bio_age_global_anc %>%
  gather(key = "global_anc", value = "proportion", CEU, YRI, PEL)

# Make a barplot of the global ancestry proportions from the melted data
bio_age_global_anc_melted %>% filter(COHORT != "NHW") %>% ggplot(aes(x = ETHNICITY, y = proportion, fill = global_anc)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ COHORT) +
  theme_minimal()

bio_age_global_anc_melted %>% filter(COHORT != "NHW") %>% ggplot(aes(x = COHORT, y = proportion, fill = global_anc)) +
  geom_boxplot() +
  theme_minimal()

bio_age_global_anc_melted %>%
  filter(COHORT == "REAAADI" & STATUS == "CONTROL") %>%
  ggplot(aes(x = global_anc, y = proportion, fill = global_anc)) +
  geom_boxplot() +
  xlab("Ancestry") +
  ylab("Global Ancestry Proportion") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

bio_age_global_anc_melted %>%
  filter(COHORT == "PRADI" & STATUS == "CONTROL") %>%
  ggplot(aes(x = global_anc, y = proportion, fill = global_anc)) +
  geom_boxplot() +
  xlab("Ancestry") +
  ylab("Global Ancestry Proportion") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

bio_age_global_anc_melted %>%
  filter(COHORT == "CuADI" & STATUS == "CONTROL") %>%
  ggplot(aes(x = global_anc, y = proportion, fill = global_anc)) +
  geom_boxplot() +
  xlab("Ancestry") +
  ylab("Global Ancestry Proportion") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

bio_age_global_anc_melted %>%
  filter(COHORT == "PERUVIAN" & STATUS == "CONTROL") %>%
  ggplot(aes(x = global_anc, y = proportion, fill = global_anc)) +
  geom_boxplot() +
  xlab("Ancestry") +
  ylab("Global Ancestry Proportion") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

# Now I'll make plots showing the correlation between predicted ages and chronological ages. I'll do this for each cohort separately.
# I'll also add the correlation for the NHW to each plot for comparison.
bio_age_global_anc %>%
  filter(COHORT == "NHW" & STATUS == "CONTROL") %>%
  ggplot(aes(x = Horvath, y = age)) +
  xlab("Horvath DNAmAge") +
  ylab("Chronological Age") +
  xlim(60,90) +
  ylim(60,90) +
  geom_point() +
  geom_smooth(method = "lm", se = F, color = "salmon") +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  ggtitle("Non-Hispanic Whites")

# ggpubr version
fig2a <- bio_age_global_anc %>%
  filter(COHORT == "NHW" & STATUS == "CONTROL") %>%
  ggscatter(
    x = "Horvath",
    y = "age",
    #add = "reg.line",          # Adds a linear regression line
    conf.int = FALSE,          # Disables confidence interval shading
    add.params = list(color = "salmon"), # Color for the regression line
    xlab = "Horvath DNAmAge",
    ylab = "Chronological Age",
    xlim = c(60, 90),
    ylim = c(60, 90),
    title = "Non-Hispanic Whites"
  ) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = 62, y = 80,
    label = "r = 0.72",
    size = 6
  )

bio_age_global_anc %>%
  filter(COHORT == "REAAADI" & STATUS == "CONTROL") %>%
  ggplot(aes(x = Horvath, y = age)) +
  xlab("Horvath DNAmAge") +
  ylab("Chronological Age") +
  xlim(60,90) +
  ylim(60,90) +
  geom_point() +
  geom_smooth(method = "lm", se = F, color = "blue") +
  geom_smooth(data = filter(bio_age_global_anc, COHORT == "NHW" & STATUS == "CONTROL"), aes(x = Horvath, y = age), method = "lm", se = F, color = "grey", linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  ggtitle("African Americans") 

# ggpubr version
bio_age_global_anc %>%
  filter(COHORT == "REAAADI" & STATUS == "CONTROL") %>%
  ggscatter(
    x = "Horvath",
    y = "age",
    add = "reg.line",          # Adds a linear regression line
    conf.int = FALSE,          # Disables confidence interval shading
    add.params = list(color = "blue"), # Color for the regression line
    xlab = "Horvath DNAmAge",
    ylab = "Chronological Age",
    xlim = c(60, 90),
    ylim = c(60, 90),
    title = "African Americans"
  ) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

bio_age_global_anc %>%
  filter(COHORT == "PRADI" & STATUS == "CONTROL") %>%
  ggplot(aes(x = Horvath, y = age)) +
  xlab("Horvath DNAmAge") +
  ylab("Chronological Age") +
  xlim(60,90) +
  ylim(65,90) +
  geom_point() +
  geom_smooth(method = "lm", se = F, color = "green") +
  geom_smooth(data = filter(bio_age_global_anc, COHORT == "NHW" & STATUS == "CONTROL"), aes(x = Horvath, y = age), method = "lm", se = F, color = "grey", linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  ggtitle("Puerto Ricans")

# ggpubr version
bio_age_global_anc %>%
  filter(COHORT == "PRADI" & STATUS == "CONTROL") %>%
  ggscatter(
    x = "Horvath",
    y = "age",
    add = "reg.line",          # Adds a linear regression line
    conf.int = FALSE,          # Disables confidence interval shading
    add.params = list(color = "green"), # Color for the regression line
    xlab = "Horvath DNAmAge",
    ylab = "Chronological Age",
    xlim = c(60, 90),
    ylim = c(60, 90),
    title = "Puerto Ricans"
  ) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

bio_age_global_anc %>%
  filter(COHORT == "CuADI" & STATUS == "CONTROL") %>%
  ggplot(aes(x = Horvath, y = age)) +
  xlab("Horvath DNAmAge") +
  ylab("Chronological Age") +
  xlim(60,90) +
  ylim(60,90) +
  geom_point() +
  geom_smooth(method = "lm", se = F, color = "red") +
  geom_smooth(data = filter(bio_age_global_anc, COHORT == "NHW" & STATUS == "CONTROL"), aes(x = Horvath, y = age), method = "lm", se = F, color = "grey", linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  ggtitle("Cubans")

# ggpubr version
bio_age_global_anc %>%
  filter(COHORT == "CuADI" & STATUS == "CONTROL") %>%
  ggscatter(
    x = "Horvath",
    y = "age",
    add = "reg.line",          # Adds a linear regression line
    conf.int = FALSE,          # Disables confidence interval shading
    add.params = list(color = "red"), # Color for the regression line
    xlab = "Horvath DNAmAge",
    ylab = "Chronological Age",
    xlim = c(60, 90),
    ylim = c(60, 90),
    title = "Cubans"
  ) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

bio_age_global_anc %>%
  filter(COHORT == "PERUVIAN" & STATUS == "CONTROL") %>%
  ggplot(aes(x = Horvath, y = age)) +
  xlab("Horvath DNAmAge") +
  ylab("Chronological Age") +
  geom_point() +
  xlim(60,90) +
  ylim(60,90) +
  geom_smooth(method = "lm", se = F, color = "orange") +
  geom_smooth(data = filter(bio_age_global_anc, COHORT == "NHW" & STATUS == "CONTROL"), aes(x = Horvath, y = age), method = "lm", se = F, color = "grey", linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  ggtitle("Peruvians")

# ggpubr version
bio_age_global_anc %>%
  filter(COHORT == "PERUVIAN" & STATUS == "CONTROL") %>%
  ggscatter(
    x = "Horvath",
    y = "age",
    add = "reg.line",          # Adds a linear regression line
    conf.int = FALSE,          # Disables confidence interval shading
    add.params = list(color = "orange"), # Color for the regression line
    xlab = "Horvath DNAmAge",
    ylab = "Chronological Age",
    xlim = c(60, 90),
    ylim = c(60, 90),
    title = "Peruvians"
  ) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

# Put each correlation plot next to the global ancestry plot for each cohort using cowplot
# African Americans
p1 <- bio_age_global_anc_melted %>%
  filter(COHORT == "REAAADI" & STATUS == "CONTROL") %>%
  ggplot(aes(x = global_anc, y = proportion, fill = global_anc)) +
  geom_boxplot() +
  xlab("Ancestry") +
  ylab("Global Ancestry") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

# ggpurb version
p1_ggpubr <- bio_age_global_anc_melted %>%
  filter(COHORT == "REAAADI" & STATUS == "CONTROL") %>%
  ggboxplot(
    x = "global_anc",
    y = "proportion",
    fill = "global_anc",
    xlab = "Ancestry",
    ylab = "Global Ancestry",
    add = "none"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

p2 <- bio_age_global_anc %>%
  filter(COHORT == "REAAADI" & STATUS == "CONTROL") %>%
  ggplot(aes(x = Horvath, y = age)) +
  xlab("Horvath DNAmAge") +
  ylab("Chronological Age") +
  xlim(60,90) +
  ylim(60,90) +
  geom_point() +
  geom_smooth(method = "lm", se = F, color = "blue") +
  geom_smooth(data = filter(bio_age_global_anc, COHORT == "NHW" & STATUS == "CONTROL"), aes(x = Horvath, y = age), method = "lm", se = F, color = "salmon", linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  ggtitle("AA")

# ggpubr version
p2_ggpubr <- bio_age_global_anc %>%
  filter(COHORT == "REAAADI" & STATUS == "CONTROL") %>%
  ggscatter(
    x = "Horvath",
    y = "age",
    #add = "reg.line",          # Adds a linear regression line
    conf.int = FALSE,          # Disables confidence interval shading
    add.params = list(color = "blue"), # Color for the regression line
    xlab = "Horvath DNAmAge",
    ylab = "Chronological Age",
    xlim = c(60, 90),
    ylim = c(60, 90),
    title = "AA"
  ) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = 65, y = 80,
    label = "r = 0.51",
    size = 6
  )

# Puerto Ricans
p3 <- bio_age_global_anc_melted %>%
  filter(COHORT == "PRADI" & STATUS == "CONTROL") %>%
  ggplot(aes(x = global_anc, y = proportion, fill = global_anc)) +
  geom_boxplot() +
  xlab("Ancestry") +
  ylab("Global Ancestry") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

# ggpubr version
p3_ggpubr <- bio_age_global_anc_melted %>%
  filter(COHORT == "PRADI" & STATUS == "CONTROL") %>%
  ggboxplot(
    x = "global_anc",
    y = "proportion",
    fill = "global_anc",
    xlab = "Ancestry",
    ylab = "Global Ancestry",
    add = "none"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

p4 <- bio_age_global_anc %>%
  filter(COHORT == "PRADI" & STATUS == "CONTROL") %>%
  ggplot(aes(x = Horvath, y = age)) +
  xlab("Horvath DNAmAge") +
  ylab("Chronological Age") +
  xlim(60,90) +
  ylim(65,90) +
  geom_point() +
  geom_smooth(method = "lm", se = F, color = "green") +
  geom_smooth(data = filter(bio_age_global_anc, COHORT == "NHW" & STATUS == "CONTROL"), aes(x = Horvath, y = age), method = "lm", se = F, color = "salmon", linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  ggtitle("Puerto Ricans")

# ggpubr version
p4_ggpubr <- bio_age_global_anc %>%
  filter(COHORT == "PRADI" & STATUS == "CONTROL") %>%
  ggscatter(
    x = "Horvath",
    y = "age",
    #add = "reg.line",          # Adds a linear regression line
    conf.int = FALSE,          # Disables confidence interval shading
    add.params = list(color = "green"), # Color for the regression line
    xlab = "Horvath DNAmAge",
    ylab = "Chronological Age",
    xlim = c(60, 90),
    ylim = c(60, 90),
    title = "PUR"
  ) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = 65, y = 80,
    label = "r = 0.45",
    size = 6
  )

# Cubans
p5 <- bio_age_global_anc_melted %>%
  filter(COHORT == "CuADI" & STATUS == "CONTROL") %>%
  ggplot(aes(x = global_anc, y = proportion, fill = global_anc)) +
  geom_boxplot() +
  xlab("Ancestry") +
  ylab("Global Ancestry") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

# ggpubr version
p5_ggpubr <- bio_age_global_anc_melted %>%
  filter(COHORT == "CuADI" & STATUS == "CONTROL") %>%
  ggboxplot(
    x = "global_anc",
    y = "proportion",
    fill = "global_anc",
    xlab = "Ancestry",
    ylab = "Global Ancestry",
    add = "none"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

p6 <- bio_age_global_anc %>%
  filter(COHORT == "CuADI" & STATUS == "CONTROL") %>%
  ggplot(aes(x = Horvath, y = age)) +
  xlab("Horvath DNAmAge") +
  ylab("Chronological Age") +
  xlim(65,90) +
  ylim(60,90) +
  geom_point() +
  geom_smooth(method = "lm", se = F, color = "red") +
  geom_smooth(data = filter(bio_age_global_anc, COHORT == "NHW" & STATUS == "CONTROL"), aes(x = Horvath, y = age), method = "lm", se = F, color = "salmon", linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  ggtitle("Cubans")

# ggpubr version
p6_ggpubr <- bio_age_global_anc %>%
  filter(COHORT == "CuADI" & STATUS == "CONTROL") %>%
  ggscatter(
    x = "Horvath",
    y = "age",
    #add = "reg.line",          # Adds a linear regression line
    conf.int = FALSE,          # Disables confidence interval shading
    add.params = list(color = "red"), # Color for the regression line
    xlab = "Horvath DNAmAge",
    ylab = "Chronological Age",
    xlim = c(60, 90),
    ylim = c(60, 90),
    title = "CUB"
  ) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = 65, y = 80,
    label = "r = 0.68",
    size = 6
  )
# Peruvians
p7 <- bio_age_global_anc_melted %>%
  filter(COHORT == "PERUVIAN" & STATUS == "CONTROL") %>%
  ggplot(aes(x = global_anc, y = proportion, fill = global_anc)) +
  geom_boxplot() +
  xlab("Ancestry") +
  ylab("Global Ancestry") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

# ggpurb version
p7_ggpubr <- bio_age_global_anc_melted %>%
  filter(COHORT == "PERUVIAN" & STATUS == "CONTROL") %>%
  ggboxplot(
    x = "global_anc",
    y = "proportion",
    fill = "global_anc",
    xlab = "Ancestry",
    ylab = "Global Ancestry",
    add = "none"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

p8 <- bio_age_global_anc %>%
  filter(COHORT == "PERUVIAN" & STATUS == "CONTROL") %>%
  ggplot(aes(x = Horvath, y = age)) +
  xlab("Horvath DNAmAge") +
  ylab("Chronological Age") +
  geom_point() +
  xlim(60,80) +
  ylim(60,80) +
  geom_smooth(method = "lm", se = F, color = "orange") +
  geom_smooth(data = filter(bio_age_global_anc, COHORT == "NHW" & STATUS == "CONTROL"), aes(x = Horvath, y = age), method = "lm", se = F, color = "salmon", linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  ggtitle("Peruvians")

# ggpubr version
p8_ggpubr <- bio_age_global_anc %>%
  filter(COHORT == "PERUVIAN" & STATUS == "CONTROL") %>%
  ggscatter(
    x = "Horvath",
    y = "age",
    #add = "reg.line",          # Adds a linear regression line
    conf.int = FALSE,          # Disables confidence interval shading
    add.params = list(color = "orange"), # Color for the regression line
    xlab = "Horvath DNAmAge",
    ylab = "Chronological Age",
    xlim = c(60, 90),
    ylim = c(60, 90),
    title = "PER"
  ) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linetype = "dashed") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  annotate(
    "text",
    x = 65, y = 80,
    label = "r = 0.72",
    size = 6
  )

# Combine the plots
plot_grid(p2, p1, p4, p3, p6, p5, p8, p7, ncol = 4)
fig2b <- plot_grid(p2_ggpubr, p1_ggpubr, p4_ggpubr, p3_ggpubr, p6_ggpubr, p5_ggpubr, p8_ggpubr, p7_ggpubr, ncol = 4)

# Make plots that show the R correlation values for the DNAmAge predictions for each population group vs the mean AFR ancestry for the group.
# Get the correlation values and mean AFR ancestry for each population group
correlation_df <- bio_age_global_anc %>%
  group_by(COHORT) %>%
  filter(STATUS == "CONTROL") %>%
  mutate(COHORT = case_when(
    COHORT == "REAAADI" ~ "AA",
    COHORT == "PRADI" ~ "PUR",
    COHORT == "CuADI" ~ "Cuban",
    COHORT == "PERUVIAN" ~ "Peruvian",
    TRUE ~ COHORT
  )) %>%
  summarize(Horvath = cor(Horvath, age),
            Hannum = cor(Hannum, age),
            PhenoAge = cor(Levine, age),
            Zhang2019_EN = cor(EN, age),
            Zhang2019_BLUP = cor(BLUP, age),
            mean_afr_ancestry = mean(YRI, na.rm = TRUE))

spearman_correlation_df <- bio_age_global_anc %>%
  group_by(COHORT) %>%
  filter(STATUS == "CONTROL") %>%
  mutate(COHORT = case_when(
    COHORT == "REAAADI" ~ "AA",
    COHORT == "PRADI" ~ "PUR",
    COHORT == "CuADI" ~ "Cuban",
    COHORT == "PERUVIAN" ~ "Peruvian",
    TRUE ~ COHORT
  )) %>%
  summarize(Horvath = cor(Horvath, age, method = "spearman"),
            Hannum = cor(Hannum, age, method = "spearman"),
            PhenoAge = cor(Levine, age, method = "spearman"),
            Zhang2019_EN = cor(EN, age, method = "spearman"),
            Zhang2019_BLUP = cor(BLUP, age, method = "spearman"),
            mean_afr_ancestry = mean(YRI, na.rm = TRUE))



# Plot them together
# correlation_df %>%
#   ggplot(aes(x = mean_afr_ancestry, y = Horvath, color = COHORT)) +
#   geom_point(size = 3) +
#   geom_text(aes(label = COHORT), hjust = 0.5, vjust = 1.5, size = 8) +
#   xlab("Mean African Ancestry") +
#   ylab("Correlation between Horvath DNAmAge and Chronological Age") +
#   theme_minimal() +
#   theme(
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 16),
#     plot.title = element_text(size = 18),
#     legend.position = "none"
#   )
# 
# # Now we'll plot the correlation values for each cohort, splitting by DNA methylation clock
# correlation_df %>%
#   gather(key = "DNAmAge", value = "correlation", Horvath, Hannum, PhenoAge) %>%
#   ggplot(aes(x = mean_afr_ancestry, y = correlation, color = COHORT)) +
#   geom_point(size = 3) +
#   geom_text(aes(label = COHORT), hjust = 0.5, vjust = 1.5, size = 8) +
#   facet_wrap(~ DNAmAge) +
#   xlab("Mean African Ancestry") +
#   ylab("Correlation between DNAmAge and Chronological Age") +
#   theme_minimal() +
#   theme(
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 16),
#     plot.title = element_text(size = 18),
#     legend.position = "none"
#   )

# Plot the correlation values for each cohort, splitting by clock
correlation_df %>%
  gather(key = "DNAmAge", value = "correlation", Horvath, Hannum, PhenoAge, Zhang2019_EN) %>%
  mutate(COHORT = factor(COHORT, levels = c("NHW", "Cuban", "Peruvian", "PUR", "AA"))) %>%
  ggplot(aes(x = COHORT, y = correlation, fill = COHORT)) +
  geom_bar(stat = "identity", position = "dodge") +
  #geom_errorbar(aes(ymin = correlation - sd(correlation), ymax = correlation + sd(correlation)), width = 0.2, position = position_dodge(width = 0.9)) +
  facet_wrap(~ DNAmAge, scales = "free_x") +
  xlab("Cohort") +
  ylab("Correlation") +
  coord_cartesian(ylim = c(0.4, NA)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.position = "none",
    strip.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # Flip x ticks
  )

spearman_correlation_df %>%
  gather(key = "DNAmAge", value = "correlation", Horvath, Hannum, PhenoAge, Zhang2019_EN) %>%
  mutate(COHORT = factor(COHORT, levels = c("NHW", "Cuban", "Peruvian", "PUR", "AA"))) %>%
  ggplot(aes(x = COHORT, y = correlation, fill = COHORT)) +
  geom_bar(stat = "identity", position = "dodge") +
  #geom_errorbar(aes(ymin = correlation - sd(correlation), ymax = correlation + sd(correlation)), width = 0.2, position = position_dodge(width = 0.9)) +
  facet_wrap(~ DNAmAge, scales = "free_x") +
  xlab("Cohort") +
  ylab("Spearman Correlation") +
  coord_cartesian(ylim = c(0.4, NA)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.position = "none",
    strip.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # Flip x ticks
  )

# Plot the correlation values for each cohort, relative to the NHW cohort, splitting by clock.
correlation_df %>%
  gather(key = "DNAmAge", value = "correlation", Horvath, Hannum, PhenoAge, Zhang2019_EN, Zhang_BLUP) %>%
  mutate(COHORT = factor(COHORT, levels = c("NHW", "Cuban", "Peruvian", "PUR", "AA"))) %>%
  group_by(DNAmAge) %>%
  mutate(NHW_correlation = correlation[COHORT == "NHW"]) %>%
  ungroup() %>%
  mutate(correlation = correlation / NHW_correlation) %>%
  ggplot(aes(x = COHORT, y = correlation, fill = COHORT)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~ DNAmAge, scales = "free_x") +
  xlab("Cohort") +
  ylab("Correlation (Relative to NHW)") +
  coord_cartesian(ylim = c(0.4, NA)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.position = "none",
    strip.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # Flip x ticks
  )

spearman_correlation_df %>%
  gather(key = "DNAmAge", value = "correlation", Horvath, Hannum, PhenoAge, Zhang2019_EN, Zhang_BLUP) %>%
  mutate(COHORT = factor(COHORT, levels = c("NHW", "Cuban", "Peruvian", "PUR", "AA"))) %>%
  group_by(DNAmAge) %>%
  mutate(NHW_correlation = correlation[COHORT == "NHW"]) %>%
  ungroup() %>%
  mutate(correlation = correlation / NHW_correlation) %>%
  ggplot(aes(x = COHORT, y = correlation, fill = COHORT)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~ DNAmAge, scales = "free_x") +
  xlab("Cohort") +
  ylab("Spearman Correlation (Relative to NHW)") +
  coord_cartesian(ylim = c(0.4, NA)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.position = "none",
    strip.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # Flip x ticks
  )


# Gather the data and calculate the NHW correlation
df <- correlation_df %>%
  gather(key = "DNAmAge", value = "correlation", Horvath, Hannum, PhenoAge, Zhang2019_EN, Zhang2019_BLUP) %>%
  mutate(COHORT = factor(COHORT, levels = c("NHW", "Cuban", "Peruvian", "PUR", "AA"))) %>%
  group_by(DNAmAge) %>%
  mutate(NHW_correlation = correlation[COHORT == "NHW"]) %>%
  ungroup()

df_2 <- spearman_correlation_df %>%
  gather(key = "DNAmAge", value = "correlation", Horvath, Hannum, PhenoAge, Zhang2019_EN) %>%
  mutate(COHORT = factor(COHORT, levels = c("NHW", "Cuban", "Peruvian", "PUR", "AA"))) %>%
  group_by(DNAmAge) %>%
  mutate(NHW_correlation = correlation[COHORT == "NHW"]) %>%
  ungroup()

# Create a new column for relative correlation
df <- df %>%
  mutate(relative_correlation = correlation - NHW_correlation)

df_2 <- df_2 %>%
  mutate(relative_correlation = correlation - NHW_correlation)

# Plot the data
ggplot(df, aes(x = COHORT, y = relative_correlation, fill = COHORT)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = ifelse(COHORT == "NHW", round(correlation, 2), "")), 
            vjust = -0.5, size = 5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ DNAmAge, scales = "free_x") +
  xlab("Cohort") +
  ylab("Pearson Correlation (Relative to NHW)") +
  coord_cartesian(ylim = c(-0.2, 0.2)) +  # Adjust the y-axis limits as needed
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.position = "none",
    strip.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)  # Flip x ticks
  )

# ggpuubr version
# We want to add the conf. intervals for the differences, so we have to do some bootstrapping and correlation comparison first using cocor.
library(cocor)
library(boot)
# Bootstrap function
foo <- function(data, indices) {
  dt <- data[indices,]
  c(
    cor(dt[, 2], dt[, 28], method = 'p'), # Horvath
    cor(dt[, 5], dt[, 28], method = 'p'), # Hannum
    cor(dt[, 8], dt[, 28], method = 'p'), # PhenoAge
    cor(dt[, 22], dt[, 28], method = 'p'), # Zhang2019_BLUP
    cor(dt[, 25], dt[, 28], method = 'p') # Zhang2019_EN
  )
}

# Create some filtered dataframes. We want to keep only the control individuals for each population group.
nhw_controls <- bio_age_global_anc %>% filter(COHORT == "NHW" & STATUS == "CONTROL")
aa_controls <- bio_age_global_anc %>% filter(COHORT == "REAAADI" & STATUS == "CONTROL")
pr_controls <- bio_age_global_anc %>% filter(COHORT == "PRADI" & STATUS == "CONTROL")
cub_controls <- bio_age_global_anc %>% filter(COHORT == "CuADI" & STATUS == "CONTROL")
per_controls <- bio_age_global_anc %>% filter(COHORT == "PERUVIAN" & STATUS == "CONTROL")

# Now leverage the "boot" function to bootstrap  these correlations 1,000 times.
set.seed(42)
nhw_cor_bootstrap <- boot(nhw_controls, foo, R = 1000)
aa_cor_bootstrap <- boot(aa_controls, foo, R = 1000)
pr_cor_bootstrap <- boot(pr_controls, foo, R = 1000)
cub_cor_bootstrap <- boot(cub_controls, foo, R = 1000)
per_cor_bootstrap <- boot(per_controls, foo, R = 1000)

compare_correlations <- function(r1_boot, r2_boot, n1, n2) {
  results <- lapply(1:length(r1_boot$t0), function(i) {
    cocor.indep.groups(
      r1.jk = r1_boot$t0[i],
      r2.hm = r2_boot$t0[i],
      n1 = n1,
      n2 = n2,
      alternative = "greater",
      alpha = 0.05,
      conf.level = 0.95,
      null.value = 0
    )
  })
  results
}

# Perform the comparisons for each group
nhw_aa_compare <- compare_correlations(aa_cor_bootstrap, nhw_cor_bootstrap, nrow(aa_controls), nrow(nhw_controls))
nhw_pr_compare <- compare_correlations(pr_cor_bootstrap, nhw_cor_bootstrap, nrow(pr_controls), nrow(nhw_controls))
nhw_cub_compare <- compare_correlations(cub_cor_bootstrap, nhw_cor_bootstrap, nrow(cub_controls), nrow(nhw_controls))
nhw_per_compare <- compare_correlations(per_cor_bootstrap, nhw_cor_bootstrap, nrow(per_controls), nrow(nhw_controls))


# Add the confidence intervals to the relative correlation df
# Add a column for the lower and upper confidence intervals, per cohort and clock
df <- df %>% mutate(lower_ci = case_when(
  COHORT == "AA" & DNAmAge == "Horvath" ~ nhw_aa_compare[[1]]@zou2007$conf.int[1],
  COHORT == "AA" & DNAmAge == "Hannum" ~ nhw_aa_compare[[2]]@zou2007$conf.int[1],
  COHORT == "AA" & DNAmAge == "PhenoAge" ~ nhw_aa_compare[[3]]@zou2007$conf.int[1],
  COHORT == "AA" & DNAmAge == "Zhang2019_EN" ~ nhw_aa_compare[[5]]@zou2007$conf.int[1],
  COHORT == "PUR" & DNAmAge == "Horvath" ~ nhw_pr_compare[[1]]@zou2007$conf.int[1],
  COHORT == "PUR" & DNAmAge == "Hannum" ~ nhw_pr_compare[[2]]@zou2007$conf.int[1],
  COHORT == "PUR" & DNAmAge == "PhenoAge" ~ nhw_pr_compare[[3]]@zou2007$conf.int[1],
  COHORT == "PUR" & DNAmAge == "Zhang2019_EN" ~ nhw_pr_compare[[5]]@zou2007$conf.int[1],
  COHORT == "Cuban" & DNAmAge == "Horvath" ~ nhw_cub_compare[[1]]@zou2007$conf.int[1],
  COHORT == "Cuban" & DNAmAge == "Hannum" ~ nhw_cub_compare[[2]]@zou2007$conf.int[1],
  COHORT == "Cuban" & DNAmAge == "PhenoAge" ~ nhw_cub_compare[[3]]@zou2007$conf.int[1],
  COHORT == "Cuban" & DNAmAge == "Zhang2019_EN" ~ nhw_cub_compare[[5]]@zou2007$conf.int[1],
  COHORT == "Peruvian" & DNAmAge == "Horvath" ~ nhw_per_compare[[1]]@zou2007$conf.int[1],
  COHORT == "Peruvian" & DNAmAge == "Hannum" ~ nhw_per_compare[[2]]@zou2007$conf.int[1],
  COHORT == "Peruvian" & DNAmAge == "PhenoAge" ~ nhw_per_compare[[3]]@zou2007$conf.int[1],
  COHORT == "Peruvian" & DNAmAge == "Zhang2019_EN" ~ nhw_per_compare[[5]]@zou2007$conf.int[1]),
  upper_ci = case_when(
    COHORT == "AA" & DNAmAge == "Horvath" ~ nhw_aa_compare[[1]]@zou2007$conf.int[2],
    COHORT == "AA" & DNAmAge == "Hannum" ~ nhw_aa_compare[[2]]@zou2007$conf.int[2],
    COHORT == "AA" & DNAmAge == "PhenoAge" ~ nhw_aa_compare[[3]]@zou2007$conf.int[2],
    COHORT == "AA" & DNAmAge == "Zhang2019_EN" ~ nhw_aa_compare[[5]]@zou2007$conf.int[2],
    COHORT == "PUR" & DNAmAge == "Horvath" ~ nhw_pr_compare[[1]]@zou2007$conf.int[2],
    COHORT == "PUR" & DNAmAge == "Hannum" ~ nhw_pr_compare[[2]]@zou2007$conf.int[2],
    COHORT == "PUR" & DNAmAge == "PhenoAge" ~ nhw_pr_compare[[3]]@zou2007$conf.int[2],
    COHORT == "PUR" & DNAmAge == "Zhang2019_EN" ~ nhw_pr_compare[[5]]@zou2007$conf.int[2],
    COHORT == "Cuban" & DNAmAge == "Horvath" ~ nhw_cub_compare[[1]]@zou2007$conf.int[2],
    COHORT == "Cuban" & DNAmAge == "Hannum" ~ nhw_cub_compare[[2]]@zou2007$conf.int[2],
    COHORT == "Cuban" & DNAmAge == "PhenoAge" ~ nhw_cub_compare[[3]]@zou2007$conf.int[2],
    COHORT == "Cuban" & DNAmAge == "Zhang2019_EN" ~ nhw_cub_compare[[5]]@zou2007$conf.int[2],
    COHORT == "Peruvian" & DNAmAge == "Horvath" ~ nhw_per_compare[[1]]@zou2007$conf.int[2],
    COHORT == "Peruvian" & DNAmAge == "Hannum" ~ nhw_per_compare[[2]]@zou2007$conf.int[2],
    COHORT == "Peruvian" & DNAmAge == "PhenoAge" ~ nhw_per_compare[[3]]@zou2007$conf.int[2],
    COHORT == "Peruvian" & DNAmAge == "Zhang2019_EN" ~ nhw_per_compare[[5]]@zou2007$conf.int[2])
)

# Remove the Zhang2019_BLUP data
df <- df %>% filter(DNAmAge != "Zhang2019_BLUP")

# Reorder the clocks
df <- df %>% mutate(DNAmAge = factor(DNAmAge, levels = c("Horvath", "Hannum", "Zhang2019_EN", "PhenoAge")))

fig2c <- ggbarplot(
  df,
  x = "COHORT",
  y = "relative_correlation",
  fill = "COHORT",
  position = position_dodge(), # Same as position = "dodge"
  ylab = "Pearson Correlation (Relative to NHW)",
  xlab = "Cohort",
  ylim = c(-0.5, 0.5), # Set y-axis limits
  lab.pos = "out", # Ensures labels are placed outside the bars
  lab.size = 5, # Label size
  lab.hjust = 0.5, # Center text horizontally
  lab.vjust = -0.5, # Adjust vertical position
  title = NULL # Add title if needed
) +
  # Add error bars using the CI values
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci),
    width = 0.2, # Adjust width as needed
    position = position_dodge(width = 0.9)
  ) +
  geom_text(
    aes(label = ifelse(COHORT == "NHW", round(NHW_correlation, 2), "")),
    vjust = -0.5, size = 5, color = "black",
    position = position_dodge(width = 0.9)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ DNAmAge, scales = "free_x") +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.position = "none",
    strip.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


plot_grid(fig2a, fig2b, fig2c, labels = c("A", "B", "C"), rel_widths = c(1,2))

plot_grid(
  plot_grid(fig2a, fig2b, ncol = 1, rel_heights = c(1, 2), labels = "AUTO"),
  fig2c,
  labels = c("A", "C"),
  rel_widths = c(2, 1)
)

########## Ensemble approach
cor(((bio_age_global_anc$Horvath + bio_age_global_anc$Hannum + bio_age_global_anc$Levine + bio_age_global_anc$EN + bio_age_global_anc$BLUP)/3), bio_age_global_anc$age)
cor(bio_age_global_anc$Horvath, bio_age_global_anc$age)
cor(bio_age_global_anc$Hannum, bio_age_global_anc$age)
cor(bio_age_global_anc$Levine, bio_age_global_anc$age)
cor(bio_age_global_anc$EN, bio_age_global_anc$age)
cor(bio_age_global_anc$BLUP, bio_age_global_anc$age)

# Now we'll plot the correlation between the ensemble DNAmAge and chronological age for all MAGENTA controls
bio_age_global_anc %>%
  filter(STATUS == "CONTROL") %>%
  ggplot(aes(x = ((Horvath + Hannum + Levine + EN + BLUP)/3), y = age)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  xlab("Ensemble DNAmAge") +
  ylab("Chronological Age") +
  ggtitle("Ensemble DNAmAge vs Chronological Age") +
  annotate("text", x = 90, y = 100, label = paste("R = ",
                                                   round(cor(((bio_age_global_anc$Horvath + bio_age_global_anc$Hannum + bio_age_global_anc$Levine + bio_age_global_anc$EN + bio_age_global_anc$BLUP)/3), bio_age_global_anc$age), 2)),
           size = 8)

# Do the same as above, but split by cohort
bio_age_global_anc %>%
  filter(STATUS == "CONTROL") %>%
  ggplot(aes(x = ((Horvath + Hannum + Levine + EN + BLUP)/3), y = age)) + 
  geom_point(aes(color = COHORT), show.legend = F) +
  geom_smooth(method = "lm", se = F) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  ) +
  xlab("Ensemble DNAmAge") +
  ylab("Chronological Age") +
  ggtitle("Ensemble DNAmAge vs Chronological Age") +
  facet_wrap(~factor(COHORT, c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")), nrow = 2, labeller = cohort_labeler)

# Get the correlation values for the ensemble DNAmAge predictions for each population group
ensemble_correlations_df <- bio_age_global_anc %>%
  group_by(COHORT) %>%
  filter(STATUS == "CONTROL") %>%
  summarize(ensemble = cor(((Horvath + Hannum + Levine + EN + BLUP)/3), age))

# Get ensemble intrinsic age accelerations for each individual, and then make a boxplot to compare (for each cohort) the age acceleration values between AD and control individuals.
cohort_labels_human_readable <- list(
  #'NHW' = 'Non-Hispanic Whites',
  'REAAADI' = 'African Americans',
  'PRADI' = 'Puerto Ricans',
  'CuADI' = 'Cubans',
  'PERUVIAN' = 'Peruvians'
)
cohort_labeler <- function(variable,value){
  return(cohort_labels_human_readable[value])
}

bio_age_global_anc <- bio_age_global_anc %>%
  group_by(COHORT) %>%
  mutate(ensemble_intrinsic_age_accel = (ageAcc2.Horvath + ageAcc2.Hannum + ageAcc2.Levine + ageAcc2.EN + ageAcc2.BLUP) / 5)

ggplot(data = bio_age_global_anc, aes(x = STATUS, y = ensemble_intrinsic_age_accel, fill = STATUS)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylab("Ensemble Intrinsic Age Acceleration") +
  theme_minimal() +
  facet_wrap(~factor(COHORT, c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")), nrow = 2, labeller = cohort_labeler)

# Test the difference in ensemble intrinsic age acceleration between AD and control groups for each cohort
wilcox.test(ensemble_intrinsic_age_accel ~ STATUS, data = filter(bio_age_global_anc, COHORT == "NHW"))
wilcox.test(ensemble_intrinsic_age_accel ~ STATUS, data = filter(bio_age_global_anc, COHORT == "REAAADI"))
wilcox.test(ensemble_intrinsic_age_accel ~ STATUS, data = filter(bio_age_global_anc, COHORT == "CuADI"))
wilcox.test(ensemble_intrinsic_age_accel ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PRADI"))
wilcox.test(ensemble_intrinsic_age_accel ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PERUVIAN"))


### Compare age acceleration values for multiple clocks
ggplot(data = bio_age_global_anc) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc.Hannum, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.Hannum), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylim(NA, 40) +
  ylab("Hannum Clock Raw Age Acceleration") + 
  theme_minimal() +
  facet_wrap(~factor(COHORT, c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")), nrow = 2, labeller = cohort_labeler)

ggplot(data = bio_age_global_anc) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.Hannum, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.Hannum), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylim(NA, 40) +
  ylab("Hannum Clock Intrinsic Age Acceleration") + 
  theme_minimal() +
  facet_wrap(~factor(COHORT, c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")), nrow = 2, labeller = cohort_labeler)

ggplot(data = bio_age_global_anc) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc.EN, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.EN), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylim(NA, 40) +
  ylab("EN Clock Raw Age Acceleration") + 
  theme_minimal() +
  facet_wrap(~factor(COHORT, c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")), nrow = 2, labeller = cohort_labeler)

ggplot(data = bio_age_global_anc) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.EN, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.EN), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylim(NA, 40) +
  ylab("EN Clock Intrinsic Age Acceleration") + 
  theme_minimal() +
  facet_wrap(~factor(COHORT, c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")), nrow = 2, labeller = cohort_labeler)


ggplot(data = bio_age_global_anc) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc.BLUP, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc.BLUP), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylim(NA, 40) +
  ylab("BLUP Clock Raw Age Acceleration") + 
  theme_minimal() +
  facet_wrap(~factor(COHORT, c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")), nrow = 2, labeller = cohort_labeler)

ggplot(data = bio_age_global_anc) +
  geom_boxplot(mapping = aes(x = STATUS, y = ageAcc2.BLUP, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = ageAcc2.BLUP), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Disease Status") +
  ylim(NA, 40) +
  ylab("BLUP Clock Intrinsic Age Acceleration") + 
  theme_minimal() +
  facet_wrap(~factor(COHORT, c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")), nrow = 2, labeller = cohort_labeler)

ggplot(data = bio_age_global_anc) +
  geom_boxplot(mapping = aes(x = STATUS, y = dunedin_preds.DunedinPACE, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = dunedin_preds.DunedinPACE), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = 2) +
  xlab("Disease Status") +
  ylab("DunedinPACE of Aging") + 
  theme_minimal() +
  facet_wrap(~factor(COHORT, c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")), nrow = 2, labeller = cohort_labeler)

ggplot(data = bio_age_global_anc) + # One more plot I forgot to make before, but DunedinPACE for AD vs Control regardless of population.
  geom_boxplot(mapping = aes(x = STATUS, y = dunedin_preds.DunedinPACE, fill = STATUS)) +
  geom_jitter(mapping = aes(x = STATUS, y = dunedin_preds.DunedinPACE), width = 0.2, alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = 2) +
  xlab("Disease Status") +
  ylab("DunedinPACE of Aging") + 
  theme_minimal()
wilcox.test(dunedin_preds.DunedinPACE ~ STATUS, data = bio_age_global_anc)

# ggpubr for Horvath Clock comparison
fig3a <- ggboxplot(
  data = bio_age_global_anc %>% 
    filter(COHORT != "NHW") %>% 
    mutate(STATUS = factor(STATUS, levels = c("AD", "CONTROL"))),  # Reorder STATUS
  x = "STATUS",
  y = "ageAcc2.Horvath",
  fill = "STATUS",
  xlab = "Disease Status",
  ylab = "Intrinsic Age Acceleration (Horvath DNAmAge)",
) +
  geom_jitter(
    mapping = aes(x = STATUS, y = ageAcc2.Horvath),
    width = 0.2, alpha = 0.7
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~factor(COHORT, c("REAAADI", "PRADI", "CuADI", "PERUVIAN")), 
             nrow = 2, labeller = cohort_labeler) +
  ylim(NA, 20) +  # Set y-axis upper limit to 20
  theme(legend.position = "none")



# Do a Mann-Whitney U test to compare the age acceleration values between the AD and control groups for each cohort
# Hannum Clock Raw Age Acceleration
wilcox.test(ageAcc.Hannum ~ STATUS, data = filter(bio_age_global_anc, COHORT == "NHW"))
wilcox.test(ageAcc.Hannum ~ STATUS, data = filter(bio_age_global_anc, COHORT == "REAAADI"))
wilcox.test(ageAcc.Hannum ~ STATUS, data = filter(bio_age_global_anc, COHORT == "CuADI"))
wilcox.test(ageAcc.Hannum ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PRADI"))
wilcox.test(ageAcc.Hannum ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PERUVIAN"))

# Hannum Clock Intrinsic Age Acceleration
wilcox.test(ageAcc2.Hannum ~ STATUS, data = filter(bio_age_global_anc, COHORT == "NHW"))
wilcox.test(ageAcc2.Hannum ~ STATUS, data = filter(bio_age_global_anc, COHORT == "REAAADI"))
wilcox.test(ageAcc2.Hannum ~ STATUS, data = filter(bio_age_global_anc, COHORT == "CuADI"))
wilcox.test(ageAcc2.Hannum ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PRADI"))
wilcox.test(ageAcc2.Hannum ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PERUVIAN"))

# EN Clock Raw Age Acceleration
wilcox.test(ageAcc.EN ~ STATUS, data = filter(bio_age_global_anc, COHORT == "NHW"))
wilcox.test(ageAcc.EN ~ STATUS, data = filter(bio_age_global_anc, COHORT == "REAAADI"))
wilcox.test(ageAcc.EN ~ STATUS, data = filter(bio_age_global_anc, COHORT == "CuADI"))
wilcox.test(ageAcc.EN ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PRADI"))
wilcox.test(ageAcc.EN ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PERUVIAN"))

# EN Clock Intrinsic Age Acceleration
wilcox.test(ageAcc2.EN ~ STATUS, data = filter(bio_age_global_anc, COHORT == "NHW"))
wilcox.test(ageAcc2.EN ~ STATUS, data = filter(bio_age_global_anc, COHORT == "REAAADI"))
wilcox.test(ageAcc2.EN ~ STATUS, data = filter(bio_age_global_anc, COHORT == "CuADI"))
wilcox.test(ageAcc2.EN ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PRADI"))
wilcox.test(ageAcc2.EN ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PERUVIAN"))

# BLUP Clock Raw Age Acceleration
wilcox.test(ageAcc.BLUP ~ STATUS, data = filter(bio_age_global_anc, COHORT == "NHW"))
wilcox.test(ageAcc.BLUP ~ STATUS, data = filter(bio_age_global_anc, COHORT == "REAAADI"))
wilcox.test(ageAcc.BLUP ~ STATUS, data = filter(bio_age_global_anc, COHORT == "CuADI"))
wilcox.test(ageAcc.BLUP ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PRADI"))
wilcox.test(ageAcc.BLUP ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PERUVIAN"))

# BLUP Clock Intrinsic Age Acceleration
wilcox.test(ageAcc2.BLUP ~ STATUS, data = filter(bio_age_global_anc, COHORT == "NHW"))
wilcox.test(ageAcc2.BLUP ~ STATUS, data = filter(bio_age_global_anc, COHORT == "REAAADI"))
wilcox.test(ageAcc2.BLUP ~ STATUS, data = filter(bio_age_global_anc, COHORT == "CuADI"))
wilcox.test(ageAcc2.BLUP ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PRADI"))
wilcox.test(ageAcc2.BLUP ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PERUVIAN"))

# DunedinPACE of Aging
wilcox.test(dunedin_preds.DunedinPACE ~ STATUS, data = filter(bio_age_global_anc, COHORT == "NHW"))
wilcox.test(dunedin_preds.DunedinPACE ~ STATUS, data = filter(bio_age_global_anc, COHORT == "REAAADI"))
wilcox.test(dunedin_preds.DunedinPACE ~ STATUS, data = filter(bio_age_global_anc, COHORT == "CuADI"))
wilcox.test(dunedin_preds.DunedinPACE ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PRADI"))
wilcox.test(dunedin_preds.DunedinPACE ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PERUVIAN"))

t.test(dunedin_preds.DunedinPACE ~ STATUS, data = filter(bio_age_global_anc, COHORT == "NHW"))
t.test(dunedin_preds.DunedinPACE ~ STATUS, data = filter(bio_age_global_anc, COHORT == "REAAADI"))
t.test(dunedin_preds.DunedinPACE ~ STATUS, data = filter(bio_age_global_anc, COHORT == "CuADI"))
t.test(dunedin_preds.DunedinPACE ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PRADI"))
t.test(dunedin_preds.DunedinPACE ~ STATUS, data = filter(bio_age_global_anc, COHORT == "PERUVIAN"))

### Heatmap section
# We'll make a heatmap showing the comparisons of age accelerations between AD and control groups for each cohort, split by clock
# Define a function to calculate t-statistic and p-value
t_test_results <- function(data, clock, cohort) {
  ad_data <- data %>% filter(COHORT == cohort & STATUS == "AD") %>% pull(clock)
  control_data <- data %>% filter(COHORT == cohort & STATUS == "CONTROL") %>% pull(clock)
  
  # Check if there are enough data points for the t-test
  if(length(ad_data) < 2 | length(control_data) < 2) {
    return(data.frame(CLOCK = clock, COHORT = cohort, t_statistic = NA, p_value = NA))
  }
  
  res <- t.test(ad_data, control_data)
  return(data.frame(CLOCK = clock, COHORT = cohort, t_statistic = res$statistic, p_value = res$p.value))
}

# List of clocks and cohorts
clocks <- c("ageAcc2.Horvath", "ageAcc2.Hannum", "ageAcc2.EN", "ageAcc2.Levine", "dunedin_preds.DunedinPACE")
cohorts <- unique(bio_age_global_anc$COHORT)

# Calculate t-statistics and p-values for each clock and cohort
results <- do.call(rbind, lapply(clocks, function(clock) {
  do.call(rbind, lapply(cohorts, function(cohort) {
    t_test_results(bio_age_global_anc, clock, cohort)
  }))
}))

# Adjust p-values for multiple comparisons (optional)
results$p_value_adj <- p.adjust(results$p_value, method = "fdr")

# Add significance markers
results$significance <- ifelse(results$p_value < 0.05, "*", "")

# Relabel cohorts and clocks
cohort_labels <- c("NHW" = "Non-Hispanic White", 
                   "REAAADI" = "African American", 
                   "CuADI" = "Cuban", 
                   "PRADI" = "Puerto Rican", 
                   "PERUVIAN" = "Peruvian")

clock_labels <- c("ageAcc2.Horvath" = "Horvath", 
                  "ageAcc2.Hannum" = "Hannum", 
                  "ageAcc2.EN" = "Zhang2019_EN",
                  "ageAcc2.Levine" = "PhenoAge", 
                  "dunedin_preds.DunedinPACE" = "DunedinPACE")

results$COHORT <- cohort_labels[results$COHORT]
results$CLOCK <- clock_labels[results$CLOCK]

# Print results for inspection
print(results)

# Create a matrix of t-statistics for the heatmap
t_stat_matrix <- results %>%
  select(CLOCK, COHORT, t_statistic) %>%
  pivot_wider(names_from = COHORT, values_from = t_statistic) %>%
  column_to_rownames("CLOCK") %>%
  as.matrix()

# Create a matrix of significance markers for the heatmap
significance_matrix <- results %>%
  select(CLOCK, COHORT, significance) %>%
  pivot_wider(names_from = COHORT, values_from = significance) %>%
  column_to_rownames("CLOCK") %>%
  as.matrix()

# Print matrices for inspection
print(t_stat_matrix)
print(significance_matrix)


# Define colors for the heatmap
breaks <- seq(min(t_stat_matrix, na.rm = TRUE), max(t_stat_matrix, na.rm = TRUE), length.out = 100)
colors <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

# Create the heatmap
pheatmap(t_stat_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colors,
         breaks = breaks,
         display_numbers = significance_matrix,
         number_color = "black",
         main = "Heatmap of T-Statistics for Age Acceleration Clocks",
         fontsize_number = 25,
         legend = TRUE,
         annotation_legend = TRUE)

# Now the same, but coloring for the difference in median intrinsic age acceleration between AD and control groups
# Define a function to calculate median differences and p-values
median_diff_results <- function(data, clock, cohort) {
  ad_data <- data %>% filter(COHORT == cohort & STATUS == "AD") %>% pull(clock)
  control_data <- data %>% filter(COHORT == cohort & STATUS == "CONTROL") %>% pull(clock)
  
  # Check if there are enough data points
  if(length(ad_data) < 2 | length(control_data) < 2) {
    return(data.frame(CLOCK = clock, COHORT = cohort, median_diff = NA, p_value = NA))
  }
  
  median_diff <- median(ad_data) - median(control_data)
  res <- t.test(ad_data, control_data) # Non-parametric test for median difference
  return(data.frame(CLOCK = clock, COHORT = cohort, median_diff = median_diff, p_value = res$p.value))
}

# List of clocks and cohorts with relabeled names
clocks <- c("ageAcc2.Horvath", "ageAcc2.Hannum", "ageAcc2.EN", "ageAcc2.Levine", "dunedin_preds.DunedinPACE")
cohorts <- c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")

# Calculate median differences and p-values for each clock and cohort
results <- do.call(rbind, lapply(clocks, function(clock) {
  do.call(rbind, lapply(cohorts, function(cohort) {
    median_diff_results(bio_age_global_anc, clock, cohort)
  }))
}))

# Adjust p-values for multiple comparisons (optional)
results$p_value_adj <- p.adjust(results$p_value, method = "bonferroni")

# Add significance markers
results$significance <- ifelse(results$p_value < 0.05, "*", "")

# Relabel cohorts and clocks
cohort_labels <- c("NHW" = "Non-Hispanic White", 
                   "REAAADI" = "African American", 
                   "CuADI" = "Cuban", 
                   "PRADI" = "Puerto Rican", 
                   "PERUVIAN" = "Peruvian")

clock_labels <- c("ageAcc2.Horvath" = "Horvath", 
                  "ageAcc2.Hannum" = "Hannum", 
                  "ageAcc2.Levine" = "PhenoAge", 
                  "ageAcc2.EN" = "Zhang_EN",
                  "dunedin_preds.DunedinPACE" = "DunedinPACE")

results$COHORT <- cohort_labels[results$COHORT]
results$CLOCK <- clock_labels[results$CLOCK]

# Create a matrix of median differences for the heatmap
median_diff_matrix <- results %>%
  select(CLOCK, COHORT, median_diff) %>%
  pivot_wider(names_from = COHORT, values_from = median_diff) %>%
  column_to_rownames("CLOCK") %>%
  as.matrix()

# Create a matrix of significance markers for the heatmap
significance_matrix <- results %>%
  select(CLOCK, COHORT, significance) %>%
  pivot_wider(names_from = COHORT, values_from = significance) %>%
  column_to_rownames("CLOCK") %>%
  as.matrix()

# Print matrices for inspection
print(median_diff_matrix)
print(significance_matrix)


# Define the range for the color scale to be symmetric around 0
max_abs_value <- max(abs(median_diff_matrix), na.rm = TRUE)
breaks <- seq(-0.5, 0.5, length.out = 100)
colors <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

# Create the heatmap
fig3b <- pheatmap(median_diff_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colors,
         breaks = breaks,
         display_numbers = significance_matrix,
         number_color = "black",
         main = "Heatmap of Median Differences in Intrinsic Age Acceleration",
         fontsize_number = 25,
         legend = TRUE,
         annotation_legend = TRUE)


# Do the same approach but dividing the median age acceleration for AD by the median age acceleration for controls
# Define a function to calculate median ratios and p-values
median_ratio_results <- function(data, clock, cohort) {
  ad_data <- data %>% filter(COHORT == cohort & STATUS == "AD") %>% pull(clock)
  control_data <- data %>% filter(COHORT == cohort & STATUS == "CONTROL") %>% pull(clock)
  
  # Check if there are enough data points
  if(length(ad_data) < 2 | length(control_data) < 2) {
    return(data.frame(CLOCK = clock, COHORT = cohort, median_ratio = NA, p_value = NA))
  }
  
  median_ratio <- median(ad_data) / median(control_data)
  res <- wilcox.test(ad_data, control_data) # Non-parametric test for median ratio
  return(data.frame(CLOCK = clock, COHORT = cohort, median_ratio = median_ratio, p_value = res$p.value))
}

# List of clocks and cohorts with relabeled names
clocks <- c("ageAcc2.Horvath", "ageAcc2.Hannum", "ageAcc2.BLUP", "ageAcc2.Levine", "dunedin_preds.DunedinPACE")
cohorts <- c("NHW", "REAAADI", "CuADI", "PRADI", "PERUVIAN")

# Calculate median ratios and p-values for each clock and cohort
results <- do.call(rbind, lapply(clocks, function(clock) {
  do.call(rbind, lapply(cohorts, function(cohort) {
    median_ratio_results(bio_age_global_anc, clock, cohort)
  }))
}))

# Adjust p-values for multiple comparisons (optional)
results$p_value_adj <- p.adjust(results$p_value, method = "bonferroni")

# Add significance markers
results$significance <- ifelse(results$p_value < 0.05, "*", "")

# Relabel cohorts and clocks
cohort_labels <- c("NHW" = "Non-Hispanic White", 
                   "REAAADI" = "African American", 
                   "CuADI" = "Cuban", 
                   "PRADI" = "Puerto Rican", 
                   "PERUVIAN" = "Peruvian")

clock_labels <- c("ageAcc2.Horvath" = "Horvath",
                  "ageAcc2.Hannum" = "Hannum", 
                  "ageAcc2.BLUP" = "Zhang_BLUP",
                  "ageAcc2.Levine" = "PhenoAge", 
                  "dunedin_preds.DunedinPACE" = "DunedinPACE")

results$COHORT <- cohort_labels[results$COHORT]
results$CLOCK <- clock_labels[results$CLOCK]

# Create a matrix of median ratios for the heatmap
median_ratio_matrix <- results %>%
  select(CLOCK, COHORT, median_ratio) %>%
  pivot_wider(names_from = COHORT, values_from = median_ratio) %>%
  column_to_rownames("CLOCK") %>%
  as.matrix()

# Create a matrix of significance markers for the heatmap
significance_matrix <- results %>%
  select(CLOCK, COHORT, significance) %>%
  pivot_wider(names_from = COHORT, values_from = significance) %>%
  column_to_rownames("CLOCK") %>%
  as.matrix()

# Print matrices for inspection
print(median_ratio_matrix)
print(significance_matrix)


# Define the range for the color scale to be symmetric around 1
max_abs_value <- max(abs(median_ratio_matrix - 1), na.rm = TRUE)
breaks <- seq(-2, 2, length.out = 100)
colors <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

# Create the heatmap
pheatmap(median_ratio_matrix - 1,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colors,
         breaks = breaks,
         display_numbers = significance_matrix,
         number_color = "black",
         main = "Heatmap of Median Ratios of Intrinsic Age Acceleration",
         fontsize_number = 25,
         legend = TRUE,
         annotation_legend = TRUE)




