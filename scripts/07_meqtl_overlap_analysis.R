#### This is a script to overlap the AFR vs EUR differentially methylated CpG sites with Horvath Clock CpG sites that are affected by meQTL

# Load libraries
library(tidyverse)
library(methylclock)

# Set the directory
setwd("~/Desktop/Capra Lab/Thesis Project/Aim_1/")

# Load the AFR vs EUR differentially methylated sites defined by Oge
load("~/Desktop/Capra Lab/Thesis Project/Aim_1/ancestry_CpGs.RData")

# Load the CpG sites for each clock that are affectd by meQTL
combined_meqtls <- readRDS("all_meqtl_affected_cpgs.rds")
all_meqtl_cpgs <- combined_meqtls$CpG


# Load the MAGENTA results
all_results <- readRDS("all_results.rds")

# Pass individual cohort data
aa_meqtl <- all_results$results$AFR$meqtl_data
pr_meqtl <- all_results$results$PR$meqtl_data
cub_meqtl <- all_results$results$CUB$meqtl_data
per_meqtl <- all_results$results$PER$meqtl_data

# Remove all_results object from memory to alleviate RAM
rm(all_results)

##### Filter AFR vs EUR CpGs
# CpGs used in Horvath clock
load_DNAm_Clocks_data()
horvath_cpgs <- unique(coefHorvath$CpGmarker)

# AFR vs EUR filtered by Horvath CpGs
#afr_vs_eur_filtered <- AFR_vs_EUR[rownames(AFR_vs_EUR) %in% horvath_cpgs, ]

# Apply significance filter
afr_vs_eur_filtered <- subset(AFR_vs_EUR, adj.P.Val < 0.05)

# Final AFR vs EUR CpG list
afr_dm_cpgs <- rownames(afr_vs_eur_filtered)

##### Gather meQTL CpGs
# All meQTL CpGs across cohorts
horvath_meqtl_cpgs <- unique(c(
  aa_meqtl$CpG, 
  pr_meqtl$PROBE, 
  cub_meqtl$PROBE, 
  per_meqtl$PROBE
))


# AFR-differentiated meQTL CpGs
afr_meqtl_cpgs <- unique(c(
  subset(aa_meqtl, african_differentiated == TRUE)$CpG,
  subset(pr_meqtl, african_differentiated == TRUE)$PROBE,
  subset(cub_meqtl, african_differentiated == TRUE)$PROBE,
  subset(per_meqtl, african_differentiated == TRUE)$PROBE
))


######## Proportion plots
#### Proportion of clock CpGs that are differentially methylated in AFR individuals relative to EUR individuals
library(ggpubr)
library(stringr)

# 1) load tidy data
df_long <- readr::read_csv("clock_cpgs_overlap_tidy.csv")

# keep the original Clock ordering as it appears in the file
original_levels <- unique(df_long$Clock)

# 2) parse totals from the Clock strings like "Horvath(n=353)"
df_totals <- df_long %>%
  distinct(Clock) %>%
  mutate(ParsedTotal = as.numeric(str_match(Clock, "n\\s*=\\s*(\\d+)")[,2]))

# 3) fallback totals: if parsing failed, use the number of unique CpGs present in the tidy file
counts_from_data <- df_long %>%
  group_by(Clock) %>%
  summarise(CountRows = n_distinct(CpGs), .groups = "drop")

df_totals <- df_totals %>%
  left_join(counts_from_data, by = "Clock") %>%
  mutate(Total = if_else(!is.na(ParsedTotal), ParsedTotal, CountRows)) %>%
  select(Clock, Total)

# 4) count WB AFR overlaps (unique CpGs)
wb_counts <- df_long %>%
  filter(Category == "WB AFR ancestry-CpGs") %>%
  group_by(Clock) %>%
  summarise(WB_AFR = n_distinct(CpGs), .groups = "drop")

# 5) summary table with percentage (using the hard-coded Total)
df_summary <- df_totals %>%
  left_join(wb_counts, by = "Clock") %>%
  mutate(WB_AFR = replace_na(WB_AFR, 0),
         Percentage = 100 * WB_AFR / Total,
         Clock = factor(Clock, levels = original_levels),
         label_pct = sprintf("%.1f%%", Percentage),
         label_counts = paste0(WB_AFR, "/", Total))
df_summary <- df_summary %>%
  mutate(Clock = str_replace(Clock, "ZhangEN\\(n=514\\)", "EN"),
         Clock = str_replace(Clock, "Horvath\\(n=353\\)", "Horvath"),
         Clock = str_replace(Clock, "Hannum\\(n=71\\)", "Hannum"),
         Clock = str_replace(Clock, "PhenoAge\\(n=513\\)", "PhenoAge"),
         Clock = str_replace(Clock, "DunedinPACE\\(n=173\\)", "DunedinPACE"))

# Set custom factor levels for plotting order
df_summary <- df_summary %>%
  mutate(Clock = factor(Clock,
                        levels = c("Horvath",
                                   "Hannum",
                                   "EN",
                                   "PhenoAge",
                                   "DunedinPACE")))

# 6) plotting: percent above each bar, counts (WB_AFR/Total) above that
ymax <- max(df_summary$Percentage, na.rm = TRUE)
if (is.na(ymax) || ymax == 0) ymax <- 1
ylim_max <- ymax * 1.25 + 5
offset <- ymax * 0.06 + 1


p <- ggbarplot(df_summary,
               x = "Clock",
               y = "Percentage",
               color = "black",
               palette = "jco",
               sort.val = "none") +
  geom_text(aes(label = label_pct, y = Percentage),
            vjust = -0.5, size = 5) +      # bigger percent labels
  geom_text(aes(label = label_counts, y = Percentage + offset),
            vjust = -0.5, size = 5) +      # bigger counts
  scale_y_continuous(limits = c(0, ylim_max), expand = c(0, 0)) +
  labs(y = "Percent of Clock CpGs with Differential Methylation in African Ancestry Individuals", x = "Clock") +
  theme_pubr(base_size = 16) +             # increase all text sizes
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))


# show + save
print(p)
saveRDS(p, "plots_RDS/fig5a.rds")

#### Proportion of clock CpGs that are associated with environmental exposures
# 1) load tidy data
df_long <- readr::read_csv("clock_cpgs_overlap_tidy.csv")

# keep the original Clock ordering as it appears in the file
original_levels <- unique(df_long$Clock)

# 2) parse totals from the Clock strings like "Horvath(n=353)"
df_totals <- df_long %>%
  distinct(Clock) %>%
  mutate(ParsedTotal = as.numeric(str_match(Clock, "n\\s*=\\s*(\\d+)")[,2]))

# 3) fallback totals: if parsing failed, use the number of unique CpGs present in the tidy file
counts_from_data <- df_long %>%
  group_by(Clock) %>%
  summarise(CountRows = n_distinct(CpGs), .groups = "drop")

df_totals <- df_totals %>%
  left_join(counts_from_data, by = "Clock") %>%
  mutate(Total = if_else(!is.na(ParsedTotal), ParsedTotal, CountRows)) %>%
  select(Clock, Total)

# 4) count all three categories (unique CpGs)
category_counts <- df_long %>%
  group_by(Clock, Category) %>%
  summarise(Count = n_distinct(CpGs), .groups = "drop") %>%
  pivot_wider(names_from = Category, values_from = Count, values_fill = 0)

# 5) create summary with all categories and their proportions
df_summary <- df_totals %>%
  left_join(category_counts, by = "Clock") %>%
  # Replace NAs with 0 for any missing categories
  mutate(across(c("All exposure-associated CpGs", 
                  "WB AFR ancestry-CpGs", 
                  "Overlapping ancestry and exposure associated-CpGs"), 
                ~replace_na(.x, 0))) %>%
  # Calculate proportions (as percentages)
  mutate(
    Prop_Exposure = 100 * `All exposure-associated CpGs` / Total,
    Prop_WB_AFR = 100 * `WB AFR ancestry-CpGs` / Total,
    Prop_Overlap = 100 * `Overlapping ancestry and exposure associated-CpGs` / Total
  )

# Clean up clock names
df_summary <- df_summary %>%
  mutate(Clock = str_replace(Clock, "ZhangEN\\(n=514\\)", "EN"),
         Clock = str_replace(Clock, "Horvath\\(n=353\\)", "Horvath"),
         Clock = str_replace(Clock, "Hannum\\(n=71\\)", "Hannum"),
         Clock = str_replace(Clock, "PhenoAge\\(n=513\\)", "PhenoAge"),
         Clock = str_replace(Clock, "DunedinPACE\\(n=173\\)", "DunedinPACE"))

# Set custom factor levels for plotting order
df_summary <- df_summary %>%
  mutate(Clock = factor(Clock,
                        levels = c("Horvath", "Hannum", "EN", "PhenoAge", "DunedinPACE")))

# 6) Reshape data for stacked bar plot
df_plot <- df_summary %>%
  select(Clock, Prop_Exposure, Prop_WB_AFR, Prop_Overlap) %>%
  pivot_longer(cols = starts_with("Prop_"), 
               names_to = "Category", 
               values_to = "Percentage") %>%
  mutate(Category = case_when(
    Category == "Prop_Exposure" ~ "All exposure-associated CpGs",
    Category == "Prop_WB_AFR" ~ "WB AFR ancestry-CpGs",
    Category == "Prop_Overlap" ~ "Overlapping ancestry and exposure associated-CpGs"
  )) %>%
  mutate(Category = factor(Category, 
                           levels = c("All exposure-associated CpGs",
                                      "WB AFR ancestry-CpGs", 
                                      "Overlapping ancestry and exposure associated-CpGs")))

# 7) Create stacked bar plot
p <- ggplot(df_plot, aes(x = Clock, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  scale_fill_manual(values = c("All exposure-associated CpGs" = "#E64B35FF",
                               "WB AFR ancestry-CpGs" = "#4DBBD5FF", 
                               "Overlapping ancestry and exposure associated-CpGs" = "#00A087FF")) +
  labs(y = "Percentage (%)", 
       x = "Clock",
       fill = "CpG Category") +
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_legend(nrow = 3)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 40))

# Alternative version with percentage labels on each segment
p_with_labels <- p +
  geom_text(data = df_plot %>% 
              filter(Percentage > 2), # Only show labels for segments > 2%
            aes(label = sprintf("%.1f%%", Percentage)),
            position = position_stack(vjust = 0.5),
            size = 6, color = "white", fontface = "bold")

# Show both versions
print(p)
print(p_with_labels)
saveRDS(p_with_labels, "plots_RDS/fig5f.rds")

# Optional: Create a summary table showing the actual counts and percentages
summary_table <- df_summary %>%
  select(Clock, Total, 
         `All exposure-associated CpGs`, Prop_Exposure,
         `WB AFR ancestry-CpGs`, Prop_WB_AFR,
         `Overlapping ancestry and exposure associated-CpGs`, Prop_Overlap) %>%
  mutate(across(starts_with("Prop_"), ~round(.x, 1))) %>%
  rename("Total CpGs" = Total,
         "Exposure Count" = `All exposure-associated CpGs`,
         "Exposure %" = Prop_Exposure,
         "WB AFR Count" = `WB AFR ancestry-CpGs`, 
         "WB AFR %" = Prop_WB_AFR,
         "Overlap Count" = `Overlapping ancestry and exposure associated-CpGs`,
         "Overlap %" = Prop_Overlap)

print("Summary Table:")
print(summary_table)




#### Proportion of clock CpGs that are associated with environmental exposures
# 1) load tidy data
df_long <- readr::read_csv("clock_cpgs_overlap_tidy.csv")

# keep the original Clock ordering as it appears in the file
original_levels <- unique(df_long$Clock)

# 2) parse totals from the Clock strings like "Horvath(n=353)"
df_totals <- df_long %>%
  distinct(Clock) %>%
  mutate(ParsedTotal = as.numeric(str_match(Clock, "n\\s*=\\s*(\\d+)")[,2]))

# 3) fallback totals: if parsing failed, use the number of unique CpGs present in the tidy file
counts_from_data <- df_long %>%
  group_by(Clock) %>%
  summarise(CountRows = n_distinct(CpGs), .groups = "drop")

df_totals <- df_totals %>%
  left_join(counts_from_data, by = "Clock") %>%
  mutate(Total = if_else(!is.na(ParsedTotal), ParsedTotal, CountRows)) %>%
  select(Clock, Total)

# 4) count all three categories (unique CpGs)
category_counts <- df_long %>%
  group_by(Clock, Category) %>%
  summarise(Count = n_distinct(CpGs), .groups = "drop") %>%
  pivot_wider(names_from = Category, values_from = Count, values_fill = 0)

# 5) compute summary with proportions
df_summary <- df_totals %>%
  left_join(category_counts, by = "Clock") %>%
  mutate(`All exposure-associated CpGs` = replace_na(`All exposure-associated CpGs`, 0)) %>%
  mutate(Prop_Exposure = 100 * `All exposure-associated CpGs` / Total)

# Clean up clock names
df_summary <- df_summary %>%
  mutate(Clock = str_replace(Clock, "ZhangEN\\(n=514\\)", "EN"),
         Clock = str_replace(Clock, "Horvath\\(n=353\\)", "Horvath"),
         Clock = str_replace(Clock, "Hannum\\(n=71\\)", "Hannum"),
         Clock = str_replace(Clock, "PhenoAge\\(n=513\\)", "PhenoAge"),
         Clock = str_replace(Clock, "DunedinPACE\\(n=173\\)", "DunedinPACE"),
         Clock = factor(Clock, 
                        levels = c("Horvath", "Hannum", "EN", "PhenoAge", "DunedinPACE")))

# 6) plot only exposure-associated CpGs
p_exposure <- ggplot(df_summary, aes(x = Clock, y = Prop_Exposure, fill = Clock)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", Prop_Exposure)), 
            vjust = -0.5, size = 6, fontface = "bold") +
  scale_fill_manual(values = c("Horvath" = "#E64B35FF",
                               "Hannum" = "#E64B35FF",
                               "EN" = "#E64B35FF",
                               "PhenoAge" = "#E64B35FF",
                               "DunedinPACE" = "#E64B35FF")) +
  labs(y = "Percentage (%)", 
       x = "Clock",
       title = "Proportion of clock CpGs that are exposure associated") +
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 20))

# Print plot
print(p_exposure)



# Separate code chunk for Horvath clock meQTL vs Environmental stacked bar plot
# STEP 1: Use your existing meQTL-affected CpGs vector
# Assumes you have a vector called 'meqtl_cpgs' in your environment
# STEP 2: Filter data for specific clocks
# Extract unique CpGs from each clock's coefficient dataframe
horvath_all_cpgs <- unique(coefHorvath$CpGmarker)
hannum_all_cpgs <- unique(coefHannum$CpGmarker)
en_all_cpgs <- unique(coefEN$CpGmarker)  # adjust name if needed
phenoage_all_cpgs <- unique(coefLevine$CpGmarker)  # adjust name if needed

# Check the counts
length(horvath_all_cpgs)
length(hannum_all_cpgs)
length(en_all_cpgs)
length(phenoage_all_cpgs)

# If you still need to filter df_long by these CpGs:
horvath_data <- df_long %>%
  filter(CpGs %in% horvath_all_cpgs)

hannum_data <- df_long %>%
  filter(CpGs %in% hannum_all_cpgs)

en_data <- df_long %>%
  filter(CpGs %in% en_all_cpgs)

phenoage_data <- df_long %>%
  filter(CpGs %in% phenoage_all_cpgs)



# Get environmentally associated CpGs for each clock
horvath_env_cpgs <- horvath_data %>%
  filter(Category == "All exposure-associated CpGs") %>%
  pull(CpGs) %>%
  unique()
hannum_env_cpgs <- hannum_data %>%
  filter(Category == "All exposure-associated CpGs") %>%
  pull(CpGs) %>%
  unique()
en_env_cpgs <- en_data %>%
  filter(Category == "All exposure-associated CpGs") %>%
  pull(CpGs) %>%
  unique()
phenoage_env_cpgs <- phenoage_data %>%
  filter(Category == "All exposure-associated CpGs") %>%
  pull(CpGs) %>%
  unique()

# STEP 3: Get total from the parsed clock name
horvath_total <- df_totals %>%
  filter(str_detect(Clock, "Horvath")) %>%
  pull(Total)
hannum_total <- df_totals %>%
  filter(str_detect(Clock, "Hannum")) %>%
  pull(Total)
en_total <- df_totals %>%
  filter(str_detect(Clock, "ZhangEN")) %>%
  pull(Total)
phenoage_total <- df_totals %>%
  filter(str_detect(Clock, "PhenoAge")) %>%
  pull(Total)

# STEP 4: Classify CpGs into categories (NEED TO ADJUST THE MEQTL CPGS DF HERE TO BE THE ONES I JUST GENERATED TODAY (09/29))
# Horvath clock
horvath_meqtl_only <- setdiff(all_meqtl_cpgs, horvath_env_cpgs)
horvath_afr_meqtl_only <- setdiff(afr_meqtl_cpgs, horvath_env_cpgs)
horvath_env_only <- setdiff(horvath_env_cpgs, all_meqtl_cpgs)
horvath_both <- intersect(all_meqtl_cpgs, horvath_env_cpgs)
horvath_afr_both <- intersect(afr_meqtl_cpgs, horvath_env_cpgs)
horvath_neither <- setdiff(horvath_all_cpgs, union(all_meqtl_cpgs, horvath_env_cpgs))

# Hannum clock
hannum_meqtl_only <- setdiff(all_meqtl_cpgs, hannum_env_cpgs)
hannum_afr_meqtl_only <- setdiff(afr_meqtl_cpgs, hannum_env_cpgs)
hannum_env_only <- setdiff(hannum_env_cpgs, all_meqtl_cpgs)
hannum_both <- intersect(all_meqtl_cpgs, hannum_env_cpgs)
hannum_afr_both <- intersect(afr_meqtl_cpgs, hannum_env_cpgs)
hannum_neither <- setdiff(hannum_all_cpgs, union(all_meqtl_cpgs, hannum_env_cpgs))

# EN clock
en_meqtl_only <- setdiff(all_meqtl_cpgs, en_env_cpgs)
en_afr_meqtl_only <- setdiff(afr_meqtl_cpgs, en_env_cpgs)
en_env_only <- setdiff(en_env_cpgs, all_meqtl_cpgs)
en_both <- intersect(all_meqtl_cpgs, en_env_cpgs)
en_afr_both <- intersect(afr_meqtl_cpgs, en_env_cpgs)
en_neither <- setdiff(en_all_cpgs, union(all_meqtl_cpgs, en_env_cpgs))

# PhenoAge clock
phenoage_meqtl_only <- setdiff(all_meqtl_cpgs, phenoage_env_cpgs)
phenoage_afr_meqtl_only <- setdiff(afr_meqtl_cpgs, phenoage_env_cpgs)
phenoage_env_only <- setdiff(phenoage_env_cpgs, all_meqtl_cpgs)
phenoage_both <- intersect(all_meqtl_cpgs, phenoage_env_cpgs)
phenoage_afr_both <- intersect(afr_meqtl_cpgs, phenoage_env_cpgs)
phenoage_neither <- setdiff(phenoage_all_cpgs, union(all_meqtl_cpgs, phenoage_env_cpgs))

# Create summary dataframe
horvath_summary <- data.frame(
  Clock = "Horvath",
  meQTL_only = length(horvath_meqtl_only),
  afr_meqtl_only = length(horvath_afr_meqtl_only),
  Environmental_only = length(horvath_env_only),
  afr_Both = length(horvath_afr_both),
  Both = length(horvath_both),
  Total = 353
)

# Hannum summary
hannum_summary <- data.frame(
  Clock = "Hannum",
  meQTL_only = length(hannum_meqtl_only),
  afr_meqtl_only = length(hannum_afr_meqtl_only),
  Environmental_only = length(hannum_env_only),
  afr_Both = length(hannum_afr_both),
  Both = length(hannum_both),
  Total = 71
)

# EN summary
en_summary <- data.frame(
  Clock = "EN",
  meQTL_only = length(en_meqtl_only),
  afr_meqtl_only = length(en_afr_meqtl_only),
  Environmental_only = length(en_env_only),
  afr_Both = length(en_afr_both),
  Both = length(en_both),
  Total = 514
)

# PhenoAge summary
phenoage_summary <- data.frame(
  Clock = "PhenoAge",
  meQTL_only = length(phenoage_meqtl_only),
  afr_meqtl_only = length(phenoage_afr_meqtl_only),
  Environmental_only = length(phenoage_env_only),
  afr_Both = length(phenoage_afr_both),
  Both = length(phenoage_both),
  Total = 513
)

# STEP 5: Calculate proportions and prepare plot data
# Function to prepare data for any clock summary table
prep_plot_data <- function(summary_df) {
  summary_df %>%
    mutate(
      Prop_meQTL_only = (253 / Total) * 100,
      Prop_afr_meqtl_only = 100 * (afr_meqtl_only / Total),
      Prop_Environmental_only = 100 * (Environmental_only / Total),
      Prop_afr_Both = 100 * (afr_Both / Total),
      Prop_Both = 100 * (Both / Total)
    ) %>%
    select(Clock, Prop_meQTL_only, Prop_afr_meqtl_only, 
           Prop_Environmental_only, Prop_Both, Prop_afr_Both) %>%
    pivot_longer(cols = starts_with("Prop_"), 
                 names_to = "Category", 
                 values_to = "Percentage") %>%
    mutate(Category = case_when(
      Category == "Prop_meQTL_only" ~ "meQTL-affected only",
      Category == "Prop_afr_meqtl_only" ~ "AFR meQTL-affected only",
      Category == "Prop_Environmental_only" ~ "Environmentally associated only",
      Category == "Prop_afr_Both" ~ "Both AFR meQTL-affected and environmentally associated",
      Category == "Prop_Both" ~ "Both meQTL-affected and environmentally associated"
    )) %>%
    mutate(Category = factor(Category, 
                             levels = c("meQTL-affected only",
                                        "AFR meQTL-affected only",
                                        "Environmentally associated only",
                                        "Both AFR meQTL-affected and environmentally associated",
                                        "Both meQTL-affected and environmentally associated"))) %>%
    # Keep only the three categories of interest
    filter(Category %in% c("AFR meQTL-affected only", 
                           "Environmentally associated only", 
                           "Both AFR meQTL-affected and environmentally associated"))
}

# Prepare datasets for all clocks
horvath_plot_data  <- prep_plot_data(horvath_summary)
hannum_plot_data   <- prep_plot_data(hannum_summary)
en_plot_data       <- prep_plot_data(en_summary)
phenoage_plot_data <- prep_plot_data(phenoage_summary)

# Combine all into one dataframe
all_clocks_plot_data <- bind_rows(
  horvath_plot_data,
  hannum_plot_data,
  en_plot_data,
  phenoage_plot_data
)
# STEP 6: Create the stacked bar plot
# Single stacked horizontal barplot
all_clocks_meqtl_plot <- ggbarplot(all_clocks_plot_data,
                                   x = "Clock",
                                   y = "Percentage",
                                   fill = "Category",
                                   color = "black",
                                   palette = c("#4DBBD5FF", "#00A087FF", "#E64B35FF"),
                                   width = 0.6,
                                   orientation = "horiz") +
  geom_text(data = all_clocks_plot_data %>% filter(Percentage >= 0.5),  # only keep ≥ 0.5%
            aes(x = Clock, y = Percentage, 
                label = sprintf("%.1f%%", Percentage), 
                fill = Category),
            position = position_stack(vjust = 0.5),
            size = 4.5, color = "white", fontface = "bold", show.legend = FALSE) +
  labs(x = "Percentage (%)", 
       y = "",
       fill = "Association Type") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(nrow = 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20))


print(all_clocks_meqtl_plot)

saveRDS(all_clocks_meqtl_plot, "plots_RDS/fig5g.rds")




####### Same thing but for AFR differentially methylated clock CpGs
# Filter down to CpG sites that are differentially methylated in AFR individuals relative to EUR individuals (based on whole blood samples)
horvath_all_cpgs <- unique(horvath_data$CpGs)
hannum_all_cpgs <- unique(hannum_data$CpGs)
en_all_cpgs <- unique(en_data$CpGs)
phenoage_all_cpgs <- unique(phenoage_data$CpGs)

horvath_diff_cpgs <- horvath_all_cpgs[horvath_all_cpgs %in% afr_dm_cpgs]
hannum_diff_cpgs <- hannum_all_cpgs[hannum_all_cpgs %in% afr_dm_cpgs]
en_diff_cpgs <- en_all_cpgs[hannum_all_cpgs %in% afr_dm_cpgs]
phenoage_diff_cpgs <- phenoage_all_cpgs[phenoage_all_cpgs %in% afr_dm_cpgs]


# Get environmentally associated CpGs for each clock
horvath_env_cpgs <- horvath_data %>%
  filter(Category == "All exposure-associated CpGs") %>%
  pull(CpGs) %>%
  unique()
hannum_env_cpgs <- hannum_data %>%
  filter(Category == "All exposure-associated CpGs") %>%
  pull(CpGs) %>%
  unique()
en_env_cpgs <- en_data %>%
  filter(Category == "All exposure-associated CpGs") %>%
  pull(CpGs) %>%
  unique()
phenoage_env_cpgs <- phenoage_data %>%
  filter(Category == "All exposure-associated CpGs") %>%
  pull(CpGs) %>%
  unique()

# STEP 3: Get total from the parsed clock name
horvath_total <- df_totals %>%
  filter(str_detect(Clock, "Horvath")) %>%
  pull(Total)
hannum_total <- df_totals %>%
  filter(str_detect(Clock, "Hannum")) %>%
  pull(Total)
en_total <- df_totals %>%
  filter(str_detect(Clock, "ZhangEN")) %>%
  pull(Total)
phenoage_total <- df_totals %>%
  filter(str_detect(Clock, "PhenoAge")) %>%
  pull(Total)

# STEP 4: Classify CpGs into categories (NEED TO ADJUST THE MEQTL CPGS DF HERE TO BE THE ONES I JUST GENERATED TODAY (09/29))
# Horvath clock
horvath_meqtl_only <- setdiff(all_meqtl_cpgs, horvath_env_cpgs)
horvath_afr_meqtl_only <- setdiff(afr_meqtl_cpgs, horvath_env_cpgs)
horvath_env_only <- setdiff(horvath_env_cpgs, all_meqtl_cpgs)
horvath_both <- intersect(all_meqtl_cpgs, horvath_env_cpgs)
horvath_afr_both <- intersect(afr_meqtl_cpgs, horvath_env_cpgs)
horvath_neither <- setdiff(horvath_all_cpgs, union(all_meqtl_cpgs, horvath_env_cpgs))

# Hannum clock
hannum_meqtl_only <- setdiff(all_meqtl_cpgs, hannum_env_cpgs)
hannum_afr_meqtl_only <- setdiff(afr_meqtl_cpgs, hannum_env_cpgs)
hannum_env_only <- setdiff(hannum_env_cpgs, all_meqtl_cpgs)
hannum_both <- intersect(all_meqtl_cpgs, hannum_env_cpgs)
hannum_afr_both <- intersect(afr_meqtl_cpgs, hannum_env_cpgs)
hannum_neither <- setdiff(hannum_all_cpgs, union(all_meqtl_cpgs, hannum_env_cpgs))

# EN clock
en_meqtl_only <- setdiff(all_meqtl_cpgs, en_env_cpgs)
en_afr_meqtl_only <- setdiff(afr_meqtl_cpgs, en_env_cpgs)
en_env_only <- setdiff(en_env_cpgs, all_meqtl_cpgs)
en_both <- intersect(all_meqtl_cpgs, en_env_cpgs)
en_afr_both <- intersect(afr_meqtl_cpgs, en_env_cpgs)
en_neither <- setdiff(en_all_cpgs, union(all_meqtl_cpgs, en_env_cpgs))

# PhenoAge clock
phenoage_meqtl_only <- setdiff(all_meqtl_cpgs, phenoage_env_cpgs)
phenoage_afr_meqtl_only <- setdiff(afr_meqtl_cpgs, phenoage_env_cpgs)
phenoage_env_only <- setdiff(phenoage_env_cpgs, all_meqtl_cpgs)
phenoage_both <- intersect(all_meqtl_cpgs, phenoage_env_cpgs)
phenoage_afr_both <- intersect(afr_meqtl_cpgs, phenoage_env_cpgs)
phenoage_neither <- setdiff(phenoage_all_cpgs, union(all_meqtl_cpgs, phenoage_env_cpgs))

# Create summary dataframe
horvath_summary <- data.frame(
  Clock = "Horvath",
  meQTL_only = length(horvath_meqtl_only),
  afr_meqtl_only = length(horvath_afr_meqtl_only),
  Environmental_only = length(horvath_env_only),
  afr_Both = length(horvath_afr_both),
  Both = length(horvath_both),
  Total = horvath_total
)

# Hannum summary
hannum_summary <- data.frame(
  Clock = "Hannum",
  meQTL_only = length(hannum_meqtl_only),
  afr_meqtl_only = length(hannum_afr_meqtl_only),
  Environmental_only = length(hannum_env_only),
  afr_Both = length(hannum_afr_both),
  Both = length(hannum_both),
  Total = hannum_total
)

# EN summary
en_summary <- data.frame(
  Clock = "EN",
  meQTL_only = length(en_meqtl_only),
  afr_meqtl_only = length(en_afr_meqtl_only),
  Environmental_only = length(en_env_only),
  afr_Both = length(en_afr_both),
  Both = length(en_both),
  Total = en_total
)

# PhenoAge summary
phenoage_summary <- data.frame(
  Clock = "PhenoAge",
  meQTL_only = length(phenoage_meqtl_only),
  afr_meqtl_only = length(phenoage_afr_meqtl_only),
  Environmental_only = length(phenoage_env_only),
  afr_Both = length(phenoage_afr_both),
  Both = length(phenoage_both),
  Total = phenoage_total
)

# STEP 5: Calculate proportions and prepare plot data
# Function to prepare data for any clock summary table
prep_plot_data <- function(summary_df) {
  summary_df %>%
    mutate(
      Prop_meQTL_only = (253 / Total) * 100,
      Prop_afr_meqtl_only = 100 * (afr_meqtl_only / Total),
      Prop_Environmental_only = 100 * (Environmental_only / Total),
      Prop_afr_Both = 100 * (afr_Both / Total),
      Prop_Both = 100 * (Both / Total)
    ) %>%
    select(Clock, Prop_meQTL_only, Prop_afr_meqtl_only, 
           Prop_Environmental_only, Prop_Both, Prop_afr_Both) %>%
    pivot_longer(cols = starts_with("Prop_"), 
                 names_to = "Category", 
                 values_to = "Percentage") %>%
    mutate(Category = case_when(
      Category == "Prop_meQTL_only" ~ "meQTL-affected only",
      Category == "Prop_afr_meqtl_only" ~ "AFR meQTL-affected only",
      Category == "Prop_Environmental_only" ~ "Environmentally associated only",
      Category == "Prop_afr_Both" ~ "Both AFR meQTL-affected and environmentally associated",
      Category == "Prop_Both" ~ "Both meQTL-affected and environmentally associated"
    )) %>%
    mutate(Category = factor(Category, 
                             levels = c("meQTL-affected only",
                                        "AFR meQTL-affected only",
                                        "Environmentally associated only",
                                        "Both AFR meQTL-affected and environmentally associated",
                                        "Both meQTL-affected and environmentally associated"))) %>%
    # Keep only the three categories of interest
    filter(Category %in% c("AFR meQTL-affected only", 
                           "Environmentally associated only", 
                           "Both AFR meQTL-affected and environmentally associated"))
}

# Prepare datasets for all clocks
horvath_plot_data  <- prep_plot_data(horvath_summary)
hannum_plot_data   <- prep_plot_data(hannum_summary)
en_plot_data       <- prep_plot_data(en_summary)
phenoage_plot_data <- prep_plot_data(phenoage_summary)

# Combine all into one dataframe
all_clocks_plot_data <- bind_rows(
  horvath_plot_data,
  hannum_plot_data,
  en_plot_data,
  phenoage_plot_data
)
# STEP 6: Create the stacked bar plot
# Single stacked horizontal barplot
all_clocks_meqtl_plot <- ggbarplot(all_clocks_plot_data,
                                   x = "Clock",
                                   y = "Percentage",
                                   fill = "Category",
                                   color = "black",
                                   palette = c("#4DBBD5FF", "#00A087FF", "#E64B35FF"),
                                   width = 0.6,
                                   orientation = "horiz") +
  geom_text(data = all_clocks_plot_data %>% filter(Percentage >= 0.5),  # only keep ≥ 0.5%
            aes(x = Clock, y = Percentage, 
                label = sprintf("%.1f%%", Percentage), 
                fill = Category),
            position = position_stack(vjust = 0.5),
            size = 4.5, color = "white", fontface = "bold", show.legend = FALSE) +
  labs(x = "Percentage (%)", 
       y = "",
       fill = "Association Type") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(nrow = 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20))


print(all_clocks_meqtl_plot)

