##### This is a script to collect the results of multiple methylation clocks applied to three external datasets, and to generate figures for the manuscript.

# Set the directory (uncomment/adjust if not running from repository root)
# setwd("~/Desktop/Capra Lab/Thesis Project/Aim_1/")

# Load some libraries
library(tidyverse)
library(patchwork)
library(cocor)
library(boot)
library(ggpubr)

# Load results
swedes <- read_csv("data/replication/horvath_swedes_dnamage_all_results.csv", col_select = -c(1))
genoa <- read_csv("data/replication/genoa_epic_dnamage_all_results.csv", col_select = -c(1))
grady <- read_csv("data/replication/aa_grady_dnamage_results.csv", col_select = -c(1))

# Remove swedes with NA values in age column
swedes <- swedes[!is.na(swedes$age), ]

# Calculate the median absolute error for the Horvath predictions, relative to age
# For Swedes dataset
absolute_errors_swedes <- abs(swedes$Horvath - swedes$age)
median_absolute_error_swedes <- median(absolute_errors_swedes, na.rm = T)

# For GENOA dataset
absolute_errors_genoa <- abs(genoa$Horvath - genoa$age)
median_absolute_error_genoa <- median(absolute_errors_genoa)

# For Grady dataset
absolute_errors_grady <- abs(grady$Horvath - grady$age)
median_absolute_error_grady <- median(absolute_errors_grady)


plot_age_comparison <- function(df, xvar, yvar, xlab, ylab,
                                title = NULL, zoomed = FALSE, show_title = TRUE) {
  # pull variables
  x <- df[[xvar]]
  y <- df[[yvar]]

  # compute stats
  r_val <- round(cor(x, y, use = "complete.obs"), 2)
  mae_val <- round(median(abs(x - y), na.rm = TRUE), 2)
  mse_val <- round(mean((x - y)^2, na.rm = TRUE), 2)

  # aligned label (tight spacing)
  label_padded <- sprintf(
    "%-4s %.2f\n%-4s %.2f\n%-4s %.2f",
    "r =",   r_val,
    "MAE=",  mae_val,
    "MSE=",  mse_val
  )

  # limits and annotation position
  if (zoomed) {
    xlim_vals <- c(45, 100)
    ylim_vals <- c(45, 100)
  } else {
    xlim_vals <- c(0, 100)
    ylim_vals <- c(0, 100)
  }
  x_pos <- xlim_vals[1] + 0.02 * diff(xlim_vals)
  y_pos <- ylim_vals[2] - 0.02 * diff(ylim_vals)

  p <- ggplot(df, aes(x = !!sym(xvar), y = !!sym(yvar))) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_abline(
      slope = 1, intercept = 0,
      color = "gray", linetype = "dashed", size = 0.8
    ) +
    annotate("text",
      x = x_pos, y = y_pos,
      label = label_padded,
      hjust = 0, vjust = 1,
      family = "Helvetica", size = 6
    ) +
    labs(x = xlab, y = ylab, title = if (show_title) title else NULL) +
    xlim(xlim_vals) +
    ylim(ylim_vals) +
    theme_pubr(base_family = "Helvetica", base_size = 16) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5)
    )

  return(p)
}



# helper to strip x or y axis labels/ticks
strip_x <- theme(axis.title.x = element_blank())

strip_y <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)

# --- FULL DATA (row 1) ---
p_swedes_full <- plot_age_comparison(swedes, "Horvath", "age",
  "Horvath DNAm Age", "Chronological Age",
  "Swedish Whites",
  zoomed = FALSE, show_title = TRUE
)

p_genoa_full <- plot_age_comparison(genoa, "Horvath", "age",
  "Horvath DNAm Age", "Chronological Age",
  "GENOA Study African Americans",
  zoomed = FALSE, show_title = TRUE
) +
  strip_y # remove y labels/ticks


p_grady_full <- plot_age_comparison(grady, "Horvath", "age",
  "Horvath DNAm Age", "Chronological Age",
  "Grady Project African Americans",
  zoomed = FALSE, show_title = TRUE
) +
  strip_y # remove y labels/ticks


# --- AGE ≥ 55 SUBSET (row 2) ---
swedes_55 <- subset(swedes, age >= 55)
genoa_55 <- subset(genoa, age >= 55)
grady_55 <- subset(grady, age >= 55)

p_swedes_55 <- plot_age_comparison(swedes_55, "Horvath", "age",
  "Horvath DNAm Age", "Chronological Age",
  title = "Swedish Whites (>= 55 yrs old)", zoomed = TRUE, show_title = T
)

p_genoa_55 <- plot_age_comparison(genoa_55, "Horvath", "age",
  "Horvath DNAm Age", "Chronological Age",
  title = "GENOA Study African Americans (>= 55 yrs old)", zoomed = TRUE, show_title = T
) +
  strip_y

p_grady_55 <- plot_age_comparison(grady_55, "Horvath", "age",
  "Horvath DNAm Age", "Chronological Age",
  title = "Grady Project African Americans (>= 55 yrs old)", zoomed = TRUE, show_title = T
) +
  strip_y

# --- Combine into 2x3 grid ---
final_grid <- (p_swedes_full | p_genoa_full | p_grady_full) /
  (p_swedes_55 | p_genoa_55 | p_grady_55)

final_grid
# saveRDS(final_grid, "plots_RDS/fig3.rds")


####### Relative accuracy analysis
# Get MAGENTA results
magenta_bio_age <- read.csv("data/bio_age_estimates_magenta_age_diff_metadata.csv", row.names = 1)

magenta_bio_age <- magenta_bio_age %>%
  group_by(COHORT) %>%
  filter(STATUS == "CONTROL") %>%
  mutate(COHORT_LABEL = case_when(
    COHORT == "NHW" ~ "White",
    COHORT == "REAAADI" ~ "AA",
    COHORT == "PRADI" ~ "PUR",
    COHORT == "CuADI" ~ "Cuban",
    COHORT == "PERUVIAN" ~ "Peruvian",
    TRUE ~ COHORT
  ))

nhw_controls <- magenta_bio_age %>% filter(COHORT == "NHW" & STATUS == "CONTROL")
aa_controls <- magenta_bio_age %>% filter(COHORT == "REAAADI" & STATUS == "CONTROL")
pr_controls <- magenta_bio_age %>% filter(COHORT == "PRADI" & STATUS == "CONTROL")
cub_controls <- magenta_bio_age %>% filter(COHORT == "CuADI" & STATUS == "CONTROL")
per_controls <- magenta_bio_age %>% filter(COHORT == "PERUVIAN" & STATUS == "CONTROL")

# Bootstrap function
foo <- function(data, indices) {
  dt <- data[indices, ]
  c(
    cor(dt$Horvath, dt$age, method = "p") # Horvath
  )
}

# Now leverage the "boot" function to bootstrap  these correlations 1,000 times.
set.seed(42)
swedes_cor_bootstrap <- boot(swedes_55, foo, R = 1000)
genoa_cor_bootstrap <- boot(genoa_55, foo, R = 1000)
grady_cor_bootstrap <- boot(grady_55, foo, R = 1000)

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
nhw_swedes_compare <- compare_correlations(swedes_cor_bootstrap, nhw_cor_bootstrap, nrow(swedes_55), nrow(nhw_controls))
nhw_genoa_compare <- compare_correlations(genoa_cor_bootstrap, nhw_cor_bootstrap, nrow(genoa_55), nrow(nhw_controls))
nhw_grady_compare <- compare_correlations(grady_cor_bootstrap, nhw_cor_bootstrap, nrow(grady_55), nrow(nhw_controls))

# Put everything into a tidy dataframe for convenience
correlation_df <- magenta_bio_age %>%
  group_by(COHORT) %>%
  filter(STATUS == "CONTROL")

# Add a cohort label column to the replication datasets to be able to merge them here
swedes_55$COHORT_LABEL <- "Swedish Whites"
genoa_55$COHORT_LABEL <- "GENOA AA"
grady_55$COHORT_LABEL <- "Grady AA"

correlation_df <- correlation_df %>%
  bind_rows(swedes_55) %>%
  bind_rows(genoa_55) %>%
  bind_rows(grady_55)


correlation_df <- correlation_df %>%
  mutate(COHORT_LABEL = factor(COHORT_LABEL, levels = c("White", "Cuban", "Peruvian", "PUR", "AA", "GENOA AA", "Grady AA", "Swedish Whites")))


# Step 1: reshape to long
long_df <- correlation_df %>%
  pivot_longer(
    cols = c(Horvath, Hannum, Levine, EN), # add all your models
    names_to = "DNAmAge",
    values_to = "pred_age"
  )

# Step 2: compute correlation within cohort × clock
cor_df <- long_df %>%
  group_by(COHORT_LABEL, DNAmAge) %>%
  summarise(
    correlation = cor(pred_age, age, use = "pairwise.complete.obs"),
    .groups = "drop"
  )

# Step 3: center relative to White baseline
cor_df <- cor_df %>%
  group_by(DNAmAge) %>%
  mutate(relative_cor = correlation - correlation[COHORT_LABEL == "White"]) %>%
  ungroup()



####### Plot relative correlations
cor_df <- cor_df %>%
  mutate(DNAmAge = recode(DNAmAge,
    "EN" = "Zhang2019_EN",
    "Levine" = "PhenoAge"
  )) %>%
  mutate(DNAmAge = factor(DNAmAge,
    levels = c("Horvath", "Hannum", "Zhang2019_EN", "PhenoAge")
  ))


# ggplot(cor_df, aes(x = COHORT_LABEL, y = relative_cor, fill = COHORT_LABEL)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_text(aes(label = ifelse(COHORT_LABEL == "White", round(correlation, 2), "")),
#             vjust = -0.5, size = 5, color = "black") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   facet_wrap(~ DNAmAge, scales = "free_x") +
#   xlab("Cohort") +
#   ylab("Pearson Correlation (Relative to White)") +
#   coord_cartesian(ylim = c(-0.3, 0.3)) +  # Adjust the y-axis limits as needed
#   theme_classic() +
#   theme(
#     axis.title = element_text(size = 16),
#     axis.text = element_text(size = 14),
#     plot.title = element_text(size = 18),
#     legend.position = "none",
#     strip.text = element_text(size = 18),
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)  # Flip x ticks
#   )
library(ggpubr)
p <- ggbarplot(
  cor_df,
  x = "COHORT_LABEL",
  y = "relative_cor",
  fill = "COHORT_LABEL",
  color = "black",
  add = "none",
  position = position_dodge()
) +
  geom_text(
    data = subset(cor_df, COHORT_LABEL == "White"),
    aes(label = round(correlation, 2), y = relative_cor),
    vjust = -0.5, size = 5, color = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~DNAmAge, scales = "free_x") +
  xlab("Cohort") +
  ylab("Pearson Correlation (Relative to White)") +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.position = "none",
    strip.text = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

p
custom_colors <- c(
  "White" = "black",
  "Swedish Whites" = "black",
  "AA" = "#E41A1C", # vivid red
  "Grady AA" = "#E41A1C", # same family
  "GENOA AA" = "#E41A1C",
  "PUR" = "#FF7F00", # bright orange
  "Cuban" = "grey60",
  "Peruvian" = "grey60"
)

p + scale_fill_manual(values = custom_colors)




# # 1) Relabel, set factor order for cohorts, AND set factor order for facets
# cor_df <- cor_df %>%
#   mutate(
#     COHORT_LABEL = recode(
#       COHORT_LABEL,
#       "White" = "White\n(MAGENTA)",
#       "Cuban" = "Cuban\n(MAGENTA)",
#       "Peruvian" = "Peruvian\n(MAGENTA)",
#       "PUR" = "Puerto Rican\n(MAGENTA)",
#       "AA" = "African American\n(MAGENTA)",
#       "GENOA AA" = "African American\n(GENOA)",
#       "Grady AA" = "African American\n(Grady)",
#       "Swedish Whites" = "White\n(Swedish)"
#     ),
#     COHORT_LABEL = factor(
#       COHORT_LABEL,
#       levels = c(
#         "White\n(MAGENTA)", "White\n(Swedish)",
#         "Cuban\n(MAGENTA)", "Peruvian\n(MAGENTA)",
#         "Puerto Rican\n(MAGENTA)",
#         "African American\n(MAGENTA)",
#         "African American\n(GENOA)",
#         "African American\n(Grady)"
#       )
#     ),
#     DNAmAge = factor(DNAmAge, levels = c("Horvath", "Hannum", "Zhang2019_EN", "PhenoAge"))
#   )
#
# # 2) Color palette
# custom_colors <- c(
#   "White\n(MAGENTA)" = "grey50",
#   "White\n(Swedish)" = "grey50",
#   "Cuban\n(MAGENTA)" = "grey70",
#   "Peruvian\n(MAGENTA)" = "grey70",
#   "Puerto Rican\n(MAGENTA)" = "#FDBE85",
#   "African American\n(MAGENTA)" = "#E6550D",
#   "African American\n(GENOA)" = "#D94801",
#   "African American\n(Grady)" = "#A63603"
# )
#
# # 3) Base plot
# p <- ggbarplot(
#   cor_df,
#   x = "COHORT_LABEL",
#   y = "relative_cor",
#   fill = "COHORT_LABEL",
#   color = "black",
#   add = "none",
#   position = position_dodge()
# ) +
#   geom_text(
#     data = subset(cor_df, COHORT_LABEL == "White\n(MAGENTA)"),
#     aes(label = round(correlation, 2), y = relative_cor),
#     vjust = -0.5, size = 6, color = "black"
#   ) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   facet_wrap(~ DNAmAge, scales = "free_x") +
#   xlab("") +
#   ylab("Pearson Correlation (Relative to White)") +
#   coord_cartesian(ylim = c(-0.3, 0.3)) +
#   scale_fill_manual(values = custom_colors) +
#   theme_classic() +
#   theme(
#     axis.title = element_text(size = 16),
#     axis.text = element_text(size = 14),
#     plot.title = element_text(size = 18),
#     legend.position = "none",
#     strip.text = element_text(size = 18),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
#
# # 4) Build dataframe with asterisk positions
# sig_map <- list(
#   "Zhang2019_EN" = c("African American\n(MAGENTA)", "African American\n(Grady)"),
#   "Hannum"       = c("African American\n(Grady)"),
#   "Horvath"      = c("Puerto Rican\n(MAGENTA)", "African American\n(MAGENTA)", "African American\n(Grady)"),
#   "PhenoAge"     = c()
# )
#
# # --- Build dataframe with single asterisks slightly above bars ---
# asterisk_df <- purrr::map_df(names(sig_map), function(clock) {
#   targets <- sig_map[[clock]]
#   if (length(targets) == 0) return(NULL)
#
#   purrr::map_df(targets, function(grp) {
#     bar_height <- max(cor_df$relative_cor[cor_df$DNAmAge == clock & cor_df$COHORT_LABEL == grp], na.rm = TRUE)
#     data.frame(
#       DNAmAge = clock,
#       COHORT_LABEL = grp,
#       y = 0.05,   # leave extra room
#       label = "*"
#     )
#   })
# }) %>%
#   mutate(DNAmAge = factor(DNAmAge, levels = c("Horvath", "Hannum", "Zhang2019_EN", "PhenoAge")))
#
# # --- Add them to the plot ---
# p_final <- p +
#   geom_text(
#     data = asterisk_df,
#     aes(x = COHORT_LABEL, y = y, label = label),
#     inherit.aes = FALSE,
#     size = 8,
#     vjust = 0
#   )
#
# print(p_final)

# Example: single column figure
# ggsave(
#   "plots_RDS/figure4_clock_performances.pdf",
#   plot = p_final,
#   width = 264/25.4,   # convert mm to inches
#   height = 150/25.4, # pick a good height
#   units = "in"
# )

# # 1) Keep your original cohort recoding, factor orders intact
# cor_df <- cor_df %>%
#   mutate(
#     COHORT_LABEL = recode(
#       COHORT_LABEL,
#       "White" = "White\n(MAGENTA)",
#       "Cuban" = "Cuban\n(MAGENTA)",
#       "Peruvian" = "Peruvian\n(MAGENTA)",
#       "PUR" = "Puerto Rican\n(MAGENTA)",
#       "AA" = "African American\n(MAGENTA)",
#       "GENOA AA" = "African American\n(GENOA)",
#       "Grady AA" = "African American\n(Grady)",
#       "Swedish Whites" = "White\n(Swedish)"
#     ),
#     COHORT_LABEL = factor(
#       COHORT_LABEL,
#       levels = c(
#         "White\n(MAGENTA)", "White\n(Swedish)",
#         "Cuban\n(MAGENTA)", "Peruvian\n(MAGENTA)",
#         "Puerto Rican\n(MAGENTA)",
#         "African American\n(MAGENTA)",
#         "African American\n(GENOA)",
#         "African American\n(Grady)"
#       )
#     ),
#     DNAmAge = factor(DNAmAge, levels = c("Horvath", "Hannum", "Zhang2019_EN", "PhenoAge"))
#   )
#
# # 2) Color palette stays with the full labels
# custom_colors <- c(
#   "White\n(MAGENTA)" = "grey50",
#   "White\n(Swedish)" = "grey50",
#   "Cuban\n(MAGENTA)" = "grey70",
#   "Peruvian\n(MAGENTA)" = "grey70",
#   "Puerto Rican\n(MAGENTA)" = "#FDBE85",
#   "African American\n(MAGENTA)" = "#E6550D",
#   "African American\n(GENOA)" = "#D94801",
#   "African American\n(Grady)" = "#A63603"
# )
#
# # 3) Function for relabeling x-axis text
# label_fun <- function(x) {
#   # Remove only the MAGENTA part (and its preceding newline if present)
#   gsub("\n?\\(MAGENTA\\)", "", x)
# }
#
# # 4) Base plot
# p <- ggbarplot(
#   cor_df,
#   x = "COHORT_LABEL",
#   y = "relative_cor",
#   fill = "COHORT_LABEL",
#   color = "black",
#   add = "none",
#   position = position_dodge()
# ) +
#   geom_text(
#     data = subset(cor_df, COHORT_LABEL == "White\n(MAGENTA)"),
#     aes(label = round(correlation, 2), y = relative_cor),
#     vjust = -0.5, size = 6, color = "black"
#   ) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   facet_wrap(~ DNAmAge, scales = "free_x") +
#   xlab("") +
#   ylab("Pearson Correlation (Relative to White)") +
#   coord_cartesian(ylim = c(-0.3, 0.3)) +
#   scale_fill_manual(values = custom_colors) +
#   scale_x_discrete(labels = label_fun) +   # <<< relabeling happens here
#   theme_classic() +
#   theme(
#     axis.title = element_text(size = 16),
#     axis.text = element_text(size = 14),
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # reduced tick size
#     plot.title = element_text(size = 18),
#     legend.position = "none",
#     strip.text = element_text(size = 18)
#   )
#
# # 5) Build dataframe with asterisk positions (no changes needed)
# sig_map <- list(
#   "Zhang2019_EN" = c("African American\n(MAGENTA)", "African American\n(Grady)"),
#   "Hannum"       = c("African American\n(Grady)"),
#   "Horvath"      = c("Puerto Rican\n(MAGENTA)", "African American\n(MAGENTA)", "African American\n(Grady)"),
#   "PhenoAge"     = c()
# )
#
# asterisk_df <- purrr::map_df(names(sig_map), function(clock) {
#   targets <- sig_map[[clock]]
#   if (length(targets) == 0) return(NULL)
#
#   purrr::map_df(targets, function(grp) {
#     bar_height <- max(cor_df$relative_cor[cor_df$DNAmAge == clock & cor_df$COHORT_LABEL == grp], na.rm = TRUE)
#     data.frame(
#       DNAmAge = clock,
#       COHORT_LABEL = grp,
#       y = 0.05,
#       label = "*"
#     )
#   })
# }) %>%
#   mutate(DNAmAge = factor(DNAmAge, levels = c("Horvath", "Hannum", "Zhang2019_EN", "PhenoAge")))
#
# # 6) Add them to the plot
# p_final <- p +
#   geom_text(
#     data = asterisk_df,
#     aes(x = COHORT_LABEL, y = y, label = label),
#     inherit.aes = FALSE,
#     size = 8,
#     vjust = 0
#   )
#
# print(p_final)


cor_df <- cor_df %>%
  mutate(DNAmAge = recode(DNAmAge,
    "EN" = "Zhang2019_EN",
    "Levine" = "PhenoAge"
  )) %>%
  filter(DNAmAge != "PhenoAge") %>% # Remove PhenoAge
  mutate(DNAmAge = factor(DNAmAge,
    levels = c("Horvath", "Hannum", "Zhang2019_EN")
  ))

# 1) Cohort recoding and updated factor order (MAGENTA first, externals last)
cor_df <- cor_df %>%
  mutate(
    COHORT_LABEL = recode(
      COHORT_LABEL,
      "White" = "White (MAGENTA)",
      "Cuban" = "Cuban (MAGENTA)",
      "Peruvian" = "Peruvian (MAGENTA)",
      "PUR" = "Puerto Rican (MAGENTA)",
      "AA" = "African American (MAGENTA)",
      "GENOA AA" = "African American (GENOA)",
      "Grady AA" = "African American (Grady)",
      "Swedish Whites" = "White (Swedish)"
    )
  )

# Create a combined factor to allow reordering within facets
cor_df <- cor_df %>%
  mutate(
    COHORT_FACET = paste(COHORT_LABEL, DNAmAge, sep = "___"),
    # Create a sorting value: White (MAGENTA) gets a high value to be at the top
    sort_val = ifelse(COHORT_LABEL == "White (MAGENTA)", 1000, relative_cor),
    COHORT_FACET = reorder(COHORT_FACET, sort_val)
  )

# 2) Color palette (unchanged)
custom_colors <- c(
  "White (MAGENTA)" = "grey50",
  "White (Swedish)" = "grey50",
  "Cuban (MAGENTA)" = "grey70",
  "Peruvian (MAGENTA)" = "grey70",
  "Puerto Rican (MAGENTA)" = "#FDBE85",
  "African American (MAGENTA)" = "#E6550D",
  "African American (GENOA)" = "#D94801",
  "African American (Grady)" = "#A63603"
)

# 3) Function for relabeling axis text
label_fun <- function(x) {
  # Remove the facet suffix and the (MAGENTA) tag
  x <- gsub("___.*", "", x)
  gsub("\\s*\\(MAGENTA\\)", "", x)
}

# 4) Asterisk dataframe
sig_map <- list(
  "Zhang2019_EN" = c("African American (MAGENTA)", "African American (Grady)"),
  "Hannum"       = c("African American (Grady)"),
  "Horvath"      = c("Puerto Rican (MAGENTA)", "African American (MAGENTA)", "African American (Grady)")
  # PhenoAge removed
)

asterisk_df <- purrr::map_df(names(sig_map), function(clock) {
  targets <- sig_map[[clock]]
  if (length(targets) == 0) {
    return(NULL)
  }

  purrr::map_df(targets, function(grp) {
    # Get the correlation height for this bar
    bar_height <- cor_df$relative_cor[
      cor_df$DNAmAge == clock & cor_df$COHORT_LABEL == grp
    ]

    # Add small offset beyond the bar end
    offset <- ifelse(bar_height >= 0, 0.02, -0.02)

    data.frame(
      DNAmAge = clock,
      COHORT_LABEL = grp,
      y = bar_height + offset,
      label = "*"
    )
  })
}) %>%
  mutate(
    DNAmAge = factor(DNAmAge, levels = c("Horvath", "Hannum", "Zhang2019_EN")),
    COHORT_FACET = paste(COHORT_LABEL, DNAmAge, sep = "___")
  )


# 5) Base plot (flipped orientation, single-line labels)
p <- ggbarplot(
  cor_df,
  x = "COHORT_FACET", # Use the combined factor
  y = "relative_cor",
  fill = "COHORT_LABEL",
  color = "black",
  add = "none",
  position = position_dodge()
) +
  geom_text(
    data = subset(cor_df, COHORT_LABEL == "White (MAGENTA)"),
    aes(label = round(correlation, 2), y = relative_cor),
    hjust = -0.3, size = 6, color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~DNAmAge, ncol = 1, scales = "free") + # Free scales for both axes
  ylab("Pearson Correlation (Relative to White)") +
  xlab("Cohort") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(labels = label_fun) +
  coord_flip() + # Removed fixed ylim
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, lineheight = 0.8),
    plot.title = element_text(size = 18),
    legend.position = "none",
    strip.text = element_text(size = 18),
    panel.spacing.y = unit(0, "cm")
  )


p_final <- p +
  geom_text(
    data = asterisk_df,
    aes(x = COHORT_FACET, y = y, label = label), # Use COHORT_FACET here too
    inherit.aes = FALSE,
    size = 8,
    hjust = -0.2
  )

print(p_final)

# --- Absolute Correlation Plot ---
# Reuse the sorted factor from before.
# Note: Sorting by relative_cor is equivalent to sorting by correlation within each facet,
# and White (MAGENTA) is forced to top by the sort_val logic.

# Recalculate asterisk positions for absolute values
asterisk_df_abs <- purrr::map_df(names(sig_map), function(clock) {
  targets <- sig_map[[clock]]
  if (length(targets) == 0) {
    return(NULL)
  }

  purrr::map_df(targets, function(grp) {
    # Get the correlation height for this bar
    bar_height <- cor_df$correlation[
      cor_df$DNAmAge == clock & cor_df$COHORT_LABEL == grp
    ]

    # Add small offset beyond the bar end
    offset <- ifelse(bar_height >= 0, 0.02, -0.02)

    data.frame(
      DNAmAge = clock,
      COHORT_LABEL = grp,
      y = bar_height + offset,
      label = "*"
    )
  })
}) %>%
  mutate(
    DNAmAge = factor(DNAmAge, levels = c("Horvath", "Hannum", "Zhang2019_EN")),
    COHORT_FACET = paste(COHORT_LABEL, DNAmAge, sep = "___")
  )

# Extract White (MAGENTA) correlations for the reference line
white_cor_df <- cor_df %>%
  filter(COHORT_LABEL == "White (MAGENTA)")

p_abs <- ggbarplot(
  cor_df,
  x = "COHORT_FACET",
  y = "correlation",
  fill = "COHORT_LABEL",
  color = "black",
  add = "none",
  position = position_dodge()
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # Add reference line for White (MAGENTA) correlation
  geom_hline(
    data = white_cor_df,
    aes(yintercept = correlation),
    linetype = "dashed",
    color = "darkgrey",
    size = 1
  ) +
  facet_wrap(~DNAmAge, ncol = 1, scales = "free") +
  ylab("Pearson Correlation") +
  xlab("Cohort") +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(labels = label_fun) +
  coord_flip(ylim = c(0, 1)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14, lineheight = 0.8),
    plot.title = element_text(size = 18),
    legend.position = "none",
    strip.text = element_text(size = 18),
    panel.spacing.y = unit(0, "cm")
  ) +
  geom_text(
    data = asterisk_df_abs,
    aes(x = COHORT_FACET, y = y, label = label),
    inherit.aes = FALSE,
    size = 8,
    hjust = -0.2
  )

print(p_abs)
