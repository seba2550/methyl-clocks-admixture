# Establish working directory (uncomment/adjust if not running from repository root)
# setwd("/Users/sgonzalez/Desktop/Capra Lab/Thesis Project/Aim_1")

# Load our libraries as needed
library(methylclock)
library(tidyverse)
library(broom)
library(readxl)
library(UpSetR)
library(ComplexUpset)
library(UpSetR)
library(ComplexUpset)
library(eulerr) # For proportional Venn/Euler diagrams





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
    TRUE ~ COHORT # Keep the original value if no match is found
  ))

# Modify the rownames to remove everything after the underscore
rownames(normalized_combined) <- gsub("_.*", "", rownames(normalized_combined))

# # Get biological age estimates by using the DNAm Age clocks
# bio_age_estimates_magenta_age_diff <- DNAmAge(normalized_combined, cell.count = F, age = sample_metadata$AGE_OF_EXAM, normalize = F)
#
# # Join the metadata to the age estimates, using sample IDs as keys
# bio_age_estimates_magenta_age_diff_metadata <- bio_age_estimates_magenta_age_diff %>% left_join(sample_metadata, join_by("id" == "Beta_ID"))

# Read in the results of the DNAmAge calculations along with the corresponding metadata
bio_age_estimates_magenta_age_diff_metadata <- read.csv("bio_age_estimates_magenta_age_diff_metadata.csv", row.names = 1)
load_DNAm_Clocks_data()

clock_cpgs <- list(
  Horvath = coefHorvath,
  Hannum = coefHannum,
  EN = coefEN,
  Levine = coefLevine
)

#### One clock to try this out first
# Choose which clock to analyze
clock_name <- "Horvath"

# Get CpGs for that clock
cpgs <- clock_cpgs[[clock_name]]
cpgs <- intersect(cpgs$CpGmarker, rownames(normalized_combined))

# Subset methylation matrix
meth_mat <- normalized_combined[cpgs, , drop = FALSE]

# Make sure sample order matches between dataframes
common_ids <- intersect(
  colnames(meth_mat),
  bio_age_estimates_magenta_age_diff_metadata$id
)

meth_mat <- meth_mat[, common_ids]
meta <- bio_age_estimates_magenta_age_diff_metadata %>%
  filter(id %in% common_ids)

# Loop through CpGs and regress ClockError ~ Methylation
results <- map_dfr(cpgs, function(cpg) {
  df <- tibble(
    ClockError = meta[[paste0("ageAcc.", clock_name)]],
    Methylation = as.numeric(meth_mat[cpg, ])
  )

  tidy(lm(ClockError ~ Methylation, data = df)) %>%
    mutate(CpG = cpg, Clock = clock_name)
})

results_clean_horvath <- results %>%
  filter(term == "Methylation") %>%
  select(CpG, estimate, std.error, statistic, p.value, Clock)

##### Hannum
clock_name <- "Hannum"

cpgs <- clock_cpgs[[clock_name]]
cpgs <- intersect(cpgs$CpGmarker, rownames(normalized_combined))

meth_mat <- normalized_combined[cpgs, , drop = FALSE]

common_ids <- intersect(
  colnames(meth_mat),
  bio_age_estimates_magenta_age_diff_metadata$id
)

meth_mat <- meth_mat[, common_ids]
meta <- bio_age_estimates_magenta_age_diff_metadata %>%
  filter(id %in% common_ids)

results <- map_dfr(cpgs, function(cpg) {
  df <- tibble(
    ClockError = meta[[paste0("ageAcc.", clock_name)]],
    Methylation = as.numeric(meth_mat[cpg, ])
  )

  tidy(lm(ClockError ~ Methylation, data = df)) %>%
    mutate(CpG = cpg, Clock = clock_name)
})

results_clean_hannum <- results %>%
  filter(term == "Methylation") %>%
  select(CpG, estimate, std.error, statistic, p.value, Clock)

###### EN
clock_name <- "EN"

cpgs <- clock_cpgs[[clock_name]]
cpgs <- intersect(cpgs$CpGmarker, rownames(normalized_combined))

meth_mat <- normalized_combined[cpgs, , drop = FALSE]

common_ids <- intersect(
  colnames(meth_mat),
  bio_age_estimates_magenta_age_diff_metadata$id
)

meth_mat <- meth_mat[, common_ids]
meta <- bio_age_estimates_magenta_age_diff_metadata %>%
  filter(id %in% common_ids)

results <- map_dfr(cpgs, function(cpg) {
  df <- tibble(
    ClockError = meta[[paste0("ageAcc.", clock_name)]],
    Methylation = as.numeric(meth_mat[cpg, ])
  )

  tidy(lm(ClockError ~ Methylation, data = df)) %>%
    mutate(CpG = cpg, Clock = clock_name)
})

results_clean_en <- results %>%
  filter(term == "Methylation") %>%
  select(CpG, estimate, std.error, statistic, p.value, Clock)

###### PhenoAge
clock_name <- "Levine"

cpgs <- clock_cpgs[[clock_name]]
cpgs <- intersect(cpgs$CpGmarker, rownames(normalized_combined))

meth_mat <- normalized_combined[cpgs, , drop = FALSE]

common_ids <- intersect(
  colnames(meth_mat),
  bio_age_estimates_magenta_age_diff_metadata$id
)

meth_mat <- meth_mat[, common_ids]
meta <- bio_age_estimates_magenta_age_diff_metadata %>%
  filter(id %in% common_ids)

results <- map_dfr(cpgs, function(cpg) {
  df <- tibble(
    ClockError = meta[[paste0("ageAcc.", clock_name)]],
    Methylation = as.numeric(meth_mat[cpg, ])
  )

  tidy(lm(ClockError ~ Methylation, data = df)) %>%
    mutate(CpG = cpg, Clock = clock_name)
})

results_clean_phenoage <- results %>%
  filter(term == "Methylation") %>%
  select(CpG, estimate, std.error, statistic, p.value, Clock)

####### Combine and adjust p-values
results_all <- bind_rows(
  results_clean_horvath,
  results_clean_hannum,
  results_clean_en,
  results_clean_phenoage
)

# Ensure the Clock order is fixed
results_all <- results_all %>%
  mutate(Clock = factor(Clock, levels = c("Horvath", "Hannum", "EN", "Levine")))

# OPTION 1: Adjust within each clock separately (recommended)
results_all <- results_all %>%
  group_by(Clock) %>%
  mutate(
    p.adjusted = p.adjust(p.value, method = "bonferroni"),
    p.fdr = p.adjust(p.value, method = "fdr")
  ) %>%
  ungroup()

# OPTION 2: Adjust across all clocks together (more conservative)
# Uncomment if you want to correct across all tests globally
# results_all <- results_all %>%
#   mutate(
#     p.adjusted = p.adjust(p.value, method = "bonferroni"),
#     p.fdr = p.adjust(p.value, method = "fdr")
#   )

####### Visualization
# Volcano plot
ggplot(results_all, aes(x = estimate, y = -log10(p.adjusted), color = Clock)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  facet_wrap(~Clock, scales = "fixed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Effect size (estimate)",
    y = expression(-log[10](adjusted ~ p ~ value)),
    title = "Volcano plots of CpG associations by clock (Bonferroni-adjusted)"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14)
  )
ggplot(results_all, aes(x = estimate, y = -log10(p.fdr), color = Clock)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  facet_wrap(~Clock, scales = "fixed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Effect size (estimate)",
    y = expression(-log[10](adjusted ~ p ~ value)),
    title = "Volcano plots of CpG associations by clock (FDR-adjusted)"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14)
  )

# Manhattan-style plot
results_all %>%
  group_by(Clock) %>%
  mutate(index = row_number()) %>%
  ggplot(aes(x = index, y = -log10(p.adjusted), color = Clock)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.7) +
  facet_wrap(~Clock, scales = "free_x") +
  theme_classic(base_size = 14) +
  labs(
    x = "Clock CpGs (index)",
    y = expression(-log[10](adjusted ~ p ~ value)),
    title = "CpG significance distribution across clocks (Bonferroni-adjusted)"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14)
  )
results_all %>%
  filter(Clock != "Levine") %>%
  group_by(Clock) %>%
  arrange(p.fdr) %>%
  mutate(index = row_number()) %>%
  ggplot(aes(x = index, y = -log10(p.fdr), color = Clock)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.7) +
  facet_wrap(~Clock, scales = "free_x") +
  coord_cartesian(ylim = c(0, 8)) +
  theme_classic(base_size = 14) +
  labs(
    x = "Clock CpGs (index)",
    y = expression(-log[10](adjusted ~ p ~ value))
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14)
  )


# Top hits (now using adjusted p-values)
top_hits <- results_all %>%
  group_by(Clock) %>%
  arrange(p.adjusted) %>%
  slice_head(n = 20)

ggplot(top_hits, aes(x = reorder(CpG, -log10(p.adjusted)), y = -log10(p.adjusted), fill = Clock)) +
  geom_col() +
  facet_wrap(~Clock, scales = "free_x") +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    x = "Top CpGs",
    y = expression(-log[10](adjusted ~ p ~ value)),
    title = "Top 20 most significant CpGs per clock (Bonferroni-adjusted)"
  ) +
  theme(axis.text.y = element_text(size = 6))

top_hits <- results_all %>%
  group_by(Clock) %>%
  arrange(p.fdr) %>%
  slice_head(n = 20)

ggplot(top_hits, aes(x = reorder(CpG, -log10(p.fdr)), y = -log10(p.fdr), fill = Clock)) +
  geom_col() +
  facet_wrap(~Clock, scales = "free_x") +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    x = "Top CpGs",
    y = expression(-log[10](adjusted ~ p ~ value)),
    title = "Top 20 most significant CpGs per clock (FDR-adjusted)"
  ) +
  theme(axis.text.y = element_text(size = 6))

# Get significant CpGs (adjusted p < 0.05)
sig_cpgs <- results_all %>%
  filter(p.fdr < 0.05) %>%
  select(CpG, Clock, p.value, p.fdr, estimate)

# Summary of significant hits per clock
sig_cpgs %>%
  group_by(Clock) %>%
  summarise(
    n_significant = n(),
    .groups = "drop"
  )

# List of significant CpGs per clock for downstream analysis
cpg_list <- sig_cpgs %>%
  group_by(Clock) %>%
  summarise(CpGs = list(unique(CpG)), .groups = "drop") %>%
  deframe()

clock_names <- names(cpg_list)

# Save error-associated CpGs for use in other scripts
# saveRDS(sig_cpgs, "error_associated_cpgs.rds")

# --- 1) Pairwise overlaps (all pairs) ---
pairwise_df <- as.data.frame(t(combn(clock_names, 2)), stringsAsFactors = FALSE) %>%
  setNames(c("Clock1", "Clock2")) %>%
  rowwise() %>%
  mutate(
    Shared_CpGs = list(intersect(cpg_list[[Clock1]], cpg_list[[Clock2]])),
    Num_Shared = length(Shared_CpGs)
  ) %>%
  ungroup()

# view counts
pairwise_df %>% select(Clock1, Clock2, Num_Shared)

# to inspect actual CpGs for the first row:
pairwise_df$Shared_CpGs[[1]]

# --- 2) Multi-way overlaps: shared across 3, 4, ... k clocks ---
multi_overlap_list <- lapply(2:length(clock_names), function(k) {
  cmbs <- combn(clock_names, k, simplify = FALSE)
  purrr::map_df(cmbs, function(cm) {
    shared <- Reduce(intersect, cpg_list[cm])
    data.frame(
      Clocks = paste(cm, collapse = "|"),
      K = k,
      Num_Shared = length(shared),
      Shared_CpGs = I(list(shared)),
      stringsAsFactors = FALSE
    )
  })
}) %>% bind_rows()

# view multi-way overlaps
multi_overlap_list %>% arrange(desc(Num_Shared))

# --- 3) (Optional) UpSet plot to visualise overlaps ---
# Quick diagnostics function
diagnose_upset_df <- function(df) {
  cat("Column classes:\n")
  print(sapply(df, class))
  cat("\nAny NA?:", any(is.na(df)), "\n")
  cat("Number of rows:", nrow(df), " Number of cols:", ncol(df), "\n")
  cat("Column sums (should be integers >=0):\n")
  print(colSums(df))
}
# Build presence/absence df (NOT converting CpG to rownames — ComplexUpset uses columns)
pa_df <- sig_cpgs %>%
  distinct(CpG, Clock) %>%
  mutate(present = 1L) %>%
  pivot_wider(
    names_from = Clock,
    values_from = present,
    values_fill = 0
  )

# Diagnostics
diagnose_upset_df(pa_df %>% select(-CpG)) # exclude CpG for column diagnostics

# Build the upset plot. Provide the set column names explicitly:
set_cols <- setdiff(names(pa_df), "CpG")

# Example: show top 10 intersections by size
ComplexUpset::upset(
  pa_df,
  intersect = set_cols,
  themes = list(default = theme_classic()),
  base_annotations = list(
    "Intersection size" = intersection_size(counts = TRUE)
  ),
  min_size = 1
) + theme(
  panel.background = element_blank(),
  axis.text.x = element_blank(), # removes x-axis text (group labels)
  axis.ticks.x = element_blank(), # removes x-axis ticks
  strip.text = element_blank(), # removes 'group' labels
  plot.title = element_text(hjust = 0.5) # centers title
) +
  ggtitle("CpG Overlaps Across Clocks")



# Count how many clocks share each CpG
shared_cpgs <- sig_cpgs %>%
  group_by(CpG) %>%
  summarise(n_clocks = n_distinct(Clock)) %>%
  count(n_clocks)

print(shared_cpgs)


##### Connecting the dots with the AFR differentially methylated CpGs (defined by Oge's analysis)
# Load the AFR vs EUR differentially methylated sites defined by Oge
load("~/Desktop/Capra Lab/Thesis Project/Aim_1/ancestry_CpGs.RData")

# Apply significance filter
afr_vs_eur_filtered <- subset(AFR_vs_EUR, adj.P.Val < 0.05)

# Final AFR vs EUR CpG list
afr_dm_cpgs <- rownames(afr_vs_eur_filtered)

sig_cpgs_afr_dm <- sig_cpgs[sig_cpgs$CpG %in% afr_dm_cpgs, ]
sig_cpgs_afr_dm[duplicated(sig_cpgs_afr_dm$CpG), ]


######## meQTL stuff

# Load the MAGENTA results
all_results <- readRDS("all_results.rds")

# Pass individual cohort data
aa_meqtl <- all_results$results$AFR$meqtl_data
pr_meqtl <- all_results$results$PR$meqtl_data
cub_meqtl <- all_results$results$CUB$meqtl_data
per_meqtl <- all_results$results$PER$meqtl_data

# Remove all_results object from memory to alleviate RAM
rm(all_results)


##### Now on a per-clock basis
# --- Step 1: Harmonize CpG column names for meQTL sets ---
aa_cpgs <- aa_meqtl %>% pull(CpG)
pr_cpgs <- pr_meqtl %>% pull(PROBE)
cub_cpgs <- cub_meqtl %>% pull(PROBE)
per_cpgs <- per_meqtl %>% pull(PROBE)

# --- Step 2: Collect all CpGs across all clocks ---
clocks <- unique(sig_cpgs$Clock)

# Get all significant CpGs across all clocks
all_sig_cpgs <- sig_cpgs %>% pull(CpG)
all_sig_afr_cpgs <- sig_cpgs_afr_dm %>% pull(CpG)

# Combine all unique CpGs
all_cpgs <- unique(c(
  aa_cpgs, pr_cpgs, cub_cpgs, per_cpgs,
  all_sig_cpgs, all_sig_afr_cpgs
))

# --- Step 3: Build combined upset dataframe ---
upset_df <- data.frame(
  CpG = all_cpgs,
  AA_meQTL = all_cpgs %in% aa_cpgs,
  PR_meQTL = all_cpgs %in% pr_cpgs,
  CUB_meQTL = all_cpgs %in% cub_cpgs,
  PER_meQTL = all_cpgs %in% per_cpgs,
  Sig_CpGs = all_cpgs %in% all_sig_cpgs,
  Sig_CpGs_AFR = all_cpgs %in% all_sig_afr_cpgs
)

# Add clock membership columns
for (clk in clocks) {
  clock_cpgs <- sig_cpgs %>%
    filter(Clock == clk) %>%
    pull(CpG)

  upset_df[[clk]] <- upset_df$CpG %in% clock_cpgs
}

# --- Step 4: Create combined UpSet plot ---
# Define which sets to include in the intersection
set_columns <- c(
  "Sig_CpGs_AFR", "Sig_CpGs", "AA_meQTL", "PR_meQTL",
  "CUB_meQTL", "PER_meQTL", "Horvath", "Hannum", "EN", "Levine"
)

upset(
  upset_df,
  intersect = set_columns,
  base_annotations = list(
    "Intersection size" = intersection_size(counts = TRUE)
  ),
  themes = upset_modify_themes(
    list("intersections_matrix" = theme_classic())
  )
) +
  ggtitle("Combined UpSet Plot: All Clocks")



#########################################
# --- Step 1: Harmonize CpG column names for meQTL sets ---
aa_cpgs <- aa_meqtl %>% pull(CpG)
pr_cpgs <- pr_meqtl %>% pull(PROBE)
cub_cpgs <- cub_meqtl %>% pull(PROBE)
per_cpgs <- per_meqtl %>% pull(PROBE)

# Combine all meQTL CpGs (any ancestry)
all_meqtl_cpgs <- unique(c(aa_cpgs, pr_cpgs, cub_cpgs, per_cpgs))

# Identify African-differentiated meQTL CpGs
# Assuming your data has an 'african_differentiated' column
afr_diff_meqtl_cpgs <- unique(c(
  aa_meqtl %>% filter(african_differentiated == TRUE) %>% pull(CpG),
  pr_meqtl %>% filter(african_differentiated == TRUE) %>% pull(PROBE),
  cub_meqtl %>% filter(african_differentiated == TRUE) %>% pull(PROBE),
  per_meqtl %>% filter(african_differentiated == TRUE) %>% pull(PROBE)
))

# --- Step 2: Create upset plot for each clock ---
clocks <- unique(sig_cpgs$Clock)

# Relabel "Levine" to "PhenoAge"
clocks <- gsub("Levine", "PhenoAge", clocks)

# Exclude PhenoAge
clocks <- clocks[clocks != "PhenoAge"]

sig_cpgs <- sig_cpgs %>%
  mutate(Clock = gsub("Levine", "PhenoAge", Clock))

# Store plots in a list
venn_plots <- list()

# --- Pre-calculate maximum scale for consistent sizing ---
# We need to find the clock with the largest union of sets to determine the scale
max_total_count <- 0

for (clk in clocks) {
  # Get CpGs significant for this specific clock
  clock_sig_cpgs <- sig_cpgs %>%
    filter(Clock == clk) %>%
    pull(CpG)

  # Determine the clock definition object
  if (clk == "Horvath") {
    clock_def_cpgs <- coefHorvath$CpGmarker
  } else if (clk == "Hannum") {
    clock_def_cpgs <- coefHannum$CpGmarker
  } else if (clk == "EN") {
    clock_def_cpgs <- coefEN$CpGmarker
  } else if (clk == "PhenoAge") {
    clock_def_cpgs <- coefLevine$CpGmarker
  }

  # Remove intercept if present
  clock_def_cpgs <- clock_def_cpgs[clock_def_cpgs != "(Intercept)"]

  # Filter meQTL lists to this clock's universe
  clock_meqtl_cpgs <- intersect(all_meqtl_cpgs, clock_def_cpgs)

  # Union size (2-way)
  union_size <- length(union(clock_sig_cpgs, clock_meqtl_cpgs))

  # Also consider 4-way unions if this clock is in the 4-way list
  if (clk != "PhenoAge") {
    # For 4-way, we add AFR meQTL and AFR vs EUR to the union consideration
    # We just need the max possible extent, so a union of all 4 sets is safe
    set_afr_meqtl <- intersect(afr_diff_meqtl_cpgs, clock_def_cpgs)
    set_afr_dm <- intersect(afr_dm_cpgs, clock_def_cpgs)

    union_size_4way <- length(unique(c(clock_sig_cpgs, clock_meqtl_cpgs, set_afr_meqtl, set_afr_dm)))
    if (union_size_4way > union_size) {
      union_size <- union_size_4way
    }
  }

  if (union_size > max_total_count) {
    max_total_count <- union_size
  }
}

# Calculate anchor "universe" size
# We want every plot to represent the same total universe size to ensure relative scaling.
# We will effectively "pad" each plot with an invisible anchor set so that (Real Data + Anchor) = max_total_count
# For the largest plot (Horvath), Anchor = 0.
anchor_universe <- max_total_count



for (clk in clocks) {
  # Get CpGs significant for this specific clock
  clock_sig_cpgs <- sig_cpgs %>%
    filter(Clock == clk) %>%
    pull(CpG)

  # Determine the clock definition object
  if (clk == "Horvath") {
    clock_def_cpgs <- coefHorvath$CpGmarker
  } else if (clk == "Hannum") {
    clock_def_cpgs <- coefHannum$CpGmarker
  } else if (clk == "EN") {
    clock_def_cpgs <- coefEN$CpGmarker
  } else if (clk == "PhenoAge") {
    clock_def_cpgs <- coefLevine$CpGmarker
  }

  # Remove intercept if present
  clock_def_cpgs <- clock_def_cpgs[clock_def_cpgs != "(Intercept)"]

  # Filter meQTL lists to this clock's universe
  clock_meqtl_cpgs <- intersect(all_meqtl_cpgs, clock_def_cpgs)

  # Calculate set sizes and overlap for Eulerr
  # Set A: Error-associated
  # Set B: meQTL-affected
  # Overlap: Intersection

  inter <- length(intersect(clock_sig_cpgs, clock_meqtl_cpgs))
  only_clock_sig <- length(setdiff(clock_sig_cpgs, clock_meqtl_cpgs))
  only_meqtl <- length(setdiff(clock_meqtl_cpgs, clock_sig_cpgs))

  # Calculate current total size for this plot
  current_total <- only_clock_sig + only_meqtl + inter

  # Determine required anchor size to reach the global universe size
  needed_anchor <- anchor_universe - current_total

  # Ensure non-negative (shouldn't happen if max_total_count is correct, but safety first)
  if (needed_anchor < 0) needed_anchor <- 0

  # Fit Euler model
  if (needed_anchor > 0) {
    # Add Dummy Anchor
    euler_input <- c(
      "Error-associated clock CpGs" = only_clock_sig,
      "meQTL affected CpGs" = only_meqtl,
      "Error-associated clock CpGs&meQTL affected CpGs" = inter,
      "zz_Anchor" = needed_anchor
    )

    euler_fit <- euler(euler_input)

    # Plot with transparent anchor and legend
    p <- plot(euler_fit,
      quantities = list(cex = 1.8),
      fills = c("lightblue", "lightpink", "transparent"),
      edges = c("black", "black", "transparent"),
      labels = FALSE,
      alpha = 0.6,
      main = paste0(clk, " Clock"),
      legend = list(labels = c("Error-associated clock CpGs", "meQTL affected CpGs"))
    )
  } else {
    # No Anchor needed (e.g. for the largest clock)
    euler_input <- c(
      "Error-associated clock CpGs" = only_clock_sig,
      "meQTL affected CpGs" = only_meqtl,
      "Error-associated clock CpGs&meQTL affected CpGs" = inter
    )

    euler_fit <- euler(euler_input)

    # Plot with legend
    p <- plot(euler_fit,
      quantities = list(cex = 1.8),
      fills = c("lightblue", "lightpink"),
      labels = FALSE,
      alpha = 0.6,
      main = paste0(clk, " Clock"),
      legend = list(labels = c("Error-associated clock CpGs", "meQTL affected CpGs"))
    )
  }

  # Capture as grob to freeze dimensions
  venn_plots[[clk]] <- grid::grid.grabExpr(print(p))
}



# Optional: Save plots
# Assemble using gridExtra for robust sizing
library(gridExtra)

# Layout matching "AABB / #CC#" (Horvath, Hannum / EN)
# Note: grid.arrange doesn't support string layouts directly like patchwork
# We'll use a layout matrix
# 1=Horvath, 2=Hannum, 3=EN
# Structure:
# 1 1 2 2
# 0 3 3 0
layout_mat <- rbind(
  c(1, 1, 2, 2),
  c(NA, 3, 3, NA)
)

# Filter plots to only include Horvath, Hannum, EN for this specific combined plot
# (If we run this loop for clocks including PhenoAge, ensure only the desired ones are in the list passed)
plot_list_2way <- venn_plots[c("Horvath", "Hannum", "EN")]

# Draw
grid.arrange(grobs = plot_list_2way, layout_matrix = layout_mat)

################ Combine with Oge's results
# Load the AFR vs EUR differentially methylated sites defined by Oge
load("~/Desktop/Capra Lab/Thesis Project/Aim_1/ancestry_CpGs.RData")

# Apply significance filter
afr_vs_eur_filtered <- subset(AFR_vs_EUR, adj.P.Val < 0.05)

# Final AFR vs EUR CpG list
afr_dm_cpgs <- rownames(afr_vs_eur_filtered)

# --- 4-Way Venn Diagrams: Error, meQTL, AFR meQTL, AFR vs EUR ---
clocks_4way <- unique(sig_cpgs$Clock)
clocks_4way <- clocks_4way[clocks_4way != "PhenoAge"]

venn_plots_4way <- list()

for (clk in clocks_4way) {
  # Get CpGs significant for this specific clock
  clock_sig_cpgs <- sig_cpgs %>%
    filter(Clock == clk) %>%
    pull(CpG)

  # Get clock definition CpGs
  if (clk == "Horvath") {
    clock_def_cpgs <- coefHorvath$CpGmarker
  } else if (clk == "Hannum") {
    clock_def_cpgs <- coefHannum$CpGmarker
  } else if (clk == "EN") {
    clock_def_cpgs <- coefEN$CpGmarker
  } else {
    next
  }

  # Remove intercept
  clock_def_cpgs <- clock_def_cpgs[clock_def_cpgs != "(Intercept)"]

  # Define sets (restricted to clock definition)
  set_error <- clock_sig_cpgs
  set_meqtl <- intersect(all_meqtl_cpgs, clock_def_cpgs)
  set_afr_meqtl <- intersect(afr_diff_meqtl_cpgs, clock_def_cpgs)
  set_afr_dm <- intersect(afr_dm_cpgs, clock_def_cpgs)

  # Prepare list for eulerr
  euler_sets <- list(
    "Error" = set_error,
    "meQTL" = set_meqtl,
    "AFR meQTL" = set_afr_meqtl,
    "AFR vs EUR" = set_afr_dm
  )

  fit <- euler(euler_sets)

  # Plot normally
  p <- plot(fit,
    quantities = list(type = c("counts"), fontsize = 8, col = "black"),
    fills = c("lightgreen", "wheat", "lightblue", "lightpink"),
    labels = FALSE,
    alpha = 0.6,
    main = paste0(clk, " Clock"),
    legend = list(labels = c("Error", "meQTL", "AFR\nmeQTL", "AFR vs\nEUR"), side = "right", fontface = "plain")
  )

  # Capture as grob
  venn_plots_4way[[clk]] <- grid::grid.grabExpr(print(p))
}

# Combine plots using grid.arrange
# Same layout: Horvath, Hannum top; EN bottom center
plot_list_4way <- venn_plots_4way[c("Horvath", "Hannum", "EN")]

grid.arrange(grobs = plot_list_4way, layout_matrix = layout_mat)

# --- Print all intersection regions underlying the 4-way Venn diagrams ---
cat("\n========== 4-Way Venn Diagram: Full Intersection Breakdown ==========\n")

for (clk in clocks_4way) {
  # Get CpGs significant for this specific clock
  clock_sig_cpgs <- sig_cpgs %>%
    filter(Clock == clk) %>%
    pull(CpG)

  # Get clock definition CpGs
  if (clk == "Horvath") {
    clock_def_cpgs <- coefHorvath$CpGmarker
  } else if (clk == "Hannum") {
    clock_def_cpgs <- coefHannum$CpGmarker
  } else if (clk == "EN") {
    clock_def_cpgs <- coefEN$CpGmarker
  } else {
    next
  }

  clock_def_cpgs <- clock_def_cpgs[clock_def_cpgs != "(Intercept)"]

  # Define the 4 sets (same as the Venn loop)
  sets <- list(
    Error       = clock_sig_cpgs,
    meQTL       = intersect(all_meqtl_cpgs, clock_def_cpgs),
    AFR_meQTL   = intersect(afr_diff_meqtl_cpgs, clock_def_cpgs),
    AFR_vs_EUR  = intersect(afr_dm_cpgs, clock_def_cpgs)
  )

  set_names <- names(sets)

  # Universe: union of all 4 sets
  universe <- unique(unlist(sets))

  # Enumerate all 2^4 = 16 regions
  n_sets <- length(sets)
  regions <- list()

  for (i in 0:(2^n_sets - 1)) {
    membership <- as.logical(intToBits(i)[1:n_sets])
    if (!any(membership)) next # skip the empty region (in none of the sets)

    # CpGs that are IN every TRUE set and NOT IN every FALSE set
    cpgs_in_region <- universe
    for (j in seq_along(sets)) {
      if (membership[j]) {
        cpgs_in_region <- intersect(cpgs_in_region, sets[[j]])
      } else {
        cpgs_in_region <- setdiff(cpgs_in_region, sets[[j]])
      }
    }

    if (length(cpgs_in_region) > 0) {
      label <- paste(set_names[membership], collapse = " ∩ ")
      only_label <- paste0(
        ifelse(membership, set_names, paste0("!", set_names)),
        collapse = " & "
      )
      regions[[length(regions) + 1]] <- data.frame(
        Region = only_label,
        Count = length(cpgs_in_region),
        stringsAsFactors = FALSE
      )
    }
  }

  region_df <- do.call(rbind, regions)
  region_df <- region_df[order(-region_df$Count), ]

  cat(paste0("\n--- ", clk, " Clock ---\n"))
  cat(sprintf("  %-60s %s\n", "Region", "Count"))
  cat(paste0("  ", strrep("-", 66), "\n"))
  for (r in seq_len(nrow(region_df))) {
    cat(sprintf("  %-60s %d\n", region_df$Region[r], region_df$Count[r]))
  }

  # Print per-set totals
  cat("\n  Set totals:\n")
  for (s in set_names) {
    cat(sprintf("    %-20s %d\n", s, length(sets[[s]])))
  }
  cat(sprintf("    %-20s %d\n", "Universe (union)", length(universe)))
  cat("\n")
}



##################### Enrichment analysis
# --- ANY meQTL Enrichment Analysis ---

# Horvath Clock
clock_cpgs <- sig_cpgs %>%
  filter(Clock == "Horvath") %>%
  pull(CpG)
all_clock_cpgs <- c(coefHorvath$CpGmarker)
all_clock_cpgs <- all_clock_cpgs[-1]

error_with_meqtl <- sum(clock_cpgs %in% all_meqtl_cpgs)
error_no_meqtl <- sum(!clock_cpgs %in% all_meqtl_cpgs)
no_error_with_meqtl <- sum(all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% all_meqtl_cpgs)
no_error_no_meqtl <- sum(!all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% all_meqtl_cpgs)

contingency_horvath_any <- matrix(
  c(
    error_with_meqtl, error_no_meqtl,
    no_error_with_meqtl, no_error_no_meqtl
  ),
  nrow = 2, byrow = TRUE
)
fisher.test(contingency_horvath_any)

# Hannum Clock
clock_cpgs <- sig_cpgs %>%
  filter(Clock == "Hannum") %>%
  pull(CpG)
all_clock_cpgs <- c(coefHannum$CpGmarker)
all_clock_cpgs <- all_clock_cpgs[-1]

error_with_meqtl <- sum(clock_cpgs %in% all_meqtl_cpgs)
error_no_meqtl <- sum(!clock_cpgs %in% all_meqtl_cpgs)
no_error_with_meqtl <- sum(all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% all_meqtl_cpgs)
no_error_no_meqtl <- sum(!all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% all_meqtl_cpgs)

contingency_hannum_any <- matrix(
  c(
    error_with_meqtl, error_no_meqtl,
    no_error_with_meqtl, no_error_no_meqtl
  ),
  nrow = 2, byrow = TRUE
)
fisher.test(contingency_hannum_any)

# EN Clock
clock_cpgs <- sig_cpgs %>%
  filter(Clock == "EN") %>%
  pull(CpG)
all_clock_cpgs <- c(coefEN$CpGmarker)
all_clock_cpgs <- all_clock_cpgs[-1]

error_with_meqtl <- sum(clock_cpgs %in% all_meqtl_cpgs)
error_no_meqtl <- sum(!clock_cpgs %in% all_meqtl_cpgs)
no_error_with_meqtl <- sum(all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% all_meqtl_cpgs)
no_error_no_meqtl <- sum(!all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% all_meqtl_cpgs)

contingency_en_any <- matrix(
  c(
    error_with_meqtl, error_no_meqtl,
    no_error_with_meqtl, no_error_no_meqtl
  ),
  nrow = 2, byrow = TRUE
)
fisher.test(contingency_en_any)

# PhenoAge Clock
clock_cpgs <- sig_cpgs %>%
  filter(Clock == "PhenoAge") %>%
  pull(CpG)
all_clock_cpgs <- c(coefLevine$CpGmarker)
all_clock_cpgs <- all_clock_cpgs[-1]

error_with_meqtl <- sum(clock_cpgs %in% all_meqtl_cpgs)
error_no_meqtl <- sum(!clock_cpgs %in% all_meqtl_cpgs)
no_error_with_meqtl <- sum(all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% all_meqtl_cpgs)
no_error_no_meqtl <- sum(!all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% all_meqtl_cpgs)

contingency_phenoage_any <- matrix(
  c(
    error_with_meqtl, error_no_meqtl,
    no_error_with_meqtl, no_error_no_meqtl
  ),
  nrow = 2, byrow = TRUE
)
fisher.test(contingency_phenoage_any)



# --- AFR-differentiated meQTL Enrichment Analysis ---

# Horvath Clock
clock_cpgs <- sig_cpgs %>%
  filter(Clock == "Horvath") %>%
  pull(CpG)
all_clock_cpgs <- c(coefHorvath$CpGmarker)
all_clock_cpgs <- all_clock_cpgs[-1]

error_with_afr_meqtl <- sum(clock_cpgs %in% afr_diff_meqtl_cpgs)
error_no_afr_meqtl <- sum(!clock_cpgs %in% afr_diff_meqtl_cpgs)
no_error_with_afr_meqtl <- sum(all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% afr_diff_meqtl_cpgs)
no_error_no_afr_meqtl <- sum(!all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% afr_diff_meqtl_cpgs)

contingency_horvath_afr <- matrix(
  c(
    error_with_afr_meqtl, error_no_afr_meqtl,
    no_error_with_afr_meqtl, no_error_no_afr_meqtl
  ),
  nrow = 2, byrow = TRUE
)
fisher.test(contingency_horvath_afr)

# Hannum Clock
clock_cpgs <- sig_cpgs %>%
  filter(Clock == "Hannum") %>%
  pull(CpG)
all_clock_cpgs <- c(coefHannum$CpGmarker)
all_clock_cpgs <- all_clock_cpgs[-1]

error_with_afr_meqtl <- sum(clock_cpgs %in% afr_diff_meqtl_cpgs)
error_no_afr_meqtl <- sum(!clock_cpgs %in% afr_diff_meqtl_cpgs)
no_error_with_afr_meqtl <- sum(all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% afr_diff_meqtl_cpgs)
no_error_no_afr_meqtl <- sum(!all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% afr_diff_meqtl_cpgs)

contingency_hannum_afr <- matrix(
  c(
    error_with_afr_meqtl, error_no_afr_meqtl,
    no_error_with_afr_meqtl, no_error_no_afr_meqtl
  ),
  nrow = 2, byrow = TRUE
)
fisher.test(contingency_hannum_afr)

# EN Clock
clock_cpgs <- sig_cpgs %>%
  filter(Clock == "EN") %>%
  pull(CpG)
all_clock_cpgs <- c(coefEN$CpGmarker)
all_clock_cpgs <- all_clock_cpgs[-1]

error_with_afr_meqtl <- sum(clock_cpgs %in% afr_diff_meqtl_cpgs)
error_no_afr_meqtl <- sum(!clock_cpgs %in% afr_diff_meqtl_cpgs)
no_error_with_afr_meqtl <- sum(all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% afr_diff_meqtl_cpgs)
no_error_no_afr_meqtl <- sum(!all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% afr_diff_meqtl_cpgs)

contingency_en_afr <- matrix(
  c(
    error_with_afr_meqtl, error_no_afr_meqtl,
    no_error_with_afr_meqtl, no_error_no_afr_meqtl
  ),
  nrow = 2, byrow = TRUE
)
fisher.test(contingency_en_afr)

# PhenoAge Clock
clock_cpgs <- sig_cpgs %>%
  filter(Clock == "PhenoAge") %>%
  pull(CpG)
all_clock_cpgs <- c(coefLevine$CpGmarker)
all_clock_cpgs <- all_clock_cpgs[-1]

error_with_afr_meqtl <- sum(clock_cpgs %in% afr_diff_meqtl_cpgs)
error_no_afr_meqtl <- sum(!clock_cpgs %in% afr_diff_meqtl_cpgs)
no_error_with_afr_meqtl <- sum(all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% afr_diff_meqtl_cpgs)
no_error_no_afr_meqtl <- sum(!all_clock_cpgs[!all_clock_cpgs %in% clock_cpgs] %in% afr_diff_meqtl_cpgs)

contingency_phenoage_afr <- matrix(
  c(
    error_with_afr_meqtl, error_no_afr_meqtl,
    no_error_with_afr_meqtl, no_error_no_afr_meqtl
  ),
  nrow = 2, byrow = TRUE
)
fisher.test(contingency_phenoage_afr)






# --- Calculate OR and CI from contingency tables ---

# Function to calculate OR, CI, and p-value from 2x2 table
calculate_or_stats <- function(cont_table) {
  # Perform Fisher's exact test
  fisher_result <- fisher.test(cont_table)

  # Extract OR and CI
  or <- fisher_result$estimate
  ci_lower <- fisher_result$conf.int[1]
  ci_upper <- fisher_result$conf.int[2]
  p_value <- fisher_result$p.value

  return(list(OR = or, CI_lower = ci_lower, CI_upper = ci_upper, p_value = p_value))
}

# --- ANY meQTL contingency tables ---
contingency_horvath_any <- matrix(c(42, 14, 196, 101), nrow = 2, byrow = TRUE)
contingency_hannum_any <- matrix(c(3, 16, 3, 49), nrow = 2, byrow = TRUE)
contingency_en_any <- matrix(c(1, 51, 6, 456), nrow = 2, byrow = TRUE)
contingency_phenoage_any <- matrix(c(11, 89, 22, 391), nrow = 2, byrow = TRUE)

# Calculate statistics for ANY meQTL
horvath_any <- calculate_or_stats(contingency_horvath_any)
hannum_any <- calculate_or_stats(contingency_hannum_any)
en_any <- calculate_or_stats(contingency_en_any)
phenoage_any <- calculate_or_stats(contingency_phenoage_any)

# Create ANY meQTL data frame
any_meqtl_results <- data.frame(
  Clock = c("Horvath", "Hannum", "EN", "PhenoAge"),
  OR = c(horvath_any$OR, hannum_any$OR, en_any$OR, phenoage_any$OR),
  CI_lower = c(horvath_any$CI_lower, hannum_any$CI_lower, en_any$CI_lower, phenoage_any$CI_lower),
  CI_upper = c(horvath_any$CI_upper, hannum_any$CI_upper, en_any$CI_upper, phenoage_any$CI_upper),
  p_value = c(horvath_any$p_value, hannum_any$p_value, en_any$p_value, phenoage_any$p_value),
  Type = "Any meQTL"
)

# --- AFR-differentiated meQTL contingency tables ---
contingency_horvath_afr <- matrix(c(9, 47, 43, 254), nrow = 2, byrow = TRUE)
contingency_hannum_afr <- matrix(c(0, 19, 1, 51), nrow = 2, byrow = TRUE)
contingency_en_afr <- matrix(c(0, 52, 3, 459), nrow = 2, byrow = TRUE)
contingency_phenoage_afr <- matrix(c(2, 98, 6, 407), nrow = 2, byrow = TRUE)

# Calculate statistics for AFR-differentiated meQTL
horvath_afr <- calculate_or_stats(contingency_horvath_afr)
hannum_afr <- calculate_or_stats(contingency_hannum_afr)
en_afr <- calculate_or_stats(contingency_en_afr)
phenoage_afr <- calculate_or_stats(contingency_phenoage_afr)

# Create AFR-differentiated meQTL data frame
afr_meqtl_results <- data.frame(
  Clock = c("Horvath", "Hannum", "EN", "PhenoAge"),
  OR = c(horvath_afr$OR, hannum_afr$OR, en_afr$OR, phenoage_afr$OR),
  CI_lower = c(horvath_afr$CI_lower, hannum_afr$CI_lower, en_afr$CI_lower, phenoage_afr$CI_lower),
  CI_upper = c(horvath_afr$CI_upper, hannum_afr$CI_upper, en_afr$CI_upper, phenoage_afr$CI_upper),
  p_value = c(horvath_afr$p_value, hannum_afr$p_value, en_afr$p_value, phenoage_afr$p_value),
  Type = "AFR-differentiated meQTL"
)

# Handle any infinite values for plotting
afr_meqtl_results$OR[is.infinite(afr_meqtl_results$OR)] <- NA
afr_meqtl_results$CI_upper[is.infinite(afr_meqtl_results$CI_upper)] <- NA

# Combine datasets
combined_results <- rbind(any_meqtl_results, afr_meqtl_results)

# Reorder clocks
combined_results$Clock <- factor(combined_results$Clock,
  levels = c("Horvath", "Hannum", "EN", "PhenoAge")
)

# Print results table
print(combined_results)

# --- Create forest plot ---
p <- ggplot(combined_results, aes(x = Clock, y = OR, color = Type)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  scale_color_manual(values = c(
    "Any meQTL" = "#2166ac",
    "AFR-differentiated meQTL" = "#b2182b"
  )) +
  scale_y_log10(limits = c(0.01, 100)) +
  labs(
    x = "Epigenetic Clock",
    y = "Odds Ratio (log scale)",
    color = "meQTL Type"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 11)
  )

print(p)
