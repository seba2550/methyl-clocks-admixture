# ==============================================================================
# 05_magenta_age_acceleration.R
# Computes epigenetic age predictions and age acceleration (raw and residual)
# for EUR, AFR-Am, Combined 50/50, and 3 reference clocks on the MAGENTA cohort.
# Evaluates associations with AD status overall and stratified by cohort.
# MEMORY OPTIMIZED VERSION: Subsets beta matrix to clock CpGs immediately.
# ==============================================================================

library(dplyr)
library(tidyr)
library(readxl)
library(methylclock)
library(ggplot2)

cat("===========================================================\n")
cat("  Analyzing Epigenetic Age Acceleration in MAGENTA Cohort\n")
cat("===========================================================\n\n")

proj_dir <- "/Users/sgonzalez/Desktop/Capra Lab/Thesis Project"
aim1_dir <- file.path(proj_dir, "Aim_1")
results_dir <- file.path(proj_dir, "results")
evaluation_dir <- file.path(results_dir, "evaluation")
figures_dir <- file.path(results_dir, "figures/magenta_age_acceleration")

# Ensure output directories exist
dir.create(evaluation_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Helper Functions ----

# Standard Horvath inverse transformation
inverse_transform_age <- function(t_age, adult_age = 20) {
    ifelse(t_age <= 0, exp(t_age + log(adult_age + 1)) - 1,
        t_age * (adult_age + 1) + adult_age
    )
}

# Helper to predict age from coefficients
predict_custom_clock <- function(coefs_df, beta_test) {
    intercept <- coefs_df$coefficient[coefs_df$cpg == "(Intercept)"]
    if (length(intercept) == 0) intercept <- 0
    cpg_coefs <- coefs_df[coefs_df$cpg != "(Intercept)", ]
    model_cpgs <- cpg_coefs$cpg
    
    available <- intersect(model_cpgs, colnames(beta_test))
    
    # Initialize samples x CpGs matrix with 0
    X <- matrix(0, nrow = nrow(beta_test), ncol = length(model_cpgs))
    colnames(X) <- model_cpgs
    
    if (length(available) > 0) {
        X[, available] <- as.matrix(beta_test[, available])
    }
    
    preds_t <- as.numeric(X %*% cpg_coefs$coefficient) + intercept
    return(inverse_transform_age(preds_t))
}

# ==========================================
# Step 1: Prepare CpG lists to save memory
# ==========================================

cat("Loading reference clock databases for CpG list...\n")
load_DNAm_Clocks_data()

# Load Coefficients for Custom Clocks to retrieve CpG names
cat("Loading custom clock coefficients to retrieve CpG names...\n")
coefs_lc <- read.csv(file.path(results_dir, "learning_curve_v2_coefs.csv"), stringsAsFactors = FALSE)
coefs_eur <- coefs_lc %>% filter(Clock_Type == "EUR", N == 1650, Iteration == 29)
coefs_afr <- coefs_lc %>% filter(Clock_Type == "AFR-Am", N == 1650, Iteration == 1)

coefs_comp <- read.csv(file.path(results_dir, "composition_v2_coefs.csv"), stringsAsFactors = FALSE)
coefs_comb <- coefs_comp %>% filter(Ratio == "50/50", Iteration == 3)

# Build unified CpG needed set
cpgs_needed <- unique(c(
    coefHorvath$CpGmarker,
    coefHannum$CpGmarker,
    coefEN$CpGmarker,
    coefs_eur$cpg,
    coefs_afr$cpg,
    coefs_comb$cpg
))
# Exclude intercept terms
cpgs_needed <- cpgs_needed[!cpgs_needed %in% c("(Intercept)", "Intercept")]

cat(sprintf("Total distinct CpGs needed for all 6 clocks: %d\n", length(cpgs_needed)))

# ==========================================
# Step 2: Load and subset MAGENTA data
# ==========================================

cat("\nLoading MAGENTA beta matrix (RDS)...\n")
beta_magenta <- readRDS(file.path(aim1_dir, "betaMatrices/normalizedBetas/beta_QGCDPB_combined.rds"))

# Clean CpG names immediately on the rows
rownames(beta_magenta) <- gsub("_.*", "", rownames(beta_magenta))

# Perform row subsetting immediately to drop 99.9% of rows and save memory
cat("Filtering beta matrix to target CpGs immediately...\n")
available_cpgs <- intersect(rownames(beta_magenta), cpgs_needed)
beta_magenta <- beta_magenta[available_cpgs, , drop = FALSE]

cat(sprintf("Filtered beta matrix dimensions: %d CpGs x %d samples\n", nrow(beta_magenta), ncol(beta_magenta)))

# Force garbage collection
gc(verbose = FALSE)

# Load MAGENTA metadata
cat("Loading MAGENTA metadata...\n")
meta_magenta <- read_xlsx(file.path(aim1_dir, "ADmethy_pheno.xlsx"))

# Align samples between metadata and beta matrix
meta_magenta <- meta_magenta %>%
    filter(Beta_ID %in% colnames(beta_magenta)) %>%
    filter(!is.na(AGE_OF_EXAM)) %>%
    filter(!is.na(STATUS))

cat(sprintf("Aligned metadata: %d samples\n", nrow(meta_magenta)))

# Extract aligned beta matrix subset and transpose for custom clocks (samples x CpGs)
beta_subset <- beta_magenta[, meta_magenta$Beta_ID, drop = FALSE] # CpGs x samples
beta_test <- t(beta_subset) # samples x CpGs

# Free full beta_magenta matrix
rm(beta_magenta)
gc(verbose = FALSE)

# ==========================================
# Step 3: Run Predictions
# ==========================================

cat("Predicting custom clocks (EUR, AFR-Am, Combined 50/50)...\n")
meta_magenta$pred_EUR <- predict_custom_clock(coefs_eur, beta_test)
meta_magenta$pred_AFR_Am <- predict_custom_clock(coefs_afr, beta_test)
meta_magenta$pred_Combined <- predict_custom_clock(coefs_comb, beta_test)

# Free transposed matrix to save memory
rm(beta_test)
gc(verbose = FALSE)

# Run Predictions for Reference Clocks using methylclock
cat("Predicting reference clocks (Horvath, Hannum, Zhang EN) via DNAmAge...\n")
preds_ref <- DNAmAge(
    beta_subset, 
    clocks = c("Horvath", "Hannum", "EN"), 
    cell.count = FALSE, 
    normalize = FALSE
)

# Merge reference predictions with metadata
meta_magenta <- meta_magenta %>%
    left_join(
        preds_ref %>% select(id, Horvath, Hannum, EN), 
        by = c("Beta_ID" = "id")
    ) %>%
    rename(
        pred_Horvath = Horvath,
        pred_Hannum = Hannum,
        pred_Zhang_EN = EN
    )

# Free remaining beta matrix subset
rm(beta_subset, preds_ref)
gc(verbose = FALSE)

# ==========================================
# Step 4: Compute Epigenetic Age Acceleration
# ==========================================

cat("\nComputing raw delta and residual age acceleration...\n")
clocks <- c("EUR", "AFR_Am", "Combined", "Horvath", "Hannum", "Zhang_EN")

for (clock in clocks) {
    pred_col <- paste0("pred_", clock)
    raw_acc_col <- paste0("raw_acc_", clock)
    res_acc_col <- paste0("res_acc_", clock)
    
    # 1. Raw Delta Age
    meta_magenta[[raw_acc_col]] <- meta_magenta[[pred_col]] - meta_magenta$AGE_OF_EXAM
    
    # 2. Residual Age Acceleration (residuals of Pred_Age on Chronological Age)
    fit <- lm(as.formula(paste(pred_col, "~ AGE_OF_EXAM")), data = meta_magenta, na.action = na.exclude)
    meta_magenta[[res_acc_col]] <- residuals(fit)
}

# Save full predictions to file for trace
write.csv(meta_magenta, file.path(evaluation_dir, "magenta_all_predicted_ages_and_acceleration.csv"), row.names = FALSE)
cat(sprintf("All predictions and acceleration scores exported to:\n  %s\n\n", 
            file.path(evaluation_dir, "magenta_all_predicted_ages_and_acceleration.csv")))

# ==========================================
# Step 5: Overall Association with AD Status
# ==========================================

cat("Evaluating overall associations with AD status...\n")
overall_results <- list()

for (clock in clocks) {
    for (acc_type in c("raw_acc", "res_acc")) {
        acc_col <- paste0(acc_type, "_", clock)
        
        # Filter NAs
        df_sub <- meta_magenta %>% filter(!is.na(.[[acc_col]]))
        
        # Calculate summary statistics
        stats <- df_sub %>%
            group_by(STATUS) %>%
            summarise(
                mean_val = mean(.data[[acc_col]]),
                median_val = median(.data[[acc_col]]),
                sd_val = sd(.data[[acc_col]]),
                n_val = n(),
                .groups = "drop"
            )
        
        mean_ad <- stats$mean_val[stats$STATUS == "AD"]
        mean_control <- stats$mean_val[stats$STATUS == "CONTROL"]
        median_ad <- stats$median_val[stats$STATUS == "AD"]
        median_control <- stats$median_val[stats$STATUS == "CONTROL"]
        sd_ad <- stats$sd_val[stats$STATUS == "AD"]
        sd_control <- stats$sd_val[stats$STATUS == "CONTROL"]
        
        # Two-sample t-test (Welch)
        ttest_res <- t.test(as.formula(paste(acc_col, "~ STATUS")), data = df_sub)
        p_ttest <- ttest_res$p.value
        
        # Linear Model: Age_Acceleration ~ STATUS + SEX
        df_sub$STATUS <- factor(df_sub$STATUS, levels = c("CONTROL", "AD"))
        df_sub$SEX <- factor(df_sub$SEX)
        
        lm_fit <- lm(as.formula(paste(acc_col, "~ STATUS + SEX")), data = df_sub)
        lm_summary <- summary(lm_fit)
        coef_names <- rownames(lm_summary$coefficients)
        
        # Extract STATUS AD coefficient
        status_row <- grep("STATUSAD", coef_names, value = TRUE)
        if (length(status_row) > 0) {
            status_beta <- lm_summary$coefficients[status_row, "Estimate"]
            status_se   <- lm_summary$coefficients[status_row, "Std. Error"]
            status_p    <- lm_summary$coefficients[status_row, "Pr(>|t|)"]
        } else {
            status_beta <- NA
            status_se   <- NA
            status_p    <- NA
        }
        
        # Extract SEX coefficient
        sex_row <- grep("SEX", coef_names, value = TRUE)
        if (length(sex_row) > 0) {
            sex_beta <- lm_summary$coefficients[sex_row, "Estimate"]
            sex_p    <- lm_summary$coefficients[sex_row, "Pr(>|t|)"]
        } else {
            sex_beta <- NA
            sex_p    <- NA
        }
        
        overall_results[[length(overall_results) + 1]] <- data.frame(
            Clock = clock,
            Acceleration_Type = ifelse(acc_type == "raw_acc", "Raw_Delta", "Residual"),
            Mean_AD = mean_ad,
            Mean_Control = mean_control,
            Median_AD = median_ad,
            Median_Control = median_control,
            SD_AD = sd_ad,
            SD_Control = sd_control,
            T_Test_P = p_ttest,
            LM_STATUS_AD_Beta = status_beta,
            LM_STATUS_AD_SE = status_se,
            LM_STATUS_AD_P = status_p,
            LM_SEX_Beta = sex_beta,
            LM_SEX_P = sex_p,
            stringsAsFactors = FALSE
        )
    }
}

overall_results_df <- do.call(rbind, overall_results)
write.csv(overall_results_df, file.path(evaluation_dir, "magenta_overall_age_acceleration.csv"), row.names = FALSE)
cat(sprintf("Overall metrics exported to:\n  %s\n\n", 
            file.path(evaluation_dir, "magenta_overall_age_acceleration.csv")))

# ==========================================
# Step 5b: Overall Associations (Excluding Cubans/CuADI)
# ==========================================
cat("Evaluating overall associations with AD status (Excluding Cubans/CuADI)...\n")
overall_nocubans_results <- list()

for (clock in clocks) {
    for (acc_type in c("raw_acc", "res_acc")) {
        acc_col <- paste0(acc_type, "_", clock)
        
        df_sub <- meta_magenta %>% filter(COHORT != "CuADI", !is.na(.data[[acc_col]]))
        
        stats <- df_sub %>%
            group_by(STATUS) %>%
            summarise(
                mean_val = mean(.data[[acc_col]]),
                median_val = median(.data[[acc_col]]),
                sd_val = sd(.data[[acc_col]]),
                n_val = n(),
                .groups = "drop"
            )
        
        mean_ad <- stats$mean_val[stats$STATUS == "AD"]
        mean_control <- stats$mean_val[stats$STATUS == "CONTROL"]
        median_ad <- stats$median_val[stats$STATUS == "AD"]
        median_control <- stats$median_val[stats$STATUS == "CONTROL"]
        sd_ad <- stats$sd_val[stats$STATUS == "AD"]
        sd_control <- stats$sd_val[stats$STATUS == "CONTROL"]
        
        # Two-sample t-test (Welch)
        ttest_res <- t.test(as.formula(paste(acc_col, "~ STATUS")), data = df_sub)
        p_ttest <- ttest_res$p.value
        
        # Linear Model: Age_Acceleration ~ STATUS + SEX
        df_sub$STATUS <- factor(df_sub$STATUS, levels = c("CONTROL", "AD"))
        df_sub$SEX <- factor(df_sub$SEX)
        
        lm_fit <- lm(as.formula(paste(acc_col, "~ STATUS + SEX")), data = df_sub)
        lm_summary <- summary(lm_fit)
        coef_names <- rownames(lm_summary$coefficients)
        
        status_row <- grep("STATUSAD", coef_names, value = TRUE)
        if (length(status_row) > 0) {
            status_beta <- lm_summary$coefficients[status_row, "Estimate"]
            status_se   <- lm_summary$coefficients[status_row, "Std. Error"]
            status_p    <- lm_summary$coefficients[status_row, "Pr(>|t|)"]
        } else {
            status_beta <- NA; status_se <- NA; status_p <- NA
        }
        
        sex_row <- grep("SEX", coef_names, value = TRUE)
        if (length(sex_row) > 0) {
            sex_beta <- lm_summary$coefficients[sex_row, "Estimate"]
            sex_p    <- lm_summary$coefficients[sex_row, "Pr(>|t|)"]
        } else {
            sex_beta <- NA; sex_p <- NA
        }
        
        overall_nocubans_results[[length(overall_nocubans_results) + 1]] <- data.frame(
            Clock = clock,
            Acceleration_Type = ifelse(acc_type == "raw_acc", "Raw_Delta", "Residual"),
            Mean_AD = mean_ad,
            Mean_Control = mean_control,
            Median_AD = median_ad,
            Median_Control = median_control,
            SD_AD = sd_ad,
            SD_Control = sd_control,
            T_Test_P = p_ttest,
            LM_STATUS_AD_Beta = status_beta,
            LM_STATUS_AD_SE = status_se,
            LM_STATUS_AD_P = status_p,
            LM_SEX_Beta = sex_beta,
            LM_SEX_P = sex_p,
            stringsAsFactors = FALSE
        )
    }
}

overall_nocubans_df <- do.call(rbind, overall_nocubans_results)
write.csv(overall_nocubans_df, file.path(evaluation_dir, "magenta_overall_age_acceleration_no_cubans.csv"), row.names = FALSE)
cat(sprintf("Overall metrics (excluding Cubans) exported to:\n  %s\n\n", 
            file.path(evaluation_dir, "magenta_overall_age_acceleration_no_cubans.csv")))

# ==========================================
# Step 6: Cohort-Stratified Association
# ==========================================

cat("Evaluating cohort-stratified associations...\n")
stratified_results <- list()

for (clock in clocks) {
    for (acc_type in c("raw_acc", "res_acc")) {
        acc_col <- paste0(acc_type, "_", clock)
        
        # Prepare subset data
        df_sub <- meta_magenta %>% filter(!is.na(.[[acc_col]]))
        df_sub$STATUS <- factor(df_sub$STATUS, levels = c("CONTROL", "AD"))
        df_sub$SEX <- factor(df_sub$SEX)
        df_sub$COHORT <- factor(df_sub$COHORT)
        
        # Combined model with STATUS * COHORT interaction to test if association differs by cohort
        m_null <- lm(as.formula(paste(acc_col, "~ STATUS + COHORT + SEX")), data = df_sub)
        m_int  <- lm(as.formula(paste(acc_col, "~ STATUS * COHORT + SEX")), data = df_sub)
        anova_res <- anova(m_null, m_int)
        interaction_p <- anova_res[2, "Pr(>F)"]
        
        # Process each cohort subset
        cohorts_list <- sort(unique(meta_magenta$COHORT))
        for (cohort_val in cohorts_list) {
            df_cohort <- df_sub %>% filter(COHORT == cohort_val)
            
            # Summary statistics
            stats_cohort <- df_cohort %>%
                group_by(STATUS) %>%
                summarise(
                    mean_val = mean(.data[[acc_col]]),
                    median_val = median(.data[[acc_col]]),
                    sd_val = sd(.data[[acc_col]]),
                    n_val = n(),
                    .groups = "drop"
                )
            
            mean_ad <- stats_cohort$mean_val[stats_cohort$STATUS == "AD"]
            mean_control <- stats_cohort$mean_val[stats_cohort$STATUS == "CONTROL"]
            median_ad <- stats_cohort$median_val[stats_cohort$STATUS == "AD"]
            median_control <- stats_cohort$median_val[stats_cohort$STATUS == "CONTROL"]
            sd_ad <- stats_cohort$sd_val[stats_cohort$STATUS == "AD"]
            sd_control <- stats_cohort$sd_val[stats_cohort$STATUS == "CONTROL"]
            
            # Two-sample t-test (Welch)
            ttest_res <- t.test(as.formula(paste(acc_col, "~ STATUS")), data = df_cohort)
            p_ttest <- ttest_res$p.value
            
            # Linear model on cohort subset
            lm_cohort_fit <- lm(as.formula(paste(acc_col, "~ STATUS + SEX")), data = df_cohort)
            lm_cohort_summary <- summary(lm_cohort_fit)
            coef_names_c <- rownames(lm_cohort_summary$coefficients)
            
            # Extract STATUS AD coefficient
            status_row_c <- grep("STATUSAD", coef_names_c, value = TRUE)
            if (length(status_row_c) > 0) {
                status_beta_c <- lm_cohort_summary$coefficients[status_row_c, "Estimate"]
                status_se_c   <- lm_cohort_summary$coefficients[status_row_c, "Std. Error"]
                status_p_c    <- lm_cohort_summary$coefficients[status_row_c, "Pr(>|t|)"]
            } else {
                status_beta_c <- NA
                status_se_c   <- NA
                status_p_c    <- NA
            }
            
            # Extract SEX coefficient
            sex_row_c <- grep("SEX", coef_names_c, value = TRUE)
            if (length(sex_row_c) > 0) {
                sex_beta_c <- lm_cohort_summary$coefficients[sex_row_c, "Estimate"]
                sex_p_c    <- lm_cohort_summary$coefficients[sex_row_c, "Pr(>|t|)"]
            } else {
                sex_beta_c <- NA
                sex_p_c    <- NA
            }
            
            stratified_results[[length(stratified_results) + 1]] <- data.frame(
                Cohort = cohort_val,
                Clock = clock,
                Acceleration_Type = ifelse(acc_type == "raw_acc", "Raw_Delta", "Residual"),
                Mean_AD = mean_ad,
                Mean_Control = mean_control,
                Median_AD = median_ad,
                Median_Control = median_control,
                SD_AD = sd_ad,
                SD_Control = sd_control,
                T_Test_P = p_ttest,
                LM_STATUS_AD_Beta = status_beta_c,
                LM_STATUS_AD_SE = status_se_c,
                LM_STATUS_AD_P = status_p_c,
                LM_SEX_Beta = sex_beta_c,
                LM_SEX_P = sex_p_c,
                Combined_Interaction_P = interaction_p,
                stringsAsFactors = FALSE
            )
        }
    }
}

stratified_results_df <- do.call(rbind, stratified_results)
write.csv(stratified_results_df, file.path(evaluation_dir, "magenta_stratified_age_acceleration.csv"), row.names = FALSE)
cat(sprintf("Cohort-stratified metrics exported to:\n  %s\n\n", 
            file.path(evaluation_dir, "magenta_stratified_age_acceleration.csv")))


# ==========================================
# Step 7: Visualizations
# ==========================================

cat("Generating visual plots...\n")

# Color palette definition: Slate Teal for Control, Terracotta Coral for AD
status_colors <- c("CONTROL" = "#2b7b8a", "AD" = "#d95f02")

# Create standard clean ggplot theme
custom_theme <- theme_bw() +
    theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b = 10)),
        strip.background = element_rect(fill = "gray95", color = "gray80"),
        strip.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold", size = 11),
        axis.text = element_text(size = 9)
    )

# Format clock names for titles
clock_title_map <- c(
    "EUR" = "EUR Clock (Custom)",
    "AFR_Am" = "AFR-Am Clock (Custom)",
    "Combined" = "Combined 50/50 Clock (Custom)",
    "Horvath" = "Horvath Clock (Reference)",
    "Hannum" = "Hannum Clock (Reference)",
    "Zhang_EN" = "Zhang EN Clock (Reference)"
)

# Helper to generate significance brackets for boxplots using Linear Model p-values (adjusting for SEX)
get_sig_brackets <- function(data, acc_col) {
    cohorts_list <- sort(unique(data$COHORT))
    brackets <- list()
    labels <- list()
    
    for (cohort_val in cohorts_list) {
        df_cohort <- data %>% filter(COHORT == cohort_val, !is.na(.data[[acc_col]]))
        if (nrow(df_cohort) < 3) next
        
        # Calculate Linear Model p-value (adjusting for SEX)
        df_cohort$STATUS <- factor(df_cohort$STATUS, levels = c("CONTROL", "AD"))
        df_cohort$SEX <- factor(df_cohort$SEX)
        
        if (length(unique(df_cohort$SEX)) > 1) {
            lm_fit <- tryCatch(
                lm(as.formula(paste(acc_col, "~ STATUS + SEX")), data = df_cohort),
                error = function(e) NULL
            )
        } else {
            lm_fit <- tryCatch(
                lm(as.formula(paste(acc_col, "~ STATUS")), data = df_cohort),
                error = function(e) NULL
            )
        }
        
        if (is.null(lm_fit)) next
        
        lm_summary <- summary(lm_fit)
        coef_names <- rownames(lm_summary$coefficients)
        status_row <- grep("STATUSAD", coef_names, value = TRUE)
        
        if (length(status_row) > 0) {
            p_val <- lm_summary$coefficients[status_row, "Pr(>|t|)"]
        } else {
            p_val <- NA
        }
        
        if (!is.na(p_val) && p_val < 0.05) {
            if (p_val < 0.001) {
                annot <- "***"
            } else if (p_val < 0.01) {
                annot <- "**"
            } else {
                annot <- "*"
            }
            annot_text <- sprintf("%s\n(p = %.3f)", annot, p_val)
            
            # Find max y value in this cohort to position the bracket
            y_cohort_max <- max(df_cohort[[acc_col]], na.rm = TRUE)
            y_cohort_min <- min(df_cohort[[acc_col]], na.rm = TRUE)
            y_range <- y_cohort_max - y_cohort_min
            if (y_range == 0) y_range <- 1
            
            y_pos <- y_cohort_max + 0.05 * y_range
            tick_height <- 0.03 * y_range
            
            bracket_coords <- data.frame(
                COHORT = cohort_val,
                x = c(1, 1, 2, 2),
                y = c(y_pos - tick_height, y_pos, y_pos, y_pos - tick_height),
                group = paste0(cohort_val, "_", acc_col),
                stringsAsFactors = FALSE
            )
            
            label_coords <- data.frame(
                COHORT = cohort_val,
                x = 1.5,
                y = y_pos + 0.02 * y_range,
                label = annot_text,
                stringsAsFactors = FALSE
            )
            
            brackets[[length(brackets) + 1]] <- bracket_coords
            labels[[length(labels) + 1]] <- label_coords
        }
    }
    
    return(list(
        path = if (length(brackets) > 0) do.call(rbind, brackets) else NULL,
        text = if (length(labels) > 0) do.call(rbind, labels) else NULL
    ))
}

# Render and save plots
for (clock in clocks) {
    clock_title <- clock_title_map[clock]
    
    # 1. Faceted Boxplots
    res_acc_col <- paste0("res_acc_", clock)
    
    # Prepare non-NA data for plotting
    plot_df <- meta_magenta %>% filter(!is.na(.[[res_acc_col]]))
    plot_df$STATUS <- factor(plot_df$STATUS, levels = c("AD", "CONTROL"))
    
    # Compute significance brackets
    brackets_data <- get_sig_brackets(plot_df, res_acc_col)
    
    p_box <- ggplot(plot_df, aes(x = STATUS, y = .data[[res_acc_col]], fill = STATUS)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "gray20", width = 0.5) +
        geom_jitter(aes(color = STATUS), width = 0.2, alpha = 0.4, size = 1.2) +
        facet_wrap(~ COHORT, nrow = 1) +
        scale_fill_manual(values = status_colors) +
        scale_color_manual(values = status_colors) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
        labs(
            title = paste("Residual Age Acceleration by AD Status -", clock_title),
            x = "Status",
            y = "Residual Age Acceleration (Years)",
            fill = "AD Status",
            color = "AD Status"
        ) +
        custom_theme
        
    # Add brackets if present
    if (!is.null(brackets_data$path)) {
        p_box <- p_box + 
            geom_path(data = brackets_data$path, aes(x = x, y = y, group = group), inherit.aes = FALSE, color = "gray30", linewidth = 0.5)
    }
    if (!is.null(brackets_data$text)) {
        p_box <- p_box + 
            geom_text(data = brackets_data$text, aes(x = x, y = y, label = label), inherit.aes = FALSE, color = "gray20", fontface = "bold", size = 3, vjust = 0, lineheight = 0.8)
    }
    
    box_filename <- file.path(figures_dir, paste0("boxplots_", clock, ".png"))
    ggsave(box_filename, plot = p_box, width = 10, height = 4.5, dpi = 300)
    
    # 2. Faceted Scatterplots
    pred_col <- paste0("pred_", clock)
    
    p_scatter <- ggplot(plot_df, aes(x = AGE_OF_EXAM, y = .data[[pred_col]], color = STATUS)) +
        geom_point(alpha = 0.5, size = 1.5) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.8) +
        geom_smooth(method = "lm", se = TRUE, alpha = 0.15, formula = y ~ x) +
        facet_wrap(~ COHORT, nrow = 1) +
        scale_color_manual(values = status_colors) +
        labs(
            title = paste("Chronological Age vs. Predicted Age -", clock_title),
            x = "Chronological Age (AGE_OF_EXAM)",
            y = "Predicted Epigenetic Age (Years)",
            color = "AD Status"
        ) +
        custom_theme
    
    scatter_filename <- file.path(figures_dir, paste0("scatterplots_", clock, ".png"))
    ggsave(scatter_filename, plot = p_scatter, width = 10, height = 4.5, dpi = 300)
}

cat("\nAnalysis and visualization generation completed successfully!\n")
cat("===========================================================\n")
