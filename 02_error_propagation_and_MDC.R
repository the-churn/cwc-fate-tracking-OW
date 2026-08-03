# ==============================================================================
# Script Name:  02_error_propagation_and_MDC.R
# Description:  Models empirical TagLab measurement error using repeat CNN 
#               annotations (CNN Error.xlsx), classifies ecological state 
#               transitions (persistence, mortality, recruitment), propagates 
#               dynamic measurement error to calculate 95% Minimum Detectable 
#               Change (MDC95), and exports Master_CWC_Tracked_MDC.xlsx for 
#               downstream plotting and statistical analyses (Script 03).
# Dependencies: readxl, readr, dplyr, writexl
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(readxl)
library(readr)
library(dplyr)
library(writexl)

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & FILE PATHS
# ------------------------------------------------------------------------------
# Base path (adjust to match your directory structure if necessary)
base_dir <- "D:/PhD_Data(Large)/06_Final_Analysis/Master_Dataset"

cnn_error_path    <- file.path(base_dir, "Error_Propagation/CNN Error.xlsx")
input_csv_path    <- file.path(base_dir, "master_CWC_tracked.csv")
output_xlsx_path  <- file.path(base_dir, "Error_Propagation/Master_CWC_Tracked_MDC.xlsx")
summary_xlsx_path <- file.path(base_dir, "Error_Propagation/Master_CWC_Species_Summary.xlsx")

# Sampling duration (2015 to 2022)
time_interval_years <- 7

# ------------------------------------------------------------------------------
# 3. MODEL EMPIRICAL MEASUREMENT ERROR (TAGLAB REPEAT SEGMENTATION)
# ------------------------------------------------------------------------------
cat("Loading CNN error calibration data...\n")
cnn_df <- read_excel(cnn_error_path)

# Standardize lower-case column names from raw repeat segmentation file
colony_stats <- cnn_df %>%
  group_by(colony_id, species, timepoint) %>%
  summarise(
    mean_area = mean(area, na.rm = TRUE),
    sd_area   = sd(area, na.rm = TRUE),
    cv_area   = sd_area / mean_area,
    .groups   = "drop"
  )

# Fit linear empirical error model: SD(area) ~ Mean(area) + Species
model_final <- lm(sd_area ~ mean_area + species, data = colony_stats)

cat("Empirical Error Model Summary:\n")
print(summary(model_final))

# ------------------------------------------------------------------------------
# 4. LOAD TRACKED DATASET & CLASSIFY STATE TRANSITIONS
# ------------------------------------------------------------------------------
cat("\nProcessing tracked master dataset from Script 01...\n")
df <- read_csv(input_csv_path, show_col_types = FALSE)

# Map initial/final area and establish presence/absence rules
df <- df %>%
  mutate(
    A1 = Area_2015,
    A2 = Area_2022,
    
    present_2015 = !is.na(A1) & A1 > 0,
    present_2022 = !is.na(A2) & A2 > 0,
    
    # Classify state transitions (exact naming expected by Script 03)
    state_transition = case_when(
      present_2015 & present_2022  ~ "persistence",
      present_2015 & !present_2022 ~ "mortality",
      !present_2015 & present_2022 ~ "recruitment",
      TRUE                         ~ "unknown"
    )
  )

# ------------------------------------------------------------------------------
# 5. DYNAMIC ERROR PROPAGATION & DETECTABLE CHANGE (MDC95)
# ------------------------------------------------------------------------------
# Prepare prediction data frames matching model terms (mean_area & species)
pred_df_2015 <- data.frame(mean_area = df$A1, species = df$Species)
pred_df_2022 <- data.frame(mean_area = df$A2, species = df$Species)

df <- df %>%
  mutate(
    # Predict dynamic measurement uncertainty (sigma) per colony based on area & species
    sigma1 = predict(model_final, newdata = pred_df_2015),
    sigma2 = predict(model_final, newdata = pred_df_2022),
    
    # Total change in planar area over 7-year monitoring window (cm²)
    change_area = A2 - A1,
    
    # Root Sum of Squares (quadrature sum) for error propagation
    sigma_change = sqrt(sigma1^2 + sigma2^2),
    
    # Absolute 95% Minimum Detectable Change threshold (cm²)
    MDC95 = 1.96 * sigma_change,
    
    # Binary detectability boolean (evaluated for persistent colonies)
    detectable_change = case_when(
      state_transition == "persistence" ~ abs(change_area) > MDC95,
      TRUE                              ~ NA
    ),
    
    # ANNUALIZED METRIC ALIGNMENT FOR SCRIPT 03 ---------------------------------
    Annual_Area_Change = change_area / time_interval_years,
    MDC95_Annual       = MDC95 / time_interval_years
  )

# ------------------------------------------------------------------------------
# 6. SPECIES-LEVEL ECOLOGICAL SUMMARY
# ------------------------------------------------------------------------------
species_summary <- df %>%
  group_by(Species) %>%
  summarise(
    n                  = n(),
    mean_change        = mean(change_area[state_transition == "persistence"], na.rm = TRUE),
    mean_MDC           = mean(MDC95[state_transition == "persistence"], na.rm = TRUE),
    percent_detectable = mean(detectable_change == TRUE, na.rm = TRUE) * 100,
    mortality_rate     = sum(state_transition == "mortality") / n,
    recruitment_rate   = sum(state_transition == "recruitment") / n,
    .groups            = "drop"
  )

cat("\nSpecies Demographic & Error Summary:\n")
print(species_summary)

# ------------------------------------------------------------------------------
# 7. EXPORT DATASETS
# ------------------------------------------------------------------------------
# Ensure directory exists before writing
dir.create(dirname(output_xlsx_path), recursive = TRUE, showWarnings = FALSE)

write_xlsx(df, output_xlsx_path)
write_xlsx(species_summary, summary_xlsx_path)

cat("\n======================================================================\n")
cat("Script 02 execution complete!\n")
cat("Master MDC file saved to:\n ", output_xlsx_path, "\n")
cat("Species Summary saved to:\n ", summary_xlsx_path, "\n")
cat("======================================================================\n")
