# ==============================================================================
# Script Name:  04_CWC_Unfiltered_Growth_Vital_Rate_Diagnostics.R
# Description:  Models growth vital rates across ALL persistent colonies (including
#                growth, stability, and shrinkage) using non-linear GNLS with 
#                asymptotic curves (SSasymp). Implements Zuur's protocol for 
#                variance structure selection, generates normalized residual
#                diagnostics, and exports fitted model object for Script 05.
# Dependencies: readxl, dplyr, nlme, bbmle, openxlsx, ggplot2, patchwork, ragg
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(nlme)
library(bbmle)
library(openxlsx)
library(ggplot2)
library(patchwork)
library(ragg)

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & FILE PATHS
# ------------------------------------------------------------------------------
base_dir <- "D:/PhD_Data(Large)/Submission_Dataset"

input_xlsx_path <- file.path(base_dir, "Master_CWC_Tracked_MDC.xlsx")
output_tab_path <- file.path(base_dir, "tables", "Script04_Formal_Vital_Rate_Models.xlsx")
output_fig_path <- file.path(base_dir, "Figures", "Script04_Vital_Rate_Diagnostics.png")
output_rds_path <- file.path(base_dir, "models", "m_growth_gnls_SSasymp.rds")

FONT_FAMILY  <- "sans"
group_levels <- c("Madrepora oculata", "Desmophyllum pertusum", "Primnoa msp.")

species_palette <- c(
  "Madrepora oculata"     = "#722082",
  "Desmophyllum pertusum" = "#B63679",
  "Primnoa msp."          = "#fb9f3a"
)

# ------------------------------------------------------------------------------
# 3. LOAD DATASET & PREPARE UNFILTERED PERSISTENCE DATASET
# ------------------------------------------------------------------------------
cat("Loading master dataset...\n")
df_raw <- read_excel(input_xlsx_path)

# Standardize species factor levels
df_clean <- df_raw %>%
  filter(Species != "Coral Recruit") %>%
  mutate(
    Group = case_when(
      Species == "Madrepora oculata" ~ "Madrepora oculata",
      Species %in% c("D. pertusum", "Desmophyllum pertusum") ~ "Desmophyllum pertusum",
      Species %in% c("Primnoa msp.5", "Primnoa msp.1") ~ "Primnoa msp.",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%
  mutate(fSpecies = factor(Group, levels = group_levels))

# UNFILTERED GROWTH DATASET: Includes ALL surviving colonies (Positive, Zero, and Negative RGR)
df_growth <- df_clean %>%
  filter(state_transition == "persistence", !is.na(RGR_Planar), !is.na(A1))

cat("Total persistent colonies included in growth model:", nrow(df_growth), "\n")

# ------------------------------------------------------------------------------
# 4. ZUUR PROTOCOL: VARIANCE STRUCTURE SELECTION (ML)
# ------------------------------------------------------------------------------
cat("\nEvaluating variance structures for continuous Growth (Zuur Protocol)...\n")

# Fit baseline GLS model (Maximum Likelihood for variance structure comparison)
m_0        <- gls(RGR_Planar ~ A1 * fSpecies, data = df_growth, method = "ML")

# Candidate Variance Structures
m_ident    <- update(m_0, weights = varIdent(form = ~ 1 | fSpecies))
m_power    <- update(m_0, weights = varPower(form = ~ A1))
m_power_sp <- update(m_0, weights = varPower(form = ~ A1 | fSpecies))
m_exp      <- update(m_0, weights = varExp(form = ~ A1))
m_exp_sp   <- update(m_0, weights = varExp(form = ~ A1 | fSpecies))
m_cp       <- update(m_0, weights = varConstPower(form = ~ A1))
m_cp_sp    <- update(m_0, weights = varConstPower(form = ~ A1 | fSpecies))

# AIC Model Comparison Matrix
model_comparison <- AIC(m_0, m_ident, m_power, m_power_sp, 
                        m_exp, m_exp_sp, m_cp, m_cp_sp)

rownames(model_comparison) <- c("Baseline (Homoscedastic)", "VarIdent", "VarPower_Global", 
                                "VarPower_Sp", "VarExp_Global", "VarExp_Sp", 
                                "VarConstPower_Global", "VarConstPower_Sp")

model_comparison <- model_comparison %>%
  tibble::rownames_to_column("Model") %>%
  mutate(dAIC = AIC - min(AIC)) %>%
  arrange(AIC)

print(model_comparison)

# ------------------------------------------------------------------------------
# 5. ASYMPTOTIC MODEL FITTING (WINNING MODEL: VarConstPower_Global)
# ------------------------------------------------------------------------------
cat("\nFitting winning non-linear Asymptotic (SSasymp) GNLS model using VarConstPower_Global...\n")

# Shared starting values via nlsList
m_list <- nlsList(
  RGR_Planar ~ SSasymp(A1, Asym, R0, lrc) | fSpecies, 
  data = df_growth, 
  na.action = na.omit
)

start_vec <- c(coef(m_list)$Asym, coef(m_list)$R0, coef(m_list)$lrc)

# Fit Final GNLS Model
m_growth_gnls <- gnls(
  RGR_Planar ~ SSasymp(A1, Asym, R0, lrc),
  params  = list(Asym ~ fSpecies, R0 ~ fSpecies, lrc ~ fSpecies),
  data    = df_growth,
  weights = varConstPower(form = ~ A1),
  start   = start_vec,
  control = gnlsControl(maxIter = 500, minScale = 1e-7, tolerance = 1e-4),
  na.action = na.omit
)

cat("\n=== GNLS Model Summary (VarConstPower_Global) ===\n")
print(summary(m_growth_gnls))

# ------------------------------------------------------------------------------
# 6. SAVE FITTED MODEL OBJECT (FOR SCRIPT 05)
# ------------------------------------------------------------------------------
cat("\nSaving fitted GNLS model object to RDS for Script 05...\n")
dir.create(dirname(output_rds_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(m_growth_gnls, file = output_rds_path)
cat(" Model object saved to:\n  ", output_rds_path, "\n")

# ------------------------------------------------------------------------------
# 7. MODEL DIAGNOSTICS: BASELINE (m_0) VS WINNING MODEL
# ------------------------------------------------------------------------------
cat("\nExtracting residuals for Baseline vs. Winning Model comparison...\n")

# Combine residuals into a single plotting data frame
df_diag_m0 <- df_growth %>%
  mutate(
    Resid = residuals(m_0, type = "pearson"),
    Model = "1. Baseline Model (Homoscedastic m_0)"
  )

df_diag_win <- df_growth %>%
  mutate(
    Resid = residuals(m_growth_gnls, type = "normalized"),
    Model = "2. Winning Model (VarConstPower_Global)"
  )

df_diag_all <- bind_rows(df_diag_m0, df_diag_win) %>%
  mutate(Model = factor(Model, levels = c(
    "1. Baseline Model (Homoscedastic m_0)", 
    "2. Winning Model (VarConstPower_Global)"
  )))

base_diag_theme <- theme_classic(base_family = FONT_FAMILY) +
  theme(
    plot.title   = element_text(face = "bold", size = 11),
    axis.text    = element_text(size = 9, color = "black"),
    axis.title   = element_text(size = 10, face = "bold"),
    strip.text   = element_text(face = "bold", size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# Panel A: Heteroscedasticity Comparison (2 columns)
p_hetero_comp <- ggplot(df_diag_all, aes(x = A1, y = Resid, color = fSpecies)) +
  geom_point(alpha = 0.55, size = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "red", linewidth = 0.5) +
  scale_x_log10() +
  scale_color_manual(values = species_palette, name = "Species") +
  facet_wrap(~ Model, ncol = 2) +
  labs(
    title = "a) Variance Stabilization: Baseline (m_0) vs. Winning Model",
    x = expression(paste("Initial Planar Area A1 (cm"^2*", log scale)")),
    y = "Residuals (Pearson / Normalized)"
  ) +
  base_diag_theme

# Diagnostic Panel B: Q-Q Plots for Winning Model
p_qq_win <- ggplot(df_diag_win, aes(sample = Resid, color = fSpecies)) +
  stat_qq(alpha = 0.55, size = 1.6) +
  stat_qq_line(color = "black", linetype = "dashed") +
  facet_wrap(~ fSpecies, scales = "free_y") +
  scale_color_manual(values = species_palette, name = "Species") + # Matches Panel A name
  labs(
    title = "b) Q-Q Normal Plots of Normalized Residuals (Winning Model)",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  base_diag_theme +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold.italic")
  )

# Combine into composite figure
diag_composite <- (p_hetero_comp) / (p_qq_win) + 
  plot_layout(guides = "collect", heights = c(1, 1)) & 
  theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# 8. EXPORT SUMMARY TABLES & HIGH-RES DIAGNOSTIC GRAPHIC
# ------------------------------------------------------------------------------
dir.create(dirname(output_tab_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_fig_path), recursive = TRUE, showWarnings = FALSE)

cat("\nSaving model comparison table and GNLS coefficients to Excel...\n")
wb <- createWorkbook()
addWorksheet(wb, "Growth_Variance_Selection")
addWorksheet(wb, "GNLS_Coefficients")

writeData(wb, "Growth_Variance_Selection", model_comparison)
writeData(wb, "GNLS_Coefficients", as.data.frame(summary(m_growth_gnls)$tTable) %>% tibble::rownames_to_column("Parameter"))

saveWorkbook(wb, output_tab_path, overwrite = TRUE)

cat("Exporting high-resolution diagnostic composite image...\n")
ggsave(
  filename = output_fig_path,
  plot     = diag_composite,
  width    = 210,
  height   = 180,
  units    = "mm",
  dpi      = 600,
  device   = ragg::agg_png,
  bg       = "white"
)

cat("\n======================================================================\n")
cat("Script 04 Execution Complete!\n")
cat("Excel table saved to:  ", output_tab_path, "\n")
cat("Diagnostic figure saved to:", output_fig_path, "\n")
cat("Fitted RDS model saved to: ", output_rds_path, "\n")
cat("======================================================================\n")
